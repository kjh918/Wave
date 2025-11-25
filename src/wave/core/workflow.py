
# flexible_workflow.py
from __future__ import annotations
import os
import re
import yaml
import shlex
import subprocess
import importlib
import pkgutil
import sys
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Any, List, Optional, Iterable

# allow local imports
sys.path.append(os.path.dirname(__file__))

from executor import SunGridExecutor, BashExecutor
from wave.core.task import Task
from wave.core.task_registry import TaskRegistry
from wave.utils.workflow_utils import (
    autoload_tasks,
    skip_finished_tasks,
    parse_placeholders,
    tpl_to_regex,
    render_value,
    tpl_to_glob,
    sh_join,
    instantiate_task,
)
from wave.utils.log import Logger


@dataclass
class Workflow:
    config_path: Optional[Path] = None
    config_dict: Optional[Dict[str, Any]] = None

    cfg: Dict[str, Any] = field(init=False)
    work_dir: Path = field(init=False)
    input_tpl: str = field(init=False)
    workflow: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        # load config
        if self.config_dict is not None:
            self.cfg = self.config_dict
        elif self.config_path is not None:
            self.cfg = yaml.safe_load(Path(self.config_path).read_text())
        else:
            raise ValueError("Workflow: either config_path or config_dict must be provided.")

        # basic params
        self.w_params = self.cfg["WORK_PARAMETERS"]
        self.input_tpl = self.w_params["input_path"]
        self.work_dir = Path(self.w_params["work_dir_path"]).resolve()
        self.work_dir.mkdir(parents=True, exist_ok=True)

        self.workflow = {
            "SETTING": self.cfg.get("SETTING", {}),
            "WORK_PARAMETERS": self.w_params,
            "SAMPLES": {},
        }

        self.max_threads_total = int(self.cfg["SETTING"]['MAX_THREADS_TOTAL'])
        self.max_samples= int(self.cfg["SETTING"]['MAX_SAMPLES'])
        
        autoload_tasks("wave.tasks")

        self.workflow_logdir = Path(self.w_params["work_dir_path"]) / 'log'
        # Logger(self.workflow_logdir)

    ## util parer 에서 fastq 파일 parsing하도록 ## 
    # --------------------------
    # sample discovery (template)
    # --------------------------
    def discover_samples(self) -> Dict[str, Dict[str, Any]]:
        tpl = self.input_tpl
        placeholders = parse_placeholders(tpl)
        if "sample_id" not in placeholders:
            raise ValueError("input_path template must include {sample_id}")

        rx = tpl_to_regex(tpl, placeholders)
        glob_pat = tpl_to_glob(tpl, placeholders)

        files = [Path(p) for p in Path("/").glob(glob_pat.lstrip("/"))]
        if not files:
            print(f"[WAVE] No files matched {glob_pat}")
            return {}

        samples: Dict[str, Dict[str, Any]] = {}
        for f in files:
            m = rx.match(str(f))
            if not m:
                continue
            g = m.groupdict()
            sid = g.pop("sample_id")
            s = samples.setdefault(sid, {"inputs": {}})
            s_inputs = s["inputs"]
            s_inputs.setdefault("raw_dir", str(f.parent))

            if not g:
                s_inputs["single"] = str(f)
                continue

            for key, val in g.items():
                bucket = s_inputs.setdefault(key, {})
                bucket[str(val)] = str(f)

        return samples

    # --------------------------
    # normalize legacy TASK_LIST
    # --------------------------
    def _set_params(self, sample_id, task_work_dir, inputs):
        setted_inputs = {}
        for key, value in inputs.items():
            value = str(value)
            for w_key, w_value in self.w_params.items():
                value = value.replace('{' + str(w_key) + '}', w_value)
            value = value.replace('{sample_id}', sample_id).replace('{WORK_DIR}', str(task_work_dir))
            setted_inputs[key] = value
        return setted_inputs

    ## 코드 정리 필요 ## 
    def _normalize_tasklist_legacy(self, task_list_raw: Any, sample_id: str) -> List[Dict[str, Any]]:
        if not isinstance(task_list_raw, list):
            raise TypeError("TASK_LIST must be a list for legacy form.")

        norm: List[Dict[str, Any]] = []
        for idx, item in enumerate(task_list_raw, 1):
            name, spec = next(iter(item.items()))
            ttype = spec.get("TOOL")
            func = spec.get("FUNC")
            threads = spec.get("THREADS") or 1

            if not isinstance(item, dict) or len(item) != 1:
                raise ValueError(f"Invalid TASK_LIST entry: {item}")
            if not isinstance(spec, dict):
                raise ValueError(f"Task spec must be a dict: {item}")
            if not ttype:
                raise ValueError(f"TASK '{name}' missing 'TOOL'")

            task_work_dir = str(spec.get("WORK_DIR", '') or self.work_dir)
            inp = spec.get("INPUT", {}) or {}
            outp = spec.get("OUTPUT", {}) or {}
            params = spec.get("PARAMS", {}) or {}
            task_work_dir = task_work_dir.replace('{work_dir_path}', str(self.work_dir)).replace('{sample_id}', sample_id)

            inputs = self._set_params(sample_id, task_work_dir, inp)
            outputs = self._set_params(sample_id, task_work_dir, outp)
            params = self._set_params(sample_id, task_work_dir, params)

            norm.append({
                "name": name,
                "tool": ttype,
                "func": func,
                "threads": threads,
                "workdir": Path(task_work_dir),
                "inputs": inputs,
                "outputs": outputs,
                "params": params,
            })

        return norm

    # --------------------------
    # build -> write per-sample master JSON
    # --------------------------
    def build(self) -> Dict[str, Any]:
        samples = self.discover_samples()
        self.workflow["SAMPLES"] = samples
        if not samples:
            return {"samples": {}, "masters": {}}

        masters: Dict[str, Path] = {}
        for sid in samples.keys():
            tasks_norm = self._normalize_tasklist_legacy(self.cfg.get("TASK_LIST", {}), sample_id=sid)

            sid_root = self.work_dir / sid
            sid_root.mkdir(parents=True, exist_ok=True)

            master_json = sid_root / f"workflow_{sid}.json"
            total_task_dict: Dict[str, Any] = {}

            for _task in tasks_norm:
                tdir: Path = _task['workdir']
                tdir.mkdir(parents=True, exist_ok=True)

                # resolve task class
                TaskCls = self._resolve_task_class(tool=_task["tool"], func=_task["func"])

                # instantiate task (compat wrapper)
                task: Task = instantiate_task(
                    TaskCls,
                    name=_task["name"],
                    tool=_task["tool"],
                    func=_task["func"],
                    threads=_task["threads"],
                    workdir=tdir,
                    inputs=_task.get("inputs", {}),
                    outputs=_task.get("outputs", {}),
                    params=_task.get("params", {}),
                )

                # obtain canonical shell line for first command
                cmd_lines = task.to_sh()
                if not cmd_lines:
                    raise RuntimeError(f"Task {task.name}.to_sh() returned empty list")
                # convert first line to shell string if needed
                first = cmd_lines[0]
                if isinstance(first, (list, tuple)):
                    cmd_str = sh_join(first)
                else:
                    cmd_str = str(first)

                total_task_dict[task.name] = {
                    'user': self.workflow['SETTING'].get('User'),
                    'node': self.workflow['SETTING'].get('Node'),
                    'job_id': f'{sid}_{task.name}',
                    'threads': int(task.threads or 1),
                    'inputs': _task.get("inputs", {}),
                    'outputs': _task.get("outputs", {}),
                    'workdir': str(tdir),
                    'cmd': cmd_str,
                }
                
            # write master json
            master_json.write_text(json.dumps(total_task_dict, indent=4))
            masters[sid] = master_json

        return {"samples": samples, "masters": masters}

    # --------------------------
    # run
    # --------------------------
    def run(self, executor: str = 'SGE', max_jobs: Optional[int] = None, max_threads_total: Optional[int] = None) -> Dict[str, Any]:
        plan = self.build()
        masters = plan.get("masters", {})
        if not masters:
            print("[WAVE] No masters to run.")
            return plan

        # state for throttling
        running_jobs = 0
        used_threads = 0

        for sid, json_file in masters.items():
            task_map: Dict[str, Any] = json.loads(Path(json_file).read_text())

            prev_qid: Optional[str] = None
            for task_name, meta in task_map.items():
                workdir = Path(meta["workdir"])
                workdir.mkdir(parents=True, exist_ok=True)

                done_flag = workdir / ".done"
                if done_flag.exists():
                    print(f"[WAVE] Skip {sid}:{task_name} (done flag found)")
                    continue

                # throttle by max_jobs / max_threads_total if provided
                if max_jobs:
                    while running_jobs >= max_jobs:
                        print("[WAVE] max_jobs reached, sleeping...")
                        import time
                        time.sleep(5)
                if max_threads_total and meta.get("threads"):
                    want = int(meta["threads"])
                    while used_threads + want > max_threads_total:
                        print("[WAVE] max_threads_total would be exceeded, sleeping...")
                        import time
                        time.sleep(5)

                if executor.upper() == 'SGE':
                    _executor = SunGridExecutor(logdir=workdir / "log")
                    script_path = _executor.make_script(cmd=meta["cmd"], job_id=meta["job_id"])
                    print(f"[WAVE] created script: {script_path}")
                    qid = _executor.qsub_sh(
                        node=meta.get("node"),
                        script_path=str(script_path),
                        threads=int(meta.get("threads", 1)),
                        job_id=meta.get("job_id"),
                        hold_jid=prev_qid,
                        finish_and_run=True,
                    )
                    print(f"[WAVE] qsub {sid}:{task_name} -> {qid}")
                    prev_qid = qid
                    used_threads += int(meta.get("threads", 1))

                else:
                    # local bash execution
                    _executor = BashExecutor(logdir=workdir / "log")
                    script_path = _executor.make_script(cmd=meta["cmd"], job_id=meta["job_id"])
                    print(f"[WAVE] local script: {script_path}")

                    _executor.run(
                        cmd = meta["cmd"],
                        job_id = f'{sid}_{task_name}',
                        workdir = meta["workdir"],
                        outputs = meta.get("outputs", {})
                    )
                    print(f"[WAVE] local run finished {sid}:{task_name}")
                    used_threads += int(meta.get("threads", 1))


        return plan

    ## 코드정리필요 ## 
    # --------------------------
    # helper: resolve task class
    # --------------------------
    def _resolve_task_class(self, tool: str, func: Optional[str] = None):
        key = tool if not func else f"{tool}.{func}"

        # 1) try registry
        try:
            return TaskRegistry.get(key)
        except KeyError:
            pass

        # 2) try import module (tasks follow convention wave.tasks.<tool>[.<func>].main)
        if func:
            module_path = f"wave.tasks.{tool}.{func}.main"
        else:
            module_path = f"wave.tasks.{tool}.main"

        try:
            importlib.import_module(module_path)
        except ModuleNotFoundError as e:
            raise KeyError(f"Task module not found (tried import): {module_path}") from e
        except Exception as e:
            # import error inside module -> re-raise with context
            raise RuntimeError(f"Error importing task module '{module_path}': {e}")

        # 3) try registry again
        try:
            return TaskRegistry.get(key)
        except KeyError as e:
            # provide helpful debug info
            registered = list(getattr(TaskRegistry, "_REG", []) if hasattr(TaskRegistry, "_REG") else [])
            raise KeyError(f"Unknown task TYPE: {key}. Registered: {registered}") from e


# --------------------------
# CLI
# --------------------------
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser("wave-flex")
    ap.add_argument("--config", required=True, help="config.yaml")
    ap.add_argument("--run", action="store_true", help="actually run (not dry-run)")
    ap.add_argument("--executor", default="SGE", help="SGE or BASH")
    ap.add_argument("--max-jobs", type=int, default=None, help="max concurrent jobs")
    ap.add_argument("--max-threads", type=int, default=None, help="max total threads across jobs")
    args = ap.parse_args()

    wf = Workflow(Path(args.config))
    plan = wf.build()
    print("[WAVE] Build complete. Masters:")
    for sid, p in plan.get("masters", {}).items():
        print(f"  - {sid}: {p}")
    if args.run:
        wf.run(executor=args.executor, max_jobs=args.max_jobs, max_threads_total=args.max_threads)