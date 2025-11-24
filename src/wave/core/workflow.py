# flexible_workflow.py
from __future__ import annotations
import os, re, yaml, shlex, subprocess, importlib, pkgutil, sys, json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Any, List, Optional, Iterable

# 베이스 Task & 레지스트리
sys.path.append(os.path.dirname(__file__))

from executor import SunGridExecutor
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


@dataclass
class Workflow:
    # 설정을 파일/딕셔너리로 모두 받을 수 있게
    config_path: Optional[Path] = None
    config_dict: Optional[Dict[str, Any]] = None

    # 내부 상태
    cfg: Dict[str, Any] = field(init=False)
    work_dir: Path = field(init=False)
    input_tpl: str = field(init=False)
    workflow: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        # 1) config 로드
        if self.config_dict is not None:
            self.cfg = self.config_dict
        elif self.config_path is not None:
            self.cfg = yaml.safe_load(Path(self.config_path).read_text())
        else:
            raise ValueError("Workflow: either config_path or config_dict must be provided.")

        # 2) 기본 설정
        self.w_params = self.cfg["WORK_PARAMETERS"]
        self.input_tpl = self.w_params["input_path"]
        self.work_dir = Path(self.w_params["work_dir_path"]).resolve()
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # 3) 스켈레톤
        self.workflow = {
            "SETTING": self.cfg.get("SETTING", {}),
            "WORK_PARAMETERS": self.w_params,
            "SAMPLES": {},             # sample_id -> {...}
        }

        # 동시에 제출할 전체 job 수 / 전체 스레드 수
        # self.max_jobs = int(setting.get("MaxJobs", 0) or 0)               # 0 → 제한 없음
        # self.max_threads_total = int(setting.get("MaxThreadsTotal", 0) or 0)


        # 4) 태스크 모듈 자동 임포트 → 레지스트리 채우기
        autoload_tasks("wave.tasks")

    # --------------------------
    # 샘플 탐색 (input_path 템플릿 기반)
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
    # 레거시 TASK_LIST 정규화 (요청 스키마 지원)
    # --------------------------

    def _set_params(self, sample_id, task_work_dir, inputs):
        setted_inputs = {}
        for key, value in inputs.items():
            value = str(value)

            for w_key, w_value in self.w_params.items():
                value = value.replace('{' + str(w_key) + '}',w_value).replace('{sample_id}', sample_id).replace('{WORK_DIR}', str(task_work_dir))
                setted_inputs[key] = value

        return setted_inputs
    
    def _normalize_tasklist_legacy(self, task_list_raw: Any, sample_id: str) -> List[Dict[str, Any]]:
        """
        TASK_LIST:
          - <name>:
              TOOL: <type>
              FUNC: <type>
              THREADS: <int>
              WORK_DIR: <Path>
              INPUT:  { read1: "...", read2: "..." }
              OUTPUT: { dirname: "..." }
              PARAMS: { ... }
        → [{name, type, inputs, outputs, params}]
        """
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
                raise ValueError(f"TASK '{name}' missing 'tool'")
            
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

    def build(self) -> Dict[str, Any]:
        samples = self.discover_samples()
        self.workflow["SAMPLES"] = samples

        if not samples:
            return {"samples": {}, "masters": {}}
        count = 0 
        masters: Dict[str, Path] = {}
        for sid, _ in samples.items():
            # 1) 레거시 TASK_LIST 정규화
            
            tasks_norm = self._normalize_tasklist_legacy(self.cfg.get("TASK_LIST", {}), sample_id=sid)

            # 샘플 작업 루트
            sid_root = self.work_dir / sid
            sid_root.mkdir(parents=True, exist_ok=True)

            # 마스터 스크립트
            master_json = sid_root / f"workflow_{sid}.json"
            total_task_dict = {}
            with open(master_json, 'w') as handle:
            
                qids: list[str] = []
                prev_qid: Optional[str] = None

                for _task in tasks_norm:
                    tdir = _task['workdir']
                    tdir.mkdir(parents=True, exist_ok=True)
                    TaskCls = self._resolve_task_class(tool=_task["tool"], func=_task["func"])
                    
                    task = instantiate_task(
                        TaskCls,
                        name = _task["name"],
                        tool = _task["tool"],
                        func = _task["func"],
                        threads = _task["threads"],
                        workdir = tdir,
                        inputs = _task.get("inputs", {}),
                        outputs = _task.get("outputs", {}),
                        params = _task.get("params", {}),
                    )
                    
                    total_task_dict[task.name] = {
                        'user': self.workflow['SETTING']['User'],
                        'node': self.workflow['SETTING']['Node'],
                        'job_id': f'{sid}_{task.name}',
                        'threads': task.threads,
                        'inputs' : _task.get("inputs", {}),
                        'outputs' : _task.get("outputs", {}),
                        'workdir': str(tdir),
                        'cmd': task.to_sh()[0]
                    }

                    workdir = Path(_task["workdir"])
                    workdir.mkdir(parents=True, exist_ok=True)

                    # .done 있으면 스킵(출력검증은 executor/flags에서 수행하는 게 베스트)
                    done_flag = workdir / ".done"

                    skip_finished_tasks(_task.get("outputs", {}), done_flag)
                    
                    sample_executor = SunGridExecutor(logdir=Path(tdir) / 'log')
                json.dump(total_task_dict, handle, indent=4)
            masters[sid] = master_json

        return {"samples": samples, "masters": masters}

    def run(self, executor='SGE') -> Dict[str, Any]:
        plan = self.build()
        masters = plan.get("masters", {})
        if not masters:
            print("[WAVE] No masters to run.")
            return plan

        for sid, json_file in masters.items():
            with open(json_file, "r") as handle:
                task_map: Dict[str, Any] = json.load(handle)

            prev_qid: Optional[str] = None
            for task_name, meta in task_map.items():

                workdir = Path(meta["workdir"])
                workdir.mkdir(parents=True, exist_ok=True)

                # .done 있으면 스킵(출력검증은 executor/flags에서 수행하는 게 베스트)
                done_flag = workdir / ".done"

                if done_flag.exists():
                    print(f"[WAVE] Skip {sid}:{task_name} (done flag found)")
                    continue

                executor = SunGridExecutor(logdir=workdir / "log")

                script_path = executor.make_script(
                    cmd=meta["cmd"], 
                    job_id=task_name,
                    )
                
                qid = executor.qsub_sh(
                    node=meta["node"],
                    script_path=str(script_path),
                    threads=int(meta["threads"]),
                    job_id=meta["job_id"],
                    hold_jid=prev_qid,
                    finish_and_run=True

                )
                # prev_qid = qid
                # print(f"[WAVE] qsub {sid}:{task_name} -> {qid}")

        return plan
    # --------------------------
    # 내부 헬퍼
    # --------------------------
    def _resolve_task_class(self, tool: str, func: Optional[str] = None):
        key = tool if not func else f"{tool}.{func}"  # fastqc / gatk4.baserecalibrator
        
        # 1) 먼저 레지스트리 조회
        try:
            return TaskRegistry.get(key)
        except KeyError:
            pass

        # 2) 모듈 임포트 시도
        if func:
            module_path = f"wave.tasks.{tool}.{func}.main"
        else:
            module_path = f"wave.tasks.{tool}.main"

        try:
            import importlib
            importlib.import_module(module_path)
        except ModuleNotFoundError as e:
            raise KeyError(f"Task module not found: {module_path}") from e

        # 3) 다시 레지스트리 조회
        return TaskRegistry.get(key)

# --------------------------
# CLI
# --------------------------
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser("wave-flex")
    ap.add_argument("--config", required=True)
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()

    wf = Workflow(Path(args.config))
    wf.run(run=args.run)