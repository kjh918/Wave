# flexible_workflow.py
from __future__ import annotations
import os, re, yaml, shlex, subprocess, importlib, pkgutil, sys, json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Any, List, Optional, Iterable

# 베이스 Task & 레지스트리
sys.path.append(os.path.dirname(__file__))

from src.executor import SunGridExecutor
from src.tasks.task import Task
from src.tasks.task import TaskRegistry  


# --------------------------
# 유틸: 태스크 오토로더 (전역 함수)
# --------------------------
def _autoload_tasks(package_root: str = "src.tasks") -> None:
    """
    src.tasks 하위 모든 서브모듈을 import해서
    @TaskRegistry.register 가 실행되도록 만든다.
    """
    pkg = importlib.import_module(package_root)
    for m in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
        # _ 로 시작하는 내부 모듈은 스킵
        if any(part.startswith("_") for part in m.name.split(".")):
            continue
        try:
            importlib.import_module(m.name)
        except Exception as e:
            print(f"[WAVE] task autoload skip {m.name}: {e}")


# --------------------------
# 유틸: 플레이스홀더 렌더/정규식/글롭
# --------------------------
def _parse_placeholders(tpl: str) -> List[str]:
    return re.findall(r"\{([a-zA-Z0-9_]+)\}", tpl)

def _tpl_to_regex(tpl: str, placeholders: List[str]) -> re.Pattern:
    esc = re.escape(tpl)
    for ph in placeholders:
        esc = esc.replace(re.escape(f"{{{ph}}}"), rf"(?P<{ph}>[^/]+)")
    return re.compile(rf"^{esc}$")

def _tpl_to_glob(tpl: str, placeholders: List[str]) -> str:
    pat = tpl
    for ph in placeholders:
        pat = pat.replace(f"{{{ph}}}", "*")
    return pat

def _render_value(val: Any, ctx: Dict[str, Any]) -> Any:
    if isinstance(val, str):
        return re.sub(r"\{([a-zA-Z0-9_]+)\}", lambda m: str(ctx.get(m.group(1), "")), val)
    if isinstance(val, list):
        return [_render_value(v, ctx) for v in val]
    if isinstance(val, dict):
        return {k: _render_value(v, ctx) for k, v in val.items()}
    return val

def _sh_join(argv: Iterable[str]) -> str:
    # 토큰 리스트를 안전하게 하나의 쉘 라인으로
    return shlex.join(list(map(str, argv)))


# --------------------------
# 핵심: Task 인스턴스 안전 생성 (여러 시그니처 지원)
# --------------------------
def _instantiate_task(TaskCls, *, name: str, workdir: Path,
                      inputs: Dict[str, Any], outputs: Dict[str, Any],
                      params: Dict[str, Any]) -> Task:
    """
    팀마다 Task 생성자 시그니처가 달라 발생하는 호환 이슈를 흡수.
    아래 순서로 시도:
      1) (name, workdir, inputs, outputs, params)
      2) (name, params, workdir)
      3) (name, workdir, params)
      4) (name=name, workdir=..., inputs=..., outputs=..., params=...)
    """
    # 1) 완전형
    try:
        return TaskCls(name=name, workdir=workdir, inputs=inputs, outputs=outputs, params=params)
    except TypeError:
        pass
    # 2) (name, params, workdir)
    try:
        return TaskCls(name, params, workdir)
    except TypeError:
        pass
    # 3) (name, workdir, params)
    try:
        return TaskCls(name, workdir, params)
    except TypeError:
        pass
    # 4) 키워드만
    try:
        return TaskCls(name=name, params=params, workdir=workdir)
    except TypeError as e:
        raise TypeError(f"Cannot instantiate task '{name}' with supported signatures: {e}")


# --------------------------
# 메인 Workflow
# --------------------------
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
        w_params = self.cfg["WORK_PARAMETERS"]
        self.input_tpl = w_params["input_path"]
        self.work_dir = Path(w_params["work_dir_path"]).resolve()
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # 3) 스켈레톤
        self.workflow = {
            "SETTING": self.cfg.get("SETTING", {}),
            "WORK_PARAMETERS": w_params,
            "SAMPLES": {},             # sample_id -> {...}
        }

        # 4) 태스크 모듈 자동 임포트 → 레지스트리 채우기
        _autoload_tasks("src.tasks")

    # --------------------------
    # 샘플 탐색 (input_path 템플릿 기반)
    # --------------------------
    def discover_samples(self) -> Dict[str, Dict[str, Any]]:
        tpl = self.input_tpl
        placeholders = _parse_placeholders(tpl)
        if "sample_id" not in placeholders:
            raise ValueError("input_path template must include {sample_id}")

        rx = _tpl_to_regex(tpl, placeholders)
        glob_pat = _tpl_to_glob(tpl, placeholders)

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
            setted_inputs[key] = value.replace('{work_dir_path}', str(self.work_dir)).replace('{sample_id}', sample_id).replace('{WORK_DIR}', str(task_work_dir))
        return setted_inputs
    
    def _normalize_tasklist_legacy(self, task_list_raw: Any, sample_id: str) -> List[Dict[str, Any]]:
        """
        TASK_LIST:
          - <name>:
              TOOL: <type>
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

            # task_dir_name = spec.get("WORK_DIR")
            # print(task_dir_name)
            # exit()
            ttype = spec.get("TOOL")

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
                "workdir": Path(task_work_dir),
                "type": ttype,
                "inputs": inputs,
                "outputs": outputs,
                "params": params,
            })
            
        return norm

    # --------------------------
    # 스크립트 빌드 (샘플 × 태스크) → task.sh / master.sh
    # --------------------------
    def build(self) -> Dict[str, Any]:
        samples = self.discover_samples()
        self.workflow["SAMPLES"] = samples

        if not samples:
            return {"samples": {}, "masters": {}}

        masters: Dict[str, Path] = {}
        for sid, _ in samples.items():
            # 1) 레거시 TASK_LIST 정규화
            
            tasks_norm = self._normalize_tasklist_legacy(self.cfg.get("TASK_LIST", {}), sample_id=sid)

            # 샘플 작업 루트
            sid_root = self.work_dir / sid
            sid_root.mkdir(parents=True, exist_ok=True)

            # 마스터 스크립트
            master_json = sid_root / f"workflow_{sid}.json"
            
            

            for _task in tasks_norm:
                tdir = _task['workdir']
                tdir.mkdir(parents=True, exist_ok=True)

                TaskCls = self._resolve_task_class(_task["type"])
                task = _instantiate_task(
                    TaskCls,
                    name = _task["name"],
                    workdir = tdir,
                    inputs = _task.get("inputs", {}),
                    outputs = _task.get("outputs", {}),
                    params = _task.get("params", {}),
                )

                task_cmd = list(self._to_shell_lines(task.to_sh()))

                executor = SunGridExecutor(
                    logdir = _task['workdir'] / 'qlog'
                    )
                # print(task_cmd)
                # exit()
                executor.run(
                    node=self.workflow['SETTING']['Node'], 
                    cmd=task_cmd[0], 
                    threads = _task.get("params", {})['threads'],
                    job_id = f'{sid}_{_task["name"]}'
                    )
                # print(tdir)
            # with open(master_json, "w") as json_file:
            #     json.dump(student_data, json_file)



            master_sh = sid_root / f"workflow_{sid}.sh"
            with master_sh.open("w") as mf:
                mf.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                for idx, t in enumerate(tasks_norm, 1):
                    tdir = t['workdir']
                    tdir.mkdir(parents=True, exist_ok=True)
                    # Task 인스턴스 생성
                    TaskCls = self._resolve_task_class(t["type"])
                    task = _instantiate_task(
                        TaskCls,
                        name=t["name"],
                        workdir=tdir,
                        inputs=t.get("inputs", {}),
                        outputs=t.get("outputs", {}),
                        params=t.get("params", {}),
                    )
                    # 커맨드 라인 생성
                    lines = list(self._to_shell_lines(task.to_sh()))

                    # task.sh 저장
                    t_sh = tdir / f"{sid}_{t['name']}.sh"
                    t_sh.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(lines) + "\n")
                    t_sh.chmod(0o755)

                    # 마스터에 실행 라인 추가
                    mf.write(f"bash {t_sh.as_posix()}\n")

            master_sh.chmod(0o755)
            masters[sid] = master_sh

        return {"samples": samples, "masters": masters}

    # --------------------------
    # 실행
    # --------------------------
    def run(self, run: bool = False):
        plan = self.build()
        masters = plan.get("masters", {})
        if not run:
            print("[WAVE] Dry-run: generated scripts")
            for sid, sh in masters.items():
                print(f"  - {sid}: {sh}")
            return plan

        for sid, sh in masters.items():
            print(f"[WAVE] Running {sid}: {sh}")
            ret = subprocess.run(["bash", str(sh)]).returncode
            if ret != 0:
                raise RuntimeError(f"Workflow for {sid} failed with code {ret}")
        return plan

    # --------------------------
    # 내부 헬퍼
    # --------------------------
    def _resolve_task_class(self, type_name: str):
        try:
            return TaskRegistry.get(type_name)
        except Exception:
            # 컨벤션: src.tasks.<type>.<type> 모듈 임포트 시도
            try:
                importlib.import_module(f"src.tasks.{type_name}.{type_name}")
            except ModuleNotFoundError:
                pass
            # 재시도 (여전히 없으면 KeyError 발생)
            return TaskRegistry.get(type_name)

    def _to_shell_lines(self, lines: Iterable[str]) -> Iterable[str]:
        """
        task.to_sh()가 이미 문자열 라인을 내면 그대로,
        만약 argv 토큰 리스트를 낸다면 shlex.join으로 변환.
        """
        for line in lines:
            if isinstance(line, (list, tuple)):
                yield _sh_join(line)
            else:
                yield str(line)


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