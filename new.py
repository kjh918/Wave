# workflow.py
from __future__ import annotations
import re, yaml, shlex, subprocess as sp, importlib, pkgutil, json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Any, List, Optional, Iterable

from src.tasks.task import TaskRegistry  # Task 클래스는 각 모듈에서 import됨


# --------------------------
# 유틸: 플레이스홀더 처리 & 쉘 join
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

def _sh_join(argv: Iterable[str]) -> str:
    return shlex.join(list(map(str, argv)))

def _render_with_ctx(val: Any, ctx: Dict[str, Any]) -> Any:
    if isinstance(val, str):
        return re.sub(r"\{([a-zA-Z0-9_]+)\}", lambda m: str(ctx.get(m.group(1), "")), val)
    if isinstance(val, list):
        return [_render_with_ctx(v, ctx) for v in val]
    if isinstance(val, dict):
        return {k: _render_with_ctx(v, ctx) for k, v in val.items()}
    return val

def _autoload_tasks(package_root: str = "src.tasks") -> None:
    """
    src.tasks 하위 모든 모듈을 import하여, 각 모듈 내의
    @TaskRegistry.register 가 실행되도록 한다.
    """
    try:
        pkg = importlib.import_module(package_root)
    except Exception as e:
        print(f"[WAVE] cannot import {package_root}: {e}")
        return

    for m in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
        # 언더스코어 시작 서브모듈은 스킵(내부/보조)
        if any(part.startswith("_") for part in m.name.split(".")):
            continue
        try:
            importlib.import_module(m.name)
        except Exception as e:
            # 개별 태스크 모듈 문제는 워크플로 전체를 막지 않기
            print(f"[WAVE] task autoload skip {m.name}: {e}")


# --------------------------
# 워크플로우
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
    # 입력/출력/파라미터 렌더
    # --------------------------
    def _render_block(self, block: Dict[str, Any], sample_id: str, work_dir_for_task: Path) -> Dict[str, Any]:
        ctx = {
            "sample_id": sample_id,
            "work_dir_path": str(self.work_dir),
            "WORK_DIR": str(work_dir_for_task),
        }
        return _render_with_ctx(block or {}, ctx)

    # --------------------------
    # TASK_LIST 정규화 (요청 포맷 지원)
    # --------------------------
    def _normalize_tasklist(self, task_list_raw: Any, sample_id: str) -> List[Dict[str, Any]]:
        """
        TASK_LIST:
          - <TaskName>:
              TOOL: <tool>       # 예: gatk, gatk3, bwa, samtools ...
              FUNC: <func>       # 예: markduplicates, baserecalibrator, index ...
              WORK_DIR: "<Path with {work_dir_path}/{sample_id}>"
              THREADs: <int>
              INPUT:  {...}
              OUTPUT: {...}
              PARAMS: {...}
        → [{name, tool, func, type, workdir(Path), threads, inputs, outputs, params}]
        """
        if not isinstance(task_list_raw, list):
            raise TypeError("TASK_LIST must be a list")

        norm: List[Dict[str, Any]] = []
        for item in task_list_raw:
            if not isinstance(item, dict) or len(item) != 1:
                raise ValueError(f"Invalid TASK_LIST entry: {item}")
            name, spec = next(iter(item.items()))
            if not isinstance(spec, dict):
                raise ValueError(f"Task spec must be a dict: {item}")

            tool = str(spec.get("TOOL") or "").strip()
            func = str(spec.get("FUNC") or "").strip()
            if not tool or not func:
                raise ValueError(f"TASK '{name}' must have TOOL and FUNC")

            # WORK_DIR 렌더
            raw_workdir = spec.get("WORK_DIR") or "{work_dir_path}/{sample_id}/" + name
            rendered_workdir = _render_with_ctx(raw_workdir, {
                "sample_id": sample_id,
                "work_dir_path": str(self.work_dir),
                "WORK_DIR": "",  # 아직 미정 (이 값은 여기서는 필요 없음)
            })
            workdir_path = Path(rendered_workdir).resolve()

            # 블록 렌더
            inputs  = self._render_block(spec.get("INPUT", {}),  sample_id, workdir_path)
            outputs = self._render_block(spec.get("OUTPUT", {}), sample_id, workdir_path)
            params  = self._render_block(spec.get("PARAMS", {}), sample_id, workdir_path)

            # THREADs 정규화
            threads = spec.get("THREADs") or spec.get("threads") or params.get("threads") or 1
            try:
                threads = int(threads)
            except Exception:
                raise ValueError(f"TASK '{name}' THREADs must be an integer (got: {threads})")

            norm.append({
                "name": name,
                "tool": tool,
                "func": func,
                "type": f"{tool}.{func}",
                "workdir": workdir_path,
                "threads": threads,
                "inputs": inputs,
                "outputs": outputs,
                "params": params,
            })
        return norm

    # --------------------------
    # Task 클래스 해석
    # --------------------------
    def _resolve_task_class(self, type_name: str):
        """
        type_name 예: "gatk.markduplicates", "samtools.sort"
        우선 TaskRegistry, 실패 시 동적 import 한 뒤 재시도
        """
        try:
            return TaskRegistry.get(type_name)
        except Exception:
            # 컨벤션: src.tasks.<tool>.<func>.main 모듈 임포트 시도
            try:
                importlib.import_module(f"src.tasks.{type_name}.main")
            except ModuleNotFoundError:
                # 혹시 모듈이 다른 파일명인 경우도 있을 수 있어 보정
                parts = type_name.split(".")
                if len(parts) == 2:
                    try:
                        importlib.import_module(f"src.tasks.{parts[0]}.{parts[1]}")
                    except ModuleNotFoundError:
                        pass
            # 재시도
            return TaskRegistry.get(type_name)

    # --------------------------
    # 쉘 라인 보정
    # --------------------------
    def _to_shell_lines(self, lines: Iterable[str]) -> Iterable[str]:
        for line in lines:
            if isinstance(line, (list, tuple)):
                yield _sh_join(line)
            else:
                yield str(line)

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
            tasks_norm = self._normalize_tasklist(self.cfg.get("TASK_LIST", []), sample_id=sid)

            # 샘플 작업 루트
            sid_root = self.work_dir / sid
            sid_root.mkdir(parents=True, exist_ok=True)

            # JSON 요약(선택)
            summary_json = sid_root / f"workflow_{sid}.json"
            summary_json.write_text(json.dumps(tasks_norm, indent=2, ensure_ascii=False))

            # 마스터 스크립트
            master_sh = sid_root / f"workflow_{sid}.sh"
            with master_sh.open("w") as mf:
                mf.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                for t in tasks_norm:
                    tdir: Path = t["workdir"]
                    tdir.mkdir(parents=True, exist_ok=True)

                    TaskCls = self._resolve_task_class(t["type"])
                    # ✅ 새로운 Task 시그니처에 맞춰 생성
                    task = TaskCls(
                        name=t["name"],
                        tool=t["tool"],
                        func=t["func"],
                        threads=t["threads"],
                        workdir=tdir,
                        inputs=t.get("inputs", {}),
                        outputs=t.get("outputs", {}),
                        params=t.get("params", {}),
                    )

                    # 쉘 라인 생성
                    lines = list(self._to_shell_lines(task.to_sh()))

                    # 태스크 스크립트 저장
                    t_sh = tdir / f"{sid}_{t['name']}.sh"
                    t_sh.write_text(
                        "#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(lines) + "\n"
                    )
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
            ret = sp.run(["bash", str(sh)]).returncode
            if ret != 0:
                raise RuntimeError(f"Workflow for {sid} failed with code {ret}")
        return plan


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