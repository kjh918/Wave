# flexible_workflow.py
from __future__ import annotations
import re
import subprocess as sp
from pathlib import Path
from typing import Dict, Any, List, Optional
import yaml

from tasks.task import TaskRegistry  # Task 구현은 외부에서 import되어 등록돼 있다고 가정


class Workflow:
    def __init__(self, config_path: Path):
        self.config_path = Path(config_path)
        self.cfg = yaml.safe_load(Path(config_path).read_text())

        wp = self.cfg.get("WORK_PARAMETERS", {})
        self.input_tpl: str = wp["input_path"]  # 예: /.../[sample_id]_R[read_info].fastq.gz
        self.work_dir = Path(wp["work_dir_path"]).resolve()
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # TASK_LIST: list[ {name,type,params}, ... ] 또는 dict[name->{type,params}]
        self.tasks_spec = self._normalize_tasks(self.cfg.get("TASK_LIST", []))

    # --------------------------
    # 유틸: 템플릿 파서
    # --------------------------
    def _parse_placeholders(self, tpl: str) -> List[str]:
        """[name] 패턴 추출"""
        return re.findall(r"\[([a-zA-Z0-9_]+)\]", tpl)

    def _tpl_to_regex(self, tpl: str, placeholders: List[str]) -> re.Pattern:
        """플레이스홀더를 정규식 그룹으로 변환"""
        esc = re.escape(tpl)
        for ph in placeholders:
            esc = esc.replace(re.escape(f"[{ph}]"), rf"(?P<{ph}>[^/]+)")
        return re.compile(rf"^{esc}$")

    def _tpl_to_glob(self, tpl: str, placeholders: List[str]) -> str:
        pat = tpl
        for ph in placeholders:
            pat = pat.replace(f"[{ph}]", "*")
        return pat

    # --------------------------
    # 샘플 탐색
    # --------------------------
    def discover_samples(self) -> Dict[str, Dict[str, Any]]:
        tpl = self.input_tpl
        placeholders = self._parse_placeholders(tpl)
        if "sample_id" not in placeholders:
            raise ValueError("input_path template must include [sample_id]")

        rx = self._tpl_to_regex(tpl, placeholders)
        glob_pat = self._tpl_to_glob(tpl, placeholders)

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
            samples.setdefault(sid, {"_files": []})
            samples[sid]["_files"].append(str(f))
            # 나머지 placeholder는 key별로 계층화
            for k, v in g.items():
                samples[sid].setdefault(k, {})
                samples[sid][k][v] = str(f)
        return samples

    # --------------------------
    # TASK_LIST 정규화/검증
    # --------------------------
    def _normalize_tasks(self, spec_raw: Any) -> List[Dict[str, Any]]:
        if isinstance(spec_raw, dict):
            items = []
            for name, spec in spec_raw.items():
                it = {"name": name, **spec}
                items.append(it)
            # dict일 경우 순서는 보장 안 될 수 있으니 필요시 name 기준 정렬 옵션 달 수 있음
            return items
        elif isinstance(spec_raw, list):
            # 각 원소: {name,type,params}
            return spec_raw
        else:
            raise TypeError("TASK_LIST must be list or dict")

    # --------------------------
    # set_up_tasks: Task 클래스 로드 & 파라미터 템플릿 준비
    # --------------------------
    def set_up_tasks(self) -> List[Dict[str, Any]]:
        tasks: List[Dict[str, Any]] = []
        for idx, spec in enumerate(self.tasks_spec, 1):
            if "type" not in spec or "name" not in spec:
                raise ValueError(f"TASK_LIST item missing 'name' or 'type': {spec}")
            # TaskRegistry에 타입이 등록되어 있어야 함
            _ = TaskRegistry.get(spec["type"])  # 미존재 시 KeyError 발생
            tasks.append({
                "order": idx,
                "name": spec["name"],
                "type": spec["type"],
                "params": spec.get("params", {}),
            })
        return tasks

    # --------------------------
    # 내부: 문자열 템플릿 전개
    # --------------------------
    def _fmt(self, val: Any, ctx: Dict[str, Any]) -> Any:
        if isinstance(val, str):
            try:
                return val.format(**ctx)
            except Exception:
                return val
        if isinstance(val, dict):
            return {k: self._fmt(v, ctx) for k, v in val.items()}
        if isinstance(val, list):
            return [self._fmt(v, ctx) for v in val]
        return val

    # --------------------------
    # build_workflow: 각 샘플 × 태스크에 대해 스크립트 작성
    # --------------------------
    def build_workflow(self, tasks: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
        tasks = tasks or self.set_up_tasks()
        samples = self.discover_samples()
        if not samples:
            return {"samples": {}, "tasks": tasks}

        master_scripts: Dict[str, Path] = {}
        for sid, info in samples.items():
            ctx = {"sample_id": sid, **{k: v for k, v in info.items() if not k.startswith("_")}}
            sid_dir = self.work_dir / sid
            sid_dir.mkdir(parents=True, exist_ok=True)

            # 샘플용 마스터 스크립트
            master_path = sid_dir / f"workflow_{sid}.sh"
            with master_path.open("w") as mf:
                mf.write("#!/usr/bin/env bash
set -euo pipefail
")
                for t in tasks:
                    TaskCls = TaskRegistry.get(t["type"])
                    # 파라미터 템플릿 전개 (예: "{sample_id}")
                    params = self._fmt(t.get("params", {}), ctx)
                    # 태스크 전용 작업 디렉토리
                    tdir = sid_dir / f"{t['order']:02d}_{t['name']}"
                    tdir.mkdir(parents=True, exist_ok=True)
                    # Task 인스턴스화
                    task = TaskCls(name=t["name"], params=params, workdir=tdir)
                    # 쉘 라인 생성
                    lines = list(task.to_sh())
                    # 태스크 스크립트 작성
                    t_sh = tdir / "task.sh"
                    t_sh.write_text("#!/usr/bin/env bash
set -euo pipefail
" + "
".join(lines) + "
")
                    t_sh.chmod(0o755)
                    # 마스터에 실행라인 추가
                    mf.write(f"bash {t_sh.as_posix()}
")
            master_path.chmod(0o755)
            master_scripts[sid] = master_path

        return {"samples": samples, "tasks": tasks, "masters": master_scripts}

    # --------------------------
    # run: 빌드된 마스터 스크립트 실행
    # --------------------------
    def run(self, run: bool = False):
        plan = self.build_workflow()
        if not run:
            print("[WAVE] Dry-run (scripts generated only)")
            for sid, mp in plan.get("masters", {}).items():
                print(f"  - {sid}: {mp}")
            return plan

        masters = plan.get("masters", {})
        for sid, mp in masters.items():
            print(f"[WAVE] Running {sid}: {mp}")
            ret = sp.run(["bash", str(mp)]).returncode
            if ret != 0:
                raise RuntimeError(f"Workflow for {sid} failed with code {ret}")
        return plan


# --------------------------
# 실행 예시 (선택)
# --------------------------
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser("wave-flex")
    ap.add_argument("--config", required=True)
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()

    wf = Workflow(Path(args.config))
    wf.run(run=args.run)