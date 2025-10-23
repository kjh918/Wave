# flexible_workflow.py
from __future__ import annotations
import re
import subprocess as sp
from pathlib import Path
from typing import Dict, Any, List, Optional
import yaml


class Workflow:
    def __init__(self, config_path: Path):
        self.config_path = Path(config_path)
        self.cfg = yaml.safe_load(Path(config_path).read_text())

        wp = self.cfg.get("WORK_PARAMETERS", {})
        self.input_tpl: str = wp["input_path"]  # 예: /.../[sample_id]_R[read_info].fastq.gz
        self.work_dir = Path(wp["work_dir_path"])

        self.tasks = self.cfg.get("TASK_LIST", {})

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
            if ph == "sample_id":
                esc = esc.replace(re.escape(f"[{ph}]"), r"(?P<sample_id>[^/]+)")
            else:
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

        # 절대경로 보정 후 스캔
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

            samples.setdefault(sid, {})
            # 남은 placeholder 들을 계층적으로 정리
            if len(g) == 1:
                # 단일 추가 변수
                key, val = next(iter(g.items()))
                samples[sid].setdefault(key, {})
                samples[sid][key][val] = f
            elif len(g) > 1:
                # 여러 변수가 있을 때 복합 key
                key = "_".join([f"{k}:{v}" for k, v in g.items()])
                samples[sid][key] = f
            else:
                samples[sid]["file"] = f

        return samples

    # --------------------------
    # 실행(또는 드라이런)
    # --------------------------
    def run(self, run: bool = False):
        samples = self.discover_samples()
        if not samples:
            print("[WAVE] No samples discovered.")
            return

        print(f"[WAVE] Found {len(samples)} samples:")
        for sid, info in samples.items():
            print(f"  ▶ {sid}")
            for key, val in info.items():
                print(f"     {key}: {val}")

        for sid, sample_info in samples.items():
            for task_name, params in self.tasks.items():
                cmd = self._build_task_cmd(task_name, sid, sample_info, params)
                print(f"\n[{sid}] {task_name} → {cmd}")
                if run:
                    sp.call(cmd, shell=True)

    # --------------------------
    # Task Command 생성 (간단히 예시)
    # --------------------------
    def _build_task_cmd(self, task_name: str, sid: str, sample_info: Dict[str, Any], params: Dict[str, Any]):
        work_dir = self.work_dir / sid / task_name
        work_dir.mkdir(parents=True, exist_ok=True)

        # 가장 기본적으로 read_info 변수가 있을 때 처리
        reads = sample_info.get("read_info", {})
        inputs = " ".join(str(v) for v in reads.values()) if reads else " ".join(
            str(v) for v in sample_info.values() if isinstance(v, Path)
        )

        threads = params.get("threads", 2)
        cmd = (
            f"fastqc --threads {threads} "
            f"--outdir {work_dir} "
            f"{inputs}"
        )
        return cmd


# --------------------------
# 실행 예시
# --------------------------
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser("wave-flex")
    ap.add_argument("--config", required=True)
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()

    wf = Workflow(Path(args.config))
    wf.run(run=args.run)
