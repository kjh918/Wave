# flexible_workflow.py
from __future__ import annotations
import re, yaml, shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Any, List, Optional

from .task import Task

@dataclass
class Workflow:

    config_path: Optional[Path] = None           # 파일 경로 or None
    config_dict: Optional[Dict[str, Any]] = None # dict directly 전달 가능

    cfg: Dict[str, Any] = field(init=False)
    work_dir: Path = field(init=False)
    input_tpl: str = field(init=False)
    workflow: Dict[str, Any] = field(default_factory=dict)
    _tasks_templates: List[TaskTemplate] = field(default_factory=list)

    def __post_init__(self):
        # --- 1. config 읽기 ---
        if self.config_dict:
            self.cfg = self.config_dict
        elif self.config_path:
            self.cfg = yaml.safe_load(Path(self.config_path).read_text())
        else:
            raise ValueError("Workflow: either config_path or config_dict must be provided.")

        # --- 2. 기본 설정 로드 ---
        w_params = self.cfg["WORK_PARAMETERS"]
        self.input_tpl = w_params["input_path"]
        self.work_dir = Path(w_params["work_dir_path"])
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # --- 3. workflow skeleton ---
        self.workflow = {
            "SETTING": self.cfg.get("SETTING", {}),
            "WORK_PARAMETERS": w_params,
            "SAMPLES": {},
        }

        # --- 4. TASK_LIST 로드 ---
        self._tasks_templates = [
            Task(
                name=t["name"],
                type=t["type"],
                params=t.get("params", {}) or {}
            )
            for t in (self.cfg.get("TASK_LIST", []) or [])
        ]

    # --------------------------
    # 유틸: 템플릿 파서
    # --------------------------
    def _parse_placeholders(self, tpl: str) -> List[str]:
        """{name} 패턴 추출"""
        # {...} 안의 식별자만 뽑기
        return re.findall(r"\{([a-zA-Z0-9_]+)\}", tpl)


    def _tpl_to_regex(self, tpl: str, placeholders: List[str]) -> re.Pattern:
        """플레이스홀더를 정규식 그룹으로 변환 (예: {sample_id}, {read_info})"""
        esc = re.escape(tpl)
        for ph in placeholders:
            if ph == "sample_id":
                esc = esc.replace(re.escape(f"{{{ph}}}"), r"(?P<sample_id>[^/]+)")
            else:
                esc = esc.replace(re.escape(f"{{{ph}}}"), rf"(?P<{ph}>[^/]+)")
        return re.compile(rf"^{esc}$")


    def _tpl_to_glob(self, tpl: str, placeholders: List[str]) -> str:
        """{name} → * 로 치환해 glob 패턴 생성"""
        pat = tpl
        for ph in placeholders:
            pat = pat.replace(f"{{{ph}}}", "*")
        return pat


    def discover_samples(self) -> Dict[str, Dict[str, Any]]:
        tpl = self.input_tpl
        placeholders = self._parse_placeholders(tpl)
        if "sample_id" not in placeholders:
            raise ValueError("input_path template must include {sample_id}")

        rx = self._tpl_to_regex(tpl, placeholders)
        glob_pat = self._tpl_to_glob(tpl, placeholders)

        # 절대경로 기준 glob (tpl이 절대경로라고 가정)
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

            # 초기 구조 생성
            s = samples.setdefault(sid, {"inputs": {}})
            s_inputs = s["inputs"]

            # raw_dir 등록(최초 값 유지)
            s_inputs.setdefault("raw_dir", str(f.parent))

            if len(g) == 0:
                # 플레이스홀더가 sample_id 말고는 없는 템플릿
                s_inputs["single"] = str(f)
                continue

            # 여러 플레이스홀더가 있을 수 있으니 전부 키별로 정리
            for key, val in g.items():
                # 예: key == "read" / val == "1" or "2"
                bucket = s_inputs.setdefault(key, {})
                # 동일 키/값이 여러 파일에 매치될 수 있으니 마지막 것으로 갱신 (필요시 list로 누적)
                bucket[str(val)] = str(f)

        return samples
        # placeholder 렌더
    def _render(self, s: str, ctx: Dict[str, Any]) -> str:
        return re.sub(r"\{([a-zA-Z0-9_]+)\}", lambda m: str(ctx.get(m.group(1), "")), s or "")

    def _render_list(self, xs: List[str], ctx: Dict[str, Any]) -> List[str]:
        return [self._render(x, ctx) for x in (xs or [])]

    # discover_samples()는 기존 그대로 (raw_dir/read dict 구성)

    def _build_task_script(self, sample_id: str, tpl: TaskTemplate) -> Dict[str, Any]:
        # 클래스 기본값/옵션/입력 스키마 꺼냄
        defaults = getattr(tpl.cls, "DEFAULTS", {}) or {}
        inputs_schema = getattr(tpl.cls, "INPUTS", {}) or {}

        ctx = {"sample_id": sample_id, "work_dir_path": str(self.work_dir)}
        # 사용자 파라미터 플레이스홀더 치환
        up = {}
        for k, v in (tpl.user_params or {}).items():
            up[k] = self._render_list(v, ctx) if isinstance(v, list) else self._render(v, ctx)

        # 출력 디렉토리: params.output가 있으면 우선
        out_dir = up.get("output") or str(Path(self.work_dir) / sample_id / tpl.name)
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        # 입력 파일 리스트 해석 (절대/상대/파일명만 → raw_dir 기준 완성)
        raw_dir = Path(self.workflow["SAMPLES"][sample_id]["inputs"]["raw_dir"])
        in_list = up.get("input", [])
        inputs_abs = []
        for token in in_list:
            p = Path(token)
            if p.is_absolute(): inputs_abs.append(str(p))
            else: inputs_abs.append(str(raw_dir / token))

        # 클래스 기본값과 merge한 params
        params = dict(defaults)
        params.update({
            "inputs": inputs_abs,
            "outDir": out_dir,
            "sample_id": sample_id,
            "SeqID": sample_id,
            "RawFastqDir": str(raw_dir),
        })
        for k, v in up.items():
            if k not in {"input", "output"}:
                params[k] = v

        # Task 인스턴스 생성 및 쉘 빌드
        task = tpl.cls(workdir=out_dir, params=params)
        script = "\n".join(list(task.to_sh()))
        return {"name": tpl.name, "tool": tpl.tool, "func": tpl.func,
                "version": tpl.version or "", "workdir": out_dir,
                "script": script, "threads": int(params.get("threads", 1))}
                
    def materialize(self) -> None:
        for sid, sdata in self.workflow["SAMPLES"].items():
            concrete = []
            for tpl in self._tasks_templates:
                concrete.append(self._build_task_script(sid, tpl))
            sdata["tasks"] = concrete
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
