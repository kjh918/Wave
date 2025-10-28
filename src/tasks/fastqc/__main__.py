from __future__ import annotations
from typing import Dict, Any, List
from pathlib import Path

from src.tasks.task import Task
from ._func import build_fastqc_cmd


class FastQCRunner(Task):
    """
    - inputs:  {"fastq_r1": "/abs/S1_R1.fastq.gz", "fastq_r2": "/abs/S1_R2.fastq.gz"}
    - outputs: {"dir": "/abs/.../00_Rawdata_FastQC"}
    - params:  {"threads": 2, "extract": True, "image": "...sif", "binds": ["/storage","/data"]}
    """

    # 메타(선택): 워크플로에서 스펙 확인용
    TYPE = "fastqc"
    INPUT_SPEC  = "map"    # 입력은 dict 형태
    OUTPUT_SPEC = "dir"    # 출력은 작업 디렉토리 중심
    DEFAULTS: Dict[str, Any] = {
        "threads": 2,
        "extract": True,
        "image": "/storage/images/fastqc-0.12.1.sif",
        "binds": ["/storage", "/data"],
        "fastqc_bin": "fastqc",
        "singularity_bin": "singularity",
    }

    # (선택) 기대하는 키를 명시하고 싶다면:
    INPUT_KEYS  = ("fastq_r1", "fastq_r2")
    OUTPUT_KEYS = ("dir",)

    def to_sh(self) -> List[str]:
        
        # --- 입력 검증 & 정리 ---
        if not isinstance(self.inputs, dict):
            raise TypeError("[FastQCRunner] inputs must be a dict")
        try:
            r1 = str(self.inputs["fastq_r1"])
            r2 = str(self.inputs["fastq_r2"])
        except KeyError as e:
            raise KeyError(f"[FastQCRunner] missing input key: {e}")

        # --- 출력 디렉토리 ---
        if not isinstance(self.outputs, dict) or "dir" not in self.outputs:
            raise TypeError("[FastQCRunner] outputs must be a dict containing 'dir'")
        out_dir = str(self.outputs["dir"])
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        # --- 파라미터 병합(기본값 -> 사용자값 우선) ---
        p: Dict[str, Any] = dict(self.DEFAULTS)
        p.update(self.params or {})

        # --- 커맨드 생성 ---
        return list(
            build_fastqc_cmd(
                inputs=[r1, r2],
                out_dir=out_dir,
                threads=int(p["threads"]),
                extract=bool(p["extract"]),
                image=p.get("image"),
                binds=p.get("binds"),
                fastqc_bin=str(p.get("fastqc_bin", "fastqc")),
                singularity_bin=str(p.get("singularity_bin", "singularity")),
            )
        )


# 워크플로 로더가 찾는 심볼
TASK_CLASS = FastQCRunner
