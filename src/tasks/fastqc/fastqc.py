from __future__ import annotations
from typing import Iterable, Dict, Any
from tasks.task import Task, TaskRegistry
from ._func import build_fastqc_cmds

@TaskRegistry.register
class FastQCRunner(Task):
    
    TYPE = "fastqc"

    # 모든 파라미터를 입력으로 받는 형태 (r2는 페어가 아닐 수 있으므로 optional)
    INPUTS   = {
        "r1": "R1 fastq", 
        "r2": "R2 fastq (optional)", 
        "threads": "int", 
        "outdir": "output dir"
        }
    OUTPUTS  = {}  # 보고서 위치를 강제하지 않고 fastqc 기본 규칙을 따름
    DEFAULTS = {}  # "모두 받는 형태" 요청에 따라 기본값 없음
    OPTIONAL = {"r2"}

    def __call__(self) -> Iterable[str]:
        p: Dict[str, Any] = self.params
        return build_fastqc_cmds(
            r1=str(p["r1"]),
            r2=str(p["r2"]) if p.get("r2") else None,
            outdir=str(p["outdir"]),
            threads=int(p["threads"]),
        )