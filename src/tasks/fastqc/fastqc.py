# tasks/fastqc/fastqc.py
from __future__ import annotations
from typing import Iterable, Dict, Any, List
from pathlib import Path
from tasks.task import Task, TaskRegistry
from ._func import build_fastqc_cmd

@TaskRegistry.register
class FastQCTask(Task):
    TYPE = "fastqc"
    INPUTS = {
        "RawFastqDir": "Input FASTQ dir",
        "SeqID": "Sample ID",
        "threads": "int",
        "extract": "bool",
        "image": "(optional) Singularity image path",
        "binds": "(optional) bind list or comma-string",
        "qcResDir": "Output dir (optional)",
    }
    DEFAULTS = {"threads": 2, "extract": True}
    OPTIONAL = {"image", "binds", "qcResDir"}

    def to_sh(self) -> Iterable[str]:
        p: Dict[str, Any] = self.params
        seqid = str(p["SeqID"])

        # ✅ output 디렉토리 자동 설정 (없을 경우)
        qcResDir = p.get("qcResDir") or str(Path(self.workdir) / f"{seqid}_fastqc")
        Path(qcResDir).mkdir(parents=True, exist_ok=True)

        return build_fastqc_cmd(
            qcResDir=qcResDir,
            RawFastqDir=str(p["RawFastqDir"]),
            SeqID=seqid,
            threads=int(p.get("threads", 2)),
            extract=bool(p.get("extract", True)),
        )