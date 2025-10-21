from __future__ import annotations
from pathlib import Path
from typing import List, Sequence
import shlex


# 로컬 바이너리 버전
def build_fastq_cmd(*, qcResDir: str, RawFastqDir: str, SeqID: str, threads: int = 2, extract: bool = True, exe: str = "fastqc") -> List[str]:
    out = Path(qcResDir)
    extract_flag = ["--extract"] if extract else []
    cmd: List[str] = [
        exe,
        *extract_flag,
        "--threads", str(int(threads)),
        "--outdir", qcResDir,
        f"{RawFastqDir}/{SeqID}_R1.fastq.gz",
        f"{RawFastqDir}/{SeqID}_R2.fastq.gz",
    ]
    return cmd