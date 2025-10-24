# tasks/fastqc/_func.py
from __future__ import annotations
from pathlib import Path
from typing import List, Optional, Iterable


def build_fastqc_cmd(
        qcResDir: str,
        RawFastqDir: str,
        SeqID: str,
        threads: int = 2,
        extract: bool = True,
        image: Optional[str] = None,
        binds: Optional[List[str] | str] = None,
    ) -> Iterable[str]:
    """
    Build a shell command for running FastQC.

    Parameters
    ----------
    qcResDir : str
        Output directory for FastQC results
    RawFastqDir : str
        Directory containing FASTQ(.gz) files
    SeqID : str
        Sample ID (used to filter FASTQs)
    threads : int
        Number of threads
    extract : bool
        Whether to extract FastQC reports (.zip)
    image : str, optional
        Singularity image path (if None, run directly)
    binds : list[str] or str, optional
        Singularity bind directories, e.g. ["/storage", "/data"] or "/storage,/data"

    Returns
    -------
    List[str]
        A list of shell command lines to execute
    """

    raw_dir = Path(RawFastqDir)
    if not raw_dir.exists():
        raise FileNotFoundError(f"[FastQC] RawFastqDir not found: {raw_dir}")

    # 🔹 SeqID에 해당하는 FASTQ 파일 자동 탐색
    fastqs = sorted(raw_dir.glob(f"*{SeqID}*.fastq*"))
    if not fastqs:
        raise FileNotFoundError(f"[FastQC] No FASTQ files found for SeqID '{SeqID}' in {raw_dir}")

    fastq_str = " ".join(str(f) for f in fastqs)

    # 🔹 extract 옵션
    extract_flag = "--extract" if extract else ""

    # 🔹 Singularity 실행 여부
    if image:
        # binds: list or comma string
        if isinstance(binds, list):
            bind_opt = "-B " + ",".join(binds)
        elif isinstance(binds, str):
            bind_opt = "-B " + binds
        else:
            bind_opt = ""

        cmd = (
            f"singularity exec {bind_opt} {image} "
            f"fastqc {extract_flag} "
            f"--threads {threads} "
            f"--outdir {qcResDir} "
            f"{fastq_str}"
        )
    else:
        cmd = (
            f"fastqc {extract_flag} "
            f"--threads {threads} "
            f"--outdir {qcResDir} "
            f"{fastq_str}"
        )

    # 반환 형태는 Iterable[str] — 여러 줄 스크립트로 쓸 수 있게
    return [cmd]
