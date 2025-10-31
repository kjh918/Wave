# src/tasks/fastqc/_func.py
from __future__ import annotations
from typing import List, Optional, Sequence


def _norm_binds(binds: Optional[Sequence[str] | str]) -> str:
    """
    binds가 리스트면 콤마-조인하여 단일 -B 인자로 사용,
    문자열이면 그대로 사용. 없으면 빈 문자열 반환.
    """
    if binds is None:
        return ""
    if isinstance(binds, (list, tuple)):
        joined = ",".join(str(b) for b in binds if str(b).strip())
        return f"-B {joined}" if joined else ""
    s = str(binds).strip().strip('[').strip(']').replace(' ','')
    return f"-B {s}" if s else ""


def build_fastqc_cmd(
        *,
        inputs: Sequence[str],
        out_dir: str,
        threads: int = 2,
        extract: bool = True,
        image: Optional[str] = None,
        fastqc_bin: str = "fastqc",
        singularity_bin: str = "singularity",
        binds: Optional[Sequence[str] | str] = None,
    ) -> List[str]:
    """
    Build full fastqc command (supports singularity or local execution).

    Parameters
    ----------
    inputs : list[str]
        One or two FASTQ files.
    out_dir : str
        Output directory path.
    threads : int
        Number of threads.
    extract : bool
        Whether to use --extract option.
    image : str, optional
        Singularity image path. If None, run locally.
    binds : list[str], optional
        List of paths to bind when using singularity.
    fastqc_bin : str
        FastQC binary name (default: 'fastqc')
    singularity_bin : str
        Singularity binary name (default: 'singularity')

    Returns
    -------
    list[str]
        Single command line string ready to write into .sh
    """
    if not inputs or len(inputs) == 0:
        raise ValueError("[build_fastqc_cmd] no input FASTQs provided")

    # 기본 FastQC argv
    cmd: List[str] = [fastqc_bin]
    if extract:
        cmd.append("--extract")
    cmd += ["--threads", str(int(threads)), "--outdir", str(out_dir)]
    cmd += list(map(str, inputs))
    # --- Singularity 사용 여부 ---
    # print(_norm_binds(binds))
    # exit()
    if image:
        cmd = [
            singularity_bin,
            "exec",
            _norm_binds(binds),
            image,
            *cmd,  # fastqc 부분 그대로 삽입
        ]
    # 최종 문자열
    return [" ".join(cmd)]