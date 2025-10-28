from __future__ import annotations
from pathlib import Path
from typing import Iterable, Optional, List, Union
import shlex


def build_fastqc_cmd(
        *,
        inputs: List[str],
        out_dir: str,
        threads: int = 2,
        extract: bool = True,
        image: Optional[str] = "/storage/images/fastqc-0.12.1.sif",
        binds: Optional[Union[List[str], str]] = ("/storage,/data"),
        fastqc_bin: str = "fastqc",
        singularity_bin: str = "singularity",
    ) -> Iterable[str]:
    """
    FastQC 실행 커맨드 생성:
      - inputs: R1/R2 경로 리스트
      - out_dir: 결과 디렉토리
      - image가 있으면 singularity exec, 없으면 로컬 바이너리 사용
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    if not inputs:
        raise ValueError("build_fastqc_cmd: inputs is empty")

    files = " ".join(shlex.quote(x) for x in inputs)
    opt_extract = "--extract " if extract else ""
    base = (
        f"{shlex.quote(fastqc_bin)} "
        f"{opt_extract}"
        f"--threads {threads} "
        f"--outdir {shlex.quote(out_dir)} "
        f"{files}"
    )

    if image:
        if isinstance(binds, list):
            bind_opt = "-B " + ",".join(binds) if binds else ""
        elif isinstance(binds, str):
            bind_opt = f"-B {binds}" if binds else ""
        else:
            bind_opt = ""
        cmd = f"{shlex.quote(singularity_bin)} exec {bind_opt} {shlex.quote(image)} {base}"
    else:
        cmd = base

    return [
        f"mkdir -p {shlex.quote(out_dir)}",
        cmd,
    ]
