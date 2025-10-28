# tasks/fastp/__main__.py
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess

from ._func import build_fastp_cmd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m tasks.fastp",
        description="fastp task CLI (generate and/or run a single fastp command)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # 필수
    p.add_argument("--rawdir", dest="RawFastqDir", required=True, help="Input FASTQ directory")
    p.add_argument("--seqid", dest="SeqID", required=True, help="Sample ID (expects <SeqID>_R1/2.fastq.gz)")
    p.add_argument("--outdir", dest="TrimFastqDir", required=True, help="Output directory for trimmed FASTQs")

    # 품질/스레드
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--length-required", dest="length_required", type=int, default=100)
    p.add_argument("--average-qual", dest="average_qual", type=int, default=10)
    p.add_argument("--qualified-quality-phred", dest="qualified_quality_phred", type=int, default=15)

    # PicoplexGold 옵션
    p.add_argument("--picoplex-gold", dest="picoplex_gold", action="store_true",
                   help="Enable Picoplex Gold trimming (defaults trim_front1/2=14 unless overridden)")
    p.add_argument("--trim-front1", dest="trim_front1", type=int, default=14)
    p.add_argument("--trim-front2", dest="trim_front2", type=int, default=14)

    # 컨테이너
    p.add_argument("--image", default="/storage/images/fastp-0.23.4.sif",
                   help="Singularity image path; set with --no-image to use local fastp")
    p.add_argument("--no-image", dest="no_image", action="store_true",
                   help="Use local fastp binary (ignore --image)")
    p.add_argument("--bind", dest="binds", action="append",
                   help="Bind path(s) for singularity, repeatable (e.g., --bind /storage --bind /data)")

    # 동작
    p.add_argument("--emit", help="Write a .sh script to this path")
    p.add_argument("--run", action="store_true", help="Run immediately after generation")
    return p


def main():
    ap = build_parser()
    a = ap.parse_args()

    # 출력 디렉토리 준비
    outdir = Path(a.TrimFastqDir)
    outdir.mkdir(parents=True, exist_ok=True)

    image = None if a.no_image else a.image
    binds = a.binds  # None | list[str]

    cmds = build_fastp_cmd(
        RawFastqDir=a.RawFastqDir,
        SeqID=a.SeqID,
        TrimFastqDir=a.TrimFastqDir,
        threads=a.threads,
        picoplex_gold=a.picoplex_gold,
        trim_front1=a.trim_front1,
        trim_front2=a.trim_front2,
        length_required=a.length_required,
        average_qual=a.average_qual,
        qualified_quality_phred=a.qualified_quality_phred,
        image=image,
        binds=binds,
    )

    line = cmds[0]  # build_fastp_cmd returns a list with single command string

    if a.emit:
        out = Path(a.emit)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + line + "\n")
        out.chmod(0o755)
        print(f"[fastp] wrote script: {out}")

    if a.run:
        print(f"[fastp] running: {line}")
        subprocess.run(["bash", "-lc", line], check=True)
    else:
        print(line)


if __name__ == "__main__":
    main()