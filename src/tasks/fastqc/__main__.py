# src/tasks/fastqc/__main__.py
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess

from ._func import build_fastqc_cmd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m src.tasks.fastqc",
        description="FastQC task CLI (generate and/or run a single FastQC command)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--r1", required=True, help="R1 FASTQ path")
    p.add_argument("--r2", help="R2 FASTQ path (optional)")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--no-extract", dest="extract", action="store_false", help="Disable --extract")

    # container & bins
    p.add_argument("--image", help="Singularity image path")
    p.add_argument("--bind", dest="binds", action="append", default=[], help="Bind path (repeatable)")
    p.add_argument("--fastqc-bin", default="fastqc")
    p.add_argument("--singularity-bin", default="singularity")

    # behavior
    p.add_argument("--emit", help="Write a .sh script to this path (optional)")
    p.add_argument("--run", action="store_true", help="Run immediately after generation")
    return p


def main():
    ap = build_parser()
    a = ap.parse_args()

    inputs = [a.r1] + ([a.r2] if a.r2 else [])
    cmds = build_fastqc_cmd(
        inputs=inputs,
        out_dir=a.outdir,
        threads=a.threads,
        extract=a.extract,
        image=a.image,
        binds=a.binds,
        fastqc_bin=a.fastqc_bin,
        singularity_bin=a.singularity_bin,
    )

    line = cmds[0]  # we only produce a single command line

    if a.emit:
        out = Path(a.emit)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + line + "\n")
        out.chmod(0o755)
        print(f"[fastqc] wrote script: {out}")

    if a.run:
        print(f"[fastqc] running: {line}")
        subprocess.run(["bash", "-lc", line], check=True)
    else:
        print(line)


if __name__ == "__main__":
    main()