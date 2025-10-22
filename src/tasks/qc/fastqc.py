#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

def build_fastqc_cmd(
        fastq_r1: str,
        fastq_r2: str | None,
        out_dir: str,
        threads: int = 4,
        singularity_bin: str = "/storage/apps/singularity/bin/singularity",
        fastqc_sif: str = "/storage/images/fastqc-0.12.1.sif",
        fastqc_exe: str = "/usr/local/bin/fastqc",
        bind_dir_list: list[str] | None = None,
        extract: bool = True,
        extra_args: list[str] | None = None,
    ) -> str:
    """
    Build a FastQC command string (supports single/paired FASTQ).
    Returns a shell command string.

    Parameters
    ----------
    fastq_r1: str           R1 FASTQ path
    fastq_r2: str|None      R2 FASTQ path (None for single-end)
    out_dir: str            Output directory for FastQC results
    threads: int            Thread count for FastQC
    singularity_bin: str    Path to singularity binary
    fastqc_sif: str         Path to FastQC .sif
    fastqc_exe: str         FastQC executable inside the container
    bind_dir_list: list[str]|None  Directories to bind (default ['/storage','/data'])
    extract: bool           --extract option (unzip reports)
    extra_args: list[str]|None     Additional FastQC args (e.g. ['--nogroup'])
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    bind_dir_list = bind_dir_list or ["/storage", "/data"]

    parts = [
        singularity_bin, "exec",
        "-B", ",".join(bind_dir_list),
        fastqc_sif,
        fastqc_exe,
    ]
    if extract:
        parts.append("--extract")
    parts += ["--outdir", out_dir, "--threads", str(threads)]

    if extra_args:
        parts += extra_args

    parts.append(fastq_r1)
    if fastq_r2:
        parts.append(fastq_r2)

    return " ".join(map(str, parts))


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build (and print) a FastQC command. Importable as a function too."
    )
    p.add_argument("--r1", required=True, help="R1 FASTQ path")
    p.add_argument("--r2", default=None, help="R2 FASTQ path (optional for paired-end)")
    p.add_argument("--out-dir", required=True, help="FastQC output directory")
    p.add_argument("--threads", type=int, default=4, help="Threads for FastQC")
    p.add_argument("--singularity", default="/storage/apps/singularity/bin/singularity",
                   help="Path to singularity binary")
    p.add_argument("--fastqc-sif", default="/storage/images/fastqc-0.12.1.sif",
                   help="Path to FastQC SIF image")
    p.add_argument("--fastqc", default="fastqc",
                   help="FastQC executable inside container")
    p.add_argument("--bind", nargs="*", default=["/storage", "/data"],
                   help="Bind directories (space-separated)")
    p.add_argument("--no-extract", action="store_true", help="Do not unzip FastQC reports")
    p.add_argument("--extra", nargs="*", default=None,
                   help="Extra FastQC args, e.g. --nogroup --kmers 7")
    return p.parse_args()


def main():
    args = _parse_args()
    cmd = build_fastqc_cmd(
        fastq_r1=args.fastq_r1,
        fastq_r2=args.fastq_r2,
        out_dir=args.out_dir,
        threads=args.threads,
        singularity_bin=args.singularity_bin,
        fastqc_sif=args.fastqc_sif,
        fastqc_exe=args.fastqc_exe,
        bind_dir_list=args.bind,
        extract=not args.no_extract,
        extra_args=args.extra,
    )
    print(cmd)


if __name__ == "__main__":
    main()
