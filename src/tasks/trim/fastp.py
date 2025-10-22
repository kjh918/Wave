#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path


def build_fastp_cmd(
        sample_id: str,
        fastq_r1: str,
        fastq_r2: str,
        out_dir: str,
        threads: int = 8,
        adapter_sequence: str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        adapter_sequence_r2: str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        length_required: int = 100,
        average_qual: int = 10,
        qualified_quality_phred: int = 15,
        trim_front1: int = 0,
        trim_front2: int = 0,
        bind_dir_list: list[str] | None = None,
        singularity: str = "/storage/apps/singularity/bin/singularity",
        fastp_sif: str = "/storage/images/fastp-0.23.4.sif",
        fastp_exe: str = "/usr/bin/fastp",
    ) -> str:
    """
        Build Fastp command (paired-end mode).
        Returns a shell command string.

        Parameters
        ----------
        fastq_r1, fastq_r2: input FASTQ.gz paths
        out_r1, out_r2: trimmed FASTQ.gz outputs
        json_path, html_path: output report paths
        threads: number of threads
        adapter_sequence(_r2): adapter sequences for trimming
        bind_dir_list: singularity bind paths
        singularity: path to singularity executable
        fastp_sif: path to fastp container image
        fastp_exe: path to fastp inside container
    """

    Path(out_r1).parent.mkdir(parents=True, exist_ok=True)
    Path(out_r2).parent.mkdir(parents=True, exist_ok=True)
    Path(json_path).parent.mkdir(parents=True, exist_ok=True)
    Path(html_path).parent.mkdir(parents=True, exist_ok=True)

    out_r1 = str(Path(args.out_dir) / f"{sample_id}_R1.trim.fastq.gz")
    out_r2 = str(Path(args.out_dir) / f"{sample_id}_R2.trim.fastq.gz")
    json_path = str(Path(args.out_dir) / f"{sample_id}.fastp.json")
    html_path = str(Path(args.out_dir) / f"{sample_id}.fastp.html")

    bind_dir_list = bind_dir_list or ["/storage", "/data"]

    cmd = (
        f"{singularity} exec "
        f"-B {','.join(bind_dir_list)} "
        f"{fastp_sif} "
        f"{fastp_exe} "
        f"--thread {threads} "
        f"--in1 {fastq_r1} "
        f"--in2 {fastq_r2} "
        f"--out1 {out_r1} "
        f"--out2 {out_r2} "
        f"--json {json_path} "
        f"--html {html_path} "
        f"--trim_poly_g --detect_adapter_for_pe "
        f"--adapter_sequence {adapter_sequence} "
        f"--adapter_sequence_r2 {adapter_sequence_r2} "
        f"--length_required {length_required} "
        f"--average_qual {average_qual} "
        f"--qualified_quality_phred {qualified_quality_phred} "
        f"--trim_front1 {trim_front1} "
        f"--trim_front2 {trim_front2}"
    )

    return cmd


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a fastp command for paired-end read trimming (importable too)."
    )
    parser.add_argument("--sample_id", required=True, help="sample_id")
    parser.add_argument("--fastq-r1", required=True, help="R1 FASTQ path")
    parser.add_argument("--fastq-r2", required=True, help="R2 FASTQ path")
    parser.add_argument("--out-dir", required=True, help="Output DIR")
    # parser.add_argument("--out-r2", required=True, help="Output trimmed R2 path")
    # parser.add_argument("--json", required=True, help="Output JSON report path")
    # parser.add_argument("--html", required=True, help="Output HTML report path")

    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--adapter-seq", default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")
    parser.add_argument("--adapter-seq-r2", default="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
    parser.add_argument("--length-required", type=int, default=100)
    parser.add_argument("--average-qual", type=int, default=10)
    parser.add_argument("--qualified-quality-phred", type=int, default=15)
    parser.add_argument("--trim-front1", type=int, default=0)
    parser.add_argument("--trim-front2", type=int, default=0)

    parser.add_argument("--singularity", default="/storage/apps/singularity/bin/singularity")
    parser.add_argument("--fastp-sif", default="/storage/images/fastp-0.23.4.sif")
    parser.add_argument("--fastp", default="fastp")
    parser.add_argument("--bind", nargs="*", default=["/storage", "/data"],
                        help="Bind directories for singularity")
    return parser.parse_args()


def main():
    args = _parse_args()

    cmd = build_fastp_cmd(
        sample_i=args.sample_id
        fastq_r1=args.fastq_r1,
        fastq_r2=args.fastq_r2,
        output_dir=args.out_dir,
        threads=args.threads,
        adapter_sequence=args.adapter_seq,
        adapter_sequence_r2=args.adapter_seq_r2,
        length_required=args.length_required,
        average_qual=args.average_qual,
        qualified_quality_phred=args.qualified_quality_phred,
        trim_front1=args.trim_front1,
        trim_front2=args.trim_front2,
        bind_dir_list=args.bind,
        singularity=args.singularity,
        fastp_sif=args.fastp_sif,
        fastp_exe=args.fastp_exe,
    )
    return cmd


if __name__ == "__main__":
    main()
