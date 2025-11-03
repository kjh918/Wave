#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

def build_haplotypecaller_cmd(
        chromosome: str,
        cmd_gatk: str,
        ref_fasta: str,
        input_bam: str,
        dbsnp_vcf: str,
        out_gvcf: str,
        log_path: str,
        gender: str = "UNKNOWN",
        tmp_dir: str = "/tmp",
        threads: int = 4,
        java_xmx: str = "32g",
    ) -> str:
    """
    Build GATK HaplotypeCaller command for a given chromosome and gender.
    - Handles ploidy rules automatically.
    - Returns a shell command string.
    """

    java_opts = (
        f'--java-options "-XX:ParallelGCThreads={threads} '
        f'-Xmx{java_xmx} -Djava.io.tmpdir={tmp_dir}"'
    )

    base = (
        f'{cmd_gatk} {java_opts} HaplotypeCaller '
        f'-R {ref_fasta} -I {input_bam} -L {chromosome} '
        f'-stand-call-conf 30 --dbsnp {dbsnp_vcf} '
        f'-O {out_gvcf} -ERC GVCF'
    )

    # --- Ploidy rules ---
    if chromosome == "chrX":
        if gender.upper() == "MALE":
            cmd = (
                f'{cmd_gatk} {java_opts} HaplotypeCaller '
                f'-ploidy 1 -R {ref_fasta} -I {input_bam} -L {chromosome} '
                f'-stand-call-conf 30 --dbsnp {dbsnp_vcf} '
                f'-O {out_gvcf} -ERC GVCF'
            )
        else:
            cmd = base
    elif chromosome in ["chrY", "chrM"]:
        cmd = (
            f'{cmd_gatk} {java_opts} HaplotypeCaller '
            f'-ploidy 1 -R {ref_fasta} -I {input_bam} -L {chromosome} '
            f'-stand-call-conf 30 --dbsnp {dbsnp_vcf} '
            f'-O {out_gvcf} -ERC GVCF'
        )
    else:
        cmd = base

    return cmd


def main():
    parser = argparse.ArgumentParser(
        description="Build HaplotypeCaller command by chromosome and gender."
    )

    parser.add_argument("--chromosome", required=True, help="Chromosome name (e.g. chr1, chrX, chrY, chrM)")
    parser.add_argument("--cmd_gatk", required=True, help="Path or singularity command for gatk (e.g. singularity exec ... gatk)")
    parser.add_argument("--ref_fasta", required=True, help="Reference FASTA file path")
    parser.add_argument("--input_bam", required=True, help="Input BAM file path")
    parser.add_argument("--dbsnp_vcf", required=True, help="dbSNP VCF file path")
    parser.add_argument("--out_gvcf", required=True, help="Output GVCF file path (e.g. sample.chr1.g.vcf.gz)")
    parser.add_argument("--log_path", required=True, help="Log file path (for 2> redirection)")
    parser.add_argument("--gender", default="UNKNOWN", help="Gender type: MALE / FEMALE / UNKNOWN")
    parser.add_argument("--tmp_dir", default="/tmp", help="Temporary directory for Java I/O")
    parser.add_argument("--threads", type=int, default=4, help="Number of GC threads for Java")
    parser.add_argument("--java_xmx", default="32g", help="Java Xmx memory setting (e.g. 16g, 32g)")

    args = parser.parse_args()

    cmd = build_haplotypecaller_cmd(
        chromosome=args.chromosome,
        cmd_gatk=args.cmd_gatk,
        ref_fasta=args.ref_fasta,
        input_bam=args.input_bam,
        dbsnp_vcf=args.dbsnp_vcf,
        out_gvcf=args.out_gvcf,
        log_path=args.log_path,
        gender=args.gender,
        tmp_dir=args.tmp_dir,
        threads=args.threads,
        java_xmx=args.java_xmx,
    )

    print(cmd)


if __name__ == "__main__":
    main()
