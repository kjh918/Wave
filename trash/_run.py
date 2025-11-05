import os
import sys
import argparse
import subprocess
from glob import glob
from typing import Dict, Any, List, Optional
import pyaml 
from pathlib import Path

## 
from src.parser import Parser
from src.align.bwa import Mapper
from src.align.bwa_after_mapping import AfterMapping
from executor import SungridUtils
## 

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="run",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--config", dest="config", help="config.yaml")
    return p

# def qsub_job():

def build_haplotypecaller_cmd(
        chromosome: str,
        cmd_gatk: str,
        ref_fasta: str,
        input_bam: str,
        dbsnp_vcf: str,
        out_gvcf: str,
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

def build_genotypegvcfs_cmd(
        cmd_gatk: str,
        ref_fasta: str,
        input_gvcf: str,
        out_vcf: str,
        tmp_dir: str = "/tmp",
        threads: int = 4,
        java_xmx: str = "32g",
    ) -> str:
    """
    Build GATK GenotypeGVCFs command.
    - Uses the requested template:
      '%s --java-options "-XX:ParallelGCThreads=4 -Xmx32g -Djava.io.tmpdir=%s" GenotypeGVCFs -R %s -V %s -O %s'
    """

    cmd = (
        f'{cmd_gatk} --java-options "-XX:ParallelGCThreads={threads} '
        f'-Xmx{java_xmx} -Djava.io.tmpdir={tmp_dir}" '
        f'GenotypeGVCFs -R {ref_fasta} -V {input_gvcf} -O {out_vcf}'
    )

    return cmd
def main():

    parser = Parser()
    rawdata_path = '/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Resources/Rawdata/PicoPLEXGold/251103'
    work_dir_path = '/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results'
    
    analysis_data_dict = parser(rawdata_path, f'{work_dir_path}/input_metadata.tsv', 'fastq')
    
    BwaIndex='/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa'
    ReferenceFasta='/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa'
    KnownSnp='/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz'
    KnownIndel1='/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.known_indels.vcf.gz'
    KnownIndel2='/storage/references_and_index/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
    

    chrom_list = [f'chr{i}' for i in range(1,23)]  # +  ['chrX','chrY','chrM']

    for chrom in chrom_list: 
        # chrom = 'chr22's
        input_gvcfs = glob(f'{work_dir_path}/*/09_HaplotypeCall_BAM/*{chrom}.gvcf.gz')
        input_gvcfs.sort()
        input_gvcfs = ' -V '.join(input_gvcfs)
        output_gvcfs = f'{work_dir_path}/ADO/MergedVCF/merged_gvcf.{chrom}.gvcf.gz'
        cmd = f'/storage/apps/gatk-4.4.0.0/gatk CombineGVCFs -R {ReferenceFasta} -V {input_gvcfs} -O {output_gvcfs}'

        input_gvcfs = glob(f'{work_dir_path}/*/09_HaplotypeCall_BAM/*{chrom}.gvcf.gz')
        input_gvcfs.sort()
        input_gvcfs = ' '.join(input_gvcfs)

        
        output_gvcfs = f'{work_dir_path}/ADO/MergedVCF/merged_gvcf.bcftools.{chrom}.gvcf.gz'
        
        cmd = f'/storage/home/jhkim/Apps/Python-3.11.13/python /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/classify_ado_bias_imbalance.py \
            {output_gvcfs} \
            --control cbNIPT_25_03_01 \
            --min-dp-control 10 \
            --min-dp-other 2 \
            --min-gq 20 \
            --ab-het-low 0.3 --ab-het-high 0.7 \
            --only-biallelic-snp \
            --site-pass-only \
            --out-summary {work_dir_path}/ADO/summary_per_sample.{chrom}.tsv \
            --out-discordant {work_dir_path}/ADO/discordant_sites.{chrom}.tsv'

        
        SungridUtils.run_sungrid(
            'jhkim', 'all.q@ngsnode1', cmd, f'{work_dir_path}/ADO/MergedVCF/log', 
            qjob_id = f'cal_ado_{chrom}', threads = 10, memory = None, random_jobid=False
        )
        
    


if __name__ == '__main__':

    main()

