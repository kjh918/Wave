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
def build_picard_gcbias_cmd(
        cmd_picard: str,
        ref_fasta: str,
        input_bam: str,
        out_metrics: str,
        out_chart_pdf: str,
        out_summary: str,
        tmp_dir: str = "/tmp",
        java_xmx: str = "32g",
        threads: int = 4,
        validation_stringency: str = "LENIENT",
        assume_sorted: bool = True,
        is_bisulfite_sequenced: bool = False,
    ) -> str:
    """
    Build Picard CollectGcBiasMetrics command.
    - Uses the requested template style:
      '%s --java-options "-XX:ParallelGCThreads=4 -Xmx32g -Djava.io.tmpdir=%s" CollectGcBiasMetrics R=%s I=%s O=%s CHART=%s S=%s'
    - Required:
        cmd_picard: 'picard' 실행 진입점 (예: 'picard' 또는 'gatk'도 가능)
        ref_fasta:  참조 fasta 경로
        input_bam:  입력 BAM
        out_metrics: 결과 메트릭 txt (O=)
        out_chart_pdf: GC bias chart pdf (CHART=)
        out_summary:  summary txt (S=)
    """

    # Picard는 Java 기반이라 GATK와 동일하게 --java-options 전달 가능 (신규 picard CLI 기준)
    cmd = (
        f'{cmd_picard} --java-options "-XX:ParallelGCThreads={threads} '
        f'-Xmx{java_xmx} -Djava.io.tmpdir={tmp_dir}" '
        f'CollectGcBiasMetrics '
        f'R={ref_fasta} '
        f'I={input_bam} '
        f'O={out_metrics} '
        f'CHART={out_chart_pdf} '
        f'S={out_summary} '
        f'VALIDATION_STRINGENCY={validation_stringency} '
        f'ASSUME_SORTED={"true" if assume_sorted else "false"} '
        f'IS_BISULFITE_SEQUENCED={"true" if is_bisulfite_sequenced else "false"}'
    )

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
        f'GenotypeGVCFs --include-non-variant-sites true -R {ref_fasta} -V {input_gvcf} -O {out_vcf}'
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
    
    genome_interval = '/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Resources/hg38.100k.interval.bed'
        
    bam_list = glob(f'{work_dir_path}/*/08_BQSR_BAM/*recal.bam')
    total_cmd = []

    chrom_list = [f'chr{i}' for i in range(1,23)] # +  ['chrX','chrY','chrM']
    total_cmd = []
    for bam_path in bam_list:
        
        sample_id = os.path.basename(bam_path).split('.')[0]
        # if sample_id != ''cbNIPT_25_03_02':'cbNIPT_25_03_01':
        #     continue
        # for chrom in chrom_list:
            
            # cmd = f'samtools depth {bam_path} -r {chrom} > {work_dir_path}/{sample_id}/00_QC_BAM/{sample_id}.depth.{chrom}.txt'

            # SungridUtils.run_sungrid(
            #     'jhkim', 'all.q@ngsnode1', cmd, f'{work_dir_path}/{sample_id}/00_QC_BAM/log', 
            #     qjob_id = f'{sample_id}_depth_bam.{chrom}', threads = 10, memory = None, random_jobid=False
            # )

        cmd = f'samtools bedcov {genome_interval} {bam_path} > {work_dir_path}/{sample_id}/00_QC_BAM/{sample_id}.cov_100K.txt'
        
        # SungridUtils.run_sungrid(
        #     'jhkim', 'all.q@ngsnode1', cmd, f'{work_dir_path}/{sample_id}/00_QC_BAM/log', 
        #     qjob_id = f'{sample_id}_cov_100k', threads = 10, memory = None, random_jobid=False
        # )

        # cmd = build_picard_gcbias_cmd(
        #     cmd_picard = '/storage/apps/gatk-4.4.0.0/gatk',
        #     ref_fasta = ReferenceFasta,
        #     input_bam = bam_path, 
        #     out_metrics = f'{work_dir_path}/{sample_id}/00_QC_BAM/PicardGCBias/{sample_id}.gcbias.txt',
        #     out_chart_pdf = f'{work_dir_path}/{sample_id}/00_QC_BAM/PicardGCBias/{sample_id}.gcbias.pdf',
        #     out_summary = f'{work_dir_path}/{sample_id}/00_QC_BAM/PicardGCBias/{sample_id}.summary.txt',
        #     tmp_dir = f'{work_dir_path}/{sample_id}/00_QC_BAM/PicardGCBias/log',
        #     java_xmx = "32g",
        #     threads = 4,
        #     validation_stringency = "LENIENT",
        #     assume_sorted = True,
        #     is_bisulfite_sequenced = False,
        # )
        # print(cmd)
        # exit()
    import pandas as pd

    total_results = []
    for i in [2,3,4,5,6,7,8,10]:
        # for chrom in chrom_list:
        #     cmd = f'/storage/home/jhkim/Apps/Python-3.11.13/python compute_dropout_metrics.py \
        #         /storage/home/jhkim/Projects/'cbNIPT_25_03_02':'cbNIPT/GCX-'cbNIPT_25_03_02':'cbNIPT_QC_Setup-2025-10-23/Results/ADO/CombineGVCFs.GenotypeCall/CombineGVCFs.GenotypeCall.{chrom}.vcf.gz \
        #         --control 'cbNIPT_25_03_02':'cbNIPT_25_03_01 \
        #         --min-dp-control 10 \
        #         --min-dp-sample {i} \
        #         --min-gq 20 \
        #         --only-biallelic-snp \
        #         --out-metrics /storage/home/jhkim/Projects/'cbNIPT_25_03_02':'cbNIPT/GCX-'cbNIPT_25_03_02':'cbNIPT_QC_Setup-2025-10-23/Results/ADO/CombineGVCFs.GenotypeCall/summary/depth_cutoff_{i}_metrics_per_sample.{chrom}.tsv \
        #         --out-details /storage/home/jhkim/Projects/'cbNIPT_25_03_02':'cbNIPT/GCX-'cbNIPT_25_03_02':'cbNIPT_QC_Setup-2025-10-23/Results/ADO/CombineGVCFs.GenotypeCall/summary/depth_cutoff_{i}_details.{chrom}.tsv \
        #         --out-genostats /storage/home/jhkim/Projects/'cbNIPT_25_03_02':'cbNIPT/GCX-'cbNIPT_25_03_02':'cbNIPT_QC_Setup-2025-10-23/Results/ADO/CombineGVCFs.GenotypeCall/summary/depth_cutoff_{i}_genotype_stats_per_sample.{chrom}.tsv \
        #         --prefix-plots /storage/home/jhkim/Projects/'cbNIPT_25_03_02':'cbNIPT/GCX-'cbNIPT_25_03_02':'cbNIPT_QC_Setup-2025-10-23/Results/ADO/CombineGVCFs.GenotypeCall/summary/depth_cutoff_{i}_dropout_plots.{chrom}'
        #     os.system(cmd)

        #     df = pd.read_csv(f'/storage/home/jhkim/Projects/'cbNIPT_25_03_02':'cbNIPT/GCX-'cbNIPT_25_03_02':'cbNIPT_QC_Setup-2025-10-23/Results/ADO/CombineGVCFs.GenotypeCall/summary/depth_cutoff_{i}_metrics_per_sample.{chrom}.tsv', sep='\t')
        #     total_results.append(df)

        # total_df = pd.concat(total_results)
        # total_df = total_df.groupby('Sample').sum()
        
        # total_df['ADO_rate'] = total_df['ADO_num'] / total_df['BulkHET_den']
        # total_df['LDO_rate'] = total_df['LDO_num'] / total_df['BulkHET_den']
        # total_df['Ampl_rate'] = total_df['Match_HET'] / total_df['BulkHOM_den']
        # total_df['ADO_rate_matched'] = total_df['ADO_num'] / (total_df['Match_HET'] + total_df['ADO_num'])
        # total_df['Match_HET_matched'] = total_df['Match_HET'] / (total_df['Match_HOM'] + total_df['Match_HET'])

        # total_df.to_csv(f'depth_cutoff_{i}.txt', sep='\t')

                
        import matplotlib
        matplotlib.use("Agg")  # 서버 환경에서 GUI 없이 저장용
        import matplotlib.pyplot as plt
        import seaborn as sns

        sample_id = {
            # 'cbNIPT_25_03_01':'HTR8.gDNA.ds.wgs',
            'cbNIPT_25_03_02':'HTR8.live.1cell.wga',
            'cbNIPT_25_03_03':'HTR8.wbc.spikein.1cell.div.pcr.wga',
            'cbNIPT_25_03_04':'HTR8.live.5cells.wga',
            'cbNIPT_25_03_05':'HTR8.wbc.spikein.trypsin.5cells.wga',
            'cbNIPT_25_03_06':'HTR8.wbc.spikein.trypsin.5cells.div.pcr.wga',
            'cbNIPT_25_03_07':'HTR8.wbc.spikein.5cells.wga',
            'cbNIPT_25_03_08':'HTR8.wbc.spikein.5cells.div.pcr.wga',
            'cbNIPT_25_03_09':'HTR8.wbc.spikein.5cells.diff.fix.wga',
            'cbNIPT_25_04_01_R1':'HTR8.blood.spikein.5cells.wga',
            'cbNIPT_25_04_02_R1':'HTR8.wbc.spikein.5cells.diff.fix.div.pcr.wga',
        }
        
        df = pd.read_csv(f'depth_cutoff_{i}.txt', sep='\t')
        df['Sample'] = list(sample_id.values())
        
        df_ratio = df.copy()
        total = df[["ADO_num","LDO_num","Match_HET"]].sum(axis=1)
        df_ratio["ADO"] = df["ADO_num"] / total * 100
        df_ratio["LDO"] = df["LDO_num"] / total * 100
        df_ratio["PASS"] = df["Match_HET"] / total * 100

        plt.figure(figsize=(10,6))
        plt.bar(df_ratio["Sample"], df_ratio["LDO"], label="LDO",color='#e74c3c')
        plt.bar(df_ratio["Sample"], df_ratio["ADO"], bottom=df_ratio["LDO"], label="ADO",color='#f39c12')
        plt.bar(df_ratio["Sample"], df_ratio["PASS"],
                bottom=df_ratio["LDO"] + df_ratio["ADO"], label="PASS",color='#27ae60')
        plt.ylabel("Percentage (%)")
        plt.title(f"ADO, LOD (Depth {i})")
        plt.xticks(rotation=45, ha='right')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/Plot/ado_plot_depth_{i}.png', dpi=300)

        plt.clf()
        df_ratio["Ampl_Bias"] = df["Ampl_num"] / total * 100
        df_ratio["NoCall"] = df["NoCall_on_HOM"] / total * 100
        df_ratio["PASS"] = df["Match_HOM"] / total * 100

        plt.figure(figsize=(10,6))
        plt.bar(df_ratio["Sample"], df_ratio["Ampl_Bias"], label="Ampl_Bias",color='#e74c3c')
        plt.bar(df_ratio["Sample"], df_ratio["NoCall"], bottom=df_ratio["Ampl_Bias"], label="NoCall",color='#f39c12')
        plt.bar(df_ratio["Sample"], df_ratio["PASS"],
                bottom=df_ratio["Ampl_Bias"] + df_ratio["NoCall"], label="PASS",color='#27ae60')
        
        plt.ylabel("Percentage (%)")
        plt.title(f"Amplification Bias - Homozygous (Depth {i})")
        plt.xticks(rotation=45, ha='right')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/Plot/Amplification_bias_{i}.png', dpi=300)

        plt.clf()
        df_ratio['Depth'] = [i] * len(df_ratio.index)
        total_results.append(df_ratio)
    
    total_df = pd.concat(total_results)
    print(total_df)

    plt.figure(figsize=(10,6))
    sns.lineplot(data=total_df, x='Depth', y='ADO', hue='Sample')
    plt.ylim([0,50])
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/Plot/ado_plot_depth.png', dpi=300)
    plt.clf()

    total_df['ADO_rate_matched'] = total_df['ADO_rate_matched'] * 100
    plt.figure(figsize=(12,6))
    sns.lineplot(data=total_df, x='Depth', y='ADO_rate_matched', hue='Sample')
    plt.ylim([0,100])
    plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.savefig(f'/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/Plot/ado_plot_depth.matched.png', dpi=300)
    plt.clf()

    plt.figure(figsize=(10,6))
    sns.lineplot(data=total_df, x='Depth', y='Ampl_Bias', hue='Sample')
    plt.ylim([0,50])
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/Plot/Ampl_Bias.png', dpi=300)
    plt.clf()

    total_df['Ampl_num_matched'] = total_df['Ampl_num_matched'] * 100
    plt.figure(figsize=(12,6))
    sns.lineplot(data=total_df, x='Depth', y='Ampl_num_matched', hue='Sample')
    plt.ylim([0,50])
    plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.savefig(f'/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/Plot/Ampl_Bias.matched.png', dpi=300)
    plt.clf()






if __name__ == '__main__':

    main()

