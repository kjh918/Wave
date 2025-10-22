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
from src.job_manager import JobManager
## 

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="run",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--config", dest="config", help="config.yaml")
    return p

# def qsub_job():


def main():
    

    # job_manager = JobManager(
    #     queue_node = 'all.q@ngsnode1', 
    #     scheduler = "SGE"
    # )

    # parser = Parser()
    # rawdata_path = '/storage/home/jhkim/Projects/cbNIPT/GCX-MakeInsilicoUPDData-2025-10/Resources/DownsampledData'
    # work_dir_path = '/storage/home/jhkim/Projects/cbNIPT/GCX-MakeInsilicoUPDData-2025-10/Results'
    
    # analysis_data_dict = parser(rawdata_path, f'{work_dir_path}/input_metadata.tsv', 'fastq')
    
    # BwaIndex='/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa'
    # ReferenceFasta='/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa'
    # KnownSnp='/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz'
    # KnownIndel1='/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.known_indels.vcf.gz'
    # KnownIndel2='/storage/references_and_index/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
    
    # total_cmd = []
    # for sample_id, reads in analysis_data_dict.items():
    #     read1, read2 = reads[0], reads[1]
        
    # #     mapper = Mapper(
    # #         seq_id = sample_id,
    # #         forward_read = read1,
    # #         reverse_read = read2,
    # #         bam_dir = f'{work_dir_path}/{sample_id}/unmapped_bam',
    # #         reference_fasta = ReferenceFasta,
    # #         bwa_index=BwaIndex
    # #     )
    # #     cmd_list = mapper.unmapped_bam()
    # #     total_cmd.append(cmd_list[0])

    # # # job_manager(
    # #     job_name = 'unmapped_bam',
    # #     cmd_list = total_cmd,
    # #     script_dir = Path(work_dir_path),
    # #     check_status = False,
    # #     max_procs = 90,
    # #     poll_interval_sec = 15
    # # )

        
    # # total_cmd = []
    # # for sample_id, reads in analysis_data_dict.items():
    # #     read1, read2 = reads[0], reads[1]
        
    # #     mapper = Mapper(
    # #         seq_id = sample_id,
    # #         forward_read = read1,
    # #         reverse_read = read2,
    # #         bam_dir = f'{work_dir_path}/{sample_id}/aligned_bam',
    # #         reference_fasta = ReferenceFasta,
    # #         bwa_index=BwaIndex
    # #     )
    # #     cmd_list = mapper.mapped_bam()
    # #     total_cmd.append(cmd_list[0])
        
    # # job_manager(
    # #     job_name = 'aligned_bam',
    # #     cmd_list = total_cmd,
    # #     script_dir = Path(work_dir_path),
    # #     check_status = False,
    # #     max_procs = 90,
    # #     poll_interval_sec = 15
    # # )

    # # total_cmd = []
    # # for sample_id, reads in analysis_data_dict.items():
    # #     read1, read2 = reads[0], reads[1]
        
    # #     mapper = Mapper(
    # #         seq_id = sample_id,
    # #         forward_read = read1,
    # #         reverse_read = read2,
    # #         bam_dir = f'{work_dir_path}/{sample_id}/merged_bam',
    # #         reference_fasta = ReferenceFasta,
    # #         bwa_index=BwaIndex
    # #     )
    # #     cmd_list = mapper.merged_bam_files()
    # #     total_cmd.append(cmd_list[0])

    # # job_manager(
    # #     job_name = 'merged_bam',
    # #     cmd_list = total_cmd,
    # #     script_dir = Path(work_dir_path),
    # #     check_status = False,
    # #     max_procs = 90,
    # #     poll_interval_sec = 15
    # # )
    
    # bam_list = glob(f'{work_dir_path}/*/merged_bam/*primary.bam')
    
    # total_cmd = []
    # for bam_path in bam_list:

    #     sample_id = os.path.basename(bam_path).split('.')[0]
    #     af_mapping = AfterMapping(
    #         seq_id = sample_id,
    #         bam_dir = f'{work_dir_path}/{sample_id}/merged_bam',
    #         qc_dir = f'{work_dir_path}/{sample_id}/qc_bam',
    #         reference_fasta = ReferenceFasta,
    #     )
    #     cmd_list = af_mapping.build_sort_index_cmds()
    #     total_cmd.append(cmd_list[0])
    
    # # job_manager(
    # #     job_name = 'merged_bam',
    # #     cmd_list = total_cmd,
    # #     script_dir = Path(work_dir_path),
    # #     check_status = False,
    # #     max_procs = 90,
    # #     poll_interval_sec = 15
    # # )

    # total_cmd = []
    # for bam_path in bam_list:

    #     sample_id = os.path.basename(bam_path).split('.')[0]
    #     af_mapping = AfterMapping(
    #         seq_id = sample_id,
    #         bam_dir = f'{work_dir_path}/{sample_id}/merged_bam',
    #         qc_dir = f'{work_dir_path}/{sample_id}/qc_bam',
    #         reference_fasta = ReferenceFasta,
    #     )
    #     cmd_list = af_mapping.build_mark_duplicates_cmd(

    #     )
    #     total_cmd.append(cmd_list[0])
    #     print(cmd_list)
    
    # # job_manager(
    # #     job_name = 'merged_bam',
    # #     cmd_list = total_cmd,
    # #     script_dir = Path(work_dir_path),
    # #     check_status = False,
    # #     max_procs = 90,
    # #     poll_interval_sec = 15 
    # # )
    
    # total_cmd = []
    # for bam_path in bam_list:

    #     sample_id = os.path.basename(bam_path).split('.')[0]
    #     af_mapping = AfterMapping(
    #         seq_id = sample_id,
    #         bam_dir = f'{work_dir_path}/{sample_id}/merged_bam',
    #         qc_dir = f'{work_dir_path}/{sample_id}/qc_bam',
    #         reference_fasta = ReferenceFasta,
    #     )
    #     cmd_list = af_mapping.build_bqsr_cmds(

    #     )
    #     total_cmd.append(cmd_list[0])
    #     print(cmd_list)
    
    # # job_manager(
    # #     job_name = 'merged_bam',
    # #     cmd_list = total_cmd,
    # #     script_dir = Path(work_dir_path),
    # #     check_status = False,
    # #     max_procs = 90,
    # #     poll_interval_sec = 15
    # # )
    
    # total_cmd = []
    # for bam_path in bam_list:

    #     sample_id = os.path.basename(bam_path).split('.')[0]
    #     af_mapping = AfterMapping(
    #         seq_id = sample_id,
    #         bam_dir = f'{work_dir_path}/{sample_id}/merged_bam',
    #         qc_dir = f'{work_dir_path}/{sample_id}/qc_bam',
    #         reference_fasta = ReferenceFasta,
    #     )
    #     cmd_list = af_mapping.build_cmds_split_bam_by_chrom(

    #     )
    #     total_cmd.append(cmd_list[0])
    
    # job_manager(
    #     job_name = 'filtered_bam',
    #     cmd_list = total_cmd,
    #     script_dir = Path(work_dir_path),
    #     check_status = False,
    #     max_procs = 90,
    #     poll_interval_sec = 15
    # )

if __name__ == '__main__':

    main()

