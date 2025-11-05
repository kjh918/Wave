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
    
    fastq_list = glob(f'{work_dir_path}/*25_04*/01_Trimdata_FastP/*_R1.trimed.fastq.gz')

    for read1 in fastq_list:
        sample_id = os.path.basename(read1).split('_R1')[0]
        read2 = read1.replace('R1','R2')

        mapper = Mapper(
            seq_id = sample_id,
            forward_read = read1,
            reverse_read = read2,
            bam_dir = f'{work_dir_path}/{sample_id}/04_Unalign_BAM',
            reference_fasta = ReferenceFasta,
            bwa_index=BwaIndex
        )
        # cmd_list = mapper.unmapped_bam()
        # print(cmd_list)
        # SungridUtils.run_sungrid(
        #     'jhkim', 'all.q@ngsnode1', cmd_list[0]['cmd'], f'{work_dir_path}/{sample_id}/04_Unalign_BAM', 
        #     qjob_id = f'{sample_id}_unmapped_bam', threads = 10, memory = None, random_jobid=False
        # )

        mapper = Mapper(
            seq_id = sample_id,
            forward_read = read1,
            reverse_read = read2,
            bam_dir = f'{work_dir_path}/{sample_id}/05_Align_BAM',
            reference_fasta = ReferenceFasta,
            bwa_index=BwaIndex
        )
        cmd_list = mapper.mapped_bam()
        # SungridUtils.run_sungrid(
        #     'jhkim', 'all.q@ngsnode1', cmd_list[0]['cmd'], f'{work_dir_path}/{sample_id}/05_Align_BAM', 
        #     qjob_id = f'{sample_id}_mapped_bam', threads = 10, memory = None, random_jobid=False
        # )
        
        mapper = Mapper(
            seq_id = sample_id,
            forward_read = read1,
            reverse_read = read2,
            bam_dir = f'{work_dir_path}/{sample_id}/06_Merged_BAM',
            reference_fasta = ReferenceFasta,
            bwa_index=BwaIndex
        )
        cmd_list = mapper.merged_bam_files()
        # print(cmd_list)

        # SungridUtils.run_sungrid(
        #     'jhkim', 'all.q@ngsnode1', cmd_list[0]['cmd'], f'{work_dir_path}/{sample_id}/06_Merged_BAM', 
        #     qjob_id = f'{sample_id}_merged_bam', threads = 10, memory = None, random_jobid=False
        # )
    
    bam_list = glob(f'{work_dir_path}/*/06_Merged_BAM/*primary.bam')
    chrom_list = [f'chr{i}' for i in range(1,23)] +  ['chrX','chrY','chrM']
    total_cmd = []
    for bam_path in bam_list:
        
        sample_id = os.path.basename(bam_path).split('.')[0]
        if sample_id.find('25_04') > 0:        

            af_mapping = AfterMapping(
                seq_id = sample_id,
                bam_dir = f'{work_dir_path}/{sample_id}/06_Merged_BAM',
                qc_dir = f'{work_dir_path}/{sample_id}/06_Merged_BAM',
                reference_fasta = ReferenceFasta,
            )
            cmd_list = af_mapping.build_mark_duplicates_cmd()
            
            # total_cmd.append(cmd_list[0])
            # os.system(cmd_list[0]['cmd'])
            # print(cmd_list)
            # exit()
            # SungridUtils.run_sungrid(
            #     'jhkim', 'all.q@ngsnode1', cmd_list[0]['cmd'], f'{work_dir_path}/{sample_id}/06_Merged_BAM', 
            #     qjob_id = f'{sample_id}_dedup_bam', threads = 10, memory = None, random_jobid=False
            # )
        
            af_mapping = AfterMapping(
                seq_id = sample_id,
                bam_dir = f'{work_dir_path}/{sample_id}/06_Merged_BAM',
                qc_dir = f'{work_dir_path}/{sample_id}/07_Realign_BAM',
                reference_fasta = ReferenceFasta,
            )

            cmd_list = af_mapping.build_local_realign_cmds()
            # total_cmd.append(cmd_list[0])
            # # print(cmd_list)
            # SungridUtils.run_sungrid(
            #     'jhkim', 'all.q@ngsnode1', cmd_list[0]['cmd'], f'{work_dir_path}/{sample_id}/07_Realign_BAM', 
            #     qjob_id = f'{sample_id}_realign_bam', threads = 10, memory = None, random_jobid=False
            # )

            af_mapping = AfterMapping(
                seq_id = sample_id,
                bam_dir = f'{work_dir_path}/{sample_id}/07_Realign_BAM',
                qc_dir = f'{work_dir_path}/{sample_id}/08_BQSR_BAM',
                reference_fasta = ReferenceFasta,
            )

            cmd_list = af_mapping.build_bqsr_cmds()
            # total_cmd.append(cmd_list[0])
            # print(cmd_list)
            SungridUtils.run_sungrid(
                'jhkim', 'all.q@ngsnode1', cmd_list[0]['cmd'], f'{work_dir_path}/{sample_id}/08_BQSR_BAM', 
                qjob_id = f'{sample_id}_bqsr_bam', threads = 10, memory = None, random_jobid=False
            )



if __name__ == '__main__':

    main()

