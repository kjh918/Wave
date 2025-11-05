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

from sklearn.preprocessing import MinMaxScaler
import pandas as pd 
## 

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

    chrom_list = [f'chr{i}' for i in range(1,23)] +  ['chrX','chrY','chrM']
    total_output = []
    for bam_path in bam_list:
        
        sample_id = os.path.basename(bam_path).split('.')[0]
    
        cmd = f'samtools bedcov {genome_interval} {bam_path} > {work_dir_path}/{sample_id}/00_QC_BAM/{sample_id}.cov_100K.txt'
        
        total_output.append(f'{work_dir_path}/{sample_id}/00_QC_BAM/{sample_id}.cov_100K.txt')
        # SungridUtils.run_sungrid(
        #     'jhkim', 'all.q@ngsnode1', cmd, f'{work_dir_path}/{sample_id}/00_QC_BAM/log', 
        #     qjob_id = f'{sample_id}_cov_100k', threads = 10, memory = None, random_jobid=False
        # )
    for output_path in total_output:

        df = pd.read_csv(output_path, sep='\t')

        # scaler = MinMaxScaler()
        # scaler.fit()


        print(df)

        exit()
    
     

if __name__ == '__main__':

    main()

