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
from src.job_manager import JobManager
from src.fastqc import FastQCRunner
from src.fastp import FastpRunner
from src.seqkit import SeqkitRunner
# from src.align.bwa import Mapper
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

    job_manager = JobManager(
        queue_node = 'all.q@ngsnode1', 
        scheduler = "SGE"
    )

    parser = Parser()
    fastqc_runner   =   FastQCRunner()
    fastp_runner    =   FastpRunner()
    seqkit_runner   =   SeqkitRunner()

    rawdata_path = '/storage/home/jhkim/Projects/cbNIPT/GCX-MakeInsilicoUPDData-2025-10/Resources/Rawdata'
    work_dir_path = '/storage/home/jhkim/Projects/cbNIPT/GCX-MakeInsilicoUPDData-2025-10/Results'
    
    analysis_data_dict = parser(rawdata_path, f'{work_dir_path}/input_metadata.tsv', 'fastq')

    # for sample_id, reads in analysis_data_dict.items():
    #     read1, read2 = reads[0], reads[1]

    # cmd_list = fastqc_runner(Path(rawdata_path), Path(work_dir_path), 'fastqc', 10)
    # job_manager(
    #     job_name = 'fastqc',
    #     cmd_list = cmd_list,
    #     script_dir = Path(work_dir_path),
    #     check_status = False,
    #     max_procs = 50,
    #     poll_interval_sec = 15
    # )
    
    # cmd_list = fastp_runner(Path(rawdata_path), Path(work_dir_path), 'fastq_trimmed', 10)
    # job_manager(
    #     job_name = 'fastq_trimmed',
    #     cmd_list = cmd_list,
    #     script_dir = Path(work_dir_path),
    #     check_status = False,
    #     max_procs = 50,
    #     poll_interval_sec = 15
    # )


    cmd_list = seqkit_runner(
        fastq_dir = f'{Path(work_dir_path)}', 
        target_read_count = 30_000_000, 
        seqkit_dir = Path(work_dir_path),
        job_name = 'fastq_downsampled',
        threads = 1)
    job_manager(
        job_name = 'fastq_downsampled',
        cmd_list = cmd_list,
        script_dir = Path(work_dir_path),
        check_status = False,
        max_procs = 3,
        poll_interval_sec = 15
    )

if __name__ == '__main__':

    main()

