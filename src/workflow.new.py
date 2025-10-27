
import os
import sys
import argparse
import subprocess
from glob import glob
from typing import Dict, Any, List, Optional
import pyaml 
from pathlib import Path 

sys.path.append('/storage/home/jhkim/scripts/Wave')
from src.workflow import Workflow 
from src.executor import SunGridExecutor 


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="run",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--config", dest="config", help="config.yaml")
    return p

def main():
    parser = build_argparser()
    args = parser.parse_args()
    workflow = Workflow(Path(args.config))

    sample_dict = workflow.discover_samples()
    tasks       = workflow.build_task_list()

    sge = SunGridExecutor   (
        logdir=str(Path(workflow.work_dir))
    )
    print(tasks)
    for sample_id, params in sample_dict.items():
        print(sample_id)
        print(params)
        print(workflow)

    #     for task in tasks:

    #         workflow.append_task(
    #             sample_id, 
    #             task.name,
    #             task.params
    #             )

    # workflow.save_workflow(
    #     logdir=str(Path(workflow.work_dir))
    #     )

    # sge.run(workflow.workflow)






        # r1 = params["read"].get("1")
        # r2 = params["read"].get("2")

        # # 3️⃣ FastQC Task 인스턴스 생성
        # task = FastQCTask(
        #     name='fastqc',
        #     workdir=str(Path(workflow.work_dir)),
        #     params={
        #         "RawFastqDir": str(Path(r1).parent),
        #         "SeqID": sample_id,
        #         "threads": 4,
        #         "extract": True,
        #         "image": "/storage/images/fastqc-0.12.1.sif",   # optional
        #         "binds": ["/storage", "/data"],
        #         "qcResDir":str(Path(workflow.work_dir) / f'{sample_id}/00_Rawdata_Fastqc')
        #     }
        # )
        
        # cmd_lines = list(task.to_sh())
        # cmd_script = "\n".join(cmd_lines)
        # sge = SunGridExecutor   (
        #     logdir=str(Path(workflow.work_dir) / f'{sample_id}/00_Rawdata_Fastqc')
        # )
        # # sge.run(
        # #     node="all.q@ngsnode1",
        # #     job_id=f"{sample_id}_fastqc",
        # #     cmd=cmd_script,
        # #     threads=4
        # # )

        
        # task = FastpTask(
        #     name='fastp',
        #     workdir=str(Path(workflow.work_dir)),
        #     params={
        #         "RawFastqDir": str(Path(r1).parent),
        #         "SeqID": sample_id,
        #         "threads": 4,
        #         "PicoplexGold": True,
        #         "image": "/storage/images/fastp-0.23.4.sif",   # optional
        #         "binds": ["/storage", "/data"],
        #         "TrimFastqDir":str(Path(workflow.work_dir) / f'{sample_id}/01_Trim_Fastp')
        #     }
        # )  

        # # 4️⃣ 명령어 생성
        # cmd_lines = list(task.to_sh())
        # cmd_script = "\n".join(cmd_lines)
        # sge = SunGridExecutor   (
        #     logdir=str(Path(workflow.work_dir) / f'{sample_id}/01_Trim_Fastp')
        # )
        
        # # sge.run(
        # #     node="all.q@ngsnode1",
        # #     job_id=f"{sample_id}_fastp",
        # #     cmd=cmd_script,
        # #     threads=4
        # # )

        # task = JellyfishTask(
        #     name='jellyfish',
        #     workdir=str(Path(workflow.work_dir)),
        #     params={
        #         "RawFastqDir": str(Path(workflow.work_dir) / f'{sample_id}/01_Trim_Fastp'),
        #         "SeqID": sample_id,
        #         "threads": 4,
        #         "outDir": str(Path(workflow.work_dir) / f'{sample_id}/02_Count_k-mer_JellyFish'),
        #     }
        # )  

        # # 4️⃣ 명령어 생성
        # # cmd_lines = list(task.to_sh())
        # # cmd_script = "\n".join(cmd_lines)
        # # sge = SunGridExecutor   (
        # #     logdir=str(Path(workflow.work_dir) / f'{sample_id}/02_Count_k-mer_JellyFish')
        # # )

        # cmd_lines = list(task.to_sh())
        # # cmd_script = "\n".join(cmd_lines)
        # # sge = SunGridExecutor   (
        # #     logdir=str(Path(workflow.work_dir) / f'{sample_id}/02_Count_k-mer_JellyFish')
        # # )
        
        
        # # print(cmd_script)
        # # sge.run(
        # #     node="all.q@ngsnode1",
        # #     job_id=f"{sample_id}_jellyfish_histo",
        # #     cmd=cmd_script,
        # #     threads=8
        # # )




if __name__ == '__main__':

    main()

