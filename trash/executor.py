import os, sys, optparse, time, subprocess, string, random, shutil, configparser, argparse, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from glob import glob
import json

class SungridUtils():

    def __init__(self) -> None:
        pass

    def check_node_info(node=None):
        nodeList = []
        sqlList = subprocess.check_output('qconf -sql',shell=True).split()
        selList = subproceqsss.check_output('qconf -sel',shell=True).split()
        for sql in sqlList :
            sql = sql.decode()
            for sel in selList :
                sel = sel.decode()
                nodeList.append("%s@%s"%(sql,sel))
        if node in nodeList :
            #logFileWrite(logFilePath, 'This pipeline runs on %s.\n'%(NODE))
            qnode = node
        #else :
        #    logFileWrite(logFilePath, 'No node was specified. SGE automatically allocates nodes.\n')
    
    def run_sungrid(user, node, cmd, qlog_path, qjob_id = None, threads = 1, memory = None, random_jobid=False):

        if random_jobid==True:
            qjob_id  =   "".join([random.choice(string.ascii_letters) for _ in range(10)])
        else:
            qjob_id  =   qjob_id
        os.makedirs(qlog_path, exist_ok=True)
        script_path     =   f'{qlog_path}/{qjob_id}.sh'
        script_out_log  =   f'{qlog_path}/{qjob_id}.stdout'
        script_err_log  =   f'{qlog_path}/{qjob_id}.stderr'

        with open(script_path, 'w') as handle:
            handle.write(f"#!/usr/bin/bash\necho -n \"$JOB_ID\"\n{cmd}\n")
        if memory == None:
            qcmd = f'qsub -N {qjob_id} -q {node} -o {script_out_log} -e {script_err_log} -pe smp {threads} -V -cwd {script_path}'
        else:
            qcmd = f'qsub -N {qjob_id} -q {node} -o {script_out_log} -e {script_err_log} -l h_vmem={memory}G -pe smp {threads} -V -cwd {script_path}'
        x = subprocess.check_output(qcmd, shell=True, encoding ='utf-8')
        job_id = x.split(' ')[2]
        return job_id
    
    def check_completed_job(jobid_list, qstat_path):
        job_count  =   len(jobid_list)
        continue_job_count =  0
        while(True):
            os.system(f'qstat > {qstat_path}')
            running_job = []
            with open(qstat_path,'r') as handle:
                for line in handle:
                    run_jobid   =   line.strip('\n').split(' ')[0]
                    if run_jobid.startswith('-'):
                        continue
                    elif run_jobid.startswith('job'):
                        continue
                    else:
                        running_job.append(run_jobid)
            if len(running_job) == 0:
                break
            else:
                for r_job in running_job:
                    if r_job in jobid_list:
                        continue_job_count += 1
            if continue_job_count == 0:
                break
            else:
                continue_job_count = 0
            time.sleep(60)
        os.system(f'rm {qstat_path}')