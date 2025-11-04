import os, sys, subprocess, string, random, time
from pathlib import Path


from src.utils.flags import skip_if_done, flag_on_complete


class Executor:
    """Base Executor — 공통 기능 정의"""

    def __init__(self, logdir: str | Path = "./qlog"):
        self.logdir = Path(logdir)
        self.logdir.mkdir(parents=True, exist_ok=True)

    def _make_script(self, cmd: str, job_id: str) -> Path:
        """명령어를 실행 가능한 bash 스크립트로 저장"""
        script_path = self.logdir / f"{job_id}.sh"
        with open(script_path, "w") as f:
            f.write(f"#!/usr/bin/env bash\nset -euo pipefail\n{cmd}\n")
        script_path.chmod(0o755)
        return script_path


# ---------------------------------------------------------------------
# ✅ 1. BashExecutor — 로컬에서 직접 실행
# ---------------------------------------------------------------------
class BashExecutor(Executor):
    """로컬 bash 실행용 executor"""

    def run(self, cmd: str, job_id: str | None = None) -> int:
        job_id = job_id or "".join(random.choice(string.ascii_letters) for _ in range(10))
        script_path = self._make_script(cmd, job_id)

        stdout_path = self.logdir / f"{job_id}.stdout"
        stderr_path = self.logdir / f"{job_id}.stderr"

        with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
            process = subprocess.Popen(["bash", str(script_path)], stdout=out, stderr=err)
            ret = process.wait()

        if ret != 0:
            raise RuntimeError(f"BashExecutor: job {job_id} failed with code {ret}")
        return ret


# ---------------------------------------------------------------------
# ✅ 2. SunGridExecutor — SGE(qsub) 기반 실행
# ---------------------------------------------------------------------
class SunGridExecutor(Executor):
    """SGE (qsub) 기반 실행기"""

    @staticmethod
    def check_node_info(node=None):
        """사용 가능한 qsub queue/node 확인"""
        nodeList = []
        sqlList = subprocess.check_output('qconf -sql', shell=True).split()
        selList = subprocess.check_output('qconf -sel', shell=True).split()
        for sql in sqlList:
            sql = sql.decode()
            for sel in selList:
                sel = sel.decode()
                nodeList.append(f"{sql}@{sel}")
        return node in nodeList

    def qsub_sh(self, node: str, script_path: str, threads: int = 1, memory: int | None = None,
            job_id: str | None = None, random_jobid: bool = False) -> str:
        """SGE qsub 제출"""
        job_id = (
            "".join(random.choice(string.ascii_letters) for _ in range(10))
            if random_jobid else job_id
        )

        stdout_path = script_path.replace(".sh",".stdout")
        stderr_path = script_path.replace(".sh",".stderr")

        if memory is None:
            qcmd = (
                f"qsub -N {job_id} -q {node} "
                f"-o {stdout_path} -e {stderr_path} "
                f"-pe smp {threads} -V -cwd {script_path}"
            )
        else:
            qcmd = (
                f"qsub -N {job_id} -q {node} "
                f"-o {stdout_path} -e {stderr_path} "
                f"-l h_vmem={memory}G -pe smp {threads} -V -cwd {script_path}"
            )

        x = subprocess.check_output(qcmd, shell=True, encoding='utf-8')
        qsub_id = x.split()[2]
        return qsub_id
    
    @skip_if_done(flag_name=".done", require_outputs_ok=True)  # 실행 전 스킵
    @flag_on_complete(flag_name=".done", fail_flag=".failed")  # 실행 후 플래그
    def run(self, node: str, cmd: str, threads: int = 1, memory: int | None = None,
            job_id: str | None = None, random_jobid: bool = False) -> str:
        """SGE qsub 제출"""
        job_id = (
            "".join(random.choice(string.ascii_letters) for _ in range(10))
            if random_jobid else job_id
        )
        script_path = self._make_script(cmd, job_id)
        stdout_path = self.logdir / f"{job_id}.stdout"
        stderr_path = self.logdir / f"{job_id}.stderr"

        if memory is None:
            qcmd = (
                f"qsub -N {job_id} -q {node} "
                f"-o {stdout_path} -e {stderr_path} "
                f"-pe smp {threads} -V -cwd {script_path}"
            )
        else:
            qcmd = (
                f"qsub -N {job_id} -q {node} "
                f"-o {stdout_path} -e {stderr_path} "
                f"-l h_vmem={memory}G -pe smp {threads} -V -cwd {script_path}"
            )

        x = subprocess.check_output(qcmd, shell=True, encoding='utf-8')
        qsub_id = x.split()[2]
        return qsub_id

    @staticmethod
    def check_completed_job(jobid_list, qstat_path: str = "./qstat.tmp"):
        """SGE에서 jobid_list의 모든 작업이 완료될 때까지 대기"""
        while True:
            os.system(f"qstat > {qstat_path}")
            running_jobs = []
            with open(qstat_path, "r") as f:
                for line in f:
                    parts = line.strip().split()
                    if not parts or parts[0].startswith(("job", "-")):
                        continue
                    running_jobs.append(parts[0])

            if all(j not in running_jobs for j in jobid_list):
                break
            time.sleep(60)
        os.remove(qstat_path)

    def run_workflow(task_dict):

        for key, value in task_dict:
            
            pass
            