import os, sys, subprocess, string, random, time
from pathlib import Path
from multiprocessing import Pool, get_context
from wave.utils.flags import auto_flag_on_complete


def run_local_shell(*, cmd: str, task, cwd=None) -> bool:
    import subprocess
    rc = subprocess.run(cmd, shell=True, cwd=cwd).returncode
    if rc != 0:
        raise RuntimeError(f"failed: {cmd}")
    return True

class Executor:
    """Base Executor — 공통 기능 정의"""

    def __init__(self, logdir: str | Path = "./qlog"):
        self.logdir = Path(logdir)
        self.logdir.mkdir(parents=True, exist_ok=True)

    def make_script(self, cmd: str, job_id: str) -> Path:
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

    @auto_flag_on_complete(outputs_arg="outputs", workdir_arg="workdir")
    def run(self, cmd: str, job_id: str | None = None, *, task=None, workdir: str | Path = None, outputs: dict | None = None) -> int:
        """
        - task: optional Task 객체 (task.workdir, task.outputs 사용)
        - workdir/outputs: explicit fallback (if task is None)
        """
        # prefer task
        if task is not None:
            self.current_task = task  # so decorators can pick it up if needed
            workdir = getattr(task, "workdir", workdir)
            outputs = getattr(task, "outputs", outputs)

        job_id = job_id or "".join(random.choice(string.ascii_letters) for _ in range(10))
        script_path = self.make_script(cmd, job_id)  # make_script should accept optional workdir
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
# src/executor.py (발췌/추가)
from pathlib import Path
import os, subprocess, shlex, time

class SunGridExecutor(Executor):
    def __init__(
            self,
            logdir: Path,                      # ✅ 중앙 로그 루트
            *,
            run_id: str = None,                  # ✅ 세션 ID (없으면 자동)
            user: str | None = None,
            max_threads_total: int | None = None,
            max_concurrent_jobs: int | None = None,
            poll_sec: int = 15,
            name_prefix: str | None = None,
            link_in_taskdir: bool = True,        # ✅ taskdir에 심볼릭 링크 생성
            layout: str = "by-sample",           # ✅ flat|by-sample|by-sample-task
        ):
        self.user = user or os.getenv("USER", "unknown")
        self.logdir = Path(logdir)
        self.run_id = run_id or time.strftime("%Y%m%d_%H%M%S")
        self.session_dir = self.logdir / self.run_id
        self.session_dir.mkdir(parents=True, exist_ok=True)
        self.link_in_taskdir = link_in_taskdir
        self.layout = layout

    def _layout_dir(self, sid: str | None, order: int | None, taskname: str) -> Path:
        if self.layout == "flat" or sid is None:
            return self.session_dir
        if self.layout == "by-sample":
            d = self.session_dir / sid
        else:  # by-sample-task
            d = self.session_dir / sid / f"{order:02d}_{taskname}" if order is not None else self.session_dir / sid / taskname
        d.mkdir(parents=True, exist_ok=True)
        return d

    def log_paths(self, *, sid: str | None, order: int | None, taskname: str):
        d = self._layout_dir(sid, order, taskname)
        base = f"{sid + '_' if sid else ''}{(order and f'{order:02d}-') or ''}{taskname}"
        return (d / f"{base}.stdout", d / f"{base}.stderr")

    def _link_into_taskdir(self, taskdir: Path, stdout_p: Path, stderr_p: Path):
        if not self.link_in_taskdir:
            return
        try:
            (taskdir / "stdout").unlink(missing_ok=True)
            (taskdir / "stderr").unlink(missing_ok=True)
            (taskdir / "stdout").symlink_to(stdout_p)
            (taskdir / "stderr").symlink_to(stderr_p)
        except Exception:
            pass

    @auto_flag_on_complete(outputs_arg="outputs", workdir_arg="workdir")
    def qsub_sh(self, *, node: str, script_path: str, threads: int, job_id: str,
                memory_gb: int | None = None,
                hold_jid: str | list[str] | None = None,   # ⬅️ 추가
                finish_and_run: bool = True) -> str:
        # (제출 전 용량확인 로직 유지)
        
        stdout_path = script_path.replace(".sh",".stdout")
        stderr_path = script_path.replace(".sh",".stderr")
        qstat_path  = script_path.replace(".sh",".qstat")

        Path(stdout_path).parent.mkdir(parents=True, exist_ok=True)
        Path(stderr_path).parent.mkdir(parents=True, exist_ok=True)

        args = [
            "qsub", "-N", job_id, "-q", node,
            "-o", str(stdout_path), "-e", str(stderr_path),
            "-pe", "smp", str(int(threads)), "-V", "-cwd"
        ]
        if memory_gb is not None: args += ["-l", f"h_vmem={int(memory_gb)}G"]
        if hold_jid: args += ["-hold_jid", str(hold_jid)]
        args.append(script_path)
        out = subprocess.check_output(args, text=True)
        jid = next((p for p in reversed(out.split()) if p.isdigit()), None) or out.split()[-1]

        if finish_and_run:
            self.check_completed_job([jid], qstat_path = qstat_path)
            
        return jid
    
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
    