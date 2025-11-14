import os, sys, subprocess, string, random, time
from pathlib import Path
from multiprocessing import Pool, get_context
from src.utils.flags import skip_if_done, flag_on_complete


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

    def _make_script(self, cmd: str, job_id: str) -> Path:
        """명령어를 실행 가능한 bash 스크립트로 저장"""
        script_path = self.logdir / f"{job_id}.sh"
        with open(script_path, "w") as f:
            f.write(f"#!/usr/bin/env bash\nset -euo pipefail\n{cmd}\n")
        script_path.chmod(0o755)
        return script_path

    def make_script(self, *, cmd: str, job_id: str, workdir: Path, outputs: list[str]) -> Path:
        # (기존 그대로) .done/.failed 생성까지 포함
        script_path = self.session_dir / f"{job_id}.sh"
        q_outputs = " ".join(shlex.quote(p) for p in outputs) if outputs else ""
        script_path.write_text(f"""#!/usr/bin/env bash
set -euo pipefail
cd {shlex.quote(str(workdir))}
{cmd}
missing=0
files=({q_outputs})
if [ ${{#files[@]}} -gt 0 ]; then
  for f in "${{files[@]}}"; do
    if [ ! -s "$f" ]; then missing=1; fi
  done
fi
if [ $missing -eq 0 ]; then echo OK > .done; else echo FAILED > .failed; exit 1; fi
""")
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
# src/executor.py (발췌/추가)
from pathlib import Path
import os, subprocess, shlex, time

class SunGridExecutor:
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


    def qsub_sh(self, *, node: str, script_path: str, threads: int, job_id: str,
                stdout_path: Path, stderr_path: Path, hold_jid: str | None = None, memory_gb: int | None = None) -> str:
        # (제출 전 용량확인 로직 유지)
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        stderr_path.parent.mkdir(parents=True, exist_ok=True)
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
        # (활성잡 등록 로직 유지)
        return jid
    
    
# class SunGridExecutor(Executor):
#     """SGE (qsub) 기반 실행기"""

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
        # (나머지 필드: 스레드/잡 제한 로직은 기존 유지)

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
    
    @skip_if_done(flag_name=".done", require_outputs_ok=True)  # 실행 전 스킵
    @flag_on_complete(flag_name=".done", fail_flag=".failed")  # 실행 후 플래그
    def qsub_sh(
            self,
            node: str,
            script_path: str,
            threads: int = 1,
            memory: int | None = None,
            job_id: str | None = None,
            random_jobid: bool = False,
            hold_jid: str | list[str] | None = None,   # ⬅️ 추가
        ) -> str:
        """SGE qsub 제출 (+ 의존성 체인 지원)"""
        job_id = (
            "".join(random.choice(string.ascii_letters) for _ in range(10))
            if random_jobid else job_id
        )

        stdout_path = script_path.replace(".sh",".stdout")
        stderr_path = script_path.replace(".sh",".stderr")

        hold_opt = ""
        if hold_jid:
            if isinstance(hold_jid, (list, tuple)):
                hold_opt = f"-hold_jid {','.join(hold_jid)} "
            else:
                hold_opt = f"-hold_jid {hold_jid} "

        if memory is None:
            qcmd = (
                f"qsub -N {job_id} -q {node} "
                f"{hold_opt}"
                f"-o {stdout_path} -e {stderr_path} "
                f"-pe smp {threads} -V -cwd {script_path}"
            )
        else:
            qcmd = (
                f"qsub -N {job_id} -q {node} "
                f"{hold_opt}"
                f"-o {stdout_path} -e {stderr_path} "
                f"-l h_vmem={memory}G -pe smp {threads} -V -cwd {script_path}"
            )

        x = subprocess.check_output(qcmd, shell=True, encoding='utf-8')
        qsub_id = x.split()[2]
        # qsub_id = ''
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

            