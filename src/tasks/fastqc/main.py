# src/tasks/fastqc/fastqc.py
from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional
from pathlib import Path

from src.tasks.task import Task
from src.tasks.task_registry import register_task
from src.tasks.utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    join_argv_lines,
)

@register_task("fastqc")
class FastQCRunner(Task):
    TYPE = "fastqc"
    DEFAULTS: Dict[str, Any] = {
        "threads": 2,
        "extract": True,
        "image": None,                    # None이면 로컬 fastqc
        "binds": None,                    # list[str] | str | None
        "fastqc_bin": "fastqc",
        "singularity_bin": "singularity",
    }

    def _build_fastqc_cmd(
            self,
            inputs: List[str],
            out_dir: str,
            threads: int = 2,
            extract: bool = True,
            image: Optional[str] = None,
            binds: Optional[List[str]] = None,
            fastqc_bin: str = "fastqc",
            singularity_bin: str = "singularity",
        ) -> List[Sequence[str] | str]:
        """
        반환: argv 토큰 리스트 또는 문자열 라인(상위에서 join_argv_lines로 통일)
        """
        # 공통 argv
        argv: List[str] = [fastqc_bin]
        if extract:
            argv.append("--extract")
        argv += ["--threads", str(int(threads)), "--outdir", out_dir]
        argv += list(map(str, inputs))  # read1, (read2)

        if image:
            return [singularity_exec_cmd(
                image=str(image),
                argv=argv,
                binds=binds,
                singularity_bin=str(singularity_bin),
            )]
        else:
            return [argv]

    def to_sh(self) -> List[str]:
        r1 = self.inputs.get("read1")
        if not r1:
            raise ValueError("[fastqc] INPUT.read1 is required")
        r2 = self.inputs.get("read2")

        # outputs.dir 우선, 없으면 workdir
        out_dir = ensure_dir(self.outputs.get("dir") or self.workdir)

        p: Dict[str, Any] = {**self.DEFAULTS, **(self.params or {})}
        binds = normalize_binds(p.get("binds"))
        threads = int(self.threads or p.get("threads", 2))

        inputs = [str(r1)] + ([str(r2)] if r2 else [])

        lines = self._build_fastqc_cmd(
            inputs=inputs,
            out_dir=out_dir,
            threads=threads,
            extract=bool(p.get("extract", True)),
            image=p.get("image"),
            binds=binds,
            fastqc_bin=p.get("fastqc_bin", "fastqc"),
            singularity_bin=p.get("singularity_bin", "singularity"),
        )

        return join_argv_lines(lines)