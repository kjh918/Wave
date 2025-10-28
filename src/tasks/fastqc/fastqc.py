# src/tasks/fastqc/fastqc.py
from __future__ import annotations
from typing import Dict, Any, List
from pathlib import Path

from src.tasks.task import Task, TaskRegistry
from ._func import build_fastqc_cmd

@TaskRegistry.register   # ✅ 등록 필수
class FastQCRunner(Task):
    TYPE = "fastqc"
    DEFAULTS: Dict[str, Any] = {
        "threads": 2,
        "extract": True,
        "image": "/storage/images/fastqc-0.12.1.sif",
        "binds": ["/storage", "/data"],
        "fastqc_bin": "fastqc",
        "singularity_bin": "singularity",
    }

    def to_sh(self) -> List[str]:
        inputs = [self.inputs.get("fastq_r1")]
        if self.inputs.get("fastq_r2"):
            inputs.append(self.inputs["fastq_r2"])

        out_dir = self.outputs.get("dir") or str(self.workdir)
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        p = self.params
        return build_fastqc_cmd(
            inputs=inputs,
            out_dir=out_dir,
            threads=int(p["threads"]),
            extract=bool(p["extract"]),
            image=p.get("image"),
            binds=p.get("binds"),
            fastqc_bin=p.get("fastqc_bin", "fastqc"),
            singularity_bin=p.get("singularity_bin", "singularity"),
        )