# src/tasks/gatk/markduplicates/main.py
from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional
import os, shlex

from src.tasks.task import Task
from src.tasks.task_registry import register_task
from src.tasks.utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    to_sh_from_builder,
)

@register_task("gatk4.markduplicates")
class GatkMarkDuplicatesTask(Task):
    TYPE = "gatk.markduplicates"

    INPUTS = {
        "bam": {"type": "path", "required": True, "desc": "Sorted BAM input"},
    }
    OUTPUTS = {
        "bam": {"type": "path", "required": False, "desc": "Dedup BAM output"},
        "metrics": {"type": "path", "required": False, "desc": "Metrics file"},
        "dir": {"type": "dir", "required": False, "desc": "Output directory (default: workdir)"},
    }
    DEFAULTS: Dict[str, Any] = {
        "gatk_bin": "gatk",
        "image": "/storage/images/gatk-4.4.0.0.sif",
        "binds": ["/storage", "/data"],
        "singularity_bin": "singularity",
        "parallel_gc_threads": 14,
        "xmx_gb": 16,
        "create_index": True,
        "remove_sequencing_duplicates": True,
    }

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id: Optional[str]=None) -> List[Sequence[str] | str]:
        in_bam = inputs["bam"]
        out_dir = ensure_dir(outputs.get("dir") or workdir)
        base = sample_id or os.path.splitext(os.path.basename(in_bam))[0]
        out_bam = outputs.get("bam") or os.path.join(out_dir, f"{base}.sorted.dedup.bam")
        metrics = outputs.get("metrics") or os.path.join(out_dir, f"{base}.mark.duplicates.metrics.txt")

        gatk_bin = str(params.get("gatk_bin", "gatk"))
        xmx = int(params.get("xmx_gb", 16))
        pgc = int(params.get("parallel_gc_threads", 14))
        create_index = str(params.get("create_index", True)).lower()
        remove_seq_dup = str(params.get("remove_sequencing_duplicates", True)).lower()

        argv = [
            gatk_bin, "MarkDuplicates",
            "--java-options", f"-XX:ParallelGCThreads={pgc} -Xmx{ xmx }g",
            "--INPUT", in_bam,
            "--OUTPUT", out_bam,
            "--METRICS_FILE", metrics,
            "--CREATE_INDEX", create_index,
            "--REMOVE_SEQUENCING_DUPLICATES", remove_seq_dup,
        ]

        image = params.get("image")
        if image:
            cmd = singularity_exec_cmd(
                image=str(image),
                argv=argv,
                binds=normalize_binds(params.get("binds")),
                singularity_bin=str(params.get("singularity_bin", "singularity")),
            )
            return [" ".join(map(shlex.quote, cmd))]
        else:
            return [" ".join(map(shlex.quote, argv))]

    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"dir": str(self.workdir)},
            params=p,
            threads=int(self.threads or p.get("threads", 1)),
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id"),
            ensure_output_dir_key="dir",
        )