# src/tasks/gatk/applybqsr/main.py
from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional
import os, shlex


from wave.core.task import Task
from wave.core.task_registry import register_task
from wave.utils.task_utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    to_sh_from_builder,
)

@register_task("[task].[func]")
class GatkApplyBQSRTask(Task):
    TYPE = "[task].[func]"

    INPUTS = {
        "bam": {"type": "path", "required": True, "desc": "BAM (post-dedup or realigned)"},
        "bqsr_table": {"type": "path", "required": True, "desc": "BQSR table"},
    }
    OUTPUTS = {
        "bam": {"type": "path", "required": False, "desc": "Recalibrated BAM"},
    }
    DEFAULTS: Dict[str, Any] = {
        "gatk_bin": "gatk",
        "image": "/storage/images/gatk-4.4.0.0.sif",
        "binds": ["/storage", "/data"],
        "singularity_bin": "singularity",
        "parallel_gc_threads": 14,
        "xmx_gb": 16,
    }

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id: Optional[str]=None) -> List[Sequence[str] | str]:
        
        bam = inputs["bam"]; table = inputs["bqsr_table"]
        out_bam = outputs.get("bam") or os.path.join(out_dir, f"{base}.recal.bam")

        gatk_bin = str(params.get("gatk_bin", "gatk"))
        xmx = int(params.get("xmx_gb", 16))
        pgc = int(params.get("parallel_gc_threads", 14))

        argv = [
            
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
            threads=int(self.threads or 1),
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id"),
            ensure_output_dir_key="dir",
        )