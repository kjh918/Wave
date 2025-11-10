# src/tasks/gatk/applybqsr/main.py
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

@register_task("gatk4.collectgcbiasmetrics")
class GatkApplyBQSRTask(Task):
    TYPE = "gatk4.collectgcbiasmetrics"

    INPUTS = {
        "bam": {"type": "path", "required": True, "desc": "BAM (post-dedup or realigned)"},
    }
    OUTPUTS = {
        "metrics": {"type": "path", "required": False, "desc": "Recalibrated BAM"},
        "summary": {"type": "path", "required": False, "desc": "Recalibrated BAM"},
        "pdf": {"type": "path", "required": False, "desc": "Recalibrated BAM"},
    }
    DEFAULTS: Dict[str, Any] = {
        "gatk_bin": "gatk",
        "image": "/storage/images/gatk-4.4.0.0.sif",
        "binds": ["/storage", "/data"],
        "validation_stringency": 'LENIENT',
        "is_bisulfite_sequenced": 'false',
        "assume_sorted": 'true',
        "singularity_bin": "singularity",
        "parallel_gc_threads": 14,
        "xmx_gb": 16,
    }

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id: Optional[str]=None) -> List[Sequence[str] | str]:
        
        bam = inputs["bam"]
        ref_fasta = inputs["reference"]
    
        out_metrics = outputs.get("metrics") or os.path.join(out_dir, f"{base}.gcbias.metrics.txt")
        out_chart_pdf = outputs.get("pdf") or os.path.join(out_dir, f"{base}.gcbias.pdf")
        out_summary = outputs.get("summary") or os.path.join(out_dir, f"{base}.gcbias.summary.txt")

        gatk_bin = str(params.get("gatk_bin", "gatk"))
        xmx = int(params.get("xmx_gb", 16))
        pgc = int(params.get("parallel_gc_threads", 14))
        validation_stringency = str(params.get("validation_stringency", 'LENIENT'))
        assume_sorted = bool(params.get("assume_sorted", True))
        is_bisulfite_sequenced = bool(params.get("is_bisulfite_sequenced", False))

        argv = [
            gatk_bin,
            "--java-options", f"-XX:ParallelGCThreads={pgc} -Xmx{xmx}g",
            'CollectGcBiasMetrics ',
            f'R={ref_fasta} ',
            f'I={bam} ',
            f'O={out_metrics} ',
            f'CHART={out_chart_pdf} ',
            f'S={out_summary} ',
            f'VALIDATION_STRINGENCY={validation_stringency} ',
            f'ASSUME_SORTED={"true" if assume_sorted else "false"} ',
            f'IS_BISULFITE_SEQUENCED={"true" if is_bisulfite_sequenced else "false"}'
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