# src/tasks/gatk/baserecalibrator/main.py
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

@register_task("gatk4.baserecalibrator")
class GatkBaseRecalibratorTask(Task):
    TYPE = "gatk4.baserecalibrator"

    INPUTS = {
        "bam": {"type": "path", "required": True, "desc": "Input BAM"},
    }
    OUTPUTS = {
        "table": {"type": "path", "required": False, "desc": "BQSR table"}
    }
    DEFAULTS: Dict[str, Any] = {
        "gatk_bin": "gatk",
        "known_snp": "gatk",
        "known_indel1": "gatk",
        "known_indel2": "gatk",
        "image": "/storage/images/gatk-4.4.0.0.sif",
        "binds": ["/storage", "/data"],
        "singularity_bin": "singularity",
        "parallel_gc_threads": 14,
        "xmx_gb": 16,
    }

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id: Optional[str]=None) -> List[Sequence[str] | str]:
        bam = inputs["bam"]
        ref = params["reference"]
        
        ks_list = [
            params["known_snp"],
            params["known_indel1"],
            params["known_indel2"],
        ]
        out_dir = ensure_dir(outputs.get("dir") or workdir)
        base = sample_id or os.path.splitext(os.path.basename(bam))[0]
        out_tbl = outputs.get("table") or os.path.join(out_dir, f"{base}.recal.table.txt")

        gatk_bin = str(params.get("gatk_bin", "gatk"))
        xmx = int(params.get("xmx_gb", 16))
        pgc = int(params.get("parallel_gc_threads", 14))

        argv = [
            gatk_bin, "BaseRecalibrator",
            "--java-options", f"-XX:ParallelGCThreads={pgc} -Xmx{xmx}g",
            "--input", bam,
            "--reference", ref,
            "--output", out_tbl,
        ]
        for vcf in ks_list:
            argv += ["--known-sites", str(vcf)]

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