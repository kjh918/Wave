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

@register_task("trtools.prancstr")
class PrancsSTRTask(Task):
    TYPE = "trtools.prancstr"

    INPUTS = {
        "vcf": {"type": "path", "required": True, "desc": "VCF (Use only HipSTR)"}
    }
    OUTPUTS = {
        "table": {"type": "path", "required": False, "desc": " tab-delimited file with one row summarizing evidence of mosaicism for each call analyzed"},
    }
    DEFAULTS: Dict[str, Any] = {
        "trtools": "/storage/apps/trtools-0.7/trtools",
        "region": "/storage/home/kangsm/myDB/STR_references/hg38.annotated.markers.trtools.bed",
        "image": "",
        "binds": ["/storage", "/data"],
        "stutter_in": False,
        "singularity_bin": "singularity"
    }

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id: Optional[str]=None) -> List[Sequence[str] | str]:
        
        ## INPUT ##
        bams = inputs["bams"]
        fasta = inputs["fasta"]
        region = inputs["region"]
        
        ## OUTPUT ##
        out_str_vcf = outputs.get("vcf") or os.path.join(out_dir, f"{base}.recal.bam")

        ## PARAMS ##
        trtools = str(params.get("trtools", "trtools"))
        
        
        xmx = int(params.get("xmx_gb", 16))
        pgc = int(params.get("parallel_gc_threads", 14))
        image = params.get("image")

        argv = [
            trtools, 
            '--bams',f'{bams}',
            '--fasta',f'{fasta}',
            '--bed',f'{region}',
            '--str-vcf',f'{out_str_vcf}'
        ]

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