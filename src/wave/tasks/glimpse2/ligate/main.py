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

@register_task("glimpse2.chunk")
class Glimpse2ChunkTask(Task):
    TYPE = "glimpse2.chunk"

    INPUTS = {
        "vcf": {"type": "path", "required": True, "desc": "VCF or BCF (Use only HipSTR)"}
    }
    OUTPUTS = {
        "table": {"type": "path", "required": False, "desc": " tab-delimited file with one row summarizing evidence of mosaicism for each call analyzed"},
    }
    DEFAULTS: Dict[str, Any] = {
        "glimpse2_chunk": "/bin/GLIMPSE2_chunk",
        "region": "/storage/home/kangsm/myDB/STR_references/hg38.annotated.markers.trtools.bed",
        "image": "/storage/home/jhkim/Apps/GLIMPSE/glimpse2.sif",
        "binds": ["/storage", "/data"],
        "stutter_in": False,
        "singularity_bin": "singularity"
    }

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id: Optional[str]=None) -> List[Sequence[str] | str]:
        
        ## INPUT ##
        vcf = inputs["vcf"]
        map_file = inputs["map_file"]
        region = inputs["region"]
        
        ## OUTPUT ##
        output_txt = outputs.get("txt")

        ## PARAMS ##
        glimpse2_chunk = str(params.get("trtools", "trtools"))
        
        image = params.get("image")

        argv = [
            glimpse2_chunk, 
            f"--input",f"{vcf}",
            f"--map",f"{map_file}",
            f"--region",f"{region}",
            f"--output",f"{output_txt}",
            f"-T",f"{threads}",
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