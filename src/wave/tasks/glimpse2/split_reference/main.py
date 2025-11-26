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

@register_task("glimpse2.split_reference")
class Glimpse2SpiltReferenceTask(Task):
    TYPE = "glimpse2.split_reference"

    INPUTS = {
        "vcf": {"type": "path", "required": True, "desc": "VCF or BCF (Use only HipSTR)"}
    }
    OUTPUTS = {
        "table": {"type": "path", "required": False, "desc": " tab-delimited file with one row summarizing evidence of mosaicism for each call analyzed"},
    }
    DEFAULTS: Dict[str, Any] = {
        "glimpse2_split_reference": "/bin/GLIMPSE2_split_reference",
        "region": "/storage/home/kangsm/myDB/STR_references/hg38.annotated.markers.trtools.bed",
        "image": "/storage/home/jhkim/Apps/GLIMPSE/glimpse2.sif",
        "binds": ["/storage", "/data"],
        "stutter_in": False,
        "singularity_bin": "singularity"
    }

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id: Optional[str]=None) -> List[Sequence[str] | str]:
        
        ## INPUT ##
        map_file = inputs["map_file"]
        reference = inputs["reference"]
        chunks_file = inputs["chunk"]
        
        ## OUTPUT ##
        output_reference = outputs.get("reference")
        ## PARAMS ##
        glimpse2_reference = str(params.get("glimpse2_split_reference", "/bin/GLIMPSE2_split_reference"))
        
        image = params.get("image")
        
        argv = [f"""
REF={reference}
MAP={map_file}

while IFS="" read -r LINE || [ -n "$LINE" ];
do
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    ID=$(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
"""]
        if image: 
            cmd = singularity_exec_cmd(
                image=str(image),
                argv=[],
                binds=normalize_binds(params.get("binds")),
                singularity_bin=str(params.get("singularity_bin", "singularity")),
            )
            argv.append(' '.join(cmd) + f""" {glimpse2_reference} \\
        --reference $REF \\
        --map $MAP \\
        --input-region $IRG \\
        --output-region $ORG \\
        --output {output_reference} \\
        --threads {threads}
done < {chunks_file}
""")        
            return [" ".join(argv)]
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