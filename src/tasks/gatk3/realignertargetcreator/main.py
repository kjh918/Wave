# src/tasks/gatk3/realignertargetcreator/main.py
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

@register_task("gatk3.realignertargetcreator")
class Gatk3RealignerTargetCreatorTask(Task):
    TYPE = "gatk3.realignertargetcreator"

    INPUTS = {
        "bam": {"type": "path", "required": True, "desc": "Dedup BAM"},
        "reference": {"type": "path", "required": True, "desc": "Reference FASTA"},
        "known_indel1": {"type": "path", "required": True, "desc": "Known indel VCF 1"},
        "known_indel2": {"type": "path", "required": True, "desc": "Known indel VCF 2"},
    }
    OUTPUTS = {
        "intervals": {"type": "path", "required": False, "desc": "Intervals output"},
        "dir": {"type": "dir", "required": False, "desc": "Output directory"},
    }
    DEFAULTS: Dict[str, Any] = {
        "image": "/storage/images/gatk-3.8-1.sif",
        "binds": ["/storage", "/data"],
        "singularity_bin": "singularity",
        "xmx_gb": 16,
        "java_bin": "java",
        "jar": "/usr/GenomeAnalysisTK.jar",
        "nt": 8,
    }

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id: Optional[str]=None) -> List[Sequence[str] | str]:
        bam = inputs["bam"]
        ref = params["reference"]
        k1 = params["known_indel1"]; k2 = params["known_indel2"]
        out_dir = ensure_dir(outputs.get("dir") or workdir)
        base = sample_id or os.path.splitext(os.path.basename(bam))[0]
        intervals = outputs.get("intervals") or os.path.join(out_dir, f"{base}.realign.target.intervals")

        xmx = int(params.get("xmx_gb", 16))
        nt = int(params.get("nt", threads or 8))
        argv = [
            str(params.get("java_bin", "java")),
            f"-Xmx{xmx}g",
            "-jar", str(params.get("jar", "/usr/GenomeAnalysisTK.jar")),
            "-T", "RealignerTargetCreator",
            "-R", ref,
            "-I", bam,
            "-o", intervals,
            "-known", k1,
            "-known", k2,
            "-nt", str(nt),
        ]

        cmd = singularity_exec_cmd(
            image=str(params.get("image")),
            argv=argv,
            binds=normalize_binds(params.get("binds")),
            singularity_bin=str(params.get("singularity_bin", "singularity")),
        )
        return [" ".join(map(shlex.quote, cmd))]

    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"dir": str(self.workdir)},
            params=p,
            threads=int(self.threads or p.get("nt", 8)),
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id"),
            ensure_output_dir_key="dir",
        )