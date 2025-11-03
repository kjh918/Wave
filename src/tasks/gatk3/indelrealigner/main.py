# src/tasks/gatk3/indelrealigner/main.py
from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional
import os, shlex
from src.tasks.task import Task, TaskRegistry
from src.tasks.util import ensure_dir, normalize_binds, singularity_exec_cmd, to_sh_from_builder

@TaskRegistry.register
class Gatk3IndelRealignerTask(Task):
    TYPE = "gatk3.indelrealigner"

    INPUTS = {
        "bam": {"type": "path", "required": True, "desc": "Dedup BAM"},
        "reference": {"type": "path", "required": True, "desc": "Reference FASTA"},
        "intervals": {"type": "path", "required": True, "desc": "Intervals from RealignerTargetCreator"},
        "known_indel1": {"type": "path", "required": True, "desc": "Known indel VCF 1"},
        "known_indel2": {"type": "path", "required": True, "desc": "Known indel VCF 2"},
    }
    OUTPUTS = {
        "bam": {"type": "path", "required": False, "desc": "Realigned BAM"},
        "dir": {"type": "dir", "required": False, "desc": "Output directory"},
    }
    DEFAULTS: Dict[str, Any] = {
        "image": "/storage/images/gatk-3.8-1.sif",
        "binds": ["/storage", "/data"],
        "singularity_bin": "singularity",
        "xmx_gb": 16,
        "java_bin": "java",
        "jar": "/usr/GenomeAnalysisTK.jar",
    }

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id: Optional[str]=None) -> List[Sequence[str] | str]:
        bam = inputs["bam"]; ref = inputs["reference"]; intervals = inputs["intervals"]
        k1 = inputs["known_indel1"]; k2 = inputs["known_indel2"]
        out_dir = ensure_dir(outputs.get("dir") or workdir)
        base = sample_id or os.path.splitext(os.path.basename(bam))[0]
        out_bam = outputs.get("bam") or os.path.join(out_dir, f"{base}.sorted.dedup.realign.bam")

        xmx = int(params.get("xmx_gb", 16))
        argv = [
            str(params.get("java_bin", "java")),
            f"-Xmx{xmx}g",
            "-jar", str(params.get("jar", "/usr/GenomeAnalysisTK.jar")),
            "-T", "IndelRealigner",
            "-R", ref,
            "-targetIntervals", intervals,
            "-known", k1,
            "-known", k2,
            "-I", bam,
            "-o", out_bam,
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
            threads=int(self.threads or 1),
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id"),
            ensure_output_dir_key="dir",
        )