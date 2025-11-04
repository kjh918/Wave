from __future__ import annotations
from typing import Dict, Any, List, Sequence
from pathlib import Path
import os, shlex

from src.tasks.task import Task
from src.tasks.task_registry import register_task
from src.tasks.utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    to_sh_from_builder,
)

@register_task("picard.mergebamalignment")
class PicardMergeBamAlignmentTask(Task):
    """
    Picard MergeBamAlignment Task

    Example:
        java -XX:ParallelGCThreads=14 -Xmx16384m -jar /storage/apps/bin/picard.jar MergeBamAlignment \
            --UNMAPPED_BAM ${BamDir}/${SeqID}.fastqtosam.bam \
            --ALIGNED_BAM ${BamDir}/${SeqID}.bwa.mem.sam \
            --REFERENCE_SEQUENCE ${ReferenceFasta} \
            --OUTPUT ${BamDir}/${SeqID}.primary.bam \
            --CREATE_INDEX true --MAX_INSERTIONS_OR_DELETIONS -1 \
            --CLIP_ADAPTERS false --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
            --ATTRIBUTES_TO_RETAIN XS \
            --EXPECTED_ORIENTATIONS FR --EXPECTED_ORIENTATIONS RF
    """

    TYPE = "picard.mergebamalignment"

    INPUTS = {
        "unmapped_bam": {"type": "path", "required": True, "desc": "UNMAPPED_BAM"},
        "aligned_bam": {"type": "path", "required": True, "desc": "ALIGNED_BAM (SAM from BWA MEM)"},
        "reference": {"type": "path", "required": True, "desc": "REFERENCE FASTA"},
    }
    OUTPUTS = {
        "bam": {"type": "path", "required": False, "desc": "Merged BAM"},
        "dir": {"type": "dir", "required": False, "desc": "Output directory"},
    }

    DEFAULTS = {
        "java_bin": "java",
        "jar": "/storage/apps/bin/picard.jar",
        "xmx_gb": 16,
        "parallel_gc_threads": 14,
        "create_index": True,
        "max_indel": -1,
        "clip_adapters": False,
        "primary_strategy": "MostDistant",
        "retain_attributes": ["XS"],
        "expected_orientations": ["FR", "RF"],
        "image": None,
        "binds": None,
        "singularity_bin": "singularity",
    }

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id=None) -> List[Sequence[str] | str]:
        unmapped = inputs["unmapped_bam"]
        aligned = inputs["aligned_bam"]
        ref = inputs["reference"]

        out_dir = ensure_dir(outputs.get("dir") or workdir)
        out_bam = outputs.get("bam") or os.path.join(out_dir, f"{sample_id}.primary.bam")

        cmd = [
            str(params.get("java_bin", "java")),
            f"-XX:ParallelGCThreads={params.get('parallel_gc_threads', 14)}",
            f"-Xmx{params.get('xmx_gb', 16)}g",
            "-jar", str(params.get("jar", "/storage/apps/bin/picard.jar")),
            "MergeBamAlignment",
            "--UNMAPPED_BAM", unmapped,
            "--ALIGNED_BAM", aligned,
            "--REFERENCE_SEQUENCE", ref,
            "--OUTPUT", out_bam,
            "--CREATE_INDEX", str(params.get("create_index", True)).lower(),
            "--MAX_INSERTIONS_OR_DELETIONS", str(params.get("max_indel", -1)),
            "--CLIP_ADAPTERS", str(params.get("clip_adapters", False)).lower(),
            "--PRIMARY_ALIGNMENT_STRATEGY", params.get("primary_strategy", "MostDistant"),
        ]
        for a in params.get("retain_attributes", ["XS"]):
            cmd += ["--ATTRIBUTES_TO_RETAIN", a]
        for eo in params.get("expected_orientations", ["FR", "RF"]):
            cmd += ["--EXPECTED_ORIENTATIONS", eo]

        image = params.get("image")
        if image:
            cmd_line = singularity_exec_cmd(
                image=image,
                argv=cmd,
                binds=normalize_binds(params.get("binds")),
                singularity_bin=params.get("singularity_bin", "singularity"),
            )
        else:
            cmd_line = " ".join(map(shlex.quote, cmd))

        return [cmd_line]

    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs,
            outputs=self.outputs,
            params=p,
            threads=1,
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id") or None,
            ensure_output_dir_key="dir",
        )