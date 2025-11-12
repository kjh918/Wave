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

@register_task("samtools.flagstats")
class SamtoolsFlagstatsTask(Task):
    """
    Index BAM only (expects a sorted BAM).
    Inputs:
      bam: sorted BAM path
    Outputs:
      dir  : base directory (optional; if output index path omitted, use dir/basename.bam.{bai|csi})
      index: index path (optional)
    Params:
      index_type (bai|csi), threads, samtools_bin, image, binds, singularity_bin
    """
    TYPE = "samtools.flagstats"

    INPUTS = {
        "bam": {"type": "path", "required": True, "desc": "Sorted BAM to index"},
    }
    OUTPUTS = {
        "txt":   {"type": "path",  "required": False, "desc": "Output directory (default: same directory as BAM or workdir)"}
    }
    DEFAULTS: Dict[str, Any] = {
        "threads": 8,
        "samtools_bin": "samtools",
        "image": None,
        "binds": None,
        "singularity_bin": "singularity",
    }

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id: Optional[str] = None) -> List[Sequence[str] | str]:
        samtools = str(params.get("samtools_bin", "samtools"))
        image = params.get("image")
        binds = normalize_binds(params.get("binds"))
        singularity_bin = str(params.get("singularity_bin", "singularity"))

        bam = inputs["bam"]
        bam_dir = os.path.dirname(bam) or workdir
        out_dir = ensure_dir(outputs.get("dir") or bam_dir)
        argv = [
            java_bin,
            f"-XX:ParallelGCThreads={pgc}",
            f"-Xmx{xmx_gb}g",
            "-jar", jar_path,
            "MergeBamAlignment",
            "--UNMAPPED_BAM", unmapped_bam,
            "--ALIGNED_BAM", aligned_bam,
            "--REFERENCE_SEQUENCE", reference,
            "--OUTPUT", out_bam,
            "--CREATE_INDEX", str(params.get("create_index", True)).lower(),
            "--MAX_INSERTIONS_OR_DELETIONS", str(params.get("max_indels", -1)),
            "--CLIP_ADAPTERS", str(params.get("clip_adapters", False)).lower(),
            "--PRIMARY_ALIGNMENT_STRATEGY", str(params.get("primary_strategy", "MostDistant")),
        ]
        
        # 컨테이너 래핑
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
        th = int(self.threads or p.get("threads", 8))
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"dir": str(self.workdir)},
            params=p,
            threads=th,
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id") or None,
            ensure_output_dir_key="dir",
        )