from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional
import os, shlex

from src.wave.core.task import Task
from src.wave.core.task_registry import register_task
from src.wave.utils.task_utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    to_sh_from_builder,
)

@register_task("samtools.depth")
class SamtoolsDepthTask(Task):
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
    TYPE = "samtools.depth"

    INPUTS = {
        "bam": {"type": "path", "required": True, "desc": "Sorted BAM to index"},
    }
    OUTPUTS = {
        "bed": {"type": "path", "required": False, "desc": "Index path (.bai or .csi)"},
    }
    DEFAULTS: Dict[str, Any] = {
        "threads": 8,
        "all_position": "true",
        "bed_format": "true",     
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
        out_depth = outputs.get("depth") or os.path.join(out_dir, f"{base}.depth.txt")

        all_position = bool(params.get("all_position", True))
        bed_format = bool(params.get("bed_format", True))

        argv = [
            samtools,
            "depth",
            f"{'-a' if all_position else ''}",
            f"-f", f"{bam}",
            f">",f"{out_depth}"
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