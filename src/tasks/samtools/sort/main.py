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

@register_task("samtools.sort")
class SamtoolsSortTask(Task):
    """
    Sort BAM (no indexing here).
    Inputs:
      bam_in: path to input BAM
    Outputs:
      dir  : output directory (optional; defaults to workdir)
      bam  : sorted BAM path (optional; auto-generated if omitted)
    Params:
      threads, mem_per_thread, tmp_dir, samtools_bin, image, binds, singularity_bin
    """
    TYPE = "samtools.sort"

    INPUTS = {
        "bam_in": {"type": "path", "required": True, "desc": "Input BAM to sort"},
    }
    OUTPUTS = {
        "dir": {"type": "dir", "required": False, "desc": "Output directory (default: workdir)"},
        "bam": {"type": "path", "required": False, "desc": "Sorted BAM path"},
    }
    DEFAULTS: Dict[str, Any] = {
        "threads": 8,
        "mem_per_thread": "1G",
        "tmp_dir": None,                # samtools -T (prefix)
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

        bam_in = inputs["bam_in"]
        out_dir = ensure_dir(outputs.get("dir") or workdir)

        if outputs.get("bam"):
            bam_sorted = outputs["bam"]
        else:
            base = sample_id or os.path.splitext(os.path.basename(bam_in))[0]
            bam_sorted = os.path.join(out_dir, f"{base}.sorted.bam")

        sort_argv = [
            samtools, "sort",
            "-@", str(threads),
            "-m", str(params.get("mem_per_thread", "1G")),
            "-o", bam_sorted,
        ]
        if params.get("tmp_dir"):
            sort_argv += ["-T", str(params["tmp_dir"])]
        sort_argv += [bam_in]

        guard = f'if [ ! -s {shlex.quote(bam_sorted)} ] || [ {shlex.quote(bam_sorted)} -ot {shlex.quote(bam_in)} ]; then '
        def wrap(argv: List[str]) -> str:
            if image:
                wrapped = singularity_exec_cmd(image=image, argv=argv, binds=binds, singularity_bin=singularity_bin)
                return " ".join(map(shlex.quote, wrapped))
            return " ".join(map(shlex.quote, argv))

        return [f"{guard}{wrap(sort_argv)}; fi"]

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