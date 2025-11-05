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

@register_task("samtools.index")
class SamtoolsIndexTask(Task):
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
    TYPE = "samtools.index"

    INPUTS = {
        "bam": {"type": "path", "required": True, "desc": "Sorted BAM to index"},
    }
    OUTPUTS = {
        "dir":   {"type": "dir",  "required": False, "desc": "Output directory (default: same directory as BAM or workdir)"},
        "index": {"type": "path", "required": False, "desc": "Index path (.bai or .csi)"},
    }
    DEFAULTS: Dict[str, Any] = {
        "threads": 8,
        "index_type": "bai",         # bai or csi
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

        idx_type = str(params.get("index_type", "bai")).lower()
        if outputs.get("index"):
            index_path = outputs["index"]
        else:
            index_path = os.path.join(out_dir, os.path.basename(bam) + (".bai" if idx_type == "bai" else ".csi"))

        idx_argv = [samtools, "index", "-@", str(threads)]
        idx_argv += ["-b" if idx_type == "bai" else "-c", bam, index_path]

        guard = f'if [ ! -s {shlex.quote(index_path)} ] || [ {shlex.quote(index_path)} -ot {shlex.quote(bam)} ]; then '
        def wrap(argv: List[str]) -> str:
            if image:
                wrapped = singularity_exec_cmd(image=image, argv=argv, binds=binds, singularity_bin=singularity_bin)
                return " ".join(map(shlex.quote, wrapped))
            return " ".join(map(shlex.quote, argv))

        return [f"{guard}{wrap(idx_argv)}; fi"]

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