from __future__ import annotations
from typing import Dict, Any, List, Sequence
from pathlib import Path
import os, shlex

from src.tasks.task import Task, TaskRegistry
from src.tasks.util import ensure_dir, normalize_binds, singularity_exec_cmd, to_sh_from_builder


@TaskRegistry.register
class BwaMemTask(Task):
    """
    BWA MEM alignment task (Singularity supported)

    Runs:
      singularity exec -B /storage,/data /storage/images/bwa-0.7.17.sif \
          bwa mem -M -t {threads} -Y -L 50,50 \
          -R "@RG\\tID:{RG_ID}\\tPL:{PL}\\tLB:{LB}\\tSM:{SM}\\tCN:{CN}" {BwaIndex} \
          {TrimFastqDir}/{SeqID}.trimmed_R1.fastq.gz {TrimFastqDir}/{SeqID}.trimmed_R2.fastq.gz \
          > {BamDir}/{SeqID}.bwa.mem.sam
    """

    TYPE = "bwa.mem"

    INPUTS = {
        "r1": {"type": "path", "required": True, "desc": "Trimmed R1 FASTQ"},
        "r2": {"type": "path", "required": True, "desc": "Trimmed R2 FASTQ"},
    }
    OUTPUTS = {
        "sam": {"type": "path", "required": False, "desc": "Output SAM file"},
        "dir": {"type": "dir", "required": False, "desc": "Work/output directory"},
    }

    DEFAULTS: Dict[str, Any] = {
        "threads": 8,
        "image": "/storage/images/bwa-0.7.17.sif",
        "binds": ["/storage", "/data"],
        "singularity_bin": "singularity",
        "bwa_bin": "bwa",
        "index": None,
        "rg_id": None,
        "rg_pl": "ILLUMINA",
        "rg_lb": None,
        "rg_sm": None,
        "rg_cn": None,
    }

    def _build_cmd(
        self, *, inputs: Dict[str, Any], outputs: Dict[str, Any], params: Dict[str, Any],
        threads: int, workdir: str, sample_id: str | None = None
    ) -> List[Sequence[str] | str]:
        r1 = inputs.get("r1")
        r2 = inputs.get("r2")
        out_dir = ensure_dir(outputs.get("dir") or workdir)
        sam_out = outputs.get("sam") or os.path.join(out_dir, f"{sample_id}.bwa.mem.sam")

        image = params.get("image")
        binds = normalize_binds(params.get("binds"))
        singularity_bin = params.get("singularity_bin", "singularity")
        bwa_bin = params.get("bwa_bin", "bwa")
        index = params.get("index")
        if not index:
            raise ValueError("BWA index path (--PARAMS.index) is required")

        rg_line = (
            f"@RG\\tID:{params.get('rg_id', sample_id)}"
            f"\\tPL:{params.get('rg_pl', 'ILLUMINA')}"
            f"\\tLB:{params.get('rg_lb', sample_id)}"
            f"\\tSM:{params.get('rg_sm', sample_id)}"
            f"\\tCN:{params.get('rg_cn', 'CENTER')}"
        )

        core = [
            bwa_bin, "mem", "-M",
            "-t", str(threads),
            "-Y", "-L", "50,50",
            "-R", rg_line,
            index, str(r1), str(r2)
        ]

        redirect = f"> {sam_out}"
        if image:
            cmd = singularity_exec_cmd(
                image=image,
                argv=core,
                binds=binds,
                singularity_bin=singularity_bin
            ) + f" {redirect}"
        else:
            cmd = " ".join(map(shlex.quote, core)) + f" {redirect}"

        return [cmd]

    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        th = int(p.get("threads", 1))
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs,
            outputs=self.outputs,
            params=p,
            threads=th,
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id") or None,
            ensure_output_dir_key="dir",
        )