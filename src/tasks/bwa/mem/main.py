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

@register_task("bwa.mem")
class BwaMemTask(Task):
    """
    BWA MEM alignment task (Singularity supported)

    Runs:
      singularity exec -B /storage,/data /storage/images/bwa-0.7.17.sif \
          bwa mem -M -t {threads} -Y -L 50,50 \
          -R "@RG\\tID:{read_group_id}\\tPL:{PL}\\tLB:{LB}\\tSM:{SM}\\tCN:{CN}" {reference_genome} \
          {TrimFastqDir}/{SeqID}.trimmed_read1.fastq.gz {TrimFastqDir}/{SeqID}.trimmed_read2.fastq.gz \
          > {BamDir}/{SeqID}.bwa.mem.sam
    """

    TYPE = "bwa.mem"

    INPUTS = {
        "read1": {"type": "path", "required": True, "desc": "Trimmed read1 FASTQ"},
        "read2": {"type": "path", "required": True, "desc": "Trimmed read2 FASTQ"},
    }
    OUTPUTS = {
        "sam": {"type": "path", "required": False, "desc": "Output SAM file"},
    }

    DEFAULTS: Dict[str, Any] = {
        "threads": 8,
        "image": "/storage/images/bwa-0.7.17.sif",
        "binds": ["/storage", "/data"],
        "singularity_bin": "singularity",
        "bwa_bin": "bwa",
        "reference_genome": None,
        "read_group_id": None,
        "platform": "ILLUMINA",
        "library_name": None,
        "sample_id": None,
        "center": None,
    }

    def _build_cmd(
            self, *, inputs: Dict[str, Any], outputs: Dict[str, Any], params: Dict[str, Any],
            threads: int, workdir: str, sample_id: str | None = None
        ) -> List[Sequence[str] | str]:
        read1 = inputs.get("read1")
        read2 = inputs.get("read2")
        sam_out = outputs.get("sam") or os.path.join(out_dir, f"{sample_id}.bwa.mem.sam")

        image = params.get("image")
        binds = normalize_binds(params.get("binds"))
        singularity_bin = params.get("singularity_bin", "singularity")
        bwa_bin = params.get("bwa_bin", "bwa")
        index = params.get("index")
        if not index:
            raise ValueError("BWA index path (--PARAMS.index) is required")

        rg_line = (
            f"@RG\\tID:{params.get('read_group_id', sample_id)}"
            f"\\tPL:{params.get('platform', 'ILLUMINA')}"
            f"\\tLB:{params.get('library_name', sample_id)}"
            f"\\tSM:{params.get('sample_id', sample_id)}"
            f"\\tCN:{params.get('center', 'CENTER')}"
        )

        core = [
            bwa_bin, "mem", "-M",
            "-t", str(threads),
            "-Y", "-L", f'{params.get("read_length")},{params.get("read_length")}',
            "-R", rg_line,
            index, str(read1), str(read2)
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