# src/tasks/fastp/fastp.py
from __future__ import annotations
from typing import Dict, Any, List, Optional, Sequence
from pathlib import Path
import os

from src.tasks.task import Task
from src.tasks.task_registry import register_task
from src.tasks.utils import (   # ← utils 로 고정
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    to_sh_from_builder,
)

def _prefix_from_read1(r1: str) -> str:
    name = os.path.basename(r1)
    for tag in ["_R1.fastq.gz", "_R1.fq.gz", ".R1.fastq.gz", ".R1.fq.gz",
                "_R1.fastq", "_R1.fq", ".R1.fastq", ".R1.fq"]:
        if name.endswith(tag):
            return name[:-len(tag)]
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"]:
        if name.endswith(ext):
            return name[:-len(ext)]
    return os.path.splitext(name)[0]

# @TaskRegistry.register('fastp')
@register_task('fastp')
class FastpTask(Task):
    TYPE = "fastp"

    INPUTS: Dict[str, Any] = {
        "read1": {"type": "path", "required": True,  "desc": "R1 FASTQ(.gz) absolute path"},
        "read2": {"type": "path", "required": False, "desc": "R2 FASTQ(.gz) absolute path (optional)"},
    }
    OUTPUTS: Dict[str, Any] = {
        "dir":   {"type": "dir",  "required": False, "desc": "Result directory for trimmed FASTQs and reports"},
        "read1": {"type": "path", "required": False, "desc": "Trimmed R1 FASTQ path"},
        "read2": {"type": "path", "required": False, "desc": "Trimmed R2 FASTQ path (if paired)"},
        "json":  {"type": "path", "required": False, "desc": "fastp JSON report"},
        "html":  {"type": "path", "required": False, "desc": "fastp HTML report"},
    }

    DEFAULTS: Dict[str, Any] = {
        "threads": 4,
        "PicoplexGold": False,
        "trim_front1": 0,
        "trim_front2": 0,
        "length_required": 100,
        "average_qual": 10,
        "qualified_quality_phred": 15,
        "adapter_sequence": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        "adapter_sequence_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "image": None,
        "binds": None,
        "fastp_bin": "fastp",
        "singularity_bin": "singularity",
    }

    def _build_fastp_cmd(
        self,
        *,
        inputs: Dict[str, Any],
        outputs: Dict[str, Any],
        params: Dict[str, Any],
        threads: int,
        workdir: str,
        sample_id: Optional[str] = None,
    ) -> List[Sequence[str] | str]:
        r1 = inputs.get("read1")
        if not r1:
            raise ValueError("[fastp] INPUT.read1 is required")
        r2 = inputs.get("read2")

        base_dir = ensure_dir(outputs.get("dir") or workdir)
        prefix = _prefix_from_read1(str(r1))

        out_r1 = outputs.get("read1") or os.path.join(base_dir, f"{prefix}.trimmed_R1.fastq.gz")
        out_r2 = outputs.get("read2") or (os.path.join(base_dir, f"{prefix}.trimmed_R2.fastq.gz") if r2 else None)
        out_json = outputs.get("json") or os.path.join(base_dir, f"{prefix}.fastp.json")
        out_html = outputs.get("html") or os.path.join(base_dir, f"{prefix}.fastp.html")

        picoplex_gold = bool(params.get("PicoplexGold", self.DEFAULTS["PicoplexGold"]))
        trim_front1 = int(params.get("trim_front1", 14 if picoplex_gold else 0))
        trim_front2 = int(params.get("trim_front2", 14 if picoplex_gold else 0))

        length_required = int(params.get("length_required", self.DEFAULTS["length_required"]))
        average_qual = int(params.get("average_qual", self.DEFAULTS["average_qual"]))
        qualified_quality_phred = int(params.get("qualified_quality_phred", self.DEFAULTS["qualified_quality_phred"]))

        adapter_sequence = str(params.get("adapter_sequence", self.DEFAULTS["adapter_sequence"]))
        adapter_sequence_r2 = str(params.get("adapter_sequence_r2", self.DEFAULTS["adapter_sequence_r2"]))

        fastp_bin = str(params.get("fastp_bin", "fastp"))
        image = params.get("image")
        binds = normalize_binds(params.get("binds"))
        singularity_bin = str(params.get("singularity_bin", "singularity"))

        argv: List[str] = [fastp_bin, "--thread", str(int(threads))]
        argv += ["--in1", str(r1), "--out1", out_r1]
        if r2:
            argv += ["--in2", str(r2)]
            if out_r2:
                argv += ["--out2", out_r2]
            argv += ["--detect_adapter_for_pe"]
            if adapter_sequence_r2:
                argv += ["--adapter_sequence_r2", adapter_sequence_r2]

        argv += [
            "--trim_poly_g",
            "--json", out_json,
            "--html", out_html,
            "--length_required", str(length_required),
            "--average_qual", str(average_qual),
            "--qualified_quality_phred", str(qualified_quality_phred),
        ]
        if adapter_sequence:
            argv += ["--adapter_sequence", adapter_sequence]
        if trim_front1 > 0:
            argv += ["--trim_front1", str(trim_front1)]
        if r2 and trim_front2 > 0:
            argv += ["--trim_front2", str(trim_front2)]

        if image:
            return [singularity_exec_cmd(
                image=str(image),
                argv=argv,
                binds=binds,
                singularity_bin=singularity_bin,
            )]
        else:
            return [argv]

    def to_sh(self) -> List[str]:
        p: Dict[str, Any] = {**self.DEFAULTS, **(self.params or {})}
        th = int(self.threads or p.get("threads", 4))
        return to_sh_from_builder(
            builder=self._build_fastp_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"dir": str(self.workdir)},
            params=p,
            threads=th,
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id") or None,
            ensure_output_dir_key="dir",
        )