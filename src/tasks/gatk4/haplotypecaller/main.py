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

@register_task("gatk4.haplotypecaller")
class GatkHaplotypeCallerTask(Task):
    """
    Gatk HaplotypeCallerTask Task

    Example:
        gatk -XX:ParallelGCThreads=14 -Xmx16384m -jar HaplotypeCaller \
            -R $ReferenceGenome \
            -I ${BamDir}/${SeqID}.bwa.mem.sam \
            -L ${region} \
            -O $df  \
            -ERC GVCF
            -stand-call-conf 30
            -plodiy 2
    """

    TYPE = "gatk4.haplotypecaller"

    INPUTS = {
        "bam": {"type": "path", "required": True, "desc": "BAM"},
        "reference": {"type": "path", "required": True, "desc": "REFERENCE FASTA"},
        "known_snp": {"type": "path", "required": True, "desc": "KNOWN VCF"},
    }
    OUTPUTS = {
        "gvcf": {"type": "path", "required": False, "desc": "gvcf"},
    }

    DEFAULTS = {
        "gatk_bin": "gatk",
        "xmx_gb": 16,
        "parallel_gc_threads": 4,
        "gender": 'UNKNOWN',
        "tmp_dir": 'tmp',
        "singularity_bin": "singularity",
    }

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id=None) -> List[Sequence[str] | str]:
        bam = inputs["bam"]
        ref = inputs["reference"]
        gvcf = outputs["aligned_bam"]

        out_dir = ensure_dir(workdir)
        out_bam = outputs.get("gvcf") or os.path.join(out_dir, f"{sample_id}.gvcf.gz")

        cmd = [
            str(params.get("gatk_bin", "gatk")),
            f"-XX:ParallelGCThreads={params.get('parallel_gc_threads', 4)}",
            f"-Xmx{params.get('xmx_gb', 16)}g",
            "HaplotypeCaller ",
            f'-R {ref} -I {bam} '
            f'-stand-call-conf 30 --dbsnp {dbsnp_vcf} '
            f'-O {out_gvcf} -ERC GVCF'
        ]

        if params['region'].startswith('ch'):
            if params['region'] == 'chrX':
                if params['gender'].upper() == 'MALE':
                    cmd += ['-ploidy 1', '-L chrX']
            elif params['region'] in ['chrY','chrM']:
                cmd += ['-ploidy 1', f'-L {params['region']}']
            else:
                cmd += [f'-L {params['region']}']
        
        
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