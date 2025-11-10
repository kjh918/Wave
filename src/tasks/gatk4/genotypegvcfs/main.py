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

@register_task("gatk4.genotypegvcfs")
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

    TYPE = "gatk4.genotypegvcfs"

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
        "all_sites": False
    }
    

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id=None) -> List[Sequence[str] | str]:
        gvcf = inputs["gvcf"]
        ref = inputs["reference"]
        out_dir = ensure_dir(workdir)
        out_vcf = outputs.get("vcf") or os.path.join(out_dir, f"{sample_id}.vcf.gz")

        total_cmd_list = []

        if params['region'] == 'chrom':
            
            chrom_list = [f'chr{i}' for i in range(1,23)] + ['chrX','chrY','chrM']
            
            for chrom in chrom_list:
                _chrom_input_gvcf = gvcf.replace('{chrom}',chrom)
                _chrom_output_gvcf = out_vcf.replace('{chrom}',chrom)
                cmd = [
                    str(params.get("gatk_bin", "gatk")),
                    f"-XX:ParallelGCThreads={params.get('parallel_gc_threads', 4)}",
                    f"-Xmx{params.get('xmx_gb', 16)}g",
                    "GenotypeGVCFs",
                    f'-R {ref} ',
                    f'-L {chrom} ',
                    f'-V {_chrom_input_gvcf} ',
                    f'-O {_chrom_output_gvcf} ',
                ]

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
                total_cmd_list.append(" ".join(cmd_line))

            return ["\n".join(total_cmd_list)]
        else:
            cmd = [
                    str(params.get("gatk_bin", "gatk")),
                    f"-XX:ParallelGCThreads={params.get('parallel_gc_threads', 4)}",
                    f"-Xmx{params.get('xmx_gb', 16)}g",
                    "GenotypeGVCFs",
                    f'-R {ref} ',
                    f'-V {gvcf} ',
                    f'-O {vcf} ',
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