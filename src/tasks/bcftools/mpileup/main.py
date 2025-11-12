from __future__ import annotations
from typing import Dict, Any, List, Sequence
from pathlib import Path
import os, shlex
from glob import glob

from src.tasks.task import Task
from src.tasks.task_registry import register_task
from src.tasks.utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    to_sh_from_builder,
)


def safe_quote(arg):
    # 리다이렉션 문자나 파이프(|) 등은 그대로 출력
    if arg in {">", ">>", "<", "|", "2>", "&>", "&&", "||"}:
        return arg
    return shlex.quote(arg)

@register_task("bcftools.mpileup")
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

        bam_path = inputs.get("bam_path")
        reference = inputs.get("reference")
        ploidy_file = inputs.get("ploidy_file")

        sample_id = str(params.get("sample_id", "bcftools"))
        out_vcf = outputs.get("vcf") 
        
        bcftools = str(params.get("bcftools", "bcftools"))
        min_map_qual = str(params.get("min_map_qual", 20))
        min_base_qual = str(params.get("min_base_qual", 20))
        threads = str(params.get("min_base_qual", 4))
        singularity_bin = str(params.get("singularity_bin", 'singularity'))
        total_cmd = []

        for bam in glob(bam_path):
            image = params.get("image")
            binds = normalize_binds(params.get("reference"))

            sample_id = os.path.basename(bam).replace('.bam','.vcf.gz')
            core = [
                bcftools, 'mpileup',
                '-f', reference,
                '-a', 'AD,DP,SP',
                '-q', min_map_qual,
                '-Q', min_base_qual,
                '--threads',threads,
                '-Ou', f'{bam}','|',
                bcftools, 'call',
                '--ploidy-file', ploidy_file,
                '-mv','-Oz','-o', f'{out_vcf}/{sample_id}'
            ]

            # if image:
            #     cmd = singularity_exec_cmd(
            #         image=image,
            #         argv=core,
            #         binds=binds,
            #         singularity_bin=singularity_bin
            #     ) 
            #     # print(cmd)
            #     cmd = " ".join(map(safe_quote, cmd))
            #     total_cmd.append(cmd)
            # else:
            cmd = " ".join(map(safe_quote, core))#  + f" {redirect}"
            total_cmd.append(cmd)

        return ['\n'.join(total_cmd)]

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