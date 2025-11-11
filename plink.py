# src/tasks/plink2/compare/main.py
from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional, Union
import os, shlex, glob

from src.tasks.task import Task
from src.tasks.task_registry import register_task
from src.tasks.utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    to_sh_from_builder,
)

@register_task("plink2.compare")
class Plink2CompareTask(Task):
    """
    여러 샘플 imputed VCF(.vcf.gz)들을 병합하고, PGEN 변환 후
    관련도/결측/PCA 리포트를 산출.

    INPUTS:
      vcf_list : list[str] | str(glob)  — 병합할 VCF들의 리스트 또는 글롭 패턴 (필수)

    OUTPUTS:
      merged_vcf   : 병합된 cohort VCF.gz (선택; 기본 {workdir}/cohort.merged.vcf.gz)
      pgen_prefix  : PGEN 세트 prefix (선택; 기본 {workdir}/cohort)
      out_dir      : 결과 디렉토리 (선택; 기본 workdir)

    PARAMS:
      threads      : int (기본 8)
      bcftools_bin : "bcftools"
      plink2_bin   : "plink2"
      maf          : float (선택; PGEN 변환 시 필터)
      geno         : float (선택; PGEN 변환 시 필터, sample-based는 --mind로 후처리 가능)
      memory_mb    : int (선택; plink2 --memory)
      image        : Singularity 이미지(optional)
      binds        : list[str] | str | None
      singularity_bin : "singularity"
    """

    TYPE = "plink2.compare"

    INPUTS: Dict[str, Any] = {
        "vcf_list": {"type": "any", "required": True, "desc": "list[str] or glob pattern"},
    }
    OUTPUTS: Dict[str, Any] = {
        "merged_vcf":  {"type": "path", "required": False},
        "pgen_prefix": {"type": "path", "required": False},
        "out_dir":     {"type": "dir",  "required": False},
    }
    DEFAULTS: Dict[str, Any] = {
        "threads": 8,
        "bcftools_bin": "bcftools",
        "plink2_bin": "plink2",
        "maf": None,
        "geno": None,
        "memory_mb": None,
        "image": None,
        "binds": None,
        "singularity_bin": "singularity",
    }

    @staticmethod
    def _expand_vcf_list(vl: Union[List[str], str]) -> List[str]:
        if isinstance(vl, list):
            return [str(x) for x in vl]
        # glob 패턴
        return sorted(map(str, glob.glob(vl)))

    def _build_cmd(
        self, *, inputs, outputs, params, threads, workdir, sample_id: Optional[str] = None
    ) -> List[Sequence[str] | str]:

        vcfs = self._expand_vcf_list(inputs.get("vcf_list"))
        if not vcfs:
            raise ValueError("[plink2.compare] inputs.vcf_list resolved to empty")

        out_dir = ensure_dir(outputs.get("out_dir") or workdir)
        merged_vcf = outputs.get("merged_vcf") or os.path.join(out_dir, "cohort.merged.vcf.gz")
        pfx = outputs.get("pgen_prefix") or os.path.join(out_dir, "cohort")

        bcftools_bin = str(params.get("bcftools_bin", "bcftools"))
        plink2_bin   = str(params.get("plink2_bin", "plink2"))
        maf = params.get("maf")
        geno = params.get("geno")
        mem = params.get("memory_mb")
        image = params.get("image")
        binds = normalize_binds(params.get("binds"))
        sing = str(params.get("singularity_bin", "singularity"))

        # 1) bcftools merge
        merge_argv = [bcftools_bin, "merge", "-m", "none", "-O", "z", "-o", merged_vcf, "--threads", str(int(threads))]
        merge_argv += vcfs
        # 2) index
        index_argv = [bcftools_bin, "index", "-f", merged_vcf]
        # 3) plink2 pgen 변환 (+선택필터)
        pgen_argv = [plink2_bin, "--vcf", merged_vcf, "--make-pgen", "--out", pfx, "--threads", str(int(threads))]
        if maf is not None:
            pgen_argv += ["--maf", str(maf)]
        if geno is not None:
            pgen_argv += ["--geno", str(geno)]
        if mem is not None:
            pgen_argv += ["--memory", str(int(mem))]

        # 4) 관련도/결측/PCA
        genome_argv = [plink2_bin, "--pfile", pfx, "--genome", "full", "--out", pfx + ".related", "--threads", str(int(threads))]
        miss_argv   = [plink2_bin, "--pfile", pfx, "--missing", "sample-only", "--out", pfx + ".missing", "--threads", str(int(threads))]
        pca_argv    = [plink2_bin, "--pfile", pfx, "--pca", "10", "--out", pfx + ".pca", "--threads", str(int(threads))]

        if image:
            def wrap(a): 
                return " ".join(map(shlex.quote, singularity_exec_cmd(image=str(image), argv=a, binds=binds, singularity_bin=sing)))
            return [
                wrap(merge_argv),
                wrap(index_argv),
                wrap(pgen_argv),
                wrap(genome_argv),
                wrap(miss_argv),
                wrap(pca_argv),
            ]
        else:
            def j(a): return " ".join(map(shlex.quote, a))
            return [j(merge_argv), j(index_argv), j(pgen_argv), j(genome_argv), j(miss_argv), j(pca_argv)]

    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"out_dir": str(self.workdir)},
            params=p,
            threads=int(self.threads or p["threads"]),
            workdir=str(self.workdir),
            sample_id=None,
            ensure_output_dir_key="out_dir",
        )