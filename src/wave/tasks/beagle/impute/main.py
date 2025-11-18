from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional
import os, shlex

from wave.core.task import Task
from wave.core.task_registry import register_task
from wave.utils.task_utils import (
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

@register_task("beagle.impute")
class BeagleImputeTask(Task):
    """
    Beagle 5.5 기반 임퓨테이션/위상화(필요 시) 태스크.

    INPUTS:
      input_vcf : 대상 샘플 VCF/BCF(.tbi/.csi) — 보통 low-pass genotypes (필수)
      ref_vcf   : 레퍼런스 패널 VCF/BCF(.tbi/.csi) (권장/보통 필수)
      genetic_map : Beagle map 파일(.map) (선택, 제공 시 정확도 향상)
      chrom     : "chr1" 같은 단일 염색체 이름 (선택)
      interval  : "chr1:1-5,000,000" 같은 구간 (선택; chrom 대신 interval 우선)
                  ※ interval이 있으면 그 구간만 처리

    OUTPUTS:
      vcf       : 결과 VCF.gz 경로 (미지정 시 {workdir}/{base}.beagle.vcf.gz)
      dir       : 출력 디렉토리(선택; vcf 미지정 시 기준)

    PARAMS:
      threads       : int (기본 8) → nthreads
      xmx_gb        : int (기본 16) → -Xmx
      java_bin      : "java" (기본)
      beagle_jar    : "/storage/apps/bin/beagle.5.5.jar" (기본)
      impute        : bool (기본 True)   # Beagle impute 파이프라인 on/off
      gp            : bool (기본 True)   # gp 출력
      ne            : int (선택)         # 유효집단크기
      err           : float (선택)       # 에러율
      seed          : int (선택)
      image         : Singularity 이미지(optional)
      binds         : list[str] | str | None
      singularity_bin : "singularity" (기본)

    출력은 `out=<prefix>` 규칙을 사용하므로, `outputs.vcf`가 지정되면 prefix를 자동 계산해 사용.
    """

    TYPE = "beagle.impute"

    # (문서/검증용) 스키마 간단 정의
    INPUTS: Dict[str, Any] = {
        "input_vcf":    {"type": "path",   "required": True,  "desc": "Target sample VCF/BCF"},
        "ref_vcf":      {"type": "path",   "required": False, "desc": "Reference panel VCF/BCF"},
        "genetic_map":  {"type": "path",   "required": False, "desc": "Beagle genetic map"},
        "chrom":        {"type": "string", "required": False, "desc": "Chromosome name"},
        "interval":     {"type": "string", "required": False, "desc": "Interval string (chr:start-end)"},
    }
    OUTPUTS: Dict[str, Any] = {
        "vcf": {"type": "path", "required": False, "desc": "Output VCF.gz"},
        "dir": {"type": "dir",  "required": False, "desc": "Base output directory"},
    }

    DEFAULTS: Dict[str, Any] = {
        "threads": 8,
        "xmx_gb": 16,
        "java_bin": "java",
        "beagle_jar": "/storage/apps/bin/beagle.5.5.jar",
        "impute": True,
        "gp": True,
        "ne": 100000,
        "err": None, #  If no err parameter is specified, the err parameter will be set equal 𝜃/(2(𝜃 + 𝐻)) where 𝐻 is the number of haplotypes and 𝜃 = 1/(0.5 + ln𝐻)
        "seed": 7890,
        "image": None,
        "binds": None,
        "singularity_bin": "singularity",
    }

    @staticmethod
    def _prefix_from_output_path(out_vcf: str) -> str:
        # out=<prefix> 규칙에 맞추기 위해 .vcf.gz / .vcf 제거
        if out_vcf.endswith(".vcf.gz"):
            return out_vcf[:-7]
        if out_vcf.endswith(".vcf"):
            return out_vcf[:-4]
        return out_vcf

    def _build_cmd(
        self, *,
        inputs: Dict[str, Any],
        outputs: Dict[str, Any],
        params: Dict[str, Any],
        threads: int,
        workdir: str,
        sample_id: Optional[str] = None,
    ) -> List[Sequence[str] | str]:

        gt = inputs.get("vcf")
        ref = inputs.get("ref")                # 없는 경우도 허용(단독 refine 용도)
        gmap = inputs.get("map")           # 없는 경우도 허용(단독 refine 용도)
        impute_vcf = outputs.get("impute_vcf")

        # 파라미터 정리
        java_bin   = str(params.get("java_bin", "java"))
        beagle_jar = str(params.get("beagle_jar", "/storage/apps/bin/beagle.5.5.jar"))
        nthreads   = int(params.get("threads", 8))
        xmx_gb     = int(params.get("xmx_gb", 16))
        impute     = bool(params.get("impute", True))
        gp         = bool(params.get("gp", True))
        chrom       = bool(params.get("chrom", True))
        ne         = params.get("ne")
        err        = params.get("err")
        seed       = params.get("seed")

        # Beagle argv 조립

        total_cmd = []
        if chrom:
            chrom_list = [f'chr{i}' for i in range(1,23)] + ['chrX','chrY'] # ,'chrM']
            for chrom in chrom_list:
                
                gt_chrom = gt.replace('{chrom}',chrom)
                ref_chrom = ref.replace('{chrom}',chrom)          # 없는 경우도 허용(단독 refine 용도)
                gmap_chrom = gmap.replace('{chrom}',chrom)

                impute_vcf_chrom = impute_vcf.replace('{chrom}',chrom)

                argv: List[str] = [
                    java_bin, f"-Xmx{xmx_gb}g",
                    "-jar", beagle_jar,
                    f"gt={gt_chrom}",
                    f"out={impute_vcf_chrom}",
                    f"nthreads={nthreads}",
                    f"impute={'true' if impute else 'false'}",
                    f"gp={'true' if gp else 'false'}",
                ]
                if ref:   argv += [f"ref={ref_chrom}"]
                if gmap:  argv += [f"map={gmap_chrom}"]
                if ne is not None:   argv += [f"ne={ne}"]
                if err is not None:  argv += [f"err={err}"]
                if seed is not None: argv += [f"seed={seed}"]
                cmd = " ".join(map(safe_quote, argv))
                total_cmd.append(cmd)

        return ['\n'.join(total_cmd)]


    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"dir": str(self.workdir)},
            params=p,
            threads=int(self.threads or p["threads"]),
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id"),
            ensure_output_dir_key="dir",
        )