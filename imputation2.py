
# src/tasks/glimpse2/sample/main.py
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

@register_task("glimpse2.sample")
class Glimpse2SampleTask(Task):
    """
    GLIMPSE2_sample: ligate 결과(phasor/posterior VCF)를 best-guess GT/DS VCF로 변환.

    INPUTS:
      ligated_vcf : GLIMPSE2_ligate 산출물(.vcf.gz/.bcf) (필수)

    OUTPUTS:
      vcf  : 최종 genotype/dosage VCF.gz (미지정 시 {workdir}/{base}.imputed.vcf.gz)
      dir  : 출력 디렉토리(선택)

    PARAMS:
      threads        : int (기본 4)  # GLIMPSE2_sample이 내부 병렬을 쓰면 전달, 아니면 무시
      format         : "vcf"|"bcf" (기본 "vcf")
      write_dosage   : bool (기본 True)  # DS 필드 쓰기(지원되는 빌드 기준)
      sample_bin     : 실행파일명 (기본 "GLIMPSE2_sample")
      image          : Singularity 이미지(optional)
      binds          : list[str] | str | None
      singularity_bin : "singularity" (기본)
    """

    TYPE = "glimpse2.sample"

    INPUTS: Dict[str, Any] = {
        "ligated_vcf": {"type": "path", "required": True, "desc": "GLIMPSE2_ligate VCF/BCF"},
    }
    OUTPUTS: Dict[str, Any] = {
        "vcf": {"type": "path", "required": False, "desc": "Output imputed VCF.gz"},
        "dir": {"type": "dir",  "required": False, "desc": "Base output directory"},
    }
    DEFAULTS: Dict[str, Any] = {
        "threads": 4,
        "format": "vcf",
        "write_dosage": True,
        "sample_bin": "GLIMPSE2_sample",
        "image": None,
        "binds": None,
        "singularity_bin": "singularity",
    }

    @staticmethod
    def _default_out_path(workdir: str, ligated_vcf: str) -> str:
        base = os.path.basename(str(ligated_vcf)).replace(".vcf.gz","").replace(".bcf","").replace(".vcf","")
        return os.path.join(workdir, f"{base}.imputed.vcf.gz")

    def _build_cmd(
        self, *, inputs, outputs, params, threads, workdir, sample_id: Optional[str] = None
    ) -> List[Sequence[str] | str]:
        lig = inputs.get("ligated_vcf")
        if not lig:
            raise ValueError("[glimpse2.sample] inputs.ligated_vcf is required")

        out_dir = ensure_dir(outputs.get("dir") or workdir)
        out_vcf = outputs.get("vcf") or self._default_out_path(out_dir, lig)

        fmt = str(params.get("format", "vcf")).lower()  # "vcf"|"bcf"
        write_dose = bool(params.get("write_dosage", True))
        sample_bin = str(params.get("sample_bin", "GLIMPSE2_sample"))

        argv = [sample_bin, "--input", lig, "--output", out_vcf]
        if fmt in ("vcf", "bcf"):
            argv += ["--format", fmt]
        if write_dose:
            argv += ["--dosage"]

        image = params.get("image")
        if image:
            cmd = singularity_exec_cmd(
                image=str(image),
                argv=argv,
                binds=normalize_binds(params.get("binds")),
                singularity_bin=str(params.get("singularity_bin", "singularity")),
            )
            # index
            idx = singularity_exec_cmd(
                image=str(image),
                argv=["bcftools", "index", "-f", out_vcf],
                binds=normalize_binds(params.get("binds")),
                singularity_bin=str(params.get("singularity_bin", "singularity")),
            )
            return [" ".join(map(shlex.quote, cmd)), " ".join(map(shlex.quote, idx))]
        else:
            return [" ".join(map(shlex.quote, argv)),
                    " ".join(["bcftools", "index", "-f", shlex.quote(out_vcf)])]

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


# src/tasks/beagle/impute/main.py
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
        "ne": None,
        "err": None,
        "seed": None,
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

        gt = inputs.get("input_vcf")
        if not gt:
            raise ValueError("[beagle.impute] inputs.input_vcf is required")

        ref = inputs.get("ref_vcf")                # 없는 경우도 허용(단독 refine 용도)
        gmap = inputs.get("genetic_map")

        # 영역 제한
        interval = inputs.get("interval")
        chrom    = inputs.get("chrom")

        # 출력 경로
        out_dir = ensure_dir(outputs.get("dir") or workdir)
        out_vcf = outputs.get("vcf")
        if not out_vcf:
            base = sample_id or os.path.basename(str(gt)).replace(".vcf.gz","").replace(".bcf","").replace(".vcf","")
            out_vcf = os.path.join(out_dir, f"{base}.beagle.vcf.gz")
        prefix = self._prefix_from_output_path(out_vcf)

        # 파라미터 정리
        java_bin   = str(params.get("java_bin", "java"))
        beagle_jar = str(params.get("beagle_jar", "/storage/apps/bin/beagle.5.5.jar"))
        nthreads   = int(params.get("threads", 8))
        xmx_gb     = int(params.get("xmx_gb", 16))
        impute     = bool(params.get("impute", True))
        gp         = bool(params.get("gp", True))
        ne         = params.get("ne")
        err        = params.get("err")
        seed       = params.get("seed")

        # Beagle argv 조립
        argv: List[str] = [
            java_bin, f"-Xmx{xmx_gb}g",
            "-jar", beagle_jar,
            f"gt={gt}",
            f"out={prefix}",
            f"nthreads={nthreads}",
            f"impute={'true' if impute else 'false'}",
            f"gp={'true' if gp else 'false'}",
        ]
        if ref:   argv += [f"ref={ref}"]
        if gmap:  argv += [f"map={gmap}"]
        if interval: argv += [f"interval={interval}"]
        elif chrom:  argv += [f"chrom={chrom}"]
        if ne is not None:   argv += [f"ne={ne}"]
        if err is not None:  argv += [f"err={err}"]
        if seed is not None: argv += [f"seed={seed}"]

        # 컨테이너 래핑
        image = params.get("image")
        if image:
            cmd = singularity_exec_cmd(
                image=str(image),
                argv=argv,
                binds=normalize_binds(params.get("binds")),
                singularity_bin=str(params.get("singularity_bin", "singularity")),
            )
            run_line = shlex.join(cmd)
            # Beagle는 out=<prefix>로 파일 생성 → prefix.vcf.gz 가 최종 산출물
            idx_line = shlex.join(singularity_exec_cmd(
                image=str(image),
                argv=["bcftools", "index", "-f", f"{prefix}.vcf.gz"],
                binds=normalize_binds(params.get("binds")),
                singularity_bin=str(params.get("singularity_bin", "singularity")),
            ))
        else:
            run_line = shlex.join(argv)
            idx_line = " ".join(["bcftools", "index", "-f", shlex.quote(f"{prefix}.vcf.gz")])

        return [run_line, idx_line]

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