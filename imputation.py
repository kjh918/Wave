# src/tasks/glimpse2/chunk/main.py
from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional, Union
import os, shlex

from src.tasks.task import Task
from src.tasks.task_registry import register_task
from src.tasks.utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    to_sh_from_builder,
)

@register_task("glimpse2.chunk")
class Glimpse2ChunkTask(Task):
    """
    GLIMPSE2 공식 chunk 생성기(GLIMPSE2_chunk)로 chunks.tsv를 만든다.

    INPUTS:
      reference_fai : 참조 FASTA의 .fai (필수)
      (optional) reference_fasta : .dict 체크 등의 추가 로직에 활용 가능 (현재는 미사용)

    OUTPUTS:
      chunks_tsv : 출력 TSV 경로 (기본: {workdir}/chunks.tsv)
      dir        : 출력 디렉토리(옵션; chunks_tsv 미지정 시 기준)

    PARAMS:
      window_size : int (기본 2_000_000)
      step_size   : int (기본 1_800_000)
      buffer_size : int | None (있으면 --buffer-size 전달)
      region      : str | None (예: "chr20" 또는 "chr20:1-10,000,000")
      regions     : list[str] | None (여러 개 region을 각각 --region으로 추가)
      image       : str | None (Singularity 이미지; 없으면 로컬 바이너리 사용)
      binds       : list[str] | str | None (Singularity 바인딩)
      singularity_bin : str (기본 "singularity")
      chunk_bin   : str (기본 "GLIMPSE2_chunk" — 로컬 실행 시)
    """

    TYPE = "glimpse2.chunk"

    INPUTS: Dict[str, Any] = {
        "reference_fai": {"type": "path", "required": True, "desc": "reference FASTA .fai index"},
        "reference_fasta": {"type": "path", "required": False, "desc": "reference FASTA (optional)"},
    }

    OUTPUTS: Dict[str, Any] = {
        "chunks_tsv": {"type": "path", "required": False, "desc": "output chunks TSV"},
        "dir": {"type": "dir", "required": False, "desc": "base output directory"},
    }

    DEFAULTS: Dict[str, Any] = {
        "window_size": 2_000_000,
        "step_size":   1_800_000,
        "buffer_size": None,
        "region": None,          # 단일 region
        "regions": None,         # 여러 region
        "image": None,           # None이면 로컬 GLIMPSE2_chunk 실행
        "binds": None,           # list | str | None
        "singularity_bin": "singularity",
        "chunk_bin": "GLIMPSE2_chunk",
    }

    # --- 실제 커맨드 빌더 (stdout을 파일로 redirect 해야 하므로 문자열 라인으로 반환)
    def _build_cmd(
        self,
        *,
        inputs: Dict[str, Any],
        outputs: Dict[str, Any],
        params: Dict[str, Any],
        threads: int,           # 인터페이스 통일용 (사용 안 함)
        workdir: str,
        sample_id: Optional[str] = None,   # 인터페이스 통일용
    ) -> List[Sequence[str] | str]:
        fai = inputs.get("reference_fai")
        if not fai:
            raise ValueError("[glimpse2.chunk] INPUTS.reference_fai is required")

        # 출력 경로 확정
        base_dir = ensure_dir(outputs.get("dir") or workdir)
        out_tsv  = outputs.get("chunks_tsv") or os.path.join(base_dir, "chunks.tsv")

        # 옵션 파라미터
        window  = int(params.get("window_size", self.DEFAULTS["window_size"]))
        step    = int(params.get("step_size", self.DEFAULTS["step_size"]))
        buffer_ = params.get("buffer_size", self.DEFAULTS["buffer_size"])
        region  = params.get("region", None)
        regions = params.get("regions", None)

        # 바이너리/컨테이너
        image = params.get("image")
        binds = normalize_binds(params.get("binds"))
        singularity_bin = str(params.get("singularity_bin", "singularity"))
        chunk_bin = str(params.get("chunk_bin", "GLIMPSE2_chunk"))

        # argv 구성 (GLIMPSE2_chunk는 stdout으로 쓰므로, 나중에 '>' 리다이렉션)
        argv: List[str] = [chunk_bin, "--input", str(fai), "--window-size", str(window), "--step-size", str(step)]
        if buffer_ is not None:
            argv += ["--buffer-size", str(int(buffer_))]

        # region/regions 추가
        def _as_list(x: Union[str, List[str], None]) -> List[str]:
            if x is None:
                return []
            if isinstance(x, list):
                return [str(v) for v in x]
            return [str(x)]

        all_regions = _as_list(region) + _as_list(regions)
        for r in all_regions:
            if r:
                argv += ["--region", r]

        # 컨테이너 래핑
        if image:
            cmd_tokens = singularity_exec_cmd(
                image=str(image),
                argv=argv,
                binds=binds,
                singularity_bin=singularity_bin,
            )
        else:
            cmd_tokens = argv

        # stdout 리다이렉션을 포함한 단일 문자열 라인으로 반환
        line = f"{shlex.join(cmd_tokens)} > {shlex.quote(out_tsv)}"
        return [line]

# src/tasks/glimpse2/glchunk/main.py
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

@register_task("glimpse2.glchunk")
class Glimpse2GLChunkTask(Task):
    """
    한 chunk(region)에 대해 GL(VCF.gz)을 생성한다.

    INPUTS:
      bam            : BAM/CRAM (필수)
      reference_fasta: FASTA (필수, .fai 필요)
      chunks_tsv     : chunks.tsv 경로 (선택; region을 여기서 찾음)
      region         : "chr:start-end" (선택; chunks_tsv 대신 직접 지정 가능)
      chunk_id       : chunks.tsv의 4번째 컬럼 값 (선택; region 미지정 시 요구)

    OUTPUTS:
      vcf            : 출력 GL VCF.gz 경로 (미지정 시 {workdir}/{chunk_id}.vcf.gz)
      dir            : 출력 디렉토리 (미지정 시 workdir)

    PARAMS:
      bcftools_bin   : 기본 "bcftools"
      image          : Singularity 이미지 경로(optional)
      binds          : list[str] | str | None
      singularity_bin: 기본 "singularity"
      mpileup_opts   : list[str] 추가 옵션 예: ["-a","FORMAT/DP"]
      call_opts      : list[str] 추가 옵션 예: ["-Aim"]
    """

    TYPE = "glimpse2.glchunk"

    DEFAULTS: Dict[str, Any] = {
        "bcftools_bin": "bcftools",
        "image": None,
        "binds": None,
        "singularity_bin": "singularity",
        "mpileup_opts": [],
        "call_opts": ["-Aim"],
    }

    @staticmethod
    def _region_from_chunks(chunks_tsv: str, chunk_id: str) -> str:
        with open(chunks_tsv) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                chrom, s, e, cid = line.strip().split()[:4]
                if cid == chunk_id:
                    return f"{chrom}:{s}-{e}"
        raise ValueError(f"[glimpse2.glchunk] chunk_id not found in {chunks_tsv}: {chunk_id}")

    def _build_cmd(
        self, *, inputs: Dict[str, Any], outputs: Dict[str, Any], params: Dict[str, Any],
        threads: int, workdir: str, sample_id: Optional[str] = None
    ) -> List[Sequence[str] | str]:
        bam = inputs.get("bam")
        ref = inputs.get("reference_fasta")
        if not bam or not ref:
            raise ValueError("[glimpse2.glchunk] inputs.bam & inputs.reference_fasta are required")

        # region 결정
        region = inputs.get("region")
        if not region:
            chunks = inputs.get("chunks_tsv")
            cid = inputs.get("chunk_id")
            if not (chunks and cid):
                raise ValueError("[glimpse2.glchunk] require either inputs.region or (inputs.chunks_tsv + inputs.chunk_id)")
            region = self._region_from_chunks(str(chunks), str(cid))

        # 출력 경로
        out_dir = ensure_dir(outputs.get("dir") or workdir)
        out_vcf = outputs.get("vcf")
        if not out_vcf:
            cid = inputs.get("chunk_id") or region.replace(":", "_").replace("-", "_")
            out_vcf = os.path.join(out_dir, f"{cid}.vcf.gz")

        bcftools = str(params.get("bcftools_bin", "bcftools"))
        mp_opts  = list(map(str, params.get("mpileup_opts", [])))
        call_opts= list(map(str, params.get("call_opts", ["-Aim"])))

        # 파이프 구조: mpileup -> call
        mp = [bcftools, "mpileup", "-f", ref, "-r", region, "-Ou", bam, *mp_opts]
        cl = [bcftools, "call", *call_opts, "-Oz", "-o", out_vcf]

        image = params.get("image")
        if image:
            mp = singularity_exec_cmd(image=str(image), argv=mp, binds=normalize_binds(params.get("binds")),
                                      singularity_bin=str(params.get("singularity_bin","singularity")))
            cl = singularity_exec_cmd(image=str(image), argv=cl, binds=normalize_binds(params.get("binds")),
                                      singularity_bin=str(params.get("singularity_bin","singularity")))

        # 단일 라인으로 파이프 연결 + 인덱스
        line1 = f"{shlex.join(mp)} | {shlex.join(cl)}"
        line2 = f"{shlex.join((image and singularity_exec_cmd(image=str(image), argv=[bcftools,'index','-f',out_vcf], binds=normalize_binds(params.get('binds')), singularity_bin=str(params.get('singularity_bin','singularity'))) or [bcftools,'index','-f',out_vcf]))}"
        return [line1, line2]

    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"dir": str(self.workdir)},
            params=p,
            threads=int(self.threads or 1),
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id"),
            ensure_output_dir_key="dir",
        )

# src/tasks/glimpse2/phasechunk/main.py
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

@register_task("glimpse2.phasechunk")
class Glimpse2PhaseChunkTask(Task):
    """
    한 chunk(region)에 대해 GLIMPSE2_phase 실행.

    INPUTS:
      gl_vcf        : per-chunk GL VCF(.vcf.gz) (필수)
      chunks_tsv    : chunks.tsv (선택; region 추출용)
      region        : "chr:start-end" (선택; 직접 지정 가능)
      chunk_id      : chunks.tsv의 4번째 컬럼 값 (선택; region 미지정 시 요구)
      panel_vcf     : 참조 패널 VCF/BCF(.tbi/.csi) (필수)

    OUTPUTS:
      phased_vcf    : 출력 VCF.gz (미지정 시 {workdir}/{chunk_id}.phased.vcf.gz)
      dir           : 출력 디렉토리

    PARAMS:
      threads       : int (기본 8)
      image         : Singularity 이미지(optional)
      binds         : list | str | None
      singularity_bin: 기본 "singularity"
      phase_bin     : "GLIMPSE2_phase" (로컬 바이너리 이름)
    """

    TYPE = "glimpse2.phasechunk"

    DEFAULTS: Dict[str, Any] = {
        "threads": 8,
        "image": None,
        "binds": None,
        "singularity_bin": "singularity",
        "phase_bin": "GLIMPSE2_phase",
    }

    @staticmethod
    def _region_from_chunks(chunks_tsv: str, chunk_id: str) -> str:
        with open(chunks_tsv) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                chrom, s, e, cid = line.strip().split()[:4]
                if cid == chunk_id:
                    return f"{chrom}:{s}-{e}"
        raise ValueError(f"[glimpse2.phasechunk] chunk_id not found in {chunks_tsv}: {chunk_id}")

    def _build_cmd(
        self, *, inputs: Dict[str, Any], outputs: Dict[str, Any], params: Dict[str, Any],
        threads: int, workdir: str, sample_id: Optional[str] = None
    ) -> List[Sequence[str] | str]:
        gl_vcf = inputs.get("gl_vcf")
        panel  = inputs.get("panel_vcf")
        if not gl_vcf or not panel:
            raise ValueError("[glimpse2.phasechunk] inputs.gl_vcf & inputs.panel_vcf are required")

        # region 결정
        region = inputs.get("region")
        if not region:
            chunks = inputs.get("chunks_tsv")
            cid = inputs.get("chunk_id")
            if not (chunks and cid):
                raise ValueError("[glimpse2.phasechunk] require either inputs.region or (inputs.chunks_tsv + inputs.chunk_id)")
            region = self._region_from_chunks(str(chunks), str(cid))

        out_dir = ensure_dir(outputs.get("dir") or workdir)
        out_vcf = outputs.get("phased_vcf")
        if not out_vcf:
            cid = inputs.get("chunk_id") or os.path.basename(gl_vcf).replace(".vcf.gz","")
            out_vcf = os.path.join(out_dir, f"{cid}.phased.vcf.gz")

        phase_bin = str(params.get("phase_bin", "GLIMPSE2_phase"))
        argv = [phase_bin, "--input", str(gl_vcf), "--reference", str(panel),
                "--region", str(region), "--output", out_vcf, "--threads", str(int(params.get("threads", 8)))]

        image = params.get("image")
        if image:
            cmd = singularity_exec_cmd(
                image=str(image), argv=argv,
                binds=normalize_binds(params.get("binds")),
                singularity_bin=str(params.get("singularity_bin","singularity")),
            )
            line1 = shlex.join(cmd)
            idx   = shlex.join(singularity_exec_cmd(
                        image=str(image),
                        argv=["bcftools","index","-f",out_vcf],
                        binds=normalize_binds(params.get("binds")),
                        singularity_bin=str(params.get("singularity_bin","singularity")),
                    ))
        else:
            line1 = shlex.join(argv)
            idx   = " ".join(["bcftools","index","-f",shlex.quote(out_vcf)])

        return [line1, idx]

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

# src/tasks/glimpse2/ligate/main.py
from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional
import os, glob, shlex

from src.tasks.task import Task
from src.tasks.task_registry import register_task
from src.tasks.utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    to_sh_from_builder,
)

@register_task("glimpse2.ligate")
class Glimpse2LigateTask(Task):
    """
    per-chunk phased VCF들을 ligate하여 하나의 VCF로 합친다.

    INPUTS:
      phased_list   : 파일 리스트(.list) 경로 (선택)
      phased_glob   : 글롭 패턴 (선택; 예: "/work/S1/phased/*.vcf.gz")
      # 둘 중 하나는 필요. 둘 다 있으면 phased_list 우선.

    OUTPUTS:
      output_vcf    : 결과 VCF.gz (필수 또는 dir 제공 시 기본 파일명 생성)
      dir           : 출력 디렉토리(선택; output_vcf 미지정 시 기준)

    PARAMS:
      image         : Singularity 이미지(optional)
      binds         : list | str | None
      singularity_bin: 기본 "singularity"
      ligate_bin    : "GLIMPSE2_ligate" (로컬 바이너리)
    """

    TYPE = "glimpse2.ligate"

    DEFAULTS: Dict[str, Any] = {
        "image": None,
        "binds": None,
        "singularity_bin": "singularity",
        "ligate_bin": "GLIMPSE2_ligate",
    }

    def _build_cmd(
        self, *, inputs: Dict[str, Any], outputs: Dict[str, Any], params: Dict[str, Any],
        threads: int, workdir: str, sample_id: Optional[str] = None
    ) -> List[Sequence[str] | str]:
        # 입력 목록 준비
        list_path = inputs.get("phased_list")
        if not list_path:
            gpat = inputs.get("phased_glob")
            if not gpat:
                raise ValueError("[glimpse2.ligate] require inputs.phased_list or inputs.phased_glob")
            files = sorted(glob.glob(gpat))
            if not files:
                raise ValueError(f"[glimpse2.ligate] no files matched: {gpat}")
            list_path = os.path.join(workdir, "phased.list")
            with open(list_path, "w") as fh:
                for f in files: fh.write(f"{f}\n")

        out_dir = ensure_dir(outputs.get("dir") or workdir)
        out_vcf = outputs.get("output_vcf") or os.path.join(out_dir, "phased.ligated.vcf.gz")

        lig_bin = str(params.get("ligate_bin", "GLIMPSE2_ligate"))
        argv = [lig_bin, "--input", list_path, "--output", out_vcf]

        image = params.get("image")
        if image:
            cmd = singularity_exec_cmd(
                image=str(image), argv=argv,
                binds=normalize_binds(params.get("binds")),
                singularity_bin=str(params.get("singularity_bin","singularity")),
            )
            line1 = shlex.join(cmd)
            idx   = shlex.join(singularity_exec_cmd(
                        image=str(image),
                        argv=["bcftools","index","-f",out_vcf],
                        binds=normalize_binds(params.get("binds")),
                        singularity_bin=str(params.get("singularity_bin","singularity")),
                    ))
        else:
            line1 = shlex.join(argv)
            idx   = " ".join(["bcftools","index","-f",shlex.quote(out_vcf)])

        return [line1, idx]

    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"dir": str(self.workdir)},
            params=p,
            threads=int(self.threads or 1),
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id"),
            ensure_output_dir_key="dir",
        )



    def to_sh(self) -> List[str]:
        # threads는 의미 없지만 인터페이스 통일을 위해 params/threads 병합
        p: Dict[str, Any] = {**self.DEFAULTS, **(self.params or {})}
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"dir": str(self.workdir)},
            params=p,
            threads=int(self.threads or 1),
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id"),
            ensure_output_dir_key="dir",
        )


#!/usr/bin/env bash
set -euo pipefail

# ===== 사용자/워크플로 변수 =====
BAM="/abs/path/sample.recal.bam"
REF_FASTA="/abs/path/ref.fa"
REGION="chr20"                           # 혹은 "chr20:1000000-2000000"
OUT_DIR="/abs/work/S1/71_Beagle"
THREADS=12

# 로컬 바이너리
BIN_BCFTOOLS="bcftools"
BIN_BEAGLE="beagle"                      # 또는 'java -jar /path/beagle.5.4.jar'
MAP_FILE=""                              # optional (있으면 경로 입력)
REF_BREF3=""                             # optional (있으면 경로 입력)

# ===== 준비 =====
mkdir -p "${OUT_DIR}"/{gl,beagle}

# 1) GL 생성
SAFE_REGION="${REGION//:/_}"
GL_VCF="${OUT_DIR}/gl/gl.${SAFE_REGION}.vcf.gz"
if [[ ! -s "${GL_VCF}" ]]; then
  echo "[GL] ${REGION} → ${GL_VCF}"
  ${BIN_BCFTOOLS} mpileup -f "${REF_FASTA}" -r "${REGION}" -Ou "${BAM}" \
    | ${BIN_BCFTOOLS} call -Aim -Oz -o "${GL_VCF}"
  ${BIN_BCFTOOLS} index -f "${GL_VCF}"
fi

# 2) Beagle imputation
OUT_PREFIX="${OUT_DIR}/beagle/imputed.${SAFE_REGION}"
CMD=( ${BIN_BEAGLE} "gl=${GL_VCF}" "out=${OUT_PREFIX}" "nthreads=${THREADS}" )
[[ -n "${REF_BREF3}" ]] && CMD+=( "ref=${REF_BREF3}" )
[[ -n "${MAP_FILE}"  ]] && CMD+=( "map=${MAP_FILE}" )

echo "[BEAGLE] ${CMD[*]}"
"${CMD[@]}"

FINAL_VCF="${OUT_PREFIX}.vcf.gz"
echo "[DONE] Beagle imputation → ${FINAL_VCF}"