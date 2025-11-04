from __future__ import annotations
from typing import Dict, Any, List, Optional, Sequence
from pathlib import Path
import os

from src.tasks.task import Task
from src.tasks.task_registry import register_task
from src.tasks.utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    to_sh_from_builder,
)


@register_task("picard.fastqtosam")   # ← 여기서 키 등록
class PicardFastqToSamTask(Task):
    """
    Picard FastqToSam 실행 태스크

    INPUTS:
      read1: R1 FASTQ(.gz) absolute path (required)
      read2: R2 FASTQ(.gz) absolute path (required)

    OUTPUTS:
      bam:   output BAM path (optional; default: {workdir}/{sample_id}.fastqtosam.bam)
      dir:   (optional) 결과 디렉토리 지정 시 bam이 없을 때 기본 dir로 생성

    PARAMS:
      sample_name:         SAMPLE_NAME (default: sample_id 또는 read1에서 추정)
      read_group_id:       READ_GROUP_NAME
      platform:            PLATFORM
      library_name:        LIBRARY_NAME
      sequencing_center:   SEQUENCING_CENTER
      tmp_dir:             TMP_DIR (default: {workdir}/tmp)
      java_bin:            "java" 기본
      jar:                 "/storage/apps/bin/picard.jar" 기본
      xmx_gb:              16  → -Xmx16g
      parallel_gc_threads: 14  → -XX:ParallelGCThreads=14
      image:               Singularity 이미지(optional)
      binds:               Singularity 바인드(optional)
      singularity_bin:     "singularity" 기본
    """

    TYPE = "picard.fastqtosam"

    # (문서/검증용) 스키마
    INPUTS: Dict[str, Any] = {
        "read1": {"type": "path", "required": True,  "desc": "R1 FASTQ(.gz)"},
        "read2": {"type": "path", "required": True,  "desc": "R2 FASTQ(.gz)"},
    }
    OUTPUTS: Dict[str, Any] = {
        "bam": {"type": "path", "required": False, "desc": "Output BAM path"},
        "dir": {"type": "dir",  "required": False, "desc": "Base directory for outputs (optional)"},
    }

    DEFAULTS: Dict[str, Any] = {
        "java_bin": "java",
        "jar": "/storage/apps/bin/picard.jar",
        "xmx_gb": 16,
        "parallel_gc_threads": 14,
        "tmp_dir": None,  # None이면 workdir/tmp
        "sample_name": None,
        "read_group_id": None,
        "platform": None,
        "library_name": None,
        "sequencing_center": None,
        "image": None,
        "binds": None,
        "singularity_bin": "singularity",
    }

    # ---- 내부 헬퍼
    @staticmethod
    def _infer_sample_name(read1_path: str) -> str:
        """R1 파일명에서 샘플 prefix 추정 (S1_R1.fastq.gz → S1 등)."""
        base = os.path.basename(read1_path)
        for tag in ["_R1.fastq.gz", "_R1.fq.gz", ".R1.fastq.gz", ".R1.fq.gz",
                    "_R1.fastq", "_R1.fq", ".R1.fastq", ".R1.fq"]:
            if base.endswith(tag):
                return base[:-len(tag)]
        # 일반 확장자 제거
        for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"]:
            if base.endswith(ext):
                return base[:-len(ext)]
        return os.path.splitext(base)[0]

    # ---- 순수 커맨드 빌더
    def _build_cmd(
        self,
        *,
        inputs: Dict[str, Any],
        outputs: Dict[str, Any],
        params: Dict[str, Any],
        threads: int,          # not used directly by picard here, but kept for consistency
        workdir: str,
        sample_id: Optional[str] = None,
    ) -> List[Sequence[str] | str]:

        r1 = inputs.get("read1")
        r2 = inputs.get("read2")
        if not r1 or not r2:
            raise ValueError("[picard.fastqtosam] inputs.read1/read2 are required")

        # 출력 경로 결정
        base_dir = outputs.get("dir") or workdir
        base_dir = ensure_dir(base_dir)
        bam_out = outputs.get("bam") or os.path.join(base_dir, f"{sample_id or self._infer_sample_name(str(r1))}.fastqtosam.bam")

        # TMP 디렉토리
        tmp_dir = params.get("tmp_dir")
        if not tmp_dir:
            tmp_dir = os.path.join(base_dir, "tmp")
        tmp_dir = ensure_dir(tmp_dir)

        # 필드 정리
        sample_name = params.get("sample_name") or sample_id or self._infer_sample_name(str(r1))
        read_group_id = params.get("read_group_id")
        platform = params.get("platform")
        library_name = params.get("library_name")
        sequencing_center = params.get("sequencing_center")

        java_bin = str(params.get("java_bin", "java"))
        jar_path = str(params.get("jar", "/storage/apps/bin/picard.jar"))
        xmx_gb = int(params.get("xmx_gb", 16))
        pgc = int(params.get("parallel_gc_threads", 14))

        core = [
            java_bin,
            f"-XX:ParallelGCThreads={pgc}",
            f"-Xmx{xmx_gb}g",
            "-jar", jar_path,
            "FastqToSam",
            "--FASTQ", str(r1),
            "--FASTQ2", str(r2),
            "--SAMPLE_NAME", str(sample_name),
            "--OUTPUT", bam_out,
        ]

        # 선택적 RG/PL/LB/SC
        if read_group_id:
            core += ["--READ_GROUP_NAME", str(read_group_id)]
        if platform:
            core += ["--PLATFORM", str(platform)]
        if library_name:
            core += ["--LIBRARY_NAME", str(library_name)]
        if sequencing_center:
            core += ["--SEQUENCING_CENTER", str(sequencing_center)]
        if tmp_dir:
            core += ["--TMP_DIR", tmp_dir]

        # 컨테이너 래핑
        image = params.get("image")
        if image:
            cmd = singularity_exec_cmd(
                image=str(image),
                argv=core,
                binds=normalize_binds(params.get("binds")),
                singularity_bin=str(params.get("singularity_bin", "singularity")),
            )
            return [cmd]
        else:
            return [core]

    # ---- to_sh: 유틸로 보일러플레이트 처리
    def to_sh(self) -> List[str]:
        p: Dict[str, Any] = {**self.DEFAULTS, **(self.params or {})}
        th = int(self.threads or p.get("threads", 1))  # picard는 threads 직접 쓰지 않지만 인터페이스 통일
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"dir": str(self.workdir)},
            params=p,
            threads=th,
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id") or None,
            ensure_output_dir_key="dir",
        )