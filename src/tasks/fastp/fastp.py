# tasks/fastp/fastp.py
from __future__ import annotations
from typing import Iterable, Dict, Any, Sequence
from pathlib import Path
import shlex

from src.tasks.task import Task, TaskRegistry
from ._func import build_fastp_cmd


@TaskRegistry.register
class FastpTask(Task):
    """
    Paired-end fastp 실행 Task

    Params (self.params):
      - RawFastqDir:   입력 FASTQ 디렉토리 (필수)
      - SeqID:         샘플 ID (필수; <SeqID>_R1/2.fastq.gz 형태를 탐색)
      - threads:       스레드 수 (기본 4)
      - PicoplexGold:  Picoplex Gold 프로토콜 여부 (bool). True면 --trim_front1/2 14 기본 적용
      - trim_front1/2: PicoplexGold=True일 때 기본 14. 직접 지정하면 우선
      - length_required / average_qual / qualified_quality_phred: 품질/길이 필터
      - TrimFastqDir:  트리밍 FASTQ 출력 디렉토리 (옵션; 없으면 workdir/<SeqID>_fastp)
      - image:         Singularity 이미지 경로 (옵션; 없으면 로컬 fastp 사용)
      - binds:         Singularity 바인드 리스트 또는 콤마 문자열 (옵션)
    """

    TYPE = "fastp"

    INPUTS = {
        "threads": "int",
        "PicoplexGold": "bool",
        "trim_front1": "int (optional; default 14 when PicoplexGold)",
        "trim_front2": "int (optional; default 14 when PicoplexGold)",
        "length_required": "int",
        "average_qual": "int",
        "qualified_quality_phred": "int",
        "TrimFastqDir": "Output dir for trimmed FASTQs (optional)",
        "image": "(optional) Singularity image path",
        "binds": "(optional) bind list or comma-string",
    }

    DEFAULTS = {
        "threads": 4,
        "PicoplexGold": False,
        "trim_front1": 14,
        "trim_front2": 14,
        "length_required": 100,
        "average_qual": 10,
        "qualified_quality_phred": 15,
    }

    def to_sh(self) -> Iterable[str]:
        p: Dict[str, Any] = self.Params


        # 출력 디렉토리 자동 설정
        trim_dir: str = p.get("TrimFastqDir") or str(Path(self.workdir) / f"{seqid}_fastp")
        Path(trim_dir).mkdir(parents=True, exist_ok=True)

        # PicoplexGold/trim 프론트 설정
        picoplex_gold: bool = bool(p.get("PicoplexGold", self.DEFAULTS["PicoplexGold"]))
        tf1_default = int(self.DEFAULTS["trim_front1"])
        tf2_default = int(self.DEFAULTS["trim_front2"])
        trim_front1 = int(p.get("trim_front1", tf1_default))
        trim_front2 = int(p.get("trim_front2", tf2_default))

        # binds 정규화
        binds = self._normalize_binds(p.get("binds"))

        # _func 빌더 호출 (이 함수는 argv 리스트 또는 문자열 라인 리스트를 반환할 수 있음)
        lines = build_fastp_cmd(
            RawFastqDir=str(p["RawFastqDir"]),
            SeqID=seqid,
            TrimFastqDir=trim_dir,
            threads=int(p.get("threads", self.DEFAULTS["threads"])),
            picoplex_gold=picoplex_gold,
            trim_front1=trim_front1,
            trim_front2=trim_front2,
            length_required=int(p.get("length_required", self.DEFAULTS["length_required"])),
            average_qual=int(p.get("average_qual", self.DEFAULTS["average_qual"])),
            qualified_quality_phred=int(p.get("qualified_quality_phred", self.DEFAULTS["qualified_quality_phred"])),
            image=p.get("image"),
            binds=binds,
        )

        # 문자열 라인으로 통일하여 반환
        return list(self._join_lines(lines))
    def to_sh(self) -> List[str]:
        inputs = [self.inputs.get("read1")]
        if self.inputs.get("read2"):
            inputs.append(self.inputs["read2"])

        out_dir = self.WORK_DIR
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        p = self.PARAMS
        return build_fastqc_cmd(
            inputs=inputs,
            out_dir=out_dir,
            threads=int(p["threads"]),
            extract=bool(p["extract"]),
            image=p.get("image"),
            binds=p.get("binds"),
            fastqc_bin=p.get("fastqc_bin", "fastqc"),
            singularity_bin=p.get("singularity_bin", "singularity"),
        )