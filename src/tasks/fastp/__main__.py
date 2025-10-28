# tasks/fastp/fastp.py
from __future__ import annotations
from typing import Iterable, Dict, Any, Optional
from pathlib import Path

import sys
import os 
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from task_registry import TaskRegistry
from task import Task
from ._func import build_fastp_cmd

TASK_CLASS = 'FastpTask'

# @register_task("fastp")
class FastpTask(Task):
    """
    Paired-end fastp 실행 Task

    필수/옵션 파라미터 정의:
      - RawFastqDir:   입력 FASTQ 디렉토리 (필수)
      - SeqID:         샘플 ID (필수; <SeqID>_R1/2.fastq.gz 형태를 탐색)
      - threads:       스레드 수 (기본 4)
      - PicoplexGold:  Picoplex Gold 프로토콜 여부(Yes/No 또는 bool). Yes면 --trim_front1/2 14 적용
      - trim_front1/2: PicoplexGold=True일 때 기본 14. 직접 지정 가능
      - length_required / average_qual / qualified_quality_phred: 품질/길이 필터
      - TrimFastqDir:  트리밍 FASTQ 출력 디렉토리 (옵션; 없으면 workdir/<SeqID>_fastp)
      - image:         Singularity 이미지 경로 (옵션; 없으면 로컬 fastp 사용)
      - binds:         Singularity 바인드 리스트 또는 콤마 문자열 (옵션)
    """

    TYPE = "fastp"

    INPUTS = {
        "RawFastqDir": "Input FASTQ dir",
        "SeqID": "Sample ID",
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

    OPTIONAL = {
        "TrimFastqDir",
        "image",
        "binds",
        "trim_front1",
        "trim_front2",
    }

    def to_sh(self) -> Iterable[str]:
        p: Dict[str, Any] = self.params
        seqid = str(p["SeqID"])

        # 출력 디렉토리 자동 설정
        trim_dir: str = p.get("TrimFastqDir") or str(Path(self.workdir) / f"{seqid}_fastp")
        Path(trim_dir).mkdir(parents=True, exist_ok=True)
        
        # PicoplexGold 여부 정규화
        picoplex_gold: bool = p.get("PicoplexGold")

        # trim_front는 PicoplexGold일 때만 강제 적용(사용자가 직접 값 넣으면 그 값을 우선)
        tf1_default = int(self.DEFAULTS["trim_front1"])
        tf2_default = int(self.DEFAULTS["trim_front2"])
        trim_front1 = int(p.get("trim_front1", tf1_default))
        trim_front2 = int(p.get("trim_front2", tf2_default))

        return build_fastp_cmd(
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
            binds=p.get("binds"),
        )
