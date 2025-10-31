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
        "trim_front1": 0,
        "trim_front2": 0,
        "length_required": 100,
        "average_qual": 10,
        "qualified_quality_phred": 15,
        "adapter_sequence": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        "adapter_sequence_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "image": "/storage/images/fastp-0.23.4.sif",
        "binds": ["/storage","/data"]
    }

    def to_sh(self) -> Iterable[str]:
        p: Dict[str, Any] = self.params

        # 출력 디렉토리 자동 설정
        trim_dir: str = self.workdir 

        Path(trim_dir).mkdir(parents=True, exist_ok=True)

        # PicoplexGold/trim 프론트 설정
        picoplex_gold: bool = bool(p.get("PicoplexGold", self.DEFAULTS["PicoplexGold"]))
        tf1_default = int(self.DEFAULTS["trim_front1"])
        tf2_default = int(self.DEFAULTS["trim_front2"])
        trim_front1 = int(p.get("trim_front1", tf1_default))
        trim_front2 = int(p.get("trim_front2", tf2_default))

        # _func 빌더 호출 (이 함수는 argv 리스트 또는 문자열 라인 리스트를 반환할 수 있음)
        lines = build_fastp_cmd(
            read1 = self.inputs.get("read1"),
            read2 = self.inputs.get("read2"),
            trim_read1 = self.outputs.get("read1"),
            trim_read2 = self.outputs.get("read2"),
            fastp_json = self.outputs.get("json"),
            fastp_html = self.outputs.get("html"),
            threads=int(p.get("threads", self.DEFAULTS["threads"])),
            picoplex_gold=picoplex_gold,
            trim_front1=trim_front1,
            trim_front2=trim_front2,    
            adapter_sequence = p.get("adapter_sequence", self.DEFAULTS["adapter_sequence"]),
            adapter_sequence_r2 = p.get("adapter_sequence_r2", self.DEFAULTS["adapter_sequence_r2"]),
            length_required=int(p.get("length_required", self.DEFAULTS["length_required"])),
            average_qual=int(p.get("average_qual", self.DEFAULTS["average_qual"])),
            qualified_quality_phred=int(p.get("qualified_quality_phred", self.DEFAULTS["qualified_quality_phred"])),
            image=p.get("image"),
            binds=p.get("binds"),
        )
        return lines
