from __future__ import annotations
from typing import Iterable, Dict, Any, Optional, List
from pathlib import Path

from tasks.task import Task
from tasks.task_registry import register_task
from ._func import build_jellyfish_count_cmd


@register_task("jellyfish_count")
class JellyfishCountTask(Task):
    """
    jellyfish count Task

    파라미터:
      - (입력 둘 중 하나 사용)
        * inputs: [R1, R2, ...]          # 직접 지정
        * RawFastqDir + SeqID            # 자동 탐색 (R1/R2)
      - k_mer: int                       # k-mer length (default 21)
      - size: str                        # hash size, e.g. "100M"
      - threads: int                     # CPU threads
      - canonical: bool                  # -C (both strands canonical)
      - outDir: str                      # 출력 디렉토리 (기본: workdir/<SeqID>_jellyfish)
      - outPrefix: str                   # 결과 파일 접두어 (기본: SeqID)
      - image: str (optional)            # Singularity 이미지 경로
      - binds: list[str] | str (optional)# Singularity -B
    """

    TYPE = "jellyfish_count"

    INPUTS = {
        "inputs": "list[str] optional; direct sequence files",
        "RawFastqDir": "Input FASTQ dir (auto-detect R1/R2 when inputs omitted)",
        "SeqID": "Sample ID",
        "k_mer": "int k-mer length",
        "size": "str hash size (e.g., 100M)",
        "threads": "int",
        "canonical": "bool",
        "outDir": "output directory",
        "outPrefix": "output prefix (without extension)",
        "image": "optional singularity image",
        "binds": "optional bind list or comma-string",
    }

    DEFAULTS = {
        "k_mer": 21,
        "size": "100M",
        "threads": 8,
        "canonical": True,
    }

    OPTIONAL = {"inputs", "RawFastqDir", "SeqID", "image", "binds", "outDir", "outPrefix"}

    def to_sh(self) -> Iterable[str]:
        p: Dict[str, Any] = self.params

        seqid = str(p.get("SeqID", self.params.get("sample_id", "sample")))
        out_dir = p.get("outDir") or str(Path(self.workdir) / f"{seqid}_jellyfish")
        out_prefix = p.get("outPrefix") or seqid

        # inputs가 주어졌으면 그대로 사용, 아니면 RawFastqDir/SeqID로 자동 탐색
        inputs: Optional[List[str]] = p.get("inputs")

        return build_jellyfish_count_cmd(
            inputs=inputs,
            RawFastqDir=p.get("RawFastqDir"),
            SeqID=seqid if not inputs else None,  # inputs가 있으면 자동탐색 미사용
            k_mer=int(p.get("k_mer", self.DEFAULTS["k_mer"])),
            size=str(p.get("size", self.DEFAULTS["size"])),
            threads=int(p.get("threads", self.DEFAULTS["threads"])),
            canonical=bool(p.get("canonical", self.DEFAULTS["canonical"])),
            outDir=out_dir,
            outPrefix=str(out_prefix),
            image=p.get("image"),
            binds=p.get("binds"),
        )
