# tasks/fastp/_func.py
from __future__ import annotations
from pathlib import Path
from typing import Iterable, Optional, List
import shlex

def build_fastp_cmd(
        RawFastqDir: str,
        SeqID: str,
        TrimFastqDir: str,
        threads: int = 4,
        picoplex_gold: bool = False,          # True면 --trim_front1/2 14 적용
        trim_front1: int = 14,                # picoplex_gold일 때 기본값
        trim_front2: int = 14,
        length_required: int = 100,
        average_qual: int = 10,
        qualified_quality_phred: int = 15,
        image: str = "/storage/images/fastp-0.23.4.sif",  # None이면 로컬 fastp 사용
        binds: Optional[List[str] | str] = ("/storage,/data"),
    ) -> Iterable[str]:
    """
    fastp 명령을 작성하여 반환한다 (Iterable[str]).
    쉘 스니펫의 동작을 그대로 재현:
      - 공통: --trim_poly_g, --detect_adapter_for_pe, 어댑터 시퀀스, 길이/품질 필터
      - PicoplexGold=True: --trim_front1/2 14 추가
    """
    raw = Path(RawFastqDir)
    if not raw.exists():
        raise FileNotFoundError(f"[fastp] RawFastqDir not found: {raw}")

    r1 = raw / f"{SeqID}_R1.fastq.gz"
    r2 = raw / f"{SeqID}_R2.fastq.gz"
    if not r1.exists():
        raise FileNotFoundError(f"[fastp] R1 not found: {r1}")
    if not r2.exists():
        raise FileNotFoundError(f"[fastp] R2 not found: {r2}")

    out_dir = Path(TrimFastqDir)
    out1 = out_dir / f"{SeqID}.trimmed_R1.fastq.gz"
    out2 = out_dir / f"{SeqID}.trimmed_R2.fastq.gz"
    json_out = out_dir / f"{SeqID}.fastp.json"
    html_out = out_dir / f"{SeqID}.fastp.html"

    # 공통 fastp 옵션
    opts = [
        f"--thread {threads}",
        f"--in1 {shlex.quote(str(r1))}",
        f"--in2 {shlex.quote(str(r2))}",
        f"--out1 {shlex.quote(str(out1))}",
        f"--out2 {shlex.quote(str(out2))}",
        f"--json {shlex.quote(str(json_out))}",
        f"--html {shlex.quote(str(html_out))}",
        "--trim_poly_g",
        "--detect_adapter_for_pe",
        '--adapter_sequence "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"',
        '--adapter_sequence_r2 "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"',
        f"--length_required {length_required}",
        f"--average_qual {average_qual}",
        f"--qualified_quality_phred {qualified_quality_phred}",
    ]

    # PicoplexGold일 때 전처리 길이 트리밍
    if picoplex_gold:
        opts += [f"--trim_front1 {trim_front1}", f"--trim_front2 {trim_front2}"]

    base_cmd = "fastp " + " ".join(opts)

    # Singularity 래핑
    if image:
        if isinstance(binds, list):
            bind_opt = "-B " + ",".join(binds)
        else:
            bind_opt = f"-B {binds}" if binds else ""
        cmd = f"singularity exec {bind_opt} {image} {base_cmd}"
    else:
        cmd = base_cmd
        
    return [cmd]
