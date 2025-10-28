# tasks/fastp/_func.py
from __future__ import annotations
from pathlib import Path
from typing import Iterable, Optional, Sequence, List
import shlex


def _norm_binds(binds: Optional[Sequence[str] | str]) -> str:
    """
    binds가 리스트면 콤마-조인하여 단일 -B 인자로 사용,
    문자열이면 그대로 사용. 없으면 빈 문자열 반환.
    """
    if binds is None:
        return ""
    if isinstance(binds, (list, tuple)):
        joined = ",".join(str(b) for b in binds if str(b).strip())
        return f"-B {joined}" if joined else ""
    s = str(binds).strip()
    return f"-B {s}" if s else ""


def build_fastp_cmd(
    *,
    RawFastqDir: str,
    SeqID: str,
    TrimFastqDir: str,
    threads: int = 4,
    picoplex_gold: bool = False,     # True면 --trim_front1/2 적용
    trim_front1: int = 14,
    trim_front2: int = 14,
    length_required: int = 100,
    average_qual: int = 10,
    qualified_quality_phred: int = 15,
    image: Optional[str] = None,     # None이면 로컬 fastp 사용, 경로 주면 singularity exec
    binds: Optional[Sequence[str] | str] = None,
) -> Iterable[str]:
    """
    fastp 실행 커맨드를 문자열 한 줄로 만들어 리스트에 담아 반환합니다.
    상위 레이어에서 그대로 파일에 쓰거나 subprocess로 실행하기 좋습니다.
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
    # 상위에서 디렉토리를 만들어 주지만, 안전하게 경로만 계산
    out1 = out_dir / f"{SeqID}.trimmed_R1.fastq.gz"
    out2 = out_dir / f"{SeqID}.trimmed_R2.fastq.gz"
    json_out = out_dir / f"{SeqID}.fastp.json"
    html_out = out_dir / f"{SeqID}.fastp.html"

    # 공통 fastp 옵션
    opts: List[str] = [
        "--thread", str(int(threads)),
        "--in1", shlex.quote(str(r1)),
        "--in2", shlex.quote(str(r2)),
        "--out1", shlex.quote(str(out1)),
        "--out2", shlex.quote(str(out2)),
        "--json", shlex.quote(str(json_out)),
        "--html", shlex.quote(str(html_out)),
        "--trim_poly_g",
        "--detect_adapter_for_pe",
        "--adapter_sequence", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        "--adapter_sequence_r2", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "--length_required", str(int(length_required)),
        "--average_qual", str(int(average_qual)),
        "--qualified_quality_phred", str(int(qualified_quality_phred)),
    ]

    # PicoplexGold일 때 read 앞부분 트리밍
    if picoplex_gold:
        opts += ["--trim_front1", str(int(trim_front1)),
                 "--trim_front2", str(int(trim_front2))]

    if image:
        # Singularity 래핑
        bflag = _norm_binds(binds)
        # singularity exec [-B ...] <image> fastp <opts...>
        tokens = ["singularity", "exec"]
        if bflag:
            tokens += shlex.split(bflag)
        tokens += [image, "fastp"] + opts
    else:
        # 로컬 fastp
        tokens = ["fastp"] + opts

    # 문자열 한 줄로 반환 (상위가 그대로 파일/실행에 사용)
    line = " ".join(shlex.quote(t) if " " in t and not t.startswith("--") else t for t in tokens)
    return [line]