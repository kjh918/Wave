from __future__ import annotations
from pathlib import Path
from typing import Iterable, List, Optional
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import shlex
import matplotlib.pyplot as plt


def _cat_cmd_for(p: Path) -> str:
    s = str(p)
    if s.endswith(".gz"):
        return f"zcat {shlex.quote(s)}"
    if s.endswith(".bz2"):
        return f"bzcat {shlex.quote(s)}"
    return f"cat {shlex.quote(s)}"


def _find_pair(raw_dir: Path, seqid: str) -> List[Path]:
    # 흔한 네이밍: *_R1*.fastq* / *_R2*.fastq*
    r1s = sorted(raw_dir.glob(f"*{seqid}*R1*.fastq*"))
    r2s = sorted(raw_dir.glob(f"*{seqid}*R2*.fastq*"))
    if not r1s:
        raise FileNotFoundError(f"[jellyfish] R1 not found for {seqid} under {raw_dir}")
    if not r2s:
        raise FileNotFoundError(f"[jellyfish] R2 not found for {seqid} under {raw_dir}")
    return [r1s[0], r2s[0]]

def build_jellyfish_count_cmd(
        *,
        # 입력: 둘 중 하나 사용
        inputs: Optional[List[str]] = None,      # 명시 입력(절대/상대 경로 가능)
        RawFastqDir: Optional[str] = None,       # 자동 탐색용
        SeqID: Optional[str] = None,             # 자동 탐색용
        # 옵션
        jellyfish: str = '/storage/home/jhkim/Apps/jellyfish-2.3.1/bin/jellyfish',
        k_mer: int = 21,
        size: str = "100M",
        threads: int = 8,
        canonical: bool = True,                  # -C
        outDir: str,
        outPrefix: str,
        image: Optional[str] = None,             # Singularity 이미지(.sif), 없으면 로컬 실행
        binds: Optional[List[str] | str] = None, # Singularity -B
    ) -> Iterable[str]:
    """
    jellyfish count 커맨드 생성.
    - 입력이 주어지지 않으면 RawFastqDir/SeqID 기반으로 R1/R2 자동 탐색
    - process substitution(<(... )>) 사용 → /bin/bash -lc 로 래핑
    - out: {outDir}/{outPrefix}.jf
    """
    out_dir = Path(outDir)
    out_dir.mkdir(parents=True, exist_ok=True)
    jf_out = out_dir / f"{outPrefix}.jf"

    # 입력 결정
    input_paths: List[Path]
    if inputs:
        input_paths = [Path(x) for x in inputs]
    else:
        if not RawFastqDir or not SeqID:
            raise ValueError("[jellyfish] Either 'inputs' or both 'RawFastqDir' and 'SeqID' must be provided")
        input_paths = _find_pair(Path(RawFastqDir), SeqID)

    # 각 입력을 process substitution으로 감싸기
    proc_inputs: List[str] = [f" <({_cat_cmd_for(p)})" for p in input_paths]
    proc_inputs_str = " ".join(proc_inputs)

    # jellyfish count 옵션 구성
    opts = [
        f"-m {k_mer}",
        f"-s {size}",
        f"-t {threads}",
        f"-o {shlex.quote(str(jf_out))}",
    ]
    if canonical:
        opts.append("-C")

    # 실제 jellyfish 커맨드 (bash 필요)
    jf_cmd = f"{jellyfish} count {' '.join(opts)} {proc_inputs_str}"

    # Singularity 래핑
    if image:
        if isinstance(binds, list):
            bind_opt = "-B " + ",".join(binds)
        elif isinstance(binds, str):
            bind_opt = f"-B {binds}" if binds else ""
        else:
            bind_opt = ""
        cmd = f"singularity exec {bind_opt} {shlex.quote(image)} {bash_cmd}"
    else:
        cmd = jf_cmd

    return [
        cmd
    ]

def plot_histogram(df: pd.DataFrame, y_percent: np.ndarray, P: float | None, peak_B: float | None, title_suffix: str = ""):
    plt.figure(figsize=(10, 6))
    plt.plot(df["B"].values, y_percent, lw=1.5, label="Histogram (percent)")
    if peak_B is not None:
        plt.axvline(peak_B, linestyle=":", linewidth=2, label=f"Single-copy peak ≈ {peak_B:.1f}")
    if P is not None:
        plt.axhline(0, alpha=0.0)  # just to keep style consistent
    plt.xlabel("k-mer bin size (B)")
    plt.ylabel("Frequency (%)")
    plt.title(f"k-mer Frequency Distribution (Jellyfish histo){title_suffix}")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.show()


def build_draw_histo_plot_cmd(
        *,
        # 입력: 둘 중 하나 사용
        inputs: Optional[List[str]] = None,      # 명시 입력(절대/상대 경로 가능)
        RawFastqDir: Optional[str] = None,       # 자동 탐색용
        SeqID: Optional[str] = None,             # 자동 탐색용
        # 옵션
        jellyfish: str = '/storage/home/jhkim/Apps/jellyfish-2.3.1/bin/jellyfish',
        k_mer: int = 21,
        size: str = "100M",
        threads: int = 8,
        canonical: bool = True,                  # -C
        outDir: str,
        outPrefix: str,
        min_bin:int = 1,
        manual_P = 4,
        image: Optional[str] = None,             # Singularity 이미지(.sif), 없으면 로컬 실행
        binds: Optional[List[str] | str] = None, # Singularity -B
    ) -> Iterable[str]:
    """
    jellyfish count 커맨드 생성.
    - 입력이 주어지지 않으면 RawFastqDir/SeqID 기반으로 R1/R2 자동 탐색
    - process substitution(<(... )>) 사용 → /bin/bash -lc 로 래핑
    - out: {outDir}/{outPrefix}.jf
    """
    rc_correction = False 

    out_dir = Path(outDir)
    out_dir.mkdir(parents=True, exist_ok=True)
    input_jf = out_dir / f"{outPrefix}.jf"
    out_histo = out_dir / f"{outPrefix}.histo"
    
    out_histo_plot = out_dir / f"{outPrefix}.histo.png"
    out_log_txt = out_dir / f"{outPrefix}.histo.log.txt"
    
    # --- 1) 읽기 ---
    df = pd.read_csv(out_histo, sep=' ', header=None, names=["B","F"])
    df = df[(df["B"]>0) & (df["F"]>0)].copy()
    df.sort_values("B", inplace=True)

    # --- 2) 그림용 퍼센트 (계산엔 사용 안 함) ---
    df["percent"] = df["F"] / df["F"].sum() * 100

    # --- 3) 컷오프 적용해 계산용 테이블 생성 ---
    calc = df[df["B"] >= min_bin].copy()

    # --- 4) P(피크 깊이) 결정 ---
    if manual_P is not None:
        P = float(manual_P)
    else:
        # 에러 구간 이후 첫 피크를 자동 탐지 (퍼센트 기반, prominence 살짝 사용)
        y = calc["percent"].to_numpy()
        peaks, _ = find_peaks(y, prominence=max(y)*0.02 if len(y) else 0.0, distance=3)
        if len(peaks) == 0:
            raise RuntimeError("True peak을 찾지 못했습니다. manual_P를 지정하세요.")
        # coverage가 가장 작은 피크를 선택
        idx = peaks[np.argmin(calc.iloc[peaks]["B"].values)]
        P = float(calc.iloc[idx]["B"])

    # --- 5) Σ(B·F) 계산 (원시 F 사용) ---
    numerator = float((calc["B"] * calc["F"]).sum())
    if rc_correction:
        numerator /= 2.0  # -C 미사용 보정 (근사)

    E_bp = numerator / P  # E = Σ(B·F) / P

    raw_numer = float((df["B"] * df["F"]).sum())
    
    with open(out_log_txt, 'w') as handle:
        handle.write(f"# SampleID : {SeqID}\n")
        handle.write(f"# Histo_Path : {out_histo}\n")
        handle.write("# === Genome size estimation ===\n")
        handle.write(f"Total Σ(Bin * Count) (raw)      : {raw_numer:,.0f}\n")
        handle.write(f"Total Σ(Bin * Count) (filtered) : {numerator:,.0f}   (min_bin >= {min_bin}{', RC/2' if rc_correction else ''})\n")
        handle.write(f"Peak depth P            : {P:.2f}\n")
        handle.write(f"Estimated genome size E : {E_bp/1e9:.3f} Gb  (E = Σ(B·F)/P)")
    # --- 6) 리포트 ---
    plt.figure(figsize=(10,6))
    plt.plot(df["B"], df["percent"], lw=1.5, color='steelblue', label='k-mer frequency (%)')
    plt.axvline(P, color='red', linestyle='--', lw=1.5, label=f'Peak depth (P={P})')
    plt.xlabel("k-mer coverage (B)")
    plt.ylabel("Frequency (%)")
    plt.title("k-mer Frequency Distribution")
    plt.xlim([0,50])
    # plt.xticks(list(range(0,5,55)))
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(out_histo_plot, dpi=300)
    plt.clf()

    return [
        1
    ]

    
def build_jellyfish_histo_cmd(
        *,
        # 입력: 둘 중 하나 사용
        inputs: Optional[List[str]] = None,      # 명시 입력(절대/상대 경로 가능)
        RawFastqDir: Optional[str] = None,       # 자동 탐색용
        SeqID: Optional[str] = None,             # 자동 탐색용
        # 옵션
        jellyfish: str = '/storage/home/jhkim/Apps/jellyfish-2.3.1/bin/jellyfish',
        k_mer: int = 21,
        size: str = "100M",
        threads: int = 8,
        canonical: bool = True,                  # -C
        outDir: str,
        outPrefix: str,
        image: Optional[str] = None,             # Singularity 이미지(.sif), 없으면 로컬 실행
        binds: Optional[List[str] | str] = None, # Singularity -B
    ) -> Iterable[str]:
    """
    jellyfish count 커맨드 생성.
    - 입력이 주어지지 않으면 RawFastqDir/SeqID 기반으로 R1/R2 자동 탐색
    - process substitution(<(... )>) 사용 → /bin/bash -lc 로 래핑
    - out: {outDir}/{outPrefix}.jf
    """
    out_dir = Path(outDir)
    out_dir.mkdir(parents=True, exist_ok=True)
    input_jf = out_dir / f"{outPrefix}.jf"
    out_histo = out_dir / f"{outPrefix}.histo"
    
    out_histo_plot = out_dir / f"{outPrefix}.histo.png"


    # 실제 jellyfish 커맨드 (bash{ 필요)
    jf_cmd = f"{jellyfish} histo {input_jf} > {out_histo}"

    # Singularity 래핑
    if image:
        if isinstance(binds, list):
            bind_opt = "-B " + ",".join(binds)
        elif isinstance(binds, str):
            bind_opt = f"-B {binds}" if binds else ""
        else:
            bind_opt = ""
        cmd = f"singularity exec {bind_opt} {shlex.quote(image)} {bash_cmd}"
    else:
        cmd = jf_cmd

    return [
        cmd
    ]