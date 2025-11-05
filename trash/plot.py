#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import entropy

def read_bed(path: str) -> pd.DataFrame:
    """BED( chr, start, end, depth ) 읽기. .gz도 자동 인식."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "depth"],
        usecols=[0, 1, 2, 3],
        dtype={0: str, 1: np.int64, 2: np.int64, 3: np.float64},
        compression="infer",
        engine="c",
    )
    # 유효값만
    df = df[np.isfinite(df["depth"])]
    return df

def align_windows(bulk: pd.DataFrame, sample: pd.DataFrame, chrom: str | None) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """공통 윈도우( chr,start,end ) 기준 inner join."""
    if chrom:
        bulk = bulk[bulk["chr"] == chrom]
        sample = sample[sample["chr"] == chrom]

    key_cols = ["chr", "start", "end"]
    merged = pd.merge(
        bulk[key_cols + ["depth"]].rename(columns={"depth": "bulk_depth"}),
        sample[key_cols + ["depth"]].rename(columns={"depth": "sample_depth"}),
        on=key_cols,
        how="inner",
        validate="one_to_one",
    ).sort_values(key_cols, kind="mergesort")

    if merged.empty:
        raise ValueError("공통 윈도우가 없습니다. (chr/start/end가 일치하는 행이 없음)")

    return merged["bulk_depth"].to_numpy(), merged["sample_depth"].to_numpy(), merged[key_cols]

def normalize_to_prob(arr: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """depth를 확률분포로 (합=1)."""
    x = np.clip(arr.astype(float), eps, None)
    s = x.sum()
    if s <= 0:
        raise ValueError("합계가 0입니다. depth 값을 확인하세요.")
    return x / s

def compute_kl_and_plot(bulk_depth: np.ndarray,
                        sample_depth: np.ndarray,
                        title_suffix: str = "",
                        bins: int = 30,
                        out_png: str | None = None,
                        log_base: float | None = None) -> float:
    """KL(sample‖bulk) 계산 후 히스토그램 플롯."""
    p = normalize_to_prob(sample_depth)
    q = normalize_to_prob(bulk_depth)

    kl_val = float(entropy(p, q, base=log_base))  # base=None → nats, base=2 → bits

    # 플롯 (규정: matplotlib, 한 그림당 하나의 차트, 색상 미지정)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(p, bins=bins, alpha=0.6, density=True, label="Sample (normalized)")
    ax.hist(q, bins=bins, alpha=0.6, density=True, label="Bulk (normalized)")
    ax.set_xlabel("Normalized 100 kb window read depth")
    ax.set_ylabel("Density")
    unit = "bits" if log_base == 2 else ("nats" if log_base is None else f"log base {log_base}")
    ax.set_title(f"Normalized depth distribution — KL(sample‖bulk) = {kl_val:.4g} {unit} {title_suffix}".strip())
    ax.legend()
    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
    else:
        plt.show()
    plt.close(fig)

    return kl_val

def main():
    ap = argparse.ArgumentParser(description="BED 기반 bulk vs sample KL(sample‖bulk) + 히스토그램 플롯")
    ap.add_argument("--bulk", required=True, help="bulk BED(.bed/.bed.gz)")
    ap.add_argument("--sample", required=True, help="sample BED(.bed/.bed.gz)")
    ap.add_argument("--chrom", default=None, help="특정 염색체만 계산 (예: chr1). 생략 시 전체 공통 윈도우 사용")
    ap.add_argument("--bins", type=int, default=100000, help="히스토그램 bin 수 (기본 30)")
    ap.add_argument("--out", default=None, help="결과 그림 PNG 경로 (미지정 시 화면 표시)")
    ap.add_argument("--log2", action="store_true", help="KL 단위를 bits로(로그 밑 2)")
    args = ap.parse_args()

    bulk = read_bed(args.bulk)
    sample = read_bed(args.sample)

    bulk_depth, sample_depth, keys = align_windows(bulk, sample, args.chrom)
    suffix = f"({args.chrom})" if args.chrom else ""

    kl_val = compute_kl_and_plot(
        bulk_depth=bulk_depth,
        sample_depth=sample_depth,
        title_suffix=suffix,
        bins=args.bins,
        out_png=args.out,
        log_base=2 if args.log2 else None,
    )

    # 요약 표준 출력
    print(f"Windows used: {len(bulk_depth)}")
    if args.chrom:
        print(f"Chromosome: {args.chrom}")
    print(f"KL(sample‖bulk) = {kl_val:.6g} {'bits' if args.log2 else 'nats'}")
    if args.out:
        print(f"Figure saved to: {args.out}")

if __name__ == "__main__":
    main()
