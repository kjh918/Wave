#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from glob import glob
from itertools import combinations
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# ======================
# 파라미터 및 유틸
# ======================

def parse_args():
    p = argparse.ArgumentParser(
        description="Compare normalized depth≥threshold patterns for CASE vs CONTROL per chromosome."
    )
    p.add_argument("--case-glob", required=True,
                   help="Glob template for CASE files. Use '{chrom}' token for chromosome. "
                        "e.g. '/data/case/*/00_QC_BAM/*depth.{chrom}.txt'")
    p.add_argument("--control-glob", required=True,
                   help="Glob template for CONTROL files. Use '{chrom}' token. "
                        "e.g. '/data/control/*/00_QC_BAM/*depth.{chrom}.txt'")
    p.add_argument("-o", "--outdir", default="depth_pattern_case_control_out",
                   help="Output directory")
    p.add_argument("--chroms", default="1-22",
                   help="Chromosomes: '1-22' or '1,2,3,X,Y'")
    p.add_argument("--min-depth", type=int, default=2,
                   help="Depth threshold (default: 2)")
    p.add_argument("--norm-scale", type=float, default=1e8,
                   help="Normalization scale (default: 1e8)")
    p.add_argument("--bin-size", type=int, default=1000,
                   help="Uniform bin size for reindex (default: 1000bp)")
    p.add_argument("--png-dpi", type=int, default=200,
                   help="DPI for figures (default: 200)")
    p.add_argument("--no-pca", action="store_true",
                   help="Skip PCA")
    p.add_argument("--include-empty", action="store_true",
                   help="Keep positions even if all-NaN/0 after merge")
    p.add_argument("--smooth-window", type=int, default=0,
                   help="If >0, apply rolling mean to group-average curves (window bins)")
    return p.parse_args()

def resolve_chroms(spec: str):
    spec = spec.strip()
    chroms = []
    if "-" in spec and "," not in spec:
        a, b = spec.split("-", 1)
        a = a.strip(); b = b.strip()
        if a.isdigit() and b.isdigit():
            chroms = [str(i) for i in range(int(a), int(b) + 1)]
        else:
            raise ValueError(f"Invalid chrom spec: {spec}")
    else:
        parts = [x.strip() for x in spec.split(",")]
        for x in parts:
            if x.isdigit() or x in ["X", "Y", "M"]:
                chroms.append(x)
            else:
                raise ValueError(f"Invalid chromosome token: {x}")
    return [f"chr{c}" for c in chroms]

# ======================
# 데이터 로드 / 정규화
# ======================

def read_depth_table(path, min_depth, norm_scale):
    # samtools depth: chrom, pos, depth (TSV, no header)
    df = pd.read_csv(path, sep="\t", header=None)
    if df.shape[1] < 3:
        raise ValueError(f"Unexpected columns in {path}")
    df.columns = ["chrom", "pos", "depth"]

    total_depth = df["depth"].sum()
    if total_depth <= 0:
        return pd.DataFrame(), set()

    # filter
    df = df[df["depth"] >= min_depth].copy()
    if df.empty:
        sample = os.path.basename(path).split(".")[0]
        return pd.DataFrame(), set()

    # normalize
    df["normalized_dp"] = df["depth"] / total_depth * norm_scale
    df.index = df["pos"].astype(int)

    sample_name = os.path.basename(path).split(".")[0]
    out_df = df[["normalized_dp"]].rename(columns={"normalized_dp": sample_name})
    present_positions = set(df.index.tolist())
    return out_df, present_positions

def build_matrix_from_glob(glob_template, chrom, min_depth, norm_scale):
    paths = sorted(glob(glob_template.format(chrom=chrom)))
    merged = []
    pos_sets = {}
    for p in paths:
        try:
            df, pos_set = read_depth_table(p, min_depth, norm_scale)
        except Exception as e:
            print(f"Skip {p}: {e}", file=sys.stderr)
            continue
        if df.empty:
            sample = os.path.basename(p).split(".")[0]
            pos_sets[sample] = set()
            continue
        merged.append(df)
        pos_sets[df.columns[0]] = pos_set

    if merged:
        mat = pd.concat(merged, axis=1, join="outer")
    else:
        mat = pd.DataFrame()

    return mat, pos_sets

# ======================
# 결측 구간 시각화용 그리드
# ======================

def reindex_to_grid(df: pd.DataFrame, bin_size: int = 1000):
    if df.empty:
        return df
    start = int(df.index.min() // bin_size * bin_size)
    end = int(np.ceil(df.index.max() / bin_size) * bin_size)
    grid = np.arange(start, end + 1, bin_size, dtype=int)
    return df.reindex(grid)

def _contiguous_spans(mask: np.ndarray):
    if mask.size == 0:
        return []
    diff = np.diff(mask.astype(int))
    starts = np.where(diff == 1)[0] + 1
    ends   = np.where(diff == -1)[0] + 1
    if mask[0]:
        starts = np.r_[0, starts]
    if mask[-1]:
        ends = np.r_[ends, mask.size]
    return list(zip(starts, ends))  # [start, end)

def shade_missing_spans(ax, xs: np.ndarray, nan_mask: np.ndarray,
                        color="#d3d3d3", alpha=0.45):
    for s, e in _contiguous_spans(nan_mask):
        x0 = xs[s]
        x1 = xs[e - 1] if e - 1 < len(xs) else xs[-1]
        ax.axvspan(x0, x1, color=color, alpha=alpha, lw=0)

# ======================
# 통계/시각화
# ======================

def jaccard_from_sets(sample2pos):
    samples = list(sample2pos.keys())
    n = len(samples)
    jac = pd.DataFrame(np.eye(n), index=samples, columns=samples, dtype=float)
    for a, b in combinations(samples, 2):
        sa, sb = sample2pos[a], sample2pos[b]
        u = len(sa | sb)
        inter = len(sa & sb)
        v = 0.0 if u == 0 else inter / u
        jac.loc[a, b] = v
        jac.loc[b, a] = v
    return jac

def plot_heatmap(mat, title, out_png, dpi=200, vmin=None, vmax=None, cmap="viridis"):
    plt.figure(figsize=(max(6, 0.45*len(mat.columns)), max(4, 0.45*len(mat.index))))
    im = plt.imshow(mat.values, aspect='auto', interpolation='nearest',
                    vmin=vmin, vmax=vmax, cmap=cmap)
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.xticks(range(len(mat.columns)), mat.columns, rotation=90, fontsize=8)
    plt.yticks(range(len(mat.index)), mat.index, fontsize=8)
    plt.title(title, fontsize=12)
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi, bbox_inches='tight')
    plt.close()

def plot_scatter_with_missing(zdf, title, out_png, dpi=200):
    xs = zdf.index.values
    num_samples = zdf.shape[1]
    fig, axes = plt.subplots(num_samples, 1, figsize=(10, 2*max(1, num_samples)), sharex=True)
    if num_samples == 1:
        axes = [axes]
    for i, col in enumerate(zdf.columns):
        ax = axes[i]
        y = zdf[col].values
        nan_mask = np.isnan(y)
        shade_missing_spans(ax, xs, nan_mask, color="#d3d3d3", alpha=0.45)
        valid = ~nan_mask
        ax.scatter(xs[valid], y[valid], s=8, alpha=0.7, edgecolors='none')
        ax.axhline(0, color='gray', lw=0.8, linestyle='--')
        ax.set_ylabel(col, rotation=0, labelpad=40, fontsize=8)
        ax.grid(True, linestyle=':', alpha=0.35)
        if len(xs) > 0:
            ax.set_xlim(xs[0], xs[-1])
    axes[-1].set_xlabel("Genomic position")
    fig.suptitle(title, fontsize=12)
    plt.tight_layout(rect=[0,0,1,0.96])
    plt.savefig(out_png, dpi=dpi, bbox_inches='tight')
    plt.close(fig)

def plot_group_means(case_df, ctrl_df, chrom, out_png, dpi=200, smooth_window=0):
    # 같은 그리드 가정(메인에서 reindex 후 전달)
    xs = case_df.index.values
    case_mean = case_df.mean(axis=1, skipna=True)
    ctrl_mean = ctrl_df.mean(axis=1, skipna=True)
    if smooth_window and smooth_window > 1:
        case_mean = case_mean.rolling(smooth_window, min_periods=1, center=True).mean()
        ctrl_mean = ctrl_mean.rolling(smooth_window, min_periods=1, center=True).mean()

    plt.figure(figsize=(12,4))
    plt.plot(xs, case_mean.values, label="CASE mean")
    plt.plot(xs, ctrl_mean.values, label="CONTROL mean")
    # 결측(두 그룹 모두 NaN) 음영
    both_nan = case_df.isna().all(axis=1) & ctrl_df.isna().all(axis=1)
    shade_missing_spans(plt.gca(), xs, both_nan.values, color="#cccccc", alpha=0.5)
    plt.xlabel("Genomic position")
    plt.ylabel("Normalized depth")
    plt.title(f"{chrom} | Group mean profiles")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi, bbox_inches='tight')
    plt.close()

def plot_group_diff(case_df, ctrl_df, chrom, out_png, dpi=200, smooth_window=0):
    xs = case_df.index.values
    case_mean = case_df.mean(axis=1, skipna=True)
    ctrl_mean = ctrl_df.mean(axis=1, skipna=True)
    diff = case_mean - ctrl_mean
    if smooth_window and smooth_window > 1:
        diff = diff.rolling(smooth_window, min_periods=1, center=True).mean()

    plt.figure(figsize=(12,3.2))
    plt.plot(xs, diff.values)
    # 결측 음영: 두 그룹 모두 결측
    both_nan = case_df.isna().all(axis=1) & ctrl_df.isna().all(axis=1)
    shade_missing_spans(plt.gca(), xs, both_nan.values, color="#cccccc", alpha=0.5)
    plt.axhline(0, color='gray', lw=0.8, linestyle='--')
    plt.xlabel("Genomic position")
    plt.ylabel("CASE - CONTROL")
    plt.title(f"{chrom} | Group mean difference")
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi, bbox_inches='tight')
    plt.close()

def run_pca(zdf, chrom, out_png, dpi=200):
    X = zdf.fillna(0).T.values
    samples = list(zdf.columns)
    if X.shape[0] < 2 or min(X.shape) < 2:
        return
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(X)
    var = pca.explained_variance_ratio_ * 100.0
    plt.figure(figsize=(6,5))
    plt.scatter(coords[:,0], coords[:,1], s=40, alpha=0.85)
    for i, s in enumerate(samples):
        plt.text(coords[i,0], coords[i,1], s, fontsize=8, ha='left', va='bottom')
    plt.xlabel(f"PC1 ({var[0]:.1f}%)")
    plt.ylabel(f"PC2 ({var[1]:.1f}%)")
    plt.title(f"{chrom} | PCA (normalized depth)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi, bbox_inches='tight')
    plt.close()

# ======================
# 메인
# ======================

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    chrom_list = resolve_chroms(args.chroms)

    for chrom in chrom_list:
        print(f"[{chrom}] processing...", file=sys.stderr)

        # 1) 각 그룹 매트릭스 생성
        case_mat, case_pos = build_matrix_from_glob(
            args.case_glob, chrom, args.min_depth, args.norm_scale
        )
        ctrl_mat, ctrl_pos = build_matrix_from_glob(
            args.control_glob, chrom, args.min_depth, args.norm_scale
        )

        if case_mat.empty and ctrl_mat.empty:
            print(f"[{chrom}] no files or all empty after filtering.", file=sys.stderr)
            continue

        # 2) 하나의 그리드로 정렬/결측 보존
        #   - 두 그룹의 인덱스 범위를 합친 뒤 동일 그리드로 맞춤
        all_idx = []
        if not case_mat.empty: all_idx.append(case_mat.index)
        if not ctrl_mat.empty: all_idx.append(ctrl_mat.index)
        combined_index = pd.Index(sorted(set(np.concatenate(all_idx)))) if all_idx else pd.Index([])

        # reindex to uniform bin grid
        if len(combined_index) > 0:
            start = int(combined_index.min() // args.bin_size * args.bin_size)
            end = int(np.ceil(combined_index.max() / args.bin_size) * args.bin_size)
            grid = np.arange(start, end + 1, args.bin_size, dtype=int)
        else:
            grid = np.array([], dtype=int)

        case_df = case_mat.reindex(grid) if not case_mat.empty else pd.DataFrame(index=grid)
        ctrl_df = ctrl_mat.reindex(grid) if not ctrl_mat.empty else pd.DataFrame(index=grid)

        # 3) 비어있는 전-NaN/0 행 제거(옵션)
        def drop_all_empty(df):
            if df.empty:
                return df
            mask = ~(df.isna().all(axis=1) | (df.fillna(0).sum(axis=1) == 0))
            return df.loc[mask]
        if not args.include_empty:
            case_df = drop_all_empty(case_df)
            ctrl_df = drop_all_empty(ctrl_df)

            # grid 재동기화
            grid = case_df.index.union(ctrl_df.index)
            case_df = case_df.reindex(grid)
            ctrl_df = ctrl_df.reindex(grid)

        # 4) 출력 폴더
        chr_out = os.path.join(args.outdir, chrom)
        os.makedirs(chr_out, exist_ok=True)

        # 5) 매트릭스 저장
        case_df.sort_index().to_csv(os.path.join(chr_out, f"{chrom}.CASE.normalized_matrix.csv"))
        ctrl_df.sort_index().to_csv(os.path.join(chr_out, f"{chrom}.CONTROL.normalized_matrix.csv"))

        # 6) 샘플 스캐터(결측 음영 포함)
        if not case_df.empty and case_df.shape[1] > 0:
            plot_scatter_with_missing(
                case_df.sort_index(),
                title=f"{chrom} | CASE Normalized_depth (bin={args.bin_size}bp)",
                out_png=os.path.join(chr_out, f"{chrom}.CASE.scatter.png"),
                dpi=args.png_dpi
            )
        if not ctrl_df.empty and ctrl_df.shape[1] > 0:
            plot_scatter_with_missing(
                ctrl_df.sort_index(),
                title=f"{chrom} | CONTROL Normalized_depth (bin={args.bin_size}bp)",
                out_png=os.path.join(chr_out, f"{chrom}.CONTROL.scatter.png"),
                dpi=args.png_dpi
            )

        # 7) 상관/자카드
        if case_df.shape[1] >= 2:
            corr = case_df.fillna(0).corr(method="pearson")
            corr.to_csv(os.path.join(chr_out, f"{chrom}.CASE.pearson_corr.csv"))
            plot_heatmap(corr, f"{chrom} | CASE Pearson corr",
                         os.path.join(chr_out, f"{chrom}.CASE.pearson_corr.png"),
                         dpi=args.png_dpi, vmin=-1, vmax=1, cmap="coolwarm")

        if ctrl_df.shape[1] >= 2:
            corr = ctrl_df.fillna(0).corr(method="pearson")
            corr.to_csv(os.path.join(chr_out, f"{chrom}.CONTROL.pearson_corr.csv"))
            plot_heatmap(corr, f"{chrom} | CONTROL Pearson corr",
                         os.path.join(chr_out, f"{chrom}.CONTROL.pearson_corr.png"),
                         dpi=args.png_dpi, vmin=-1, vmax=1, cmap="coolwarm")

        # 자카드(존재/부재)
        if len(case_pos) >= 2:
            jac = jaccard_from_sets(case_pos)
            jac.to_csv(os.path.join(chr_out, f"{chrom}.CASE.jaccard.csv"))
            plot_heatmap(jac, f"{chrom} | CASE Jaccard (depth≥{args.min_depth})",
                         os.path.join(chr_out, f"{chrom}.CASE.jaccard.png"),
                         dpi=args.png_dpi, vmin=0, vmax=1, cmap="viridis")
        if len(ctrl_pos) >= 2:
            jac = jaccard_from_sets(ctrl_pos)
            jac.to_csv(os.path.join(chr_out, f"{chrom}.CONTROL.jaccard.csv"))
            plot_heatmap(jac, f"{chrom} | CONTROL Jaccard (depth≥{args.min_depth})",
                         os.path.join(chr_out, f"{chrom}.CONTROL.jaccard.png"),
                         dpi=args.png_dpi, vmin=0, vmax=1, cmap="viridis")

        # 8) PCA
        if not args.no_pca:
            if case_df.shape[1] >= 2:
                run_pca(case_df.sort_index(), chrom,
                        os.path.join(chr_out, f"{chrom}.CASE.pca.png"), dpi=args.png_dpi)
            if ctrl_df.shape[1] >= 2:
                run_pca(ctrl_df.sort_index(), chrom,
                        os.path.join(chr_out, f"{chrom}.CONTROL.pca.png"), dpi=args.png_dpi)

        # 9) 그룹 평균/차이
        if not case_df.empty and not ctrl_df.empty:
            plot_group_means(case_df.sort_index(), ctrl_df.sort_index(), chrom,
                             os.path.join(chr_out, f"{chrom}.GroupMeans.png"),
                             dpi=args.png_dpi, smooth_window=args.smooth_window)
            plot_group_diff(case_df.sort_index(), ctrl_df.sort_index(), chrom,
                            os.path.join(chr_out, f"{chrom}.GroupDiff.png"),
                            dpi=args.png_dpi, smooth_window=args.smooth_window)

        print(f"[{chrom}] Done -> {chr_out}", file=sys.stderr)

if __name__ == "__main__":
    main()
