#!/usr/bin/env python3
import sys, argparse, pysam, re, csv, math
from collections import defaultdict
from typing import Optional, Tuple, List

# ---------- helpers ----------
import matplotlib
matplotlib.use("Agg")  # 서버 환경에서 GUI 없이 저장용
import matplotlib.pyplot as plt

BASE_RE = re.compile(r"^[ACGT]$", re.I)

def is_base(x: Optional[str]) -> bool:
    return x is not None and BASE_RE.match(x) and len(x) == 1

def is_biallelic_snp(rec) -> bool:
    return is_base(rec.ref) and rec.alts and len(rec.alts) == 1 and is_base(rec.alts[0])

def norm_gt(gt: Optional[Tuple[int, ...]]) -> Optional[Tuple[int, int]]:
    if not gt or any(a is None for a in gt) or len(gt) != 2:
        return None
    a, b = sorted(gt)
    return (a, b)

def zygosity(gt01: Optional[Tuple[int,int]]) -> str:
    if gt01 is None: return "NO_CALL"
    a,b = gt01
    if a==0 and b==0: return "HOM_REF"
    if a!=b:          return "HET"
    if a==b and a>0:  return "HOM_ALT"
    return "NO_CALL"

def callable_sample(sample, min_dp: int, min_gq: int) -> bool:
    dp = sample.get("DP")
    gq = sample.get("GQ")
    gt = norm_gt(sample.get("GT"))
    if gt is None: return False
    if min_dp and (dp is None or dp < min_dp): return False
    if min_gq and (gq is None or gq < min_gq): return False
    return True

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(
        description="ADO / LDO / Amplification per single-cell vs control + genotype stats + plots"
    )
    ap.add_argument("vcf", help="joint-called all-sites VCF (.vcf.gz) with index")
    ap.add_argument("--control", required=True, help="control sample name (cell line gDNA)")
    ap.add_argument("--min-dp-control", type=int, default=15)
    ap.add_argument("--min-dp-sample", type=int, default=10)
    ap.add_argument("--min-gq", type=int, default=0)
    ap.add_argument("--only-biallelic-snp", action="store_true")
    ap.add_argument("--out-metrics", default="metrics_per_sample.tsv")
    ap.add_argument("--out-details", default="details.tsv")
    ap.add_argument("--out-genostats", default="genotype_stats_per_sample.tsv")
    ap.add_argument("--prefix-plots", default="dropout_plots")
    ap.add_argument("--max-records", type=int, default=0)
    args = ap.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    samples = list(vcf.header.samples)
    if args.control not in samples:
        sys.exit(f"[ERROR] Control '{args.control}' not in VCF: {samples}")
    sc_samples = [s for s in samples if s != args.control]

    # Metrics per single-cell
    M = {s: defaultdict(int) for s in sc_samples}   # bulk_het_den, bulk_hom_den, ado_num, ldo_num, amp_num, match_het, match_hom, other_mismatch, nocall_on_hom
    # Genotype distribution per sample (including control)
    G = {s: defaultdict(int) for s in samples}      # HOM_REF/HET/HOM_ALT/NO_CALL counts
    DPsum = defaultdict(int)

    details: List[dict] = []
    processed = 0

    for rec in vcf.fetch():
        if args.only_biallelic_snp and not is_biallelic_snp(rec):
            continue

        ctrl = rec.samples[args.control]
        if not callable_sample(ctrl, args.min_dp_control, args.min_gq):
            continue

        ctrl_gt = norm_gt(ctrl.get("GT"))
        ctrl_zy = zygosity(ctrl_gt)

        # genotype distribution for all samples at this locus
        for s in samples:
            sm = rec.samples[s]
            gt = norm_gt(sm.get("GT"))
            zy = zygosity(gt)
            G[s][zy] += 1
            dp = sm.get("DP")
            if dp:
                DPsum[s] += dp

        # metrics for sc samples
        for s in sc_samples:
            sm = rec.samples[s]

            if ctrl_zy == "HET":
                M[s]["bulk_het_den"] += 1

                if not callable_sample(sm, args.min_dp_sample, args.min_gq):
                    M[s]["ldo_num"] += 1
                    details.append({
                        "CHROM": rec.chrom, "POS": rec.pos, "REF": rec.ref,
                        "ALT": rec.alts[0] if rec.alts else ".",
                        "class": "LDO",
                        "control_GT": "/".join(map(str, ctrl.get("GT"))) if ctrl.get("GT") else "./.",
                        "sample": s, "sample_GT": "./.",
                        "note": "single-cell not callable"
                    })
                    continue

                sm_gt = norm_gt(sm.get("GT"))
                sm_zy = zygosity(sm_gt)

                if sm_zy in ("HOM_REF","HOM_ALT"):
                    M[s]["ado_num"] += 1
                    details.append({
                        "CHROM": rec.chrom, "POS": rec.pos, "REF": rec.ref,
                        "ALT": rec.alts[0] if rec.alts else ".",
                        "class": "ADO",
                        "control_GT": "/".join(map(str, ctrl.get("GT"))),
                        "sample": s, "sample_GT": "/".join(map(str, sm.get("GT")))
                    })
                elif sm_zy == "HET":
                    M[s]["match_het"] += 1
                else:
                    M[s]["other_mismatch"] += 1

            elif ctrl_zy in ("HOM_REF","HOM_ALT"):
                M[s]["bulk_hom_den"] += 1

                if not callable_sample(sm, args.min_dp_sample, args.min_gq):
                    M[s]["nocall_on_hom"] += 1
                    continue

                sm_gt = norm_gt(sm.get("GT"))
                sm_zy = zygosity(sm_gt)
                if sm_zy == "HET":
                    M[s]["amp_num"] += 1
                    details.append({
                        "CHROM": rec.chrom, "POS": rec.pos, "REF": rec.ref,
                        "ALT": rec.alts[0] if rec.alts else ".",
                        "class": "Amplification_Error",
                        "control_GT": "/".join(map(str, ctrl.get("GT"))),
                        "sample": s, "sample_GT": "/".join(map(str, sm.get("GT")))
                    })
                elif sm_zy in ("HOM_REF","HOM_ALT"):
                    M[s]["match_hom"] += 1
                else:
                    M[s]["other_mismatch"] += 1

        processed += 1
        if args.max_records and processed >= args.max_records:
            break

    # ----- write TSVs -----
    with open(args.out_metrics, "w") as w:
        w.write("Sample\tBulkHET_den\tADO_num\tLDO_num\tADO_rate\tLDO_rate\t"
                "BulkHOM_den\tAmpl_num\tAmpl_rate\tMatch_HET\tMatch_HOM\tNoCall_on_HOM\tOther_Mismatch\n")
        for s in sc_samples:
            bhet = M[s]["bulk_het_den"]
            bhom = M[s]["bulk_hom_den"]
            ado  = M[s]["ado_num"]
            ldo  = M[s]["ldo_num"]
            amp  = M[s]["amp_num"]
            r_ado = (ado / bhet) if bhet > 0 else 0.0
            r_ldo = (ldo / bhet) if bhet > 0 else 0.0
            r_amp = (amp / bhom) if bhom > 0 else 0.0
            w.write(
                f"{s}\t{bhet}\t{ado}\t{ldo}\t{r_ado:.6f}\t{r_ldo:.6f}\t"
                f"{bhom}\t{amp}\t{r_amp:.6f}\t"
                f"{M[s]['match_het']}\t{M[s]['match_hom']}\t{M[s]['nocall_on_hom']}\t{M[s]['other_mismatch']}\n"
            )

    with open(args.out_genostats, "w") as w:
        w.write("Sample\tHOM_REF\tHET\tHOM_ALT\tNO_CALL\tCallable\tMean_DP\n")
        for s in samples:
            hr = G[s]["HOM_REF"]
            ht = G[s]["HET"]
            ha = G[s]["HOM_ALT"]
            nc = G[s]["NO_CALL"]
            callable_sites = hr + ht + ha
            mean_dp = (DPsum[s]/callable_sites) if callable_sites>0 else 0.0
            w.write(f"{s}\t{hr}\t{ht}\t{ha}\t{nc}\t{callable_sites}\t{mean_dp:.2f}\n")

    with open(args.out_details, "w", newline="") as w:
        cols = ["CHROM","POS","REF","ALT","class","control_GT","sample","sample_GT","note"]
        dw = csv.DictWriter(w, fieldnames=cols, delimiter="\t")
        dw.writeheader()
        for r in details:
            dw.writerow(r)

    # ----- plots -----
    # load back small tables for plotting without pandas, keep it simple
    def read_metrics(path):
        rows = []
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                parts = line.rstrip("\n").split("\t")
                rows.append(dict(zip(header, parts)))
        return rows

    def read_genostats(path):
        rows = []
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                parts = line.rstrip("\n").split("\t")
                rows.append(dict(zip(header, parts)))
        return rows

    met = read_metrics(args.out_metrics)
    gst = read_genostats(args.out_genostats)

    # Figure 1: ADO / LDO / Ampl rates per single-cell
    samples_x = [r["Sample"] for r in met]
    ado_rates = [float(r["ADO_rate"]) for r in met]
    ldo_rates = [float(r["LDO_rate"]) for r in met]
    amp_rates = [float(r["Ampl_rate"]) for r in met]

    x = list(range(len(samples_x)))
    width = 0.25
    fig1 = plt.figure(figsize=(max(8, len(x)*0.6), 5))
    plt.bar([i - width for i in x], ado_rates, width, label="ADO rate")
    plt.bar(x, ldo_rates, width, label="LDO rate")
    plt.bar([i + width for i in x], amp_rates, width, label="Amplification rate")
    plt.xticks(x, samples_x, rotation=45, ha="right")
    plt.ylabel("Rate")
    plt.title("ADO / LDO / Amplification rates per single-cell")
    plt.legend()
    plt.tight_layout()
    fig1.savefig(f"{args.prefix_plots}.rates.png", dpi=200)

    # Figure 2: Bulk denominators per single-cell (bars)
    bhet = [int(r["BulkHET_den"]) for r in met]
    bhom = [int(r["BulkHOM_den"]) for r in met]
    fig2 = plt.figure(figsize=(max(8, len(x)*0.6), 4.5))
    plt.bar([i - width/2 for i in x], bhet, width, label="Bulk HET denom")
    plt.bar([i + width/2 for i in x], bhom, width, label="Bulk HOM denom")
    plt.xticks(x, samples_x, rotation=45, ha="right")
    plt.ylabel("# sites")
    plt.title("Denominators per single-cell (bulk-based)")
    plt.legend()
    plt.tight_layout()
    fig2.savefig(f"{args.prefix_plots}.denoms.png", dpi=200)

    # Figure 3: Genotype distribution per sample (stacked bars)
    gs_samples = [r["Sample"] for r in gst]
    hr = [int(r["HOM_REF"]) for r in gst]
    ht = [int(r["HET"]) for r in gst]
    ha = [int(r["HOM_ALT"]) for r in gst]
    nc = [int(r["NO_CALL"]) for r in gst]

    xpos = list(range(len(gs_samples)))
    fig3 = plt.figure(figsize=(max(10, len(xpos)*0.6), 6))
    b1 = plt.bar(xpos, hr, label="HOM_REF")
    b2 = plt.bar(xpos, ht, bottom=hr, label="HET")
    b3 = plt.bar(xpos, ha, bottom=[hr[i]+ht[i] for i in range(len(hr))], label="HOM_ALT")
    b4 = plt.bar(xpos, nc, bottom=[hr[i]+ht[i]+ha[i] for i in range(len(hr))], label="NO_CALL")
    plt.xticks(xpos, gs_samples, rotation=45, ha="right")
    plt.ylabel("# sites")
    plt.title("Genotype distribution per sample")
    plt.legend()
    plt.tight_layout()
    fig3.savefig(f"{args.prefix_plots}.genotypes.png", dpi=200)

    # Figure 4: Mean DP per sample (bar)
    mean_dp = [float(r["Mean_DP"]) for r in gst]
    fig4 = plt.figure(figsize=(max(10, len(xpos)*0.6), 4.5))
    plt.bar(xpos, mean_dp)
    plt.xticks(xpos, gs_samples, rotation=45, ha="right")
    plt.ylabel("Mean DP over callable sites")
    plt.title("Mean depth per sample")
    plt.tight_layout()
    fig4.savefig(f"{args.prefix_plots}.meanDP.png", dpi=200)

    print(f"[PLOTS] Saved: "
          f"{args.prefix_plots}.rates.png, "
          f"{args.prefix_plots}.denoms.png, "
          f"{args.prefix_plots}.genotypes.png, "
          f"{args.prefix_plots}.meanDP.png")

if __name__ == "__main__":
    main()
