#!/usr/bin/env python3
import sys
import re
import csv
import argparse
from typing import Optional, List, Tuple, Dict
import pysam

BASE_RE = re.compile(r"^[ACGT]$", re.I)

def is_base(x: Optional[str]) -> bool:
    return x is not None and BASE_RE.match(x or "") and len(x) == 1

def pick_single_snv_alt_index(rec: pysam.libcbcf.VariantRecord) -> Optional[int]:
    """
    심볼릭 ALT(<NON_REF> 등), INDEL, 다중 ALT 제거하고
    '실제 염기 ALT'가 정확히 1개인 경우 그 ALT의 1-based index 반환.
    (AD에서 ref=0, alt=해당 인덱스)
    """
    if not is_base(rec.ref):
        return None
    if rec.alts is None:
        return None
    real = [i for i, a in enumerate(rec.alts, start=1) if is_base(a)]
    return real[0] if len(real) == 1 else None

def norm_gt(gt: Optional[Tuple[int, ...]], allowed_alt_idx: int) -> Optional[Tuple[int, int]]:
    """GT를 (0,0)/(0,1)/(1,1)로 정규화. 허용 ALT 외 인덱스가 끼면 None."""
    if not gt or any(a is None for a in gt) or len(gt) != 2:
        return None
    a, b = gt
    allowed = {0, allowed_alt_idx}
    if a not in allowed or b not in allowed:
        return None
    a = 0 if a == 0 else 1
    b = 0 if b == 0 else 1
    return tuple(sorted((a, b)))

def zygosity01(gt01: Optional[Tuple[int,int]]) -> str:
    if gt01 is None: return "NA"
    a,b = gt01
    if a==0 and b==0: return "HOM_REF"
    if a!=b:          return "HET"
    if a==b and a==1: return "HOM_ALT"
    return "NA"

def callable_sample(sm, min_dp: int, min_gq: int, allowed_alt_idx: int) -> bool:
    gt01 = norm_gt(sm.get("GT"), allowed_alt_idx)
    if gt01 is None: return False
    dp = sm.get("DP")
    if min_dp and (dp is None or dp < min_dp): return False
    gq = sm.get("GQ")
    if min_gq and (gq is None or gq < min_gq): return False
    return True

def get_ad_pair(sm, alt_idx_1based: int) -> Tuple[Optional[int], Optional[int]]:
    """
    FORMAT/AD에서 ref 카운트(AD[0])와 선택 ALT 카운트(AD[alt_idx])만 반환.
    없거나 파싱 불가하면 (None, None).
    """
    ad = sm.get("AD")
    if ad is None: return (None, None)
    try:
        refc = int(ad[0]) if len(ad) > 0 else None
        altc = int(ad[alt_idx_1based]) if len(ad) > alt_idx_1based else None
        return (refc, altc)
    except Exception:
        return (None, None)

def frac(refc: Optional[int], altc: Optional[int]) -> Optional[float]:
    if refc is None or altc is None: return None
    tot = refc + altc
    if tot <= 0: return None
    return altc / tot

def main():
    ap = argparse.ArgumentParser(
        description="Make a table of per-sample allele counts & fractions at sites where CONTROL is heterozygous."
    )
    ap.add_argument("vcf", help="joint-called all-sites VCF (.vcf.gz) with index")
    ap.add_argument("--control", required=True, help="control sample name")
    ap.add_argument("--samples", nargs="*", help="subset of samples to include (default: all including control)")
    ap.add_argument("--min-dp-control", type=int, default=10)
    ap.add_argument("--min-gq-control", type=int, default=0)
    ap.add_argument("--min-dp-sample", type=int, default=0)
    ap.add_argument("--min-gq-sample", type=int, default=0)
    ap.add_argument("--only-biallelic-snp", action="store_true", help="require exactly one ALT and both REF/ALT are single bases")
    ap.add_argument("--out", default="controlHET_alleleFractions.tsv")
    ap.add_argument("--format", choices=["wide", "long"], default="wide",
                    help="wide: 한 행에 모든 샘플 컬럼 / long: 샘플별 행 분리")
    ap.add_argument("--max-records", type=int, default=0)
    args = ap.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    all_samples = list(vcf.header.samples)
    if args.control not in all_samples:
        sys.exit(f"[ERROR] Control '{args.control}' not found. Available: {all_samples}")

    if args.samples:
        # 보존: control은 항상 포함, 지정 목록과 교집합
        include = [s for s in args.samples if s in all_samples]
        if args.control not in include:
            include = [args.control] + include
        samples = include
    else:
        samples = all_samples

    # wide 포맷일 때 헤더 구성
    if args.format == "wide":
        header = ["CHROM", "POS", "REF", "ALT"]
        for s in samples:
            header += [f"{s}:ref", f"{s}:alt", f"{s}:dp", f"{s}:alt_frac"]
    else:  # long
        header = ["CHROM", "POS", "REF", "ALT", "SAMPLE", "REF_COUNT", "ALT_COUNT", "DP", "ALT_FRAC"]

    outf = open(args.out, "w", newline="")
    writer = csv.writer(outf, delimiter="\t")
    writer.writerow(header)

    processed = 0
    kept = 0

    for rec in vcf.fetch():
        # 1) SNV ALT 1개 고르기
        alt_idx = pick_single_snv_alt_index(rec)
        if alt_idx is None:
            continue
        if args.only_biallelic_snp and (rec.alts is None or len(rec.alts) != 1):
            continue

        # 2) control이 callable + HET인 자리만 유지
        ctrl = rec.samples[args.control]
        if not callable_sample(ctrl, args.min_dp_control, args.min_gq_control, alt_idx):
            continue
        ctrl_gt01 = norm_gt(ctrl.get("GT"), alt_idx)
        if ctrl_gt01 is None or zygosity01(ctrl_gt01) != "HET":
            continue

        chrom, pos, ref, alt = rec.chrom, rec.pos, rec.ref, rec.alts[alt_idx-1]
        kept += 1

        if args.format == "wide":
            row: List[Optional[str]] = [chrom, pos, ref, alt]
            for s in samples:
                sm = rec.samples[s]
                # 샘플 callable 체크(컨트롤과 달리: 표에는 남기되, 미달이면 빈칸 처리)
                dp = sm.get("DP")
                gt01 = norm_gt(sm.get("GT"), alt_idx)
                callable_ok = callable_sample(sm, args.min_dp_sample, args.min_gq_sample, alt_idx)

                refc, altc = get_ad_pair(sm, alt_idx)
                af = frac(refc, altc)

                row += [
                    "" if refc is None or not callable_ok else str(refc),
                    "" if altc is None or not callable_ok else str(altc),
                    "" if dp is None   or not callable_ok else str(dp),
                    "" if af is None   or not callable_ok else f"{af:.4f}",
                ]
            writer.writerow(row)

        else:  # long
            for s in samples:
                sm = rec.samples[s]
                dp = sm.get("DP")
                callable_ok = callable_sample(sm, args.min_dp_sample, args.min_gq_sample, alt_idx)
                refc, altc = get_ad_pair(sm, alt_idx)
                af = frac(refc, altc)
                writer.writerow([
                    chrom, pos, ref, alt, s,
                    "" if refc is None or not callable_ok else refc,
                    "" if altc is None or not callable_ok else altc,
                    "" if dp   is None or not callable_ok else dp,
                    "" if af   is None or not callable_ok else f"{af:.4f}",
                ])

        processed += 1
        if args.max_records and processed >= args.max_records:
            break

    outf.close()
    print(f"[DONE] Wrote: {args.out}")
    print(f"[INFO] Control: {args.control}")
    print(f"[INFO] Samples: {', '.join(samples)}")
    print(f"[INFO] Sites kept (control HET): {kept}")

if __name__ == "__main__":
    main()
