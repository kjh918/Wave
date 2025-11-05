#!/usr/bin/env python3
import argparse
import gzip
import sys

def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize per-sample genotype counts and called variant counts from a multi-sample VCF(.gz)."
    )
    p.add_argument("vcf", help="Input VCF(.vcf.gz) file")
    p.add_argument("-o", "--out", help="Output TSV file (default: stdout)")
    p.add_argument("--include-filtered", action="store_true",
                   help="Include FILTER!=PASS variants as well (default: included). (Flag kept for clarity)")
    return p.parse_args()

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def classify_gt(gt):
    """
    Normalize and classify GT.
    Returns one of: 'missing', 'hom_ref', 'het', 'hom_alt'
    Also returns booleans: called (not missing), non_ref (has any ALT allele)
    """
    if gt is None or gt == "" or gt == ".":
        return "missing", False, False

    # GT could have extra fields like "0/1:35:..." – caller passes only the GT token
    # Normalize phased/unphased
    sep = "/" if "/" in gt else ("|" if "|" in gt else None)
    if sep is None:
        # single allele or malformed; treat '.' as missing
        if "." in gt:
            return "missing", False, False
        alleles = [gt]
    else:
        alleles = gt.split(sep)

    # Missing if any allele is '.'
    if any(a == "." or a == "" for a in alleles):
        return "missing", False, False

    # Called if we reached here
    called = True

    # Determine non_ref: any allele not equal to '0'
    non_ref = any(a != "0" for a in alleles)

    # Hom/Het classification
    uniq = set(alleles)
    if uniq == {"0"}:
        return "hom_ref", called, False
    if len(uniq) == 1 and "0" not in uniq:
        # e.g., {'1'} or {'2'}
        return "hom_alt", called, True
    # Mixed (e.g., 0/1, 1/2)
    return "het", called, True

def main():
    args = parse_args()

    with open_maybe_gzip(args.vcf) as fh:
        samples = []
        gt_index = None

        # stats per sample
        stats = {}
        total_variants = 0

        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                # fixed VCF columns: 0..8, then samples start at 9
                samples = header[9:]
                # initialize stats
                for s in samples:
                    stats[s] = {
                        "hom_ref": 0,
                        "het": 0,
                        "hom_alt": 0,
                        "missing": 0,
                        "called_variants": 0,
                        "non_ref_variants": 0,
                    }
                continue

            # data line
            total_variants += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                # no samples
                continue

            fmt = fields[8].split(":")
            try:
                gt_index = fmt.index("GT")
            except ValueError:
                # no GT field; count all samples as missing for this site
                for s in samples:
                    stats[s]["missing"] += 1
                continue

            # Iterate samples
            for i, s in enumerate(samples, start=9):
                sample_field = fields[i]
                tokens = sample_field.split(":")
                gt = tokens[gt_index] if gt_index < len(tokens) else "."
                cls, called, non_ref = classify_gt(gt)

                stats[s][cls] += 1
                if called:
                    stats[s]["called_variants"] += 1
                if non_ref:
                    stats[s]["non_ref_variants"] += 1

        # Output
        outfh = open(args.out, "w") if args.out else sys.stdout
        with outfh:
            header_cols = [
                "sample",
                "total_variants",
                "called_variants",
                "non_ref_variants",
                "hom_ref",
                "het",
                "hom_alt",
                "missing",
            ]
            print("\t".join(header_cols), file=outfh)
            for s in samples:
                row = [
                    s,
                    str(total_variants),
                    str(stats[s]["called_variants"]),
                    str(stats[s]["non_ref_variants"]),
                    str(stats[s]["hom_ref"]),
                    str(stats[s]["het"]),
                    str(stats[s]["hom_alt"]),
                    str(stats[s]["missing"]),
                ]
                print("\t".join(row), file=outfh)

if __name__ == "__main__":
    main()
