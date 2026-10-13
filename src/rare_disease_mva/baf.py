"""Per-chromosome allele balance at heterozygous PASS SNVs.

Mosaic aneuploidy of a chromosome in a fraction f of cells shifts het-SNP B-allele
fractions away from 0.5 (trisomy: 1/2 -> (1+f)/(2+f) or 1/(2+f)), which widens the
|BAF - 0.5| distribution relative to the other autosomes. Also reports mean DP.
"""
import sys
from collections import defaultdict

import numpy as np
from cyvcf2 import VCF


def main(vcf_path, min_dp=25):
    baf, dp = defaultdict(list), defaultdict(list)
    for v in VCF(vcf_path):
        if v.FILTER is not None or not v.is_snp:
            continue
        if v.gt_types[0] != 1:  # het
            continue
        ad = v.format("AD")[0]
        tot = int(ad[0]) + int(ad[1])
        if tot < min_dp:
            continue
        baf[v.CHROM].append(ad[1] / tot)
        dp[v.CHROM].append(tot)
    autos = [f"chr{i}" for i in range(1, 23)]
    dev_all = np.concatenate([np.abs(np.array(baf[c]) - 0.5) for c in autos])
    print("chrom\tn_het\tmean_dp\tmedian_baf\tmean_absdev\tz_vs_autosomes")
    mu = np.mean([np.mean(np.abs(np.array(baf[c]) - 0.5)) for c in autos])
    sd = np.std([np.mean(np.abs(np.array(baf[c]) - 0.5)) for c in autos])
    for c in autos + ["chrX"]:
        b = np.array(baf[c])
        if len(b) == 0:
            continue
        m = np.mean(np.abs(b - 0.5))
        print(f"{c}\t{len(b)}\t{np.mean(dp[c]):.1f}\t{np.median(b):.3f}\t{m:.4f}\t{(m - mu) / sd:+.2f}")
    print(f"# autosome pooled mean |BAF-0.5| = {dev_all.mean():.4f}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1])
