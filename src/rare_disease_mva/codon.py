"""Minimal CDS-level consequence for SNVs on a + or - strand transcript (sanity check vs VEP)."""
import gzip
import re
import sys

import pysam
from Bio.Seq import Seq


def cds_exons(gtf_gz, transcript_id):
    out, strand = [], None
    with gzip.open(gtf_gz, "rt") as fh:
        for line in fh:
            if transcript_id not in line:
                continue
            f = line.split("\t")
            if f[2] == "CDS" and re.search(f'transcript_id "{transcript_id}', f[8]):
                out.append((f[0], int(f[3]), int(f[4])))
                strand = f[6]
    return sorted(out, key=lambda x: x[1]), strand


def snv_consequence(fasta, exons, strand, pos, alt):
    chrom = exons[0][0]
    seq = "".join(fasta.fetch(chrom, s - 1, e) for _, s, e in exons)
    coords = [p for _, s, e in exons for p in range(s, e + 1)]
    if strand == "-":
        seq, coords = str(Seq(seq).reverse_complement()), coords[::-1]
        alt = str(Seq(alt).reverse_complement())
    i = coords.index(pos)
    mut = seq[:i] + alt + seq[i + 1:]
    aa_i = i // 3
    ref_aa = str(Seq(seq[aa_i * 3:aa_i * 3 + 3]).translate())
    alt_aa = str(Seq(mut[aa_i * 3:aa_i * 3 + 3]).translate())
    return f"c.{i + 1}{seq[i]}>{alt} p.{ref_aa}{aa_i + 1}{alt_aa} (CDS len {len(seq)}, {len(seq)//3 - 1} aa)"


if __name__ == "__main__":
    fa_path, gtf, tid, *variants = sys.argv[1:]
    exons, strand = cds_exons(gtf, tid)
    fa = pysam.FastaFile(fa_path)
    for v in variants:
        pos, alt = v.split(":")
        print(v, snv_consequence(fa, exons, strand, int(pos), alt))
