"""Stage 1: conservation of BUBR1 Asn1002 (and known MVA1 residues) across orthologs.

Pulls BUB1B orthologues with protein sequences from Ensembl Compara, aligns each
to the human sequence (MAFFT, all together), and reports the residue observed at
each position of interest per species, plus the BUB1 paralog. Writes
reports/track2/tables/s1_conservation.tsv and s1_conservation_summary.tsv.

Usage: python -I scripts/track2/s1_conservation.py
"""

import collections
import json
import subprocess
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data/work/track2/cache"
OUT = ROOT / "reports/track2/tables"
MAFFT = Path.home() / "micromamba/envs/bio/bin/mafft"

BUB1B = "ENSG00000156970"
BUB1 = "ENSG00000169679"
POSITIONS = {737: "L737 (nonsense)", 814: "R814 (MVA1)", 844: "L844 (MVA1)",
             909: "I909 (MVA1)", 921: "Q921 (MVA1)", 1002: "N1002 (patient)",
             1008: "T1008 (PLK1 site)", 1012: "L1012 (MVA1)"}


def get_json(url, dest):
    if not dest.exists():
        req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
        dest.write_bytes(urllib.request.urlopen(req, timeout=120).read())
    return json.loads(dest.read_text())


def human_seq(gene):
    d = get_json(f"https://rest.ensembl.org/sequence/id/{gene}?type=protein;multiple_sequences=1;content-type=application/json",
                 CACHE / f"seq_{gene}.json")
    return max(d, key=lambda r: len(r["seq"]))


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    homol = get_json(
        f"https://rest.ensembl.org/homology/id/human/{BUB1B}?type=orthologues;sequence=protein;"
        "aligned=0;content-type=application/json", CACHE / "compara_BUB1B_orthologues.json")
    seqs = {}
    uniprot = json.loads((CACHE / "uniprot_O60566.json").read_text())["sequence"]["value"]
    seqs["homo_sapiens|BUB1B"] = uniprot
    for h in homol["data"][0]["homologies"]:
        t = h["target"]
        if h["type"] != "ortholog_one2one":
            continue
        seqs[f"{t['species']}|{t['protein_id']}"] = t["seq"].replace("-", "")
    para = human_seq(BUB1)
    seqs["homo_sapiens|BUB1_paralog"] = para["seq"]

    fasta = CACHE / "bub1b_orthologs.fa"
    fasta.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    aln_txt = subprocess.run([str(MAFFT), "--auto", "--quiet", str(fasta)],
                             capture_output=True, text=True, check=True).stdout
    aln, name = {}, None
    for line in aln_txt.splitlines():
        if line.startswith(">"):
            name = line[1:]
            aln[name] = ""
        else:
            aln[name] += line.strip()

    ref = aln["homo_sapiens|BUB1B"]
    col_of, n = {}, 0
    for i, c in enumerate(ref):
        if c != "-":
            n += 1
            col_of[n] = i
    assert uniprot[1001] == "N" and ref[col_of[1002]] == "N"

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "s1_conservation.tsv", "w") as fh:
        fh.write("sequence\t" + "\t".join(POSITIONS.values()) + "\n")
        for k, s in aln.items():
            fh.write(k + "\t" + "\t".join(s[col_of[p]] for p in POSITIONS) + "\n")

    orth = [k for k in aln if not k.endswith("BUB1_paralog")]
    with open(OUT / "s1_conservation_summary.tsv", "w") as fh:
        fh.write("position\tn_orthologs\tidentity_to_human\tresidues_seen\tlysine_seen\tBUB1_paralog\n")
        for p, label in POSITIONS.items():
            col = [aln[k][col_of[p]] for k in orth]
            obs = [c for c in col if c != "-"]
            h = uniprot[p - 1]
            cnt = collections.Counter(obs)
            fh.write(f"{label}\t{len(obs)}\t{cnt[h] / len(obs):.2f}\t"
                     f"{','.join(f'{a}:{b}' for a, b in cnt.most_common())}\t"
                     f"{'yes' if cnt['K'] else 'no'}\t{aln['homo_sapiens|BUB1_paralog'][col_of[p]]}\n")
    print((OUT / "s1_conservation_summary.tsv").read_text())
    # Local window for the report
    lo, hi = col_of[990], col_of[1020]
    for k in ["homo_sapiens|BUB1B"] + [k for k in aln if any(s in k for s in
              ("mus_musculus", "gallus", "xenopus", "danio", "callorhinchus", "petromyzon"))] + ["homo_sapiens|BUB1_paralog"]:
        print(f"{k[:40]:40s} {aln[k][lo:hi + 1]}")


if __name__ == "__main__":
    main()
