"""Stage 3: reproducible approved-drug longlist from network targets.

For each gene in config/track2/intervention_targets.tsv and each high-confidence
STRING partner of BUB1B, look up drugs with a curated mechanism of action on
that target in ChEMBL (via /mechanism) and keep molecules at max_phase 4
(approved). Open Targets known-drug evidence is added as a cross-check.

Writes reports/track2/tables/s3_longlist.tsv with one row per (drug, target).

Usage: python -I scripts/track2/s3_drugs.py
"""

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data/work/track2/cache/drugs"
OUT = ROOT / "reports/track2/tables"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
OT = "https://api.platform.opentargets.org/api/v4/graphql"
UA = {"User-Agent": "rare-disease-mva/0.1", "Content-Type": "application/json"}


def get_json(url, dest, data=None):
    if not dest.exists():
        for i in range(4):
            try:
                req = urllib.request.Request(url, data=data, headers=UA)
                dest.write_bytes(urllib.request.urlopen(req, timeout=120).read())
                break
            except urllib.error.HTTPError as e:
                if 400 <= e.code < 500 and e.code != 429:
                    raise
                time.sleep(3 * (i + 1))
            except Exception:
                time.sleep(3 * (i + 1))
        else:
            raise RuntimeError(url)
        time.sleep(0.3)
    return json.loads(dest.read_text())


def chembl_targets(gene):
    q = urllib.parse.quote(gene)
    d = get_json(f"{CHEMBL}/target.json?target_components__target_component_synonyms__component_synonym__iexact={q}"
                 f"&organism=Homo%20sapiens&limit=100", CACHE / f"target_{gene}.json")
    # single protein plus complexes/families that contain the gene
    return [t["target_chembl_id"] for t in d["targets"]
            if any(s["component_synonym"].upper() == gene and s["syn_type"] == "GENE_SYMBOL"
                   for c in t["target_components"] for s in c["target_component_synonyms"])]


def chembl_mechanisms(tid):
    return get_json(f"{CHEMBL}/mechanism.json?target_chembl_id={tid}&limit=1000",
                    CACHE / f"mech_{tid}.json")["mechanisms"]


def molecules(ids):
    out = {}
    ids = sorted(set(ids))
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        d = get_json(f"{CHEMBL}/molecule.json?molecule_chembl_id__in={','.join(chunk)}&limit=50",
                     CACHE / f"mol_{chunk[0]}_{len(chunk)}.json")
        for m in d["molecules"]:
            out[m["molecule_chembl_id"]] = m
    return out


def ensembl_id(gene):
    d = get_json(f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{gene}?content-type=application/json",
                 CACHE / f"ens_{gene}.json")
    return next((x["id"] for x in d if x["id"].startswith("ENSG")), None)


def ot_drugs(ensg):
    # Platform schema as of 2026 (knownDrugs was replaced by drugAndClinicalCandidates)
    query = ("query($id:String!){target(ensemblId:$id){approvedSymbol drugAndClinicalCandidates{rows{"
             "maxClinicalStage drug{id name}}}}}")
    body = json.dumps({"query": query, "variables": {"id": ensg}}).encode()
    d = get_json(OT, CACHE / f"ot_{ensg}.json", data=body)
    t = (d.get("data") or {}).get("target") or {}
    return [r for r in ((t.get("drugAndClinicalCandidates") or {}).get("rows") or []) if r.get("drug")]


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(ROOT / "config/track2/intervention_targets.tsv") as fh:
        targets = {r["gene"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with open(OUT / "s2_string_partners.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            targets.setdefault(r["partner"], dict(gene=r["partner"], node="BUB1B interactor (STRING)",
                                                  desired_direction="context", rationale=f"STRING {r['combined_score']}"))
    targets.setdefault("BUB1B", dict(gene="BUB1B", node="BUBR1 abundance", desired_direction="increase",
                                     rationale="Index gene"))

    rows, mech_mols = [], []
    for gene, t in targets.items():
        for tid in chembl_targets(gene):
            for m in chembl_mechanisms(tid):
                mech_mols.append(m["molecule_chembl_id"])
                rows.append(dict(gene=gene, target_chembl_id=tid, molecule=m["molecule_chembl_id"],
                                 action=m["action_type"], mechanism=m["mechanism_of_action"]))
    mols = molecules(mech_mols)

    ot = {}
    for gene in targets:
        try:
            ensg = ensembl_id(gene)
        except RuntimeError:
            ensg = None  # Ensembl lookup failed; Open Targets cross-check skipped for this gene
        if ensg:
            for r in ot_drugs(ensg):
                ot.setdefault((gene, r["drug"]["name"].upper()), r)

    out, seen = [], set()
    for r in rows:
        m = mols.get(r["molecule"], {})
        if str(m.get("max_phase")) not in ("4", "4.0"):
            continue
        name = (m.get("pref_name") or r["molecule"]).upper()
        key = (name, r["gene"], r["action"])
        if key in seen:
            continue
        seen.add(key)
        t = targets[r["gene"]]
        out.append(dict(drug=name, molecule=r["molecule"], gene=r["gene"], node=t["node"],
                        desired_direction=t["desired_direction"], action=r["action"],
                        mechanism=r["mechanism"],
                        first_approval=m.get("first_approval") or "",
                        open_targets_confirms="yes" if (r["gene"], name) in ot else "no"))

    def flag(r):
        d, a = r["desired_direction"], r["action"]
        if d.startswith("avoid") and ((d == "avoid_decrease" and a in ("INHIBITOR", "ANTAGONIST", "NEGATIVE MODULATOR"))
                                      or (d == "avoid_increase")):
            return "AVOID"
        if d == "decrease" and a in ("INHIBITOR", "ANTAGONIST", "NEGATIVE ALLOSTERIC MODULATOR", "NEGATIVE MODULATOR", "BLOCKER"):
            return "direction_match"
        if d == "increase" and a in ("ACTIVATOR", "AGONIST", "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE MODULATOR"):
            return "direction_match"
        if d == "context":
            return "context_only"
        return "direction_mismatch"

    for r in out:
        r["direction_flag"] = flag(r)
    out.sort(key=lambda r: (r["direction_flag"], r["node"], r["drug"]))
    with open(OUT / "s3_longlist.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    from collections import Counter
    print(Counter(r["direction_flag"] for r in out))
    print(f"{len(out)} approved drug-target rows; {len({r['drug'] for r in out})} unique drugs")


if __name__ == "__main__":
    main()
