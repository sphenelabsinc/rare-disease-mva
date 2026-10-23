"""Stage 2: evidence-annotated causal network.

Combines (a) the hand-curated causal edges in config/track2/causal_edges.tsv
(each with a PMID, an evidence tier and a demonstrated/inferred/speculative
status) with (b) context from STRING (high-confidence BUB1B partners) and
Reactome (pathways containing BUB1B). Writes:

  reports/track2/tables/s2_string_partners.tsv
  reports/track2/tables/s2_reactome_pathways.tsv
  reports/track2/tables/s2_network_edges.tsv   (curated + STRING, one schema)
  reports/track2/s2_network.mmd                (Mermaid diagram of curated edges)

Usage: python -I scripts/track2/s2_network.py
"""

import csv
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data/work/track2/cache"
OUT = ROOT / "reports/track2"
EDGES = ROOT / "config/track2/causal_edges.tsv"

STRING_URL = ("https://string-db.org/api/json/interaction_partners?identifiers=BUB1B"
              "&species=9606&required_score=700&limit=40&caller_identity=rare-disease-mva")
REACTOME_URL = "https://reactome.org/ContentService/data/mapping/UniProt/O60566/pathways?species=9606"

SCHEMA = ["source", "target", "relation", "module", "status", "tier", "pmids", "note", "origin"]


def fetch_json(url, dest):
    if not dest.exists():
        req = urllib.request.Request(url, headers={"User-Agent": "rare-disease-mva/0.1"})
        dest.write_bytes(urllib.request.urlopen(req, timeout=120).read())
    return json.loads(dest.read_text())


def main():
    tables = OUT / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    partners = fetch_json(STRING_URL, CACHE / "string_bub1b.json")
    with open(tables / "s2_string_partners.tsv", "w") as fh:
        fh.write("partner\tcombined_score\texperimental\tdatabase\ttextmining\n")
        for p in sorted(partners, key=lambda p: -p["score"]):
            fh.write(f"{p['preferredName_B']}\t{p['score']:.3f}\t{p['escore']:.3f}\t{p['dscore']:.3f}\t{p['tscore']:.3f}\n")

    pathways = fetch_json(REACTOME_URL, CACHE / "reactome_bub1b.json")
    with open(tables / "s2_reactome_pathways.tsv", "w") as fh:
        fh.write("stId\tname\n")
        for p in pathways:
            fh.write(f"{p['stId']}\t{p['displayName']}\n")

    with open(EDGES) as fh:
        curated = list(csv.DictReader(fh, delimiter="\t"))
    for e in curated:
        e["origin"] = "curated"
        missing = [k for k in SCHEMA if k != "origin" and not e.get(k)]
        assert not [k for k in missing if k not in ("note",)], (e, missing)

    string_edges = [dict(source="BUB1B", target=p["preferredName_B"], relation="interacts",
                         module="primary_context", status="database", tier="E4",
                         pmids="", note=f"STRING combined {p['score']:.2f}", origin="STRING")
                    for p in partners]

    with open(tables / "s2_network_edges.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=SCHEMA, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(curated + string_edges)

    # Mermaid diagram: solid = demonstrated, dashed = inferred, dotted = speculative
    arrow = {"demonstrated": "-->", "inferred": "-.->", "speculative": "-..->"}
    ids = {}

    def nid(name):
        return ids.setdefault(name, f"n{len(ids)}")

    lines = ["flowchart TD"]
    for e in curated:
        a, b = nid(e["source"]), nid(e["target"])
        lines.append(f'  {a}["{e["source"]}"] {arrow[e["status"]]}|{e["relation"]} {e["tier"]}| {b}["{e["target"]}"]')
    (OUT / "s2_network.mmd").write_text("\n".join(lines) + "\n")

    by_status = {}
    for e in curated:
        by_status[e["status"]] = by_status.get(e["status"], 0) + 1
    print(f"STRING partners: {len(partners)}; Reactome pathways: {len(pathways)}; "
          f"curated edges: {len(curated)} {by_status}")


if __name__ == "__main__":
    main()
