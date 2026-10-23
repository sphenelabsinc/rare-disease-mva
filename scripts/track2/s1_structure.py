"""Stage 1: structural context of BUBR1 p.Asn1002Lys in the AlphaFold model.

Downloads the AlphaFold model of O60566 and reports, for the patient residue and
the known MVA1 pseudokinase-domain missense residues: pLDDT, relative solvent
accessibility, and neighbouring residues. Writes
reports/track2/tables/s1_structure.tsv.

Usage: python -I scripts/track2/s1_structure.py
"""

import json
import urllib.request
from pathlib import Path

from Bio.PDB import NeighborSearch, PDBParser
from Bio.PDB.SASA import ShrakeRupley

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data/work/track2/cache"
OUT = ROOT / "reports/track2/tables"

UNIPROT = "O60566"
PATIENT = 1002
# MVA1 missense in the C-terminal domain (UniProt natural variants)
MVA_MISSENSE = {814: "R814H", 844: "L844F", 909: "I909T", 921: "Q921H", 1012: "L1012P"}
PLK1_SITE = 1008  # T1008, PLK1 phosphosite

# Max ASA (Tien et al. 2013, theoretical) for relative accessibility
MAX_ASA = dict(ALA=129, ARG=274, ASN=195, ASP=193, CYS=167, GLN=225, GLU=223,
               GLY=104, HIS=224, ILE=197, LEU=201, LYS=236, MET=224, PHE=240,
               PRO=159, SER=155, THR=172, TRP=285, TYR=263, VAL=174)


def fetch(url, dest):
    if not dest.exists():
        urllib.request.urlretrieve(url, dest)
    return dest


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    meta_path = fetch(f"https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT}",
                      CACHE / f"afdb_{UNIPROT}.json")
    meta = json.loads(meta_path.read_text())[0]
    pdb_path = fetch(meta["pdbUrl"], CACHE / f"AF-{UNIPROT}.pdb")

    model = PDBParser(QUIET=True).get_structure("bubr1", pdb_path)[0]
    chain = model["A"]
    ShrakeRupley().compute(model, level="R")
    atoms = [a for a in model.get_atoms()]
    ns = NeighborSearch(atoms)

    def describe(resnum, label):
        res = chain[resnum]
        plddt = res["CA"].get_bfactor()
        rsa = res.sasa / MAX_ASA[res.get_resname()]
        near = set()
        for atom in res:
            for other in ns.search(atom.coord, 4.5, level="R"):
                n = other.id[1]
                if abs(n - resnum) > 2:
                    near.add(n)
        return dict(residue=f"{res.get_resname()}{resnum}", label=label,
                    plddt=round(plddt, 1), rsa=round(rsa, 2),
                    burial="buried" if rsa < 0.2 else "intermediate" if rsa < 0.4 else "exposed",
                    contacts_4p5A=",".join(str(n) for n in sorted(near)))

    rows = [describe(PATIENT, "patient p.Asn1002Lys")]
    rows += [describe(k, f"MVA1 {v}") for k, v in MVA_MISSENSE.items()]
    rows.append(describe(PLK1_SITE, "PLK1 site T1008"))

    ca = lambda n: chain[n]["CA"]
    dist = {f"N1002-{v}": round(ca(PATIENT) - ca(k), 1) for k, v in MVA_MISSENSE.items()}
    dist["N1002-T1008"] = round(ca(PATIENT) - ca(PLK1_SITE), 1)

    # Local pLDDT over the C-terminal domain, for context
    dom = [chain[n]["CA"].get_bfactor() for n in range(766, 1051) if n in chain]
    with open(OUT / "s1_structure.tsv", "w") as fh:
        cols = list(rows[0])
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(OUT / "s1_structure_distances.tsv", "w") as fh:
        fh.write("pair\tca_distance_A\n")
        for k, v in dist.items():
            fh.write(f"{k}\t{v}\n")
    print(f"AlphaFold model version {meta.get('latestVersion')}; "
          f"domain 766-1050 mean pLDDT {sum(dom) / len(dom):.1f}")
    for r in rows:
        print(r)
    print(dist)


if __name__ == "__main__":
    main()
