"""Literature retrieval: run fixed PubMed queries and cache records.

Each query in QUERIES is run through NCBI E-utilities (esearch, then efetch for
abstracts). Records are cached as JSON under data/work/track2/cache/pubmed/, and
a deduplicated index is written to reports/track2/tables/s0_pubmed_hits.tsv
(query id, PMID, year, journal, first author, title). The query list is fixed
so retrieval is reproducible apart from new publications.

Usage: python -I scripts/track2/s0_pubmed.py
"""

import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data/work/track2/cache/pubmed"
OUT = ROOT / "reports/track2/tables"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

QUERIES = {
    "mva_bub1b": '("mosaic variegated aneuploidy" OR "premature chromatid separation") AND (BUB1B OR BUBR1)',
    "mva_review": '"mosaic variegated aneuploidy" AND review[pt]',
    "mva_cancer_rx": '"mosaic variegated aneuploidy" AND (chemotherapy OR rhabdomyosarcoma OR Wilms OR leukemia OR treatment)',
    "allele_models": '(BubR1 OR BUB1B) AND (L1012P OR "allelic" OR "allele combination") AND mice',
    "pseudokinase": '(BubR1 OR BUBR1) AND (pseudokinase OR "kinase domain") AND (KARD OR PP2A OR "B56" OR stability OR conformation)',
    "kard_pp2a": '(BubR1 OR BUBR1) AND KARD',
    "missense_function": '(BUB1B OR BUBR1) AND (missense OR mutant OR mutation) AND (stability OR degradation OR "protein level") AND (MVA OR aneuploidy)',
    "sirt2": '(SIRT2 OR "sirtuin 2") AND (BubR1 OR BUB1B)',
    "nad": '(NAD OR "nicotinamide mononucleotide" OR "nicotinamide riboside") AND (BubR1 OR BUB1B)',
    "senescence": '(BubR1 OR BUB1B) AND (senescence OR senolytic OR p16Ink4a OR INK-ATTAC)',
    "overexpression": '(BubR1 OR BUB1B) AND (overexpression OR "high levels") AND (lifespan OR aging OR tumor)',
    "mtor": '(BubR1 OR BUB1B OR "mosaic variegated aneuploidy") AND (mTOR OR rapamycin OR sirolimus OR everolimus)',
    "aneuploid_selective": '(aneuploid OR aneuploidy) AND ("selectively" OR "selective") AND (kill OR inhibit OR target) AND (AICAR OR 17-AAG OR chloroquine OR "energy stress" OR drug)',
    "cin_cgas": '("chromosomal instability" OR aneuploidy OR micronuclei) AND (cGAS OR STING) AND (senescence OR inflammation)',
    "cin_taxane": '(BubR1 OR BUB1B OR "spindle assembly checkpoint") AND (vincristine OR paclitaxel OR taxol OR "microtubule poison") AND (sensitivity OR resistance)',
    "metformin_exercise": '(BubR1 OR BUB1B) AND (metformin OR exercise OR "dietary restriction" OR "caloric restriction")',
    "p53_aneuploidy": '(BubR1 OR BUB1B) AND (p53 OR p21) AND (aneuploidy OR senescence OR apoptosis)',
    "cep57_trip13": '("mosaic variegated aneuploidy") AND (CEP57 OR TRIP13)',
    "bubr1_acetylation": '(BubR1 OR BUB1B) AND (acetylation OR ubiquitination OR degradation) AND (PCAF OR APC OR mitosis)',
    "aneuploidy_proteotoxic": '(aneuploidy OR aneuploid) AND (proteotoxic OR proteasome OR HSP90 OR autophagy) AND (tolerance OR stress)',
    # Targeted lookups for foundational papers the broad queries missed
    "baker2004": 'Baker[au] AND BubR1 AND 2004[dp] AND progeroid',
    "baker2013_ncb": '"Increased expression of BubR1 protects against aneuploidy and cancer"',
    "tang2011": 'Tang[au] AND aneuploidy-selective AND 2011[dp]',
    "bubr1_dasatinib": 'BubR1 AND dasatinib',
    "mva_vincristine": '"mosaic variegated aneuploidy" AND vincristine',
}
RETMAX = 60


def get(url, dest=None, tries=4):
    for i in range(tries):
        try:
            data = urllib.request.urlopen(url, timeout=60).read()
            if dest:
                dest.write_bytes(data)
            return data
        except Exception:
            time.sleep(2 * (i + 1))
    raise RuntimeError(url)


def parse(xml_bytes):
    recs = []
    for art in ET.fromstring(xml_bytes).iter("PubmedArticle"):
        pmid = art.findtext(".//PMID")
        title = "".join(art.find(".//ArticleTitle").itertext())
        abstract = " ".join("".join(t.itertext()) for t in art.iter("AbstractText"))
        year = art.findtext(".//JournalIssue/PubDate/Year") or (art.findtext(".//JournalIssue/PubDate/MedlineDate") or "")[:4]
        journal = art.findtext(".//Journal/ISOAbbreviation") or ""
        au = art.find(".//AuthorList/Author")
        first = (au.findtext("LastName") or au.findtext("CollectiveName") or "") if au is not None else ""
        types = [t.text for t in art.iter("PublicationType")]
        recs.append(dict(pmid=pmid, year=year, journal=journal, first_author=first,
                         title=title, abstract=abstract, types=types))
    return recs


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    seen = {}
    for qid, q in QUERIES.items():
        dest = CACHE / f"{qid}.json"
        if not dest.exists():
            term = urllib.parse.quote(q)
            ids = json.loads(get(f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax={RETMAX}&sort=relevance&term={term}"))["esearchresult"]["idlist"]
            recs = parse(get(f"{EUTILS}/efetch.fcgi?db=pubmed&retmode=xml&id={','.join(ids)}")) if ids else []
            dest.write_text(json.dumps(dict(query=q, records=recs), indent=1))
            time.sleep(0.4)
        recs = json.loads(dest.read_text())["records"]
        print(f"{qid}: {len(recs)}")
        for r in recs:
            seen.setdefault(r["pmid"], dict(r, queries=[]))["queries"].append(qid)
    with open(OUT / "s0_pubmed_hits.tsv", "w") as fh:
        fh.write("pmid\tyear\tjournal\tfirst_author\tqueries\ttitle\n")
        for r in sorted(seen.values(), key=lambda r: r["year"], reverse=True):
            fh.write(f"{r['pmid']}\t{r['year']}\t{r['journal']}\t{r['first_author']}\t{','.join(r['queries'])}\t{r['title']}\n")
    (CACHE / "all_records.json").write_text(json.dumps(seen, indent=1))
    print(f"{len(seen)} unique records")


if __name__ == "__main__":
    main()
