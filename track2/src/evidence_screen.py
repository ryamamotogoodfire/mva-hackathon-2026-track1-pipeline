#!/usr/bin/env python3
"""Systematic candidate screen machinery for Track 2.

Two evidence lanes, both public and reproducible:

1. Literature lane (Europe PMC REST). Fixed queries per candidate class; the
   top hits by citation count are stored with PMID/DOI/title/journal/year so
   every retained mechanistic claim can carry a checkable citation.

2. Regulatory lane. For each named agent we record
   - ChEMBL max trial phase + preferred name (EBI REST)
   - DailyMed label presence (US, SPL services)
   ChEMBL/DailyMed absence is reported as absence; final approval claims in
   the report are always stated against a named jurisdiction and verified on
   the regulator's own page (EMA EPAR / FDA label) during curation.

Outputs results/evidence_raw.json (raw API payloads, for audit) and
results/literature_hits.json (curated query -> hit lists used by curation).
"""
import json
import time
import urllib.parse
import urllib.request

from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
RESULTS.mkdir(exist_ok=True)

QUERIES = {
    "q_bubr1_hypomorph": '((TITLE_ABS:"BubR1" OR TITLE_ABS:"BUB1B") AND (TITLE_ABS:"hypomorph" OR TITLE_ABS:"insufficiency")) AND PUB_YEAR:[2000 TO 2026]',
    "q_bubr1_kinase_activity": '(TITLE_ABS:"BubR1" AND (TITLE_ABS:"kinase-dead" OR TITLE_ABS:"kinase activity" OR TITLE_ABS:"kinase-dead"))',
    "q_nmd_50nt_rule": '(TITLE_ABS:"nonsense-mediated mRNA decay" AND (TITLE_ABS:"50 nucleotides" OR TITLE_ABS:"exon junction complex" OR TITLE_ABS:"rule"))',
    "q_ataluren_dmd": '(TITLE_ABS:"ataluren" AND (TITLE_ABS:"Duchenne" OR TITLE_ABS:"nonsense"))',
    "q_aminoglycoside_readthrough": '((TITLE_ABS:"gentamicin" OR TITLE_ABS:"aminoglycoside") AND TITLE_ABS:"readthrough" AND (TITLE_ABS:"cystic fibrosis" OR TITLE_ABS:"CFTR" OR TITLE_ABS:"nonsense"))',
    "q_amlexanox_nmd": '(TITLE_ABS:"amlexanox" AND (TITLE_ABS:"nonsense" OR TITLE_ABS:"nonsense-mediated mRNA decay" OR TITLE_ABS:"readthrough"))',
    "q_aneuploidy_metabolic_fitness": '(TITLE_ABS:"aneuploid" AND (TITLE_ABS:"metabolic" OR TITLE_ABS:"energy stress" OR TITLE_ABS:"fitness"))',
    "q_aneuploidy_proteotoxic": '(TITLE_ABS:"aneuploid" AND (TITLE_ABS:"proteotoxic" OR TITLE_ABS:"proteostasis" OR TITLE_ABS:"Hsp90" OR TITLE_ABS:"autophagy"))',
    "q_selective_killing_aneuploid": '(TITLE_ABS:"aneuploid" AND (TITLE_ABS:"selective" OR TITLE_ABS:"synthetic lethal" OR TITLE_ABS:"target") AND TITLE_ABS:"cell")',
    "q_r_mva_syndrome_review": '(TITLE_ABS:"mosaic variegated aneuploidy")',
    "q_everolimus_tsc_pediatric": '(TITLE_ABS:"everolimus" AND TITLE_ABS:"tuberous sclerosis" AND (TITLE_ABS:"pediatric" OR TITLE_ABS:"children" OR TITLE_ABS:"subependymal"))',
    "q_metformin_pediatric": '(TITLE_ABS:"metformin" AND (TITLE_ABS:"children" OR TITLE_ABS:"pediatric" OR TITLE_ABS:"adolescents") AND TITLE_ABS:"type 2 diabetes")',
    "q_metformin_cin": '(TITLE_ABS:"metformin" AND (TITLE_ABS:"aneuploid" OR TITLE_ABS:"chromosomal instability"))',
    "q_aspirin_lynch_capp2": '(TITLE_ABS:"aspirin" AND (TITLE_ABS:"Lynch" OR TITLE_ABS:"CAPP2"))',
    "q_mps1_inhibitor_trial": '((TITLE_ABS:"Mps1" OR TITLE_ABS:"TTK") AND TITLE_ABS:"inhibitor" AND (TITLE_ABS:"trial" OR TITLE_ABS:"phase 1" OR TITLE_ABS:"clinical"))',
    "q_mva_surveillance": '((TITLE_ABS:"mosaic variegated aneuploidy" OR TITLE_ABS:"BUB1B") AND (TITLE_ABS:"surveillance" OR TITLE_ABS:"screening" OR TITLE_ABS:"management" OR TITLE_ABS:"tumour" OR TITLE_ABS:"tumor"))',
    "q_readthrough_context": '(TITLE_ABS:"stop codon" AND (TITLE_ABS:"readthrough" OR TITLE_ABS:"context") AND (TITLE_ABS:"UGA" OR TITLE_ABS:"tetranucleotide" OR TITLE_ABS:"termination"))',
    "q_hsp90_aneuploidy_buffer": '(TITLE_ABS:"Hsp90" AND TITLE_ABS:"aneuploid")',
}

DRUGS = [
    "ataluren", "gentamicin", "amikacin", "amlexanox",
    "sirolimus", "everolimus", "temsirolimus",
    "metformin", "phenformin", "aspirin", "bortezomib",
    "reversine", "barasertib",
]

PMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CHEMBL_SEARCH = "https://www.ebi.ac.uk/chembl/api/data/search.json"
CHEMBL_MOLECULE = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
DAILYMED = "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json"


def get(url, tries=4):
    req = urllib.request.Request(url, headers={"User-Agent": "silico-track2/1.0", "Accept": "application/json"})
    last = None
    for a in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode("utf-8", errors="replace"))
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 + 3 * a)
    return {"_error": str(last)}


raw = {"queries": {}, "chEMBL": {}, "dailymed": {}}
lit = {}

for key, q in QUERIES.items():
    url = PMC + "?" + urllib.parse.urlencode({
        "query": q, "format": "json", "pageSize": 6, "sort": "CITED desc",
        "resultType": "core",
    })
    js = get(url)
    raw["queries"][key] = js
    hits = []
    for r in js.get("resultList", {}).get("result", []):
        hits.append({
            "pmid": r.get("pmid"), "doi": r.get("doi"),
            "pmcid": r.get("pmcid"),
            "title": r.get("title"), "journal": r.get("journalInfo", {}).get("journal", {}).get("title"),
            "year": r.get("pubYear"), "citedByCount": r.get("citedByCount"),
            "authors": ", ".join(a.get("fullName", "") for a in r.get("authorList", {}).get("author", [])[:4]),
        })
    lit[key] = {"query": q, "hitCount": js.get("hitCount"), "top": hits}
    print(f"{key}: {js.get('hitCount')} hits")

for d in DRUGS:
    js = get(CHEMBL_SEARCH + "?" + urllib.parse.urlencode({"q": d, "limit": 5}))
    raw["chEMBL"][d] = js
    dl = get(DAILYMED + "?" + urllib.parse.urlencode({"drug_name": d, "pagesize": 3}))
    raw["dailymed"][d] = dl
    n_labels = len(dl.get("data", [])) if isinstance(dl, dict) else 0
    print(f"{d}: dailymed labels={n_labels}")

RESULTS.joinpath("evidence_raw.json").write_text(json.dumps(raw, indent=2) + "\n")
RESULTS.joinpath("literature_hits.json").write_text(json.dumps(lit, indent=2) + "\n")
print("OK wrote results/evidence_raw.json results/literature_hits.json")
