#!/usr/bin/env python3
"""Regulatory status for every named agent (ChEMBL REST).

Resolution order per name: exact preferred-name match, then synonym substring
match. For each resolved molecule we record max_phase, first_approval,
withdrawn_flag, black-box flag, molecule type, and the named regulator cross
references (EMA EPAR, DailyMed) ChEMBL carries. Missing molecule or null
max_phase means 'not approved / not registered' for grading purposes; the
final report checks the regulator's own page for each positive claim.

Output: results/regulatory_status.json
"""
import json
import time
import urllib.parse
import urllib.request

from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"

DRUGS = [
    "ataluren", "gentamicin", "amikacin", "amlexanox",
    "sirolimus", "everolimus", "temsirolimus",
    "metformin", "phenformin", "aspirin", "bortezomib",
    "reversine", "barasertib", "cycloheximide",
]

BASE = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"


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
    return {"molecules": [], "_error": str(last)}


def summarize(m):
    return {
        "chembl_id": m.get("molecule_chembl_id"),
        "pref_name": m.get("pref_name"),
        "max_phase": m.get("max_phase"),
        "first_approval": m.get("first_approval"),
        "withdrawn_flag": m.get("withdrawn_flag"),
        "black_box_warning": m.get("black_box_warning"),
        "molecule_type": m.get("molecule_type"),
        "availability_type": m.get("availability_type"),
        "cross_refs": [{"src": x.get("xref_src"), "id": x.get("xref_id"), "name": x.get("xref_name")}
                       for x in (m.get("cross_references") or [])],
        "atc": m.get("atc_classifications"),
    }


out = {}
for d in DRUGS:
    js = get(BASE + "?" + urllib.parse.urlencode({"pref_name__iexact": d}))
    ms = js.get("molecules", [])
    via = "pref_name"
    if not ms:
        js = get(BASE + "?" + urllib.parse.urlencode({"molecule_synonyms__molecule_synonym__icontains": d, "limit": 8}))
        ms = js.get("molecules", [])
        via = "synonym"
    recs = [summarize(m) for m in ms]
    out[d] = {"resolved_via": via, "n": len(recs), "entries": recs}
    best = recs[0] if recs else None
    if best:
        print(f"{d:16s} via={via:9s} n={len(recs):2d} best={best['pref_name']} phase={best['max_phase']} appr={best['first_approval']} withdrawn={best['withdrawn_flag']}")
    else:
        print(f"{d:16s} NOT FOUND")

RESULTS.joinpath("regulatory_status.json").write_text(json.dumps(out, indent=2) + "\n")
print("OK wrote results/regulatory_status.json")
