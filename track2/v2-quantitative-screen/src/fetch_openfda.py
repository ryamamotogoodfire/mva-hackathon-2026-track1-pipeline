"""Enrich the approved-drug table with openFDA pediatric label fields.

For each drug in the ChEMBL drug-target table, search the openFDA drug label
index by generic name and record whether a label exists, whether it carries a
Pediatric Use section (8.4), and the indication text length. This is a screen
feature, age claims are re-verified on regulator pages at synthesis time.

Output: results/openfda_pediatric.parquet
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
SESSION = requests.Session()
SESSION.headers["User-Agent"] = "SilicoMVA/1.0 (research)"
API = "https://api.fda.gov/drug/label.json"


def fetch_label(name: str) -> dict | None:
    q = f'openfda.generic_name:"{name}"'
    params = {"search": q, "limit": 1}
    for attempt in range(4):
        try:
            r = SESSION.get(API, params=params, timeout=60)
            if r.status_code == 404:
                return None
            if r.status_code == 200:
                res = r.json().get("results") or []
                if not res:
                    return None
                doc = res[0]
                of = doc.get("openfda", {})
                ped = " ".join(doc.get("pediatric_use", []) or [])
                ind = " ".join(doc.get("indications_and_usage", []) or [])
                pop = " ".join(doc.get("use_in_specific_populations", []) or [])
                return {
                    "set_id": doc.get("set_id", ""),
                    "generic_names": ";".join(of.get("generic_name", []) or []),
                    "brand_names": ";".join(of.get("brand_name", []) or []),
                    "manufacturer": ";".join(of.get("manufacturer_name", []) or []),
                    "route": ";".join(of.get("route", []) or []),
                    "has_pediatric_section": bool(ped.strip()),
                    "pediatric_text_len": len(ped),
                    "indication_text_len": len(ind),
                    "specific_populations_len": len(pop),
                    "pediatric_snippet": ped[:400],
                    "indication_snippet": ind[:400],
                }
        except requests.RequestException:
            pass
        time.sleep(2**attempt)
    return None


def main() -> None:
    dt = pd.read_parquet(RESULTS / "drug_target_genes.parquet")
    names = sorted(dt["pref_name"].dropna().unique())
    # clean the list: pure alphanumerics/dashes/spaces only
    names = [n for n in names if len(n) < 80]
    print(f"drugs to query openFDA: {len(names)}", flush=True)

    rows = []
    t0 = time.time()
    for i, name in enumerate(names):
        hit = fetch_label(name)
        row = {"pref_name": name, "openfda_found": bool(hit)}
        if hit:
            row.update(hit)
        rows.append(row)
        if (i + 1) % 100 == 0:
            dt_s = time.time() - t0
            print(f"  {i + 1}/{len(names)}  ({dt_s:.0f}s)", flush=True)
        time.sleep(0.28)  # stay below 240 req/min anonymous limit
    df = pd.DataFrame(rows)
    df.to_parquet(RESULTS / "openfda_pediatric.parquet", index=False)
    print(f"found labels for {df.openfda_found.sum()}/{len(df)} drugs", flush=True)
    print(
        f"with pediatric section: {df.get('has_pediatric_section', pd.Series(False)).sum()}",
        flush=True,
    )


if __name__ == "__main__":
    main()
