"""Fetch approved small-molecule structures from ChEMBL for the docking screen.

Writes results/approved_ligands.parquet with molecule_chembl_id, pref_name, smiles,
mw, inchikey, and the dedup key. Screen filter: max_phase=4 small molecules with a
canonical SMILES and 150 <= MW <= 600 (the user-specified docking window).
"""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

BASE = "https://www.ebi.ac.uk/chembl/api/data"
OUT = Path(__file__).resolve().parent.parent / "results"


def get(url: str, tries: int = 4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def main() -> None:
    rows = []
    url = (f"{BASE}/molecule.json?max_phase=4&limit=1000"
           "&only=molecule_chembl_id,pref_name,molecule_type,max_phase,first_approval,"
           "withdrawn_flag,molecule_structures,molecule_properties")
    page = 0
    while url:
        d = get(url)
        for m in d["molecules"]:
            st = m.get("molecule_structures") or {}
            pr = m.get("molecule_properties") or {}
            rows.append({
                "molecule_chembl_id": m["molecule_chembl_id"],
                "pref_name": m.get("pref_name"),
                "molecule_type": m.get("molecule_type"),
                "first_approval": m.get("first_approval"),
                "withdrawn_flag": m.get("withdrawn_flag"),
                "smiles": st.get("canonical_smiles"),
                "inchikey": st.get("standard_inchi_key"),
                "mw": pd.to_numeric(pr.get("full_mwt"), errors="coerce"),
                "heavy_atoms": pd.to_numeric(pr.get("heavy_atoms"), errors="coerce"),
                "alogp": pd.to_numeric(pr.get("alogp"), errors="coerce"),
            })
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
        page += 1
        print(f"page {page}: {len(rows)} molecules", flush=True)

    df = pd.DataFrame(rows)
    df.to_parquet(OUT / "approved_molecules_full.parquet", index=False)
    sel = df[(df.molecule_type == "Small molecule") & df.smiles.notna() & df.mw.between(150, 600)].copy()
    sel["stem"] = sel.inchikey.fillna("").str.split("-").str[0]
    sel = sel.sort_values(["stem", "first_approval"], na_position="last")
    sel = sel[~(sel.stem.duplicated() & sel.stem.ne(""))]
    sel = sel[sel.smiles.str.len() < 300]
    sel.to_parquet(OUT / "approved_ligands.parquet", index=False)
    print(f"all approved molecules: {len(df)}; docking set (small molecule, 150<=MW<=600, dedup by InChIKey block): {len(sel)}")
    print(sel[["molecule_chembl_id", "pref_name", "mw"]].head(8).to_string())


if __name__ == "__main__":
    main()
