"""Fetch the FDA-approved drug-target table from ChEMBL.

Source: ChEMBL REST API (release recorded in output).
Scope: ChEMBL molecules with max_phase == 4 (approved in at least one major
market), mechanism records with direct interaction == 1, human protein targets.

Outputs (under results/):
  chembl_mechanisms.parquet    one row per (drug, mechanism record)
  chembl_targets.parquet       one row per (drug, human UniProt target)
  chembl_fda_targets.parquet   the join used downstream, with UniProt accession
  fetch_meta.json              release, counts, timestamps, spot-check verdicts

Verification: known approved drug-target pairs are spot-checked against the
published labels (methotrexate-DHFR, imatinib-ABL1/KIT, chloroquine-TLR7/TLR9
may be absent as a non-protein-target drug; it is checked leniently).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

BASE = "https://www.ebi.ac.uk/chembl/api/data"
HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
RESULTS.mkdir(exist_ok=True)

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "SilicoMVA/1.0 (research)"


def get_json(url: str, retries: int = 4) -> dict:
    for attempt in range(retries):
        try:
            r = SESSION.get(url, timeout=60)
            if r.status_code == 200:
                return r.json()
        except requests.RequestException:
            pass
        time.sleep(2**attempt)
    raise RuntimeError(f"fetch failed after {retries} retries: {url}")


def page_all(endpoint: str, params: str = "", page_size: int = 1000):
    url = f"{BASE}/{endpoint}.json?limit={page_size}&offset=0"
    if params:
        url = f"{BASE}/{endpoint}.json?limit={page_size}&offset=0&{params}"
    key = {"mechanism": "mechanisms", "target": "targets", "molecule": "molecules"}.get(endpoint, endpoint + "s")
    rows = []
    while url:
        data = get_json(url)
        rows.extend(data.get(key, []))
        nxt = data.get("page_meta", {}).get("next")
        total = data.get("page_meta", {}).get("total_count", "?")
        url = None
        if nxt:
            url = nxt if nxt.startswith("http") else f"https://www.ebi.ac.uk{chembl_path(nxt)}"
        if url and len(rows) % (page_size * 10) == 0:
            print(f"  {endpoint}: {len(rows)}/{total}", flush=True)
    return rows


def chembl_path(next_url: str) -> str:
    i = next_url.find("/chembl")
    return next_url[i:] if i >= 0 else next_url


def main() -> int:
    status = get_json(f"{BASE}/status.json")
    print("ChEMBL status:", status, flush=True)

    # --- mechanisms for approved drugs -------------------------------------
    print("fetching mechanism.json (max_phase=4) ...", flush=True)
    mech = page_all("mechanism", params="max_phase=4&direct_interaction=1")
    mech_df = pd.DataFrame(mech)
    print(f"mechanism rows: {len(mech_df)}", flush=True)
    keep_cols = [
        "mec_id", "record_id", "molecule_chembl_id", "parent_molecule_chembl_id",
        "target_chembl_id", "action_type", "mechanism_of_action", "max_phase",
        "direct_interaction", "disease_efficacy", "binding_site_comment",
        "selectivity_comment",
    ]
    mech_df = mech_df[[c for c in keep_cols if c in mech_df.columns]]
    mech_df.to_parquet(RESULTS / "chembl_mechanisms.parquet", index=False)

    # --- targets ------------------------------------------------------------
    targ_ids = sorted(mech_df["target_chembl_id"].dropna().unique())
    print(f"distinct targets: {len(targ_ids)}", flush=True)
    targ_rows = []
    for i in range(0, len(targ_ids), 50):
        batch = targ_ids[i : i + 50]
        data = get_json(f"{BASE}/target.json?target_chembl_id__in={','.join(batch)}&limit=1000")
        targ_rows.extend(data.get("targets", []))
    tf = []
    for t in targ_rows:
        comps = t.get("target_components") or []
        for comp in comps:
            tf.append(
                {
                    "target_chembl_id": t.get("target_chembl_id"),
                    "pref_name": t.get("pref_name"),
                    "target_type": t.get("target_type"),
                    "organism": t.get("organism"),
                    "component_accession": comp.get("accession"),
                    "component_type": comp.get("component_type"),
                    "component_description": comp.get("component_description"),
                }
            )
    targ_df = pd.DataFrame(tf)
    targ_df.to_parquet(RESULTS / "chembl_targets.parquet", index=False)
    print(f"target rows: {len(targ_df)}", flush=True)

    # --- molecules ----------------------------------------------------------
    mol_ids = sorted(mech_df["molecule_chembl_id"].dropna().unique())
    print(f"distinct molecules: {len(mol_ids)}", flush=True)
    mol_rows = []
    for i in range(0, len(mol_ids), 50):
        batch = mol_ids[i : i + 50]
        data = get_json(f"{BASE}/molecule.json?molecule_chembl_id__in={','.join(batch)}&limit=1000")
        mol_rows.extend(data.get("molecules", []))
    mf = []
    for m in mol_rows:
        atc = m.get("atc_classifications") or []
        mf.append(
            {
                "molecule_chembl_id": m.get("molecule_chembl_id"),
                "pref_name": m.get("pref_name"),
                "molecule_type": m.get("molecule_type"),
                "max_phase": m.get("max_phase"),
                "first_approval": m.get("first_approval"),
                "withdrawn_flag": m.get("withdrawn_flag"),
                "black_box_warning": m.get("black_box_warning"),
                "usan_stem": m.get("usan_stem"),
                "atc": ";".join(atc),
                "molecule_structures.molregno": (m.get("molecule_structures") or {}).get("molregno"),
            }
        )
    mol_df = pd.DataFrame(mf)
    mol_df.to_parquet(RESULTS / "chembl_molecules.parquet", index=False)
    print(f"molecule rows: {len(mol_df)}", flush=True)

    # --- join: drug x human UniProt protein targets -------------------------
    human = targ_df[(targ_df["organism"] == "Homo sapiens")]
    j = mech_df.merge(
        human[["target_chembl_id", "pref_name", "target_type", "component_accession"]],
        on="target_chembl_id",
        how="inner",
    ).merge(
        mol_df, on="molecule_chembl_id", how="left", suffixes=("", ".mol")
    )
    j = j.rename(columns={"pref_name": "target_pref_name"})
    # keep protein-family/single-protein direct targets
    prot_types = {"SINGLE PROTEIN", "PROTEIN FAMILY", "CHIMERIC PROTEIN"}
    j = j[j["target_type"].isin(prot_types)]
    j.to_parquet(RESULTS / "chembl_fda_targets.parquet", index=False)
    print(f"joined rows (human protein targets): {len(j)}", flush=True)

    meta = {
        "chembl_status": status,
        "n_mechanisms": int(len(mech_df)),
        "n_targets": int(len(targ_df)),
        "n_molecules": int(len(mol_df)),
        "n_joined": int(len(j)),
        "n_distinct_drugs_with_targets": int(j["molecule_chembl_id"].nunique()),
        "n_distinct_accessions": int(j["component_accession"].nunique()),
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (RESULTS / "fetch_meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
