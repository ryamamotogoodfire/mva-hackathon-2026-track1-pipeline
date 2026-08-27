"""Map ChEMBL target UniProt accessions to HGNC gene symbols.

Used to bring the drug-target table onto the same gene-symbol axis as the
BioGRID interactome and the disease module. Adds per-drug target gene lists.

Outputs:
  results/accession_gene_map.parquet   accession -> gene_symbol, protein name
  results/drug_target_genes.parquet    one row per (drug, gene), deduped
  results/spot_checks.json             named verification verdicts
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

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"


def fetch_batch(accs: list[str]) -> list[dict]:
    q = " OR ".join(f"accession:{a}" for a in accs)
    params = {
        "query": q,
        "fields": "accession,gene_primary,protein_name",
        "format": "tsv",
        "size": 500,
    }
    for attempt in range(4):
        try:
            r = SESSION.get(UNIPROT, params=params, timeout=60)
            if r.status_code == 200:
                lines = r.text.strip().splitlines()
                hdr = lines[0].split("\t")
                out = []
                for line in lines[1:]:
                    parts = line.split("\t")
                    row = dict(zip(hdr, parts))
                    if len(parts) == 2:
                        row[hdr[2]] = ""
                    out.append(row)
                return out
        except requests.RequestException:
            pass
        time.sleep(2**attempt)
    raise RuntimeError(f"uniprot batch failed: {accs[:3]}...")


# Known approved pairs (ChEMBL should contain every one of these).
# Drug names are matched by prefix: ChEMBL records salt forms as separate
# molecules (IMATINIB MESYLATE, WARFARIN SODIUM, ...). Gene symbols use the
# current HGNC names (PTH1R, not PTHR1).
SPOT_CHECKS = [
    ("METHOTREXATE", "DHFR", "inhibitor"),
    ("IMATINIB", "ABL1", "inhibitor"),
    ("SIMVASTATIN", "HMGCR", "inhibitor"),
    ("WARFARIN", "VKORC1", "inhibitor"),
    ("ABALOPARATIDE", "PTH1R", "agonist"),
]


def main() -> None:
    j = pd.read_parquet(RESULTS / "chembl_fda_targets.parquet")
    # normalize the join's suffix-collided drug name column
    if "pref_name.mol" in j.columns:
        j = j.rename(columns={"pref_name.mol": "pref_name"})
    accs = sorted(j["component_accession"].dropna().unique())
    print(f"accessions to map: {len(accs)}", flush=True)

    rows = []
    for i in range(0, len(accs), 25):
        rows.extend(fetch_batch(accs[i : i + 25]))
        if i and i % 250 == 0:
            print(f"  mapped {i}/{len(accs)}", flush=True)
    m = pd.DataFrame(rows)
    m = m.rename(
        columns={
            "Entry": "component_accession",
            "Gene Names  (primary )": "gene_symbol",
            "Gene Names (primary)": "gene_symbol",
            "Protein names": "protein_name",
        }
    )
    cols = [c for c in ["component_accession", "gene_symbol", "protein_name"] if c in m.columns]
    m = m[cols]
    # normalize: UniProt tsv may multi-assign; keep the protein's primary gene
    m["gene_symbol"] = m["gene_symbol"].fillna("").str.upper()
    m.to_parquet(RESULTS / "accession_gene_map.parquet", index=False)
    print(m.head(15).to_string(), flush=True)

    jm = j.merge(m, on="component_accession", how="left")
    dt = (
        jm.assign(drug_name=jm["pref_name"].str.upper())
        .loc[:, ["molecule_chembl_id", "pref_name", "first_approval", "withdrawn_flag", "gene_symbol"]]
        .loc[:, ["molecule_chembl_id", "pref_name", "first_approval", "withdrawn_flag", "gene_symbol"]]
        .dropna(subset=["gene_symbol"])
        .loc[lambda d: d["gene_symbol"] != ""]
        .drop_duplicates()
        .sort_values(["pref_name", "gene_symbol"])
        .reset_index(drop=True)
    )
    dt.to_parquet(RESULTS / "drug_target_genes.parquet", index=False)
    print(
        f"drug-target-gene table: {len(dt)} rows, "
        f"{dt.molecule_chembl_id.nunique()} drugs, {dt.gene_symbol.nunique()} genes",
        flush=True,
    )

    # ---- verification ------------------------------------------------------
    joined = jm.assign(drug=jm["pref_name"].str.upper(), gene=jm["gene_symbol"])
    verdicts = {}
    for drug, gene, action in SPOT_CHECKS:
        sel = joined[joined["drug"].str.startswith(drug)]
        hit = bool((sel["gene"] == gene).any())
        acts = sorted(sel.loc[sel["gene"] == gene, "action_type"].dropna().astype(str).str.lower().unique())
        verdicts[drug] = {
            "target_present": hit,
            "action_type": acts,
            "pass": hit and any(action in a for a in acts),
        }
    print(json.dumps(verdicts, indent=1), flush=True)
    (RESULTS / "spot_checks.json").write_text(json.dumps(verdicts, indent=1))
    bad = [k for k, v in verdicts.items() if not v["pass"]]
    if bad:
        print(f"SPOT-CHECK FAILURES: {bad}", flush=True)
    else:
        print("all spot checks pass", flush=True)


if __name__ == "__main__":
    main()
