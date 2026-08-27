"""Join the screen's three lanes into one per-drug table for synthesis.

Lanes (one row per ChEMBL molecule):
  - prox_z: Barabasi/Guney network-proximity z-score vs the MVA module
    (lower = nearer).
  - kg_pct_mva1 / kg_sim_mva1 / kg_direct_mva1_edge: Bioteque disgenet-inferred
    embedding cosine percentile to DOID:0080141 (higher = nearer), with direct
    KG edges to MVA1 flagged.
  - l1000_n_sig, l1000_tau_left_mass(frac of signatures with tau < -30),
    l1000_tau_min: Phase I BUB1B-low reversal screen (instrument-limited).
  - openfda_found / has_pediatric_section.

IDs joined on molecule_chembl_id; L1000 joined via two-tier InChIKey mapping.

Output: results/screen_joined.parquet + results/screen_joined_summary.json
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"


def main() -> None:
    dt = pd.read_parquet(RESULTS / "drug_target_genes.parquet")
    drugs = (
        dt[["molecule_chembl_id", "pref_name", "first_approval", "withdrawn_flag"]]
        .drop_duplicates("molecule_chembl_id")
    )
    targets = dt.groupby("molecule_chembl_id").gene_symbol.agg(sorted).rename("targets").reset_index()
    drugs = drugs.merge(targets, on="molecule_chembl_id", how="left")

    prox = pd.read_parquet(RESULTS / "proximity_scores.parquet")[
        ["molecule_chembl_id", "prox_obs", "prox_z"]
    ]
    drugs = drugs.merge(prox, on="molecule_chembl_id", how="left")

    kg = pd.read_parquet(RESULTS / "bioteque_scores.parquet")
    drugs = drugs.merge(
        kg[["molecule_chembl_id", "kg_sim_mva1", "kg_pct_mva1", "kg_sim_mva_parent", "kg_pct_mva_parent", "match_tier", "kg_direct_mva1_edge"]],
        on="molecule_chembl_id", how="left",
    )

    pa = pd.read_parquet(RESULTS / "openfda_pediatric.parquet")
    drugs = drugs.merge(
        pa[["pref_name", "openfda_found", "has_pediatric_section", "pediatric_text_len"]],
        on="pref_name", how="left",
    )

    # ---- L1000 per-drug collapse, mapped by InChIKey ------------------------
    tau = pd.read_parquet(RESULTS / "l1000/reversal_tau.parquet")
    ik = pd.read_parquet(RESULTS / "molecule_inchikeys.parquet")
    ik_exact = ik.set_index("standard_inchi_key").molecule_chembl_id.dropna().to_dict()
    ik_first = ik.dropna(subset=["inchi_key_prefix"]).groupby("inchi_key_prefix").molecule_chembl_id.first().to_dict()

    # name tiers (case-insensitive exact, then salt-stripped stem)
    name_exact = ik.set_index(ik.pref_name.str.upper()).molecule_chembl_id.to_dict()
    def stem(name: str) -> str:
        toks = name.upper().split()
        cut = {"HYDROCHLORIDE", "DIHYDROCHLORIDE", "SODIUM", "POTASSIUM", "SULFATE", "MESYLATE", "TOSYLATE", "FOs", "FUMARATE", "ANHYDROUS", "MONOHYDRATE", "DIPOTASSIUM", "PAMOATE"}
        toks = [w for w in toks if w not in cut]
        return " ".join(toks) if toks else name.upper()
    name_stem = {}
    for nm, mid in ik.set_index(ik.pref_name.str.upper()).molecule_chembl_id.items():
        name_stem.setdefault(stem(nm), mid)

    def map_pert(row) -> tuple[str | None, str]:
        key = row["inchi_key"]
        pref = row["inchi_key_prefix"]
        if isinstance(key, str) and key in ik_exact:
            return ik_exact[key], "exact"
        if isinstance(pref, str) and pref in ik_first:
            return ik_first[pref], "first_block"
        nm = str(row["pert_iname"]).upper()
        if nm in name_exact:
            return name_exact[nm], "name"
        st = stem(nm)
        if st in name_stem:
            return name_stem[st], "name_stem"
        return None, "absent"

    mp = tau[["pert_id", "pert_iname", "inchi_key", "inchi_key_prefix"]].drop_duplicates("pert_id").copy()
    mp["pair"] = mp.apply(map_pert, axis=1)
    mp["molecule_chembl_id"] = mp.pair.str[0]
    mp["map_tier"] = mp.pair.str[1]
    tau2 = tau.merge(mp[["pert_id", "molecule_chembl_id", "map_tier"]], on="pert_id", how="left")

    strong_rev = tau2[tau2.tau < -30]
    per_drug = (
        tau2.dropna(subset=["molecule_chembl_id"])
        .groupby("molecule_chembl_id")
        .agg(
            l1000_n_sig=("tau", "size"),
            l1000_pert_ids=("pert_id", "nunique"),
            l1000_tau_min=("tau", "min"),
            l1000_tau_med=("tau", "median"),
        )
        .reset_index()
    )
    sr = (
        strong_rev.dropna(subset=["molecule_chembl_id"])
        .groupby("molecule_chembl_id")
        .size()
        .rename("l1000_n_strong_rev")
        .reset_index()
    )
    per_drug = per_drug.merge(sr, on="molecule_chembl_id", how="left")
    per_drug["l1000_n_strong_rev"] = per_drug.l1000_n_strong_rev.fillna(0).astype(int)
    per_drug["l1000_frac_strong_rev"] = per_drug.l1000_n_strong_rev / per_drug.l1000_n_sig
    # background fraction for calibration signal
    bg = len(tau2[tau2.tau < -30]) / len(tau2)
    drugs = drugs.merge(per_drug, on="molecule_chembl_id", how="left")

    drugs.to_parquet(RESULTS / "screen_joined.parquet", index=False)

    summ = {
        "n_drugs": int(len(drugs)),
        "n_with_prox_z": int(drugs.prox_z.notna().sum()),
        "n_with_kg": int(drugs.kg_sim_mva1.notna().sum()),
        "n_with_l1000": int(drugs.l1000_n_sig.notna().sum()),
        "l1000_background_frac_tau_lt_-30": float(bg),
        "n_openfda_found": int(drugs.openfda_found.fillna(False).sum()),
        "n_pediatric_section": int(drugs.has_pediatric_section.fillna(False).sum()),
    }
    (RESULTS / "screen_joined_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1), flush=True)


if __name__ == "__main__":
    main()
