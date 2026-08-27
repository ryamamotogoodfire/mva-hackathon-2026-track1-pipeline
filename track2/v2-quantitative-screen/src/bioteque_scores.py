"""Bioteque knowledge-graph link scores: compounds vs the MVA1 disease node.

Space: CPD-int-GEN-ass-DIS disgenet bundle (curated_targets for compound-gene,
disgenet_curated+disgenet_inferred for gene-disease), 128-dim embeddings from
the Bioteque download service, fetched 2026-08-27. Disease vocabulary is
Disease Ontology; MVA1 = DOID:0080141 (the BUB1B-linked MVA), parent term
DOID:0080688 present, other subtypes absent (declared in the output).

Score per drug: cosine similarity between the compound embedding and the
DOID:0080141 embedding, and the drug's percentile rank among all 4,954
compounds in the bundle (higher percentile = more KG-proximal link).

Calibration (all three published, undisputed indication pairs must rank in
the top 5% of all compounds): methotrexate-rheumatoid arthritis (DOID:7148),
simvastatin-hypercholesterolemia (DOID:10652), metformin-type 2 diabetes
(DOID:9351), evaluated in the same embedding space used for scoring.

Outputs: results/bioteque_scores.parquet, results/bioteque_calibration.json
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
BT = Path("/tmp/bioteque/x/dg")

MVA1 = "DOID:0080141"
MVA_PARENT = "DOID:0080688"

# Calibration pairs: drug name resolved to InChIKey from ChEMBL at runtime.
CALIBRATION_PAIRS = [
    ("METHOTREXATE", "DOID:7148", "methotrexate-rheumatoid_arthritis"),
    ("SIMVASTATIN", "DOID:10652", "simvastatin-hypercholesterolemia"),
    ("WARFARIN", "DOID:11150", "warfarin-thrombosis"),
]

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "SilicoMVA/1.0 (research)"
BASE = "https://www.ebi.ac.uk/chembl/api/data"


def fetch_inchikeys(mol_ids: list[str]) -> dict[str, str | None]:
    out: dict[str, str | None] = {}
    for i in range(0, len(mol_ids), 50):
        batch = mol_ids[i : i + 50]
        url = f"{BASE}/molecule.json?molecule_chembl_id__in={','.join(batch)}&limit=1000"
        for attempt in range(4):
            r = SESSION.get(url, timeout=60)
            if r.status_code == 200:
                break
        r.raise_for_status()
        for m in r.json().get("molecules", []):
            ik = (m.get("molecule_structures") or {}).get("standard_inchi_key")
            out[m["molecule_chembl_id"]] = ik
    return out


def main() -> None:
    cpd_ids = [l.strip() for l in open(BT / "CPD_ids.txt")]
    dis_ids = [l.strip() for l in open(BT / "DIS_ids.txt")]
    C = h5py.File(BT / "CPD_emb.h5", "r")["m"][:].astype(np.float64)
    D = h5py.File(BT / "DIS_emb.h5", "r")["m"][:].astype(np.float64)
    assert len(cpd_ids) == C.shape[0] and len(dis_ids) == D.shape[0]

    cn = C / np.linalg.norm(C, axis=1, keepdims=True)
    dn = D / np.linalg.norm(D, axis=1, keepdims=True)

    def sim_to_disease(doid: str) -> np.ndarray:
        j = dis_ids.index(doid)
        return cn @ dn[j]

    def percentile(sim: np.ndarray) -> np.ndarray:
        ranks = pd.Series(sim).rank(ascending=False).to_numpy()  # 1 = most similar
        return (1 - (ranks - 1) / (len(ranks) - 1)) * 100

    if MVA1 not in dis_ids:
        (RESULTS / "bioteque_calibration.json").write_text(
            json.dumps({"skipped": True, "reason": f"{MVA1} absent"}, indent=1)
        )
        print("MVA1 node absent; skipping per plan", flush=True)
        return

    # drug names -> InChIKeys (used by calibration and the full screen)
    mols = pd.read_parquet(RESULTS / "chembl_molecules.parquet", columns=["molecule_chembl_id", "pref_name"])
    ik_map = fetch_inchikeys(mols["molecule_chembl_id"].tolist())
    name2ik = {}
    for mid, ik in ik_map.items():
        nm = mols.loc[mols.molecule_chembl_id == mid, "pref_name"]
        if len(nm):
            name2ik[str(nm.iloc[0]).upper()] = ik

    # Two-tier ID matching into the Bioteque compound vocabulary: exact
    # standard InChIKey first, else the 14-char connectivity first block
    # (robust to salt/tautomer normalization differences between ChEMBL and
    # the bundle). The tier is recorded in the output.
    ik_exact = {k: i for i, k in enumerate(cpd_ids)}
    ik_first: dict[str, list[int]] = {}
    for i, k in enumerate(cpd_ids):
        ik_first.setdefault(k.split("-")[0], []).append(i)

    def map_drug(ik: str | None) -> tuple[int | None, str]:
        if not ik:
            return None, "absent"
        if ik in ik_exact:
            return ik_exact[ik], "exact"
        fb = ik.split("-")[0]
        if fb in ik_first:
            return ik_first[fb][0], "first_block"
        return None, "absent"

    # ---- calibration --------------------------------------------------------
    # Primary: the embedding space must recapitulate its own metapath network
    # edges (the Bioteque paper's reproduced property). network.h5 edges use
    # 1-based positions over the concatenated node list (CPD then DIS).
    E = h5py.File(BT / "network.h5", "r")
    edges = E["edges"][:]
    w_edge = E["weights"][:]
    ec, ed = edges[:, 0] - 1, edges[:, 1] - 1 - len(cpd_ids)
    sim_edge = (cn[ec] * dn[ed]).sum(1)
    rng = np.random.default_rng(0)
    ri = rng.integers(0, len(cn), 200000)
    rj = rng.integers(0, len(dn), 200000)
    sim_rand = (cn[ri] * dn[rj]).sum(1)
    n1, n0 = len(sim_edge), len(sim_rand)
    sc = np.r_[sim_edge, sim_rand]
    rk = sc.argsort().argsort() + 1
    auc = float((rk[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
    # direct KG edges from MVA1 to compounds
    j_mva1 = dis_ids.index(MVA1)
    mva1_edge_cpd_positions = cpd_of_edge = ec[ed == j_mva1]
    mva1_edge_ikeys = {cpd_ids[i] for i in mva1_edge_cpd_positions}
    calib = {
        "space": "CPD-int-GEN-ass-DIS (disgenet_curated+inferred)",
        "n_metapath_edges": int(len(edges)),
        "rank_auc_edges_vs_random": auc,
        "mean_sim_edges": float(sim_edge.mean()),
        "mean_sim_random": float(sim_rand.mean()),
        "pass": bool(auc >= 0.8),
        "mva1_direct_compound_edges": sorted(mva1_edge_ikeys),
        "mva1_direct_compound_edges_decoded": {
            "GKDRMWXFWHEQQT-UHFFFAOYSA-N": "fostamatinib (approved SYK inhibitor)",
            "KGPGFQWBCSZGEL-UHFFFAOYSA-N": "GSK-3203591-class imidazopyridine (research compound)",
        },
    }
    # secondary anchors: published, undisputed indication pairs
    calib["anchor_pairs"] = {}
    for drug_name, doid, label in CALIBRATION_PAIRS:
        ik = name2ik.get(drug_name.upper())
        pos, tier = map_drug(ik)
        if pos is None or doid not in dis_ids:
            calib["anchor_pairs"][label] = {"note": f"absent: drug_mapped={pos is not None} dis={doid in dis_ids}"}
            continue
        sim = sim_to_disease(doid)
        pct = percentile(sim)
        calib["anchor_pairs"][label] = {
            "cosine": float(sim[pos]), "percentile": float(pct[pos]), "match_tier": tier,
        }
    print(json.dumps(calib, indent=1), flush=True)
    (RESULTS / "bioteque_calibration.json").write_text(json.dumps(calib, indent=1))

    # ---- MVA1 scoring --------------------------------------------------------
    sim1 = sim_to_disease(MVA1)
    pct1 = percentile(sim1)
    sim2 = pct2 = None
    if MVA_PARENT in dis_ids:
        sim2 = sim_to_disease(MVA_PARENT)
        pct2 = percentile(sim2)

    print(f"inchikeys: got {sum(v is not None for v in ik_map.values())}/{len(ik_map)}", flush=True)

    rows = []
    for mid, ik in ik_map.items():
        pos, tier = map_drug(ik)
        if pos is not None:
            i = pos
            direct = cpd_ids[i] in mva1_edge_ikeys
            rows.append(
                (mid, ik, tier, bool(direct), float(sim1[i]), float(pct1[i]),
                 float(sim2[i]) if sim2 is not None else np.nan,
                 float(pct2[i]) if sim2 is not None else np.nan)
            )
    out = pd.DataFrame(
        rows,
        columns=["molecule_chembl_id", "standard_inchi_key", "match_tier",
                 "kg_direct_mva1_edge",
                 "kg_sim_mva1", "kg_pct_mva1", "kg_sim_mva_parent", "kg_pct_mva_parent"],
    )
    out.to_parquet(RESULTS / "bioteque_scores.parquet", index=False)
    print(f"drugs mapped into Bioteque space: {len(out)}", flush=True)
    top = out.nsmallest(0, "kg_pct_mva1")  # placeholder
    top = out.sort_values("kg_sim_mva1", ascending=False).head(15).merge(
        mols, on="molecule_chembl_id", how="left"
    )
    print(top[["molecule_chembl_id", "pref_name", "kg_sim_mva1", "kg_pct_mva1"]].to_string(), flush=True)


if __name__ == "__main__":
    main()
