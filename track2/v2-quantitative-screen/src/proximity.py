"""Barabasi/Guney-style network proximity z-scores, drug targets vs a disease module.

Method (Guney et al., Nat Commun 2016;7:10331):
  p(S,T) = (1/|T|) * sum over t in T of min over s in S of shortest path d(s,t)
  z      = (p - mean(null)) / sd(null), null = 1000 random target sets of the
           same size and degree distribution.

Disease module (declared before scoring):
  MVA founder genes: BUB1B, BUB1, CEP57, TRIP13
  core spindle-checkpoint: MAD1L1, MAD2L1, BUB3, TTK, CDC20, CENPE, AURKB, PLK1

Sanity checks (run before any new score is trusted):
  1. ivacaftor against the cystic fibrosis module {CFTR}: mechanism-direct,
     published-proximal case; expect distance 0, strongly negative z.
  2. imatinib against {BCR, ABL1}: same class of mechanism-direct check.
  3. simvastatin against the MVA module: off-axis drug, expect |z| small.

Outputs:
  results/proximity_scores.parquet
  results/proximity_sanity.json
"""

from __future__ import annotations

import json
from collections import deque
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"

MVA_MODULE = [
    "BUB1B", "BUB1", "CEP57", "TRIP13",
    "MAD1L1", "MAD2L1", "BUB3", "TTK", "CDC20", "CENPE", "AURKB", "PLK1",
]
N_NULL = 1000
RNG = np.random.default_rng(42)


def distances_to_module(G: nx.Graph, module: list[str]) -> dict[str, float]:
    """Multi-source BFS: distance from the module to every node."""
    dist: dict[str, float] = {}
    q = deque()
    for s in module:
        if s in G:
            dist[s] = 0.0
            q.append(s)
    while q:
        u = q.popleft()
        for v in G.neighbors(u):
            if v not in dist:
                dist[v] = dist[u] + 1.0
                q.append(v)
    return dist


def degree_pools(nodes_sorted_idx: np.ndarray, degrees: np.ndarray):
    """For each distinct target degree, sample pool = exact-degree nodes; when
    fewer than 50 nodes share the degree, widen to a +-100-rank window on the
    degree-sorted order, following the Guney toolbox's matched-sampling rule."""
    dord = np.sort(degrees)
    pools = {}
    for d in np.unique(degrees):
        idx = np.where(degrees == d)[0]
        if len(idx) < 50:
            r_lo = max(np.searchsorted(dord, d, side="left") - 100, 0)
            r_hi = min(np.searchsorted(dord, d, side="right") + 100, len(dord) - 1)
            idx = np.where((degrees >= dord[r_lo]) & (degrees <= dord[r_hi]))[0]
        pools[int(d)] = nodes_sorted_idx[idx]
    return pools


def z_for_drug(targets: list[str], deg_of: dict[str, int],
               pools: dict[float, np.ndarray], dist_arr: np.ndarray, node_pos: dict[str, int],
               n_null: int = N_NULL, rng: np.random.Generator = RNG) -> tuple[float, float, float, float]:
    """Return (p_obs, z, mu_null, sd_null). Distance mean over targets vs null sets."""
    present = [t for t in targets if t in deg_of and t in node_pos]
    if not present:
        return float("nan"), float("nan"), float("nan"), float("nan")
    pos = np.array([node_pos[t] for t in present])
    obs = float(dist_arr[pos].mean())
    null_means = np.empty(n_null)
    for i in range(n_null):
        vals = []
        for t in present:
            pool = pools[deg_of[t]]
            vals.append(dist_arr[pool[rng.integers(len(pool))]])
        null_means[i] = np.mean(vals)
    mu, sd = float(null_means.mean()), float(null_means.std(ddof=1))
    z = (obs - mu) / sd if sd > 0 else float("nan")
    return obs, z, mu, sd


def main() -> None:
    edges = pd.read_parquet(RESULTS / "biogrid_edges.parquet")
    G = nx.Graph()
    G.add_edges_from(edges[["a", "b"]].itertuples(index=False))
    n = G.number_of_nodes()
    node_pos = {node: i for i, node in enumerate(G.nodes())}
    deg_of = dict(G.degree())
    degrees = np.array([deg_of[nod] for nod in node_pos], dtype=int)

    # ---- sanity-check drugs and modules -----------------------------------
    all_dist_cache = {}

    def dist_to(module: list[str]) -> np.ndarray:
        key = tuple(module)
        if key not in all_dist_cache:
            d = distances_to_module(G, list(module))
            # nodes unreachable -> network diameter + 1
            arr = np.full(n, 10.0)
            for nod, dval in d.items():
                arr[node_pos[nod]] = dval
            all_dist_cache[key] = arr
        return all_dist_cache[key]

    pools = degree_pools(np.arange(n), degrees)

    def score(drug_targets: list[str], module: list[str]):
        arr = dist_to(module)
        return z_for_drug(drug_targets, deg_of, pools, arr, node_pos)

    sanity = {}
    # ivacaftor targets CFTR only
    obs, z, mu, sd = score(["CFTR"], ["CFTR"])
    sanity["ivacaftor_vs_CFCmodule(CFTR-only)"] = {"p_obs": obs, "z": z, "mu_null": mu, "sd_null": sd, "pass": obs == 0 and z < -1.5}
    # imatinib targets ABL1,BCR,KIT,PDGFRB vs CML module {BCR,ABL1}: two targets
    # sit inside the module, so the published-proximal expectation is z << 0
    obs, z, mu, sd = score(["ABL1", "BCR", "KIT", "PDGFRB"], ["BCR", "ABL1"])
    sanity["imatinib_vs_cml_module(BCR-ABL1)"] = {"p_obs": obs, "z": z, "mu_null": mu, "sd_null": sd, "pass": z < -1.5}
    # simvastatin vs MVA module: expect not extreme
    obs, z, mu, sd = score(["HMGCR"], MVA_MODULE)
    sanity["simvastatin_vs_MVA_module"] = {"p_obs": obs, "z": z, "mu_null": mu, "sd_null": sd, "pass": abs(z) < 1.5}
    print(json.dumps(sanity, indent=1), flush=True)
    (RESULTS / "proximity_sanity.json").write_text(json.dumps(sanity, indent=1))

    # ---- full screen: MVA module ------------------------------------------
    missing = [g for g in MVA_MODULE if g not in G]
    print("module genes missing from network:", missing, flush=True)
    arr = dist_to(MVA_MODULE)
    dt = pd.read_parquet(RESULTS / "drug_target_genes.parquet")
    drug_targets = dt.groupby(["molecule_chembl_id", "pref_name"])["gene_symbol"].agg(sorted).reset_index()
    drug_targets["target_list"] = drug_targets["gene_symbol"]
    rows = []
    for r in drug_targets.itertuples(index=False):
        targets = [t for t in r.target_list if t in deg_of]
        if not targets:
            rows.append((r.molecule_chembl_id, r.pref_name, len(r.target_list), 0, None, None, None, None))
            continue
        obs, z, mu, sd = z_for_drug(targets, deg_of, pools, arr, node_pos)
        rows.append((r.molecule_chembl_id, r.pref_name, len(r.target_list), len(targets), obs, z, mu, sd))
    out = pd.DataFrame(rows, columns=[
        "molecule_chembl_id", "pref_name", "n_targets_named", "n_targets_in_network",
        "prox_obs", "prox_z", "prox_mu_null", "prox_sd_null",
    ])
    out.to_parquet(RESULTS / "proximity_scores.parquet", index=False)
    print(out.sort_values("prox_z").head(25).to_string(), flush=True)
    print("drugs scored:", int(out.prox_z.notna().sum()), flush=True)


if __name__ == "__main__":
    main()
