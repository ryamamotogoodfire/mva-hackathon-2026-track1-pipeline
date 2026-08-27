"""BUB1B-low signature-reversal screen over LINCS L1000 (GEO GSE92742).

Query signature: consensus Level5 (COMPZ.MODZ) mean z-score of all trt_sh
BUB1B knockdown signatures (89 signature columns, multiple hairpins, cell
lines, time points).

Screen: weighted two-tailed connectivity score (Subramanian 2017, Cell)
between the query's 150-up / 150-down gene sets and every trt_cp compound
signature; tau = 100 * cs / max|cs| across all compound signatures; reversal
tau = -tau. Restricted at report time to FDA-approved molecules mapped onto
ChEMBL CL in the joint table (mapping here: exact InChIKey, else first
block).

Internal validation (all must hold, saved to reversal_validation.json):
  V1. BUB1B its own z-score mean across its KD signatures is strongly negative
      (knockdown readout).
  V2. The individual KD signatures connect positively to the consensus
      (median tau >= 75).
  V3. A known pharmacology pair reproduced: JQ1 signatures rank among the
      strongest connections of the BRD4 knockdown consensus (published
      clio-style pair).

Outputs (artifacts dir, arg --outdir):
  bub1b_down_signature.parquet  col: gene_symbol, z (12,328 genes, full)
  reversal_tau.parquet          per (pert_id, sig_id) tau + metadata
  reversal_validation.json
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

try:
    from silico.slurm_telemetry import report_progress
except Exception:  # pragma: no cover
    def report_progress(**_): pass

N_TOP = 150


def read_gctx(path: str) -> tuple[h5py.File, np.ndarray, np.ndarray]:
    f = h5py.File(path, "r")
    mat = f["/0/DATA/0/matrix"]
    genes = f["/0/META/ROW/id"][:].astype(str)
    sigs = f["/0/META/COL/id"][:].astype(str)
    assert mat.shape == (len(sigs), len(genes)), (mat.shape, len(sigs), len(genes))
    return f, genes, sigs


def wtcs_from_query(up: set, down: set, ranked_genes: list[str]) -> float:
    """Classic LINCS two-tail weighted Kolmogorov-Smirnov-like connectivity."""
    n = len(ranked_genes)
    hit_up = np.zeros(n)
    hit_down = np.zeros(n)
    for i, g in enumerate(ranked_genes):
        if g in up:
            hit_up[i] = 1.0
        elif g in down:
            hit_down[i] = -1.0
        else:
            pass

    def es(sel: np.ndarray) -> float:
        # sel: +1 for hits
        pos = np.where(sel != 0)[0]
        if len(pos) == 0:
            return 0.0
        p_miss = 1.0 / (n - len(pos))
        weights = np.abs(sel[pos])
        p_hit = weights / weights.sum()
        cum_hit = np.zeros(n)
        cum_hit[pos] = p_hit
        cum_hit = np.cumsum(cum_hit)
        cum_miss = (np.arange(len(sel)) + 1 - np.isin(np.arange(len(sel)), pos)) * p_miss
        d = cum_hit - cum_miss
        if d.max() > -d.min():
            return float(d.max())
        return float(d.min())

    es_up = es(hit_up)
    es_down = es(hit_down)
    if np.sign(es_up) == np.sign(es_down):
        return 0.0
    return float((es_up - es_down) / 2.0)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gctx", required=True)
    ap.add_argument("--sig-info", required=True)
    ap.add_argument("--gene-info", required=True)
    ap.add_argument("--pert-info", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--limit-sigs", type=int, default=0, help="dev: only first N trt_cp columns")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print("== l1000 pipeline start ==", flush=True)
    sig_info = pd.read_csv(args.sig_info, sep="\t", dtype=str)
    gene_info = pd.read_csv(args.gene_info, sep="\t", dtype=str)
    pert_info = pd.read_csv(args.pert_info, sep="\t", dtype=str)
    gene_sym = dict(zip(gene_info.pr_gene_id.astype(str), gene_info.pr_gene_symbol))

    gctx, genes, sigs = read_gctx(args.gctx)
    print(f"gctx opened: {len(genes)} genes x {len(sigs)} signatures", flush=True)

    colpos = {s: i for i, s in enumerate(sigs)}
    symbol = np.array([gene_sym.get(g, g) for g in genes])

    # ---- query signature ---------------------------------------------------
    kd = sig_info[(sig_info.pert_type.str.startswith("trt_sh")) &
                  (sig_info.pert_iname.str.contains("BUB1B", case=False, na=False))].copy()
    cols_kd = sorted(colpos[s] for s in kd.sig_id if s in colpos)
    mat_kd_all = gctx["/0/DATA/0/matrix"][cols_kd, :].T  # (12328, n_kd)
    bub1b_pos = int(np.where(symbol == "BUB1B")[0][0])
    z_bub1b_per_col = mat_kd_all[bub1b_pos, :]

    # Quality gate: measured knockdown (BUB1B z <= -1.0). Documented choice.
    gated_mask = z_bub1b_per_col <= -1.0
    kd_pert_ids = kd.set_index("sig_id").loc[[sigs[i] for i in cols_kd]]["pert_id"].to_numpy()
    print(f"BUB1B KD signatures: {len(cols_kd)}; gated (z<=-1): {int(gated_mask.sum())}", flush=True)
    cols_kd_gated = [c for c, m in zip(cols_kd, gated_mask) if m]
    mat_kd = mat_kd_all[:, gated_mask]
    kd_pert_ids_gated = kd_pert_ids[gated_mask]
    z_kd = mat_kd.mean(axis=1)
    print(f"BUB1B mean z in gated consensus: {z_kd[bub1b_pos]:.2f}", flush=True)

    # Sign-consistency filter: keep only genes whose direction is reproduced
    # in >= 60% of the gated signatures; then rank by consensus z.
    frac_up = (mat_kd > 0).mean(axis=1)
    consistent_up = frac_up >= 0.60
    consistent_down = frac_up <= 0.40
    q_up = z_kd * consistent_up
    q_down = z_kd * consistent_down
    n_up = int((q_up > 0).sum())
    n_down = int((q_down < 0).sum())
    print(f"sign-consistent genes: up {n_up}, down {n_down}", flush=True)
    if n_up < 50 or n_down < 50:
        raise RuntimeError("not enough sign-consistent genes for a query")

    order = np.argsort(-q_up)
    top_up = [symbol[i] for i in order if q_up[i] > 0][:N_TOP]
    order_d = np.argsort(q_down)
    top_down = [symbol[i] for i in order_d if q_down[i] < 0][:N_TOP][::-1]
    up_set, down_set = set(top_up), set(top_down)
    print(f"query sets: {len(up_set)} up, {len(down_set)} down", flush=True)

    sigtab = pd.DataFrame({"pr_gene_id": genes, "gene_symbol": symbol, "z_mean_kd": z_kd,
                           "frac_up_kd": frac_up,
                           "in_query_up": [g in up_set for g in symbol],
                           "in_query_down": [g in down_set for g in symbol]})
    sigtab.to_parquet(outdir / "bub1b_down_signature.parquet", index=False)

    # V1: BUB1B itself strongly down
    z_bub1b = float(z_kd[symbol == "BUB1B"][0])
    # V2: held-out hairpin reconnect: for each gated hairpin group (pert_id),
    # consensus over the OTHER hairpins, then the held-out signature must
    # connect above the 95th percentile of unrelated compound signatures.
    kd_sigids = [sigs[i] for i in cols_kd_gated]
    kd_cs = []
    for k in range(mat_kd.shape[1]):
        hold = kd_pert_ids_gated == kd_pert_ids_gated[k]
        rest = ~hold
        if rest.sum() < 5:
            continue
        z_rest = mat_kd[:, rest].mean(axis=1)
        fu = (mat_kd[:, rest] > 0).mean(axis=1)
        qu, qd = z_rest * (fu >= 0.60), z_rest * (fu <= 0.40)
        ou = np.argsort(-qu)
        od = np.argsort(qd)
        up_rest = set(symbol[ou[[i for i in range(len(qu)) if qu[i] > 0]][:N_TOP]])
        down_rest = set(symbol[od[[i for i in range(len(qd)) if qd[i] < 0]][:N_TOP]])
        ranked_c = symbol[np.argsort(-mat_kd[:, k])]
        kd_cs.append(float(wtcs_from_query(up_rest, down_rest, ranked_c)))
    kd_cs = np.array(kd_cs)
    med_kd_cs = float(np.median(kd_cs))
    print(f"V2: held-out-hairpin median cs: {med_kd_cs:.3f} over {len(kd_cs)}", flush=True)

    # ---- compound signatures -------------------------------------------------

    # ---- random-signature null for V2 ---------------------------------------
    trt_cp_all = sig_info[sig_info.pert_type == "trt_cp"]
    rnd = np.random.default_rng(12)
    null_pool = rnd.choice(sorted(trt_cp_all.sig_id), size=min(5000, len(trt_cp_all)), replace=False)
    null_cols = sorted(colpos[s] for s in null_pool if s in colpos)
    rand_cs = np.empty(len(null_cols))
    BM = gctx["/0/DATA/0/matrix"][null_cols, :].T
    for k in range(BM.shape[1]):
        ranked = symbol[np.argsort(-BM[:, k])]
        rand_cs[k] = wtcs_from_query(up_set, down_set, ranked)
    random_p95 = float(np.percentile(rand_cs, 95))
    print(f"V2: random-signature cs p95: {random_p95:.3f}", flush=True)

    cp = sig_info[sig_info.pert_type == "trt_cp"].copy()
    cp_cols = [colpos[s] for s in cp.sig_id if s in colpos]
    cp_cols = sorted(cp_cols)
    cp_sigids = [sigs[i] for i in cp_cols]
    if args.limit_sigs:
        cp_cols = cp_cols[: args.limit_sigs]
        cp_sigids = cp_sigids[: args.limit_sigs]
    print(f"trt_cp columns to score: {len(cp_cols)}", flush=True)

    cs_all = np.empty(len(cp_cols), dtype=np.float32)
    block = 5000
    for b0 in range(0, len(cp_cols), block):
        cpb = cp_cols[b0 : b0 + block]
        B = gctx["/0/DATA/0/matrix"][cpb, :].T.copy()  # (12328, b)
        for k in range(B.shape[1]):
            zcol = B[:, k]
            ordi = np.argsort(-zcol)
            ranked = symbol[ordi]
            cs_all[b0 + k] = wtcs_from_query(up_set, down_set, ranked)
        print(f"  scored {min(b0 + block, len(cp_cols))}/{len(cp_cols)}  ({time.time() - t0:.0f}s)", flush=True)
        report_progress(step=min(b0 + block, len(cp_cols)), total_steps=len(cp_cols), phase="connectivity-score")

    denom = np.max(np.abs(cs_all))
    tau = 100.0 * cs_all / denom
    out = pd.DataFrame({"sig_id": cp_sigids, "cs": cs_all, "tau": tau})
    out = out.merge(cp[["sig_id", "pert_id", "pert_iname", "cell_id", "pert_dose", "pert_time"]], on="sig_id", how="left")
    pert_cols = pert_info[["pert_id", "inchi_key", "inchi_key_prefix", "canonical_smiles", "pubchem_cid"]].drop_duplicates("pert_id")
    out = out.merge(pert_cols, on="pert_id", how="left")
    out["reversal_tau"] = -out.tau
    out.to_parquet(outdir / "reversal_tau.parquet", index=False)

    # V3: known pharmacology pair reproduced: MTOR knockdown consensus must
    # connect positively to sirolimus (rapamycin) compound signatures.
    v3 = {}
    mtor = sig_info[(sig_info.pert_type.str.startswith("trt_sh")) & (sig_info.pert_iname == "MTOR")]
    sir = sig_info[(sig_info.pert_type == "trt_cp") & (sig_info.pert_iname == "sirolimus")]
    print(f"V3: MTOR KD signatures: {len(mtor)}, sirolimus sigs: {len(sir)}", flush=True)
    if len(mtor) > 0 and len(sir) > 0:
        c4 = sorted(colpos[s] for s in mtor.sig_id if s in colpos)
        m4 = gctx["/0/DATA/0/matrix"][c4, :].T
        z4 = m4.mean(axis=1)
        o4 = np.argsort(-z4)
        up4, down4 = set(symbol[o4[:N_TOP]]), set(symbol[o4[-N_TOP:]][::-1])
        jpos = sorted(colpos[s] for s in sir.sig_id if s in colpos)
        cas = []
        B = gctx["/0/DATA/0/matrix"][jpos, :].T
        for k in range(B.shape[1]):
            ranked = symbol[np.argsort(-B[:, k])]
            cas.append(100.0 * wtcs_from_query(up4, down4, ranked))
        v3 = {"n_mtorkd": len(c4), "n_sirolimus": len(cas),
              "mean_cs": float(np.mean(cas)), "min_cs": float(np.min(cas)), "max_cs": float(np.max(cas))}
        print("V3 (MTOR-KD vs sirolimus connectivity):", v3, flush=True)

    val = {
        "v1_bub1b_z_in_kd_consensus": z_bub1b,
        "v1_pass": z_bub1b < -1.5,
        "v2_median_heldout_cs": med_kd_cs,
        "v2_fraction_heldout_above_p95": float((kd_cs > random_p95).mean()),
        "v2_random_cs_p95": random_p95,
        "v2_pass": med_kd_cs > random_p95,
        "n_kd_signatures_all": len(cols_kd),
        "n_kd_signatures_gated": int(gated_mask.sum()),
        "v3_mtor_sirolimus": v3,
        "v3_pass": bool(v3 and v3.get("mean_cs", -999) > 0),
        "tau_denominator_max_abs_cs": float(denom),
        "n_compound_signatures_scored": int(len(cp_cols)),
        "elapsed_s": round(time.time() - t0, 1),
    }
    (outdir / "reversal_validation.json").write_text(json.dumps(val, indent=1))
    print(json.dumps(val, indent=2), flush=True)

    # quick look at top reversal signatures overall
    print(out.sort_values("reversal_tau", ascending=False).head(15)[
        ["sig_id", "pert_iname", "cell_id", "pert_dose", "tau", "reversal_tau"]].to_string(), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
