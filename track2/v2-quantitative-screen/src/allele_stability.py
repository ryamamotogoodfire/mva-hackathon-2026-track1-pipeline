"""Allele stability ticket: RaSP saturation-mutagenesis ddG on AF-O60566-F1.

Runs the KULL-Centre RaSP pipeline (Blaabjerg et al., eLife 2023) exactly as
their maintained Colab recipe, on the AlphaFold v6 BUBR1 model. Reports ddG
(kcal/mol, positive = destabilizing) for every single substitution, with
focal extraction for p.Asn1002Lys and the experimentally validated adjacent
instability allele p.Leu1012Pro (Suijkerbuijk et al., Cancer Res 2010).

Inputs (artifacts): rasp_assets/{bubr1_v6.pdb, colab_additional/, ds_models/,
cavity_models/, PrismData.py, cavity_model.py, extract_environments.py}
Deps come from image torch/pandas + session-installed openmm/biopython.

Output: <outdir>/rasp_bubr1_variants.csv and <outdir>/rasp_focal.json
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")


def sh(cmd: list[str], env: dict | None = None) -> None:
    print("+", " ".join(cmd), flush=True)
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit(r.returncode)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assets", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--pdbid", default="O60566")
    args = ap.parse_args()
    assets = Path(args.assets)
    outdir = Path(args.outdir)
    raw = outdir / "raw"
    cleaned = outdir / "cleaned"
    parsed = outdir / "parsed"
    for d in (raw, cleaned, parsed):
        d.mkdir(parents=True, exist_ok=True)

    pdbid = args.pdbid
    vendor = assets / "vendor"
    if not vendor.exists():
        vendor = Path(__file__).resolve().parent.parent / "vendor" / "rasp"
    src_pdb = assets / "bubr1_v6.pdb"

    # 1) chain select + tidy with pdb-tools (same ops as the notebook cell)
    raw_pdb = raw / f"{pdbid}_uniquechain.pdb"
    sh(["bash", "-c",
        f"pdb_selchain -A {src_pdb} | pdb_delhetatm | pdb_tidy > {raw_pdb}"],
       env=dict(os.environ))
    assert raw_pdb.stat().st_size > 1000, "raw pdb empty"

    reduce_exe = vendor / "reduce"
    # 2) clean with pdbfixer + reduce
    sh([sys.executable, str(vendor / "clean_pdb.py"),
        "--pdb_file_in", str(raw_pdb), "--out_dir", str(cleaned),
        "--reduce_exe", str(reduce_exe)], env=dict(os.environ))
    cpdb = cleaned / f"{pdbid}_uniquechain_clean.pdb"
    assert cpdb.exists(), f"missing cleaned pdb {cpdb}"

    # 3) extract residue environments
    sh([sys.executable, str(vendor / "extract_environments.py"),
        "--pdb_in", str(cpdb), "--out_dir", str(parsed)], env=dict(os.environ))

    # 4) predictions
    sys.path.insert(0, str(vendor))
    os.environ["RASP_DS_MODELS"] = str(assets / "ds_models")
    import pandas as pd  # noqa
    import torch
    from Bio.PDB.Polypeptide import index_to_one, one_to_index
    from cavity_model import CavityModel, DownstreamModel, ResidueEnvironmentsDataset
    from helpers import (cavity_to_prism, ds_pred, get_seq_from_variant,
                         init_lin_weights)

    pdb_filenames_ds = sorted(str(p) for p in parsed.glob("*coord*"))
    dataset_structure = ResidueEnvironmentsDataset(pdb_filenames_ds, transformer=None)
    resenv_dataset = {}
    for resenv in dataset_structure:
        key = (f"--{pdbid}--{resenv.chain_id}--{resenv.pdb_residue_number}--{index_to_one(resenv.restype_index)}--")
        resenv_dataset[key] = resenv
    df_no = pd.DataFrame.from_dict(resenv_dataset, orient="index", columns=["resenv"])
    df_no.reset_index(inplace=True)
    df_no["index"] = df_no["index"].astype(str)
    res_info = pd.DataFrame(df_no["index"].str.split("--").tolist(),
                            columns=["blank", "pdb_id", "chain_id", "pos", "wt_AA", "blank2"])
    df_no["pdbid"] = res_info["pdb_id"]
    df_no["chainid"] = res_info["chain_id"]
    df_no["variant"] = res_info["wt_AA"] + res_info["pos"] + "X"
    aa_list = AA_LIST
    df_structure = pd.DataFrame(df_no.values.repeat(20, axis=0), columns=df_no.columns)
    for i in range(0, len(df_structure), 20):
        for j in range(20):
            df_structure.iloc[i + j, :]["variant"] = df_structure.iloc[i + j, :]["variant"][:-1] + aa_list[j]
    df_structure.drop(columns="index", inplace=True)

    npz = vendor / "pdb_frequencies.npz"
    if not npz.exists():
        npz = assets / "colab_additional" / "pdb_frequencies.npz"
    pdb_nlfs = -np.log(np.load(npz)["frequencies"])
    df_structure["wt_idx"] = df_structure.apply(lambda row: one_to_index(row["variant"][0]), axis=1)
    df_structure["mt_idx"] = df_structure.apply(lambda row: one_to_index(row["variant"][-1]), axis=1)
    df_structure["wt_nlf"] = df_structure.apply(lambda row: pdb_nlfs[row["wt_idx"]], axis=1)
    df_structure["mt_nlf"] = df_structure.apply(lambda row: pdb_nlfs[row["mt_idx"]], axis=1)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print("device:", device, flush=True)
    # notebook uses /content/output/cavity_models/<best>. Remap locally:
    best_path = assets / "cavity_models" / "cavity_model_15.pt"
    cavity_model_net = CavityModel(get_latent=True).to(device)
    cavity_model_net.load_state_dict(torch.load(best_path, map_location=device, weights_only=False))
    cavity_model_net.eval()
    ds_model_net = DownstreamModel().to(device)
    ds_model_net.apply(init_lin_weights)
    ds_model_net.eval()

    # neutered cuda noise if cpu
    df_ml = ds_pred(cavity_model_net, ds_model_net, df_structure,
                    "predictions", 10, device)
    df_total = df_structure.merge(df_ml, on=["pdbid", "chainid", "variant"], how="outer")
    df_total = df_total.drop("resenv", axis=1)
    print("rows:", len(df_total), flush=True)

    df_total = df_total.assign(pos=df_total["variant"].str[1:-1].astype(int))
    outfile = outdir / "rasp_bubr1_variants.csv"
    df_total.to_csv(outfile, index=False)
    print("wrote", outfile, flush=True)

    ddg_col = [c for c in df_total.columns if "ddg" in c.lower() or "score" in c.lower()]
    print("prediction columns:", [c for c in df_total.columns], flush=True)

    # sanity self-check: variance of ddG for WT->WT self substitutions must be ~0
    if "ddG_ml" in df_total.columns:
        col = "ddG_ml"
    elif "score_ml" in df_total.columns:
        col = "score_ml"
    else:
        col = ddg_col[0]
    self_subs = df_total[df_total.variant.str[0] == df_total.variant.str[-1]]
    foc_sub = df_total[df_total.pos.isin([1002, 1012])]
    out = {
        "ddg_col": col,
        "device": device,
        "n_variants": int(len(df_total)),
        "self_substitutions_ddG_range": [float(self_subs[col].min()), float(self_subs[col].max())],
        "n_positions": int(df_total.pos.nunique()),
        "pos_1002": foc_sub[foc_sub.pos == 1002][["variant", col]].to_dict("records"),
        "pos_1012": foc_sub[foc_sub.pos == 1012][["variant", col]].to_dict("records"),
    }
    (outdir / "rasp_focal.json").write_text(json.dumps(out, indent=1))
    for v in ["N1002K", "L1012P"]:
        sel = df_total[df_total.variant == v]
        if len(sel):
            print(v, "RaSP ddG_ml =", float(sel[col].iloc[0]), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
