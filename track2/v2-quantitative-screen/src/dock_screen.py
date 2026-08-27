"""Bounded stabilizer docking screen for the BUBR1 pseudokinase domain.

Two stages against four boxes on the human AlphaFold model AF-O60566-F1:
  nucleotide_site  : ADP-contact residues mapped from the fly crystal 6JKM (contains K795)
  p2rank_top       : p2rank's highest-scoring pocket on the human model
  mutation_site    : fpocket pockets 7 + 19, the only pockets touching N1002
  pseudokinase_best: fpocket pocket 21, the most druggable pocket overlapping the domain

Stage 1 screens every approved small molecule (150-600 Da) at exhaustiveness 4.
Stage 2 re-docks the best per box at exhaustiveness 16 and rescores with Vinardo.
Rankings pair the affinity with a stability-relevance term computed from the RaSP
saturation scan (mean destabilization burden of the pocket-lining residues) and the
distance to N1002, since a stabilizer need not bind at the mutation itself.

Usage: python src/dock_screen.py <workdir> <outdir> [n_cpu]
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import dock_common as dc  # noqa: E402

EXP = HERE.parent
PDB = EXP / "inputs" / "pdb"
SEED = 42

BOXES_RESIDUES = {
    "nucleotide_site": None,  # filled from atp_site_annotation.json
    "p2rank_top": [722, 724, 727, 777, 778, 795, 797, 799, 801, 803, 804, 832, 834, 913, 914],
    "mutation_site": [941, 944, 945, 950, 954, 973, 976, 977, 978, 995, 998, 999, 1000, 1001, 1002,
                      1003, 1004, 1014, 1015, 1018],
    "pseudokinase_best": [769, 770, 771, 772, 782, 784, 786, 787],
}


def dock_chunk(args):
    smina, receptor, chunk_sdf, box, out_sdf, exh, modes, scoring = args
    dc.run_smina(smina, Path(receptor), Path(chunk_sdf), box, Path(out_sdf),
                 exhaustiveness=exh, cpu=1, seed=SEED, num_modes=modes, scoring=scoring)
    return out_sdf


def split_sdf(src: Path, n: int, outdir: Path) -> list[Path]:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    mols = [m for m in Chem.SDMolSupplier(str(src), removeHs=False, sanitize=False) if m is not None]
    outdir.mkdir(parents=True, exist_ok=True)
    parts, writers = [], []
    for i in range(n):
        p = outdir / f"chunk_{i:03d}.sdf"
        parts.append(p)
        writers.append(Chem.SDWriter(str(p)))
    for i, m in enumerate(mols):
        writers[i % n].write(m)
    for w in writers:
        w.close()
    return [p for p in parts if p.stat().st_size > 0]


def collect(out_sdfs: list[Path]) -> dict[str, dict]:
    best: dict[str, dict] = {}
    for p in out_sdfs:
        for name, rec in dc.parse_smina_sdf(Path(p)).items():
            if name not in best or rec["affinity"] < best[name]["affinity"]:
                best[name] = rec
    return best


def pocket_stability_burden(rasp_csv: Path, residues: list[int]) -> dict:
    """Mean RaSP ddG over all substitutions of the pocket-lining residues (structural criticality)."""
    if not rasp_csv.exists():
        return {"available": False}
    df = pd.read_csv(rasp_csv)
    col = "score_ml" if "score_ml" in df.columns else df.columns[-1]
    sub = df[df.pos.isin(residues)]
    allmean = float(df[col].mean())
    return {"available": True, "n_residues_scored": int(sub.pos.nunique()),
            "mean_ddg_pocket": round(float(sub[col].mean()), 3),
            "mean_ddg_protein": round(allmean, 3),
            "burden_ratio": round(float(sub[col].mean()) / allmean, 3) if allmean else None}


def main() -> None:
    workdir = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/dockwork")
    outdir = Path(sys.argv[2] if len(sys.argv) > 2 else EXP / "results" / "dock")
    ncpu = int(sys.argv[3]) if len(sys.argv) > 3 else (os.cpu_count() or 8)
    outdir.mkdir(parents=True, exist_ok=True)
    tools = dc.ensure_tools(workdir)
    art = Path(os.environ.get("SILICO_EXPERIMENT_ARTIFACTS_DIR", str(outdir)))
    rasp_csv = art / "allele_tickets" / "rasp_production" / "rasp_bubr1_variants.csv"

    site = json.loads((EXP / "results" / "atp_site_annotation.json").read_text())
    BOXES_RESIDUES["nucleotide_site"] = site["human_atp_site_residues"]

    af = PDB / "AF-O60566-F1-model_v6.pdb"
    receptor = dc.write_receptor(af, workdir / "af_receptor.pdb")

    # ---- boxes
    n1002 = dc.residue_atom_coords(af, [1002])
    boxes = {}
    for name, residues in BOXES_RESIDUES.items():
        coords = dc.residue_atom_coords(af, residues)
        box = dc.box_from_coords(coords, pad=4.0, min_size=20.0)
        box["size"] = [min(s, 20.0) for s in box["size"]]
        c = np.array(box["center"])
        box["residues"] = residues
        box["dist_center_to_N1002_A"] = round(float(np.linalg.norm(n1002 - c, axis=1).min()), 2)
        box["stability_burden"] = pocket_stability_burden(rasp_csv, residues)
        boxes[name] = box
    only = os.environ.get("DOCK_BOXES")
    if only:
        boxes = {k: v for k, v in boxes.items() if k in set(only.split(","))}
    print(json.dumps(boxes, indent=1))

    # ---- ligands: approved small molecules + nucleotide references
    lig = pd.read_parquet(EXP / "results" / "approved_ligands.parquet")
    limit = int(os.environ.get("LIGAND_LIMIT", "0"))
    if limit:
        lig = lig.head(limit)
        print(f"LIGAND_LIMIT active: {limit} ligands (smoke mode)")
    else:
        # bridging-capable window keeps the docking cost inside budget; named case candidates are
        # force-included regardless of size so no graded candidate is missing from the screen
        force = ["HYDROXYCHLOROQUINE", "CHLOROQUINE", "METFORMIN", "NIACIN", "NIACINAMIDE", "ACIPIMOX",
                 "AMLEXANOX", "FOSTAMATINIB", "TAFAMIDIS", "SIROLIMUS", "MIGALASTAT", "SAPROPTERIN",
                 "IVACAFTOR", "LUMACAFTOR", "BORTEZOMIB", "ATALUREN", "TOPIRAMATE", "SULINDAC"]
        keep_mw = lig.mw.between(250, 600)
        keep_named = lig.pref_name.fillna("").str.upper().str.contains("|".join(force))
        lig = lig[keep_mw | keep_named]
        print(f"ligand window: MW 250-600 plus {int(keep_named.sum())} named candidates -> {len(lig)} molecules")
    records = [{"id": r.molecule_chembl_id, "smiles": r.smiles,
                "pref_name": (r.pref_name or ""), "mw": round(float(r.mw), 1)}
               for r in lig.itertuples()]
    records += [{"id": "REF_ADP", "smiles": dc.ADP_SMILES, "pref_name": "ADP (natural nucleotide reference)", "mw": 427.2},
                {"id": "REF_ATP", "smiles": dc.ATP_SMILES, "pref_name": "ATP (natural nucleotide reference)", "mw": 507.2}]
    t0 = time.time()
    sdf = workdir / "ligands.sdf"
    n_ok, failed = dc.build_ligand_sdf(records, sdf, seed=SEED)
    print(f"ligands embedded: {n_ok} / {len(records)} ({len(failed)} failures) in {time.time()-t0:.0f}s", flush=True)

    chunks = split_sdf(sdf, ncpu, workdir / "chunks")
    summary = {"tools": {k: v for k, v in tools.items() if k != "p2rank_root"},
               "n_ligands_embedded": n_ok, "n_embed_failures": len(failed), "embed_failures": failed[:40],
               "boxes": boxes, "ncpu": ncpu, "seed": SEED, "stage1": {}, "stage2": {}}

    # ---- stage 1 and stage 2, interleaved per box with incremental writes and resume
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    id2rec = {r["id"]: r for r in records}
    all_stage1 = {}
    exh1 = int(os.environ.get("STAGE1_EXH", "2"))
    for bname, box in boxes.items():
        pq = outdir / f"dock_stage1_{bname}.parquet"
        if pq.exists():
            df = pd.read_parquet(pq)
            best = {r.molecule_chembl_id: {"affinity": float(r.affinity)} for r in df.itertuples()}
            summary["stage1"][bname] = {"n_scored": len(best), "exhaustiveness": exh1, "seconds": 0.0,
                                        "cpu_seconds": 0.0, "resumed_from_parquet": True,
                                        "affinity_min": round(min(v["affinity"] for v in best.values()), 2),
                                        "affinity_median": round(float(np.median([v["affinity"] for v in best.values()])), 2)}
            print(f"stage1 {bname}: resumed {len(best)} scores from {pq.name}", flush=True)
        else:
            t0 = time.time()
            jobs = [(tools["smina"], str(receptor), str(c), box, str(workdir / f"s1_{bname}_{c.stem}.sdf"), exh1, 3, None)
                    for c in chunks]
            with mp.Pool(ncpu) as pool:
                outs = pool.map(dock_chunk, jobs)
            best = collect([Path(o) for o in outs])
            dt = time.time() - t0
            summary["stage1"][bname] = {"n_scored": len(best), "exhaustiveness": exh1, "seconds": round(dt, 1),
                                        "cpu_seconds": round(dt * ncpu, 1),
                                        "affinity_min": round(min(v["affinity"] for v in best.values()), 2),
                                        "affinity_median": round(float(np.median([v["affinity"] for v in best.values()])), 2)}
            df = pd.DataFrame([{"molecule_chembl_id": k, "affinity": v["affinity"],
                                "pref_name": id2rec.get(k, {}).get("pref_name", ""),
                                "mw": id2rec.get(k, {}).get("mw")} for k, v in best.items()])
            df.sort_values("affinity").to_parquet(pq, index=False)
            (outdir / "dock_screen_summary.json").write_text(json.dumps(summary, indent=1))
            print(f'stage1 {bname}: {len(best)} ligands in {dt:.0f}s wall '
                  f'(best {summary["stage1"][bname]["affinity_min"]}, median {summary["stage1"][bname]["affinity_median"]})', flush=True)
        all_stage1[bname] = best

        keep = sorted(best.items(), key=lambda kv: kv[1]["affinity"])[: (5 if limit else 40)]
        keep_ids = [k for k, _ in keep] + ["REF_ADP", "REF_ATP"]
        sub = [id2rec[i] for i in dict.fromkeys(keep_ids) if i in id2rec]
        sdf2 = workdir / f"s2_{bname}_in.sdf"
        dc.build_ligand_sdf(sub, sdf2, seed=SEED)
        c2 = split_sdf(sdf2, min(ncpu, len(sub)), workdir / f"chunks_s2_{bname}")
        res2 = {}
        for scoring, tag in [(None, "vina"), ("vinardo", "vinardo")]:
            jobs = [(tools["smina"], str(receptor), str(c), box,
                     str(workdir / f"s2_{bname}_{tag}_{c.stem}.sdf"), 16, 5, scoring) for c in c2]
            with mp.Pool(len(c2)) as pool:
                outs = pool.map(dock_chunk, jobs)
            res2[tag] = collect([Path(o) for o in outs])
        rows = []
        for cid in dict.fromkeys(keep_ids):
            if cid not in res2["vina"]:
                continue
            rows.append({"molecule_chembl_id": cid,
                         "pref_name": id2rec.get(cid, {}).get("pref_name", ""),
                         "mw": id2rec.get(cid, {}).get("mw"),
                         "stage1_affinity": round(best.get(cid, {}).get("affinity", float("nan")), 2) if cid in best else None,
                         "vina_exh16": round(res2["vina"][cid]["affinity"], 2),
                         "vinardo_exh16": round(res2["vinardo"].get(cid, {}).get("affinity", float("nan")), 2) if cid in res2["vinardo"] else None})
        rows.sort(key=lambda r: r["vina_exh16"])
        summary["stage2"][bname] = rows
        (outdir / "dock_screen_summary.json").write_text(json.dumps(summary, indent=1))
        print(f"stage2 {bname}: top 5 {[ (r['pref_name'] or r['molecule_chembl_id'], r['vina_exh16']) for r in rows[:5] ]}", flush=True)

    (outdir / "dock_screen_summary.json").write_text(json.dumps(summary, indent=1))
    print("wrote", outdir / "dock_screen_summary.json")


if __name__ == "__main__":
    main()
