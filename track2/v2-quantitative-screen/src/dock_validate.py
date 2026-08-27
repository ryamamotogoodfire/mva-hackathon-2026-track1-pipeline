"""Docking instrument checks before the bounded stabilizer screen.

Checks, in order:
  1. p2rank (ML pocket predictor, AlphaFold model) on human AF-O60566-F1 and on the fly crystal 6JKM.
  2. ADP redock into its own crystal (6JKM): does the engine recover the experimental pose?
  3. ADP / ATP docked into the human pseudokinase site box (mapped from the fly structure).
  4. Throughput timing on a handful of approved drugs, to size the batch screen.

Usage: python src/dock_validate.py <workdir> <outjson>
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import dock_common as dc  # noqa: E402

EXP = HERE.parent
PDB = EXP / "inputs" / "pdb"


def rmsd_to_reference(pose_mol, ref_pdb_lines: list[str], template_smiles: str) -> float:
    """Best heavy-atom RMSD (symmetry-aware) between a docked pose and the crystal ligand."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    ref_block = "\n".join(ref_pdb_lines) + "\nEND\n"
    ref = Chem.MolFromPDBBlock(ref_block, removeHs=True, sanitize=False)
    tmpl = Chem.MolFromSmiles(template_smiles)
    try:
        ref = AllChem.AssignBondOrdersFromTemplate(tmpl, ref)
    except Exception:
        pass
    pose = Chem.RemoveHs(Chem.Mol(pose_mol))
    try:
        pose = AllChem.AssignBondOrdersFromTemplate(tmpl, pose)
    except Exception:
        pass
    try:
        return float(rdMolAlign.CalcRMS(pose, ref))
    except Exception:
        # fall back: nearest-neighbour heavy atom deviation (no re-alignment; both are in the crystal frame)
        p = pose.GetConformer().GetPositions()
        r = ref.GetConformer().GetPositions()
        d = np.linalg.norm(p[:, None, :] - r[None, :, :], axis=2)
        return float(np.sqrt((d.min(axis=1) ** 2).mean()))


def p2rank(tools: dict, pdb: Path, outdir: Path, config: str) -> list[dict]:
    env = dict(os.environ)
    env["PATH"] = tools["java_bin"] + os.pathsep + env["PATH"]
    env["JAVA_HOME"] = str(Path(tools["java_bin"]).parent)
    cmd = f'"{tools["p2rank"]}" predict -f "{pdb}" -c {config} -o "{outdir}" -threads 4'
    subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, env=env)
    pred = next(outdir.glob("*predictions.csv"))
    rows = []
    lines = [l.strip() for l in pred.read_text().splitlines() if l.strip()]
    header = [h.strip() for h in lines[0].split(",")]
    for line in lines[1:]:
        vals = [v.strip() for v in line.split(",")]
        rec = dict(zip(header, vals))
        rows.append({
            "rank": int(rec["rank"]), "name": rec["name"], "score": float(rec["score"]),
            "probability": float(rec.get("probability", "nan")),
            "center": [float(rec["center_x"]), float(rec["center_y"]), float(rec["center_z"])],
            "residue_ids": rec.get("residue_ids", ""),
        })
    return rows


def main() -> None:
    workdir = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/dockwork")
    outjson = Path(sys.argv[2] if len(sys.argv) > 2 else EXP / "results" / "dock_validation.json")
    workdir.mkdir(parents=True, exist_ok=True)
    tools = dc.ensure_tools(workdir)
    res: dict = {"tools": {k: v for k, v in tools.items() if k != "p2rank_root"}}

    site = json.loads((EXP / "results" / "atp_site_annotation.json").read_text())
    human_site = site["human_atp_site_residues"]

    af = PDB / "AF-O60566-F1-model_v6.pdb"
    jkm = PDB / "6jkm.pdb"

    # ---- 1. ML pocket prediction
    t0 = time.time()
    res["p2rank_af_human"] = p2rank(tools, af, workdir / "p2rank_af", "alphafold")[:10]
    res["p2rank_6jkm_fly"] = p2rank(tools, jkm, workdir / "p2rank_jkm", "default")[:10]
    res["p2rank_seconds"] = round(time.time() - t0, 1)

    # p2rank pocket closest to the mapped human ATP site, and to N1002
    af_site = dc.residue_atom_coords(af, human_site)
    n1002 = dc.residue_atom_coords(af, [1002])
    for p in res["p2rank_af_human"]:
        c = np.array(p["center"])
        p["dist_to_atp_site_A"] = round(float(np.linalg.norm(af_site - c, axis=1).min()), 2)
        p["dist_to_N1002_A"] = round(float(np.linalg.norm(n1002 - c, axis=1).min()), 2)

    # ---- 2. ADP redock into its own crystal
    prot, het = dc.read_pdb_atoms(jkm)
    adp_lines = het["ADP"]
    rec_jkm = dc.write_receptor(jkm, workdir / "6jkm_receptor.pdb")
    adp_ref_pdb = dc.het_to_pdb(adp_lines, workdir / "adp_ref.pdb")
    adp_coords = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in adp_lines])
    box_cryst = dc.box_from_coords(adp_coords, pad=6.0, min_size=20.0)
    n_ok, failed = dc.build_ligand_sdf([{"id": "ADP", "smiles": dc.ADP_SMILES, "name": "ADP"}],
                                      workdir / "adp.sdf")
    t0 = time.time()
    dc.run_smina(tools["smina"], rec_jkm, workdir / "adp.sdf", box_cryst, workdir / "adp_redock.sdf",
                 exhaustiveness=16, cpu=8, num_modes=9)
    redock_s = time.time() - t0
    poses = dc.poses_by_name(workdir / "adp_redock.sdf").get("ADP", [])
    scored = []
    for m in poses:
        try:
            aff = float(m.GetProp("minimizedAffinity"))
        except Exception:
            continue
        scored.append({"affinity": aff, "rmsd_to_crystal_A": round(rmsd_to_reference(m, adp_lines, dc.ADP_SMILES), 2)})
    res["adp_redock_6jkm"] = {
        "box": box_cryst, "seconds": round(redock_s, 1), "n_poses": len(scored),
        "top_pose": scored[0] if scored else None,
        "best_rmsd_pose": min(scored, key=lambda s: s["rmsd_to_crystal_A"]) if scored else None,
        "all_poses": scored,
        "embed_failures": failed,
    }

    # ---- 3. nucleotide docking into the human pseudokinase site
    rec_af = dc.write_receptor(af, workdir / "af_receptor.pdb")
    box_human = dc.box_from_coords(af_site, pad=5.0, min_size=20.0)
    dc.build_ligand_sdf([{"id": "ADP", "smiles": dc.ADP_SMILES}, {"id": "ATP", "smiles": dc.ATP_SMILES}],
                        workdir / "nucleotides.sdf")
    dc.run_smina(tools["smina"], rec_af, workdir / "nucleotides.sdf", box_human,
                 workdir / "nucleotides_af.sdf", exhaustiveness=16, cpu=8, num_modes=5)
    res["nucleotide_docking_human_af"] = {
        "box": box_human,
        "scores": {k: round(v["affinity"], 2) for k, v in dc.parse_smina_sdf(workdir / "nucleotides_af.sdf").items()},
    }
    # same nucleotides into the fly crystal site, as the cross-species reference
    res["nucleotide_docking_fly_crystal"] = {
        "scores": {k: round(v["affinity"], 2) for k, v in dc.parse_smina_sdf(workdir / "adp_redock.sdf").items()}
    }

    # ---- 4. throughput timing on a few approved drugs
    probe = [
        {"id": "imatinib", "smiles": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"},
        {"id": "hydroxychloroquine", "smiles": "CCN(CCO)CCCC(C)Nc1ccnc2cc(Cl)ccc12"},
        {"id": "tafamidis", "smiles": "OC(=O)c1cc2cc(Cl)cc(Cl)c2o1"},
        {"id": "sirolimus_fragment_ref", "smiles": "CC1CCC2CC(OC)C(=O)OC2C1"},
        {"id": "niacin", "smiles": "OC(=O)c1cccnc1"},
    ]
    dc.build_ligand_sdf(probe, workdir / "probe.sdf")
    t0 = time.time()
    dc.run_smina(tools["smina"], rec_af, workdir / "probe.sdf", box_human, workdir / "probe_out.sdf",
                 exhaustiveness=8, cpu=1, num_modes=3)
    dt = time.time() - t0
    res["throughput"] = {"n_ligands": len(probe), "exhaustiveness": 8, "cpu": 1,
                         "seconds_total": round(dt, 1), "seconds_per_ligand_1cpu": round(dt / len(probe), 2),
                         "probe_scores": {k: round(v["affinity"], 2) for k, v in dc.parse_smina_sdf(workdir / "probe_out.sdf").items()}}

    outjson.parent.mkdir(parents=True, exist_ok=True)
    outjson.write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ["tools", "adp_redock_6jkm", "nucleotide_docking_human_af", "throughput"]}, indent=1)[:4000])
    print("\ntop p2rank pockets (human AF):")
    for p in res["p2rank_af_human"][:6]:
        print(f'  rank {p["rank"]} score {p["score"]:.2f} prob {p["probability"]:.2f} '
              f'd(ATP site)={p["dist_to_atp_site_A"]} d(N1002)={p["dist_to_N1002_A"]}')
    print("wrote", outjson)


if __name__ == "__main__":
    main()
