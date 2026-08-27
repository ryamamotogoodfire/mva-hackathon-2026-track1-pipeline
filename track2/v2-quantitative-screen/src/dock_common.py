"""Shared helpers for the bounded BUBR1 stabilizer docking screen.

Toolchain, all job-time installed / downloaded (no baked deps touched):
  smina (static, AutoDock Vina 1.1.2 scoring) as the docking engine
  p2rank 2.4.2 (ML pocket predictor, AlphaFold-specific model) via jdk4py's JRE
  RDKit for ligand 3D generation
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

SMINA_URL = "https://sourceforge.net/projects/smina/files/smina.static/download"
P2RANK_URL = "https://github.com/rdk/p2rank/releases/download/2.4.2/p2rank_2.4.2.tar.gz"
# p2rank 2.4.2 bundles Groovy 4.0.21, which refuses Java 25 class files, so pin a Temurin 21 JRE
JRE_URL = "https://api.adoptium.net/v3/binary/latest/21/ga/linux/x64/jre/hotspot/normal/eclipse"

AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
       "MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}

ADP_SMILES = "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"
ATP_SMILES = "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"


def sh(cmd: str, **kw) -> str:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, **kw).stdout


def ensure_tools(workdir: Path) -> dict:
    """Download smina + p2rank once into workdir; return tool paths."""
    workdir.mkdir(parents=True, exist_ok=True)
    smina = workdir / "smina.static"
    if not smina.exists() or smina.stat().st_size < 1_000_000:
        sh(f'curl -sL -o "{smina}" "{SMINA_URL}"')
        smina.chmod(0o755)
    ver = subprocess.run([str(smina), "--version"], capture_output=True, text=True).stdout.strip()

    p2root = workdir / "p2rank_2.4.2"
    if not p2root.exists():
        tgz = workdir / "p2rank.tar.gz"
        if not tgz.exists():
            sh(f'curl -sL -o "{tgz}" "{P2RANK_URL}"')
        sh(f'tar xzf "{tgz}" -C "{workdir}"')

    jres = sorted(workdir.glob("jdk-21*/bin"))
    if not jres:
        tgz = workdir / "jre21.tar.gz"
        if not tgz.exists():
            sh(f'curl -sL -o "{tgz}" "{JRE_URL}"')
        sh(f'tar xzf "{tgz}" -C "{workdir}"')
        jres = sorted(workdir.glob("jdk-21*/bin"))
    java_bin = jres[0]
    java_ver = subprocess.run([str(java_bin / "java"), "-version"], capture_output=True, text=True).stderr.splitlines()[0]
    return {"smina": str(smina), "smina_version": ver, "p2rank": str(p2root / "prank"),
            "java_bin": str(java_bin), "java_version": java_ver, "p2rank_root": str(p2root)}


# ---------------------------------------------------------------- structures

def read_pdb_atoms(path: Path):
    """Return (protein_records, hetero_by_resname) as lists of raw PDB lines."""
    prot, het = [], {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("ATOM"):
            prot.append(line)
        elif line.startswith("HETATM"):
            het.setdefault(line[17:20].strip(), []).append(line)
    return prot, het


def write_receptor(src: Path, dest: Path) -> Path:
    """Protein-only receptor PDB (drops HETATM/waters/ligands)."""
    prot, _ = read_pdb_atoms(src)
    dest.write_text("\n".join(prot) + "\nEND\n")
    return dest


def residue_atom_coords(pdb: Path, resnums: list[int]) -> np.ndarray:
    want = set(resnums)
    out = []
    for line in Path(pdb).read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            num = int(line[22:26])
        except ValueError:
            continue
        if num in want:
            out.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(out)


def box_from_coords(coords: np.ndarray, pad: float = 5.0, min_size: float = 18.0) -> dict:
    lo, hi = coords.min(axis=0) - pad, coords.max(axis=0) + pad
    center = (lo + hi) / 2.0
    size = np.maximum(hi - lo, min_size)
    return {"center": [round(float(c), 3) for c in center],
            "size": [round(float(s), 3) for s in size]}


def het_to_pdb(lines: list[str], dest: Path) -> Path:
    dest.write_text("\n".join(lines) + "\nEND\n")
    return dest


# ---------------------------------------------------------------- ligands

def build_ligand_sdf(records: list[dict], dest: Path, seed: int = 42, max_confs: int = 1) -> tuple[int, list[str]]:
    """3D-embed SMILES records [{name, smiles, id}] into one multi-mol SDF."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    ok, failed = 0, []
    w = Chem.SDWriter(str(dest))
    for rec in records:
        try:
            m = Chem.MolFromSmiles(rec["smiles"])
            if m is None:
                failed.append(rec["id"]); continue
            m = Chem.AddHs(m)
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            if AllChem.EmbedMolecule(m, params) != 0:
                if AllChem.EmbedMolecule(m, useRandomCoords=True, randomSeed=seed) != 0:
                    failed.append(rec["id"]); continue
            try:
                AllChem.MMFFOptimizeMolecule(m, maxIters=400)
            except Exception:
                pass
            m.SetProp("_Name", rec["id"])
            for k, v in rec.items():
                if k != "smiles":
                    m.SetProp(str(k), str(v))
            w.write(m)
            ok += 1
        except Exception:
            failed.append(rec["id"])
    w.close()
    return ok, failed


# ---------------------------------------------------------------- docking

def run_smina(smina: str, receptor: Path, ligands: Path, box: dict, out_sdf: Path,
              exhaustiveness: int = 8, cpu: int = 8, seed: int = 42,
              scoring: str | None = None, num_modes: int = 3, extra: str = "") -> str:
    cx, cy, cz = box["center"]
    sx, sy, sz = box["size"]
    cmd = (f'"{smina}" -r "{receptor}" -l "{ligands}" -o "{out_sdf}" '
           f'--center_x {cx} --center_y {cy} --center_z {cz} '
           f'--size_x {sx} --size_y {sy} --size_z {sz} '
           f'--exhaustiveness {exhaustiveness} --num_modes {num_modes} --seed {seed} --cpu {cpu} --quiet {extra}')
    if scoring:
        cmd += f" --scoring {scoring}"
    subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    return cmd


def parse_smina_sdf(path: Path) -> dict[str, dict]:
    """Best (most negative) minimizedAffinity per ligand name from a smina output SDF."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    best: dict[str, dict] = {}
    supp = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    for m in supp:
        if m is None:
            continue
        name = m.GetProp("_Name") if m.HasProp("_Name") else None
        if not name:
            continue
        try:
            aff = float(m.GetProp("minimizedAffinity"))
        except Exception:
            continue
        rec = best.get(name)
        if rec is None or aff < rec["affinity"]:
            props = {k: m.GetProp(k) for k in m.GetPropNames() if k not in ("minimizedAffinity",)}
            best[name] = {"affinity": aff, "props": props}
    return best


def poses_by_name(path: Path):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    out = {}
    for m in Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False):
        if m is None:
            continue
        nm = m.GetProp("_Name") if m.HasProp("_Name") else "?"
        out.setdefault(nm, []).append(m)
    return out
