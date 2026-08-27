"""Persist molecule_chembl_id -> standard InChIKey mapping for joins."""
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
BASE = "https://www.ebi.ac.uk/chembl/api/data"
S = requests.Session()
S.headers["User-Agent"] = "SilicoMVA/1.0 (research)"


def main() -> None:
    mols = pd.read_parquet(RESULTS / "chembl_molecules.parquet", columns=["molecule_chembl_id"])
    ids = mols.molecule_chembl_id.tolist()
    rows = []
    for i in range(0, len(ids), 50):
        batch = ids[i:i + 50]
        for attempt in range(4):
            r = S.get(f"{BASE}/molecule.json?molecule_chembl_id__in={','.join(batch)}&limit=1000", timeout=60)
            if r.status_code == 200:
                break
        r.raise_for_status()
        for m in r.json().get("molecules", []):
            rows.append({
                "molecule_chembl_id": m["molecule_chembl_id"],
                "pref_name": m.get("pref_name"),
                "standard_inchi_key": (m.get("molecule_structures") or {}).get("standard_inchi_key"),
                "inchi_key_prefix": ((m.get("molecule_structures") or {}).get("standard_inchi_key") or "").split("-")[0] or None,
            })
    df = pd.DataFrame(rows)
    df.to_parquet(RESULTS / "molecule_inchikeys.parquet", index=False)
    print(len(df), "molecules,", df.standard_inchi_key.notna().sum(), "with inchikey")


if __name__ == "__main__":
    main()
