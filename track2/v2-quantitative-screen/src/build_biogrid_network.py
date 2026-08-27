"""Build the human physical PPI network from BioGRID 4.4.236 tab3.

Input: /tmp/biogrid_organism_44236.zip (BIOGRID-ORGANISM-4.4.236.tab3.zip)
Filters: organism 9606 both interactors, EXPERIMENTAL_SYSTEM_TYPE == Physical,
gene-symbol axis, self interactions removed, largest connected component kept.

Outputs:
  results/biogrid_edges.parquet    columns a, b (gene symbols), biogrid_id
  results/biogrid_network_stats.json
"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import networkx as nx
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
ZIP = Path("/tmp/biogrid_organism_44236.zip")


def main() -> None:
    with zipfile.ZipFile(ZIP) as z:
        name = [n for n in z.namelist() if "Homo_sapiens" in n and n.endswith(".txt")][0]
        print("reading", name, flush=True)
        df = pd.read_csv(
            z.open(name),
            sep="\t",
            low_memory=False,
            dtype=str,
            usecols=[
                "#BioGRID Interaction ID",
                "Official Symbol Interactor A",
                "Official Symbol Interactor B",
                "Experimental System Type",
                "Organism ID Interactor A",
                "Organism ID Interactor B",
            ],
        )
    print("raw rows:", len(df), flush=True)

    df = df[
        (df["Organism ID Interactor A"] == "9606")
        & (df["Organism ID Interactor B"] == "9606")
        & (df["Experimental System Type"].str.lower() == "physical")
    ]
    print("human physical rows:", len(df), flush=True)

    df = df[["#BioGRID Interaction ID", "Official Symbol Interactor A", "Official Symbol Interactor B"]]
    df.columns = ["biogrid_id", "a", "b"]
    df = df.dropna(subset=["a", "b"])
    df["a"] = df["a"].str.upper().str.strip()
    df["b"] = df["b"].str.upper().str.strip()
    df = df[df["a"] != df["b"]]
    df = df.drop_duplicates(subset=["a", "b"])

    G = nx.Graph()
    G.add_edges_from(df[["a", "b"]].itertuples(index=False))
    lcc = max(nx.connected_components(G), key=len)
    Gl = G.subgraph(lcc).copy()
    keep = set(Gl.nodes())
    out = df[df.a.isin(keep) & df.b.isin(keep)].reset_index(drop=True)
    out.to_parquet(RESULTS / "biogrid_edges.parquet", index=False)

    deg = dict(Gl.degree())
    stats = {
        "source": "BIOGRID-ORGANISM-4.4.236.tab3.zip",
        "n_raw_rows": int(len(df)),
        "n_edges_lcc": int(Gl.number_of_edges()),
        "n_nodes_lcc": int(Gl.number_of_nodes()),
        "n_components": int(nx.number_connected_components(G)),
        "lcc_share": round(len(Gl) / G.number_of_nodes(), 4),
        "degree_mean": round(sum(deg.values()) / len(deg), 2),
        "degree_max": max(deg.values()),
    }
    print(json.dumps(stats, indent=1), flush=True)
    (RESULTS / "biogrid_network_stats.json").write_text(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
