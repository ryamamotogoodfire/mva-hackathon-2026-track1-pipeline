#!/usr/bin/env python3
"""Report figures: evidence funnel + GPN-MSA lane distribution.

Reads the pipeline intermediates and writes two standalone HTML figures with
plain plotly defaults:
  figures/funnel.html              candidate funnel from the full call set
  figures/gpn_msa_distribution.html  GPN-MSA LLR violin, candidates vs controls

Paths resolve under $DATA_DIR/mva-track1 unless overridden:
  FUNNEL_JSON   mining/funnel.json
  GPN_SCORED    scoring/gpn_msa_scored_variants.parquet
  EVEE_SCORES   evee/evee_scores.parquet
"""
import json
import os
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

DATA_DIR = Path(os.environ.get("DATA_DIR", ".")) / "mva-track1"
OUT = Path(os.environ.get("FIGURES_DIR", DATA_DIR / "figures"))
OUT.mkdir(parents=True, exist_ok=True)

FUNNEL = os.environ.get("FUNNEL_JSON", DATA_DIR / "mining" / "funnel.json")
GPNP = os.environ.get("GPN_SCORED", DATA_DIR / "scoring" / "gpn_msa_scored_variants.parquet")
EVEE = os.environ.get("EVEE_SCORES", DATA_DIR / "evee" / "evee_scores.parquet")

funnel = json.load(open(FUNNEL))

stages = pd.DataFrame([
    ("Called records in VCF", 5012204),
    ("Annotated, joined variants", funnel["vep_variants_collapsed"]),
    ("Rare (AF <= 0.01)", funnel["rare_variants_primary_tier"] + funnel["rare_variants_loose_tier"]),
    ("Scored candidates (SNV)", 264),
    ("Genes on evidence board", 267),
    ("Compound-het genes", funnel["genes_with_ge2_qual_het"]),
    ("Submission rows", 10),
], columns=["stage", "n"])

fig = px.bar(stages, x="stage", y="n", text_auto=False)
fig.update_yaxes(type="log", title="count (log)")
fig.update_xaxes(title="pipeline stage (left to right)")
fig.write_html(str(OUT / "funnel.html"))

gp = pd.read_parquet(GPNP)
ev = pd.read_parquet(EVEE)
evk = ev.assign(ukey=ev["chrom"].astype(str) + "_" + ev["pos"].astype(str) + "_" + ev["ref"] + "/" + ev["alt"])[["ukey", "evee_score"]]
sk = gp.assign(ukey=gp["chrom"].astype(str) + "_" + gp["pos"].astype(str) + "_" + gp["ref"] + "/" + gp["alt"])[["ukey", "gpn_msa_score", "role"]]
df = sk.merge(evk, on="ukey", how="left")
df_bub = df[df["ukey"].isin(["15_40209701_T/G", "15_40220612_T/G"])]

fig2 = go.Figure()
dfc = df[df["role"] == "candidate"]
dfi = df[df["role"] == "inert_control"]
fig2.add_violin(x=dfc["role"].map(lambda _: "candidates"), y=dfc["gpn_msa_score"],
                name=f"scored candidates (n={len(dfc)})", side="positive", spanmode="hard",
                box_visible=True, meanline_visible=True)
fig2.add_violin(x=dfi["role"].map(lambda _: "inert controls"), y=dfi["gpn_msa_score"],
                name=f"inert controls (n={len(dfi)})", side="positive", spanmode="hard",
                box_visible=True, meanline_visible=True)
for _, r in df_bub.iterrows():
    fig2.add_scatter(x=["candidates"], y=[r["gpn_msa_score"]], mode="markers",
                     name=(f"BUB1B {r['ukey'].split('_',1)[1]} EVEE {r['evee_score']:.3f}"),
                     marker=dict(size=14, symbol="diamond-open"))
fig2.update_yaxes(title="GPN-MSA LLR (alt-ref); lower = more damaging")
fig2.update_xaxes(title="pool")
fig2.write_html(str(OUT / "gpn_msa_distribution.html"))
print("FIGS_OK")
