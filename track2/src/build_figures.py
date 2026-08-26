#!/usr/bin/env python3
"""Figures for the Track 2 results page. Renders from saved results JSON.

F1 transcript_geometry: BUB1B CDS strip with exons, the PTC, the last
    exon-exon junction, and the 748 nt NMD bracket.
F2 dosage_ladder: residual-BUBR1 levels from mouse genetics with outcomes.
F3 verdict_matrix: candidate class x grading axis with the per-axis result.
"""
import json
from pathlib import Path

import plotly.graph_objects as go
from silico_figures import apply_theme, save_figure_bundle

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
allele = json.loads(RESULTS.joinpath("allele_map.json").read_text())
cands = json.loads(RESULTS.joinpath("candidate_grades.json").read_text())
syn = json.loads(RESULTS.joinpath("recommendation_synthesis.json").read_text())

# ---------------------------------------------------------------- F1 -------
blocks = allele["cds_blocks"]
cds_len = allele["cds_length_nt"]
ptc = allele["alleles"]["p.Leu737Ter"]["cdna_position"]
miss = allele["alleles"]["p.Asn1002Lys"]["cdna_position"]
last_j = allele["nmd_rule"]["last_exon_junction_cds"]

fig = go.Figure()
y0, y1 = 0.0, 1.0
cum = 0
for i, b in enumerate(blocks):
    x0 = cum
    x1 = cum + b["cds_block_len"]
    cum = x1
    fig.add_shape(type="rect", x0=x0, x1=x1, y0=y0, y1=y1,
                  fillcolor="#2E5EAA" if b["rank"] % 2 else "#7FA6D9",
                  line=dict(color="white", width=1.5))
    if b["rank"] in (1, 17, 23):
        fig.add_annotation(x=(x0 + x1) / 2, y=-0.22, text=f"exon {b['rank']}",
                           showarrow=False, font=dict(size=12))
fig.add_shape(type="line", x0=ptc, x1=ptc, y0=-0.12, y1=1.12,
              line=dict(color="#D64545", width=3))
fig.add_annotation(x=ptc, y=1.24, text="PTC c.2210 (exon 17)", showarrow=False,
                   font=dict(size=12, color="#D64545"))
fig.add_shape(type="line", x0=last_j, x1=last_j, y0=-0.12, y1=1.12,
              line=dict(color="#1D7874", width=2, dash="dot"))
fig.add_annotation(x=last_j, y=1.24, text="last exon junction c.2957", showarrow=False,
                   font=dict(size=12, color="#1D7874"))
fig.add_shape(type="line", x0=miss, x1=miss, y0=-0.12, y1=1.12,
              line=dict(color="#8C564B", width=2))
fig.add_annotation(x=miss - 160, y=1.38, text="p.Asn1002Lys (exon 23, escapes NMD)",
                   showarrow=False, font=dict(size=12, color="#8C564B"))
# NMD bracket
fig.add_shape(type="line", x0=ptc, x1=last_j, y0=0.5, y1=0.5,
              line=dict(color="black", width=1.5, dash="dash"))
fig.add_annotation(x=(ptc + last_j) / 2, y=0.62, text="748 nucleotides to last junction (rule says about 50)",
                   showarrow=False, font=dict(size=12))
fig.update_xaxes(title_text="Coding-sequence position (nucleotides)", range=[0, cds_len + 30])
fig.update_yaxes(visible=False, range=[-0.45, 1.55])
fig.update_layout(title="BUB1B coding exons and the nonsense-mediated decay arithmetic")
apply_theme(fig, height=380)
save_figure_bundle(fig, "transcript_geometry", root=HERE.parent / "figures",
                   data={"cds_blocks": blocks, "ptc": ptc, "last_junction": last_j,
                          "missense": miss, "cds_length": cds_len},
                   alt="Strip diagram of 23 BUB1B coding exons. The premature stop in exon 17 sits 748 nucleotides before the last exon junction, far past the 50-nucleotide decay rule.")

# ---------------------------------------------------------------- F2 -------
rungs = [
    {"label": "about 10 percent (compound hypomorph mice)", "x": 10,
     "outcome": "severe aneuploidy, progeroid disease"},
    {"label": "about 50 percent (heterozygotes)", "x": 50,
     "outcome": "more carcinogen-driven tumors"},
    {"label": "above 100 percent (overexpression)", "x": 150,
     "outcome": "protected from aneuploidy-driven cancer"},
]
fig = go.Figure()
fig.add_trace(go.Bar(
    x=[r["label"] for r in rungs], y=[r["x"] for r in rungs],
    marker_color=["#D64545", "#E29A00", "#1D7874"],
    hovertext=[r["outcome"] for r in rungs], hoverinfo="text+y",
))
fig.add_shape(type="line", x0=-0.5, x1=2.5, y0=50, y1=50, line=dict(color="black", width=1.5, dash="dash"))
fig.add_annotation(x=0.15, y=56, text="heterozygote level, the minimal protective boundary",
                   showarrow=False, font=dict(size=12))
fig.add_annotation(x=1.9, y=6, text="readthrough drugs reach single digits",
                   showarrow=False, font=dict(size=12, color="#D64545"))
fig.update_xaxes(title_text="Mouse genetic state", type="category")
fig.update_yaxes(title_text="Residual BUBR1 level, percent of normal")
fig.update_layout(title="Dosage threshold for BUBR1 from mouse genetics")
apply_theme(fig, height=420)
save_figure_bundle(fig, "dosage_ladder", root=HERE.parent / "figures",
                   data={"rungs": rungs},
                   alt="Bar chart of residual BUBR1 levels in mice. Ten percent causes severe disease, half causes more tumors, and overexpression protects. Readthrough drugs reach only single digits.")

# ---------------------------------------------------------------- F3 -------
rows = []
for c in cands:
    if c["verdict"] == "comparator":
        continue
    rows.append({
        "class": c["class"].split(",")[0],
        "mechanism": c["mechanism_tie"],
        "regulatory": c["regulatory"].split(" (")[0],
        "pediatric": c["pediatric_experience"].split(" (")[0],
        "verdict": c["verdict"],
    })

def score(v):
    vl = v.lower()
    if any(k in vl for k in ("wrong-direction", "rejected", "none", "contraindicated", "lapsed", "harm")):
        return 0, "fail"
    if any(k in vl for k in ("strong", "direct", "approved-current", "labeled", "primary")):
        return 2, "pass"
    return 1, "mid"

axes = ["mechanism", "regulatory", "pediatric"]
z, text = [], []
for r in rows:
    zz, tt = [], []
    for a in axes:
        s, lab = score(r[a])
        zz.append(s)
        tt.append(r[a])
    zs, zlab = score(r["verdict"].replace("fallback-research", "mid").replace("primary", "pass").replace("rejected", "fail"))
    verdict_num = {"primary": 2, "fallback-research": 1, "rejected": 0}[r["verdict"]]
    zz.append(verdict_num)
    tt.append(r["verdict"])
    z.append(zz)
    text.append(tt)

colorscale = [[0, "#E8B4B0"], [0.5, "#F2E1B8"], [1, "#BADDC5"]]
fig = go.Figure(go.Heatmap(
    z=z, x=["Mechanism tie", "Regulatory", "Pediatric", "Verdict"],
    y=[r["class"] for r in rows], text=text, texttemplate="%{text}",
    colorscale=colorscale, zmin=0, zmax=2, showscale=False,
    xgap=3, ygap=3,
))
fig.update_yaxes(autorange="reversed", title_text="Candidate class")
fig.update_xaxes(side="top", title_text="Grading axis")
fig.update_layout(title="Candidate grading summary. Green passes, amber partial, red fails")
apply_theme(fig, height=60 * len(rows) + 140)
save_figure_bundle(fig, "verdict_matrix", root=HERE.parent / "figures",
                   data={"rows": rows, "z": z, "text": text},
                   alt="Matrix of candidate class against mechanism, regulatory and pediatric axes plus verdict. Only metformin passes all axes.")

print("OK figures written")
