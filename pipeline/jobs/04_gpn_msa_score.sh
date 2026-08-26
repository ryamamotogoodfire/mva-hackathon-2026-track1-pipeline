#!/bin/bash
# MVA Track 1 step 4: score surviving candidate variants (plus matched inert
# controls) with GPN-MSA through the maintained upstream CLI, on one GPU.
#
# Environment (see README):
#   DATA_DIR      pipeline working directory with mining/score_inputs.parquet
#   PIPELINE_DIR  this repository's pipeline/ directory (default: script dir)
#   GPN installed from the public repository pinned to
#     git+https://github.com/songlab-cal/gpn.git@690557d949309cf4f4234554888bb5421c49aede
#   plus transformers==5.15.0 safetensors==0.8.0 jaxtyping==0.3.11 cyclopts "zarr<3"
#   and a CUDA build of torch. The alignment store streams over https.
#
# Inputs  : $DATA_DIR/mva-track1/mining/score_inputs.parquet (chrom, pos, ref, alt, role)
# Outputs : $DATA_DIR/mva-track1/scoring/gpn_msa_scored_variants.parquet + metadata
set -euo pipefail

: "${DATA_DIR:?set DATA_DIR to the pipeline working directory}"
PIPELINE_DIR="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

ART="$DATA_DIR/mva-track1"
IN="$ART/mining/score_inputs.parquet"
OUTD="$ART/scoring"
mkdir -p "$OUTD"

cd "$PIPELINE_DIR"

# sort by locus for chunk-local reads of the streamed alignment store
python - "$IN" "$OUTD/score_inputs_sorted.parquet" <<'EOF'
import pandas as pd, sys
df = pd.read_parquet(sys.argv[1])
def chrom_key(c):
    c = str(c)
    return int(c) if c.isdigit() else (23 if c == "X" else 24 if c == "Y" else 25)
df["ck"] = df["chrom"].map(chrom_key)
df = df.sort_values(["ck", "pos"]).drop(columns="ck").reset_index(drop=True)
df.to_parquet(sys.argv[2])
print(f"sorted {len(df)} variants; roles: {df['role'].value_counts().to_dict() if 'role' in df else 'n/a'}")
EOF

STORE_URL="zip:///::https://huggingface.co/datasets/songlab/multiz100way/resolve/main/89.zarr.zip"

echo "[04] scoring with gpn msa vep"
python -m gpn.cli msa vep \
  --input-path "$OUTD/score_inputs_sorted.parquet" \
  --msa-path "$STORE_URL" \
  --window-size 512 \
  --model-path songlab/gpn-msa-sapiens \
  --output-path "$OUTD/gpn_msa_scores.parquet" \
  --checkpoint-batch-size 500 \
  --per-device-eval-batch-size 32

python - "$OUTD" "$OUTD/score_inputs_sorted.parquet" <<'EOF'
import json, sys
import pandas as pd
outd, in_path = sys.argv[1], sys.argv[2]
scores = pd.read_parquet(f"{outd}/gpn_msa_scores.parquet")
inp = pd.read_parquet(f"{outd}/score_inputs_sorted.parquet")
meta = {
    "n_input": int(len(inp)),
    "n_scored": int(len(scores)),
    "score_summary": {k: float(v) for k, v in scores["score"].describe().items()},
}
# join back side-by-side (CLI preserves input order)
inp["gpn_msa_score"] = scores["score"].values
for role, g in inp.groupby("role"):
    meta[f"n_{role}"] = int(len(g))
    meta[f"mean_{role}"] = float(g["gpn_msa_score"].mean())
inp.drop(columns=["score"], errors="ignore").to_parquet(f"{outd}/gpn_msa_scored_variants.parquet")
json.dump(meta, open(f"{outd}/scoring_meta.json", "w"), indent=1)
print(json.dumps(meta, indent=1))
n_in, n_out = len(inp), len(scores)
assert n_in == n_out, f"row mismatch: {n_in} in vs {n_out} out"
print("[04] SCORING_META_OK")
EOF
echo "[04] DONE"
