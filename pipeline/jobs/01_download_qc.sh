#!/bin/bash
# MVA Track 1 step 1: download the gated proband VCF + index + phenotype doc
# from the Hugging Face dataset repo into $DATA_DIR, then run VCF QC.
#
# Required environment:
#   DATA_DIR  Local working directory for all pipeline inputs and outputs.
#   HF_TOKEN  Hugging Face access token accepted on the gated dataset
#             SageBio/mva-hackathon-2026-data.
set -euo pipefail

: "${DATA_DIR:?set DATA_DIR to a local working directory}"
: "${HF_TOKEN:?set HF_TOKEN to a Hugging Face token accepted on the gated dataset}"

OUT="$DATA_DIR/mva-track1/vcf"
mkdir -p "$OUT"
BASE="https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data/resolve/main"

for f in WGS_EX2312012_HGWCNDSX7.vcf.gz WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi Challenge_Clinical_Phenotype_1.docx README.md; do
  echo "[01] downloading $f"
  curl -fSL --retry 5 --retry-all-errors -C - \
       -H "Authorization: Bearer ${HF_TOKEN}" \
       -o "$OUT/$f" "$BASE/$f"
done

echo "[01] sizes:"
ls -la "$OUT"

gzip -t "$OUT/WGS_EX2312012_HGWCNDSX7.vcf.gz" && echo "[01] gzip integrity OK"

python3 - <<'EOF'
import gzip, json, re, os

vcf = os.path.join(os.environ["DATA_DIR"],
                   "mva-track1", "vcf", "WGS_EX2312012_HGWCNDSX7.vcf.gz")

qc = {"file": vcf}
header_meta = {}
contigs = []
samples = []
n_records = 0
filter_counts = {}
n_multiallelic = 0
n_star = 0
n_snv = 0
n_indel = 0
n_symbolic = 0
chrom_counts = {}
first_records = []

with gzip.open(vcf, "rt") as fh:
    for line in fh:
        if line.startswith("##"):
            m = re.match(r"##([^=]+)=(.*)", line.strip())
            if m:
                k, v = m.group(1), m.group(2)
                if k == "contig":
                    contigs.append(v.strip("<>"))
                elif k in ("reference", "fileformat", "source", "ALT", "phasing"):
                    header_meta.setdefault(k, []).append(v)
                else:
                    header_meta.setdefault(k, []).append(v)
            continue
        if line.startswith("#"):
            samples = line.strip().split("\t")[9:]
            continue
        n_records += 1
        f = line.split("\t", 7)
        chrom, pos, vid, ref, alt, qual, filt = f[:7]
        chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        filter_counts[filt] = filter_counts.get(filt, 0) + 1
        if len(first_records) < 3:
            first_records.append(line.strip())
        if alt == "*":
            n_star += 1
            continue
        alts = alt.split(",")
        if len(alts) > 1:
            n_multiallelic += 1
        if any(a.startswith("<") for a in alts):
            n_symbolic += 1
        elif len(ref) == 1 and all(len(a) == 1 for a in alts):
            n_snv += 1
        else:
            n_indel += 1

qc["n_records"] = n_records
qc["n_contigs_header"] = len(contigs)
qc["contigs"] = contigs
qc["samples"] = samples
qc["filter_counts"] = filter_counts
qc["n_multiallelic"] = n_multiallelic
qc["n_star_allele"] = n_star
qc["n_snv_like"] = n_snv
qc["n_indel_like"] = n_indel
qc["n_symbolic"] = n_symbolic
qc["header_reference"] = header_meta.get("reference")
qc["header_fileformat"] = header_meta.get("fileformat")
qc["chrom_prefix_chr"] = all(c.startswith("chr") for c in chrom_counts) if chrom_counts else None
qc["chrom_name_example"] = sorted(chrom_counts)[:5]
qc["per_chrom_counts"] = {k: chrom_counts[k] for k in sorted(chrom_counts)}
qc["first_records"] = first_records

out = os.path.join(os.environ["DATA_DIR"], "mva-track1", "qc_vcf.json")
with open(out, "w") as fh:
    json.dump(qc, fh, indent=2)
print(json.dumps({k: (v if k not in ("contigs", "first_records", "per_chrom_counts") else "...") for k, v in qc.items()}, indent=2))
print("per-chrom sample:", dict(sorted(chrom_counts.items())[:6]))
print("[01] QC written to", out)
EOF
echo "[01] DONE"
