#!/usr/bin/env python3
"""Track 1 submission CSV writer.

Takes a ranked candidate table (one row per proposed finding; compound-het
pairs carry both variants) and emits the exact challenge submission schema:
  proband_id, chrom_1, pos_1, ref_1, alt_1, chrom_2, pos_2, ref_2, alt_2,
  epcr, finding_type, notes
Coordinates must be GRCh38 with 'chr'-prefixed chromosomes per the template
example. Enforces: at most 10 rows, epcr in (0, 1] and non-increasing with
rank order.
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

COLUMNS = [
    "proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
    "chrom_2", "pos_2", "ref_2", "alt_2",
    "epcr", "finding_type", "notes",
]


def format_row(row, proband: str) -> dict:
    def c(cstr):
        c = str(cstr)
        return c if c.startswith("chr") else "chr" + c

    def b(x):
        return "" if pd.isna(x) or x in ("", "-") else x

    out = {
        "proband_id": proband,
        "chrom_1": c(row["chrom1"]),
        "pos_1": int(row["pos1"]),
        "ref_1": row["ref1"],
        "alt_1": row["alt1"],
        "chrom_2": "" if b(row.get("chrom2", "")) == "" else c(row["chrom2"]),
        "pos_2": b(row.get("pos2", "")),
        "ref_2": b(row.get("ref2", "")),
        "alt_2": b(row.get("alt2", "")),
        "epcr": row["epcr"],
        "finding_type": row["finding_type"],
        "notes": row.get("notes", "") or "",
    }
    return out


def validate(df: pd.DataFrame) -> None:
    errs = []
    if not (1 <= len(df) <= 10):
        errs.append(f"row count {len(df)} outside 1..10")
    if not ((df["epcr"] > 0) & (df["epcr"] <= 1)).all():
        errs.append("epcr outside (0, 1]")
    if (df["epcr"].diff().fillna(0) > 1e-12).any():
        errs.append("epcr not non-increasing with rank")
    if not df["finding_type"].isin(["primary", "secondary"]).all():
        errs.append("finding_type must be primary or secondary")
    for i, r in df.iterrows():
        alts2 = r["alt_2"]
        fields2 = [r["chrom_2"], r["pos_2"], r["ref_2"], alts2]
        if (alts2 == "") and any(str(x) != "" for x in fields2):
            errs.append(f"row {i}: partial second-variant fields")
        if r["chrom_2"] == "" and alts2 != "":
            errs.append(f"row {i}: second variant alt without chrom")
    if errs:
        raise SystemExit("submission validation failed:\n" + "\n".join(errs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ranked_parquet", help="ranked candidates (paired columns chrom1..alt2)")
    ap.add_argument("out_csv")
    ap.add_argument("--proband", default="PROBAND01")
    ap.add_argument("--max-rows", type=int, default=10)
    args = ap.parse_args()

    df = pd.read_parquet(args.ranked_parquet).head(args.max_rows)
    if (df["alt1"] == "").any():
        raise SystemExit("empty alt1 in ranked table")
    out = pd.DataFrame([format_row(r, args.proband) for _, r in df.iterrows()])[COLUMNS]
    validate(out)
    out.to_csv(args.out_csv, index=False)
    print(f"SUBMISSION_OK rows={len(out)} -> {args.out_csv}")
    print(out.to_string())


if __name__ == "__main__":
    main()
