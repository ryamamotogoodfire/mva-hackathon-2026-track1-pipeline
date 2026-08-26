#!/usr/bin/env python3
"""Small end-to-end check of the maintained `gpn msa vep` CLI against the
streaming alignment store + GPU, using synthetic SNVs whose REF alleles are
taken from the store's own human row (so the REF check can only fail on real
bugs). Prints per-variant scores and wall time."""
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

STORE_URL = "zip:///::https://huggingface.co/datasets/songlab/multiz100way/resolve/main/89.zarr.zip"
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/gpn_cli_check")
OUT.mkdir(parents=True, exist_ok=True)


def synth_variants():
    from gpn.msa.data import GenomeMSA

    gm = GenomeMSA(STORE_URL)
    bases = np.array(
        [
            gm.get_msa("6", p - 1, p)[0, 0].decode()
            for p in range(31575600, 31575600 + 40)
        ],
        dtype="U1",
    )
    rows = []
    others = {"A": ["C", "G", "T"], "C": ["A", "G", "T"], "G": ["A", "C", "T"], "T": ["A", "C", "G"]}
    for i, b in enumerate(bases):
        if b.upper() in "ACGT":
            rows.append(
                dict(chrom="6", pos=31575600 + i, ref=b.upper(), alt=others[b.upper()][i % 3])
            )
    # a few on chr2 and chrX to test chunk switching + sex chromosomes
    for chrom, pos in [("2", 71536710), ("2", 241695590), ("X", 153296333), ("X", 33043000)]:
        b = gm.get_msa(chrom, pos - 1, pos)[0, 0].decode().upper()
        if b in "ACGT":
            rows.append(dict(chrom=chrom, pos=pos, ref=b, alt=others[b][0]))
    df = pd.DataFrame(rows)
    df.to_parquet(OUT / "variants.parquet")
    print(f"{len(df)} synthetic variants written")
    return df


def main():
    df = synth_variants()
    env = dict(os.environ)
    cmd = [
        sys.executable, "-m", "gpn.cli", "msa", "vep",
        "--input-path", str(OUT / "variants.parquet"),
        "--msa-path", STORE_URL,
        "--window-size", "512",
        "--model-path", "songlab/gpn-msa-sapiens",
        "--output-path", str(OUT / "scores.parquet"),
        "--per-device-eval-batch-size", "16",
    ]
    t0 = time.time()
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    dt = time.time() - t0
    print(r.stdout[-2000:])
    print(r.stderr[-2000:])
    r.check_returncode()
    print(f"CLI wall time for {len(df)} variants: {dt:.1f}s")
    out = pd.read_parquet(OUT / "scores.parquet")
    print(out.describe())
    print(out.head(10).to_string())
    print("CLI_CHECK_OK")


if __name__ == "__main__":
    main()
