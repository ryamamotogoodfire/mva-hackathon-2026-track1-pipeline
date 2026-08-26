#!/usr/bin/env python3
"""Validate the streaming whole-genome alignment store read path.

Pulls the exact chr6 test-fixture window from the remote multiz100way 89-species
store (zip over https), tokenizes it, and compares element-wise against the
shipped fixture tokens. A mismatch would mean wrong species order or wrong store;
equality means candidate windows can be read directly from the network without
materializing the 36 GB genome-wide store.
"""
import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

FIXTURE_URL = "https://raw.githubusercontent.com/songlab-cal/gpn/690557d949309cf4f4234554888bb5421c49aede/tests/fixtures/hg38_chr6_31575665_31575793_multiz100way.npz"
BASELINE_URL = "https://raw.githubusercontent.com/songlab-cal/gpn/690557d949309cf4f4234554888bb5421c49aede/tests/fixtures/published_model_baseline.json"
STORE_URL = "zip:///::https://huggingface.co/datasets/songlab/multiz100way/resolve/main/89.zarr.zip"


def main():
    from gpn.data import Tokenizer
    from gpn.msa.data import GenomeMSA

    t0 = time.time()
    gm = GenomeMSA(STORE_URL)
    open_dt = time.time() - t0
    print(f"store opened in {open_dt:.1f}s; chroms: {list(gm.data.keys())[:30]}")

    t0 = time.time()
    arr = gm.get_msa("6", 31575665, 31575793)
    slice_dt = time.time() - t0
    print(f"one 128bp window read in {slice_dt:.2f}s; shape {arr.shape} dtype {arr.dtype}")

    tok = Tokenizer()(arr)
    fix = urllib.request.urlretrieve(FIXTURE_URL, "/tmp/fixture.npz")
    with open(fix[0] if isinstance(fix, tuple) else fix, "rb") as fh:
        with np.load(fh) as aln:
            fixture_tokens = aln["gpn_msa_tokens"]
    print("fixture shape", fixture_tokens.shape)
    same = tok.shape == fixture_tokens.shape and (tok == fixture_tokens).all()
    print("STORE_MATCHES_FIXTURE:", same)

    # timing probe: 20 successive reads around the same area
    t0 = time.time()
    for i in range(20):
        _ = gm.get_msa("6", 31575665 + i * 512, 31575793 + i * 512)
    dt = time.time() - t0
    print(f"20 streaming reads took {dt:.1f}s ({dt/20:.2f}s/window)")

    base = json.load(urllib.request.urlopen(BASELINE_URL))
    species = base["alignment_fixture"]["gpn_msa_species"]
    print("n species expected:", len(species), species[:3])


if __name__ == "__main__":
    main()
