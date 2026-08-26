#!/usr/bin/env python3
"""Reference check for the GPN-MSA scoring stack.

Runs the maintained gpn msa forward pass on the published chr6 test fixture and
compares the masked-position logits to the shipped baseline (rtol=atol=1e-4).
Never loads the whole-genome store: this validates model weights, auto-class
registration, and numeric behavior only.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import torch

DNA_VOCAB = "-ACGT?"
MASK_TOKEN_ID = DNA_VOCAB.index("?")
NUCLEOTIDE_INDICES = [DNA_VOCAB.index(n) for n in "ACGT"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assets", type=Path, required=True)
    ap.add_argument("--revision", type=str, default=None)
    ap.add_argument(
        "--fixture-url-base",
        type=str,
        default="https://raw.githubusercontent.com/songlab-cal/gpn/690557d949309cf4f4234554888bb5421c49aede/tests/fixtures",
    )
    args = ap.parse_args()

    def _asset(name: str) -> Path:
        p = args.assets / name
        if not p.exists():
            import urllib.request

            p.parent.mkdir(parents=True, exist_ok=True)
            urllib.request.urlretrieve(f"{args.fixture_url_base}/{name}", p)
        return p

    baseline = json.loads(_asset("published_model_baseline.json").read_text())
    record = baseline["models"]["gpn_msa"]
    revision = args.revision or record["revision"]
    expected = record["expected"]

    from gpn import register_auto_classes
    from transformers import AutoModelForMaskedLM

    register_auto_classes("msa")
    fixture = _asset("hg38_chr6_31575665_31575793_multiz100way.npz")
    with open(fixture, "rb") as fh, np.load(fh) as aln:
        tokens = aln["gpn_msa_tokens"].astype(np.int64)
    msa = torch.from_numpy(tokens).unsqueeze(0)
    input_ids = msa[:, :, 0].clone()
    aux_features = msa[:, :, 1:]
    position = expected["input"]["masked_position_zero_based"]
    input_ids[0, position] = MASK_TOKEN_ID

    model = AutoModelForMaskedLM.from_pretrained(record["model_id"], revision=revision).eval()
    with torch.inference_mode():
        out = model(input_ids=input_ids, aux_features=aux_features).logits
    actual = out[0, position, NUCLEOTIDE_INDICES].detach().cpu().to(torch.float32)
    exp = torch.tensor(expected["logits"], dtype=torch.float32)
    torch.testing.assert_close(actual, exp, rtol=1e-4, atol=1e-4)
    print("GPN_MSA_BASELINE_OK", model.config.model_type, "n_params", sum(p.numel() for p in model.parameters()))
    print("logits", actual.tolist())


if __name__ == "__main__":
    main()
