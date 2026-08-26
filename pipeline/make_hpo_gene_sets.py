#!/usr/bin/env python3
"""Build the HPO term-to-gene index that step 3 (mine candidates) consumes.

The index is derived entirely from the public HPO annotation file
genes_to_phenotype.txt. The only case-specific input is the list of phenotype
term IDs, which the reader extracts from the gated challenge phenotype document
(see README). No case data ships with this repository.

Usage:
  python3 make_hpo_gene_sets.py --terms hpo_terms.txt
    --genes-to-phenotype genes_to_phenotype.txt      # downloaded if missing
    --out "$DATA_DIR/mva-track1/hpo/hpo_gene_sets.json"

hpo_terms.txt holds one HPO identifier per line, for example HP:0004322.
Blank lines and lines starting with # are ignored.
"""
import argparse
import json
import os
import urllib.request

HPO_URL = "https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--terms", required=True, help="text file, one HP:XXXXXXX id per line")
    ap.add_argument("--genes-to-phenotype", default="genes_to_phenotype.txt")
    ap.add_argument("--out", default=os.path.join(os.environ.get("DATA_DIR", "."),
                                                  "mva-track1", "hpo", "hpo_gene_sets.json"))
    args = ap.parse_args()

    terms = []
    for line in open(args.terms):
        line = line.strip()
        if line and not line.startswith("#"):
            terms.append(line)
    if not terms:
        raise SystemExit("no HPO terms found in --terms file")

    g2p = args.genes_to_phenotype
    if not os.path.exists(g2p):
        print(f"[hpo] downloading {HPO_URL}")
        urllib.request.urlretrieve(HPO_URL, g2p)

    wanted = set(terms)
    names = {}
    exact = {t: set() for t in terms}
    with open(g2p) as fh:
        for line in fh:
            if line.startswith("ncbi_gene_id") or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            symbol, hpo_id, hpo_name = parts[1], parts[2], parts[3]
            if hpo_id in wanted:
                exact[hpo_id].add(symbol)
                names[hpo_id] = hpo_name

    missing = [t for t in terms if not exact[t]]
    if missing:
        raise SystemExit(f"no genes_to_phenotype rows for terms: {missing}")

    out = {
        "names": {t: names[t] for t in terms},
        "exact": {t: sorted(exact[t]) for t in terms},
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=False)
    total = sum(len(v) for v in exact.values())
    print(f"[hpo] wrote {args.out}: {len(terms)} terms, {total} gene mappings")


if __name__ == "__main__":
    main()
