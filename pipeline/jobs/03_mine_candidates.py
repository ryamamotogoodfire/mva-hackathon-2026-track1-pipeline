#!/usr/bin/env python3
"""MVA Track 1 step 3: mine candidate causal variants from the VEP annotation.

  - rsID index stores every ALT key per rsID and the join selects the exact
    allele match (no sibling-ALT collapse at multiallelic records)
  - lane SNV filters exclude alt/ref '-' explicitly
  - compound-het evidence pool spans primary + loose AF tiers (tier recorded)
  - candidates require FILTER PASS and execution on canonical chromosomes
    (1..22, X, Y, M excluded deliberately: mtDNA not in plan scope);
    unobserved population frequency flagged explicitly ("no population data")
    rather than promoted as extreme rarity
  - join failures are written to join_unresolved.tsv, with a count of how
    many carried any coding consequence

Inputs (all overridable; defaults resolve under $DATA_DIR/mva-track1):
  --vep-tab  VEP tab output from step 2 (default vep_out/proband_vep.tab.bgz)
  --vcf      proband VCF from step 1
  --hpo      HPO gene-set JSON from make_hpo_gene_sets.py
             (default $DATA_DIR/mva-track1/hpo/hpo_gene_sets.json)
  --out      output directory (default $DATA_DIR/mva-track1/mining)
"""
import argparse
import gzip
import json
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd

AP = argparse.ArgumentParser()
AP.add_argument("--vep-tab", default=None)
AP.add_argument("--vcf", default=None)
AP.add_argument("--out", default=None)
AP.add_argument("--hpo", default=None)
ARGS = AP.parse_args()

if "DATA_DIR" not in os.environ:
    raise SystemExit("set DATA_DIR to the pipeline working directory")
ART = os.path.join(os.environ["DATA_DIR"], "mva-track1")
VEP_TAB = ARGS.vep_tab or os.path.join(ART, "vep_out", "proband_vep.tab.bgz")
VCF = ARGS.vcf or os.path.join(ART, "vcf", "WGS_EX2312012_HGWCNDSX7.vcf.gz")
OUT = ARGS.out or os.path.join(ART, "mining")
os.makedirs(OUT, exist_ok=True)

HPO_JSON = ARGS.hpo or os.path.join(ART, "hpo", "hpo_gene_sets.json")
HPO_SETS = json.load(open(HPO_JSON))
HPO_EXACT = HPO_SETS["exact"]
HPO_GENE_UNION = set().union(*[set(v) for v in HPO_EXACT.values()]) if HPO_EXACT else set()

MVA_PANEL = ["BUB1B", "CEP57", "TRIP13"]

PRIMARY_AF = 0.001
LOOSE_AF = 0.01

HIGH_CONSEQUENCES = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
}
CODING_CONSEQUENCES = HIGH_CONSEQUENCES | {
    "missense_variant", "splice_region_variant", "inframe_deletion", "inframe_insertion",
    "initiator_codon_variant", "start_retained_variant", "stop_retained_variant",
    "protein_altering_variant", "coding_sequence_variant", "splice_polypyrimidine_tract_variant",
    "splice_donor_5th_base_variant", "splice_donor_region_variant",
}
SPLICE_ONLY = {c for c in CODING_CONSEQUENCES if c.startswith("splice_")}

AF_ALLPOP = ["AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "gnomADe_AF", "gnomADg_AF"]

DAMAGING = ["ptv", "missense_consensus", "splice_rare", "coding_other"]
QUALIFYING = DAMAGING + ["missense_any"]

CANONICAL_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}


def minrep(pos, ref, alt):
    while ref and alt and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    return pos, ref or "-", alt or "-"


def ukey(chrom, pos, ref, alt):
    p, r, a = minrep(pos, ref, alt)
    return f"{chrom}_{p}_{r}/{a}"


def key_parts(k):
    """chrom, pos, ref, alt from a ukey, robust to '_' inside contig names."""
    i = k.rfind("/")
    alt = k[i + 1:]
    pre = k[:i]
    j = pre.rfind("_")
    ref = pre[j + 1:]
    cp = pre[:j]
    j2 = cp.rfind("_")
    return cp[:j2], int(cp[j2 + 1:]), ref, alt


def parse_vcf_genotypes():
    gtmap = {}
    rsindex = defaultdict(list)
    with gzip.open(VCF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 10)
            if len(f) < 10:
                continue
            chrom, pos, vid, ref, alts, qual, filt, info, fmt, sample = f[:10]
            fmtk = fmt.split(":")
            rec = dict(zip(fmtk, sample.split(":")))
            gt = rec.get("GT", "")
            if gt in ("0/0", "0|0", "./."):
                continue
            ad = rec.get("AD", "")
            dp = rec.get("DP", "")
            ads = ad.split(",") if ad else []
            def toi(i0, k):
                try:
                    return int(i0[k])
                except Exception:
                    return -1
            for i, alt in enumerate(alts.split(","), start=1):
                gt_norm = gt.replace("|", "/")
                if gt_norm.count(str(i)) == 0:
                    continue
                k = ukey(chrom, int(pos), ref, alt)
                gtmap[k] = dict(
                    zygosity="hom_alt" if gt_norm.count(str(i)) >= 2 else "het",
                    gt=gt_norm, ad_ref=toi(ads, 0), ad_alt=toi(ads, i),
                    dp=int(dp) if dp.isdigit() else -1, vcf_filter=filt)
                if vid and vid != ".":
                    if k not in rsindex[vid]:
                        rsindex[vid].append(k)
    return gtmap, rsindex


def to_float(x):
    try:
        return float(x)
    except Exception:
        return np.nan


def strrep(x):
    return "-" if x in (None, "", "-") else str(x)


def parse():
    print("[03] parsing VCF genotypes", flush=True)
    gt, rsindex = parse_vcf_genotypes()
    n_rs_multi = sum(1 for v in rsindex.values() if len(v) > 1)
    print(f"[03] {len(gt)} carrying alleles; rs-index {len(rsindex)}; rs multi-key {n_rs_multi}")

    variant = {}
    rows_seen = [0]
    join_stats = {"positional_exact": 0, "rs_exactalt": 0, "rs_ambiguous": 0,
                  "location_fallback": 0, "unresolved": 0}
    unresolved_rows = []

    def collapse(rows):
        def rtag(row):
            if row.get("MANE_SELECT") not in ("-", "", None):
                return 0
            if row.get("CANONICAL") == "YES":
                return 1
            return 2
        rows = sorted(rows, key=rtag)
        rep = rows[0]
        cons = sorted({c for r in rows for c in r["Consequence"].split(",")})
        lof_hc = any(r.get("LoF") == "HC" for r in rows)
        lof_any = any(r.get("LoF") in ("HC", "LC") for r in rows)
        votes = 0
        if "deleterious" in str(rep.get("SIFT", "-")):
            votes += 1
        if "probably_damaging" in str(rep.get("PolyPhen", "-")):
            votes += 1
        if rep.get("am_class") == "likely_pathogenic":
            votes += 1
        afs = [to_float(rep.get(c)) for c in AF_ALLPOP]
        afs = [v for v in afs if v == v]
        max_af = max(afs) if afs else 0.0
        max_af_vep = to_float(rep.get("MAX_AF"))
        if max_af_vep != max_af_vep:
            max_af_vep = 0.0
        return dict(
            gene=strrep(rep.get("SYMBOL")), gene_id=strrep(rep.get("Gene")),
            feature=strrep(rep.get("Feature")), canonical=strrep(rep.get("CANONICAL")),
            mane=strrep(rep.get("MANE_SELECT")),
            consequence=",".join(cons), impact=rep.get("IMPACT", "-"),
            hgvsc=strrep(rep.get("HGVSc")), hgvsp=strrep(rep.get("HGVSp")),
            max_af_contemp=max_af, max_af_vep=max_af_vep,
            af_observed=bool(afs),
            clin_sig=strrep(rep.get("ClinVar_CLNSIG")), clin_sig_cache=strrep(rep.get("CLIN_SIG")),
            existing=strrep(rep.get("Existing_variation")),
            lof_hc=lof_hc, lof_any=lof_any,
            lof_filter=strrep(rep.get("LoF_filter")), lof_flags=strrep(rep.get("LoF_flags")),
            am_class=strrep(rep.get("am_class")), am_score=to_float(rep.get("am_pathogenicity")),
            sift=strrep(rep.get("SIFT")), polyphen=strrep(rep.get("PolyPhen")),
            consensus_votes=votes,
        )

    pending = []
    last_key = [None]
    consequences_hist = defaultdict(int)

    def flush():
        if last_key[0] is not None and pending:
            variant[last_key[0]] = collapse(pending)
        pending.clear()

    def resolve(row):
        vk_raw = row["Uploaded_variation"]
        allele = row.get("Allele", "")
        if vk_raw in gt:
            return vk_raw, "positional_exact"
        if vk_raw in rsindex:
            cands = rsindex[vk_raw]
            matched = [k for k in cands if key_parts(k)[3] == allele]
            if len(matched) == 1:
                return matched[0], "rs_exactalt"
            if len(matched) > 1:
                join_stats["rs_ambiguous"] += 1
                return matched[0], "rs_ambiguous"
        loc = row.get("Location", "")
        m = re.match(r"^(\S+?):(\d+)(?:-(\d+))?$", loc)
        if m:
            chrom, p1, p2 = m.group(1), int(m.group(2)), int(m.group(3) or m.group(2))
            for pp in sorted({p1, p2, p1 - 1, p1 + 1}):
                if pp < 1:
                    continue
                for kk in gt_by_pos.get((chrom, pp), ()):
                    if key_parts(kk)[3] == allele:
                        return kk, "location_fallback"
        return None, "unresolved"

    global gt_by_pos
    gt_by_pos = defaultdict(list)
    for k in gt:
        c, p, r, a = key_parts(k)
        gt_by_pos[(c, p)].append(k)

    open_fn = open
    with open(VEP_TAB, "rb") as probe:
        if probe.read(2) == b"\x1f\x8b":
            open_fn = gzip.open
    with open_fn(VEP_TAB, "rt") as fh:
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            if header is None and line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if header is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            rows_seen[0] += 1
            row = dict(zip(header, parts))
            if row.get("Feature_type") != "Transcript":
                continue
            vk, style = resolve(row)
            join_stats[style] += 1
            if vk is None:
                unresolved_rows.append({
                    "uploaded_variation": row["Uploaded_variation"],
                    "location": row.get("Location", "-"),
                    "allele": row.get("Allele", "-"),
                    "consequence": row.get("Consequence", "-"),
                    "symbol": row.get("SYMBOL", "-"),
                })
                flush()
                last_key[0] = None
                continue
            if vk != last_key[0]:
                flush()
                last_key[0] = vk
            for c in row["Consequence"].split(","):
                consequences_hist[c] += 1
            pending.append(row)
        flush()

    if unresolved_rows:
        pd.DataFrame(unresolved_rows).to_csv(os.path.join(OUT, "join_unresolved.tsv"),
                                             sep="\t", index=False)

    joined = 0
    for k, rec in variant.items():
        g = gt.get(k)
        rec.update(g if g else dict(zygosity="unknown", gt="-/-", ad_ref=-1, ad_alt=-1, dp=-1, vcf_filter="-"))
        rec["has_gt"] = g is not None
        joined += g is not None
    print(f"[03] rows {rows_seen[0]}; collapsed {len(variant)}; joined {joined}; join styles {join_stats}")
    return variant, rows_seen[0], joined, join_stats, n_rs_multi, unresolved_rows, consequences_hist


gt_by_pos = {}


def classify(rec):
    cons = rec["consequence"].split(",")
    rec["is_ptv"] = any(c in HIGH_CONSEQUENCES for c in cons) or rec["lof_hc"]
    rec["is_coding"] = any(c in CODING_CONSEQUENCES for c in cons)
    rec["is_splice_only"] = rec["is_coding"] and not rec["is_ptv"] and all(c in SPLICE_ONLY for c in cons)
    rec["is_missense"] = "missense_variant" in cons
    rec["missense_damaging"] = rec["is_missense"] and rec["consensus_votes"] >= 2
    rec["missense_any"] = rec["is_missense"] and rec["consensus_votes"] >= 1
    rec["max_af"] = max(rec["max_af_contemp"], rec["max_af_vep"])
    rec["af_tier"] = "primary" if rec["max_af"] <= PRIMARY_AF else ("loose" if rec["max_af"] <= LOOSE_AF else "common")
    if not rec["af_observed"] and rec["max_af"] == 0.0:
        rec["af_label"] = "no population data"
    else:
        rec["af_label"] = f"{rec['max_af']:.2e}"
    rec["hpo_terms"] = sorted([t for t, gs in HPO_EXACT.items() if rec["gene"] in gs])
    rec["in_hpo"] = rec["gene"] in HPO_GENE_UNION
    rec["panel"] = rec["gene"] in MVA_PANEL
    if rec["af_tier"] == "common":
        rec["dmg_class"] = "common"
    elif rec["is_ptv"]:
        rec["dmg_class"] = "ptv"
    elif rec["missense_damaging"]:
        rec["dmg_class"] = "missense_consensus"
    elif rec["is_splice_only"] and rec["af_tier"] == "primary":
        rec["dmg_class"] = "splice_rare"
    elif rec["is_coding"] and rec["af_tier"] == "primary":
        rec["dmg_class"] = "coding_other"
    else:
        rec["dmg_class"] = "low"


def main():
    variant, rows_seen, joined, join_stats, n_rs_multi, unresolved_rows, consequences_hist = parse()
    n_unres = len(unresolved_rows)
    n_unres_coding = sum(any(c in CODING_CONSEQUENCES for c in r["consequence"].split(",")) for r in unresolved_rows)

    for k, rec in variant.items():
        classify(rec)

    decomp = pd.DataFrame([(k, *key_parts(k)) for k in variant],
                          columns=["ukey", "chrom", "pos", "ref", "alt"])
    df = decomp.merge(pd.DataFrame(list(variant.values())).assign(ukey=list(variant)), on="ukey")
    df.to_parquet(os.path.join(OUT, "variants_all_classes.parquet"), index=False)

    # quality/board masks
    df["pass_vcf"] = df["vcf_filter"] == "PASS"
    df["canonical"] = df["chrom"].isin(CANONICAL_CHROMS)
    in_board = df["pass_vcf"] & df["canonical"]

    het = df[in_board & (df["zygosity"] == "het") & df["dmg_class"].isin(QUALIFYING)
             & df["af_tier"].isin(["primary", "loose"])]
    homa = df[in_board & (df["zygosity"] == "hom_alt") & df["dmg_class"].isin(QUALIFYING)
              & df["af_tier"].isin(["primary", "loose"])]

    genes = defaultdict(lambda: {"het": [], "hom": []})
    for g, ser in het.groupby("gene")["ukey"]:
        if g not in ("-", ""):
            genes[g]["het"] = list(ser)
    for g, ser in homa.groupby("gene")["ukey"]:
        if g not in ("-", ""):
            genes[g]["hom"] = list(ser)

    gsum = []
    for gene, mem in genes.items():
        hpo_terms = sorted([t for t, g in HPO_EXACT.items() if gene in g])
        gsum.append(dict(
            gene=gene, n_qual_het=len(mem["het"]), n_qual_hom=len(mem["hom"]),
            compound_het=len(mem["het"]) >= 2,
            panel=gene in MVA_PANEL,
            n_hpo_terms=len(hpo_terms), hpo_terms=",".join(hpo_terms),
            members="|".join(mem["het"] + mem["hom"]),
        ))
    gsf = pd.DataFrame(gsum) if gsum else pd.DataFrame(columns=["gene"])
    if len(gsf):
        gsf = gsf.sort_values(
            ["panel", "n_hpo_terms", "compound_het", "n_qual_het"],
            ascending=[False, False, False, False])
    gsf.to_parquet(os.path.join(OUT, "genes_summary.parquet"), index=False)

    funnel = {
        "vcf_carried_alleles": joined,
        "vep_transcript_rows": rows_seen,
        "vep_variants_collapsed": len(variant),
        "join_row_rate": round((rows_seen - n_unres) / max(rows_seen, 1), 5),
        "join_styles": join_stats,
        "rs_multikey_records": n_rs_multi,
        "unresolved_rows_written": n_unres,
        "unresolved_with_coding_consequence": int(n_unres_coding),
        "variants_excluded_nonpass": int((~df["pass_vcf"]).sum()),
        "variants_excluded_noncanonical": int((~df["canonical"]).sum()),
        "rare_variants_primary_tier": int((df["af_tier"] == "primary").sum()),
        "rare_variants_loose_tier": int((df["af_tier"] == "loose").sum()),
        "no_population_data": int((df["af_label"] == "no population data").sum()),
        "rare_ptv": int(((df["af_tier"] == "primary") & df["is_ptv"]).sum()),
        "rare_missense_consensus": int(((df["af_tier"] == "primary") & df["missense_damaging"]).sum()),
        "rare_missense_any": int(((df["af_tier"] == "primary") & df["missense_any"]).sum()),
        "rare_splice": int(((df["af_tier"] == "primary") & df["is_splice_only"]).sum()),
        "rare_coding_other": int(((df["af_tier"] == "primary") & (df["dmg_class"] == "coding_other")).sum()),
        "genes_with_ge2_qual_het": int(gsf["compound_het"].sum()) if len(gsf) else 0,
        "genes_with_qual_hom": int((gsf["n_qual_hom"] > 0).sum()) if len(gsf) else 0,
        "panel_gene_rows": {g: int((df["gene"] == g).sum()) for g in MVA_PANEL},
        "panel_gene_rare_primary": {g: int(((df["gene"] == g) & (df["af_tier"] == "primary")).sum()) for g in MVA_PANEL},
        "panel_gene_damaging_in_board": {
            g: int(((df["gene"] == g) & in_board & df["dmg_class"].isin(DAMAGING)).sum()) for g in MVA_PANEL
        },
        "consequences_top": dict(sorted(consequences_hist.items(), key=lambda kv: -kv[1])[:30]),
    }
    with open(os.path.join(OUT, "funnel.json"), "w") as f:
        json.dump(funnel, f, indent=1)

    # model-lane inputs: SNV-only (true SNVs, no '-'), canonical contigs, PASS,
    # qualifying damage class, primary OR loose tier
    is_snv = (df["ref"].isin(list("ACGT")) & df["alt"].isin(list("ACGT")))
    lane_pool = df[in_board & is_snv & df["dmg_class"].isin(QUALIFYING)
                   & df["af_tier"].isin(["primary", "loose"])].copy()
    inert = df[in_board & is_snv & (df["consequence"] == "synonymous_variant")
               & (df["af_tier"] == "primary")]
    inert_sel = []
    for chrom, gdf in lane_pool.groupby("chrom"):
        cand = inert[inert["chrom"] == chrom]
        n = min(len(gdf) * 2, len(cand), 500)
        if n > 0:
            inert_sel.append(cand.sample(n, random_state=42))
    inert_df = pd.concat(inert_sel) if inert_sel else inert.iloc[:0]
    lane_pool["role"] = "candidate"
    inert_df["role"] = "inert_control"
    cols = ["chrom", "pos", "ref", "alt", "gene", "consequence", "dmg_class",
            "af_tier", "af_label", "vcf_filter", "dp", "ad_ref", "ad_alt", "role"]
    out_df = pd.concat([lane_pool[cols], inert_df[cols]], ignore_index=True)
    out_df.to_parquet(os.path.join(OUT, "score_inputs.parquet"), index=False)
    print(f"[03] lane inputs: {len(out_df)} ({len(lane_pool)} candidates, {len(inert_df)} inert)")

    # candidates_ranked: genes of interest with member evidence
    goi = gsf[(gsf["panel"]) | (gsf["compound_het"]) | (gsf["n_hpo_terms"] > 0) | (gsf["n_qual_hom"] > 0)] if len(gsf) else gsf
    card_rows = []
    for _, grow in goi.iterrows():
        for mem in grow["members"].split("|"):
            rec = variant.get(mem)
            if rec:
                r = {"ukey": mem, "gene": grow["gene"], "compound_het": bool(grow["compound_het"]),
                     "panel": bool(grow["panel"]), "hpo_terms": grow["hpo_terms"]}
                r.update({k: rec[k] for k in ["consequence", "impact", "hgvsc", "hgvsp", "max_af",
                                              "af_tier", "af_label", "dmg_class", "consensus_votes",
                                              "lof_hc", "lof_any", "lof_filter", "lof_flags",
                                              "am_class", "am_score", "sift", "polyphen", "clin_sig",
                                              "existing", "zygosity", "gt", "ad_ref", "ad_alt", "dp",
                                              "vcf_filter"]})
                card_rows.append(r)
    cd = pd.DataFrame(card_rows)
    cd.to_parquet(os.path.join(OUT, "candidates_ranked.parquet"), index=False)
    print(f"[03] candidate cards: {len(cd)} across {len(goi)} genes")
    print(json.dumps(funnel, indent=1))
    print("[03] MINE_OK")


if __name__ == "__main__":
    main()
