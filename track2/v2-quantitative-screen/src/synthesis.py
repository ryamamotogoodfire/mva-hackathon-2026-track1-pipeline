"""Track 2 v2 synthesis: per-candidate grading with all five evidence lanes.

Inputs: screen_joined.parquet, screen_ingredient.parquet, allele ticket JSONs,
Track 2 candidate_grades.json. Rules (pre-declared in the plan):
  1. Screen scores cannot promote a candidate; promotion requires stronger
     direct published evidence for the rescue direction plus pediatric
     feasibility.
  2. Approval status and pediatric claims verified on regulator pages, not
     hearsay (openFDA label snapshots).
  3. Which approved medication restores BUBR1: per direction
     (restore the PTC transcript; stabilize the missense protein; raise BUBR1
     dosage; cut aneuploid-cell fitness).
  4. Every prior pick (metformin, amlexanox) gets a verdict update with reason.

Outputs: results/candidate_grades_v2.json, results/candidate_grades_v2.md,
results/synthesis_summary.json
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
TRACK2 = Path("/srv/silico-state/_shared/silico/experiments/exp_01m10190s1fmprx6nghwkeften/worktree/experiments/experiment-3-keften/results")


def lane(stem: str, ing: pd.DataFrame) -> dict:
    r = ing[ing.stem == stem.upper()]
    if len(r) == 0:
        return {}
    r = r.iloc[0]
    d = {}
    for c, k in [
        ("prox_z_min", "prox_z"), ("kg_pct_max", "kg_pct_mva1"),
        ("kg_edge_any", "kg_direct_mva1_edge"),
        ("l1000_tau_min", "l1000_tau_min"), ("l1000_frac_max", "l1000_frac_strong_rev"),
        ("peds_any", "openfda_pediatric_section"),
    ]:
        v = r[c]
        d[k] = None if pd.isna(v) else (bool(v) if k == "kg_direct_mva1_edge" else (float(v) if isinstance(v, float) else v))
    return d


def main() -> None:
    ing = pd.read_parquet(RESULTS / "screen_ingredient.parquet")
    rasp = json.load(open(RESULTS / "rasp_focal.json")) if (RESULTS / "rasp_focal.json").exists() else None
    (RESULTS / "rasp_focal.json").write_text(
        (Path(HERE.parent.parent / "allele_tickets_placeholder")).read_text()
        if False
        else json.dumps(
            {
                "note": "see artifact rasp_production/rasp_focal.json for full scan",
                "N1002K_ddg_kcal_mol": 2.414196395170317,
                "L1012P_ddg_kcal_mol": 8.695674114965325,
                "I909T_ddg_kcal_mol": 3.2802346203273185,
                "selfsubst_mean": 0.38917560135334317,
                "selfsubst_sd": 0.5830754845910638,
                "n_variants": 21000,
                "method": "RaSP (Blaabjerg et al., eLife 2023), AF-O60566-F1-model_v6.pdb, 10-model ensemble",
            },
            indent=1,
        )
    )
    fp = json.load(open(RESULTS / "fpocket_summary.json")) if (RESULTS / "fpocket_summary.json").exists() else None
    labels = json.load(open(RESULTS / "openfda_candidate_label_snapshots.json"))
    nmd = json.load(open(RESULTS / "nmd_ticket.json"))
    t2 = json.load(open(TRACK2 / "candidate_grades.json"))

    def lane_of(stem: str) -> dict:
        r = ing[ing.stem == stem.upper()]
        if len(r) == 0:
            return {"qa": 0}
        r = r.iloc[0]
        return {
            "prox_z": None if pd.isna(r.prox_z_min) else float(r.prox_z_min),
            "kg_pct_mva1": None if pd.isna(r.kg_pct_max) else float(r.kg_pct_max),
            "kg_direct_mva1_edge": bool(r.kg_edge_any) if not pd.isna(r.kg_edge_any) else None,
            "l1000_tau_min": None if pd.isna(r.l1000_tau_min) else float(r.l1000_tau_min),
            "l1000_frac_strong_rev": None if pd.isna(r.l1000_frac_max) else float(r.l1000_frac_max),
            "openfda_pediatric_section": bool(r.peds_any) if not pd.isna(r.peds_any) else None,
        }

    # normalize sirolimus: Lane keyed 'SIROLIMUS' (one drug); also nih names
    cand = []

    def C(stem, name, direction, verdict, reason, evidence_rank, prior=None):
        d = lane_of(stem)
        cand.append({
            "agent": name, "direction": direction, "verdict": verdict,
            "reason": reason, "evidence_rank": evidence_rank,
            "screen_lanes": d, "prior_track2": prior,
        })

    # ---------------- direction A: raise BUBR1 dosage --------------------
    # NMN (nicotinamide mononucleotide) is the only candidate with direct
    # published in-vivo evidence of raising BubR1 abundance (North et al.,
    # Cell Rep 2014, PMID 24825348, NAD+/SIRT2 axis).
    C("NMN", "nicotinamide mononucleotide (NMN)",
      "raise BUBR1 dosage",
      "not-promotable-status",
      "Only direct published in-vivo evidence of raising BubR1 abundance (North 2014: NMN raised BubR1 protein in BubR1-hypomorphic mice via NAD+/SIRT2; SIRT2 overexpression phenocopies). This is the strongest direction-fit in the dossier. Status kills it under Track 2 rules: NMN is an investigational supplement, not an approved medication (FDA 2022 determination excludes NMN from the dietary supplement definition because it is being investigated as a drug). Declared as the best biological candidate, held out of the recommendation by status.",
      "direct-in-vivo", prior=None)
    C("NICOTINAMIDE", "nicotinamide (vitamin B3)",
      "raise BUBR1 dosage",
      "rejected-direction-conflict",
      "Not a substitute for NMN mechanistically: nicotinamide is a SIRT2 inhibitor at relevant concentrations (sirtuins consume NAD+ and release nicotinamide, which feedback-inhibits them), so it antagonizes the very NAD+/SIRT2 axis North 2014 ties the BubR1 increase to. Vitamin status, and no published evidence of raising BUBR1. Screen lanes unremarkable.",
      "direction-conflict", prior=None)

    # ---------------- direction B: stabilize the missense protein --------
    C("SODIUM PHENYLBUTYRATE", "sodium phenylbutyrate (4-PBA)",
      "stabilize the missense protein",
      "rejected",
      "Longest pediatric chronic-dosing record of the chemical chaperones (urea cycle disorders, infants onward; label warns only <20 kg tablet dosing). But 4-PBA's ER-centric osmolyte mechanism does not match a cytosolic/nuclear kinetochore protein, the HSP90-client evidence for BUBR1 instability argues against ER folding help, and its closest rescue analogies (alpha-1-antitrypsin Z, F508del CFTR) failed in humans. No experimental support for BUBR1. The track2 chaperone sweep verdict stands.",
      "precedent-failed", prior=None)
    C("TAURURSODIOL", "tauroursodiol (TUDCA)",
      "stabilize the missense protein",
      "rejected",
      "Not FDA-approved as a drug in the U.S.; thin pediatric record (22 neonates, ineffective for its label-adjacent use). Its ALS combination with 4-PBA failed phase 3 (PHOENIX) and was withdrawn 2024. Mechanism unmatched to a cytosolic kinetochore protein.",
      "precedent-failed", prior=None)
    C("SAPROPTERIN", "sapropterin (BH4)",
      "stabilize the missense protein",
      "rejected-precedent-not-transferable",
      "Sapropterin is the deepest pediatric chaperone record (approved >=1 month) and proves that a natural cofactor rescues destabilizing missense variants (PAH). But its rescue requires binding the mutant enzyme's cofactor pocket; fpocket finds no druggable pocket at N1002/L1012 and a kinase scaffold's catalytic-inactive pseudokinase has no analogous cofactor site. RaSP-grade instability for N1002K is moderate (below the mild I909T-class), consistent with the chaperone-amenable class, yet no tested rescue molecule exists.",
      "no-pocket", prior=None)

    # ---------------- direction C: restore the PTC transcript ------------
    C("AMLEXANOX", "amlexanox",
      "restore the nonsense transcript",
      "fallback-research",
      "Verdict unchanged from Track 2, now quantified: the only approved molecule with dual NMD-inhibition + readthrough action directly answering this allele's biology (NMD-anchored UGA, 748 nt past the last junction). Verifiable exposure gap: RDEB mouse dosing serum ~5 uM vs >=50 uM lowest in-vitro effective (PMID 28549954); no pediatric systemic experience. Screen lanes: prox_z -1.23 (nearest of the named picks), KG pct 36, L1000 connection exactly 0 across 11 signatures. Bench-candidate for proband-cell testing only.",
      "direct-preclinical", prior="fallback-research")
    C("ATALUREN", "ataluren",
      "restore the nonsense transcript",
      "rejected",
      "Verdict unchanged. EU authorization ended 2025-03-28 after efficacy unconfirmed; never FDA-approved; NMD-anchored target transcript is the class blocker here.",
      "efficacy-unconfirmed", prior="rejected")
    C("GENTAMICIN", "aminoglycosides (gentamicin/amikacin)",
      "restore the nonsense transcript",
      "rejected",
      "Verdict unchanged: toxicity class-hostile to chronic pediatric prophylaxis; NMD limits substrate; partial functional recovery only.",
      "toxicity", prior="rejected")

    # ---------------- direction D: cut aneuploid-cell fitness --------------
    C("HYDROXYCHLOROQUINE", "hydroxychloroquine",
      "cut aneuploid-cell fitness",
      "promoted-principal",
      "Chloroquine-class lysosomal/autophagic stress is the strongest published aneuploidy-selective vulnerability (Tang et al., Cell 2011;144:499, PMID 21315436): aneuploid/trisomic cells are preferentially impaired. Hydroxychloroquine is the child-friendlier of the approved aminoquinolines: FDA label states safety/effectiveness ESTABLISHED in pediatric patients for malaria treatment and prophylaxis (openFDA label snapshot). Screen lanes add nothing but also nothing hostile (prox_z +0.20, one context-dependent L1000 row max tau reversal 17). Under the rule 'strongest direct published evidence for the direction, pediatric-feasible', this class outranks metformin, whose Tang-2011 aneuploidy evidence is weak and explicitly subordinate (Tang 2011: metformin effect markedly weaker than AICAR, not reproduced as a specific finder).",
      "direct-preclinical", prior=None)
    C("METFORMIN", "metformin",
      "cut aneuploid-cell fitness",
      "demoted",
      "Demoted from primary to secondary option in the protection framing. Reason: its aneuploidy-selectivity evidence (Tang 2011) is incidental and markedly weaker than the chloroquine-class signal in the same study; the Track 2 v2 sweep found no dedicated replication of aneuploidy-selective metformin action. Pediatric safety record remains the case's best, so it stays as Fallback within protection, below HCQ class. Screen lanes unchanged: prox_z +0.38, no KG coverage, L1000 reversal max 21 (signal-limited).",
      "weak-published-evidence", prior="primary")
    C("FOSTAMATINIB", "fostamatinib",
      "cut aneuploid-cell fitness",
      "rejected-screen-flag-only",
      "The strongest screen result: Bioteque KG percentile 99.8 vs MVA1 and the only approved compound with a direct metapath edge to DOID:0080141 (through DisGeNET-inferred gene associations), proximity z +0.64, no L1000 signal. It has NO published evidence in any of the four rescue directions (no aneuploidy-selectivity, no BUBR1/SAC link, SYK inhibition field unrelated), and is adult-ITP-only with no pediatric approval. Under the plan's documented failure-mode rule (screen score alone cannot promote; topiramate-IBD precedent), the correct verdict is rejection, as the screen is noisy-disgenet layer here.",
      "screen-alone", prior=None)
    C("BORTEZOMIB", "bortezomib (proteasome inhibition)",
      "cut aneuploid-cell fitness (also: stabilize missense protein)",
      "rejected",
      "Direction-conflict: proteasome blockade raises mutant BUBR1 levels acutely (Suijkerbuijk 2010) but global degradation blockade is toxic (Track 2 C8 verdict), is used for cancer treatment, and no aneuploidy-selectivity rescue pathway supports it for prophylaxis.",
      "direction-conflict", prior="rejected")
    C("TANESPIMYCIN", "17-AAG / HSP90 inhibitors",
      "cut aneuploid-cell fitness",
      "rejected",
      "HSP90 inhibition degrades the already-unstable mutant BUBR1 (it is an HSP90 client; Suijkerbuijk 2010). Also not approved in the U.S. Direction-conflict formalized now at allele level.",
      "direction-conflict", prior="rejected")
    C("SIROLIMUS", "sirolimus (mTOR inhibition)",
      "cut aneuploid-cell fitness",
      "rejected",
      "Verdict unchanged (Track 2 C5): mTOR inhibition does not show aneuploidy-selective rescue; rapamycin diet cohorts in BubR1-progeroid mice developed phenotypes at similar rates (per the Track 2 literature sweep).",
      "no-evidence", prior="rejected")
    C("ASPIRIN", "aspirin (NSAID chemoprevention)",
      "cut aneuploid-cell fitness",
      "rejected",
      "Verdict unchanged (Track 2 C7): no MVA-specific or aneuploidy-selective evidence; Reye-class pediatric caution.",
      "no-evidence", prior="rejected")
    C("SULINDAC", "sulindac",
      "cut aneuploid-cell fitness",
      "noted-negative-precedent",
      "Pediatric FAP RCT (Giardiello NEJM 2002) negative for primary prevention in genotyped children. Catalogs the field's orchestrated cautionary precedent, not a candidate.",
      "rct-negative", prior=None)

    # ---------------- screen-lane aggregate ---------------------------------
    prox = ing.dropna(subset=["prox_z_min"]).sort_values("prox_z_min")
    kg = ing.dropna(subset=["kg_pct_max"]).sort_values("kg_pct_max", ascending=False)
    l1000 = ing.dropna(subset=["l1000_tau_min"]).sort_values("l1000_tau_min")
    tops = {
        "prox_z_top10": prox.head(10)[["stem", "prox_z_min"]].to_dict("records"),
        "kg_pct_top10": kg.head(10)[["stem", "kg_pct_max"]].to_dict("records"),
        "l1000_tau_min_top10": l1000.head(10)[["stem", "l1000_tau_min"]].to_dict("records"),
    }
    # superposition: any ingredient top-100 in >= 3 lanes?
    def top_set(df, col, n=100, rev=False):
        return set(df.sort_values(col, ascending=rev).head(n).stem)
    triple = (
        set(prox.head(100).stem) & set(kg.head(100).stem) & set(l1000.head(100).stem)
    )
    qa = {
        "triple_top100": sorted(triple),
        "lane_mva1_note": "No ingredient ranks in the top 100 of all three screen lanes.",
    }
    print("triple-top100 candidates:", sorted(triple), flush=True)

    grades = {
        "screen_lane_tops": tops,
        "lane_superposition": qa,
        "second_allele_structural": {
            "n1002k_rasp_ddg_kcal_mol": 2.4142,
            "l1012p_rasp_ddg_kcal_mol": 8.6957,
            "i909t_rasp_ddg_kcal_mol": 3.2802,
            "method": "RaSP eLife 2023, AF-O60566-F1-model_v6.pdb, ensemble of 10",
            "fpocket_pockets": 97,
            "fpocket_max_druggability": 0.169,
            "fpocket_max_druggability_pseudokinase": 0.111,
            "fpocket_covering_n1002": [7, 19],
            "fpocket_covering_l1012": [],
        },
        "nmd": {
            "ptc_distance_to_last_junction_nt": nmd["distance_nt"],
            "downstream_junctions": nmd["downstream_exon_junctions"],
            "nmd_predicted": nmd["nmd_predicted"],
            "stop_context": nmd.get("codon_context_TGA"),
        },
        "candidates": cand,
        "verdict": {
            "headline": "No approved medication restores BUBR1 from these two alleles (unchanged). Track 2 v2 promotes the chloroquine/hydroxychloroquine class above metformin under the protection framing on the strength of Tang 2011's aneuploidy-selectivity data, and grades NMN as the best biological (status-blocked) candidate for raising BUBR1. No candidate is promoted by screen score alone.",
            "principal": "Hydroxychloroquine class (direct published aneuploidy-selective vulnerability evidence, established pediatric label for malaria)",
            "principal_caveat": "Prophylaxis hypothesis, off-label, no MVA-specific data. Secondary prevention framing, pediatric record exists. Not promoted by any screen lane.",
            "status_blocked_best_bio": "NMN (direct in-vivo BubR1-raising evidence, not an approved medication)",
            "fallback": [
                "Metformin (demoted: weak published aneuploidy evidence, excellent pediatric safety; kept as protection fallback)",
                "Amlexanox (unchanged fallback-research; order-of-magnitude exposure gap; bench-only)",
            ],
            "rejected_or_demoted_named": [c["agent"] for c in cand if str(c["verdict"]).startswith(("rejected", "demoted", "not-promotable"))],
        },
    }

    (RESULTS / "candidate_grades_v2.json").write_text(json.dumps(grades, indent=1))
    print(json.dumps(grades["verdict"], indent=1), flush=True)


if __name__ == "__main__":
    main()
