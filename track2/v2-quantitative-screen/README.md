# Track 2 v2: quantitative full-FDA screen follow-on (2026-08)

This subdirectory holds the quantitative re-grade of the Track 2 drug-repurposing list, produced by the follow-on Silico experiment (the "Quantitative full-FDA screen for BUBR1 rescue" run). It supersedes no Track 2 coordinate, allele fact, or class rejection; when a prior pick changes rank, the change is stated with its reason.

## What the screen did

Three orthogonal quantitative lanes over all FDA-approved drugs with mechanism targets in ChEMBL 37 (1,811 drugs, 650 human targets, 3,927 drug-gene rows), plus per-allele structural tickets:

1. **Network proximity**: Guney/Barabasi-degree-matched z-scores over BioGRID 4.4.236 human PPI (19,908 proteins, 849,005 edges) against a declared MVA module (BUB1B, BUB1, CEP57, TRIP13, MAD1L1, MAD2L1, BUB3, TTK, CDC20, CENPE, AURKB, PLK1). Calibration: ivacaftor-CFTR z=-3.41, imatinib-BCR-ABL1 z=-10.66, simvastatin-MVA z=+0.37.
2. **Knowledge graph**: Bioteque embeddings (CPD-int-GEN-ass-DIS, disgenet bundle). MVA1 = DOID:0080141 scored by cosine; the bundle recovers its own edges at AUC 0.917.
3. **LINCS reversal**: BUB1B-knockdown proxy signature (58 quality-gated trt_sh signatures, GEO GSE92742) scored against 205,034 compound signatures, tau-style. Instrument checks: own knockdown readout passes (z=-2.95), the published MTOR-sirolimus pair reconnects (+19.1 mean cs), but held-out hairpins do not (0.101 vs null p95 0.184), so this lane is instrument-limited and no verdict rests on it.
4. **Allele tickets**: RaSP stability ddG for p.Asn1002Lys on AF-O60566-F1 (+2.41 kcal/mol; validated anchors L1012P +8.70, I909T +3.28), fpocket pocket map (no druggable pocket at the locus), and the independently recomputed NMD margin (748 nt, six downstream junctions, exon 17/23).

The pre-declared rule decides: a candidate promotes only with the strongest direct published evidence for its rescue direction plus pediatric feasibility, never via screen score alone.

## What the screen says (Q1)

No approved medication is promoted by the quantitative screen itself. **No ingredient ranks in the top 100 of all three lanes.** The screen's strongest beat, fostamatinib (knowledge-graph percentile 99.8, the only approved compound with a direct metapath edge to MVA1), is rejected: no published evidence in any rescue direction and no pediatric approval. This is the field's documented screen-only failure pattern and the plan's stated refutation watch.

## Verdict updates (Q2)

- **Hydroxychloroquine/chloroquine class promoted above metformin** as the secondary-prevention hypothesis: tang-2011's aneuploidy-selective vulnerability data is the strongest published route for that direction (chloroquine stronger than metformin in the same study); hydroxychloroquine's FDA label states safety and effectiveness *established* in pediatric patients for malaria treatment and prophylaxis. The same label limits the film-coated tablets to children of at least 31 kg (tablets cannot be crushed or divided), so the pediatric leg narrows for small children; the class route under 31 kg is chloroquine liquid 16.67 mg/ml, the compound that carried the Tang 2011 aneuploidy-selectivity evidence. This is a hypothesis for secondary prevention in this MVA case, off-label for MVA.
- **Metformin demoted** from primary to protection fallback: the aneuploidy-selectivity evidence is incidental and weak (subordinate to AICAR in Tang 2011, no dedicated replication). Pediatric safety record remains the case's best.
- **NMN** graded the best biological candidate for raising BUBR1: North 2014 raised BubR1 protein in young and aged wild-type mice via NMN, while the paper's BubR1-hypomorphic cohorts tested SIRT2 overexpression. NMN is **status-blocked**: FDA reversed its 2022 drug-preclusion exclusion on 2025-09-29 (response to citizen petition FDA-2023-P-0872), so NMN is a lawful dietary supplement, and a lawful supplement is still not an FDA-approved medication. Nicotinamide is a direction-conflict (SIRT2 inhibitor) and rejected.
- **Amlexanox unchanged** bench-only (dual NMD+readthrough mechanism answers the allele's biology; exposure-vs-potency gap ~10x; no L1000 reversal signal for it).
- **Other Track 2 verdicts hold** (ataluren, aminoglycosides, rapalogs, aspirin, proteostasis inhibitors, checkpoint-kinase inhibitors rejected; HSP90-inhibition direction-conflict for an unstable BUBR1 client now explicit).

## Two additions that close gaps the lane screen left (2026-08-27)

**Approved NAD+ precursors.** The screen's drug universe is target-annotated pharmacology: a drug enters only through a ChEMBL 37 mechanism record with a resolvable human protein target. Niacin (CHEMBL573) and nicotinamide (CHEMBL1140) carry zero mechanism records, and acipimox (CHEMBL345714) carries one record with a null target, so the whole NAD+ precursor class was invisible to all three lanes — a coverage gap, not a negative result, and it falls on exactly the class the North 2014 NAD+/SIRT2-to-BubR1 mechanism points to. The class was therefore graded separately on regulator facts and mechanism tie (`results/nad_precursor_grades.{json,md}`). Verdicts: niacin is a genuine FDA prescription NAD+ precursor (NIACOR ANDA040378 plus extended-release ANDAs) and is still rejected for this case, because its labels state pediatric safety and effectiveness are not established and it has no BubR1 or aneuploidy data; nicotinamide is a direction-conflict through SIRT2 inhibition; acipimox has no FDA approval at all; nicotinamide riboside and NMN stay status-blocked as supplements. The dosage-raise axis is therefore not empty of approved medicines, and nothing on it is promotable.

**Pseudokinase stabilizer screen.** A bounded computation asked whether any approved small molecule could act as a tafamidis-style chaperone on the BUBR1 pseudokinase domain, binding away from the N1002K mutation. The published structures 6JKK and 6JKM are *Drosophila* BubR1, so the ADP-contact site of the 6JKM ADP+Mg complex was mapped onto human numbering (validated at K795, the kinase-dead anchor). p2rank 2.4.2 ranks that site first on the human AlphaFold model and puts its nearest pocket 25.05 angstrom from N1002, agreeing with fpocket that no pocket sits at the mutation. 2,240 ligands (approved small molecules plus ADP/ATP references) were docked with smina against three boxes (nucleotide site, the two fpocket pockets touching N1002, the most druggable domain pocket), then the top 40 per box were re-docked at exhaustiveness 16 with Vinardo rescoring (35.5 CPU-hours total). The result is an honest negative bounded by a declared instrument limit: redocking ADP into its own crystal misses the experimental pose by 6.03 angstrom, ADP reaches only the 83rd percentile in its own pocket with 372 ordinary approved drugs outscoring it, known chaperones land at both extremes of the ranking (migalastat bottom 6 percent, lumacaftor above the 99th percentile), and no molecule separates from background. A stabilizer route for this allele would need an experimentally discovered site. Docking inputs, scores, and code live under `results/dock/` and `src/dock_*.py`.

## Files

- `results/` — candidate_grades_v2.{json,md}, screen tables (drug_target_genes, proximity_scores, bioteque_scores, screen_joined/ingredient), nmd_ticket, rasp_focal, fpocket_summary, openfda snapshots, nad_precursor_grades.{json,md}, atp_site_annotation.json.
- `results/dock/` — dock_screen_summary.json (stage-1/2 scores for all boxes), dock_validation.json (p2rank + redock controls), dock_named_candidates.json (percentiles for every graded candidate), stage-1 score parquets per box.
- `src/` — the full v2 pipeline: lane code, allele tickets, the docking screen (dock_common/validate/screen), the NAD+ grading, and the report fold script.
- `report_summary.json`, `claims_manifest.json` — machine-readable verdict and claim-to-source map used for the results page.
- `report.md` — plain reading of the v2 grading.

Full pipeline code: the Silico experiment worktree (experiment-4-m81y59) carries src/ and vendor/ (RaSP replay) with provenance notes; the external inputs (ChEMBL, openFDA, BioGRID, Bioteque, GEO, AlphaFold) are fetched reproducibly at run time.
