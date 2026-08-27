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

- **Hydroxychloroquine/chloroquine class promoted above metformin** as the secondary-prevention hypothesis: tang-2011's aneuploidy-selective vulnerability data is the strongest published route for that direction (chloroquine stronger than metformin in the same study); hydroxychloroquine's FDA label states safety and effectiveness *established* in pediatric patients for malaria treatment and prophylaxis. This is a hypothesis for secondary prevention in this MVA case, off-label for MVA.
- **Metformin demoted** from primary to protection fallback: the aneuploidy-selectivity evidence is incidental and weak (subordinate to AICAR in Tang 2011, no dedicated replication). Pediatric safety record remains the case's best.
- **NMN** graded the best biological candidate for raising BUBR1 (North 2014, in vivo BubR1 raise via NAD+/SIRT2) and **status-blocked**: NMN is not an approved medication. Nicotinamide is a direction-conflict (SIRT2 inhibitor) and rejected.
- **Amlexanox unchanged** bench-only (dual NMD+readthrough mechanism answers the allele's biology; exposure-vs-potency gap ~10x; no L1000 reversal signal for it).
- **Other Track 2 verdicts hold** (ataluren, aminoglycosides, rapalogs, aspirin, proteostasis inhibitors, checkpoint-kinase inhibitors rejected; HSP90-inhibition direction-conflict for an unstable BUBR1 client now explicit).

## Files

- `results/` — candidate_grades_v2.{json,md}, screen tables (drug_target_genes, proximity_scores, bioteque_scores, screen_joined/ingredient), nmd_ticket, rasp_focal, fpocket_summary, openfda snapshots.
- `report_summary.json`, `claims_manifest.json` — machine-readable verdict and claim-to-source map used for the results page.
- `report.md` — plain reading of the v2 grading.

Full pipeline code: the Silico experiment worktree (experiment-4-m81y59) carries src/ and vendor/ (RaSP replay) with provenance notes; the external inputs (ChEMBL, openFDA, BioGRID, Bioteque, GEO, AlphaFold) are fetched reproducibly at run time.
