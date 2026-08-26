# Track 2. Drug repurposing for the BUB1B compound-heterozygous case

This chapter adds the Track 2 analysis to the repository. It builds on the Track 1 top-ranked pair without changing any Track 1 coordinate or result.

## What this is

- Silico-EVEE_track2_report.md, the full Track 2 submission report with complete citations and the methods-description-form answers. Named per the Track 2 file convention (team-name prefix).
- results/allele_map.json, the verified mechanism map of the two alleles.
- results/candidate_grades.json and results/candidate_grades.md, the graded candidate screen.
- results/recommendation_synthesis.json, the falsifiable recommendation set.
- src/, the scripts that produced every file above, CPU-only, rerunnable in under an hour.

## How to reproduce

From this directory, with a Python 3.11+ environment that has only the standard library for steps 1 to 4, and plotly plus the silico-figures package for step 5:

1. python src/allele_map.py, writes results/allele_map.json after querying Ensembl REST, UniProt, and AlphaFold DB.
2. python src/evidence_screen.py, writes results/evidence_raw.json and results/literature_hits.json from Europe PMC, ChEMBL, and DailyMed.
3. python src/regulatory_status.py, writes results/regulatory_status.json from ChEMBL.
4. python src/build_candidates.py and python src/build_synthesis.py, write the curated grading and synthesis tables.
5. python src/build_figures.py, writes figure bundles used by the results page.

## Headline findings

- The stop-gain transcript p.Leu737Ter is destroyed by nonsense-mediated decay 748 nucleotides before the last exon junction, so no approved readthrough medicine can act on it. UGA is the readthrough-friendliest stop class and the TGAA context is permissive, so NMD, not the stop context, is what disables the direct-readthrough class with drugs alone. That is exactly why amlexanox, an approved NMD-inhibiting readthrough promoter, is the bench-grade fallback.
- The missense allele p.Asn1002Lys is a last-exon change inside the pseudokinase domain and escapes decay, class-consistent with instability-driven abundance loss.
- No approved medication restores BUBR1 from these alleles. Metformin is proposed as a secondary-prevention hypothesis. Amlexanox is named at bench grade only.
- Recorded rejections with reasons, ataluren, aminoglycosides, sirolimus and everolimus, checkpoint kinase inhibitors, aspirin, bortezomib and Hsp90 inhibitors.
