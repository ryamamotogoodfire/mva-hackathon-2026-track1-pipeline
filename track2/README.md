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


## Track 2 v2 (2026-08 follow-on)

The quantitative re-grade of this candidate list lives in [v2-quantitative-screen/](v2-quantitative-screen/) and in the v2 deliverables below. Headlines: no approved medication restores BUBR1 (unchanged); the hydroxychloroquine/chloroquine class is promoted above metformin as the secondary-prevention hypothesis (Tang 2011 aneuploidy-selective evidence + established pediatric malaria label); NMN is the best biological but status-blocked route to raise BUBR1 (North 2014, EMBO J, PMID 24825348); every verdict update is recorded with its reason in [v2-quantitative-screen/results/candidate_grades_v2.md](v2-quantitative-screen/results/candidate_grades_v2.md).


## Track 2 v2 deliverables (2026-08-27)

- [Silico-EVEE_track2_report_v2.md](Silico-EVEE_track2_report_v2.md) supersedes the v1 report: same verified mechanism chain, candidate grading rebuilt on the quantitative full-FDA screen, all verdict updates with reasons, revised falsifiable statements, and the v2 answers to the official methods form.
- [Silico-EVEE_track2_methods_description_form_v2.xlsx](Silico-EVEE_track2_methods_description_form_v2.xlsx) is the official methods form, Track 2 sheet, with every answer revised to v2 (sha256 f0507f523860e5b58553efe8432b8ee17d1c60865030ffe476fef5e6bbb65b3e, recorded in the experiment's verification note).
- [Silico-EVEE_track2_pitch_video_outline_v2.md](Silico-EVEE_track2_pitch_video_outline_v2.md) is the refreshed three-minute pitch outline.
- [v2-quantitative-screen/](v2-quantitative-screen/) holds the screen README plus every graded writeup and table the report cites.
