# Methods description, MVA Hackathon 2026, Track 1 (Variant Prediction)

Team name: Silico EVEE.

## Approach in detail

We built an automated rare-variant prioritization pipeline over the proband's whole-genome VCF (5,012,204 records, GRCh38) and the eight provided HPO terms.

1. **Annotation.** Every variant-allele unit was annotated with Ensembl VEP release 112: transcript consequence, SIFT, PolyPhen, LOFTEE loss-of-function confidence, ClinVar, AlphaMissense (release 151), and population allele frequencies from gnomAD exomes/genomes and 1000 Genomes. 3,237,738 variant-allele units were recovered after joining annotations to genotypes, with join failures audited separately.
2. **Rare-damaging filtration.** Primary tier requires allele frequency at or below 0.001 in every population reference; a looser 0.01 tier is reported in parallel. Damaging classes are protein-truncating variants (VEP high impact or LOFTEE high confidence), missense with 2-of-3 consensus among SIFT, PolyPhen, and AlphaMissense, splice-only, and other coding. All candidates additionally require FILTER=PASS on a canonical chromosome. 319 variants survived.
3. **Inheritance-pattern mining.** Candidates are grouped by gene into compound-heterozygous pairs (17 genes) and homozygous-alternate patterns (15 genes). The known chromosomal-instability disease panel (CEP57, BUB1B, TRIP13 and related spindle-checkpoint genes) is checked openly as one panel among the genome-wide search.
4. **Phenotype matching.** Surviving genes are matched against the proband's HPO terms using the public HPO gene-to-phenotype table.
5. **Sequence-model lanes.** Two genomic foundation models score every surviving SNV. GPN-MSA (songlab/gpn-msa-sapiens) contributes reference-minus-alternate log-likelihood over 512 bp windows. EVEE (Goodfire's public released covariance probe over Evo 2 7B block-27 embeddings, 16,384 bp windows) contributes a calibrated pathogenicity probability.
6. **Ranking.** A transparent weighted-evidence score combines HPO term matches, compound-het pairing, disease-panel membership, ClinVar precedent in the gene, homozygous pattern, per-variant damage class scaled by rarity, and averaged model-lane percentiles. Pairs outrank singletons.

## Automated output or manual review

The ranking is automated end to end, with one disclosed manual override. The purely mechanical total placed a homozygous 45 bp deletion in PEX5 (31.0) ahead of the compound-heterozygous BUB1B pair (27.4). We promoted the BUB1B pair to rank 1 because BUB1B is the known MVA1 gene and matches the phenotype cluster (rhabdomyosarcoma predisposition, growth restriction), whereas homozygous PEX5 loss causes a severe neonatal peroxisome-biogenesis spectrum without cancer predisposition. The swap is documented on the PEX5 row of the submission file and in this report. No other manual curation was applied.

## Data

Publicly available data only: gnomAD (v2.1.1 exomes, v3.1.2 genomes), 1000 Genomes, ClinVar (August 2026 weekly release), AlphaMissense hg38, LOFTEE, Ensembl gene annotation release 112, GRCh38 reference sequence, the multiz100way alignment store, the HPO gene-to-phenotype table, songlab GPN-MSA public weights, and Goodfire's publicly released EVEE probe weights. No proprietary data.

## Compound heterozygous pairs

Yes. The pipeline natively outputs compound-heterozygous candidate pairs, ranked ahead of singletons, and encodes them in one submission row using the second-variant columns. Note that with proband-only VCF data, trans-phase (which parent carried each allele) cannot be verified; the pair is inferred from two qualifying heterozygous variants in one gene.

## Secondary findings

Secondary rows in the submission carry evidence-based scores with finding_type=secondary, including a homozygous protein-truncating deletion in PEX5 that matched four HPO terms mechanically. Secondary findings are separated from the primary candidate so judges can review them independently; they do not affect automated scoring.

## Run time and cost

Total wall-clock approximately 6 hours on cloud single-node infrastructure, including dataset download, annotation (about 5 hours of CPU time), candidate mining, and both model lanes (about 3 GPU-hours on one B200-class GPU). Monetary cost is on the order of tens of dollars of cloud compute.

## Method abstract (about 300 words)

We prioritized causal-variant candidates for a pediatric MVA case from a single whole-genome VCF plus eight HPO terms. All 5 million variant records were annotated with VEP (consequence, SIFT, PolyPhen, LOFTEE, ClinVar, AlphaMissense, gnomAD and 1000 Genomes frequencies), filtered to 319 rare damaging candidates, grouped by gene into compound-heterozygous and homozygous patterns, matched to the phenotype, and scored by two genomic foundation models (GPN-MSA zero-shot sequence likelihood and Goodfire's released EVEE probe). A transparent weighted-evidence score produced the ranking, with one disclosed manual override promoting the BUB1B compound-heterozygous pair over a PEX5 homozygous deletion on syndrome specificity.

**Headline result.** Rank 1 is a compound-heterozygous pair in BUB1B (the known MVA1 gene): a ClinVar-pathogenic stop-gain (p.Leu737Ter, heterozygous allele depth 21/25) plus a missense (p.Asn1002Lys) never observed in population databases with deleterious predictions from SIFT, PolyPhen, and elevated EVEE pathogenicity. Three proband HPO terms map to BUB1B exactly.

**Strengths.** Fully automated and auditable point system; inheritance-pattern awareness; independent corroboration from two sequence models; open disease-gene panel checked inside a genome-wide search so the pipeline could, in principle, have missed BUB1B.

**Limitations.** No parental data, so compound-het phasing is inferred rather than proven. Indels carry annotation evidence only (sequence-model lanes scored SNVs). EVEE windows were reduced from 65 kb to 16 kb for cost. The ranking weights are a case-specific heuristic, not a calibrated published method, and manual review influenced the rank-1 decision as disclosed.

## Supporting materials

- Predictions file: `Silico-EVEE_bub1b-comphet.csv`.
- Code: public GitHub repository link TO FILL (pipeline source is packaged with this submission package).
