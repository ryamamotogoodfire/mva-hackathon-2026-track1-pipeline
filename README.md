# MVA Hackathon 2026, Track 1. Rare-variant prioritization on the proband whole genome

This repository holds the pipeline behind our Track 1 submission for the Rare Disease, Real Kid MVA Hackathon 2026. The pipeline searches a pediatric whole-genome VCF for candidate causal variants and ranks them. Our submitted rank 1 is a compound-heterozygous pair in BUB1B, the gene behind mosaic variegated aneuploidy type 1. The pair combines a ClinVar-pathogenic stop gain with a previously unobserved damaging missense.

The repository reproduces the full analysis from the gated challenge dataset to the ranked submission CSV. It contains no case data. The proband VCF, phenotype document, and every derived variant table stay outside this repository.

## Glossary

- HPO, Human Phenotype Ontology. The vocabulary of clinical signs used to describe the proband.
- VEP, Ensembl Variant Effect Predictor. The annotation tool that labels each variant-allele unit with consequence, frequencies, and deleteriousness predictors.
- PTV, protein-truncating variant. A stop gained, frameshift, or similar allele that shortens the protein.
- Compound heterozygous. Two different damaging alleles in the two copies of one gene. With proband-only data the trans phase is inferred from co-occurrence, since parental reads are unavailable.
- AF, allele frequency. The highest frequency observed across population references decides the rarity tier.
- GPN-MSA. A genomic foundation model that scores a variant from multiple-sequence alignment context. The score is the alternate-minus-reference log-likelihood ratio, and more negative means more damaging under the model.
- EVEE. A released supervised probe over Evo 2 foundation-model embeddings that returns a pathogenicity probability per variant.
- LOFTEE. A VEP plugin that grades loss-of-function confidence.
- EPCR, estimated probability of causal relationship. The confidence field on (0, 1] required by the challenge submission schema, where 1 is most likely causal.

## Repository contents

- `pipeline/` holds the runnable source. It covers download and QC, VEP annotation, candidate mining, the GPN-MSA scoring lane, the combine-and-rank step, the submission writer, and the report figures.
- `submission/TEAMNAME_track1_report.md` is the methods report for the challenge upload. Rename TEAMNAME with your team name before uploading.
- `submission/TEAMNAME_bub1b-comphet.csv` is the validated submission file in challenge schema. Rename TEAMNAME with your team name. It matches `pipeline/track1_submission_template.csv` in format.
- `LICENSE` is the CC BY 4.0 text required by the hackathon for released submissions.

The EVEE lane is the one stage not runnable from this repository alone. It ran with Goodfire's publicly released code at github.com/goodfire-ai/evee-manuscript, pinned to commit e5e4a43d60bb217819144c09ab509ed44f49cd5f. The exact reproduction parameters are in the EVEE lane section below. Everything else runs from `pipeline/` with public dependencies only.

## External inputs the reader must obtain

1. The gated challenge dataset SageBio/mva-hackathon-2026-data on Hugging Face. Access requires accepting the dataset terms. Set `HF_TOKEN` to a token for an accepted account. The pipeline uses the proband VCF, its index, and the phenotype document. The phenotype document contributes only the list of HPO term identifiers.
2. GRCh38 reference sequence. Step 2 downloads Homo_sapiens.GRCh38.dna.primary_assembly from Ensembl release 112 automatically.
3. Ensembl VEP release 112 cache for GRCh38, downloaded by step 2.
4. The public HPO gene-to-phenotype table genes_to_phenotype.txt, downloaded automatically by `pipeline/make_hpo_gene_sets.py` from purl.obolibrary.org.
5. GPN-MSA weights from the public Hugging Face model songlab/gpn-msa-sapiens, plus the streamed multiz100way 89-species alignment store songlab/multiz100way.
6. The public EVEE probe weights and demo code from github.com/goodfire-ai/evee-manuscript at the pinned commit above.
7. Public annotation databases pulled by step 2. ClinVar weekly VCF for GRCh38, AlphaMissense hg38, LOFTEE at commit a46b502a with its GERP bigwig and human ancestor FASTA, gnomAD and 1000 Genomes frequencies through the VEP cache, and SIFT with PolyPhen through VEP.

## Environment

Annotation and mining run on one CPU node. Step 2 used 32 cores and about five hours for 5 million records. The two sequence-model lanes each used one B200-class GPU, about 3 GPU-hours combined.

Software. Python 3.12, plus the core package list in `pipeline/requirements.txt`. Shell tools are curl, git, gzip, perl, and tabix with bgzip from htslib. VEP ran as the public container ensemblorg/ensembl-vep release_112.0. The GPN-MSA lane needs a CUDA build of PyTorch.

## Reproduction steps

All steps assume a fresh working directory.

```bash
export DATA_DIR=/path/to/workdir          # receives every input and output
export HF_TOKEN=...                       # token accepted on the gated dataset
pip install -r pipeline/requirements.txt
```

1. Download and QC. Fetches the proband VCF and phenotype document, checks gzip integrity, and writes a record-count and header QC report.

   ```bash
   bash pipeline/jobs/01_download_qc.sh
   ```

2. VEP annotation. Downloads the VEP cache, reference FASTA, ClinVar, AlphaMissense, and LOFTEE inputs, runs a 30,000-record sanity slice with a column and population check, then annotates the full VCF with SIFT, PolyPhen, LOFTEE, AlphaMissense, ClinVar, gnomAD, and 1000 Genomes.

   ```bash
   bash pipeline/jobs/02_vep_annotate.sh
   ```

3. HPO gene sets. Extract the HPO term identifiers from the gated phenotype document into a text file with one ID per line, then build the term-to-gene index from the public HPO table. This recreates the input that candidate matching uses, with no case data in the repository.

   ```bash
   python3 pipeline/make_hpo_gene_sets.py --terms hpo_terms.txt
   ```

4. Candidate mining. Joins VEP rows to genotypes with an allele-exact key, audits the join, assigns rarity tiers at maximum AF 0.001 and 0.01, grades damage classes, requires FILTER PASS on canonical chromosomes, crosses the HPO index, groups genes into compound-heterozygous and homozygous patterns, and writes the model-lane scoring inputs plus inert synonymous controls.

   ```bash
   python3 pipeline/jobs/03_mine_candidates.py
   ```

5. GPN-MSA lane, one GPU. Installs the public gpn package at the pinned commit and scores every surviving SNV over 512-bp windows from the streamed alignment store.

   ```bash
   pip install "git+https://github.com/songlab-cal/gpn.git@690557d949309cf4f4234554888bb5421c49aede"
   pip install transformers==5.15.0 safetensors==0.8.0 jaxtyping==0.3.11 cyclopts "zarr<3"
   bash pipeline/jobs/04_gpn_msa_score.sh
   ```

   Three optional checks validate this lane before a long run. `pipeline/gpn_msa_baseline_check.py` reproduces the published baseline logits on the upstream chr6 fixture to 1e-4 tolerance. `pipeline/gpn_stream_check.py` verifies that streamed store reads match the fixture tokens exactly. `pipeline/gpn_cli_check.py` scores synthetic variants end to end through the command line interface.

6. EVEE lane, one GPU, external code. Described in the next section.

7. Combine and rank. Merges annotation evidence, gene-level pattern bonuses, and the model-lane percentiles into one transparent weighted score, pairs compound-heterozygous genes first, applies the documented diagnostic override, and writes the ranked board plus the submission input table.

   ```bash
   python3 pipeline/jobs/05_combine_rank.py
   ```

8. Submission writer. Emits the exact challenge schema, with at most 10 rows, chromosome names prefixed with chr, monotone non-increasing EPCR values, and complete second-variant fields on every pair row.

   ```bash
   python3 pipeline/write_submission.py \
     "$DATA_DIR/mva-track1/ranking/submission_input.parquet" \
     submission/TEAMNAME_bub1b-comphet.csv
   ```

9. Report builder. `pipeline/plot_report_figs.py` regenerates the two report figures, the candidate funnel and the GPN-MSA candidate-versus-control distribution, into `figures/`. The methods narrative for the challenge upload is `submission/TEAMNAME_track1_report.md`.

## The EVEE lane

EVEE scoring used the public code and released probe weights from github.com/goodfire-ai/evee-manuscript at commit e5e4a43d60bb217819144c09ab509ed44f49cd5f, run by us during the analysis. The scoring recipe per variant builds reference and alternate windows of 16,384 bp centered on the variant, in both reading directions. The paper used 65,536 bp windows, and the reduction is recorded in the methods report. Each sequence passes through Evo 2 7B, and the block-27 hidden states are captured. For each direction, the 256 positions with the largest cosine divergence between variant and reference feed the released covariance probe, which returns a pathogenicity probability. The probe weights come from that repository at artifacts/samples/probe.safetensors, 4,279,352 bytes, byte-verified in our run. The lane writes `evee/evee_scores.parquet` with one `evee_score` per candidate SNV, which step 7 joins as an optional input. Before production scoring, the released probe was reproduced on the repository's 8-variant holdout set with 8 of 8 agreements.

## Validation performed in our run

- The 30,000-record sanity slice used the same image and flags as the full run and passed every column and plugin check.
- The variant-to-genotype join joined 99.04% of transcript rows. Unresolved rows are written to `join_unresolved.tsv`, and none survive rarity and quality filters into the candidate lanes.
- The GPN-MSA stack reproduced the published baseline logits to four decimal places, and streamed alignment-store reads matched the published fixture tokens exactly.
- The EVEE probe matched the released holdout predictions 8 of 8.
- The submission CSV passed schema validation, 10 rows, chr-prefixed chromosomes, monotone EPCR, and pair rows with complete second-variant fields.

## Acknowledgement

> This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

Reuse of the underlying individual-level data must follow the dataset terms. The data are gated and must be requested from SageBio/mva-hackathon-2026-data on Hugging Face. Publications must cite the dataset using the reference on the Hackathon Synapse page at the time of publication, and must not include information that could re-identify the data subject or their family beyond what the family has already made public. This repository itself contains no case data and no derived case tables. Released under CC BY 4.0, see `LICENSE`.
