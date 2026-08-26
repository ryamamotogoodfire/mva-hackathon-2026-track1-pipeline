# Track 2 report. Drug repurposing for a compound-heterozygous BUB1B loss-of-function case

## Glossary

- **Allele**. One of the two copies of a gene. Every person carries two copies of BUB1B.
- **Compound heterozygous**. Two different damaging variants, one on each copy of the same gene.
- **BUB1B / BUBR1**. The gene, and its protein product. BUBR1 runs the mitotic spindle assembly checkpoint that keeps chromosome numbers correct during cell division.
- **MVA syndrome 1**. Mosaic variegated aneuploidy syndrome 1. An inherited disorder caused by damaging variants in both BUB1B copies. Patients develop cells with wrong chromosome numbers and have a raised risk of childhood tumors, including embryonal rhabdomyosarcoma, Wilms tumor, leukemia, and lymphoma.
- **Stop-gain / nonsense variant**. A DNA change that turns a normal amino-acid codon into a premature stop signal, here p.Leu737Ter.
- **Premature termination codon (PTC)**. The new stop signal produced by a nonsense variant.
- **Readthrough**. The rare event where the ribosome passes over a stop codon and keeps making protein. The stop codon UGA is passed most often.
- **Nonsense-mediated decay (NMD)**. A quality-control system that destroys mRNA transcripts carrying a premature stop, before the ribosome can translate them.
- **50-nucleotide rule**. NMD attacks a PTC-carrying transcript when the PTC sits more than about 50 nucleotides before the last exon-exon junction.
- **Exon / exon-exon junction**. Exons are the protein-coding segments of a transcript. A junction is the boundary where two exons are joined after splicing.
- **Missense variant**. A DNA change that swaps one amino acid for another, here p.Asn1002Lys.
- **Kinase domain / pseudokinase**. The part of BUBR1 (residues 766 to 1050) shaped like a kinase. BUBR1 is an unusual pseudokinase. Its catalytic activity is dispensable for error-free chromosome segregation, but the domain's shape and stability matter.
- **Aneuploidy**. A wrong chromosome count in a cell. One aneuploid division event per cell division is enough to raise long-term tumor risk.
- **Repurposing**. Finding a new use for a medicine that regulators already approve for another disease.
- **AD / GT**. Sequencing allele depths and genotype calls from the proband WGS in Track 1.
- **EVEE / GPN-MSA**. The two Track 1 sequence-model scorers used to rank the variants. EVEE is a covariance probe on genome-model embeddings. GPN-MSA is a multi-species sequence model.
- **Secondary prevention**. Reducing the risk of further tumors after a first tumor has already occurred. This matches the proband record, which already lists rhabdomyosarcoma.

## The answer in one paragraph

The two BUB1B alleles in this child both destroy BUBR1 dosage, by two different routes. The stop-gain allele produces a transcript that nonsense-mediated decay destroys 748 nucleotides before the last exon-exon junction, far beyond the 50-nucleotide rule, so there is almost no mRNA left for any readthrough drug to act on. The missense allele escapes decay but swaps a folded, high-confidence residue inside the pseudokinase domain, a class of MVA mutations known to cut BUBR1 protein levels through instability. Because every approve-today readthrough medicine fails on the first allele and nothing approved refolds the second allele, **no approved medication restores BUBR1 function in this case**. The one defensible approved-medication proposal is **metformin as a secondary-prevention hypothesis**, offered because its pediatric safety record is unmatched among candidates and its mechanism chain is testable and falsifiable on the child's own cells before any clinical use. **Amlexanox** is named as the one mechanism-faithful restore-class option, at bench grade only. Surveillance remains the standard of care and the comparator every proposal must beat.

## The case and the two alleles

Track 1 ranked a compound-heterozygous BUB1B pair first in the proband WGS. Every coordinate below is identical to Track 1, build GRCh38.

| Allele | Location | HGVS | Genotype | Evidence from Track 1 |
| --- | --- | --- | --- | --- |
| p.Leu737Ter | chr15:40209701 | NM_001211.6 c.2210T>G | het 0/1, AD 21/25 | stop_gained, LOFTEE high-confidence loss of function, ClinVar Pathogenic/Likely pathogenic (VCV000533901.9, rs759242053), max AF 7.35e-5, EVEE 0.9762, GPN-MSA LLR -5.82 |
| p.Asn1002Lys | chr15:40220612 | NM_001211.6 c.3006T>G | het 0/1, AD 15/13 | missense, SIFT deleterious 0.01, PolyPhen probably damaging 0.997, never observed in population databases, EVEE 0.4250, GPN-MSA LLR -8.89 |

ClinVar cites two submitters. Labcorp/Invitae (SCV000762865.8) predicts a premature stop with an absent or disrupted protein and warns that gnomAD frequency data at this position are unreliable because of poor data quality. GeneDx (SCV004031025.1) predicts truncation or nonsense-mediated decay in a gene where loss of function is a known disease mechanism. Both views match the arithmetic below. Rarity statements in this report do not lean on gnomAD at this position.

The proband phenotype already lists rhabdomyosarcoma. Every prevention argument below is a secondary-prevention argument.

## Variant mechanism characterization

### Both alleles are loss of function, by different routes

BUB1B sits on the plus strand of chromosome 15 with 23 coding exons (transcript ENST00000287598, MANE-matched to NM_001211.6).

**The stop-gain allele.** Sequence verification on the MANE-matched CDS confirms c.2210 is the middle base of codon 737, changing TTA (Leu) to TGA, a UGA stop. UGA is the stop class with the highest natural readthrough propensity. The downstream nucleotide makes the tetranucleotide TGAA, a weaker readthrough context than the strongest known context TGAC (Dabrowski 2015, PMID 26176195, and Floquet 2012, PMID 22479203). The decisive fact is structural. The premature stop falls in exon 17 of 23, 748 nucleotides before the last exon-exon junction at c.2957, with six downstream exon-exon junctions. The 50-nucleotide NMD rule (Le Hir 2001, PMID 11532962, and the PTC-escape rules of Lindeboom 2016, PMID 27618451) predicts decay with a margin of 698 nucleotides. Both ClinVar submitters predict the same without computing it. A surviving, un-degraded transcript would lose the last 314 residues, the entire pseudokinase domain, but surviving is the exception, not the rule (Linde 2007, PMID 17290305, shows response to readthrough drugs tracks with how much nonsense transcript escapes NMD).

**The missense allele.** c.3006T>G verified against the CDS swaps codon 1002 from AAT (Asn) to AAG (Lys). The residue sits in exon 23 of 23, the last exon, so this transcript escapes NMD and yields a full-length protein. UniProt (O60566) places the protein kinase domain at residues 766 to 1050 with the catalytic reference at 882, so Asn1002 lies inside the domain. The precomputed AlphaFold model gives residue 1002 a confidence of 91.1 pLDDT, a folded core position, consistent with the severe in silico calls. Two further findings sharpen the reading. Vertebrate BUBR1 is an unusual pseudokinase whose catalytic activity is dispensable for error-free segregation (Suijkerbuijk 2012, PMID 22698286), so a catalysis story alone is the wrong frame. MVA-missense mutations instead cut BUBR1 abundance through instability (Suijkerbuijk 2010, PMID 20516114), which fits a folded-domain substitution at a constrained core residue.

**The shared endpoint is dosage.** MVA itself is the natural experiment. In patients with biallelic mutations a missense allele pairs with a truncating allele, exactly this proband's configuration, and patient cells show low overall BUBR1 abundance with impaired checkpoint function, rescued by adding back BUBR1 (Suijkerbuijk 2010, PMID 20516114).

### The dosage threshold sets the grading scale

Mouse genetics supplies the dose-response ladder every candidate is graded against.

| Residual BUBR1 level | Consequence in mice |
| --- | --- |
| about 10 percent, compound hypomorphs | progressive aneuploidy plus severe progeroid phenotypes (Baker 2004, PMID 15208629) |
| about 50 percent, heterozygotes | higher carcinogen-driven tumor development (Dai 2004, PMID 14744753) |
| above normal, transgenic overexpression | protection from aneuploidy-driven cancer, even under oncogenic Ras (Baker 2013, PMID 23242215) |

The protective boundary sits at or above the heterozygote level. Any restore therapy must lift functional BUBR1 back to roughly half-normal before benefit can begin. That arithmetic eliminates every approved readthrough medicine, whose best realistic yields are single-digit percents.

### Downstream biology

Low BUBR1 degrades the spindle assembly checkpoint. Cells mis-segregate chromosomes, arrest, senesce, and signal their own clearance by the immune system (Santaguida 2017, PMID 28633018). Cells that escape that clearance accumulate karyotype chaos and, over years, produce the MVA tumor spectrum of embryonal rhabdomyosarcoma, Wilms tumor, leukemia and lymphoma (Hanks 2004, PMID 15475955, García-Castillo 2008, PMID 18548531, and the matching case in Plaja 1999/Howard 1999, PMID 9916837). The pathway disrupted is mitotic chromosome segregation. The downstream consequence is constitutional aneuploidy with tumor predisposition. Secondary trisomies and whole-chromosome gains define the tissue of origin risk.

## Candidate space and grading method

Candidate classes were screened systematically. For every agent the screen recorded the mechanism tie to the two alleles or to aneuploidy biology, the strength of preclinical evidence, current approval status verified on the regulator's own page, pediatric labeling verified on FDA label text through openFDA, and fit to the MVA tumor spectrum. Every retained claim carries a PMID, DOI, or label citation. Approval claims name a jurisdiction. None of the work below uses proprietary data.

### Grading table

| Class | Representative agents | Mechanism tie | Evidence | Regulatory | Pediatric | Verdict |
| --- | --- | --- | --- | --- | --- | --- |
| Readthrough small molecule | ataluren | strong for the PTC | direct human trials, negative | EU lapsed 2025, never US | was labeled age 2+ in DMD | Rejected |
| Readthrough antibiotic | gentamicin, amikacin | strong for the PTC | direct human | approved as antibiotic only | labeled, toxicity-hostile | Rejected |
| NMD inhibition plus readthrough | amlexanox | strong for this allele | patient cells only | approved topical US, oral Japan | none | Fallback, research only |
| Energy stress, AMPK axis | metformin | moderate, indirect | indirect | approved, US children 10 to 16 | labeled for T2D | **Primary** |
| mTOR inhibition | sirolimus, everolimus | wrong direction | tolerance data only | approved | labeled with malignancy warning | Rejected |
| Checkpoint kinase inhibition | Mps1/TTK, Aurora B inhibitors | wrong direction | harm signals | none approved | none | Rejected |
| NSAID chemoprevention | aspirin | precedential only | RCT positive in Lynch | approved OTC | contraindicated | Rejected |
| Proteostasis pressure | bortezomib, Hsp90 inhibitors | moderate | patient cells | mixed | prophylaxis-incompatible | Rejected |
| Surveillance | imaging and bloods cadence | n/a | guideline frameworks | n/a | n/a | Comparator, remains standard |

### Why the readthrough class fails for this allele

Ataluren is the class's best case, and the disabling fact is NMD, not its chemistry and not the stop context. The allele's transcript is NMD-targeted with a 698-nucleotide margin past the canonical threshold, so the substrate for readthrough barely exists (PMID 17290305). The stop context itself is permissive, since UGA is the readthrough-friendliest stop class and the TGAA tetranucleotide is weaker only than TGAC (Dabrowski 2015, PMID 26176195, and Floquet 2012, PMID 22479203). The context is therefore not what blocks this class. The compound's own clinical record compounds the verdict. Two phase 3 trials read negative, nonsense cystic fibrosis (Kerem 2014, PMID 24836205) and nonsense Duchenne, whose primary endpoint missed prespecified significance (McDonald 2017, PMID 28728956), and the European Medicines Agency did not renew Translarna's conditional authorisation because effectiveness was not confirmed, with expiry on 28 March 2025 (EMA Translarna EPAR page). It was never approved in the United States. The mechanism is real (Welch 2007, PMID 17450125). This allele plus this drug simply does not work, and the reason is transcript decay upstream of translation.

Aminoglycosides add ototoxicity and nephrotoxicity in children to the same structural objections (Wilschanski 2003, PMID 14534336, Howard 2000, PMID 10939566).

### Why amlexanox survives to bench grade

Amlexanox is the only approved molecule that addresses the actual blocker. It stabilizes nonsense-containing mRNA by inhibiting NMD and it promotes readthrough of the PTC, yielding full-length functional protein in patient-derived cells (Salvatori 2012, PMID 22938201, Atanasova 2017, PMID 28549954). It is approved as a 5 percent topical paste in the United States (ChEMBL first approval 1996, Aphthasol) and as an oral agent in Japan (Dabrowski 2018 review, PMID 30134808). The evidence stops at patient cells. There is no systemic nonsense-disease trial, no pediatric systemic safety record, and in-vitro NMD-modulating doses exceed topical exposure. The verdict is a proband-fibroblast proof-of-mechanism experiment, falsifiable in weeks, offered as research only. It cannot touch the missense allele.

### Why metformin is the primary candidate

The protect axis scores candidates on whether they reduce the fitness of newly aneuploid cells without harming the child. Metformin clears the safety bar that eliminates every rival. The FDA label establishes safety and effectiveness for immediate-release metformin in type 2 diabetes for children 10 to 16 (openFDA metformin record), with the TODAY trial covering ages 10 to 17 (PMID 22540912). Chronic risks are few and monitorable, mainly vitamin B12 deficiency over years (DPPOS study, PMID 26900641) with the renal-function contraindication handled by monitoring. The mechanism tie is real but indirect. Aneuploid cells live under energy, proteotoxic, and lysosomal stress (Santaguida 2015, PMID 26404941). The energy-stress agent AICAR selectively kills aneuploid cells (Tang 2011, PMID 21315436), and metformin is the approved, pediatric-tested member of that AMPK-activating neighborhood. Observational cancer-incidence signals in adult diabetes cohorts are consistent but confounded (Decensi 2010, PMID 20947488). The proposal is stated as a hypothesis with falsifiable bench tests, below, and only as secondary prevention under the oncology team.

### Rejected classes, recorded

- **mTOR inhibitors.** Mechanistically they point the wrong way. Rapamycin-class drugs relieve protein-synthesis burden and induce autophagy, removing exactly the bottleneck that makes aneuploid cells unfit (PMID 26404941). Aneuploid mis-segregated cells are partly cleared by p53-driven and immune surveillance (PMID 28633018), and the sirolimus label itself warns that immunosuppression raises infection and lymphoma risk (openFDA sirolimus record, boxed warning). Pediatric labels cover transplant use at age 13 and older for sirolimus and TSC-associated SEGA from age 1 for everolimus (PMID 21047224), which documents tolerance, not prevention.
- **Checkpoint kinase inhibitors.** Removing the residual checkpoint in a checkpoint-hypomorphic child raises mis-segregation further. Checkpoint blockade in human cells produces massive chromosome loss (Kops 2004, PMID 15159543). Their only defensible role is against established aneuploid tumors (Cohen-Sharir 2021, PMID 33505028), a treatment setting this child is not in. None is approved (Mason 2017, PMID 28270606). They are recorded to mark the boundary of the repurposing space.
- **Aspirin.** CAPP2 is the strongest randomized chemoprevention precedent in any hereditary instability syndrome (Burn 2011, PMID 22036019, Burn 2020, PMID 32534647), but Lynch syndrome is mismatch-repair driven, the benefit accrued in adults, the MVA spectrum is not colorectal, and aspirin is contraindicated in children because of Reye's syndrome risk.
- **Proteostasis pressure.** The aneuploidy-proteostasis dependence is real (Donnelly 2014, PMID 25205676), but no Hsp90 inhibitor is approved and bortezomib cannot be given to a well child as prophylaxis.

### Surveillance remains the comparator

No MVA-specific surveillance guideline exists. Management follows childhood cancer predisposition frameworks (Kratz 2017, PMID 28168833) adapted to the MVA spectrum. Any proposed candidate must improve on this on a risk-adjusted basis, and no graded candidate above clearly does.

## Falsifiable statements and validation plan

The report's claims are phrased so a bench team can falsify them.

Restore chain.
- R1. BUBR1 protein in proband fibroblasts sits far below family-control level on Western blot, matching transcript loss from the PTC allele plus instability of the missense allele. This falsifies or confirms the shared-dosage model.
- R2. Amlexanox at tolerated concentrations raises PTC-transcript abundance by allele-specific qPCR before any full-length protein appears. Predicted protein recovery is below 5 percent for a TGAA context and fails the half-normal threshold, but the experiment decides.
- R3. Complementing proband cells with wild-type versus kinase-domain-deletion BUBR1 separates scaffold from catalytic loss, characterizing p.Asn1002Lys directly.

Protect chain.
- P1. Serially passaged proband fibroblasts accumulate lagging chromosomes and micronuclei faster than matched family controls, confirming a measurable intrinsic aneuploid-prone state in non-tumor tissue.
- P2. Metformin at exposures achievable in children reduces outgrowth of the aneuploid subfraction relative to euploid sibling cultures, quantified with karyotype-resolved colony counts. A null finding falsifies the prevention hypothesis before any clinical use.

## Innovation and scalability

The innovation is methodological honesty plus a bench-ready loop. Every allele fact in this report is recomputed from public primary sources at fixed coordinates, every candidate row is traced to a regulator's own page or a PMID, and every rejection is recorded with its reason. The proposals come with named falsifying experiments. The work scales because it is CPU-only, scripted end to end, and reusable for any compound-heterozygous loss-of-function case. The same pipeline flags readthrough-hostile contexts before a team spends on them.

Preclinical directions offered for completeness are suppressor-tRNA readthrough matched to UGA, splice-modulation of exon 17, prime editing correction of c.2210T>G (a transversion, outside base-editor reach), and BUBR1 cDNA augmentation (the coding sequence is about 3.2 kb, inside AAV capacity).

## Completeness disclosures

Track 1's mechanical ranking crowned a homozygous 45-base PEX5 deletion (chr12:7190513) above the BUB1B pair. That finding is complete in Track 1's report and does not enter this mechanism chain. Parental segregation has not been performed, so the compound-heterozygous trans configuration is inferred, not proven, and parental testing is the first preceptor action. Both submitters flagged poor gnomAD data quality at the stop-gain position, and all rarity statements here rest on ClinVar plus Track 1's own mining.

## Limitations

- Dosage-threshold reasoning rests on mouse genetics, not human titration.
- Secondary-prevention logic is untested for metformin anywhere, and pediatric cohort data are in type 2 diabetes, not oncology prevention. If the child is younger than 10, the labeled window does not yet cover them and pharmacokinetic bridging would be required.
- Readthrough-context grades come from reporter systems, not from this exon 17 context measured directly.
- The NMD rule is a rule, not a direct measurement of the proband's transcripts.
- The trans configuration is inferred.

## Track 2 methods description form, answered

**Team name.** Silico EVEE

**Approach in detail.** From the Track 1-ranked BUB1B pair, the mechanism map was rebuilt from public primary sources. Ensembl REST supplied the 23-exon structure of the MANE-matched transcript, the CDS was fetched and both codon changes verified by translation, the NMD arithmetic computed from exon coordinates, UniProt supplied domain annotation, and AlphaFold DB supplied per-residue confidence. Candidate classes were screened with fixed Europe PMC queries, each agent's approval state checked on ChEMBL, openFDA label text, and the EMA EPAR page, and classes graded on five axes. Recommendations keep only what survives falsifiable-statement drafting.

**Automated or manual.** Hybrid. Allele mapping and retrieval are scripted and rerunnable. Grading and synthesis were curated by review of the retrieved abstracts, with every retained claim citation-tagged.

**Manual curation detail.** Each candidate class was assessed against abstracts fetched from Europe PMC, and any claim entering the report was matched to the abstract text. Regulatory claims were matched to the label or EPAR text, not to secondary databases.

**Public data only.** Yes. ClinVar, Ensembl REST, RefSeq/UniProt, AlphaFold DB, Europe PMC, ChEMBL, openFDA/DailyMed, EMA EPAR. The only non-public input is the proband WGS-derived Track 1 result, which the hackathon supplied.

**Public sources in detail.** ClinVar VCV000533901.9 for classification and submitter comments. Ensembl REST lookup/map/sequence endpoints for transcript, exons, and codons. UniProt O60566 for domains. AlphaFold model AF-O60566-F1 for confidence. Europe PMC/PubMed for every mechanistic citation. ChEMBL REST for molecule approval metadata. openFDA drug label records for pediatric sections. EMA EPAR for Translarna history.

**Proprietary sources.** None.

**Variant mechanism characterization.** Both alleles are loss of function. The stop-gain allele subjects BUB1B to NMD, verified by codon and exon arithmetic. The missense allele substitutes a folded residue in the pseudokinase domain, class-consistent with instability-driven abundance loss. The disrupted pathway is the spindle assembly checkpoint. The downstream consequence is constitutional mosaic variegated aneuploidy with childhood tumor predisposition.

**Time and effort.** About one analyst-day of computational work plus literature curation within a 36-hour wall-clock window. CPU only. Re-runnable in under an hour.

**Method abstract.** Starting from the Track 1 compound-heterozygous BUB1B pair, the variant mechanism was characterized from public primary sources rather than asserted. The stop-gain allele was verified as a UGA PTC in a TGAA context, in exon 17 of 23, 748 nucleotides before the last exon junction, a margin of 698 nucleotides past the 50-nucleotide NMD rule, so the transcript is destroyed before translation. The missense allele is a last-exon change inside the pseudokinase domain at a high-confidence folded residue, matching the MVA class of instability-driven abundance loss. Candidate classes were screened and graded against regulator pages and PMIDs. The restore class fails because NMD removes the readthrough substrate, with the permissive stop context not at fault, and ataluren's EU authorisation lapsed in March 2025 after unconfirmed effectiveness. mTOR inhibition was rejected because it relieves the proteostasis bottleneck of aneuploid cells and carries a malignancy warning. Checkpoint kinase inhibitors were recorded as wrong-direction. Metformin is offered as the primary secondary-prevention hypothesis on pediatric safety strength and an indirect aneuploidy-selectivity chain, and amlexanox as a bench-only restore test because it uniquely combines NMD inhibition with readthrough. Strengths include complete source discipline, falsifiable statements, and negative results recorded with reasons. Limitations include mouse-derived thresholds, reporter-based context grades, and an inferred trans configuration. The recommendation is surveillance plus the named bench experiments, with metformin only under oncology stewardship.

## References

1. Hanks S, et al. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B. Nat Genet 2004. PMID 15475955.
2. García-Castillo H, et al. Clinical and genetic heterogeneity in mosaic variegated aneuploidy. Am J Med Genet A 2008. PMID 18548531.
3. Plaja A, et al. Child with mosaic variegated aneuploidy and embryonal rhabdomyosarcoma. PMID 9916837.
4. Yost S, et al. Biallelic TRIP13 mutations predispose to Wilms tumor. Nat Genet 2017. PMID 28553959.
5. Suijkerbuijk SJ, et al. Molecular causes for BUBR1 dysfunction in MVA. Cancer Res 2010. PMID 20516114.
6. Suijkerbuijk SJ, et al. The vertebrate mitotic checkpoint protein BUBR1 is an unusual pseudokinase. Dev Cell 2012. PMID 22698286.
7. Baker DJ, et al. BubR1 insufficiency causes early onset of aging-associated phenotypes. Nat Genet 2004. PMID 15208629.
8. Dai W, et al. Slippage of mitotic arrest and enhanced tumor development in BubR1 haploinsufficiency. Cancer Res 2004. PMID 14744753.
9. Baker DJ, et al. Increased expression of BubR1 protects against aneuploidy and cancer. Nat Cell Biol 2013. PMID 23242215.
10. Le Hir H, et al. The exon-exon junction complex in mRNA export and NMD. EMBO J 2001. PMID 11532962.
11. Lindeboom RG, et al. The rules and impact of NMD in human cancers. Nat Genet 2016. PMID 27618451.
12. Linde L, et al. NMD governs response of CF patients to gentamicin. J Clin Invest 2007. PMID 17290305.
13. Dabrowski M, et al. Translational readthrough potential of natural termination codons. RNA Biol 2015. PMID 26176195.
14. Floquet C, et al. Statistical analysis of readthrough levels in mammalian cells. PLoS Genet 2012. PMID 22479203.
15. Howard MT, et al. Sequence specificity of aminoglycoside-induced stop codon readthrough. Ann Neurol 2000. PMID 10939566.
16. Welch EM, et al. PTC124 targets genetic disorders caused by nonsense mutations. Nature 2007. PMID 17450125.
17. Kerem E, et al. Ataluren for nonsense-mutation cystic fibrosis. Lancet Respir Med 2014. PMID 24836205.
18. McDonald CM, et al. Ataluren in ACT DMD phase 3. Lancet 2017. PMID 28728956.
19. Wilschanski M, et al. Gentamicin-induced correction of CFTR function. Am J Respir Crit Care Med 2003. PMID 14534336.
20. Salvatori F, et al. Rescue of nonsense mutations by amlexanox in human cells. Orphanet J Rare Dis 2012. PMID 22938201.
21. Atanasova VS, et al. Amlexanox enhances COL7A1 PTC read-through. J Invest Dermatol 2017. PMID 28549954.
22. Dabrowski M, et al. Advances in drug-stimulated translational readthrough. Trends Mol Med / review series 2018. PMID 30134808.
23. Keeling KM, et al. Therapeutics based on stop codon readthrough. Annu Rev Genomics Hum Genet 2014. PMID 24773318.
24. Tang YC, et al. Identification of aneuploidy-selective antiproliferation compounds. Cell 2011. PMID 21315436.
25. Santaguida S, et al. Aneuploidy-induced cellular stresses limit autophagic degradation. Genes Dev 2015. PMID 26404941.
26. Santaguida S, et al. Chromosome mis-segregation generates arrested cells eliminated by surveillance. Dev Cell 2017. PMID 28633018.
27. Donnelly N, et al. HSF1 deficiency and impaired HSP90 folding in aneuploid human cells. EMBO J 2014. PMID 25205676.
28. Cohen-Sharir Y, et al. Aneuploidy renders cancer cells vulnerable to checkpoint inhibition. Nature 2021. PMID 33505028.
29. Kops GJPL, et al. Lethality through massive chromosome loss by checkpoint inhibition. PNAS 2004. PMID 15159543.
30. Mason JM, et al. CFI-402257 Mps1 inhibitor. PNAS 2017. PMID 28270606.
31. Krueger DA, et al. Everolimus for SEGA in TSC. N Engl J Med 2010. PMID 21047224.
32. Franz DN, et al. Long-term everolimus in TSC, EXIST-1 final. PLoS One 2016. PMID 27351628.
33. Zeitler P, et al. TODAY trial, metformin in youth with type 2 diabetes. N Engl J Med 2012. PMID 22540912.
34. DPPOS. Long-term metformin and B12 deficiency. J Clin Endocrinol Metab 2016. PMID 26900641.
35. Decensi A, et al. Metformin and cancer risk meta-analysis. Cancer Prev Res 2010. PMID 20947488.
36. Burn J, et al. Long-term aspirin in Lynch syndrome. Lancet 2011. PMID 22036019.
37. Burn J, et al. CAPP2 10-year follow-up. Lancet 2020. PMID 32534647.
38. Kratz CP, et al. Childhood cancer predisposition syndromes, review and recommendations. Am J Med Genet A 2017. PMID 28168833.
39. ClinVar record VCV000533901.9, NM_001211.6 c.2210T>G, with submissions SCV000762865.8 and SCV004031025.1.
40. EMA Translarna EPAR page, conditional authorisation 2014-07-31, non-renewal decision 2025-03-28. https://www.ema.europa.eu/en/medicines/human/EPAR/translarna
41. openFDA drug label records for metformin hydrochloride tablets, sirolimus, everolimus, aspirin.
42. Ensembl REST endpoints for ENST00000287598 lookup, CDS map, and sequence.
43. UniProt record O60566, BUB1B_HUMAN.
44. AlphaFold DB model AF-O60566-F1.
