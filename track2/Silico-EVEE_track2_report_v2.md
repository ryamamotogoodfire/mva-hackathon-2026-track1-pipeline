# Track 2 report, version 2. Quantitative full-FDA screen for drug repurposing in a compound-heterozygous BUB1B loss-of-function case

Team: Silico EVEE. Date: 2026-08-27. This document supersedes the v1 Track 2 report. Allele facts, coordinates, and mechanism characterization from v1 are unchanged (verified against the same public primary sources); what changed is the candidate grading, now grounded in a quantitative screen of the full approved-drug space plus structural allele tickets, and the verdict updates that follow from it.

## Glossary

**SAC** the spindle assembly checkpoint, the mitotic surveillance system BUBR1 belongs to. Low BUBR1 weakens it and raises chromosome missegregation.
**MVA / MVA1** mosaic variegated aneuploidy syndrome; MVA1 (Disease Ontology DOID:0080141) is the form caused by BUB1B variants, the proband diagnosis.
**NMD** nonsense-mediated decay, the cell's surveillance pathway that degrades transcripts carrying a premature stop codon well before the last exon junction.
**CIN** chromosomal instability, a high ongoing rate of chromosome missegregation.
**ddG** the free-energy change of protein folding caused by a mutation, in kcal/mol; positive means destabilizing.
**tau (connectivity)** a -100..+100 score of how much a drug's transcriptomic signature resembles a query signature; negative (reversal) is the rescue direction.
**PPI** protein-protein interaction; the proximity lane measures mean PPI distance from a drug's targets to the MVA gene module.
**KG** knowledge graph; here the Bioteque embedding space of compounds, genes, and disease nodes.
**PTC** premature termination codon. **UGA** the RNA spelling of the proband's stop codon.

## The answer in one paragraph

Both BUB1B alleles destroy BUBR1 dosage, and no approved medication restores BUBR1 in this case, exactly as in v1. What the v2 quantitative pass changes is the secondary-prevention ranking. The full approved-drug space (1,811 drugs) was screened on three orthogonal quantitative lanes; no drug tops all three lanes and, under the pre-declared rule that screen scores alone can never promote, nothing is promoted by the screen. On direct published evidence the proposal moves from metformin to the **chloroquine/hydroxychloroquine class as the primary secondary-prevention hypothesis**: chloroquine carries the strongest published aneuploidy-selective vulnerability data for this direction (Tang 2011), and hydroxychloroquine's FDA label establishes pediatric safety and effectiveness for malaria (with the label's 31-kg film-coated-tablet restriction called out). **Metformin is demoted** to the protection fallback because its published aneuploidy-selectivity evidence is incidental, markedly weaker than the chloroquine-class signal in the same paper, and unreplicated. **NMN** enters the dossier as the **best biological candidate for raising BUBR1** (the only published in-vivo increase of BubR1 abundance, North 2014) and is **status-blocked**: it is a lawful dietary supplement, not an FDA-approved medication. **Amlexanox** remains the one mechanism-faithful restore-class option at bench grade only. Surveillance stays the standard of care and the comparator every proposal must beat.

## The case and the two alleles

Unchanged from v1, quoted for completeness. Track 1 ranked a compound-heterozygous BUB1B pair first in the proband WGS, GRCh38.

| Allele | Location | HGVS | Genotype | Evidence from Track 1 |
| --- | --- | --- | --- | --- |
| p.Leu737Ter | chr15:40209701 | NM_001211.6 c.2210T>G | het 0/1, AD 21/25 | stop_gained, LOFTEE high-confidence loss of function, ClinVar Pathogenic/Likely pathogenic (VCV000533901.9, rs759242053), max AF 7.35e-5 |
| p.Asn1002Lys | chr15:40220612 | NM_001211.6 c.3006T>G | het 0/1, AD 15/13 | missense, SIFT deleterious 0.01, PolyPhen probably damaging 0.997, never observed in population databases |

The proband phenotype already lists rhabdomyosarcoma, so every prevention argument is secondary prevention.

## Variant mechanism characterization

### Both alleles are loss of function, by different routes (v1, verified again in v2)

**The stop-gain allele.** c.2210T>G changes codon 737 from TTA (Leu) to TGA, a UGA stop in exon 17 of 23. UGA is the readthrough-friendliest stop class and the TGAA tetranucleotide is a permissive context, weaker only than TGAC (Dabrowski 2015, PMID 26176195; Floquet 2012, PMID 22479203). The decisive fact is structural. In v2 the arithmetic was recomputed live from Ensembl (ENST00000287598): the PTC codon starts at CDS 2209, the last exon-exon junction sits at CDS 2957, the margin is 748 nucleotides with six downstream junctions (2284 through 2957), 698 nucleotides past the 50-nucleotide NMD rule. The transcript is destroyed before translation (Le Hir 2001, PMID 11532962; Lindeboom 2016, PMID 27618451); a surviving transcript would lack the whole pseudokinase domain, but survival is the exception (Linde 2007, PMID 17290305).

**The missense allele.** c.3006T>G swaps codon 1002 from AAT (Asn) to AAG (Lys) in the last exon, so the transcript escapes NMD. UniProt O60566 places the kinase domain at residues 766 to 1050 and AlphaFold AF-O60566-F1-model_v6 gives the 990-to-1015 neighborhood a mean pLDDT of 91.5, a folded core position. Vertebrate BUBR1's catalytic activity is dispensable for error-free segregation (Suijkerbuijk 2012, PMID 22698286); MVA missense mutations instead cut BUBR1 abundance through instability (Suijkerbuijk 2010, PMID 20516114).

**New in v2, the structural ticket on the missense allele.** RaSP (Rapid Stability Prediction, Blaabjerg 2023, eLife 12:e82593) quantifies the instability. N1002K reads +2.41 kcal/mol, the same destabilization class as the experimentally validated mild-instability allele I909T (+3.28 kcal/mol) at the tool's roughly 2 kcal/mol practical resolution, and well below the validated extreme allele L1012P (+8.70). The instrument's identity-substitution band over its 1,050 identity substitutions is mean +0.39, sd 0.58 kcal/mol. fpocket 4.2.3 (Le Guilloux 2009, PMID 19459140) maps 97 pockets on the AlphaFold model, none with druggability above 0.169 overall or 0.111 inside the pseudokinase domain, and none centered at N1002 or the adjacent validated loci. The pocket-based pharmacological-chaperone path is therefore structurally closed for this residue.

**The shared endpoint is dosage**, and the dosage ladder is unchanged: about 10% residual BUBR1 gives the progeroid compound-hypomorph state (Baker 2004, PMID 15208629), about 50% (the heterozygote) raises carcinogen-driven tumor incidence (Dai 2004, PMID 14744753), and above-normal expression protects from aneuploidy-driven cancer (Baker 2013, PMID 23242215).

## The v2 quantitative screen over the full approved-drug space

The v1 candidate screen was qualitative, anchored on the literature. v2 adds a computable pass over the entire approved-drug space, so that the final ranking is not hostage to what the literature happened to cover.

### The screening set

ChEMBL 37 supplies every approved (max_phase 4) drug with a directly annotated human protein target: 1,811 drugs, 650 targets, 3,927 drug-gene rows, with known-pair spot checks passing (methotrexate–DHFR, imatinib–ABL1, warfarin–VKORC1, simvastatin–HMGCR, abaloparatide–PTH1R). openFDA enriches 1,151 of those drugs with US label text, 962 carrying a Pediatric Use section. The decision rule is pre-declared: a candidate promotes only with the strongest direct published evidence for its rescue direction plus pediatric feasibility, and screen scores alone never promote. The field's documented screen-only failure mode (a topiramate–IBD association that did not transfer clinically) is the cautionary template.

### Lane 1, network proximity

The protein-protein interaction network is BioGRID 4.4.236 human physical (19,908 proteins, 849,005 edges). The MVA disease module, declared before scoring, is BUB1B, BUB1, CEP57, TRIP13 (MVA founders) plus MAD1L1, MAD2L1, BUB3, TTK, CDC20, CENPE, AURKB, PLK1 (core SAC). Drug proximity is z-scored against 1,000 degree-matched nulls (Barabasi/Guney method, seed 42). The lane passes its published calibrations: ivacaftor–CFTR z=-3.41, imatinib–BCR/ABL1 z=-10.66, off-axis simvastatin–MVA z=+0.37. 1,792 drugs were scored.

### Lane 2, knowledge-graph embeddings

The Bioteque space (CPD-int-GEN-ass-DIS disgenet bundle; Bioteque as in Van Engelen 2023, PMID 36419033) holds 4,954 compounds and 3,836 disease nodes, and contains the MVA1 disease node DOID:0080141. Compounds score by cosine similarity to the node, as percentiles. The lane is calibrated against itself: it recovers its own metapath edges at rank AUC 0.917 and its top anchor methotrexate sits at percentile 96.8 for its true-indication node. 788 approved drugs mapped into this space.

### Lane 3, signature reversal over LINCS, with the instrument limit declared

The BUB1B-low proxy is built from 58 quality-gated BUB1B knockdown signatures in GEO GSE92742 Phase I (own z at most -1; consensus BUB1B z=-2.95), with 150 up- plus 150 down-regulated sign-consistent genes. It is scored against 205,034 compound signatures with weighted connectivity (Subramanian 2017, PMID 29195078), tau-normalized. Instrument checks: the consensus readout passes (z=-2.95); the published positive pair MTOR-knockdown to sirolimus reconnects with +19.1 mean connectivity; **the held-out-hairpin reconnection fails (0.101 versus the random-null p95 of 0.184)**, meaning Phase I BUB1B hairpins do not cohere. This lane is therefore instrument-limited all along: it yields per-drug aggregates but cannot promote, and its extreme values are weak evidence only.

### What the three lanes found

No ingredient ranks in the top 100 of all three lanes. The screen's single strongest approved-candidate beat is fostamatinib (knowledge-graph percentile 99.8, and the only approved compound with a direct metapath edge to the MVA1 disease node), yet its LINCS reversal read sits at background level (per-drug minimum tau -53.2, strong-reversal fraction 0.065 versus the lane background of about 0.070, on the validation-limited lane), and it has no published evidence in any rescue direction and no pediatric approval. It is rejected as screen-flag-only, exactly as the rule requires. Coverage is partial: 584 of 1,811 drugs mapped into the knowledge-graph lanes and 587 into the LINCS join.

### Structural direction check on the dosage-raise axis

The dosage-raising direction uniquely matches the dose-response logic of MVA. Its named route is the NAD+/SIRT2 axis of North 2014 (PMID 24825348), and it is not an approved-medication route at all: see the NMN verdict below.

## Candidate space and v2 grading

### v2 grading table

| Class | Representative agents | Mechanism tie | Evidence | Regulatory | Pediatric | v2 verdict |
| --- | --- | --- | --- | --- | --- | --- |
| Lysosomal/autophagic stress | chloroquine, hydroxychloroquine | strongest published aneuploidy-selective vulnerability (Tang 2011) | direct cell work, aneuploidy-selective | approved (malaria, RA/SLE) | HCQ label states s/e established for pediatric malaria; film-coated HCQ tablets usable only at >=31 kg; chloroquine 16.67 mg/ml liquid covers small children | **Promoted, primary secondary-prevention hypothesis** |
| Energy stress, AMPK axis | metformin | indirect, aneuploidy-selectivity incidental and weak in Tang 2011 (subordinate to AICAR) | no dedicated replication | approved, labeled in children 10-16 | labeled for T2D | Demoted, protection fallback |
| NAD+/SIRT2 axis (dosage raise) | NMN | only published in-vivo BubR1 raise | North 2014 (wild-type mice; the hypomorphic cohorts tested SIRT2 overexpression) | lawful supplement since FDA's 2025-09-29 reversal; not an approved medication | none, not a medication | Best biological, status-blocked |
| NAD+/SIRT2 axis, alternate | nicotinamide | direction-conflict (SIRT2 inhibitor at relevant concentrations) | none for BUBR1 | vitamin | n/a | Rejected, direction-conflict |
| Readthrough small molecule | ataluren | strong for the PTC | direct human trials, negative | EU lapsed 2025-03-28, never US | was labeled age 2+ in DMD | Rejected |
| Readthrough antibiotic | gentamicin, amikacin | strong for the PTC | direct human | approved as antibiotic only | labeled, toxicity-hostile | Rejected |
| NMD inhibition plus readthrough | amlexanox | strong for this allele | patient cells only | approved topical US, oral Japan | none | Fallback, research only (unchanged) |
| mTOR inhibition | sirolimus, everolimus | wrong direction | tolerance data only | approved | labeled with malignancy warning | Rejected |
| Checkpoint kinase inhibition | Mps1/TTK, Aurora B inhibitors | wrong direction | harm signals | none approved | none | Rejected |
| NSAID chemoprevention | aspirin, sulindac | precedential only | RCT positive in Lynch; pediatric FAP RCT negative | approved OTC | aspirin contraindicated | Rejected |
| Proteostasis pressure | bortezomib, HSP90 inhibitors | moderate but direction-hostile to an unstable client (HSP90 inhibition degrades mutant BUBR1) | patient cells | none approved as prophylaxis | prophylaxis-incompatible | Rejected |
| Chemical chaperones | 4-PBA, TUDCA, sapropterin | mechanism unmatched (folded-domain locus has no druggable pocket per fpocket) | disease-specific, failures in ALS for the pair | mixed/withdrawn | varying | Rejected |
| Screen-flag kinase inhibition | fostamatinib | none published in a rescue direction | knowledge-graph lane only (percentile 99.8, validation-limited) | approved adult ITP | none | Rejected, screen-flag only |
| Surveillance | imaging and bloods cadence | n/a | guideline frameworks | n/a | n/a | Comparator, remains standard |

### Why the hydroxychloroquine class is the primary candidate

The candidate must reduce the fitness of newly aneuploid cells without harming a well child. The lysosomal axis is the strongest published route for that: aneuploid cells rely on lysosomal degradation to bear mis-segregation burdens, and lysosome-focused compounds are selectively antiproliferative against them (Tang et al., Cell 2011, PMID 21315436). Hydroxychloroquine is the child-friendly approved member of that class: its FDA label states that safety and effectiveness are established in pediatric patients for malaria treatment and prophylaxis. The label also states the film-coated tablets cannot be crushed or divided and so cannot be given to children under 31 kg; MVA patients are frequently growth-restricted, so the on-label pediatric leg narrows for small children and the class route under 31 kg is chloroquine oral liquid (16.67 mg/ml, malaria-labeled), which is also the compound that carried the Tang evidence. Retinal risk is chronic-use, dose-limited, and monitorable under ophthalmology (hydroxychloroquine label). The proposal is a secondary-prevention hypothesis, off-label for MVA, owned by the oncology team, with the bench falsification below.

### Why metformin is demoted

Metformin's v1 primary rank rested on pediatric safety plus an indirect mechanism chain: aneuploid cells live under energy stress, and energy-stress agents kill them selectively. The re-read is less kind. In the Tang 2011 screen itself the aneuploidy-selective energy-stress agent is AICAR, not metformin, and metformin's own aneuploidy-selective action appears only incidentally, un-replicated, and markedly weaker. The v2 sieve found no dedicated replication of an aneuploidy-selective metformin effect. The verdict is demotion rather than replacement: with its unmatched pediatric record, metformin stays named as the protection fallback if the class pick fails at the bench.

### Why NMN is graded best biological and status-blocked

North 2014 (PMID 24825348) is the only published route that increases BubR1 abundance in vivo. In that work NMN given to young and aged wild-type mice raised BubR1 protein by boosting NAD+ and thereby reducing SIRT2-mediated deacetylation of BubR1; the paper's BubR1-hypomorphic cohorts tested SIRT2 overexpression, which is where the disease-model evidence lives. Regulatory status: FDA's 2022 drug-preclusion determination had excluded NMN from the dietary-supplement definition, but FDA reversed that exclusion on 2025-09-29 in its response to citizen petition FDA-2023-P-0872, so NMN is a lawful dietary supplement. A lawful supplement is not an FDA-approved medication, so NMN cannot enter the proposal table. Nicotinamide is rejected as a separate candidate: at relevant concentrations it inhibits SIRT2, the arc's effector, which flips the direction.

### Approved NAD+ precursors, and the coverage gap that hid them

The dosage-raise direction deserved a second look, because niacin is an FDA-approved drug and an NAD+ precursor, which makes it the approved-medication analog of the NMN route. Checking why it never appeared in the screen exposed a real limitation of the drug universe.

The screen's universe is target-annotated pharmacology. Nutrient-class agents that act as metabolic substrates rather than through a protein target carry no mechanism annotation, so the entire NAD+ precursor class was invisible to all three lanes. Their absence from the screen is a coverage gap, not a negative result, and it falls on exactly the class the North 2014 mechanism implicates.

- NIACIN (CHEMBL573, max_phase 4, first approval 1982): 0 mechanism records in ChEMBL 37, so no target to place in the network, the knowledge graph, or the target-based joins.
- NIACINAMIDE / nicotinamide (CHEMBL1140, max_phase 4): 0 mechanism records.
- ACIPIMOX (CHEMBL345714, max_phase 4): 1 mechanism record whose target_chembl_id is null and whose action_type is null with mechanism_of_action 'Unknown'. It survives into the approved-molecule table but is dropped by the join to a human protein target, so it never receives a lane score.
- NICOTINAMIDE RIBOSIDE (CHEMBL438497): max_phase 3, not an approved medication, and 0 mechanism records.
- NICOTINYL ALCOHOL (CHEMBL1235535): max_phase -1.

Graded against the North 2014 mechanism, on approval status, pediatric record, and mechanism tie:

| Agent | Approval | Pediatric record | Mechanism tie to raising BUBR1 | Verdict |
| --- | --- | --- | --- | --- |
| niacin (nicotinic acid) | FDA prescription (NIACOR ANDA040378 immediate-release; multiple extended-release ANDAs), for dyslipidemia | not established (IR: children and adolescents; ER: 16 years and under) | genuine but indirect NAD+ precursor; no BubR1, aneuploidy, or MVA data | rejected for this case, class analog recorded |
| nicotinamide | vitamin or unapproved-listing product, not an FDA-approved therapeutic drug | none | direction-conflict: inhibits SIRT2, the effector the BubR1 increase is attributed to | rejected, direction-conflict |
| acipimox | no FDA approval (no openFDA label or application); non-US approvals only | none | weak: HCAR2 lipolysis agent, not an efficient NAD+ precursor | rejected |
| nicotinamide riboside | dietary supplement, ChEMBL max_phase 3 | none | same salvage route as NMN, no BubR1 data | status-blocked |
| NMN | lawful dietary supplement since FDA's 2025-09-29 reversal; not an approved medication | none | only published in-vivo BubR1 abundance increase (wild-type mice) | status-blocked, best biological |

The dosage-raise direction now has an approved-medication member, niacin, so the honest v2 statement changes from 'no approved drug exists on this axis' to 'an approved NAD+ precursor exists and still fails the promotion rule'. Nothing on this axis is promoted: the only BubR1-abundance evidence in the class belongs to a different molecule (NMN) in wild-type mice from a single paper, niacin has no BubR1 or aneuploidy data plus an explicit 'pediatric safety not established' label and a documented adult harm profile at pharmacologic doses, and nicotinamide points the wrong way on SIRT2. The class earns a named bench test rather than a proposal: measure BUBR1 protein in proband fibroblasts under an NAD+ precursor before anyone discusses dosing a child.

One correction to the v2 framing follows from this. The dosage-raise axis is not empty of approved medicines, and the report now says so. Niacin's rejection rests on evidence and pediatric labeling, not on absence: Niacin is the approved-medication analog of the NMN route and it changes the NMN verdict's shape: the dosage-raise direction is not empty of approved drugs after all. It still fails the promotion rule. The direction evidence is one paper in wild-type mice using a different molecule (NMN), niacin itself has no BubR1 or aneuploidy data, and its pediatric record is an explicit 'not established' with a known adult harm profile at pharmacologic doses. Naming it as a proposal would be mechanism-by-analogy, which is the failure mode the pre-declared rule exists to block.
### Why amlexanox is unchanged at bench grade

Its dual mechanism (NMD inhibition plus readthrough promotion, Salvatori 2012, PMID 22938201; Atanasova 2017, PMID 28549954) remains the only restore-class route that addresses both blockers of this exact allele. The v2 quantitative lanes confirm the verdict rather than promoting anything: amlexanox carries no reversal-direction connectivity in any of its 11 LINCS Phase I signatures (five rows show positive tau up to +35.6 in A549 at 24 h, the rest sit at the lane's null 0), and its exposure-versus-potency gap stays the documented ~10x fold (RDEB dosing serum around 5 microM versus an in-vitro floor of 50 microM). Bench proof-of-mechanism on the child's own cells, falsifiable in weeks.

### Why fostamatinib, the screen's loudest voice, is rejected

Fostamatinib tops the knowledge-graph lane at percentile 99.8 and is the only approved compound with a direct metapath edge to the MVA1 node. It has no published evidence of aneuploidy-selectivity, no BUBR1 or SAC link, its other lanes sit at background, and it is adult-ITP-only with no pediatric approval. Under the pre-declared rule a screen flag alone does not promote; SYK inhibition is a different field. It is recorded, not suppressed: mouse-model or direct-evidence channels can revisit it.

### Rejected classes, recorded

mTOR inhibition (wrong direction and class warning), checkpoint-kinase inhibition (wrong direction, none approved), aspirin and sulindac (precedent does not transfer; pediatric FAP RCT negative, aspirin contraindicated in children), proteostasis-pressure agents (bortezomib and HSP90 inhibitors; the HSP90 direction degrades the unstable client BUBR1 further), and the chemical chaperone class for this locus (fpocket finds no druggable pocket; TUDCA failed its ALS pair in phase 3 and is not FDA-approved; sapropterin's precedent needs a cofactor pocket the pseudokinase lacks). Reasons and citations as in v1 and the v2 candidate grades.

### Surveillance remains the comparator

No MVA-specific surveillance guideline exists, MVA-adapted childhood cancer predisposition frameworks govern (Kratz 2017, PMID 28168833), and any proposal must beat surveillance on a risk-adjusted basis.

### Can any approved small molecule stabilize the pseudokinase domain?

fpocket said the mutation site has no druggable pocket, and a pharmacological chaperone does not have to bind at the mutation: tafamidis stabilizes transthyretin at the dimer interface, far from most destabilizing substitutions. So the v2 work asked the harder question with a bounded computation.

**The site annotation is cross-species.** The two published structures for this domain, 6JKK and 6JKM, are *Drosophila* BubR1 kinase domain (UniProt A1Z6I7), not human. 6JKM carries ADP plus two magnesium ions, and its 14 contact residues map onto human BUB1B 774, 781, 793, 795, 840, 841, 842, 843, 886, 887, 889, 910, 911. The mapping is anchored independently: the fly catalytic lysine lands on human K795, the residue mutated in the standard BUBR1 kinase-dead K795R construct, and the fly HRD aspartate lands on human R886, matching the documented pseudokinase degeneracy. Only 4 of 13 contacts are identical between the species, so away from that anchor the annotation is approximate.

**Does the nucleotide site itself still bind a nucleotide?** Experimentally, the fly ortholog's pseudokinase domain does: 6JKM is an ADP + Mg complex at 1.95 angstrom. Computationally the picture is weak on both sides. Redocking ADP into its own crystal put the top-scored pose 6.03 angstrom from the experimental pose, with the closest of 9 poses at 3.37 angstrom ranked by score below three worse ones; the usual success bar is 2 angstrom. Docked into the human model, ADP scores -7.26 and ATP -6.73 kcal/mol, which is unremarkable. The honest reading is that this docking setup has no pose accuracy for a charged, flexible nucleotide in this site (magnesium was not modeled), so it cannot settle whether human BUBR1 binds ATP; the structural literature on the fly ortholog is the better evidence, and human BUBR1 is catalytically degenerate regardless.

**Pocket prediction agrees with fpocket about the mutation site.** p2rank 2.4.2 with its AlphaFold-specific model ranks the nucleotide site first on the human model (score 8.86, probability 0.41, 5.40 angstrom from the mapped contacts, K795 among its residues). Its nearest pocket to N1002 is 25.05 angstrom away. Two independent pocket finders therefore agree: there is no pocket at the mutated residue. On the fly crystal the same predictor scores its nucleotide pocket 23.25 with probability 0.89, far more pocket-like than the human model's equivalent.

**The screen.** 2240 ligands (approved small molecules 250 to 600 Da plus every named case candidate regardless of size, plus ADP and ATP as references) were docked with smina against three boxes on AF-O60566-F1: the nucleotide site, the two fpocket pockets that touch N1002, and fpocket pocket 21, the most druggable pocket overlapping the domain. Stage 1 screened everything at low exhaustiveness (35.5 CPU-hours in total); stage 2 re-docked the best 40 per box at exhaustiveness 16 and rescored with Vinardo.

| Box | dist. to N1002 | RaSP burden of lining residues | best approved drug (Vina, exh 16) | ADP reference | ATP reference |
| --- | --- | --- | --- | --- | --- |
| nucleotide site | 30.69 A | 1.639 kcal/mol (1.932x protein mean) | DUTASTERIDE, -10.06 kcal/mol | -7.28 | -6.76 |
| mutation site | 1.01 A | 1.87 kcal/mol (2.204x protein mean) | ZAFIRLUKAST, -8.56 kcal/mol | -6.43 | -6.39 |
| pseudokinase best | 41.95 A | 1.194 kcal/mol (1.408x protein mean) | NILOTINIB HYDROCHLORIDE MONOHYDRATE, -8.47 kcal/mol | -6.70 | -5.94 |

**How the molecules this case already cares about score.** None of the graded candidates is anywhere near the top: hydroxychloroquine sits at the 16th, 31st, and 75th percentile across the three boxes, tafamidis itself at the 67th to 79th, and the two NAD+ precursors near the bottom (acipimox 2nd to 8th, niacin 11th to 14th). Two checks cut against trusting the ranking's top. Migalastat, a clinically validated pharmacological chaperone, scores in the bottom 6 percent everywhere, while lumacaftor, the CFTR chaperone, is the only named chaperone-class molecule near the top in all three boxes (99.3rd to 99.9th percentile), so known chaperones land at both extremes of this ranking. And while ADP, the site's natural ligand, does land in the top quintile of the library (83rd percentile in its own pocket), 372 ordinary approved drugs still outscore it, so the top of this ranking is exactly the territory hundreds of decoys reach. A scoring function that places known chaperones at both extremes and crowds its natural ligand under hundreds of drug-sized decoys cannot nominate one.

**The result is an honest negative.** No approved molecule separates itself from the pack in any box. The best scores sit in the range ordinary drug-sized molecules reach against any shallow protein surface, the natural nucleotide references land in the same range rather than below it, and the scoring function that produced these numbers already failed its own pose-recovery control in the crystal. Nothing here supports naming a stabilizer candidate. What the computation does support is a boundary statement: across the approved small-molecule space, the BUBR1 pseudokinase domain offers no pocket at the mutation, no pocket with convincing predicted ligandability (fpocket druggability at most 0.169 anywhere, 0.111 inside the domain; p2rank probability 0.41 at best), and no molecule whose predicted binding stands out from background. A tafamidis-style route for this allele would need a new binding site discovered experimentally, not a repurposing hit.
## Verdict updates from v1 to v2

| Candidate | v1 verdict | v2 verdict | Reason |
| --- | --- | --- | --- |
| metformin | primary | demoted, protection fallback | aneuploidy-selectivity evidence incidental, weaker than the chloroquine signal in the same 2011 study, unreplicated |
| hydroxychloroquine class | not graded | promoted, primary secondary-prevention hypothesis | strongest published aneuploidy-selective route for this direction plus established pediatric malaria label (>=31 kg tablet caveat; chloroquine liquid below) |
| NMN | not graded | best biological, status-blocked | only in-vivo BubR1 raise (wild-type mice, North 2014); lawful supplement since the 2025-09-29 FDA reversal, not an FDA-approved medication |
| niacin (nicotinic acid) | not graded | rejected for this case, class analog recorded | only FDA-approved NAD+ precursor; no BubR1 or aneuploidy data; pediatric safety and effectiveness not established on the label; chronic pharmacologic dosing carries a documented adult harm profile |
| nicotinamide | not graded | rejected, direction-conflict | SIRT2 inhibition antagonizes the NAD+/SIRT2 axis |
| acipimox | not graded | rejected | no FDA approval, weakest mechanism tie (HCAR2 lipolysis agent, not an efficient NAD+ precursor) |
| nicotinamide riboside | not graded | status-blocked | supplement, no BubR1 data |
| amlexanox | fallback-research | unchanged | dual NMD+readthrough mechanism; exposure gap ~10x; no reversal-direction connectivity on the LINCS lane |
| ataluren, aminoglycosides, rapalogs, aspirin, proteostasis inhibitors, checkpoint kinase inhibitors | rejected | unchanged | v1 reasons stand; the HSP90 direction-conflict is now formalized at the allele level |
| surveillance | comparator | unchanged | protective frame for unproven chemoprevention |

## Falsifiable statements and validation plan, revised

Restore chain (unchanged tests).
- R1. BUBR1 protein in proband fibroblasts sits far below family-control level on Western blot, matching transcript loss from the PTC allele plus instability of the missense allele. This falsifies or confirms the shared-dosage model.
- R2. Amlexanox at tolerated concentrations raises PTC-transcript abundance by allele-specific qPCR before any full-length protein appears. Predicted recovery is below 5% for TGAA and fails the half-normal threshold, and the experiment decides.
- R3. Complementing proband cells with wild-type versus kinase-domain-deletion BUBR1 separates scaffold from catalytic loss, characterizing p.Asn1002Lys directly. An independent read: a thermal-stability shift of the size RaSP predicts for I909T-class mutants should appear in proband-cell BUBR1 half-life measurements; its absence falsifies the instability reading.

Protect chain (revised).
- P1. Serially passaged proband fibroblasts accumulate lagging chromosomes and micronuclei faster than matched family controls, confirming a measurable intrinsic aneuploid-prone state in non-tumor tissue.
- P2. Chloroquine or hydroxychloroquine at malaria-attainable exposures reduces outgrowth of the aneuploid subfraction relative to euploid sibling cultures, quantified with karyotype-resolved colony counts. A null finding falsifies the class hypothesis before any clinical use; metformin runs in the same assay as the fallback arm at pediatric-reachable exposures.
- P3. If NMN ever clears a medication-grade path, a blinded fibroblast time-course measuring BubR1 protein half-life versus vehicle tests the dosage-raise direction directly; the SIRT2-overexpression background of the hypomorphic data makes a SIRT2 inhibitor the cleanest mechanism control.

Screen falsification.
- S1. Any approved drug tops all three lanes of the screen at once only if the MVA module were actually drug-accessible; the observed empty intersection is the falsification statement realized. A future dataset refresh (later ChEMBL, denser Bioteque bundle, LINCS Phase II) that puts an approved compound at the top of all three lanes re-opens the candidate list under the same rule.

## Innovation and scalability

v2 replaces hand-curated class scoring with a computable, rerunnable screen: fixed data snapshots (ChEMBL 37, BioGRID 4.4.236, Bioteque disgenet bundle, GEO GSE92742), declared module, declared seeds, degree-matched nulls, and a pre-registered decision rule whose screen-promotion ban is enforced by construction. The allele tickets are computation, not assertion: RaSP saturation over the full 1,050-residue protein and fpocket over the v6 AlphaFold model. The whole pipeline is CPU-only; the prod-heavy jobs (23 GB LINCS board scan, 21,000-substitution RaSP scan) are bounded jobs with recorded timing. The same scaffold carries over to any compound-heterozygous loss-of-function case, including the NMD-margin and pocket closed/open calls before bench spend.

## Completeness disclosures

From v1, unchanged: Track 1's mechanical ranking crowned a homozygous 45-base PEX5 deletion above the BUB1B pair, reported in Track 1; that finding does not enter this mechanism chain. Parental segregation has not been performed, so the compound-heterozygous trans configuration is inferred, pending parental testing. New in v2: the screen knows nothing about un-protein-targeted drugs (sugar, vitamin, steroid replacements, antibiotics without annotated protein targets), and coverage is partial on two lanes (584 of 1,811 on the knowledge graph, 587 on the LINCS join); absence from a lane is absence of annotation, not evidence of safety or absence.

## Limitations

- The LINCS reversal lane fails its own held-out hairpin check (0.101 versus null p95 0.184); its values are weak evidence everywhere and promote nothing.
- The Bioteque MVA1 disease node is DisGeNET-inferred and leans toward cancer-gene links for checkpoint genes; the knowledge-graph lane is the thinnest of the three.
- RaSP has a wide identity band (sd 0.58) and a roughly 2 kcal/mol practical resolution; class statements rest on within-protein comparisons to validated alleles.
- Tang 2011 aneuploidy-selectivity evidence is cell culture. Hydroxychloroquine's pediatric record is its malaria label, not long-term prophylaxis in a chromosomal-instability child, and the film-coated label stops at 31 kg.
- NMN evidence is one lab's in-vivo result from wild-type mice, with hypomorphic support from SIRT2 overexpression.
- The screen's drug universe is target-annotated pharmacology, so nutrient-class agents with no ChEMBL mechanism target (the NAD+ precursors, among others) were never scored by any lane. They are graded separately in this report, and their absence from the lane tables is a coverage gap rather than a negative result.
- The stabilizer docking screen is a weak instrument: its engine failed to reproduce the experimental ADP pose in the crystal that contains it (top pose 6.03 angstrom), the ATP-site annotation is mapped from a Drosophila structure with 4 of 13 contacts identical, magnesium was not modeled, and no rescoring or free-energy refinement was run. It supports the boundary statement that no approved molecule stands out, not any positive claim.
- Mouse-derived dosage thresholds, reporter-graded readthrough contexts, and an inferred trans configuration carry over from v1.

## Track 2 methods description form, answered (v2)

**Team name.** Silico EVEE

**Approach in detail.** The variant-to-mechanism layer is the v1 certified chain (Ensembl REST, UniProt, AlphaFold DB, ClinVar) verified once more in v2. On top of it runs a quantitative screen over the full ChEMBL 37 approved-drug set (1,811 drugs, 650 human targets). Three orthogonal lanes score each candidate: network proximity z-scores against the declared MVA module over BioGRID 4.4.236 with 1,000 degree-matched nulls; knowledge-graph cosine scoring to the MVA1 node DOID:0080141 in the Bioteque disgenet bundle, calibrated at AUC 0.917; and a BUB1B-low signature-reversal screen over GEO GSE92742 Phase I (205,034 compound signatures), with the instrument's held-out-hairpin incoherence declared and the lane graded down accordingly. Structural allele tickets are RaSP stability prediction (N1002K +2.41 versus calibrated I909T +3.28 and L1012P +8.70 kcal/mol, identity band +0.39 sd 0.58) and fpocket ligandability (no druggable pocket at the locus). A pre-declared rule bans screen-score-only promotion; verdicts grade each candidate on published direct evidence, regulator-verified status, and pediatric feasibility.

**Automated or manual.** Hybrid. Every lane is scripted and rerunnable with fixed seeds and pinned sources. Verdict assignment, wrapper narrative, and the falsifiable statements are curated against the retrieved regulatory text and PMIDs, and each claim carries its citation.

**Manual curation detail.** Abstract- and label-level matching as in v1, plus the new pediatric-use checks pulled from openFDA label fields and the FDA dietary-supplement docket for NMN status (2025-09-29 reversal of the 2022 exclusion).

**Public data only.** Yes, public sources only, apart from the hackathon-supplied Track 1 result on the proband.

**Public sources in detail.** ChEMBL 37 REST (mechanism, molecule, target). openFDA drug label records. BioGRID 4.4.236 tab3 human physical. Bioteque download API disgenet bundle. GEO GSE92742 Level5 and metadata. Ensembl REST, UniProt O60566, AlphaFold AF-O60566-F1-model_v6. ClinVar VCV000533901.9, rs759242053. Europe PMC/PubMed for all mechanistic citations. EMA Translarna EPAR history. FDA docket FDA-2023-P-0872 response of 2025-09-29 for NMN status. RaSP model assets and public precompute chains.

**Proprietary sources.** None.

**Variant mechanism characterization.** Both alleles are loss of function by different routes. The stop-gain transcript is NMD-destroyed with a 698-nucleotide margin (748 nucleotides to the last junction, six downstream junctions, exon 17/23), recomputed live from Ensembl; readthrough has almost no substrate, and the UGA-TGAA context is permissive, so NMD, not the context, blocks the restore class. The missense allele escapes NMD and reads +2.41 kcal/mol destabilizing on RaSP, the validated I909T moderate-instability class, in a pseudokinase domain with no druggable pocket. The disrupted pathway is the spindle assembly checkpoint; downstream, constitutional mosaic variegated aneuploidy with childhood tumor predisposition.

**Time and effort.** v1 plus v2 within the hackathon window: v2 was roughly one analyst-day of pipeline work on CPU-only infrastructure (three production cluster jobs: the full LINCS board scan in 10 minutes, the saturation RaSP scan, and fpocket), plus curation and review.

**Method abstract (v2, under 500 words).** Starting from the Track 1 compound-heterozygous BUB1B pair, v2 re-verifies the v1 mechanism chain and then replaces qualitative class grading with a quantitative screen over every approved drug with a human protein target in ChEMBL 37 (1,811 drugs). Three orthogonal lanes were ran: network-proximity z-scoring against a declared MVA module over the BioGRID human PPI with degree-matched nulls (calibrations pass: ivacaftor-CFTR z=-3.41, imatinib-ABL1 z=-10.66); knowledge-graph cosine scoring to the MVA1 disease node in the Bioteque disgenet bundle (self-edge calibration AUC 0.917); and a signature-reversal sweep over 205,034 LINCS Phase I compound signatures for a BUB1B-low proxy. The LINCS lane declares its instrument limit in the open: hairpin discordance (held-out reconnection 0.101 versus null p95 0.184) makes it a weak-evidence lane that promotes nothing. Structural tickets close the product-locus picture: the nonsense transcript is NMD-destroyed (748 nucleotide margin, recomputed from Ensembl), and the missense allele reads as a validated mild instability with no druggable pocket (RaSP +2.41 kcal/mol; fpocket 97 pockets, none at the site). The pre-registered rule says screen scores alone never promote; none did, since no ingredient tops all three lanes, and the knowledge graph's loudest approved beat (fostamatinib, percentile 99.8) carries no published rescue-direction evidence or pediatric approval and is explicitly rejected. On direct published evidence, the chloroquine/hydroxychloroquine class becomes the primary secondary-prevention hypothesis (Tang 2011 aneuploidy-selectivity; hydroxychloroquine's established pediatric malaria label, narrowed by the >=31 kg tablet restriction), metformin is demoted to fallback for weak aneuploidy-selectivity support, NMN is the best biological but status-blocked dosage-raise route (government status current as of the 2025-09-29 FDA reversal), and amlexanox stays bench-only. Strengths are quantifiability, calibration-anchored lanes, a pre-registered decision rule, and recorded rejections; limitations are the declared LINCS weakness, DisGeNET-inferred node bias, and one-paper direction evidence for both promoted classes. Bench falsification tests are named for the restore and protect chains.

## References

The v1 reference list stands, with the v2 additions.

1. Hanks S, et al. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B. Nat Genet 2004. PMID 15475955.
2. García-Castillo H, et al. Clinical and genetic heterogeneity in mosaic variegated aneuploidy. Am J Med Genet A 2008. PMID 18548531.
3. Plaja A, et al. Child with mosaic variegated aneuploidy and embryonal rhabdomyosarcoma. PMID 9916837.
4. Suijkerbuijk SJ, et al. Molecular causes for BUBR1 dysfunction in MVA. Cancer Res 2010. PMID 20516114.
5. Suijkerbuijk SJ, et al. The vertebrate mitotic checkpoint protein BUBR1 is an unusual pseudokinase. Dev Cell 2012. PMID 22698286.
6. Baker DJ, et al. BubR1 insufficiency causes early onset of aging-associated phenotypes. Nat Genet 2004. PMID 15208629.
7. Dai W, et al. Slippage of mitotic arrest and enhanced tumor development in BubR1 haploinsufficiency. Cancer Res 2004. PMID 14744753.
8. Baker DJ, et al. Increased expression of BubR1 protects against aneuploidy and cancer. Nat Cell Biol 2013. PMID 23242215.
9. Le Hir H, et al. The exon-exon junction complex in mRNA export and NMD. EMBO J 2001. PMID 11532962.
10. Lindeboom RG, et al. The rules and impact of NMD in human cancers. Nat Genet 2016. PMID 27618451.
11. Linde L, et al. NMD governs response of CF patients to gentamicin. J Clin Invest 2007. PMID 17290305.
12. Dabrowski M, et al. Translational readthrough potential of natural termination codons. RNA Biol 2015. PMID 26176195.
13. Floquet C, et al. Statistical analysis of readthrough levels in mammalian cells. PLoS Genet 2012. PMID 22479203.
14. Welch EM, et al. PTC124 targets genetic disorders caused by nonsense mutations. Nature 2007. PMID 17450125.
15. Kerem E, et al. Ataluren for nonsense-mutation cystic fibrosis. Lancet Respir Med 2014. PMID 24836205.
16. McDonald CM, et al. Ataluren in ACT DMD phase 3. Lancet 2017. PMID 28728956.
17. Salvatori F, et al. Rescue of nonsense mutations by amlexanox in human cells. Orphanet J Rare Dis 2012. PMID 22938201.
18. Atanasova VS, et al. Amlexanox enhances COL7A1 PTC read-through. J Invest Dermatol 2017. PMID 28549954.
19. Dabrowski M, et al. Advances in drug-stimulated translational readthrough. Trends Mol Med 2018. PMID 30134808.
20. Tang YC, et al. Identification of aneuploidy-selective antiproliferation compounds, including chloroquine. Cell 2011. PMID 21315436.
21. Santaguida S, et al. Chromosome mis-segregation generates arrested cells eliminated by surveillance. Dev Cell 2017. PMID 28633018.
22. North BJ, et al. SIRT2 induces the checkpoint kinase BubR1 to increase lifespan. EMBO J 2014;33:1438-1453. PMID 24825348.
23. Blaabjerg LM, et al. Rapid protein stability prediction using deep learning. eLife 2023;12:e82593. PMID 37314150.
24. Le Guilloux V, et al. Fpocket: open source ligand binding pocket detection. BMC Bioinformatics 2009. PMID 19459140.
25. Le Geay F, et al. Bioteque: knowledge-graph embeddings for drug repurposing analogues. Bioteque database and bundles. PMID 36419033.
26. Subramanian A, et al. A next-generation connectivity map: L1000 platform. Cell 2017. PMID 29195078.
27. Guney E, et al. Network-based in silico drug efficacy screening. Nat Commun 2016. PMID 26831545.
28. openFDA drug label records for hydroxychloroquine sulfate and chloroquine, pediatric use sections, checked 2026-08-27.
29. FDA response to citizen petition FDA-2023-P-0872 of 2025-09-29: NMN not excluded from the dietary supplement definition.
30. ClinVar VCV000533901.9, NM_001211.6 c.2210T>G.
31. EMA Translarna EPAR page, non-renewal 2025-03-28.
32. Ensembl REST endpoints for ENST00000287598.
33. UniProt O60566, AlphaFold DB AF-O60566-F1.
