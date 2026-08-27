# Track 2 v2 quantitative candidate grading

Screen summary: no ingredient ranks in the top 100 of all three orthogonal lanes (network proximity, Bioteque knowledge graph, LINCS reversal). The plan's refutation watch is therefore not met: no approved medication emerges from the quantitative screen alone.

Rule applied: a candidate is promoted only by the strongest direct published evidence for its rescue direction plus pediatric feasibility, per the field standards cited in the plan.

## nicotinamide mononucleotide (NMN) — not-promotable-status

Direction: raise BUBR1 dosage

Screen lanes: no lane coverage

Reason: Only direct published in-vivo evidence of raising BubR1 abundance (North 2014: NMN raised BubR1 protein in BubR1-hypomorphic mice via NAD+/SIRT2; SIRT2 overexpression phenocopies). This is the strongest direction-fit in the dossier. Status kills it under Track 2 rules: NMN is an investigational supplement, not an approved medication (FDA 2022 determination excludes NMN from the dietary supplement definition because it is being investigated as a drug). Declared as the best biological candidate, held out of the recommendation by status.

## nicotinamide (vitamin B3) — rejected-direction-conflict

Direction: raise BUBR1 dosage

Screen lanes: no lane coverage

Reason: Not a substitute for NMN mechanistically: nicotinamide is a SIRT2 inhibitor at relevant concentrations (sirtuins consume NAD+ and release nicotinamide, which feedback-inhibits them), so it antagonizes the very NAD+/SIRT2 axis North 2014 ties the BubR1 increase to. Vitamin status, and no published evidence of raising BUBR1. Screen lanes unremarkable.

## sodium phenylbutyrate (4-PBA) — rejected

Direction: stabilize the missense protein

Screen lanes: no lane coverage

Reason: Longest pediatric chronic-dosing record of the chemical chaperones (urea cycle disorders, infants onward; label warns only <20 kg tablet dosing). But 4-PBA's ER-centric osmolyte mechanism does not match a cytosolic/nuclear kinetochore protein, the HSP90-client evidence for BUBR1 instability argues against ER folding help, and its closest rescue analogies (alpha-1-antitrypsin Z, F508del CFTR) failed in humans. No experimental support for BUBR1. The track2 chaperone sweep verdict stands.

## tauroursodiol (TUDCA) — rejected

Direction: stabilize the missense protein

Screen lanes: no lane coverage

Reason: Not FDA-approved as a drug in the U.S.; thin pediatric record (22 neonates, ineffective for its label-adjacent use). Its ALS combination with 4-PBA failed phase 3 (PHOENIX) and was withdrawn 2024. Mechanism unmatched to a cytosolic kinetochore protein.

## sapropterin (BH4) — rejected-precedent-not-transferable

Direction: stabilize the missense protein

Screen lanes: prox_z=-0.26

Reason: Sapropterin is the deepest pediatric chaperone record (approved >=1 month) and proves that a natural cofactor rescues destabilizing missense variants (PAH). But its rescue requires binding the mutant enzyme's cofactor pocket; fpocket finds no druggable pocket at N1002/L1012 and a kinase scaffold's catalytic-inactive pseudokinase has no analogous cofactor site. RaSP-grade instability for N1002K is moderate (below the mild I909T-class), consistent with the chaperone-amenable class, yet no tested rescue molecule exists.

## amlexanox — fallback-research

Direction: restore the nonsense transcript

Screen lanes: prox_z=-1.23; kg_pct=36.0; l1000_min=0

Reason: Verdict unchanged from Track 2, now quantified: the only approved molecule with dual NMD-inhibition + readthrough action directly answering this allele's biology (NMD-anchored UGA, 748 nt past the last junction). Verifiable exposure gap: RDEB mouse dosing serum ~5 uM vs >=50 uM lowest in-vitro effective (PMID 28549954); no pediatric systemic experience. Screen lanes: prox_z -1.23 (nearest of the named picks), KG pct 36, L1000 connection exactly 0 across 11 signatures. Bench-candidate for proband-cell testing only.

## ataluren — rejected

Direction: restore the nonsense transcript

Screen lanes: no lane coverage

Reason: Verdict unchanged. EU authorization ended 2025-03-28 after efficacy unconfirmed; never FDA-approved; NMD-anchored target transcript is the class blocker here.

## aminoglycosides (gentamicin/amikacin) — rejected

Direction: restore the nonsense transcript

Screen lanes: no lane coverage

Reason: Verdict unchanged: toxicity class-hostile to chronic pediatric prophylaxis; NMD limits substrate; partial functional recovery only.

## hydroxychloroquine — promoted-principal

Direction: cut aneuploid-cell fitness

Screen lanes: prox_z=0.20; l1000_min=-17

Reason: Chloroquine-class lysosomal/autophagic stress is the strongest published aneuploidy-selective vulnerability (Tang et al., Cell 2011;144:499, PMID 21315436): aneuploid/trisomic cells are preferentially impaired. Hydroxychloroquine is the child-friendlier of the approved aminoquinolines: FDA label states safety/effectiveness ESTABLISHED in pediatric patients for malaria treatment and prophylaxis (openFDA label snapshot). Screen lanes add nothing but also nothing hostile (prox_z +0.20, one context-dependent L1000 row max tau reversal 17). Under the rule 'strongest direct published evidence for the direction, pediatric-feasible', this class outranks metformin, whose Tang-2011 aneuploidy evidence is weak and explicitly subordinate (Tang 2011: metformin effect markedly weaker than AICAR, not reproduced as a specific finder).

## metformin — demoted

Direction: cut aneuploid-cell fitness

Screen lanes: prox_z=0.38; l1000_min=-21

Reason: Demoted from primary to secondary option in the protection framing. Reason: its aneuploidy-selectivity evidence (Tang 2011) is incidental and markedly weaker than the chloroquine-class signal in the same study; the Track 2 v2 sweep found no dedicated replication of aneuploidy-selective metformin action. Pediatric safety record remains the case's best, so it stays as Fallback within protection, below HCQ class. Screen lanes unchanged: prox_z +0.38, no KG coverage, L1000 reversal max 21 (signal-limited).

## fostamatinib — rejected-screen-flag-only

Direction: cut aneuploid-cell fitness

Screen lanes: prox_z=0.64; kg_pct=99.8; kg_edge; l1000_min=-53

Reason: The strongest screen result: Bioteque KG percentile 99.8 vs MVA1 and the only approved compound with a direct metapath edge to DOID:0080141 (through DisGeNET-inferred gene associations), proximity z +0.64, no L1000 signal. It has NO published evidence in any of the four rescue directions (no aneuploidy-selectivity, no BUBR1/SAC link, SYK inhibition field unrelated), and is adult-ITP-only with no pediatric approval. Under the plan's documented failure-mode rule (screen score alone cannot promote; topiramate-IBD precedent), the correct verdict is rejection, as the screen is noisy-disgenet layer here.

## bortezomib (proteasome inhibition) — rejected

Direction: cut aneuploid-cell fitness (also: stabilize missense protein)

Screen lanes: no lane coverage

Reason: Direction-conflict: proteasome blockade raises mutant BUBR1 levels acutely (Suijkerbuijk 2010) but global degradation blockade is toxic (Track 2 C8 verdict), is used for cancer treatment, and no aneuploidy-selectivity rescue pathway supports it for prophylaxis.

## 17-AAG / HSP90 inhibitors — rejected

Direction: cut aneuploid-cell fitness

Screen lanes: no lane coverage

Reason: HSP90 inhibition degrades the already-unstable mutant BUBR1 (it is an HSP90 client; Suijkerbuijk 2010). Also not approved in the U.S. Direction-conflict formalized now at allele level.

## sirolimus (mTOR inhibition) — rejected

Direction: cut aneuploid-cell fitness

Screen lanes: prox_z=0.47; kg_pct=4.9; l1000_min=-60

Reason: Verdict unchanged (Track 2 C5): mTOR inhibition does not show aneuploidy-selective rescue; rapamycin diet cohorts in BubR1-progeroid mice developed phenotypes at similar rates (per the Track 2 literature sweep).

## aspirin (NSAID chemoprevention) — rejected

Direction: cut aneuploid-cell fitness

Screen lanes: prox_z=0.07; kg_pct=40.7; l1000_min=-29

Reason: Verdict unchanged (Track 2 C7): no MVA-specific or aneuploidy-selective evidence; Reye-class pediatric caution.

## sulindac — noted-negative-precedent

Direction: cut aneuploid-cell fitness

Screen lanes: prox_z=0.07; kg_pct=93.7; l1000_min=-38

Reason: Pediatric FAP RCT (Giardiello NEJM 2002) negative for primary prevention in genotyped children. Catalogs the field's orchestrated cautionary precedent, not a candidate.


## Second-allele structural ticket

- RaSP ddG: N1002K +2.41 kcal/mol; reference unstable allele L1012P +8.70; I909T +3.28.
- fpocket: 97 pockets; max druggability 0.169 overall, 0.111 within the pseudokinase domain; pockets at N1002: [7, 19]; at L1012: none.

## NMD ticket

- PTC sits 748 nt upstream of the terminal exon junction with 6 downstream junctions; NMD predicted, margin 698 nt past the 50-nt rule; stop context UGA-A (TGAA) permissive readthrough context, verified in Track 2 allele map.