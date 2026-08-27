# Approved NAD+ precursors graded against the BubR1 dosage-raise direction

## Mechanism being graded against

North et al., EMBO J 2014 (PMID 24825348): raising NAD+ lowers SIRT2-mediated deacetylation of BubR1, which raises BubR1 protein abundance. NMN raised BubR1 in young and aged WILD-TYPE mice. The BubR1-hypomorphic cohorts in the same paper tested SIRT2 overexpression, not NMN, so the disease-model evidence belongs to SIRT2 rather than to any NAD+ precursor. One laboratory, one paper.

## Why this class was missing from the 1,811-drug screen

The screen's universe is target-annotated pharmacology. Nutrient-class agents that act as metabolic substrates rather than through a protein target carry no mechanism annotation, so the entire NAD+ precursor class was invisible to all three lanes. Their absence from the screen is a coverage gap, not a negative result, and it falls on exactly the class the North 2014 mechanism implicates.

- NIACIN (CHEMBL573, max_phase 4, first approval 1982): 0 mechanism records in ChEMBL 37, so no target to place in the network, the knowledge graph, or the target-based joins.
- NIACINAMIDE / nicotinamide (CHEMBL1140, max_phase 4): 0 mechanism records.
- ACIPIMOX (CHEMBL345714, max_phase 4): 1 mechanism record whose target_chembl_id is null and whose action_type is null with mechanism_of_action 'Unknown'. It survives into the approved-molecule table but is dropped by the join to a human protein target, so it never receives a lane score.
- NICOTINAMIDE RIBOSIDE (CHEMBL438497): max_phase 3, not an approved medication, and 0 mechanism records.
- NICOTINYL ALCOHOL (CHEMBL1235535): max_phase -1.

## Candidate grades

### niacin (nicotinic acid) — rejected-for-this-case, class-analog noted

**Approval.** FDA-approved prescription drug. Immediate-release NIACOR (ANDA040378, tablet) and multiple extended-release niacin tablet applications (for example ANDA090892, ANDA204178, ANDA201273). Indicated for dyslipidemia, not for any checkpoint or aneuploidy indication.

**Pediatric record.** Not established. The immediate-release label states safety and effectiveness in children and adolescents have not been established; the extended-release label states safety and effectiveness of niacin therapy in pediatric patients (16 years or younger) have not been established.

**Mechanism tie.** Genuine but indirect: niacin is a NAD+ precursor through the Preiss-Handler route, so it raises NAD+ and can in principle drive the same NAD+/SIRT2 axis North 2014 links to BubR1 abundance. No published experiment tests niacin against BubR1 abundance, BUB1B dosage, aneuploidy rate, or any MVA endpoint.

**Caveats.**
- The pharmacologic niacin doses used in dyslipidemia cause flushing, hepatotoxicity risk (especially sustained-release forms), insulin resistance, and hyperuricemia, and AIM-HIGH plus HPS2-THRIVE showed no cardiovascular benefit with added harm in adults, so chronic pediatric dosing has no supporting safety base.
- Whether an oral NAD+ precursor raises NAD+ in the tissues that matter for chromosome segregation in a child is untested.
- The rescue direction needs a sustained abundance increase across the roughly 50 percent dosage boundary from the mouse ladder; no precursor has been shown to move BUBR1 that far in any system.

**Verdict reason.** Niacin is the approved-medication analog of the NMN route and it changes the NMN verdict's shape: the dosage-raise direction is not empty of approved drugs after all. It still fails the promotion rule. The direction evidence is one paper in wild-type mice using a different molecule (NMN), niacin itself has no BubR1 or aneuploidy data, and its pediatric record is an explicit 'not established' with a known adult harm profile at pharmacologic doses. Naming it as a proposal would be mechanism-by-analogy, which is the failure mode the pre-declared rule exists to block.

### nicotinamide (niacinamide, vitamin B3 amide) — rejected-direction-conflict

**Approval.** Marketed as a vitamin ingredient and in unapproved-drug or supplement listings; not an FDA-approved drug product for a therapeutic indication (openFDA returns a listing with no indications and no application number). ONTRAC used it as an over-the-counter skin-cancer chemoprevention agent in adults.

**Pediatric record.** No FDA-approved pediatric indication or labeled pediatric dosing as a therapeutic drug.

**Mechanism tie.** Direction-conflicting. Nicotinamide is a NAD+ precursor through the salvage route, but it is also a product-feedback inhibitor of sirtuins including SIRT2 at relevant concentrations, and the North 2014 arc raises BubR1 by reducing SIRT2-mediated deacetylation, so nicotinamide pushes the effector in the wrong direction.

**Caveats.**
- The net effect of a substrate that is also an inhibitor of the pathway's effector is not predictable from first principles and has not been measured for BubR1.

**Verdict reason.** Unchanged from the v2 grading: the same molecule that supplies NAD+ inhibits SIRT2, the effector the BubR1 increase is attributed to, and there is no BubR1 evidence for it.

### acipimox — rejected

**Approval.** Not FDA-approved: openFDA returns no label and no application for acipimox. ChEMBL records max_phase 4 on the strength of non-US approvals (European national approvals, ATC C10AD06) as a nicotinic acid analog for dyslipidemia.

**Pediatric record.** No FDA pediatric labeling, since there is no FDA approval; no pediatric use record in the case-relevant sense.

**Mechanism tie.** Weak. Acipimox is a nicotinic acid analog acting on HCAR2 for lipolysis suppression; it is not an efficient NAD+ precursor, and no BubR1, SIRT2, or aneuploidy data exist.

**Caveats.**
- Its ChEMBL max_phase 4 flag with a null mechanism target is exactly why it entered the approved-molecule table and then vanished before scoring.

**Verdict reason.** No US approval, no pediatric record, and the weakest mechanism tie of the three: an HCAR2 lipolysis agent is not a route to raising BUBR1 abundance.

### nicotinamide riboside (NR) — not-promotable-status

**Approval.** Not an approved medication. Marketed as a dietary supplement with FDA new-dietary-ingredient notifications; ChEMBL max_phase 3 reflects trial activity, not approval.

**Pediatric record.** None.

**Mechanism tie.** Same NAD+ salvage route as NMN, one step removed; no BubR1 or aneuploidy data.

**Caveats.**
- Supplement status carries the same blocking logic as NMN.

**Verdict reason.** Closest supplement analog to NMN and blocked for the same reason: a supplement is not an approved medication, and there is no BubR1 evidence specific to it.

### nicotinamide mononucleotide (NMN) — not-promotable-status (best biological)

**Approval.** Lawful dietary supplement since FDA's 2025-09-29 response to citizen petition FDA-2023-P-0872 reversed the 2022 drug-preclusion exclusion. Not an FDA-approved medication.

**Pediatric record.** None.

**Mechanism tie.** The only molecule in this class with a published in-vivo BubR1 abundance increase (North 2014, wild-type mice).

**Caveats.**
- Wild-type mice, not a BubR1-deficient model; the hypomorphic cohorts in that paper tested SIRT2 overexpression.
- One laboratory, one paper, no replication, no human data on BubR1.

**Verdict reason.** Unchanged: strongest direction fit in the dossier, blocked by medication status.

## Direction verdict

The dosage-raise direction now has an approved-medication member, niacin, so the honest v2 statement changes from 'no approved drug exists on this axis' to 'an approved NAD+ precursor exists and still fails the promotion rule'. Nothing on this axis is promoted: the only BubR1-abundance evidence in the class belongs to a different molecule (NMN) in wild-type mice from a single paper, niacin has no BubR1 or aneuploidy data plus an explicit 'pediatric safety not established' label and a documented adult harm profile at pharmacologic doses, and nicotinamide points the wrong way on SIRT2. The class earns a named bench test rather than a proposal: measure BUBR1 protein in proband fibroblasts under an NAD+ precursor before anyone discusses dosing a child.
