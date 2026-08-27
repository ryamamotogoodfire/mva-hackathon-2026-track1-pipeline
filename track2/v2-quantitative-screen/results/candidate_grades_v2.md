# Track 2 v2 candidate grades (quantitative full-FDA screen)

Regenerated 2026-08-27 after final-review fixes. Decision rule: a candidate promotes only with the strongest direct published evidence for its direction plus pediatric feasibility, never via screen score alone.

## Verdict summary

- headline: No approved medication restores BUBR1 from these two alleles (unchanged). Track 2 v2 promotes the chloroquine/hydroxychloroquine class above metformin under the protection framing on the strength of Tang 2011's aneuploidy-selectivity data, and grades NMN as the best biological (status-blocked) candidate for raising BUBR1. No candidate is promoted by screen score alone.

- principal: Hydroxychloroquine class (direct published aneuploidy-selective vulnerability evidence, established pediatric label for malaria)

- principal caveat: Prophylaxis hypothesis, off-label for MVA. HCQ film-coated tablets cannot be crushed/divided and so are labelled only for children >=31 kg; the class route for smaller children is chloroquine liquid 16.67 mg/ml. Not promoted by any screen lane.

- status blocked best bio: NMN (in-vivo BubR1-raising evidence in wild-type mice, North 2014; lawful dietary supplement since FDA's 2025-09-29 reversal of the 2022 exclusion, still not an FDA-approved medication)

- fallback: Metformin (demoted: weak published aneuploidy evidence, excellent pediatric safety; kept as protection fallback); Amlexanox (unchanged fallback-research; order-of-magnitude exposure gap; bench-only)

- rejected or demoted named: nicotinamide mononucleotide (NMN); nicotinamide (vitamin B3); sodium phenylbutyrate (4-PBA); tauroursodiol (TUDCA); sapropterin (BH4); ataluren; aminoglycosides (gentamicin/amikacin); metformin; fostamatinib; bortezomib (proteasome inhibition); 17-AAG / HSP90 inhibitors; sirolimus (mTOR inhibition); aspirin (NSAID chemoprevention)


## Screen lane tops (percentiles/z-scores, none promoted)

- prox_z_top10: [{'stem': 'BRENTUXIMAB VEDOTIN', 'prox_z_min': -9.477298234105637}, {'stem': 'MITAPIVAT', 'prox_z_min': -5.099117570717508}, {'stem': 'NAFTOPIDIL', 'prox_z_min': -3.751947435631389}, {'stem': 'DESERPIDINE', 'prox_z_min': -3.34423785476875}, {'stem': 'TETRABENAZINE', 'prox_z_min': -3.258416650731985}, {'stem': 'VALBENAZINE', 'prox_z_min': -3.237863280110767}, {'stem': 'DENOSUMAB', 'prox_z_min': -3.178207041713928}, {'stem': 'ZILEUTON', 'prox_z_min': -3.1213675662874065}, {'stem': 'DEUTETRABENAZINE', 'prox_z_min': -3.0152792724570463}, {'stem': 'RESERPINE', 'prox_z_min': -2.9656563602445285}]
- kg_pct_top10: [{'stem': 'FOSTAMATINIB', 'kg_pct_max': 99.83848172824551}, {'stem': 'ALPELISIB', 'kg_pct_max': 99.59620432061377}, {'stem': 'GEFITINIB', 'kg_pct_max': 99.47506561679789}, {'stem': 'VANDETANIB', 'kg_pct_max': 97.85988289925298}, {'stem': 'RILMENIDINE', 'kg_pct_max': 97.63779527559055}, {'stem': 'TOVORAFENIB', 'kg_pct_max': 97.59741570765192}, {'stem': 'TENIPOSIDE', 'kg_pct_max': 97.23399959620433}, {'stem': 'ABIRATERONE', 'kg_pct_max': 97.1936200282657}, {'stem': 'TUCATINIB', 'kg_pct_max': 97.13305067635777}, {'stem': 'MOXONIDINE', 'kg_pct_max': 96.87058348475671}]
- l1000_tau_min_top10: [{'stem': 'NIFEDIPINE', 'l1000_tau_min': -69.9043960571289}, {'stem': 'ISOTRETINOIN', 'l1000_tau_min': -67.40926361083984}, {'stem': 'TICLOPIDINE', 'l1000_tau_min': -64.25770568847656}, {'stem': 'CLOPIDOGREL', 'l1000_tau_min': -63.575584411621094}, {'stem': 'PSEUDOEPHEDRINE', 'l1000_tau_min': -63.08785629272461}, {'stem': 'LEVONORGESTREL', 'l1000_tau_min': -63.072998046875}, {'stem': 'FLUPHENAZINE', 'l1000_tau_min': -63.0550537109375}, {'stem': 'TRAMADOL', 'l1000_tau_min': -62.888912200927734}, {'stem': 'VEMURAFENIB', 'l1000_tau_min': -62.82879638671875}, {'stem': 'VALPROIC ACID', 'l1000_tau_min': -62.821292877197266}]

Triple-lane top-100 superposition: []

## Candidates

### nicotinamide mononucleotide (NMN) — not-promotable-status

Only direct published in-vivo evidence of raising BubR1 abundance. North 2014 (EMBO J 33(13):1438-1453; PMID 24825348) raised BubR1 protein in young and aged wild-type mice via NMN (NAD+ raise, reducing SIRT2 deacetylation of BubR1), and the same paper tested SIRT2 overexpression in the BubR1-hypomorphic cohorts, which is where the disease-model evidence lives. Status: FDA's 2022 drug-preclusion determination had excluded NMN from the dietary-supplement definition, but FDA reversed that exclusion on 2025-09-29 in its response to citizen petition FDA-2023-P-0872, so NMN is a lawful dietary supplement. A lawful supplement is not an FDA-approved medication, so NMN remains status-blocked as a candidate medication. Screen lanes unremarkable: Bioteque percentile 42, proximity z -0.9, no L1000 beat.

Regulatory status note: Lawful dietary supplement since FDA's 2025-09-29 response to citizen petition FDA-2023-P-0872 (reversed the 2022 drug-preclusion exclusion); not an FDA-approved medication

Screen lanes: qa=0

### nicotinamide (vitamin B3) — rejected-direction-conflict

Not a substitute for NMN mechanistically: nicotinamide is a SIRT2 inhibitor at relevant concentrations (sirtuins consume NAD+ and release nicotinamide, which feedback-inhibits them), so it antagonizes the very NAD+/SIRT2 axis North 2014 ties the BubR1 increase to. Vitamin status, and no published evidence of raising BUBR1. Screen lanes unremarkable.

Screen lanes: qa=0

### sodium phenylbutyrate (4-PBA) — rejected

Longest pediatric chronic-dosing record of the chemical chaperones (urea cycle disorders, infants onward; label warns only <20 kg tablet dosing). But 4-PBA's ER-centric osmolyte mechanism does not match a cytosolic/nuclear kinetochore protein, the HSP90-client evidence for BUBR1 instability argues against ER folding help, and its closest rescue analogies (alpha-1-antitrypsin Z, F508del CFTR) failed in humans. No experimental support for BUBR1. The track2 chaperone sweep verdict stands.

Screen lanes: qa=0

### tauroursodiol (TUDCA) — rejected

Not FDA-approved as a drug in the U.S.; thin pediatric record (22 neonates, ineffective for its label-adjacent use). Its ALS combination with 4-PBA failed phase 3 (PHOENIX) and was withdrawn 2024. Mechanism unmatched to a cytosolic kinetochore protein.

Screen lanes: qa=0

### sapropterin (BH4) — rejected-precedent-not-transferable

Sapropterin is the deepest pediatric chaperone record (approved >=1 month) and proves that a natural cofactor rescues destabilizing missense variants (PAH). But its rescue requires binding the mutant enzyme's cofactor pocket; fpocket finds no druggable pocket at N1002/L1012 and a kinase scaffold's catalytic-inactive pseudokinase has no analogous cofactor site. RaSP-grade instability for N1002K is moderate (below the mild I909T-class), consistent with the chaperone-amenable class, yet no tested rescue molecule exists.

Screen lanes: prox_z=-0.255312258620436; openfda_pediatric_section=True

### amlexanox — fallback-research

Verdict unchanged from Track 2, now quantified: the only approved molecule with dual NMD-inhibition + readthrough action directly answering this allele's biology (NMD-anchored UGA, 748 nt past the last junction). Verifiable exposure gap: RDEB mouse dosing serum ~5 uM vs >=50 uM lowest in-vitro effective (PMID 28549954); no pediatric systemic experience. Screen lanes: prox_z -1.23 (nearest of the named picks), Bioteque percentile 36 (no top beat). L1000 lane: no reversal-direction connectivity in any of its 11 Phase I signatures; five rows have positive tau (peak +35.6 at A549 24 h) and the rest are exactly 0, the lane's null default across the whole board, so the lane offers no reversal support for it either.

Screen lanes: prox_z=-1.2330347086898545; kg_pct_mva1=36.03876438522108; kg_direct_mva1_edge=False; l1000_tau_min=0.0; l1000_frac_strong_rev=0.0

### ataluren — rejected

Verdict unchanged. EU authorization ended 2025-03-28 after efficacy unconfirmed; never FDA-approved; NMD-anchored target transcript is the class blocker here.

Screen lanes: qa=0

### aminoglycosides (gentamicin/amikacin) — rejected

Verdict unchanged: toxicity class-hostile to chronic pediatric prophylaxis; NMD limits substrate; partial functional recovery only.

Screen lanes: qa=0

### hydroxychloroquine — promoted-principal

Chloroquine-class lysosomal/autophagic stress is the strongest published aneuploidy-selective vulnerability (Tang et al., Cell 2011;144:499, PMID 21315436): aneuploid/trisomic cells are preferentially impaired. Hydroxychloroquine is the child-friendlier of the approved aminoquinolines: the FDA label states safety and effectiveness ESTABLISHED in pediatric patients for malaria treatment and prophylaxis (openFDA label text for hydroxychloroquine sulfate, checked 2026-08-27). Pediatric feasibility caveat: the same label restricts the film-coated 200 mg tablets to children of at least 31 kg because they cannot be crushed or divided, so on-label malaria dosing does not reach small children, and MVA patients are often growth-restricted. For smaller children the class route is the chloroquine 16.67 mg/ml liquid (openFDA label, malaria), which is also the compound that carried the Tang 2011 aneuploidy-selectivity evidence. Grading as promoted-principal for the class with this scope note, because the decidion rule's pediatric leg is narrowed, not closed.

Pediatric formulation limitation: Film-coated HCQ tablets cannot be crushed/divided and so are labelled only for children >=31 kg (openFDA pediatric_use). The class route for under-31-kg children is chloroquine liquid 16.67 mg/ml (openFDA label).

Screen lanes: prox_z=0.2047237416077082; l1000_tau_min=-17.248668670654297; l1000_frac_strong_rev=0.0; openfda_pediatric_section=True

### metformin — demoted

Demoted from primary to secondary option in the protection framing. Reason: its aneuploidy-selectivity evidence (Tang 2011) is incidental and markedly weaker than the chloroquine-class signal in the same study; the Track 2 v2 sweep found no dedicated replication of aneuploidy-selective metformin action. Pediatric safety record remains the case's best, so it stays as Fallback within protection, below HCQ class. Screen lanes unchanged: prox_z +0.38, no KG coverage, L1000 reversal max 21 (signal-limited).

Screen lanes: prox_z=0.3795004544999581; l1000_tau_min=-21.13161277770996; l1000_frac_strong_rev=0.0; openfda_pediatric_section=True

### fostamatinib — rejected-screen-flag-only

The strongest screen result: Bioteque KG percentile 99.8 vs MVA1 and the only approved compound with a direct metapath edge to DOID:0080141 (through DisGeNET-inferred gene associations), proximity z +0.64. Its LINCS reversal read sits at the lane's background level (tau_min -53.2; strong-reversal fraction 0.065 versus a lane background of ~0.070) on a lane that failed its held-out-hairpin check, so the lane cannot score it either way. It has no published evidence in any of the four rescue directions (no aneuploidy-selectivity, no BUBR1/SAC link, SYK inhibition field unrelated) and is adult-ITP-only with no pediatric approval. Rejected under the pre-declared rule: screen-flag-only does not promote (the field's documented topiramate-IBD failure mode).

Screen lanes: prox_z=0.6418712217698518; kg_pct_mva1=99.83848172824551; kg_direct_mva1_edge=True; l1000_tau_min=-53.204505920410156; l1000_frac_strong_rev=0.06542056074766354; openfda_pediatric_section=True

### bortezomib (proteasome inhibition) — rejected

Direction-conflict: proteasome blockade raises mutant BUBR1 levels acutely (Suijkerbuijk 2010) but global degradation blockade is toxic (Track 2 C8 verdict), is used for cancer treatment, and no aneuploidy-selectivity rescue pathway supports it for prophylaxis.

Screen lanes: qa=0

### 17-AAG / HSP90 inhibitors — rejected

HSP90 inhibition degrades the already-unstable mutant BUBR1 (it is an HSP90 client; Suijkerbuijk 2010). Also not approved in the U.S. Direction-conflict formalized now at allele level.

Screen lanes: qa=0

### sirolimus (mTOR inhibition) — rejected

Verdict unchanged (Track 2 C5): mTOR inhibition does not show aneuploidy-selective rescue; rapamycin diet cohorts in BubR1-progeroid mice developed phenotypes at similar rates (per the Track 2 literature sweep).

Screen lanes: prox_z=0.4651112416973863; kg_pct_mva1=4.92630728851201; kg_direct_mva1_edge=False; l1000_tau_min=-60.225215911865234; l1000_frac_strong_rev=0.04042553191489362; openfda_pediatric_section=True

### aspirin (NSAID chemoprevention) — rejected

Verdict unchanged (Track 2 C7): no MVA-specific or aneuploidy-selective evidence; Reye-class pediatric caution.

Screen lanes: prox_z=0.06616931598844276; kg_pct_mva1=40.66222491419342; kg_direct_mva1_edge=False; l1000_tau_min=-28.748361587524414; l1000_frac_strong_rev=0.0; openfda_pediatric_section=False

### sulindac — noted-negative-precedent

Pediatric FAP RCT (Giardiello NEJM 2002) negative for primary prevention in genotyped children. Catalogs the field's orchestrated cautionary precedent, not a candidate.

Screen lanes: prox_z=0.07420175715129616; kg_pct_mva1=93.66040783363619; kg_direct_mva1_edge=False; l1000_tau_min=-37.98935317993164; l1000_frac_strong_rev=0.07142857142857142; openfda_pediatric_section=True
