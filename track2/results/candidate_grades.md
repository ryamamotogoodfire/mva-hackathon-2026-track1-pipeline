# Track 2 candidate screen, graded

Q1 restore axis, Q2 protect axis, plus explicit wrong-direction rejections and the surveillance comparator.
Grades, mechanism_tie / preclinical_evidence / regulatory / pediatric_experience / spectrum_fit.

## C1-ataluren, nonsense readthrough, selective small-molecule
- verdict: **rejected**
- reason: NMD-dominant allele plus a context weaker than TGAC plus a failed-efficacy regulatory record removes ataluren as a credible candidate.
- agent: ataluren (Translarna/PTC124)
- mechanism tie: strong. evidence: direct-human-trial-negative. regulatory: approved-lapsed. pediatric: labeled (was, age >= 2y, DMD indication). spectrum: none
- note: Selectively induces ribosomal readthrough of premature, not normal, termination codons (Welch 2007). The target allele carries a UGA PTC, the class ataluren was optimized against.
  - blocker: The PTC transcript is NMD-targeted: 748 nt upstream of the last exon-exon junction with six downstream junctions (allele_map.json), far past the ~50 nt rule (Le Hir 2001; Lindeboom 2016). NMD degrades the transcript before translation, so the readthrough substrate is largely absent (Linde 2007).
  - blocker: The tetranucleotide is UGA-A (TGAA), a permissive readthrough context, though weaker than TGAC, which gives the strongest basal and drug-induced readthrough (Dabrowski 2015; Floquet 2012). This window is not the class blocker; the NMD substrate loss is.
  - blocker: Clinical efficacy record is negative to mixed: phase 3 in nonsense-mutation CF missed its primary endpoint (Kerem 2014); ACT DMD did not meet its primary endpoint at the prespecified significance level (McDonald 2017).
  - blocker: Regulatory status today: the EU conditional marketing authorisation was not renewed and expired 2025-03-28 because effectiveness was not confirmed (EMA EPAR page). Never approved in the US.
- citations: 17450125 -- PTC124 readthrough mechanism and nonsense-targeting selectivity; 28728956 -- ACT DMD phase 3 outcome; 24836205 -- nonsense-CF phase 3 primary endpoint missed; 17290305 -- NMD governs response to readthrough treatment; 26176195 -- stop-codon context ranking, UGA-C best; 22479203 -- readthrough level determinants; https://www.ema.europa.eu/en/medicines/human/EPAR/translarna -- EU authorisation history and 2025-03-28 expiry after non-renewal

## C2-aminoglycosides, nonsense readthrough, aminoglycoside antibiotics
- verdict: **rejected**
- reason: Toxicity, weak context-matched efficacy, and the NMD substrate problem; keep as a lab control only.
- agents: gentamicin, amikacin, G418 (research only)
- mechanism tie: strong. evidence: direct-human. regulatory: approved-current (antibiotic indications only, not readthrough). pediatric: labeled (antibiotic); toxicity profile hostile to chronic use. spectrum: none
- note: Aminoglycosides bind the decoding center and force near-cognate tRNA acceptance at stop codons; efficacy depends on stop context (Howard 2000). UGA-A is a permissive context for this class, weaker than TGAC.
  - blocker: Same NMD problem as ataluren: little transcript survives to be read through (Linde 2007).
  - blocker: Class-limiting pediatric toxicity: ototoxicity and nephrotoxicity. Chronic dosing in a child for a prophylaxis hypothesis is not acceptable without a dramatic efficacy lead.
  - blocker: Clinical rescue levels are small: in CFTR stop-mutation patients, gentamicin gave partial functional correction, not protein-level normalization (Wilschanski 2003).
  - blocker: Eliminating even 5-20% function recovery is far below the level this allele requires to matter (see dosage-threshold evidence: Dai 2004, Baker 2004, Baker 2013).
- citations: 10939566 -- context-specificity of aminoglycoside readthrough; 14534336 -- gentamicin partial CFTR functional correction in patients; 17290305 -- NMD limits readthrough substrate

## C3-amlexanox, dual NMD inhibition plus PTC readthrough, repurposed small molecule
- verdict: **fallback-research**
- reason: Only approved-class candidate whose mechanism answers the allele's actual blocker (NMD). Recommendation: proband fibroblast/patient-cell proof-of-mechanism study, not a clinical proposal. Falsifiable at the bench within weeks.
- agent: amlexanox (Aphthasol; Solfa-A in Japan)
- mechanism tie: strong. evidence: direct-human-cells. regulatory: approved-current (topical, US; oral, Japan). pediatric: none. spectrum: none
- note: Uniquely among approved agents, amlexanox both stabilizes nonsense-containing mRNAs (inhibits NMD) and promotes PTC readthrough, producing full-length functional protein in patient-derived cells (Salvatori 2012; Atanasova 2017). This is the only approved molecule whose dual action directly addresses the blocking biology of this allele (an NMD-targeted UGA PTC).
  - blocker: Evidence is early preclinical: patient-cell studies in three genes (TP53/CFTR/DMD-class reporters in Salvatori 2012) and COL7A1 (Atanasova 2017); no animal efficacy study for a systemic nonsense disease; no human trial for any nonsense indication.
  - blocker: Systemic exposure: marketed as a 5% topical paste in the US and as an oral anti-allergy drug in Japan (Dabrowski 2018); pediatric systemic safety is not established.
  - blocker: Dose required for NMD modulation in vitro sits above typical systemic exposure from topical dosing.
  - blocker: The second allele (p.Asn1002Lys) is a last-exon kinase-domain missense, unreachable by this mechanism.
- citations: 22938201 -- amlexanox stabilizes nonsense mRNAs and yields full-length functional protein in human patient cells; 28549954 -- full-length type VII collagen from PTC alleles ex vivo; 30134808 -- drug-stimulated translational readthrough review and amlexanox approval landscape

## C4-metformin, metabolic stress, AMPK-pathway biguanide
- verdict: **primary**
- reason: Only candidate with an unimpeachable pediatric safety record plus a defensible, falsifiable mechanism chain. Presented as a prophylaxis hypothesis with explicit preclinical tests, not a treatment recommendation.
- agent: metformin (immediate-release)
- mechanism tie: moderate. evidence: indirect. regulatory: approved-current (US >=10y T2D; worldwide). pediatric: labeled (10-16y T2D). spectrum: indirect (general tumor-prevention logic, no MVA-specific data)
- note: Aneuploid cells suffer energy, proteotoxic, and lysosomal stress; AICAR, the pharmacological energy-stress inducer in the same AMPK-family axis, selectively kills aneuploid cells (Tang 2011; Santaguida 2015; reviews 2020s). Metformin is the approved pediatric-safe agent in that mechanistic neighborhood. The tie is indirect: no study has shown aneuploidy-selective fitness reduction with metformin itself.
  - blocker: No direct aneuploidy-selectivity evidence for metformin; the AICAR composite in Tang 2011 combined three agents.
  - blocker: A germline carrier's normal tissues are not yet aneuploid; the drug can only act by reducing the fitness of emerging aneuploid clones, a prevention scenario untested in any organism.
  - blocker: Observational diabetes-cohort data cannot be transferred to germline-LoF prevention.
  - support: Pediatric safety is the strongest of any candidate class: FDA label establishes safety and effectiveness for type 2 diabetes in children 10-16 years (openFDA label record); the TODAY trial studied late teens (Zeitler 2012).
  - support: Chronic-use risk profile is mild and monitorable: vitamin B12 deficiency on long-term use (DPPOS, PMID 26900641); boxed warning for lactic acidosis is a renal-failure-context contraindication, managed by renal monitoring.
  - support: Population-level cancer-prevention signals in diabetic cohorts are consistent but confounded (Decensi 2010 meta-analysis; Pollak 2012 review).
- citations: 21315436 -- energy stress (AICAR) selectively antiproliferates aneuploid cells; 26404941 -- aneuploidy-induced lysosomal/TFEB stress state; 22540912 -- pediatric metformin trial (TODAY); 26900641 -- long-term B12 deficiency risk; 20947488 -- cancer-incidence meta-analysis in metformin users; openFDA: metformin hydrochloride tablets, sections 1 and 8.4 -- pediatric 10-16y labeling and boxed warning

## C5-rapalogs, mTOR inhibition
- verdict: **rejected**
- reason: Mechanism points the wrong way on aneuploid-cell fitness, and the class's own label warns of malignancy from immunosuppression in the highest-surveillance-need population.
- agents: sirolimus, everolimus, temsirolimus
- mechanism tie: wrong-direction for aneuploidy selectivity. evidence: direct-animal (tolerance), none for the aneuploidy-prevention claim. regulatory: approved-current. pediatric: labeled (sirolimus >=13y transplant; everolimus >=1y TSC-SEGA). spectrum: none
- note: mTOR inhibition relieves translation burden and induces autophagy. Aneuploid cells are protein-synthesis-stressed and autophagy-limited (Santaguida 2015), so mTOR inhibition is expected to relieve exactly the bottleneck that disadvantages aneuploid cells, improving their survival fitness rather than reducing it. The pediatric TSC program (everolimus, EXIST-1) shows the class can be given to young children, but that experience concerns tolerance, not prevention of aneuploidy-driven tumors.
  - blocker: Sirolimus carries a boxed warning: immunosuppression increases susceptibility to infection and development of lymphoma and other malignancies (openFDA label) - the specific harm an MVA patient can least tolerate, since residual immune surveillance may clear aneuploid clones (Santaguida 2017, PMID 28633018, shows p53/intrinsic clearance of these cells).
  - blocker: Pediatric sirolimus labeling holds only for renal-transplant use at age >= 13y; everolimus SEGA labeling starts at 1y (EXIST-1), an oncology-tolerance precedent rather than prevention data.
- citations: 26404941 -- aneuploid cells are autophagy/lysosome-limited; 28633018 -- elimination of mis-segregated cells by p53-dependent surveillance; 21047224 -- everolimus pediatric TSC-SEGA efficacy/tolerance; openFDA: sirolimus, boxed warning and section 8.4 -- malignancy-from-immunosuppression warning and >=13y labeling

## C6-checkpoint-kinase-inhibitors, checkpoint-adjacent kinase inhibition (Mps1/TTK, Aurora B)
- verdict: **rejected**
- reason: Documented wrong-direction candidates; listing them explicitly marks the boundary of the repurposing space.
- agents: CFI-402257, BAY 1217389, barasertib, reversine (research only)
- mechanism tie: wrong-direction. evidence: direct-animal, as harm. regulatory: none approved. pediatric: none. spectrum: none for prevention; theoretical for treatment of any future tumor
- note: These drugs remove residual checkpoint signaling. In an MVA background that is already hypomorphic for a checkpoint protein, checkpoint inhibition raises chromosome mis-segregation: checkpoint blockade causes massive chromosome loss (Kops 2004). The only defensible use is against established, highly aneuploid tumors, where checkpoint inhibition is exploitative (Cohen-Sharir 2021) - a treatment setting, not prevention, and a population this patient does not yet belong to. Mps1 inhibitors remain investigational (Mason 2017); none is approved.
- citations: 15159543 -- checkpoint inhibition causes massive chromosome loss; 33505028 -- aneuploid cancer cells vulnerable to checkpoint inhibition (treatment-only setting); 28270606 -- Mps1 inhibitor CFI-402257 remains investigational (trial-stage)

## C7-aspirin, cancer chemoprevention in a hereditary instability syndrome (NSAID)
- verdict: **rejected**
- reason: Only syndrome-precedent value. Pediatric safety and spectrum mismatch rule it out.
- agent: aspirin
- mechanism tie: weak-for-MVA (CIN-spectrum precedent only). evidence: direct-human-trial-positive in a different syndrome, adult-only. regulatory: approved-current (OTC). pediatric: labeled against. spectrum: none for MVA tumors (RMS/leukemia/Wilms, not colorectal)
- note: CAPP2 randomized 600 mg daily aspirin in Lynch syndrome and showed reduced colorectal cancer incidence at >=2 years of use (Burn 2011, Burn 2020). This is the strongest randomized chemoprevention evidence in any hereditary cancer-instability syndrome, but Lynch is mismatch-repair-driven, not chromosome-instable, and the benefit accrued in adults.
  - blocker: Pediatric contraindication class-wide: aspirin is not given to children/teenagers because of Reye's syndrome risk (standard US label language), and the MVA risk window is childhood.
- citations: 22036019 -- CAPP2 first long-term aspirin report; 32534647 -- CAPP2 10-year follow-up, primary endpoint

## C8-proteostasis-pressure, proteotoxic-stress amplification (Hsp90 inhibitors, proteasome inhibitors)
- verdict: **rejected**
- reason: Cannot be given as long-term prophylaxis to a child; retained only as preclinical confirmation that the aneuploidy-proteostasis axis is druggable.
- agents: 17-AAG/tanespimycin (investigational), bortezomib
- mechanism tie: moderate. evidence: direct-human-cells. regulatory: approved-current for bortezomib (oncology); Hsp90 inhibitors none. pediatric: studied (bortezomib in pediatric ALL), prophylaxis-incompatible. spectrum: indirect
- note: Aneuploid human cells are HSF1/Hsp90-dependent (Donnelly 2014), and Hsp90 inhibition plus energy stress was aneuploidy-selective in Tang 2011. However no aneuploidy-selective Hsp90 inhibitor is approved; approved proteasome inhibition (bortezomib) is an oncology agent intolerable as prophylaxis in a child.
- citations: 25205676 -- HSP90/HSF1 dependence of aneuploid human cells; 21315436 -- 17-AAG in the aneuploidy-selective screen

## C9-surveillance, tumor surveillance plus symptom-directed management
- verdict: **comparator**
- mechanism tie: n/a. evidence: n/a. regulatory: n/a. pediatric: n/a. spectrum: direct
- note: No syndrome-specific surveillance guideline exists for MVA1. Practice follows the childhood cancer predisposition frameworks (Kratz 2017) and the MVA tumor spectrum: Wilms tumor, embryonal rhabdomyosarcoma, leukemias/lymphomas (Hanks 2004; García-Castillo 2008; case series). Surveillance is what any candidate must beat on a risk-adjusted basis.
- citations: 28168833 -- childhood cancer predisposition screening recommendations framework; 15475955 -- MVA1 description and tumor spectrum; 18548531 -- MVA clinical heterogeneity delineation; 9916837 -- MVA plus embryonal rhabdomyosarcoma case matching this proband; 28553959 -- TRIP13 MVA2 Wilms predisposition for spectrum context

