# New Directions Exploration: From a Confounded H&E HER2 Signal to a Cheap HER2 "Targetability" Marker

Status: strategic exploration / brainstorm, recorded 2026-06-06. This is hypothesis-generating
direction-setting, NOT a validated result. It captures a working session on where the project
should go after the TCGA H&E HER2-low-vs-zero direction was shown to be confounded. No experiments
in this document have been run yet. Decisions here are still open — the author will choose a
direction later.

For the existing project state see `00_start_here.md`, `PROJECT_SUMMARY.md` (top level), and
`advisor_brief.md`.

## 1. Why the current direction is stuck (diagnosis)

The flagship question so far — "can image features from H&E predict clinical HER2 status?" — is
stuck, and the reason is upstream of execution. By the project's own controls the TCGA-BRCA
HER2-low-vs-zero signal is most parsimoniously TCGA acquisition/batch plus tumor grade and immune
context, not tumor-cell HER2: two generic foundation models (H-Optimus-0, Virchow2) reproduce the
separation and collapse the same way under source-site holdout, external BCNB validation is modest
(~0.60 AUC), and the within-slide RNA audit shows the GigaTIME virtual channels are only weakly
marker-specific.

The execution has been excellent; the problem is **problem choice**. The original question is (a)
crowded — Valieris et al. 2024 already published H&E-based HER2-low prediction; (b) confound-prone
— demonstrated repeatedly here; and (c) **not aligned with the lab's actual edge**. Generic imaging
foundation models are a commodity; the lab's differentiated asset is HER2 isoform / alternative-splicing
biology. The fix is to repoint the same computational-pathology skill set at a question only this lab
can ask.

## 2. Grounding in the lab's own work (Galante Lab, Hospital Sírio-Libanês)

PI: Pedro A. F. Galante, Centro de Oncologia Molecular, Instituto Sírio-Libanês de Ensino e Pesquisa,
São Paulo. Approx. metrics (Google Scholar, retrieved 2026-06-06): ~4,967 citations, h-index 38,
i10-index 90, 100+ papers.

The lab's five signature research lines:

1. **Alternative splicing & isoform biology** (the methodological core): intron-retention transcriptome
   (RNA 2004), SPLOOCE splicing-variant portal (RNA Biol 2012), Reboot isoform-prognosis tool
   (NAR Cancer 2021), STK3/MST2 splicing isoform (Oncogene 2024), culminating in **Guardia et al.,
   "Alternative splicing generates HER2 isoform diversity underlying antibody-drug conjugate resistance
   in breast cancer," Genome Research 2025** (PMID 40664477, DOI 10.1101/gr.280304.124).
2. **Retrocopies / mobile elements** (the lab's most distinctive niche): RCPedia (Bioinformatics
   2013/2024), retrocopy copy-number polymorphism (PLoS Genet 2013), sideRETRO pipeline (Bioinformatics
   2021), endogenous-retrovirus-in-cancer reviews.
3. **Cancer cell surfaceome** (directly relevant to HER2 / ADC targeting): human cell surfaceome
   (PNAS 2009), surfaceome of breast tumors (2013), SurfaceomeDB (2012).
4. **RNA-binding proteins in brain tumors**: SNRPB (Genome Biology 2016), SERBP1 (Genome Biology 2020),
   Musashi1 networks.
5. **Translational Brazilian cancer genomics**: TP53 R337H founder variant (Lancet Reg Health Am 2025),
   1,171 admixed genomes (Nat Commun 2022), TMB / PD-1-blockade biomarkers (Oncotarget 2015), rectal-cancer
   intratumoral heterogeneity (Annals of Surgery 2017).

The Guardia 2025 paper is the anchor. It expanded the catalog of ERBB2 protein-coding isoforms from
13 to 90, profiled them across 561 primary breast cancers plus mass spectrometry, and showed that
trastuzumab/ADC-resistant cells shift expression toward isoforms **lacking the antibody-binding domain**.
Resistance is not "less HER2" — it is HER2 the drug cannot grab. Threads 1, 3, and 5 all converge on the
direction below; the lab has never combined them with computational pathology, which is the white space.

## 3. Three candidate directions considered

| # | Direction | Uses lab's edge | Confound risk | Feasibility | Impact |
|---|-----------|---|---|---|---|
| A | Isoform-resolved "targetable HER2" score → predicts ADC response/resistance | core | low | high (needs transcript-level data) | high |
| B | Can morphology / spatial see the isoform-resistance state? (image ↔ isoform bridge) | core | medium | medium | medium-high |
| C | A quantitative molecular HER2 measure to replace the irreproducible IHC 0/low/ultralow scale | strong | n/a (framing) | high | high |

Working lean (author, 2026-06-06): pursue the **isoform-targetability biology as the engine and ground
truth** (A), keep the **imaging/spatial bridge** (B) as the distinctive, uncrowded contribution, and use
the **broken-IHC-scale narrative** (C) as the framing/impact story. Direction is still open.

Clinical motivation: T-DXd (trastuzumab deruxtecan) now benefits HER2-low (DESTINY-Breast04) and even
HER2-ultralow (DESTINY-Breast06, 2024), so the IHC 0/low/ultralow ordinal is now therapeutically
load-bearing — while being one of the least reproducible calls in pathology (the same 0-vs-1+ noise that
confounded the TCGA result).

## 4. Competitive landscape scan (2024–2026)

A web scan on 2026-06-06 shows the *pure-biology* lane of thread 1 is becoming an active subfield, which
matters for problem choice:

- A medRxiv patient preprint (2024-11-25, doi 10.1101/2024.11.25.24317569) investigates alternative-splicing
  profiles as a mechanism of resistance to anti-HER2 therapies — overlaps the patient clinical-validation
  move. (Verify whether this is a competitor or an adjacent/own group.)
- ERBB2 i14e / ESRP1-2 splicing work (reviewed in IJMS 2025, doi 10.3390/ijms262411918) already paired a
  mechanism with a therapeutic ASO that re-sensitizes tumors to trastuzumab/ADCs.
- A gallbladder-cancer full-length transcriptome (PMC11825701) shows ERBB2-splicing→trastuzumab-resistance
  generalizing to other cancers via long-read atlases.
- npj Breast Cancer 2025 (s41523-025-00868-y) profiles T-DXd resistance multi-omically; PMC12631751 covers
  T-DXd resistance via loss of HER2 expression/binding. A 2025 review on splicing-mediated resistance to
  antibody therapies already exists.

Implication: the splicing→ADC-resistance *concept* is no longer white space, and the lab's head start
(Guardia 2025) is months, not years. **None of these competitors bring imaging, spatial, or computational
pathology, or a cheap deployable assay.** That is where the defensible contribution remains.

## 5. The concrete idea: a cheap HER2 "targetability" marker

Goal (parameter deliberately fixed): a **cheap, deployable** marker that classifies HER2 along an axis
**orthogonal to IHC (amount) and grade (morphology)** — namely *targetability*: how much of a tumor's HER2
still carries the drug-binding epitope.

### 5.1 The mismatch the idea exploits (verified 2026-06-06)

- Standard clinical HER2 IHC antibodies (Ventana 4B5; Dako HercepTest) bind the **intracellular domain**
  (4B5 reported epitope ≈ aa 1231–1250). They score HER2 regardless of whether the extracellular
  drug-binding region is intact.
- Trastuzumab / T-DXd bind **extracellular subdomain IV** (reported ≈ aa 489–630).
- The clinic already has one proof case: **p95HER2**, a truncated form lacking the ECD — it is
  IHC-positive but trastuzumab cannot bind it. Current IHC cannot distinguish intact from truncated HER2.

So the "distinct axis from IHC/grade" is not hypothetical — it is a known blind spot of the current test.
Guardia's 90-isoform landscape expands that blind spot from one proteolytic variant (p95) to a whole
splice-isoform space.

### 5.2 The unifying readout

> targetability index = (HER2 carrying the drug-binding epitope) ÷ (total HER2 made)

IHC measures only the denominator's amount; grade measures morphology; neither sees this ratio.

### 5.3 Three cheap ways to measure it (ranked)

| | Marker | How | Cost | Deployable | Role |
|---|--------|-----|---|---|---|
| A | **Dual-epitope IHC ratio** | ECD-IV stain (trastuzumab epitope / biotinylated trastuzumab) ÷ standard ICD (4B5); positive ICD with low ECD-IV = "HER2 present, untargetable" | ~2 routine stains | any path lab | **Flagship** |
| B | **Junction RT-qPCR panel** | small primer set spanning Guardia's resistance-associated splice junctions; ratio of antibody-blind to canonical junctions | cheap, FFPE-compatible | molecular labs | Molecular confirmatory |
| C | **Reuse-existing-RNA proxy** | exon/junction ratio computable from RNA-seq tumors already have | "free" where data exists | analysis, not assay | Retrospective discovery |

Recommended flagship: **A** — cheapest, most clinically legible, fits existing workflow (one extra
epitope-specific stain + a ratio). B is the molecular confirmatory; C mines it retrospectively for free.

### 5.4 Prior art and what is genuinely new

Not virgin territory: p95HER2 assays (VeraTag/HERmark era) and a trastuzumab-epitope-specific stain in
gastric cancer (PMC6995606) already exist. What is new and defensible here:

1. **Splice-isoform basis** (Guardia's catalog), not just proteolytic p95 — a much larger, unmapped fraction.
2. **T-DXd specifically.** The bystander payload was assumed to make T-DXd forgiving of low/heterogeneous
   HER2, but the antibody must still bind ECD-IV to internalize; whether epitope integrity predicts T-DXd
   response is genuinely open.
3. A **cheap ratio index validated against the isoform ground truth and real ADC outcomes** — which the
   lab has and the p95 papers did not.

### 5.5 First experiment — the orthogonality gate (cheap, near-term)

The cheapest possible first test uses assets already in hand (the trastuzumab/ADC-sensitive vs. resistant
cell-line pair from the Guardia work, plus the isoform pipeline):

1. On the cell-line pair, measure the targetability index three ways (ECD-IV÷ICD staining; junction
   RT-qPCR; isoform-seq ground truth) and confirm resistant lines drop on the index **while staying
   IHC-positive**. That is the whole thesis on one plate.
2. On the 561-sample cohort, plot targetable-fraction against IHC score and against grade.
   - **Go condition:** wide *spread within categories* (targetability varies among IHC-2+ tumors and among
     grade-2 tumors) → a real orthogonal axis and a target for the cheap marker and the imaging bridge.
   - **No-go condition:** targetable-fraction is flat within IHC/grade categories → no new axis; stop
     cheaply.

Bake in the project's hard-won discipline from day one: proper cross-validation, orthogonality to grade,
and single-scanner / leave-site-out wherever imaging enters.

### 5.6 Honest risks

- **Reagent reality:** does a reliable ECD-IV-specific antibody perform on FFPE, or should biotinylated
  trastuzumab be used as the stain? (Feasibility gate for modality A.)
- **Baseline vs. acquired:** if antibody-blind isoforms mainly emerge under treatment pressure (Guardia's
  resistant lines were selected), the marker may be strongest on on-treatment biopsies (monitoring) rather
  than at baseline (selection). Shapes the use-case; does not kill it.
- **Bystander caveat:** an untargetable cell may still die from a neighbor's payload, so the index may
  predict degree/heterogeneity of T-DXd response rather than all-or-none.

## 6. Open decisions and next steps

Author will decide direction later. The concrete forks are:

1. **Use-case:** predictive biomarker for ADC response vs. quantitative replacement for the irreproducible
   0/low/ultralow IHC scale (different validation designs; the targetability axis fits the predictive use
   most naturally).
2. **Lead modality:** dual-epitope IHC ratio (needs an ECD-IV antibody / biotinylated trastuzumab) vs.
   junction RT-qPCR.
3. **Data in hand:** confirm access to the sensitive/resistant cell-line pair, transcript-level RNA
   (long-read / re-quantifiable BAMs, not just gene-level STAR counts), and the ADC-treated cohort with
   outcomes (size, which ADCs, what endpoint).
4. If greenlit, spec the orthogonality-gate experiment in full (exact index definition, controls, go/no-go
   thresholds).

## 7. References

Bibliographic details retrieved 2026-06-06 via PubMed and web search; verify before citing in a manuscript.

- Guardia GDA, et al. Alternative splicing generates HER2 isoform diversity underlying antibody-drug
  conjugate resistance in breast cancer. Genome Research 2025;35(9):1942–1958. PMID 40664477.
  https://doi.org/10.1101/gr.280304.124
- Valieris R, et al. Weakly-supervised deep learning models enable HER2-low prediction from H&E stained
  slides. Breast Cancer Research 2024. https://pmc.ncbi.nlm.nih.gov/articles/PMC11331614/
- Alternative splicing profiles as a mechanism of resistance to anti-HER2 therapies (medRxiv 2024).
  https://doi.org/10.1101/2024.11.25.24317569
- Alternative splicing-mediated resistance to antibody-based therapies (review, IJMS 2025).
  https://doi.org/10.3390/ijms262411918
- Full-length transcriptome atlas of gallbladder cancer: ERBB2 alternative splicing and trastuzumab
  resistance. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11825701/
- Mechanisms of resistance to trastuzumab deruxtecan by multi-omic profiling. npj Breast Cancer 2025.
  https://www.nature.com/articles/s41523-025-00868-y
- Trastuzumab deruxtecan resistance via loss of HER2 expression and binding.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12631751/
- Trastuzumab-specific epitope evaluation as a predictive/prognostic biomarker (gastric cancer).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6995606/
- HercepTest mAb pharmDx (Dako Omnis) vs Ventana PATHWAY anti-HER-2/neu (4B5) epitope-detection comparison.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11537468/
- Paratope plasticity determines anti-HER2 (4B5) antibody specificity.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12599367/
- Galante PAF, Google Scholar profile. https://scholar.google.com/citations?user=-5x8qCwAAAAJ
- Galante Lab publications. https://www.bioinfo.mochsl.org.br/
