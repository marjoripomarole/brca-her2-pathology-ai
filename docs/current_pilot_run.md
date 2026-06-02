# Current Pilot Run

Status: Current project status summary. For the cleanest latest-results presentation, use `docs/clinical_her2_high_trust_tile128_results.md`; for the history, use `docs/paper_proposal_process_log.md`.

## Run Status

The current workspace contains both the original 30-slide clinical HER2 pilot and the expanded 60-slide clinical HER2 run. The original pilot processed:

- HER2-positive cases processed: 10
- HER2-low cases processed: 10
- HER2-zero cases processed: 10
- Slides processed: 30
- Initial GigaTIME tiles per slide: 64 random tissue tiles
- Robustness GigaTIME tiles per slide: up to 256 random tissue tiles
- Initial total tile predictions: about 1,920
- Robustness total tile predictions: about 7,600, depending on available tissue tiles per slide
- Device used in the current run: Apple MPS

The HER2 group labels come from `data/tcga_brca/clinical_her2_labels.csv`, which was built from TCGA-BRCA clinical HER2 IHC/ISH fields. The selected 30-case cohort is recorded in `data/tcga_brca/clinical_her2_cohort_cases.csv` and summarized in `docs/clinical_her2_cohort_selection.md`.

## Expanded 20/20/20 Run

An expanded clinical HER2 run is now complete:

- HER2-positive cases processed: 20
- HER2-low cases processed: 20
- HER2-zero cases processed: 20
- Slides processed: 60
- Tile sampling: up to 256 random tissue tiles per slide
- Total tile predictions: 15,225
- Matched STAR-count RNA-seq expression downloaded for all 60 selected cases

The expanded run strengthened the HER2-low versus HER2-zero finding. Several all-tissue and QC-cellular virtual immune/myeloid/checkpoint channels now pass within-view BH correction for the HER2-low versus HER2-zero pairwise comparison. The HER2-low versus HER2-zero classifier also remained at about 0.80 balanced accuracy in the expanded 40-case binary comparison.

The three-group interpretation became more nuanced: HER2-low is often the lowest virtual immune/checkpoint group, while HER2-positive becomes highest for several broader virtual immune programs. RNA validation remains weak, so this is still hypothesis-generating.

See `docs/clinical_her2_expanded20_results.md`.

## Laptop-Realistic 61/61/61 Cohort Ready For Next Run

A larger balanced cohort has now been prepared and downloaded for scaling beyond the 60-slide analysis:

- HER2-positive cases selected: 61
- HER2-low cases selected: 61
- HER2-zero cases selected: 61
- Slides selected: 183
- Slides downloaded locally: 183 of 183
- Slide storage footprint: about 31 GB under `data/tcga_brca/slides/`
- Remaining local disk after download: about 144 GiB free

This is the largest balanced clinical HER2 cohort available from the current label table because HER2-zero is the limiting group with 61 candidate cases. The selection used one diagnostic primary-tumor H&E whole-slide image per case, prioritized direct clinical labels, and favored already-local or smaller slide files for laptop feasibility.

This larger cohort became the basis for the current strict high-trust 171-slide analysis. The expanded 20/20/20 run is now historical and mainly useful as a 256-tile comparison point.

HER2 label and slide-metadata QC has also been run on this 183-slide cohort. After checking Guardia et al., Genome Research 2025, PMID 40664477, "Alternative splicing generates HER2 isoform diversity underlying antibody-drug conjugate resistance in breast cancer," we kept the relevant female-patient TCGA sample-selection principle and added our own H&E slide-specific checks. The current trustworthy-slide list finds 171 strict high label+slide trust slides, 9 slides needing review before primary analysis, and 3 male-patient slides excluded from primary analysis. See `docs/clinical_her2_trustworthy_slide_list.md`.

## Strict High-Trust 171-Slide Tile128 GigaTIME Analysis

The high-trust cohort was processed with a laptop-safe GigaTIME run, then filtered to the strict 171-slide primary analysis set:

- HER2-positive strict high-trust slides analyzed: 53
- HER2-low high-trust slides processed: 57
- HER2-zero high-trust slides processed: 61
- Total strict primary-analysis slides: 171
- Tile sampling: 128 random tissue tiles per slide
- Total primary-analysis tile predictions: 21,888
- Device: Apple MPS
- Heatmaps: disabled for runtime/storage practicality
- Runner change: `scripts/run_gigatime_tcga_brca.py` now supports `--resume`

This is now the largest and cleanest analyzed result in the project. The main finding is again HER2-low versus HER2-zero separation: HER2-low is lower than HER2-zero for `CD68`, `PD-L1`, `CD11c`, `PD-1`, `CD4`, `CD3`, and `CK`, with BH-corrected pairwise q values around 0.002-0.003 in the all-sampled-tissue view.

The parameter/settings robustness check compared the overlapping slides from the earlier 60-slide 256-tile run with the current strict high-trust 128-tile analysis. There were 56 overlapping slide IDs, including all 20 HER2-low and all 20 HER2-zero expanded-run slides. All 8 tested key channels had the same HER2-low versus HER2-zero direction across runs, and 7 of 8 had HER2-low lower than HER2-zero in both runs.

The classifier result is useful but not diagnostic:

- Best HER2-low versus HER2-zero balanced accuracy: 0.727
- Best HER2-low versus HER2-zero macro AUC: 0.787
- QC-cellular view preserved the signal with balanced accuracy 0.719
- HER2-positive versus HER2-negative classification remained weak, with best balanced accuracy about 0.574
- Three-class HER2 classification remained weak/moderate, around 0.50-0.52 balanced accuracy

The ER/PR and HER2-detail sensitivity analysis strengthens the main interpretation:

- In all sampled tissue, 7 of 8 tested key channels remain significant after ER/PR adjustment.
- The signal remains visible across HER2-low IHC `1+`, HER2-low IHC `2+`/ISH-negative, HER2-zero IHC `0`/ISH-negative, and HER2-zero IHC `0`/ISH-not-evaluated subgroups.
- The strict CK-enriched view still weakens, keeping the interpretation focused on tissue context rather than pure tumor-cell HER2 phenotype.

The case-level driver analysis adds a review layer:

- 118 HER2-low/HER2-zero slides were scored using the significant low-versus-zero virtual channels.
- 71 slides matched the expected HER2-low/HER2-zero direction in at least 3 of 4 cleanup views.
- 47 slides showed the opposite profile in at least 2 cleanup views.
- 37 cases were misclassified by the best low-vs-zero classifier in at least 2 cleanup views.
- These opposite-profile and classifier-error cases are now the highest-priority manual pathology/QC review list.

The case-driver visual QC adds an important caveat:

- 8 representative case-driver cases were rendered as H&E plus virtual mIF panels.
- Low-like selected tiles had high tissue fraction but very low virtual CK, CD68, PD-L1, and CD11c.
- Several low-like tiles looked stromal/collagen-rich rather than clearly tumor-rich.
- This means the HER2-low versus HER2-zero signal may partly reflect tissue composition unless pathologist review or tumor-rich tile restriction supports it.

The tissue-composition sensitivity analysis strengthens that caveat:

- HER2-low has more low-marker tiles than HER2-zero: 0.349 versus 0.180, BH q = 0.000265.
- HER2-low has fewer high-marker and absolute CK-high QC tiles than HER2-zero.
- The case-driver score is strongly correlated with marker/tissue composition.
- After adjusting for low-marker tile fraction, most HER2-low versus HER2-zero channel effects collapse.
- Current best interpretation: GigaTIME detects a HER2-low versus HER2-zero tissue-context association, not yet a tumor-cell HER2 biology or diagnostic signal.

The tumor-rich proxy sensitivity analysis adds the next layer:

- Fixed-count virtual CK-rich tile views weaken individual channel-test significance.
- Top 8 CK tiles per slide has no q<0.05 individual channel differences, but the multichannel low-vs-zero classifier still reaches balanced accuracy 0.727 and macro AUC 0.755.
- Absolute CK-high QC tiles reach low-vs-zero classifier balanced accuracy 0.761 and macro AUC 0.782, but fewer slides qualify and HER2-zero retains more CK-high tiles than HER2-low.
- Best interpretation: there may be a multichannel image pattern in CK-rich proxy regions, but this remains GigaTIME-derived and needs real tumor-rich/pathologist validation.

The classifier permutation sanity check asks whether those low-vs-zero classifiers still look good after shuffling labels:

- Observed repeated-CV balanced accuracy across selected views: 0.693-0.729.
- Shuffled-label null mean balanced accuracy: 0.482-0.488.
- All tested views beat the shuffled-label null with empirical p = 0.0099 and BH q = 0.0099.
- Best interpretation: the classifier signal is not obviously random, but this is post-hoc exploratory evidence and still needs tumor-rich/pathologist validation.

The nested model-selection check is stricter:

- The model chooses the best feature set inside each outer training fold.
- Observed nested balanced accuracy across proxy views: 0.672-0.721.
- Fully nested shuffled-label null mean balanced accuracy: 0.492-0.514.
- All tested views beat the nested null with empirical p = 0.0323 and BH q = 0.0323.
- Best interpretation: the low-vs-zero classifier signal survives a stronger internal validation check, but this still is not external validation or clinical diagnosis.

The clinical/source-site covariate sensitivity check adds the strongest caveat so far:

- HER2-low and HER2-zero are highly imbalanced by TCGA source site and slide size.
- HER2-low mean slide file size is about 100 MB; HER2-zero mean slide file size is about 279 MB.
- Source-site group imbalance has chi-square p = 7.97e-10.
- In the top 8 CK proxy view, slide-size covariates alone reach balanced accuracy 0.879, source-site covariates alone reach 0.878, and source-site plus slide-size reaches 0.897.
- In the same view, GigaTIME mean channels reach balanced accuracy 0.745.
- Best interpretation: the current classifier signal is not safe to present as independent HER2 biology. It required source-site/slide-size matched sensitivity, and after that matched check it still needs external validation.

The matched low-versus-zero sensitivity check has now been run:

- Exact source-site matching gives 12 HER2-low/HER2-zero pairs.
- Slide-size matching gives 14 pairs with log-size caliper 0.25 and 20 pairs with log-size caliper 0.50.
- In exact source-site pairs, GigaTIME mean channels reach balanced accuracy 0.708 and AUC 0.750 in the top 8 CK proxy view.
- In the slide-size matched subsets, GigaTIME mean channels remain modestly above chance, with balanced accuracy 0.679 and 0.675.
- However, source-site and slide-size baselines remain competitive or stronger in the larger matched subsets.
- Paired channel tests in the matched subsets do not reach BH q < 0.05.
- Best interpretation: matching keeps the signal worth studying, but it does not make the TCGA result safe as independent HER2 biology. External/site-balanced validation and tumor-rich/pathologist review are still necessary.

The source-site held-out generalization check has now been run:

- In the top 8 CK proxy view, GigaTIME mean channels drop from balanced accuracy 0.745 under repeated stratified cross-validation to 0.669 when entire TCGA source sites are held out.
- In the same view, slide-size covariates remain strong under leave-source-site-out validation, with balanced accuracy 0.882 and AUC 0.915.
- GigaTIME mean-channel performance drops under source-site holdout across every tested proxy view.
- Best interpretation: the low-vs-zero classifier contains some above-chance image signal, but it is not source-independent. Slide-size/acquisition structure remains too strong for independent HER2 biology claims.

The HER2 isoform validation feasibility audit has also been run:

- The current local RNA files are GDC STAR augmented gene-count TSVs, not transcript-level isoform quantification.
- Local RNA data include 110 STAR gene-count cases, 56 strict high-trust cases with local STAR counts, and 40 HER2-low/HER2-zero high-trust cases with local STAR counts.
- There are no local BAM, FASTQ, junction-count, isoform, or transcript-level files under `data/tcga_brca`.
- Best interpretation: current data can support ERBB2 gene-level and RNA-program context, but cannot validate Guardia-style HER2 isoform biology. For that, we need sample-level isoform labels from the paper workflow or appropriate RNA-seq read/junction data.

See `docs/clinical_her2_high_trust_tile128_results.md`.

Key files:

- `data/tcga_brca/clinical_her2_laptop_balanced61_cases.csv`
- `data/tcga_brca/clinical_her2_laptop_balanced61_slides_files.csv`
- `data/tcga_brca/clinical_her2_laptop_balanced61_slide_manifest.tsv`
- `data/tcga_brca/clinical_her2_laptop_balanced61_slide_download_status.json`
- `docs/clinical_her2_laptop_balanced61_selection.md`
- `docs/clinical_her2_trustworthy_slide_list.md`
- `docs/assets/clinical_her2_trustworthy_slide_list/trustworthy_slides.csv`
- `docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/slide_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/tile_scores.csv`
- `docs/clinical_her2_high_trust_tile128_results.md`
- `docs/clinical_her2_high_trust_tile128_gigatime_data_cleanup.md`
- `docs/clinical_her2_high_trust_tile128_cleaned_classifier_comparison.md`
- `docs/clinical_her2_high_trust_tile128_case_driver_analysis.md`
- `docs/clinical_her2_high_trust_tile128_case_driver_visual_qc.md`
- `docs/clinical_her2_high_trust_tile128_tissue_composition_sensitivity.md`
- `docs/clinical_her2_high_trust_tile128_tumor_proxy_sensitivity.md`
- `docs/clinical_her2_high_trust_tile128_classifier_permutation_sanity.md`
- `docs/clinical_her2_high_trust_tile128_nested_classifier_model_selection.md`
- `docs/clinical_her2_high_trust_tile128_clinical_covariate_sensitivity.md`
- `docs/clinical_her2_high_trust_tile128_matched_low_zero_sensitivity.md`
- `docs/clinical_her2_high_trust_tile128_source_site_generalization.md`
- `docs/her2_isoform_validation_feasibility.md`
- `docs/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity.md`
- `docs/clinical_her2_high_trust_tile128_vs_expanded20_tile256_agreement.md`

## Main Outputs

- `data/tcga_brca/clinical_her2_labels.csv`
- `data/tcga_brca/clinical_her2_cohort_cases.csv`
- `data/tcga_brca/clinical_her2_cohort_slides_files.csv`
- `data/tcga_brca/clinical_her2_cohort_slide_download_status.json`
- `results/gigatime_tcga_brca_clinical_her2/slide_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2/tile_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2/heatmaps/`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_summary.md`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_channel_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_pairwise_tests.csv`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/*.png`
- `results/gigatime_tcga_brca_clinical_her2/rna_validation/gigatime_rna_signature_correlations.csv`
- `results/gigatime_tcga_brca_clinical_her2/rna_validation/rna_validation_summary.md`
- `results/gigatime_tcga_brca_clinical_her2_tile256/slide_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/clinical_summary/clinical_her2_summary.md`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_validation/rna_validation_summary.md`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/rna_program_validation_summary.md`
- `results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/classifier_baseline_summary.md`
- `docs/assets/clinical_her2_visual_qc/clinical_her2_visual_qc_selected_cases.csv`
- `docs/assets/clinical_her2_visual_qc/*_he_vs_virtual_mif_qc.png`
- `docs/assets/clinical_her2_tile256/`
- `docs/assets/clinical_her2_rna_program_validation/`
- `docs/assets/clinical_her2_classifier_baseline/`
- `docs/assets/clinical_her2_visual_qc_tile256/`
- `data/tcga_brca/clinical_her2_cohort_expanded20_cases.csv`
- `data/tcga_brca/erbb2_expression_expanded20.csv`
- `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/`
- `docs/clinical_her2_expanded20_results.md`
- `docs/assets/clinical_her2_expanded20_gigatime_cleanup/`
- `docs/assets/clinical_her2_expanded20_cleaned_classifier/`

The earlier ERBB2-high versus ERBB2-low pilot outputs are still present under `results/gigatime_tcga_brca_extremes/`, and the documentation-facing virtual mIF images are still under:

- `docs/assets/virtual_mif_channels/`
- `docs/assets/virtual_mif_composites/`

## Current Finding

The full 30-slide clinical HER2 pilot suggests a possible HER2-zero versus HER2-low immune-microenvironment signal.

The strongest three-group virtual-channel differences were:

| Channel | Kruskal p | Highest mean group | Lowest mean group |
|---|---:|---|---|
| CD68 | 0.0242 | HER2-zero | HER2-low |
| PD-L1 | 0.0423 | HER2-zero | HER2-low |
| CD11c | 0.0494 | HER2-zero | HER2-low |

Pairwise HER2-low versus HER2-zero comparisons were strongest for CD68, CD11c, PD-L1, CD4, and Ki67, but none remained significant after Benjamini-Hochberg correction.

## 256-Tile Robustness Finding

The same 30 selected slides were rerun with up to 256 random tissue tiles per slide. The main HER2-zero versus HER2-low pattern persisted and became slightly stronger for the leading channels:

| Channel | 64-tile p | 256-tile p | 64 max-min | 256 max-min | 256 direction |
|---|---:|---:|---:|---:|---|
| CD68 | 0.0242 | 0.0167 | 0.00913 | 0.01044 | HER2-zero > HER2-low |
| PD-L1 | 0.0423 | 0.0211 | 0.01749 | 0.02061 | HER2-zero > HER2-low |
| CD11c | 0.0494 | 0.0384 | 0.00450 | 0.00504 | HER2-zero > HER2-low |

The top 256-tile pairwise comparisons were again HER2-low versus HER2-zero. CD68, PD-L1, and CD11c had BH q values of 0.1133, improved from the 64-tile run but still not FDR-significant.

## RNA Validation Check

The first indirect validation layer compared GigaTIME virtual channels with matched RNA-seq marker signatures from the same 30 cases.

Result:

- No channel had an FDR-significant GigaTIME-versus-RNA signature correlation.
- `Ki67` had the strongest positive trend, Spearman rho 0.294.
- `CD68`, `PD-L1`, and `CD11c` did not show strong positive correlations with their matching RNA signatures.

Interpretation: the clinical HER2 virtual immune signal is interesting and robust to denser tile sampling, but not yet validated. It needs pathologist review and stronger orthogonal validation before making biological claims.

## Broader RNA Program Validation

The next validation step tested broader RNA immune and tissue programs instead of only marker-level RNA signatures.

Result:

- No broad RNA immune program showed an FDR-significant HER2-group difference.
- The virtual myeloid/checkpoint composite retained the HER2-zero > HER2-low direction, but did not pass FDR correction: p 0.0176, BH q 0.0878.
- The strongest FDR-significant virtual-vs-RNA associations were negative correlations with endothelial RNA signal:
  - Virtual T cell/checkpoint vs endothelial RNA: Spearman rho -0.585, BH q 0.0309.
  - Virtual all immune/checkpoint vs endothelial RNA: Spearman rho -0.556, BH q 0.0320.

Interpretation: the virtual signal is reproducible inside GigaTIME, but broader RNA validation still does not confirm it. This strengthens the need for pathologist review, tissue-composition checks, and external validation.

## First HER2 Classifier Baseline

The first diagnostic-model style classifier is now complete. It used slide-level GigaTIME features from the 256-tile run with leave-one-out cross-validation.

Best GigaTIME/H&E results:

| Task | Best GigaTIME feature set | Accuracy | Balanced accuracy | Macro AUC |
|---|---|---:|---:|---:|
| HER2-low vs HER2-zero | Mean + fraction channels | 0.800 | 0.800 | 0.870 |
| HER2-positive vs HER2-negative | Mean + fraction channels | 0.533 | 0.475 | 0.430 |
| Three-class HER2 group | Mean + fraction channels | 0.333 | 0.333 | 0.555 |

Interpretation: the classifier result is promising only for HER2-low versus HER2-zero in this tiny pilot. It is not reliable for HER2-positive detection or full three-class diagnosis.

## Pre-Classifier GigaTIME Data Cleanup

We then returned to the tile-level GigaTIME data to create cleaner feature views before retraining the classifier. The goal was to test whether the signal depends on all sampled tissue, cellular tissue, or more tumor/epithelial-enriched tiles.

Cleanup views:

| Cleanup view | Median retained tiles | Median retained fraction | Median DAPI | Median CK |
|---|---:|---:|---:|---:|
| All sampled tissue | 256.0 | 1.000 | 0.324 | 0.231 |
| QC cellular tissue | 190.5 | 0.744 | 0.360 | 0.249 |
| CK-enriched top 50% | 96.0 | 0.375 | 0.450 | 0.359 |
| CK-enriched top 25% | 48.0 | 0.188 | 0.493 | 0.431 |

Main cleanup interpretation:

- The HER2-zero > HER2-low CD68/PD-L1/CD11c signal persisted after cellular-tissue filtering.
- The signal weakened under stricter CK-enriched tile selection, especially in the top 25% CK view.
- This suggests the original GigaTIME signal is not simply blank-tile artifact, but it may depend partly on broader tissue context rather than only tumor-rich epithelial regions.

See `docs/clinical_her2_gigatime_data_cleanup.md`.

## Cleaned-View Classifier Comparison

The classifier was rerun separately on each cleaned GigaTIME feature view.

HER2-low versus HER2-zero result:

| Cleanup view | Best feature set | Accuracy | Balanced accuracy | Macro AUC |
|---|---|---:|---:|---:|
| All sampled tissue | Mean + fraction channels | 0.800 | 0.800 | 0.870 |
| QC cellular tissue | Mean + fraction channels | 0.800 | 0.800 | 0.900 |
| CK-enriched top 50% | Interpretable means | 0.650 | 0.650 | 0.670 |
| CK-enriched top 25% | Interpretable means | 0.650 | 0.650 | 0.630 |

Interpretation: cellular-tissue cleanup preserved the HER2-low versus HER2-zero classifier signal, which argues against blank/background tissue as the sole explanation. The signal weakened in CK-enriched views, suggesting the current GigaTIME signal may depend more on broader tissue or microenvironment context than on a purely epithelial tumor-cell HER2 phenotype.

HER2-positive versus HER2-negative performance remained weak. CK-enriched top 25% reached balanced accuracy 0.550, but sensitivity was only 0.200, so this is not useful for clinical HER2-positive detection.

See `docs/clinical_her2_cleaned_classifier_comparison.md`.

## Reframed Scientific Direction

The next proposal angle should be broader than "can the model classify HER2?" A stronger scientific question is whether image-derived GigaTIME features predict or associate with HER2-related biological states.

The especially interesting hypotheses are:

- HER2-low versus HER2-zero tumors with hidden or alternate ERBB2 transcript/isoform expression.
- HER2-positive tumors with image-derived states associated with trastuzumab or antibody-drug conjugate resistance.
- Tumors with preserved HER2 pathway signaling but reduced antibody targetability.

The current pilot does not prove any of these. It only suggests that the HER2-low versus HER2-zero boundary may have an image-derived tissue-context signal worth validating. We should say that image AI predicts or associates with HER2 isoform/state hypotheses, not that image AI detects HER2 isoforms.

See `docs/her2_isoform_state_hypothesis.md`.

## Local ERBB2 RNA Update

The latest local ERBB2 validation extracted gene-level ERBB2 TPM from all 110 local STAR count cases. In the strict high-trust overlap, ERBB2 RNA strongly supports broad HER2-positive status but weakly separates HER2-low from HER2-zero:

| Task | N cases | ERBB2-only AUC | Best-threshold balanced accuracy |
|---|---:|---:|---:|
| HER2-positive vs non-positive | 56 | 0.905 | 0.831 |
| HER2-low vs HER2-zero | 40 | 0.605 | 0.625 |

This helps the paper framing: the low/zero GigaTIME signal is not simply a strong gene-level ERBB2 RNA split. It still does not validate HER2 isoforms because these are gene-level STAR count files, not transcript-level isoform or junction files.

See `docs/clinical_her2_high_trust_tile128_local_erbb2_validation.md`.

## Within-Source-Site Update

We also restricted HER2-low versus HER2-zero analysis to TCGA source sites that contain both groups. Only four source sites qualify: A1, A2, A8, and AO. This leaves 51 cases total: 12 HER2-low and 39 HER2-zero.

In this mixed-site subset:

- Site-fixed channel tests keep 7 channel/view effects at BH q < 0.05.
- In the top 8 CK proxy view, GigaTIME all mean channels reach 0.667 repeated-CV balanced accuracy and 0.628 leave-mixed-source-site-out balanced accuracy.
- GigaTIME key mean channels are much weaker, around 0.505 repeated-CV and 0.490 leave-site-out balanced accuracy.

Interpretation: some all-channel GigaTIME signal remains inside mixed source sites, but the subset is too small and imbalanced to prove source-independent HER2 biology.

See `docs/clinical_her2_high_trust_tile128_within_source_site_low_zero.md`.

## Visual QC Check

The first visual QC pass selected the top `CD68` + `PD-L1` + `CD11c` case from each HER2 group:

| Clinical HER2 group | Selected case | Combined signal |
|---|---|---:|
| HER2-positive | TCGA-A2-A0EQ | 0.115 |
| HER2-low | TCGA-A2-A04Q | 0.086 |
| HER2-zero | TCGA-A2-A0T2 | 0.126 |

The high-scoring tiles were tissue-containing and cellular rather than obvious blank regions. This supports continued investigation, but it does not validate the virtual marker biology. The selected HER2-positive case also had visually plausible high-signal tiles, so the current result should be framed as a slide-level pilot trend rather than a clean visual separation.

In the 256-tile visual QC pass, the same representative cases were selected. The HER2-zero case `TCGA-A2-A0T2` had a higher combined `CD68` + `PD-L1` + `CD11c` signal than the selected HER2-positive and HER2-low representatives.

## Commands

The clinical HER2 run used:

```bash
conda run -n gigatime-tcga python scripts/run_gigatime_tcga_brca.py \
  --slide-table data/tcga_brca/clinical_her2_cohort_slides_files.csv \
  --missing-slide-policy skip \
  --out-dir results/gigatime_tcga_brca_clinical_her2 \
  --tile-limit 64 \
  --tile-order random \
  --batch-size 16 \
  --device auto \
  --save-tile-csv
```

The clinical HER2 summary used:

```bash
conda run -n gigatime-tcga python scripts/summarize_clinical_her2_gigatime.py \
  --slide-scores results/gigatime_tcga_brca_clinical_her2/slide_scores.csv \
  --cohort data/tcga_brca/clinical_her2_cohort_cases.csv \
  --out-dir results/gigatime_tcga_brca_clinical_her2/clinical_summary
```

The 256-tile robustness run used:

```bash
conda run -n gigatime-tcga python scripts/run_gigatime_tcga_brca.py \
  --slide-table data/tcga_brca/clinical_her2_cohort_slides_files.csv \
  --missing-slide-policy skip \
  --out-dir results/gigatime_tcga_brca_clinical_her2_tile256 \
  --tile-limit 256 \
  --tile-order random \
  --random-seed 42 \
  --batch-size 16 \
  --device auto \
  --save-tile-csv
```

## Caveat

This is still a pilot, not a definitive biological result. It is stronger than the first ERBB2-expression proof-of-work because it uses clinical HER2 groups and now includes both the original balanced 10/10/10 run and an expanded balanced 20/20/20 run. The expanded run strengthens the HER2-low versus HER2-zero image-signal argument, but marker-level and broader RNA-program validation remain weak. The classifier baselines are useful but not clinically reliable. The next scientific step is pathologist review, tumor-rich tile selection, stronger tissue QC, tumor-purity or immune-deconvolution adjustment, transcript/isoform-aware HER2 validation if available, and ideally an external dataset with real mIF or therapy-response data.
