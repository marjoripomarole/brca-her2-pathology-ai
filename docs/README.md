# Documentation Guide

Status: Start here. This is the navigation map for the docs folder.

Project: `brca-her2-pathology-ai`.

This folder contains both current summaries and historical analysis reports for the BRCA HER2 pathology-AI project. Start here if you are new to the project.

## Read First

1. `clinical_her2_high_trust_tile128_results.md`
   - Best current presentation summary.
   - Covers the strict high-trust 171-slide analysis set: 53 HER2-positive, 57 HER2-low, 61 HER2-zero.
   - Use this for the latest results.

2. `clinical_her2_expanded20_results.md`
   - Previous expanded 60-slide run: 20 HER2-positive, 20 HER2-low, 20 HER2-zero.
   - Useful as a denser 256-tile comparison point.

3. `advisor_brief.md`
   - Short advisor-facing summary.
   - Good for meetings and high-level discussion.

4. `plain_language_methodology.md`
   - Explains the project for someone without genetics, pathology, or AI background.
   - Best teaching document.

5. `paper_proposal_process_log.md`
   - The living history/process log.
   - This is the document that records what we did over time and why.

6. `her2_isoform_state_hypothesis.md`
   - The paper-proposal framing around HER2 state, isoform hypotheses, targetability, and careful language.

7. `tcga_her2_label_quality_assessment.md`
   - Explains how trustworthy the TCGA-BRCA HER2 labels are, what assumptions we make, and which sensitivity checks should be run.

8. `clinical_her2_trustworthy_slide_list.md`
   - Lists the 61/61/61 cohort slides by trust level after HER2 label checks, slide metadata checks, file-size checks, OpenSlide readability, and available GigaTIME tissue QC.

9. `clinical_her2_high_trust_tile128_case_driver_analysis.md`
   - New case-level HER2-low versus HER2-zero driver analysis.
   - Use this to find the strongest label-consistent examples and the highest-priority manual QC/pathology review cases.

10. `clinical_her2_high_trust_tile128_case_driver_visual_qc.md`
   - Small H&E plus virtual mIF visual QC set for the case-driver shortlist.
   - Important caveat: low-like tiles can be stromal/collagen-rich rather than clearly tumor-rich.

11. `clinical_her2_high_trust_tile128_tissue_composition_sensitivity.md`
   - Quantifies the tissue-composition caveat across all HER2-low/HER2-zero slides.
   - Shows that HER2-low has more low-marker tiles and that low-marker adjustment removes most channel effects.

12. `clinical_her2_high_trust_tile128_tumor_proxy_sensitivity.md`
   - Tests stricter virtual CK/tumor-rich proxy tile filters.
   - Shows a nuanced result: individual channel tests weaken, but exploratory low-vs-zero classifier performance remains above chance.

13. `clinical_her2_high_trust_tile128_classifier_permutation_sanity.md`
   - Shuffled-label sanity check for the selected HER2-low versus HER2-zero classifiers.
   - Shows the classifiers beat a random-label null, while remaining post-hoc and exploratory.

14. `clinical_her2_high_trust_tile128_nested_classifier_model_selection.md`
   - Stricter nested model-selection check for the HER2-low versus HER2-zero classifiers.
   - Shows the signal remains above shuffled-label null when feature-set selection is repeated inside validation folds.

15. `clinical_her2_high_trust_tile128_clinical_covariate_sensitivity.md`
   - Clinical/source-site/slide-size confounder sensitivity.
   - Important caution: source-site and slide-size covariates classify HER2-low versus HER2-zero better than GigaTIME features.

16. `clinical_her2_high_trust_tile128_matched_low_zero_sensitivity.md`
   - Matched HER2-low/HER2-zero sensitivity after the confounder check.
   - Important caution: matching keeps modest GigaTIME signal but does not remove TCGA source-site/slide-size concern.

17. `clinical_her2_high_trust_tile128_source_site_generalization.md`
   - Leave-TCGA-source-site-out classifier validation.
   - Important caution: GigaTIME performance drops under source-site holdout while slide-size covariates remain very strong.

18. `clinical_her2_high_trust_tile128_within_source_site_low_zero.md`
   - Restricts HER2-low versus HER2-zero analysis to TCGA source sites that contain both classes.
   - Useful stress test: some signal remains, but only 12 HER2-low and 39 HER2-zero cases qualify.

19. `clinical_her2_high_trust_tile128_local_erbb2_validation.md`
   - Expanded local STAR gene-level ERBB2 validation.
   - Useful sanity check: ERBB2 strongly supports broad HER2-positive labels, but weakly separates HER2-low from HER2-zero.

20. `her2_isoform_validation_feasibility.md`
   - Audits whether the current local RNA files can validate the Guardia et al. HER2 isoform biology.
   - Important caution: local STAR gene-count files support gene-level RNA context, not direct isoform or junction validation.

## Current Latest Results

The latest top-line result is the strict high-trust 171-slide clinical HER2 analysis:

- 171 TCGA-BRCA diagnostic H&E slides total after HER2 label, slide-readability, file-integrity, and female-patient filtering. Only the female-patient TCGA sample-selection principle is borrowed from Guardia et al., Genome Research 2025, PMID 40664477; the H&E slide QC is our project-specific check.
- 53 HER2-positive, 57 HER2-low, 61 HER2-zero.
- 128 random tissue tiles per slide.
- 21,888 primary-analysis tile predictions after filtering the existing 174-slide inference output.
- All slides are high label+slide trust after HER2 label QC, female-patient filtering, OpenSlide readability, and file-size checks.

The strongest current finding:

- GigaTIME/H&E features continue to separate HER2-low from HER2-zero.
- HER2-low is lower than HER2-zero for multiple virtual immune/myeloid/checkpoint and tissue-context channels, including `CD68`, `PD-L1`, `CD11c`, `PD-1`, `CD4`, `CD3`, and `CK`.
- The signal is directionally robust across overlapping slides from the earlier 60-slide 256-tile run and the current 171-slide 128-tile analysis set: 8 of 8 key channels have the same HER2-low versus HER2-zero direction, and 7 of 8 have HER2-low lower than HER2-zero in both runs.
- The all-sampled-tissue HER2-low versus HER2-zero signal mostly survives ER/PR adjustment: 7 of 8 tested key channels remain BH q < 0.05 after ER/PR adjustment.
- The HER2-low versus HER2-zero classifier reaches balanced accuracy about 0.727 and macro AUC about 0.787 with all sampled tissue features.
- A new case-level driver check shows that 71 of 118 HER2-low/HER2-zero slides match the expected low-versus-zero profile in at least 3 of 4 cleanup views, while 47 slides show the opposite profile in at least 2 views and need manual review.
- A small visual QC set suggests the low-like signal may partly reflect stromal/collagen-rich tissue composition; this needs pathologist review and tumor-rich tile restriction before stronger biological claims.
- The quantified tissue-composition analysis strengthens the caveat: HER2-low has a higher fraction of low-marker tiles than HER2-zero, and adjusting for low-marker tile fraction removes most low-versus-zero channel effects.
- The tumor-rich proxy analysis adds a useful counterbalance: strict fixed-count CK-rich views weaken univariate channel significance, but low-vs-zero classifier balanced accuracy remains around 0.708-0.761 across virtual tumor-rich proxy views.
- The classifier permutation sanity check adds a second counterbalance: all tested low-vs-zero classifier views beat shuffled-label null tests with empirical p = 0.0099 and BH q = 0.0099.
- The nested model-selection check is stricter: when the feature set is selected inside each training fold, low-vs-zero balanced accuracy remains about 0.672-0.721 across proxy views and still beats nested shuffled-label null tests.
- The clinical/source-site covariate check is the strongest current caveat: slide-size-only and source-site-only baselines classify HER2-low versus HER2-zero better than GigaTIME features. The classifier result may therefore be confounded by TCGA cohort construction or acquisition/site differences.
- The matched low-vs-zero sensitivity check reduces obvious source/size imbalance, but it does not rescue the classifier as independent biology: GigaTIME remains modestly above chance while site/slide-size baselines remain competitive or stronger, and paired channel tests do not reach FDR significance.
- The source-site generalization check is another strong caution: GigaTIME mean channels drop from 0.745 to 0.669 balanced accuracy in the top 8 CK proxy view when entire TCGA source sites are held out, while slide-size covariates remain about 0.882.
- The within-source-site sensitivity check is a small but useful stress test: only A1, A2, A8, and AO contain both low and zero cases, giving 12 HER2-low and 39 HER2-zero cases. Some site-fixed channel effects and all-channel classifiers remain above chance, but the subset is too small and imbalanced to prove source-independent biology.
- The expanded local ERBB2 validation adds RNA context: ERBB2-only AUC is 0.905 for HER2-positive versus non-positive, but only 0.605 for HER2-low versus HER2-zero, with low/zero pairwise p/q 0.262/0.262. This supports the clinical HER2-positive label but does not explain the low/zero image signal as a strong ERBB2 expression split.
- The HER2 isoform feasibility audit shows that our local RNA files are gene-level STAR count/TPM files, not transcript-level isoform or junction-count data. We need sample-level isoform labels or read/junction data before testing the Guardia-style HER2 isoform hypothesis.
- RNA validation remains weak, so this is still hypothesis-generating and not clinical diagnosis.

## Larger Cohort Status

For laptop-realistic scaling, the project downloaded a larger balanced clinical HER2 cohort:

- 183 TCGA-BRCA diagnostic H&E slides total.
- 61 HER2-positive, 61 HER2-low, and 61 HER2-zero cases.
- This is the largest balanced set possible with the current clinical HER2 labels because HER2-zero has 61 candidate cases.
- All 183 selected slide files are present locally under `data/tcga_brca/slides/`.
- The cohort definition is in `clinical_her2_laptop_balanced61_selection.md`.

HER2 label, slide-file, OpenSlide, and Guardia et al.-aligned sex/gender QC marked 171 of these slides as strict high label+slide trust. Those 171 high-trust slides have now replaced the 60-slide result as the current primary analysis.

## Current Detailed Reports

- `clinical_her2_high_trust_tile128_results.md`: current best larger strict high-trust 171-slide findings summary.
- `clinical_her2_expanded20_results.md`: previous expanded 60-slide findings summary.
- `clinical_her2_cohort_expanded20_selection.md`: previous 20/20/20 cohort selection.
- `clinical_her2_expanded20_gigatime_data_cleanup.md`: previous expanded cleanup/tile-filtering report.
- `clinical_her2_expanded20_cleaned_classifier_comparison.md`: previous expanded cleaned classifier comparison.
- `clinical_her2_laptop_balanced61_selection.md`: larger 61/61/61 downloaded cohort selection.
- `clinical_her2_trustworthy_slide_list.md`: trustworthy-slide list for the 183 downloaded slides.
- `clinical_her2_high_trust_tile128_gigatime_data_cleanup.md`: high-trust cleanup/tile-filtering report.
- `clinical_her2_high_trust_tile128_cleaned_classifier_comparison.md`: high-trust cleaned classifier comparison.
- `clinical_her2_high_trust_tile128_case_driver_analysis.md`: case-level driver score, view-stability check, classifier-error review list, and manual-review shortlist.
- `clinical_her2_high_trust_tile128_case_driver_visual_qc.md`: H&E plus virtual mIF visual QC panels for representative label-consistent and opposite-profile cases.
- `clinical_her2_high_trust_tile128_tissue_composition_sensitivity.md`: low-marker/CK-high tissue-composition sensitivity analysis.
- `clinical_her2_high_trust_tile128_tumor_proxy_sensitivity.md`: virtual tumor-rich proxy statistics and classifier sensitivity.
- `clinical_her2_high_trust_tile128_classifier_permutation_sanity.md`: shuffled-label sanity check for the selected low-vs-zero classifiers.
- `clinical_her2_high_trust_tile128_nested_classifier_model_selection.md`: nested feature-set selection and shuffled-label null check for the low-vs-zero classifiers.
- `clinical_her2_high_trust_tile128_clinical_covariate_sensitivity.md`: clinical/source-site/slide-size confounder sensitivity.
- `clinical_her2_high_trust_tile128_matched_low_zero_sensitivity.md`: matched source-site/slide-size sensitivity after the confounder finding.
- `clinical_her2_high_trust_tile128_source_site_generalization.md`: leave-source-site-out classifier validation after the confounder finding.
- `clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity.md`: ER/PR-adjusted and HER2 IHC/ISH subgroup sensitivity analysis.
- `clinical_her2_high_trust_tile128_vs_expanded20_tile256_agreement.md`: parameter/settings robustness check comparing the 128-tile high-trust run against the overlapping 256-tile expanded run.
- `her2_isoform_validation_feasibility.md`: practical feasibility audit for testing HER2 isoform/state biology with current local files.

## Historical 30-Slide Reports

These are still useful because they show how the project developed, but they should not be cited as the latest result:

- `clinical_her2_cohort_selection.md`: original 10/10/10 cohort.
- `clinical_her2_gigatime_run.md`: original 30-slide GigaTIME run.
- `clinical_her2_tile_sampling_robustness.md`: 30-slide 256-tile robustness check.
- `clinical_her2_rna_validation.md`: 30-slide marker-level RNA validation.
- `clinical_her2_rna_program_validation.md`: 30-slide broader RNA program validation.
- `clinical_her2_classifier_baseline.md`: 30-slide first classifier baseline.
- `clinical_her2_gigatime_data_cleanup.md`: 30-slide cleanup report.
- `clinical_her2_cleaned_classifier_comparison.md`: 30-slide cleaned classifier report.
- `clinical_her2_visual_qc.md`: initial 30-slide visual QC.

## Visual Explanation Files

- `virtual_mif_channel_outputs.md`: explains virtual mIF-style images and channel visualizations.
- `assets/`: tracked figures used by the markdown reports.

## How To Use This Folder

For a new researcher:

1. Read `clinical_her2_high_trust_tile128_results.md`.
2. Read `plain_language_methodology.md` if the biology or workflow is unfamiliar.
3. Read `paper_proposal_process_log.md` to understand the history.
4. Use historical reports only when you need details from a specific earlier analysis.

For a presentation:

1. Lead with `clinical_her2_high_trust_tile128_results.md`.
2. Use `advisor_brief.md` for a concise narrative.
3. Use `her2_isoform_state_hypothesis.md` for careful biological framing.
4. Use `her2_isoform_validation_feasibility.md` if someone asks whether the current data can validate HER2 isoforms.
5. Use `clinical_her2_high_trust_tile128_case_driver_analysis.md` if someone asks which cases should be inspected next.
