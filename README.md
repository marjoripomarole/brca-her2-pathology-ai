# BRCA HER2 Pathology AI

Repository slug: `brca-her2-pathology-ai`

This repository is a computational pathology workspace for studying HER2-related breast cancer states in TCGA-BRCA diagnostic H&E whole-slide images. GigaTIME is the current first image model used here, but the project is intentionally broader than one model: it includes slide download, HER2 label cleanup, virtual mIF feature extraction, RNA/ERBB2 context checks, visual QC, confounder analyses, and exploratory classifiers.

The goal is to generate image-derived features from TCGA-BRCA diagnostic H&E whole-slide images, join those features to clinical HER2 labels and ERBB2 RNA expression, and produce advisor-ready summary tables and figures for HER2-positive, HER2-low, and HER2-zero research questions.

## What This Contains

- `external/GigaTIME/`: cloned official GigaTIME code from `prov-gigatime/GigaTIME`, currently used as the virtual mIF feature generator.
- `scripts/gdc_query_tcga_brca.py`: queries GDC for TCGA-BRCA diagnostic slides and STAR-count RNA-seq files, writes GDC manifests, and can extract ERBB2 expression.
- `scripts/build_tcga_brca_clinical_her2_labels.py`: queries the GDC TCGA-BRCA clinical supplement and builds reproducible clinical HER2-positive/HER2-low/HER2-zero labels.
- `scripts/select_clinical_her2_cohort.py`: selects a balanced clinical HER2-positive/HER2-low/HER2-zero cohort and writes a slide manifest for the next GigaTIME run.
- `scripts/download_clinical_her2_cohort_slides.py`: downloads the selected clinical HER2 cohort slide files from GDC by file ID and writes a resumable status JSON.
- `scripts/build_tcga_her2_trustworthy_slide_list.py`: checks the 61/61/61 cohort for HER2 label trust, slide metadata, file-size integrity, and OpenSlide readability, then writes high-trust slide lists.
- `scripts/run_gigatime_tcga_brca.py`: tiles TCGA-BRCA `.svs` slides, runs the official GigaTIME model, and aggregates virtual mIF channels per slide.
- `scripts/summarize_her2_gigatime.py`: joins GigaTIME slide scores with ERBB2 expression and makes HER2-high/HER2-low summary figures.
- `scripts/summarize_clinical_her2_gigatime.py`: compares GigaTIME virtual mIF outputs across clinical HER2-positive/HER2-low/HER2-zero groups.
- `scripts/validate_gigatime_with_rna_signatures.py`: compares GigaTIME virtual channels with matched RNA-seq marker signatures as an indirect validation check.
- `scripts/validate_gigatime_with_rna_programs.py`: compares GigaTIME virtual composite programs with broader RNA immune and tissue programs.
- `scripts/download_selected_star_counts.py`: downloads STAR-count RNA-seq files for a selected case list and extracts ERBB2 expression.
- `scripts/cleanup_gigatime_tile_features.py`: builds pre-classifier cleaned GigaTIME feature views from cellular and CK-enriched tile subsets.
- `scripts/train_her2_classifier_baseline.py`: trains first slide-level HER2 classifier baselines from GigaTIME features with leave-one-out cross-validation.
- `scripts/train_her2_cleaned_classifier_comparison.py`: reruns HER2 classifiers across all-tissue, cellular-tissue, and CK-enriched GigaTIME feature views.
- `scripts/analyze_high_trust_her2_sensitivity.py`: runs ER/PR-adjusted and HER2 IHC/ISH subgroup sensitivity checks for the high-trust result.
- `scripts/compare_gigatime_run_agreement.py`: compares overlapping slide-level GigaTIME outputs across the earlier 256-tile expanded run and current 128-tile high-trust run.
- `scripts/analyze_high_trust_case_drivers.py`: ranks case-level HER2-low versus HER2-zero driver scores, view stability, classifier errors, and manual-review priorities.
- `scripts/render_case_driver_visual_qc.py`: renders H&E plus virtual mIF panels for representative case-driver and opposite-profile review cases.
- `scripts/analyze_tissue_composition_sensitivity.py`: quantifies whether HER2-low versus HER2-zero GigaTIME signals track low-marker/stromal-like tissue composition.
- `scripts/analyze_tumor_proxy_sensitivity.py`: tests whether HER2-low versus HER2-zero statistics and classifiers survive stricter virtual tumor-rich proxy tile filters.
- `scripts/analyze_classifier_permutation_sanity.py`: runs shuffled-label sanity checks for the selected HER2-low versus HER2-zero classifiers.
- `scripts/analyze_nested_classifier_model_selection.py`: repeats low-vs-zero classifier feature-set selection inside cross-validation and compares it with a shuffled-label null.
- `scripts/analyze_clinical_covariate_sensitivity.py`: tests whether HER2-low versus HER2-zero findings are confounded by clinical, TCGA source-site, or slide-size covariates.
- `scripts/analyze_matched_low_zero_sensitivity.py`: builds source-site and slide-size matched HER2-low/HER2-zero sensitivity subsets after the confounder check.
- `scripts/analyze_source_site_generalization.py`: compares ordinary cross-validation with leave-TCGA-source-site-out validation for HER2-low versus HER2-zero classifiers.
- `scripts/analyze_within_source_site_low_zero.py`: restricts HER2-low/HER2-zero analysis to TCGA source sites containing both classes and runs site-fixed channel tests plus mixed-site classifiers.
- `scripts/audit_her2_isoform_validation_feasibility.py`: audits whether the current local RNA files can support Guardia-style HER2 isoform validation.
- `scripts/analyze_local_erbb2_expression_validation.py`: extracts ERBB2 gene-level TPM from all local STAR count files and checks whether clinical HER2/GigaTIME findings track gene-level ERBB2 expression.
- `scripts/run_hoptimus_tcga_brca.py`: extracts H-Optimus/H0-mini tile embeddings from TCGA-BRCA diagnostic H&E slides and aggregates them to slide-level embedding features.
- `scripts/render_virtual_mif_channel_images.py`: renders all-channel virtual mIF figures from GigaTIME tile and slide predictions.
- `scripts/render_virtual_mif_composites.py`: reruns GigaTIME on selected tiles and renders fluorescence-style virtual mIF composites from the full predicted channel maps.
- `scripts/render_clinical_her2_visual_qc.py`: renders clinical HER2 visual QC panels for cases driving high virtual `CD68`/`PD-L1`/`CD11c` signal.
- `scripts/build_clinical_her2_findings_report.py`: builds a simple display notebook and HTML report for the current clinical HER2 findings.
- `docs/virtual_mif_channel_outputs.md`: explains the generated virtual mIF channel images and how to interpret them.
- `docs/README.md`: start-here guide that separates current summaries from historical 30-slide reports.
- `docs/plain_language_methodology.md`: detailed non-specialist explanation of the study background, methodology, outputs, and current limitations.
- `docs/paper_proposal_process_log.md`: living process log for turning the pilot into a paper or grant proposal.
- `docs/clinical_her2_cohort_selection.md`: selected 30-case clinical HER2 pilot cohort and selection counts.
- `docs/clinical_her2_gigatime_run.md`: selected-cohort GigaTIME run status and full 30-slide clinical HER2 summary.
- `docs/clinical_her2_rna_validation.md`: first RNA-seq validation check for the clinical HER2 GigaTIME pilot.
- `docs/clinical_her2_visual_qc.md`: first visual/spatial QC pass for the clinical HER2 virtual immune-channel signal.
- `docs/clinical_her2_tile_sampling_robustness.md`: 256-tile robustness check showing whether the 64-tile HER2-zero versus HER2-low signal persists with denser sampling.
- `docs/clinical_her2_rna_program_validation.md`: broader RNA immune/tissue program validation after the 256-tile robustness run.
- `docs/clinical_her2_gigatime_data_cleanup.md`: pre-classifier tile cleanup using cellular tissue and virtual CK-enriched GigaTIME views.
- `docs/clinical_her2_classifier_baseline.md`: first diagnostic-model style classifier baseline for HER2-positive/negative, HER2-low/zero, and three-class HER2 prediction.
- `docs/clinical_her2_cleaned_classifier_comparison.md`: classifier comparison after GigaTIME tile cleanup and CK-enriched feature selection.
- `docs/clinical_her2_expanded20_results.md`: presentation-oriented summary of the expanded 20/20/20 clinical HER2 run.
- `docs/tcga_her2_label_quality_assessment.md`: source-backed assessment of TCGA-BRCA HER2 label quality, issues, assumptions, and checks.
- `docs/clinical_her2_trustworthy_slide_list.md`: trustworthy-slide list for the 183-slide 61/61/61 downloaded cohort.
- `docs/clinical_her2_high_trust_tile128_results.md`: current strict high-trust 171-slide GigaTIME/HER2 result.
- `docs/clinical_her2_high_trust_tile128_gigatime_data_cleanup.md`: current high-trust tile-cleanup sensitivity analysis.
- `docs/clinical_her2_high_trust_tile128_cleaned_classifier_comparison.md`: current high-trust cleaned classifier comparison.
- `docs/clinical_her2_high_trust_tile128_case_driver_analysis.md`: case-level HER2-low versus HER2-zero driver analysis and manual-review shortlist.
- `docs/clinical_her2_high_trust_tile128_case_driver_visual_qc.md`: small H&E plus virtual mIF visual QC set for the case-driver shortlist.
- `docs/clinical_her2_high_trust_tile128_tissue_composition_sensitivity.md`: quantifies the tissue-composition caveat for the HER2-low versus HER2-zero result.
- `docs/clinical_her2_high_trust_tile128_tumor_proxy_sensitivity.md`: stricter virtual CK/tumor-rich proxy sensitivity analysis for the HER2-low versus HER2-zero result.
- `docs/clinical_her2_high_trust_tile128_classifier_permutation_sanity.md`: shuffled-label sanity check for the selected HER2-low versus HER2-zero classifiers.
- `docs/clinical_her2_high_trust_tile128_nested_classifier_model_selection.md`: stricter nested model-selection classifier check for HER2-low versus HER2-zero.
- `docs/clinical_her2_high_trust_tile128_clinical_covariate_sensitivity.md`: clinical/source-site/slide-size confounder sensitivity for the HER2-low versus HER2-zero result.
- `docs/clinical_her2_high_trust_tile128_matched_low_zero_sensitivity.md`: source-site and slide-size matched sensitivity check for the HER2-low versus HER2-zero result.
- `docs/clinical_her2_high_trust_tile128_source_site_generalization.md`: leave-source-site-out validation for HER2-low versus HER2-zero classifiers.
- `docs/clinical_her2_high_trust_tile128_within_source_site_low_zero.md`: within-source-site low-vs-zero sensitivity using only TCGA source sites with both HER2-low and HER2-zero cases.
- `docs/clinical_her2_high_trust_tile128_local_erbb2_validation.md`: expanded local STAR gene-level ERBB2 validation for clinical HER2 labels and GigaTIME low/zero findings.
- `docs/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity.md`: ER/PR-adjusted and HER2 IHC/ISH subgroup sensitivity analysis for the current high-trust result.
- `docs/clinical_her2_high_trust_tile128_vs_expanded20_tile256_agreement.md`: channel-agreement and HER2-low/zero direction robustness across 128-tile and 256-tile runs.
- `docs/her2_isoform_validation_feasibility.md`: audit showing that current local STAR gene-count files support gene-level RNA context but not direct HER2 isoform validation.
- `docs/her2_isoform_state_hypothesis.md`: sharper paper-proposal framing around HER2 state, isoform hypotheses, targetability, and language guardrails.
- `docs/hoptimus_embedding_baseline.md`: first-run plan and commands for the H-Optimus/H0-mini embedding baseline.
- `docs/advisor_brief.md`: concise project framing and discussion points.
- `docs/current_pilot_run.md`: current two-case run status and advisor-facing caveats.
- `configs/tcga_brca_her2.yaml`: default paths and pilot settings.
- `notebooks/clinical_her2_findings_simple.ipynb` and `notebooks/clinical_her2_findings_simple.html`: simple presentation-ready summary of the findings so far.

## Requirements

GigaTIME requires access to the gated Hugging Face model. Accept the terms on the model card, then set a read-only token:

```bash
export HF_TOKEN=<huggingface_read_token>
```

Create the working environment:

```bash
conda env create -f envs/gigatime-tcga.yml
conda activate gigatime-tcga
```

The workflow is research-only and not for clinical decision-making, matching the GigaTIME model license notice.

## 1. Query TCGA-BRCA and Extract ERBB2

Start with a pilot subset before running all BRCA slides:

```bash
python scripts/gdc_query_tcga_brca.py \
  --out-dir data/tcga_brca \
  --case-limit 25 \
  --download-expression
```

This writes:

- `data/tcga_brca/tcga_brca_diagnostic_slides_manifest.tsv`
- `data/tcga_brca/tcga_brca_star_counts_manifest.tsv`
- `data/tcga_brca/erbb2_expression.csv`
- `data/tcga_brca/file_metadata_*.json`

To download slide files, either use the manifest with the GDC Data Transfer Tool:

```bash
gdc-client download \
  -m data/tcga_brca/tcga_brca_diagnostic_slides_manifest.tsv \
  -d data/tcga_brca/slides
```

or download a very small pilot directly:

```bash
python scripts/gdc_query_tcga_brca.py \
  --out-dir data/tcga_brca \
  --case-limit 5 \
  --download-slides \
  --max-slide-downloads 5 \
  --slide-download-order smallest
```

To pull one specific case, add `--slide-case-id TCGA-3C-AALI`.

## 2. Build Clinical HER2 Labels

For analyses that compare HER2-positive, HER2-low, and HER2-zero disease, build labels from the TCGA-BRCA clinical supplement:

```bash
conda run -n gigatime-tcga python scripts/build_tcga_brca_clinical_her2_labels.py
```

This writes:

- `data/tcga_brca/clinical_her2_labels.csv`
- `data/tcga_brca/clinical_her2_labels_metadata.json`
- `data/tcga_brca/clinical/nationwidechildrens.org_clinical_patient_brca.txt`

The label rules are:

- `HER2-positive`: IHC `3+`, ISH positive, or positive IHC receptor status when detailed score/ISH are missing.
- `HER2-low`: IHC `1+` with no positive ISH, or IHC `2+` with ISH negative.
- `HER2-zero`: IHC `0` with no positive ISH.
- `HER2-unknown`: missing, not evaluated, equivocal without definitive ISH, or otherwise incomplete HER2 fields.

## 3. Select a Balanced Clinical HER2 Cohort

After clinical labels are available, select a balanced 10/10/10 pilot cohort:

```bash
conda run -n gigatime-tcga python scripts/select_clinical_her2_cohort.py
```

This writes:

- `data/tcga_brca/clinical_her2_cohort_cases.csv`
- `data/tcga_brca/clinical_her2_cohort_slides_files.csv`
- `data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv`
- `data/tcga_brca/clinical_her2_cohort_summary.json`

The default selector chooses 10 cases per clinical group, prioritizing direct clinical HER2 labels, already-downloaded slides, smaller slide files, and deterministic case IDs.

To download the selected slide files with the GDC Data Transfer Tool:

```bash
gdc-client download \
  -m data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv \
  -d data/tcga_brca/slides
```

If `gdc-client` is not installed, use the project downloader:

```bash
conda run -n gigatime-tcga python scripts/download_clinical_her2_cohort_slides.py \
  --only-missing
```

This downloads each selected slide by GDC file ID into `data/tcga_brca/slides/<case>/` and writes:

- `data/tcga_brca/clinical_her2_cohort_slide_download_status.json`

## 4. Run GigaTIME on TCGA-BRCA Slides

```bash
python scripts/run_gigatime_tcga_brca.py \
  --slides-dir data/tcga_brca/slides \
  --out-dir results/gigatime_tcga_brca \
  --tile-limit 512 \
  --tile-order random \
  --batch-size 16 \
  --device auto \
  --save-tile-csv
```

Key output:

- `results/gigatime_tcga_brca/slide_scores.csv`
- `results/gigatime_tcga_brca/tile_scores.csv`
- `results/gigatime_tcga_brca/heatmaps/*.png`

For the first advisor meeting, `--tile-limit 512` is enough to demonstrate the pipeline. Increase or remove it for the full run.

## 4b. Extract H-Optimus/H0-Mini Embeddings

H-Optimus is the next generic H&E foundation-model baseline after the GigaTIME virtual mIF analysis. The first practical run should use the laptop-friendlier `h0-mini` preset, then optionally repeat with full `H-optimus-0` if the hardware and Hugging Face access are ready.

Both presets are gated on Hugging Face. Accept the model terms and set a read token:

```bash
export HF_TOKEN=<huggingface_read_token>
```

Check runtime and slide discovery:

```bash
conda run -n gigatime-tcga python scripts/run_hoptimus_tcga_brca.py --dry-run
```

Inspect one slide's extraction geometry without loading the model or running inference:

```bash
conda run -n gigatime-tcga python scripts/run_hoptimus_tcga_brca.py \
  --dry-run \
  --inspect-slide \
  --max-slides 1
```

Tiny smoke run:

```bash
conda run -n gigatime-tcga python scripts/run_hoptimus_tcga_brca.py \
  --model-preset h0-mini \
  --out-dir results/hoptimus_tcga_brca_high_trust_tile224_smoke \
  --max-slides 2 \
  --tile-limit 16 \
  --batch-size 4
```

Scaled high-trust run:

```bash
conda run -n gigatime-tcga python scripts/run_hoptimus_tcga_brca.py \
  --model-preset h0-mini \
  --out-dir results/hoptimus_tcga_brca_high_trust_tile224 \
  --tile-limit 64 \
  --batch-size 8 \
  --resume
```

See `docs/hoptimus_embedding_baseline.md` for the full H-Optimus-0 command and analysis plan.

To run GigaTIME only on the selected clinical HER2 cohort and skip slides that have not been downloaded yet:

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

The first robustness rerun used the same selected slides with denser sampling:

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

## 5. Summarize Results

```bash
python scripts/summarize_her2_gigatime.py \
  --slide-scores results/gigatime_tcga_brca/slide_scores.csv \
  --expression data/tcga_brca/erbb2_expression.csv \
  --out-dir results/gigatime_tcga_brca/advisor_summary
```

This writes joined data, channel-level HER2-high versus HER2-low summaries, figures, and `advisor_summary.md`.

For the clinical HER2-positive/HER2-low/HER2-zero cohort:

```bash
conda run -n gigatime-tcga python scripts/summarize_clinical_her2_gigatime.py \
  --slide-scores results/gigatime_tcga_brca_clinical_her2/slide_scores.csv \
  --cohort data/tcga_brca/clinical_her2_cohort_cases.csv \
  --out-dir results/gigatime_tcga_brca_clinical_her2/clinical_summary
```

To run the first indirect RNA-seq validation layer:

```bash
conda run -n gigatime-tcga python scripts/validate_gigatime_with_rna_signatures.py
```

This compares GigaTIME channels such as `CD68`, `PD-L1`, `CD11c`, and `Ki67` with simple matched RNA marker signatures from the available STAR-count files.

To run broader RNA program validation after the 256-tile clinical HER2 rerun:

```bash
conda run -n gigatime-tcga python scripts/validate_gigatime_with_rna_programs.py
```

This compares virtual composite programs such as myeloid/checkpoint and T-cell/checkpoint with broader RNA programs such as cytotoxic T-cell, checkpoint/IFNG, myeloid/macrophage, B-cell, stromal, endothelial, epithelial, and proliferation signatures.

To run the first slide-level HER2 classifier baseline:

```bash
conda run -n gigatime-tcga python scripts/train_her2_classifier_baseline.py
```

This trains regularized logistic and nearest-centroid baselines with leave-one-out cross-validation. It reports HER2-positive versus negative, HER2-low versus zero, and full three-class HER2 prediction performance.

To build cleaned pre-classifier GigaTIME feature views from the 256-tile output:

```bash
conda run -n gigatime-tcga python scripts/cleanup_gigatime_tile_features.py
```

This creates cellular-tissue and virtual CK-enriched slide feature tables under `results/gigatime_tcga_brca_clinical_her2_tile256/gigatime_cleanup/` and tracked figures under `docs/assets/clinical_her2_gigatime_cleanup/`.

To rerun HER2 classifiers across those cleaned feature views:

```bash
conda run -n gigatime-tcga python scripts/train_her2_cleaned_classifier_comparison.py
```

This writes cleaned-view classifier metrics under `results/gigatime_tcga_brca_clinical_her2_tile256/cleaned_classifier_comparison/` and tracked figures under `docs/assets/clinical_her2_cleaned_classifier/`.

## 6. Render All Virtual mIF Channel Images

```bash
conda run -n gigatime-tcga python scripts/render_virtual_mif_channel_images.py
```

This writes documentation-facing figures to `docs/assets/virtual_mif_channels/`, including all-channel group means, a slide-by-channel activation matrix, and HER2-high/HER2-low reference grids for the 23 GigaTIME virtual mIF channels. See `docs/virtual_mif_channel_outputs.md` for interpretation.

To create fluorescence-style virtual mIF images that look closer to real multiplex immunofluorescence panels:

```bash
conda run -n gigatime-tcga python scripts/render_virtual_mif_composites.py
```

This writes H&E-versus-virtual-mIF panels and marker-composite montages to `docs/assets/virtual_mif_composites/`. These are still GigaTIME predictions, not experimental mIF data.

To render the clinical HER2 visual QC panels for cases driving high virtual `CD68`, `PD-L1`, and `CD11c`:

```bash
conda run -n gigatime-tcga python scripts/render_clinical_her2_visual_qc.py
```

This writes tracked QC panels and selected-case tables to `docs/assets/clinical_her2_visual_qc/`.

To rebuild the simple presentation notebook and HTML report:

```bash
conda run -n gigatime-tcga python scripts/build_clinical_her2_findings_report.py
```

This writes:

- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

To rerun the shuffled-label classifier sanity check:

```bash
conda run -n gigatime-tcga python scripts/analyze_classifier_permutation_sanity.py
```

This writes the observed-versus-null classifier summary to `docs/clinical_her2_high_trust_tile128_classifier_permutation_sanity.md` and the permutation figures to `docs/assets/clinical_her2_high_trust_tile128_classifier_permutation/`.

To rerun the stricter nested classifier model-selection check:

```bash
conda run -n gigatime-tcga python scripts/analyze_nested_classifier_model_selection.py
```

This writes the nested observed-versus-null classifier summary to `docs/clinical_her2_high_trust_tile128_nested_classifier_model_selection.md` and figures to `docs/assets/clinical_her2_high_trust_tile128_nested_classifier/`.

To rerun the clinical/source-site/slide-size confounder check:

```bash
conda run -n gigatime-tcga python scripts/analyze_clinical_covariate_sensitivity.py
```

This writes the confounder report to `docs/clinical_her2_high_trust_tile128_clinical_covariate_sensitivity.md` and figures to `docs/assets/clinical_her2_high_trust_tile128_clinical_covariates/`.

To rerun the matched HER2-low versus HER2-zero sensitivity check:

```bash
conda run -n gigatime-tcga python scripts/analyze_matched_low_zero_sensitivity.py
```

This writes the matched-subset report to `docs/clinical_her2_high_trust_tile128_matched_low_zero_sensitivity.md` and figures to `docs/assets/clinical_her2_high_trust_tile128_matched_low_zero/`.

To rerun the source-site held-out generalization check:

```bash
conda run -n gigatime-tcga python scripts/analyze_source_site_generalization.py
```

This writes the source-site generalization report to `docs/clinical_her2_high_trust_tile128_source_site_generalization.md` and figures to `docs/assets/clinical_her2_high_trust_tile128_source_site_generalization/`.

To rerun the within-source-site low-vs-zero sensitivity check:

```bash
conda run -n gigatime-tcga python scripts/analyze_within_source_site_low_zero.py
```

This writes the mixed-site sensitivity report to `docs/clinical_her2_high_trust_tile128_within_source_site_low_zero.md`, figures to `docs/assets/clinical_her2_high_trust_tile128_within_source_site/`, and machine-readable outputs to `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/within_source_site_low_zero/`.

To rerun the HER2 isoform validation feasibility audit:

```bash
python3 scripts/audit_her2_isoform_validation_feasibility.py
```

This writes the feasibility report to `docs/her2_isoform_validation_feasibility.md` and machine-readable outputs to `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/her2_isoform_validation_feasibility/`.

To rerun the expanded local gene-level ERBB2 validation:

```bash
conda run -n gigatime-tcga python scripts/analyze_local_erbb2_expression_validation.py
```

This writes the ERBB2 validation report to `docs/clinical_her2_high_trust_tile128_local_erbb2_validation.md`, figures to `docs/assets/clinical_her2_high_trust_tile128_local_erbb2_validation/`, and machine-readable outputs to `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/local_erbb2_expression_validation/`.

## Notes for the Advisor Discussion

Current top-line result: the strict high-trust 171-slide tile128 analysis is the result to present. It uses 53 HER2-positive, 57 HER2-low, and 61 HER2-zero TCGA-BRCA diagnostic H&E slides after HER2 label, slide-integrity, OpenSlide, and female-patient primary-tumor filtering. The female-patient filter follows the relevant TCGA sample-selection principle from Guardia et al., Genome Research 2025, PMID 40664477; the remaining H&E file and slide QC is project-specific. Earlier bullets describe the precursor 30-slide and expanded 60-slide runs and why we scaled further.

Larger cohort status: a laptop-realistic 61/61/61 clinical HER2 cohort is downloaded locally, with 183 TCGA-BRCA diagnostic H&E slides total. This is the largest balanced cohort available from the current labels because HER2-zero has 61 candidate cases. After label, slide-integrity, OpenSlide, and female-patient QC, 171 strict high-trust slides are used as the current primary analysis set. The raw GigaTIME inference output contains the previous 174-slide set, and the current summaries filter that output to the 171 trustworthy slides.

- HER2 is represented here by `ERBB2` RNA expression from TCGA-BRCA STAR-count files.
- GigaTIME outputs virtual mIF maps for 23 channels, including immune markers such as `CD3`, `CD8`, `CD4`, `CD20`, `CD68`, `PD-1`, and `PD-L1`.
- The first deliverable is a replication/adaptation pilot, not a new model: run the released model on TCGA-BRCA H&E slides and ask whether virtual TIME signatures differ across clinical HER2 groups.
- The initial run should be treated as exploratory until tissue QC, slide-level aggregation, and HER2 clinical annotations are reviewed.
- The current clinical HER2 pilot has processed 30 selected slides: 10 HER2-positive, 10 HER2-low, and 10 HER2-zero. The strongest pilot signal is higher GigaTIME-predicted CD68, PD-L1, and CD11c in HER2-zero versus HER2-low, but these are hypothesis-generating and not FDR-significant after pairwise correction.
- The 256-tile robustness rerun reproduced the same HER2-zero greater than HER2-low direction for CD68, PD-L1, and CD11c. The leading pairwise q values improved to about 0.113 but remained above 0.05.
- The first RNA-seq validation check did not strongly confirm the virtual immune-channel signal; correlations between matched RNA marker signatures and GigaTIME channels were weak and not FDR-significant.
- Broader RNA program validation also did not positively confirm the virtual immune/checkpoint signal. The strongest FDR-significant associations were negative correlations between virtual immune/checkpoint programs and endothelial RNA signal.
- The first classifier baseline suggests possible GigaTIME signal for HER2-low versus HER2-zero, but not reliable HER2-positive/negative or three-class diagnosis. This is not clinically usable.
- The pre-classifier cleanup shows that the HER2-zero greater than HER2-low CD68/PD-L1/CD11c signal persists after cellular-tissue filtering, but weakens under strict CK-enriched tile selection. This suggests the original signal may depend partly on broader tissue context, not only tumor-rich tiles.
- The cleaned-view classifier comparison preserves HER2-low versus HER2-zero balanced accuracy at 0.800 after cellular-tissue filtering, but drops to 0.650 in CK-enriched views. This supports a microenvironment/tissue-context interpretation more than a purely tumor-epithelial HER2 classifier.
- The first visual QC pass found that high virtual CD68/PD-L1/CD11c tiles were tissue-containing and cellular rather than obvious blank background, but this still does not validate the virtual marker biology.
- The expanded 20/20/20 run processed 60 slides and strengthened the HER2-low versus HER2-zero signal. Several all-tissue or QC-cellular pairwise differences now pass within-view BH correction for `CD3`, `CD4`, `CD11c`, `CD68`, and QC-cellular `PD-L1`.
- In the expanded run, the best HER2-low versus HER2-zero GigaTIME/H&E classifier remained around balanced accuracy 0.800 and macro AUC 0.820. HER2-positive classification from GigaTIME/H&E remained weak.
- RNA marker and RNA program validation remain weak in the expanded run, so the strongest current claim is still hypothesis-generating image-derived HER2-state association, not clinical diagnosis.
- The strict high-trust 171-slide analysis uses only slides passing HER2 label, female-patient, slide-integrity, and readability checks. HER2-low remained lower than HER2-zero for multiple virtual immune/myeloid/checkpoint/tissue-context channels, including `CD68`, `PD-L1`, `CD11c`, `PD-1`, `CD4`, `CD3`, and `CK`, with all-sampled-tissue BH q values around 0.0016-0.0020 for the HER2-low versus HER2-zero pairwise tests.
- In the high-trust run, the best HER2-low versus HER2-zero GigaTIME/H&E classifier reached balanced accuracy 0.727 and macro AUC 0.787. HER2-positive classification remained weak, so this supports a HER2-low/HER2-zero tissue-context association rather than a diagnostic HER2 model.
- ER/PR-adjusted sensitivity checks support the same direction: in all sampled tissue, 7 of 8 tested key channels remain significant after ER/PR adjustment. The HER2-low versus HER2-zero signal is also visible across the main IHC/ISH detail subgroups.
- Parameter robustness supports the same interpretation: among 56 overlapping slides between the expanded 256-tile run and strict high-trust 128-tile analysis, all 8 tested key channels keep the same HER2-low versus HER2-zero direction and 7 of 8 keep HER2-low lower than HER2-zero.
- Tumor-rich proxy sensitivity is nuanced: fixed-count CK-rich tile views weaken individual low-versus-zero channel tests, but the multichannel low-versus-zero classifier remains around 0.708-0.761 balanced accuracy across the virtual tumor-rich proxy views. This keeps the signal worth pursuing, but it still needs real tumor-rich/pathologist validation.
- Classifier permutation sanity checks support that the low-versus-zero classifiers are not just fitting random labels: all selected views beat shuffled-label null distributions with empirical p = 0.0099 and BH q = 0.0099. This is still post-hoc exploratory evidence, not clinical validation.
- Nested classifier model-selection checks reduce feature-set selection bias: when the feature set is selected inside each training fold, low-vs-zero balanced accuracy remains about 0.672-0.721 across proxy views and beats nested shuffled-label null tests. This is stronger internal classifier evidence, but still not external validation.
- Clinical/source-site covariate sensitivity adds a major caveat: slide-size-only and source-site-only baselines classify HER2-low versus HER2-zero better than GigaTIME features. The follow-up matched sensitivity keeps modest GigaTIME performance in source-site/slide-size matched subsets, but source-site/slide-size baselines remain competitive or stronger and paired channel tests are not FDR-significant. This means the current TCGA classifier signal is still not safe to present as independent HER2 biology.
- Source-site held-out validation strengthens that caution: in the top 8 CK proxy view, GigaTIME mean channels drop from 0.745 balanced accuracy under repeated stratified CV to 0.669 when entire TCGA source sites are held out, while slide-size covariates remain about 0.882 balanced accuracy. The classifier signal is therefore not source-independent HER2 biology yet.
- Within-source-site sensitivity adds a smaller stress test: only four TCGA source sites contain both HER2-low and HER2-zero cases, with 12 low and 39 zero cases total. Site-fixed channel tests keep 7 channel/view effects at BH q < 0.05, and all-channel GigaTIME classifiers retain some above-chance mixed-site performance, but the subset is too small and imbalanced to rescue the current classifier as source-independent biology.
- Expanded local ERBB2 validation adds a useful RNA sanity check: ERBB2 gene expression strongly validates broad HER2-positive status, with ERBB2-only AUC 0.905 for HER2-positive versus non-positive, but weakly separates HER2-low from HER2-zero, with AUC 0.605 and pairwise p/q 0.262/0.262. The low/zero GigaTIME signal is therefore not simply a strong gene-level ERBB2 expression split.
- HER2 isoform feasibility audit adds a second guardrail: the current local RNA files are GDC STAR augmented gene-count TSVs. They support ERBB2 gene expression and RNA-program context, but not direct HER2 isoform quantification, SUPPA2 PSI, or rMATS junction confirmation. Isoform-state analysis requires sample-level isoform labels from the Guardia/Galante workflow or appropriate RNA-seq read/junction data.
- The sharper paper angle is inspired by Guardia et al.'s HER2 isoform/ADC-resistance biology: ask whether image-derived features predict or associate with HER2-related biological states, such as ERBB2 isoform/transcript context, signaling, or targetability. We should not claim that image AI detects HER2 isoforms without transcript-level or protein-level validation.
