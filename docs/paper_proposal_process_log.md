# Paper Proposal Process Log

Status: Living history/process log. This is the document that records the sequence of decisions, analyses, and findings over time.

Last updated: 2026-06-02

This document keeps a running record of the research process for a future paper or grant proposal. It is written to preserve both the scientific reasoning and the concrete computational steps used in this BRCA HER2 pathology-AI project.

## Working Project Title

BRCA HER2 Pathology AI: evaluating image-derived tumor microenvironment features across the HER2 axis in TCGA breast cancer.

## Central Research Question

Can computational pathology methods, starting with the existing histopathology model GigaTIME, infer biologically interpretable virtual multiplex immunofluorescence features from TCGA-BRCA H&E slides, and do those predicted features differ across HER2-related breast cancer states?

The current project began with an ERBB2 RNA-expression pilot. The proposed next version should move toward clinically meaningful HER2 groups:

- HER2-positive
- HER2-low
- HER2-zero

This distinction matters because HER2-low and HER2-zero are defined by clinical protein-level testing, not simply by the amount of ERBB2 RNA measured in sequencing data.

## Plain-Language Framing

Breast cancer slides stained with H&E show tissue structure, but they do not directly show all immune-cell markers or molecular markers. GigaTIME is a released AI model that takes H&E image tiles and predicts virtual multiplex immunofluorescence-like marker maps. In this project, those predictions are used as research features that may describe the tumor immune microenvironment.

The study asks whether these predicted immune and tumor marker patterns vary along the HER2 axis. HER2 is clinically important in breast cancer because it affects biology and treatment options. The project is not trying to diagnose HER2 from H&E. Instead, it asks whether a previously released AI model produces tissue-microenvironment signals that correlate with known HER2 biology.

## Data Sources

### TCGA-BRCA

The data source is TCGA-BRCA, the breast invasive carcinoma project from The Cancer Genome Atlas. The files are accessed through the Genomic Data Commons (GDC).

The current workflow uses:

- Diagnostic H&E whole-slide images in `.svs` format.
- RNA-seq STAR-count files used to extract ERBB2 expression.
- Public clinical supplement files used to investigate HER2 IHC/ISH status.

Official references:

- GDC API overview: https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/
- GDC clinical supplement description: https://docs.gdc.cancer.gov/Encyclopedia/pages/Clinical_Supplement/

### GigaTIME

The project uses the official released GigaTIME implementation and model weights. The model predicts 23 virtual mIF channels from H&E pathology image tiles:

`DAPI`, `TRITC`, `Cy5`, `PD-1`, `CD14`, `CD4`, `T-bet`, `CD34`, `CD68`, `CD16`, `CD11c`, `CD138`, `CD20`, `CD3`, `CD8`, `PD-L1`, `CK`, `Ki67`, `Tryptase`, `Actin-D`, `Caspase3-D`, `PHH3-B`, and `Transgelin`.

These outputs are model predictions, not laboratory-measured mIF.

## Methods Completed So Far

### 1. Queried TCGA-BRCA Files From GDC

The script `scripts/gdc_query_tcga_brca.py` was created to query GDC for:

- TCGA-BRCA diagnostic H&E slide images.
- TCGA-BRCA STAR-count RNA-seq expression files.

Main local outputs:

- `data/tcga_brca/tcga_brca_diagnostic_slides_manifest.tsv`
- `data/tcga_brca/tcga_brca_diagnostic_slides_files.csv`
- `data/tcga_brca/tcga_brca_star_counts_manifest.tsv`
- `data/tcga_brca/tcga_brca_star_counts_files.csv`
- `data/tcga_brca/file_metadata_slides.json`
- `data/tcga_brca/file_metadata_star_counts.json`

### 2. Extracted ERBB2 Expression

The workflow downloaded selected TCGA-BRCA STAR-count files and extracted ERBB2 expression from each file.

The ERBB2 gene was identified as:

- Gene symbol: `ERBB2`
- Ensembl gene ID: `ENSG00000141736`

Main local output:

- `data/tcga_brca/erbb2_expression.csv`

Current status:

- ERBB2 expression was extracted for 80 TCGA-BRCA cases.
- These expression values were used as the first HER2-biology proxy.

### 3. Selected ERBB2-Extreme Cases

The script `scripts/select_her2_extremes.py` selected the top and bottom ERBB2 TPM cases from the current expression pilot.

Main local output:

- `data/tcga_brca/her2_extreme_cases.csv`

Current selected cohort:

- 20 ERBB2-high cases.
- 20 ERBB2-low cases.

Important caveat: these labels are expression-based research labels. They do not mean confirmed clinical HER2-positive, HER2-low, or HER2-zero status.

### 4. Downloaded and Processed a Slide Subset

The workflow attempted to download H&E whole-slide images for the selected ERBB2-extreme cases. Slide downloads were slow and occasionally unstable, so the current processed subset is smaller than the target cohort.

Current processed pilot:

- 12 TCGA-BRCA diagnostic slides.
- 7 ERBB2-high slides.
- 5 ERBB2-low slides.
- 64 random tissue tiles per slide.
- 768 total tile predictions.
- CPU inference.

Main local outputs:

- `results/gigatime_tcga_brca_extremes/slide_scores.csv`
- `results/gigatime_tcga_brca_extremes/tile_scores.csv`
- `results/gigatime_tcga_brca_extremes/heatmaps/`

### 5. Ran GigaTIME on H&E Tiles

The script `scripts/run_gigatime_tcga_brca.py` tiles each slide, filters for tissue-containing tiles, runs GigaTIME, and aggregates virtual mIF channel predictions.

For each processed slide, the workflow stores:

- Tile-level channel activations.
- Slide-level mean activations.
- Slide-level fraction-positive or thresholded summaries for each channel.
- Heatmap-style visual outputs.

### 6. Summarized Virtual mIF Features by ERBB2 Group

The script `scripts/summarize_her2_gigatime.py` joins slide-level GigaTIME predictions to ERBB2 expression and compares ERBB2-high versus ERBB2-low groups.

Main local outputs:

- `results/gigatime_tcga_brca_extremes/advisor_summary/joined_slide_her2_gigatime.csv`
- `results/gigatime_tcga_brca_extremes/advisor_summary/her2_group_channel_summary.csv`
- `results/gigatime_tcga_brca_extremes/advisor_summary/advisor_summary.md`
- Summary figures in `results/gigatime_tcga_brca_extremes/advisor_summary/`

This analysis is exploratory because the current processed subset contains only 12 slides.

### 7. Generated Virtual mIF Channel Figures

The script `scripts/render_virtual_mif_channel_images.py` renders documentation-facing all-channel GigaTIME figures.

Main local outputs:

- `docs/assets/virtual_mif_channels/virtual_mif_all_channel_group_means.png`
- `docs/assets/virtual_mif_channels/virtual_mif_slide_channel_matrix.png`
- `docs/assets/virtual_mif_channels/her2_high_reference_all_virtual_mif_channels.png`
- `docs/assets/virtual_mif_channels/her2_low_reference_all_virtual_mif_channels.png`

These figures show predicted channels across groups, slides, and spatial tile positions.

### 8. Generated Fluorescence-Style Virtual mIF Composites

The script `scripts/render_virtual_mif_composites.py` reruns GigaTIME on selected tiles and combines predicted channel maps into black-background fluorescence-style composites.

Main local outputs:

- `docs/assets/virtual_mif_composites/her2_high_immune_checkpoint_virtual_mif_montage.png`
- `docs/assets/virtual_mif_composites/her2_low_immune_checkpoint_virtual_mif_montage.png`
- `docs/assets/virtual_mif_composites/her2_high_immune_checkpoint_he_vs_virtual_mif.png`
- `docs/assets/virtual_mif_composites/her2_low_immune_checkpoint_he_vs_virtual_mif.png`
- Additional tumor/proliferation and myeloid/B-cell virtual panels.

These images are closer in appearance to real mIF images than the dot-grid figures, but they remain virtual predictions from H&E.

### 9. Built a Reproducible Clinical HER2 Label Table

The script `scripts/build_tcga_brca_clinical_her2_labels.py` now queries the GDC TCGA-BRCA clinical supplement, downloads the patient-level BCR Biotab, extracts HER2 IHC/ISH fields, and writes a traceable clinical HER2 label table.

Command:

```bash
conda run -n gigatime-tcga python scripts/build_tcga_brca_clinical_her2_labels.py
```

Main local outputs:

- `data/tcga_brca/clinical_her2_labels.csv`
- `data/tcga_brca/clinical_her2_labels_metadata.json`
- `data/tcga_brca/clinical/nationwidechildrens.org_clinical_patient_brca.txt`

The generated label table includes one row per TCGA-BRCA clinical case, the raw HER2 IHC/ISH values, the assigned clinical HER2 group, and the exact rule used to assign that group.

### 10. Selected a Balanced Clinical HER2 Pilot Cohort

The script `scripts/select_clinical_her2_cohort.py` now joins the clinical HER2 label table with ERBB2 expression and slide metadata, then selects a deterministic balanced cohort for the next GigaTIME run.

Command:

```bash
conda run -n gigatime-tcga python scripts/select_clinical_her2_cohort.py
```

Main local outputs:

- `data/tcga_brca/clinical_her2_cohort_cases.csv`
- `data/tcga_brca/clinical_her2_cohort_slides_files.csv`
- `data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv`
- `data/tcga_brca/clinical_her2_cohort_summary.json`
- `docs/clinical_her2_cohort_selection.md`

Selection priority:

- Clinical HER2 group must be one of HER2-positive, HER2-low, or HER2-zero.
- Direct clinical HER2 labels are preferred over inferred labels.
- Already-downloaded slides are preferred.
- Smaller slide files are preferred to make the next pilot more practical.
- Case IDs are used for deterministic tie-breaking.

Selected cohort after downloading the selected slides:

| Cohort group | Selected cases | Slides now downloaded |
|---|---:|---:|
| HER2-positive | 10 | 10 |
| HER2-low | 10 | 10 |
| HER2-zero | 10 | 10 |

This gives the first clean 30-case clinical HER2 pilot cohort for running GigaTIME across HER2-positive, HER2-low, and HER2-zero disease.

### 11. Ran the Availability-Limited Selected-Slide Clinical HER2 GigaTIME Pilot

The script `scripts/run_gigatime_tcga_brca.py` now supports `--slide-table`, allowing GigaTIME to process only the selected clinical HER2 cohort slides instead of scanning every slide under `data/tcga_brca/slides`.

The first selected-cohort run used the local slides already available:

- 8 selected slides processed.
- 4 HER2-positive slides.
- 3 HER2-low slides.
- 1 HER2-zero slide.
- 22 selected slides still missing locally.

Command:

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

The clinical three-group summary is generated by `scripts/summarize_clinical_her2_gigatime.py`.

Historical preliminary result:

- No definitive group difference should be claimed yet because HER2-zero has only 1 processed slide.
- The strongest availability-limited three-group signal among summarized channels was CD4, but it was not statistically significant in this small subset.

### 12. Downloaded the Remaining Clinical HER2 Slides and Ran the Full 30-Slide Pilot

The selected cohort was completed locally by downloading the 22 missing selected slides with:

```bash
conda run -n gigatime-tcga python scripts/download_clinical_her2_cohort_slides.py \
  --only-missing
```

Main local output:

- `data/tcga_brca/clinical_her2_cohort_slide_download_status.json`

The full selected clinical HER2 pilot was then rerun with:

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

Full clinical HER2 pilot status:

- 30 selected slides processed.
- 10 HER2-positive slides.
- 10 HER2-low slides.
- 10 HER2-zero slides.
- 64 random tissue tiles per slide.

The full clinical three-group summary was generated by:

```bash
conda run -n gigatime-tcga python scripts/summarize_clinical_her2_gigatime.py \
  --slide-scores results/gigatime_tcga_brca_clinical_her2/slide_scores.csv \
  --cohort data/tcga_brca/clinical_her2_cohort_cases.csv \
  --out-dir results/gigatime_tcga_brca_clinical_her2/clinical_summary
```

Current full-pilot result:

- The top unadjusted three-group differences were CD68, PD-L1, and CD11c.
- For these channels, HER2-zero had the highest mean virtual signal and HER2-low had the lowest mean virtual signal.
- HER2-positive was generally intermediate rather than clearly separated from HER2-low.
- Pairwise HER2-low versus HER2-zero tests were strongest for CD68, CD11c, PD-L1, CD4, and Ki67.
- No pairwise comparison remained significant after Benjamini-Hochberg correction, so the result is hypothesis-generating rather than definitive.

### 13. Compared GigaTIME Virtual Channels With RNA-Seq Marker Signatures

Because matched real mIF is not currently available for the TCGA slides in this project, the first indirect validation layer compared GigaTIME virtual-channel scores with matched bulk RNA-seq marker signatures.

Command:

```bash
conda run -n gigatime-tcga python scripts/validate_gigatime_with_rna_signatures.py
```

Main local outputs:

- `results/gigatime_tcga_brca_clinical_her2/rna_validation/case_rna_signatures.csv`
- `results/gigatime_tcga_brca_clinical_her2/rna_validation/joined_gigatime_rna_signatures.csv`
- `results/gigatime_tcga_brca_clinical_her2/rna_validation/gigatime_rna_signature_correlations.csv`
- `results/gigatime_tcga_brca_clinical_her2/rna_validation/gigatime_rna_group_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2/rna_validation/rna_validation_summary.md`

Tracked documentation:

- `docs/clinical_her2_rna_validation.md`
- `docs/assets/clinical_her2_rna_validation/gigatime_rna_correlation_heatmap.png`
- `docs/assets/clinical_her2_rna_validation/top_gigatime_rna_signature_scatter.png`

RNA validation result:

- All 30 clinical HER2 pilot cases had matched RNA-seq files available locally.
- No GigaTIME channel had an FDR-significant correlation with its RNA marker signature.
- `Ki67` had the strongest positive trend, Spearman rho 0.294, but this was not significant after correction.
- The main virtual immune-signal channels, CD68, PD-L1, and CD11c, did not show strong positive RNA-signature correlations.

Interpretation:

- The HER2-zero versus HER2-low GigaTIME immune/checkpoint signal is not yet validated by bulk RNA-seq.
- This does not prove the virtual signal is wrong, because bulk RNA-seq and H&E tile-level virtual mIF measure different tissue layers.
- The result raises the bar for the next step: visual QC, more tile sampling per slide, and stronger orthogonal validation are needed before making biological claims.

### 14. Rendered Visual QC Panels for High Virtual Immune-Channel Cases

The first visual/spatial QC pass selected the top case from each clinical HER2 group by combined GigaTIME signal:

```text
mean_CD68 + mean_PD-L1 + mean_CD11c
```

Command:

```bash
conda run -n gigatime-tcga python scripts/render_clinical_her2_visual_qc.py
```

Selected cases:

| Clinical HER2 group | Selected case | Combined signal | mean CD68 | mean PD-L1 | mean CD11c |
|---|---|---:|---:|---:|---:|
| HER2-positive | TCGA-A2-A0EQ | 0.115 | 0.029 | 0.072 | 0.014 |
| HER2-low | TCGA-A2-A04Q | 0.086 | 0.018 | 0.058 | 0.010 |
| HER2-zero | TCGA-A2-A0T2 | 0.126 | 0.037 | 0.069 | 0.021 |

Tracked outputs:

- `docs/clinical_her2_visual_qc.md`
- `docs/assets/clinical_her2_visual_qc/clinical_her2_visual_qc_selected_cases.csv`
- `docs/assets/clinical_her2_visual_qc/clinical_her2_visual_qc_manifest.csv`
- `docs/assets/clinical_her2_visual_qc/*_he_vs_virtual_mif_qc.png`
- `docs/assets/clinical_her2_visual_qc/*_sampled_tile_overlay.png`

Visual QC result:

- The high-scoring tiles were tissue-containing and cellular rather than obvious blank background.
- The HER2-zero selected case had the highest combined slide-level signal among the selected group representatives.
- The selected HER2-positive case also had visually plausible high-signal tiles, so the pattern is not unique to HER2-zero at the tile level.
- This supports continued investigation but does not validate the virtual marker biology.

### 15. Built a Simple Display Notebook and HTML Report

To make the current findings easier to present, a simple display notebook and HTML report were generated from the current clinical HER2 result tables and tracked figure assets.

Command:

```bash
conda run -n gigatime-tcga python scripts/build_clinical_her2_findings_report.py
```

Tracked outputs:

- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`
- `docs/assets/clinical_her2_findings/clinical_her2_channel_boxplots.png`
- `docs/assets/clinical_her2_findings/clinical_her2_group_mean_heatmap.png`
- `docs/assets/clinical_her2_findings/erbb2_tpm_by_clinical_her2_group.png`

The report is intentionally simple. It emphasizes:

- The balanced 10/10/10 clinical HER2 design.
- The leading HER2-zero versus HER2-low virtual immune/checkpoint signal.
- The weak RNA validation result.
- The visual QC result.
- The correct cautious proposal language.

See `docs/clinical_her2_gigatime_run.md` for the exact commands, local output paths, and current pilot table.

### 16. Ran a 256-Tile Robustness Check on the Same Clinical HER2 Cohort

The next robustness step was completed by rerunning the same 30 selected slides with up to 256 random tissue tiles per slide.

Command:

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

Then the same clinical summary, RNA validation, visual QC, and display-report steps were repeated.

Main local outputs:

- `results/gigatime_tcga_brca_clinical_her2_tile256/slide_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/tile_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/clinical_summary/clinical_her2_summary.md`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_validation/rna_validation_summary.md`
- `docs/assets/clinical_her2_tile256/`
- `docs/assets/clinical_her2_visual_qc_tile256/`
- `docs/clinical_her2_tile_sampling_robustness.md`

Robustness result:

| Channel | 64-tile p | 256-tile p | 64 max-min | 256 max-min | Direction |
|---|---:|---:|---:|---:|---|
| CD68 | 0.0242 | 0.0167 | 0.00913 | 0.01044 | HER2-zero > HER2-low |
| PD-L1 | 0.0423 | 0.0211 | 0.01749 | 0.02061 | HER2-zero > HER2-low |
| CD11c | 0.0494 | 0.0384 | 0.00450 | 0.00504 | HER2-zero > HER2-low |

The leading pairwise q values improved to 0.1133 for CD68, PD-L1, and CD11c, but remained above 0.05. RNA validation remained weak and no channel had an FDR-significant correlation with matched RNA marker signatures.

Interpretation:

- The HER2-zero greater than HER2-low virtual immune/checkpoint signal is now more robust to tile sampling.
- The result is still not biologically validated.
- The next proposal step should emphasize pathologist review and orthogonal validation rather than simply claiming a HER2 biology discovery.

### 17. Tested Broader RNA Immune and Tissue Programs

After the marker-level RNA validation remained weak, the next validation step compared GigaTIME virtual composite programs with broader RNA programs.

Command:

```bash
conda run -n gigatime-tcga python scripts/validate_gigatime_with_rna_programs.py
```

Main local outputs:

- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/case_rna_programs.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/case_virtual_programs.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/joined_virtual_rna_programs.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/virtual_rna_program_correlations.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/rna_program_group_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/virtual_program_group_summary.csv`
- `docs/clinical_her2_rna_program_validation.md`
- `docs/assets/clinical_her2_rna_program_validation/`

Virtual programs tested:

- Myeloid/checkpoint: `CD68`, `CD11c`, `PD-L1`
- T cell/checkpoint: `CD3`, `CD4`, `CD8`, `PD-1`
- All immune/checkpoint: `CD3`, `CD4`, `CD8`, `CD20`, `CD68`, `CD11c`, `PD-1`, `PD-L1`
- Proliferation: `Ki67`
- Epithelial: `CK`

RNA programs tested:

- T cell/cytotoxic
- Checkpoint/IFNG
- Myeloid/macrophage
- Dendritic/APC
- B cell
- Proliferation
- Epithelial/tumor
- Stromal/fibroblast
- Endothelial

Main result:

- The virtual myeloid/checkpoint composite retained the HER2-zero greater than HER2-low direction, but did not pass FDR correction: p 0.0176, BH q 0.0878.
- No broad RNA immune program showed an FDR-significant HER2-group difference.
- The strongest FDR-significant virtual-vs-RNA associations were negative correlations with endothelial RNA signal:
  - Virtual T cell/checkpoint versus endothelial RNA: Spearman rho -0.585, BH q 0.0309.
  - Virtual all immune/checkpoint versus endothelial RNA: Spearman rho -0.556, BH q 0.0320.

Interpretation:

- The virtual signal is reproducible within GigaTIME and stable across tile sampling.
- The signal is still not validated by orthogonal RNA evidence.
- The endothelial negative correlations raise a tissue-composition concern that should be reviewed before any strong biological claim.

### 18. Trained a First Slide-Level HER2 Classifier Baseline

The next methodological step moved beyond group-average comparisons. A first classifier baseline was trained to ask whether slide-level GigaTIME features can predict held-out clinical HER2 labels.

Command:

```bash
conda run -n gigatime-tcga python scripts/train_her2_classifier_baseline.py
```

Main local outputs:

- `results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/classifier_crossval_predictions.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/classifier_metrics.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/classifier_confusion_matrices.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/classifier_baseline_summary.md`
- `docs/clinical_her2_classifier_baseline.md`
- `docs/assets/clinical_her2_classifier_baseline/`

Model setup:

- Input features: slide-level GigaTIME virtual channel means, thresholded fraction-positive channel summaries, interpretable marker subsets, and composite virtual immune/tumor programs from the 256-tile run.
- Reference feature: ERBB2 RNA TPM, included only as a non-H&E benchmark.
- Output labels: clinical HER2-positive, HER2-low, and HER2-zero groups from TCGA IHC/ISH fields.
- Models: regularized logistic classifier and nearest-centroid baseline.
- Evaluation: leave-one-out cross-validation, with accuracy, balanced accuracy, macro AUC, sensitivity, specificity, and confusion matrices.

Main GigaTIME/H&E result:

| Task | Best GigaTIME/H&E feature set | Accuracy | Balanced accuracy | Macro AUC |
|---|---|---:|---:|---:|
| HER2-low vs HER2-zero | GigaTIME mean + fraction channels | 0.800 | 0.800 | 0.870 |
| HER2-positive vs HER2-negative | GigaTIME mean + fraction channels | 0.533 | 0.475 | 0.430 |
| Three-class HER2 group | GigaTIME mean + fraction channels | 0.333 | 0.333 | 0.555 |

The ERBB2 RNA reference classified HER2-positive versus HER2-negative much better than GigaTIME/H&E features, with balanced accuracy 0.850 and macro AUC 0.800. This is useful because it shows the clinical labels contain molecular HER2 signal, but the current H&E-derived features are not capturing the HER2-positive diagnostic signal reliably.

Interpretation:

- The first classifier framework now works end to end.
- GigaTIME features look most promising for the subtle HER2-low versus HER2-zero comparison.
- GigaTIME features do not yet reliably detect HER2-positive disease.
- Full three-class clinical HER2 prediction is at chance in this tiny pilot.
- This is a feasibility and failure-mode analysis, not a diagnostic model.

### 19. Returned to Pre-Classifier GigaTIME Data Cleanup

After the first classifier baseline, we returned to the GigaTIME tile-level data to ask whether the input features were too broad. The classifier baseline averaged all sampled tissue tiles, which can mix tumor, stroma, immune regions, normal tissue, and other non-tumor context. Because HER2 is clinically assessed in tumor cells, the next cleanup step created more biologically focused feature views before retraining any classifier.

Command:

```bash
conda run -n gigatime-tcga python scripts/cleanup_gigatime_tile_features.py
```

Main local outputs:

- `results/gigatime_tcga_brca_clinical_her2_tile256/gigatime_cleanup/tile_qc_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/gigatime_cleanup/cleaned_slide_features.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/gigatime_cleanup/filter_retention_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/gigatime_cleanup/cleanup_channel_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/gigatime_cleanup/cleanup_pairwise_tests.csv`
- `docs/clinical_her2_gigatime_data_cleanup.md`
- `docs/assets/clinical_her2_gigatime_cleanup/`

Cleanup views:

- All sampled tissue tiles from the original 256-tile run.
- QC cellular tissue: tissue fraction at least 0.70 and virtual DAPI mean at least 0.05.
- CK-enriched top 50%: the top half of virtual CK tiles within each slide after QC.
- CK-enriched top 25%: the top quarter of virtual CK tiles within each slide after QC.

Important caveat: virtual DAPI and virtual CK are still GigaTIME predictions from H&E, not real stains or pathologist-annotated tumor masks. These views are tumor-enriched research feature views, not confirmed tumor regions.

Tile retention:

| Cleanup view | Median retained tiles | Median retained fraction | Median DAPI | Median CK |
|---|---:|---:|---:|---:|
| All sampled tissue | 256.0 | 1.000 | 0.324 | 0.231 |
| QC cellular tissue | 190.5 | 0.744 | 0.360 | 0.249 |
| CK-enriched top 50% | 96.0 | 0.375 | 0.450 | 0.359 |
| CK-enriched top 25% | 48.0 | 0.188 | 0.493 | 0.431 |

Main result:

- The HER2-zero greater than HER2-low CD68/PD-L1/CD11c signal persisted after cellular-tissue QC and became slightly stronger by mean difference.
- The same signal weakened under stricter CK-enriched tile selection, especially in the top 25% CK view.
- This suggests the original signal is not only blank-tile artifact, but it may depend partly on broader tissue context rather than only tumor-rich epithelial regions.

Interpretation:

- The cleaned feature tables are now ready for a second classifier run.
- If the HER2-low versus HER2-zero classifier remains strong using QC cellular tissue but weakens with CK-enriched views, that would argue the model is learning tissue-context or microenvironment signal rather than a tumor-cell HER2 phenotype.
- If CK-enriched views improve classification, that would support the idea that tumor-region GigaTIME features are more relevant for HER2 prediction.

### 20. Reran HER2 Classifiers Across Cleaned GigaTIME Views

The second classifier step used the cleaned feature views from step 19 and reran the same leave-one-out classifier evaluation separately for each view.

Command:

```bash
conda run -n gigatime-tcga python scripts/train_her2_cleaned_classifier_comparison.py
```

Main local outputs:

- `results/gigatime_tcga_brca_clinical_her2_tile256/cleaned_classifier_comparison/cleaned_classifier_predictions.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/cleaned_classifier_comparison/cleaned_classifier_metrics.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/cleaned_classifier_comparison/cleaned_classifier_confusion_matrices.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/cleaned_classifier_comparison/cleaned_classifier_best_h_e_metrics.csv`
- `docs/clinical_her2_cleaned_classifier_comparison.md`
- `docs/assets/clinical_her2_cleaned_classifier/`

HER2-low versus HER2-zero result:

| Cleanup view | Best feature set | Accuracy | Balanced accuracy | Macro AUC |
|---|---|---:|---:|---:|
| All sampled tissue | Mean + fraction channels | 0.800 | 0.800 | 0.870 |
| QC cellular tissue | Mean + fraction channels | 0.800 | 0.800 | 0.900 |
| CK-enriched top 50% | Interpretable means | 0.650 | 0.650 | 0.670 |
| CK-enriched top 25% | Interpretable means | 0.650 | 0.650 | 0.630 |

HER2-positive versus HER2-negative remained weak. The CK-enriched top 25% view reached balanced accuracy 0.550, but sensitivity was only 0.200, so this is not clinically useful for HER2-positive detection.

Interpretation:

- Cellular-tissue cleanup preserved HER2-low versus HER2-zero performance, arguing against blank/background artifact as the sole explanation.
- Strict CK enrichment weakened the HER2-low versus HER2-zero classifier, suggesting that the current GigaTIME signal may depend more on broader tissue or microenvironment context than on purely epithelial tumor-cell features.
- Full three-class prediction remained near chance across views.
- The next useful analysis is visual inspection of cases whose predictions change between all-tissue/QC-cellular views and CK-enriched views.

### 21. Expanded the Clinical HER2 Cohort to 20/20/20

We then expanded the balanced clinical HER2 run from 30 slides to 60 slides:

- 20 HER2-positive cases.
- 20 HER2-low cases.
- 20 HER2-zero cases.
- Up to 256 random tissue tiles per slide.
- 15,225 total GigaTIME tile predictions.
- STAR-count RNA-seq expression downloaded for all 60 selected cases.

New helper added:

```bash
conda run -n gigatime-tcga python scripts/download_selected_star_counts.py \
  --star-counts data/tcga_brca_full_query/tcga_brca_star_counts_files.csv \
  --cases data/tcga_brca/clinical_her2_cohort_expanded20_cases.csv \
  --expression-out data/tcga_brca/erbb2_expression_expanded20.csv \
  --status-out data/tcga_brca/clinical_her2_cohort_expanded20_star_counts_download_status.json
```

Main expanded-run outputs:

- `docs/clinical_her2_expanded20_results.md`
- `docs/clinical_her2_cohort_expanded20_selection.md`
- `docs/clinical_her2_expanded20_gigatime_data_cleanup.md`
- `docs/clinical_her2_expanded20_cleaned_classifier_comparison.md`
- `results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/`
- `docs/assets/clinical_her2_expanded20_gigatime_cleanup/`
- `docs/assets/clinical_her2_expanded20_cleaned_classifier/`

Strongest updated finding:

- The HER2-low versus HER2-zero image-derived signal persisted and became stronger after expansion.
- Several HER2-low versus HER2-zero pairwise channel differences now pass within-view BH correction in all-tissue or QC-cellular views:
  - `CD4`: q 0.0252 all tissue; q 0.0200 QC-cellular.
  - `CD3`: q 0.0252 all tissue; q 0.0200 QC-cellular.
  - `CD11c`: q 0.0252 all tissue; q 0.0206 QC-cellular.
  - `CD68`: q 0.0326 all tissue; q 0.0320 QC-cellular.
  - `PD-L1`: q 0.0320 in QC-cellular tissue.

Expanded classifier result:

| Input view | Best feature set | N | Balanced accuracy | Macro AUC |
|---|---|---:|---:|---:|
| All sampled tissue | Mean + fraction channels | 40 | 0.800 | 0.820 |
| QC cellular tissue | Mean + fraction channels | 40 | 0.775 | 0.820 |
| CK-enriched top 50% | Mean channels | 40 | 0.750 | 0.807 |
| CK-enriched top 25% | Mean + fraction channels | 40 | 0.800 | 0.820 |

Expanded interpretation:

- The strongest current presentation result is now the expanded HER2-low versus HER2-zero signal, not only the original 30-slide pilot.
- HER2-low often appears to be the lowest virtual immune/checkpoint group.
- HER2-positive becomes highest for several broader virtual immune programs, so the three-group pattern is not simply "HER2-zero is highest for every immune marker."
- RNA marker and RNA program validation remain weak, so the result is still hypothesis-generating and should be described as an image-derived HER2-state association, not clinical diagnosis.
- GigaTIME/H&E still does not reliably classify HER2-positive disease; ERBB2 RNA performs better for that task, as expected.

## Initial Biological Findings From the ERBB2-Extreme Pilot

The current processed dataset is too small for strong claims. The main result so far is that the workflow is feasible and produces interpretable tables and figures.

Current working interpretation:

- GigaTIME can be run on TCGA-BRCA H&E slides.
- Its 23 predicted virtual mIF channels can be aggregated by slide and compared with ERBB2 expression.
- The current figures are useful for advisor discussion and proposal planning.
- The current analysis should not yet be framed as a definitive HER2 biological result.

The correct proposal language is "exploratory pilot" rather than "validated biomarker analysis."

## Investigation of Clinical HER2-Zero, HER2-Low, and HER2-Positive Labels

The advisor-facing research direction requires true clinical HER2 groups, not only ERBB2 RNA extremes.

HER2-low is generally defined as:

- IHC `1+`, or
- IHC `2+` with ISH negative.

The FDA describes HER2-low breast cancer using this IHC/ISH definition in the context of trastuzumab deruxtecan eligibility. Reference: https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-fam-trastuzumab-deruxtecan-nxki-unresectable-or-metastatic-hr-positive-her2-low-or-her2

### Clinical Fields Found in TCGA-BRCA

The TCGA-BRCA GDC clinical patient Biotab contains HER2-related fields, including:

- `lab_proc_her2_neu_immunohistochemistry_receptor_status`
- `her2_erbb_pos_finding_cell_percent_category`
- `her2_immunohistochemistry_level_result`
- `pos_finding_her2_erbb2_other_measurement_scale_text`
- `lab_procedure_her2_neu_in_situ_hybrid_outcome_type`
- `her2_neu_chromosone_17_signal_ratio_value`

This means TCGA-BRCA can support a better clinical HER2 grouping than the original ERBB2-expression-only pilot.

### Implemented Clinical HER2 Mapping

The implemented mapping is:

- `HER2-positive`: IHC `3+`, ISH positive, or receptor status positive when detailed fields are missing.
- `HER2-low`: IHC `1+`, or IHC `2+` with ISH negative.
- `HER2-zero`: IHC `0` with no positive ISH.
- `HER2-unknown`: missing, not evaluated, contradictory, or incomplete HER2 data.

This mapping should be documented carefully in the methods section because TCGA clinical supplement fields are not always complete.

### Counts From the Clinical HER2 Label Table

Using the implemented mapping above:

| Dataset | HER2-positive | HER2-low | HER2-zero | HER2-unknown |
|---|---:|---:|---:|---:|
| TCGA-BRCA clinical rows | 174 | 407 | 61 | 455 |
| Cases with current slide metadata and ERBB2 expression | 13 | 35 | 10 | 22 |
| Current 40 ERBB2-extreme selected cases | 13 | 10 | 8 | 9 |
| Current 12 GigaTIME-processed slides | 6 | 3 | 1 | 2 |
| Current 30-slide clinical HER2 GigaTIME pilot | 10 | 10 | 10 | 0 |
| Expanded 60-slide clinical HER2 GigaTIME run | 20 | 20 | 20 | 0 |

Interpretation:

- TCGA-BRCA appears useful for clinical HER2 grouping.
- The 30-slide clinical HER2 pilot now supports a first balanced three-group comparison.
- The expanded 60-slide run strengthens the HER2-low versus HER2-zero image-signal argument, but still requires molecular, pathology, and external validation.

## Multiplex Immunofluorescence Comparison Question

We investigated whether TCGA has real multiplex immunofluorescence results matched to the H&E slides used here. The current working conclusion is that this project does not yet have matched real mIF for those TCGA slides.

Therefore:

- We cannot currently validate GigaTIME virtual mIF predictions by direct TCGA matched mIF comparison.
- The virtual mIF images should be presented as model-generated predictions, not as experimental ground truth.
- Trustworthiness needs to be assessed using indirect and external validation strategies.

Potential validation strategies:

- Compare GigaTIME-predicted immune channels with bulk expression signatures from RNA-seq.
- Compare predicted epithelial/tumor channels with pathology/tumor-purity annotations where available.
- Compare HER2 group trends with known breast cancer biology.
- Use an external dataset with paired H&E and real mIF, if available.
- Perform manual pathology review of selected H&E and virtual mIF panels.
- Check whether tile-level high-signal regions correspond to plausible tissue structures.

## Proposed Next Analyses

### Analysis 1: Re-select Cases by Clinical HER2 Group

This is now completed for both a first balanced 10/10/10 clinical HER2 pilot and an expanded balanced 20/20/20 run. Instead of selecting only ERBB2-high and ERBB2-low cases, the workflow now creates clinical HER2 cohorts:

- HER2-positive cases.
- HER2-low cases.
- HER2-zero cases.

The selection balances slide availability and chooses one primary-tumor slide per case whenever possible.

### Analysis 2: Run GigaTIME on a Larger Clinical HER2 Cohort

The first balanced 30-slide clinical HER2 pilot and the expanded balanced 60-slide clinical HER2 run are now complete. The proposal should now target validation rather than another simple expansion.

Completed expanded run:

- 20 HER2-positive slides.
- 20 HER2-low slides.
- 20 HER2-zero slides.

Possible future target, only if time and compute permit:

- Larger balanced cohort with more HER2-zero cases.
- More complete whole-slide sampling.
- External validation cohort with H&E plus real mIF, IHC/ISH, or treatment-response data.

### Analysis 3: Compare Virtual mIF Features Across Clinical HER2 Groups

For each GigaTIME channel, compare slide-level mean activations across HER2 groups.

Candidate statistical tests:

- Kruskal-Wallis test for three-group comparison.
- Pairwise Wilcoxon rank-sum tests for HER2-positive versus HER2-low, HER2-low versus HER2-zero, and HER2-positive versus HER2-zero.
- Benjamini-Hochberg FDR correction across channels.
- Effect sizes with confidence intervals where possible.

Primary endpoints should be pre-specified before scaling:

- Immune checkpoint channels: `PD-1`, `PD-L1`.
- T-cell channels: `CD3`, `CD8`, `CD4`.
- Myeloid/macrophage channels: `CD68`, `CD11c`, `CD14`.
- Tumor/proliferation channels: `CK`, `Ki67`.

### Analysis 4: Compare Virtual mIF With RNA-Seq Immune Signatures

Because matched real mIF is not currently available, use RNA-seq as an indirect validation layer.

Examples:

- Compare virtual `CD3` or `CD8` channels with T-cell gene signatures.
- Compare virtual `PD-L1` with `CD274` expression.
- Compare virtual macrophage-like channels with macrophage-related genes or signatures.
- Compare virtual `Ki67` with proliferation-related genes such as `MKI67`.

This does not prove the virtual mIF is correct, but it helps determine whether predicted tissue signals are directionally consistent with molecular data.

The first implementation of this RNA validation layer is complete for simple marker signatures and was repeated after the 256-tile rerun. It did not strongly validate the current virtual immune-channel pattern, so future validation should consider richer immune signatures, tumor purity adjustment, external data, and visual review.

### Analysis 5: Visual QC and Trustworthiness Review

For selected cases from each HER2 group:

- Show source H&E tiles.
- Show GigaTIME virtual mIF composites.
- Show all-channel spatial maps.
- Flag artifacts, blank tissue, necrosis, folds, staining variation, and suspicious tile predictions.

This step is important because model outputs can look polished while still being wrong. The proposal should explicitly state that image-level QC is part of the methodology.

The first 256-tile visual QC repeated the same representative cases and again showed tissue-containing high-signal tiles. The next visual step should be human review by an advisor/pathologist, not only automated figure generation.

### Analysis 6: Train and Evaluate HER2 Classifiers

This is now implemented as a first baseline. The current baseline uses slide-level aggregate features, which is useful for feasibility but not sufficient for a final diagnostic model.

The next classifier versions should:

- Restrict inputs to tumor-rich tiles rather than all tissue tiles.
- Add tile distribution features such as percentiles, maximum signal, and spatial heterogeneity.
- Add H&E tile embeddings if available from GigaTIME or another pathology foundation model.
- Use multiple-instance learning once tile-level features or embeddings are organized reliably.
- Use nested cross-validation or a separate held-out test set before reporting tuned model performance.
- Report confusion matrices, AUC, sensitivity, specificity, and calibrated probabilities for each HER2 task.

## Paper Proposal Structure

### Background

HER2 is a clinically important axis in breast cancer. The emergence of HER2-low as a therapeutically relevant category creates a need to understand whether tumor morphology and microenvironmental patterns differ across HER2-positive, HER2-low, and HER2-zero disease.

### Gap

Traditional H&E slides are widely available, but they do not directly provide multiplex immune-marker information. Real mIF is informative but expensive and not routinely available for large public cohorts. Virtual mIF models may provide a scalable way to generate hypotheses about immune microenvironment differences from existing pathology slides.

The sharper biology gap is that HER2 categories are often treated as if they only represent "how much HER2" is present. For a stronger paper, we should ask whether image-derived features associate with HER2-related biological states, such as ERBB2 transcript/isoform context, preserved signaling, altered antibody targetability, or treatment-resistance hypotheses.

### HER2 Isoform/State Hypothesis

The most interesting future direction is not to claim that image AI sees HER2 isoforms directly. The safer and more scientifically useful question is whether image-derived features predict or associate with HER2 isoform/state hypotheses.

Potential high-impact hypotheses:

- HER2-low versus HER2-zero tumors may differ in hidden or alternate ERBB2 transcript/isoform expression.
- Some HER2-positive tumors may have image-derived tissue states associated with trastuzumab or antibody-drug conjugate resistance.
- Some tumors may preserve HER2 pathway signaling while having reduced antibody targetability.

Language guardrails:

- Use: "associated with," "predicts," "stratifies," and "image-derived correlate of HER2 state."
- Avoid: "detects HER2 isoforms," "diagnoses isoforms," "directly measures targetability," or "proves therapy resistance."

Validation required:

- ERBB2 transcript-level or isoform-aware quantification if available.
- Protein-level or antibody-based validation using IHC, ISH, real mIF, proteomics, or similar data.
- External therapy-response cohorts for trastuzumab or ADC resistance questions, because TCGA-BRCA alone is not enough to test treatment resistance.

See `docs/her2_isoform_state_hypothesis.md` for the current working version of this framing.

### Objective

Evaluate whether GigaTIME-derived virtual mIF features from TCGA-BRCA H&E slides differ across clinical HER2 groups and whether those predicted features associate with molecular, clinical, or HER2-related state annotations.

### Methods Overview

1. Retrieve TCGA-BRCA H&E slides, RNA-seq data, and clinical HER2 supplement fields from GDC.
2. Assign clinical HER2-positive, HER2-low, HER2-zero, or unknown labels using IHC/ISH fields.
3. Run GigaTIME on tissue-containing H&E tiles from selected diagnostic slides.
4. Aggregate virtual mIF channel predictions at tile and slide levels.
5. Compare GigaTIME channels across HER2 clinical groups.
6. Perform indirect validation against RNA-seq immune and proliferation signatures.
7. Generate visual QC panels and virtual mIF composites for interpretability.
8. Train cross-validated slide-level classifier baselines and evaluate diagnostic-model failure modes.
9. If data permits, test whether image-derived features associate with ERBB2 transcript-level, isoform-aware, targetability-related, or therapy-response evidence.

### Expected Contribution

This study would not claim that GigaTIME diagnoses HER2 status or detects HER2 isoforms. Instead, it would evaluate whether a released virtual mIF model can produce biologically interpretable tissue and immune-context features from public breast cancer H&E slides, whether those features vary across clinically meaningful HER2 categories, and whether they can be developed into hypotheses about HER2-related biological state.

## 2026-06-02: TCGA HER2 Label And Slide Trustworthiness QC

We reviewed TCGA-BRCA HER2 label quality against GDC/TCGA clinical documentation and, after advisor clarification, the relevant Galante-lab paper by Guardia et al. (Genome Research 2025, PMID `40664477`). The useful methodological lesson for this image project is that TCGA-derived analyses should use explicit traceable inputs, restrict the primary cohort to appropriate primary tumor samples, account for hormone-receptor context, exclude male TCGA-BRCA samples for consistency with the HER2 isoform paper, and rely on sensitivity analyses rather than a single fragile model run.

We then created a reproducible trustworthiness script:

```text
scripts/build_tcga_her2_trustworthy_slide_list.py
```

The script checks:

- Direct versus inferred clinical HER2 label rules.
- IHC/ISH discordance flags.
- HER2-low subgroups: IHC `1+` versus IHC `2+`/ISH-negative.
- HER2-zero subgroups: IHC `0`/ISH-negative versus IHC `0`/ISH-not-evaluated.
- ER/PR and ERBB2 RNA context as possible confounders.
- Primary-tumor slide metadata.
- SVS file existence and exact file-size match.
- OpenSlide readability.
- Existing GigaTIME/tile-level tissue QC where available.

Outputs:

- `docs/clinical_her2_trustworthy_slide_list.md`
- `docs/assets/clinical_her2_trustworthy_slide_list/trustworthy_slides.csv`
- `docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv`
- `data/tcga_brca/clinical_her2_laptop_balanced61_trustworthy_slides.csv`
- `data/tcga_brca/clinical_her2_laptop_balanced61_high_trust_slides.csv`

Main result for the 183 downloaded 61/61/61 cohort:

| Trust category | Slides |
|---|---:|
| High label+slide trust | 174 |
| Review before primary analysis | 9 |
| Exclude from primary analysis | 0 |

High-trust slides by group:

| HER2 group | High-trust slides | Review slides |
|---|---:|---:|
| HER2-positive | 56 | 5 |
| HER2-low | 57 | 4 |
| HER2-zero | 61 | 0 |

All 183 slides were present locally, matched expected file size, and could be opened with OpenSlide. The nine review slides are flagged because of HER2 label ambiguity or IHC/ISH discordance, not because of file integrity.

## 2026-06-02: High-Trust 174-Slide GigaTIME Run

We then processed the high-trust subset with GigaTIME using a laptop-safe 128-tile run.

Run:

```text
conda run -n gigatime-tcga python scripts/run_gigatime_tcga_brca.py \
  --slide-table docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv \
  --slide-path-column slide_local_path \
  --missing-slide-policy skip \
  --out-dir results/gigatime_tcga_brca_clinical_her2_high_trust_tile128 \
  --tile-limit 128 \
  --batch-size 16 \
  --tile-order random \
  --random-seed 2026 \
  --save-tile-csv \
  --heatmap-channels "" \
  --resume
```

The GigaTIME runner now has a `--resume` option so large runs can skip slides already present in `slide_scores.csv`.

Outputs:

- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/slide_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/tile_scores.csv`
- `docs/clinical_her2_high_trust_tile128_results.md`
- `docs/clinical_her2_high_trust_tile128_gigatime_data_cleanup.md`
- `docs/clinical_her2_high_trust_tile128_cleaned_classifier_comparison.md`

Run size:

| Group | Slides |
|---|---:|
| HER2-positive | 56 |
| HER2-low | 57 |
| HER2-zero | 61 |
| Total | 174 |

Tile predictions:

- 128 random tissue tiles per slide.
- 22,272 total tile predictions.
- All 174 slides reached the 128-tile cap.

The high-trust run strengthened the main scientific pattern. HER2-low was lower than HER2-zero for several virtual immune/myeloid/checkpoint and tissue-context channels:

| Channel | HER2-low minus HER2-zero | Mann-Whitney p | BH q |
|---|---:|---:|---:|
| CD68 | -0.00537 | 0.000371 | 0.00223 |
| CK | -0.06377 | 0.000129 | 0.00223 |
| PD-L1 | -0.01301 | 0.000302 | 0.00223 |
| PD-1 | -0.03948 | 0.000225 | 0.00223 |
| CD11c | -0.00325 | 0.000272 | 0.00223 |
| CD4 | -0.02379 | 0.000615 | 0.00263 |
| CD3 | -0.02433 | 0.000778 | 0.00292 |

Cleanup/tile-filtering result:

- The HER2-low versus HER2-zero signal persisted after cellular-tissue QC.
- The signal partially persisted after CK enrichment.
- The signal weakened in the strictest CK-enriched top 25% view, suggesting the current image-derived signal may reflect broader tissue/microenvironment context more than a purely epithelial tumor-cell phenotype.

Classifier result:

| Task | Best view/features | Balanced accuracy | Macro AUC |
|---|---|---:|---:|
| HER2-low vs HER2-zero | All sampled tissue, mean + fraction channels | 0.727 | 0.787 |
| HER2-low vs HER2-zero | QC-cellular tissue, mean channels | 0.719 | 0.741 |
| HER2-positive vs HER2-negative | CK-enriched top 25%, mean + fraction channels | 0.574 | 0.600 |
| Three-class HER2 | CK-enriched top 25%, mean + fraction channels | 0.518 | 0.689 |

Interpretation:

- The strongest current paper-ready result is not diagnostic HER2 prediction.
- The stronger result is that high-trust TCGA-BRCA H&E/GigaTIME features reproducibly associate with the HER2-low versus HER2-zero boundary.
- HER2-positive classification remains weak, so we should not present GigaTIME as a reliable HER2 diagnostic model.

## 2026-06-02: ER/PR and HER2-Detail Sensitivity Checks

After the pre-sex-filter high-trust 174-slide run, later refined to the strict 171-slide analysis set, we tested two obvious alternative explanations for the HER2-low versus HER2-zero signal:

1. The signal might be explained by hormone-receptor imbalance between HER2-low and HER2-zero cases.
2. The signal might be driven by only one HER2 IHC/ISH detail subgroup.

Script:

```text
conda run -n gigatime-tcga python scripts/analyze_high_trust_her2_sensitivity.py
```

Outputs:

- `docs/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity.md`
- `docs/assets/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity/erpr_adjusted_low_zero_q_heatmap.png`
- `docs/assets/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity/her2_detail_subgroup_boxplots.png`
- `docs/assets/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity/low_zero_erpr_adjusted_tests.csv`
- `docs/assets/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity/low_zero_erpr_stratified_tests.csv`
- `docs/assets/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity/her2_detail_subgroup_tests.csv`
- `docs/assets/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity/her2_detail_subgroup_contrasts.csv`
- `docs/assets/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity/sensitivity_summary.json`

Input counts:

| Group | Slides | ER positive | ER negative | PR positive | PR negative |
|---|---:|---:|---:|---:|---:|
| HER2-positive | 56 | 46 | 10 | 32 | 24 |
| HER2-low | 57 | 45 | 12 | 41 | 16 |
| HER2-zero | 61 | 43 | 18 | 37 | 24 |

ER/PR-adjusted result:

| Cleanup view | Unadjusted q<0.05 channels | ER/PR adjusted q<0.05 channels | ER/PR+ERBB2 adjusted q<0.05 channels |
|---|---:|---:|---:|
| All sampled tissue | 8 | 7 | 4 |
| QC cellular tissue | 7 | 4 | 3 |
| CK-enriched top 50% | 7 | 4 | 3 |
| CK-enriched top 25% | 7 | 1 | 2 |

All-sampled-tissue ER/PR-adjusted channels that remained significant:

- `CK`
- `PD-1`
- `CD4`
- `CD3`
- `PD-L1`
- `CD68`
- `CD11c`

`Ki67` did not remain significant after ER/PR adjustment.

HER2-detail subgroup result:

- HER2-low IHC `1+` cases remained lower than HER2-zero overall.
- HER2-low IHC `2+`/ISH-negative cases remained lower than HER2-zero overall.
- HER2-low overall remained lower than HER2-zero IHC `0`/ISH-negative.
- HER2-low overall remained lower than HER2-zero IHC `0`/ISH-not-evaluated.

Interpretation:

- The main all-sampled-tissue signal is not obviously explained by ER/PR imbalance.
- The signal is not obviously driven by only one HER2-low or HER2-zero detail subgroup.
- The strict CK-enriched view still weakens, so the best current biological framing remains a broader tissue-context association rather than a purely epithelial HER2 phenotype.
- This is still a sensitivity check, not causal proof and not clinical validation.

## 2026-06-02: 128-Tile Versus 256-Tile Run Agreement

We then compared the current high-trust 128-tile run against the earlier expanded 60-slide 256-tile run on overlapping slides.

Script:

```text
conda run -n gigatime-tcga python scripts/compare_gigatime_run_agreement.py
```

Outputs:

- `docs/clinical_her2_high_trust_tile128_vs_expanded20_tile256_agreement.md`
- `docs/assets/clinical_her2_high_trust_tile128_vs_expanded20_tile256/tile128_vs_tile256_channel_correlation_heatmap.png`
- `docs/assets/clinical_her2_high_trust_tile128_vs_expanded20_tile256/key_channel_tile128_vs_tile256_scatter.png`
- `docs/assets/clinical_her2_high_trust_tile128_vs_expanded20_tile256/low_zero_delta_tile128_vs_tile256.png`
- `docs/assets/clinical_her2_high_trust_tile128_vs_expanded20_tile256/run_channel_agreement.csv`
- `docs/assets/clinical_her2_high_trust_tile128_vs_expanded20_tile256/low_zero_direction_comparison.csv`
- `docs/assets/clinical_her2_high_trust_tile128_vs_expanded20_tile256/overlap_slide_scores.csv`
- `docs/assets/clinical_her2_high_trust_tile128_vs_expanded20_tile256/run_agreement_summary.json`

Overlap:

| Quantity | Count |
|---|---:|
| Expanded tile256 run slides | 60 |
| High-trust tile128 run slides | 174 |
| Overlapping slide IDs | 58 |
| Overlapping HER2-low slides | 20 |
| Overlapping HER2-zero slides | 20 |
| Overlapping HER2-positive slides | 18 |

The two missing expanded-run cases were HER2-positive review/excluded cases, so the HER2-low and HER2-zero comparison has complete overlap.

Key-channel agreement:

| Channel | Spearman rho | Reference low-zero delta | High-trust low-zero delta | Both low lower than zero |
|---|---:|---:|---:|---|
| PD-1 | 0.988 | -0.05030 | -0.04810 | yes |
| CD68 | 0.988 | -0.00699 | -0.00632 | yes |
| PD-L1 | 0.986 | -0.01447 | -0.01303 | yes |
| CD3 | 0.985 | -0.03711 | -0.03387 | yes |
| Ki67 | 0.984 | 0.00181 | 0.00198 | no |
| CK | 0.982 | -0.06973 | -0.07023 | yes |
| CD4 | 0.977 | -0.03719 | -0.03365 | yes |
| CD11c | 0.970 | -0.00508 | -0.00427 | yes |

Interpretation:

- Slide-level GigaTIME channel agreement is very high across the two run settings.
- All 8 tested key channels preserve the same HER2-low versus HER2-zero direction across runs.
- 7 of 8 tested key channels have HER2-low lower than HER2-zero in both runs.
- This supports that the HER2-low versus HER2-zero signal is not a single-run or tile-count artifact.
- The comparison is not a perfect randomized parameter experiment because the tile samples differ and the high-trust run excludes two HER2-positive review cases, but it is strong enough to support run-level robustness for the low-versus-zero boundary.

## Current Limitations to State Clearly

- The project now includes a strict high-trust 171-slide TCGA-BRCA clinical HER2 analysis, but this is still small for a clinical model and still lacks external validation.
- The first full clinical HER2 run used 64 random tissue tiles per slide; the 256-tile reruns improve sampling robustness but are still not exhaustive whole-slide analysis.
- The high-trust run used 128 random tissue tiles per slide, which is larger in cases but still not exhaustive whole-slide analysis.
- The earlier ERBB2 RNA-expression extreme comparison should not be treated as the clinical HER2 result.
- Clinical HER2 fields in TCGA are incomplete for many cases.
- TCGA clinical supplement files may contain missing, not evaluated, or inconsistent fields.
- No matched real mIF validation data is currently present in this project.
- GigaTIME predictions are research features, not clinical measurements.
- The classifier baselines should not be interpreted as diagnostic performance.
- The current data do not prove HER2 isoform state, targetability, or therapy resistance. Those are future hypotheses requiring molecular, protein-level, and treatment-response validation.
- Whole-slide sampling, tile quality, tumor purity, and batch/stain variation need stronger QC.

## Reproducibility Checklist

Current workflow scripts:

- `scripts/gdc_query_tcga_brca.py`
- `scripts/build_tcga_brca_clinical_her2_labels.py`
- `scripts/select_clinical_her2_cohort.py`
- `scripts/download_clinical_her2_cohort_slides.py`
- `scripts/select_her2_extremes.py`
- `scripts/run_gigatime_tcga_brca.py`
- `scripts/summarize_her2_gigatime.py`
- `scripts/summarize_clinical_her2_gigatime.py`
- `scripts/validate_gigatime_with_rna_signatures.py`
- `scripts/validate_gigatime_with_rna_programs.py`
- `scripts/cleanup_gigatime_tile_features.py`
- `scripts/train_her2_classifier_baseline.py`
- `scripts/train_her2_cleaned_classifier_comparison.py`
- `scripts/render_clinical_her2_visual_qc.py`
- `scripts/build_clinical_her2_findings_report.py`
- `scripts/render_he_slide_images.py`
- `scripts/render_virtual_mif_channel_images.py`
- `scripts/render_virtual_mif_composites.py`
- `scripts/build_tcga_her2_trustworthy_slide_list.py`
- `scripts/analyze_high_trust_case_drivers.py`
- `scripts/render_case_driver_visual_qc.py`
- `scripts/analyze_tissue_composition_sensitivity.py`

Current documentation:

- `README.md`
- `docs/current_pilot_run.md`
- `docs/advisor_brief.md`
- `docs/plain_language_methodology.md`
- `docs/virtual_mif_channel_outputs.md`
- `docs/paper_proposal_process_log.md`
- `docs/clinical_her2_cohort_selection.md`
- `docs/clinical_her2_gigatime_run.md`
- `docs/clinical_her2_rna_validation.md`
- `docs/clinical_her2_visual_qc.md`
- `docs/clinical_her2_tile_sampling_robustness.md`
- `docs/clinical_her2_rna_program_validation.md`
- `docs/clinical_her2_gigatime_data_cleanup.md`
- `docs/clinical_her2_classifier_baseline.md`
- `docs/clinical_her2_cleaned_classifier_comparison.md`
- `docs/tcga_her2_label_quality_assessment.md`
- `docs/clinical_her2_trustworthy_slide_list.md`
- `docs/clinical_her2_high_trust_tile128_results.md`
- `docs/clinical_her2_high_trust_tile128_case_driver_analysis.md`
- `docs/clinical_her2_high_trust_tile128_case_driver_visual_qc.md`
- `docs/clinical_her2_high_trust_tile128_tissue_composition_sensitivity.md`
- `docs/clinical_her2_high_trust_tile128_gigatime_data_cleanup.md`
- `docs/clinical_her2_high_trust_tile128_cleaned_classifier_comparison.md`
- `docs/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity.md`
- `docs/clinical_her2_high_trust_tile128_vs_expanded20_tile256_agreement.md`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

Current key result files:

- `data/tcga_brca/erbb2_expression.csv`
- `data/tcga_brca/clinical_her2_labels.csv`
- `data/tcga_brca/clinical_her2_labels_metadata.json`
- `data/tcga_brca/clinical_her2_laptop_balanced61_trustworthy_slides.csv`
- `data/tcga_brca/clinical_her2_laptop_balanced61_high_trust_slides.csv`
- `data/tcga_brca/clinical_her2_cohort_cases.csv`
- `data/tcga_brca/clinical_her2_cohort_slides_files.csv`
- `data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv`
- `data/tcga_brca/clinical_her2_cohort_slide_download_status.json`
- `data/tcga_brca/clinical_her2_cohort_summary.json`
- `data/tcga_brca/her2_extreme_cases.csv`
- `results/gigatime_tcga_brca_extremes/slide_scores.csv`
- `results/gigatime_tcga_brca_extremes/tile_scores.csv`
- `results/gigatime_tcga_brca_extremes/advisor_summary/joined_slide_her2_gigatime.csv`
- `results/gigatime_tcga_brca_extremes/advisor_summary/her2_group_channel_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2/slide_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_summary.md`
- `results/gigatime_tcga_brca_clinical_her2/rna_validation/gigatime_rna_signature_correlations.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/clinical_summary/clinical_her2_summary.md`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_validation/gigatime_rna_signature_correlations.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation/virtual_rna_program_correlations.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/gigatime_cleanup/cleaned_slide_features.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/gigatime_cleanup/cleanup_channel_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/classifier_metrics.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/classifier_baseline/classifier_crossval_predictions.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/cleaned_classifier_comparison/cleaned_classifier_best_h_e_metrics.csv`
- `results/gigatime_tcga_brca_clinical_her2_tile256/cleaned_classifier_comparison/cleaned_classifier_metrics.csv`
- `docs/assets/clinical_her2_visual_qc/clinical_her2_visual_qc_selected_cases.csv`
- `docs/assets/clinical_her2_visual_qc_tile256/clinical_her2_visual_qc_selected_cases.csv`
- `docs/assets/clinical_her2_findings/clinical_her2_group_mean_heatmap.png`
- `docs/assets/clinical_her2_tile256/clinical_her2_group_mean_heatmap.png`
- `docs/assets/clinical_her2_rna_program_validation/virtual_rna_program_correlation_heatmap.png`
- `docs/assets/clinical_her2_gigatime_cleanup/cleanup_key_channel_heatmap.png`
- `docs/assets/clinical_her2_classifier_baseline/classifier_balanced_accuracy.png`
- `docs/assets/clinical_her2_cleaned_classifier/cleaned_classifier_best_by_view.png`

## 2026-06-02: Guardia et al. Paper Correction And Strict Trustworthy Slide List

The advisor clarified that the intended Galante-lab reference was Guardia et al., Genome Research 2025, PMID `40664477`: "Alternative splicing generates HER2 isoform diversity underlying antibody-drug conjugate resistance in breast cancer."

We corrected the project documentation away from the earlier wrong-reference framing and used only the HER2 isoform paper's relevant TCGA sample-selection principles. Important boundary: Guardia et al. is an RNA-seq isoform/ADC-resistance paper, not an H&E slide-QC or virtual mIF validation paper. Our HER2 label, file-integrity, OpenSlide, and tile-quality checks are project-specific.

- Use TCGA-BRCA primary tumor material.
- Exclude male TCGA-BRCA samples from the strict primary analysis.
- Keep HER2 groups tied to IHC/FISH/ISH clinical annotation.
- Stratify or adjust by hormone-receptor context where possible.
- Treat transcript/isoform evidence as validation or hypothesis support, not as a substitute for clinical HER2 IHC/ISH labels.

This changed the primary trustworthy slide list:

| Trust category | Slides |
|---|---:|
| Strict high label+slide trust | 171 |
| Review before primary analysis | 9 |
| Exclude from primary analysis | 3 |

The three newly excluded slides are male HER2-positive TCGA-BRCA cases:

| Case | HER2 group | Reason |
|---|---|---|
| TCGA-A1-A0SM | HER2-positive | Male patient |
| TCGA-E2-A14W | HER2-positive | Male patient |
| TCGA-EW-A1PD | HER2-positive | Male patient |

The strict primary analysis now uses:

| HER2 group | Slides |
|---|---:|
| HER2-positive | 53 |
| HER2-low | 57 |
| HER2-zero | 61 |

The raw GigaTIME output directory still contains the earlier 174-slide inference run. We did not need to rerun image inference; instead, we filtered the existing GigaTIME outputs to the 171 strict trustworthy slides and regenerated the clinical summary, cleanup analysis, ER/PR and HER2-detail sensitivity analysis, overlap/run-agreement check, and cleaned classifier comparison.

Updated strongest result:

- HER2-low remains lower than HER2-zero for `CD68`, `CK`, `PD-L1`, `PD-1`, `CD11c`, `CD4`, and `CD3`.
- All-sampled-tissue HER2-low versus HER2-zero BH q values are about `0.0016-0.0020` for those key channels.
- The ER/PR-adjusted sensitivity result remains: 7 of 8 all-sampled-tissue key channels stay significant after ER/PR adjustment.
- The overlap check now uses 56 matched slides between the 60-slide 256-tile run and the strict 171-slide 128-tile analysis. All 8 tested key channels still keep the same HER2-low versus HER2-zero direction, and 7 of 8 keep HER2-low lower than HER2-zero.
- The cleaned HER2-low versus HER2-zero classifier remains exploratory but nontrivial: all-sampled-tissue balanced accuracy `0.727`, macro AUC `0.787`; QC-cellular balanced accuracy `0.719`, macro AUC `0.741`.

Updated key files:

- `docs/clinical_her2_trustworthy_slide_list.md`
- `docs/assets/clinical_her2_trustworthy_slide_list/trustworthy_slides.csv`
- `docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv`
- `docs/clinical_her2_high_trust_tile128_results.md`
- `docs/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity.md`
- `docs/clinical_her2_high_trust_tile128_vs_expanded20_tile256_agreement.md`
- `docs/clinical_her2_high_trust_tile128_cleaned_classifier_comparison.md`

## 2026-06-02: Case-Level Driver Analysis For The HER2-Low Versus HER2-Zero Signal

After regenerating the strict 171-slide analysis, we added a case-level driver analysis to test whether the most interesting result is broad or dominated by a few unusual slides.

Method:

- Use the strict 171-slide high-trust GigaTIME output.
- Restrict the driver score to HER2-low and HER2-zero slides.
- Within each cleanup view, select only the virtual channels that significantly separate HER2-low from HER2-zero.
- Standardize those channels and orient them so higher values are HER2-zero-like and lower values are HER2-low-like.
- Check whether each slide keeps the expected direction across all sampled tissue, QC-cellular tissue, CK-enriched top 50%, and CK-enriched top 25%.
- Join the driver score to the best low-versus-zero classifier predictions to create a manual review list.

Current result:

| Case-level check | Result |
|---|---:|
| HER2-low/HER2-zero slides scored | 118 |
| Slides matching expected direction in at least 3 of 4 cleanup views | 71 |
| Slides matching expected direction in all 4 cleanup views | 63 |
| Slides with opposite profile in at least 2 cleanup views | 47 |
| Cases misclassified by the best low-vs-zero classifier in at least 2 cleanup views | 37 |

All-sampled-tissue case score:

| Group | N | Mean zero-like score | Median zero-like score |
|---|---:|---:|---:|
| HER2-low | 57 | -0.219 | -0.488 |
| HER2-zero | 61 | 0.205 | 0.123 |

Interpretation:

- The HER2-low versus HER2-zero result is not only a group-average table; it now has a slide-level score and a stability check.
- The result is promising but imperfect. A substantial subset of slides behaves opposite to the expected direction across multiple cleanup views.
- Those opposite-profile and classifier-error cases should be reviewed in H&E and virtual mIF-like overlays before making stronger biological claims.
- This analysis strengthens the paper path because it gives a concrete pathologist/QC shortlist rather than hiding model errors.

Updated key files:

- `scripts/analyze_high_trust_case_drivers.py`
- `docs/clinical_her2_high_trust_tile128_case_driver_analysis.md`
- `docs/assets/clinical_her2_high_trust_tile128_case_drivers/`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/case_driver_analysis/`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## 2026-06-02: Case-Driver Visual QC And Tissue-Composition Caveat

We then rendered a small H&E plus virtual mIF visual QC set from the case-driver shortlist.

Cases rendered:

- 2 stable label-consistent HER2-low cases.
- 2 stable label-consistent HER2-zero cases.
- 4 opposite-profile/manual-review cases.
- 4 selected tiles per case.

The visual panels were generated by rerunning GigaTIME on selected tile coordinates so that the report could show actual spatial virtual mIF-like maps, not only tile-level CSV means.

Selected-tile summary:

| Review category | Group | Tiles | Mean tissue | Mean zero-like tile score | Mean CK | Mean CD68 | Mean PD-L1 | Mean CD11c |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Label-consistent HER2-low | HER2-low | 8 | 0.987 | -0.605 | 0.0003 | 0.0002 | 0.0006 | 0.0002 |
| Label-consistent HER2-zero | HER2-zero | 8 | 0.955 | 1.675 | 0.0487 | 0.1300 | 0.3153 | 0.1173 |
| Opposite-profile manual review | HER2-low | 8 | 0.923 | 2.800 | 0.0684 | 0.1786 | 0.3879 | 0.1483 |
| Opposite-profile manual review | HER2-zero | 8 | 0.986 | -0.618 | 0.0012 | 0.0003 | 0.0011 | 0.0003 |

Interpretation:

- This is a serious tissue-composition caveat.
- The low-like selected tiles are tissue-containing, but they can be very low for virtual CK, CD68, PD-L1, and CD11c.
- Visual inspection suggests some low-like tiles may be stromal/collagen-rich rather than clearly tumor-rich.
- The current HER2-low versus HER2-zero signal may therefore partly reflect broader tissue composition, not a purely tumor-cell HER2 phenotype.
- This does not make the analysis useless; it makes the next step clear: pathologist review and tumor-rich tile restriction.

Updated key files:

- `scripts/render_case_driver_visual_qc.py`
- `docs/clinical_her2_high_trust_tile128_case_driver_visual_qc.md`
- `docs/assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/case_driver_visual_qc/`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## 2026-06-02: Quantified Tissue-Composition Sensitivity

The case-driver visual QC raised a concern that the HER2-low versus HER2-zero signal could be driven by stromal/collagen-rich or low-marker tissue context. We quantified that concern across all HER2-low and HER2-zero slides.

Definitions:

- Marker burden = mean virtual `CK`, `CD68`, `PD-L1`, `CD11c`, `CD3`, `CD4`, `CD20`, and `Ki67` per tile.
- Low-marker tile = marker burden in the bottom quartile across strict high-trust tiles.
- Absolute CK-high tile = QC-cellular tile with virtual CK in the top quartile across QC tiles.

Main result:

| Metric | HER2-low mean | HER2-zero mean | Low-zero delta | BH q |
|---|---:|---:|---:|---:|
| Fraction low-marker tiles | 0.349 | 0.180 | 0.169 | 0.000265 |
| Fraction very-low-marker tiles | 0.150 | 0.060 | 0.089 | 0.000265 |
| Mean marker burden | 0.0479 | 0.0669 | -0.0190 | 0.000288 |
| Mean CK | 0.167 | 0.231 | -0.0638 | 0.000288 |
| Fraction high-marker tiles | 0.168 | 0.287 | -0.119 | 0.00139 |
| Fraction high-CK QC tiles | 0.145 | 0.218 | -0.0728 | 0.00192 |

The case-driver score also tracks tissue composition strongly:

- Spearman rho with mean marker burden: 0.980.
- Spearman rho with low-marker tile fraction: -0.782.

The most important covariate result:

- After adjusting for low-marker tile fraction, most HER2-low versus HER2-zero channel effects collapse.
- CD68, PD-L1, PD-1, CD11c, CD4, CD3, and CK all lose FDR significance after this adjustment.

Interpretation:

- The strongest honest result is now a GigaTIME-derived HER2-low versus HER2-zero tissue-context association.
- We should not present this as tumor-cell HER2 biology, HER2 isoform detection, or diagnostic HER2 prediction.
- This is still meaningful: it may indicate that HER2-low and HER2-zero TCGA tumors differ in broader tissue composition or microenvironment patterns visible through H&E-derived virtual channels.
- The next required validation is tumor-rich/pathologist-approved tile restriction.

Updated key files:

- `scripts/analyze_tissue_composition_sensitivity.py`
- `docs/clinical_her2_high_trust_tile128_tissue_composition_sensitivity.md`
- `docs/assets/clinical_her2_high_trust_tile128_tissue_composition/`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/tissue_composition_sensitivity/`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## Next Immediate Step

The next step is not another download. The 30-slide clinical HER2 pilot, 256-tile robustness check, broader RNA-program validation, first classifier baseline, pre-classifier GigaTIME cleanup, cleaned-view classifier comparison, expanded 20/20/20 run, strict 171-slide high-trust analysis, and ER/PR/subgroup sensitivity checks are complete. The next scientific step is trustworthiness review of the cases driving model behavior:

- Ask an advisor/pathologist to review whether the H&E regions driving high virtual CD68, PD-L1, CD11c, CK, CD3, and CD4 are biologically plausible.
- Inspect the highest-priority opposite-profile and classifier-error cases from `clinical_her2_high_trust_tile128_case_driver_analysis.md`.
- Use `clinical_her2_high_trust_tile128_case_driver_visual_qc.md` to decide whether low-like driver tiles are tumor-rich or mostly stromal/collagen-rich.
- Use `clinical_her2_high_trust_tile128_tissue_composition_sensitivity.md` to frame the result as tissue-context sensitive.
- If feasible, rerun the key HER2-low versus HER2-zero statistics using tumor-rich/pathologist-approved tile subsets.
- Add tile distribution features and, if available, GigaTIME/pathology embeddings.
- Adjust for tumor purity or immune deconvolution if available.
- Check whether endothelial/stromal/tissue-composition differences explain part of the virtual signal.
- Search for an external dataset with paired H&E and real mIF for direct validation.

## 2026-06-02: Tumor-Rich Proxy Sensitivity

After the tissue-composition caveat, we tested stricter virtual tumor-rich proxy views using GigaTIME-derived virtual DAPI, CK, and marker-burden features. These are not pathologist tumor annotations, but they help answer whether the HER2-low versus HER2-zero signal completely disappears when we push toward CK-rich/epithelial-like tiles.

Proxy views tested:

- QC cellular tissue.
- CK top 25% within slide.
- Top 16 CK tiles per slide.
- Top 8 CK tiles per slide.
- Top 16 CK tiles after removing low-marker tiles.
- Absolute CK-high QC tiles, using the global QC-tile CK top quartile.

Main low-versus-zero results:

| Proxy view | Classifier balanced accuracy | Macro AUC | q<0.05 low-zero channels |
|---|---:|---:|---|
| QC cellular tissue | 0.719 | 0.741 | CD68, PD-L1, PD-1, CD11c, CD4, CD3, CD20, CK |
| CK top 25% within slide | 0.708 | 0.744 | CD68, PD-L1, PD-1, CD4, CD3, CK |
| Top 16 CK tiles per slide | 0.711 | 0.766 | CD68, PD-1, CK |
| Top 8 CK tiles per slide | 0.727 | 0.755 | none |
| Top 16 CK, non-low-marker | 0.708 | 0.761 | none |
| Absolute CK-high QC tiles | 0.761 | 0.782 | PD-1, CD11c, CD4, CD3 |

Interpretation:

- The univariate marker-channel tests weaken under the strictest fixed-count CK-rich proxy filters.
- The multichannel low-versus-zero classifier remains above chance across proxy views.
- This suggests the signal is not only blank/background artifact and may contain a multichannel CK-rich/tumor-proxy pattern.
- However, because the proxy itself is GigaTIME-derived and not pathologist-confirmed tumor annotation, this still does not prove tumor-cell HER2 biology.
- The most honest next step is pathologist-approved tumor-rich tile selection or a validated tumor segmentation model.

Updated key files:

- `scripts/analyze_tumor_proxy_sensitivity.py`
- `docs/clinical_her2_high_trust_tile128_tumor_proxy_sensitivity.md`
- `docs/assets/clinical_her2_high_trust_tile128_tumor_proxy_sensitivity/`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/tumor_proxy_sensitivity/`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## 2026-06-02: Classifier Permutation Sanity Check

After the tumor-rich proxy analysis, we added a shuffled-label sanity check for the selected HER2-low versus HER2-zero classifiers.

Question:

- If the HER2-low/HER2-zero labels are randomly shuffled, do the same classifiers still perform well?
- If they do, the classifier result could be mostly modeling artifact.
- If they do not, the classifier has at least some real label-associated structure.

Method:

- Task: HER2-low versus HER2-zero.
- Model: regularized logistic regression using the same selected GigaTIME/H&E feature set per view.
- Evaluation: repeated stratified 5-fold cross-validation with 3 repeats.
- Null model: 100 shuffled-label permutations per view using the same folds and feature columns.

Main result:

| View | LOOCV balanced accuracy | Repeated-CV balanced accuracy | Null mean | Null 95% | Empirical p | BH q | Repeated-CV AUC |
|---|---:|---:|---:|---:|---:|---:|---:|
| QC cellular tissue | 0.719 | 0.705 | 0.484 | 0.566 | 0.0099 | 0.0099 | 0.744 |
| CK top 25% within slide | 0.708 | 0.693 | 0.485 | 0.557 | 0.0099 | 0.0099 | 0.731 |
| Top 16 CK tiles per slide | 0.711 | 0.705 | 0.485 | 0.561 | 0.0099 | 0.0099 | 0.763 |
| Top 8 CK tiles per slide | 0.727 | 0.716 | 0.488 | 0.586 | 0.0099 | 0.0099 | 0.741 |
| Top 16 CK, non-low-marker | 0.708 | 0.710 | 0.482 | 0.567 | 0.0099 | 0.0099 | 0.767 |
| Absolute CK-high QC tiles | 0.761 | 0.729 | 0.484 | 0.575 | 0.0099 | 0.0099 | 0.764 |

Interpretation:

- The classifier result is not obviously random: every selected view beats its shuffled-label null distribution.
- This strengthens the trustworthiness story for the HER2-low versus HER2-zero classifier.
- It is still a post-hoc sanity check, not a fully nested model-selection permutation test.
- It does not solve the tissue-composition caveat, validate real mIF, prove tumor-cell HER2 biology, or detect HER2 isoforms.

Updated key files:

- `scripts/analyze_classifier_permutation_sanity.py`
- `docs/clinical_her2_high_trust_tile128_classifier_permutation_sanity.md`
- `docs/assets/clinical_her2_high_trust_tile128_classifier_permutation/`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/classifier_permutation_sanity/`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## 2026-06-02: Nested Classifier Model-Selection Check

The previous permutation check tested selected low-versus-zero feature sets. To reduce the concern that we were reporting a winner chosen after looking at the whole dataset, we added a stricter nested model-selection check.

Question:

- Does the HER2-low versus HER2-zero classifier still work if feature-set selection happens inside cross-validation?
- Does it still beat a shuffled-label null when feature-set selection is also repeated inside each shuffled-label run?

Method:

- Task: HER2-low versus HER2-zero.
- Model: regularized logistic regression.
- Outer evaluation: repeated stratified 5-fold cross-validation with 3 repeats.
- Inner feature-set selection: stratified 4-fold cross-validation inside each outer training fold only.
- Candidate feature sets: GigaTIME mean channels, mean+fraction channels, interpretable marker means, interpretable distribution features, and virtual programs.
- Null: 30 shuffled-label runs per view, each with nested feature-set selection repeated.

Main result:

| View | Nested balanced accuracy | Nested AUC | Null mean | Null 95% | Empirical p | BH q | Most selected feature set |
|---|---:|---:|---:|---:|---:|---:|---|
| QC cellular tissue | 0.674 | 0.717 | 0.498 | 0.577 | 0.0323 | 0.0323 | Mean channels |
| CK top 25% within slide | 0.672 | 0.706 | 0.500 | 0.578 | 0.0323 | 0.0323 | Mean channels |
| Top 16 CK tiles per slide | 0.706 | 0.731 | 0.506 | 0.566 | 0.0323 | 0.0323 | Mean channels |
| Top 8 CK tiles per slide | 0.721 | 0.739 | 0.514 | 0.565 | 0.0323 | 0.0323 | Mean channels |
| Top 16 CK, non-low-marker | 0.676 | 0.717 | 0.492 | 0.548 | 0.0323 | 0.0323 | Mean channels |
| Absolute CK-high QC tiles | 0.717 | 0.766 | 0.505 | 0.553 | 0.0323 | 0.0323 | Mean channels |

Interpretation:

- The low-versus-zero classifier signal survives nested feature-set selection.
- This reduces the concern that the signal is only a feature-set selection artifact.
- The most commonly selected feature set is usually GigaTIME mean channels, suggesting the useful information is in the core virtual marker means rather than a fragile exotic feature set.
- This is still internal validation. It does not prove clinical diagnosis, real mIF validity, HER2 isoform biology, or therapy-response prediction.

Updated key files:

- `scripts/analyze_nested_classifier_model_selection.py`
- `docs/clinical_her2_high_trust_tile128_nested_classifier_model_selection.md`
- `docs/assets/clinical_her2_high_trust_tile128_nested_classifier/`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/nested_classifier_model_selection/`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## 2026-06-02: Clinical, Source-Site, And Slide-Size Confounder Sensitivity

After the nested classifier check, we asked whether the low-versus-zero signal could be explained by non-image metadata.

Question:

- Are HER2-low and HER2-zero slides balanced for ordinary clinical covariates, TCGA source site, and slide size?
- Can non-image covariates classify HER2-low versus HER2-zero?
- Do channel-level low-versus-zero effects survive after clinical/source-site/slide-size adjustment?

Main covariate imbalance:

| Covariate | HER2-low mean | HER2-zero mean | Low-zero delta | p |
|---|---:|---:|---:|---:|
| Slide file size MB | 100.486 | 278.961 | -178.475 | 1.45e-15 |
| Slide width | 61509.123 | 102679.066 | -41169.943 | 2.58e-12 |
| Slide height | 23300.702 | 31759.115 | -8458.413 | 6.49e-07 |
| Mean marker burden | 0.056 | 0.076 | -0.019 | 5.79e-04 |
| Mean virtual DAPI | 0.306 | 0.384 | -0.078 | 1.63e-04 |
| Mean virtual CK | 0.191 | 0.245 | -0.054 | 0.0017 |

Categorical imbalance:

- Histology group chi-square p = 0.0413.
- TCGA source-site group chi-square p = 7.97e-10.
- ER, PR, and broad stage group were not the major imbalance signals in this check.

Classifier result in the top 8 CK proxy view:

| Feature set | Balanced accuracy | AUC |
|---|---:|---:|
| Clinical covariates | 0.536 | 0.512 |
| Slide-size covariates | 0.879 | 0.921 |
| Source-site covariates | 0.878 | 0.925 |
| Source-site + slide-size covariates | 0.897 | 0.965 |
| GigaTIME mean channels | 0.745 | 0.751 |
| GigaTIME + clinical + site/slide | 0.890 | 0.952 |

Channel-level adjustment:

- In QC-cellular tissue, several key channels were significant before broad site/slide adjustment.
- After clinical plus source-site/slide-size adjustment, no key QC-cellular channel remained significant at q < 0.05.
- In the stricter top 8 CK and absolute CK-high views, no key channel was significant at q < 0.05 even before full site/slide adjustment.

Interpretation:

- This is the strongest caveat so far.
- The current TCGA HER2-low/HER2-zero classifier signal is not safe to present as independent HER2 biology.
- Slide size and source site classify HER2-low versus HER2-zero better than GigaTIME image features, suggesting a major cohort-construction/acquisition confounding risk.
- The result is still useful scientifically because it identifies a concrete failure mode and next step: source-site/slide-size matched sensitivity analysis.

Updated key files:

- `scripts/analyze_clinical_covariate_sensitivity.py`
- `docs/clinical_her2_high_trust_tile128_clinical_covariate_sensitivity.md`
- `docs/assets/clinical_her2_high_trust_tile128_clinical_covariates/`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/clinical_covariate_sensitivity/`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## 2026-06-02: Matched HER2-Low Versus HER2-Zero Sensitivity

After identifying source-site and slide-size confounding, we ran a matched sensitivity analysis.

Question:

- Does the HER2-low versus HER2-zero GigaTIME signal survive when slides are matched by source site or slide size?
- Do non-image source-site/slide-size baselines still perform well after matching?
- Do paired HER2-low minus HER2-zero channel tests remain significant?

Matched subsets:

| Matched subset | Pairs | Same-source-site pairs | Median abs log-size diff | Median abs MB diff |
|---|---:|---:|---:|---:|
| Exact source-site, nearest size | 12 | 12 | 0.150 | 17.4 |
| Slide-size matched, caliper 0.25 | 14 | 2 | 0.017 | 2.0 |
| Slide-size matched, caliper 0.50 | 20 | 2 | 0.031 | 3.1 |

Top 8 CK proxy view classifier result:

| Matched subset | Slide-size BA | Source-site BA | Site+size BA | GigaTIME BA | GigaTIME AUC |
|---|---:|---:|---:|---:|---:|
| Exact source-site, nearest size | 0.667 | 0.500 | 0.625 | 0.708 | 0.750 |
| Slide-size matched, caliper 0.25 | 0.607 | 0.821 | 0.786 | 0.679 | 0.694 |
| Slide-size matched, caliper 0.50 | 0.650 | 0.750 | 0.850 | 0.675 | 0.623 |

Paired channel result:

- No paired top 8 CK proxy channel test reached BH q < 0.05.
- Exact source-site subset strongest paired result: `CD20`, q = 0.347.
- Strict slide-size subset strongest paired result: `CK`, q = 0.706.
- Wider slide-size subset strongest paired result: `CK`, q = 0.138.

Interpretation:

- Matching keeps a modest GigaTIME signal, especially in exact source-site pairs.
- Matching does not eliminate the confounder concern.
- In larger matched subsets, non-image source-site/slide-size baselines remain competitive or stronger than GigaTIME.
- The safest conclusion is that GigaTIME may contain a HER2-low/HER2-zero tissue-context signal, but TCGA alone is not clean enough to claim independent HER2 biology.

Updated key files:

- `scripts/analyze_matched_low_zero_sensitivity.py`
- `docs/clinical_her2_high_trust_tile128_matched_low_zero_sensitivity.md`
- `docs/assets/clinical_her2_high_trust_tile128_matched_low_zero/`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/matched_low_zero_sensitivity/`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## 2026-06-02: Source-Site Held-Out Generalization

After matching, we ran a stricter classifier robustness check: leave one TCGA source site out during validation.

Question:

- Does the HER2-low versus HER2-zero classifier generalize to held-out source sites?
- Does GigaTIME performance drop compared with ordinary repeated stratified cross-validation?
- Do slide-size covariates remain predictive even when source sites are held out?

Main source-site imbalance:

| TSS | HER2-low | HER2-zero | Cases | Both classes |
|---|---:|---:|---:|---|
| AO | 2 | 17 | 19 | yes |
| A2 | 7 | 12 | 19 | yes |
| BH | 0 | 12 | 12 | no |
| A8 | 2 | 8 | 10 | yes |
| A7 | 10 | 0 | 10 | no |
| AN | 0 | 9 | 9 | no |

Top 8 CK proxy view:

| Feature set | Repeated CV BA | Leave-source-site-out BA | Leave-source-site-out AUC |
|---|---:|---:|---:|
| Slide-size covariates | 0.879 | 0.882 | 0.915 |
| Tissue/QC covariates | 0.581 | 0.507 | 0.478 |
| GigaTIME mean channels | 0.745 | 0.669 | 0.679 |
| GigaTIME + slide-size | 0.857 | 0.837 | 0.894 |
| GigaTIME + tissue/QC | 0.734 | 0.668 | 0.683 |

GigaTIME mean-channel performance dropped under source-site holdout across every tested feature view:

- QC cellular tissue: 0.705 to 0.617.
- CK top 25% within slide: 0.675 to 0.614.
- Top 16 CK tiles per slide: 0.755 to 0.601.
- Top 8 CK tiles per slide: 0.745 to 0.669.
- Top 16 CK, non-low-marker: 0.702 to 0.597.
- Absolute CK-high QC tiles: 0.704 to 0.643.

Interpretation:

- This is a strong classifier caveat.
- GigaTIME retains some above-chance low-vs-zero signal, but it is weaker when source sites are held out.
- Slide-size covariates remain extremely strong under source-site holdout, so the current TCGA signal is still dominated by technical/acquisition structure.
- The classifier should be presented as hypothesis-generating internal evidence, not source-independent HER2 biology.

Updated key files:

- `scripts/analyze_source_site_generalization.py`
- `docs/clinical_her2_high_trust_tile128_source_site_generalization.md`
- `docs/assets/clinical_her2_high_trust_tile128_source_site_generalization/`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/source_site_generalization/`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## 2026-06-02: HER2 Isoform Validation Feasibility Audit

After correcting the advisor paper to Guardia et al., Genome Research 2025, PMID 40664477, we checked whether the current local RNA files can directly validate HER2 isoform biology.

Question:

- Do our local TCGA RNA files contain transcript-level isoform quantification?
- Do we have RNA-seq reads, BAM files, or junction-count outputs that could reproduce kallisto/SUPPA2/rMATS-style analysis?
- Can we test whether GigaTIME features associate with HER2 isoform states right now?

Main local audit:

| Item | Count / value |
|---|---:|
| Strict high-trust cases | 171 |
| Local STAR gene-count cases | 110 |
| High-trust cases with local STAR gene counts | 56 |
| Low/zero high-trust cases with local STAR gene counts | 40 |
| Local BAM files under `data/tcga_brca` | 0 |
| Local FASTQ files under `data/tcga_brca` | 0 |
| Local junction files under `data/tcga_brca` | 0 |
| Local isoform files under `data/tcga_brca` | 0 |
| Expression file has `transcript_id` column | False |
| Expression file has junction columns | False |

Interpretation:

- The current local RNA files support gene-level ERBB2 expression and broad RNA-program context.
- They do not support direct HER2 isoform quantification.
- We cannot compute kallisto transcript TPM, SUPPA2 PSI, rMATS junction confirmation, p95/Delta16 isoform states, or antibody-binding-domain loss from the current local files.
- To test the isoform hypothesis directly, we need sample-level HER2 isoform labels from the Guardia/Galante workflow or appropriate RNA-seq read/junction data.
- The correct language remains: image AI predicts or associates with HER2 isoform/state hypotheses; it does not detect HER2 isoforms from the current data.

Updated key files:

- `scripts/audit_her2_isoform_validation_feasibility.py`
- `docs/her2_isoform_validation_feasibility.md`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/her2_isoform_validation_feasibility/`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## 2026-06-02: Expanded Local ERBB2 Gene-Level Validation

We extracted ERBB2 gene-level TPM from all local GDC STAR augmented gene-count files and joined it to the strict high-trust GigaTIME/HER2 cohort.

Question:

- Do the clinical HER2 labels look plausible against gene-level ERBB2 RNA?
- Does ERBB2 RNA explain the HER2-low versus HER2-zero GigaTIME signal?
- Can this serve as HER2 isoform validation?

Main local audit:

| Item | Count / value |
|---|---:|
| Local STAR ERBB2 cases | 110 |
| Strict high-trust cases with local ERBB2 | 56 |
| HER2-positive high-trust cases with local ERBB2 | 16 |
| HER2-low/HER2-zero high-trust cases with local ERBB2 | 40 |

Clinical HER2 group ERBB2 medians:

| Group | N | Median ERBB2 TPM |
|---|---:|---:|
| HER2-positive | 16 | 778.7 |
| HER2-low | 20 | 83.4 |
| HER2-zero | 20 | 62.7 |

Classifier-style reference checks using ERBB2 TPM alone:

| Task | N | AUC | Best-threshold balanced accuracy |
|---|---:|---:|---:|
| HER2-positive vs non-positive | 56 | 0.905 | 0.831 |
| HER2-low vs HER2-zero | 40 | 0.605 | 0.625 |

Low/zero pairwise ERBB2 test:

- HER2-low median TPM 83.4 versus HER2-zero median TPM 62.7.
- AUC 0.605.
- Mann-Whitney p/q 0.262/0.262.

Interpretation:

- Gene-level ERBB2 RNA strongly validates broad HER2-positive status as a sanity check.
- Gene-level ERBB2 RNA weakly separates HER2-low from HER2-zero.
- Therefore, the GigaTIME low/zero signal is not simply a strong ERBB2 expression split.
- This still does not validate Guardia-style HER2 isoforms because the local files are gene-level STAR count/TPM files, not transcript-level isoform, PSI, or junction files.
- After adjusting low-vs-zero GigaTIME channel tests for log ERBB2 TPM in the small RNA-overlap subset, 15 channel/view effects remain BH q < 0.05, but this does not remove the stronger source-site, slide-size, and tissue-composition caveats.

Updated key files:

- `scripts/analyze_local_erbb2_expression_validation.py`
- `docs/clinical_her2_high_trust_tile128_local_erbb2_validation.md`
- `docs/assets/clinical_her2_high_trust_tile128_local_erbb2_validation/`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/local_erbb2_expression_validation/`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## 2026-06-02: Within-Source-Site HER2-Low Versus HER2-Zero Sensitivity

After source-site and slide-size covariate analyses showed major confounding risk, we ran a stricter sensitivity check restricted to TCGA source sites that contain both HER2-low and HER2-zero cases.

Question:

- Does any GigaTIME low/zero signal remain when we avoid comparing source sites that contain only one HER2 group?
- How many cases are actually available for this within-source-site comparison?
- Does this rescue the classifier as source-independent HER2 biology?

Mixed source sites:

| TSS | HER2-low | HER2-zero | Cases |
|---|---:|---:|---:|
| A2 | 7 | 12 | 19 |
| AO | 2 | 17 | 19 |
| A8 | 2 | 8 | 10 |
| A1 | 1 | 2 | 3 |

Main result:

- Only 4 source sites qualify.
- The mixed-source-site subset has 51 cases total: 12 HER2-low and 39 HER2-zero.
- Site-fixed channel models retain 7 channel/view effects with BH q < 0.05.
- In the top 8 CK proxy view, GigaTIME all mean channels reach 0.667 repeated-CV balanced accuracy and 0.628 leave-mixed-source-site-out balanced accuracy.
- GigaTIME key mean channels weaken strongly: 0.505 repeated-CV balanced accuracy and 0.490 leave-mixed-source-site-out balanced accuracy.
- Source-site one-hot in repeated CV is not useful in this mixed-site subset: balanced accuracy 0.483.

Interpretation:

- This partially supports continued investigation because some all-channel GigaTIME signal remains inside mixed source sites.
- It does not solve the confounding problem. The subset is very small, HER2-low is underrepresented, and classifier specificity is low.
- The result should be described as a stress test that keeps the hypothesis alive, not proof of source-independent HER2 biology.

Updated key files:

- `scripts/analyze_within_source_site_low_zero.py`
- `docs/clinical_her2_high_trust_tile128_within_source_site_low_zero.md`
- `docs/assets/clinical_her2_high_trust_tile128_within_source_site/`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/within_source_site_low_zero/`
- `notebooks/clinical_her2_findings_simple.ipynb`
- `notebooks/clinical_her2_findings_simple.html`

## 2026-06-04: Generic-Embedding Confound Controls (H-Optimus-0 + Virchow2), Repo Cleanup, And External-Validation Scouting

This session had one scientific goal: stop adding GigaTIME-internal sensitivity analyses to the HER2-low versus HER2-zero result and instead run the single most decisive control, then decide where the project goes next. It is recorded here in full so future work understands what was explored and why.

### Framing: what the prior evidence already implied

Coming into the session, the accumulated sensitivity analyses (clinical/source-site covariates, matched subsets, leave-source-site-out, within-source-site, tissue composition) all pointed the same way: the HER2-low versus HER2-zero signal is real inside TCGA but heavily confounded by slide-size and TCGA source-site acquisition structure. Key prior numbers: HER2-zero slides are ~2.8x larger files than HER2-low (279 vs 100 MB, p ~ 1.5e-15); slide-size-only and source-site-only covariates classify low/zero better than GigaTIME; GigaTIME drops under source-site holdout; adjusting for low-marker tile fraction collapses the channel effects.

The open question those analyses did not resolve: is the GigaTIME "virtual immune/myeloid/checkpoint" interpretation actually required to explain the signal, or would any generic morphology representation do the same? That is what a generic-embedding control answers directly, and it had been listed as "planned" but never run.

### Repo housekeeping

We removed an unrelated GATK germline variant-calling scaffold left over from the repo's initial template: `data/reference/` (hg38 FASTA + BWA/GATK index, ~8.3 GB), `data/known_sites/` (dbSNP138/Mills/known-indels VCFs, ~10 GB), and `data/raw/NA12878.bam` (Genome-in-a-Bottle germline sample). These were never git-tracked and unreferenced by the pathology workflow; removal reclaimed ~19 GB locally (51 GB -> 33 GB on disk). The README data-hygiene note was updated to record this. Commit `0885de7`. (A separate in-flight DeepSpot one-tile smoke update was also committed, `42ba8a3`.)

### Generic-embedding control 1: H-Optimus-0

Method: extracted `bioptimus/H-optimus-0` embeddings (1536-d, mean-pooled over 128 random tissue tiles per slide) for the same 171 strict high-trust slides used in the GigaTIME primary run, via `scripts/run_hoptimus_tcga_brca.py` (the existing runner already targeted the high-trust list). Then `scripts/analyze_hoptimus_embedding_control.py` classified HER2-low versus HER2-zero (118 slides: 57 low, 61 zero) on identical cross-validation folds as the GigaTIME analyses, comparing slide-size, source-site, GigaTIME-channel, and embedding feature sets, with PCA fit inside each training fold to avoid leakage.

The control asks three questions: (1) does a generic embedding with no immune interpretation separate low/zero? (2) does it collapse under leave-source-site-out the way GigaTIME does? (3) does the portable slide-size baseline still beat it?

Result, all three "yes":

| Feature set | Repeated-CV balanced accuracy | Leave-source-site-out balanced accuracy |
|---|---:|---:|
| Slide-size covariates | 0.888 | 0.882 |
| Source-site covariates | 0.873 | 0.500 |
| GigaTIME mean channels | 0.710 | 0.617 |
| H-Optimus-0 embedding | 0.726 | 0.586 |

H-Optimus-0 separates low/zero at 0.726 (beats shuffled-label null, mean 0.488, empirical p = 0.005; stable across PCA 10/20/30 components = 0.724/0.726/0.705), slightly exceeding GigaTIME's 0.710. It collapses under source-site holdout (0.726 -> 0.586) exactly as GigaTIME does (0.710 -> 0.617), while slide-size stays portable at 0.882. The source-site 0.500 under holdout is the expected degenerate case (one-hot site identity cannot predict an unseen held-out site).

Conclusion: the GigaTIME-specific virtual-immune framing is not required. A generic morphology embedding reproduces the separation and the same source-site collapse. Committed `8d21b02`, wired into RUN_REGISTRY, the high-trust results doc, the advisor brief, and the model-experiments/script maps.

### Generic-embedding control 2: Virchow2 (independent replication)

To make the control hard to dismiss, we ran a second, architecturally distinct foundation model. Built `scripts/run_virchow2_tcga_brca.py` (cohort runner reusing the register-token-aware embedding recipe from `scripts/run_virchow2_one_slide_smoke.py`: Virchow2 returns 261 tokens = 1 class + 4 register + 256 patch; the embedding concatenates the class token with the mean of the patch tokens after the register tokens, giving 2560-d). Generalized `scripts/analyze_hoptimus_embedding_control.py` to accept any model via `--model-label`/`--model-id` with a dynamic embedding dimension; verified by regression that the H-Optimus path reproduces byte-identically.

Virchow2 (`paige-ai/Virchow2`, 2560-d, same 171 slides, 128 tiles) replicated H-Optimus-0:

| Feature set (low vs zero, n=118) | Repeated-CV BA | Leave-source-site-out BA |
|---|---:|---:|
| Slide-size covariates | 0.888 | 0.882 |
| GigaTIME mean channels | 0.710 | 0.617 |
| H-Optimus-0 embedding | 0.726 | 0.586 |
| Virchow2 embedding | 0.693 | 0.551 |

Virchow2: 0.693 repeated-CV (beats null, p = 0.005; PCA 10/20/30 = 0.693/0.693/0.659), collapsing to 0.551 under source-site holdout, beaten by portable slide-size 0.882. Committed `b58d795`, wired into the same spine docs (results-doc section now shows the consolidated two-model table).

Bottom line of the two controls: three independent image representations (GigaTIME virtual channels, H-Optimus-0, Virchow2) all show the same three-part signature: separate low/zero at ~0.69-0.73, collapse under source-site holdout, and lose to a 3-number portable slide-size baseline. This is strong, replicated evidence that the low-versus-zero axis is generic morphology/tissue-composition tracking TCGA acquisition structure, not GigaTIME-specific virtual immune biology. It is an internal TCGA control, not external validation.

### Engineering notes for future runs

- Both embedding runners write `slide_embeddings.csv` incrementally per slide and support `--resume` (keyed on slide_id). Two accidental mid-run stops this session were recovered with zero loss by re-launching the identical command with `--resume`.
- Do not pass `--save-tile-csv` for full-cohort embedding runs: the runner rewrites the entire growing tile CSV after every slide (quadratic I/O). The slide-mean `slide_embeddings.csv` is all the control needs.
- Measured throughput on Apple MPS: H-Optimus-0 (ViT-g, 1.1B) ~0.21 s/tile, full 171-slide x ~123-tile run ~82 min; Virchow2 (ViT-H, 632M) faster. H-Optimus-0 and Virchow2 weights are gated on Hugging Face but loaded from local cache without a shell token (a stored `~/.cache/huggingface/token` plus accepted licenses).
- The control script is now model-agnostic: a third model drops in via `--embeddings/--model-label/--model-id`.

### Decision: would pulling more TCGA slides help? No.

We explicitly evaluated whether expanding the TCGA cohort would strengthen the analysis. It would not, for the question that matters:

- HER2-zero is the binding limit. The TCGA-BRCA clinical label table has exactly 61 HER2-zero cases total (verified live: positive 174, low 407, zero 61, unknown 455, of 1,097 rows), and the high-trust cohort already uses all 61. The 455 "unknown" cases cannot be cleanly recovered into HER2-zero (most are receptor-negative without an IHC score, so 0 vs 1+ is indeterminate). The low/zero comparison is therefore already N-maxed.
- The bottleneck is confounding, not power. More TCGA-BRCA slides carry the same acquisition structure, so they tighten confidence intervals around a biased estimate (can make a confounded signal look more significant). The within-site rescue does not scale either: only 4 TCGA sites contain both classes (51 cases, 12 low / 39 zero), and more slides cannot manufacture site-balanced low/zero that TCGA does not have.

What would actually help is variation independent of HER2 status: an external, single-scanner/single-institution cohort with H&E + real HER2 IHC/ISH (ideally with IHC-0-vs-1+ granularity), real IHC/mIF to ground the virtual channels, or pathologist tumor-region annotation.

### External-validation scouting

We scouted external cohorts (web/literature search). Full shortlist, attributes, access notes, the HER2-low-vs-zero granularity caveat, and prior-work context are recorded in `docs/external_validation_candidates.md`. Headline points:

- The exact comparison (HER2-low IHC 1+/2+ISH- vs HER2-zero IHC 0) is the hard part externally: most public H&E+HER2 sets are binary +/-, and even with the IHC slide, pathologist interobserver agreement on 0 vs 1+ is low.
- Closest prior work: Valieris et al., Breast Cancer Research 2024 (PMC11331614), weakly-supervised HER2-low from H&E across ACCCC (private single-institution, Brazil), HEROHE, and TCGA. They saw external performance drop and noted TCGA technical diversity affects classification, but did not do the leave-site-out / embedding-control confound analysis this project did. Our method is more rigorous than the published SOTA on this exact task; they are a natural citation and possible collaborator (ACCCC is the cleanest fit cohort).
- Best candidates: ACCCC (single-institution, has neg/low/high, private -> request), BCNB (1,058 WSI, single scanner, free with registration, HER2 0-vs-1+ separability to confirm), ACROBAT (4,212 WSI, single Swedish source, paired H&E + HER2-IHC, public), HEROHE (single scanner, binary HER2 only), Yale HER2-TUMOR-ROIS on TCIA (H&E + trastuzumab response), IMPRESS (HER2+/TNBC + neoadjuvant response + real multiplex IHC, which could validate GigaTIME's virtual immune channels against measured markers).

### Commits this session

- `42ba8a3` Record real H-Optimus one-tile DeepSpot smoke result (in-flight work).
- `0885de7` Note removal of unrelated GATK genomics scaffold (~19 GB reclaimed).
- `8d21b02` Add generic H-Optimus-0 embedding control for HER2-low vs zero.
- `b58d795` Add Virchow2 as a second generic-embedding control.

### Current honest status and next steps

The TCGA-internal evidence is now exhausted: the low-versus-zero signal has been shown to be acquisition-confounded via covariate baselines, matched subsets, source-site holdout, within-site restriction, and two orthogonal foundation-model embedding controls. The two scientifically meaningful paths forward are: (1) external/site-controlled validation with real HER2 IHC/ISH (see `docs/external_validation_candidates.md`), and (2) writing a cautionary-methods paper, for which the two-model embedding control is a strong centerpiece figure and Valieris et al. 2024 is the key comparison. Pulling more TCGA slides is not a productive next step.

Updated key files:

- `scripts/run_hoptimus_tcga_brca.py` (used for full H-Optimus-0 extraction)
- `scripts/run_virchow2_one_slide_smoke.py`, `scripts/run_virchow2_tcga_brca.py`
- `scripts/analyze_hoptimus_embedding_control.py` (generalized to any embedding model)
- `docs/clinical_her2_high_trust_tile128_hoptimus_embedding_control.md`
- `docs/clinical_her2_high_trust_tile128_virchow2_embedding_control.md`
- `docs/external_validation_candidates.md`
- `docs/clinical_her2_high_trust_tile128_results.md`, `docs/advisor_brief.md`, `docs/RUN_REGISTRY.md` (spine updates)
- `results/hoptimus_tcga_brca_high_trust_tile128/`, `results/virchow2_tcga_brca_high_trust_tile128/`
