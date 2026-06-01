# Clinical HER2 GigaTIME Run

This document records the selected clinical HER2-positive / HER2-low / HER2-zero GigaTIME pilot run.

## Run Status

The selected cohort contains 30 TCGA-BRCA cases, and all selected diagnostic H&E slide files are now available locally:

| Clinical HER2 group | Selected cases | Slides available locally | Slides processed |
|---|---:|---:|---:|
| HER2-positive | 10 | 10 | 10 |
| HER2-low | 10 | 10 | 10 |
| HER2-zero | 10 | 10 | 10 |

The earlier availability-limited pilot processed only 8 of 30 selected slides. The current full pilot processed all 30 selected slides using 64 random tissue tiles per slide, followed by a 256-tile robustness rerun on the same 30 slides.

## Slide Download Command

The 22 missing selected slides were downloaded with the project Python downloader because `gdc-client` was not installed in the local environment:

```bash
conda run -n gigatime-tcga python scripts/download_clinical_her2_cohort_slides.py \
  --only-missing
```

This wrote a local status file:

- `data/tcga_brca/clinical_her2_cohort_slide_download_status.json`

## GigaTIME Command

The selected-slide run uses the clinical HER2 cohort slide table and skips missing local slides if any are absent:

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

Regenerated local outputs:

- `results/gigatime_tcga_brca_clinical_her2/slide_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2/tile_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2/heatmaps/`

For full-cohort statistics, use the regenerated `slide_scores.csv` and `clinical_summary/` outputs. The heatmap directory can contain files from earlier pilot runs if it is not cleaned before rerunning.

## Clinical HER2 Summary Command

```bash
conda run -n gigatime-tcga python scripts/summarize_clinical_her2_gigatime.py \
  --slide-scores results/gigatime_tcga_brca_clinical_her2/slide_scores.csv \
  --cohort data/tcga_brca/clinical_her2_cohort_cases.csv \
  --out-dir results/gigatime_tcga_brca_clinical_her2/clinical_summary
```

Regenerated local outputs:

- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/joined_slide_clinical_her2_gigatime.csv`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_channel_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_pairwise_tests.csv`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_summary.md`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_channel_boxplots.png`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_group_mean_heatmap.png`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/erbb2_tpm_by_clinical_her2_group.png`

## Full 30-Slide Pilot Findings

The joined analysis includes 30 slides from 30 unique TCGA cases:

- 10 HER2-positive cases.
- 10 HER2-low cases.
- 10 HER2-zero cases.

The strongest three-group differences were:

| Channel | Kruskal p | BH q | Highest mean group | Lowest mean group | Max-min mean |
|---|---:|---:|---|---|---:|
| CD68 | 0.0242 | 0.1647 | HER2-zero | HER2-low | 0.00913 |
| PD-L1 | 0.0423 | 0.1647 | HER2-zero | HER2-low | 0.01749 |
| CD11c | 0.0494 | 0.1647 | HER2-zero | HER2-low | 0.00450 |
| CD4 | 0.0794 | 0.1722 | HER2-zero | HER2-low | 0.02739 |
| Ki67 | 0.0920 | 0.1722 | HER2-zero | HER2-low | 0.00416 |

The repeated pattern is that HER2-zero had the highest predicted mean signal for several immune and checkpoint-related virtual channels, while HER2-low had the lowest mean signal. HER2-positive was usually intermediate for these channels.

Top pairwise unadjusted tests were all HER2-low versus HER2-zero:

| Channel | Comparison | Delta mean | Mann-Whitney p | BH q |
|---|---|---:|---:|---:|
| CD68 | HER2-low vs HER2-zero | -0.00913 | 0.00911 | 0.2113 |
| CD11c | HER2-low vs HER2-zero | -0.00450 | 0.01726 | 0.2113 |
| PD-L1 | HER2-low vs HER2-zero | -0.01749 | 0.02113 | 0.2113 |
| CD4 | HER2-low vs HER2-zero | -0.02739 | 0.03121 | 0.2258 |
| Ki67 | HER2-low vs HER2-zero | -0.00416 | 0.03764 | 0.2258 |

No pairwise result remained significant after Benjamini-Hochberg correction. These should be treated as pilot signals for follow-up, not final biological claims.

## Interpretation

The current result suggests a possible virtual microenvironment pattern: in this selected TCGA-BRCA pilot, HER2-zero slides showed higher GigaTIME-predicted macrophage/myeloid, checkpoint, and broader immune-channel signals than HER2-low slides. HER2-positive slides did not show a clear separation from HER2-low across the top channels.

This is useful because it gives the project a specific next hypothesis:

> In TCGA-BRCA, GigaTIME-predicted immune microenvironment features may separate clinically HER2-zero from HER2-low tumors more clearly than they separate HER2-positive from HER2-low tumors.

## Guardrails

- Clinical HER2 labels are based on TCGA IHC/ISH clinical supplement fields.
- GigaTIME outputs are virtual mIF research features, not experimental mIF measurements.
- The cohort is balanced but small: 10 cases per HER2 group.
- The first run used 64 random tissue tiles per slide. The 256-tile rerun improves sampling robustness, but it is still not exhaustive whole-slide analysis.
- The current findings are hypothesis-generating and need validation with either real mIF, orthogonal immune deconvolution, pathology review, or a larger TCGA-BRCA cohort.

## 256-Tile Robustness Update

The next robustness step has now been run on the same 30 selected slides using up to 256 random tissue tiles per slide:

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

The same clinical summary, RNA validation, and visual QC steps were repeated for `results/gigatime_tcga_brca_clinical_her2_tile256/`.

The 256-tile run reproduced the main 64-tile pattern:

| Channel | 64-tile Kruskal p | 256-tile Kruskal p | 64-tile max-min | 256-tile max-min | 256-tile direction |
|---|---:|---:|---:|---:|---|
| CD68 | 0.0242 | 0.0167 | 0.00913 | 0.01044 | HER2-zero > HER2-low |
| PD-L1 | 0.0423 | 0.0211 | 0.01749 | 0.02061 | HER2-zero > HER2-low |
| CD11c | 0.0494 | 0.0384 | 0.00450 | 0.00504 | HER2-zero > HER2-low |
| CD3 | 0.1876 | 0.0686 | 0.02661 | 0.03289 | HER2-zero > HER2-low |
| CD4 | 0.0794 | 0.0819 | 0.02739 | 0.03277 | HER2-zero > HER2-low |

The strongest 256-tile pairwise tests were again HER2-low versus HER2-zero. For CD68, PD-L1, and CD11c, the leading BH q values improved to 0.1133 but remained above 0.05.

The 256-tile RNA validation still did not show FDR-significant correlations between GigaTIME virtual channels and matched RNA marker signatures. The robustness result therefore strengthens the sampling argument but does not validate the marker biology.

## Next Step

The next analysis step is validation. The 256-tile rerun supports that the HER2-zero versus HER2-low signal is not simply a 64-tile sampling accident. However, RNA-seq validation remains weak, and visual QC still does not establish biological validity.

1. Ask an advisor/pathologist to review the source H&E tiles and virtual mIF composites for cases driving high `CD68`, `PD-L1`, and `CD11c`.
2. Compare GigaTIME predictions with richer orthogonal immune evidence, such as RNA-derived immune signatures, tumor purity estimates, or published TCGA immune subtype annotations.
3. Consider a 512-tile or more exhaustive run if compute time allows.
4. Expand beyond the 30 selected cases if more clinical HER2-zero cases can be reliably included.

See `docs/clinical_her2_rna_validation.md` for the first RNA-seq marker-signature comparison.
See `docs/clinical_her2_visual_qc.md` for the first H&E-versus-virtual-mIF visual QC pass.
See `docs/clinical_her2_tile_sampling_robustness.md` for the 256-tile robustness check.
