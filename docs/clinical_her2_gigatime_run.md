# Clinical HER2 GigaTIME Run

This document records the first GigaTIME run using the selected clinical HER2-positive / HER2-low / HER2-zero cohort.

## Run Status

The selected cohort contains 30 TCGA-BRCA cases:

- 10 HER2-positive
- 10 HER2-low
- 10 HER2-zero

At the time of this run, 8 of the 30 selected slide files were already downloaded locally:

| Clinical HER2 group | Selected cases | Slides available locally | Slides processed |
|---|---:|---:|---:|
| HER2-positive | 10 | 4 | 4 |
| HER2-low | 10 | 3 | 3 |
| HER2-zero | 10 | 1 | 1 |

The remaining 22 selected slides still need to be downloaded before the full 30-case clinical HER2 analysis can be run.

## GigaTIME Command

The selected-slide run uses the clinical HER2 cohort slide table and skips missing local slides:

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

Local outputs:

- `results/gigatime_tcga_brca_clinical_her2/slide_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2/tile_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2/heatmaps/`

## Clinical HER2 Summary Command

```bash
conda run -n gigatime-tcga python scripts/summarize_clinical_her2_gigatime.py \
  --slide-scores results/gigatime_tcga_brca_clinical_her2/slide_scores.csv \
  --cohort data/tcga_brca/clinical_her2_cohort_cases.csv \
  --out-dir results/gigatime_tcga_brca_clinical_her2/clinical_summary
```

Local outputs:

- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/joined_slide_clinical_her2_gigatime.csv`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_channel_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_pairwise_tests.csv`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_summary.md`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_channel_boxplots.png`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_group_mean_heatmap.png`
- `results/gigatime_tcga_brca_clinical_her2/clinical_summary/erbb2_tpm_by_clinical_her2_group.png`

## Preliminary Findings

This is an availability-limited pilot, not the final comparison. The current joined set contains only 8 slides, and HER2-zero is represented by only 1 slide.

Current top three-group differences among the summarized channels:

| Channel | Kruskal p | Highest mean group | Lowest mean group | Max-min mean |
|---|---:|---|---|---:|
| CD4 | 0.2966 | HER2-zero | HER2-low | 0.04051 |
| CD68 | 0.4140 | HER2-zero | HER2-low | 0.00909 |
| CD11c | 0.4994 | HER2-zero | HER2-low | 0.00661 |
| CK | 0.5778 | HER2-positive | HER2-low | 0.10538 |
| Ki67 | 0.5778 | HER2-positive | HER2-low | 0.00283 |

None of these should be interpreted as a definitive biological result yet. The value of this run is that the clinical HER2 grouping, selected-slide GigaTIME path, and three-group summary workflow are now reproducible.

## Next Step

Download the 22 missing selected slides using:

```bash
gdc-client download \
  -m data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv \
  -d data/tcga_brca/slides
```

Then rerun the GigaTIME and clinical summary commands above. The intended next full pilot will have 10 HER2-positive, 10 HER2-low, and 10 HER2-zero cases.
