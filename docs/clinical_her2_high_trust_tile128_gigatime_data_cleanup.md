# Clinical HER2 GigaTIME Data Cleanup

This cleanup step goes back before classifier training. It asks whether the GigaTIME input features should be aggregated from all sampled tissue tiles, or from cleaner tile subsets that are more cellular and more epithelial/tumor enriched.

## Cleanup Views

- `all_sampled_tissue`: all tissue tiles selected by the original GigaTIME run.
- `qc_cellular_tissue`: tiles with H&E tissue fraction >= 0.70 and virtual DAPI mean >= 0.05.
- `ck_enriched_top50`: the top 50% virtual CK tiles within each slide after cellular-tissue QC.
- `ck_enriched_top25`: the top 25% virtual CK tiles within each slide after cellular-tissue QC.

Important: CK and DAPI are still GigaTIME virtual predictions from H&E, not laboratory stains. These filters create tumor-enriched research feature views, not confirmed tumor masks.

## Tile Retention

| Cleanup view | Median retained tiles | Median retained fraction | Median DAPI | Median CK |
| --- | --- | --- | --- | --- |
| All sampled tissue | 128.0 | 1.000 | 0.322 | 0.190 |
| QC cellular tissue | 96.0 | 0.750 | 0.359 | 0.213 |
| CK-enriched top 50% | 49.0 | 0.383 | 0.452 | 0.320 |
| CK-enriched top 25% | 25.0 | 0.195 | 0.492 | 0.409 |

![Tile retention by cleanup view](assets/clinical_her2_high_trust_tile128_gigatime_cleanup/cleanup_retained_tiles_by_filter.png)

## Tile-Level Cleanup Map

This plot shows how the cleanup rules move from broad tissue tiles toward cellular, CK-enriched tiles.

![CK and DAPI tile distribution](assets/clinical_her2_high_trust_tile128_gigatime_cleanup/cleanup_ck_dapi_distribution.png)

## Group Means After Cleanup

![Cleaned GigaTIME group mean heatmap](assets/clinical_her2_high_trust_tile128_gigatime_cleanup/cleanup_key_channel_heatmap.png)

## Top Three-Group Signals Across Cleanup Views

| Cleanup view | Channel | Kruskal p | BH q within view | Highest group | Lowest group | Max-min mean |
| --- | --- | --- | --- | --- | --- | --- |
| All sampled tissue | PD-L1 | 0.0002 | 0.0009 | HER2-positive | HER2-low | 0.0291 |
| All sampled tissue | CD11c | 0.0002 | 0.0009 | HER2-positive | HER2-low | 0.0079 |
| All sampled tissue | CD4 | 0.0004 | 0.0010 | HER2-positive | HER2-low | 0.0329 |
| All sampled tissue | CD3 | 0.0004 | 0.0010 | HER2-positive | HER2-low | 0.0328 |
| All sampled tissue | CD68 | 0.0007 | 0.0012 | HER2-positive | HER2-low | 0.0155 |
| QC cellular tissue | PD-L1 | 0.0007 | 0.0035 | HER2-positive | HER2-low | 0.0320 |
| QC cellular tissue | CD4 | 0.0011 | 0.0035 | HER2-positive | HER2-low | 0.0375 |
| QC cellular tissue | CD3 | 0.0014 | 0.0035 | HER2-positive | HER2-low | 0.0373 |
| QC cellular tissue | CD11c | 0.0016 | 0.0035 | HER2-positive | HER2-low | 0.0090 |
| QC cellular tissue | CD68 | 0.0038 | 0.0058 | HER2-positive | HER2-low | 0.0177 |
| CK-enriched top 50% | CD4 | 0.0019 | 0.0091 | HER2-positive | HER2-low | 0.0335 |
| CK-enriched top 50% | CD3 | 0.0020 | 0.0091 | HER2-positive | HER2-low | 0.0329 |
| CK-enriched top 50% | PD-L1 | 0.0035 | 0.0097 | HER2-positive | HER2-low | 0.0340 |
| CK-enriched top 50% | CK | 0.0043 | 0.0097 | HER2-zero | HER2-low | 0.0701 |
| CK-enriched top 50% | CD11c | 0.0061 | 0.0109 | HER2-positive | HER2-low | 0.0094 |
| CK-enriched top 25% | CD4 | 0.0052 | 0.0264 | HER2-positive | HER2-low | 0.0284 |
| CK-enriched top 25% | CD3 | 0.0109 | 0.0264 | HER2-positive | HER2-low | 0.0277 |
| CK-enriched top 25% | PD-L1 | 0.0113 | 0.0264 | HER2-positive | HER2-low | 0.0328 |
| CK-enriched top 25% | CK | 0.0117 | 0.0264 | HER2-zero | HER2-low | 0.0689 |
| CK-enriched top 25% | CD68 | 0.0165 | 0.0297 | HER2-positive | HER2-low | 0.0188 |

## HER2-Low Versus HER2-Zero Focus

Negative delta means HER2-low is lower than HER2-zero.

| Cleanup view | Channel | HER2-low minus HER2-zero | Mann-Whitney p | BH q within view |
| --- | --- | --- | --- | --- |
| All sampled tissue | CD68 | -0.0054 | 0.0004 | 0.0020 |
| All sampled tissue | PD-L1 | -0.0130 | 0.0003 | 0.0020 |
| All sampled tissue | CD11c | -0.0033 | 0.0003 | 0.0020 |
| QC cellular tissue | CD68 | -0.0053 | 0.0019 | 0.0057 |
| QC cellular tissue | PD-L1 | -0.0135 | 0.0005 | 0.0057 |
| QC cellular tissue | CD11c | -0.0034 | 0.0012 | 0.0057 |
| CK-enriched top 50% | CD68 | -0.0051 | 0.0039 | 0.0150 |
| CK-enriched top 50% | PD-L1 | -0.0104 | 0.0026 | 0.0128 |
| CK-enriched top 50% | CD11c | -0.0023 | 0.0049 | 0.0156 |
| CK-enriched top 25% | CD68 | -0.0045 | 0.0070 | 0.0337 |
| CK-enriched top 25% | PD-L1 | -0.0063 | 0.0147 | 0.0462 |
| CK-enriched top 25% | CD11c | -0.0015 | 0.0339 | 0.0708 |

![Key channel boxplots after cleanup](assets/clinical_her2_high_trust_tile128_gigatime_cleanup/cleanup_key_channel_boxplots.png)

## Interpretation

This cleanup does not validate the virtual markers, but it makes the next classifier input more biologically defensible. The original baseline averaged all sampled tissue tiles. The cleaned tables let us rerun summaries or classifiers using cellular tissue tiles and CK-enriched tumor-context tiles.

If HER2-low versus HER2-zero signal remains after CK enrichment, it is less likely to be explained only by blank tissue or broad non-cellular sampling. If the signal disappears, then the original classifier may have been leaning on non-tumor tissue composition rather than tumor-region biology.

## Outputs

- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/gigatime_cleanup/tile_qc_scores.csv`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/gigatime_cleanup/cleaned_slide_features.csv`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/gigatime_cleanup/filter_retention_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/gigatime_cleanup/cleanup_channel_summary.csv`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/gigatime_cleanup/cleanup_pairwise_tests.csv`

## Next Step

Use `cleaned_slide_features.csv` to rerun the classifier separately for each cleanup view, especially `qc_cellular_tissue`, `ck_enriched_top50`, and `ck_enriched_top25`. The comparison should show whether tumor-enriched GigaTIME features improve HER2 prediction or expose the current model as tissue-composition driven.
