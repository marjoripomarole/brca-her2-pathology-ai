# Case-Driver Visual QC Panels

This report renders a small visual QC set from the strict high-trust HER2-low versus HER2-zero case-driver analysis.

The goal is not to validate the whole cohort visually. The goal is to inspect representative label-consistent cases and opposite-profile/manual-review cases, then decide whether the GigaTIME signal appears tissue-plausible or artifact-prone.

Important: the fluorescence-style panels are GigaTIME virtual predictions from H&E tiles, not real multiplex immunofluorescence.

## Cases Rendered

| Case | Group | Review category | HER2 detail | Zero-like score | Expected views | Opposite views | Selection reason |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TCGA-B6-A409 | HER2-low | label_consistent_her2_low | HER2-low_IHC1_ISH-negative | -1.272 | 4 | 0 | Lowest HER2-zero-like score among stable HER2-low cases |
| TCGA-E2-A108 | HER2-low | label_consistent_her2_low | HER2-low_IHC1_ISH-not-evaluated | -1.272 | 4 | 0 | Lowest HER2-zero-like score among stable HER2-low cases |
| TCGA-AO-A128 | HER2-zero | label_consistent_her2_zero | HER2-zero_IHC0_ISH-negative | 2.913 | 4 | 0 | Highest HER2-zero-like score among stable HER2-zero cases |
| TCGA-BH-A0H0 | HER2-zero | label_consistent_her2_zero | HER2-zero_IHC0_ISH-not-evaluated | 1.761 | 4 | 0 | Highest HER2-zero-like score among stable HER2-zero cases |
| TCGA-A2-A0EW | HER2-zero | opposite_profile_manual_review | HER2-zero_IHC0_ISH-negative | -1.021 | 0 | 4 | Opposite profile and classifier-error priority case |
| TCGA-A7-A13E | HER2-low | opposite_profile_manual_review | HER2-low_IHC2_ISH-negative | 0.940 | 0 | 4 | Opposite profile and classifier-error priority case |
| TCGA-AO-A03N | HER2-zero | opposite_profile_manual_review | HER2-zero_IHC0_ISH-not-evaluated | -1.026 | 0 | 4 | Opposite profile and classifier-error priority case |
| TCGA-AO-A0JG | HER2-low | opposite_profile_manual_review | HER2-low_IHC1_ISH-not-evaluated | 0.299 | 0 | 4 | Opposite profile and classifier-error priority case |

## Signal Channels Used For Tile Selection

| Channel | Low-zero delta | BH q |
| --- | --- | --- |
| CD11c | -0.0033 | 0.0020 |
| CD68 | -0.0054 | 0.0020 |
| CK | -0.0638 | 0.0020 |
| PD-L1 | -0.0130 | 0.0020 |
| CD4 | -0.0238 | 0.0020 |
| CD3 | -0.0243 | 0.0021 |
| CD20 | -0.0182 | 0.0166 |
| Ki67 | -5.34e-04 | 0.0350 |

## Selected Tile Summary

This summary is the key visual-QC caveat. Low-like selected tiles can be tissue-containing but nearly blank for virtual CK, CD68, PD-L1, and CD11c, so they may represent stromal/collagen-rich tissue context rather than tumor-rich biology.

| Review category | Group | Tiles | Mean tissue | Mean zero-like tile score | Mean CK | Mean CD68 | Mean PD-L1 | Mean CD11c |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| label consistent her2 low | HER2-low | 8 | 0.987 | -0.605 | 3.08e-04 | 1.62e-04 | 6.32e-04 | 1.74e-04 |
| label consistent her2 zero | HER2-zero | 8 | 0.955 | 1.675 | 0.0487 | 0.1300 | 0.3153 | 0.1173 |
| opposite profile manual review | HER2-low | 8 | 0.923 | 2.800 | 0.0684 | 0.1786 | 0.3879 | 0.1483 |
| opposite profile manual review | HER2-zero | 8 | 0.986 | -0.618 | 0.0012 | 2.57e-04 | 0.0011 | 3.47e-04 |

## Visual Panels

### TCGA-B6-A409 | HER2-low | label_consistent_her2_low

![TCGA-B6-A409 H&E and virtual mIF panel](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/label_consistent_her2_low_TCGA-B6-A409_he_vs_virtual_mif_qc.png)

![TCGA-B6-A409 tile score overlay](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/label_consistent_her2_low_TCGA-B6-A409_tile_driver_overlay.png)

### TCGA-E2-A108 | HER2-low | label_consistent_her2_low

![TCGA-E2-A108 H&E and virtual mIF panel](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/label_consistent_her2_low_TCGA-E2-A108_he_vs_virtual_mif_qc.png)

![TCGA-E2-A108 tile score overlay](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/label_consistent_her2_low_TCGA-E2-A108_tile_driver_overlay.png)

### TCGA-AO-A128 | HER2-zero | label_consistent_her2_zero

![TCGA-AO-A128 H&E and virtual mIF panel](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/label_consistent_her2_zero_TCGA-AO-A128_he_vs_virtual_mif_qc.png)

![TCGA-AO-A128 tile score overlay](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/label_consistent_her2_zero_TCGA-AO-A128_tile_driver_overlay.png)

### TCGA-BH-A0H0 | HER2-zero | label_consistent_her2_zero

![TCGA-BH-A0H0 H&E and virtual mIF panel](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/label_consistent_her2_zero_TCGA-BH-A0H0_he_vs_virtual_mif_qc.png)

![TCGA-BH-A0H0 tile score overlay](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/label_consistent_her2_zero_TCGA-BH-A0H0_tile_driver_overlay.png)

### TCGA-A2-A0EW | HER2-zero | opposite_profile_manual_review

![TCGA-A2-A0EW H&E and virtual mIF panel](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/opposite_profile_manual_review_TCGA-A2-A0EW_he_vs_virtual_mif_qc.png)

![TCGA-A2-A0EW tile score overlay](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/opposite_profile_manual_review_TCGA-A2-A0EW_tile_driver_overlay.png)

### TCGA-A7-A13E | HER2-low | opposite_profile_manual_review

![TCGA-A7-A13E H&E and virtual mIF panel](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/opposite_profile_manual_review_TCGA-A7-A13E_he_vs_virtual_mif_qc.png)

![TCGA-A7-A13E tile score overlay](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/opposite_profile_manual_review_TCGA-A7-A13E_tile_driver_overlay.png)

### TCGA-AO-A03N | HER2-zero | opposite_profile_manual_review

![TCGA-AO-A03N H&E and virtual mIF panel](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/opposite_profile_manual_review_TCGA-AO-A03N_he_vs_virtual_mif_qc.png)

![TCGA-AO-A03N tile score overlay](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/opposite_profile_manual_review_TCGA-AO-A03N_tile_driver_overlay.png)

### TCGA-AO-A0JG | HER2-low | opposite_profile_manual_review

![TCGA-AO-A0JG H&E and virtual mIF panel](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/opposite_profile_manual_review_TCGA-AO-A0JG_he_vs_virtual_mif_qc.png)

![TCGA-AO-A0JG tile score overlay](assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/opposite_profile_manual_review_TCGA-AO-A0JG_tile_driver_overlay.png)

## How To Read This

- Label-consistent HER2-low cases should generally show low HER2-zero-like tile/channel scores.
- Label-consistent HER2-zero cases should generally show higher HER2-zero-like tile/channel scores.
- Opposite-profile cases are the most important QC cases. They may reflect label noise, slide artifact, non-tumor tissue sampling, or real biological exceptions.
- If the selected H&E tiles are blank, folded, necrotic, mostly stroma, or not tumor-rich, the case should be flagged before making biological claims.

## Output Files

- `docs/clinical_her2_high_trust_tile128_case_driver_visual_qc.md`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/case_driver_visual_qc/case_driver_visual_qc_manifest.csv`
- `results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/case_driver_visual_qc/case_driver_visual_qc_selected_tiles.csv`
- `docs/assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/`
