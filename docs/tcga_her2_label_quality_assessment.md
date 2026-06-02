# TCGA-BRCA HER2 Label Quality Assessment

Status: Working quality-control note for the TCGA-BRCA HER2 labels used in this project.

## Bottom Line

The TCGA-BRCA HER2 labels are usable for exploratory research, especially for asking whether H&E-derived GigaTIME features associate with broad HER2 clinical states. They should not be treated as perfect diagnostic ground truth.

The main reason is that TCGA gives clinical supplement fields for HER2 IHC and FISH/ISH, but those fields are incomplete, sometimes discordant, and were collected across many institutions and historical testing periods. The safest language is:

> We derived research HER2 groups from TCGA-BRCA clinical HER2 IHC/ISH fields using transparent rules.

Not:

> TCGA provides perfectly validated HER2-zero, HER2-low, and HER2-positive diagnostic labels.

## What The Online Sources Support

GDC clinical data are standardized around the GDC data model and data dictionary, but TCGA can also include project-specific clinical supplement files such as BCR Biotab files.

Useful sources:

- GDC Clinical Data overview: https://docs.gdc.cancer.gov/Encyclopedia/pages/Clinical_Data/
- GDC Data Dictionary overview: https://docs.gdc.cancer.gov/Data_Dictionary/
- TCGA-BRCA publication page: https://gdc.cancer.gov/about-data/publications/brca_2012
- TCGA-BRCA Breast Enrollment Form: https://gdc.cancer.gov/system/files/public/file/Breast%20Enrollment%20Form.pdf
- CAP/ASCO HER2 Testing in Breast Cancer 2023 guideline page: https://www.cap.org/protocols-and-guidelines/cap-guidelines/current-cap-guidelines/recommendations-for-human-epidermal-growth-factor-2-testing-in-breast-cancer
- ASCO/CAP 2018 HER2 focused update: https://ascopubs.org/doi/10.1200/JCO.2018.77.8738

The TCGA-BRCA enrollment form explicitly collected HER2/ERBB2 IHC status, IHC percent category, IHC intensity score, HER2/ERBB2 FISH status, HER2 copy number, centromere 17 copy number, and number of cells counted. That means our labels are based on real clinical annotation fields, not on image model output.

However, CAP/ASCO guidance also makes clear that HER2 testing depends on proper tissue handling, assay performance, IHC/ISH interpretation, and reporting. The 2023 CAP/ASCO update does not endorse "HER2-low" as a formal interpretive diagnostic category in the same way as HER2-positive versus HER2-negative. For this project, HER2-low should therefore be described as a research grouping based on IHC 1+ or IHC 2+/ISH-negative clinical fields.

## What We Did Locally

The local label builder is:

```text
scripts/build_tcga_brca_clinical_her2_labels.py
```

It uses the GDC API to download the TCGA-BRCA patient-level BCR Biotab:

```text
nationwidechildrens.org_clinical_patient_brca.txt
```

Local metadata:

- GDC project: `TCGA-BRCA`
- Source file ID: `8162d394-8b64-4da2-9f5b-d164c54b9608`
- Source file name: `nationwidechildrens.org_clinical_patient_brca.txt`
- Local label table: `data/tcga_brca/clinical_her2_labels.csv`
- Local metadata file: `data/tcga_brca/clinical_her2_labels_metadata.json`

## Current Label Counts

The current TCGA-BRCA clinical label table contains 1,097 patient rows:

| Group | Count |
|---|---:|
| HER2-low | 407 |
| HER2-positive | 174 |
| HER2-zero | 61 |
| HER2-unknown | 455 |

This immediately tells us that HER2-zero is the limiting group. That is why the largest balanced local cohort is 61/61/61.

## Current Label Rules

The project currently assigns labels using these rules:

| Assigned group | Rule |
|---|---|
| HER2-positive | IHC score `3+`, ISH positive, or positive HER2 receptor status when detailed score/ISH are missing |
| HER2-low | IHC score `1+` with no positive ISH, or IHC score `2+` with ISH negative |
| HER2-zero | IHC score `0` with no positive ISH |
| HER2-unknown | Missing, not evaluated, indeterminate, equivocal without definitive ISH, or incomplete fields |

This is a conservative and traceable labeling approach, but it is still an operational research definition.

## Issues We Might Encounter

1. Missing HER2 detail is common.

455 of 1,097 cases are currently `HER2-unknown`. Many have negative receptor status but no IHC score, which means we cannot distinguish HER2-zero from HER2-low.

2. HER2-low is a research grouping here, not a clean historical TCGA diagnostic category.

Our HER2-low group uses the modern research definition of IHC `1+` or IHC `2+` with ISH negative. TCGA was not originally built around HER2-low versus HER2-zero as a therapeutic category.

3. Some HER2-positive labels are weaker than others.

145 HER2-positive cases are direct by IHC `3+` or ISH positive. Another 29 are inferred from a positive HER2 receptor-status field when detailed IHC score and ISH fields are missing. Those 29 should be tested in a sensitivity analysis.

4. Discordant IHC/ISH combinations exist.

Current local flags:

| Flag | Count | Current handling |
|---|---:|---|
| IHC `1+` with ISH positive | 3 | Labeled HER2-positive because ISH positive takes priority |
| IHC `3+` with ISH negative | 5 | Labeled HER2-positive because IHC `3+` takes priority |
| IHC `2+` with missing/not-evaluated ISH | 22 | Labeled HER2-unknown |
| ISH equivocal | 5 | Labeled HER2-unknown |
| ISH indeterminate | 4 | Mostly labeled HER2-unknown |

The discordant positive cases are not necessarily wrong, but they should be reviewed or excluded in strict sensitivity analyses.

5. Patient-level HER2 labels are not always slide-level truth.

The HER2 label is attached to the patient/case clinical record. The H&E slide is a diagnostic whole-slide image from the case. That does not guarantee the exact same block, region, or tumor subclone was used for HER2 testing.

6. Clinical HER2 status is not the same as ERBB2 RNA.

ERBB2 RNA expression is useful for biological validation, but clinical HER2 status is based on protein expression and/or gene amplification. A slide can be clinically HER2-low or HER2-zero while still showing variable ERBB2 RNA.

7. Clinical HER2 status is not the same as PAM50 HER2-enriched subtype.

The TCGA-BRCA molecular paper describes integrated molecular classes and HER2-related signaling, but molecular HER2-enriched subtype should not be treated as identical to clinical HER2-positive disease.

## Checks We Should Run Before Claiming Strong Results

### Required Checks

1. Reproduce the label table from the GDC source file and preserve the source file ID.

2. Report label counts by group, rule, and confidence:

```text
clinical_her2_group
clinical_her2_group_rule
clinical_her2_group_confidence
```

3. Run a direct-only sensitivity analysis:

- Keep only direct HER2-positive, HER2-low, and HER2-zero labels.
- Exclude the 29 inferred HER2-positive cases.
- Compare whether the GigaTIME findings survive.

Status: done for the current strict high-trust 171-slide analysis set. The high-trust list keeps only direct clinical HER2 labels.

4. Run a discordance-excluded sensitivity analysis:

- Exclude IHC `1+`/ISH-positive cases.
- Exclude IHC `3+`/ISH-negative cases.
- Exclude all ISH equivocal or indeterminate cases.

Status: done for the current strict high-trust 171-slide analysis set. Discordant/review cases were excluded from the primary high-trust list.

5. Split HER2-low into its two biological/clinical components:

- IHC `1+`
- IHC `2+` and ISH negative

Status: done for the current high-trust run in `docs/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity.md`. Both HER2-low IHC `1+` and HER2-low IHC `2+`/ISH-negative cases remained lower than HER2-zero overall for key GigaTIME channels.

6. Split HER2-zero by ISH availability:

- IHC `0` and ISH negative
- IHC `0` and ISH not evaluated

Status: done for the current high-trust run in `docs/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity.md`. HER2-low remained lower than both HER2-zero IHC `0`/ISH-negative and HER2-zero IHC `0`/ISH-not-evaluated subgroups for key GigaTIME channels.

7. Check selection bias across groups:

- ER status
- PR status
- tumor type/histology if available
- stage/grade if available
- tissue amount and tumor content proxies
- slide size and number of usable tissue tiles

Status: partially done. ER/PR adjustment was run in `docs/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity.md`; slide/file trust and available tile QC were checked in `docs/clinical_her2_trustworthy_slide_list.md`. Tumor type, grade, and richer batch/stain proxies still need follow-up.

8. Check molecular concordance:

- ERBB2 RNA expression
- ERBB2 copy-number/amplification if available
- HER2 or phospho-HER2 RPPA if available
- PAM50 subtype if available

These are validation layers, not replacements for clinical HER2 IHC/ISH.

9. Verify slide-level tissue quality:

- Use only diagnostic primary-tumor H&E slides.
- Confirm slides are tissue-rich.
- Check that GigaTIME is sampling tumor-containing and cellular tissue, not mostly fat, blank background, necrosis, or benign stroma.

10. Validate externally if possible.

The strongest future validation would be an external cohort with:

- H&E slides
- real HER2 IHC scores
- ISH/FISH status
- ideally real mIF or spatial immune-marker data
- ideally treatment-response information for anti-HER2 or ADC therapy

## Completed 61/61/61 Cohort QC Checks

The first pass of these checks has now been run on the downloaded 183-slide laptop-balanced cohort. The reproducible script is:

```text
scripts/build_tcga_her2_trustworthy_slide_list.py
```

Outputs:

- `data/tcga_brca/clinical_her2_laptop_balanced61_trustworthy_slides.csv`
- `data/tcga_brca/clinical_her2_laptop_balanced61_high_trust_slides.csv`
- `docs/assets/clinical_her2_trustworthy_slide_list/trustworthy_slides.csv`
- `docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv`
- `docs/clinical_her2_trustworthy_slide_list.md`

### Slide Integrity And Metadata

All 183 selected files passed the basic local slide checks:

| Check | Passing slides |
|---|---:|
| Local slide file exists | 183 |
| Local file size matches GDC metadata | 183 |
| OpenSlide can read the SVS file | 183 |
| Primary tumor sample metadata | 183 |
| SVS slide image metadata | 183 |
| Female patient after Guardia et al.-aligned cleanup | 180 |

No selected slide is currently excluded because of file corruption, size mismatch, non-primary sample type, or unreadable SVS format.

### Label Trust

| Trust level | Slides |
|---|---:|
| High label+slide trust | 171 |
| Review before primary analysis | 9 |
| Exclude from primary analysis | 3 |

High label+slide trust means: known clinical HER2 group, direct HER2 label rule, no flagged IHC/ISH discordance, female patient, primary-tumor diagnostic SVS metadata, local file present, exact file size match, and readable by OpenSlide.

The nine review slides are still technically usable, but they should be excluded from the strictest primary sensitivity analysis:

| Group | Review slides | Main reason |
|---|---:|---|
| HER2-positive | 5 | IHC/ISH discordance, such as IHC `3+` with ISH negative or IHC `1+` with ISH positive |
| HER2-low | 4 | HER2-low group with positive receptor-status field despite low detailed IHC/ISH rule |
| HER2-zero | 0 | No review flags in selected cohort |

Three additional slides are excluded from the strict primary analysis because Guardia et al., Genome Research 2025, PMID 40664477, excluded male TCGA-BRCA samples:

| Group | Excluded slides | Main reason |
|---|---:|---|
| HER2-positive | 3 | Male TCGA-BRCA patient |
| HER2-low | 0 | None |
| HER2-zero | 0 | None |

### High-Trust Counts By Group

| Group | High-trust slides | Review slides |
|---|---:|---:|
| HER2-positive | 53 | 5 |
| HER2-low | 57 | 4 |
| HER2-zero | 61 | 0 |

### HER2-Low And HER2-Zero Subgroups

| Group | Subgroup | Slides |
|---|---|---:|
| HER2-low | IHC `1+`, ISH negative | 8 |
| HER2-low | IHC `1+`, ISH not evaluated | 30 |
| HER2-low | IHC `2+`, ISH negative | 23 |
| HER2-zero | IHC `0`, ISH negative | 18 |
| HER2-zero | IHC `0`, ISH not evaluated | 43 |

This matters because HER2-low IHC `1+` and IHC `2+`/ISH-negative may not behave identically. Likewise, HER2-zero cases with and without ISH evaluation should be tested separately if the sample size permits.

### Selection Bias Snapshot

| HER2 group | Slides | ER positive | PR positive | Median ERBB2 TPM | Median slide size MB |
|---|---:|---:|---:|---:|---:|
| HER2-positive | 61 | 49 | 33 | 550 | 143 |
| HER2-low | 61 | 49 | 45 | 83.4 | 101 |
| HER2-zero | 61 | 43 | 37 | 62.7 | 234 |

These differences are potential confounders. A model could learn ER/PR context, ERBB2 RNA level, tissue size, or tissue-composition patterns instead of clinically meaningful HER2 status.

### Available GigaTIME Tissue QC

Only 60 of the 183 downloaded slides have existing GigaTIME/tile-level QC from the expanded 20/20/20 run:

| Processed tissue QC | Slides |
|---|---:|
| Strict pass | 50 |
| Moderate pass | 7 |
| Weak/review | 3 |
| Not processed yet | 123 |

Strict processed tissue QC currently means at least 200 sampled tiles, mean tissue fraction at least 0.70, and QC-cellular retained fraction at least 0.50.

The remaining 123 downloaded slides need GigaTIME processing before we can apply tile-level tissue QC.

## Recommended Language

Use:

> TCGA-BRCA clinical supplement HER2 IHC/ISH fields were converted into research HER2-positive, HER2-low, HER2-zero, and unknown groups using prespecified rules.

Use:

> HER2-low is used here as an operational research group based on IHC 1+ or IHC 2+/ISH-negative annotation.

Avoid:

> TCGA gives perfectly reliable HER2-low and HER2-zero labels.

Avoid:

> GigaTIME detects HER2 status.

Better:

> GigaTIME-derived H&E features associate with TCGA clinical HER2 label groups, pending sensitivity analysis and external validation.
