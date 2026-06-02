# Plain-Language Methodology: BRCA HER2 Pathology AI

Status: Current plain-language explanation, updated through the strict high-trust 171-slide clinical HER2 analysis, tumor-rich proxy sensitivity check, and matched HER2-low/HER2-zero sensitivity check.

This document explains what has been done in this project so far for a reader who does not have a genetics, pathology, or computational biology background.

The short version is: this project takes public breast cancer microscope images, uses computational pathology methods to extract image-derived signals, and compares those signals across clinically defined HER2 groups. The first major image model used here is an existing artificial intelligence model called GigaTIME, which predicts virtual immune-marker patterns from small pieces of H&E slides.

This is an early research project. It is not a clinical test, not a diagnostic tool, and not a validated HER2 prediction system.

## 1. The Research Question

The current working question is:

Can computational pathology features from public TCGA breast cancer H&E pathology slides, starting with GigaTIME virtual mIF predictions, generate interpretable immune-environment signals that look different between HER2-positive, HER2-low, and HER2-zero tumors?

In plainer terms:

- We have breast cancer tissue images.
- We have clinical HER2 information from TCGA clinical records.
- We use image-AI methods, currently centered on GigaTIME, to estimate immune and tumor marker patterns from the tissue images.
- We ask whether those image-derived marker patterns differ across HER2-positive, HER2-low, and HER2-zero breast cancers.

This began as a replication/adaptation pilot using a previously released model. It now also includes exploratory classifier analyses, but those classifiers are research tools, not clinical models.

## 2. Biological Background Without Assuming Genetics Knowledge

### What is HER2?

HER2 is a protein that can be present on the surface of breast cancer cells. In some tumors, HER2 is strongly present or amplified. That can affect tumor biology and treatment options.

HER2 status is usually assessed in clinical pathology using protein or gene-copy tests:

- IHC, or immunohistochemistry, estimates how much HER2 protein is visible in the tumor tissue.
- ISH/FISH, or in situ hybridization, checks whether the HER2 gene region is amplified.

### What are HER2-positive, HER2-low, and HER2-zero?

For this project, the groups are defined from TCGA clinical HER2 fields:

- `HER2-positive`: IHC `3+`, ISH positive, or a positive HER2 receptor status when more detailed fields are missing.
- `HER2-low`: IHC `1+`, or IHC `2+` with ISH negative.
- `HER2-zero`: IHC `0` with no positive ISH evidence.
- `HER2-unknown`: missing, not evaluated, equivocal, contradictory, or incomplete fields.

This distinction matters because HER2-low and HER2-zero can be biologically and clinically different, even though both are not classic HER2-positive disease.

### What is ERBB2?

`ERBB2` is the gene that contains the instructions for making the HER2 protein.

At the start of this project, `ERBB2` RNA expression was used as a first HER2-biology proxy. That was useful for building the pipeline, but RNA expression is not the same as clinical HER2 status. The current better analysis uses clinical IHC/ISH-derived HER2 groups.

## 3. Pathology Background Without Assuming Medical Training

### What is an H&E slide?

An H&E slide is a standard pathology slide. The tissue is stained with two dyes:

- Hematoxylin stains cell nuclei, usually blue or purple.
- Eosin stains other tissue structures, usually pink.

Pathologists use H&E slides every day to look at tissue structure under a microscope.

In this project, the H&E slides are digital whole-slide images. A whole-slide image is a very large scanned image of the entire tissue section.

### Why split a slide into tiles?

Whole-slide images are huge. An AI model usually cannot process the entire image at once.

So the workflow cuts each slide into many small image patches called tiles. The GigaTIME inference script uses 256 by 256 pixel tiles by default.

Not every tile is useful. Some areas are blank background or have little tissue. The script estimates the tissue fraction of each tile and keeps tiles that pass a minimum tissue threshold.

For the first clinical HER2 pilot, the workflow processed 64 random tissue tiles per slide. The same 30 slides were then rerun with up to 256 random tissue tiles per slide to check whether the result was stable when more of each slide was sampled. The current larger strict high-trust analysis uses 128 sampled tissue tiles per slide across 171 primary-analysis slides.

## 4. What is GigaTIME?

GigaTIME is an existing released model. This project uses the official GigaTIME implementation and model weights.

The key idea is that GigaTIME takes ordinary H&E pathology image tiles as input and predicts virtual multiplex immunofluorescence-like marker maps as output.

That phrase has several parts:

- Multiplex immunofluorescence, often shortened to mIF, is a lab technique that can stain tissue for many biological markers at the same time.
- A marker is a biological signal associated with a cell type or process. For example, `CD3` and `CD8` are often used as immune-cell markers.
- "Virtual" means the model predicts marker-like signal from the H&E image computationally. The tissue was not actually stained for those markers in this project.

GigaTIME predicts 23 channels:

`DAPI`, `TRITC`, `Cy5`, `PD-1`, `CD14`, `CD4`, `T-bet`, `CD34`, `CD68`, `CD16`, `CD11c`, `CD138`, `CD20`, `CD3`, `CD8`, `PD-L1`, `CK`, `Ki67`, `Tryptase`, `Actin-D`, `Caspase3-D`, `PHH3-B`, and `Transgelin`.

For the clinical HER2 summary analysis, the project focuses on a smaller set of interpretable channels:

`CD3`, `CD8`, `CD4`, `CD20`, `CD68`, `CD11c`, `PD-1`, `PD-L1`, `CK`, and `Ki67`.

Very roughly:

- `CD3`, `CD4`, and `CD8` relate to T cells, a type of immune cell.
- `CD20` relates to B cells.
- `CD68` and `CD11c` relate to macrophage/myeloid immune biology.
- `PD-1` and `PD-L1` relate to immune checkpoint biology.
- `CK` relates to epithelial/tumor-cell structure.
- `Ki67` relates to cell proliferation.

These are model-predicted research features. They should not be interpreted as direct laboratory measurements without validation.

## 5. Data Source

The data source is TCGA-BRCA.

TCGA means The Cancer Genome Atlas. It is a large public cancer research dataset.

BRCA is TCGA's breast cancer project. In this context, BRCA means breast invasive carcinoma; it is not the same thing as the `BRCA1` or `BRCA2` genes.

This project uses three kinds of TCGA-BRCA data from the Genomic Data Commons, or GDC:

- Diagnostic H&E whole-slide images, stored as `.svs` slide files.
- RNA-seq gene expression files, used to extract `ERBB2` expression.
- Clinical supplement fields, used to assign HER2-positive, HER2-low, HER2-zero, or HER2-unknown labels.

## 6. What Has Been Done So Far

### Step 1: Query TCGA-BRCA files from GDC

The workflow queried GDC for:

- TCGA-BRCA slide images where `data_type` is `Slide Image` and `data_format` is `SVS`.
- TCGA-BRCA RNA-seq expression files where the workflow type is `STAR - Counts`.
- TCGA-BRCA clinical supplement files containing HER2-related fields.

Important local metadata files include:

```text
data/tcga_brca/tcga_brca_diagnostic_slides_manifest.tsv
data/tcga_brca/tcga_brca_diagnostic_slides_files.csv
data/tcga_brca/tcga_brca_star_counts_manifest.tsv
data/tcga_brca/tcga_brca_star_counts_files.csv
data/tcga_brca/file_metadata_slides.json
data/tcga_brca/file_metadata_star_counts.json
```

### Step 2: Extract ERBB2 expression

The GDC script downloaded selected STAR-count RNA-seq files and extracted the row corresponding to `ERBB2`.

The extracted expression table is:

```text
data/tcga_brca/erbb2_expression.csv
```

This was used first to make an ERBB2-high versus ERBB2-low pilot. That pilot proved the workflow could run, but it is not the main clinical HER2 comparison.

### Step 3: Build clinical HER2 labels

The clinical HER2 labeling script is:

```bash
scripts/build_tcga_brca_clinical_her2_labels.py
```

It downloads the TCGA-BRCA patient-level clinical supplement and extracts HER2 IHC/ISH fields.

Main outputs:

```text
data/tcga_brca/clinical_her2_labels.csv
data/tcga_brca/clinical_her2_labels_metadata.json
data/tcga_brca/clinical/nationwidechildrens.org_clinical_patient_brca.txt
```

The resulting label table found:

| Clinical HER2 group | TCGA-BRCA clinical rows |
|---|---:|
| HER2-positive | 174 |
| HER2-low | 407 |
| HER2-zero | 61 |
| HER2-unknown | 455 |

### Step 4: Select balanced clinical HER2 cohorts

The cohort selection script is:

```bash
scripts/select_clinical_her2_cohort.py
```

It joins clinical HER2 labels, ERBB2 expression, and slide metadata, then selects a balanced pilot:

| Clinical HER2 group | Selected cases |
|---|---:|
| HER2-positive | 10 |
| HER2-low | 10 |
| HER2-zero | 10 |

The selection prefers direct clinical labels, local slide availability, smaller slide files, and deterministic case IDs.

Main outputs:

```text
data/tcga_brca/clinical_her2_cohort_cases.csv
data/tcga_brca/clinical_her2_cohort_slides_files.csv
data/tcga_brca/clinical_her2_cohort_slide_manifest.tsv
data/tcga_brca/clinical_her2_cohort_summary.json
docs/clinical_her2_cohort_selection.md
```

The project later expanded from this first 30-slide pilot to a 60-slide 20/20/20 cohort, and then to a laptop-realistic 183-slide downloaded cohort: 61 HER2-positive, 61 HER2-low, and 61 HER2-zero slides. HER2-zero is the limiting group because the TCGA-BRCA label table contains only 61 usable HER2-zero clinical rows under the current rules.

### Step 4b: Clean the larger cohort before trusting it

The larger 183-slide cohort was checked for HER2 label trust and slide-file integrity.

The cleanup script is:

```bash
scripts/build_tcga_her2_trustworthy_slide_list.py
```

It checks whether each selected slide has:

- a known clinical HER2 group;
- a direct HER2 label rule;
- no flagged IHC/ISH discordance;
- primary-tumor diagnostic SVS metadata;
- a local slide file;
- a file size that matches GDC metadata;
- an SVS file that OpenSlide can read.

The result was:

| Trust category | Slides |
|---|---:|
| High label+slide trust | 171 |
| Review before primary analysis | 9 |
| Exclude for strict primary analysis | 3 |

The high-trust primary analysis therefore used:

| Clinical HER2 group | High-trust slides |
|---|---:|
| HER2-positive | 53 |
| HER2-low | 57 |
| HER2-zero | 61 |

The machine-readable trustworthy-slide lists are:

```text
docs/assets/clinical_her2_trustworthy_slide_list/trustworthy_slides.csv
docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv
```

### Step 5: Download the selected slides

The selected clinical HER2 cohort needed 30 diagnostic H&E slides. Eight were already present locally, and 22 missing slides were downloaded with:

```bash
scripts/download_clinical_her2_cohort_slides.py
```

The downloader uses GDC file IDs and writes a status file:

```text
data/tcga_brca/clinical_her2_cohort_slide_download_status.json
```

Current selected-slide status:

| Clinical HER2 group | Selected slides | Downloaded slides |
|---|---:|---:|
| HER2-positive | 10 | 10 |
| HER2-low | 10 | 10 |
| HER2-zero | 10 | 10 |

### Step 6: Run GigaTIME on the clinical HER2 cohort

The GigaTIME inference script is:

```bash
scripts/run_gigatime_tcga_brca.py
```

For each slide, the script:

1. Opens the digital pathology slide.
2. Divides the slide into 256 by 256 pixel tiles.
3. Estimates whether each tile contains enough tissue.
4. Keeps tissue-containing tiles.
5. Randomly samples tiles when a tile limit is used.
6. Normalizes the tile image.
7. Runs the tile through the GigaTIME model.
8. Gets predicted marker maps for the GigaTIME channels.
9. Summarizes each tile into marker scores.
10. Aggregates tile scores into one row per slide.

Current strict high-trust clinical HER2 analysis:

- Slides in primary analysis: 171.
- Cases in primary analysis: 171.
- Groups analyzed: 53 HER2-positive, 57 HER2-low, 61 HER2-zero.
- Tiles per slide: 128 sampled tissue tiles.
- Total primary-analysis tile predictions: 21,888.
- Device used: Apple MPS in the current local run.

The raw GigaTIME output folder still contains the earlier 174-slide inference run. The current statistical summaries filter those existing predictions to the 171 strict trustworthy slides after excluding three male HER2-positive cases. That matches the relevant female-patient TCGA sample-selection principle used by Guardia et al.; our HER2 label, file-integrity, OpenSlide, and tile-quality checks are additional project-specific H&E safeguards.

We also ran a simple "does this disappear if the labels are shuffled?" check for the HER2-low versus HER2-zero classifiers. The real labels performed clearly better than shuffled labels. That means the classifier result is not obviously random, but it still does not make the model a clinical diagnostic test and does not prove HER2 isoform biology.

Then we ran a stricter version where the computer had to choose which image features to use using only the training data, before testing on held-out cases. The result still stayed above the shuffled-label comparison. This makes the classifier result more trustworthy as a research signal, but it still needs outside validation and pathologist-reviewed tumor regions.

The newest and most important caution is that non-image information can also separate HER2-low from HER2-zero in this TCGA cohort. The HER2-zero slides are often much larger and come from different TCGA source sites than the HER2-low slides. A model using only slide size or source-site information performs better than the GigaTIME image-feature model. This means the current result may be partly caused by how TCGA slides were collected or selected, not only by tumor biology.

We then compared more similar HER2-low and HER2-zero slides by matching them on source site or slide size. This reduced some obvious imbalance, but it did not fully solve the problem. GigaTIME image features still performed modestly above chance, but slide-size and source-site information remained competitive or stronger, and individual matched marker tests were not statistically convincing after multiple-testing correction. In plain terms: the signal is interesting, but TCGA alone is not clean enough to prove it is true HER2 biology.

We also tested whether the classifier works when an entire TCGA source site is left out during training. This is like asking whether the pattern travels to a different hospital/scanner/source setting. GigaTIME performance dropped in that harder test. In one important view, the GigaTIME classifier dropped from about 74.5% balanced accuracy to about 66.9%, while slide-size information alone still reached about 88.2%. This is a strong warning that the current classifier may still be learning technical differences between slides, not only tumor biology.

Then we ran an even more local check: only compare HER2-low and HER2-zero cases inside TCGA source sites that contain both groups. This is closer to comparing patients from the same source setting. The problem is that only four source sites qualify, leaving only 12 HER2-low and 39 HER2-zero cases. Some GigaTIME signal remains, especially when using all mean channels, but the smaller key-marker classifier weakens. In plain terms: this keeps the idea alive, but it does not remove the worry that TCGA site differences are influencing the result.

We also rechecked the available RNA files for the HER2 gene itself, `ERBB2`. This was a sanity check, not an isoform test. The result was helpful: ERBB2 RNA was very good at separating HER2-positive tumors from the other tumors, which supports that the broad HER2-positive labels make sense. But ERBB2 RNA was weak at separating HER2-low from HER2-zero: the ERBB2-only AUC was about 0.605. In simple terms, the HER2-low versus HER2-zero image signal is not just because HER2-low has dramatically more ERBB2 RNA. That makes the GigaTIME signal more interesting, but it still does not prove the signal is true HER2 biology because source-site, slide-size, and tissue-composition caveats remain.

Earlier historical runs processed 30 slides at 64 and 256 tiles per slide, then 60 slides at up to 256 tiles per slide.

Main outputs:

```text
results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/slide_scores.csv
results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/tile_scores.csv
results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/clinical_summary/
results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/gigatime_cleanup/
results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/cleaned_classifier_comparison/
results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/matched_low_zero_sensitivity/
```

### Step 7: Summarize the GigaTIME predictions by clinical HER2 group

The clinical summary script is:

```bash
scripts/summarize_clinical_her2_gigatime.py
```

It combines:

- The slide-level GigaTIME output.
- The selected clinical HER2 cohort table.
- The clinical HER2 group labels.

It writes:

```text
results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/clinical_summary/joined_slide_clinical_her2_gigatime.csv
results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/clinical_summary/clinical_her2_channel_summary.csv
results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/clinical_summary/clinical_her2_pairwise_tests.csv
results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/clinical_summary/clinical_her2_summary.md
```

## 7. What the Current Results Show

The current main analysis is the strict high-trust 171-slide analysis:

- 53 HER2-positive slides.
- 57 HER2-low slides.
- 61 HER2-zero slides.
- 21,888 primary-analysis GigaTIME tile predictions.

The strongest result is HER2-low versus HER2-zero. In the all-sampled-tissue view, HER2-low was lower than HER2-zero for multiple GigaTIME-predicted immune, myeloid, checkpoint, and tissue-context channels:

| Channel | HER2-low minus HER2-zero | Mann-Whitney p | BH q |
|---|---:|---:|---:|
| CD68 | -0.00537 | 0.000371 | 0.002227 |
| CK | -0.06377 | 0.000129 | 0.002227 |
| PD-L1 | -0.01301 | 0.000302 | 0.002227 |
| PD-1 | -0.03948 | 0.000225 | 0.002227 |
| CD11c | -0.00325 | 0.000272 | 0.002227 |
| CD4 | -0.02379 | 0.000615 | 0.002634 |
| CD3 | -0.02433 | 0.000778 | 0.002918 |

Plain-language interpretation:

> In the larger high-trust TCGA-BRCA run, HER2-low tumors had lower GigaTIME-predicted immune/myeloid/checkpoint and CK-associated image signals than HER2-zero tumors. This is a stronger and more reproducible research signal than the first pilot, but it is still not a clinical diagnostic test.

The tile-cleanup analysis asked whether this pattern was only caused by bad tiles. The signal survived cellular-tissue QC and partly survived CK enrichment. It became weaker under the strictest CK top-25% filter, which suggests the pattern may involve broader tissue microenvironment context rather than only epithelial tumor-cell signal.

The cleaned classifier comparison found:

- HER2-low versus HER2-zero: best balanced accuracy 0.727, macro AUC 0.787.
- HER2-positive versus HER2-negative: weak.
- Three-class HER2-positive versus HER2-low versus HER2-zero: weak to moderate.

So the strongest current claim is an image-derived HER2-low versus HER2-zero association. The project should not claim reliable HER2 diagnosis yet.

The sensitivity analysis then asked whether this association was only caused by ER/PR status or by one HER2 label subgroup. The answer was reassuring but still cautious:

- In all sampled tissue, 7 of 8 tested key GigaTIME channels stayed statistically significant after ER/PR adjustment.
- The pattern was visible in HER2-low IHC `1+` cases and HER2-low IHC `2+`/ISH-negative cases.
- The pattern was visible against both HER2-zero IHC `0`/ISH-negative and HER2-zero IHC `0`/ISH-not-evaluated cases.
- The strictest CK-enriched view still became weaker, so the result likely involves broader tissue context, not only tumor epithelial cells.

A parameter/settings robustness check then compared the earlier 60-slide run that used 256 tiles per slide against the current strict high-trust analysis that used 128 tiles per slide. There were 56 overlapping slides. The same slide-level GigaTIME channel scores were very similar across runs, and all 8 tested key channels kept the same HER2-low versus HER2-zero direction. This makes it less likely that the main pattern is only a random tile-sampling artifact.

See:

```text
docs/clinical_her2_high_trust_tile128_results.md
docs/clinical_her2_high_trust_tile128_gigatime_data_cleanup.md
docs/clinical_her2_high_trust_tile128_cleaned_classifier_comparison.md
docs/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity.md
docs/clinical_her2_high_trust_tile128_vs_expanded20_tile256_agreement.md
```

### Historical 30-slide pilot

The full clinical HER2 pilot includes 30 joined slides:

- 10 HER2-positive.
- 10 HER2-low.
- 10 HER2-zero.

The strongest three-group differences were:

| Channel | Kruskal p | Highest mean group | Lowest mean group |
|---|---:|---|---|
| CD68 | 0.0242 | HER2-zero | HER2-low |
| PD-L1 | 0.0423 | HER2-zero | HER2-low |
| CD11c | 0.0494 | HER2-zero | HER2-low |
| CD4 | 0.0794 | HER2-zero | HER2-low |
| Ki67 | 0.0920 | HER2-zero | HER2-low |

The pattern is that HER2-zero had higher predicted mean signal than HER2-low for several immune and checkpoint-related virtual channels. HER2-positive was usually between those two groups rather than clearly separated from HER2-low.

The strongest pairwise comparisons were HER2-low versus HER2-zero:

| Channel | Direction | Mann-Whitney p | BH q |
|---|---|---:|---:|
| CD68 | HER2-zero higher than HER2-low | 0.0091 | 0.2113 |
| CD11c | HER2-zero higher than HER2-low | 0.0173 | 0.2113 |
| PD-L1 | HER2-zero higher than HER2-low | 0.0211 | 0.2113 |
| CD4 | HER2-zero higher than HER2-low | 0.0312 | 0.2258 |
| Ki67 | HER2-zero higher than HER2-low | 0.0376 | 0.2258 |

No pairwise comparison remained statistically significant after multiple-testing correction. This means the result is a signal worth following, not proof of a biological conclusion.

Plain-language interpretation:

> In this first 30-slide pilot, GigaTIME predicted more immune/checkpoint-like signal in HER2-zero tumors than in HER2-low tumors, especially for CD68, PD-L1, and CD11c. The sample is small, so this should be treated as a hypothesis for validation.

### 256-tile robustness check

The same 30 slides were rerun with more tile sampling. This was done because a whole-slide image is large, and 64 tiles might not represent the whole tissue section.

The denser run gave the same main result:

| Channel | 64-tile p | 256-tile p | 64 max-min | 256 max-min | Direction |
|---|---:|---:|---:|---:|---|
| CD68 | 0.0242 | 0.0167 | 0.00913 | 0.01044 | HER2-zero > HER2-low |
| PD-L1 | 0.0423 | 0.0211 | 0.01749 | 0.02061 | HER2-zero > HER2-low |
| CD11c | 0.0494 | 0.0384 | 0.00450 | 0.00504 | HER2-zero > HER2-low |

In plain language: when the model looked at more pieces of each slide, the same HER2-zero versus HER2-low signal was still there. That makes the finding more trustworthy than before, but it still does not prove the biology is real.

The most important caution is that matched RNA-seq validation was still weak after the 256-tile run. So the result is more robust to tile sampling, but still not validated by an independent biological measurement.

### Broader RNA program validation

The next validation step asked a broader question:

> If GigaTIME predicts immune-like signal, do the same cases also show higher RNA programs for immune biology?

Instead of only checking one marker at a time, this step tested larger RNA programs for T cells, checkpoint/IFNG biology, macrophages, B cells, proliferation, epithelial/tumor signal, stroma, and endothelial tissue.

The answer was still cautious:

- The GigaTIME virtual myeloid/checkpoint program still showed HER2-zero higher than HER2-low.
- Broad RNA immune programs did not show a matching significant HER2-zero higher than HER2-low pattern.
- The strongest significant virtual-vs-RNA associations were negative correlations with endothelial RNA signal.

In plain language: the AI signal is repeatable inside the image model, but the RNA data still does not confirm that it represents real immune-marker biology.

### First classifier baseline

The next step was to stop only comparing group averages and train a simple classifier.

The classifier asked:

> Can GigaTIME features from a slide predict the HER2 group for a slide the model did not train on?

The answer was mixed:

- For HER2-low versus HER2-zero, GigaTIME features did fairly well in this tiny pilot: balanced accuracy 0.800.
- For HER2-positive versus HER2-negative, GigaTIME features did not do well: balanced accuracy 0.475.
- For the full three-class problem, GigaTIME features were at chance: balanced accuracy 0.333.

In plain language: the current GigaTIME features may contain some signal for separating HER2-low from HER2-zero, but they are not yet reliable for clinical HER2 diagnosis.

### Expanded 20/20/20 clinical HER2 run

After the first classifier and cleanup analyses, the cohort was expanded from 30 slides to 60 slides:

- 20 HER2-positive.
- 20 HER2-low.
- 20 HER2-zero.
- Up to 256 tissue tiles per slide.
- 15,225 total GigaTIME tile predictions.
- Matched STAR-count RNA-seq expression downloaded for all 60 selected cases.

This made the main finding stronger. In the expanded run, several HER2-low versus HER2-zero differences passed multiple-testing correction in all-tissue or QC-cellular views:

- `CD3`
- `CD4`
- `CD11c`
- `CD68`
- `PD-L1` in the QC-cellular view

The HER2-low versus HER2-zero classifier also held up:

- All sampled tissue: balanced accuracy 0.800, macro AUC 0.820.
- QC-cellular tissue: balanced accuracy 0.775, macro AUC 0.820.
- CK-enriched top 25%: balanced accuracy 0.800, macro AUC 0.820.

Plain-language interpretation:

> The larger run supports the idea that HER2-low and HER2-zero tumors may have different image-derived tissue-context patterns. The result is stronger than before, but it is still not a clinical diagnostic model and still needs biological validation.

The expanded run also changed the three-group story. HER2-low often looked like the lowest immune/checkpoint group, while HER2-positive became highest for several broader virtual immune programs. So the result should not be simplified to "HER2-zero is always highest."

## 8. What the Output Tables Mean

### `slide_scores.csv`

This table has one row per processed slide.

Each row contains:

- The slide file path.
- The slide ID.
- The TCGA case ID.
- The number of tiles analyzed.
- The average tissue fraction.
- For each GigaTIME channel, the average predicted activation across tiles.
- For each GigaTIME channel, thresholded summaries.

This is the main slide-level model output.

### `tile_scores.csv`

This table has one row per tile.

Each row contains:

- The tile coordinate on the slide.
- The estimated tissue fraction.
- The slide ID.
- The TCGA case ID.
- The GigaTIME marker scores for that tile.

This table is more detailed than `slide_scores.csv` and is useful for heatmaps or spatial inspection.

### `joined_slide_clinical_her2_gigatime.csv`

This table joins model output to clinical HER2 labels.

Each row represents a processed slide with:

- GigaTIME slide-level marker scores.
- The matching TCGA case ID.
- The clinical HER2 group.
- ER, PR, HER2 IHC/ISH, and ERBB2 expression context where available.

This is the table used for the clinical HER2 comparison.

### `clinical_her2_channel_summary.csv`

This table summarizes three-group differences by marker channel.

For each marker, it asks:

- Which HER2 group has the highest average virtual marker activation?
- Which HER2 group has the lowest average virtual marker activation?
- Is the three-group difference suggestive by Kruskal-Wallis testing?
- How large is the difference between the highest and lowest group means?

### `clinical_her2_pairwise_tests.csv`

This table compares pairs of HER2 groups:

- HER2-positive versus HER2-low.
- HER2-positive versus HER2-zero.
- HER2-low versus HER2-zero.

It includes Mann-Whitney p values and Benjamini-Hochberg corrected q values.

## 9. Visual Outputs

The project also generated documentation-facing virtual mIF images:

```text
docs/assets/virtual_mif_channels/
docs/assets/virtual_mif_composites/
```

The files under `docs/assets/virtual_mif_channels/` show all 23 predicted channels, including group-level channel means and slide-by-channel matrices.

The files under `docs/assets/virtual_mif_composites/` look closer to real multiplex immunofluorescence images. They are made by rerunning GigaTIME on selected H&E tiles, keeping the full predicted channel maps, and compositing marker colors on a black background.

These images are still virtual predictions, not real mIF measurements.

## 10. What This Study Has Not Done Yet

The current project has not yet:

- Processed a large external validation cohort beyond the local TCGA-BRCA runs.
- Used more exhaustive whole-slide sampling.
- Performed pathologist-confirmed tumor-region review on every sampled tile.
- Validated GigaTIME predictions against real multiplex immunofluorescence staining in these TCGA slides.
- Trained or fine-tuned a new model.
- Produced a clinically deployable classifier.

The project has now performed indirect RNA-seq validation using simple marker-gene signatures and broader RNA immune/tissue programs. These checks did not strongly confirm the current virtual immune-channel pattern, so they should be treated as cautionary evidence rather than completed validation.

These limitations are important. The current goal is proof of workflow and early biological exploration, not a final scientific claim.

## 11. Why This Is Still Useful

This pilot is useful because it proves several practical pieces:

- TCGA-BRCA H&E slides can be queried and downloaded through GDC.
- TCGA-BRCA clinical HER2 IHC/ISH fields can be used to create HER2-positive, HER2-low, and HER2-zero groups.
- Balanced 10/10/10 and 20/20/20 clinical HER2 cohorts can be selected reproducibly.
- A larger 61/61/61 downloaded cohort can be cleaned into a 171-slide strict high-trust primary analysis set.
- The released GigaTIME model can be run on TCGA-BRCA pathology tiles.
- The GigaTIME outputs can be aggregated into slide-level marker features.
- Those marker features can be compared across clinical HER2 groups.
- The results can be summarized in tables, plots, and visual examples.

In other words, the technical pipeline now works end to end for the clinically meaningful HER2 grouping.

## 12. Current Scientific Interpretation

The safest interpretation is:

This is an exploratory feasibility run showing that an existing H&E-to-virtual-mIF model can be applied to TCGA-BRCA breast cancer slides and connected to clinical HER2 labels.

The current strict high-trust 171-slide result supports a specific hypothesis: HER2-low and HER2-zero tumors may differ in image-derived tissue-context and immune/checkpoint patterns.

RNA-seq validation still did not strongly support that signal. GigaTIME virtual channel predictions did not show strong positive correlations with matched RNA marker signatures, and broader RNA programs did not confirm the virtual immune/checkpoint pattern. This does not automatically invalidate the model, but it means the project should be careful and should not claim that the virtual mIF signal has been validated.

The first visual QC check looked at the H&E tiles driving high virtual `CD68`, `PD-L1`, and `CD11c` predictions. Those tiles contained real cellular tissue rather than obvious blank background. That supports continuing the analysis, but it still does not prove that the virtual markers are biologically correct.

The newest tumor-rich proxy check asked a more specific question: if we keep only virtual CK-rich, cellular, epithelial-like tiles, does the HER2-low versus HER2-zero pattern disappear? The answer is mixed. Individual marker differences become weaker in the strictest fixed-count CK-rich tile views, which keeps the tissue-composition warning alive. But a multichannel HER2-low versus HER2-zero classifier still performs above chance, with balanced accuracy around 0.708 to 0.761 across virtual tumor-rich proxy views.

Plain-language interpretation:

> There may still be a useful image pattern in the CK-rich/tumor-like regions, but we have not proven that yet because these are virtual GigaTIME filters, not pathologist-confirmed tumor regions.

The newest matched-source/slide-size check asks whether the HER2-low versus HER2-zero signal remains when the slides being compared are more similar in obvious technical ways. The result is cautious: GigaTIME still shows a modest signal, but non-image site/slide-size baselines remain strong and individual matched marker differences are not FDR-significant. This means the project should now move toward external or site-balanced validation rather than making strong claims from TCGA alone.

The source-site holdout check strengthens that caution. When entire TCGA source sites are held out, GigaTIME performance drops, but slide-size information remains very predictive. This means the current classifier is not yet robust enough to be described as source-independent HER2 biology.

We also checked whether the local RNA files can directly test the HER2 isoform idea from the Guardia et al. paper. They cannot yet. The files we have are gene-level STAR count/TPM files. They can tell us about total ERBB2 RNA expression and broad RNA programs, but they do not tell us which HER2 transcript isoforms are present. For that, we would need transcript-level isoform labels or RNA-seq read/junction data.

The 256-tile rerun, expanded 20/20/20 run, and strict high-trust 171-slide analysis suggest the HER2-low versus HER2-zero pattern remains when more tissue tiles and more cases are included. The next scientific step is to test whether the pattern becomes more trustworthy when:

- Even more tiles per slide are sampled or more complete whole-slide coverage is used.
- More cases are included from a source-balanced or external dataset if reliable HER2-zero cases are available.
- A human reviews representative H&E tiles and virtual mIF composites for plausibility.
- Richer validation layers are added, such as tumor purity adjustment, immune deconvolution, published immune subtype annotations, or an external dataset with real mIF.
- Classifier inputs are improved by focusing on tumor-rich tiles instead of all tissue tiles.

## 13. Reproducible Workflow Summary

The workflow is organized around these scripts:

```text
scripts/gdc_query_tcga_brca.py
scripts/build_tcga_brca_clinical_her2_labels.py
scripts/select_clinical_her2_cohort.py
scripts/build_tcga_her2_trustworthy_slide_list.py
scripts/download_clinical_her2_cohort_slides.py
scripts/run_gigatime_tcga_brca.py
scripts/summarize_clinical_her2_gigatime.py
scripts/cleanup_gigatime_tile_features.py
scripts/train_her2_cleaned_classifier_comparison.py
scripts/render_he_slide_images.py
scripts/render_clinical_her2_visual_qc.py
scripts/render_virtual_mif_channel_images.py
scripts/render_virtual_mif_composites.py
```

The current clinical HER2 flow is:

```text
Query GDC files
  -> extract ERBB2 expression
  -> build clinical HER2 labels from IHC/ISH fields
  -> select balanced HER2-positive / HER2-low / HER2-zero cases
  -> download selected H&E slides
  -> clean labels and slide files into a high-trust primary cohort
  -> run GigaTIME on slide tiles
  -> aggregate tile predictions to slide scores
  -> join slide scores to clinical HER2 groups
  -> summarize and visualize group differences
  -> run tile-cleanup sensitivity analyses and baseline classifiers
  -> use visual QC, RNA validation, and external validation as next checks
```

## 14. Practical Notes for a New Reader

If you are reading the project for the first time, start with:

1. `README.md` for commands and file locations.
2. `docs/clinical_her2_high_trust_tile128_results.md` for the current main result.
3. `docs/clinical_her2_trustworthy_slide_list.md` for the trustworthy slide list.
4. `docs/tcga_her2_label_quality_assessment.md` for label assumptions and caveats.
5. `docs/clinical_her2_high_trust_tile128_gigatime_data_cleanup.md` for tile-cleanup sensitivity analyses.
6. `docs/clinical_her2_high_trust_tile128_cleaned_classifier_comparison.md` for the current classifier comparison.
7. This document for the conceptual explanation.
8. `docs/advisor_brief.md` for the short advisor-facing summary.

The most important caution is that GigaTIME outputs are predicted virtual mIF research features. They are not real multiplex immunofluorescence measurements and should be validated before making biological or clinical claims.
