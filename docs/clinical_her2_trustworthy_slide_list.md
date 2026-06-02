# Trustworthy TCGA-BRCA HER2 Slide List

Status: Generated QC/trustworthiness list for the 61/61/61 laptop-balanced TCGA-BRCA HER2 cohort.

## Bottom Line

- Slides checked: 183
- High label+slide trust: 171
- Review before primary analysis: 9
- Exclude from primary analysis: 3
- Output CSV: `data/tcga_brca/clinical_her2_laptop_balanced61_trustworthy_slides.csv`
- Tracked CSV for GitHub: `docs/assets/clinical_her2_trustworthy_slide_list/trustworthy_slides.csv`
- High-trust-only CSV: `data/tcga_brca/clinical_her2_laptop_balanced61_high_trust_slides.csv`
- Tracked high-trust-only CSV for GitHub: `docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv`

High label+slide trust means: known clinical HER2 group, direct HER2 label rule, no flagged IHC/ISH discordance, female patient, primary-tumor diagnostic SVS metadata, local file present, exact file size match, and readable by OpenSlide when that check is available.

This does not mean the slide has pathologist-confirmed tumor-rich sampled tiles. For that, we still need GigaTIME/tile-level QC or human review.

## What We Borrowed From Guardia et al.

The paper the advisor pointed to is Guardia et al., Genome Research 2025, PMID 40664477: `Alternative splicing generates HER2 isoform diversity underlying antibody-drug conjugate resistance in breast cancer`.

Guardia et al. profiled HER2 isoforms in 561 TCGA-BRCA primary breast tumors from female patients and cell-line models. For the TCGA part, they selected and stratified samples using technical considerations such as library preparation and mapped reads, excluded male samples, and grouped tumors by hormone-receptor status plus HER2-high, HER2-low, or HER2-zero status from IHC and/or FISH.

We use this paper as a TCGA cohort-cleanup and HER2 biology reference, not as validation that GigaTIME virtual mIF channels are real measured mIF.

For our H&E/GigaTIME project, the practical translation is: use only traceable direct HER2 labels for primary analysis, exclude or review IHC/ISH discordance, keep female primary-tumor diagnostic slides, stratify or adjust by ER and PR, separate HER2-low and HER2-zero detail subgroups, check slide integrity and tissue quality, and treat RNA/isoform evidence as validation rather than as a substitute for clinical IHC/ISH.

## Group Counts

| HER2 group | Value | Slides |
|---|---|---:|
| HER2-low | high_label_and_slide_trust | 57 |
| HER2-low | review_before_primary_analysis | 4 |
| HER2-positive | exclude_from_primary_analysis | 3 |
| HER2-positive | high_label_and_slide_trust | 53 |
| HER2-positive | review_before_primary_analysis | 5 |
| HER2-zero | high_label_and_slide_trust | 61 |

## HER2 Detail Subgroups

| HER2 group | Value | Slides |
|---|---|---:|
| HER2-low | HER2-low_IHC1_ISH-negative | 8 |
| HER2-low | HER2-low_IHC1_ISH-not-evaluated | 30 |
| HER2-low | HER2-low_IHC2_ISH-negative | 23 |
| HER2-positive | HER2-positive_IHC1_ISH-positive_discordant | 1 |
| HER2-positive | HER2-positive_IHC2_ISH-positive | 15 |
| HER2-positive | HER2-positive_IHC3 | 31 |
| HER2-positive | HER2-positive_IHC3_ISH-negative_discordant | 4 |
| HER2-positive | HER2-positive_ISH-positive_IHC-missing | 10 |
| HER2-zero | HER2-zero_IHC0_ISH-negative | 18 |
| HER2-zero | HER2-zero_IHC0_ISH-not-evaluated | 43 |

## Selection Bias Snapshot

| HER2 group | Slides | Female | ER positive | PR positive | ERBB2 TPM median | Slide size median MB |
|---|---:|---:|---:|---:|---:|---:|
| HER2-positive | 61 | 58 | 49 | 33 | 550 | 143 |
| HER2-low | 61 | 61 | 49 | 45 | 83.4 | 101 |
| HER2-zero | 61 | 61 | 43 | 37 | 62.7 | 234 |

These ER/PR and ERBB2 RNA differences are important because a classifier could accidentally learn hormone-receptor or molecular-subtype context instead of HER2-low versus HER2-zero biology.

## Clinical Context Checks

Patient sex/gender:

| HER2 group | Value | Slides |
|---|---|---:|
| HER2-low | FEMALE | 61 |
| HER2-positive | FEMALE | 58 |
| HER2-positive | MALE | 3 |
| HER2-zero | FEMALE | 61 |

Histology:

| HER2 group | Value | Slides |
|---|---|---:|
| HER2-low | Infiltrating Ductal Carcinoma | 37 |
| HER2-low | Infiltrating Lobular Carcinoma | 19 |
| HER2-low | Mixed Histology (please specify) | 1 |
| HER2-low | Other, specify | 4 |
| HER2-positive | Infiltrating Ductal Carcinoma | 46 |
| HER2-positive | Infiltrating Lobular Carcinoma | 13 |
| HER2-positive | Mucinous Carcinoma | 1 |
| HER2-positive | Other, specify | 1 |
| HER2-zero | Infiltrating Ductal Carcinoma | 49 |
| HER2-zero | Infiltrating Lobular Carcinoma | 7 |
| HER2-zero | Other, specify | 5 |

Pathologic stage:

| HER2 group | Value | Slides |
|---|---|---:|
| HER2-low | Stage I | 3 |
| HER2-low | Stage IA | 6 |
| HER2-low | Stage IB | 1 |
| HER2-low | Stage IIA | 24 |
| HER2-low | Stage IIB | 13 |
| HER2-low | Stage III | 2 |
| HER2-low | Stage IIIA | 8 |
| HER2-low | Stage IIIC | 1 |
| HER2-low | Stage IV | 3 |
| HER2-positive | Stage I | 2 |
| HER2-positive | Stage IA | 4 |
| HER2-positive | Stage IIA | 24 |
| HER2-positive | Stage IIB | 16 |
| HER2-positive | Stage IIIA | 7 |
| HER2-positive | Stage IIIB | 1 |
| HER2-positive | Stage IIIC | 6 |
| HER2-positive | Stage IV | 1 |
| HER2-zero | Stage I | 3 |
| HER2-zero | Stage IA | 7 |
| HER2-zero | Stage IIA | 22 |
| HER2-zero | Stage IIB | 16 |
| HER2-zero | Stage IIIA | 6 |
| HER2-zero | Stage IIIB | 3 |
| HER2-zero | Stage IIIC | 2 |
| HER2-zero | Stage IV | 1 |
| HER2-zero | Stage X | 1 |

## Available GigaTIME Tissue QC

- Slides with existing GigaTIME tissue QC: 60
- Processed slides passing strict tissue QC: 50

Strict processed tissue QC currently means at least 200 sampled tiles, mean tissue fraction at least 0.70, and QC-cellular retained fraction at least 0.50.

## High-Trust Slides By Group

The full machine-readable list is in the CSV above. The table below lists the high-trust slide IDs selected for primary analysis.

| HER2 group | Case | Detail subgroup | Slide ID | ER | PR | ERBB2 TPM | Processed tissue QC |
|---|---|---|---|---|---|---:|---|
| HER2-low | TCGA-5L-AAT0 | HER2-low_IHC1_ISH-not-evaluated | TCGA-5L-AAT0-01A-01-TS1.E7894281-DAB5-4603-B42A-A4C71FBE3D37 | Positive | Positive | 156.7722 | strict_pass |
| HER2-low | TCGA-5T-A9QA | HER2-low_IHC2_ISH-negative | TCGA-5T-A9QA-01A-01-TSA.2224E235-416C-4426-BAE4-55961186FB4E | Positive | Negative | 216.2315 | strict_pass |
| HER2-low | TCGA-A1-A0SJ | HER2-low_IHC2_ISH-negative | TCGA-A1-A0SJ-01A-01-BSA.f9a1ae6c-0cdf-4fdb-b397-17a1e2fc8b74 | Positive | Positive | 157.4375 | moderate_pass |
| HER2-low | TCGA-A2-A04Q | HER2-low_IHC2_ISH-negative | TCGA-A2-A04Q-01A-02-TSB.abf5f4de-662d-4479-a556-2cf0b59d4ee3 | Negative | Negative | 35.8567 | strict_pass |
| HER2-low | TCGA-A2-A04T | HER2-low_IHC2_ISH-negative | TCGA-A2-A04T-01A-02-TSB.93063466-272a-42e9-b4cb-ccf262b7e4c4 | Negative | Negative | 74.48 | strict_pass |
| HER2-low | TCGA-A2-A0CT | HER2-low_IHC2_ISH-negative | TCGA-A2-A0CT-01A-02-BSB.0aa72ee2-8d99-4e03-8f40-7a55d0acb7ad | Positive | Negative | 257.8515 | strict_pass |
| HER2-low | TCGA-A2-A0EN | HER2-low_IHC2_ISH-negative | TCGA-A2-A0EN-01A-01-TSA.d57caa16-813c-41c9-9b9a-3e906e576bb0 | Positive | Positive | 92.4114 | strict_pass |
| HER2-low | TCGA-A2-A0ES | HER2-low_IHC1_ISH-negative | TCGA-A2-A0ES-01A-01-TSA.9e47cc5c-6fcd-4515-ac20-1d81d3c58100 | Positive | Positive | 151.7257 | strict_pass |
| HER2-low | TCGA-A2-A0SV | HER2-low_IHC2_ISH-negative | TCGA-A2-A0SV-01A-01-BSA.94a1dce7-dd72-40a0-8e19-a688d46281f5 | Positive | Positive | 110.2184 | strict_pass |
| HER2-low | TCGA-A2-A0T6 | HER2-low_IHC1_ISH-not-evaluated | TCGA-A2-A0T6-01A-01-BSA.b82fe5ae-dafc-45d6-9410-b6ef20725739 | Positive | Positive | 252.815 | strict_pass |
| HER2-low | TCGA-A7-A0DB | HER2-low_IHC1_ISH-not-evaluated | TCGA-A7-A0DB-01C-02-BS2.97C8D20F-C401-43EA-8374-E1ABA485B0EB | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-A7-A0DC | HER2-low_IHC1_ISH-not-evaluated | TCGA-A7-A0DC-01B-04-BS4.7F68E084-12C4-4591-AD2C-8DADE9DF09F1 | Positive | Negative | nan | not_processed_yet |
| HER2-low | TCGA-A7-A13D | HER2-low_IHC2_ISH-negative | TCGA-A7-A13D-01B-04-BS4.CEB8AD38-4ED8-4D93-B532-945AB128C9DA | Negative | Positive | nan | not_processed_yet |
| HER2-low | TCGA-A7-A13E | HER2-low_IHC2_ISH-negative | TCGA-A7-A13E-01B-03-BS3.62D23E20-E3E7-46B0-BAE1-12204C0F04B7 | Positive | Negative | nan | not_processed_yet |
| HER2-low | TCGA-A7-A13G | HER2-low_IHC1_ISH-not-evaluated | TCGA-A7-A13G-01B-04-BS4.AB30BD66-3064-4733-AA25-C9E13F861211 | Positive | Positive | 54.1363 | strict_pass |
| HER2-low | TCGA-A7-A26E | HER2-low_IHC1_ISH-not-evaluated | TCGA-A7-A26E-01A-01-TSA.c3b5e8a1-130a-41b4-b05f-06f86abd3e21 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-A7-A26F | HER2-low_IHC2_ISH-negative | TCGA-A7-A26F-01B-04-BS4.5818E4B9-7CBD-4241-AD3A-22061417AB30 | Negative | Negative | 32.4939 | strict_pass |
| HER2-low | TCGA-A7-A26I | HER2-low_IHC2_ISH-negative | TCGA-A7-A26I-01B-06-BS6.E89976C5-6194-4D3A-82D3-1A8D1C3D88EF | Negative | Negative | 52.735 | strict_pass |
| HER2-low | TCGA-A7-A26J | HER2-low_IHC1_ISH-not-evaluated | TCGA-A7-A26J-01B-02-BS2.2BDFB544-F62A-402C-9D97-DE2B6766DEDC | Positive | Positive | 70.9602 | moderate_pass |
| HER2-low | TCGA-A7-A6VX | HER2-low_IHC1_ISH-not-evaluated | TCGA-A7-A6VX-01Z-00-DX2.9EE94B59-6A2C-4507-AA4F-DC6402F2B74F | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-A8-A07F | HER2-low_IHC1_ISH-negative | TCGA-A8-A07F-01A-01-TS1.b4547802-0778-49f5-8465-4e3de9ffbf41 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-A8-A08O | HER2-low_IHC1_ISH-negative | TCGA-A8-A08O-01A-02-BS2.6b4babf5-86ae-4b21-ba69-9ea213e76c0d | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-AC-A3OD | HER2-low_IHC1_ISH-not-evaluated | TCGA-AC-A3OD-01A-01-TS1.4915B413-202C-4D3E-A360-E1142F365438 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-AC-A3QP | HER2-low_IHC2_ISH-negative | TCGA-AC-A3QP-01A-01-TS1.51D476D1-143B-4530-A00B-B7764B7DAA0D | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-AO-A0J7 | HER2-low_IHC1_ISH-not-evaluated | TCGA-AO-A0J7-01A-01-BSA.24323379-2daa-45cf-b3e8-f8e4da05a90a | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-AO-A0JG | HER2-low_IHC1_ISH-not-evaluated | TCGA-AO-A0JG-01A-03-BSC.6b4d1364-f805-4583-8aac-7382e95a2cc7 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-AR-A0U2 | HER2-low_IHC1_ISH-not-evaluated | TCGA-AR-A0U2-01A-01-BSA.1bd5ddeb-f22c-4147-b905-0f1316e5903a | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-AR-A24N | HER2-low_IHC2_ISH-negative | TCGA-AR-A24N-01A-01-TSA.b617b299-5c2d-45c3-b018-6b8cf0e84811 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-AR-A2LH | HER2-low_IHC2_ISH-negative | TCGA-AR-A2LH-01A-03-TSC.F95ACB2C-4C81-4104-A40D-D4DD52A1CC2F | Negative | Negative | nan | not_processed_yet |
| HER2-low | TCGA-AR-A2LK | HER2-low_IHC1_ISH-negative | TCGA-AR-A2LK-01A-01-TSA.53C36847-0443-432D-97FF-02949CCBA21B | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-AR-A2LO | HER2-low_IHC2_ISH-negative | TCGA-AR-A2LO-01A-03-TSC.62917FFC-65A4-46A2-AED6-1E9B368023E2 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-B6-A3ZX | HER2-low_IHC1_ISH-not-evaluated | TCGA-B6-A3ZX-01A-01-TS1.68F70986-83F0-4B65-AFF5-496B0C60270F | Negative | Negative | 45.4756 | strict_pass |
| HER2-low | TCGA-B6-A400 | HER2-low_IHC1_ISH-not-evaluated | TCGA-B6-A400-01A-01-TS1.825AC7BC-8A50-49C2-8346-8F39A55D117B | Negative | Negative | nan | not_processed_yet |
| HER2-low | TCGA-B6-A402 | HER2-low_IHC1_ISH-not-evaluated | TCGA-B6-A402-01A-01-TS1.45F6B344-FEE5-4F83-AB34-4F1BA99E445D | Negative | Negative | nan | not_processed_yet |
| HER2-low | TCGA-B6-A409 | HER2-low_IHC1_ISH-negative | TCGA-B6-A409-01A-01-TS1.02A10D0A-2B31-4165-94E3-4E6CA784A092 | Negative | Negative | 53.3979 | moderate_pass |
| HER2-low | TCGA-B6-A40B | HER2-low_IHC1_ISH-not-evaluated | TCGA-B6-A40B-01A-01-TS1.5B2B3C90-F8E8-4853-B1F8-1CEE52253550 | Positive | Positive | 193.373 | strict_pass |
| HER2-low | TCGA-B6-A40C | HER2-low_IHC1_ISH-not-evaluated | TCGA-B6-A40C-01A-01-TS1.8B901E88-9348-4988-9F69-C130D7D8E78E | Positive | Positive | 140.743 | strict_pass |
| HER2-low | TCGA-C8-A1HM | HER2-low_IHC1_ISH-not-evaluated | TCGA-C8-A1HM-01A-01-TSA.8a4cc75b-80ec-45e3-8b1a-a5e66aea10b5 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-D8-A1XR | HER2-low_IHC1_ISH-not-evaluated | TCGA-D8-A1XR-01A-01-TS1.087d05ee-9a03-4433-9060-98a80832867a | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-D8-A27G | HER2-low_IHC2_ISH-negative | TCGA-D8-A27G-01A-01-TS1.59dc5ba9-df66-4674-82ba-c791b36adfe9 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-D8-A27K | HER2-low_IHC1_ISH-not-evaluated | TCGA-D8-A27K-01A-01-TS1.b8d6df65-efe5-4317-b903-b93ec05aaa5f | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-D8-A27P | HER2-low_IHC1_ISH-not-evaluated | TCGA-D8-A27P-01A-01-TSA.df327c6c-1fde-4c29-ad33-de4b8420170f | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-D8-A73U | HER2-low_IHC1_ISH-negative | TCGA-D8-A73U-01A-01-TS1.112044FF-39FE-4588-9517-31F7821746D0 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-E2-A106 | HER2-low_IHC2_ISH-negative | TCGA-E2-A106-01A-01-TS1.c5fe7052-cb83-49bd-8c1d-a82067606921 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-E2-A108 | HER2-low_IHC1_ISH-not-evaluated | TCGA-E2-A108-01A-01-BSA.285fd0a0-82be-4775-bce9-85cd5497e382 | Positive | Positive | 48.3077 | moderate_pass |
| HER2-low | TCGA-E2-A10C | HER2-low_IHC1_ISH-not-evaluated | TCGA-E2-A10C-01A-02-TSB.4d05f801-0494-42ad-9365-9b59bf957f4a | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-E2-A10F | HER2-low_IHC2_ISH-negative | TCGA-E2-A10F-01A-01-TS1.2e848f27-cde2-4da2-baaf-f4de637f891c | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-E2-A14T | HER2-low_IHC1_ISH-negative | TCGA-E2-A14T-01A-01-TSA.3e947c63-817a-4b85-8065-30adf00015b5 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-E2-A1LA | HER2-low_IHC2_ISH-negative | TCGA-E2-A1LA-01A-01-TS1.f300df8b-2423-47e1-8a95-24048aaeac80 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-E2-A1LK | HER2-low_IHC2_ISH-negative | TCGA-E2-A1LK-01A-02-TSB.5f4a6e3e-8c5c-43b3-9053-44a4d0f9f954 | Negative | Negative | nan | not_processed_yet |
| HER2-low | TCGA-EW-A1PA | HER2-low_IHC1_ISH-not-evaluated | TCGA-EW-A1PA-01A-01-TSA.b6516a72-91fb-4dd8-8da0-2a9989a48e14 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-GM-A5PX | HER2-low_IHC1_ISH-negative | TCGA-GM-A5PX-01A-01-TS1.B6353D86-B0CD-464C-B976-752100D5CAE8 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-LL-A442 | HER2-low_IHC2_ISH-negative | TCGA-LL-A442-01A-01-TSA.EB28ADD6-530F-423A-B6C4-998490972A34 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-LL-A5YO | HER2-low_IHC1_ISH-not-evaluated | TCGA-LL-A5YO-01A-02-TS2.5819C502-EEA4-4D98-B4EC-E184B374ECFB | Negative | Negative | nan | not_processed_yet |
| HER2-low | TCGA-LL-A8F5 | HER2-low_IHC1_ISH-not-evaluated | TCGA-LL-A8F5-01A-01-TS1.B5357CA4-8C24-4F27-9054-06B801C3D8DF | Positive | Negative | nan | not_processed_yet |
| HER2-low | TCGA-S3-A6ZF | HER2-low_IHC2_ISH-negative | TCGA-S3-A6ZF-01A-03-TS3.2DE8CF7F-0109-42C3-8C05-4D708BBE19F6 | Positive | Positive | nan | not_processed_yet |
| HER2-low | TCGA-WT-AB44 | HER2-low_IHC1_ISH-not-evaluated | TCGA-WT-AB44-01A-01-TS1.B6C0EEDB-E5B9-4B0D-8599-23879A0419EB | Positive | Positive | 65.4996 | strict_pass |
| HER2-positive | TCGA-A2-A04W | HER2-positive_ISH-positive_IHC-missing | TCGA-A2-A04W-01A-03-TSC.83981bb3-9dfd-4632-91c0-00068ad599f3 | Negative | Negative | 1133.2294 | strict_pass |
| HER2-positive | TCGA-A2-A04X | HER2-positive_IHC3 | TCGA-A2-A04X-01A-02-BS2.e8b8e62b-18f3-47f1-8b53-94659913afda | Positive | Positive | 877.1314 | weak_or_review |
| HER2-positive | TCGA-A2-A0CX | HER2-positive_IHC3 | TCGA-A2-A0CX-01A-02-TSB.eae9814e-e25c-4dfe-9564-e348c74699ba | Positive | Negative | 1465.7897 | strict_pass |
| HER2-positive | TCGA-A2-A0D1 | HER2-positive_IHC3 | TCGA-A2-A0D1-01A-01-MSA.3bc9e5b7-75ab-4eea-8c50-7f67dd11b2c7 | Negative | Negative | 1122.8745 | strict_pass |
| HER2-positive | TCGA-A2-A0EY | HER2-positive_IHC3 | TCGA-A2-A0EY-01A-01-TSA.b3cf2e58-b037-4496-b85e-6fefe2d42e01 | Positive | Negative | 1680.5636 | strict_pass |
| HER2-positive | TCGA-A2-A0SY | HER2-positive_ISH-positive_IHC-missing | TCGA-A2-A0SY-01A-03-BSC.111f2ea3-8315-412c-807d-5582f160cb42 | Positive | Positive | 420.0192 | moderate_pass |
| HER2-positive | TCGA-A2-A0T1 | HER2-positive_IHC3 | TCGA-A2-A0T1-01A-02-TSB.63e971bf-39ec-41c9-986e-a1937437dea8 | Negative | Negative | 1699.1491 | strict_pass |
| HER2-positive | TCGA-A2-A1G1 | HER2-positive_IHC2_ISH-positive | TCGA-A2-A1G1-01A-02-TSB.084536c6-00fe-428c-b744-47d391f0a0ef | Negative | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-A2-A3XV | HER2-positive_IHC2_ISH-positive | TCGA-A2-A3XV-01A-02-TSB.FF8434E6-B703-43FE-AC0A-AE53131F1EC6 | Positive | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-A7-A2KD | HER2-positive_IHC3 | TCGA-A7-A2KD-01A-03-TSC.CC1E46B0-7920-4F91-AE61-A9A78BF724DB | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-A7-A4SF | HER2-positive_IHC3 | TCGA-A7-A4SF-01A-01-TSA.A9D65EF2-4CD0-41BC-98D0-FF6495369BF3 | Positive | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-A8-A06U | HER2-positive_IHC2_ISH-positive | TCGA-A8-A06U-01A-01-TS1.DFD8B445-C7A5-4247-A39F-222591C6D7E2 | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-A8-A075 | HER2-positive_IHC2_ISH-positive | TCGA-A8-A075-01A-01-TS1.3658536c-c029-4685-979f-11e914e53e4e | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-A8-A07P | HER2-positive_IHC2_ISH-positive | TCGA-A8-A07P-01A-01-BS1.4bb95702-b8db-403d-9f30-1519eb9c7b11 | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-A8-A08B | HER2-positive_IHC3 | TCGA-A8-A08B-01A-01-BS1.cd256d6e-99c5-4b32-8bbc-7c1811004b8b | Positive | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-A8-A08C | HER2-positive_IHC2_ISH-positive | TCGA-A8-A08C-01A-01-TS1.fc1246d4-246a-482f-a819-b80f25b95f87 | Positive | Positive | 157.259 | strict_pass |
| HER2-positive | TCGA-A8-A08G | HER2-positive_IHC3 | TCGA-A8-A08G-01A-01-TS1.00c53463-683f-443f-b693-5673b3b564ff | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-A8-A08H | HER2-positive_IHC3 | TCGA-A8-A08H-01A-02-BS2.26653396-4d2c-4614-9cb1-244e7d4fdfae | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-A8-A08T | HER2-positive_IHC2_ISH-positive | TCGA-A8-A08T-01A-02-BS2.f401c61c-3713-4bb9-a709-71423c478c8b | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-A8-A09G | HER2-positive_IHC2_ISH-positive | TCGA-A8-A09G-01A-01-BS1.d16e6e1e-268a-4973-97c5-8aa9cdef622b | Positive | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-A8-A09N | HER2-positive_IHC3 | TCGA-A8-A09N-01A-01-BS1.48ca5ab6-429b-420a-9a89-7605aee874ab | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-A8-A0A7 | HER2-positive_IHC3 | TCGA-A8-A0A7-01A-01-BS1.4a4d1511-7d0f-4482-bd6c-7a1ebbfe2ce6 | Negative | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-A8-A0AB | HER2-positive_IHC2_ISH-positive | TCGA-A8-A0AB-01Z-00-DX1.103ED338-A0F9-403B-A10C-49840BD60EB8 | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-AC-A23C | HER2-positive_ISH-positive_IHC-missing | TCGA-AC-A23C-01Z-00-DX1.0E67C785-83D3-49AF-B600-FB5B909AE6ED | Positive | Positive | 680.2476 | strict_pass |
| HER2-positive | TCGA-AC-A23G | HER2-positive_ISH-positive_IHC-missing | TCGA-AC-A23G-01Z-00-DX1.2F0326F7-6B77-4B3F-B4FA-59ADB785AA07 | Positive | Positive | 118.5879 | strict_pass |
| HER2-positive | TCGA-AC-A23H | HER2-positive_ISH-positive_IHC-missing | TCGA-AC-A23H-01Z-00-DX1.8E0AE339-1047-4CA5-BFC5-37A3B10FD8B5 | Positive | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-AC-A2FB | HER2-positive_IHC3 | TCGA-AC-A2FB-01Z-00-DX1.A4D93E32-BBD7-45E4-8ACF-3724B059ECBC | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-AC-A3QQ | HER2-positive_IHC3 | TCGA-AC-A3QQ-01B-06-BS6.C20D50B8-6EC0-4BF4-BE23-D5664C18B702 | Positive | Positive | 264.2144 | weak_or_review |
| HER2-positive | TCGA-AC-A3W5 | HER2-positive_ISH-positive_IHC-missing | TCGA-AC-A3W5-01Z-00-DX1.522F702C-776E-45F1-84D3-1648DF04137C | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-AC-A3W6 | HER2-positive_IHC3 | TCGA-AC-A3W6-01Z-00-DX1.88CC534C-F032-4E5D-9CC4-4BB50AA46880 | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-AN-A04C | HER2-positive_IHC3 | TCGA-AN-A04C-01A-01-BS1.4f8165fd-0791-4fa6-97ea-66b2b3eb0730 | Negative | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-AN-A0AJ | HER2-positive_IHC3 | TCGA-AN-A0AJ-01A-01-TSA.c21f9dd1-b47c-406f-9ad6-6969bd6f3246 | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-AO-A0JM | HER2-positive_IHC3 | TCGA-AO-A0JM-01A-02-TSB.d19f822d-0237-41cb-be13-4d3023dcb507 | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-AQ-A04L | HER2-positive_IHC3 | TCGA-AQ-A04L-01B-02-BSB.229b081b-f757-4082-85fb-303bf19e611c | Positive | Negative | 100.2928 | strict_pass |
| HER2-positive | TCGA-AQ-A0Y5 | HER2-positive_IHC2_ISH-positive | TCGA-AQ-A0Y5-01A-01-BSA.0402a8fd-a121-4a86-927e-073a2e4a7e7e | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-AQ-A1H2 | HER2-positive_IHC2_ISH-positive | TCGA-AQ-A1H2-01A-01-TSA.529445fe-410c-43db-9fe3-e5ac13387eec | Positive | Positive | 191.2381 | strict_pass |
| HER2-positive | TCGA-AR-A0TX | HER2-positive_IHC3 | TCGA-AR-A0TX-01A-01-BSA.2d59aca1-69fd-4105-830a-a687eddd27c9 | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-AR-A24U | HER2-positive_IHC3 | TCGA-AR-A24U-01A-01-TSA.e0a5c6e6-f748-47c5-8ede-f526fef07b4f | Negative | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-AR-A255 | HER2-positive_IHC3 | TCGA-AR-A255-01A-01-TSA.85414c8e-dc54-4ae2-80dd-c58db439c5c7 | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-AR-A2LJ | HER2-positive_IHC3 | TCGA-AR-A2LJ-01A-01-TSA.D42CBC19-1BB0-44B0-A70D-4F872BA53E2C | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-C8-A12P | HER2-positive_IHC3 | TCGA-C8-A12P-01A-01-TSA.a635ee87-3392-4ce7-93a3-87fc5680e4e9 | Negative | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-C8-A3M8 | HER2-positive_IHC3 | TCGA-C8-A3M8-01A-01-TSA.F27E2918-AAA7-41F0-90D3-56322773E678 | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-D8-A1J9 | HER2-positive_IHC3 | TCGA-D8-A1J9-01A-01-TSA.a7772c0a-357e-4d4c-b74c-1c14d13ab744 | Positive | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-D8-A27N | HER2-positive_IHC3 | TCGA-D8-A27N-01A-01-TS1.f4371b61-6f8f-48c6-9493-d9c96402634f | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-D8-A27W | HER2-positive_IHC2_ISH-positive | TCGA-D8-A27W-01A-01-TSA.73251955-9d38-4c32-bce3-567a50db25a4 | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-E2-A1LE | HER2-positive_IHC3 | TCGA-E2-A1LE-01A-01-TSA.3D0F1C87-1DDC-4BF0-A586-401299A7286D | Negative | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-EW-A1OZ | HER2-positive_IHC2_ISH-positive | TCGA-EW-A1OZ-01A-01-TSA.5b179149-f5d2-4b79-8988-5c2fc2268cdc | Positive | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-EW-A2FR | HER2-positive_IHC2_ISH-positive | TCGA-EW-A2FR-01A-01-TSA.516FAE11-477F-423F-82CB-9041FA0D8111 | Negative | Negative | nan | not_processed_yet |
| HER2-positive | TCGA-LL-A7T0 | HER2-positive_IHC3 | TCGA-LL-A7T0-01A-03-TS3.07761D74-57D6-4859-AD0B-BD2F8C68D905 | Positive | Positive | nan | not_processed_yet |
| HER2-positive | TCGA-OL-A5RY | HER2-positive_ISH-positive_IHC-missing | TCGA-OL-A5RY-01Z-00-DX1.AE4E9D74-FC1C-4C1E-AE6D-5DF38899BBA6 | Positive | Negative | 1875.9017 | strict_pass |
| HER2-positive | TCGA-OL-A5RZ | HER2-positive_ISH-positive_IHC-missing | TCGA-OL-A5RZ-01Z-00-DX1.6394C05E-1C34-4F4B-8859-F5E961E7EFF9 | Positive | Negative | 3532.8051 | weak_or_review |
| HER2-positive | TCGA-OL-A5S0 | HER2-positive_ISH-positive_IHC-missing | TCGA-OL-A5S0-01Z-00-DX1.49A7AC9D-C186-406C-BA67-2D73DE82E13B | Positive | Negative | 100.044 | strict_pass |
| HER2-positive | TCGA-S3-AA14 | HER2-positive_IHC3 | TCGA-S3-AA14-01A-01-TS1.8A693383-EEAF-4F83-9FA1-8F59AB8B619B | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-A1-A0SK | HER2-zero_IHC0_ISH-not-evaluated | TCGA-A1-A0SK-01A-01-BSA.29fd9144-1886-47ad-b075-7800cb57c9f5 | Negative | Negative | 6.7185 | strict_pass |
| HER2-zero | TCGA-A1-A0SP | HER2-zero_IHC0_ISH-not-evaluated | TCGA-A1-A0SP-01A-01-TS1.0cda5515-8bec-4df9-a50e-f4a649aba4bf | Negative | Negative | 51.0204 | strict_pass |
| HER2-zero | TCGA-A2-A0CM | HER2-zero_IHC0_ISH-not-evaluated | TCGA-A2-A0CM-01A-03-BSC.c240a44d-0583-496e-bdef-5cdb5e0ee167 | Negative | Negative | 42.1815 | strict_pass |
| HER2-zero | TCGA-A2-A0D0 | HER2-zero_IHC0_ISH-negative | TCGA-A2-A0D0-01A-01-TSA.231084a7-765a-481e-9e40-2c9b131b38ea | Negative | Negative | 38.3901 | strict_pass |
| HER2-zero | TCGA-A2-A0D2 | HER2-zero_IHC0_ISH-negative | TCGA-A2-A0D2-01A-02-BSB.8b0210e1-3e64-4e14-8a7f-785bc507f7d8 | Negative | Negative | 61.5892 | strict_pass |
| HER2-zero | TCGA-A2-A0EU | HER2-zero_IHC0_ISH-negative | TCGA-A2-A0EU-01A-02-TSB.648d1299-d273-4ae2-8926-cfad39ff735a | Positive | Positive | 117.0857 | strict_pass |
| HER2-zero | TCGA-A2-A0EV | HER2-zero_IHC0_ISH-negative | TCGA-A2-A0EV-01A-01-TSA.49658afd-55bd-4adb-a3c2-06a4e633efc5 | Positive | Positive | 226.5908 | strict_pass |
| HER2-zero | TCGA-A2-A0EW | HER2-zero_IHC0_ISH-negative | TCGA-A2-A0EW-01A-02-TSB.6ebf62a3-15ec-43dd-a5de-36f0a2effb41 | Positive | Positive | 141.1141 | strict_pass |
| HER2-zero | TCGA-A2-A0T0 | HER2-zero_IHC0_ISH-not-evaluated | TCGA-A2-A0T0-01A-02-TSB.802f97da-e482-4c7a-aa5e-a692871ad1ea | Negative | Negative | 54.4828 | strict_pass |
| HER2-zero | TCGA-A2-A0T2 | HER2-zero_IHC0_ISH-not-evaluated | TCGA-A2-A0T2-01A-01-BSA.3613219b-687b-419f-ac8f-33b668b72446 | Negative | Negative | 40.535 | strict_pass |
| HER2-zero | TCGA-A2-A0YE | HER2-zero_IHC0_ISH-not-evaluated | TCGA-A2-A0YE-01A-01-BS1.7203b6fc-3a45-41e1-9104-14255fe4eb6e | Negative | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-A2-A0YF | HER2-zero_IHC0_ISH-not-evaluated | TCGA-A2-A0YF-01A-02-TSB.65b05e09-46e2-47f8-8ef6-ccc5914752c5 | Positive | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-A2-A0YH | HER2-zero_IHC0_ISH-not-evaluated | TCGA-A2-A0YH-01A-01-TSA.77998a85-49c4-4bdf-8635-ec9516aad084 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-A2-A0YJ | HER2-zero_IHC0_ISH-not-evaluated | TCGA-A2-A0YJ-01A-01-BSA.d8a11cff-f5fd-4d70-8655-c8db88f7147e | Positive | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-A8-A06N | HER2-zero_IHC0_ISH-negative | TCGA-A8-A06N-01A-01-BS1.12bf6cf9-767e-469a-ab07-accd6b3b134e | Positive | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-A8-A07L | HER2-zero_IHC0_ISH-negative | TCGA-A8-A07L-01A-01-TS1.ef6b6a88-8ca0-4d4f-90fe-4996e6a12e0f | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-A8-A086 | HER2-zero_IHC0_ISH-negative | TCGA-A8-A086-01A-01-BSA.f75d10e3-5a89-4b4f-9b40-6bca77731bb2 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-A8-A09C | HER2-zero_IHC0_ISH-negative | TCGA-A8-A09C-01A-01-TS1.60187ffe-746b-4a98-a6b2-06a6f884f04a | Positive | Positive | 80.3384 | strict_pass |
| HER2-zero | TCGA-A8-A09K | HER2-zero_IHC0_ISH-negative | TCGA-A8-A09K-01A-01-TS1.8443b0b9-57e1-4bd2-957d-033a523b4caf | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-A8-A09V | HER2-zero_IHC0_ISH-negative | TCGA-A8-A09V-01A-01-BS1.2d168522-0774-4867-a685-672e2fa1d0d0 | Positive | Positive | 148.1888 | strict_pass |
| HER2-zero | TCGA-A8-A09W | HER2-zero_IHC0_ISH-negative | TCGA-A8-A09W-01A-01-BS1.35bdda00-d682-4d57-8ce9-d1e7234b0db8 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-A8-A0A2 | HER2-zero_IHC0_ISH-negative | TCGA-A8-A0A2-01A-01-TS1.124134ed-c4d2-4a08-937c-1def71bfcc42 | Positive | Positive | 126.1578 | strict_pass |
| HER2-zero | TCGA-AN-A03Y | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AN-A03Y-01A-01-BS1.b7b9c538-7d46-47c8-8071-9b6080f6f330 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-AN-A049 | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AN-A049-01A-01-BS1.01c2a132-3509-42c7-b003-8811f0947bc1 | Positive | Positive | 279.4578 | strict_pass |
| HER2-zero | TCGA-AN-A04A | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AN-A04A-01A-01-BS1.647f0482-49a8-4794-b9c4-5941b14fd1af | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-AN-A04D | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AN-A04D-01A-02-BS2.5c72d68a-9c0e-4551-bc1f-21f91c168acc | Negative | Negative | 63.8607 | strict_pass |
| HER2-zero | TCGA-AN-A0AL | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AN-A0AL-01A-01-BSA.56e46756-aa0f-4c59-a189-ee4c9601b3ea | Negative | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-AN-A0AR | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AN-A0AR-01A-01-TSA.8a0c5187-915e-4608-9eed-9ea98f1775ad | Negative | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-AN-A0AS | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AN-A0AS-01A-01-BSA.1c59c3dd-eebc-421c-ad6f-953c0b924b13 | Positive | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-AN-A0FF | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AN-A0FF-01A-01-BSA.ef124d17-e63f-4c55-a50a-1cdea1b4cbde | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-AN-A0FY | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AN-A0FY-01A-01-BSA.288f78df-2abe-4b01-aa91-ad2088788fd9 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-AO-A03N | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AO-A03N-01B-01-BSA.e6642355-b3e1-4cc4-81ca-129b9ddf9240 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-AO-A03P | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AO-A03P-01A-01-BS1.b5b1b2eb-6a49-4d06-a855-7263aaa78ac1 | Positive | Positive | 101.0977 | strict_pass |
| HER2-zero | TCGA-AO-A03R | HER2-zero_IHC0_ISH-negative | TCGA-AO-A03R-01A-01-BS1.7a857be4-4080-4cf0-88b2-794b20a06a97 | Positive | Positive | 14.5836 | moderate_pass |
| HER2-zero | TCGA-AO-A03T | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AO-A03T-01A-02-MS2.910c2a57-801f-422c-ae31-17dace94b475 | Positive | Positive | 29.2507 | moderate_pass |
| HER2-zero | TCGA-AO-A03U | HER2-zero_IHC0_ISH-negative | TCGA-AO-A03U-01B-02-BSB.dcb167f4-c3ab-4dcc-8f40-41c4ce453847 | Negative | Negative | 61.0963 | strict_pass |
| HER2-zero | TCGA-AO-A03V | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AO-A03V-01A-01-TS1.9b67fe65-977d-4b18-bf14-79378230f289 | Positive | Positive | 109.2495 | strict_pass |
| HER2-zero | TCGA-AO-A0J6 | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AO-A0J6-01A-01-TSA.d62f4c33-ae94-4321-b833-53663942e846 | Negative | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-AO-A0JB | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AO-A0JB-01A-01-BSA.705ff788-9aea-491f-8eb4-d19e21c984ba | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-AO-A0JC | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AO-A0JC-01A-01-BSA.2a81630f-81ee-4c33-ba26-0e825cfa5f72 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-AO-A0JJ | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AO-A0JJ-01A-01-BSA.132f03bb-3b31-485c-b814-0f4bf204440e | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-AO-A124 | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AO-A124-01A-01-TSA.ace4fa54-5ef0-4cb9-85be-8c7680889278 | Negative | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-AO-A125 | HER2-zero_IHC0_ISH-negative | TCGA-AO-A125-01A-01-BS1.e881fb71-5e5e-40ec-8284-3ee696de7239 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-AO-A128 | HER2-zero_IHC0_ISH-negative | TCGA-AO-A128-01A-01-BSA.9019856f-fe50-42af-8588-b73395fc0165 | Negative | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-AO-A129 | HER2-zero_IHC0_ISH-negative | TCGA-AO-A129-01A-02-BSB.a781dca2-08b1-46a7-a7c9-7c2944546e58 | Negative | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-AO-A12B | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AO-A12B-01A-01-TSA.635646b4-dc81-45e2-92b1-8f087abc5240 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-AO-A12E | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AO-A12E-01A-01-BSA.9b55ad59-32de-4238-9958-b54e3a43eb2e | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-AO-A12H | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AO-A12H-01A-01-BSA.eb4b965d-5a91-4c8b-a61d-177b895d2fe5 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-AQ-A04J | HER2-zero_IHC0_ISH-not-evaluated | TCGA-AQ-A04J-01A-01-BS1.83b89953-2edf-4a7c-bb98-c308fadffc3d | Negative | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-BH-A0AY | HER2-zero_IHC0_ISH-not-evaluated | TCGA-BH-A0AY-01A-02-BSB.581eb097-e917-4666-9344-f62b088cb024 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-BH-A0B9 | HER2-zero_IHC0_ISH-not-evaluated | TCGA-BH-A0B9-01A-01-BSA.f1673d45-8491-4821-87e3-62e7cbcb2994 | Negative | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-BH-A0BC | HER2-zero_IHC0_ISH-not-evaluated | TCGA-BH-A0BC-01A-02-TSB.a7afb5f2-2676-428b-95c4-5b8764004820 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-BH-A0BM | HER2-zero_IHC0_ISH-not-evaluated | TCGA-BH-A0BM-01A-01-BSA.0c2b78c5-2579-4a8a-81cb-e2f6cbb4a019 | Positive | Negative | nan | not_processed_yet |
| HER2-zero | TCGA-BH-A0BV | HER2-zero_IHC0_ISH-not-evaluated | TCGA-BH-A0BV-01A-01-BSA.1d27d74a-405e-4a7f-befa-6d5f7b1fa8e9 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-BH-A0DH | HER2-zero_IHC0_ISH-not-evaluated | TCGA-BH-A0DH-01A-01-BSA.8ac4fdc4-cc2c-4c7f-b0b9-e55e4a6d225d | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-BH-A0DK | HER2-zero_IHC0_ISH-not-evaluated | TCGA-BH-A0DK-01A-02-BSB.d2b81f36-c819-44a9-8171-5297fd30cd73 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-BH-A0DP | HER2-zero_IHC0_ISH-not-evaluated | TCGA-BH-A0DP-01A-02-BSB.ad74deb9-5256-4688-9a5c-c4db6be91e67 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-BH-A0DQ | HER2-zero_IHC0_ISH-not-evaluated | TCGA-BH-A0DQ-01A-01-MSA.3e981cdd-f9e6-454a-9746-16e29a072b17 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-BH-A0H0 | HER2-zero_IHC0_ISH-not-evaluated | TCGA-BH-A0H0-01A-01-TSA.70fed29e-bd3c-4cb6-92c4-0c16bd6048d6 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-BH-A0HF | HER2-zero_IHC0_ISH-not-evaluated | TCGA-BH-A0HF-01A-01-BSA.e917388d-c4db-4345-b7de-172bc1894173 | Positive | Positive | nan | not_processed_yet |
| HER2-zero | TCGA-BH-A0HK | HER2-zero_IHC0_ISH-not-evaluated | TCGA-BH-A0HK-01A-01-BSA.9503c4b2-4f84-4b5e-80a1-2e29605606b2 | Positive | Negative | nan | not_processed_yet |

## Review Before Primary Analysis

| HER2 group | Case | Slide ID | Reason |
|---|---|---|---|
| HER2-low | TCGA-AC-A8OS | TCGA-AC-A8OS-01A-01-TS1.31C88D4F-7621-4B24-83CC-8FD69EC72933 | HER2 field flag: low_group_with_positive_receptor_status |
| HER2-low | TCGA-AN-A03X | TCGA-AN-A03X-01A-02-BS2.88b66ddf-fc7b-4b5f-904d-7eafc9654a1a | HER2 field flag: low_group_with_positive_receptor_status |
| HER2-low | TCGA-AN-A0XW | TCGA-AN-A0XW-01Z-00-DX1.811E11E7-FA67-46BB-9BC6-1FD0106B789D | HER2 field flag: low_group_with_positive_receptor_status |
| HER2-low | TCGA-E9-A295 | TCGA-E9-A295-01A-01-TSA.14a39b7e-5d34-43ea-971c-09a8db7694e8 | HER2 field flag: low_group_with_positive_receptor_status |
| HER2-positive | TCGA-A2-A04U | TCGA-A2-A04U-01A-01-BSA.65a9cf7d-8f03-45b6-8751-5a535db11426 | HER2 field flag: low_or_zero_ihc_with_positive_ish |
| HER2-positive | TCGA-A2-A0EQ | TCGA-A2-A0EQ-01A-01-BSA.fa7a323e-2939-4941-a56f-3b7352b16050 | HER2 field flag: ihc3_with_negative_ish |
| HER2-positive | TCGA-A7-A425 | TCGA-A7-A425-01A-01-TSA.D8ECCE2B-4CC3-4459-8C20-9F21DE0F9241 | HER2 field flag: ihc3_with_negative_ish |
| HER2-positive | TCGA-A7-A4SC | TCGA-A7-A4SC-01A-01-TS1.039AAA5C-0866-45EF-9FD2-36F9A6AA58AA | HER2 field flag: ihc3_with_negative_ish |
| HER2-positive | TCGA-LL-A5YL | TCGA-LL-A5YL-01A-01-TS1.FC54A71A-76C6-4C15-BBD2-5BDC859BA356 | HER2 field flag: ihc3_with_negative_ish |

## Excluded From Primary Analysis

| HER2 group | Case | Slide ID | Reason |
|---|---|---|---|
| HER2-positive | TCGA-A1-A0SM | TCGA-A1-A0SM-01A-01-BS1.f31192b0-feeb-4274-bd1e-53b1a036c888 | patient gender is not FEMALE |
| HER2-positive | TCGA-E2-A14W | TCGA-E2-A14W-01A-01-TSA.292e68b2-d49a-46e0-aca5-a587709df55a | patient gender is not FEMALE |
| HER2-positive | TCGA-EW-A1PD | TCGA-EW-A1PD-01A-01-TSA.fc5cdd95-913b-4466-bcbb-6f482bfb0764 | patient gender is not FEMALE |

## Sources Checked

- Guardia et al., PMID 40664477: https://pubmed.ncbi.nlm.nih.gov/40664477/
- Guardia et al., Genome Research abstract: https://genome.cshlp.org/content/35/9/1942.abstract
- Guardia et al., supplemental material: https://genome.cshlp.org/content/35/9/1942/suppl/DC1
- Guardia et al., medRxiv full-text preprint used to verify methods language: https://www.medrxiv.org/content/10.1101/2024.11.25.24317569v1.full-text
- GDC Clinical Data overview: https://docs.gdc.cancer.gov/Encyclopedia/pages/Clinical_Data/
- GDC Data Dictionary overview: https://docs.gdc.cancer.gov/Data_Dictionary/
- TCGA-BRCA Breast Enrollment Form: https://gdc.cancer.gov/system/files/public/file/Breast%20Enrollment%20Form.pdf
- CAP/ASCO HER2 testing guideline page: https://www.cap.org/protocols-and-guidelines/cap-guidelines/current-cap-guidelines/recommendations-for-human-epidermal-growth-factor-2-testing-in-breast-cancer

## Recommended Use

- Primary analysis: use `label_slide_trust == high_label_and_slide_trust`.
- Strict sensitivity analysis: exclude `processed_tissue_qc != strict_pass` after GigaTIME is run on the full 183-slide cohort.
- HER2-low sensitivity analysis: run IHC `1+` and IHC `2+`/ISH-negative separately.
- HER2-zero sensitivity analysis: run IHC `0`/ISH-negative and IHC `0`/ISH-not-evaluated separately.
- Discordance sensitivity analysis: exclude all cases with non-empty `discordance_flags`.
