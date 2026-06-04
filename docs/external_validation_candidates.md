# External Validation Candidate Cohorts

Status: working reference for external (non-TCGA) validation of the HER2-low versus HER2-zero image signal. Compiled 2026-06-04 from a web/literature scout. Numbers and access details should be re-verified against the primary sources before committing to any cohort.

## Why External Data Is Now The Bottleneck

The HER2-low versus HER2-zero signal in this project has been shown, repeatedly and from multiple angles, to be confounded by TCGA acquisition structure (slide size and source site) rather than tumor-cell HER2 biology. As of 2026-06-04 this is supported by clinical/source-site covariate baselines, matched subsets, leave-source-site-out validation, within-source-site restriction, and two independent generic foundation-model embedding controls (H-Optimus-0 and Virchow2) that reproduce the separation and the same source-site collapse. See `docs/clinical_her2_high_trust_tile128_results.md`.

TCGA-internal evidence is therefore exhausted. Pulling more TCGA-BRCA slides does not help: HER2-zero is capped at 61 cases in all of TCGA-BRCA (already fully used), and more data along a confounded axis tightens confidence intervals around a biased estimate. The only data that can change the conclusion is variation independent of HER2 status, i.e. an external cohort, ideally single-scanner / single-institution so the slide-size/source-site confound is removed by construction.

## The Hard Constraint: HER2-Low-vs-Zero Granularity

The exact comparison this project cares about is HER2-low (IHC 1+, or 2+/ISH-negative) versus HER2-zero (IHC 0). This is the difficult attribute to source externally, for two reasons:

1. Most public H&E + HER2 datasets label only binary HER2-positive versus HER2-negative. HER2-low is a post-2019 clinical category (driven by trastuzumab deruxtecan / DESTINY-Breast04), so older public datasets predate it and collapse 0 and 1+ into "negative."
2. Even with the IHC slide in hand, pathologist interobserver agreement on IHC 0 versus 1+ is low. The ground truth on this exact boundary is intrinsically noisy.

Practical implication: a clean external reproduction of the low-versus-zero signal would be a strong positive result; failure to reproduce would not be surprising and would be consistent with the confound interpretation. A cohort is only fully useful for the primary question if it exposes IHC score 0 separately from 1+.

## Prior Work (Closest Published Analogue)

Valieris et al., "Weakly-supervised deep learning models enable HER2-low prediction from H&E stained slides," Breast Cancer Research 2024 (PMC11331614). This is essentially the published version of this project's question. They predicted HER2-low from H&E across three cohorts: ACCCC (private, single-institution, A.C. Camargo Cancer Center, Brazil, 546 slides / 504 patients, Leica Aperio AT2), HEROHE (public), and TCGA-BRCA (535 slides). They observed external performance drop and explicitly noted that TCGA aggregates many institutions with varied protocols affecting classification, but they did not run leave-site-out or generic-embedding confound controls.

Takeaway: this project's confound analysis is more rigorous than the published state of the art on this exact task. Valieris et al. is the key citation for a cautionary-methods paper, and the ACCCC group is a natural collaborator (their cohort is the cleanest fit for the primary question).

## Candidate Shortlist

| Cohort | What it is | Single-source? | HER2 granularity | Access |
|---|---|---|---|---|
| ACCCC (A.C. Camargo, Brazil) | 546 H&E WSI / 504 pts, Leica Aperio AT2, 0.25 um/px | Yes — one institution, one scanner | neg / low / high (the exact split) | Private; request from Valieris et al. |
| BCNB (Early Breast Cancer Core-Needle Biopsy) | 1,058 core-biopsy H&E WSI, China, iScan Coreo 200x | Yes — one institution, one scanner | Clinical HER2 present; IHC-0-vs-1+ separability to confirm | Free, registration, non-commercial |
| ACROBAT | 4,212 WSI / ~1,153 pts, Swedish routine diagnostics; paired H&E + IHC (ER/PR/HER2/Ki67) consecutive sections | Yes — one source | HER2 as stained IHC slide; score likely needs deriving | Public (grand-challenge, CC) |
| HEROHE | 509 cases, single scanner (3DHistech Pannoramic 1000, Ipatimup) | Yes — one scanner | Binary positive/negative only | Public |
| Yale "HER2-TUMOR-ROIS" (TCIA) | H&E + HER2 status + trastuzumab response + tumor ROI annotations | Yes — single institution | HER2-positive focus, not low/zero | Public (TCIA) |
| IMPRESS | 126 WSI (62 HER2+, 64 TNBC), neoadjuvant chemo response, multiplex IHC (PD-L1/CD8/CD163) | Yes — single-source | HER2+ vs TNBC, not low/zero | Public |

## Recommended Path

1. Contact the Valieris / ACCCC group. They have the single best-fit cohort (single-institution, single-scanner, with neg/low/high labels) and have hit the same wall. A collaboration or data request is the cleanest possible external test of the confound finding, and this project's leave-site-out / embedding-control method is exactly what their analysis lacked.
2. Pull BCNB now (free) and check whether its HER2 field encodes IHC score 0 separately from 1+. If yes, it is an immediate single-scanner external low-versus-zero test; if binary, it still gives a clean HER2-positive-vs-negative reproducibility check under controlled acquisition.
3. Use ACROBAT as the strongest "does anything survive single-scanner" stress test: 4,000+ WSIs from one source, with paired HER2-IHC slides to derive status.
4. Use IMPRESS multiplex IHC (real PD-L1/CD8/CD163) to validate GigaTIME's virtual immune channels against measured immune markers, closing the RNA-validation gap that never closed.

## Open Items To Verify Before Committing

- BCNB: confirm exact HER2 value categories (is IHC score 0 separable from 1+, so HER2-zero is recoverable?).
- ACROBAT: confirm whether per-case HER2 IHC scores are tabulated in metadata, or only the stained IHC slides are provided (requiring score derivation).
- ACCCC: confirm data-sharing terms / whether the Valieris group will collaborate or release.

## Sources

- Valieris et al., Breast Cancer Research 2024: https://pmc.ncbi.nlm.nih.gov/articles/PMC11331614/
- ACROBAT (Scientific Data 2023): https://www.nature.com/articles/s41597-023-02422-6 ; challenge: https://acrobat.grand-challenge.org/
- BCNB: https://bupt-ai-cz.github.io/BCNB/ ; challenge: https://bcnb.grand-challenge.org/
- HEROHE (J. Imaging 2022): https://www.mdpi.com/2313-433X/8/8/213 ; challenge: https://ecdp2020.grand-challenge.org/
- Yale HER2-TUMOR-ROIS (TCIA): https://www.cancerimagingarchive.net/collection/her2-tumor-rois/
- IMPRESS (npj Precision Oncology 2023): https://www.nature.com/articles/s41698-023-00352-5
- Breast H&E WSI dataset scoping review: https://arxiv.org/html/2306.01546v2
