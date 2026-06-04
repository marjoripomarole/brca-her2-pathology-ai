# BCNB Patient-Level Patch Embedding Control (H-Optimus-0)

Status: BCNB external-cohort patch analysis for HER2-low versus HER2-zero.

## Method

- Cohort: 781 BCNB patients with precomputed patch embeddings (654 HER2-low, 127 HER2-zero).
- Embedding input: patient-level mean of capped precomputed 256x256 H&E patches from `paper_patches.zip`.
- Model: `bioptimus/H-optimus-0`, 1536-d patient embedding.
- Classifier: class-balanced regularized logistic regression with repeated stratified 5-fold CV (5 repeats).
- Embedding dimensionality reduction: PCA fit inside each training fold only (20 components).
- Sanity: 200 shuffled-label permutations for the embedding.

## Results

| Feature set | Features | PCA | Balanced accuracy | AUC | Sensitivity | Specificity |
| --- | --- | --- | --- | --- | --- | --- |
| Patch-count / tissue QC | 4 |  | 0.559 | 0.562 | 0.578 | 0.540 |
| Grade only | 2 |  | 0.595 | 0.604 | 0.517 | 0.673 |
| ER/PR only | 4 |  | 0.606 | 0.570 | 0.354 | 0.858 |
| Clinical covariates | 13 |  | 0.643 | 0.627 | 0.532 | 0.753 |
| H-Optimus-0 embedding (PCA) | 1536 | 20 | 0.597 | 0.640 | 0.539 | 0.655 |
| H-Optimus-0 + clinical covariates | 1553 | 20 | 0.595 | 0.641 | 0.532 | 0.658 |

![BCNB patch embedding control](assets/bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero/bcnb_patch_embedding_balanced_accuracy.png)

## Embedding PCA Robustness

| PCA components | Balanced accuracy | AUC |
| --- | --- | --- |
| 5 | 0.573 | 0.590 |
| 10 | 0.590 | 0.611 |
| 20 | 0.597 | 0.640 |
| 30 | 0.626 | 0.680 |
| 50 | 0.610 | 0.663 |

## Shuffled-Label Sanity

| Metric | Observed | Null mean | Null 95% | Empirical p |
| --- | --- | --- | --- | --- |
| Balanced accuracy | 0.597 | 0.499 | 0.540 | 0.0050 |
| AUC | 0.640 | 0.495 | 0.553 | 0.0050 |

## Interpretation

- H-Optimus-0 patch embeddings reach balanced accuracy 0.597 and AUC 0.640 versus 0.643 and AUC 0.627 for clinical covariates.
- Interpret this as external-cohort effect-size evidence, not just a p-value: a statistically non-null but small signal is not a strong image classifier.
- This is patient-level analysis, not patch-level analysis; patch-level splits would leak patient identity and overweight patients with many patches.
- Because these are precomputed tumor-region patches, this does not test whole-slide slide-size or tissue-area confounding. Full WSIs remain the stronger input if the patch signal is interesting.

## Output Files

- `docs/bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero.md`
- `results/bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero/bcnb_patch_embedding_metrics.csv`
- `results/bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero/bcnb_patch_embedding_pca_robustness.csv`
- `results/bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero/bcnb_patch_embedding_permutation.csv`
- `docs/assets/bcnb_patch_embedding_control_hoptimus0_hash_capped10_low_zero/`
