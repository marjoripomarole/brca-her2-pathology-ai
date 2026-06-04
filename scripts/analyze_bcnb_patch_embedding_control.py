#!/usr/bin/env python3
"""Patient-level BCNB HER2-low versus zero analysis from patch embeddings."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from analyze_classifier_permutation_sanity import fmt, make_repeated_stratified_folds, markdown_table, metric_dict
from train_her2_classifier_baseline import softmax, standardize_train_test


NEGATIVE_CLASS = "HER2-low"
POSITIVE_CLASS = "HER2-zero"
LOW_ZERO_GROUPS = [NEGATIVE_CLASS, POSITIVE_CLASS]

FEATURE_SETS = [
    ("patch_qc", "Patch-count / tissue QC"),
    ("grade_only", "Grade only"),
    ("er_pr_only", "ER/PR only"),
    ("clinical_covariates", "Clinical covariates"),
    ("embedding", "{model} embedding (PCA)"),
    ("embedding_plus_clinical", "{model} + clinical covariates"),
]
EMBEDDING_FEATURE_SETS = {"embedding", "embedding_plus_clinical"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", required=True, help="Patient-level BCNB patch embeddings CSV.")
    parser.add_argument("--model-label", default="H-Optimus-0")
    parser.add_argument("--model-id", default="")
    parser.add_argument("--out-dir", type=Path, default=Path("results/bcnb_patch_embedding_control"))
    parser.add_argument("--asset-dir", type=Path, default=Path("docs/assets/bcnb_patch_embedding_control"))
    parser.add_argument("--out-markdown", type=Path, default=Path("docs/bcnb_patch_embedding_control.md"))
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--permutations", type=int, default=200)
    parser.add_argument("--pca-components", type=int, default=20)
    parser.add_argument("--pca-grid", default="5,10,20,30,50")
    parser.add_argument("--seed", type=int, default=20260604)
    parser.add_argument("--l2-penalty", type=float, default=1.0)
    return parser.parse_args()


def require_libs(mpl_config_dir: Path):
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
        from scipy import optimize, stats
    except ModuleNotFoundError as exc:
        raise SystemExit(f"Missing Python package: {exc.name}. Use `conda run -n gigatime-tcga ...`.") from exc
    sns.set_theme(style="whitegrid", context="notebook")
    return pd, plt, sns, optimize, stats


def embedding_columns(frame) -> list[str]:
    return sorted(col for col in frame.columns if col.startswith("embedding_"))


def labeled(template: str, model_label: str) -> str:
    return template.format(model=model_label) if "{model}" in template else template


def normalize_inputs(pd, rows):
    rows = rows.copy()
    for col in ["grade", "ki67", "n_manifest_patches", "n_used_patches", "mean_tissue_fraction", "min_tissue_fraction"]:
        if col in rows.columns:
            rows[col] = pd.to_numeric(rows[col], errors="coerce")
    for col in ["ER", "PR", "molecular_subtype", "aln_status"]:
        if col in rows.columns:
            rows[col] = rows[col].fillna("Unknown").astype(str).replace({"": "Unknown", "nan": "Unknown"})
    rows["grade_missing"] = rows["grade"].isna().astype(float) if "grade" in rows.columns else 1.0
    return rows


def design_matrix(pd, rows, feature_set: str, embed_cols: list[str]):
    pieces = []
    if feature_set in {"patch_qc", "embedding_plus_clinical"}:
        cols = [col for col in ["n_manifest_patches", "n_used_patches", "mean_tissue_fraction", "min_tissue_fraction"] if col in rows]
        if cols:
            pieces.append(rows[cols].astype(float).reset_index(drop=True))
    if feature_set in {"grade_only", "clinical_covariates", "embedding_plus_clinical"}:
        cols = [col for col in ["grade", "grade_missing"] if col in rows]
        if cols:
            pieces.append(rows[cols].astype(float).reset_index(drop=True))
    if feature_set in {"er_pr_only", "clinical_covariates", "embedding_plus_clinical"}:
        cat_cols = [col for col in ["ER", "PR"] if col in rows]
        if cat_cols:
            pieces.append(pd.get_dummies(rows[cat_cols].astype(str), prefix=cat_cols, dummy_na=False).reset_index(drop=True))
    if feature_set in {"clinical_covariates", "embedding_plus_clinical"}:
        numeric = [col for col in ["ki67"] if col in rows]
        if numeric:
            pieces.append(rows[numeric].astype(float).reset_index(drop=True))
        cat_cols = [col for col in ["molecular_subtype", "aln_status"] if col in rows]
        if cat_cols:
            pieces.append(pd.get_dummies(rows[cat_cols].astype(str), prefix=cat_cols, dummy_na=False).reset_index(drop=True))
    if feature_set in {"embedding", "embedding_plus_clinical"}:
        pieces.append(rows[embed_cols].astype(float).reset_index(drop=True))
    if not pieces:
        return pd.DataFrame(index=rows.index)
    return pd.concat(pieces, axis=1)


def pca_fit_transform(x_train: np.ndarray, x_test: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    k = int(min(k, x_train.shape[0], x_train.shape[1]))
    mean = x_train.mean(axis=0)
    x_train_c = x_train - mean
    x_test_c = x_test - mean
    _u, _s, vt = np.linalg.svd(x_train_c, full_matrices=False)
    components = vt[:k]
    return x_train_c @ components.T, x_test_c @ components.T


def fit_predict_balanced_logistic(optimize, x_train, y_train, x_test, classes, l2_penalty: float):
    n_features = x_train.shape[1]
    n_classes = len(classes)
    x_aug = np.c_[np.ones(x_train.shape[0]), x_train]
    x_test_aug = np.c_[np.ones(x_test.shape[0]), x_test]
    sample_weights = np.zeros(len(y_train), dtype=float)
    for class_idx in range(n_classes):
        n_class = int(np.sum(y_train == class_idx))
        if n_class:
            sample_weights[y_train == class_idx] = len(y_train) / (n_classes * n_class)
    weight_norm = float(sample_weights.sum())

    def unpack(params):
        return params.reshape(n_features + 1, n_classes)

    def objective(params):
        weights = unpack(params)
        logits = x_aug @ weights
        probs = softmax(logits)
        log_likelihood = -np.log(probs[np.arange(len(y_train)), y_train] + 1e-12)
        nll = float((sample_weights * log_likelihood).sum() / weight_norm)
        penalty = 0.5 * l2_penalty * np.sum(weights[1:] ** 2) / len(y_train)
        grad_logits = probs
        grad_logits[np.arange(len(y_train)), y_train] -= 1.0
        grad_logits *= sample_weights[:, None] / weight_norm
        grad = x_aug.T @ grad_logits
        grad[1:] += l2_penalty * weights[1:] / len(y_train)
        return nll + penalty, grad.ravel()

    initial = np.zeros((n_features + 1, n_classes), dtype=float).ravel()
    result = optimize.minimize(objective, initial, jac=True, method="L-BFGS-B", options={"maxiter": 500})
    weights = unpack(result.x)
    return softmax(x_test_aug @ weights)


def prepare_fold_matrices(x: np.ndarray, folds, use_pca: bool, pca_k: int):
    prepared = []
    for train_idx, test_idx in folds:
        x_train, x_test = standardize_train_test(x[train_idx], x[test_idx])
        if use_pca and pca_k and x_train.shape[1] > pca_k:
            x_train, x_test = pca_fit_transform(x_train, x_test, pca_k)
        prepared.append((train_idx, test_idx, x_train, x_test))
    return prepared


def evaluate_prepared(optimize, stats, prepared, y: np.ndarray, l2_penalty: float):
    true_labels: list[int] = []
    pred_labels: list[int] = []
    positive_probs: list[float] = []
    classes = [NEGATIVE_CLASS, POSITIVE_CLASS]
    for train_idx, test_idx, x_train, x_test in prepared:
        if len(np.unique(y[train_idx])) < 2:
            continue
        probs = fit_predict_balanced_logistic(optimize, x_train, y[train_idx], x_test, classes, l2_penalty)
        true_labels.extend(y[test_idx].tolist())
        pred_labels.extend(np.argmax(probs, axis=1).astype(int).tolist())
        positive_probs.extend(probs[:, 1].astype(float).tolist())
    metrics = metric_dict(
        stats,
        np.array(true_labels, dtype=int),
        np.array(pred_labels, dtype=int),
        np.array(positive_probs, dtype=float),
    )
    metrics["n_cv_predictions"] = len(true_labels)
    return metrics


def evaluate(optimize, stats, x: np.ndarray, y: np.ndarray, folds, use_pca: bool, pca_k: int, l2_penalty: float):
    prepared = prepare_fold_matrices(x, folds, use_pca, pca_k)
    return evaluate_prepared(optimize, stats, prepared, y, l2_penalty)


def run_analysis(pd, optimize, stats, args):
    rng = np.random.default_rng(args.seed)
    rows = pd.read_csv(args.embeddings)
    rows = rows.loc[rows["clinical_her2_group"].isin(LOW_ZERO_GROUPS)].copy()
    rows = normalize_inputs(pd, rows).reset_index(drop=True)
    embed_cols = embedding_columns(rows)
    if not embed_cols:
        raise SystemExit(f"No embedding_* columns found in {args.embeddings}")

    y = rows["clinical_her2_group"].map({NEGATIVE_CLASS: 0, POSITIVE_CLASS: 1}).to_numpy(dtype=int)
    n_low = int((rows["clinical_her2_group"] == NEGATIVE_CLASS).sum())
    n_zero = int((rows["clinical_her2_group"] == POSITIVE_CLASS).sum())
    if n_low < args.folds or n_zero < args.folds:
        raise SystemExit(f"Too few low/zero patients for {args.folds}-fold CV: low={n_low}, zero={n_zero}")

    folds = make_repeated_stratified_folds(y, args.folds, args.repeats, rng)
    metric_rows = []
    designs = {}
    for feature_set, label_template in FEATURE_SETS:
        design = design_matrix(pd, rows, feature_set, embed_cols)
        if design.empty:
            continue
        x = design.to_numpy(dtype=float)
        designs[feature_set] = x
        use_pca = feature_set in EMBEDDING_FEATURE_SETS
        metrics = evaluate(optimize, stats, x, y, folds, use_pca, args.pca_components, args.l2_penalty)
        metric_rows.append(
            {
                "feature_set": feature_set,
                "feature_set_label": labeled(label_template, args.model_label),
                "n_patients": len(rows),
                "n_low": n_low,
                "n_zero": n_zero,
                "n_features": int(design.shape[1]),
                "pca_components": args.pca_components if use_pca else "",
                **metrics,
            }
        )
    metrics = pd.DataFrame(metric_rows)

    pca_grid = [int(k) for k in str(args.pca_grid).split(",") if k.strip()]
    pca_rows = []
    embed_x = designs["embedding"]
    for k in pca_grid:
        m = evaluate(optimize, stats, embed_x, y, folds, True, k, args.l2_penalty)
        pca_rows.append({"pca_components": k, **m})
    pca_robustness = pd.DataFrame(pca_rows)

    observed = metrics.loc[metrics["feature_set"] == "embedding"].iloc[0]
    observed_ba = float(observed["balanced_accuracy"])
    observed_auc = float(observed["macro_auc_ovr"])
    embedding_permutation_folds = prepare_fold_matrices(embed_x, folds, True, args.pca_components)
    null_ba = []
    null_auc = []
    for _ in range(args.permutations):
        y_perm = y.copy()
        rng.shuffle(y_perm)
        m = evaluate_prepared(optimize, stats, embedding_permutation_folds, y_perm, args.l2_penalty)
        null_ba.append(m["balanced_accuracy"])
        null_auc.append(m["macro_auc_ovr"])
    null_arr = np.array(null_ba, dtype=float)
    null_auc_arr = np.array(null_auc, dtype=float)
    permutation = {
        "feature_set": "embedding",
        "pca_components": args.pca_components,
        "observed_repeated_cv_balanced_accuracy": observed_ba,
        "observed_repeated_cv_auc": observed_auc,
        "n_permutations": args.permutations,
        "null_balanced_accuracy_mean": float(np.nanmean(null_arr)),
        "null_balanced_accuracy_sd": float(np.nanstd(null_arr, ddof=1)),
        "null_balanced_accuracy_p95": float(np.nanquantile(null_arr, 0.95)),
        "empirical_p_balanced_accuracy": float((1 + np.sum(null_arr >= observed_ba)) / (1 + args.permutations)),
        "null_auc_mean": float(np.nanmean(null_auc_arr)),
        "null_auc_sd": float(np.nanstd(null_auc_arr, ddof=1)),
        "null_auc_p95": float(np.nanquantile(null_auc_arr, 0.95)),
        "empirical_p_auc": float((1 + np.sum(null_auc_arr >= observed_auc)) / (1 + args.permutations)),
    }
    null_df = pd.DataFrame({"null_balanced_accuracy": null_arr, "null_auc": null_auc_arr})
    return rows, metrics, pca_robustness, permutation, null_df, len(embed_cols)


def metric_table(metrics) -> list[list[str]]:
    order = {name: idx for idx, (name, _) in enumerate(FEATURE_SETS)}
    selected = metrics.copy()
    selected["_order"] = selected["feature_set"].map(order)
    selected = selected.sort_values("_order")
    rows = []
    for _, row in selected.iterrows():
        rows.append(
            [
                row["feature_set_label"],
                str(int(row["n_features"])),
                str(row["pca_components"]),
                fmt(row["balanced_accuracy"], 3),
                fmt(row["macro_auc_ovr"], 3),
                fmt(row["sensitivity"], 3),
                fmt(row["specificity"], 3),
            ]
        )
    return rows


def plot_metrics(plt, sns, metrics, asset_dir: Path, model_label: str) -> None:
    order = [labeled(label, model_label) for _, label in FEATURE_SETS]
    plt.figure(figsize=(10.5, 5.2))
    sns.barplot(data=metrics, x="feature_set_label", y="balanced_accuracy", order=order)
    plt.axhline(0.5, color="#374151", linestyle="--", linewidth=1)
    plt.ylim(0, 1)
    plt.xlabel("Feature set")
    plt.ylabel("Balanced accuracy (HER2-low vs HER2-zero)")
    plt.title(f"BCNB Patient-Level Patch Embedding Control ({model_label})")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(asset_dir / "bcnb_patch_embedding_balanced_accuracy.png", dpi=180)
    plt.close()


def write_markdown(path: Path, asset_dir: Path, args, rows, metrics, pca_robustness, permutation, embedding_dim: int) -> None:
    n_low = int((rows["clinical_her2_group"] == NEGATIVE_CLASS).sum())
    n_zero = int((rows["clinical_her2_group"] == POSITIVE_CLASS).sum())
    pca_rows = [
        [str(int(row["pca_components"])), fmt(row["balanced_accuracy"], 3), fmt(row["macro_auc_ovr"], 3)]
        for _, row in pca_robustness.iterrows()
    ]
    emb = metrics.loc[metrics["feature_set"] == "embedding"].iloc[0]
    clinical = metrics.loc[metrics["feature_set"] == "clinical_covariates"].iloc[0]
    lines = [
        f"# BCNB Patient-Level Patch Embedding Control ({args.model_label})",
        "",
        "Status: BCNB external-cohort patch analysis for HER2-low versus HER2-zero.",
        "",
        "## Method",
        "",
        f"- Cohort: {len(rows)} BCNB patients with precomputed patch embeddings ({n_low} HER2-low, {n_zero} HER2-zero).",
        f"- Embedding input: patient-level mean of capped precomputed 256x256 H&E patches from `paper_patches.zip`.",
        f"- Model: `{args.model_id or args.model_label}`, {embedding_dim}-d patient embedding.",
        f"- Classifier: class-balanced regularized logistic regression with repeated stratified {args.folds}-fold CV ({args.repeats} repeats).",
        f"- Embedding dimensionality reduction: PCA fit inside each training fold only ({args.pca_components} components).",
        f"- Sanity: {args.permutations} shuffled-label permutations for the embedding.",
        "",
        "## Results",
        "",
        markdown_table(
            ["Feature set", "Features", "PCA", "Balanced accuracy", "AUC", "Sensitivity", "Specificity"],
            metric_table(metrics),
        ),
        "",
        f"![BCNB patch embedding control](assets/{asset_dir.name}/bcnb_patch_embedding_balanced_accuracy.png)",
        "",
        "## Embedding PCA Robustness",
        "",
        markdown_table(["PCA components", "Balanced accuracy", "AUC"], pca_rows),
        "",
        "## Shuffled-Label Sanity",
        "",
        markdown_table(
            ["Metric", "Observed", "Null mean", "Null 95%", "Empirical p"],
            [[
                "Balanced accuracy",
                fmt(permutation["observed_repeated_cv_balanced_accuracy"], 3),
                fmt(permutation["null_balanced_accuracy_mean"], 3),
                fmt(permutation["null_balanced_accuracy_p95"], 3),
                fmt(permutation["empirical_p_balanced_accuracy"], 4),
            ], [
                "AUC",
                fmt(permutation["observed_repeated_cv_auc"], 3),
                fmt(permutation["null_auc_mean"], 3),
                fmt(permutation["null_auc_p95"], 3),
                fmt(permutation["empirical_p_auc"], 4),
            ]],
        ),
        "",
        "## Interpretation",
        "",
        f"- {args.model_label} patch embeddings reach balanced accuracy {fmt(emb['balanced_accuracy'], 3)} "
        f"and AUC {fmt(emb['macro_auc_ovr'], 3)} versus {fmt(clinical['balanced_accuracy'], 3)} "
        f"and AUC {fmt(clinical['macro_auc_ovr'], 3)} for clinical covariates.",
        "- Interpret this as external-cohort effect-size evidence, not just a p-value: a statistically non-null but small signal is not a strong image classifier.",
        "- This is patient-level analysis, not patch-level analysis; patch-level splits would leak patient identity and overweight patients with many patches.",
        "- Because these are precomputed tumor-region patches, this does not test whole-slide slide-size or tissue-area confounding. Full WSIs remain the stronger input if the patch signal is interesting.",
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_embedding_metrics.csv'}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_embedding_pca_robustness.csv'}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_embedding_permutation.csv'}`",
        f"- `{asset_dir}/`",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, optimize, stats = require_libs(args.out_dir / ".matplotlib")
    rows, metrics, pca_robustness, permutation, null_df, embedding_dim = run_analysis(pd, optimize, stats, args)

    metrics.to_csv(args.out_dir / "bcnb_patch_embedding_metrics.csv", index=False)
    pca_robustness.to_csv(args.out_dir / "bcnb_patch_embedding_pca_robustness.csv", index=False)
    pd.DataFrame([permutation]).to_csv(args.out_dir / "bcnb_patch_embedding_permutation.csv", index=False)
    null_df.to_csv(args.out_dir / "bcnb_patch_embedding_permutation_null.csv", index=False)
    (args.out_dir / "bcnb_patch_embedding_metadata.json").write_text(
        json.dumps(
            {
                "task": "bcnb_her2_low_vs_zero_patch_embedding_control",
                "model_label": args.model_label,
                "model_id": args.model_id,
                "embedding_dimensions": embedding_dim,
                "folds": args.folds,
                "repeats": args.repeats,
                "permutations": args.permutations,
                "pca_components": args.pca_components,
                "seed": args.seed,
                "n_patients": int(len(rows)),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    plot_metrics(plt, sns, metrics, args.asset_dir, args.model_label)
    write_markdown(args.out_markdown, args.asset_dir, args, rows, metrics, pca_robustness, permutation, embedding_dim)
    print(f"Wrote BCNB patch embedding control outputs to {args.out_dir}")
    print(f"Wrote BCNB patch embedding markdown to {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
