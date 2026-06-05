#!/usr/bin/env python3
"""Compare BCNB HER2-low versus zero signal across patch embedding models."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from analyze_bcnb_patch_embedding_control import (
    LOW_ZERO_GROUPS,
    NEGATIVE_CLASS,
    POSITIVE_CLASS,
    design_matrix,
    embedding_columns,
    fit_predict_balanced_logistic,
    normalize_inputs,
    pca_fit_transform,
)
from analyze_classifier_permutation_sanity import fmt, make_repeated_stratified_folds, markdown_table, metric_dict
from train_her2_classifier_baseline import standardize_train_test


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hoptimus-embeddings",
        required=True,
        help="Patient-level H-Optimus BCNB patch embeddings CSV.",
    )
    parser.add_argument(
        "--virchow2-embeddings",
        required=True,
        help="Patient-level Virchow2 BCNB patch embeddings CSV.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero"),
    )
    parser.add_argument(
        "--asset-dir",
        type=Path,
        default=Path("docs/assets/bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero"),
    )
    parser.add_argument(
        "--out-markdown",
        type=Path,
        default=Path("docs/bcnb_patch_model_comparison_hoptimus0_virchow2_hash_capped10_low_zero.md"),
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--permutations", type=int, default=200)
    parser.add_argument("--pca-components", type=int, default=20)
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


def load_aligned_embeddings(pd, hoptimus_path: Path, virchow_path: Path):
    hoptimus = pd.read_csv(hoptimus_path, low_memory=False)
    virchow = pd.read_csv(virchow_path, low_memory=False)
    hoptimus = hoptimus.loc[hoptimus["clinical_her2_group"].isin(LOW_ZERO_GROUPS)].copy()
    virchow = virchow.loc[virchow["clinical_her2_group"].isin(LOW_ZERO_GROUPS)].copy()
    hoptimus = hoptimus.sort_values("patient_id").reset_index(drop=True)
    virchow = virchow.sort_values("patient_id").reset_index(drop=True)
    if not hoptimus["patient_id"].equals(virchow["patient_id"]):
        raise SystemExit("H-Optimus and Virchow2 embedding files do not contain the same ordered patient IDs.")

    check_cols = [
        "clinical_her2_group",
        "her2_status",
        "her2_ihc",
        "grade",
        "ER",
        "PR",
        "ki67",
        "molecular_subtype",
        "aln_status",
        "n_manifest_patches",
        "n_used_patches",
        "mean_tissue_fraction",
        "min_tissue_fraction",
    ]
    for col in check_cols:
        if col in hoptimus.columns and col in virchow.columns:
            left = hoptimus[col].fillna("__NA__").astype(str)
            right = virchow[col].fillna("__NA__").astype(str)
            if not left.equals(right):
                raise SystemExit(f"Column mismatch between embedding files: {col}")

    h_cols = embedding_columns(hoptimus)
    v_cols = embedding_columns(virchow)
    if not h_cols or not v_cols:
        raise SystemExit("Both input files must contain embedding_* columns.")
    meta_cols = [col for col in hoptimus.columns if not col.startswith("embedding_")]
    rows = hoptimus[meta_cols].copy()
    return rows, hoptimus[h_cols].to_numpy(dtype=float), virchow[v_cols].to_numpy(dtype=float), len(h_cols), len(v_cols)


def prepare_blocks(blocks: list[tuple[str, np.ndarray, int | None]], folds):
    prepared = []
    for train_idx, test_idx in folds:
        train_parts = []
        test_parts = []
        for _name, values, pca_k in blocks:
            x_train, x_test = standardize_train_test(values[train_idx], values[test_idx])
            if pca_k and x_train.shape[1] > pca_k:
                x_train, x_test = pca_fit_transform(x_train, x_test, pca_k)
            train_parts.append(x_train)
            test_parts.append(x_test)
        prepared.append((train_idx, test_idx, np.concatenate(train_parts, axis=1), np.concatenate(test_parts, axis=1)))
    return prepared


def evaluate_prepared_with_oof(optimize, stats, prepared, y: np.ndarray, patient_ids: np.ndarray, l2_penalty: float):
    true_labels = []
    pred_labels = []
    positive_probs = []
    oof_rows = []
    classes = [NEGATIVE_CLASS, POSITIVE_CLASS]
    for fold_id, (train_idx, test_idx, x_train, x_test) in enumerate(prepared):
        if len(np.unique(y[train_idx])) < 2:
            continue
        probs = fit_predict_balanced_logistic(optimize, x_train, y[train_idx], x_test, classes, l2_penalty)
        pred = np.argmax(probs, axis=1).astype(int)
        positive = probs[:, 1].astype(float)
        true_labels.extend(y[test_idx].tolist())
        pred_labels.extend(pred.tolist())
        positive_probs.extend(positive.tolist())
        for row_idx, patient_id in enumerate(patient_ids[test_idx]):
            oof_rows.append(
                {
                    "fold_id": int(fold_id),
                    "patient_id": patient_id,
                    "y_true": int(y[test_idx][row_idx]),
                    "prob_her2_zero": float(positive[row_idx]),
                    "predicted_label": POSITIVE_CLASS if int(pred[row_idx]) == 1 else NEGATIVE_CLASS,
                }
            )
    metrics = metric_dict(
        stats,
        np.array(true_labels, dtype=int),
        np.array(pred_labels, dtype=int),
        np.array(positive_probs, dtype=float),
    )
    metrics["n_cv_predictions"] = len(true_labels)
    return metrics, oof_rows


def permutation_sanity(rng, optimize, stats, prepared, y: np.ndarray, observed_ba: float, observed_auc: float, args):
    null_ba = []
    null_auc = []
    patient_ids = np.arange(len(y))
    for _ in range(args.permutations):
        y_perm = y.copy()
        rng.shuffle(y_perm)
        metrics, _ = evaluate_prepared_with_oof(optimize, stats, prepared, y_perm, patient_ids, args.l2_penalty)
        null_ba.append(metrics["balanced_accuracy"])
        null_auc.append(metrics["macro_auc_ovr"])
    null_ba = np.array(null_ba, dtype=float)
    null_auc = np.array(null_auc, dtype=float)
    summary = {
        "feature_set": "dual_model_block_pca",
        "observed_repeated_cv_balanced_accuracy": observed_ba,
        "observed_repeated_cv_auc": observed_auc,
        "n_permutations": args.permutations,
        "null_balanced_accuracy_mean": float(np.nanmean(null_ba)),
        "null_balanced_accuracy_sd": float(np.nanstd(null_ba, ddof=1)),
        "null_balanced_accuracy_p95": float(np.nanquantile(null_ba, 0.95)),
        "empirical_p_balanced_accuracy": float((1 + np.sum(null_ba >= observed_ba)) / (1 + args.permutations)),
        "null_auc_mean": float(np.nanmean(null_auc)),
        "null_auc_sd": float(np.nanstd(null_auc, ddof=1)),
        "null_auc_p95": float(np.nanquantile(null_auc, 0.95)),
        "empirical_p_auc": float((1 + np.sum(null_auc >= observed_auc)) / (1 + args.permutations)),
    }
    null_df = {
        "null_balanced_accuracy": null_ba.tolist(),
        "null_auc": null_auc.tolist(),
    }
    return summary, null_df


def ensemble_from_oof(pd, stats, h_oof, v_oof):
    h_frame = pd.DataFrame(h_oof).rename(
        columns={"prob_her2_zero": "hoptimus_prob_her2_zero", "predicted_label": "hoptimus_predicted_label"}
    )
    v_frame = pd.DataFrame(v_oof).rename(
        columns={"prob_her2_zero": "virchow2_prob_her2_zero", "predicted_label": "virchow2_predicted_label"}
    )
    merged = h_frame.merge(
        v_frame[["fold_id", "patient_id", "virchow2_prob_her2_zero", "virchow2_predicted_label"]],
        on=["fold_id", "patient_id"],
        how="inner",
        validate="one_to_one",
    )
    merged["ensemble_prob_her2_zero"] = (
        merged["hoptimus_prob_her2_zero"] + merged["virchow2_prob_her2_zero"]
    ) / 2.0
    merged["ensemble_pred"] = np.where(merged["ensemble_prob_her2_zero"] >= 0.5, 1, 0)
    metrics = metric_dict(
        stats,
        merged["y_true"].to_numpy(dtype=int),
        merged["ensemble_pred"].to_numpy(dtype=int),
        merged["ensemble_prob_her2_zero"].to_numpy(dtype=float),
    )
    metrics["n_cv_predictions"] = len(merged)
    return metrics, merged


def agreement_summary(stats, oof_predictions):
    patient = (
        oof_predictions.groupby("patient_id", as_index=False)
        .agg(
            y_true=("y_true", "first"),
            hoptimus_prob_her2_zero=("hoptimus_prob_her2_zero", "mean"),
            virchow2_prob_her2_zero=("virchow2_prob_her2_zero", "mean"),
            ensemble_prob_her2_zero=("ensemble_prob_her2_zero", "mean"),
        )
        .copy()
    )
    patient["clinical_her2_group"] = np.where(patient["y_true"] == 1, POSITIVE_CLASS, NEGATIVE_CLASS)
    patient["hoptimus_pred"] = np.where(patient["hoptimus_prob_her2_zero"] >= 0.5, 1, 0)
    patient["virchow2_pred"] = np.where(patient["virchow2_prob_her2_zero"] >= 0.5, 1, 0)
    oof_pearson = stats.pearsonr(
        oof_predictions["hoptimus_prob_her2_zero"], oof_predictions["virchow2_prob_her2_zero"]
    )
    oof_spearman = stats.spearmanr(
        oof_predictions["hoptimus_prob_her2_zero"], oof_predictions["virchow2_prob_her2_zero"]
    )
    patient_pearson = stats.pearsonr(patient["hoptimus_prob_her2_zero"], patient["virchow2_prob_her2_zero"])
    patient_spearman = stats.spearmanr(patient["hoptimus_prob_her2_zero"], patient["virchow2_prob_her2_zero"])
    summary = {
        "n_patients": int(len(patient)),
        "n_oof_predictions": int(len(oof_predictions)),
        "oof_probability_pearson_r": float(oof_pearson.statistic),
        "oof_probability_pearson_p": float(oof_pearson.pvalue),
        "oof_probability_spearman_rho": float(oof_spearman.statistic),
        "oof_probability_spearman_p": float(oof_spearman.pvalue),
        "patient_mean_probability_pearson_r": float(patient_pearson.statistic),
        "patient_mean_probability_pearson_p": float(patient_pearson.pvalue),
        "patient_mean_probability_spearman_rho": float(patient_spearman.statistic),
        "patient_mean_probability_spearman_p": float(patient_spearman.pvalue),
        "patient_threshold_disagreement_fraction": float((patient["hoptimus_pred"] != patient["virchow2_pred"]).mean()),
        "patient_mean_absolute_probability_difference": float(
            np.mean(np.abs(patient["hoptimus_prob_her2_zero"] - patient["virchow2_prob_her2_zero"]))
        ),
    }
    return summary, patient


def metric_rows(metrics_by_name: dict[str, dict], feature_counts: dict[str, int]) -> list[list[str]]:
    rows = []
    for name, metrics in metrics_by_name.items():
        rows.append(
            [
                name,
                str(feature_counts.get(name, "")),
                fmt(metrics["balanced_accuracy"], 3),
                fmt(metrics["macro_auc_ovr"], 3),
                fmt(metrics["sensitivity"], 3),
                fmt(metrics["specificity"], 3),
            ]
        )
    return rows


def plot_metrics(plt, sns, pd, metrics_by_name: dict[str, dict], asset_dir: Path) -> None:
    plot_df = pd.DataFrame(
        [
            {
                "feature_set": name,
                "Balanced accuracy": metrics["balanced_accuracy"],
                "AUC": metrics["macro_auc_ovr"],
            }
            for name, metrics in metrics_by_name.items()
        ]
    )
    long_df = plot_df.melt(id_vars="feature_set", value_vars=["Balanced accuracy", "AUC"], var_name="metric")
    fig, axis = plt.subplots(figsize=(10.8, 5.2))
    sns.barplot(data=long_df, x="feature_set", y="value", hue="metric", ax=axis)
    axis.axhline(0.5, color="#374151", linestyle="--", linewidth=1)
    axis.set_ylim(0, 1)
    axis.set_xlabel("Feature set")
    axis.set_ylabel("Repeated-CV score")
    axis.set_title("BCNB H-Optimus-0 vs Virchow2 Patch Model Comparison")
    axis.tick_params(axis="x", rotation=25)
    for label in axis.get_xticklabels():
        label.set_horizontalalignment("right")
    fig.tight_layout()
    fig.savefig(asset_dir / "bcnb_patch_model_comparison_metrics.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_probability_agreement(plt, sns, patient_predictions, asset_dir: Path) -> None:
    fig, axis = plt.subplots(figsize=(6.2, 5.6))
    sns.scatterplot(
        data=patient_predictions,
        x="hoptimus_prob_her2_zero",
        y="virchow2_prob_her2_zero",
        hue="clinical_her2_group",
        alpha=0.75,
        linewidth=0,
        ax=axis,
    )
    axis.axhline(0.5, color="#374151", linestyle="--", linewidth=1)
    axis.axvline(0.5, color="#374151", linestyle="--", linewidth=1)
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_xlabel("H-Optimus-0 mean out-of-fold P(HER2-zero)")
    axis.set_ylabel("Virchow2 mean out-of-fold P(HER2-zero)")
    axis.set_title("Patient-Level Cross-Model Probability Agreement")
    fig.tight_layout()
    fig.savefig(asset_dir / "bcnb_patch_model_probability_agreement.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_markdown(path: Path, asset_dir: Path, args, rows, h_dim, v_dim, metrics_by_name, feature_counts, agreement, permutation):
    n_low = int((rows["clinical_her2_group"] == NEGATIVE_CLASS).sum())
    n_zero = int((rows["clinical_her2_group"] == POSITIVE_CLASS).sum())
    lines = [
        "# BCNB Patch Model Comparison: H-Optimus-0 Versus Virchow2",
        "",
        "Status: paired external-cohort comparison using the same BCNB patients and the same deterministic patch sample.",
        "",
        "## Method",
        "",
        f"- Cohort: {len(rows)} BCNB patients ({n_low} HER2-low, {n_zero} HER2-zero).",
        "- Inputs: patient-level means of the same 10 deterministic hash-sampled 256x256 tumor-region patches per patient.",
        f"- Models: H-Optimus-0 ({h_dim}-d) and Virchow2 ({v_dim}-d).",
        f"- Classifier: class-balanced regularized logistic regression with repeated stratified {args.folds}-fold CV ({args.repeats} repeats).",
        f"- Dimensionality control: PCA is fit inside each training fold only; dual-model rows use {args.pca_components} components per embedding model.",
        f"- Sanity: {args.permutations} shuffled-label permutations for the dual-model embedding classifier.",
        "",
        "## Results",
        "",
        markdown_table(
            ["Feature set", "Features after PCA", "Balanced accuracy", "AUC", "Sensitivity", "Specificity"],
            metric_rows(metrics_by_name, feature_counts),
        ),
        "",
        f"![BCNB patch model comparison](assets/{asset_dir.name}/bcnb_patch_model_comparison_metrics.png)",
        "",
        "## Cross-Model Agreement",
        "",
        markdown_table(
            ["Agreement measure", "Value"],
            [
                ["OOF probability Pearson r", fmt(agreement["oof_probability_pearson_r"], 3)],
                ["OOF probability Spearman rho", fmt(agreement["oof_probability_spearman_rho"], 3)],
                ["Patient-mean probability Pearson r", fmt(agreement["patient_mean_probability_pearson_r"], 3)],
                ["Patient-mean probability Spearman rho", fmt(agreement["patient_mean_probability_spearman_rho"], 3)],
                ["Patient threshold disagreement fraction", fmt(agreement["patient_threshold_disagreement_fraction"], 3)],
                ["Patient mean absolute probability difference", fmt(agreement["patient_mean_absolute_probability_difference"], 3)],
            ],
        ),
        "",
        f"![BCNB patch model probability agreement](assets/{asset_dir.name}/bcnb_patch_model_probability_agreement.png)",
        "",
        "## Dual-Model Shuffled-Label Sanity",
        "",
        markdown_table(
            ["Metric", "Observed", "Null mean", "Null 95%", "Empirical p"],
            [
                [
                    "Balanced accuracy",
                    fmt(permutation["observed_repeated_cv_balanced_accuracy"], 3),
                    fmt(permutation["null_balanced_accuracy_mean"], 3),
                    fmt(permutation["null_balanced_accuracy_p95"], 3),
                    fmt(permutation["empirical_p_balanced_accuracy"], 4),
                ],
                [
                    "AUC",
                    fmt(permutation["observed_repeated_cv_auc"], 3),
                    fmt(permutation["null_auc_mean"], 3),
                    fmt(permutation["null_auc_p95"], 3),
                    fmt(permutation["empirical_p_auc"], 4),
                ],
            ],
        ),
        "",
        "## Interpretation",
        "",
        "- H-Optimus-0 and Virchow2 produce concordant but not identical patient scores on the same BCNB patches.",
        "- The dual-model embedding does not create a large jump over either model alone, which argues against a hidden strong signal missed by one encoder.",
        "- Clinical covariates remain at least as strong by balanced accuracy, so the current BCNB result should be framed as a weak, reproducible morphology/covariate-associated signal rather than a clinically deployable HER2-low versus zero classifier.",
        "- The next escalation, if needed for a manuscript, is patch-sampling sensitivity or full-WSI processing to test whether broader tissue context changes the effect size.",
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_model_comparison_metrics.csv'}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_model_oof_predictions.csv'}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_model_patient_predictions.csv'}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_model_agreement.json'}`",
        f"- `{asset_dir}/`",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def run_analysis(pd, plt, sns, optimize, stats, args):
    rng = np.random.default_rng(args.seed)
    rows, h_x, v_x, h_dim, v_dim = load_aligned_embeddings(
        pd,
        Path(args.hoptimus_embeddings),
        Path(args.virchow2_embeddings),
    )
    rows = normalize_inputs(pd, rows).reset_index(drop=True)
    patient_ids = rows["patient_id"].to_numpy()
    y = rows["clinical_her2_group"].map({NEGATIVE_CLASS: 0, POSITIVE_CLASS: 1}).to_numpy(dtype=int)
    folds = make_repeated_stratified_folds(y, args.folds, args.repeats, rng)
    clinical_x = design_matrix(pd, rows, "clinical_covariates", []).to_numpy(dtype=float)

    feature_specs = [
        ("Clinical covariates", [("clinical", clinical_x, None)]),
        ("H-Optimus-0 embedding", [("hoptimus0", h_x, args.pca_components)]),
        ("Virchow2 embedding", [("virchow2", v_x, args.pca_components)]),
        (
            "H-Optimus-0 + Virchow2",
            [("hoptimus0", h_x, args.pca_components), ("virchow2", v_x, args.pca_components)],
        ),
        (
            "H-Optimus-0 + Virchow2 + clinical",
            [("clinical", clinical_x, None), ("hoptimus0", h_x, args.pca_components), ("virchow2", v_x, args.pca_components)],
        ),
    ]
    metrics_by_name = {}
    oof_by_name = {}
    prepared_by_name = {}
    feature_counts = {}
    for name, blocks in feature_specs:
        prepared = prepare_blocks(blocks, folds)
        metrics, oof = evaluate_prepared_with_oof(optimize, stats, prepared, y, patient_ids, args.l2_penalty)
        metrics_by_name[name] = metrics
        oof_by_name[name] = oof
        prepared_by_name[name] = prepared
        feature_counts[name] = int(prepared[0][2].shape[1]) if prepared else 0

    ensemble_metrics, oof_predictions = ensemble_from_oof(
        pd,
        stats,
        oof_by_name["H-Optimus-0 embedding"],
        oof_by_name["Virchow2 embedding"],
    )
    metrics_by_name["Average probability ensemble"] = ensemble_metrics
    feature_counts["Average probability ensemble"] = 2
    agreement, patient_predictions = agreement_summary(stats, oof_predictions)

    dual_metrics = metrics_by_name["H-Optimus-0 + Virchow2"]
    permutation, null_df = permutation_sanity(
        rng,
        optimize,
        stats,
        prepared_by_name["H-Optimus-0 + Virchow2"],
        y,
        float(dual_metrics["balanced_accuracy"]),
        float(dual_metrics["macro_auc_ovr"]),
        args,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    metrics_table = pd.DataFrame(
        [
            {
                "feature_set": name,
                "n_features_after_pca": feature_counts[name],
                **metrics,
            }
            for name, metrics in metrics_by_name.items()
        ]
    )
    metrics_table.to_csv(args.out_dir / "bcnb_patch_model_comparison_metrics.csv", index=False)
    oof_predictions.to_csv(args.out_dir / "bcnb_patch_model_oof_predictions.csv", index=False)
    patient_predictions.to_csv(args.out_dir / "bcnb_patch_model_patient_predictions.csv", index=False)
    pd.DataFrame([permutation]).to_csv(args.out_dir / "bcnb_patch_model_dual_permutation.csv", index=False)
    pd.DataFrame(null_df).to_csv(args.out_dir / "bcnb_patch_model_dual_permutation_null.csv", index=False)
    (args.out_dir / "bcnb_patch_model_agreement.json").write_text(json.dumps(agreement, indent=2) + "\n", encoding="utf-8")
    plot_metrics(plt, sns, pd, metrics_by_name, args.asset_dir)
    plot_probability_agreement(plt, sns, patient_predictions, args.asset_dir)
    write_markdown(args.out_markdown, args.asset_dir, args, rows, h_dim, v_dim, metrics_by_name, feature_counts, agreement, permutation)
    return metrics_table, agreement, permutation


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, optimize, stats = require_libs(args.out_dir / ".matplotlib")
    metrics, agreement, permutation = run_analysis(pd, plt, sns, optimize, stats, args)
    print(f"Wrote BCNB patch model comparison outputs to {args.out_dir}")
    print(f"Wrote BCNB patch model comparison markdown to {args.out_markdown}")
    print(metrics[["feature_set", "balanced_accuracy", "macro_auc_ovr"]].to_string(index=False))
    print(json.dumps(agreement, indent=2))
    print(json.dumps(permutation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
