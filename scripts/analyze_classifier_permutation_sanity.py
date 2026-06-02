#!/usr/bin/env python3
"""Permutation sanity check for HER2-low versus HER2-zero classifiers."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np

from train_her2_classifier_baseline import binary_auc, fit_predict_logistic, standardize_train_test
from train_her2_cleaned_classifier_comparison import FEATURE_LABELS, feature_sets_for


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
DEFAULT_TUMOR_PROXY_DIR = BASE_RESULT_DIR / "tumor_proxy_sensitivity"
POSITIVE_CLASS = "HER2-zero"
NEGATIVE_CLASS = "HER2-low"
LOW_ZERO_GROUPS = [NEGATIVE_CLASS, POSITIVE_CLASS]
VIEW_ORDER = [
    "qc_cellular_tissue",
    "ck_top25_within_slide",
    "ck_top16_within_slide",
    "ck_top8_within_slide",
    "ck_top16_non_low_marker",
    "absolute_ck_high_q75",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--slide-features",
        default=str(DEFAULT_TUMOR_PROXY_DIR / "tumor_proxy_slide_features.csv"),
        help="Slide-level tumor-proxy feature table.",
    )
    parser.add_argument(
        "--best-metrics",
        default=str(DEFAULT_TUMOR_PROXY_DIR / "tumor_proxy_classifier_best_h_e_metrics.csv"),
        help="Best H&E/GigaTIME low-vs-zero feature set per view.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(BASE_RESULT_DIR / "classifier_permutation_sanity"),
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_classifier_permutation",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_classifier_permutation_sanity.md",
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260602)
    parser.add_argument("--l2-penalty", type=float, default=1.0)
    return parser.parse_args()


def require_analysis_libs(mpl_config_dir: Path):
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
        from scipy import optimize, stats
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"Missing Python package: {exc.name}. Use `conda run -n gigatime-tcga ...`."
        ) from exc
    sns.set_theme(style="whitegrid", context="notebook")
    return pd, plt, sns, optimize, stats


def fmt(value: object, digits: int = 3) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "" if value is None else str(value)
    if numeric != numeric:
        return ""
    if abs(numeric) < 0.001 and numeric != 0:
        return f"{numeric:.2e}"
    return f"{numeric:.{digits}f}"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    if rows:
        lines.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    else:
        lines.append("| " + " | ".join("" for _ in headers) + " |")
    return "\n".join(lines)


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    indexed = [(idx, p) for idx, p in enumerate(p_values) if not math.isnan(p)]
    if not indexed:
        return [float("nan")] * len(p_values)
    ranked = sorted(indexed, key=lambda item: item[1])
    adjusted = [float("nan")] * len(p_values)
    previous = 1.0
    m = len(ranked)
    for rank, (idx, p_value) in reversed(list(enumerate(ranked, start=1))):
        q_value = min(previous, p_value * m / rank)
        adjusted[idx] = q_value
        previous = q_value
    return adjusted


def make_repeated_stratified_folds(labels: np.ndarray, n_folds: int, n_repeats: int, rng: np.random.Generator):
    unique = np.unique(labels)
    folds = []
    for _repeat in range(n_repeats):
        fold_indices = [[] for _ in range(n_folds)]
        for label in unique:
            label_indices = np.where(labels == label)[0].copy()
            rng.shuffle(label_indices)
            for offset, index in enumerate(label_indices):
                fold_indices[offset % n_folds].append(int(index))
        for fold in fold_indices:
            test_idx = np.array(sorted(fold), dtype=int)
            train_mask = np.ones(len(labels), dtype=bool)
            train_mask[test_idx] = False
            train_idx = np.where(train_mask)[0]
            if len(np.unique(labels[train_idx])) < len(unique):
                continue
            folds.append((train_idx, test_idx))
    return folds


def metric_dict(stats, y_true: np.ndarray, y_pred: np.ndarray, y_prob_positive: np.ndarray) -> dict[str, float]:
    low_mask = y_true == 0
    zero_mask = y_true == 1
    low_recall = float(np.mean(y_pred[low_mask] == 0)) if low_mask.any() else float("nan")
    zero_recall = float(np.mean(y_pred[zero_mask] == 1)) if zero_mask.any() else float("nan")
    specificity = low_recall
    sensitivity = zero_recall
    return {
        "accuracy": float(np.mean(y_true == y_pred)),
        "balanced_accuracy": float(np.nanmean([low_recall, zero_recall])),
        "macro_auc_ovr": binary_auc(stats, y_true.astype(int), y_prob_positive),
        "sensitivity": sensitivity,
        "specificity": specificity,
    }


def evaluate_cv(optimize, stats, x: np.ndarray, y: np.ndarray, folds, l2_penalty: float) -> dict[str, float]:
    true_labels = []
    pred_labels = []
    positive_probs = []
    classes = [NEGATIVE_CLASS, POSITIVE_CLASS]
    for train_idx, test_idx in folds:
        y_train = y[train_idx]
        if len(np.unique(y_train)) < 2:
            continue
        x_train, x_test = standardize_train_test(x[train_idx], x[test_idx])
        probs = fit_predict_logistic(optimize, x_train, y_train, x_test, classes, l2_penalty)
        true_labels.extend(y[test_idx].tolist())
        pred_labels.extend(np.argmax(probs, axis=1).astype(int).tolist())
        positive_probs.extend(probs[:, 1].astype(float).tolist())
    return metric_dict(
        stats,
        np.array(true_labels, dtype=int),
        np.array(pred_labels, dtype=int),
        np.array(positive_probs, dtype=float),
    )


def selected_low_zero_feature_sets(pd, best_metrics_path: Path) -> dict[str, str]:
    best = pd.read_csv(best_metrics_path)
    selected = best.loc[
        (best["task"] == "her2_low_vs_zero")
        & (best["model"] == "regularized_logistic")
        & (best["feature_set"] != "erbb2_rna_reference_not_h_e")
    ].copy()
    return dict(zip(selected["feature_view"], selected["feature_set"]))


def run_permutation_checks(pd, optimize, stats, args: argparse.Namespace):
    rng = np.random.default_rng(args.seed)
    slide_features = pd.read_csv(args.slide_features)
    best_metrics = pd.read_csv(args.best_metrics)
    selected_feature_sets = selected_low_zero_feature_sets(pd, Path(args.best_metrics))
    summary_rows = []
    permutation_rows = []

    for view in VIEW_ORDER:
        rows = slide_features.loc[
            (slide_features["feature_view"] == view)
            & (slide_features["clinical_her2_group"].isin(LOW_ZERO_GROUPS))
        ].copy()
        if rows.empty or view not in selected_feature_sets:
            continue
        feature_sets = feature_sets_for(rows)
        feature_set = selected_feature_sets[view]
        feature_cols = feature_sets.get(feature_set)
        if not feature_cols:
            continue
        rows = rows.dropna(subset=feature_cols + ["clinical_her2_group"]).copy()
        y = rows["clinical_her2_group"].map({NEGATIVE_CLASS: 0, POSITIVE_CLASS: 1}).to_numpy(dtype=int)
        x = rows[feature_cols].to_numpy(dtype=float)
        folds = make_repeated_stratified_folds(y, args.folds, args.repeats, rng)
        observed = evaluate_cv(optimize, stats, x, y, folds, args.l2_penalty)

        view_best = best_metrics.loc[
            (best_metrics["feature_view"] == view)
            & (best_metrics["task"] == "her2_low_vs_zero")
            & (best_metrics["feature_set"] == feature_set)
            & (best_metrics["model"] == "regularized_logistic")
        ].iloc[0]

        null_balanced = []
        null_auc = []
        for permutation_idx in range(args.permutations):
            y_perm = y.copy()
            rng.shuffle(y_perm)
            permuted = evaluate_cv(optimize, stats, x, y_perm, folds, args.l2_penalty)
            null_balanced.append(permuted["balanced_accuracy"])
            null_auc.append(permuted["macro_auc_ovr"])
            permutation_rows.append(
                {
                    "feature_view": view,
                    "feature_view_label": rows["feature_view_label"].iloc[0],
                    "feature_set": feature_set,
                    "feature_set_label": FEATURE_LABELS.get(feature_set, feature_set),
                    "permutation_index": permutation_idx,
                    "balanced_accuracy": permuted["balanced_accuracy"],
                    "macro_auc_ovr": permuted["macro_auc_ovr"],
                }
            )

        null_balanced_array = np.array(null_balanced, dtype=float)
        null_auc_array = np.array(null_auc, dtype=float)
        p_balanced = float((1 + np.sum(null_balanced_array >= observed["balanced_accuracy"])) / (1 + args.permutations))
        p_auc = float((1 + np.sum(null_auc_array >= observed["macro_auc_ovr"])) / (1 + args.permutations))
        summary_rows.append(
            {
                "feature_view": view,
                "feature_view_label": rows["feature_view_label"].iloc[0],
                "feature_set": feature_set,
                "feature_set_label": FEATURE_LABELS.get(feature_set, feature_set),
                "n_cases": int(len(rows)),
                "n_features": int(len(feature_cols)),
                "cv_folds": int(args.folds),
                "cv_repeats": int(args.repeats),
                "n_cv_predictions": int(len(rows) * args.repeats),
                "n_permutations": int(args.permutations),
                "loocv_balanced_accuracy": float(view_best["balanced_accuracy"]),
                "loocv_macro_auc_ovr": float(view_best["macro_auc_ovr"]),
                "observed_repeated_cv_accuracy": observed["accuracy"],
                "observed_repeated_cv_balanced_accuracy": observed["balanced_accuracy"],
                "observed_repeated_cv_macro_auc_ovr": observed["macro_auc_ovr"],
                "observed_repeated_cv_sensitivity": observed["sensitivity"],
                "observed_repeated_cv_specificity": observed["specificity"],
                "null_balanced_accuracy_mean": float(np.nanmean(null_balanced_array)),
                "null_balanced_accuracy_sd": float(np.nanstd(null_balanced_array, ddof=1)),
                "null_balanced_accuracy_p95": float(np.nanquantile(null_balanced_array, 0.95)),
                "empirical_p_balanced_accuracy": p_balanced,
                "null_macro_auc_mean": float(np.nanmean(null_auc_array)),
                "null_macro_auc_sd": float(np.nanstd(null_auc_array, ddof=1)),
                "null_macro_auc_p95": float(np.nanquantile(null_auc_array, 0.95)),
                "empirical_p_macro_auc": p_auc,
            }
        )

    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary["empirical_q_balanced_accuracy_bh"] = benjamini_hochberg(
            summary["empirical_p_balanced_accuracy"].tolist()
        )
        summary["empirical_q_macro_auc_bh"] = benjamini_hochberg(summary["empirical_p_macro_auc"].tolist())
    permutations = pd.DataFrame(permutation_rows)
    return summary, permutations


def plot_null_distributions(plt, sns, summary, permutations, asset_dir: Path) -> None:
    if summary.empty or permutations.empty:
        return
    plot_df = permutations.merge(
        summary[
            [
                "feature_view",
                "observed_repeated_cv_balanced_accuracy",
                "observed_repeated_cv_macro_auc_ovr",
            ]
        ],
        on="feature_view",
        how="left",
    )
    grid = sns.displot(
        data=plot_df,
        x="balanced_accuracy",
        col="feature_view_label",
        col_wrap=3,
        bins=18,
        color="#64748b",
        height=3.2,
        aspect=1.25,
    )
    for axis, (_view_label, group) in zip(grid.axes.flat, plot_df.groupby("feature_view_label", sort=False)):
        axis.axvline(group["observed_repeated_cv_balanced_accuracy"].iloc[0], color="#dc2626", linewidth=2)
        axis.axvline(0.5, color="#374151", linestyle="--", linewidth=1)
    grid.set_axis_labels("Balanced accuracy under shuffled labels", "Count")
    grid.set_titles("{col_name}")
    grid.fig.subplots_adjust(top=0.88)
    grid.fig.suptitle("HER2-Low vs HER2-Zero Classifier Permutation Null")
    grid.savefig(asset_dir / "classifier_permutation_balanced_accuracy_null.png", dpi=180)
    plt.close(grid.fig)

    grid = sns.displot(
        data=plot_df,
        x="macro_auc_ovr",
        col="feature_view_label",
        col_wrap=3,
        bins=18,
        color="#0f766e",
        height=3.2,
        aspect=1.25,
    )
    for axis, (_view_label, group) in zip(grid.axes.flat, plot_df.groupby("feature_view_label", sort=False)):
        axis.axvline(group["observed_repeated_cv_macro_auc_ovr"].iloc[0], color="#dc2626", linewidth=2)
        axis.axvline(0.5, color="#374151", linestyle="--", linewidth=1)
    grid.set_axis_labels("Macro AUC under shuffled labels", "Count")
    grid.set_titles("{col_name}")
    grid.fig.subplots_adjust(top=0.88)
    grid.fig.suptitle("HER2-Low vs HER2-Zero Classifier AUC Permutation Null")
    grid.savefig(asset_dir / "classifier_permutation_auc_null.png", dpi=180)
    plt.close(grid.fig)


def summary_table_rows(summary) -> list[list[str]]:
    rows = []
    for _, row in summary.iterrows():
        rows.append(
            [
                row["feature_view_label"],
                row["feature_set_label"],
                str(int(row["n_cases"])),
                str(int(row["n_features"])),
                fmt(row["loocv_balanced_accuracy"]),
                fmt(row["observed_repeated_cv_balanced_accuracy"]),
                fmt(row["null_balanced_accuracy_mean"]),
                fmt(row["null_balanced_accuracy_p95"]),
                fmt(row["empirical_p_balanced_accuracy"], 4),
                fmt(row["empirical_q_balanced_accuracy_bh"], 4),
                fmt(row["observed_repeated_cv_macro_auc_ovr"]),
                fmt(row["empirical_p_macro_auc"], 4),
            ]
        )
    return rows


def asset_link(asset_dir: Path, filename: str) -> str:
    return str(asset_dir / filename).replace("docs/", "")


def write_markdown(path: Path, asset_dir: Path, summary, args: argparse.Namespace) -> None:
    lines = [
        "# Classifier Permutation Sanity Check",
        "",
        "This analysis asks whether the selected HER2-low versus HER2-zero GigaTIME/H&E classifiers perform better than the same classifiers trained on shuffled labels.",
        "",
        "Important caveat: this is a post-hoc sanity check for the selected feature set in each view. It is not a fully nested model-selection permutation test, so it should be used as evidence that the signal is not obviously random, not as final clinical validation.",
        "",
        "Method:",
        "",
        f"- Task: HER2-low versus HER2-zero.",
        f"- Model: regularized logistic regression using the same selected GigaTIME/H&E feature set per view.",
        f"- Evaluation: repeated stratified {args.folds}-fold cross-validation with {args.repeats} repeats.",
        f"- Null: {args.permutations} label shuffles per view using the same folds and feature columns.",
        "",
        "## Results",
        "",
        markdown_table(
            [
                "View",
                "Feature set",
                "N",
                "Features",
                "LOOCV bal acc",
                "Repeated-CV bal acc",
                "Null mean",
                "Null 95%",
                "Empirical p",
                "BH q",
                "Repeated-CV AUC",
                "AUC p",
            ],
            summary_table_rows(summary),
        ),
        "",
        f"![Balanced accuracy permutation null]({asset_link(asset_dir, 'classifier_permutation_balanced_accuracy_null.png')})",
        "",
        f"![AUC permutation null]({asset_link(asset_dir, 'classifier_permutation_auc_null.png')})",
        "",
        "## Interpretation",
        "",
        "A low empirical p value means the observed repeated-CV classifier result is rarely matched by shuffled HER2-low/HER2-zero labels. This supports the idea that GigaTIME features contain real label-associated structure.",
        "",
        "This still does not make the model diagnostic. It does not validate real mIF biology, does not prove tumor-cell HER2 biology, and does not solve the tissue-composition caveat. It is a useful classifier trustworthiness check before discussing the result with an advisor.",
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{Path(args.out_dir) / 'classifier_permutation_summary.csv'}`",
        f"- `{Path(args.out_dir) / 'classifier_permutation_null_metrics.csv'}`",
        f"- `{asset_dir}/`",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    asset_dir = Path(args.asset_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, optimize, stats = require_analysis_libs(out_dir / ".matplotlib")

    summary, permutations = run_permutation_checks(pd, optimize, stats, args)
    summary.to_csv(out_dir / "classifier_permutation_summary.csv", index=False)
    permutations.to_csv(out_dir / "classifier_permutation_null_metrics.csv", index=False)
    metadata = {
        "task": "her2_low_vs_zero",
        "positive_class": POSITIVE_CLASS,
        "folds": args.folds,
        "repeats": args.repeats,
        "permutations": args.permutations,
        "seed": args.seed,
        "l2_penalty": args.l2_penalty,
    }
    (out_dir / "classifier_permutation_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    plot_null_distributions(plt, sns, summary, permutations, asset_dir)
    write_markdown(Path(args.out_markdown), asset_dir, summary, args)

    print(f"Wrote classifier permutation outputs to {out_dir}")
    print(f"Wrote classifier permutation figures to {asset_dir}")
    print(f"Wrote classifier permutation markdown to {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
