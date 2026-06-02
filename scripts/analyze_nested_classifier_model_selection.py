#!/usr/bin/env python3
"""Nested model-selection check for HER2-low versus HER2-zero classifiers."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np

from analyze_classifier_permutation_sanity import (
    LOW_ZERO_GROUPS,
    NEGATIVE_CLASS,
    POSITIVE_CLASS,
    VIEW_ORDER,
    benjamini_hochberg,
    fmt,
    make_repeated_stratified_folds,
    markdown_table,
    metric_dict,
    require_analysis_libs,
)
from train_her2_classifier_baseline import binary_auc, fit_predict_logistic, standardize_train_test
from train_her2_cleaned_classifier_comparison import FEATURE_LABELS, feature_sets_for


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
DEFAULT_TUMOR_PROXY_DIR = BASE_RESULT_DIR / "tumor_proxy_sensitivity"
VIEW_SORT_ORDER = {view: idx for idx, view in enumerate(VIEW_ORDER)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--slide-features",
        default=str(DEFAULT_TUMOR_PROXY_DIR / "tumor_proxy_slide_features.csv"),
        help="Slide-level tumor-proxy feature table.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(BASE_RESULT_DIR / "nested_classifier_model_selection"),
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_nested_classifier",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_nested_classifier_model_selection.md",
    )
    parser.add_argument("--outer-folds", type=int, default=5)
    parser.add_argument("--outer-repeats", type=int, default=3)
    parser.add_argument("--inner-folds", type=int, default=4)
    parser.add_argument("--permutations", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20260602)
    parser.add_argument("--l2-penalty", type=float, default=1.0)
    return parser.parse_args()


def candidate_feature_sets(rows) -> dict[str, list[str]]:
    return {
        name: cols
        for name, cols in feature_sets_for(rows).items()
        if name != "erbb2_rna_reference_not_h_e" and cols
    }


def evaluate_feature_set_cv(
    optimize,
    stats,
    x: np.ndarray,
    y: np.ndarray,
    folds,
    l2_penalty: float,
) -> dict[str, float]:
    true_labels = []
    pred_labels = []
    positive_probs = []
    classes = [NEGATIVE_CLASS, POSITIVE_CLASS]
    for train_idx, test_idx in folds:
        if len(np.unique(y[train_idx])) < 2:
            continue
        x_train, x_test = standardize_train_test(x[train_idx], x[test_idx])
        probs = fit_predict_logistic(optimize, x_train, y[train_idx], x_test, classes, l2_penalty)
        true_labels.extend(y[test_idx].tolist())
        pred_labels.extend(np.argmax(probs, axis=1).astype(int).tolist())
        positive_probs.extend(probs[:, 1].astype(float).tolist())
    return metric_dict(
        stats,
        np.array(true_labels, dtype=int),
        np.array(pred_labels, dtype=int),
        np.array(positive_probs, dtype=float),
    )


def choose_feature_set(
    optimize,
    stats,
    rows,
    feature_sets: dict[str, list[str]],
    train_idx: np.ndarray,
    y: np.ndarray,
    inner_folds: int,
    rng: np.random.Generator,
    l2_penalty: float,
) -> tuple[str, dict[str, float], list[dict[str, object]]]:
    y_train = y[train_idx]
    inner_cv = make_repeated_stratified_folds(y_train, inner_folds, 1, rng)
    candidate_rows = []
    for feature_set_name, feature_cols in feature_sets.items():
        x_train = rows.iloc[train_idx][feature_cols].to_numpy(dtype=float)
        metrics = evaluate_feature_set_cv(optimize, stats, x_train, y_train, inner_cv, l2_penalty)
        candidate_rows.append(
            {
                "feature_set": feature_set_name,
                "feature_set_label": FEATURE_LABELS.get(feature_set_name, feature_set_name),
                "n_features": len(feature_cols),
                **metrics,
            }
        )
    candidate_rows.sort(
        key=lambda row: (
            row["balanced_accuracy"],
            row["macro_auc_ovr"] if row["macro_auc_ovr"] == row["macro_auc_ovr"] else -1.0,
            -row["n_features"],
        ),
        reverse=True,
    )
    winner = candidate_rows[0]
    return str(winner["feature_set"]), winner, candidate_rows


def run_nested_cv_for_view(
    optimize,
    stats,
    rows,
    feature_sets: dict[str, list[str]],
    y: np.ndarray,
    outer_folds: int,
    outer_repeats: int,
    inner_folds: int,
    rng: np.random.Generator,
    l2_penalty: float,
    view: str,
    permutation_index: int | None = None,
):
    folds = make_repeated_stratified_folds(y, outer_folds, outer_repeats, rng)
    classes = [NEGATIVE_CLASS, POSITIVE_CLASS]
    prediction_rows = []
    selection_rows = []
    all_true = []
    all_pred = []
    all_positive_probs = []

    for fold_idx, (train_idx, test_idx) in enumerate(folds):
        selected_feature_set, inner_winner, candidate_rows = choose_feature_set(
            optimize,
            stats,
            rows,
            feature_sets,
            train_idx,
            y,
            inner_folds,
            rng,
            l2_penalty,
        )
        feature_cols = feature_sets[selected_feature_set]
        x_train = rows.iloc[train_idx][feature_cols].to_numpy(dtype=float)
        x_test = rows.iloc[test_idx][feature_cols].to_numpy(dtype=float)
        x_train_std, x_test_std = standardize_train_test(x_train, x_test)
        probs = fit_predict_logistic(optimize, x_train_std, y[train_idx], x_test_std, classes, l2_penalty)
        pred = np.argmax(probs, axis=1).astype(int)

        selection_rows.append(
            {
                "feature_view": view,
                "feature_view_label": rows["feature_view_label"].iloc[0],
                "permutation_index": permutation_index,
                "outer_fold_index": fold_idx,
                "selected_feature_set": selected_feature_set,
                "selected_feature_set_label": FEATURE_LABELS.get(selected_feature_set, selected_feature_set),
                "selected_n_features": len(feature_cols),
                "inner_balanced_accuracy": inner_winner["balanced_accuracy"],
                "inner_macro_auc_ovr": inner_winner["macro_auc_ovr"],
                "n_candidate_feature_sets": len(candidate_rows),
            }
        )

        for local_offset, sample_idx in enumerate(test_idx):
            row = rows.iloc[sample_idx]
            true_idx = int(y[sample_idx])
            pred_idx = int(pred[local_offset])
            prediction_rows.append(
                {
                    "feature_view": view,
                    "feature_view_label": row["feature_view_label"],
                    "permutation_index": permutation_index,
                    "outer_fold_index": fold_idx,
                    "case_submitter_id": row["case_submitter_id"],
                    "clinical_her2_group": row["clinical_her2_group"],
                    "true_label": classes[true_idx],
                    "predicted_label": classes[pred_idx],
                    "correct": true_idx == pred_idx,
                    "selected_feature_set": selected_feature_set,
                    "selected_feature_set_label": FEATURE_LABELS.get(selected_feature_set, selected_feature_set),
                    f"prob_{NEGATIVE_CLASS}": float(probs[local_offset, 0]),
                    f"prob_{POSITIVE_CLASS}": float(probs[local_offset, 1]),
                }
            )
        all_true.extend(y[test_idx].tolist())
        all_pred.extend(pred.tolist())
        all_positive_probs.extend(probs[:, 1].astype(float).tolist())

    metrics = metric_dict(
        stats,
        np.array(all_true, dtype=int),
        np.array(all_pred, dtype=int),
        np.array(all_positive_probs, dtype=float),
    )
    metrics["n_cv_predictions"] = len(all_true)
    metrics["n_outer_folds"] = len(folds)
    return metrics, prediction_rows, selection_rows


def run_nested_analysis(pd, optimize, stats, args: argparse.Namespace):
    rng = np.random.default_rng(args.seed)
    slide_features = pd.read_csv(args.slide_features)
    summary_rows = []
    prediction_rows = []
    selection_rows = []
    null_rows = []

    for view in VIEW_ORDER:
        rows = slide_features.loc[
            (slide_features["feature_view"] == view)
            & (slide_features["clinical_her2_group"].isin(LOW_ZERO_GROUPS))
        ].copy()
        if rows.empty:
            continue
        rows = rows.reset_index(drop=True)
        feature_sets = candidate_feature_sets(rows)
        if not feature_sets:
            continue
        y = rows["clinical_her2_group"].map({NEGATIVE_CLASS: 0, POSITIVE_CLASS: 1}).to_numpy(dtype=int)
        observed, view_predictions, view_selections = run_nested_cv_for_view(
            optimize,
            stats,
            rows,
            feature_sets,
            y,
            args.outer_folds,
            args.outer_repeats,
            args.inner_folds,
            rng,
            args.l2_penalty,
            view,
        )
        prediction_rows.extend(view_predictions)
        selection_rows.extend(view_selections)

        null_balanced = []
        null_auc = []
        for permutation_idx in range(args.permutations):
            y_perm = y.copy()
            rng.shuffle(y_perm)
            null_metrics, _null_predictions, null_selections = run_nested_cv_for_view(
                optimize,
                stats,
                rows,
                feature_sets,
                y_perm,
                args.outer_folds,
                args.outer_repeats,
                args.inner_folds,
                rng,
                args.l2_penalty,
                view,
                permutation_index=permutation_idx,
            )
            null_balanced.append(null_metrics["balanced_accuracy"])
            null_auc.append(null_metrics["macro_auc_ovr"])
            null_rows.append(
                {
                    "feature_view": view,
                    "feature_view_label": rows["feature_view_label"].iloc[0],
                    "permutation_index": permutation_idx,
                    "nested_balanced_accuracy": null_metrics["balanced_accuracy"],
                    "nested_macro_auc_ovr": null_metrics["macro_auc_ovr"],
                    "most_selected_feature_set": most_common_selected_feature(null_selections),
                }
            )

        null_balanced_array = np.array(null_balanced, dtype=float)
        null_auc_array = np.array(null_auc, dtype=float)
        selected_counts = Counter(row["selected_feature_set"] for row in view_selections)
        most_selected = selected_counts.most_common(1)[0][0]
        summary_rows.append(
            {
                "feature_view": view,
                "feature_view_label": rows["feature_view_label"].iloc[0],
                "n_cases": len(rows),
                "n_feature_sets": len(feature_sets),
                "outer_folds": args.outer_folds,
                "outer_repeats": args.outer_repeats,
                "inner_folds": args.inner_folds,
                "n_permutations": args.permutations,
                "observed_nested_accuracy": observed["accuracy"],
                "observed_nested_balanced_accuracy": observed["balanced_accuracy"],
                "observed_nested_macro_auc_ovr": observed["macro_auc_ovr"],
                "observed_nested_sensitivity": observed["sensitivity"],
                "observed_nested_specificity": observed["specificity"],
                "n_cv_predictions": observed["n_cv_predictions"],
                "n_outer_folds_total": observed["n_outer_folds"],
                "most_selected_feature_set": most_selected,
                "most_selected_feature_set_label": FEATURE_LABELS.get(most_selected, most_selected),
                "most_selected_count": selected_counts[most_selected],
                "selection_counts_json": json.dumps(dict(sorted(selected_counts.items())), sort_keys=True),
                "null_balanced_accuracy_mean": float(np.nanmean(null_balanced_array)),
                "null_balanced_accuracy_sd": float(np.nanstd(null_balanced_array, ddof=1)),
                "null_balanced_accuracy_p95": float(np.nanquantile(null_balanced_array, 0.95)),
                "empirical_p_balanced_accuracy": float(
                    (1 + np.sum(null_balanced_array >= observed["balanced_accuracy"])) / (1 + args.permutations)
                ),
                "null_macro_auc_mean": float(np.nanmean(null_auc_array)),
                "null_macro_auc_sd": float(np.nanstd(null_auc_array, ddof=1)),
                "null_macro_auc_p95": float(np.nanquantile(null_auc_array, 0.95)),
                "empirical_p_macro_auc": float(
                    (1 + np.sum(null_auc_array >= observed["macro_auc_ovr"])) / (1 + args.permutations)
                ),
            }
        )

    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary["empirical_q_balanced_accuracy_bh"] = benjamini_hochberg(
            summary["empirical_p_balanced_accuracy"].tolist()
        )
        summary["empirical_q_macro_auc_bh"] = benjamini_hochberg(summary["empirical_p_macro_auc"].tolist())
    predictions = pd.DataFrame(prediction_rows)
    selections = pd.DataFrame(selection_rows)
    null_metrics = pd.DataFrame(null_rows)
    return summary, predictions, selections, null_metrics


def most_common_selected_feature(selection_rows: list[dict[str, object]]) -> str:
    if not selection_rows:
        return ""
    counts = Counter(str(row["selected_feature_set"]) for row in selection_rows)
    return counts.most_common(1)[0][0]


def plot_nested_vs_null(plt, sns, summary, null_metrics, asset_dir: Path) -> None:
    if summary.empty or null_metrics.empty:
        return
    plot_df = summary.copy()
    plot_df["feature_view_label"] = plot_df["feature_view_label"].astype(str)
    plot_df = plot_df.sort_values(
        "feature_view", key=lambda values: values.map(lambda value: VIEW_SORT_ORDER.get(value, 99))
    )
    plt.figure(figsize=(11.0, 5.4))
    x = np.arange(len(plot_df))
    plt.bar(x - 0.18, plot_df["null_balanced_accuracy_p95"], width=0.36, color="#cbd5e1", label="Shuffled-label 95%")
    plt.bar(
        x + 0.18,
        plot_df["observed_nested_balanced_accuracy"],
        width=0.36,
        color="#0f766e",
        label="Observed nested CV",
    )
    plt.axhline(0.5, color="#374151", linestyle="--", linewidth=1)
    plt.xticks(x, plot_df["feature_view_label"], rotation=25, ha="right")
    plt.ylabel("Balanced accuracy")
    plt.xlabel("GigaTIME proxy view")
    plt.title("Nested Feature-Selection Classifier Versus Shuffled-Label Null")
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()
    plt.savefig(asset_dir / "nested_classifier_observed_vs_null.png", dpi=180)
    plt.close()

    grid = sns.displot(
        data=null_metrics,
        x="nested_balanced_accuracy",
        col="feature_view_label",
        col_wrap=3,
        bins=12,
        color="#64748b",
        height=3.2,
        aspect=1.25,
    )
    observed_lookup = dict(zip(summary["feature_view_label"], summary["observed_nested_balanced_accuracy"]))
    for axis, view_label in zip(grid.axes.flat, [axis.get_title().replace("feature_view_label = ", "") for axis in grid.axes.flat]):
        if view_label in observed_lookup:
            axis.axvline(observed_lookup[view_label], color="#dc2626", linewidth=2)
        axis.axvline(0.5, color="#374151", linestyle="--", linewidth=1)
    grid.set_axis_labels("Nested balanced accuracy under shuffled labels", "Count")
    grid.set_titles("{col_name}")
    grid.fig.subplots_adjust(top=0.88)
    grid.fig.suptitle("Nested HER2-Low vs HER2-Zero Shuffled-Label Null")
    grid.savefig(asset_dir / "nested_classifier_balanced_accuracy_null.png", dpi=180)
    plt.close(grid.fig)


def plot_selection_frequency(plt, sns, selections, asset_dir: Path) -> None:
    if selections.empty:
        return
    counts = (
        selections.groupby(["feature_view_label", "selected_feature_set_label"])
        .size()
        .reset_index(name="n_outer_folds")
    )
    plt.figure(figsize=(11.0, 5.5))
    sns.barplot(
        data=counts,
        x="feature_view_label",
        y="n_outer_folds",
        hue="selected_feature_set_label",
    )
    plt.xticks(rotation=25, ha="right")
    plt.xlabel("GigaTIME proxy view")
    plt.ylabel("Outer folds selected")
    plt.title("Feature Set Selected Inside Nested CV")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "nested_classifier_feature_selection_frequency.png", dpi=180)
    plt.close()


def summary_rows(summary) -> list[list[str]]:
    rows = []
    for _, row in summary.iterrows():
        rows.append(
            [
                row["feature_view_label"],
                str(int(row["n_cases"])),
                row["most_selected_feature_set_label"],
                str(int(row["most_selected_count"])),
                fmt(row["observed_nested_balanced_accuracy"]),
                fmt(row["observed_nested_macro_auc_ovr"]),
                fmt(row["null_balanced_accuracy_mean"]),
                fmt(row["null_balanced_accuracy_p95"]),
                fmt(row["empirical_p_balanced_accuracy"], 4),
                fmt(row["empirical_q_balanced_accuracy_bh"], 4),
            ]
        )
    return rows


def asset_link(asset_dir: Path, filename: str) -> str:
    return str(asset_dir / filename).replace("docs/", "")


def write_markdown(path: Path, asset_dir: Path, summary, args: argparse.Namespace) -> None:
    lines = [
        "# Nested Classifier Model-Selection Check",
        "",
        "This analysis is a stricter follow-up to the post-hoc classifier permutation check. Instead of fixing the previously selected feature set, it chooses the best GigaTIME/H&E feature set inside each training fold and only then evaluates the held-out fold.",
        "",
        "Method:",
        "",
        "- Task: HER2-low versus HER2-zero.",
        "- Model: regularized logistic regression.",
        f"- Outer evaluation: repeated stratified {args.outer_folds}-fold cross-validation with {args.outer_repeats} repeats.",
        f"- Inner model selection: stratified {args.inner_folds}-fold cross-validation on the outer-training set only.",
        "- Candidate feature sets: GigaTIME mean channels, mean+fraction channels, interpretable marker means, interpretable distribution features, and virtual programs when available.",
        f"- Null: {args.permutations} shuffled-label runs per view, with feature-set selection repeated inside each shuffled run.",
        "",
        "Important caveat: this is still not external clinical validation. It is a stronger internal sanity check that reduces feature-set selection bias.",
        "",
        "## Results",
        "",
        markdown_table(
            [
                "View",
                "N",
                "Most selected feature set",
                "Selected folds",
                "Nested bal acc",
                "Nested AUC",
                "Null mean",
                "Null 95%",
                "Empirical p",
                "BH q",
            ],
            summary_rows(summary),
        ),
        "",
        f"![Nested classifier observed versus null]({asset_link(asset_dir, 'nested_classifier_observed_vs_null.png')})",
        "",
        f"![Nested classifier null distributions]({asset_link(asset_dir, 'nested_classifier_balanced_accuracy_null.png')})",
        "",
        f"![Nested classifier feature selection frequency]({asset_link(asset_dir, 'nested_classifier_feature_selection_frequency.png')})",
        "",
        "## Interpretation",
        "",
        "If the nested result remains above the shuffled-label null, the low-versus-zero classifier signal is less likely to be a simple artifact of choosing the best feature set after looking at the whole data set.",
        "",
        "This still does not prove clinical diagnosis, real mIF validity, HER2 isoform biology, or treatment-response biology. It is classifier-methodology evidence that the next validation step should focus on tumor-rich/pathologist-approved regions and external molecular/protein validation.",
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{Path(args.out_dir) / 'nested_classifier_summary.csv'}`",
        f"- `{Path(args.out_dir) / 'nested_classifier_predictions.csv'}`",
        f"- `{Path(args.out_dir) / 'nested_classifier_feature_selection.csv'}`",
        f"- `{Path(args.out_dir) / 'nested_classifier_null_metrics.csv'}`",
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

    summary, predictions, selections, null_metrics = run_nested_analysis(pd, optimize, stats, args)
    summary.to_csv(out_dir / "nested_classifier_summary.csv", index=False)
    predictions.to_csv(out_dir / "nested_classifier_predictions.csv", index=False)
    selections.to_csv(out_dir / "nested_classifier_feature_selection.csv", index=False)
    null_metrics.to_csv(out_dir / "nested_classifier_null_metrics.csv", index=False)
    metadata = {
        "task": "her2_low_vs_zero",
        "positive_class": POSITIVE_CLASS,
        "outer_folds": args.outer_folds,
        "outer_repeats": args.outer_repeats,
        "inner_folds": args.inner_folds,
        "permutations": args.permutations,
        "seed": args.seed,
        "l2_penalty": args.l2_penalty,
    }
    (out_dir / "nested_classifier_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    plot_nested_vs_null(plt, sns, summary, null_metrics, asset_dir)
    plot_selection_frequency(plt, sns, selections, asset_dir)
    write_markdown(Path(args.out_markdown), asset_dir, summary, args)

    print(f"Wrote nested classifier outputs to {out_dir}")
    print(f"Wrote nested classifier figures to {asset_dir}")
    print(f"Wrote nested classifier markdown to {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
