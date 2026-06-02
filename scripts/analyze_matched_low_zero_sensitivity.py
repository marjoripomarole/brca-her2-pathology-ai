#!/usr/bin/env python3
"""Matched HER2-low versus HER2-zero sensitivity after site/slide-size confounding."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np

from analyze_classifier_permutation_sanity import (
    NEGATIVE_CLASS,
    POSITIVE_CLASS,
    VIEW_ORDER,
    benjamini_hochberg,
    fmt,
    markdown_table,
    metric_dict,
)
from analyze_clinical_covariate_sensitivity import (
    KEY_CHANNELS,
    LOW_ZERO_GROUPS,
    VIEW_LABELS,
    add_covariates,
    one_hot_design,
    require_analysis_libs,
)
from train_her2_classifier_baseline import fit_predict_logistic, standardize_train_test
from train_her2_cleaned_classifier_comparison import GIGATIME_CHANNELS


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
SUBSET_ORDER = [
    "exact_source_site_nearest_size",
    "slide_size_caliper_0.25",
    "slide_size_caliper_0.50",
]
SUBSET_LABELS = {
    "exact_source_site_nearest_size": "Exact source-site, nearest size",
    "slide_size_caliper_0.25": "Slide-size matched, caliper 0.25",
    "slide_size_caliper_0.50": "Slide-size matched, caliper 0.50",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--slide-features",
        default=str(BASE_RESULT_DIR / "tumor_proxy_sensitivity/tumor_proxy_slide_features.csv"),
    )
    parser.add_argument(
        "--high-trust-slides",
        default="docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv",
    )
    parser.add_argument(
        "--out-dir",
        default=str(BASE_RESULT_DIR / "matched_low_zero_sensitivity"),
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_matched_low_zero",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_matched_low_zero_sensitivity.md",
    )
    parser.add_argument("--seed", type=int, default=20260602)
    parser.add_argument("--l2-penalty", type=float, default=1.0)
    parser.add_argument("--min-site-count", type=int, default=5)
    return parser.parse_args()


def low_zero_metadata(metadata):
    rows = metadata.loc[metadata["clinical_her2_group"].isin(LOW_ZERO_GROUPS)].copy()
    rows["tss_code"] = rows["case_submitter_id"].astype(str).str.split("-").str[1]
    rows["slide_file_size_mb"] = rows["slide_file_size_mb"].astype(float)
    rows["log_slide_file_size_mb"] = np.log1p(rows["slide_file_size_mb"])
    return rows.reset_index(drop=True)


def greedy_pairs(low_rows, zero_rows, subset: str, exact_site: bool = False, caliper: float | None = None):
    candidates = []
    for _, low in low_rows.iterrows():
        for _, zero in zero_rows.iterrows():
            if exact_site and low["tss_code"] != zero["tss_code"]:
                continue
            distance = abs(float(low["log_slide_file_size_mb"]) - float(zero["log_slide_file_size_mb"]))
            if caliper is not None and distance > caliper:
                continue
            candidates.append((distance, str(low["case_submitter_id"]), str(zero["case_submitter_id"]), low, zero))
    candidates.sort(key=lambda item: (item[0], item[1], item[2]))
    used_low = set()
    used_zero = set()
    pairs = []
    for distance, low_case, zero_case, low, zero in candidates:
        if low_case in used_low or zero_case in used_zero:
            continue
        used_low.add(low_case)
        used_zero.add(zero_case)
        pair_index = len(pairs) + 1
        pairs.append(
            {
                "matched_subset": subset,
                "matched_subset_label": SUBSET_LABELS[subset],
                "pair_id": f"{subset}_{pair_index:03d}",
                "low_case_submitter_id": low_case,
                "zero_case_submitter_id": zero_case,
                "low_slide_id": low["slide_id"],
                "zero_slide_id": zero["slide_id"],
                "low_tss_code": low["tss_code"],
                "zero_tss_code": zero["tss_code"],
                "same_tss": low["tss_code"] == zero["tss_code"],
                "low_slide_file_size_mb": float(low["slide_file_size_mb"]),
                "zero_slide_file_size_mb": float(zero["slide_file_size_mb"]),
                "abs_log_slide_size_diff": float(distance),
                "abs_slide_file_size_mb_diff": abs(float(low["slide_file_size_mb"]) - float(zero["slide_file_size_mb"])),
            }
        )
    return pairs


def build_all_pairs(pd, metadata):
    rows = low_zero_metadata(metadata)
    low_rows = rows.loc[rows["clinical_her2_group"] == NEGATIVE_CLASS]
    zero_rows = rows.loc[rows["clinical_her2_group"] == POSITIVE_CLASS]
    pairs = []
    pairs.extend(greedy_pairs(low_rows, zero_rows, "exact_source_site_nearest_size", exact_site=True))
    pairs.extend(greedy_pairs(low_rows, zero_rows, "slide_size_caliper_0.25", caliper=0.25))
    pairs.extend(greedy_pairs(low_rows, zero_rows, "slide_size_caliper_0.50", caliper=0.50))
    return pd.DataFrame(pairs)


def matched_feature_rows(pd, features, pairs):
    rows = []
    pair_lookup = {}
    for _, pair in pairs.iterrows():
        pair_lookup[(pair["matched_subset"], pair["low_case_submitter_id"])] = (pair["pair_id"], "HER2-low")
        pair_lookup[(pair["matched_subset"], pair["zero_case_submitter_id"])] = (pair["pair_id"], "HER2-zero")
    for subset in SUBSET_ORDER:
        subset_cases = {
            case_id
            for matched_subset, case_id in pair_lookup
            if matched_subset == subset
        }
        subset_features = features.loc[features["case_submitter_id"].isin(subset_cases)].copy()
        for _, row in subset_features.iterrows():
            key = (subset, row["case_submitter_id"])
            if key not in pair_lookup:
                continue
            pair_id, role = pair_lookup[key]
            record = row.to_dict()
            record["matched_subset"] = subset
            record["matched_subset_label"] = SUBSET_LABELS[subset]
            record["pair_id"] = pair_id
            record["pair_role"] = role
            rows.append(record)
    return pd.DataFrame(rows)


def design_matrix_for_feature_set(pd, rows, feature_set: str):
    image_cols = [f"mean_{channel}" for channel in GIGATIME_CHANNELS if f"mean_{channel}" in rows.columns]
    if feature_set == "slide_size_only":
        return one_hot_design(pd, rows, [], ["log_slide_file_size_mb", "log_slide_width", "log_slide_height"])
    if feature_set == "source_site_only":
        return one_hot_design(pd, rows, ["tss_group"], [])
    if feature_set == "source_site_slide_size":
        return one_hot_design(pd, rows, ["tss_group"], ["log_slide_file_size_mb", "log_slide_width", "log_slide_height"])
    if feature_set == "gigatime_mean_channels":
        return rows[image_cols].astype(float)
    if feature_set == "gigatime_plus_source_site_slide_size":
        return pd.concat(
            [
                rows[image_cols].astype(float),
                one_hot_design(pd, rows, ["tss_group"], ["log_slide_file_size_mb", "log_slide_width", "log_slide_height"]),
            ],
            axis=1,
        )
    raise ValueError(feature_set)


def leave_pair_out_metrics(pd, optimize, stats, rows, feature_set: str, l2_penalty: float):
    design = design_matrix_for_feature_set(pd, rows, feature_set)
    y = rows["clinical_her2_group"].map({NEGATIVE_CLASS: 0, POSITIVE_CLASS: 1}).to_numpy(dtype=int)
    x = design.to_numpy(dtype=float)
    classes = [NEGATIVE_CLASS, POSITIVE_CLASS]
    true_labels = []
    pred_labels = []
    positive_probs = []
    prediction_rows = []
    pair_ids = rows["pair_id"].drop_duplicates().tolist()
    for pair_id in pair_ids:
        test_idx = np.where(rows["pair_id"].to_numpy() == pair_id)[0]
        train_idx = np.where(rows["pair_id"].to_numpy() != pair_id)[0]
        if len(np.unique(y[train_idx])) < 2:
            continue
        x_train, x_test = standardize_train_test(x[train_idx], x[test_idx])
        probs = fit_predict_logistic(optimize, x_train, y[train_idx], x_test, classes, l2_penalty)
        pred = np.argmax(probs, axis=1).astype(int)
        true_labels.extend(y[test_idx].tolist())
        pred_labels.extend(pred.tolist())
        positive_probs.extend(probs[:, 1].astype(float).tolist())
        for local_idx, row_idx in enumerate(test_idx):
            source = rows.iloc[row_idx]
            prediction_rows.append(
                {
                    "matched_subset": source["matched_subset"],
                    "matched_subset_label": source["matched_subset_label"],
                    "feature_view": source["feature_view"],
                    "feature_view_label": source["feature_view_label"],
                    "feature_set": feature_set,
                    "pair_id": pair_id,
                    "case_submitter_id": source["case_submitter_id"],
                    "clinical_her2_group": source["clinical_her2_group"],
                    "true_label": classes[int(y[row_idx])],
                    "predicted_label": classes[int(pred[local_idx])],
                    "correct": bool(y[row_idx] == pred[local_idx]),
                    f"prob_{NEGATIVE_CLASS}": float(probs[local_idx, 0]),
                    f"prob_{POSITIVE_CLASS}": float(probs[local_idx, 1]),
                }
            )
    metrics = metric_dict(
        stats,
        np.array(true_labels, dtype=int),
        np.array(pred_labels, dtype=int),
        np.array(positive_probs, dtype=float),
    )
    return metrics, prediction_rows, design.shape[1]


def run_classifiers(pd, optimize, stats, matched_features, l2_penalty: float):
    feature_sets = [
        ("slide_size_only", "Slide-size covariates"),
        ("source_site_only", "Source-site covariates"),
        ("source_site_slide_size", "Source-site + slide-size"),
        ("gigatime_mean_channels", "GigaTIME mean channels"),
        ("gigatime_plus_source_site_slide_size", "GigaTIME + source-site/slide-size"),
    ]
    metric_rows = []
    prediction_rows = []
    for subset in SUBSET_ORDER:
        for view in VIEW_ORDER:
            rows = matched_features.loc[
                (matched_features["matched_subset"] == subset) & (matched_features["feature_view"] == view)
            ].copy()
            if rows.empty or rows["pair_id"].nunique() < 4:
                continue
            rows = rows.reset_index(drop=True)
            for feature_set, feature_set_label in feature_sets:
                metrics, predictions, n_features = leave_pair_out_metrics(pd, optimize, stats, rows, feature_set, l2_penalty)
                prediction_rows.extend(predictions)
                metric_rows.append(
                    {
                        "matched_subset": subset,
                        "matched_subset_label": SUBSET_LABELS[subset],
                        "feature_view": view,
                        "feature_view_label": VIEW_LABELS.get(view, view),
                        "feature_set": feature_set,
                        "feature_set_label": feature_set_label,
                        "n_pairs": rows["pair_id"].nunique(),
                        "n_cases": len(rows),
                        "n_features": n_features,
                        "accuracy": metrics["accuracy"],
                        "balanced_accuracy": metrics["balanced_accuracy"],
                        "macro_auc_ovr": metrics["macro_auc_ovr"],
                        "sensitivity": metrics["sensitivity"],
                        "specificity": metrics["specificity"],
                    }
                )
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows)


def paired_channel_tests(pd, stats, matched_features):
    rows = []
    for subset in SUBSET_ORDER:
        for view in VIEW_ORDER:
            view_rows = matched_features.loc[
                (matched_features["matched_subset"] == subset) & (matched_features["feature_view"] == view)
            ].copy()
            if view_rows.empty:
                continue
            p_values = []
            row_indices = []
            for channel in KEY_CHANNELS:
                col = f"mean_{channel}"
                if col not in view_rows.columns:
                    continue
                pair_diffs = []
                for pair_id, pair in view_rows.groupby("pair_id"):
                    if {NEGATIVE_CLASS, POSITIVE_CLASS} - set(pair["clinical_her2_group"]):
                        continue
                    low_value = float(pair.loc[pair["clinical_her2_group"] == NEGATIVE_CLASS, col].iloc[0])
                    zero_value = float(pair.loc[pair["clinical_her2_group"] == POSITIVE_CLASS, col].iloc[0])
                    pair_diffs.append(low_value - zero_value)
                if not pair_diffs:
                    continue
                p_value = math.nan
                if any(abs(value) > 1e-12 for value in pair_diffs):
                    p_value = float(stats.wilcoxon(pair_diffs, zero_method="wilcox", alternative="two-sided").pvalue)
                row_indices.append(len(rows))
                p_values.append(p_value)
                rows.append(
                    {
                        "matched_subset": subset,
                        "matched_subset_label": SUBSET_LABELS[subset],
                        "feature_view": view,
                        "feature_view_label": VIEW_LABELS.get(view, view),
                        "channel": channel,
                        "n_pairs": len(pair_diffs),
                        "mean_low_minus_zero": float(np.mean(pair_diffs)),
                        "median_low_minus_zero": float(np.median(pair_diffs)),
                        "wilcoxon_p_value": p_value,
                    }
                )
            q_values = benjamini_hochberg(p_values)
            for idx, q_value in zip(row_indices, q_values):
                rows[idx]["q_value_bh_within_subset_view"] = q_value
    return pd.DataFrame(rows)


def pair_summary(pd, pairs):
    rows = []
    for subset, group in pairs.groupby("matched_subset", sort=False):
        rows.append(
            {
                "matched_subset": subset,
                "matched_subset_label": SUBSET_LABELS[subset],
                "n_pairs": int(len(group)),
                "n_same_tss_pairs": int(group["same_tss"].sum()),
                "median_abs_log_slide_size_diff": float(group["abs_log_slide_size_diff"].median()),
                "mean_abs_log_slide_size_diff": float(group["abs_log_slide_size_diff"].mean()),
                "median_abs_slide_file_size_mb_diff": float(group["abs_slide_file_size_mb_diff"].median()),
                "mean_abs_slide_file_size_mb_diff": float(group["abs_slide_file_size_mb_diff"].mean()),
            }
        )
    return pd.DataFrame(rows)


def plot_pair_balance(plt, sns, pairs, asset_dir: Path) -> None:
    plot_rows = []
    for _, pair in pairs.iterrows():
        plot_rows.append(
            {
                "matched_subset_label": pair["matched_subset_label"],
                "pair_id": pair["pair_id"],
                "clinical_her2_group": NEGATIVE_CLASS,
                "slide_file_size_mb": pair["low_slide_file_size_mb"],
            }
        )
        plot_rows.append(
            {
                "matched_subset_label": pair["matched_subset_label"],
                "pair_id": pair["pair_id"],
                "clinical_her2_group": POSITIVE_CLASS,
                "slide_file_size_mb": pair["zero_slide_file_size_mb"],
            }
        )
    import pandas as pd

    plot_df = pd.DataFrame(plot_rows)
    plt.figure(figsize=(10.8, 5.6))
    sns.boxplot(data=plot_df, x="matched_subset_label", y="slide_file_size_mb", hue="clinical_her2_group")
    plt.xticks(rotation=20, ha="right")
    plt.xlabel("Matched sensitivity subset")
    plt.ylabel("Slide file size MB")
    plt.title("Slide-Size Balance After Matching")
    plt.tight_layout()
    plt.savefig(asset_dir / "matched_low_zero_slide_size_balance.png", dpi=180)
    plt.close()


def plot_classifier_metrics(plt, sns, classifier, asset_dir: Path) -> None:
    focus = classifier.loc[classifier["feature_view"] == "ck_top8_within_slide"].copy()
    plt.figure(figsize=(11.5, 6.0))
    sns.barplot(
        data=focus,
        x="matched_subset_label",
        y="balanced_accuracy",
        hue="feature_set_label",
    )
    plt.axhline(0.5, color="#374151", linestyle="--", linewidth=1)
    plt.ylim(0, 1)
    plt.xlabel("Matched sensitivity subset")
    plt.ylabel("Leave-pair-out balanced accuracy")
    plt.title("Matched HER2-Low vs HER2-Zero Classifier Sensitivity")
    plt.xticks(rotation=20, ha="right")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "matched_low_zero_classifier_sensitivity.png", dpi=180)
    plt.close()


def plot_channel_q_counts(plt, sns, channel_tests, asset_dir: Path) -> None:
    counts = (
        channel_tests.assign(significant=channel_tests["q_value_bh_within_subset_view"] < 0.05)
        .groupby(["matched_subset_label", "feature_view_label"])["significant"]
        .sum()
        .reset_index(name="n_q_lt_0_05_channels")
    )
    plt.figure(figsize=(11.5, 5.6))
    sns.barplot(data=counts, x="matched_subset_label", y="n_q_lt_0_05_channels", hue="feature_view_label")
    plt.xlabel("Matched sensitivity subset")
    plt.ylabel("Key channels with paired q < 0.05")
    plt.title("Paired Channel Signals After Matching")
    plt.xticks(rotation=20, ha="right")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "matched_low_zero_channel_q_counts.png", dpi=180)
    plt.close()


def table_pair_summary(summary) -> list[list[str]]:
    return [
        [
            row["matched_subset_label"],
            str(int(row["n_pairs"])),
            str(int(row["n_same_tss_pairs"])),
            fmt(row["median_abs_log_slide_size_diff"], 3),
            fmt(row["median_abs_slide_file_size_mb_diff"], 1),
        ]
        for _, row in summary.iterrows()
    ]


def classifier_table_rows(classifier, subset: str, view: str = "ck_top8_within_slide") -> list[list[str]]:
    selected = classifier.loc[(classifier["matched_subset"] == subset) & (classifier["feature_view"] == view)].copy()
    return [
        [
            row["feature_set_label"],
            str(int(row["n_pairs"])),
            str(int(row["n_features"])),
            fmt(row["balanced_accuracy"], 3),
            fmt(row["macro_auc_ovr"], 3),
        ]
        for _, row in selected.iterrows()
    ]


def channel_table_rows(channel_tests, subset: str, view: str = "ck_top8_within_slide") -> list[list[str]]:
    selected = channel_tests.loc[(channel_tests["matched_subset"] == subset) & (channel_tests["feature_view"] == view)].copy()
    selected = selected.sort_values("q_value_bh_within_subset_view")
    return [
        [
            row["channel"],
            str(int(row["n_pairs"])),
            fmt(row["mean_low_minus_zero"], 4),
            fmt(row["wilcoxon_p_value"], 4),
            fmt(row["q_value_bh_within_subset_view"], 4),
        ]
        for _, row in selected.head(8).iterrows()
    ]


def write_markdown(path: Path, asset_dir: Path, pair_summary_df, classifier, channel_tests) -> None:
    lines = [
        "# Matched HER2-Low Versus HER2-Zero Sensitivity",
        "",
        "This analysis responds to the clinical/source-site confounder finding. It builds matched HER2-low/HER2-zero subsets and reruns the key GigaTIME checks with leave-one-pair-out classifier evaluation and paired channel tests.",
        "",
        "Important caveat: exact source-site matching leaves only a small sample. These matched analyses are sensitivity checks, not final validation.",
        "",
        "## Matched Subsets",
        "",
        markdown_table(
            ["Matched subset", "Pairs", "Same-source-site pairs", "Median abs log-size diff", "Median abs MB diff"],
            table_pair_summary(pair_summary_df),
        ),
        "",
        f"![Slide-size balance]({str(asset_dir / 'matched_low_zero_slide_size_balance.png').replace('docs/', '')})",
        "",
        "## Leave-Pair-Out Classifier Sensitivity",
        "",
    ]
    for subset in SUBSET_ORDER:
        lines.extend(
            [
                f"### {SUBSET_LABELS[subset]}",
                "",
                markdown_table(
                    ["Feature set", "Pairs", "Features", "Balanced accuracy", "AUC"],
                    classifier_table_rows(classifier, subset),
                ),
                "",
            ]
        )
    lines.extend(
        [
            f"![Classifier sensitivity]({str(asset_dir / 'matched_low_zero_classifier_sensitivity.png').replace('docs/', '')})",
            "",
            "## Paired Channel Tests",
            "",
            "The table below shows the top paired low-minus-zero channel differences in the top 8 CK proxy view for each matched subset.",
            "",
        ]
    )
    for subset in SUBSET_ORDER:
        lines.extend(
            [
                f"### {SUBSET_LABELS[subset]}",
                "",
                markdown_table(
                    ["Channel", "Pairs", "Mean low-zero", "Wilcoxon p", "BH q"],
                    channel_table_rows(channel_tests, subset),
                ),
                "",
            ]
        )
    lines.extend(
        [
            f"![Paired channel q counts]({str(asset_dir / 'matched_low_zero_channel_q_counts.png').replace('docs/', '')})",
            "",
            "## Interpretation",
            "",
            "- GigaTIME mean channels remain modestly above chance in the matched subsets.",
            "- The confounder concern is not solved if source-site or slide-size baselines remain competitive or stronger than GigaTIME.",
            "- Paired channel tests should be interpreted cautiously; lack of BH q < 0.05 means the matched channel-level evidence is weak.",
            "- The safest conclusion is that the HER2-low/HER2-zero GigaTIME signal remains worth studying, but TCGA alone is not clean enough to support an independent HER2-biology or diagnostic claim.",
            "- Because exact source-site matching is small, any surviving signal should be treated as hypothesis-generating and should motivate external/site-balanced validation plus pathologist-reviewed tumor-rich tile analysis.",
            "",
            "## Output Files",
            "",
            f"- `{path}`",
            f"- `{BASE_RESULT_DIR / 'matched_low_zero_sensitivity/matched_pairs.csv'}`",
            f"- `{BASE_RESULT_DIR / 'matched_low_zero_sensitivity/matched_slide_features.csv'}`",
            f"- `{BASE_RESULT_DIR / 'matched_low_zero_sensitivity/matched_classifier_metrics.csv'}`",
            f"- `{BASE_RESULT_DIR / 'matched_low_zero_sensitivity/matched_channel_tests.csv'}`",
            f"- `{asset_dir}/`",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    asset_dir = Path(args.asset_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, optimize, stats = require_analysis_libs(out_dir / ".matplotlib")

    features = pd.read_csv(args.slide_features)
    metadata = pd.read_csv(args.high_trust_slides)
    features = add_covariates(pd, features, metadata, args.min_site_count)
    pairs = build_all_pairs(pd, metadata)
    matched_features = matched_feature_rows(pd, features, pairs)
    pair_summary_df = pair_summary(pd, pairs)
    classifier, predictions = run_classifiers(pd, optimize, stats, matched_features, args.l2_penalty)
    channel_tests = paired_channel_tests(pd, stats, matched_features)

    pairs.to_csv(out_dir / "matched_pairs.csv", index=False)
    matched_features.to_csv(out_dir / "matched_slide_features.csv", index=False)
    pair_summary_df.to_csv(out_dir / "matched_pair_summary.csv", index=False)
    classifier.to_csv(out_dir / "matched_classifier_metrics.csv", index=False)
    predictions.to_csv(out_dir / "matched_classifier_predictions.csv", index=False)
    channel_tests.to_csv(out_dir / "matched_channel_tests.csv", index=False)
    metadata_json = {
        "task": "her2_low_vs_zero",
        "subsets": SUBSET_LABELS,
        "seed": args.seed,
        "l2_penalty": args.l2_penalty,
        "min_site_count": args.min_site_count,
    }
    (out_dir / "matched_low_zero_sensitivity_metadata.json").write_text(
        json.dumps(metadata_json, indent=2) + "\n",
        encoding="utf-8",
    )

    plot_pair_balance(plt, sns, pairs, asset_dir)
    plot_classifier_metrics(plt, sns, classifier, asset_dir)
    plot_channel_q_counts(plt, sns, channel_tests, asset_dir)
    write_markdown(Path(args.out_markdown), asset_dir, pair_summary_df, classifier, channel_tests)

    print(f"Wrote matched low-zero sensitivity outputs to {out_dir}")
    print(f"Wrote matched low-zero sensitivity figures to {asset_dir}")
    print(f"Wrote matched low-zero sensitivity markdown to {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
