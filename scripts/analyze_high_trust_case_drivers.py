#!/usr/bin/env python3
"""Rank case-level drivers of the strict HER2-low versus HER2-zero signal."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
LOW_ZERO_GROUPS = ["HER2-low", "HER2-zero"]
FEATURE_VIEW_ORDER = ["all_sampled_tissue", "qc_cellular_tissue", "ck_enriched_top50", "ck_enriched_top25"]
FEATURE_VIEW_LABELS = {
    "all_sampled_tissue": "All sampled tissue",
    "qc_cellular_tissue": "QC cellular tissue",
    "ck_enriched_top50": "CK-enriched top 50%",
    "ck_enriched_top25": "CK-enriched top 25%",
}
METADATA_COLUMNS = [
    "label_slide_trust",
    "processed_tissue_qc",
    "trust_reasons",
    "her2_detail_subgroup",
    "her2_ihc_score",
    "her2_ish_status",
    "er_status",
    "pr_status",
    "patient_gender",
    "histological_type",
    "pathologic_stage",
    "history_neoadjuvant_treatment",
    "slide_file_name",
    "slide_local_path",
    "slide_file_size_mb",
    "slide_width",
    "slide_height",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cleaned-features",
        default=str(BASE_RESULT_DIR / "gigatime_cleanup/cleaned_slide_features.csv"),
    )
    parser.add_argument(
        "--pairwise",
        default=str(BASE_RESULT_DIR / "gigatime_cleanup/cleanup_pairwise_tests.csv"),
    )
    parser.add_argument(
        "--classifier-predictions",
        default=str(BASE_RESULT_DIR / "cleaned_classifier_comparison/cleaned_classifier_predictions.csv"),
    )
    parser.add_argument(
        "--best-classifier-metrics",
        default=str(BASE_RESULT_DIR / "cleaned_classifier_comparison/cleaned_classifier_best_h_e_metrics.csv"),
    )
    parser.add_argument(
        "--high-trust-slides",
        default="docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv",
    )
    parser.add_argument(
        "--out-dir",
        default=str(BASE_RESULT_DIR / "case_driver_analysis"),
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_case_drivers",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_case_driver_analysis.md",
    )
    parser.add_argument("--q-cutoff", type=float, default=0.05)
    return parser.parse_args()


def require_analysis_libs(mpl_config_dir: Path):
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    try:
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError as exc:
        raise SystemExit(
            "Missing analysis dependency. Run this inside the project conda environment."
        ) from exc
    sns.set_theme(style="whitegrid", context="notebook")
    return np, pd, plt, sns


def fmt(value, digits: int = 3) -> str:
    if value is None:
        return ""
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
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


def asset_link(asset_dir: Path, image_name: str) -> str:
    return str(asset_dir / image_name).replace("docs/", "")


def z_col(channel: str) -> str:
    return "zero_like_z_" + channel.replace("-", "_").replace(" ", "_")


def load_inputs(pd, args: argparse.Namespace):
    features = pd.read_csv(args.cleaned_features)
    pairwise = pd.read_csv(args.pairwise)
    predictions = pd.read_csv(args.classifier_predictions)
    best_metrics = pd.read_csv(args.best_classifier_metrics)
    cohort = pd.read_csv(args.high_trust_slides)

    metadata = ["slide_id", "case_submitter_id"] + [
        col for col in METADATA_COLUMNS if col in cohort.columns and col not in features.columns
    ]
    merged = features.merge(
        cohort[metadata],
        on=["slide_id", "case_submitter_id"],
        how="left",
        validate="many_to_one",
    )
    return merged, pairwise, predictions, best_metrics


def build_signal_channels(pd, pairwise, q_cutoff: float):
    low_zero = pairwise.loc[
        pairwise["group_a"].isin(LOW_ZERO_GROUPS) & pairwise["group_b"].isin(LOW_ZERO_GROUPS)
    ].copy()
    rows = []
    for _, row in low_zero.iterrows():
        if row["group_a"] == "HER2-low":
            mean_low = row["mean_a"]
            mean_zero = row["mean_b"]
        else:
            mean_low = row["mean_b"]
            mean_zero = row["mean_a"]
        delta = float(mean_low) - float(mean_zero)
        q_value = float(row["mannwhitney_q_value_bh_within_view"])
        if q_value > q_cutoff or delta == 0:
            continue
        rows.append(
            {
                "feature_view": row["feature_view"],
                "feature_view_label": row["feature_view_label"],
                "channel": row["channel"],
                "mean_her2_low": float(mean_low),
                "mean_her2_zero": float(mean_zero),
                "delta_low_minus_zero": delta,
                "zero_higher_direction": 1 if mean_zero > mean_low else -1,
                "mannwhitney_q_value_bh_within_view": q_value,
            }
        )
    output = pd.DataFrame(rows)
    if not output.empty:
        output["view_order"] = output["feature_view"].map({view: idx for idx, view in enumerate(FEATURE_VIEW_ORDER)})
        output = output.sort_values(["view_order", "mannwhitney_q_value_bh_within_view", "channel"])
    return output


def build_case_driver_scores(np, pd, features, signal_channels):
    base = features.loc[features["clinical_her2_group"].isin(LOW_ZERO_GROUPS)].copy()
    scored_views = []
    long_rows = []

    for view in FEATURE_VIEW_ORDER:
        view_rows = base.loc[base["feature_view"] == view].copy()
        view_signals = signal_channels.loc[signal_channels["feature_view"] == view].copy()
        z_cols = []
        for _, signal in view_signals.iterrows():
            channel = signal["channel"]
            column = f"mean_{channel}"
            if column not in view_rows.columns:
                continue
            values = view_rows[column].astype(float)
            center = float(values.mean())
            scale = float(values.std(ddof=0))
            if scale == 0 or scale != scale:
                continue
            direction = float(signal["zero_higher_direction"])
            score_column = z_col(channel)
            view_rows[score_column] = ((values - center) / scale) * direction
            z_cols.append(score_column)
            for _, case_row in view_rows.iterrows():
                long_rows.append(
                    {
                        "feature_view": view,
                        "feature_view_label": FEATURE_VIEW_LABELS.get(view, view),
                        "slide_id": case_row["slide_id"],
                        "case_submitter_id": case_row["case_submitter_id"],
                        "clinical_her2_group": case_row["clinical_her2_group"],
                        "channel": channel,
                        "raw_mean_channel_value": float(case_row[column]),
                        "zero_like_channel_z": float(case_row[score_column]),
                        "delta_low_minus_zero": float(signal["delta_low_minus_zero"]),
                        "mannwhitney_q_value_bh_within_view": float(
                            signal["mannwhitney_q_value_bh_within_view"]
                        ),
                    }
                )
        if not z_cols:
            continue
        view_rows["zero_like_score"] = view_rows[z_cols].mean(axis=1)
        view_rows["n_signal_channels"] = len(z_cols)
        view_rows["expected_profile_score"] = np.where(
            view_rows["clinical_her2_group"] == "HER2-zero",
            view_rows["zero_like_score"],
            -view_rows["zero_like_score"],
        )
        view_rows["opposite_profile"] = view_rows["expected_profile_score"] < 0
        view_rows["zero_like_score_percentile"] = view_rows["zero_like_score"].rank(pct=True, method="average")
        view_rows["expected_driver_rank_within_group"] = view_rows.groupby("clinical_her2_group")[
            "expected_profile_score"
        ].rank(method="first", ascending=False)
        scored_views.append(view_rows)

    driver_scores = pd.concat(scored_views, ignore_index=True) if scored_views else pd.DataFrame()
    channel_scores = pd.DataFrame(long_rows)
    if not driver_scores.empty:
        driver_scores["view_order"] = driver_scores["feature_view"].map(
            {view: idx for idx, view in enumerate(FEATURE_VIEW_ORDER)}
        )
        driver_scores = driver_scores.sort_values(
            ["view_order", "clinical_her2_group", "expected_driver_rank_within_group"]
        )
    return driver_scores, channel_scores


def build_view_stability(pd, driver_scores):
    group_columns = ["slide_id", "case_submitter_id", "clinical_her2_group"]
    metadata_columns = [
        col
        for col in [
            "her2_detail_subgroup",
            "her2_ihc_score",
            "her2_ish_status",
            "er_status",
            "pr_status",
            "erbb2_tpm",
            "histological_type",
            "pathologic_stage",
            "processed_tissue_qc",
        ]
        if col in driver_scores.columns
    ]
    rows = []
    for keys, group in driver_scores.groupby(group_columns, sort=False):
        record = dict(zip(group_columns, keys))
        for column in metadata_columns:
            record[column] = group[column].dropna().iloc[0] if not group[column].dropna().empty else ""
        record["views_available"] = int(group["feature_view"].nunique())
        record["expected_profile_views"] = int((~group["opposite_profile"]).sum())
        record["opposite_profile_views"] = int(group["opposite_profile"].sum())
        record["mean_expected_profile_score"] = float(group["expected_profile_score"].mean())
        record["min_expected_profile_score"] = float(group["expected_profile_score"].min())
        record["mean_zero_like_score"] = float(group["zero_like_score"].mean())
        record["zero_like_score_range"] = float(group["zero_like_score"].max() - group["zero_like_score"].min())
        record["zero_like_score_sd"] = float(group["zero_like_score"].std(ddof=0))
        for view in FEATURE_VIEW_ORDER:
            view_group = group.loc[group["feature_view"] == view]
            record[f"{view}_zero_like_score"] = (
                float(view_group["zero_like_score"].iloc[0]) if not view_group.empty else float("nan")
            )
            record[f"{view}_expected_profile_score"] = (
                float(view_group["expected_profile_score"].iloc[0]) if not view_group.empty else float("nan")
            )
        if record["opposite_profile_views"] >= 3:
            record["stability_category"] = "opposite_profile_in_most_views"
        elif record["opposite_profile_views"] >= 2:
            record["stability_category"] = "mixed_or_view_sensitive"
        elif record["zero_like_score_range"] >= 2.0:
            record["stability_category"] = "large_view_shift"
        else:
            record["stability_category"] = "mostly_label_consistent"
        rows.append(record)
    output = pd.DataFrame(rows)
    if not output.empty:
        output = output.sort_values(
            ["opposite_profile_views", "zero_like_score_range", "mean_expected_profile_score"],
            ascending=[False, False, True],
        )
    return output


def build_classifier_review(pd, predictions, best_metrics, view_stability):
    best = best_metrics.loc[
        (best_metrics["task"] == "her2_low_vs_zero") & (best_metrics["model"] == "regularized_logistic")
    ][["feature_view", "feature_set"]].copy()
    preds = predictions.loc[
        (predictions["task"] == "her2_low_vs_zero") & (predictions["model"] == "regularized_logistic")
    ].merge(best, on=["feature_view", "feature_set"], how="inner", validate="many_to_one")
    preds["prob_true_label"] = preds.apply(lambda row: row.get(f"prob_{row['true_label']}", float("nan")), axis=1)
    preds["prob_predicted_label"] = preds.apply(
        lambda row: row.get(f"prob_{row['predicted_label']}", float("nan")), axis=1
    )
    preds["probability_margin"] = (preds["prob_HER2-zero"] - preds["prob_HER2-low"]).abs()
    preds["correct"] = preds["correct"].astype(str).str.lower().eq("true")
    preds["view_order"] = preds["feature_view"].map({view: idx for idx, view in enumerate(FEATURE_VIEW_ORDER)})
    preds = preds.sort_values(["case_submitter_id", "view_order"])

    summary_rows = []
    for (case_id, group), case_group in preds.groupby(["case_submitter_id", "clinical_her2_group"], sort=False):
        record = {
            "case_submitter_id": case_id,
            "clinical_her2_group": group,
            "classifier_views_available": int(case_group["feature_view"].nunique()),
            "classifier_correct_views": int(case_group["correct"].sum()),
            "classifier_incorrect_views": int((~case_group["correct"]).sum()),
            "mean_prob_true_label": float(case_group["prob_true_label"].mean()),
            "min_prob_true_label": float(case_group["prob_true_label"].min()),
            "mean_probability_margin": float(case_group["probability_margin"].mean()),
            "predictions_by_view": "; ".join(
                f"{row.feature_view_label}: {row.predicted_label}"
                for row in case_group.itertuples(index=False)
            ),
        }
        summary_rows.append(record)
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        merge_cols = [
            "slide_id",
            "case_submitter_id",
            "her2_detail_subgroup",
            "er_status",
            "pr_status",
            "erbb2_tpm",
            "opposite_profile_views",
            "expected_profile_views",
            "zero_like_score_range",
            "all_sampled_tissue_zero_like_score",
            "stability_category",
        ]
        merge_cols = [col for col in merge_cols if col in view_stability.columns]
        summary = summary.merge(view_stability[merge_cols], on="case_submitter_id", how="left")
        summary["manual_review_priority_score"] = (
            summary["classifier_incorrect_views"].fillna(0)
            + summary.get("opposite_profile_views", 0).fillna(0)
            + (summary.get("zero_like_score_range", 0).fillna(0) >= 2.0).astype(int)
        )
        summary = summary.sort_values(
            ["manual_review_priority_score", "classifier_incorrect_views", "opposite_profile_views"],
            ascending=[False, False, False],
        )
    return preds, summary


def summarize_counts(pd, signal_channels, driver_scores, view_stability, classifier_best_predictions):
    signal_counts = (
        signal_channels.groupby(["feature_view", "feature_view_label"], as_index=False)
        .agg(n_signal_channels=("channel", "nunique"), channels=("channel", lambda values: ", ".join(values)))
        .sort_values("feature_view", key=lambda col: col.map({view: idx for idx, view in enumerate(FEATURE_VIEW_ORDER)}))
    )
    stability_counts = (
        view_stability.groupby(["clinical_her2_group", "expected_profile_views"], as_index=False)
        .size()
        .rename(columns={"size": "n_slides"})
    )
    classifier_counts = (
        classifier_best_predictions.groupby(["feature_view", "feature_view_label"], as_index=False)
        .agg(
            n_cases=("case_submitter_id", "nunique"),
            n_correct=("correct", "sum"),
            n_incorrect=("correct", lambda values: int((~values).sum())),
        )
        .sort_values("feature_view", key=lambda col: col.map({view: idx for idx, view in enumerate(FEATURE_VIEW_ORDER)}))
    )

    all_view = driver_scores.loc[driver_scores["feature_view"] == "all_sampled_tissue"]
    all_view_summary = (
        all_view.groupby("clinical_her2_group")["zero_like_score"].agg(["count", "mean", "median"]).reset_index()
    )
    return {
        "n_low_zero_slides": int(all_view["slide_id"].nunique()),
        "signal_channels_by_view": signal_counts.to_dict(orient="records"),
        "stability_counts": stability_counts.to_dict(orient="records"),
        "classifier_counts_by_view": classifier_counts.to_dict(orient="records"),
        "all_sampled_tissue_zero_like_score_by_group": all_view_summary.to_dict(orient="records"),
        "n_slides_expected_profile_4_views": int((view_stability["expected_profile_views"] == 4).sum()),
        "n_slides_expected_profile_3plus_views": int((view_stability["expected_profile_views"] >= 3).sum()),
        "n_slides_opposite_profile_2plus_views": int((view_stability["opposite_profile_views"] >= 2).sum()),
        "n_slides_classifier_wrong_2plus_views": int(
            (classifier_best_predictions.groupby("case_submitter_id")["correct"].apply(lambda values: int((~values).sum())) >= 2).sum()
        ),
    }


def plot_zero_like_score_by_view(plt, sns, driver_scores, asset_dir: Path) -> None:
    plot_df = driver_scores.copy()
    plot_df["feature_view_label"] = plot_df["feature_view"].map(FEATURE_VIEW_LABELS)
    order = [FEATURE_VIEW_LABELS[view] for view in FEATURE_VIEW_ORDER]
    plt.figure(figsize=(11.5, 5.8))
    sns.boxplot(
        data=plot_df,
        x="feature_view_label",
        y="zero_like_score",
        hue="clinical_her2_group",
        order=order,
        fliersize=0,
    )
    sns.stripplot(
        data=plot_df,
        x="feature_view_label",
        y="zero_like_score",
        hue="clinical_her2_group",
        order=order,
        dodge=True,
        alpha=0.42,
        size=3,
        linewidth=0,
        legend=False,
    )
    plt.axhline(0, color="#6b7280", linewidth=1, linestyle="--")
    plt.xlabel("GigaTIME cleanup view")
    plt.ylabel("Case score oriented toward HER2-zero")
    plt.title("Case-Level HER2-Zero-Like Signal by Cleanup View")
    plt.xticks(rotation=18, ha="right")
    plt.legend(title="Clinical group", loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "case_zero_like_score_by_view.png", dpi=180)
    plt.close()


def plot_driver_heatmap(pd, plt, sns, driver_scores, signal_channels, asset_dir: Path) -> None:
    all_view = driver_scores.loc[driver_scores["feature_view"] == "all_sampled_tissue"].copy()
    all_signals = signal_channels.loc[signal_channels["feature_view"] == "all_sampled_tissue"].copy()
    channels = all_signals["channel"].tolist()
    columns = [z_col(channel) for channel in channels if z_col(channel) in all_view.columns]
    if not columns:
        return
    top_low = all_view.loc[all_view["clinical_her2_group"] == "HER2-low"].nsmallest(12, "zero_like_score")
    top_zero = all_view.loc[all_view["clinical_her2_group"] == "HER2-zero"].nlargest(12, "zero_like_score")
    selected = pd.concat([top_low, top_zero], ignore_index=True)
    selected["case_label"] = selected["case_submitter_id"] + " | " + selected["clinical_her2_group"]
    matrix = selected.set_index("case_label")[columns]
    matrix.columns = [column.replace("zero_like_z_", "").replace("_", "-") for column in columns]
    plt.figure(figsize=(9.5, 8.8))
    sns.heatmap(
        matrix,
        cmap="vlag",
        center=0,
        linewidths=0.3,
        linecolor="#f3f4f6",
        cbar_kws={"label": "Zero-like channel z-score"},
    )
    plt.title("Most Label-Consistent HER2-Low and HER2-Zero Driver Cases")
    plt.xlabel("GigaTIME virtual channel")
    plt.ylabel("Case")
    plt.tight_layout()
    plt.savefig(asset_dir / "case_driver_channel_heatmap.png", dpi=180)
    plt.close()


def plot_stability_counts(pd, plt, sns, view_stability, asset_dir: Path) -> None:
    counts = (
        view_stability.groupby(["clinical_her2_group", "expected_profile_views"], as_index=False)
        .size()
        .rename(columns={"size": "n_slides"})
    )
    all_rows = []
    for group in LOW_ZERO_GROUPS:
        for views in range(5):
            match = counts.loc[
                (counts["clinical_her2_group"] == group) & (counts["expected_profile_views"] == views)
            ]
            all_rows.append(
                {
                    "clinical_her2_group": group,
                    "expected_profile_views": views,
                    "n_slides": int(match["n_slides"].iloc[0]) if not match.empty else 0,
                }
            )
    plot_df = pd.DataFrame(all_rows)
    plt.figure(figsize=(8.4, 5.2))
    sns.barplot(
        data=plot_df,
        x="expected_profile_views",
        y="n_slides",
        hue="clinical_her2_group",
    )
    plt.xlabel("Cleanup views matching the expected HER2-low/HER2-zero profile")
    plt.ylabel("Number of slides")
    plt.title("View Stability of Case-Level Driver Scores")
    plt.tight_layout()
    plt.savefig(asset_dir / "case_driver_view_stability.png", dpi=180)
    plt.close()


def plot_classifier_review(plt, sns, classifier_best_predictions, driver_scores, asset_dir: Path) -> None:
    all_preds = classifier_best_predictions.loc[
        classifier_best_predictions["feature_view"] == "all_sampled_tissue"
    ].copy()
    all_scores = driver_scores.loc[
        driver_scores["feature_view"] == "all_sampled_tissue",
        ["case_submitter_id", "zero_like_score"],
    ].copy()
    plot_df = all_preds.merge(all_scores, on="case_submitter_id", how="left")
    plt.figure(figsize=(8.2, 5.8))
    sns.scatterplot(
        data=plot_df,
        x="zero_like_score",
        y="prob_HER2-zero",
        hue="clinical_her2_group",
        style="correct",
        s=82,
        alpha=0.88,
    )
    plt.axvline(0, color="#6b7280", linewidth=1, linestyle="--")
    plt.axhline(0.5, color="#6b7280", linewidth=1, linestyle="--")
    plt.xlabel("All-tissue case score oriented toward HER2-zero")
    plt.ylabel("Classifier probability of HER2-zero")
    plt.title("Classifier Probability Tracks the Case-Level Driver Score")
    plt.tight_layout()
    plt.savefig(asset_dir / "case_driver_classifier_probability.png", dpi=180)
    plt.close()


def format_signal_rows(signal_channels) -> list[list[str]]:
    rows = []
    for view in FEATURE_VIEW_ORDER:
        subset = signal_channels.loc[signal_channels["feature_view"] == view]
        rows.append(
            [
                FEATURE_VIEW_LABELS[view],
                str(subset["channel"].nunique()),
                ", ".join(subset["channel"].tolist()),
            ]
        )
    return rows


def format_group_score_rows(summary: dict) -> list[list[str]]:
    rows = []
    for item in summary["all_sampled_tissue_zero_like_score_by_group"]:
        rows.append(
            [
                item["clinical_her2_group"],
                str(int(item["count"])),
                fmt(item["mean"], 3),
                fmt(item["median"], 3),
            ]
        )
    return rows


def format_stability_rows(summary: dict) -> list[list[str]]:
    counts = {
        (row["clinical_her2_group"], int(row["expected_profile_views"])): int(row["n_slides"])
        for row in summary["stability_counts"]
    }
    return [
        [group] + [str(counts.get((group, views), 0)) for views in range(5)]
        for group in LOW_ZERO_GROUPS
    ]


def format_classifier_rows(summary: dict) -> list[list[str]]:
    rows = []
    for row in summary["classifier_counts_by_view"]:
        n_cases = int(row["n_cases"])
        n_correct = int(row["n_correct"])
        rows.append(
            [
                row["feature_view_label"],
                str(n_cases),
                str(n_correct),
                str(int(row["n_incorrect"])),
                f"{100 * n_correct / n_cases:.1f}%" if n_cases else "",
            ]
        )
    return rows


def format_top_driver_rows(driver_scores, limit_per_group: int = 6) -> list[list[str]]:
    all_view = driver_scores.loc[driver_scores["feature_view"] == "all_sampled_tissue"].copy()
    top_low = all_view.loc[all_view["clinical_her2_group"] == "HER2-low"].nsmallest(
        limit_per_group, "zero_like_score"
    )
    top_zero = all_view.loc[all_view["clinical_her2_group"] == "HER2-zero"].nlargest(
        limit_per_group, "zero_like_score"
    )
    selected = top_low.assign(driver_direction="low-like") 
    selected = selected.copy()
    top_zero = top_zero.copy()
    top_zero["driver_direction"] = "zero-like"
    selected = pd_concat([selected, top_zero])
    rows = []
    for _, row in selected.iterrows():
        rows.append(
            [
                row["case_submitter_id"],
                row["clinical_her2_group"],
                row.get("her2_detail_subgroup", ""),
                row["driver_direction"],
                fmt(row["zero_like_score"], 3),
                fmt(row.get("expected_profile_score", ""), 3),
                row.get("er_status", ""),
                row.get("pr_status", ""),
                fmt(row.get("erbb2_tpm", ""), 2),
            ]
        )
    return rows


def pd_concat(frames):
    import pandas as pd

    return pd.concat(frames, ignore_index=True)


def format_review_rows(classifier_review, limit: int = 12) -> list[list[str]]:
    selected = classifier_review.head(limit)
    rows = []
    for _, row in selected.iterrows():
        rows.append(
            [
                row["case_submitter_id"],
                row["clinical_her2_group"],
                row.get("her2_detail_subgroup", ""),
                str(int(row.get("classifier_incorrect_views", 0))),
                str(int(row.get("opposite_profile_views", 0))) if row.get("opposite_profile_views", "") == row.get("opposite_profile_views", "") else "",
                fmt(row.get("zero_like_score_range", ""), 3),
                fmt(row.get("all_sampled_tissue_zero_like_score", ""), 3),
                row.get("predictions_by_view", ""),
            ]
        )
    return rows


def write_markdown(path: Path, asset_dir: Path, summary: dict, signal_channels, driver_scores, classifier_review) -> None:
    lines = [
        "# Strict High-Trust HER2-Low vs HER2-Zero Case Driver Analysis",
        "",
        "This analysis asks whether the strongest GigaTIME result is broad and case-level believable, or whether it is driven by a few unusual slides.",
        "",
        "The score is built only from GigaTIME virtual channels that significantly separated HER2-low from HER2-zero within each cleanup view. Each channel is standardized and oriented so that higher values are more HER2-zero-like and lower values are more HER2-low-like.",
        "",
        "Important: these are virtual H&E-derived GigaTIME features. This analysis creates a review shortlist; it does not validate the virtual channels as real mIF or prove HER2 isoform biology.",
        "",
        "## Signal Channels Used",
        "",
        markdown_table(["Cleanup view", "N signal channels", "Channels"], format_signal_rows(signal_channels)),
        "",
        "## All-Tissue Case Score",
        "",
        markdown_table(
            ["Clinical group", "N", "Mean zero-like score", "Median zero-like score"],
            format_group_score_rows(summary),
        ),
        "",
        f"{summary['n_slides_expected_profile_3plus_views']} of {summary['n_low_zero_slides']} HER2-low/HER2-zero slides matched their expected direction in at least 3 of 4 cleanup views. {summary['n_slides_opposite_profile_2plus_views']} slides had the opposite profile in at least 2 views and should be prioritized for manual QC/pathology review.",
        "",
        f"![Case zero-like score by view]({asset_link(asset_dir, 'case_zero_like_score_by_view.png')})",
        "",
        "## View Stability",
        "",
        markdown_table(
            ["Clinical group", "0 views", "1 view", "2 views", "3 views", "4 views"],
            format_stability_rows(summary),
        ),
        "",
        f"![Case driver view stability]({asset_link(asset_dir, 'case_driver_view_stability.png')})",
        "",
        "## Most Label-Consistent Driver Cases",
        "",
        "These are not automatically the most important biological cases; they are the clearest examples of the current GigaTIME pattern.",
        "",
        markdown_table(
            [
                "Case",
                "Group",
                "HER2 detail",
                "Driver direction",
                "Zero-like score",
                "Expected-profile score",
                "ER",
                "PR",
                "ERBB2 TPM",
            ],
            format_top_driver_rows(driver_scores),
        ),
        "",
        f"![Case driver channel heatmap]({asset_link(asset_dir, 'case_driver_channel_heatmap.png')})",
        "",
        "## Classifier Error Review",
        "",
        markdown_table(
            ["Cleanup view", "N cases", "Correct", "Incorrect", "Accuracy"],
            format_classifier_rows(summary),
        ),
        "",
        f"{summary['n_slides_classifier_wrong_2plus_views']} cases were misclassified by the best low-vs-zero classifier in at least 2 cleanup views. These are useful cases to inspect because they may represent label noise, slide artifact, tumor-region sampling problems, or real biological exceptions.",
        "",
        f"![Classifier probability and driver score]({asset_link(asset_dir, 'case_driver_classifier_probability.png')})",
        "",
        "## Highest-Priority Manual Review Cases",
        "",
        markdown_table(
            [
                "Case",
                "Group",
                "HER2 detail",
                "Wrong classifier views",
                "Opposite-profile views",
                "Score range",
                "All-tissue score",
                "Predictions by view",
            ],
            format_review_rows(classifier_review),
        ),
        "",
        "## Interpretation",
        "",
        "A useful paper-facing result is beginning to emerge: the HER2-low versus HER2-zero GigaTIME signal is not just a group-average table. We now have a slide-level score, view-stability check, and classifier-error review list that can guide manual pathology review.",
        "",
        "The next highest-value validation step is to open the top label-consistent driver cases and the highest-priority manual-review cases, inspect the sampled H&E regions and virtual mIF-like overlays, and decide whether the signal is tumor-rich, immune/stroma-rich, necrosis/artifact-driven, or plausibly biological.",
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{path.parent.parent / 'results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/case_driver_analysis/case_driver_scores.csv'}`",
        f"- `{path.parent.parent / 'results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/case_driver_analysis/view_stability_by_slide.csv'}`",
        f"- `{path.parent.parent / 'results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/case_driver_analysis/low_zero_classifier_review_cases.csv'}`",
        f"- `{asset_dir}/`",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    asset_dir = Path(args.asset_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)

    np, pd, plt, sns = require_analysis_libs(out_dir / ".matplotlib")
    features, pairwise, predictions, best_metrics = load_inputs(pd, args)

    signal_channels = build_signal_channels(pd, pairwise, args.q_cutoff)
    driver_scores, channel_scores = build_case_driver_scores(np, pd, features, signal_channels)
    view_stability = build_view_stability(pd, driver_scores)
    classifier_best_predictions, classifier_review = build_classifier_review(
        pd, predictions, best_metrics, view_stability
    )
    summary = summarize_counts(pd, signal_channels, driver_scores, view_stability, classifier_best_predictions)

    signal_channels.to_csv(out_dir / "low_zero_signal_channels_by_view.csv", index=False)
    driver_scores.to_csv(out_dir / "case_driver_scores.csv", index=False)
    channel_scores.to_csv(out_dir / "case_driver_channel_z_scores.csv", index=False)
    view_stability.to_csv(out_dir / "view_stability_by_slide.csv", index=False)
    classifier_best_predictions.to_csv(out_dir / "low_zero_best_classifier_predictions_by_view.csv", index=False)
    classifier_review.to_csv(out_dir / "low_zero_classifier_review_cases.csv", index=False)
    (out_dir / "case_driver_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    plot_zero_like_score_by_view(plt, sns, driver_scores, asset_dir)
    plot_driver_heatmap(pd, plt, sns, driver_scores, signal_channels, asset_dir)
    plot_stability_counts(pd, plt, sns, view_stability, asset_dir)
    plot_classifier_review(plt, sns, classifier_best_predictions, driver_scores, asset_dir)

    write_markdown(Path(args.out_markdown), asset_dir, summary, signal_channels, driver_scores, classifier_review)
    print(f"Wrote case-driver outputs to {out_dir}")
    print(f"Wrote case-driver figures to {asset_dir}")
    print(f"Wrote case-driver markdown to {args.out_markdown}")


if __name__ == "__main__":
    main()
