#!/usr/bin/env python3
"""Test HER2 signals under stricter virtual tumor-rich proxy tile filters."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np

from train_her2_classifier_baseline import TASKS, leave_one_out_predictions, metrics_for_predictions
from train_her2_cleaned_classifier_comparison import (
    FEATURE_LABELS,
    GIGATIME_CHANNELS,
    INTERPRETABLE_CHANNELS,
    feature_sets_for,
)


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
GROUP_ORDER = ["HER2-positive", "HER2-low", "HER2-zero"]
LOW_ZERO_GROUPS = ["HER2-low", "HER2-zero"]
KEY_CHANNELS = ["CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CD20", "CK", "Ki67"]
MARKER_BURDEN_CHANNELS = ["CK", "CD68", "PD-L1", "CD11c", "CD3", "CD4", "CD20", "Ki67"]
PROXY_ORDER = [
    "qc_cellular_tissue",
    "ck_top25_within_slide",
    "ck_top16_within_slide",
    "ck_top8_within_slide",
    "ck_top16_non_low_marker",
    "absolute_ck_high_q75",
]
PROXY_LABELS = {
    "qc_cellular_tissue": "QC cellular tissue",
    "ck_top25_within_slide": "CK top 25% within slide",
    "ck_top16_within_slide": "Top 16 CK tiles per slide",
    "ck_top8_within_slide": "Top 8 CK tiles per slide",
    "ck_top16_non_low_marker": "Top 16 CK, non-low-marker",
    "absolute_ck_high_q75": "Absolute CK-high QC tiles",
}
PROXY_MIN_TILES = {
    "qc_cellular_tissue": 16,
    "ck_top25_within_slide": 8,
    "ck_top16_within_slide": 8,
    "ck_top8_within_slide": 4,
    "ck_top16_non_low_marker": 8,
    "absolute_ck_high_q75": 4,
}
TASK_ORDER = ["her2_low_vs_zero", "her2_positive_vs_negative", "her2_three_class"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tile-qc",
        default=str(BASE_RESULT_DIR / "gigatime_cleanup/tile_qc_scores.csv"),
        help="Tile-level high-trust GigaTIME table with clinical labels and cleanup flags.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(BASE_RESULT_DIR / "tumor_proxy_sensitivity"),
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_tumor_proxy_sensitivity",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_tumor_proxy_sensitivity.md",
    )
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


def cliffs_delta(values_a, values_b) -> float:
    a = list(values_a)
    b = list(values_b)
    if not a or not b:
        return float("nan")
    greater = 0
    less = 0
    for value_a in a:
        for value_b in b:
            if value_a > value_b:
                greater += 1
            elif value_a < value_b:
                less += 1
    return (greater - less) / (len(a) * len(b))


def safe_bool_series(series):
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def add_proxy_flags(tiles):
    tiles = tiles.copy()
    marker_cols = [f"mean_{channel}" for channel in MARKER_BURDEN_CHANNELS]
    tiles["marker_burden"] = tiles[marker_cols].mean(axis=1)
    qc_mask = safe_bool_series(tiles["qc_cellular_tissue"])
    qc_tiles = tiles.loc[qc_mask].copy()
    thresholds = {
        "marker_burden_q25_all_tiles": float(tiles["marker_burden"].quantile(0.25)),
        "marker_burden_q50_all_tiles": float(tiles["marker_burden"].quantile(0.50)),
        "ck_q75_qc_tiles": float(qc_tiles["mean_CK"].quantile(0.75)),
    }

    tiles["qc_cellular_tissue"] = qc_mask
    tiles["ck_top25_within_slide"] = safe_bool_series(tiles["ck_enriched_top25"])
    tiles["absolute_ck_high_q75"] = qc_mask & (tiles["mean_CK"] >= thresholds["ck_q75_qc_tiles"])
    tiles["non_low_marker"] = tiles["marker_burden"] > thresholds["marker_burden_q25_all_tiles"]

    tiles["qc_ck_rank_within_slide"] = float("nan")
    qc_rank = qc_tiles.groupby("slide_id")["mean_CK"].rank(method="first", ascending=False)
    tiles.loc[qc_rank.index, "qc_ck_rank_within_slide"] = qc_rank
    tiles["ck_top16_within_slide"] = qc_mask & (tiles["qc_ck_rank_within_slide"] <= 16)
    tiles["ck_top8_within_slide"] = qc_mask & (tiles["qc_ck_rank_within_slide"] <= 8)

    eligible = tiles.loc[qc_mask & tiles["non_low_marker"]].copy()
    tiles["non_low_marker_ck_rank_within_slide"] = float("nan")
    eligible_rank = eligible.groupby("slide_id")["mean_CK"].rank(method="first", ascending=False)
    tiles.loc[eligible_rank.index, "non_low_marker_ck_rank_within_slide"] = eligible_rank
    tiles["ck_top16_non_low_marker"] = (
        qc_mask & tiles["non_low_marker"] & (tiles["non_low_marker_ck_rank_within_slide"] <= 16)
    )
    return tiles, thresholds


def top_fraction_mean(values, fraction: float = 0.10) -> float:
    clean = sorted([float(value) for value in values if not math.isnan(float(value))], reverse=True)
    if not clean:
        return float("nan")
    n_top = max(1, math.ceil(len(clean) * fraction))
    return float(sum(clean[:n_top]) / n_top)


def aggregate_view(pd, tiles, view: str):
    rows = []
    selected = tiles.loc[tiles[view]].copy()
    total_counts = tiles.groupby("slide_id").size().to_dict()
    mean_cols = [f"mean_{channel}" for channel in GIGATIME_CHANNELS if f"mean_{channel}" in tiles.columns]
    frac_cols = [f"frac_{channel}" for channel in GIGATIME_CHANNELS if f"frac_{channel}" in tiles.columns]
    for slide_id, group in selected.groupby("slide_id", sort=False):
        min_tiles = PROXY_MIN_TILES[view]
        if len(group) < min_tiles:
            continue
        first = group.iloc[0]
        row = {
            "feature_view": view,
            "feature_view_label": PROXY_LABELS[view],
            "slide_id": slide_id,
            "case_submitter_id": first["case_submitter_id"],
            "clinical_her2_group": first["clinical_her2_group"],
            "n_tiles_total": int(total_counts.get(slide_id, len(group))),
            "n_tiles_retained": int(len(group)),
            "min_tiles_required": int(min_tiles),
            "retained_fraction": float(len(group) / total_counts.get(slide_id, len(group))),
            "mean_tissue_fraction": float(group["tissue_fraction"].mean()),
            "mean_marker_burden": float(group["marker_burden"].mean()),
            "fraction_non_low_marker": float(group["non_low_marker"].mean()),
        }
        for optional_col in [
            "clinical_her2_group_rule",
            "clinical_her2_group_confidence",
            "her2_ihc_score",
            "her2_ish_status",
            "erbb2_tpm",
            "er_status",
            "pr_status",
        ]:
            if optional_col in group.columns:
                row[optional_col] = first[optional_col]
        for col in mean_cols:
            channel = col.replace("mean_", "", 1)
            values = group[col].dropna()
            row[col] = float(values.mean())
            row[f"median_{channel}"] = float(values.median())
            row[f"p90_{channel}"] = float(values.quantile(0.90))
            row[f"std_{channel}"] = float(values.std(ddof=1)) if len(values) > 1 else float("nan")
            row[f"top10_mean_{channel}"] = top_fraction_mean(values)
        for col in frac_cols:
            row[col] = float(group[col].mean())
        row["virtual_myeloid_checkpoint_score"] = float(group[["mean_CD68", "mean_CD11c", "mean_PD-L1"]].mean(axis=1).mean())
        row["virtual_t_cell_score"] = float(group[["mean_CD3", "mean_CD4", "mean_CD8"]].mean(axis=1).mean())
        rows.append(row)
    return pd.DataFrame(rows)


def aggregate_all_views(pd, tiles):
    return pd.concat([aggregate_view(pd, tiles, view) for view in PROXY_ORDER], ignore_index=True)


def build_retention_summary(pd, tiles, slide_features):
    rows = []
    for view in PROXY_ORDER:
        selected = tiles.loc[tiles[view]].copy()
        raw_counts = selected.groupby(["slide_id", "clinical_her2_group"]).size().reset_index(name="raw_retained_tiles")
        eligible = slide_features.loc[slide_features["feature_view"] == view]
        for group in GROUP_ORDER:
            raw_group = raw_counts.loc[raw_counts["clinical_her2_group"] == group]
            eligible_group = eligible.loc[eligible["clinical_her2_group"] == group]
            rows.append(
                {
                    "feature_view": view,
                    "feature_view_label": PROXY_LABELS[view],
                    "clinical_her2_group": group,
                    "n_slides_with_any_tiles": int(raw_group["slide_id"].nunique()),
                    "n_slides_passing_min_tiles": int(eligible_group["slide_id"].nunique()),
                    "median_tiles_before_min_filter": float(raw_group["raw_retained_tiles"].median())
                    if len(raw_group)
                    else float("nan"),
                    "median_tiles_after_min_filter": float(eligible_group["n_tiles_retained"].median())
                    if len(eligible_group)
                    else float("nan"),
                    "median_retained_fraction": float(eligible_group["retained_fraction"].median())
                    if len(eligible_group)
                    else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def build_low_zero_pairwise(pd, stats, slide_features):
    rows = []
    for view in PROXY_ORDER:
        view_rows = slide_features.loc[slide_features["feature_view"] == view].copy()
        for channel in KEY_CHANNELS:
            col = f"mean_{channel}"
            if col not in view_rows.columns:
                continue
            low = view_rows.loc[view_rows["clinical_her2_group"] == "HER2-low", col].dropna()
            zero = view_rows.loc[view_rows["clinical_her2_group"] == "HER2-zero", col].dropna()
            if len(low) and len(zero):
                test = stats.mannwhitneyu(low, zero, alternative="two-sided")
                p_value = float(test.pvalue)
                cliff = cliffs_delta(low, zero)
            else:
                p_value = float("nan")
                cliff = float("nan")
            rows.append(
                {
                    "feature_view": view,
                    "feature_view_label": PROXY_LABELS[view],
                    "channel": channel,
                    "n_low": int(len(low)),
                    "n_zero": int(len(zero)),
                    "mean_low": float(low.mean()) if len(low) else float("nan"),
                    "mean_zero": float(zero.mean()) if len(zero) else float("nan"),
                    "delta_low_minus_zero": float(low.mean() - zero.mean()) if len(low) and len(zero) else float("nan"),
                    "mannwhitney_p_value": p_value,
                    "cliffs_delta": cliff,
                }
            )
    output = pd.DataFrame(rows)
    output["mannwhitney_q_value_bh_within_view"] = float("nan")
    for view, group in output.groupby("feature_view"):
        output.loc[group.index, "mannwhitney_q_value_bh_within_view"] = benjamini_hochberg(
            group["mannwhitney_p_value"].tolist()
        )
    return output


def run_classifiers(pd, optimize, stats, slide_features, l2_penalty: float):
    prediction_tables = []
    metric_tables = []
    confusion_tables = []
    feature_sets_by_view = {}
    for view in PROXY_ORDER:
        view_rows = slide_features.loc[slide_features["feature_view"] == view].copy()
        if view_rows.empty:
            continue
        view_prediction_tables = []
        feature_sets = feature_sets_for(view_rows)
        feature_sets_by_view[view] = feature_sets
        for task in TASKS:
            for feature_name, feature_cols in feature_sets.items():
                predictions = leave_one_out_predictions(pd, optimize, view_rows, task, feature_name, feature_cols, l2_penalty)
                if predictions.empty:
                    continue
                predictions["feature_view"] = view
                predictions["feature_view_label"] = PROXY_LABELS[view]
                view_prediction_tables.append(predictions)
        if not view_prediction_tables:
            continue
        view_predictions = pd.concat(view_prediction_tables, ignore_index=True)
        view_metrics, view_confusion = metrics_for_predictions(
            pd, stats, view_predictions, {task.name: task for task in TASKS}
        )
        view_metrics["feature_view"] = view
        view_metrics["feature_view_label"] = PROXY_LABELS[view]
        view_confusion["feature_view"] = view
        view_confusion["feature_view_label"] = PROXY_LABELS[view]
        prediction_tables.append(view_predictions)
        metric_tables.append(view_metrics)
        confusion_tables.append(view_confusion)
    predictions = pd.concat(prediction_tables, ignore_index=True) if prediction_tables else pd.DataFrame()
    metrics = pd.concat(metric_tables, ignore_index=True) if metric_tables else pd.DataFrame()
    confusion = pd.concat(confusion_tables, ignore_index=True) if confusion_tables else pd.DataFrame()
    return predictions, metrics, confusion, feature_sets_by_view


def best_h_e_rows(metrics):
    if metrics.empty:
        return metrics
    h_e = metrics.loc[
        (metrics["model"] == "regularized_logistic")
        & (metrics["feature_set"] != "erbb2_rna_reference_not_h_e")
    ].copy()
    h_e["task_order"] = h_e["task"].map({task: idx for idx, task in enumerate(TASK_ORDER)})
    h_e["view_order"] = h_e["feature_view"].map({view: idx for idx, view in enumerate(PROXY_ORDER)})
    return (
        h_e.sort_values(["view_order", "task_order", "balanced_accuracy"], ascending=[True, True, False])
        .groupby(["feature_view", "task"], as_index=False)
        .head(1)
        .sort_values(["view_order", "task_order"])
    )


def best_reference_rows(metrics):
    if metrics.empty:
        return metrics
    ref = metrics.loc[
        (metrics["model"] == "regularized_logistic")
        & (metrics["feature_set"] == "erbb2_rna_reference_not_h_e")
    ].copy()
    ref["task_order"] = ref["task"].map({task: idx for idx, task in enumerate(TASK_ORDER)})
    ref["view_order"] = ref["feature_view"].map({view: idx for idx, view in enumerate(PROXY_ORDER)})
    return (
        ref.sort_values(["view_order", "task_order", "balanced_accuracy"], ascending=[True, True, False])
        .groupby(["feature_view", "task"], as_index=False)
        .head(1)
        .sort_values(["view_order", "task_order"])
    )


def plot_retention(plt, sns, slide_features, asset_dir: Path) -> None:
    plt.figure(figsize=(12.2, 5.5))
    sns.boxplot(
        data=slide_features,
        x="feature_view_label",
        y="n_tiles_retained",
        hue="clinical_her2_group",
        order=[PROXY_LABELS[view] for view in PROXY_ORDER],
        hue_order=GROUP_ORDER,
        fliersize=0,
    )
    sns.stripplot(
        data=slide_features,
        x="feature_view_label",
        y="n_tiles_retained",
        hue="clinical_her2_group",
        order=[PROXY_LABELS[view] for view in PROXY_ORDER],
        hue_order=GROUP_ORDER,
        dodge=True,
        color="#111827",
        alpha=0.36,
        size=2.5,
        legend=False,
    )
    plt.xlabel("Virtual tumor-rich proxy view")
    plt.ylabel("Retained tiles per slide")
    plt.title("Tile Retention Under Virtual Tumor-Rich Proxy Filters")
    plt.xticks(rotation=24, ha="right")
    legend = plt.gca().get_legend()
    if legend:
        legend.set_title("Clinical HER2 group")
    plt.tight_layout()
    plt.savefig(asset_dir / "tumor_proxy_tile_retention.png", dpi=180)
    plt.close()


def plot_low_zero_deltas(plt, sns, pairwise, asset_dir: Path) -> None:
    plot_df = pairwise.copy()
    plot_df["significant_q05"] = plot_df["mannwhitney_q_value_bh_within_view"] < 0.05
    grid = sns.catplot(
        data=plot_df,
        x="channel",
        y="delta_low_minus_zero",
        col="feature_view_label",
        col_wrap=3,
        kind="bar",
        order=KEY_CHANNELS,
        color="#0f766e",
        sharey=False,
        height=3.4,
        aspect=1.25,
    )
    for axis in grid.axes.flat:
        axis.axhline(0, color="#374151", linewidth=1)
        axis.tick_params(axis="x", rotation=30)
    grid.set_axis_labels("GigaTIME channel", "HER2-low minus HER2-zero")
    grid.set_titles("{col_name}")
    grid.fig.subplots_adjust(top=0.88)
    grid.fig.suptitle("HER2-Low vs HER2-Zero Channel Deltas Under Tumor-Rich Proxy Filters")
    grid.savefig(asset_dir / "tumor_proxy_low_zero_channel_deltas.png", dpi=180)
    plt.close(grid.fig)


def plot_low_zero_classifier(plt, sns, metrics, asset_dir: Path) -> None:
    plot_df = metrics.loc[
        (metrics["model"] == "regularized_logistic")
        & (metrics["task"] == "her2_low_vs_zero")
        & (metrics["feature_set"] != "erbb2_rna_reference_not_h_e")
    ].copy()
    if plot_df.empty:
        return
    plot_df["feature_set_label"] = plot_df["feature_set"].map(FEATURE_LABELS).fillna(plot_df["feature_set"])
    plt.figure(figsize=(12.5, 5.8))
    sns.barplot(
        data=plot_df,
        x="feature_view_label",
        y="balanced_accuracy",
        hue="feature_set_label",
        order=[PROXY_LABELS[view] for view in PROXY_ORDER],
    )
    plt.axhline(0.5, color="#6b7280", linestyle="--", linewidth=1)
    plt.ylim(0, 1)
    plt.xlabel("Virtual tumor-rich proxy view")
    plt.ylabel("HER2-low vs HER2-zero balanced accuracy")
    plt.title("Classifier Performance Under Tumor-Rich Proxy Filters")
    plt.xticks(rotation=24, ha="right")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "tumor_proxy_low_zero_classifier_feature_sets.png", dpi=180)
    plt.close()


def plot_best_classifier(plt, sns, best_h_e, asset_dir: Path) -> None:
    plot_df = best_h_e.copy()
    if plot_df.empty:
        return
    plt.figure(figsize=(12.0, 5.6))
    sns.barplot(
        data=plot_df,
        x="feature_view_label",
        y="balanced_accuracy",
        hue="task_label",
        order=[PROXY_LABELS[view] for view in PROXY_ORDER],
    )
    plt.axhline(0.5, color="#6b7280", linestyle="--", linewidth=1)
    plt.axhline(1 / 3, color="#9ca3af", linestyle=":", linewidth=1)
    plt.ylim(0, 1)
    plt.xlabel("Virtual tumor-rich proxy view")
    plt.ylabel("Best balanced accuracy")
    plt.title("Best GigaTIME/H&E Classifier by Tumor-Rich Proxy View")
    plt.xticks(rotation=24, ha="right")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "tumor_proxy_best_classifier_by_view.png", dpi=180)
    plt.close()


def asset_link(asset_dir: Path, image_name: str) -> str:
    return str(asset_dir / image_name).replace("docs/", "")


def retention_rows(retention) -> list[list[str]]:
    rows = []
    for view in PROXY_ORDER:
        for group in GROUP_ORDER:
            match = retention.loc[(retention["feature_view"] == view) & (retention["clinical_her2_group"] == group)]
            if match.empty:
                continue
            row = match.iloc[0]
            rows.append(
                [
                    row["feature_view_label"],
                    group,
                    str(int(row["n_slides_passing_min_tiles"])),
                    fmt(row["median_tiles_after_min_filter"], 1),
                    fmt(row["median_retained_fraction"], 3),
                ]
            )
    return rows


def pairwise_focus_rows(pairwise) -> list[list[str]]:
    rows = []
    focus_channels = ["CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CK"]
    for view in PROXY_ORDER:
        view_rows = pairwise.loc[pairwise["feature_view"] == view].copy()
        for channel in focus_channels:
            match = view_rows.loc[view_rows["channel"] == channel]
            if match.empty:
                continue
            row = match.iloc[0]
            rows.append(
                [
                    row["feature_view_label"],
                    channel,
                    str(int(row["n_low"])),
                    str(int(row["n_zero"])),
                    fmt(row["mean_low"], 4),
                    fmt(row["mean_zero"], 4),
                    fmt(row["delta_low_minus_zero"], 4),
                    fmt(row["mannwhitney_q_value_bh_within_view"], 4),
                ]
            )
    return rows


def classifier_rows(best_h_e) -> list[list[str]]:
    rows = []
    for _, row in best_h_e.iterrows():
        rows.append(
            [
                row["feature_view_label"],
                row["task_label"],
                FEATURE_LABELS.get(row["feature_set"], row["feature_set"]),
                str(int(row["n_cases"])),
                fmt(row["accuracy"]),
                fmt(row["balanced_accuracy"]),
                fmt(row["macro_auc_ovr"]),
                fmt(row.get("sensitivity", float("nan"))),
                fmt(row.get("specificity", float("nan"))),
            ]
        )
    return rows


def low_zero_summary(best_h_e, pairwise) -> list[str]:
    lines = []
    for view in PROXY_ORDER:
        best = best_h_e.loc[(best_h_e["feature_view"] == view) & (best_h_e["task"] == "her2_low_vs_zero")]
        tests = pairwise.loc[pairwise["feature_view"] == view]
        significant = tests.loc[tests["mannwhitney_q_value_bh_within_view"] < 0.05]
        if best.empty:
            continue
        row = best.iloc[0]
        significant_channels = ", ".join(significant["channel"].tolist()) if len(significant) else "none"
        lines.append(
            f"- {row['feature_view_label']}: best low-vs-zero balanced accuracy {fmt(row['balanced_accuracy'])}, "
            f"macro AUC {fmt(row['macro_auc_ovr'])}; q<0.05 channels: {significant_channels}."
        )
    return lines


def write_markdown(
    path: Path,
    asset_dir: Path,
    thresholds: dict[str, float],
    retention,
    pairwise,
    best_h_e,
    best_ref,
    summary: dict[str, object],
) -> None:
    ref_rows = []
    for _, row in best_ref.iterrows():
        ref_rows.append(
            [
                row["feature_view_label"],
                row["task_label"],
                str(int(row["n_cases"])),
                fmt(row["balanced_accuracy"]),
                fmt(row["macro_auc_ovr"]),
            ]
        )

    lines = [
        "# Tumor-Rich Proxy Sensitivity for High-Trust HER2 Analysis",
        "",
        "This analysis asks whether the HER2-low versus HER2-zero GigaTIME signal survives stricter virtual tumor-rich proxy tile filters.",
        "",
        "Important: these are not pathologist tumor annotations. They are GigaTIME-derived proxies using virtual DAPI, CK, marker burden, and existing cellular-tissue cleanup. A true tumor-rich analysis still requires pathologist review or a validated tumor segmentation model.",
        "",
        "## Proxy Definitions",
        "",
        f"- Marker burden = mean of virtual `{', '.join(MARKER_BURDEN_CHANNELS)}` per tile.",
        f"- Low-marker threshold = bottom quartile of marker burden across all tiles: `{thresholds['marker_burden_q25_all_tiles']:.4f}`.",
        f"- Absolute CK-high threshold = top quartile of virtual CK among QC-cellular tiles: `{thresholds['ck_q75_qc_tiles']:.4f}`.",
        "- QC cellular tissue = existing tissue_fraction and virtual DAPI cleanup.",
        "- CK top 25% within slide = existing within-slide CK-enriched top-quarter view.",
        "- Top 16/top 8 CK tiles = fixed-count strongest virtual CK tiles within each slide after QC.",
        "- Top 16 CK, non-low-marker = fixed-count virtual CK tiles after removing low-marker tiles.",
        "",
        "## Retention",
        "",
        markdown_table(
            ["Proxy view", "Clinical group", "Slides passing min tiles", "Median retained tiles", "Median retained fraction"],
            retention_rows(retention),
        ),
        "",
        f"![Tumor proxy tile retention]({asset_link(asset_dir, 'tumor_proxy_tile_retention.png')})",
        "",
        "## HER2-Low Versus HER2-Zero Channel Tests",
        "",
        markdown_table(
            ["Proxy view", "Channel", "N low", "N zero", "Mean low", "Mean zero", "Low-zero delta", "BH q"],
            pairwise_focus_rows(pairwise),
        ),
        "",
        f"![Low-zero channel deltas]({asset_link(asset_dir, 'tumor_proxy_low_zero_channel_deltas.png')})",
        "",
        "## Classifier Sensitivity",
        "",
        "Every classifier result below is leave-one-out cross-validated. The H&E/GigaTIME models exclude the ERBB2 RNA reference.",
        "",
        markdown_table(
            [
                "Proxy view",
                "Task",
                "Best H&E/GigaTIME feature set",
                "N",
                "Accuracy",
                "Balanced accuracy",
                "Macro AUC",
                "Sensitivity",
                "Specificity",
            ],
            classifier_rows(best_h_e),
        ),
        "",
        f"![Best classifier by proxy view]({asset_link(asset_dir, 'tumor_proxy_best_classifier_by_view.png')})",
        "",
        f"![Low-zero classifier feature sets]({asset_link(asset_dir, 'tumor_proxy_low_zero_classifier_feature_sets.png')})",
        "",
        "## Low-Zero Summary",
        "",
        *low_zero_summary(best_h_e, pairwise),
        "",
        "## ERBB2 RNA Reference",
        "",
        "ERBB2 RNA is shown only as a non-H&E reference. It is not affected biologically by these tile filters and should not be used as image-derived evidence.",
        "",
        markdown_table(["Proxy view", "Task", "N", "Balanced accuracy", "Macro AUC"], ref_rows),
        "",
        "## Interpretation",
        "",
        "If HER2-low versus HER2-zero separation remains strong in fixed-count CK-high/non-low-marker views, that supports the idea that GigaTIME is capturing something closer to tumor-rich epithelial context. If it weakens, the current result should remain framed as a broader tissue-context association.",
        "",
        "This still does not validate real mIF channels, diagnose HER2 status, or detect HER2 isoforms. It is a stricter failure-mode analysis that tells us where the signal lives.",
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{summary['out_dir']}/tumor_proxy_slide_features.csv`",
        f"- `{summary['out_dir']}/tumor_proxy_low_zero_pairwise_tests.csv`",
        f"- `{summary['out_dir']}/tumor_proxy_classifier_metrics.csv`",
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

    tiles = pd.read_csv(args.tile_qc, low_memory=False)
    tiles, thresholds = add_proxy_flags(tiles)
    slide_features = aggregate_all_views(pd, tiles)
    retention = build_retention_summary(pd, tiles, slide_features)
    pairwise = build_low_zero_pairwise(pd, stats, slide_features)
    predictions, metrics, confusion, feature_sets = run_classifiers(pd, optimize, stats, slide_features, args.l2_penalty)
    best_h_e = best_h_e_rows(metrics)
    best_ref = best_reference_rows(metrics)

    by_view_group = slide_features.groupby(["feature_view", "clinical_her2_group"])["slide_id"].nunique().astype(int)
    summary = {
        "thresholds": thresholds,
        "proxy_min_tiles": PROXY_MIN_TILES,
        "n_tile_rows": int(len(tiles)),
        "n_slide_feature_rows": int(len(slide_features)),
        "n_slides_by_view": slide_features.groupby("feature_view")["slide_id"].nunique().to_dict(),
        "n_slides_by_view_group": {f"{view}|{group}": int(value) for (view, group), value in by_view_group.items()},
        "out_dir": str(out_dir),
    }

    tiles.to_csv(out_dir / "tumor_proxy_tile_flags.csv", index=False)
    slide_features.to_csv(out_dir / "tumor_proxy_slide_features.csv", index=False)
    retention.to_csv(out_dir / "tumor_proxy_retention_summary.csv", index=False)
    pairwise.to_csv(out_dir / "tumor_proxy_low_zero_pairwise_tests.csv", index=False)
    predictions.to_csv(out_dir / "tumor_proxy_classifier_predictions.csv", index=False)
    metrics.to_csv(out_dir / "tumor_proxy_classifier_metrics.csv", index=False)
    confusion.to_csv(out_dir / "tumor_proxy_classifier_confusion_matrices.csv", index=False)
    best_h_e.to_csv(out_dir / "tumor_proxy_classifier_best_h_e_metrics.csv", index=False)
    best_ref.to_csv(out_dir / "tumor_proxy_classifier_erbb2_reference_metrics.csv", index=False)
    (out_dir / "tumor_proxy_feature_sets.json").write_text(json.dumps(feature_sets, indent=2) + "\n", encoding="utf-8")
    (out_dir / "tumor_proxy_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    plot_retention(plt, sns, slide_features, asset_dir)
    plot_low_zero_deltas(plt, sns, pairwise, asset_dir)
    plot_low_zero_classifier(plt, sns, metrics, asset_dir)
    plot_best_classifier(plt, sns, best_h_e, asset_dir)
    write_markdown(
        Path(args.out_markdown),
        asset_dir,
        thresholds,
        retention,
        pairwise,
        best_h_e,
        best_ref,
        summary,
    )

    print(f"Wrote tumor-proxy sensitivity outputs to {out_dir}")
    print(f"Wrote tumor-proxy figures to {asset_dir}")
    print(f"Wrote tumor-proxy markdown to {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
