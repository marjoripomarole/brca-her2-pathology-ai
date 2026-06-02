#!/usr/bin/env python3
"""Quantify tissue-composition sensitivity for the high-trust HER2 result."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
LOW_ZERO_GROUPS = ["HER2-low", "HER2-zero"]
KEY_CHANNELS = ["CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CD20", "CK", "Ki67"]
MARKER_BURDEN_CHANNELS = ["CK", "CD68", "PD-L1", "CD11c", "CD3", "CD4", "CD20", "Ki67"]
COMPOSITION_METRICS = [
    "fraction_low_marker_q25",
    "fraction_very_low_marker_q10",
    "fraction_high_marker_q75",
    "mean_marker_burden",
    "mean_CK",
    "fraction_low_ck_q25",
    "fraction_high_ck_q75",
    "mean_tissue_fraction",
    "mean_DAPI",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tile-qc",
        default=str(BASE_RESULT_DIR / "gigatime_cleanup/tile_qc_scores.csv"),
    )
    parser.add_argument(
        "--cleaned-features",
        default=str(BASE_RESULT_DIR / "gigatime_cleanup/cleaned_slide_features.csv"),
    )
    parser.add_argument(
        "--case-driver-scores",
        default=str(BASE_RESULT_DIR / "case_driver_analysis/case_driver_scores.csv"),
    )
    parser.add_argument(
        "--out-dir",
        default=str(BASE_RESULT_DIR / "tissue_composition_sensitivity"),
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_tissue_composition",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_tissue_composition_sensitivity.md",
    )
    parser.add_argument("--min-absolute-ck-high-tiles", type=int, default=4)
    return parser.parse_args()


def require_analysis_libs(mpl_config_dir: Path):
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
        from scipy import stats
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"Missing Python package: {exc.name}. Use `conda run -n gigatime-tcga ...`."
        ) from exc
    sns.set_theme(style="whitegrid", context="notebook")
    return pd, plt, sns, stats


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


def asset_link(asset_dir: Path, image_name: str) -> str:
    return str(asset_dir / image_name).replace("docs/", "")


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    indexed = [(idx, p) for idx, p in enumerate(p_values) if not math.isnan(p)]
    if not indexed:
        return [float("nan")] * len(p_values)
    ranked = sorted(indexed, key=lambda item: item[1])
    adjusted = [float("nan")] * len(p_values)
    prev = 1.0
    m = len(ranked)
    for rank, (idx, p_value) in reversed(list(enumerate(ranked, start=1))):
        q_value = min(prev, p_value * m / rank)
        adjusted[idx] = q_value
        prev = q_value
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


def load_inputs(pd, args: argparse.Namespace):
    tile_qc = pd.read_csv(args.tile_qc, low_memory=False)
    cleaned = pd.read_csv(args.cleaned_features)
    drivers = pd.read_csv(args.case_driver_scores)
    return tile_qc, cleaned, drivers


def add_composition_flags(tile_qc):
    tiles = tile_qc.copy()
    marker_cols = [f"mean_{channel}" for channel in MARKER_BURDEN_CHANNELS]
    tiles["marker_burden"] = tiles[marker_cols].mean(axis=1)
    qc_mask = safe_bool_series(tiles["qc_cellular_tissue"]) if "qc_cellular_tissue" in tiles else tiles.index == tiles.index
    qc_tiles = tiles.loc[qc_mask]
    thresholds = {
        "marker_burden_q10": float(tiles["marker_burden"].quantile(0.10)),
        "marker_burden_q25": float(tiles["marker_burden"].quantile(0.25)),
        "marker_burden_q75": float(tiles["marker_burden"].quantile(0.75)),
        "ck_q25_qc": float(qc_tiles["mean_CK"].quantile(0.25)),
        "ck_q75_qc": float(qc_tiles["mean_CK"].quantile(0.75)),
    }
    tiles["low_marker_q25"] = tiles["marker_burden"] <= thresholds["marker_burden_q25"]
    tiles["very_low_marker_q10"] = tiles["marker_burden"] <= thresholds["marker_burden_q10"]
    tiles["high_marker_q75"] = tiles["marker_burden"] >= thresholds["marker_burden_q75"]
    tiles["low_ck_q25"] = qc_mask & (tiles["mean_CK"] <= thresholds["ck_q25_qc"])
    tiles["high_ck_q75"] = qc_mask & (tiles["mean_CK"] >= thresholds["ck_q75_qc"])
    tiles["high_tissue_low_marker_q25"] = (tiles["tissue_fraction"] >= 0.90) & tiles["low_marker_q25"]
    return tiles, thresholds


def aggregate_slide_composition(pd, tiles):
    rows = []
    total_counts = tiles.groupby("slide_id").size().to_dict()
    for slide_id, group in tiles.groupby("slide_id", sort=False):
        first = group.iloc[0]
        row = {
            "slide_id": slide_id,
            "case_submitter_id": first["case_submitter_id"],
            "clinical_her2_group": first["clinical_her2_group"],
            "n_tiles": int(total_counts[slide_id]),
            "mean_tissue_fraction": float(group["tissue_fraction"].mean()),
            "mean_DAPI": float(group["mean_DAPI"].mean()),
            "mean_CK": float(group["mean_CK"].mean()),
            "mean_marker_burden": float(group["marker_burden"].mean()),
            "median_marker_burden": float(group["marker_burden"].median()),
            "fraction_low_marker_q25": float(group["low_marker_q25"].mean()),
            "fraction_very_low_marker_q10": float(group["very_low_marker_q10"].mean()),
            "fraction_high_marker_q75": float(group["high_marker_q75"].mean()),
            "fraction_low_ck_q25": float(group["low_ck_q25"].mean()),
            "fraction_high_ck_q75": float(group["high_ck_q75"].mean()),
            "fraction_high_tissue_low_marker_q25": float(group["high_tissue_low_marker_q25"].mean()),
        }
        for channel in KEY_CHANNELS:
            column = f"mean_{channel}"
            if column in group.columns:
                row[f"all_tile_mean_{channel}"] = float(group[column].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def build_low_zero_composition_tests(pd, stats, slide_composition):
    low_zero = slide_composition.loc[slide_composition["clinical_her2_group"].isin(LOW_ZERO_GROUPS)].copy()
    rows = []
    for metric in COMPOSITION_METRICS:
        low = low_zero.loc[low_zero["clinical_her2_group"] == "HER2-low", metric].dropna()
        zero = low_zero.loc[low_zero["clinical_her2_group"] == "HER2-zero", metric].dropna()
        if len(low) and len(zero):
            test = stats.mannwhitneyu(low, zero, alternative="two-sided")
            p_value = float(test.pvalue)
            cliff = cliffs_delta(low, zero)
        else:
            p_value = float("nan")
            cliff = float("nan")
        rows.append(
            {
                "metric": metric,
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
    output["mannwhitney_q_value_bh"] = benjamini_hochberg(output["mannwhitney_p_value"].tolist())
    return output.sort_values("mannwhitney_q_value_bh")


def aggregate_absolute_ck_high(pd, tiles, min_tiles: int):
    selected = tiles.loc[tiles["high_ck_q75"]].copy()
    rows = []
    for slide_id, group in selected.groupby("slide_id", sort=False):
        if len(group) < min_tiles:
            continue
        first = group.iloc[0]
        row = {
            "feature_view": "absolute_ck_high_q75",
            "feature_view_label": "Absolute CK-high QC tiles",
            "slide_id": slide_id,
            "case_submitter_id": first["case_submitter_id"],
            "clinical_her2_group": first["clinical_her2_group"],
            "n_tiles_retained": int(len(group)),
            "mean_tissue_fraction": float(group["tissue_fraction"].mean()),
            "mean_marker_burden": float(group["marker_burden"].mean()),
        }
        for channel in KEY_CHANNELS:
            row[f"mean_{channel}"] = float(group[f"mean_{channel}"].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def build_absolute_ck_pairwise(pd, stats, absolute_ck):
    rows = []
    low_zero = absolute_ck.loc[absolute_ck["clinical_her2_group"].isin(LOW_ZERO_GROUPS)].copy()
    for channel in KEY_CHANNELS:
        column = f"mean_{channel}"
        low = low_zero.loc[low_zero["clinical_her2_group"] == "HER2-low", column].dropna()
        zero = low_zero.loc[low_zero["clinical_her2_group"] == "HER2-zero", column].dropna()
        if len(low) and len(zero):
            test = stats.mannwhitneyu(low, zero, alternative="two-sided")
            p_value = float(test.pvalue)
            cliff = cliffs_delta(low, zero)
        else:
            p_value = float("nan")
            cliff = float("nan")
        rows.append(
            {
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
    output["mannwhitney_q_value_bh"] = benjamini_hochberg(output["mannwhitney_p_value"].tolist())
    return output.sort_values("mannwhitney_q_value_bh")


def ols_group_effect(stats, frame, y_col: str, covariates: list[str]):
    rows = frame[["clinical_her2_group", y_col] + covariates].dropna().copy()
    rows = rows.loc[rows["clinical_her2_group"].isin(LOW_ZERO_GROUPS)]
    if len(rows) <= len(covariates) + 2:
        return {"n": int(len(rows)), "beta_low_vs_zero": float("nan"), "p_value": float("nan")}
    y = rows[y_col].astype(float).to_numpy()
    group_low = (rows["clinical_her2_group"] == "HER2-low").astype(float).to_numpy()
    columns = [np.ones(len(rows)), group_low]
    for covariate in covariates:
        values = rows[covariate].astype(float).to_numpy()
        std = values.std(ddof=0)
        if std > 0:
            values = (values - values.mean()) / std
        else:
            values = values - values.mean()
        columns.append(values)
    x = np.column_stack(columns)
    beta, _residuals, rank, _singular = np.linalg.lstsq(x, y, rcond=None)
    if rank < x.shape[1]:
        return {"n": int(len(rows)), "beta_low_vs_zero": float(beta[1]), "p_value": float("nan")}
    residual = y - x @ beta
    df = len(rows) - x.shape[1]
    sigma2 = float((residual @ residual) / df)
    covariance = sigma2 * np.linalg.inv(x.T @ x)
    se = math.sqrt(float(covariance[1, 1]))
    t_stat = float(beta[1] / se) if se > 0 else float("nan")
    p_value = float(2 * stats.t.sf(abs(t_stat), df)) if t_stat == t_stat else float("nan")
    return {"n": int(len(rows)), "beta_low_vs_zero": float(beta[1]), "p_value": p_value}


def build_adjusted_channel_tests(pd, stats, cleaned, slide_composition):
    all_view = cleaned.loc[cleaned["feature_view"] == "all_sampled_tissue"].copy()
    merged = all_view.merge(
        slide_composition[
            [
                "slide_id",
                "case_submitter_id",
                "fraction_low_marker_q25",
                "fraction_very_low_marker_q10",
                "fraction_high_ck_q75",
                "fraction_low_ck_q25",
                "mean_marker_burden",
                "composition_mean_CK",
            ]
        ],
        on=["slide_id", "case_submitter_id"],
        how="inner",
        validate="one_to_one",
    )
    model_specs = {
        "unadjusted": [],
        "adjusted_low_marker_fraction": ["fraction_low_marker_q25"],
        "adjusted_very_low_marker_fraction": ["fraction_very_low_marker_q10"],
        "adjusted_high_ck_fraction": ["fraction_high_ck_q75"],
        "adjusted_mean_ck": ["composition_mean_CK"],
    }
    rows = []
    for channel in KEY_CHANNELS:
        y_col = f"mean_{channel}"
        for model_name, covariates in model_specs.items():
            result = ols_group_effect(stats, merged, y_col, covariates)
            rows.append(
                {
                    "channel": channel,
                    "model": model_name,
                    "covariates": ", ".join(covariates),
                    "n": result["n"],
                    "beta_low_vs_zero": result["beta_low_vs_zero"],
                    "p_value": result["p_value"],
                }
            )
    output = pd.DataFrame(rows)
    output["q_value_bh_within_model"] = float("nan")
    for model_name, group in output.groupby("model"):
        output.loc[group.index, "q_value_bh_within_model"] = benjamini_hochberg(group["p_value"].tolist())
    return output


def build_driver_correlations(pd, stats, slide_composition, drivers):
    all_drivers = drivers.loc[drivers["feature_view"] == "all_sampled_tissue"].copy()
    merged = all_drivers[
        ["slide_id", "case_submitter_id", "clinical_her2_group", "zero_like_score", "expected_profile_score"]
    ].merge(slide_composition, on=["slide_id", "case_submitter_id", "clinical_her2_group"], how="inner")
    merged = merged.loc[merged["clinical_her2_group"].isin(LOW_ZERO_GROUPS)].copy()
    rows = []
    for metric in COMPOSITION_METRICS:
        values = merged[["zero_like_score", metric]].dropna()
        if len(values) >= 3:
            corr = stats.spearmanr(values["zero_like_score"], values[metric])
            rho = float(corr.statistic)
            p_value = float(corr.pvalue)
        else:
            rho = float("nan")
            p_value = float("nan")
        rows.append({"metric": metric, "n": int(len(values)), "spearman_rho": rho, "p_value": p_value})
    output = pd.DataFrame(rows)
    output["q_value_bh"] = benjamini_hochberg(output["p_value"].tolist())
    return output.sort_values("q_value_bh"), merged


def plot_composition_boxplots(plt, sns, slide_composition, asset_dir: Path) -> None:
    plot_df = slide_composition.loc[slide_composition["clinical_her2_group"].isin(LOW_ZERO_GROUPS)].copy()
    metrics = [
        ("fraction_low_marker_q25", "Fraction low-marker tiles"),
        ("fraction_high_marker_q75", "Fraction high-marker tiles"),
        ("mean_marker_burden", "Mean marker burden"),
        ("fraction_high_ck_q75", "Fraction absolute CK-high tiles"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.0))
    axes = axes.ravel()
    for axis, (metric, label) in zip(axes, metrics):
        sns.boxplot(data=plot_df, x="clinical_her2_group", y=metric, order=LOW_ZERO_GROUPS, ax=axis, fliersize=0)
        sns.stripplot(data=plot_df, x="clinical_her2_group", y=metric, order=LOW_ZERO_GROUPS, ax=axis, color="#111827", alpha=0.45, size=3)
        axis.set_xlabel("")
        axis.set_ylabel(label)
    fig.suptitle("HER2-Low vs HER2-Zero Tissue-Composition Metrics")
    fig.tight_layout()
    fig.savefig(asset_dir / "tissue_composition_low_zero_boxplots.png", dpi=180)
    plt.close(fig)


def plot_driver_scatter(plt, sns, driver_composition, asset_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0))
    sns.scatterplot(
        data=driver_composition,
        x="fraction_low_marker_q25",
        y="zero_like_score",
        hue="clinical_her2_group",
        s=72,
        alpha=0.86,
        ax=axes[0],
    )
    axes[0].axhline(0, color="#6b7280", linestyle="--", linewidth=1)
    axes[0].set_xlabel("Fraction low-marker tiles")
    axes[0].set_ylabel("Case zero-like score")
    sns.scatterplot(
        data=driver_composition,
        x="fraction_high_ck_q75",
        y="zero_like_score",
        hue="clinical_her2_group",
        s=72,
        alpha=0.86,
        ax=axes[1],
        legend=False,
    )
    axes[1].axhline(0, color="#6b7280", linestyle="--", linewidth=1)
    axes[1].set_xlabel("Fraction absolute CK-high tiles")
    axes[1].set_ylabel("Case zero-like score")
    fig.suptitle("Case Driver Score Versus Tissue-Composition Metrics")
    fig.tight_layout()
    fig.savefig(asset_dir / "case_driver_vs_tissue_composition.png", dpi=180)
    plt.close(fig)


def plot_adjusted_betas(plt, sns, adjusted_tests, asset_dir: Path) -> None:
    keep_models = ["unadjusted", "adjusted_low_marker_fraction", "adjusted_mean_ck"]
    plot_df = adjusted_tests.loc[adjusted_tests["model"].isin(keep_models)].copy()
    plt.figure(figsize=(11.5, 5.5))
    sns.barplot(data=plot_df, x="channel", y="beta_low_vs_zero", hue="model", order=KEY_CHANNELS)
    plt.axhline(0, color="#374151", linewidth=1)
    plt.xlabel("GigaTIME channel")
    plt.ylabel("Low vs zero beta")
    plt.title("HER2-Low vs HER2-Zero Channel Effects Before and After Composition Adjustment")
    plt.xticks(rotation=30, ha="right")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "composition_adjusted_channel_betas.png", dpi=180)
    plt.close()


def plot_absolute_ck_high(plt, sns, absolute_ck_pairwise, asset_dir: Path) -> None:
    plot_df = absolute_ck_pairwise.copy()
    plt.figure(figsize=(9.5, 4.8))
    sns.barplot(data=plot_df, x="channel", y="delta_low_minus_zero", order=KEY_CHANNELS, color="#0f766e")
    plt.axhline(0, color="#374151", linewidth=1)
    plt.xlabel("GigaTIME channel")
    plt.ylabel("Low minus zero delta")
    plt.title("Low-Zero Channel Delta in Absolute CK-High QC Tiles")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(asset_dir / "absolute_ck_high_low_zero_deltas.png", dpi=180)
    plt.close()


def composition_rows(table) -> list[list[str]]:
    rows = []
    label_map = {
        "fraction_low_marker_q25": "Fraction low-marker tiles",
        "fraction_very_low_marker_q10": "Fraction very-low-marker tiles",
        "fraction_high_marker_q75": "Fraction high-marker tiles",
        "mean_marker_burden": "Mean marker burden",
        "mean_CK": "Mean CK",
        "fraction_low_ck_q25": "Fraction low-CK QC tiles",
        "fraction_high_ck_q75": "Fraction high-CK QC tiles",
        "mean_tissue_fraction": "Mean tissue fraction",
        "mean_DAPI": "Mean DAPI",
    }
    for _, row in table.iterrows():
        rows.append(
            [
                label_map.get(row["metric"], row["metric"]),
                str(int(row["n_low"])),
                str(int(row["n_zero"])),
                fmt(row["mean_low"], 4),
                fmt(row["mean_zero"], 4),
                fmt(row["delta_low_minus_zero"], 4),
                fmt(row["mannwhitney_p_value"], 4),
                fmt(row["mannwhitney_q_value_bh"], 4),
                fmt(row["cliffs_delta"], 3),
            ]
        )
    return rows


def driver_correlation_rows(table) -> list[list[str]]:
    return [
        [row["metric"], str(int(row["n"])), fmt(row["spearman_rho"], 3), fmt(row["p_value"], 4), fmt(row["q_value_bh"], 4)]
        for _, row in table.iterrows()
    ]


def absolute_ck_rows(table) -> list[list[str]]:
    return [
        [
            row["channel"],
            str(int(row["n_low"])),
            str(int(row["n_zero"])),
            fmt(row["mean_low"], 4),
            fmt(row["mean_zero"], 4),
            fmt(row["delta_low_minus_zero"], 4),
            fmt(row["mannwhitney_p_value"], 4),
            fmt(row["mannwhitney_q_value_bh"], 4),
        ]
        for _, row in table.iterrows()
    ]


def adjusted_rows(table, model: str) -> list[list[str]]:
    subset = table.loc[table["model"] == model].copy()
    subset["channel_order"] = subset["channel"].map({channel: idx for idx, channel in enumerate(KEY_CHANNELS)})
    subset = subset.sort_values("channel_order")
    return [
        [
            row["channel"],
            str(int(row["n"])),
            fmt(row["beta_low_vs_zero"], 4),
            fmt(row["p_value"], 4),
            fmt(row["q_value_bh_within_model"], 4),
        ]
        for _, row in subset.iterrows()
    ]


def write_markdown(
    path: Path,
    asset_dir: Path,
    thresholds: dict[str, float],
    composition_tests,
    driver_correlations,
    absolute_ck_pairwise,
    adjusted_tests,
    summary: dict[str, object],
) -> None:
    lines = [
        "# Tissue-Composition Sensitivity for HER2-Low vs HER2-Zero",
        "",
        "This analysis quantifies the caveat raised by the case-driver visual QC: the HER2-low versus HER2-zero signal may partly reflect tissue composition, especially low-marker/stromal-like tiles.",
        "",
        "Definitions:",
        "",
        f"- Marker burden = mean of virtual `{', '.join(MARKER_BURDEN_CHANNELS)}` per tile.",
        f"- Low-marker tile = marker burden <= bottom quartile threshold `{thresholds['marker_burden_q25']:.4f}`.",
        f"- Very-low-marker tile = marker burden <= bottom decile threshold `{thresholds['marker_burden_q10']:.4f}`.",
        f"- Absolute CK-high tile = QC-cellular tile with virtual CK >= QC-tile top quartile threshold `{thresholds['ck_q75_qc']:.4f}`.",
        "",
        "## Main Tissue-Composition Result",
        "",
        markdown_table(
            ["Metric", "N low", "N zero", "Mean low", "Mean zero", "Low-zero delta", "p", "BH q", "Cliff"],
            composition_rows(composition_tests),
        ),
        "",
        f"![Tissue composition boxplots]({asset_link(asset_dir, 'tissue_composition_low_zero_boxplots.png')})",
        "",
        "Interpretation: if HER2-low has more low-marker tiles and fewer high-marker or CK-high tiles than HER2-zero, then the GigaTIME HER2-low versus HER2-zero signal is at least partly a tissue-composition signal.",
        "",
        "## Case-Driver Score Correlation",
        "",
        markdown_table(["Metric", "N", "Spearman rho", "p", "BH q"], driver_correlation_rows(driver_correlations)),
        "",
        f"![Case driver score versus tissue composition]({asset_link(asset_dir, 'case_driver_vs_tissue_composition.png')})",
        "",
        "Interpretation: these correlations test whether the slide-level HER2-zero-like case-driver score is tracking tissue composition. A strong positive correlation with high-marker or CK-high fraction, and a negative correlation with low-marker fraction, supports the tissue-composition caveat.",
        "",
        "## Absolute CK-High Tile Restriction",
        "",
        f"This view keeps only QC-cellular tiles with CK above the global QC-tile 75th percentile and requires at least {summary['min_absolute_ck_high_tiles']} retained tiles per slide.",
        "",
        markdown_table(
            ["Channel", "N low", "N zero", "Mean low", "Mean zero", "Low-zero delta", "p", "BH q"],
            absolute_ck_rows(absolute_ck_pairwise),
        ),
        "",
        f"![Absolute CK-high low-zero deltas]({asset_link(asset_dir, 'absolute_ck_high_low_zero_deltas.png')})",
        "",
        "Interpretation: this is a stricter virtual tumor/epithelial-enriched proxy than the per-slide CK top-25% view. If the low-zero signal weakens here, the current finding should be framed as broader tissue context rather than tumor-cell intrinsic HER2 biology.",
        "",
        "## Composition-Adjusted Channel Effects",
        "",
        "Low-vs-zero beta is from a simple OLS model where HER2-low is compared against HER2-zero. These are exploratory covariate checks, not a final causal model.",
        "",
        "Unadjusted:",
        "",
        markdown_table(["Channel", "N", "Low-vs-zero beta", "p", "BH q"], adjusted_rows(adjusted_tests, "unadjusted")),
        "",
        "Adjusted for low-marker tile fraction:",
        "",
        markdown_table(
            ["Channel", "N", "Low-vs-zero beta", "p", "BH q"],
            adjusted_rows(adjusted_tests, "adjusted_low_marker_fraction"),
        ),
        "",
        "Adjusted for mean CK:",
        "",
        markdown_table(["Channel", "N", "Low-vs-zero beta", "p", "BH q"], adjusted_rows(adjusted_tests, "adjusted_mean_ck")),
        "",
        f"![Composition-adjusted channel betas]({asset_link(asset_dir, 'composition_adjusted_channel_betas.png')})",
        "",
        "## Bottom Line",
        "",
        "- The HER2-low versus HER2-zero GigaTIME signal remains interesting.",
        "- The new tissue-composition analysis makes the caveat stronger, not weaker.",
        "- A presentable claim should emphasize tissue-context association.",
        "- The next scientific step is tumor-rich/pathologist-approved tile restriction before claiming HER2 biology.",
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{path.parent.parent / 'results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/tissue_composition_sensitivity/slide_tissue_composition_metrics.csv'}`",
        f"- `{path.parent.parent / 'results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/tissue_composition_sensitivity/low_zero_composition_tests.csv'}`",
        f"- `{path.parent.parent / 'results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/tissue_composition_sensitivity/composition_adjusted_channel_tests.csv'}`",
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
    pd, plt, sns, stats = require_analysis_libs(out_dir / ".matplotlib")

    tile_qc, cleaned, drivers = load_inputs(pd, args)
    tiles, thresholds = add_composition_flags(tile_qc)
    slide_composition = aggregate_slide_composition(pd, tiles)
    slide_composition["composition_mean_CK"] = slide_composition["mean_CK"]
    composition_tests = build_low_zero_composition_tests(pd, stats, slide_composition)
    absolute_ck = aggregate_absolute_ck_high(pd, tiles, args.min_absolute_ck_high_tiles)
    absolute_ck_pairwise = build_absolute_ck_pairwise(pd, stats, absolute_ck)
    adjusted_tests = build_adjusted_channel_tests(pd, stats, cleaned, slide_composition)
    driver_correlations, driver_composition = build_driver_correlations(pd, stats, slide_composition, drivers)

    summary = {
        "thresholds": thresholds,
        "n_tiles": int(len(tiles)),
        "n_slides": int(slide_composition["slide_id"].nunique()),
        "n_absolute_ck_high_slides": int(absolute_ck["slide_id"].nunique()),
        "n_absolute_ck_high_by_group": absolute_ck.groupby("clinical_her2_group")["slide_id"].nunique().to_dict(),
        "min_absolute_ck_high_tiles": int(args.min_absolute_ck_high_tiles),
    }

    tiles.drop(columns=[], errors="ignore").to_csv(out_dir / "tile_tissue_composition_flags.csv", index=False)
    slide_composition.to_csv(out_dir / "slide_tissue_composition_metrics.csv", index=False)
    composition_tests.to_csv(out_dir / "low_zero_composition_tests.csv", index=False)
    absolute_ck.to_csv(out_dir / "absolute_ck_high_slide_features.csv", index=False)
    absolute_ck_pairwise.to_csv(out_dir / "absolute_ck_high_low_zero_pairwise_tests.csv", index=False)
    adjusted_tests.to_csv(out_dir / "composition_adjusted_channel_tests.csv", index=False)
    driver_correlations.to_csv(out_dir / "driver_tissue_composition_correlations.csv", index=False)
    driver_composition.to_csv(out_dir / "driver_tissue_composition_joined.csv", index=False)
    (out_dir / "tissue_composition_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    plot_composition_boxplots(plt, sns, slide_composition, asset_dir)
    plot_driver_scatter(plt, sns, driver_composition, asset_dir)
    plot_adjusted_betas(plt, sns, adjusted_tests, asset_dir)
    plot_absolute_ck_high(plt, sns, absolute_ck_pairwise, asset_dir)
    write_markdown(
        Path(args.out_markdown),
        asset_dir,
        thresholds,
        composition_tests,
        driver_correlations,
        absolute_ck_pairwise,
        adjusted_tests,
        summary,
    )

    print(f"Wrote tissue-composition sensitivity outputs to {out_dir}")
    print(f"Wrote tissue-composition figures to {asset_dir}")
    print(f"Wrote tissue-composition markdown to {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
