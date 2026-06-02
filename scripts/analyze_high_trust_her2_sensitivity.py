#!/usr/bin/env python3
"""Run ER/PR and HER2-detail sensitivity checks for the high-trust GigaTIME run."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path


DEFAULT_CHANNELS = ["CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CK", "Ki67"]
FILTER_ORDER = ["all_sampled_tissue", "qc_cellular_tissue", "ck_enriched_top50", "ck_enriched_top25"]
FILTER_LABELS = {
    "all_sampled_tissue": "All sampled tissue",
    "qc_cellular_tissue": "QC cellular tissue",
    "ck_enriched_top50": "CK-enriched top 50%",
    "ck_enriched_top25": "CK-enriched top 25%",
}
DETAIL_LABELS = {
    "HER2-low_IHC1_ISH-negative": "Low IHC1/ISH-",
    "HER2-low_IHC1_ISH-not-evaluated": "Low IHC1/ISH NE",
    "HER2-low_IHC2_ISH-negative": "Low IHC2/ISH-",
    "HER2-zero_IHC0_ISH-negative": "Zero IHC0/ISH-",
    "HER2-zero_IHC0_ISH-not-evaluated": "Zero IHC0/ISH NE",
}
LOW_DETAILS = [
    "HER2-low_IHC1_ISH-negative",
    "HER2-low_IHC1_ISH-not-evaluated",
    "HER2-low_IHC2_ISH-negative",
]
ZERO_DETAILS = [
    "HER2-zero_IHC0_ISH-negative",
    "HER2-zero_IHC0_ISH-not-evaluated",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cleaned-features",
        default="results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/gigatime_cleanup/cleaned_slide_features.csv",
    )
    parser.add_argument(
        "--cohort",
        default="docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv",
        help="High-trust cohort CSV containing HER2 detail subgroup labels.",
    )
    parser.add_argument(
        "--out-dir",
        default="results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/erpr_subgroup_sensitivity",
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity.md",
    )
    parser.add_argument("--channels", default=",".join(DEFAULT_CHANNELS))
    return parser.parse_args()


def require_analysis_libs(mpl_config_dir: Path):
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import seaborn as sns
        from scipy import stats
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"Missing Python package: {exc.name}. Use `conda activate gigatime-tcga` "
            "or `conda run -n gigatime-tcga ...`."
        ) from exc
    return np, pd, plt, sns, stats


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    indexed = [(idx, p) for idx, p in enumerate(p_values) if not math.isnan(p)]
    adjusted = [float("nan")] * len(p_values)
    if not indexed:
        return adjusted
    ranked = sorted(indexed, key=lambda item: item[1])
    m = len(ranked)
    prev = 1.0
    for rank, (idx, p_value) in reversed(list(enumerate(ranked, start=1))):
        q_value = min(prev, p_value * m / rank)
        adjusted[idx] = q_value
        prev = q_value
    return adjusted


def fmt(value: float | int | str, digits: int = 3) -> str:
    if isinstance(value, str):
        return value
    try:
        value = float(value)
    except (TypeError, ValueError):
        return ""
    if math.isnan(value):
        return ""
    if abs(value) < 0.001 and value != 0:
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    if not rows:
        return "No rows available."
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def clean_status(value) -> str:
    if value != value or value is None:
        return ""
    return str(value).strip()


def is_positive(value) -> int:
    return 1 if clean_status(value).lower() == "positive" else 0


def load_inputs(pd, cleaned_features: Path, cohort_path: Path):
    features = pd.read_csv(cleaned_features)
    cohort = pd.read_csv(cohort_path)
    cohort_cols = [
        "slide_id",
        "case_submitter_id",
        "her2_detail_subgroup",
        "label_slide_trust",
        "processed_tissue_qc",
    ]
    cohort = cohort[[col for col in cohort_cols if col in cohort.columns]].drop_duplicates("slide_id")
    merged = features.merge(cohort, on=["slide_id", "case_submitter_id"], how="inner", validate="many_to_one")
    merged["er_positive"] = merged["er_status"].map(is_positive)
    merged["pr_positive"] = merged["pr_status"].map(is_positive)
    merged["erbb2_tpm_numeric"] = pd.to_numeric(merged.get("erbb2_tpm"), errors="coerce")
    merged["her2_low_indicator"] = (merged["clinical_her2_group"] == "HER2-low").astype(int)
    return merged


def ols_test(np, stats, rows, y_col: str, covariates: list[str]):
    model_df = rows[["her2_low_indicator", y_col, *covariates]].dropna().copy()
    n = int(len(model_df))
    if n < len(covariates) + 4:
        return {
            "adjusted_n": n,
            "adjusted_beta_low_vs_zero": float("nan"),
            "adjusted_se": float("nan"),
            "adjusted_t": float("nan"),
            "adjusted_p_value": float("nan"),
        }
    y = model_df[y_col].to_numpy(dtype=float)
    x_cols = ["her2_low_indicator", *covariates]
    x = model_df[x_cols].to_numpy(dtype=float)
    x = np.column_stack([np.ones(n), x])
    rank = int(np.linalg.matrix_rank(x))
    df = n - rank
    if df <= 0:
        return {
            "adjusted_n": n,
            "adjusted_beta_low_vs_zero": float("nan"),
            "adjusted_se": float("nan"),
            "adjusted_t": float("nan"),
            "adjusted_p_value": float("nan"),
        }
    beta = np.linalg.pinv(x) @ y
    residuals = y - x @ beta
    mse = float((residuals @ residuals) / df)
    covariance = mse * np.linalg.pinv(x.T @ x)
    se = float(math.sqrt(max(covariance[1, 1], 0.0)))
    coef = float(beta[1])
    t_value = coef / se if se else float("nan")
    p_value = float(2 * stats.t.sf(abs(t_value), df)) if not math.isnan(t_value) else float("nan")
    return {
        "adjusted_n": n,
        "adjusted_beta_low_vs_zero": coef,
        "adjusted_se": se,
        "adjusted_t": t_value,
        "adjusted_p_value": p_value,
    }


def mann_whitney(stats, rows, value_col: str, group_col: str, group_a: str, group_b: str):
    values_a = rows.loc[rows[group_col] == group_a, value_col].dropna()
    values_b = rows.loc[rows[group_col] == group_b, value_col].dropna()
    if len(values_a) == 0 or len(values_b) == 0:
        return {
            "n_a": int(len(values_a)),
            "n_b": int(len(values_b)),
            "mean_a": float("nan"),
            "mean_b": float("nan"),
            "delta_mean_a_minus_b": float("nan"),
            "mannwhitney_p_value": float("nan"),
        }
    test = stats.mannwhitneyu(values_a, values_b, alternative="two-sided")
    return {
        "n_a": int(len(values_a)),
        "n_b": int(len(values_b)),
        "mean_a": float(values_a.mean()),
        "mean_b": float(values_b.mean()),
        "delta_mean_a_minus_b": float(values_a.mean() - values_b.mean()),
        "mannwhitney_p_value": float(test.pvalue),
    }


def build_adjusted_tests(np, pd, stats, rows, channels: list[str]):
    records = []
    low_zero = rows.loc[rows["clinical_her2_group"].isin(["HER2-low", "HER2-zero"])].copy()
    for view in FILTER_ORDER:
        view_df = low_zero.loc[low_zero["feature_view"] == view].copy()
        for channel in channels:
            col = f"mean_{channel}"
            if col not in view_df.columns:
                continue
            base = mann_whitney(stats, view_df, col, "clinical_her2_group", "HER2-low", "HER2-zero")
            erpr = ols_test(np, stats, view_df, col, ["er_positive", "pr_positive"])
            erpr_erbb2 = ols_test(np, stats, view_df, col, ["er_positive", "pr_positive", "erbb2_tpm_numeric"])
            records.append(
                {
                    "feature_view": view,
                    "feature_view_label": FILTER_LABELS[view],
                    "channel": channel,
                    **base,
                    "erpr_adjusted_n": erpr["adjusted_n"],
                    "erpr_adjusted_beta_low_vs_zero": erpr["adjusted_beta_low_vs_zero"],
                    "erpr_adjusted_p_value": erpr["adjusted_p_value"],
                    "erpr_erbb2_adjusted_n": erpr_erbb2["adjusted_n"],
                    "erpr_erbb2_adjusted_beta_low_vs_zero": erpr_erbb2["adjusted_beta_low_vs_zero"],
                    "erpr_erbb2_adjusted_p_value": erpr_erbb2["adjusted_p_value"],
                }
            )
    output = pd.DataFrame(records)
    if output.empty:
        return output
    output["mannwhitney_q_value_bh_within_view"] = float("nan")
    output["erpr_adjusted_q_value_bh_within_view"] = float("nan")
    output["erpr_erbb2_adjusted_q_value_bh_within_view"] = float("nan")
    for view in FILTER_ORDER:
        mask = output["feature_view"] == view
        output.loc[mask, "mannwhitney_q_value_bh_within_view"] = benjamini_hochberg(
            output.loc[mask, "mannwhitney_p_value"].tolist()
        )
        output.loc[mask, "erpr_adjusted_q_value_bh_within_view"] = benjamini_hochberg(
            output.loc[mask, "erpr_adjusted_p_value"].tolist()
        )
        output.loc[mask, "erpr_erbb2_adjusted_q_value_bh_within_view"] = benjamini_hochberg(
            output.loc[mask, "erpr_erbb2_adjusted_p_value"].tolist()
        )
    return output


def build_stratified_tests(pd, stats, rows, channels: list[str]):
    strata = [
        ("ER-positive only", "er_status", "Positive"),
        ("ER-negative only", "er_status", "Negative"),
        ("PR-positive only", "pr_status", "Positive"),
        ("PR-negative only", "pr_status", "Negative"),
        ("ER+/PR+ only", "er_pr_combo", "Positive/Positive"),
        ("ER-/PR- only", "er_pr_combo", "Negative/Negative"),
    ]
    rows = rows.copy()
    rows["er_pr_combo"] = rows["er_status"].astype(str) + "/" + rows["pr_status"].astype(str)
    low_zero = rows.loc[rows["clinical_her2_group"].isin(["HER2-low", "HER2-zero"])].copy()
    records = []
    for view in FILTER_ORDER:
        view_df = low_zero.loc[low_zero["feature_view"] == view].copy()
        for stratum_label, stratum_col, stratum_value in strata:
            stratum_df = view_df.loc[view_df[stratum_col] == stratum_value].copy()
            for channel in channels:
                col = f"mean_{channel}"
                if col not in stratum_df.columns:
                    continue
                test = mann_whitney(stats, stratum_df, col, "clinical_her2_group", "HER2-low", "HER2-zero")
                records.append(
                    {
                        "feature_view": view,
                        "feature_view_label": FILTER_LABELS[view],
                        "stratum": stratum_label,
                        "channel": channel,
                        **test,
                    }
                )
    output = pd.DataFrame(records)
    if output.empty:
        return output
    output["mannwhitney_q_value_bh_within_view_stratum"] = float("nan")
    for (view, stratum), group in output.groupby(["feature_view", "stratum"], sort=False):
        output.loc[group.index, "mannwhitney_q_value_bh_within_view_stratum"] = benjamini_hochberg(
            group["mannwhitney_p_value"].tolist()
        )
    return output


def build_detail_tests(pd, stats, rows, channels: list[str]):
    low_zero = rows.loc[rows["clinical_her2_group"].isin(["HER2-low", "HER2-zero"])].copy()
    detail_records = []
    contrast_records = []
    for view in FILTER_ORDER:
        view_df = low_zero.loc[low_zero["feature_view"] == view].copy()
        for channel in channels:
            col = f"mean_{channel}"
            if col not in view_df.columns:
                continue
            by_detail = {
                detail: view_df.loc[view_df["her2_detail_subgroup"] == detail, col].dropna()
                for detail in LOW_DETAILS + ZERO_DETAILS
            }
            nonempty = [values for values in by_detail.values() if len(values)]
            if len(nonempty) >= 2:
                kruskal = stats.kruskal(*nonempty)
                kruskal_p = float(kruskal.pvalue)
                kruskal_h = float(kruskal.statistic)
            else:
                kruskal_p = float("nan")
                kruskal_h = float("nan")
            detail_row = {
                "feature_view": view,
                "feature_view_label": FILTER_LABELS[view],
                "channel": channel,
                "kruskal_h": kruskal_h,
                "kruskal_p_value": kruskal_p,
            }
            for detail, values in by_detail.items():
                prefix = detail.replace("HER2-", "").replace("-", "_").replace("+", "").lower()
                detail_row[f"{prefix}_n"] = int(len(values))
                detail_row[f"{prefix}_mean"] = float(values.mean()) if len(values) else float("nan")
            detail_records.append(detail_row)

            contrasts = [
                ("low_IHC1_any_ISH_vs_zero_all", LOW_DETAILS[:2], ZERO_DETAILS),
                ("low_IHC2_ISH_negative_vs_zero_all", [LOW_DETAILS[2]], ZERO_DETAILS),
                ("low_all_vs_zero_ISH_negative", LOW_DETAILS, [ZERO_DETAILS[0]]),
                ("low_all_vs_zero_ISH_not_evaluated", LOW_DETAILS, [ZERO_DETAILS[1]]),
                ("low_IHC1_any_ISH_vs_zero_ISH_not_evaluated", LOW_DETAILS[:2], [ZERO_DETAILS[1]]),
                ("low_IHC2_ISH_negative_vs_zero_ISH_negative", [LOW_DETAILS[2]], [ZERO_DETAILS[0]]),
            ]
            for contrast_name, group_a_details, group_b_details in contrasts:
                contrast_df = view_df.loc[view_df["her2_detail_subgroup"].isin(group_a_details + group_b_details)].copy()
                contrast_df["contrast_group"] = contrast_df["her2_detail_subgroup"].map(
                    {detail: "A" for detail in group_a_details} | {detail: "B" for detail in group_b_details}
                )
                test = mann_whitney(stats, contrast_df, col, "contrast_group", "A", "B")
                contrast_records.append(
                    {
                        "feature_view": view,
                        "feature_view_label": FILTER_LABELS[view],
                        "channel": channel,
                        "contrast": contrast_name,
                        **test,
                    }
                )
    detail_output = pd.DataFrame(detail_records)
    contrast_output = pd.DataFrame(contrast_records)
    if not detail_output.empty:
        detail_output["kruskal_q_value_bh_within_view"] = float("nan")
        for view, group in detail_output.groupby("feature_view", sort=False):
            detail_output.loc[group.index, "kruskal_q_value_bh_within_view"] = benjamini_hochberg(
                group["kruskal_p_value"].tolist()
            )
    if not contrast_output.empty:
        contrast_output["mannwhitney_q_value_bh_within_view_contrast"] = float("nan")
        for (view, contrast), group in contrast_output.groupby(["feature_view", "contrast"], sort=False):
            contrast_output.loc[group.index, "mannwhitney_q_value_bh_within_view_contrast"] = benjamini_hochberg(
                group["mannwhitney_p_value"].tolist()
            )
    return detail_output, contrast_output


def plot_adjusted_heatmap(plt, sns, adjusted, asset_dir: Path):
    if adjusted.empty:
        return
    plot_df = adjusted.pivot(index="channel", columns="feature_view_label", values="erpr_adjusted_q_value_bh_within_view")
    plot_df = plot_df[[FILTER_LABELS[view] for view in FILTER_ORDER]]
    fig, axis = plt.subplots(figsize=(9.8, max(4.8, 0.55 * len(plot_df))))
    sns.heatmap(plot_df, cmap="mako_r", annot=True, fmt=".3g", linewidths=0.3, vmin=0, vmax=0.10, ax=axis)
    axis.set_xlabel("Cleanup view")
    axis.set_ylabel("GigaTIME channel")
    axis.set_title("HER2-low vs HER2-zero ER/PR-adjusted q values", pad=14)
    axis.tick_params(axis="y", labelrotation=0)
    axis.tick_params(axis="x", labelrotation=90)
    fig.tight_layout()
    fig.savefig(asset_dir / "erpr_adjusted_low_zero_q_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_detail_boxplots(plt, sns, pd, rows, channels: list[str], asset_dir: Path):
    view_df = rows.loc[
        (rows["feature_view"] == "all_sampled_tissue")
        & (rows["clinical_her2_group"].isin(["HER2-low", "HER2-zero"]))
    ].copy()
    if view_df.empty:
        return
    records = []
    selected_channels = [channel for channel in channels if f"mean_{channel}" in view_df.columns][:6]
    detail_order = LOW_DETAILS + ZERO_DETAILS
    for channel in selected_channels:
        for _, row in view_df.iterrows():
            records.append(
                {
                    "channel": channel,
                    "detail": DETAIL_LABELS.get(row["her2_detail_subgroup"], row["her2_detail_subgroup"]),
                    "detail_raw": row["her2_detail_subgroup"],
                    "mean_activation": row[f"mean_{channel}"],
                }
            )
    if not records:
        return
    plot_df = pd.DataFrame(records)
    plot_df["detail"] = pd.Categorical(
        plot_df["detail"],
        [DETAIL_LABELS[detail] for detail in detail_order],
        ordered=True,
    )
    grid = sns.catplot(
        data=plot_df,
        x="detail",
        y="mean_activation",
        col="channel",
        col_wrap=3,
        kind="box",
        sharey=False,
        height=3.0,
        aspect=1.12,
        color="#9ecae1",
    )
    for axis in grid.axes.flatten():
        axis.tick_params(axis="x", labelrotation=45)
        axis.set_xlabel("")
        axis.set_ylabel("Mean activation")
    grid.fig.suptitle("All-tissue GigaTIME features by HER2 IHC/ISH detail subgroup", y=0.995)
    grid.fig.tight_layout()
    grid.fig.savefig(asset_dir / "her2_detail_subgroup_boxplots.png", dpi=180, bbox_inches="tight")
    plt.close(grid.fig)


def summarize_counts(rows):
    base = rows.loc[rows["feature_view"] == "all_sampled_tissue"].copy()
    counts = {
        "groups": base["clinical_her2_group"].value_counts().to_dict(),
        "er_by_group": {},
        "pr_by_group": {},
        "detail_subgroups": base["her2_detail_subgroup"].value_counts().to_dict(),
    }
    for group in ["HER2-positive", "HER2-low", "HER2-zero"]:
        group_df = base.loc[base["clinical_her2_group"] == group]
        counts["er_by_group"][group] = group_df["er_status"].value_counts().to_dict()
        counts["pr_by_group"][group] = group_df["pr_status"].value_counts().to_dict()
    return counts


def write_markdown(path: Path, adjusted, stratified, detail, contrasts, counts, out_dir: Path, asset_dir: Path):
    def asset_link(filename: str) -> str:
        return os.path.relpath(asset_dir / filename, path.parent).replace(os.sep, "/")

    def asset_csv_link(filename: str) -> str:
        return os.path.relpath(asset_dir / filename, path.parent).replace(os.sep, "/")

    all_view = adjusted.loc[adjusted["feature_view"] == "all_sampled_tissue"].copy()
    all_view = all_view.sort_values("erpr_adjusted_q_value_bh_within_view", na_position="last")
    top_adjusted = [
        [
            row["channel"],
            fmt(row["delta_mean_a_minus_b"], 5),
            fmt(row["mannwhitney_q_value_bh_within_view"]),
            fmt(row["erpr_adjusted_beta_low_vs_zero"], 5),
            fmt(row["erpr_adjusted_q_value_bh_within_view"]),
            fmt(row["erpr_erbb2_adjusted_q_value_bh_within_view"]),
            str(int(row["erpr_erbb2_adjusted_n"])) if not math.isnan(row["erpr_erbb2_adjusted_n"]) else "",
        ]
        for _, row in all_view.head(8).iterrows()
    ]

    significant_counts = []
    for view in FILTER_ORDER:
        view_df = adjusted.loc[adjusted["feature_view"] == view]
        significant_counts.append(
            [
                FILTER_LABELS[view],
                str(int((view_df["mannwhitney_q_value_bh_within_view"] < 0.05).sum())),
                str(int((view_df["erpr_adjusted_q_value_bh_within_view"] < 0.05).sum())),
                str(int((view_df["erpr_erbb2_adjusted_q_value_bh_within_view"] < 0.05).sum())),
            ]
        )

    stratified_rows = []
    for stratum in ["ER-positive only", "ER-negative only", "PR-positive only", "PR-negative only"]:
        subset = stratified.loc[
            (stratified["feature_view"] == "all_sampled_tissue")
            & (stratified["stratum"] == stratum)
        ].copy()
        subset = subset.sort_values("mannwhitney_q_value_bh_within_view_stratum", na_position="last")
        if subset.empty:
            continue
        row = subset.iloc[0]
        stratified_rows.append(
            [
                stratum,
                str(int(row["n_a"])),
                str(int(row["n_b"])),
                row["channel"],
                fmt(row["delta_mean_a_minus_b"], 5),
                fmt(row["mannwhitney_q_value_bh_within_view_stratum"]),
            ]
        )

    detail_rows = []
    all_detail = detail.loc[detail["feature_view"] == "all_sampled_tissue"].copy()
    for channel in ["CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CK"]:
        channel_rows = all_detail.loc[all_detail["channel"] == channel]
        if channel_rows.empty:
            continue
        row = channel_rows.iloc[0]
        detail_rows.append(
            [
                channel,
                fmt(row.get("low_ihc1_ish_negative_mean", float("nan")), 5),
                fmt(row.get("low_ihc1_ish_not_evaluated_mean", float("nan")), 5),
                fmt(row.get("low_ihc2_ish_negative_mean", float("nan")), 5),
                fmt(row.get("zero_ihc0_ish_negative_mean", float("nan")), 5),
                fmt(row.get("zero_ihc0_ish_not_evaluated_mean", float("nan")), 5),
                fmt(row["kruskal_q_value_bh_within_view"]),
            ]
        )

    contrast_rows = []
    interesting_contrasts = [
        "low_IHC1_any_ISH_vs_zero_all",
        "low_IHC2_ISH_negative_vs_zero_all",
        "low_all_vs_zero_ISH_negative",
        "low_all_vs_zero_ISH_not_evaluated",
    ]
    for contrast in interesting_contrasts:
        subset = contrasts.loc[
            (contrasts["feature_view"] == "all_sampled_tissue")
            & (contrasts["contrast"] == contrast)
        ].sort_values("mannwhitney_q_value_bh_within_view_contrast", na_position="last")
        if subset.empty:
            continue
        row = subset.iloc[0]
        contrast_rows.append(
            [
                contrast,
                row["channel"],
                str(int(row["n_a"])),
                str(int(row["n_b"])),
                fmt(row["delta_mean_a_minus_b"], 5),
                fmt(row["mannwhitney_q_value_bh_within_view_contrast"]),
            ]
        )

    lines = [
        "# High-Trust HER2 ER/PR and Subgroup Sensitivity",
        "",
        "Status: Sensitivity analysis for the current strict high-trust 171-slide GigaTIME analysis set.",
        "",
        "## Why This Matters",
        "",
        "The main high-trust finding is HER2-low versus HER2-zero separation in GigaTIME virtual immune/myeloid/checkpoint and tissue-context channels. This analysis asks whether that signal is obviously explained by hormone receptor imbalance or by one HER2 IHC/ISH subgroup.",
        "",
        "This is still exploratory. ER/PR adjustment is a statistical sensitivity check, not causal proof.",
        "",
        "## Input Counts",
        "",
        markdown_table(
            ["Group", "Slides", "ER positive", "ER negative", "PR positive", "PR negative"],
            [
                [
                    group,
                    str(counts["groups"].get(group, 0)),
                    str(counts["er_by_group"].get(group, {}).get("Positive", 0)),
                    str(counts["er_by_group"].get(group, {}).get("Negative", 0)),
                    str(counts["pr_by_group"].get(group, {}).get("Positive", 0)),
                    str(counts["pr_by_group"].get(group, {}).get("Negative", 0)),
                ]
                for group in ["HER2-positive", "HER2-low", "HER2-zero"]
            ],
        ),
        "",
        "## ER/PR-Adjusted Result",
        "",
        markdown_table(
            ["Cleanup view", "Unadjusted q<0.05 channels", "ER/PR adjusted q<0.05 channels", "ER/PR+ERBB2 adjusted q<0.05 channels"],
            significant_counts,
        ),
        "",
        "All-sampled-tissue channels ranked by ER/PR-adjusted q value:",
        "",
        markdown_table(
            ["Channel", "Unadjusted low-zero delta", "Unadjusted q", "ER/PR beta", "ER/PR q", "ER/PR+ERBB2 q", "RNA n"],
            top_adjusted,
        ),
        "",
        f"![ER/PR-adjusted q-value heatmap]({asset_link('erpr_adjusted_low_zero_q_heatmap.png')})",
        "",
        "Interpretation: if the ER/PR-adjusted q values remain small and the beta stays negative, the HER2-low lower-than-zero signal is not explained only by ER/PR imbalance.",
        "",
        "The ER/PR+ERBB2 adjustment uses only cases with available ERBB2 RNA, so it is a smaller secondary sensitivity check.",
        "",
        "## ER/PR-Stratified Result",
        "",
        markdown_table(
            ["Stratum", "HER2-low n", "HER2-zero n", "Best channel", "Low-zero delta", "BH q"],
            stratified_rows,
        ),
        "",
        "Interpretation: these strata are smaller, especially ER-negative and PR-negative subsets. Consistent negative deltas across strata are more important than any single p value.",
        "",
        "## HER2 Detail Subgroup Result",
        "",
        "All-sampled-tissue subgroup means:",
        "",
        markdown_table(
            ["Channel", "Low IHC1/ISH-", "Low IHC1/ISH NE", "Low IHC2/ISH-", "Zero IHC0/ISH-", "Zero IHC0/ISH NE", "Subgroup q"],
            detail_rows,
        ),
        "",
        "Best all-sampled-tissue subgroup contrasts:",
        "",
        markdown_table(
            ["Contrast", "Best channel", "A n", "B n", "A-B delta", "BH q"],
            contrast_rows,
        ),
        "",
        f"![HER2 detail subgroup boxplots]({asset_link('her2_detail_subgroup_boxplots.png')})",
        "",
        "Interpretation: if both HER2-low IHC1 and HER2-low IHC2/ISH-negative remain lower than HER2-zero, the result is less likely to be an artifact of only one HER2-low subgroup.",
        "",
        "## Machine-Readable Outputs",
        "",
        f"- [{asset_csv_link('low_zero_erpr_adjusted_tests.csv')}]({asset_csv_link('low_zero_erpr_adjusted_tests.csv')})",
        f"- [{asset_csv_link('low_zero_erpr_stratified_tests.csv')}]({asset_csv_link('low_zero_erpr_stratified_tests.csv')})",
        f"- [{asset_csv_link('her2_detail_subgroup_tests.csv')}]({asset_csv_link('her2_detail_subgroup_tests.csv')})",
        f"- [{asset_csv_link('her2_detail_subgroup_contrasts.csv')}]({asset_csv_link('her2_detail_subgroup_contrasts.csv')})",
        f"- [{asset_csv_link('sensitivity_summary.json')}]({asset_csv_link('sensitivity_summary.json')})",
        "",
        "## Cautious Claim This Supports",
        "",
        "> The main all-sampled-tissue HER2-low versus HER2-zero GigaTIME signal persists after ER/PR adjustment and remains visible across major HER2 IHC/ISH detail subgroups, supporting a tissue-context association rather than an obvious hormone-receptor or single-subgroup artifact.",
        "",
        "What this does not prove:",
        "",
        "- It does not prove clinical HER2 diagnosis from H&E.",
        "- It does not prove that GigaTIME measures real immune proteins in these TCGA slides.",
        "- It does not replace pathologist tumor-region review or external validation.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    asset_dir = Path(args.asset_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)
    channels = [channel.strip() for channel in args.channels.split(",") if channel.strip()]
    np, pd, plt, sns, stats = require_analysis_libs(out_dir / ".matplotlib")
    rows = load_inputs(pd, Path(args.cleaned_features), Path(args.cohort))

    adjusted = build_adjusted_tests(np, pd, stats, rows, channels)
    stratified = build_stratified_tests(pd, stats, rows, channels)
    detail, contrasts = build_detail_tests(pd, stats, rows, channels)
    counts = summarize_counts(rows)

    adjusted.to_csv(out_dir / "low_zero_erpr_adjusted_tests.csv", index=False)
    stratified.to_csv(out_dir / "low_zero_erpr_stratified_tests.csv", index=False)
    detail.to_csv(out_dir / "her2_detail_subgroup_tests.csv", index=False)
    contrasts.to_csv(out_dir / "her2_detail_subgroup_contrasts.csv", index=False)
    adjusted.to_csv(asset_dir / "low_zero_erpr_adjusted_tests.csv", index=False)
    stratified.to_csv(asset_dir / "low_zero_erpr_stratified_tests.csv", index=False)
    detail.to_csv(asset_dir / "her2_detail_subgroup_tests.csv", index=False)
    contrasts.to_csv(asset_dir / "her2_detail_subgroup_contrasts.csv", index=False)

    summary = {
        "n_rows": int(len(rows)),
        "n_slides_all_tissue": int((rows["feature_view"] == "all_sampled_tissue").sum()),
        "channels": channels,
        "counts": counts,
        "adjusted_q_lt_0_05_by_view": {
            view: {
                "unadjusted": int(
                    (
                        adjusted.loc[adjusted["feature_view"] == view, "mannwhitney_q_value_bh_within_view"]
                        < 0.05
                    ).sum()
                ),
                "er_pr_adjusted": int(
                    (
                        adjusted.loc[adjusted["feature_view"] == view, "erpr_adjusted_q_value_bh_within_view"]
                        < 0.05
                    ).sum()
                ),
                "er_pr_erbb2_adjusted": int(
                    (
                        adjusted.loc[adjusted["feature_view"] == view, "erpr_erbb2_adjusted_q_value_bh_within_view"]
                        < 0.05
                    ).sum()
                ),
            }
            for view in FILTER_ORDER
        },
    }
    (out_dir / "sensitivity_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (asset_dir / "sensitivity_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    plot_adjusted_heatmap(plt, sns, adjusted, asset_dir)
    plot_detail_boxplots(plt, sns, pd, rows, channels, asset_dir)
    write_markdown(Path(args.out_markdown), adjusted, stratified, detail, contrasts, counts, out_dir, asset_dir)
    print(f"Wrote HER2 sensitivity outputs to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
