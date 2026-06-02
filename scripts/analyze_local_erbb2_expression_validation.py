#!/usr/bin/env python3
"""Validate clinical HER2/GigaTIME findings against local ERBB2 gene expression."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

import numpy as np


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
GROUP_ORDER = ["HER2-positive", "HER2-low", "HER2-zero"]
KEY_CHANNELS = ["CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CD20", "CK", "Ki67"]
VIEW_LABELS = {
    "all_sampled_tissue": "All sampled tissue",
    "qc_cellular_tissue": "QC cellular tissue",
    "ck_enriched_top50": "CK-enriched top 50%",
    "ck_enriched_top25": "CK-enriched top 25%",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression-dir", default="data/tcga_brca/expression_files")
    parser.add_argument(
        "--high-trust-slides",
        default="docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv",
    )
    parser.add_argument(
        "--slide-features",
        default=str(BASE_RESULT_DIR / "gigatime_cleanup/cleaned_slide_features.csv"),
    )
    parser.add_argument(
        "--out-dir",
        default=str(BASE_RESULT_DIR / "local_erbb2_expression_validation"),
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_local_erbb2_validation",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_local_erbb2_validation.md",
    )
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
    if value is None:
        return "NA"
    try:
        f_value = float(value)
    except (TypeError, ValueError):
        return str(value)
    if math.isnan(f_value) or math.isinf(f_value):
        return "NA"
    if abs(f_value) < 0.001 and f_value != 0:
        return f"{f_value:.2e}"
    return f"{f_value:.{digits}g}"


def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    clean = [(idx, p) for idx, p in enumerate(p_values) if not math.isnan(p)]
    q_values = [math.nan] * len(p_values)
    if not clean:
        return q_values
    ordered = sorted(clean, key=lambda item: item[1])
    m = len(ordered)
    running = 1.0
    for rank_from_end, (idx, p_value) in enumerate(reversed(ordered), start=1):
        rank = m - rank_from_end + 1
        running = min(running, p_value * m / rank)
        q_values[idx] = min(running, 1.0)
    return q_values


def auc_from_scores(stats, y_true: np.ndarray, scores: np.ndarray) -> float:
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores, dtype=float)
    valid = np.isfinite(scores)
    y_true = y_true[valid]
    scores = scores[valid]
    n_pos = int(y_true.sum())
    n_neg = int((1 - y_true).sum())
    if n_pos == 0 or n_neg == 0:
        return math.nan
    ranks = stats.rankdata(scores)
    pos_rank_sum = float(ranks[y_true == 1].sum())
    u_stat = pos_rank_sum - n_pos * (n_pos + 1) / 2
    return float(u_stat / (n_pos * n_neg))


def balanced_accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    pos = y_true == 1
    neg = y_true == 0
    sensitivity = float((y_pred[pos] == 1).mean()) if pos.any() else math.nan
    specificity = float((y_pred[neg] == 0).mean()) if neg.any() else math.nan
    return float(np.nanmean([sensitivity, specificity]))


def read_erbb2_from_star_file(path: Path) -> dict[str, object]:
    case_id = path.parent.name
    row = {
        "case_submitter_id": case_id,
        "expression_file": str(path),
        "gene_id": "",
        "gene_name": "ERBB2",
        "unstranded": math.nan,
        "stranded_first": math.nan,
        "stranded_second": math.nan,
        "tpm_unstranded": math.nan,
        "fpkm_unstranded": math.nan,
        "fpkm_uq_unstranded": math.nan,
        "gene_model": "",
    }
    with path.open("r", encoding="utf-8") as handle:
        header: list[str] | None = None
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("# gene-model:"):
                row["gene_model"] = line.replace("# gene-model:", "").strip()
                continue
            if not line or line.startswith("#"):
                continue
            values = line.split("\t")
            if header is None:
                header = values
                continue
            record = dict(zip(header, values))
            if record.get("gene_name") != "ERBB2":
                continue
            for key in row:
                if key in record:
                    row[key] = record[key]
            break
    for col in [
        "unstranded",
        "stranded_first",
        "stranded_second",
        "tpm_unstranded",
        "fpkm_unstranded",
        "fpkm_uq_unstranded",
    ]:
        try:
            row[col] = float(row[col])
        except (TypeError, ValueError):
            row[col] = math.nan
    return row


def extract_erbb2_expression(pd, expression_dir: Path):
    files = sorted(expression_dir.glob("*/*.rna_seq.augmented_star_gene_counts.tsv"))
    rows = [read_erbb2_from_star_file(path) for path in files]
    expression = pd.DataFrame(rows)
    if expression.empty:
        return expression
    expression["local_star_erbb2_tpm"] = pd.to_numeric(expression["tpm_unstranded"], errors="coerce")
    expression["local_star_log1p_erbb2_tpm"] = np.log1p(expression["local_star_erbb2_tpm"])
    expression = expression.sort_values("case_submitter_id").drop_duplicates("case_submitter_id", keep="first")
    return expression


def group_summary(pd, joined):
    rows: list[dict[str, object]] = []
    for group in GROUP_ORDER:
        values = joined.loc[joined["clinical_her2_group"] == group, "local_star_erbb2_tpm"].dropna()
        rows.append(
            {
                "clinical_her2_group": group,
                "n_with_local_erbb2": int(values.shape[0]),
                "median_erbb2_tpm": float(values.median()) if len(values) else math.nan,
                "mean_erbb2_tpm": float(values.mean()) if len(values) else math.nan,
                "q25_erbb2_tpm": float(values.quantile(0.25)) if len(values) else math.nan,
                "q75_erbb2_tpm": float(values.quantile(0.75)) if len(values) else math.nan,
            }
        )
    return pd.DataFrame(rows)


def pairwise_tests(pd, stats, joined):
    comparisons = [
        ("HER2-positive", "HER2-low"),
        ("HER2-positive", "HER2-zero"),
        ("HER2-low", "HER2-zero"),
    ]
    rows: list[dict[str, object]] = []
    for group_a, group_b in comparisons:
        values_a = joined.loc[joined["clinical_her2_group"] == group_a, "local_star_erbb2_tpm"].dropna()
        values_b = joined.loc[joined["clinical_her2_group"] == group_b, "local_star_erbb2_tpm"].dropna()
        p_value = (
            float(stats.mannwhitneyu(values_a, values_b, alternative="two-sided").pvalue)
            if len(values_a) and len(values_b)
            else math.nan
        )
        auc = math.nan
        if len(values_a) and len(values_b):
            y_true = np.array([1] * len(values_a) + [0] * len(values_b))
            y_score = np.concatenate([values_a.to_numpy(), values_b.to_numpy()])
            raw_auc = auc_from_scores(stats, y_true, y_score)
            auc = max(raw_auc, 1.0 - raw_auc)
        rows.append(
            {
                "comparison": f"{group_a} vs {group_b}",
                "n_a": int(len(values_a)),
                "n_b": int(len(values_b)),
                "median_a_erbb2_tpm": float(values_a.median()) if len(values_a) else math.nan,
                "median_b_erbb2_tpm": float(values_b.median()) if len(values_b) else math.nan,
                "delta_median_a_minus_b": (
                    float(values_a.median() - values_b.median()) if len(values_a) and len(values_b) else math.nan
                ),
                "mannwhitney_p_value": p_value,
                "separation_auc_abs_direction": auc,
            }
        )
    q_values = benjamini_hochberg([float(row["mannwhitney_p_value"]) for row in rows])
    for row, q_value in zip(rows, q_values):
        row["mannwhitney_bh_q_value"] = q_value
    return pd.DataFrame(rows)


def classifier_reference_metrics(pd, stats, joined):
    tasks = [
        ("HER2-positive vs non-positive", {"HER2-positive"}, {"HER2-low", "HER2-zero"}),
        ("HER2-low vs HER2-zero", {"HER2-low"}, {"HER2-zero"}),
    ]
    rows: list[dict[str, object]] = []
    for task_name, positive_groups, negative_groups in tasks:
        subset = joined.loc[
            joined["clinical_her2_group"].isin(positive_groups | negative_groups)
            & joined["local_star_erbb2_tpm"].notna()
        ].copy()
        if subset.empty:
            continue
        y_true = subset["clinical_her2_group"].isin(positive_groups).astype(int).to_numpy()
        score = subset["local_star_erbb2_tpm"].to_numpy(dtype=float)
        if len(np.unique(y_true)) < 2:
            auc = math.nan
            ba = math.nan
            best_threshold = math.nan
        else:
            auc_raw = auc_from_scores(stats, y_true, score)
            if auc_raw < 0.5:
                score = -score
                auc = 1.0 - auc_raw
            else:
                auc = auc_raw
            thresholds = np.unique(score)
            best_ba = -1.0
            best_threshold = float(thresholds[0])
            for threshold in thresholds:
                y_pred = (score >= threshold).astype(int)
                current = balanced_accuracy(y_true, y_pred)
                if current > best_ba:
                    best_ba = current
                    best_threshold = float(threshold)
            ba = best_ba
        rows.append(
            {
                "task": task_name,
                "n_cases": int(subset.shape[0]),
                "n_positive": int(y_true.sum()),
                "n_negative": int((1 - y_true).sum()),
                "auc": auc,
                "best_threshold_balanced_accuracy": ba,
                "best_threshold_on_score": best_threshold,
            }
        )
    return pd.DataFrame(rows)


def channel_columns(features) -> list[str]:
    return [f"mean_{channel}" for channel in KEY_CHANNELS if f"mean_{channel}" in features.columns]


def spearman_correlations(pd, stats, features, erbb2, subset_name: str, groups: set[str] | None):
    rows: list[dict[str, object]] = []
    merged = features.merge(
        erbb2[["case_submitter_id", "local_star_erbb2_tpm", "local_star_log1p_erbb2_tpm"]],
        on="case_submitter_id",
        how="inner",
    )
    if groups is not None:
        merged = merged.loc[merged["clinical_her2_group"].isin(groups)].copy()
    for view in sorted(merged["feature_view"].dropna().unique()):
        view_rows = merged.loc[merged["feature_view"] == view].copy()
        for col in channel_columns(view_rows):
            valid = view_rows[[col, "local_star_log1p_erbb2_tpm"]].dropna()
            if valid.shape[0] < 5:
                rho = math.nan
                p_value = math.nan
            else:
                rho, p_value = stats.spearmanr(valid["local_star_log1p_erbb2_tpm"], valid[col])
                rho = float(rho)
                p_value = float(p_value)
            rows.append(
                {
                    "subset": subset_name,
                    "feature_view": view,
                    "feature_view_label": VIEW_LABELS.get(view, view),
                    "channel": col.replace("mean_", ""),
                    "n_cases": int(valid.shape[0]),
                    "spearman_rho_log_erbb2": rho,
                    "p_value": p_value,
                }
            )
    q_values = benjamini_hochberg([float(row["p_value"]) for row in rows])
    for row, q_value in zip(rows, q_values):
        row["bh_q_value"] = q_value
    return pd.DataFrame(rows)


def adjusted_low_zero_tests(pd, stats, features, erbb2):
    rows: list[dict[str, object]] = []
    merged = features.merge(
        erbb2[["case_submitter_id", "local_star_erbb2_tpm", "local_star_log1p_erbb2_tpm"]],
        on="case_submitter_id",
        how="inner",
    )
    merged = merged.loc[merged["clinical_her2_group"].isin({"HER2-low", "HER2-zero"})].copy()
    merged["is_her2_zero"] = (merged["clinical_her2_group"] == "HER2-zero").astype(float)
    for view in sorted(merged["feature_view"].dropna().unique()):
        view_rows = merged.loc[merged["feature_view"] == view].copy()
        for col in channel_columns(view_rows):
            valid = view_rows[[col, "is_her2_zero", "local_star_log1p_erbb2_tpm"]].dropna()
            if valid.shape[0] < 8 or valid["is_her2_zero"].nunique() < 2:
                beta_group = math.nan
                p_value = math.nan
                n_cases = int(valid.shape[0])
            else:
                y = valid[col].to_numpy(dtype=float)
                x = np.column_stack(
                    [
                        np.ones(valid.shape[0]),
                        valid["is_her2_zero"].to_numpy(dtype=float),
                        valid["local_star_log1p_erbb2_tpm"].to_numpy(dtype=float),
                    ]
                )
                beta, *_ = np.linalg.lstsq(x, y, rcond=None)
                residuals = y - x @ beta
                dof = valid.shape[0] - x.shape[1]
                if dof <= 0:
                    beta_group = float(beta[1])
                    p_value = math.nan
                else:
                    sigma2 = float((residuals @ residuals) / dof)
                    cov = sigma2 * np.linalg.pinv(x.T @ x)
                    se = math.sqrt(max(float(cov[1, 1]), 0.0))
                    t_stat = float(beta[1] / se) if se > 0 else math.nan
                    p_value = float(2 * stats.t.sf(abs(t_stat), dof)) if not math.isnan(t_stat) else math.nan
                    beta_group = float(beta[1])
                n_cases = int(valid.shape[0])
            rows.append(
                {
                    "feature_view": view,
                    "feature_view_label": VIEW_LABELS.get(view, view),
                    "channel": col.replace("mean_", ""),
                    "n_cases": n_cases,
                    "beta_her2_zero_vs_low_adjusted_log_erbb2": beta_group,
                    "p_value": p_value,
                }
            )
    q_values = benjamini_hochberg([float(row["p_value"]) for row in rows])
    for row, q_value in zip(rows, q_values):
        row["bh_q_value"] = q_value
    return pd.DataFrame(rows)


def plot_erbb2_by_group(plt, sns, joined, asset_dir: Path):
    plot_rows = joined.loc[joined["local_star_erbb2_tpm"].notna()].copy()
    plot_rows["clinical_her2_group"] = plot_rows["clinical_her2_group"].astype(str)
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.boxplot(
        data=plot_rows,
        x="clinical_her2_group",
        y="local_star_log1p_erbb2_tpm",
        order=GROUP_ORDER,
        ax=ax,
        color="#dbeafe",
        fliersize=0,
    )
    sns.stripplot(
        data=plot_rows,
        x="clinical_her2_group",
        y="local_star_log1p_erbb2_tpm",
        order=GROUP_ORDER,
        ax=ax,
        color="#1f2937",
        size=4,
        alpha=0.75,
        jitter=0.18,
    )
    ax.set_xlabel("Clinical HER2 group")
    ax.set_ylabel("log1p ERBB2 TPM from local STAR counts")
    ax.set_title("Local gene-level ERBB2 expression by clinical HER2 group")
    fig.tight_layout()
    path = asset_dir / "local_erbb2_by_clinical_her2_group.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_reference_metrics(plt, sns, reference_metrics, asset_dir: Path):
    if reference_metrics.empty:
        return
    plot_rows = reference_metrics.melt(
        id_vars=["task"],
        value_vars=["auc", "best_threshold_balanced_accuracy"],
        var_name="metric",
        value_name="value",
    )
    label_map = {
        "auc": "AUC",
        "best_threshold_balanced_accuracy": "Best threshold balanced accuracy",
    }
    plot_rows["metric"] = plot_rows["metric"].map(label_map)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    sns.barplot(data=plot_rows, x="task", y="value", hue="metric", ax=ax, palette=["#0f766e", "#1d4ed8"])
    ax.axhline(0.5, color="#6b7280", linestyle="--", linewidth=1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("")
    ax.set_ylabel("Score")
    ax.set_title("ERBB2 gene expression as a simple reference classifier")
    ax.tick_params(axis="x", rotation=15)
    fig.tight_layout()
    fig.savefig(asset_dir / "local_erbb2_reference_classifier_metrics.png", dpi=180)
    plt.close(fig)


def plot_correlation_heatmap(plt, sns, correlations, asset_dir: Path):
    subset = correlations.loc[
        (correlations["subset"] == "HER2-low/HER2-zero")
        & (correlations["feature_view"].isin(["all_sampled_tissue", "qc_cellular_tissue", "ck_enriched_top25"]))
    ].copy()
    if subset.empty:
        return
    subset["view_channel"] = subset["feature_view_label"] + " - " + subset["channel"]
    pivot = subset.pivot_table(
        index="channel",
        columns="feature_view_label",
        values="spearman_rho_log_erbb2",
        aggfunc="first",
    )
    pivot = pivot.reindex([channel for channel in KEY_CHANNELS if channel in pivot.index])
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    sns.heatmap(
        pivot,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-1,
        vmax=1,
        annot=True,
        fmt=".2f",
        linewidths=0.5,
        cbar_kws={"label": "Spearman rho"},
    )
    ax.set_title("Low/zero GigaTIME channel correlation with gene-level ERBB2")
    ax.set_xlabel("")
    ax.set_ylabel("GigaTIME channel")
    fig.tight_layout()
    fig.savefig(asset_dir / "local_erbb2_gigatime_correlation_heatmap.png", dpi=180)
    plt.close(fig)


def plot_adjusted_counts(plt, sns, adjusted_tests, asset_dir: Path):
    if adjusted_tests.empty:
        return
    summary = (
        adjusted_tests.assign(significant=adjusted_tests["bh_q_value"] < 0.05)
        .groupby("feature_view_label", as_index=False)["significant"]
        .sum()
        .rename(columns={"significant": "n_q_lt_0_05"})
    )
    fig, ax = plt.subplots(figsize=(8, 4.6))
    sns.barplot(data=summary, x="feature_view_label", y="n_q_lt_0_05", ax=ax, color="#0f766e")
    ax.set_xlabel("")
    ax.set_ylabel("Channels with BH q < 0.05")
    ax.set_title("Low-vs-zero channel effects after ERBB2 gene-expression adjustment")
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(asset_dir / "local_erbb2_adjusted_low_zero_q_counts.png", dpi=180)
    plt.close(fig)


def build_markdown(
    summary: dict[str, object],
    group_summary_df,
    pairwise_df,
    reference_metrics_df,
    correlations_df,
    adjusted_tests_df,
    args: argparse.Namespace,
) -> str:
    low_zero_pair = pairwise_df.loc[pairwise_df["comparison"] == "HER2-low vs HER2-zero"]
    low_zero_auc = float(low_zero_pair["separation_auc_abs_direction"].iloc[0]) if not low_zero_pair.empty else math.nan
    low_zero_p = float(low_zero_pair["mannwhitney_p_value"].iloc[0]) if not low_zero_pair.empty else math.nan
    low_zero_q = float(low_zero_pair["mannwhitney_bh_q_value"].iloc[0]) if not low_zero_pair.empty else math.nan
    pos_ref = reference_metrics_df.loc[reference_metrics_df["task"] == "HER2-positive vs non-positive"]
    low_zero_ref = reference_metrics_df.loc[reference_metrics_df["task"] == "HER2-low vs HER2-zero"]
    pos_auc = float(pos_ref["auc"].iloc[0]) if not pos_ref.empty else math.nan
    low_zero_ref_auc = float(low_zero_ref["auc"].iloc[0]) if not low_zero_ref.empty else math.nan

    group_rows = [
        [
            row["clinical_her2_group"],
            int(row["n_with_local_erbb2"]),
            fmt(row["median_erbb2_tpm"]),
            fmt(row["q25_erbb2_tpm"]),
            fmt(row["q75_erbb2_tpm"]),
        ]
        for _, row in group_summary_df.iterrows()
    ]
    pair_rows = [
        [
            row["comparison"],
            int(row["n_a"]),
            int(row["n_b"]),
            fmt(row["median_a_erbb2_tpm"]),
            fmt(row["median_b_erbb2_tpm"]),
            fmt(row["separation_auc_abs_direction"]),
            fmt(row["mannwhitney_p_value"]),
            fmt(row["mannwhitney_bh_q_value"]),
        ]
        for _, row in pairwise_df.iterrows()
    ]
    reference_rows = [
        [
            row["task"],
            int(row["n_cases"]),
            fmt(row["auc"]),
            fmt(row["best_threshold_balanced_accuracy"]),
        ]
        for _, row in reference_metrics_df.iterrows()
    ]

    low_zero_cor = correlations_df.loc[correlations_df["subset"] == "HER2-low/HER2-zero"].copy()
    low_zero_cor["abs_rho"] = low_zero_cor["spearman_rho_log_erbb2"].abs()
    top_cor_rows = [
        [
            row["feature_view_label"],
            row["channel"],
            int(row["n_cases"]),
            fmt(row["spearman_rho_log_erbb2"]),
            fmt(row["p_value"]),
            fmt(row["bh_q_value"]),
        ]
        for _, row in low_zero_cor.sort_values("abs_rho", ascending=False).head(10).iterrows()
    ]

    adjusted_sorted = adjusted_tests_df.sort_values("bh_q_value", na_position="last").head(12)
    adjusted_rows = [
        [
            row["feature_view_label"],
            row["channel"],
            int(row["n_cases"]),
            fmt(row["beta_her2_zero_vs_low_adjusted_log_erbb2"]),
            fmt(row["p_value"]),
            fmt(row["bh_q_value"]),
        ]
        for _, row in adjusted_sorted.iterrows()
    ]
    n_adjusted_sig = int((adjusted_tests_df["bh_q_value"] < 0.05).sum()) if not adjusted_tests_df.empty else 0

    return f"""# Local ERBB2 Gene-Level Validation

Status: expanded gene-level ERBB2 validation using all local GDC STAR augmented gene-count files currently downloaded in this workspace.

## Bottom Line

The local STAR files add gene-level ERBB2 context for {summary["n_local_erbb2_cases"]} TCGA-BRCA cases, including {summary["n_high_trust_with_local_erbb2"]} strict high-trust GigaTIME/HER2 cases and {summary["n_low_zero_high_trust_with_local_erbb2"]} HER2-low/HER2-zero high-trust cases.

This is a useful sanity check, but it is not HER2 isoform validation. These files contain gene-level ERBB2 TPM from STAR counts, not transcript-level isoform proportions, PSI, junction evidence, or antibody-binding-domain information.

Main interpretation:

- ERBB2 gene expression strongly supports the HER2-positive label as a broad molecular sanity check: the simple ERBB2-only reference classifier has AUC {fmt(pos_auc)} for HER2-positive versus non-positive.
- ERBB2 gene expression is much weaker for HER2-low versus HER2-zero: low/zero ERBB2-only AUC is {fmt(low_zero_ref_auc)}, and the pairwise low-vs-zero Mann-Whitney p/q are {fmt(low_zero_p)}/{fmt(low_zero_q)}.
- Therefore, the current GigaTIME HER2-low versus HER2-zero signal is not simply a strong gene-level ERBB2 expression separation.
- In the low/zero subset, GigaTIME channel correlations with gene-level ERBB2 are limited and should be treated as context, not validation.
- After adjusting low-vs-zero channel tests for log ERBB2 TPM in the small local RNA-overlap subset, {n_adjusted_sig} tested channel/view effects remain BH q < 0.05.

## Local ERBB2 Coverage

| Item | Count |
|---|---:|
| Local STAR ERBB2 cases | {summary["n_local_erbb2_cases"]} |
| Strict high-trust slides/cases | {summary["n_high_trust_cases"]} |
| Strict high-trust cases with local ERBB2 | {summary["n_high_trust_with_local_erbb2"]} |
| HER2-low/HER2-zero high-trust cases with local ERBB2 | {summary["n_low_zero_high_trust_with_local_erbb2"]} |

## ERBB2 Expression By Clinical HER2 Group

{markdown_table(["Clinical HER2 group", "N", "Median TPM", "Q25 TPM", "Q75 TPM"], group_rows)}

![Local ERBB2 by clinical HER2 group](assets/clinical_her2_high_trust_tile128_local_erbb2_validation/local_erbb2_by_clinical_her2_group.png)

## Pairwise ERBB2 Gene-Level Tests

{markdown_table(["Comparison", "N A", "N B", "Median A", "Median B", "AUC", "p", "BH q"], pair_rows)}

## ERBB2-Only Reference Classifier

This is not an image model. It is a sanity-check reference using only gene-level ERBB2 TPM.

{markdown_table(["Task", "N cases", "AUC", "Best-threshold balanced accuracy"], reference_rows)}

![Local ERBB2 reference classifier metrics](assets/clinical_her2_high_trust_tile128_local_erbb2_validation/local_erbb2_reference_classifier_metrics.png)

## GigaTIME Correlation With Gene-Level ERBB2

Top absolute low/zero correlations between GigaTIME mean channels and log ERBB2 TPM:

{markdown_table(["Feature view", "Channel", "N", "Spearman rho", "p", "BH q"], top_cor_rows)}

![GigaTIME ERBB2 correlation heatmap](assets/clinical_her2_high_trust_tile128_local_erbb2_validation/local_erbb2_gigatime_correlation_heatmap.png)

## Low-Vs-Zero Channel Tests Adjusted For ERBB2

These models test whether the low-vs-zero GigaTIME channel difference remains after adding log ERBB2 TPM as a covariate. This is limited by the small RNA-overlap subset and should not be overinterpreted.

{markdown_table(["Feature view", "Channel", "N", "Beta zero-vs-low adjusted", "p", "BH q"], adjusted_rows)}

![ERBB2-adjusted channel q counts](assets/clinical_her2_high_trust_tile128_local_erbb2_validation/local_erbb2_adjusted_low_zero_q_counts.png)

## What This Means For The Paper

This result helps us say something careful:

- The clinical HER2-positive labels look molecularly plausible because ERBB2 RNA is high in many HER2-positive cases.
- HER2-low versus HER2-zero is not strongly resolved by gene-level ERBB2 RNA alone in the local overlap subset.
- That makes the GigaTIME low/zero image signal more interesting, but not automatically biological. It could reflect tissue context, source-site/slide-size effects, stromal composition, immune context, or real biology.
- This analysis does not test the Guardia et al. HER2 isoform hypothesis directly. For that, we still need transcript-level isoform labels or RNA-seq reads/junction evidence.

## Output Files

- `{args.out_markdown}`
- `{args.out_dir}/local_erbb2_expression.csv`
- `{args.out_dir}/high_trust_local_erbb2_joined.csv`
- `{args.out_dir}/local_erbb2_group_summary.csv`
- `{args.out_dir}/local_erbb2_pairwise_tests.csv`
- `{args.out_dir}/local_erbb2_reference_classifier_metrics.csv`
- `{args.out_dir}/local_erbb2_gigatime_correlations.csv`
- `{args.out_dir}/local_erbb2_adjusted_low_zero_channel_tests.csv`
- `{args.out_dir}/local_erbb2_validation_summary.json`
"""


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    asset_dir = Path(args.asset_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, stats = require_analysis_libs(out_dir / ".matplotlib")

    erbb2 = extract_erbb2_expression(pd, Path(args.expression_dir))
    high_trust = pd.read_csv(args.high_trust_slides)
    features = pd.read_csv(args.slide_features)
    features = features.loc[features["feature_view"].isin(VIEW_LABELS)].copy()

    high_trust = high_trust.drop_duplicates("case_submitter_id", keep="first").copy()
    joined = high_trust.merge(
        erbb2[
            [
                "case_submitter_id",
                "gene_id",
                "local_star_erbb2_tpm",
                "local_star_log1p_erbb2_tpm",
                "fpkm_unstranded",
                "fpkm_uq_unstranded",
                "expression_file",
                "gene_model",
            ]
        ],
        on="case_submitter_id",
        how="left",
    )

    group_summary_df = group_summary(pd, joined)
    pairwise_df = pairwise_tests(pd, stats, joined)
    reference_metrics_df = classifier_reference_metrics(pd, stats, joined)
    correlations_all = spearman_correlations(pd, stats, features, erbb2, "All high-trust overlap", None)
    correlations_low_zero = spearman_correlations(
        pd, stats, features, erbb2, "HER2-low/HER2-zero", {"HER2-low", "HER2-zero"}
    )
    correlations_df = pd.concat([correlations_all, correlations_low_zero], ignore_index=True)
    adjusted_tests_df = adjusted_low_zero_tests(pd, stats, features, erbb2)

    summary = {
        "n_local_erbb2_cases": int(erbb2["case_submitter_id"].nunique()) if not erbb2.empty else 0,
        "n_high_trust_cases": int(high_trust["case_submitter_id"].nunique()),
        "n_high_trust_with_local_erbb2": int(joined["local_star_erbb2_tpm"].notna().sum()),
        "n_low_zero_high_trust_with_local_erbb2": int(
            joined.loc[
                joined["clinical_her2_group"].isin({"HER2-low", "HER2-zero"}),
                "local_star_erbb2_tpm",
            ]
            .notna()
            .sum()
        ),
        "n_positive_high_trust_with_local_erbb2": int(
            joined.loc[joined["clinical_her2_group"] == "HER2-positive", "local_star_erbb2_tpm"].notna().sum()
        ),
        "n_erbb2_adjusted_channel_tests_q_lt_0_05": int((adjusted_tests_df["bh_q_value"] < 0.05).sum())
        if not adjusted_tests_df.empty
        else 0,
    }

    erbb2.to_csv(out_dir / "local_erbb2_expression.csv", index=False)
    joined.to_csv(out_dir / "high_trust_local_erbb2_joined.csv", index=False)
    group_summary_df.to_csv(out_dir / "local_erbb2_group_summary.csv", index=False)
    pairwise_df.to_csv(out_dir / "local_erbb2_pairwise_tests.csv", index=False)
    reference_metrics_df.to_csv(out_dir / "local_erbb2_reference_classifier_metrics.csv", index=False)
    correlations_df.to_csv(out_dir / "local_erbb2_gigatime_correlations.csv", index=False)
    adjusted_tests_df.to_csv(out_dir / "local_erbb2_adjusted_low_zero_channel_tests.csv", index=False)
    with (out_dir / "local_erbb2_validation_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    plot_erbb2_by_group(plt, sns, joined, asset_dir)
    plot_reference_metrics(plt, sns, reference_metrics_df, asset_dir)
    plot_correlation_heatmap(plt, sns, correlations_df, asset_dir)
    plot_adjusted_counts(plt, sns, adjusted_tests_df, asset_dir)

    markdown = build_markdown(
        summary,
        group_summary_df,
        pairwise_df,
        reference_metrics_df,
        correlations_df,
        adjusted_tests_df,
        args,
    )
    Path(args.out_markdown).write_text(markdown, encoding="utf-8")

    print(f"Wrote {args.out_markdown}")
    print(f"Wrote {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
