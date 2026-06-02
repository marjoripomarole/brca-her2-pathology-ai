#!/usr/bin/env python3
"""Clinical, slide, and site covariate sensitivity for HER2-low versus HER2-zero."""

from __future__ import annotations

import argparse
import json
import math
import os
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
)
from train_her2_classifier_baseline import fit_predict_logistic, standardize_train_test
from train_her2_cleaned_classifier_comparison import GIGATIME_CHANNELS


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
KEY_CHANNELS = ["CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CD20", "CK", "Ki67"]
VIEW_LABELS = {
    "qc_cellular_tissue": "QC cellular tissue",
    "ck_top25_within_slide": "CK top 25% within slide",
    "ck_top16_within_slide": "Top 16 CK tiles per slide",
    "ck_top8_within_slide": "Top 8 CK tiles per slide",
    "ck_top16_non_low_marker": "Top 16 CK, non-low-marker",
    "absolute_ck_high_q75": "Absolute CK-high QC tiles",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--slide-features",
        default=str(BASE_RESULT_DIR / "tumor_proxy_sensitivity/tumor_proxy_slide_features.csv"),
        help="Slide-level tumor-proxy feature table.",
    )
    parser.add_argument(
        "--high-trust-slides",
        default="docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv",
    )
    parser.add_argument(
        "--out-dir",
        default=str(BASE_RESULT_DIR / "clinical_covariate_sensitivity"),
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_clinical_covariates",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_clinical_covariate_sensitivity.md",
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260602)
    parser.add_argument("--l2-penalty", type=float, default=1.0)
    parser.add_argument("--min-site-count", type=int, default=5)
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


def collapse_histology(value: object) -> str:
    text = str(value or "").lower()
    if "ductal" in text:
        return "Ductal"
    if "lobular" in text:
        return "Lobular"
    if text in {"", "nan", "[not available]"}:
        return "Unknown"
    return "Other"


def collapse_stage(value: object) -> str:
    text = str(value or "").upper()
    if "IV" in text:
        return "Stage IV"
    if "III" in text:
        return "Stage III"
    if "II" in text:
        return "Stage II"
    if "I" in text:
        return "Stage I"
    return "Unknown"


def tss_code(case_submitter_id: object) -> str:
    parts = str(case_submitter_id or "").split("-")
    return parts[1] if len(parts) > 1 else "Unknown"


def site_group(values, min_count: int) -> list[str]:
    counts = values.value_counts(dropna=False)
    return [str(value) if counts.get(value, 0) >= min_count else "Other" for value in values]


def add_covariates(pd, features, metadata, min_site_count: int):
    meta_cols = [
        "case_submitter_id",
        "slide_id",
        "histological_type",
        "pathologic_stage",
        "history_neoadjuvant_treatment",
        "patient_gender",
        "slide_file_size_mb",
        "slide_width",
        "slide_height",
    ]
    available = [col for col in meta_cols if col in metadata.columns]
    merged = features.merge(
        metadata[available].drop_duplicates(["case_submitter_id", "slide_id"]),
        on=["case_submitter_id", "slide_id"],
        how="left",
    )
    merged["histology_group"] = merged.get("histological_type", "").map(collapse_histology)
    merged["stage_group"] = merged.get("pathologic_stage", "").map(collapse_stage)
    merged["tss_code"] = merged["case_submitter_id"].map(tss_code)
    low_zero = merged.loc[merged["clinical_her2_group"].isin(LOW_ZERO_GROUPS)]
    tss_map = dict(zip(low_zero.index, site_group(low_zero["tss_code"], min_site_count)))
    merged["tss_group"] = [tss_map.get(idx, "Other") for idx in merged.index]
    for col in ["er_status", "pr_status", "histology_group", "stage_group", "tss_group"]:
        merged[col] = merged[col].fillna("Unknown").astype(str).replace({"nan": "Unknown", "": "Unknown"})
    for col in ["slide_file_size_mb", "slide_width", "slide_height", "erbb2_tpm"]:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")
    if "slide_file_size_mb" in merged.columns:
        merged["log_slide_file_size_mb"] = np.log1p(merged["slide_file_size_mb"])
    if "slide_width" in merged.columns:
        merged["log_slide_width"] = np.log1p(merged["slide_width"])
    if "slide_height" in merged.columns:
        merged["log_slide_height"] = np.log1p(merged["slide_height"])
    if "erbb2_tpm" in merged.columns:
        merged["log_erbb2_tpm"] = np.log1p(merged["erbb2_tpm"])
    return merged


def numeric_summary_rows(stats, rows, columns: list[tuple[str, str]]) -> list[dict[str, object]]:
    output = []
    for col, label in columns:
        if col not in rows.columns:
            continue
        low = rows.loc[rows["clinical_her2_group"] == NEGATIVE_CLASS, col].dropna().astype(float)
        zero = rows.loc[rows["clinical_her2_group"] == POSITIVE_CLASS, col].dropna().astype(float)
        p_value = float(stats.mannwhitneyu(low, zero, alternative="two-sided").pvalue) if len(low) and len(zero) else math.nan
        output.append(
            {
                "covariate": col,
                "covariate_label": label,
                "n_low": len(low),
                "n_zero": len(zero),
                "mean_low": float(low.mean()) if len(low) else math.nan,
                "mean_zero": float(zero.mean()) if len(zero) else math.nan,
                "median_low": float(low.median()) if len(low) else math.nan,
                "median_zero": float(zero.median()) if len(zero) else math.nan,
                "delta_low_minus_zero": float(low.mean() - zero.mean()) if len(low) and len(zero) else math.nan,
                "mannwhitney_p_value": p_value,
            }
        )
    return output


def categorical_balance_rows(pd, stats, rows, columns: list[tuple[str, str]]) -> list[dict[str, object]]:
    output = []
    for col, label in columns:
        if col not in rows.columns:
            continue
        table = pd.crosstab(rows["clinical_her2_group"], rows[col])
        table = table.reindex(LOW_ZERO_GROUPS).fillna(0)
        if table.shape[1] <= 1:
            p_value = math.nan
        else:
            p_value = float(stats.chi2_contingency(table.to_numpy())[1])
        output.append(
            {
                "covariate": col,
                "covariate_label": label,
                "levels": "; ".join(f"{level}: {int(table[level].sum())}" for level in table.columns),
                "chi_square_p_value": p_value,
            }
        )
    return output


def one_hot_design(pd, rows, categorical_cols: list[str], numeric_cols: list[str]):
    parts = []
    for col in numeric_cols:
        if col in rows.columns:
            parts.append(pd.to_numeric(rows[col], errors="coerce").rename(col))
    for col in categorical_cols:
        if col in rows.columns:
            dummies = pd.get_dummies(rows[col].fillna("Unknown").astype(str), prefix=col, drop_first=True, dtype=float)
            parts.append(dummies)
    if not parts:
        return pd.DataFrame(index=rows.index)
    design = pd.concat(parts, axis=1)
    return design.loc[:, ~design.columns.duplicated()].astype(float)


def design_matrix_for_set(pd, rows, feature_set: str):
    clinical_cat = ["er_status", "pr_status", "histology_group", "stage_group"]
    site_cat = clinical_cat + ["tss_group"]
    clinical_num = []
    site_num = ["log_slide_file_size_mb", "log_slide_width", "log_slide_height"]
    tissue_num = ["mean_tissue_fraction", "retained_fraction", "n_tiles_retained"]
    image_cols = [f"mean_{channel}" for channel in GIGATIME_CHANNELS if f"mean_{channel}" in rows.columns]
    if feature_set == "clinical_basic":
        return one_hot_design(pd, rows, clinical_cat, clinical_num)
    if feature_set == "slide_size_only":
        return one_hot_design(pd, rows, [], site_num)
    if feature_set == "source_site_only":
        return one_hot_design(pd, rows, ["tss_group"], [])
    if feature_set == "site_slide_only":
        return one_hot_design(pd, rows, ["tss_group"], site_num)
    if feature_set == "clinical_site_slide":
        return one_hot_design(pd, rows, site_cat, site_num)
    if feature_set == "gigatime_mean_channels":
        return rows[image_cols].astype(float)
    if feature_set == "gigatime_plus_clinical_basic":
        return pd.concat([rows[image_cols].astype(float), one_hot_design(pd, rows, clinical_cat, clinical_num)], axis=1)
    if feature_set == "gigatime_plus_clinical_site_slide":
        return pd.concat([rows[image_cols].astype(float), one_hot_design(pd, rows, site_cat, site_num)], axis=1)
    if feature_set == "gigatime_plus_clinical_site_slide_tissue":
        return pd.concat([rows[image_cols].astype(float), one_hot_design(pd, rows, site_cat, site_num + tissue_num)], axis=1)
    raise ValueError(feature_set)


def evaluate_cv(optimize, stats, x: np.ndarray, y: np.ndarray, folds, l2_penalty: float) -> dict[str, float]:
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


def run_classifier_covariate_checks(pd, optimize, stats, rows, args: argparse.Namespace):
    rng = np.random.default_rng(args.seed)
    feature_sets = [
        ("clinical_basic", "Clinical covariates"),
        ("slide_size_only", "Slide-size covariates"),
        ("source_site_only", "Source-site covariates"),
        ("site_slide_only", "Source-site + slide-size covariates"),
        ("clinical_site_slide", "Clinical + site/slide covariates"),
        ("gigatime_mean_channels", "GigaTIME mean channels"),
        ("gigatime_plus_clinical_basic", "GigaTIME + clinical covariates"),
        ("gigatime_plus_clinical_site_slide", "GigaTIME + clinical + site/slide"),
        ("gigatime_plus_clinical_site_slide_tissue", "GigaTIME + clinical + site/slide + tissue QC"),
    ]
    output = []
    for view in VIEW_ORDER:
        view_rows = rows.loc[
            (rows["feature_view"] == view) & (rows["clinical_her2_group"].isin(LOW_ZERO_GROUPS))
        ].copy()
        if view_rows.empty:
            continue
        view_rows = view_rows.reset_index(drop=True)
        y = view_rows["clinical_her2_group"].map({NEGATIVE_CLASS: 0, POSITIVE_CLASS: 1}).to_numpy(dtype=int)
        folds = make_repeated_stratified_folds(y, args.folds, args.repeats, rng)
        for feature_set, feature_set_label in feature_sets:
            design = design_matrix_for_set(pd, view_rows, feature_set)
            if design.empty:
                continue
            x = design.to_numpy(dtype=float)
            metrics = evaluate_cv(optimize, stats, x, y, folds, args.l2_penalty)
            output.append(
                {
                    "feature_view": view,
                    "feature_view_label": VIEW_LABELS.get(view, view),
                    "feature_set": feature_set,
                    "feature_set_label": feature_set_label,
                    "n_cases": len(view_rows),
                    "n_features": design.shape[1],
                    "accuracy": metrics["accuracy"],
                    "balanced_accuracy": metrics["balanced_accuracy"],
                    "macro_auc_ovr": metrics["macro_auc_ovr"],
                    "sensitivity": metrics["sensitivity"],
                    "specificity": metrics["specificity"],
                }
            )
    return pd.DataFrame(output)


def ols_group_effect(stats, y, group_low, covariates):
    y = np.asarray(y, dtype=float)
    group_low = np.asarray(group_low, dtype=float)
    covariates = np.asarray(covariates, dtype=float)
    if covariates.ndim == 1:
        covariates = covariates.reshape(-1, 1)
    x = np.column_stack([np.ones(len(y)), group_low, covariates])
    keep = np.isfinite(y) & np.all(np.isfinite(x), axis=1)
    y = y[keep]
    x = x[keep]
    if len(y) <= x.shape[1] + 1:
        return math.nan, math.nan, int(len(y)), int(x.shape[1])
    rank = np.linalg.matrix_rank(x)
    if rank < x.shape[1]:
        return math.nan, math.nan, int(len(y)), int(x.shape[1])
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residual = y - x @ beta
    dof = len(y) - x.shape[1]
    mse = float((residual @ residual) / dof)
    cov_beta = mse * np.linalg.inv(x.T @ x)
    se = math.sqrt(max(float(cov_beta[1, 1]), 0.0))
    if se == 0:
        return float(beta[1]), math.nan, int(len(y)), int(x.shape[1])
    t_stat = float(beta[1] / se)
    p_value = float(2 * stats.t.sf(abs(t_stat), dof))
    return float(beta[1]), p_value, int(len(y)), int(x.shape[1])


def covariate_design_for_adjustment(pd, rows, model: str):
    if model == "unadjusted":
        return pd.DataFrame(index=rows.index)
    if model == "er_pr":
        return one_hot_design(pd, rows, ["er_status", "pr_status"], [])
    if model == "clinical":
        return one_hot_design(pd, rows, ["er_status", "pr_status", "histology_group", "stage_group"], [])
    if model == "clinical_slide":
        return one_hot_design(
            pd,
            rows,
            ["er_status", "pr_status", "histology_group", "stage_group"],
            ["log_slide_file_size_mb", "log_slide_width", "log_slide_height"],
        )
    if model == "clinical_site_slide":
        return one_hot_design(
            pd,
            rows,
            ["er_status", "pr_status", "histology_group", "stage_group", "tss_group"],
            ["log_slide_file_size_mb", "log_slide_width", "log_slide_height"],
        )
    raise ValueError(model)


def run_adjusted_channel_tests(pd, stats, rows):
    model_labels = {
        "unadjusted": "Unadjusted",
        "er_pr": "ER/PR adjusted",
        "clinical": "ER/PR + histology + stage",
        "clinical_slide": "Clinical + slide size",
        "clinical_site_slide": "Clinical + site/slide size",
    }
    output = []
    for view in VIEW_ORDER:
        view_rows = rows.loc[
            (rows["feature_view"] == view) & (rows["clinical_her2_group"].isin(LOW_ZERO_GROUPS))
        ].copy()
        if view_rows.empty:
            continue
        group_low = (view_rows["clinical_her2_group"] == NEGATIVE_CLASS).astype(float).to_numpy()
        for model, model_label in model_labels.items():
            covariates = covariate_design_for_adjustment(pd, view_rows, model)
            if covariates.empty:
                covariate_matrix = np.empty((len(view_rows), 0))
            else:
                covariate_matrix = covariates.to_numpy(dtype=float)
            p_values = []
            row_indices = []
            for channel in KEY_CHANNELS:
                col = f"mean_{channel}"
                if col not in view_rows.columns:
                    continue
                beta, p_value, n, n_parameters = ols_group_effect(
                    stats,
                    view_rows[col].to_numpy(dtype=float),
                    group_low,
                    covariate_matrix,
                )
                row_indices.append(len(output))
                p_values.append(p_value)
                output.append(
                    {
                        "feature_view": view,
                        "feature_view_label": VIEW_LABELS.get(view, view),
                        "model": model,
                        "model_label": model_label,
                        "channel": channel,
                        "beta_low_minus_zero": beta,
                        "p_value": p_value,
                        "n_cases": n,
                        "n_parameters": n_parameters,
                        "n_covariates": covariate_matrix.shape[1],
                    }
                )
            q_values = benjamini_hochberg(p_values)
            for idx, q_value in zip(row_indices, q_values):
                output[idx]["q_value_bh_within_model_view"] = q_value
    return pd.DataFrame(output)


def plot_classifier_results(plt, sns, classifier, asset_dir: Path) -> None:
    plot_df = classifier.copy()
    plt.figure(figsize=(12.0, 6.0))
    sns.barplot(
        data=plot_df,
        x="feature_view_label",
        y="balanced_accuracy",
        hue="feature_set_label",
        order=[VIEW_LABELS[view] for view in VIEW_ORDER],
    )
    plt.axhline(0.5, color="#374151", linestyle="--", linewidth=1)
    plt.ylim(0, 1)
    plt.xlabel("GigaTIME proxy view")
    plt.ylabel("Repeated-CV balanced accuracy")
    plt.title("HER2-Low vs HER2-Zero: GigaTIME Versus Clinical/Site Covariates")
    plt.xticks(rotation=25, ha="right")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "clinical_covariate_classifier_comparison.png", dpi=180)
    plt.close()


def plot_adjusted_q_counts(plt, sns, adjusted, asset_dir: Path) -> None:
    counts = (
        adjusted.assign(significant=adjusted["q_value_bh_within_model_view"] < 0.05)
        .groupby(["feature_view_label", "model_label"])["significant"]
        .sum()
        .reset_index(name="n_q_lt_0_05_channels")
    )
    plt.figure(figsize=(11.0, 5.8))
    sns.barplot(
        data=counts,
        x="feature_view_label",
        y="n_q_lt_0_05_channels",
        hue="model_label",
    )
    plt.xlabel("GigaTIME proxy view")
    plt.ylabel("Key channels with adjusted q < 0.05")
    plt.title("Low-vs-Zero Channel Signal After Clinical/Site Adjustment")
    plt.xticks(rotation=25, ha="right")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "clinical_covariate_adjusted_q_counts.png", dpi=180)
    plt.close()


def plot_covariate_balance(plt, sns, low_zero_rows, asset_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))
    sns.countplot(data=low_zero_rows, x="histology_group", hue="clinical_her2_group", ax=axes[0, 0])
    axes[0, 0].set_title("Histology")
    sns.countplot(data=low_zero_rows, x="stage_group", hue="clinical_her2_group", ax=axes[0, 1])
    axes[0, 1].set_title("Pathologic stage")
    sns.countplot(data=low_zero_rows, x="tss_group", hue="clinical_her2_group", ax=axes[1, 0])
    axes[1, 0].set_title("TCGA source-site group")
    axes[1, 0].tick_params(axis="x", rotation=45)
    sns.boxplot(data=low_zero_rows, x="clinical_her2_group", y="slide_file_size_mb", ax=axes[1, 1])
    axes[1, 1].set_title("Slide file size")
    for axis in axes.flat:
        axis.set_xlabel("")
    fig.suptitle("HER2-Low vs HER2-Zero Covariate Balance")
    fig.tight_layout()
    fig.savefig(asset_dir / "clinical_covariate_balance.png", dpi=180)
    plt.close(fig)


def table_rows_from_balance(numeric, categorical) -> tuple[list[list[str]], list[list[str]]]:
    numeric_rows = [
        [
            row["covariate_label"],
            str(row["n_low"]),
            str(row["n_zero"]),
            fmt(row["mean_low"], 3),
            fmt(row["mean_zero"], 3),
            fmt(row["delta_low_minus_zero"], 3),
            fmt(row["mannwhitney_p_value"], 4),
        ]
        for _, row in numeric.iterrows()
    ]
    categorical_rows = [
        [
            row["covariate_label"],
            row["levels"],
            fmt(row["chi_square_p_value"], 4),
        ]
        for _, row in categorical.iterrows()
    ]
    return numeric_rows, categorical_rows


def classifier_table_rows(classifier, view: str) -> list[list[str]]:
    selected = classifier.loc[classifier["feature_view"] == view].copy()
    return [
        [
            row["feature_set_label"],
            str(int(row["n_features"])),
            fmt(row["balanced_accuracy"], 3),
            fmt(row["macro_auc_ovr"], 3),
            fmt(row["sensitivity"], 3),
            fmt(row["specificity"], 3),
        ]
        for _, row in selected.iterrows()
    ]


def adjusted_summary_rows(adjusted, view: str) -> list[list[str]]:
    selected = adjusted.loc[adjusted["feature_view"] == view]
    rows = []
    for model_label, group in selected.groupby("model_label", sort=False):
        sig = group.loc[group["q_value_bh_within_model_view"] < 0.05]
        rows.append(
            [
                model_label,
                str(len(sig)),
                ", ".join(sig["channel"].tolist()) if not sig.empty else "none",
                fmt(group["q_value_bh_within_model_view"].min(), 4),
            ]
        )
    return rows


def write_markdown(path: Path, asset_dir: Path, numeric, categorical, classifier, adjusted) -> None:
    numeric_rows, categorical_rows = table_rows_from_balance(numeric, categorical)
    focus_view = "ck_top8_within_slide"
    absolute_view = "absolute_ck_high_q75"
    lines = [
        "# Clinical/Site Covariate Sensitivity",
        "",
        "This analysis asks whether the current HER2-low versus HER2-zero GigaTIME signal could be explained by ordinary clinical, slide, or TCGA source-site covariates.",
        "",
        "Important caveat: this still uses retrospective TCGA metadata. It is a confounder sensitivity check, not external validation.",
        "",
        "## Covariate Balance",
        "",
        markdown_table(
            ["Numeric covariate", "N low", "N zero", "Mean low", "Mean zero", "Low-zero delta", "p"],
            numeric_rows,
        ),
        "",
        markdown_table(["Categorical covariate", "Level totals", "Chi-square p"], categorical_rows),
        "",
        f"![Covariate balance]({str(asset_dir / 'clinical_covariate_balance.png').replace('docs/', '')})",
        "",
        "## Classifier Sensitivity",
        "",
        "The classifier comparison asks whether GigaTIME features add signal beyond non-image clinical/site/slide covariates.",
        "",
        f"### {VIEW_LABELS[focus_view]}",
        "",
        markdown_table(
            ["Feature set", "Features", "Balanced accuracy", "AUC", "Sensitivity", "Specificity"],
            classifier_table_rows(classifier, focus_view),
        ),
        "",
        f"### {VIEW_LABELS[absolute_view]}",
        "",
        markdown_table(
            ["Feature set", "Features", "Balanced accuracy", "AUC", "Sensitivity", "Specificity"],
            classifier_table_rows(classifier, absolute_view),
        ),
        "",
        f"![Classifier comparison]({str(asset_dir / 'clinical_covariate_classifier_comparison.png').replace('docs/', '')})",
        "",
        "## Adjusted Channel Tests",
        "",
        "These models test the HER2-low minus HER2-zero coefficient for each key GigaTIME channel after adding covariates.",
        "",
        f"### {VIEW_LABELS[focus_view]}",
        "",
        markdown_table(["Adjustment model", "q<0.05 channels", "Channels", "Best q"], adjusted_summary_rows(adjusted, focus_view)),
        "",
        f"### {VIEW_LABELS[absolute_view]}",
        "",
        markdown_table(["Adjustment model", "q<0.05 channels", "Channels", "Best q"], adjusted_summary_rows(adjusted, absolute_view)),
        "",
        f"![Adjusted q counts]({str(asset_dir / 'clinical_covariate_adjusted_q_counts.png').replace('docs/', '')})",
        "",
        "## Interpretation",
        "",
        "- HER2-low and HER2-zero are not perfectly balanced for histology, stage, source site, or slide size.",
        "- The source-site and slide-size imbalance is especially important because it can reflect TCGA batch/collection differences rather than tumor biology.",
        "- If source-site or slide-size covariates alone classify well, the image result may be partly confounded by cohort construction or acquisition differences.",
        "- If GigaTIME plus covariates does not clearly improve beyond site/slide covariates alone, the image-derived classifier should not be presented as a strong independent biological model.",
        "- Any result that remains after this adjustment still needs pathologist-reviewed tumor-rich regions and external validation before biological claims.",
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{BASE_RESULT_DIR / 'clinical_covariate_sensitivity/covariate_balance_numeric.csv'}`",
        f"- `{BASE_RESULT_DIR / 'clinical_covariate_sensitivity/covariate_balance_categorical.csv'}`",
        f"- `{BASE_RESULT_DIR / 'clinical_covariate_sensitivity/covariate_classifier_metrics.csv'}`",
        f"- `{BASE_RESULT_DIR / 'clinical_covariate_sensitivity/covariate_adjusted_channel_tests.csv'}`",
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

    features = pd.read_csv(args.slide_features)
    metadata = pd.read_csv(args.high_trust_slides)
    rows = add_covariates(pd, features, metadata, args.min_site_count)
    low_zero_rows = rows.loc[(rows["feature_view"] == VIEW_ORDER[0]) & (rows["clinical_her2_group"].isin(LOW_ZERO_GROUPS))]

    numeric_covariates = [
        ("slide_file_size_mb", "Slide file size MB"),
        ("slide_width", "Slide width"),
        ("slide_height", "Slide height"),
        ("mean_tissue_fraction", "Mean tissue fraction"),
        ("retained_fraction", "Retained tile fraction"),
        ("mean_marker_burden", "Mean marker burden"),
        ("mean_DAPI", "Mean virtual DAPI"),
        ("mean_CK", "Mean virtual CK"),
        ("erbb2_tpm", "ERBB2 TPM subset"),
    ]
    categorical_covariates = [
        ("er_status", "ER status"),
        ("pr_status", "PR status"),
        ("histology_group", "Histology group"),
        ("stage_group", "Pathologic stage group"),
        ("tss_group", "TCGA source-site group"),
    ]
    numeric = pd.DataFrame(numeric_summary_rows(stats, low_zero_rows, numeric_covariates))
    categorical = pd.DataFrame(categorical_balance_rows(pd, stats, low_zero_rows, categorical_covariates))
    classifier = run_classifier_covariate_checks(pd, optimize, stats, rows, args)
    adjusted = run_adjusted_channel_tests(pd, stats, rows)

    numeric.to_csv(out_dir / "covariate_balance_numeric.csv", index=False)
    categorical.to_csv(out_dir / "covariate_balance_categorical.csv", index=False)
    classifier.to_csv(out_dir / "covariate_classifier_metrics.csv", index=False)
    adjusted.to_csv(out_dir / "covariate_adjusted_channel_tests.csv", index=False)
    metadata_json = {
        "task": "her2_low_vs_zero",
        "folds": args.folds,
        "repeats": args.repeats,
        "seed": args.seed,
        "l2_penalty": args.l2_penalty,
        "min_site_count": args.min_site_count,
    }
    (out_dir / "clinical_covariate_sensitivity_metadata.json").write_text(
        json.dumps(metadata_json, indent=2) + "\n",
        encoding="utf-8",
    )

    plot_covariate_balance(plt, sns, low_zero_rows, asset_dir)
    plot_classifier_results(plt, sns, classifier, asset_dir)
    plot_adjusted_q_counts(plt, sns, adjusted, asset_dir)
    write_markdown(Path(args.out_markdown), asset_dir, numeric, categorical, classifier, adjusted)

    print(f"Wrote clinical covariate sensitivity outputs to {out_dir}")
    print(f"Wrote clinical covariate sensitivity figures to {asset_dir}")
    print(f"Wrote clinical covariate sensitivity markdown to {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
