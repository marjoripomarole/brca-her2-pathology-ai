#!/usr/bin/env python3
"""Within-source-site sensitivity analysis for HER2-low versus HER2-zero."""

from __future__ import annotations

import argparse
import json
import math
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
from analyze_clinical_covariate_sensitivity import VIEW_LABELS, add_covariates, one_hot_design, require_analysis_libs
from train_her2_classifier_baseline import fit_predict_logistic, standardize_train_test
from train_her2_cleaned_classifier_comparison import GIGATIME_CHANNELS


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
KEY_CHANNELS = ["CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CD20", "CK", "Ki67"]
FEATURE_SETS = [
    ("slide_size_only", "Slide-size covariates"),
    ("tissue_qc_only", "Tissue/QC covariates"),
    ("source_site_one_hot", "Source-site one-hot"),
    ("gigatime_key_mean_channels", "GigaTIME key mean channels"),
    ("gigatime_mean_channels", "GigaTIME all mean channels"),
    ("gigatime_plus_slide_size", "GigaTIME + slide-size"),
]


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
        default=str(BASE_RESULT_DIR / "within_source_site_low_zero"),
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_within_source_site",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_within_source_site_low_zero.md",
    )
    parser.add_argument("--folds", type=int, default=4)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260602)
    parser.add_argument("--l2-penalty", type=float, default=1.0)
    return parser.parse_args()


def tss_code(case_submitter_id: object) -> str:
    parts = str(case_submitter_id or "").split("-")
    return parts[1] if len(parts) > 1 else "Unknown"


def mixed_site_balance(pd, rows):
    case_rows = rows.drop_duplicates(["case_submitter_id", "slide_id"]).copy()
    table = pd.crosstab(case_rows["tss_code"], case_rows["clinical_her2_group"]).reindex(columns=LOW_ZERO_GROUPS, fill_value=0)
    table["n_cases"] = table.sum(axis=1)
    table["has_both_classes"] = (table[NEGATIVE_CLASS] > 0) & (table[POSITIVE_CLASS] > 0)
    table = table.sort_values(["has_both_classes", "n_cases"], ascending=False)
    output = table.reset_index().rename(
        columns={
            NEGATIVE_CLASS: f"n_{NEGATIVE_CLASS}",
            POSITIVE_CLASS: f"n_{POSITIVE_CLASS}",
        }
    )
    output["has_both_classes"] = output["has_both_classes"].astype(int)
    return output


def channel_cols(rows) -> list[str]:
    return [f"mean_{channel}" for channel in GIGATIME_CHANNELS if f"mean_{channel}" in rows.columns]


def key_channel_cols(rows) -> list[str]:
    return [f"mean_{channel}" for channel in KEY_CHANNELS if f"mean_{channel}" in rows.columns]


def ols_group_test(pd, stats, valid, channel_col: str) -> dict[str, object]:
    y = valid[channel_col].to_numpy(dtype=float)
    is_zero = (valid["clinical_her2_group"] == POSITIVE_CLASS).astype(float).to_numpy()
    site_dummies = pd.get_dummies(valid["tss_code"].astype(str), drop_first=True, dtype=float)
    x = np.column_stack([np.ones(valid.shape[0]), is_zero, site_dummies.to_numpy(dtype=float)])
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    residuals = y - x @ beta
    dof = valid.shape[0] - x.shape[1]
    if dof <= 0:
        return {
            "beta_zero_vs_low_site_fixed": float(beta[1]),
            "p_value": math.nan,
            "dof": int(dof),
        }
    sigma2 = float((residuals @ residuals) / dof)
    cov = sigma2 * np.linalg.pinv(x.T @ x)
    se = math.sqrt(max(float(cov[1, 1]), 0.0))
    t_stat = float(beta[1] / se) if se > 0 else math.nan
    p_value = float(2 * stats.t.sf(abs(t_stat), dof)) if not math.isnan(t_stat) else math.nan
    return {
        "beta_zero_vs_low_site_fixed": float(beta[1]),
        "p_value": p_value,
        "dof": int(dof),
    }


def run_channel_tests(pd, stats, mixed_rows):
    rows_out: list[dict[str, object]] = []
    for view in VIEW_ORDER:
        view_rows = mixed_rows.loc[mixed_rows["feature_view"] == view].copy()
        if view_rows.empty:
            continue
        for col in channel_cols(view_rows):
            valid = view_rows[["clinical_her2_group", "tss_code", col]].dropna().copy()
            n_low = int((valid["clinical_her2_group"] == NEGATIVE_CLASS).sum())
            n_zero = int((valid["clinical_her2_group"] == POSITIVE_CLASS).sum())
            if n_low == 0 or n_zero == 0 or valid["tss_code"].nunique() < 2:
                test = {"beta_zero_vs_low_site_fixed": math.nan, "p_value": math.nan, "dof": math.nan}
            else:
                test = ols_group_test(pd, stats, valid, col)
            low_mean = float(valid.loc[valid["clinical_her2_group"] == NEGATIVE_CLASS, col].mean()) if n_low else math.nan
            zero_mean = float(valid.loc[valid["clinical_her2_group"] == POSITIVE_CLASS, col].mean()) if n_zero else math.nan
            rows_out.append(
                {
                    "feature_view": view,
                    "feature_view_label": VIEW_LABELS.get(view, view),
                    "channel": col.replace("mean_", ""),
                    "n_cases": int(valid.shape[0]),
                    "n_low": n_low,
                    "n_zero": n_zero,
                    "n_source_sites": int(valid["tss_code"].nunique()),
                    "mean_low": low_mean,
                    "mean_zero": zero_mean,
                    "delta_zero_minus_low_unadjusted": zero_mean - low_mean if n_low and n_zero else math.nan,
                    **test,
                }
            )
    by_view = {}
    for view in VIEW_ORDER:
        idxs = [idx for idx, row in enumerate(rows_out) if row["feature_view"] == view]
        q_values = benjamini_hochberg([float(rows_out[idx]["p_value"]) for idx in idxs])
        by_view[view] = q_values
        for idx, q_value in zip(idxs, q_values):
            rows_out[idx]["q_value_bh_within_view"] = q_value
    return pd.DataFrame(rows_out)


def run_per_site_deltas(pd, mixed_rows):
    rows_out: list[dict[str, object]] = []
    for view in VIEW_ORDER:
        view_rows = mixed_rows.loc[mixed_rows["feature_view"] == view].copy()
        for site, site_rows in view_rows.groupby("tss_code"):
            for col in key_channel_cols(site_rows):
                low = site_rows.loc[site_rows["clinical_her2_group"] == NEGATIVE_CLASS, col].dropna()
                zero = site_rows.loc[site_rows["clinical_her2_group"] == POSITIVE_CLASS, col].dropna()
                rows_out.append(
                    {
                        "feature_view": view,
                        "feature_view_label": VIEW_LABELS.get(view, view),
                        "tss_code": site,
                        "channel": col.replace("mean_", ""),
                        "n_low": int(len(low)),
                        "n_zero": int(len(zero)),
                        "mean_low": float(low.mean()) if len(low) else math.nan,
                        "mean_zero": float(zero.mean()) if len(zero) else math.nan,
                        "delta_zero_minus_low": float(zero.mean() - low.mean()) if len(low) and len(zero) else math.nan,
                    }
                )
    return pd.DataFrame(rows_out)


def design_for_feature_set(pd, rows, feature_set: str):
    image_cols = channel_cols(rows)
    key_cols = key_channel_cols(rows)
    slide_size_cols = ["log_slide_file_size_mb", "log_slide_width", "log_slide_height"]
    tissue_cols = ["n_tiles_retained", "retained_fraction", "mean_tissue_fraction", "mean_marker_burden"]
    if feature_set == "slide_size_only":
        return one_hot_design(pd, rows, [], slide_size_cols)
    if feature_set == "tissue_qc_only":
        return one_hot_design(pd, rows, [], tissue_cols)
    if feature_set == "source_site_one_hot":
        return one_hot_design(pd, rows, ["tss_code"], [])
    if feature_set == "gigatime_key_mean_channels":
        return rows[key_cols].astype(float)
    if feature_set == "gigatime_mean_channels":
        return rows[image_cols].astype(float)
    if feature_set == "gigatime_plus_slide_size":
        return pd.concat([rows[image_cols].astype(float), one_hot_design(pd, rows, [], slide_size_cols)], axis=1)
    raise ValueError(feature_set)


def leave_site_out_folds(rows, y):
    folds = []
    sites = rows["tss_code"].astype(str).to_numpy()
    for site in sorted(set(sites)):
        test_idx = np.where(sites == site)[0]
        train_idx = np.where(sites != site)[0]
        if len(test_idx) == 0 or len(np.unique(y[train_idx])) < 2:
            continue
        if len(np.unique(y[test_idx])) < 2:
            continue
        folds.append((train_idx, test_idx))
    return folds


def evaluate_folds(optimize, stats, rows, design, y, folds, feature_set: str, validation_scheme: str, l2_penalty: float):
    classes = [NEGATIVE_CLASS, POSITIVE_CLASS]
    x = design.to_numpy(dtype=float)
    true_labels = []
    pred_labels = []
    positive_probs = []
    predictions = []
    for fold_idx, (train_idx, test_idx) in enumerate(folds):
        if len(np.unique(y[train_idx])) < 2:
            continue
        x_train, x_test = standardize_train_test(x[train_idx], x[test_idx])
        probs = fit_predict_logistic(optimize, x_train, y[train_idx], x_test, classes, l2_penalty)
        pred = np.argmax(probs, axis=1).astype(int)
        true_labels.extend(y[test_idx].tolist())
        pred_labels.extend(pred.tolist())
        positive_probs.extend(probs[:, 1].astype(float).tolist())
        heldout_sites = sorted(set(rows.iloc[test_idx]["tss_code"].astype(str)))
        for local_idx, row_idx in enumerate(test_idx):
            row = rows.iloc[row_idx]
            predictions.append(
                {
                    "validation_scheme": validation_scheme,
                    "feature_view": row["feature_view"],
                    "feature_view_label": row["feature_view_label"],
                    "feature_set": feature_set,
                    "fold_index": fold_idx,
                    "heldout_tss_codes": ",".join(heldout_sites),
                    "case_submitter_id": row["case_submitter_id"],
                    "clinical_her2_group": row["clinical_her2_group"],
                    "tss_code": row["tss_code"],
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
    metrics["n_cv_predictions"] = len(true_labels)
    metrics["n_folds_run"] = len({row["fold_index"] for row in predictions})
    return metrics, predictions


def run_classifiers(pd, optimize, stats, rng, mixed_rows, args):
    metrics_rows = []
    prediction_rows = []
    for view in VIEW_ORDER:
        rows = mixed_rows.loc[mixed_rows["feature_view"] == view].copy().reset_index(drop=True)
        if rows.empty:
            continue
        y = rows["clinical_her2_group"].map({NEGATIVE_CLASS: 0, POSITIVE_CLASS: 1}).to_numpy(dtype=int)
        random_folds = make_repeated_stratified_folds(y, args.folds, args.repeats, rng)
        site_folds = leave_site_out_folds(rows, y)
        for feature_set, feature_set_label in FEATURE_SETS:
            design = design_for_feature_set(pd, rows, feature_set)
            if design.empty:
                continue
            schemes = [("repeated_stratified_cv", "Repeated stratified CV", random_folds)]
            if feature_set != "source_site_one_hot":
                schemes.append(("leave_mixed_source_site_out", "Leave mixed source site out", site_folds))
            for scheme, scheme_label, folds in schemes:
                metrics, predictions = evaluate_folds(
                    optimize,
                    stats,
                    rows,
                    design,
                    y,
                    folds,
                    feature_set,
                    scheme,
                    args.l2_penalty,
                )
                prediction_rows.extend(predictions)
                metrics_rows.append(
                    {
                        "feature_view": view,
                        "feature_view_label": VIEW_LABELS.get(view, view),
                        "feature_set": feature_set,
                        "feature_set_label": feature_set_label,
                        "validation_scheme": scheme,
                        "validation_scheme_label": scheme_label,
                        "n_cases": int(len(rows)),
                        "n_low": int((y == 0).sum()),
                        "n_zero": int((y == 1).sum()),
                        "n_source_sites": int(rows["tss_code"].nunique()),
                        "n_features": int(design.shape[1]),
                        "n_folds_available": int(len(folds)),
                        **metrics,
                    }
                )
    return pd.DataFrame(metrics_rows), pd.DataFrame(prediction_rows)


def plot_site_balance(pd, plt, sns, site_balance, asset_dir: Path):
    mixed = site_balance.loc[site_balance["has_both_classes"] == 1].copy()
    plot_rows = []
    for _, row in mixed.iterrows():
        plot_rows.append({"tss_code": row["tss_code"], "HER2 group": NEGATIVE_CLASS, "n_cases": row[f"n_{NEGATIVE_CLASS}"]})
        plot_rows.append({"tss_code": row["tss_code"], "HER2 group": POSITIVE_CLASS, "n_cases": row[f"n_{POSITIVE_CLASS}"]})
    plot_rows = pd.DataFrame(plot_rows)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    sns.barplot(data=plot_rows, x="tss_code", y="n_cases", hue="HER2 group", ax=ax, palette=["#0f766e", "#1d4ed8"])
    ax.set_xlabel("TCGA source-site code")
    ax.set_ylabel("Cases")
    ax.set_title("Mixed TCGA source sites used for within-site low-vs-zero analysis")
    fig.tight_layout()
    fig.savefig(asset_dir / "within_source_site_mixed_site_balance.png", dpi=180)
    plt.close(fig)


def plot_channel_q_counts(plt, sns, channel_tests, asset_dir: Path):
    summary = (
        channel_tests.assign(significant=channel_tests["q_value_bh_within_view"] < 0.05)
        .groupby("feature_view_label", as_index=False)["significant"]
        .sum()
        .rename(columns={"significant": "n_q_lt_0_05"})
    )
    fig, ax = plt.subplots(figsize=(9, 4.8))
    sns.barplot(data=summary, x="feature_view_label", y="n_q_lt_0_05", ax=ax, color="#0f766e")
    ax.set_xlabel("")
    ax.set_ylabel("Channels with BH q < 0.05")
    ax.set_title("Site-fixed low-vs-zero GigaTIME channel tests")
    ax.tick_params(axis="x", rotation=18)
    fig.tight_layout()
    fig.savefig(asset_dir / "within_source_site_channel_q_counts.png", dpi=180)
    plt.close(fig)


def plot_classifier_metrics(plt, sns, metrics, asset_dir: Path):
    focus = metrics.loc[
        (metrics["feature_view"] == "ck_top8_within_slide")
        & (
            metrics["feature_set"].isin(
                ["slide_size_only", "tissue_qc_only", "source_site_one_hot", "gigatime_key_mean_channels", "gigatime_mean_channels"]
            )
        )
    ].copy()
    fig, ax = plt.subplots(figsize=(11.2, 5.6))
    sns.barplot(data=focus, x="feature_set_label", y="balanced_accuracy", hue="validation_scheme_label", ax=ax)
    ax.axhline(0.5, color="#374151", linestyle="--", linewidth=1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("")
    ax.set_ylabel("Balanced accuracy")
    ax.set_title("Mixed-source-site classifier sensitivity: top 8 CK proxy view")
    ax.tick_params(axis="x", rotation=18)
    fig.tight_layout()
    fig.savefig(asset_dir / "within_source_site_classifier_sensitivity.png", dpi=180)
    plt.close(fig)


def top_channel_rows(channel_tests, view: str = "ck_top8_within_slide", limit: int = 12) -> list[list[str]]:
    rows = channel_tests.loc[channel_tests["feature_view"] == view].copy()
    rows = rows.sort_values(["q_value_bh_within_view", "p_value"], na_position="last").head(limit)
    return [
        [
            row["feature_view_label"],
            row["channel"],
            int(row["n_low"]),
            int(row["n_zero"]),
            fmt(float(row["beta_zero_vs_low_site_fixed"]), 4),
            fmt(float(row["p_value"]), 4),
            fmt(float(row["q_value_bh_within_view"]), 4),
        ]
        for _, row in rows.iterrows()
    ]


def q_count_rows(channel_tests) -> list[list[str]]:
    rows = []
    for view, group in channel_tests.groupby(["feature_view", "feature_view_label"]):
        _view_id, view_label = view
        rows.append(
            [
                view_label,
                int(group.shape[0]),
                int((group["q_value_bh_within_view"] < 0.05).sum()),
                fmt(float(group["q_value_bh_within_view"].min()), 4),
            ]
        )
    return rows


def classifier_focus_rows(metrics, view: str = "ck_top8_within_slide") -> list[list[str]]:
    focus_sets = ["slide_size_only", "tissue_qc_only", "source_site_one_hot", "gigatime_key_mean_channels", "gigatime_mean_channels"]
    rows = metrics.loc[(metrics["feature_view"] == view) & metrics["feature_set"].isin(focus_sets)].copy()
    order = {name: idx for idx, name in enumerate(focus_sets)}
    rows = rows.sort_values(["validation_scheme", "feature_set"], key=lambda series: series.map(order).fillna(series))
    return [
        [
            row["feature_set_label"],
            row["validation_scheme_label"],
            int(row["n_cases"]),
            int(row["n_features"]),
            fmt(float(row["balanced_accuracy"]), 3),
            fmt(float(row["macro_auc_ovr"]), 3),
            fmt(float(row["sensitivity"]), 3),
            fmt(float(row["specificity"]), 3),
        ]
        for _, row in rows.iterrows()
    ]


def site_balance_rows(site_balance) -> list[list[str]]:
    mixed = site_balance.loc[site_balance["has_both_classes"] == 1]
    return [
        [
            row["tss_code"],
            int(row[f"n_{NEGATIVE_CLASS}"]),
            int(row[f"n_{POSITIVE_CLASS}"]),
            int(row["n_cases"]),
        ]
        for _, row in mixed.iterrows()
    ]


def build_markdown(summary: dict[str, object], site_balance, channel_tests, classifier_metrics, args) -> str:
    return f"""# Within-Source-Site HER2-Low Versus HER2-Zero Sensitivity

Status: sensitivity check asking whether the GigaTIME HER2-low versus HER2-zero signal remains when analysis is restricted to TCGA source sites that contain both groups.

## Bottom Line

Only {summary["n_mixed_source_sites"]} TCGA source sites contain both HER2-low and HER2-zero strict high-trust cases. That mixed-site subset contains {summary["n_mixed_cases"]} cases: {summary["n_mixed_low"]} HER2-low and {summary["n_mixed_zero"]} HER2-zero.

This is a useful but very small sensitivity check. It cannot prove biology, but it helps separate two possibilities:

- If the signal persists within mixed source sites, it is less likely to be only a between-site artifact.
- If the signal collapses within mixed source sites, then the TCGA source-site confounding concern becomes even stronger.

## Mixed Source Sites

{markdown_table(["TSS", "HER2-low", "HER2-zero", "Cases"], site_balance_rows(site_balance))}

![Mixed source-site balance](assets/clinical_her2_high_trust_tile128_within_source_site/within_source_site_mixed_site_balance.png)

## Site-Fixed Channel Tests

These models test each GigaTIME channel with a source-site fixed effect: channel score ~ HER2-low/zero group + TCGA source site. Positive beta means HER2-zero is higher than HER2-low after accounting for source site.

{markdown_table(["Feature view", "Channels tested", "Channels q<0.05", "Best BH q"], q_count_rows(channel_tests))}

Top site-fixed channel tests in the top 8 CK proxy view:

{markdown_table(["Feature view", "Channel", "N low", "N zero", "Beta zero-vs-low", "p", "BH q"], top_channel_rows(channel_tests))}

![Site-fixed channel q counts](assets/clinical_her2_high_trust_tile128_within_source_site/within_source_site_channel_q_counts.png)

## Mixed-Site Classifier Sensitivity

This classifier is run only inside the mixed-source-site subset. Repeated stratified CV is still random internal validation. Leave-mixed-source-site-out is harder: one mixed source site is held out at a time.

{markdown_table(["Feature set", "Validation", "N", "Features", "Balanced accuracy", "AUC", "Sensitivity", "Specificity"], classifier_focus_rows(classifier_metrics))}

![Mixed source-site classifier sensitivity](assets/clinical_her2_high_trust_tile128_within_source_site/within_source_site_classifier_sensitivity.png)

## Interpretation

This analysis should be presented as a stress test, not as a main result. The sample is small and imbalanced because most TCGA source sites contain only HER2-low or only HER2-zero cases. A persistent signal here would support continued investigation, while a weak signal here reinforces the need for external/site-balanced validation.

The strongest conclusion remains cautious: TCGA-BRCA can generate a hypothesis about HER2-low versus HER2-zero tissue context, but it is not clean enough by itself to prove source-independent HER2 biology or clinical diagnostic performance.

## Output Files

- `{args.out_markdown}`
- `{args.out_dir}/within_source_site_summary.json`
- `{args.out_dir}/mixed_source_site_balance.csv`
- `{args.out_dir}/within_source_site_channel_tests.csv`
- `{args.out_dir}/within_source_site_per_site_channel_deltas.csv`
- `{args.out_dir}/within_source_site_classifier_metrics.csv`
- `{args.out_dir}/within_source_site_classifier_predictions.csv`
"""


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    asset_dir = Path(args.asset_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, optimize, stats = require_analysis_libs(out_dir / ".matplotlib")
    rng = np.random.default_rng(args.seed)

    features = pd.read_csv(args.slide_features)
    metadata = pd.read_csv(args.high_trust_slides)
    features = add_covariates(pd, features, metadata, min_site_count=1)
    low_zero = features.loc[features["clinical_her2_group"].isin(LOW_ZERO_GROUPS)].copy()
    low_zero["tss_code"] = low_zero["case_submitter_id"].map(tss_code)

    site_balance = mixed_site_balance(pd, low_zero)
    mixed_sites = set(site_balance.loc[site_balance["has_both_classes"] == 1, "tss_code"].astype(str))
    mixed_rows = low_zero.loc[low_zero["tss_code"].isin(mixed_sites)].copy()
    mixed_unique = mixed_rows.drop_duplicates(["case_submitter_id", "slide_id"])

    channel_tests = run_channel_tests(pd, stats, mixed_rows)
    per_site_deltas = run_per_site_deltas(pd, mixed_rows)
    classifier_metrics, classifier_predictions = run_classifiers(pd, optimize, stats, rng, mixed_rows, args)

    summary = {
        "n_all_low_zero_cases": int(low_zero.drop_duplicates(["case_submitter_id", "slide_id"]).shape[0]),
        "n_mixed_source_sites": int(len(mixed_sites)),
        "mixed_source_sites": sorted(mixed_sites),
        "n_mixed_cases": int(mixed_unique.shape[0]),
        "n_mixed_low": int((mixed_unique["clinical_her2_group"] == NEGATIVE_CLASS).sum()),
        "n_mixed_zero": int((mixed_unique["clinical_her2_group"] == POSITIVE_CLASS).sum()),
        "n_channel_tests_q_lt_0_05": int((channel_tests["q_value_bh_within_view"] < 0.05).sum()) if not channel_tests.empty else 0,
    }

    site_balance.to_csv(out_dir / "mixed_source_site_balance.csv", index=False)
    channel_tests.to_csv(out_dir / "within_source_site_channel_tests.csv", index=False)
    per_site_deltas.to_csv(out_dir / "within_source_site_per_site_channel_deltas.csv", index=False)
    classifier_metrics.to_csv(out_dir / "within_source_site_classifier_metrics.csv", index=False)
    classifier_predictions.to_csv(out_dir / "within_source_site_classifier_predictions.csv", index=False)
    with (out_dir / "within_source_site_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    plot_site_balance(pd, plt, sns, site_balance, asset_dir)
    plot_channel_q_counts(plt, sns, channel_tests, asset_dir)
    plot_classifier_metrics(plt, sns, classifier_metrics, asset_dir)

    markdown = build_markdown(summary, site_balance, channel_tests, classifier_metrics, args)
    Path(args.out_markdown).write_text(markdown, encoding="utf-8")

    print(f"Wrote {args.out_markdown}")
    print(f"Wrote {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
