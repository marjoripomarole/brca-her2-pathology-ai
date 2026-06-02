#!/usr/bin/env python3
"""Source-site held-out generalization for HER2-low versus HER2-zero classifiers."""

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
    fmt,
    make_repeated_stratified_folds,
    markdown_table,
    metric_dict,
)
from analyze_clinical_covariate_sensitivity import (
    VIEW_LABELS,
    add_covariates,
    one_hot_design,
    require_analysis_libs,
)
from train_her2_classifier_baseline import fit_predict_logistic, standardize_train_test
from train_her2_cleaned_classifier_comparison import GIGATIME_CHANNELS


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
FEATURE_SETS = [
    ("slide_size_only", "Slide-size covariates"),
    ("tissue_qc_only", "Tissue/QC covariates"),
    ("gigatime_mean_channels", "GigaTIME mean channels"),
    ("gigatime_plus_slide_size", "GigaTIME + slide-size"),
    ("gigatime_plus_tissue_qc", "GigaTIME + tissue/QC"),
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
        default=str(BASE_RESULT_DIR / "source_site_generalization"),
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_source_site_generalization",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_source_site_generalization.md",
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260602)
    parser.add_argument("--l2-penalty", type=float, default=1.0)
    parser.add_argument("--min-site-count", type=int, default=5)
    return parser.parse_args()


def design_for_feature_set(pd, rows, feature_set: str):
    image_cols = [f"mean_{channel}" for channel in GIGATIME_CHANNELS if f"mean_{channel}" in rows.columns]
    slide_size_cols = ["log_slide_file_size_mb", "log_slide_width", "log_slide_height"]
    tissue_cols = ["n_tiles_retained", "retained_fraction", "mean_tissue_fraction", "mean_marker_burden"]
    if feature_set == "slide_size_only":
        return one_hot_design(pd, rows, [], slide_size_cols)
    if feature_set == "tissue_qc_only":
        return one_hot_design(pd, rows, [], tissue_cols)
    if feature_set == "gigatime_mean_channels":
        return rows[image_cols].astype(float)
    if feature_set == "gigatime_plus_slide_size":
        return pd.concat([rows[image_cols].astype(float), one_hot_design(pd, rows, [], slide_size_cols)], axis=1)
    if feature_set == "gigatime_plus_tissue_qc":
        return pd.concat([rows[image_cols].astype(float), one_hot_design(pd, rows, [], tissue_cols)], axis=1)
    raise ValueError(feature_set)


def evaluate_folds(optimize, stats, rows, design, y, folds, feature_set: str, validation_scheme: str, l2_penalty: float):
    classes = [NEGATIVE_CLASS, POSITIVE_CLASS]
    x = design.to_numpy(dtype=float)
    true_labels = []
    pred_labels = []
    positive_probs = []
    prediction_rows = []
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
            source = rows.iloc[row_idx]
            prediction_rows.append(
                {
                    "validation_scheme": validation_scheme,
                    "feature_view": source["feature_view"],
                    "feature_view_label": source["feature_view_label"],
                    "feature_set": feature_set,
                    "fold_index": fold_idx,
                    "heldout_tss_codes": ",".join(heldout_sites),
                    "case_submitter_id": source["case_submitter_id"],
                    "clinical_her2_group": source["clinical_her2_group"],
                    "tss_code": source["tss_code"],
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
    metrics["n_folds_run"] = len({row["fold_index"] for row in prediction_rows})
    return metrics, prediction_rows


def leave_source_site_out_folds(rows, y):
    folds = []
    sites = rows["tss_code"].astype(str).to_numpy()
    for site in sorted(set(sites)):
        test_idx = np.where(sites == site)[0]
        train_idx = np.where(sites != site)[0]
        if len(test_idx) == 0 or len(np.unique(y[train_idx])) < 2:
            continue
        folds.append((train_idx, test_idx))
    return folds


def source_site_balance_rows(pd, rows):
    table = pd.crosstab(rows["tss_code"], rows["clinical_her2_group"]).reindex(columns=LOW_ZERO_GROUPS, fill_value=0)
    table["n_cases"] = table.sum(axis=1)
    table = table.sort_values(["n_cases", POSITIVE_CLASS, NEGATIVE_CLASS], ascending=False)
    output = []
    for site, row in table.iterrows():
        output.append(
            {
                "tss_code": site,
                f"n_{NEGATIVE_CLASS}": int(row.get(NEGATIVE_CLASS, 0)),
                f"n_{POSITIVE_CLASS}": int(row.get(POSITIVE_CLASS, 0)),
                "n_cases": int(row["n_cases"]),
                "has_both_classes": int(row.get(NEGATIVE_CLASS, 0) > 0 and row.get(POSITIVE_CLASS, 0) > 0),
            }
        )
    return pd.DataFrame(output)


def site_prediction_summary(pd, predictions):
    if predictions.empty:
        return predictions
    rows = []
    for (view, feature_set, site), group in predictions.groupby(["feature_view", "feature_set", "tss_code"]):
        counts = group["clinical_her2_group"].value_counts()
        rows.append(
            {
                "feature_view": view,
                "feature_view_label": group["feature_view_label"].iloc[0],
                "feature_set": feature_set,
                "feature_set_label": dict(FEATURE_SETS).get(feature_set, feature_set),
                "tss_code": site,
                f"n_{NEGATIVE_CLASS}": int(counts.get(NEGATIVE_CLASS, 0)),
                f"n_{POSITIVE_CLASS}": int(counts.get(POSITIVE_CLASS, 0)),
                "n_cases": int(len(group)),
                "accuracy": float(group["correct"].mean()),
                f"mean_prob_{POSITIVE_CLASS}": float(group[f"prob_{POSITIVE_CLASS}"].mean()),
            }
        )
    return pd.DataFrame(rows)


def run_analysis(pd, optimize, stats, args):
    rng = np.random.default_rng(args.seed)
    features = pd.read_csv(args.slide_features)
    metadata = pd.read_csv(args.high_trust_slides)
    features = add_covariates(pd, features, metadata, args.min_site_count)
    low_zero = features.loc[features["clinical_her2_group"].isin(LOW_ZERO_GROUPS)].copy()
    low_zero["tss_code"] = low_zero["case_submitter_id"].astype(str).str.split("-").str[1]

    site_balance = source_site_balance_rows(pd, low_zero.drop_duplicates(["case_submitter_id", "slide_id"]))
    metrics_rows = []
    all_predictions = []
    for view in VIEW_ORDER:
        rows = low_zero.loc[low_zero["feature_view"] == view].copy().reset_index(drop=True)
        if rows.empty:
            continue
        y = rows["clinical_her2_group"].map({NEGATIVE_CLASS: 0, POSITIVE_CLASS: 1}).to_numpy(dtype=int)
        random_folds = make_repeated_stratified_folds(y, args.folds, args.repeats, rng)
        site_folds = leave_source_site_out_folds(rows, y)
        for feature_set, feature_set_label in FEATURE_SETS:
            design = design_for_feature_set(pd, rows, feature_set)
            if design.empty:
                continue
            for scheme, scheme_label, folds in [
                ("repeated_stratified_cv", "Repeated stratified CV", random_folds),
                ("leave_source_site_out", "Leave source site out", site_folds),
            ]:
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
                all_predictions.extend(predictions)
                metrics_rows.append(
                    {
                        "feature_view": view,
                        "feature_view_label": VIEW_LABELS.get(view, view),
                        "feature_set": feature_set,
                        "feature_set_label": feature_set_label,
                        "validation_scheme": scheme,
                        "validation_scheme_label": scheme_label,
                        "n_cases": len(rows),
                        "n_features": design.shape[1],
                        "n_folds_available": len(folds),
                        **metrics,
                    }
                )
    metrics = pd.DataFrame(metrics_rows)
    predictions = pd.DataFrame(all_predictions)
    site_summary = site_prediction_summary(pd, predictions.loc[predictions["validation_scheme"] == "leave_source_site_out"])
    return site_balance, metrics, predictions, site_summary


def plot_metric_comparison(plt, sns, metrics, asset_dir: Path) -> None:
    focus = metrics.loc[
        (metrics["feature_view"] == "ck_top8_within_slide")
        & (metrics["feature_set"].isin(["slide_size_only", "tissue_qc_only", "gigatime_mean_channels", "gigatime_plus_slide_size"]))
    ].copy()
    plt.figure(figsize=(10.8, 5.6))
    sns.barplot(data=focus, x="feature_set_label", y="balanced_accuracy", hue="validation_scheme_label")
    plt.axhline(0.5, color="#374151", linestyle="--", linewidth=1)
    plt.ylim(0, 1)
    plt.xlabel("Feature set")
    plt.ylabel("Balanced accuracy")
    plt.title("HER2-Low vs HER2-Zero Generalization Across TCGA Source Sites")
    plt.xticks(rotation=20, ha="right")
    plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(asset_dir / "source_site_generalization_balanced_accuracy.png", dpi=180)
    plt.close()


def plot_view_drop(plt, sns, metrics, asset_dir: Path) -> None:
    focus = metrics.loc[metrics["feature_set"] == "gigatime_mean_channels"].copy()
    pivot = focus.pivot_table(
        index=["feature_view", "feature_view_label"],
        columns="validation_scheme",
        values="balanced_accuracy",
        aggfunc="first",
    ).reset_index()
    pivot["ba_drop_site_holdout_minus_random_cv"] = (
        pivot["leave_source_site_out"] - pivot["repeated_stratified_cv"]
    )
    plt.figure(figsize=(10.2, 5.2))
    sns.barplot(data=pivot, x="feature_view_label", y="ba_drop_site_holdout_minus_random_cv")
    plt.axhline(0, color="#374151", linewidth=1)
    plt.xlabel("Feature view")
    plt.ylabel("Source-site holdout BA minus random CV BA")
    plt.title("GigaTIME Generalization Drop Under Source-Site Holdout")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(asset_dir / "source_site_generalization_drop_by_view.png", dpi=180)
    plt.close()


def metric_table_rows(metrics, view: str = "ck_top8_within_slide"):
    selected = metrics.loc[metrics["feature_view"] == view].copy()
    order = {name: idx for idx, (name, _) in enumerate(FEATURE_SETS)}
    scheme_order = {"repeated_stratified_cv": 0, "leave_source_site_out": 1}
    selected["_feature_order"] = selected["feature_set"].map(order)
    selected["_scheme_order"] = selected["validation_scheme"].map(scheme_order)
    selected = selected.sort_values(["_feature_order", "_scheme_order"])
    rows = []
    for _, row in selected.iterrows():
        rows.append(
            [
                row["feature_set_label"],
                row["validation_scheme_label"],
                str(int(row["n_features"])),
                fmt(row["balanced_accuracy"], 3),
                fmt(row["macro_auc_ovr"], 3),
                fmt(row["sensitivity"], 3),
                fmt(row["specificity"], 3),
            ]
        )
    return rows


def site_balance_table_rows(site_balance, limit: int = 12):
    return [
        [
            row["tss_code"],
            str(int(row[f"n_{NEGATIVE_CLASS}"])),
            str(int(row[f"n_{POSITIVE_CLASS}"])),
            str(int(row["n_cases"])),
            "yes" if int(row["has_both_classes"]) else "no",
        ]
        for _, row in site_balance.head(limit).iterrows()
    ]


def write_markdown(path: Path, asset_dir: Path, site_balance, metrics) -> None:
    top8 = metrics.loc[metrics["feature_view"] == "ck_top8_within_slide"].copy()
    gigatime_random = top8.loc[
        (top8["feature_set"] == "gigatime_mean_channels")
        & (top8["validation_scheme"] == "repeated_stratified_cv")
    ].iloc[0]
    gigatime_site = top8.loc[
        (top8["feature_set"] == "gigatime_mean_channels")
        & (top8["validation_scheme"] == "leave_source_site_out")
    ].iloc[0]
    lines = [
        "# Source-Site Held-Out Generalization",
        "",
        "This analysis asks whether the HER2-low versus HER2-zero classifier travels across TCGA source sites. It compares ordinary repeated stratified cross-validation with a harsher leave-one-source-site-out validation.",
        "",
        "Important caveat: many TCGA source sites have only HER2-low or only HER2-zero cases. Leave-source-site-out validation is therefore conservative and can be unstable, but it directly tests the acquisition/source-site confounding concern.",
        "",
        "## Source-Site Balance",
        "",
        markdown_table(
            ["TSS", f"N {NEGATIVE_CLASS}", f"N {POSITIVE_CLASS}", "N cases", "Both classes"],
            site_balance_table_rows(site_balance),
        ),
        "",
        "## Top 8 CK Proxy View",
        "",
        markdown_table(
            ["Feature set", "Validation", "Features", "Balanced accuracy", "AUC", "Sensitivity", "Specificity"],
            metric_table_rows(metrics),
        ),
        "",
        f"GigaTIME mean channels drop from balanced accuracy {fmt(gigatime_random['balanced_accuracy'], 3)} under repeated stratified CV to {fmt(gigatime_site['balanced_accuracy'], 3)} under leave-source-site-out validation in the top 8 CK proxy view.",
        "",
        f"![Source-site generalization balanced accuracy]({str(asset_dir / 'source_site_generalization_balanced_accuracy.png').replace('docs/', '')})",
        "",
        f"![Source-site generalization drop by view]({str(asset_dir / 'source_site_generalization_drop_by_view.png').replace('docs/', '')})",
        "",
        "## Interpretation",
        "",
        "- GigaTIME mean channels lose performance under source-site holdout across every tested feature view.",
        "- Slide-size covariates remain very strong even when entire source sites are held out, meaning the low-versus-zero cohort still carries a portable technical/size imbalance.",
        "- Adding slide-size covariates to GigaTIME largely preserves the slide-size signal, but that does not make the image model more biologically trustworthy.",
        "- The safest conclusion is that the low-versus-zero GigaTIME classifier remains hypothesis-generating and internally interesting, but it is not yet robust evidence of source-independent HER2 biology.",
        "- The next validation step should be external/site-balanced data or pathologist-approved tumor-rich regions with stronger acquisition controls.",
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{BASE_RESULT_DIR / 'source_site_generalization/source_site_generalization_metrics.csv'}`",
        f"- `{BASE_RESULT_DIR / 'source_site_generalization/source_site_generalization_predictions.csv'}`",
        f"- `{BASE_RESULT_DIR / 'source_site_generalization/source_site_balance.csv'}`",
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
    site_balance, metrics, predictions, site_summary = run_analysis(pd, optimize, stats, args)

    site_balance.to_csv(out_dir / "source_site_balance.csv", index=False)
    metrics.to_csv(out_dir / "source_site_generalization_metrics.csv", index=False)
    predictions.to_csv(out_dir / "source_site_generalization_predictions.csv", index=False)
    site_summary.to_csv(out_dir / "source_site_prediction_summary.csv", index=False)
    metadata = {
        "task": "her2_low_vs_zero",
        "validation_schemes": ["repeated_stratified_cv", "leave_source_site_out"],
        "seed": args.seed,
        "folds": args.folds,
        "repeats": args.repeats,
        "l2_penalty": args.l2_penalty,
    }
    (out_dir / "source_site_generalization_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    plot_metric_comparison(plt, sns, metrics, asset_dir)
    plot_view_drop(plt, sns, metrics, asset_dir)
    write_markdown(Path(args.out_markdown), asset_dir, site_balance, metrics)
    print(f"Wrote source-site generalization outputs to {out_dir}")
    print(f"Wrote source-site generalization figures to {asset_dir}")
    print(f"Wrote source-site generalization markdown to {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
