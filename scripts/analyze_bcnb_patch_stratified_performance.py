#!/usr/bin/env python3
"""Slice BCNB low/zero patch-model performance by clinical covariates."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from analyze_bcnb_patch_embedding_control import (
    NEGATIVE_CLASS,
    POSITIVE_CLASS,
    design_matrix,
    normalize_inputs,
)
from analyze_bcnb_patch_model_comparison import (
    ensemble_from_oof,
    evaluate_prepared_with_oof,
    load_aligned_embeddings,
    prepare_blocks,
)
from analyze_classifier_permutation_sanity import fmt, make_repeated_stratified_folds, markdown_table, metric_dict


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hoptimus-embeddings",
        default="results/bcnb_patch_embeddings_hoptimus0_hash_capped10_low_zero/patient_embeddings.csv",
    )
    parser.add_argument(
        "--virchow2-embeddings",
        default="results/bcnb_patch_embeddings_virchow2_hash_capped10_low_zero/patient_embeddings.csv",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/bcnb_patch_stratified_performance_hoptimus0_virchow2_hash_capped10_low_zero"),
    )
    parser.add_argument(
        "--asset-dir",
        type=Path,
        default=Path("docs/assets/bcnb_patch_stratified_performance_hoptimus0_virchow2_hash_capped10_low_zero"),
    )
    parser.add_argument(
        "--out-markdown",
        type=Path,
        default=Path("docs/bcnb_patch_stratified_performance_hoptimus0_virchow2_hash_capped10_low_zero.md"),
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--pca-components", type=int, default=20)
    parser.add_argument("--min-per-class", type=int, default=15)
    parser.add_argument("--seed", type=int, default=20260604)
    parser.add_argument("--l2-penalty", type=float, default=1.0)
    return parser.parse_args()


def require_libs(mpl_config_dir: Path):
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
        from scipy import optimize, stats
    except ModuleNotFoundError as exc:
        raise SystemExit(f"Missing Python package: {exc.name}. Use `conda run -n gigatime-tcga ...`.") from exc
    sns.set_theme(style="whitegrid", context="notebook")
    return pd, plt, sns, optimize, stats


def patient_predictions_from_oof(pd, oof_rows, rows, model_label: str):
    oof = pd.DataFrame(oof_rows)
    meta_cols = [
        "patient_id",
        "clinical_her2_group",
        "grade",
        "grade_missing",
        "ER",
        "PR",
        "ki67",
        "molecular_subtype",
        "aln_status",
        "n_manifest_patches",
        "n_used_patches",
        "mean_tissue_fraction",
        "min_tissue_fraction",
    ]
    available_meta = [col for col in meta_cols if col in rows.columns]
    meta = rows[available_meta].copy()
    patient = (
        oof.groupby("patient_id", as_index=False)
        .agg(
            y_true=("y_true", "first"),
            mean_prob_her2_zero=("prob_her2_zero", "mean"),
            n_oof_predictions=("prob_her2_zero", "size"),
        )
        .merge(meta, on="patient_id", how="left", validate="one_to_one")
    )
    patient["model"] = model_label
    patient["predicted_label"] = np.where(patient["mean_prob_her2_zero"] >= 0.5, POSITIVE_CLASS, NEGATIVE_CLASS)
    patient["predicted_int"] = np.where(patient["mean_prob_her2_zero"] >= 0.5, 1, 0)
    return patient


def add_strata(pd, rows, min_per_class: int):
    stratum_rows = []

    def add(mask, family: str, label: str):
        subset = rows.loc[mask]
        n_low = int((subset["clinical_her2_group"] == NEGATIVE_CLASS).sum())
        n_zero = int((subset["clinical_her2_group"] == POSITIVE_CLASS).sum())
        if n_low >= min_per_class and n_zero >= min_per_class:
            stratum_rows.append(
                {
                    "stratum_family": family,
                    "stratum_label": label,
                    "n_patients": int(len(subset)),
                    "n_low": n_low,
                    "n_zero": n_zero,
                    "patient_ids": set(subset["patient_id"].tolist()),
                }
            )

    add(rows.index == rows.index, "All", "All low/zero patients")
    add(rows["grade"].notna(), "Grade", "Grade known")
    for grade in [1.0, 2.0, 3.0]:
        add(rows["grade"] == grade, "Grade", f"Grade {int(grade)}")
    add(rows["grade"].isna(), "Grade", "Grade missing")
    for col in ["ER", "PR", "molecular_subtype", "aln_status"]:
        if col not in rows.columns:
            continue
        for value in sorted(rows[col].dropna().astype(str).unique()):
            add(rows[col].astype(str) == value, col, f"{col}: {value}")
    if "ki67" in rows.columns:
        valid = rows["ki67"].notna()
        add(valid & (rows["ki67"] < 0.20), "Ki67", "Ki67 < 20%")
        add(valid & (rows["ki67"] >= 0.20), "Ki67", "Ki67 >= 20%")
        add(valid & (rows["ki67"] >= 0.30), "Ki67", "Ki67 >= 30%")
    for grade in [2.0, 3.0]:
        for er in ["Positive", "Negative"]:
            add((rows["grade"] == grade) & (rows["ER"] == er), "Grade x ER", f"Grade {int(grade)}, ER {er}")
    counts = pd.DataFrame([{k: v for k, v in row.items() if k != "patient_ids"} for row in stratum_rows])
    return stratum_rows, counts


def metrics_for_slice(stats, patient_predictions, patient_ids: set):
    subset = patient_predictions.loc[patient_predictions["patient_id"].isin(patient_ids)].copy()
    if subset.empty or subset["y_true"].nunique() < 2:
        return None
    metrics = metric_dict(
        stats,
        subset["y_true"].to_numpy(dtype=int),
        subset["predicted_int"].to_numpy(dtype=int),
        subset["mean_prob_her2_zero"].to_numpy(dtype=float),
    )
    metrics["n_patients"] = int(len(subset))
    metrics["n_low"] = int((subset["clinical_her2_group"] == NEGATIVE_CLASS).sum())
    metrics["n_zero"] = int((subset["clinical_her2_group"] == POSITIVE_CLASS).sum())
    metrics["mean_prob_zero_low"] = float(
        subset.loc[subset["clinical_her2_group"] == NEGATIVE_CLASS, "mean_prob_her2_zero"].mean()
    )
    metrics["mean_prob_zero_zero"] = float(
        subset.loc[subset["clinical_her2_group"] == POSITIVE_CLASS, "mean_prob_her2_zero"].mean()
    )
    metrics["delta_prob_zero_minus_low"] = metrics["mean_prob_zero_zero"] - metrics["mean_prob_zero_low"]
    return metrics


def build_stratified_metrics(pd, stats, patient_prediction_frames, stratum_rows):
    rows = []
    for predictions in patient_prediction_frames:
        model = predictions["model"].iloc[0]
        for stratum in stratum_rows:
            metrics = metrics_for_slice(stats, predictions, stratum["patient_ids"])
            if metrics is None:
                continue
            rows.append(
                {
                    "model": model,
                    "stratum_family": stratum["stratum_family"],
                    "stratum_label": stratum["stratum_label"],
                    **metrics,
                }
            )
    return pd.DataFrame(rows)


def clinical_associations(pd, stats, rows):
    records = []
    for col in ["grade", "ER", "PR", "molecular_subtype", "aln_status"]:
        tmp = rows.copy()
        tmp[col] = tmp[col].fillna("Missing").astype(str)
        for value, subset in tmp.groupby(col, dropna=False):
            n_low = int((subset["clinical_her2_group"] == NEGATIVE_CLASS).sum())
            n_zero = int((subset["clinical_her2_group"] == POSITIVE_CLASS).sum())
            if n_low + n_zero == 0:
                continue
            records.append(
                {
                    "clinical_field": col,
                    "level": value,
                    "n_low": n_low,
                    "n_zero": n_zero,
                    "zero_fraction": n_zero / (n_low + n_zero),
                }
            )
    if "ki67" in rows.columns:
        low = rows.loc[rows["clinical_her2_group"] == NEGATIVE_CLASS, "ki67"].dropna()
        zero = rows.loc[rows["clinical_her2_group"] == POSITIVE_CLASS, "ki67"].dropna()
        if len(low) and len(zero):
            test = stats.mannwhitneyu(low, zero, alternative="two-sided")
            records.append(
                {
                    "clinical_field": "ki67",
                    "level": "continuous",
                    "n_low": int(len(low)),
                    "n_zero": int(len(zero)),
                    "zero_fraction": float("nan"),
                    "low_mean": float(low.mean()),
                    "zero_mean": float(zero.mean()),
                    "mannwhitney_p_value": float(test.pvalue),
                }
            )
    return pd.DataFrame(records)


def plot_stratified_auc(plt, sns, metrics, asset_dir: Path):
    selected = metrics.loc[
        metrics["model"].isin(["Clinical covariates", "H-Optimus-0 embedding", "Virchow2 embedding", "H-Optimus-0 + Virchow2"])
        & metrics["stratum_family"].isin(["All", "Grade", "ER", "PR", "molecular_subtype"])
    ].copy()
    selected["stratum"] = selected["stratum_label"].str.replace("molecular_subtype: ", "Subtype: ", regex=False)
    pivot = selected.pivot_table(index="stratum", columns="model", values="macro_auc_ovr", aggfunc="first")
    order = [
        "All low/zero patients",
        "Grade known",
        "Grade 2",
        "Grade 3",
        "ER: Positive",
        "ER: Negative",
        "PR: Positive",
        "PR: Negative",
        "Subtype: Luminal A",
        "Subtype: Luminal B",
        "Subtype: Triple negative",
    ]
    pivot = pivot.reindex([item for item in order if item in pivot.index])
    fig, axis = plt.subplots(figsize=(9.5, max(5.2, 0.36 * len(pivot))))
    sns.heatmap(pivot, annot=True, fmt=".2f", cmap="viridis", vmin=0.45, vmax=0.75, linewidths=0.3, ax=axis)
    axis.set_xlabel("Model")
    axis.set_ylabel("Clinical stratum")
    axis.set_title("BCNB HER2-low vs zero patient-mean OOF AUC by stratum")
    fig.tight_layout()
    fig.savefig(asset_dir / "bcnb_patch_stratified_auc_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_clinical_zero_fraction(plt, sns, associations, asset_dir: Path):
    selected = associations.loc[
        associations["clinical_field"].isin(["grade", "ER", "PR", "molecular_subtype", "aln_status"])
    ].copy()
    field_labels = {
        "grade": "Grade",
        "ER": "ER",
        "PR": "PR",
        "molecular_subtype": "Subtype",
        "aln_status": "ALN status",
    }
    selected["field_label"] = selected["clinical_field"].map(field_labels).fillna(selected["clinical_field"])
    selected["label"] = selected["field_label"] + ": " + selected["level"]
    selected = selected.sort_values(["clinical_field", "zero_fraction"])
    fig, axis = plt.subplots(figsize=(8.8, max(5.0, 0.25 * len(selected))))
    sns.barplot(data=selected, y="label", x="zero_fraction", ax=axis)
    axis.set_xlabel("HER2-zero fraction within low/zero cohort")
    axis.set_ylabel("Clinical stratum")
    axis.set_title("BCNB clinical covariate imbalance")
    fig.tight_layout()
    fig.savefig(asset_dir / "bcnb_clinical_zero_fraction.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def metric_table(metrics, labels: list[str]) -> list[list[str]]:
    rows = []
    selected = metrics.loc[
        (metrics["stratum_label"].isin(labels))
        & metrics["model"].isin(["Clinical covariates", "H-Optimus-0 embedding", "Virchow2 embedding", "H-Optimus-0 + Virchow2"])
    ].copy()
    label_order = {label: index for index, label in enumerate(labels)}
    model_order = {
        "Clinical covariates": 0,
        "H-Optimus-0 embedding": 1,
        "Virchow2 embedding": 2,
        "H-Optimus-0 + Virchow2": 3,
    }
    selected["_label_order"] = selected["stratum_label"].map(label_order)
    selected["_model_order"] = selected["model"].map(model_order)
    selected = selected.sort_values(["_label_order", "_model_order"])
    for _, row in selected.iterrows():
        rows.append(
            [
                row["stratum_label"],
                row["model"],
                str(int(row["n_low"])),
                str(int(row["n_zero"])),
                fmt(row["balanced_accuracy"], 3),
                fmt(row["macro_auc_ovr"], 3),
                fmt(row["delta_prob_zero_minus_low"], 3),
            ]
        )
    return rows


def write_markdown(path: Path, asset_dir: Path, args, rows, metrics, associations):
    labels = [
        "All low/zero patients",
        "Grade 2",
        "Grade 3",
        "ER: Positive",
        "ER: Negative",
        "molecular_subtype: Luminal A",
        "molecular_subtype: Luminal B",
        "molecular_subtype: Triple negative",
    ]
    dual = metrics.loc[
        (metrics["model"] == "H-Optimus-0 + Virchow2")
        & (metrics["stratum_label"].isin(["All low/zero patients", "Grade 2", "Grade 3", "ER: Positive", "ER: Negative"]))
    ].copy()
    clinical = metrics.loc[
        (metrics["model"] == "Clinical covariates")
        & (metrics["stratum_label"].isin(["All low/zero patients", "Grade 2", "Grade 3", "ER: Positive", "ER: Negative"]))
    ].copy()
    lines = [
        "# BCNB Patch Embedding Stratified Performance",
        "",
        "Status: clinical-slice robustness check for the BCNB external patch pilots.",
        "",
        "## Method",
        "",
        f"- Cohort: {len(rows)} BCNB low/zero patients ({int((rows['clinical_her2_group'] == NEGATIVE_CLASS).sum())} HER2-low, {int((rows['clinical_her2_group'] == POSITIVE_CLASS).sum())} HER2-zero).",
        "- Inputs: the same hash-capped 10-patch patient embeddings used in the H-Optimus-0 and Virchow2 BCNB pilots.",
        f"- Classifiers: class-balanced logistic regression with repeated stratified {args.folds}-fold CV ({args.repeats} repeats); patient-level probabilities are averaged across out-of-fold repeats before slice scoring.",
        f"- Embedding PCA: {args.pca_components} components fit inside each training fold only.",
        f"- Strata are reported only when both HER2-low and HER2-zero have at least {args.min_per_class} patients.",
        "",
        "## Clinical Imbalance",
        "",
        f"![BCNB clinical zero fraction](assets/{asset_dir.name}/bcnb_clinical_zero_fraction.png)",
        "",
        "Grade, ER/PR, subtype, nodal status, and Ki67 all show some low/zero imbalance. This is exactly why pooled image performance is not sufficient evidence of a HER2-specific visual signal.",
        "",
        "## Stratified Patient-Mean OOF Performance",
        "",
        markdown_table(
            ["Stratum", "Model", "Low", "Zero", "Balanced accuracy", "AUC", "Mean P0 zero-low"],
            metric_table(metrics, labels),
        ),
        "",
        f"![BCNB stratified AUC heatmap](assets/{asset_dir.name}/bcnb_patch_stratified_auc_heatmap.png)",
        "",
        "## Interpretation",
        "",
        "- The pooled dual-model result remains modest, and slice performance is uneven rather than uniformly strong across clinically meaningful subgroups.",
        "- ER-negative and triple-negative slices are particularly important stress tests because they reduce one major receptor-status imbalance. The image models do not become a strong classifier there.",
        "- Grade 2 and grade 3 slices still show some image-readable separation, but the effect is not large enough to support a clinical HER2-low/zero classifier claim.",
        "- Overall, this strengthens the current manuscript framing: BCNB supports a weak, reproducible morphology/covariate-associated signal, not a robust standalone HER2-low versus HER2-zero detector.",
        "",
        "## Selected Dual-Model Versus Clinical Rows",
        "",
        markdown_table(
            ["Stratum", "Dual-model BA", "Dual-model AUC", "Clinical BA", "Clinical AUC"],
            [
                [
                    label,
                    fmt(
                        dual.loc[dual["stratum_label"] == label, "balanced_accuracy"].iloc[0]
                        if not dual.loc[dual["stratum_label"] == label].empty
                        else float("nan"),
                        3,
                    ),
                    fmt(
                        dual.loc[dual["stratum_label"] == label, "macro_auc_ovr"].iloc[0]
                        if not dual.loc[dual["stratum_label"] == label].empty
                        else float("nan"),
                        3,
                    ),
                    fmt(
                        clinical.loc[clinical["stratum_label"] == label, "balanced_accuracy"].iloc[0]
                        if not clinical.loc[clinical["stratum_label"] == label].empty
                        else float("nan"),
                        3,
                    ),
                    fmt(
                        clinical.loc[clinical["stratum_label"] == label, "macro_auc_ovr"].iloc[0]
                        if not clinical.loc[clinical["stratum_label"] == label].empty
                        else float("nan"),
                        3,
                    ),
                ]
                for label in ["All low/zero patients", "Grade 2", "Grade 3", "ER: Positive", "ER: Negative"]
            ],
        ),
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_stratified_patient_metrics.csv'}`",
        f"- `{Path(args.out_dir) / 'bcnb_patch_stratified_patient_predictions.csv'}`",
        f"- `{Path(args.out_dir) / 'bcnb_clinical_stratum_counts.csv'}`",
        f"- `{Path(args.out_dir) / 'bcnb_clinical_associations.csv'}`",
        f"- `{asset_dir}/`",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def run_analysis(pd, plt, sns, optimize, stats, args):
    rng = np.random.default_rng(args.seed)
    rows, h_x, v_x, _h_dim, _v_dim = load_aligned_embeddings(
        pd,
        Path(args.hoptimus_embeddings),
        Path(args.virchow2_embeddings),
    )
    rows = normalize_inputs(pd, rows).reset_index(drop=True)
    y = rows["clinical_her2_group"].map({NEGATIVE_CLASS: 0, POSITIVE_CLASS: 1}).to_numpy(dtype=int)
    patient_ids = rows["patient_id"].to_numpy()
    folds = make_repeated_stratified_folds(y, args.folds, args.repeats, rng)
    clinical_x = design_matrix(pd, rows, "clinical_covariates", []).to_numpy(dtype=float)
    specs = [
        ("Clinical covariates", [("clinical", clinical_x, None)]),
        ("H-Optimus-0 embedding", [("hoptimus0", h_x, args.pca_components)]),
        ("Virchow2 embedding", [("virchow2", v_x, args.pca_components)]),
        (
            "H-Optimus-0 + Virchow2",
            [("hoptimus0", h_x, args.pca_components), ("virchow2", v_x, args.pca_components)],
        ),
        (
            "H-Optimus-0 + Virchow2 + clinical",
            [("clinical", clinical_x, None), ("hoptimus0", h_x, args.pca_components), ("virchow2", v_x, args.pca_components)],
        ),
    ]
    oof_by_model = {}
    prediction_frames = []
    for model_label, blocks in specs:
        prepared = prepare_blocks(blocks, folds)
        _metrics, oof_rows = evaluate_prepared_with_oof(optimize, stats, prepared, y, patient_ids, args.l2_penalty)
        oof_by_model[model_label] = oof_rows
        prediction_frames.append(patient_predictions_from_oof(pd, oof_rows, rows, model_label))
    ensemble_metrics, ensemble_oof = ensemble_from_oof(
        pd,
        stats,
        oof_by_model["H-Optimus-0 embedding"],
        oof_by_model["Virchow2 embedding"],
    )
    del ensemble_metrics
    ensemble_rows = ensemble_oof.rename(
        columns={"ensemble_prob_her2_zero": "prob_her2_zero", "ensemble_pred": "predicted_int"}
    )
    ensemble_rows["predicted_label"] = np.where(ensemble_rows["predicted_int"] == 1, POSITIVE_CLASS, NEGATIVE_CLASS)
    prediction_frames.append(
        patient_predictions_from_oof(
            pd,
            ensemble_rows[["fold_id", "patient_id", "y_true", "prob_her2_zero", "predicted_label"]],
            rows,
            "Average probability ensemble",
        )
    )
    stratum_rows, stratum_counts = add_strata(pd, rows, args.min_per_class)
    metrics = build_stratified_metrics(pd, stats, prediction_frames, stratum_rows)
    patient_predictions = pd.concat(prediction_frames, axis=0, ignore_index=True)
    associations = clinical_associations(pd, stats, rows)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.out_dir / "bcnb_patch_stratified_patient_metrics.csv", index=False)
    patient_predictions.to_csv(args.out_dir / "bcnb_patch_stratified_patient_predictions.csv", index=False)
    stratum_counts.to_csv(args.out_dir / "bcnb_clinical_stratum_counts.csv", index=False)
    associations.to_csv(args.out_dir / "bcnb_clinical_associations.csv", index=False)
    (args.out_dir / "bcnb_patch_stratified_metadata.json").write_text(
        json.dumps(
            {
                "task": "bcnb_patch_stratified_performance",
                "folds": args.folds,
                "repeats": args.repeats,
                "pca_components": args.pca_components,
                "min_per_class": args.min_per_class,
                "seed": args.seed,
                "n_patients": int(len(rows)),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    plot_stratified_auc(plt, sns, metrics, args.asset_dir)
    plot_clinical_zero_fraction(plt, sns, associations, args.asset_dir)
    write_markdown(args.out_markdown, args.asset_dir, args, rows, metrics, associations)
    return metrics, associations


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.asset_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, optimize, stats = require_libs(args.out_dir / ".matplotlib")
    metrics, _associations = run_analysis(pd, plt, sns, optimize, stats, args)
    print(f"Wrote BCNB patch stratified outputs to {args.out_dir}")
    print(f"Wrote BCNB patch stratified markdown to {args.out_markdown}")
    selected = metrics.loc[
        (metrics["stratum_label"].isin(["All low/zero patients", "Grade 2", "Grade 3", "ER: Positive", "ER: Negative"]))
        & (metrics["model"].isin(["Clinical covariates", "H-Optimus-0 embedding", "Virchow2 embedding", "H-Optimus-0 + Virchow2"]))
    ]
    print(selected[["stratum_label", "model", "n_low", "n_zero", "balanced_accuracy", "macro_auc_ovr"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
