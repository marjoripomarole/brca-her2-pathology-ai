#!/usr/bin/env python3
"""Build a simple notebook and HTML report for the current HER2 findings."""

from __future__ import annotations

import argparse
import csv
import html
import json
import shutil
from pathlib import Path


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
BASE_ASSET_DIR = Path("docs/assets")

ASSET_COPIES = [
    (
        BASE_RESULT_DIR / "clinical_summary/clinical_her2_channel_boxplots.png",
        BASE_ASSET_DIR / "clinical_her2_high_trust_tile128_findings/clinical_her2_channel_boxplots.png",
    ),
    (
        BASE_RESULT_DIR / "clinical_summary/clinical_her2_group_mean_heatmap.png",
        BASE_ASSET_DIR / "clinical_her2_high_trust_tile128_findings/clinical_her2_group_mean_heatmap.png",
    ),
    (
        BASE_RESULT_DIR / "clinical_summary/erbb2_tpm_by_clinical_her2_group.png",
        BASE_ASSET_DIR / "clinical_her2_high_trust_tile128_findings/erbb2_tpm_by_clinical_her2_group.png",
    ),
]

FEATURE_VIEW_ORDER = {
    "all_sampled_tissue": 0,
    "qc_cellular_tissue": 1,
    "ck_enriched_top50": 2,
    "ck_enriched_top25": 3,
}

TUMOR_PROXY_VIEW_ORDER = {
    "qc_cellular_tissue": 0,
    "ck_top25_within_slide": 1,
    "ck_top16_within_slide": 2,
    "ck_top8_within_slide": 3,
    "ck_top16_non_low_marker": 4,
    "absolute_ck_high_q75": 5,
}

FEATURE_SET_LABELS = {
    "gigatime_mean_and_fraction_channels": "Mean + fraction channels",
    "gigatime_mean_channels": "Mean channels",
    "interpretable_distribution_features": "Distribution features",
    "interpretable_marker_means": "Interpretable marker means",
    "virtual_programs": "Virtual programs",
}

KEY_CHANNELS = ["CD68", "CK", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CD20", "Ki67"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-notebook", default="notebooks/clinical_her2_findings_simple.ipynb")
    parser.add_argument("--out-html", default="notebooks/clinical_her2_findings_simple.html")
    parser.add_argument(
        "--trust-summary",
        default="docs/assets/clinical_her2_trustworthy_slide_list/trustworthy_summary.json",
    )
    parser.add_argument(
        "--high-trust-slides",
        default="docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv",
    )
    parser.add_argument(
        "--clinical-pairwise",
        default=str(BASE_RESULT_DIR / "clinical_summary/clinical_her2_pairwise_tests.csv"),
    )
    parser.add_argument(
        "--clinical-channels",
        default=str(BASE_RESULT_DIR / "clinical_summary/clinical_her2_channel_summary.csv"),
    )
    parser.add_argument(
        "--cleanup-pairwise",
        default=str(BASE_RESULT_DIR / "gigatime_cleanup/cleanup_pairwise_tests.csv"),
    )
    parser.add_argument(
        "--cleanup-channels",
        default=str(BASE_RESULT_DIR / "gigatime_cleanup/cleanup_channel_summary.csv"),
    )
    parser.add_argument(
        "--classifier-metrics",
        default=str(BASE_RESULT_DIR / "cleaned_classifier_comparison/cleaned_classifier_best_h_e_metrics.csv"),
    )
    parser.add_argument(
        "--erpr-adjusted",
        default="docs/assets/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity/low_zero_erpr_adjusted_tests.csv",
    )
    parser.add_argument(
        "--sensitivity-summary",
        default="docs/assets/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity/sensitivity_summary.json",
    )
    parser.add_argument(
        "--run-agreement-summary",
        default="docs/assets/clinical_her2_high_trust_tile128_vs_expanded20_tile256/run_agreement_summary.json",
    )
    parser.add_argument(
        "--direction-comparison",
        default="docs/assets/clinical_her2_high_trust_tile128_vs_expanded20_tile256/low_zero_direction_comparison.csv",
    )
    parser.add_argument(
        "--case-driver-summary",
        default=str(BASE_RESULT_DIR / "case_driver_analysis/case_driver_summary.json"),
    )
    parser.add_argument(
        "--case-driver-review",
        default=str(BASE_RESULT_DIR / "case_driver_analysis/low_zero_classifier_review_cases.csv"),
    )
    parser.add_argument(
        "--visual-qc-selected-tiles",
        default=str(BASE_RESULT_DIR / "case_driver_visual_qc/case_driver_visual_qc_selected_tiles.csv"),
    )
    parser.add_argument(
        "--tissue-composition-tests",
        default=str(BASE_RESULT_DIR / "tissue_composition_sensitivity/low_zero_composition_tests.csv"),
    )
    parser.add_argument(
        "--tissue-composition-correlations",
        default=str(BASE_RESULT_DIR / "tissue_composition_sensitivity/driver_tissue_composition_correlations.csv"),
    )
    parser.add_argument(
        "--absolute-ck-pairwise",
        default=str(BASE_RESULT_DIR / "tissue_composition_sensitivity/absolute_ck_high_low_zero_pairwise_tests.csv"),
    )
    parser.add_argument(
        "--composition-adjusted-tests",
        default=str(BASE_RESULT_DIR / "tissue_composition_sensitivity/composition_adjusted_channel_tests.csv"),
    )
    parser.add_argument(
        "--tumor-proxy-pairwise",
        default=str(BASE_RESULT_DIR / "tumor_proxy_sensitivity/tumor_proxy_low_zero_pairwise_tests.csv"),
    )
    parser.add_argument(
        "--tumor-proxy-classifier",
        default=str(BASE_RESULT_DIR / "tumor_proxy_sensitivity/tumor_proxy_classifier_best_h_e_metrics.csv"),
    )
    parser.add_argument(
        "--tumor-proxy-retention",
        default=str(BASE_RESULT_DIR / "tumor_proxy_sensitivity/tumor_proxy_retention_summary.csv"),
    )
    parser.add_argument(
        "--classifier-permutation-summary",
        default=str(BASE_RESULT_DIR / "classifier_permutation_sanity/classifier_permutation_summary.csv"),
    )
    parser.add_argument(
        "--nested-classifier-summary",
        default=str(BASE_RESULT_DIR / "nested_classifier_model_selection/nested_classifier_summary.csv"),
    )
    parser.add_argument(
        "--covariate-numeric",
        default=str(BASE_RESULT_DIR / "clinical_covariate_sensitivity/covariate_balance_numeric.csv"),
    )
    parser.add_argument(
        "--covariate-categorical",
        default=str(BASE_RESULT_DIR / "clinical_covariate_sensitivity/covariate_balance_categorical.csv"),
    )
    parser.add_argument(
        "--covariate-classifier",
        default=str(BASE_RESULT_DIR / "clinical_covariate_sensitivity/covariate_classifier_metrics.csv"),
    )
    parser.add_argument(
        "--covariate-adjusted",
        default=str(BASE_RESULT_DIR / "clinical_covariate_sensitivity/covariate_adjusted_channel_tests.csv"),
    )
    parser.add_argument(
        "--matched-pair-summary",
        default=str(BASE_RESULT_DIR / "matched_low_zero_sensitivity/matched_pair_summary.csv"),
    )
    parser.add_argument(
        "--matched-classifier",
        default=str(BASE_RESULT_DIR / "matched_low_zero_sensitivity/matched_classifier_metrics.csv"),
    )
    parser.add_argument(
        "--matched-channel-tests",
        default=str(BASE_RESULT_DIR / "matched_low_zero_sensitivity/matched_channel_tests.csv"),
    )
    parser.add_argument(
        "--isoform-feasibility-summary",
        default=str(BASE_RESULT_DIR / "her2_isoform_validation_feasibility/her2_isoform_validation_feasibility_summary.json"),
    )
    parser.add_argument(
        "--isoform-feasibility-table",
        default=str(BASE_RESULT_DIR / "her2_isoform_validation_feasibility/her2_isoform_validation_feasibility_table.csv"),
    )
    parser.add_argument(
        "--source-site-generalization",
        default=str(BASE_RESULT_DIR / "source_site_generalization/source_site_generalization_metrics.csv"),
    )
    parser.add_argument(
        "--source-site-balance",
        default=str(BASE_RESULT_DIR / "source_site_generalization/source_site_balance.csv"),
    )
    parser.add_argument(
        "--within-source-site-summary",
        default=str(BASE_RESULT_DIR / "within_source_site_low_zero/within_source_site_summary.json"),
    )
    parser.add_argument(
        "--within-source-site-balance",
        default=str(BASE_RESULT_DIR / "within_source_site_low_zero/mixed_source_site_balance.csv"),
    )
    parser.add_argument(
        "--within-source-site-channel-tests",
        default=str(BASE_RESULT_DIR / "within_source_site_low_zero/within_source_site_channel_tests.csv"),
    )
    parser.add_argument(
        "--within-source-site-classifier",
        default=str(BASE_RESULT_DIR / "within_source_site_low_zero/within_source_site_classifier_metrics.csv"),
    )
    parser.add_argument(
        "--local-erbb2-summary",
        default=str(BASE_RESULT_DIR / "local_erbb2_expression_validation/local_erbb2_validation_summary.json"),
    )
    parser.add_argument(
        "--local-erbb2-group-summary",
        default=str(BASE_RESULT_DIR / "local_erbb2_expression_validation/local_erbb2_group_summary.csv"),
    )
    parser.add_argument(
        "--local-erbb2-pairwise",
        default=str(BASE_RESULT_DIR / "local_erbb2_expression_validation/local_erbb2_pairwise_tests.csv"),
    )
    parser.add_argument(
        "--local-erbb2-reference-classifier",
        default=str(BASE_RESULT_DIR / "local_erbb2_expression_validation/local_erbb2_reference_classifier_metrics.csv"),
    )
    parser.add_argument(
        "--local-erbb2-correlations",
        default=str(BASE_RESULT_DIR / "local_erbb2_expression_validation/local_erbb2_gigatime_correlations.csv"),
    )
    parser.add_argument(
        "--local-erbb2-adjusted-tests",
        default=str(BASE_RESULT_DIR / "local_erbb2_expression_validation/local_erbb2_adjusted_low_zero_channel_tests.csv"),
    )
    return parser.parse_args()


def read_rows(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_json(path: str | Path) -> dict[str, object]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def as_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def fmt(value: float | str, digits: int = 3) -> str:
    if isinstance(value, str):
        return value
    if value != value:
        return ""
    if abs(value) < 0.001 and value != 0:
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def pct(value: float) -> str:
    if value != value:
        return ""
    return f"{100 * value:.1f}%"


def feature_set_label(value: str) -> str:
    return FEATURE_SET_LABELS.get(value, value.replace("_", " "))


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


def html_table(headers: list[str], rows: list[list[str]]) -> str:
    header_html = "".join(f"<th>{html.escape(str(header))}</th>" for header in headers)
    if rows:
        body = "".join(
            "<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>"
            for row in rows
        )
    else:
        body = "<tr>" + "".join("<td></td>" for _ in headers) + "</tr>"
    return f"<table><thead><tr>{header_html}</tr></thead><tbody>{body}</tbody></table>"


def copy_assets() -> None:
    for source, destination in ASSET_COPIES:
        if not source.exists():
            raise FileNotFoundError(source)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def view_sort_key(row: dict[str, str]) -> tuple[int, float, str]:
    return (
        FEATURE_VIEW_ORDER.get(row.get("feature_view", ""), 99),
        as_float(row, "mannwhitney_q_value_bh_within_view"),
        row.get("channel", ""),
    )


def group_count_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    counts: dict[str, int] = {}
    for row in rows:
        group = row.get("clinical_her2_group", "")
        counts[group] = counts.get(group, 0) + 1
    order = ["HER2-positive", "HER2-low", "HER2-zero"]
    return [[group, str(counts.get(group, 0))] for group in order]


def low_zero_pairwise_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = [
        row
        for row in rows
        if row.get("group_a") == "HER2-low"
        and row.get("group_b") == "HER2-zero"
        and row.get("channel") in KEY_CHANNELS
        and as_float(row, "mannwhitney_q_value_bh") <= 0.05
    ]
    selected.sort(key=lambda row: as_float(row, "mannwhitney_q_value_bh"))
    return [
        [
            row["channel"],
            fmt(as_float(row, "mean_a"), 4),
            fmt(as_float(row, "mean_b"), 4),
            fmt(as_float(row, "delta_mean_a_minus_b"), 4),
            fmt(as_float(row, "mannwhitney_p_value"), 4),
            fmt(as_float(row, "mannwhitney_q_value_bh"), 4),
            fmt(as_float(row, "cliffs_delta"), 3),
        ]
        for row in selected
    ]


def three_group_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = [row for row in rows if row.get("channel") in KEY_CHANNELS]
    selected.sort(key=lambda row: as_float(row, "kruskal_q_value_bh"))
    return [
        [
            row["channel"],
            fmt(as_float(row, "kruskal_p_value"), 4),
            fmt(as_float(row, "kruskal_q_value_bh"), 4),
            row["highest_mean_group"],
            row["lowest_mean_group"],
            fmt(as_float(row, "max_minus_min_mean"), 4),
        ]
        for row in selected[:8]
    ]


def cleanup_low_zero_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = [
        row
        for row in rows
        if row.get("group_a") == "HER2-low"
        and row.get("group_b") == "HER2-zero"
        and row.get("channel") in {"CD68", "PD-L1", "CD11c", "CK", "CD4", "CD3"}
    ]
    selected.sort(key=view_sort_key)
    return [
        [
            row["feature_view_label"],
            row["channel"],
            fmt(as_float(row, "delta_mean_a_minus_b"), 4),
            fmt(as_float(row, "mannwhitney_p_value"), 4),
            fmt(as_float(row, "mannwhitney_q_value_bh_within_view"), 4),
        ]
        for row in selected
    ]


def cleanup_three_group_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = [
        row
        for row in rows
        if row.get("channel") in {"CD68", "PD-L1", "CD11c", "CK", "CD4", "CD3"}
        and as_float(row, "kruskal_q_value_bh_within_view") <= 0.05
    ]
    selected.sort(
        key=lambda row: (
            FEATURE_VIEW_ORDER.get(row.get("feature_view", ""), 99),
            as_float(row, "kruskal_q_value_bh_within_view"),
            row.get("channel", ""),
        )
    )
    return [
        [
            row["feature_view_label"],
            row["channel"],
            fmt(as_float(row, "kruskal_p_value"), 4),
            fmt(as_float(row, "kruskal_q_value_bh_within_view"), 4),
            row["highest_mean_group"],
            row["lowest_mean_group"],
            fmt(as_float(row, "max_minus_min_mean"), 4),
        ]
        for row in selected[:18]
    ]


def classifier_low_zero_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = [
        row
        for row in rows
        if row.get("task") == "her2_low_vs_zero" and row.get("model") == "regularized_logistic"
    ]
    selected.sort(key=lambda row: FEATURE_VIEW_ORDER.get(row.get("feature_view", ""), 99))
    return [
        [
            row["feature_view_label"],
            feature_set_label(row["feature_set"]),
            row["n_cases"],
            pct(as_float(row, "accuracy")),
            pct(as_float(row, "balanced_accuracy")),
            fmt(as_float(row, "macro_auc_ovr"), 3),
            pct(as_float(row, "sensitivity")),
            pct(as_float(row, "specificity")),
        ]
        for row in selected
    ]


def classifier_other_task_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = [
        row
        for row in rows
        if row.get("task") in {"her2_positive_vs_negative", "her2_three_class"}
        and row.get("model") == "regularized_logistic"
    ]
    best_by_task: dict[str, dict[str, str]] = {}
    for row in selected:
        task = row["task"]
        if task not in best_by_task or as_float(row, "balanced_accuracy") > as_float(
            best_by_task[task], "balanced_accuracy"
        ):
            best_by_task[task] = row
    ordered = [best_by_task[key] for key in ("her2_positive_vs_negative", "her2_three_class") if key in best_by_task]
    return [
        [
            row["task_label"],
            row["feature_view_label"],
            feature_set_label(row["feature_set"]),
            row["n_cases"],
            pct(as_float(row, "balanced_accuracy")),
            fmt(as_float(row, "macro_auc_ovr"), 3),
            pct(as_float(row, "sensitivity")),
            pct(as_float(row, "specificity")),
        ]
        for row in ordered
    ]


def erpr_adjusted_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = [
        row
        for row in rows
        if row.get("feature_view") == "all_sampled_tissue" and row.get("channel") in KEY_CHANNELS
    ]
    selected.sort(key=lambda row: as_float(row, "erpr_adjusted_q_value_bh_within_view"))
    return [
        [
            row["channel"],
            fmt(as_float(row, "delta_mean_a_minus_b"), 4),
            fmt(as_float(row, "mannwhitney_q_value_bh_within_view"), 4),
            fmt(as_float(row, "erpr_adjusted_beta_low_vs_zero"), 4),
            fmt(as_float(row, "erpr_adjusted_q_value_bh_within_view"), 4),
            fmt(as_float(row, "erpr_erbb2_adjusted_q_value_bh_within_view"), 4),
            row.get("erpr_erbb2_adjusted_n", ""),
        ]
        for row in selected
    ]


def case_driver_signal_rows(summary: dict[str, object]) -> list[list[str]]:
    rows = []
    order = {view: idx for idx, view in enumerate(FEATURE_VIEW_ORDER)}
    for row in sorted(summary.get("signal_channels_by_view", []), key=lambda item: order.get(item["feature_view"], 99)):
        rows.append([row["feature_view_label"], str(row["n_signal_channels"]), row["channels"]])
    return rows


def case_driver_group_score_rows(summary: dict[str, object]) -> list[list[str]]:
    return [
        [
            row["clinical_her2_group"],
            str(int(row["count"])),
            fmt(row["mean"], 3),
            fmt(row["median"], 3),
        ]
        for row in summary.get("all_sampled_tissue_zero_like_score_by_group", [])
    ]


def case_driver_stability_rows(summary: dict[str, object]) -> list[list[str]]:
    counts = {
        (row["clinical_her2_group"], int(row["expected_profile_views"])): int(row["n_slides"])
        for row in summary.get("stability_counts", [])
    }
    return [
        [group] + [str(counts.get((group, views), 0)) for views in range(5)]
        for group in ["HER2-low", "HER2-zero"]
    ]


def case_driver_review_rows(rows: list[dict[str, str]], limit: int = 8) -> list[list[str]]:
    selected = rows[:limit]
    return [
        [
            row["case_submitter_id"],
            row["clinical_her2_group"],
            row.get("her2_detail_subgroup", ""),
            row.get("classifier_incorrect_views", ""),
            row.get("opposite_profile_views", ""),
            fmt(as_float(row, "all_sampled_tissue_zero_like_score"), 3),
        ]
        for row in selected
    ]


def visual_qc_tile_summary_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        key = (row.get("review_category", ""), row.get("clinical_her2_group", ""))
        groups.setdefault(key, []).append(row)

    def mean(group_rows: list[dict[str, str]], key: str) -> float:
        values = [as_float(row, key) for row in group_rows]
        values = [value for value in values if value == value]
        return sum(values) / len(values) if values else float("nan")

    order = [
        ("label_consistent_her2_low", "HER2-low"),
        ("label_consistent_her2_zero", "HER2-zero"),
        ("opposite_profile_manual_review", "HER2-low"),
        ("opposite_profile_manual_review", "HER2-zero"),
    ]
    rows_out = []
    for key in order:
        group_rows = groups.get(key, [])
        rows_out.append(
            [
                key[0].replace("_", " "),
                key[1],
                str(len(group_rows)),
                fmt(mean(group_rows, "tissue_fraction"), 3),
                fmt(mean(group_rows, "tile_zero_like_score"), 3),
                fmt(mean(group_rows, "mean_CK"), 4),
                fmt(mean(group_rows, "mean_CD68"), 4),
                fmt(mean(group_rows, "mean_PD-L1"), 4),
                fmt(mean(group_rows, "mean_CD11c"), 4),
            ]
        )
    return rows_out


def tissue_composition_rows(rows: list[dict[str, str]], limit: int = 8) -> list[list[str]]:
    label_map = {
        "fraction_low_marker_q25": "Fraction low-marker tiles",
        "fraction_very_low_marker_q10": "Fraction very-low-marker tiles",
        "mean_marker_burden": "Mean marker burden",
        "mean_CK": "Mean CK",
        "fraction_high_marker_q75": "Fraction high-marker tiles",
        "fraction_high_ck_q75": "Fraction high-CK QC tiles",
        "fraction_low_ck_q25": "Fraction low-CK QC tiles",
        "mean_DAPI": "Mean DAPI",
    }
    selected = rows[:limit]
    return [
        [
            label_map.get(row["metric"], row["metric"]),
            fmt(as_float(row, "mean_low"), 4),
            fmt(as_float(row, "mean_zero"), 4),
            fmt(as_float(row, "delta_low_minus_zero"), 4),
            fmt(as_float(row, "mannwhitney_q_value_bh"), 4),
            fmt(as_float(row, "cliffs_delta"), 3),
        ]
        for row in selected
    ]


def tissue_correlation_rows(rows: list[dict[str, str]], limit: int = 6) -> list[list[str]]:
    return [
        [
            row["metric"],
            row["n"],
            fmt(as_float(row, "spearman_rho"), 3),
            fmt(as_float(row, "q_value_bh"), 4),
        ]
        for row in rows[:limit]
    ]


def absolute_ck_rows(rows: list[dict[str, str]], limit: int = 8) -> list[list[str]]:
    selected = rows[:limit]
    return [
        [
            row["channel"],
            row["n_low"],
            row["n_zero"],
            fmt(as_float(row, "delta_low_minus_zero"), 4),
            fmt(as_float(row, "mannwhitney_q_value_bh"), 4),
        ]
        for row in selected
    ]


def composition_adjusted_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    by_channel: dict[str, dict[str, dict[str, str]]] = {}
    for row in rows:
        by_channel.setdefault(row["channel"], {})[row["model"]] = row
    selected_channels = ["CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CK"]
    output = []
    for channel in selected_channels:
        models = by_channel.get(channel, {})
        unadjusted = models.get("unadjusted", {})
        low_marker = models.get("adjusted_low_marker_fraction", {})
        mean_ck = models.get("adjusted_mean_ck", {})
        output.append(
            [
                channel,
                fmt(as_float(unadjusted, "beta_low_vs_zero"), 4),
                fmt(as_float(unadjusted, "q_value_bh_within_model"), 4),
                fmt(as_float(low_marker, "beta_low_vs_zero"), 4),
                fmt(as_float(low_marker, "q_value_bh_within_model"), 4),
                fmt(as_float(mean_ck, "beta_low_vs_zero"), 4),
                fmt(as_float(mean_ck, "q_value_bh_within_model"), 4),
            ]
        )
    return output


def tumor_proxy_retention_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = [
        row
        for row in rows
        if row.get("feature_view")
        in {
            "ck_top16_within_slide",
            "ck_top8_within_slide",
            "ck_top16_non_low_marker",
            "absolute_ck_high_q75",
        }
    ]
    selected.sort(
        key=lambda row: (
            TUMOR_PROXY_VIEW_ORDER.get(row.get("feature_view", ""), 99),
            {"HER2-positive": 0, "HER2-low": 1, "HER2-zero": 2}.get(row.get("clinical_her2_group", ""), 99),
        )
    )
    return [
        [
            row["feature_view_label"],
            row["clinical_her2_group"],
            row["n_slides_passing_min_tiles"],
            fmt(as_float(row, "median_tiles_after_min_filter"), 1),
            fmt(as_float(row, "median_retained_fraction"), 3),
        ]
        for row in selected
    ]


def tumor_proxy_pairwise_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    focus_channels = {"CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CK"}
    selected = [
        row
        for row in rows
        if row.get("feature_view")
        in {"ck_top16_within_slide", "ck_top8_within_slide", "ck_top16_non_low_marker", "absolute_ck_high_q75"}
        and row.get("channel") in focus_channels
    ]
    selected.sort(
        key=lambda row: (
            TUMOR_PROXY_VIEW_ORDER.get(row.get("feature_view", ""), 99),
            as_float(row, "mannwhitney_q_value_bh_within_view"),
            row.get("channel", ""),
        )
    )
    return [
        [
            row["feature_view_label"],
            row["channel"],
            row["n_low"],
            row["n_zero"],
            fmt(as_float(row, "delta_low_minus_zero"), 4),
            fmt(as_float(row, "mannwhitney_q_value_bh_within_view"), 4),
        ]
        for row in selected
    ]


def tumor_proxy_classifier_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = [
        row
        for row in rows
        if row.get("task") == "her2_low_vs_zero" and row.get("model") == "regularized_logistic"
    ]
    selected.sort(key=lambda row: TUMOR_PROXY_VIEW_ORDER.get(row.get("feature_view", ""), 99))
    return [
        [
            row["feature_view_label"],
            feature_set_label(row["feature_set"]),
            row["n_cases"],
            pct(as_float(row, "accuracy")),
            pct(as_float(row, "balanced_accuracy")),
            fmt(as_float(row, "macro_auc_ovr"), 3),
            pct(as_float(row, "sensitivity")),
            pct(as_float(row, "specificity")),
        ]
        for row in selected
    ]


def classifier_permutation_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = sorted(
        rows,
        key=lambda row: TUMOR_PROXY_VIEW_ORDER.get(row.get("feature_view", ""), 99),
    )
    return [
        [
            row["feature_view_label"],
            row["feature_set_label"],
            row["n_cases"],
            row["n_features"],
            pct(as_float(row, "loocv_balanced_accuracy")),
            pct(as_float(row, "observed_repeated_cv_balanced_accuracy")),
            pct(as_float(row, "null_balanced_accuracy_mean")),
            pct(as_float(row, "null_balanced_accuracy_p95")),
            fmt(as_float(row, "empirical_p_balanced_accuracy"), 4),
            fmt(as_float(row, "empirical_q_balanced_accuracy_bh"), 4),
            fmt(as_float(row, "observed_repeated_cv_macro_auc_ovr"), 3),
        ]
        for row in selected
    ]


def nested_classifier_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = sorted(
        rows,
        key=lambda row: TUMOR_PROXY_VIEW_ORDER.get(row.get("feature_view", ""), 99),
    )
    return [
        [
            row["feature_view_label"],
            row["n_cases"],
            row["most_selected_feature_set_label"],
            row["most_selected_count"],
            pct(as_float(row, "observed_nested_balanced_accuracy")),
            fmt(as_float(row, "observed_nested_macro_auc_ovr"), 3),
            pct(as_float(row, "null_balanced_accuracy_mean")),
            pct(as_float(row, "null_balanced_accuracy_p95")),
            fmt(as_float(row, "empirical_p_balanced_accuracy"), 4),
            fmt(as_float(row, "empirical_q_balanced_accuracy_bh"), 4),
        ]
        for row in selected
    ]


def covariate_balance_rows(numeric_rows: list[dict[str, str]], categorical_rows: list[dict[str, str]]) -> list[list[str]]:
    rows = []
    numeric_focus = {"slide_file_size_mb", "slide_width", "slide_height", "mean_marker_burden", "mean_DAPI", "mean_CK"}
    for row in numeric_rows:
        if row.get("covariate") in numeric_focus:
            rows.append(
                [
                    row["covariate_label"],
                    fmt(as_float(row, "mean_low"), 3),
                    fmt(as_float(row, "mean_zero"), 3),
                    fmt(as_float(row, "delta_low_minus_zero"), 3),
                    fmt(as_float(row, "mannwhitney_p_value"), 4),
                ]
            )
    for row in categorical_rows:
        if row.get("covariate") in {"histology_group", "tss_group"}:
            rows.append(
                [
                    row["covariate_label"],
                    "",
                    "",
                    "imbalanced",
                    fmt(as_float(row, "chi_square_p_value"), 4),
                ]
            )
    return rows


def covariate_classifier_rows(rows: list[dict[str, str]], view: str = "ck_top8_within_slide") -> list[list[str]]:
    focus_sets = [
        "clinical_basic",
        "slide_size_only",
        "source_site_only",
        "site_slide_only",
        "gigatime_mean_channels",
        "gigatime_plus_clinical_site_slide",
    ]
    selected = [
        row
        for row in rows
        if row.get("feature_view") == view and row.get("feature_set") in set(focus_sets)
    ]
    selected.sort(key=lambda row: focus_sets.index(row["feature_set"]))
    return [
        [
            row["feature_set_label"],
            row["n_features"],
            pct(as_float(row, "balanced_accuracy")),
            fmt(as_float(row, "macro_auc_ovr"), 3),
        ]
        for row in selected
    ]


def covariate_adjusted_summary_rows(rows: list[dict[str, str]], view: str = "qc_cellular_tissue") -> list[list[str]]:
    selected = [row for row in rows if row.get("feature_view") == view]
    order = [
        "Unadjusted",
        "ER/PR adjusted",
        "ER/PR + histology + stage",
        "Clinical + slide size",
        "Clinical + site/slide size",
    ]
    output = []
    for label in order:
        group = [row for row in selected if row.get("model_label") == label]
        if not group:
            continue
        significant = [row["channel"] for row in group if as_float(row, "q_value_bh_within_model_view") < 0.05]
        best_q = min(as_float(row, "q_value_bh_within_model_view") for row in group)
        output.append([label, str(len(significant)), ", ".join(significant) if significant else "none", fmt(best_q, 4)])
    return output


def matched_pair_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    return [
        [
            row["matched_subset_label"],
            row["n_pairs"],
            row["n_same_tss_pairs"],
            fmt(as_float(row, "median_abs_log_slide_size_diff"), 3),
            fmt(as_float(row, "median_abs_slide_file_size_mb_diff"), 1),
        ]
        for row in rows
    ]


def matched_classifier_rows(rows: list[dict[str, str]], view: str = "ck_top8_within_slide") -> list[list[str]]:
    focus_sets = [
        "slide_size_only",
        "source_site_only",
        "source_site_slide_size",
        "gigatime_mean_channels",
        "gigatime_plus_source_site_slide_size",
    ]
    selected = [
        row
        for row in rows
        if row.get("feature_view") == view and row.get("feature_set") in set(focus_sets)
    ]
    subset_order = {
        "exact_source_site_nearest_size": 0,
        "slide_size_caliper_0.25": 1,
        "slide_size_caliper_0.50": 2,
    }
    selected.sort(
        key=lambda row: (
            subset_order.get(row.get("matched_subset", ""), 99),
            focus_sets.index(row["feature_set"]),
        )
    )
    return [
        [
            row["matched_subset_label"],
            row["feature_set_label"],
            row["n_pairs"],
            row["n_features"],
            pct(as_float(row, "balanced_accuracy")),
            fmt(as_float(row, "macro_auc_ovr"), 3),
        ]
        for row in selected
    ]


def matched_channel_rows(rows: list[dict[str, str]], view: str = "ck_top8_within_slide") -> list[list[str]]:
    focus_channels = {"CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CD20", "CK", "Ki67"}
    selected = [
        row
        for row in rows
        if row.get("feature_view") == view and row.get("channel") in focus_channels
    ]
    subset_order = {
        "exact_source_site_nearest_size": 0,
        "slide_size_caliper_0.25": 1,
        "slide_size_caliper_0.50": 2,
    }
    selected.sort(
        key=lambda row: (
            subset_order.get(row.get("matched_subset", ""), 99),
            as_float(row, "q_value_bh_within_subset_view"),
            row.get("channel", ""),
        )
    )
    output = []
    seen_counts: dict[str, int] = {}
    for row in selected:
        subset = row["matched_subset"]
        seen_counts[subset] = seen_counts.get(subset, 0) + 1
        if seen_counts[subset] > 4:
            continue
        output.append(
            [
                row["matched_subset_label"],
                row["channel"],
                row["n_pairs"],
                fmt(as_float(row, "mean_low_minus_zero"), 4),
                fmt(as_float(row, "wilcoxon_p_value"), 4),
                fmt(as_float(row, "q_value_bh_within_subset_view"), 4),
            ]
        )
    return output


def isoform_audit_rows(summary: dict[str, object]) -> list[list[str]]:
    schema = summary.get("expression_schema", {})
    return [
        ["High-trust cases", summary.get("n_high_trust_cases", "")],
        ["Local STAR gene-count cases", summary.get("n_local_star_gene_count_cases", "")],
        ["High-trust cases with local STAR gene counts", summary.get("n_high_trust_cases_with_local_star_gene_counts", "")],
        ["Low/zero high-trust cases with local STAR gene counts", summary.get("n_low_zero_high_trust_cases_with_local_star_gene_counts", "")],
        ["Local BAM files under data/tcga_brca", summary.get("local_bam_files", "")],
        ["Local FASTQ files under data/tcga_brca", summary.get("local_fastq_files", "")],
        ["Local junction files under data/tcga_brca", summary.get("local_junction_files", "")],
        ["Expression file has transcript_id", schema.get("has_transcript_id", "")],
        ["Expression file has junction columns", schema.get("has_junction_columns", "")],
    ]


def isoform_feasibility_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    return [
        [
            row["analysis"],
            row["status"],
            row["what_it_can_support"],
            row["next_action"],
        ]
        for row in rows
    ]


def local_erbb2_group_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    order = {"HER2-positive": 0, "HER2-low": 1, "HER2-zero": 2}
    selected = sorted(rows, key=lambda row: order.get(row.get("clinical_her2_group", ""), 99))
    return [
        [
            row["clinical_her2_group"],
            row["n_with_local_erbb2"],
            fmt(as_float(row, "median_erbb2_tpm")),
            fmt(as_float(row, "q25_erbb2_tpm")),
            fmt(as_float(row, "q75_erbb2_tpm")),
        ]
        for row in selected
    ]


def local_erbb2_pairwise_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    return [
        [
            row["comparison"],
            row["n_a"],
            row["n_b"],
            fmt(as_float(row, "median_a_erbb2_tpm")),
            fmt(as_float(row, "median_b_erbb2_tpm")),
            fmt(as_float(row, "separation_auc_abs_direction"), 3),
            fmt(as_float(row, "mannwhitney_p_value"), 3),
            fmt(as_float(row, "mannwhitney_bh_q_value"), 3),
        ]
        for row in rows
    ]


def local_erbb2_reference_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    return [
        [
            row["task"],
            row["n_cases"],
            fmt(as_float(row, "auc"), 3),
            pct(as_float(row, "best_threshold_balanced_accuracy")),
        ]
        for row in rows
    ]


def local_erbb2_correlation_rows(rows: list[dict[str, str]], limit: int = 8) -> list[list[str]]:
    selected = [row for row in rows if row.get("subset") == "HER2-low/HER2-zero"]
    selected.sort(key=lambda row: abs(as_float(row, "spearman_rho_log_erbb2")), reverse=True)
    return [
        [
            row["feature_view_label"],
            row["channel"],
            row["n_cases"],
            fmt(as_float(row, "spearman_rho_log_erbb2"), 3),
            fmt(as_float(row, "p_value"), 3),
            fmt(as_float(row, "bh_q_value"), 3),
        ]
        for row in selected[:limit]
    ]


def local_erbb2_adjusted_rows(rows: list[dict[str, str]], limit: int = 8) -> list[list[str]]:
    selected = sorted(rows, key=lambda row: as_float(row, "bh_q_value"))
    return [
        [
            row["feature_view_label"],
            row["channel"],
            row["n_cases"],
            fmt(as_float(row, "beta_her2_zero_vs_low_adjusted_log_erbb2"), 4),
            fmt(as_float(row, "p_value"), 3),
            fmt(as_float(row, "bh_q_value"), 3),
        ]
        for row in selected[:limit]
    ]


def source_site_balance_rows(rows: list[dict[str, str]], limit: int = 8) -> list[list[str]]:
    return [
        [
            row["tss_code"],
            row.get("n_HER2-low", ""),
            row.get("n_HER2-zero", ""),
            row.get("n_cases", ""),
            "yes" if row.get("has_both_classes") == "1" else "no",
        ]
        for row in rows[:limit]
    ]


def source_site_generalization_rows(rows: list[dict[str, str]], view: str = "ck_top8_within_slide") -> list[list[str]]:
    focus_sets = ["slide_size_only", "tissue_qc_only", "gigatime_mean_channels", "gigatime_plus_slide_size"]
    scheme_order = {"repeated_stratified_cv": 0, "leave_source_site_out": 1}
    selected = [
        row
        for row in rows
        if row.get("feature_view") == view and row.get("feature_set") in set(focus_sets)
    ]
    selected.sort(key=lambda row: (focus_sets.index(row["feature_set"]), scheme_order.get(row["validation_scheme"], 99)))
    return [
        [
            row["feature_set_label"],
            row["validation_scheme_label"],
            row["n_features"],
            pct(as_float(row, "balanced_accuracy")),
            fmt(as_float(row, "macro_auc_ovr"), 3),
        ]
        for row in selected
    ]


def source_site_gigatime_drop_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    by_key = {
        (row["feature_view"], row["validation_scheme"]): row
        for row in rows
        if row.get("feature_set") == "gigatime_mean_channels"
    }
    order = [
        "qc_cellular_tissue",
        "ck_top25_within_slide",
        "ck_top16_within_slide",
        "ck_top8_within_slide",
        "ck_top16_non_low_marker",
        "absolute_ck_high_q75",
    ]
    output = []
    for view in order:
        random_row = by_key.get((view, "repeated_stratified_cv"))
        site_row = by_key.get((view, "leave_source_site_out"))
        if not random_row or not site_row:
            continue
        random_ba = as_float(random_row, "balanced_accuracy")
        site_ba = as_float(site_row, "balanced_accuracy")
        output.append(
            [
                random_row["feature_view_label"],
                pct(random_ba),
                pct(site_ba),
                pct(site_ba - random_ba),
            ]
        )
    return output


def erpr_summary_rows(summary: dict[str, object]) -> list[list[str]]:
    counts = summary.get("adjusted_q_lt_0_05_by_view", {})
    labels = {
        "all_sampled_tissue": "All sampled tissue",
        "qc_cellular_tissue": "QC cellular tissue",
        "ck_enriched_top50": "CK-enriched top 50%",
        "ck_enriched_top25": "CK-enriched top 25%",
    }
    rows = []
    for view, label in labels.items():
        view_counts = counts.get(view, {}) if isinstance(counts, dict) else {}
        rows.append(
            [
                label,
                str(view_counts.get("unadjusted", "")),
                str(view_counts.get("er_pr_adjusted", "")),
                str(view_counts.get("er_pr_erbb2_adjusted", "")),
            ]
        )
    return rows


def direction_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = [row for row in rows if row.get("channel") in {"CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CK", "Ki67"}]
    selected.sort(key=lambda row: abs(as_float(row, "reference_low_zero_delta")), reverse=True)
    return [
        [
            row["channel"],
            fmt(as_float(row, "reference_low_zero_delta"), 4),
            fmt(as_float(row, "comparison_low_zero_delta"), 4),
            "yes" if row.get("same_direction") == "True" else "no",
            "yes" if row.get("both_low_lower_than_zero") == "True" else "no",
        ]
        for row in selected
    ]


def within_site_balance_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    selected = [row for row in rows if row.get("has_both_classes") == "1"]
    return [
        [
            row["tss_code"],
            row.get("n_HER2-low", ""),
            row.get("n_HER2-zero", ""),
            row.get("n_cases", ""),
        ]
        for row in selected
    ]


def within_site_channel_q_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["feature_view_label"], []).append(row)
    output = []
    for label, group in grouped.items():
        q_values = [as_float(row, "q_value_bh_within_view") for row in group]
        output.append(
            [
                label,
                len(group),
                sum(1 for value in q_values if value < 0.05),
                fmt(min(q_values), 4) if q_values else "",
            ]
        )
    output.sort(key=lambda row: row[0])
    return output


def within_site_top_channel_rows(rows: list[dict[str, str]], view: str = "ck_top8_within_slide") -> list[list[str]]:
    selected = [row for row in rows if row.get("feature_view") == view]
    selected.sort(key=lambda row: (as_float(row, "q_value_bh_within_view"), as_float(row, "p_value")))
    return [
        [
            row["feature_view_label"],
            row["channel"],
            row["n_low"],
            row["n_zero"],
            fmt(as_float(row, "beta_zero_vs_low_site_fixed"), 4),
            fmt(as_float(row, "p_value"), 4),
            fmt(as_float(row, "q_value_bh_within_view"), 4),
        ]
        for row in selected[:8]
    ]


def within_site_classifier_rows(rows: list[dict[str, str]], view: str = "ck_top8_within_slide") -> list[list[str]]:
    focus_sets = [
        "slide_size_only",
        "tissue_qc_only",
        "source_site_one_hot",
        "gigatime_key_mean_channels",
        "gigatime_mean_channels",
    ]
    order = {feature: idx for idx, feature in enumerate(focus_sets)}
    selected = [
        row
        for row in rows
        if row.get("feature_view") == view and row.get("feature_set") in focus_sets
    ]
    selected.sort(key=lambda row: (row.get("validation_scheme_label", ""), order.get(row.get("feature_set", ""), 99)))
    return [
        [
            row["feature_set_label"],
            row["validation_scheme_label"],
            row["n_cases"],
            row["n_features"],
            pct(as_float(row, "balanced_accuracy")),
            fmt(as_float(row, "macro_auc_ovr"), 3),
            fmt(as_float(row, "sensitivity"), 3),
            fmt(as_float(row, "specificity"), 3),
        ]
        for row in selected
    ]


def build_content(args: argparse.Namespace) -> dict[str, object]:
    copy_assets()
    high_trust_slides = read_rows(args.high_trust_slides)
    trust_summary = read_json(args.trust_summary)
    sensitivity_summary = read_json(args.sensitivity_summary)
    agreement_summary = read_json(args.run_agreement_summary)
    clinical_pairwise = read_rows(args.clinical_pairwise)
    clinical_channels = read_rows(args.clinical_channels)
    cleanup_pairwise = read_rows(args.cleanup_pairwise)
    cleanup_channels = read_rows(args.cleanup_channels)
    classifier_metrics = read_rows(args.classifier_metrics)
    erpr_adjusted = read_rows(args.erpr_adjusted)
    direction_comparison = read_rows(args.direction_comparison)
    case_driver_summary = read_json(args.case_driver_summary)
    case_driver_review = read_rows(args.case_driver_review)
    visual_qc_selected_tiles = read_rows(args.visual_qc_selected_tiles)
    tissue_composition_tests = read_rows(args.tissue_composition_tests)
    tissue_composition_correlations = read_rows(args.tissue_composition_correlations)
    absolute_ck_pairwise = read_rows(args.absolute_ck_pairwise)
    composition_adjusted_tests = read_rows(args.composition_adjusted_tests)
    tumor_proxy_pairwise = read_rows(args.tumor_proxy_pairwise)
    tumor_proxy_classifier = read_rows(args.tumor_proxy_classifier)
    tumor_proxy_retention = read_rows(args.tumor_proxy_retention)
    classifier_permutation = read_rows(args.classifier_permutation_summary)
    nested_classifier = read_rows(args.nested_classifier_summary)
    covariate_numeric = read_rows(args.covariate_numeric)
    covariate_categorical = read_rows(args.covariate_categorical)
    covariate_classifier = read_rows(args.covariate_classifier)
    covariate_adjusted = read_rows(args.covariate_adjusted)
    matched_pair_summary = read_rows(args.matched_pair_summary)
    matched_classifier = read_rows(args.matched_classifier)
    matched_channel_tests = read_rows(args.matched_channel_tests)
    isoform_feasibility_summary = read_json(args.isoform_feasibility_summary)
    isoform_feasibility_table = read_rows(args.isoform_feasibility_table)
    source_site_generalization = read_rows(args.source_site_generalization)
    source_site_balance = read_rows(args.source_site_balance)
    within_source_site_summary = read_json(args.within_source_site_summary)
    within_source_site_balance = read_rows(args.within_source_site_balance)
    within_source_site_channel_tests = read_rows(args.within_source_site_channel_tests)
    within_source_site_classifier = read_rows(args.within_source_site_classifier)
    local_erbb2_summary = read_json(args.local_erbb2_summary)
    local_erbb2_group_summary = read_rows(args.local_erbb2_group_summary)
    local_erbb2_pairwise = read_rows(args.local_erbb2_pairwise)
    local_erbb2_reference_classifier = read_rows(args.local_erbb2_reference_classifier)
    local_erbb2_correlations = read_rows(args.local_erbb2_correlations)
    local_erbb2_adjusted_tests = read_rows(args.local_erbb2_adjusted_tests)

    return {
        "trust_summary": trust_summary,
        "sensitivity_summary": sensitivity_summary,
        "agreement_summary": agreement_summary,
        "case_driver_summary": case_driver_summary,
        "group_counts": group_count_rows(high_trust_slides),
        "low_zero_pairwise": low_zero_pairwise_rows(clinical_pairwise),
        "three_group": three_group_rows(clinical_channels),
        "cleanup_low_zero": cleanup_low_zero_rows(cleanup_pairwise),
        "cleanup_three_group": cleanup_three_group_rows(cleanup_channels),
        "classifier_low_zero": classifier_low_zero_rows(classifier_metrics),
        "classifier_other_tasks": classifier_other_task_rows(classifier_metrics),
        "erpr_summary": erpr_summary_rows(sensitivity_summary),
        "erpr_adjusted": erpr_adjusted_rows(erpr_adjusted),
        "direction_comparison": direction_rows(direction_comparison),
        "case_driver_signal": case_driver_signal_rows(case_driver_summary),
        "case_driver_group_scores": case_driver_group_score_rows(case_driver_summary),
        "case_driver_stability": case_driver_stability_rows(case_driver_summary),
        "case_driver_review": case_driver_review_rows(case_driver_review),
        "visual_qc_tiles": visual_qc_tile_summary_rows(visual_qc_selected_tiles),
        "tissue_composition": tissue_composition_rows(tissue_composition_tests),
        "tissue_correlations": tissue_correlation_rows(tissue_composition_correlations),
        "absolute_ck_pairwise": absolute_ck_rows(absolute_ck_pairwise),
        "composition_adjusted": composition_adjusted_rows(composition_adjusted_tests),
        "tumor_proxy_retention": tumor_proxy_retention_rows(tumor_proxy_retention),
        "tumor_proxy_pairwise": tumor_proxy_pairwise_rows(tumor_proxy_pairwise),
        "tumor_proxy_classifier": tumor_proxy_classifier_rows(tumor_proxy_classifier),
        "classifier_permutation": classifier_permutation_rows(classifier_permutation),
        "nested_classifier": nested_classifier_rows(nested_classifier),
        "covariate_balance": covariate_balance_rows(covariate_numeric, covariate_categorical),
        "covariate_classifier": covariate_classifier_rows(covariate_classifier),
        "covariate_adjusted": covariate_adjusted_summary_rows(covariate_adjusted),
        "matched_pairs": matched_pair_rows(matched_pair_summary),
        "matched_classifier": matched_classifier_rows(matched_classifier),
        "matched_channels": matched_channel_rows(matched_channel_tests),
        "isoform_audit": isoform_audit_rows(isoform_feasibility_summary),
        "isoform_feasibility": isoform_feasibility_rows(isoform_feasibility_table),
        "source_site_balance": source_site_balance_rows(source_site_balance),
        "source_site_generalization": source_site_generalization_rows(source_site_generalization),
        "source_site_gigatime_drop": source_site_gigatime_drop_rows(source_site_generalization),
        "within_source_site_summary": within_source_site_summary,
        "within_source_site_balance": within_site_balance_rows(within_source_site_balance),
        "within_source_site_channel_q": within_site_channel_q_rows(within_source_site_channel_tests),
        "within_source_site_top_channels": within_site_top_channel_rows(within_source_site_channel_tests),
        "within_source_site_classifier": within_site_classifier_rows(within_source_site_classifier),
        "local_erbb2_summary": local_erbb2_summary,
        "local_erbb2_group_summary": local_erbb2_group_rows(local_erbb2_group_summary),
        "local_erbb2_pairwise": local_erbb2_pairwise_rows(local_erbb2_pairwise),
        "local_erbb2_reference_classifier": local_erbb2_reference_rows(local_erbb2_reference_classifier),
        "local_erbb2_correlations": local_erbb2_correlation_rows(local_erbb2_correlations),
        "local_erbb2_adjusted": local_erbb2_adjusted_rows(local_erbb2_adjusted_tests),
    }


def notebook_cell(source: str) -> dict[str, object]:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": [line + "\n" for line in source.strip().splitlines()],
    }


def md_image(path: str, label: str) -> str:
    return f"![{label}]({path})"


def build_markdown_sections(content: dict[str, object]) -> list[tuple[str, str]]:
    trust = content["trust_summary"]
    sensitivity = content["sensitivity_summary"]
    agreement = content["agreement_summary"]
    case_driver = content["case_driver_summary"]
    local_erbb2 = content["local_erbb2_summary"]
    within_site = content["within_source_site_summary"]
    n_high = trust["n_high_label_slide_trust"]
    n_review = trust["n_review_before_primary_analysis"]
    n_exclude = trust["n_exclude_from_primary_analysis"]
    n_tiles = int(sensitivity["n_slides_all_tissue"]) * 128
    n_overlap = agreement["n_overlap_slides"]
    n_low_zero = case_driver["n_low_zero_slides"]
    n_expected_3plus = case_driver["n_slides_expected_profile_3plus_views"]
    n_opposite_2plus = case_driver["n_slides_opposite_profile_2plus_views"]
    n_classifier_wrong_2plus = case_driver["n_slides_classifier_wrong_2plus_views"]
    n_local_erbb2 = local_erbb2["n_local_erbb2_cases"]
    n_high_erbb2 = local_erbb2["n_high_trust_with_local_erbb2"]
    n_low_zero_erbb2 = local_erbb2["n_low_zero_high_trust_with_local_erbb2"]
    n_erbb2_adjusted_sig = local_erbb2["n_erbb2_adjusted_channel_tests_q_lt_0_05"]
    n_mixed_sites = within_site["n_mixed_source_sites"]
    n_mixed_cases = within_site["n_mixed_cases"]
    n_mixed_low = within_site["n_mixed_low"]
    n_mixed_zero = within_site["n_mixed_zero"]

    return [
        (
            "BRCA HER2 Pathology AI Findings",
            f"""
**Simple display report for the current result.**

The current primary analysis uses **{n_high} strict high-trust TCGA-BRCA diagnostic H&E slides** after HER2 label, slide-file, OpenSlide, and female-patient primary-tumor filtering. The female-patient filter follows the relevant TCGA sample-selection principle from Guardia et al., Genome Research 2025, PMID 40664477; the rest of our slide QC is specific to this pathology-AI project. The raw GigaTIME run still contains the previous 174-slide inference output, but the current summaries filter those predictions to the strict {n_high}-slide analysis set.

**Most interesting result:** GigaTIME-derived H&E virtual mIF features reproducibly separate HER2-low from HER2-zero tumors. HER2-low is lower than HER2-zero for several virtual immune/myeloid/checkpoint and CK-associated channels.
""",
        ),
        (
            "One-Slide Takeaway",
            f"""
- Downloaded cohort checked: 183 TCGA-BRCA diagnostic slides.
- Strict high-trust primary analysis: {n_high} slides.
- Review before primary analysis: {n_review} slides.
- Excluded from primary analysis: {n_exclude} slides, all male HER2-positive cases.
- Groups: 53 HER2-positive, 57 HER2-low, 61 HER2-zero.
- GigaTIME sampling: 128 tissue tiles per primary-analysis slide, {n_tiles:,} tile predictions after filtering.
- Strongest finding: HER2-low versus HER2-zero, not HER2-positive diagnosis.
- Presentation-safe claim: image-derived GigaTIME features associate with HER2-low versus HER2-zero state differences; they do not yet provide clinical HER2 diagnosis or direct HER2 isoform detection.
""",
        ),
        (
            "Trustworthy Slide List",
            markdown_table(["Clinical HER2 group", "Strict high-trust slides"], content["group_counts"])
            + "\n\nThe strict list is in `docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv`. The three excluded cases are male HER2-positive TCGA-BRCA patients. Guardia et al. also excluded male TCGA-BRCA samples for their HER2 isoform/ADC-resistance analysis; our additional file-integrity, OpenSlide, and H&E tile checks are project-specific.",
        ),
        (
            "Strongest Pairwise Result",
            "Negative delta means HER2-low has lower mean virtual activation than HER2-zero.\n\n"
            + markdown_table(
                [
                    "Channel",
                    "HER2-low mean",
                    "HER2-zero mean",
                    "Low-zero delta",
                    "p",
                    "BH q",
                    "Cliff's delta",
                ],
                content["low_zero_pairwise"],
            )
            + "\n\nThis is the cleanest current statistical result: HER2-low is lower than HER2-zero for multiple GigaTIME virtual immune/myeloid/checkpoint and tissue-context channels after correction for multiple testing.",
        ),
        (
            "Three-Group Pattern",
            markdown_table(
                [
                    "Channel",
                    "Kruskal p",
                    "BH q",
                    "Highest group",
                    "Lowest group",
                    "Max-min mean",
                ],
                content["three_group"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_findings/clinical_her2_group_mean_heatmap.png",
                "Strict high-trust HER2 group mean heatmap",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_findings/clinical_her2_channel_boxplots.png",
                "Strict high-trust HER2 virtual channel boxplots",
            ),
        ),
        (
            "Cleanup View Result",
            "The signal is not simply a blank-tile artifact. It persists after requiring more cellular tissue and partially persists after CK enrichment, but it weakens in the strictest CK-enriched view.\n\n"
            + markdown_table(
                ["Cleanup view", "Channel", "Low-zero delta", "p", "BH q"],
                content["cleanup_low_zero"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_gigatime_cleanup/cleanup_key_channel_heatmap.png",
                "Strict high-trust cleanup key channel heatmap",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_gigatime_cleanup/cleanup_key_channel_boxplots.png",
                "Strict high-trust cleanup key channel boxplots",
            ),
        ),
        (
            "Top Cleanup Three-Group Signals",
            markdown_table(
                [
                    "Cleanup view",
                    "Channel",
                    "Kruskal p",
                    "BH q",
                    "Highest group",
                    "Lowest group",
                    "Max-min mean",
                ],
                content["cleanup_three_group"],
            ),
        ),
        (
            "Classifier Result",
            "This is diagnostic-model style evidence, but it is still exploratory and not clinical.\n\n"
            + markdown_table(
                [
                    "Input view",
                    "Best feature set",
                    "N",
                    "Accuracy",
                    "Balanced accuracy",
                    "Macro AUC",
                    "Sensitivity",
                    "Specificity",
                ],
                content["classifier_low_zero"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_cleaned_classifier/cleaned_classifier_best_by_view.png",
                "Strict high-trust cleaned classifier best model by view",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_cleaned_classifier/cleaned_classifier_low_zero_confusions.png",
                "Strict high-trust HER2-low versus HER2-zero confusion matrices",
            ),
        ),
        (
            "Where The Classifier Is Weak",
            markdown_table(
                [
                    "Task",
                    "Best input view",
                    "Best feature set",
                    "N",
                    "Balanced accuracy",
                    "Macro AUC",
                    "Sensitivity",
                    "Specificity",
                ],
                content["classifier_other_tasks"],
            )
            + "\n\nHER2-positive detection remains weak. The current GigaTIME/H&E signal is stronger for HER2-low versus HER2-zero than for classic HER2-positive diagnosis.",
        ),
        (
            "Case-Level Driver Check",
            f"The case-level check asks whether the HER2-low versus HER2-zero result is broad enough to believe or whether a few unusual slides carry the signal. A score was built from the significant low-zero virtual channels in each cleanup view. Higher scores are more HER2-zero-like; lower scores are more HER2-low-like.\n\n"
            + markdown_table(
                ["Clinical group", "N", "Mean zero-like score", "Median zero-like score"],
                content["case_driver_group_scores"],
            )
            + f"\n\n{n_expected_3plus} of {n_low_zero} HER2-low/HER2-zero slides matched the expected direction in at least 3 of 4 cleanup views. {n_opposite_2plus} slides showed the opposite profile in at least 2 views, and {n_classifier_wrong_2plus} cases were misclassified by the best low-vs-zero classifier in at least 2 cleanup views. This is promising, but it also gives us a concrete manual-review list.\n\n"
            + markdown_table(
                ["Cleanup view", "N signal channels", "Channels"],
                content["case_driver_signal"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_case_drivers/case_zero_like_score_by_view.png",
                "Case-level HER2-zero-like score by cleanup view",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_case_drivers/case_driver_view_stability.png",
                "Case-level driver view stability",
            ),
        ),
        (
            "Manual Review Shortlist",
            "These cases are the first ones to inspect in H&E and virtual mIF-like overlays. They are not failures; they are the cases most likely to reveal label noise, sampling problems, artifact, or interesting biology.\n\n"
            + markdown_table(
                [
                    "Case",
                    "Group",
                    "HER2 detail",
                    "Wrong classifier views",
                    "Opposite-profile views",
                    "All-tissue score",
                ],
                content["case_driver_review"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_case_drivers/case_driver_channel_heatmap.png",
                "Most label-consistent case-driver channel heatmap",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_case_drivers/case_driver_classifier_probability.png",
                "Classifier probability versus case-driver score",
            ),
        ),
        (
            "Visual QC Spot Check",
            "A small visual QC set was rendered from the case-driver shortlist. This is the most important cautionary update: the low-like tiles, including label-consistent HER2-low and opposite-profile HER2-zero examples, often have high tissue fraction but very low virtual CK, CD68, PD-L1, and CD11c. Visually, these can look stromal/collagen-rich rather than clearly tumor-rich. That means the current HER2-low versus HER2-zero signal may partly reflect tissue composition unless a pathologist confirms tumor-rich regions.\n\n"
            + markdown_table(
                [
                    "Review category",
                    "Group",
                    "Tiles",
                    "Mean tissue",
                    "Mean zero-like tile score",
                    "Mean CK",
                    "Mean CD68",
                    "Mean PD-L1",
                    "Mean CD11c",
                ],
                content["visual_qc_tiles"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/label_consistent_her2_zero_TCGA-AO-A128_he_vs_virtual_mif_qc.png",
                "Label-consistent HER2-zero visual QC panel",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_case_driver_visual_qc/opposite_profile_manual_review_TCGA-A2-A0EW_he_vs_virtual_mif_qc.png",
                "Opposite-profile HER2-zero visual QC panel",
            )
            + "\n\nFull visual QC report: `docs/clinical_her2_high_trust_tile128_case_driver_visual_qc.md`.",
        ),
        (
            "Tissue-Composition Sensitivity",
            "We quantified the visual QC caveat across all HER2-low and HER2-zero slides. A low-marker tile is a tissue tile whose mean virtual CK/CD68/PD-L1/CD11c/CD3/CD4/CD20/Ki67 marker burden is in the bottom quartile. HER2-low has substantially more low-marker tiles than HER2-zero.\n\n"
            + markdown_table(
                [
                    "Metric",
                    "HER2-low mean",
                    "HER2-zero mean",
                    "Low-zero delta",
                    "BH q",
                    "Cliff",
                ],
                content["tissue_composition"],
            )
            + "\n\nThe case-driver score also tracks marker/tissue composition very strongly. This is partly expected because the driver score is built from GigaTIME marker channels, but it confirms that the current signal is composition-sensitive.\n\n"
            + markdown_table(
                ["Metric", "N", "Spearman rho", "BH q"],
                content["tissue_correlations"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_tissue_composition/tissue_composition_low_zero_boxplots.png",
                "HER2-low versus HER2-zero tissue-composition boxplots",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_tissue_composition/case_driver_vs_tissue_composition.png",
                "Case-driver score versus tissue-composition metrics",
            ),
        ),
        (
            "Composition Adjustment",
            "The strongest caution is the adjusted model: after adjusting only for low-marker tile fraction, most HER2-low versus HER2-zero channel effects collapse. This does not make the project bad; it tells us the honest paper angle is tissue-context association unless tumor-rich/pathologist-approved tiles reproduce the result.\n\n"
            + markdown_table(
                [
                    "Channel",
                    "Unadjusted beta",
                    "Unadjusted q",
                    "Low-marker-adjusted beta",
                    "Low-marker-adjusted q",
                    "Mean-CK-adjusted beta",
                    "Mean-CK-adjusted q",
                ],
                content["composition_adjusted"],
            )
            + "\n\nIn a stricter absolute CK-high proxy view, some immune/T-cell channels remain different, but the signal weakens and fewer slides qualify.\n\n"
            + markdown_table(
                ["Channel", "N low", "N zero", "Low-zero delta", "BH q"],
                content["absolute_ck_pairwise"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_tissue_composition/composition_adjusted_channel_betas.png",
                "Composition-adjusted low-zero channel effects",
            )
            + "\n\nFull tissue-composition report: `docs/clinical_her2_high_trust_tile128_tissue_composition_sensitivity.md`.",
        ),
        (
            "Tumor-Rich Proxy Sensitivity",
            "We then tested stricter virtual tumor-rich proxy filters. These are still GigaTIME-derived proxies, not pathologist tumor annotations. The result is nuanced: individual low-versus-zero channel tests weaken in the strictest fixed-count CK views, but the multichannel low-vs-zero classifier remains around 0.71-0.76 balanced accuracy. That means the signal is not solved, but it is still worth pursuing with real tumor-rich review.\n\n"
            + markdown_table(
                [
                    "Proxy view",
                    "Group",
                    "Slides passing min tiles",
                    "Median retained tiles",
                    "Median retained fraction",
                ],
                content["tumor_proxy_retention"],
            )
            + "\n\nLow-versus-zero univariate channel effects under stricter proxy filters:\n\n"
            + markdown_table(
                ["Proxy view", "Channel", "N low", "N zero", "Low-zero delta", "BH q"],
                content["tumor_proxy_pairwise"],
            )
            + "\n\nLow-versus-zero classifier under the same proxy filters:\n\n"
            + markdown_table(
                [
                    "Proxy view",
                    "Best feature set",
                    "N",
                    "Accuracy",
                    "Balanced accuracy",
                    "Macro AUC",
                    "Sensitivity",
                    "Specificity",
                ],
                content["tumor_proxy_classifier"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_tumor_proxy_sensitivity/tumor_proxy_low_zero_channel_deltas.png",
                "Tumor-proxy low-zero channel deltas",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_tumor_proxy_sensitivity/tumor_proxy_low_zero_classifier_feature_sets.png",
                "Tumor-proxy low-zero classifier feature sets",
            )
            + "\n\nFull tumor-proxy report: `docs/clinical_her2_high_trust_tile128_tumor_proxy_sensitivity.md`.",
        ),
        (
            "Classifier Permutation Sanity Check",
            "We also asked a blunt question: if we shuffle the HER2-low/HER2-zero labels, do these classifiers still look good? They should not. This is not a final clinical validation because the feature views were selected after earlier analyses, but it is an important sanity check that the classifier signal is not obviously random label structure.\n\n"
            + markdown_table(
                [
                    "Proxy view",
                    "Feature set",
                    "N",
                    "Features",
                    "LOOCV balanced accuracy",
                    "Repeated-CV balanced accuracy",
                    "Null mean",
                    "Null 95%",
                    "Empirical p",
                    "BH q",
                    "Repeated-CV AUC",
                ],
                content["classifier_permutation"],
            )
            + "\n\nEvery tested view stayed above its shuffled-label null distribution with empirical p = 0.0099 and BH q = 0.0099. This supports a real label-associated image pattern, but it still does not solve the tissue-composition caveat or prove HER2 isoform biology.\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_classifier_permutation/classifier_permutation_balanced_accuracy_null.png",
                "Classifier balanced accuracy versus shuffled-label null",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_classifier_permutation/classifier_permutation_auc_null.png",
                "Classifier AUC versus shuffled-label null",
            )
            + "\n\nFull classifier permutation report: `docs/clinical_her2_high_trust_tile128_classifier_permutation_sanity.md`.",
        ),
        (
            "Nested Classifier Model Selection",
            "The previous permutation test fixed the feature set that had already been selected. This stricter check chooses the best feature set inside each training fold, then evaluates the held-out fold. In plain language: the model has to pick its recipe without seeing the test samples.\n\n"
            + markdown_table(
                [
                    "Proxy view",
                    "N",
                    "Most selected feature set",
                    "Selected folds",
                    "Nested balanced accuracy",
                    "Nested AUC",
                    "Null mean",
                    "Null 95%",
                    "Empirical p",
                    "BH q",
                ],
                content["nested_classifier"],
            )
            + "\n\nThe nested balanced accuracies stay around 67-72%, and all tested views remain above the fully nested shuffled-label null. This is stronger classifier-methodology evidence than the post-hoc permutation check, but it is still internal validation rather than clinical proof.\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_nested_classifier/nested_classifier_observed_vs_null.png",
                "Nested classifier observed versus shuffled-label null",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_nested_classifier/nested_classifier_feature_selection_frequency.png",
                "Nested classifier feature selection frequency",
            )
            + "\n\nFull nested classifier report: `docs/clinical_her2_high_trust_tile128_nested_classifier_model_selection.md`.",
        ),
        (
            "Clinical/Site Confounder Check",
            "This is the biggest cautionary result so far. HER2-low and HER2-zero slides are not balanced for TCGA source site or slide size. Non-image source-site and slide-size covariates alone can classify HER2-low versus HER2-zero better than the GigaTIME image features. That means the current TCGA signal may partly reflect cohort construction, scanner/site differences, or slide-selection artifacts.\n\n"
            + markdown_table(
                ["Covariate", "HER2-low mean", "HER2-zero mean", "Low-zero delta", "p"],
                content["covariate_balance"],
            )
            + "\n\nTop 8 CK tile classifier comparison:\n\n"
            + markdown_table(
                ["Feature set", "Features", "Balanced accuracy", "AUC"],
                content["covariate_classifier"],
            )
            + "\n\nQC-cellular channel tests after clinical/site adjustment:\n\n"
            + markdown_table(
                ["Adjustment model", "q<0.05 channels", "Channels", "Best q"],
                content["covariate_adjusted"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_clinical_covariates/clinical_covariate_balance.png",
                "HER2-low versus HER2-zero covariate balance",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_clinical_covariates/clinical_covariate_classifier_comparison.png",
                "Clinical, source-site, slide-size, and GigaTIME classifier comparison",
            )
            + "\n\nFull covariate sensitivity report: `docs/clinical_her2_high_trust_tile128_clinical_covariate_sensitivity.md`.",
        ),
        (
            "Matched Low-Versus-Zero Sensitivity",
            "After finding strong TCGA source-site and slide-size imbalance, we built matched HER2-low/HER2-zero subsets. This is the most important new check because it asks whether the image signal still exists when obvious technical/source imbalance is reduced.\n\n"
            + markdown_table(
                [
                    "Matched subset",
                    "Pairs",
                    "Same-source-site pairs",
                    "Median abs log-size diff",
                    "Median abs MB diff",
                ],
                content["matched_pairs"],
            )
            + "\n\nTop 8 CK tile leave-one-pair-out classifier comparison:\n\n"
            + markdown_table(
                ["Matched subset", "Feature set", "Pairs", "Features", "Balanced accuracy", "AUC"],
                content["matched_classifier"],
            )
            + "\n\nTop paired channel tests in the same top 8 CK proxy view:\n\n"
            + markdown_table(
                ["Matched subset", "Channel", "Pairs", "Mean low-zero", "Wilcoxon p", "BH q"],
                content["matched_channels"],
            )
            + "\n\nInterpretation: GigaTIME mean channels remain modestly above chance in the matched subsets, around 0.675-0.708 balanced accuracy. But source-site and slide-size baselines remain competitive or stronger in the larger matched subsets. No paired channel test reaches BH q < 0.05. This keeps the HER2-low/HER2-zero signal worth studying, but it makes the current TCGA result unsafe as an independent HER2 biology claim.\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_matched_low_zero/matched_low_zero_classifier_sensitivity.png",
                "Matched HER2-low versus HER2-zero classifier sensitivity",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_matched_low_zero/matched_low_zero_channel_q_counts.png",
                "Matched HER2-low versus HER2-zero paired channel q counts",
            )
            + "\n\nFull matched sensitivity report: `docs/clinical_her2_high_trust_tile128_matched_low_zero_sensitivity.md`.",
        ),
        (
            "Source-Site Generalization",
            "We then asked whether the low-vs-zero classifier travels to held-out TCGA source sites. This is stricter than ordinary random cross-validation because all cases from one source site are removed from training and used as the test set.\n\n"
            + markdown_table(
                ["TSS", "N HER2-low", "N HER2-zero", "N cases", "Both classes"],
                content["source_site_balance"],
            )
            + "\n\nTop 8 CK tile classifier comparison:\n\n"
            + markdown_table(
                ["Feature set", "Validation", "Features", "Balanced accuracy", "AUC"],
                content["source_site_generalization"],
            )
            + "\n\nGigaTIME mean-channel performance drop by feature view:\n\n"
            + markdown_table(
                ["Feature view", "Random CV BA", "Leave-site-out BA", "Difference"],
                content["source_site_gigatime_drop"],
            )
            + "\n\nInterpretation: GigaTIME mean-channel performance drops under source-site holdout across every tested view. In the top 8 CK proxy view it drops from 74.5% balanced accuracy under repeated stratified CV to 66.9% under leave-source-site-out validation. Slide-size covariates remain very strong even under source-site holdout, reaching 88.2% balanced accuracy. This is a major warning that the current TCGA classifier is not source-independent HER2 biology.\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_source_site_generalization/source_site_generalization_balanced_accuracy.png",
                "Source-site held-out classifier balanced accuracy",
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_source_site_generalization/source_site_generalization_drop_by_view.png",
                "GigaTIME source-site holdout performance drop by view",
            )
            + "\n\nFull source-site generalization report: `docs/clinical_her2_high_trust_tile128_source_site_generalization.md`.",
        ),
        (
            "Within-Source-Site Sensitivity",
            f"Because source-site confounding is so strong, we restricted the analysis to TCGA source sites that contain both HER2-low and HER2-zero cases. Only {n_mixed_sites} source sites qualify, giving {n_mixed_cases} cases total: {n_mixed_low} HER2-low and {n_mixed_zero} HER2-zero. This is small and imbalanced, so it is a stress test rather than proof.\n\n"
            + markdown_table(
                ["TSS", "HER2-low", "HER2-zero", "Cases"],
                content["within_source_site_balance"],
            )
            + "\n\nSite-fixed channel tests ask whether HER2-zero remains higher than HER2-low after adding TCGA source-site fixed effects:\n\n"
            + markdown_table(
                ["Feature view", "Channels tested", "Channels q<0.05", "Best BH q"],
                content["within_source_site_channel_q"],
            )
            + "\n\nTop 8 CK proxy view, top site-fixed channel tests:\n\n"
            + markdown_table(
                ["Feature view", "Channel", "N low", "N zero", "Beta zero-vs-low", "p", "BH q"],
                content["within_source_site_top_channels"],
            )
            + "\n\nMixed-site classifier sensitivity in the top 8 CK proxy view:\n\n"
            + markdown_table(
                ["Feature set", "Validation", "N", "Features", "Balanced accuracy", "AUC", "Sensitivity", "Specificity"],
                content["within_source_site_classifier"],
            )
            + "\n\nInterpretation: this partially supports continued investigation but does not solve the confounding problem. Some site-fixed channel effects remain in QC-cellular and CK-high views, and all-channel GigaTIME classifiers retain some above-chance performance in the mixed-site subset. But the subset is tiny, key-marker-only classifiers weaken, and specificity is low. This is evidence for a possible within-site image pattern, not enough for a source-independent HER2 biology claim.\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_within_source_site/within_source_site_classifier_sensitivity.png",
                "Within-source-site classifier sensitivity",
            )
            + "\n\nFull within-source-site report: `docs/clinical_her2_high_trust_tile128_within_source_site_low_zero.md`.",
        ),
        (
            "ER/PR Sensitivity",
            markdown_table(
                [
                    "Cleanup view",
                    "Unadjusted q<0.05 channels",
                    "ER/PR adjusted q<0.05 channels",
                    "ER/PR+ERBB2 adjusted q<0.05 channels",
                ],
                content["erpr_summary"],
            )
            + "\n\nAll-sampled-tissue channels ranked by ER/PR-adjusted q value:\n\n"
            + markdown_table(
                [
                    "Channel",
                    "Unadjusted low-zero delta",
                    "Unadjusted q",
                    "ER/PR beta",
                    "ER/PR q",
                    "ER/PR+ERBB2 q",
                    "RNA n",
                ],
                content["erpr_adjusted"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_erpr_subgroup_sensitivity/erpr_adjusted_low_zero_q_heatmap.png",
                "ER/PR-adjusted low-zero q-value heatmap",
            ),
        ),
        (
            "Run Agreement",
            f"The current strict 128-tile analysis was compared with the earlier 60-slide 256-tile run on {n_overlap} overlapping slide IDs. All 8 tested key channels kept the same HER2-low versus HER2-zero direction; 7 of 8 kept HER2-low lower than HER2-zero.\n\n"
            + markdown_table(
                [
                    "Channel",
                    "256-tile low-zero delta",
                    "128-tile low-zero delta",
                    "Same direction",
                    "Both low lower",
                ],
                content["direction_comparison"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_vs_expanded20_tile256/low_zero_delta_tile128_vs_tile256.png",
                "HER2-low versus zero delta comparison",
            ),
        ),
        (
            "Local ERBB2 Gene-Level Validation",
            f"We extracted ERBB2 gene-level TPM from all local GDC STAR augmented gene-count files. This gives ERBB2 RNA context for {n_local_erbb2} TCGA-BRCA cases, including {n_high_erbb2} strict high-trust GigaTIME/HER2 cases and {n_low_zero_erbb2} HER2-low/HER2-zero high-trust cases. This is useful, but it is not isoform validation: these files do not contain transcript-level isoform proportions, PSI, junction evidence, or antibody-binding-domain information.\n\n"
            + "Clinical HER2 group summary from local STAR ERBB2 TPM:\n\n"
            + markdown_table(
                ["Clinical HER2 group", "N", "Median TPM", "Q25 TPM", "Q75 TPM"],
                content["local_erbb2_group_summary"],
            )
            + "\n\nPairwise gene-level ERBB2 tests:\n\n"
            + markdown_table(
                ["Comparison", "N A", "N B", "Median A", "Median B", "AUC", "p", "BH q"],
                content["local_erbb2_pairwise"],
            )
            + "\n\nSimple ERBB2-only reference classifier:\n\n"
            + markdown_table(
                ["Task", "N cases", "AUC", "Best-threshold balanced accuracy"],
                content["local_erbb2_reference_classifier"],
            )
            + "\n\nInterpretation: ERBB2 RNA strongly supports the HER2-positive label as a molecular sanity check, but it weakly separates HER2-low from HER2-zero. The low/zero ERBB2-only AUC is only about 0.605, so the current GigaTIME low/zero signal is not simply a strong gene-level ERBB2 expression separation.\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_local_erbb2_validation/local_erbb2_by_clinical_her2_group.png",
                "Local ERBB2 by clinical HER2 group",
            )
            + "\n\nLow/zero GigaTIME channel correlations with gene-level ERBB2 are limited:\n\n"
            + markdown_table(
                ["Feature view", "Channel", "N", "Spearman rho", "p", "BH q"],
                content["local_erbb2_correlations"],
            )
            + f"\n\nAfter adjusting low-vs-zero channel tests for log ERBB2 TPM in the small RNA-overlap subset, {n_erbb2_adjusted_sig} channel/view effects remain BH q < 0.05. This supports the idea that the GigaTIME low/zero signal is not fully explained by gene-level ERBB2 RNA, but it does not remove the stronger source-site, slide-size, and tissue-composition caveats.\n\n"
            + markdown_table(
                ["Feature view", "Channel", "N", "Beta zero-vs-low adjusted", "p", "BH q"],
                content["local_erbb2_adjusted"],
            )
            + "\n\n"
            + md_image(
                "../docs/assets/clinical_her2_high_trust_tile128_local_erbb2_validation/local_erbb2_gigatime_correlation_heatmap.png",
                "Local ERBB2 GigaTIME correlation heatmap",
            )
            + "\n\nFull local ERBB2 validation report: `docs/clinical_her2_high_trust_tile128_local_erbb2_validation.md`.",
        ),
        (
            "HER2 Isoform Feasibility",
            "After correcting the target paper to Guardia et al., we checked whether the local project can directly validate HER2 isoform biology. The answer is no for the current local files. We have gene-level STAR count/TPM files, which can support ERBB2 expression and RNA-program context, but we do not have transcript-level isoform quantification, raw reads, BAMs, or junction-count outputs in the current TCGA workspace.\n\n"
            + markdown_table(["Audit item", "Value"], content["isoform_audit"])
            + "\n\n"
            + markdown_table(
                ["Analysis", "Status", "What it can support", "Next action"],
                content["isoform_feasibility"],
            )
            + "\n\nThe safe conclusion is that GigaTIME features can be tested for association with HER2 isoform/state hypotheses only after we obtain sample-level isoform labels or reproduce the RNA-seq isoform workflow. The current H&E/GigaTIME work does not detect HER2 isoforms.\n\nFull feasibility report: `docs/her2_isoform_validation_feasibility.md`.",
        ),
        (
            "What We Learned",
            """
1. The strongest research signal is HER2-low versus HER2-zero, not HER2-positive diagnosis.
2. HER2-low tends to show lower GigaTIME virtual immune/checkpoint, myeloid, and CK-associated activation than HER2-zero.
3. The signal survives cellular-tissue cleanup, so it is unlikely to be only blank/background artifact.
4. The signal weakens in strict CK-enriched tumor-focused views, suggesting broader tissue context may matter.
5. The signal mostly survives ER/PR adjustment and remains directionally stable across overlapping 128-tile and 256-tile runs.
6. The case-level driver check shows a real but imperfect pattern: many cases are stable, while a substantial manual-review subset may contain label noise, artifacts, sampling problems, or biological exceptions.
7. The visual QC spot check raises a serious tissue-composition caveat: low-like tiles can be stromal/collagen-rich and not clearly tumor-rich.
8. The quantified tissue-composition analysis strengthens that caveat: HER2-low has more low-marker tiles, and low-marker adjustment removes most channel effects.
9. Stricter virtual tumor-rich proxy filters weaken individual channel tests but preserve exploratory low-vs-zero classifier performance, so the signal remains worth validating with real tumor-rich/pathologist review.
10. Shuffled-label classifier checks show the low-vs-zero classifiers are not obviously random, but this remains post-hoc exploratory evidence.
11. Nested model-selection checks reduce feature-set selection bias and still show above-null low-vs-zero classifier performance.
12. Clinical/site covariate checks show a serious TCGA confounding risk: slide size and source site classify HER2-low versus HER2-zero better than GigaTIME features.
13. Matched low-vs-zero sensitivity keeps modest GigaTIME signal but does not eliminate the confounder concern; source-site/slide-size baselines remain competitive and paired channel tests are not FDR-significant.
14. Source-site held-out validation weakens GigaTIME performance while slide-size covariates stay very strong, so the current classifier is not source-independent HER2 biology.
15. Within-source-site sensitivity partially supports continued investigation: some site-fixed channel effects and all-channel classifiers remain above chance, but the mixed-site subset is tiny and imbalanced.
16. Local gene-level ERBB2 RNA strongly validates the broad HER2-positive label but weakly separates HER2-low from HER2-zero, so the low/zero GigaTIME signal is not simply a strong ERBB2 expression split.
17. The current local RNA files are gene-level STAR count/TPM files, not transcript-level isoform or junction-count data, so Guardia-style HER2 isoform validation is not available yet.
18. We should say image AI predicts or associates with HER2 state hypotheses, not that image AI directly detects HER2 isoforms.
""",
        ),
        (
            "Best Presentation Wording",
            """
In a strict high-trust female TCGA-BRCA diagnostic H&E cohort, using HER2 isoform/state biology from Guardia et al. as the motivating context, GigaTIME-derived virtual mIF features reproducibly associate with the HER2-low versus HER2-zero boundary. HER2-low tumors show lower virtual immune/myeloid/checkpoint and CK-associated signals than HER2-zero tumors, and low-vs-zero classifiers beat shuffled-label null tests even under nested feature-set selection. Local gene-level ERBB2 RNA validates the broad HER2-positive label but only weakly separates HER2-low from HER2-zero, suggesting the image signal is not just a strong ERBB2 expression split. However, clinical/site covariate, matched-subset, source-site held-out, and within-source-site sensitivity checks reveal major TCGA confounding risk: source-site and slide-size covariates classify HER2-low versus HER2-zero very well, matching does not fully remove the concern, GigaTIME performance drops when entire source sites are held out, and the mixed-site subset is small and imbalanced. These findings support a hypothesis-generating tissue-context observation and a clear validation plan requiring external/site-balanced cohorts, tumor-rich/pathologist review, and molecular validation; they do not support clinical HER2 diagnosis or direct HER2 isoform detection.
""",
        ),
        (
            "Best Next Step",
            """
The next strongest move is to validate outside the current confounded TCGA slice while doing manual pathology/QC review:

1. Open the strongest label-consistent HER2-low and HER2-zero driver cases.
2. Open the highest-priority opposite-profile and classifier-error cases.
3. Ask a pathologist to review whether high-signal H&E regions are tumor-rich, immune-rich, stromal, necrotic, or artifact-prone.
4. Rerun the key statistics using tumor-rich or pathologist-approved tile subsets before making stronger HER2-biology claims.
5. Add tumor purity, subtype, stain/batch, and immune-deconvolution covariates if available.
6. Look for external H&E plus real HER2 IHC/ISH/mIF or transcript/protein validation data.
7. If isoform data are available, test whether GigaTIME features associate with ERBB2 isoform groups rather than only clinical HER2 labels.
""",
        ),
    ]


def build_notebook(args: argparse.Namespace, content: dict[str, object]) -> None:
    nb_path = Path(args.out_notebook)
    nb_path.parent.mkdir(parents=True, exist_ok=True)
    cells = [
        notebook_cell(f"# {title}\n\n{body}" if index == 0 else f"## {title}\n\n{body}")
        for index, (title, body) in enumerate(build_markdown_sections(content))
    ]
    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    nb_path.write_text(json.dumps(notebook, indent=2) + "\n", encoding="utf-8")


def paragraphize_markdown(body: str) -> str:
    blocks = []
    for raw_block in body.strip().split("\n\n"):
        block = raw_block.strip()
        if not block:
            continue
        if block.startswith("| "):
            lines = block.splitlines()
            headers = [item.strip() for item in lines[0].strip("|").split("|")]
            rows = [[item.strip() for item in line.strip("|").split("|")] for line in lines[2:]]
            blocks.append(html_table(headers, rows))
        elif block.startswith("!["):
            alt = block.split("](", 1)[0][2:]
            src = block.split("](", 1)[1].rstrip(")")
            blocks.append(f"<img src='{html.escape(src)}' alt='{html.escape(alt)}'>")
        elif block.startswith("- "):
            items = block.splitlines()
            blocks.append("<ul>" + "".join(f"<li>{html.escape(item[2:])}</li>" for item in items) + "</ul>")
        elif block.startswith("1. "):
            items = block.splitlines()
            html_items = []
            for item in items:
                _, text = item.split(". ", 1)
                html_items.append(f"<li>{html.escape(text)}</li>")
            blocks.append("<ol>" + "".join(html_items) + "</ol>")
        else:
            escaped = html.escape(block).replace("\n", "<br>")
            escaped = escaped.replace("**", "")
            blocks.append(f"<p>{escaped}</p>")
    return "\n".join(blocks)


def section(title: str, body: str) -> str:
    return f"<section><h2>{html.escape(title)}</h2>{paragraphize_markdown(body)}</section>"


def build_html(args: argparse.Namespace, content: dict[str, object]) -> None:
    html_path = Path(args.out_html)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    sections = build_markdown_sections(content)
    title, hero_body = sections[0]
    trust = content["trust_summary"]
    sensitivity = content["sensitivity_summary"]
    agreement = content["agreement_summary"]
    n_high = trust["n_high_label_slide_trust"]
    n_tiles = int(sensitivity["n_slides_all_tissue"]) * 128
    n_overlap = agreement["n_overlap_slides"]
    css = """
:root {
  color-scheme: light;
  --ink: #17202a;
  --muted: #536171;
  --line: #d9e0e8;
  --paper: #ffffff;
  --bg: #f5f7fa;
  --accent: #0f766e;
  --accent-2: #1d4ed8;
  --warn: #b45309;
}
body {
  margin: 0;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  color: var(--ink);
  background: var(--bg);
}
main {
  max-width: 1160px;
  margin: 0 auto;
  padding: 34px 24px 72px;
}
.hero {
  background: var(--paper);
  border: 1px solid var(--line);
  border-top: 6px solid var(--accent);
  border-radius: 8px;
  padding: 30px;
}
.hero h1 {
  margin: 0 0 12px;
  font-size: clamp(28px, 4vw, 42px);
  line-height: 1.08;
  letter-spacing: 0;
}
.hero p {
  max-width: 980px;
  margin: 12px 0 0;
  color: var(--muted);
  font-size: 18px;
  line-height: 1.55;
}
.stats {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
  gap: 12px;
  margin-top: 18px;
}
.stat {
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 13px 14px;
  background: #fbfcfe;
}
.stat strong {
  display: block;
  font-size: 23px;
  color: var(--accent);
}
.stat span {
  color: var(--muted);
  font-size: 13px;
}
section {
  background: var(--paper);
  border: 1px solid var(--line);
  border-radius: 8px;
  margin-top: 18px;
  padding: 24px;
}
section:nth-of-type(2) {
  border-top: 5px solid var(--accent-2);
}
h2 {
  margin: 0 0 14px;
  font-size: 23px;
  line-height: 1.2;
  letter-spacing: 0;
}
p, li {
  font-size: 16px;
  line-height: 1.58;
}
p {
  margin: 11px 0;
}
ul, ol {
  margin: 10px 0 0;
  padding-left: 22px;
}
table {
  width: 100%;
  border-collapse: collapse;
  margin: 14px 0;
  font-size: 14px;
  line-height: 1.35;
}
th, td {
  text-align: left;
  border-bottom: 1px solid var(--line);
  padding: 9px 10px;
  vertical-align: top;
}
th {
  background: #eef2f6;
  font-weight: 650;
}
tr:nth-child(even) td {
  background: #fbfcfe;
}
img {
  display: block;
  width: 100%;
  max-width: 1040px;
  height: auto;
  margin: 16px auto;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: white;
}
.note {
  margin-top: 18px;
  padding: 14px 16px;
  border-left: 5px solid var(--warn);
  background: #fff7ed;
}
@media (max-width: 720px) {
  main { padding: 18px 12px 48px; }
  .hero, section { padding: 18px; }
  table { font-size: 12px; }
  th, td { padding: 7px 6px; }
}
"""
    hero_html = f"""
<div class="hero">
  <h1>{html.escape(title)}</h1>
  {paragraphize_markdown(hero_body)}
  <div class="stats">
    <div class="stat"><strong>{n_high}</strong><span>strict high-trust slides</span></div>
    <div class="stat"><strong>53/57/61</strong><span>HER2-positive / low / zero</span></div>
    <div class="stat"><strong>{n_tiles:,}</strong><span>primary-analysis tile predictions</span></div>
    <div class="stat"><strong>{n_overlap}</strong><span>overlap slides for 128-vs-256 check</span></div>
  </div>
  <div class="note">Use this file as the simple presentation notebook. The result is promising, but still not validated for clinical diagnosis.</div>
</div>
"""
    body_sections = "\n".join(section(title, body) for title, body in sections[1:])
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>BRCA HER2 Pathology AI Findings</title>
  <style>{css}</style>
</head>
<body>
<main>
{hero_html}
{body_sections}
</main>
</body>
</html>
"""
    html_path.write_text(html_text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    content = build_content(args)
    build_notebook(args, content)
    build_html(args, content)
    print(f"Wrote {args.out_notebook}")
    print(f"Wrote {args.out_html}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
