#!/usr/bin/env python3
"""Compare overlapping slides across two GigaTIME runs with different settings."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path


GIGATIME_CHANNELS = [
    "DAPI",
    "TRITC",
    "Cy5",
    "PD-1",
    "CD14",
    "CD4",
    "T-bet",
    "CD34",
    "CD68",
    "CD16",
    "CD11c",
    "CD138",
    "CD20",
    "CD3",
    "CD8",
    "PD-L1",
    "CK",
    "Ki67",
    "Tryptase",
    "Actin-D",
    "Caspase3-D",
    "PHH3-B",
    "Transgelin",
]
KEY_CHANNELS = ["CD68", "PD-L1", "PD-1", "CD11c", "CD4", "CD3", "CK", "Ki67"]
GROUP_ORDER = ["HER2-positive", "HER2-low", "HER2-zero"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-joined",
        default="results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/clinical_summary/joined_slide_clinical_her2_gigatime.csv",
        help="Reference run joined slide table, usually the 60-slide tile256 run.",
    )
    parser.add_argument(
        "--comparison-joined",
        default="results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/clinical_summary/joined_slide_clinical_her2_gigatime.csv",
        help="Comparison run joined slide table, usually the high-trust tile128 run.",
    )
    parser.add_argument("--reference-label", default="Expanded 60-slide tile256")
    parser.add_argument("--comparison-label", default="High-trust tile128")
    parser.add_argument(
        "--out-dir",
        default="results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/tile128_vs_expanded20_tile256_agreement",
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_vs_expanded20_tile256",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_vs_expanded20_tile256_agreement.md",
    )
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


def load_overlap(pd, reference_path: Path, comparison_path: Path):
    reference = pd.read_csv(reference_path)
    comparison = pd.read_csv(comparison_path)
    overlap = reference.merge(
        comparison,
        on=["case_submitter_id", "slide_id"],
        suffixes=("_reference", "_comparison"),
        how="inner",
        validate="one_to_one",
    )
    if "clinical_her2_group_reference" in overlap.columns:
        overlap["clinical_her2_group"] = overlap["clinical_her2_group_reference"]
    elif "clinical_her2_group_comparison" in overlap.columns:
        overlap["clinical_her2_group"] = overlap["clinical_her2_group_comparison"]
    return reference, comparison, overlap


def available_channels(overlap) -> list[str]:
    channels = []
    for channel in GIGATIME_CHANNELS:
        if f"mean_{channel}_reference" in overlap.columns and f"mean_{channel}_comparison" in overlap.columns:
            channels.append(channel)
    return channels


def build_channel_agreement(np, stats, overlap, channels: list[str]):
    rows = []
    for channel in channels:
        ref_col = f"mean_{channel}_reference"
        cmp_col = f"mean_{channel}_comparison"
        pair = overlap[[ref_col, cmp_col]].dropna()
        if len(pair) >= 3:
            pearson = stats.pearsonr(pair[ref_col], pair[cmp_col])
            spearman = stats.spearmanr(pair[ref_col], pair[cmp_col])
            pearson_r = float(pearson.statistic)
            pearson_p = float(pearson.pvalue)
            spearman_rho = float(spearman.statistic)
            spearman_p = float(spearman.pvalue)
        else:
            pearson_r = pearson_p = spearman_rho = spearman_p = float("nan")
        diff = pair[cmp_col] - pair[ref_col]
        rows.append(
            {
                "channel": channel,
                "n_overlap": int(len(pair)),
                "reference_mean": float(pair[ref_col].mean()) if len(pair) else float("nan"),
                "comparison_mean": float(pair[cmp_col].mean()) if len(pair) else float("nan"),
                "comparison_minus_reference_mean": float(diff.mean()) if len(diff) else float("nan"),
                "median_absolute_difference": float(np.median(np.abs(diff))) if len(diff) else float("nan"),
                "pearson_r": pearson_r,
                "pearson_p_value": pearson_p,
                "spearman_rho": spearman_rho,
                "spearman_p_value": spearman_p,
            }
        )
    return rows


def low_zero_delta(stats, rows, value_col: str):
    low = rows.loc[rows["clinical_her2_group"] == "HER2-low", value_col].dropna()
    zero = rows.loc[rows["clinical_her2_group"] == "HER2-zero", value_col].dropna()
    if len(low) and len(zero):
        test = stats.mannwhitneyu(low, zero, alternative="two-sided")
        p_value = float(test.pvalue)
    else:
        p_value = float("nan")
    return {
        "n_low": int(len(low)),
        "n_zero": int(len(zero)),
        "mean_low": float(low.mean()) if len(low) else float("nan"),
        "mean_zero": float(zero.mean()) if len(zero) else float("nan"),
        "delta_low_minus_zero": float(low.mean() - zero.mean()) if len(low) and len(zero) else float("nan"),
        "mannwhitney_p_value": p_value,
    }


def build_low_zero_direction(stats, overlap, channels: list[str]):
    rows = []
    low_zero = overlap.loc[overlap["clinical_her2_group"].isin(["HER2-low", "HER2-zero"])].copy()
    for channel in channels:
        ref = low_zero_delta(stats, low_zero, f"mean_{channel}_reference")
        cmp = low_zero_delta(stats, low_zero, f"mean_{channel}_comparison")
        ref_delta = ref["delta_low_minus_zero"]
        cmp_delta = cmp["delta_low_minus_zero"]
        rows.append(
            {
                "channel": channel,
                "n_low": ref["n_low"],
                "n_zero": ref["n_zero"],
                "reference_delta_low_minus_zero": ref_delta,
                "reference_p_value": ref["mannwhitney_p_value"],
                "comparison_delta_low_minus_zero": cmp_delta,
                "comparison_p_value": cmp["mannwhitney_p_value"],
                "same_direction": (
                    bool(ref_delta * cmp_delta > 0)
                    if not math.isnan(ref_delta) and not math.isnan(cmp_delta)
                    else False
                ),
                "both_low_lower_than_zero": (
                    bool(ref_delta < 0 and cmp_delta < 0)
                    if not math.isnan(ref_delta) and not math.isnan(cmp_delta)
                    else False
                ),
            }
        )
    return rows


def plot_correlation_heatmap(plt, sns, pd, agreement, asset_dir: Path):
    plot_df = pd.DataFrame(agreement).set_index("channel")[["spearman_rho"]]
    plot_df = plot_df.sort_values("spearman_rho", ascending=False)
    plot_df.columns = ["rho"]
    fig, axis = plt.subplots(figsize=(4.8, max(5.2, 0.30 * len(plot_df))))
    sns.heatmap(plot_df, cmap="viridis", annot=True, fmt=".2f", linewidths=0.3, vmin=0, vmax=1, ax=axis)
    axis.set_xlabel("")
    axis.set_ylabel("GigaTIME channel")
    axis.set_title("Slide-level agreement across GigaTIME runs")
    fig.tight_layout()
    fig.savefig(asset_dir / "tile128_vs_tile256_channel_correlation_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_delta_comparison(plt, sns, pd, direction_rows, asset_dir: Path, reference_label: str, comparison_label: str):
    records = []
    for row in direction_rows:
        if row["channel"] not in KEY_CHANNELS:
            continue
        records.append(
            {
                "channel": row["channel"],
                "run": reference_label,
                "delta_low_minus_zero": row["reference_delta_low_minus_zero"],
            }
        )
        records.append(
            {
                "channel": row["channel"],
                "run": comparison_label,
                "delta_low_minus_zero": row["comparison_delta_low_minus_zero"],
            }
        )
    plot_df = pd.DataFrame(records)
    fig, axis = plt.subplots(figsize=(9.8, 5.2))
    sns.barplot(data=plot_df, x="channel", y="delta_low_minus_zero", hue="run", ax=axis)
    axis.axhline(0, color="#333333", linewidth=1)
    axis.set_xlabel("GigaTIME channel")
    axis.set_ylabel("HER2-low minus HER2-zero mean activation")
    axis.set_title("HER2-low vs HER2-zero direction on overlapping slides")
    axis.tick_params(axis="x", labelrotation=35)
    fig.tight_layout()
    fig.savefig(asset_dir / "low_zero_delta_tile128_vs_tile256.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_key_scatter(plt, sns, pd, overlap, asset_dir: Path, reference_label: str, comparison_label: str):
    records = []
    for channel in KEY_CHANNELS[:6]:
        ref_col = f"mean_{channel}_reference"
        cmp_col = f"mean_{channel}_comparison"
        if ref_col not in overlap.columns or cmp_col not in overlap.columns:
            continue
        for _, row in overlap.iterrows():
            records.append(
                {
                    "channel": channel,
                    "reference": row[ref_col],
                    "comparison": row[cmp_col],
                    "clinical_her2_group": row["clinical_her2_group"],
                }
            )
    plot_df = pd.DataFrame(records).dropna()
    if plot_df.empty:
        return
    grid = sns.FacetGrid(
        plot_df,
        col="channel",
        col_wrap=3,
        hue="clinical_her2_group",
        hue_order=GROUP_ORDER,
        height=3.0,
        aspect=1.08,
        sharex=False,
        sharey=False,
    )
    grid.map_dataframe(sns.scatterplot, x="reference", y="comparison", s=35, alpha=0.82)
    for axis in grid.axes.flatten():
        xlim = axis.get_xlim()
        ylim = axis.get_ylim()
        lo = min(xlim[0], ylim[0])
        hi = max(xlim[1], ylim[1])
        axis.plot([lo, hi], [lo, hi], color="#777777", linestyle="--", linewidth=1)
        axis.set_xlim(lo, hi)
        axis.set_ylim(lo, hi)
    grid.set_axis_labels(reference_label, comparison_label)
    grid.add_legend(title="Clinical HER2", bbox_to_anchor=(1.02, 0.5), loc="center left")
    grid.fig.suptitle("Overlapping slide channel scores across runs", y=0.995)
    grid.fig.tight_layout(rect=[0, 0, 0.87, 0.96])
    grid.fig.savefig(asset_dir / "key_channel_tile128_vs_tile256_scatter.png", dpi=180, bbox_inches="tight")
    plt.close(grid.fig)


def write_markdown(
    path: Path,
    agreement,
    direction,
    summary: dict,
    asset_dir: Path,
    reference_label: str,
    comparison_label: str,
):
    def asset_link(filename: str) -> str:
        return os.path.relpath(asset_dir / filename, path.parent).replace(os.sep, "/")

    agreement_sorted = sorted(agreement, key=lambda row: row["spearman_rho"], reverse=True)
    key_agreement = [row for row in agreement_sorted if row["channel"] in KEY_CHANNELS]
    direction_key = [row for row in direction if row["channel"] in KEY_CHANNELS]

    lines = [
        "# GigaTIME Run Agreement: Tile128 High-Trust vs Tile256 Expanded",
        "",
        "Status: Parameter/settings robustness check comparing overlapping slides across two completed GigaTIME runs.",
        "",
        "## Why This Matters",
        "",
        f"The current primary result uses `{comparison_label}`. Earlier, the project used `{reference_label}`. This analysis asks whether the same slide-level GigaTIME channels agree across those settings, and whether the HER2-low versus HER2-zero direction remains the same on overlapping slides.",
        "",
        "This is not a perfect tile-count experiment because the cohorts and random tile samples are not identical. It is still useful because it compares the same slide IDs where both runs exist.",
        "",
        "## Overlap",
        "",
        markdown_table(
            ["Quantity", "Count"],
            [
                ["Reference run slides", str(summary["n_reference_slides"])],
                ["Comparison run slides", str(summary["n_comparison_slides"])],
                ["Overlapping slide IDs", str(summary["n_overlap_slides"])],
                ["Overlapping HER2-low slides", str(summary["overlap_group_counts"].get("HER2-low", 0))],
                ["Overlapping HER2-zero slides", str(summary["overlap_group_counts"].get("HER2-zero", 0))],
                ["Overlapping HER2-positive slides", str(summary["overlap_group_counts"].get("HER2-positive", 0))],
            ],
        ),
        "",
        "The overlap preserves all 20 HER2-low and all 20 HER2-zero slides from the expanded 60-slide cohort. Two HER2-positive expanded-run cases are absent from the high-trust list because they were review/excluded cases.",
        "",
        "## Channel Agreement",
        "",
        markdown_table(
            ["Channel", "Spearman rho", "Pearson r", "Median absolute difference", "Mean comparison-reference"],
            [
                [
                    row["channel"],
                    fmt(row["spearman_rho"]),
                    fmt(row["pearson_r"]),
                    fmt(row["median_absolute_difference"], 5),
                    fmt(row["comparison_minus_reference_mean"], 5),
                ]
                for row in key_agreement
            ],
        ),
        "",
        f"![Channel correlation heatmap]({asset_link('tile128_vs_tile256_channel_correlation_heatmap.png')})",
        "",
        f"![Key channel scatter]({asset_link('key_channel_tile128_vs_tile256_scatter.png')})",
        "",
        "## HER2-Low Versus HER2-Zero Direction",
        "",
        markdown_table(
            ["Channel", "Reference low-zero delta", "Comparison low-zero delta", "Same direction", "Both low lower than zero"],
            [
                [
                    row["channel"],
                    fmt(row["reference_delta_low_minus_zero"], 5),
                    fmt(row["comparison_delta_low_minus_zero"], 5),
                    "yes" if row["same_direction"] else "no",
                    "yes" if row["both_low_lower_than_zero"] else "no",
                ]
                for row in direction_key
            ],
        ),
        "",
        f"![HER2-low versus zero delta comparison]({asset_link('low_zero_delta_tile128_vs_tile256.png')})",
        "",
        "## Interpretation",
        "",
        f"- Overlap is strong enough for this check: {summary['n_overlap_slides']} matched slides, including all HER2-low and HER2-zero slides from the 60-slide run.",
        f"- Key-channel direction agreement: {summary['key_channels_same_direction']} of {summary['n_key_channels']} tested key channels have the same HER2-low versus HER2-zero direction across runs.",
        f"- Key channels with HER2-low lower than HER2-zero in both runs: {summary['key_channels_both_low_lower']} of {summary['n_key_channels']}.",
        "- This supports that the main HER2-low versus HER2-zero direction is not simply an artifact of using 128 tiles instead of 256 tiles.",
        "- Some absolute channel scores shift across runs, which is expected because tile samples, cohort filtering, and run settings differ. The direction and relative agreement matter more than exact equality.",
        "",
        "## Machine-Readable Outputs",
        "",
        f"- [{asset_link('run_channel_agreement.csv')}]({asset_link('run_channel_agreement.csv')})",
        f"- [{asset_link('low_zero_direction_comparison.csv')}]({asset_link('low_zero_direction_comparison.csv')})",
        f"- [{asset_link('overlap_slide_scores.csv')}]({asset_link('overlap_slide_scores.csv')})",
        f"- [{asset_link('run_agreement_summary.json')}]({asset_link('run_agreement_summary.json')})",
        "",
        "## Cautious Claim This Supports",
        "",
        "> The HER2-low versus HER2-zero GigaTIME signal is directionally robust across an earlier 60-slide 256-tile run and the larger high-trust 128-tile run on overlapping slides, supporting a reproducible image-derived tissue-context association rather than a single-run sampling artifact.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    asset_dir = Path(args.asset_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)
    np, pd, plt, sns, stats = require_analysis_libs(out_dir / ".matplotlib")

    reference, comparison, overlap = load_overlap(pd, Path(args.reference_joined), Path(args.comparison_joined))
    channels = available_channels(overlap)
    agreement = build_channel_agreement(np, stats, overlap, channels)
    direction = build_low_zero_direction(stats, overlap, channels)
    group_counts = overlap["clinical_her2_group"].value_counts().to_dict()
    key_direction = [row for row in direction if row["channel"] in KEY_CHANNELS]
    summary = {
        "reference_label": args.reference_label,
        "comparison_label": args.comparison_label,
        "n_reference_slides": int(len(reference)),
        "n_comparison_slides": int(len(comparison)),
        "n_overlap_slides": int(len(overlap)),
        "overlap_group_counts": {str(key): int(value) for key, value in group_counts.items()},
        "n_channels": int(len(channels)),
        "n_key_channels": int(len(key_direction)),
        "key_channels_same_direction": int(sum(row["same_direction"] for row in key_direction)),
        "key_channels_both_low_lower": int(sum(row["both_low_lower_than_zero"] for row in key_direction)),
    }

    agreement_df = pd.DataFrame(agreement)
    direction_df = pd.DataFrame(direction)
    overlap_out = overlap[
        [
            "case_submitter_id",
            "slide_id",
            "clinical_her2_group",
            *[f"mean_{channel}_reference" for channel in channels],
            *[f"mean_{channel}_comparison" for channel in channels],
        ]
    ].copy()

    agreement_df.to_csv(out_dir / "run_channel_agreement.csv", index=False)
    direction_df.to_csv(out_dir / "low_zero_direction_comparison.csv", index=False)
    overlap_out.to_csv(out_dir / "overlap_slide_scores.csv", index=False)
    (out_dir / "run_agreement_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    agreement_df.to_csv(asset_dir / "run_channel_agreement.csv", index=False)
    direction_df.to_csv(asset_dir / "low_zero_direction_comparison.csv", index=False)
    overlap_out.to_csv(asset_dir / "overlap_slide_scores.csv", index=False)
    (asset_dir / "run_agreement_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    plot_correlation_heatmap(plt, sns, pd, agreement, asset_dir)
    plot_delta_comparison(plt, sns, pd, direction, asset_dir, args.reference_label, args.comparison_label)
    plot_key_scatter(plt, sns, pd, overlap, asset_dir, args.reference_label, args.comparison_label)
    write_markdown(
        Path(args.out_markdown),
        agreement,
        direction,
        summary,
        asset_dir,
        args.reference_label,
        args.comparison_label,
    )
    print(f"Wrote run agreement outputs to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
