#!/usr/bin/env python3
"""Render H&E and virtual mIF QC panels for case-driver review cases."""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from run_gigatime_tcga_brca import (
    GIGATIME_CHANNELS,
    import_runtime,
    load_model,
    preprocess_tile,
    resolve_device,
)


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")
FEATURE_VIEW = "all_sampled_tissue"

PANELS = {
    "immune_checkpoint": [
        ("DAPI", "#2D5BFF"),
        ("CK", "#FF3B30"),
        ("CD3", "#00E676"),
        ("CD4", "#30D158"),
        ("PD-L1", "#FF2DCE"),
        ("PD-1", "#FFD60A"),
    ],
    "myeloid_context": [
        ("DAPI", "#2D5BFF"),
        ("CK", "#FF3B30"),
        ("CD68", "#FF9500"),
        ("CD11c", "#64D2FF"),
        ("PD-L1", "#FF2DCE"),
        ("CD20", "#BF5AF2"),
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case-driver-scores",
        default=str(BASE_RESULT_DIR / "case_driver_analysis/case_driver_scores.csv"),
    )
    parser.add_argument(
        "--signal-channels",
        default=str(BASE_RESULT_DIR / "case_driver_analysis/low_zero_signal_channels_by_view.csv"),
    )
    parser.add_argument(
        "--classifier-review",
        default=str(BASE_RESULT_DIR / "case_driver_analysis/low_zero_classifier_review_cases.csv"),
    )
    parser.add_argument(
        "--tile-qc",
        default=str(BASE_RESULT_DIR / "gigatime_cleanup/tile_qc_scores.csv"),
    )
    parser.add_argument(
        "--high-trust-slides",
        default="docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv",
    )
    parser.add_argument(
        "--out-dir",
        default=str(BASE_RESULT_DIR / "case_driver_visual_qc"),
    )
    parser.add_argument(
        "--asset-dir",
        default="docs/assets/clinical_her2_high_trust_tile128_case_driver_visual_qc",
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/clinical_her2_high_trust_tile128_case_driver_visual_qc.md",
    )
    parser.add_argument("--gigatime-repo", default="external/GigaTIME")
    parser.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"])
    parser.add_argument("--tile-size", type=int, default=256)
    parser.add_argument("--tiles-per-case", type=int, default=4)
    parser.add_argument("--case-count-per-category", type=int, default=2)
    parser.add_argument("--manual-review-count", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--thumbnail-width", type=int, default=1500)
    return parser.parse_args()


def require_runtime(mpl_config_dir: Path):
    cache_dir = Path(tempfile.gettempdir()) / "gigatime_tcga_case_driver_visual_qc_mplconfig"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import pandas as pd

    return pd, plt, mpatches


def slug(value: str) -> str:
    return (
        value.lower()
        .replace(" ", "_")
        .replace("-", "_")
        .replace("/", "_")
        .replace("|", "_")
        .replace("+", "plus")
    )


def hex_to_rgb(hex_color: str) -> np.ndarray:
    clean = hex_color.lstrip("#")
    return np.array([int(clean[i : i + 2], 16) for i in (0, 2, 4)], dtype=np.float32) / 255.0


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


def load_inputs(pd, args: argparse.Namespace):
    driver_scores = pd.read_csv(args.case_driver_scores)
    signal_channels = pd.read_csv(args.signal_channels)
    classifier_review = pd.read_csv(args.classifier_review)
    tile_qc = pd.read_csv(args.tile_qc)
    high_trust = pd.read_csv(args.high_trust_slides)
    return driver_scores, signal_channels, classifier_review, tile_qc, high_trust


def selected_case_rows(pd, driver_scores, classifier_review, high_trust, args: argparse.Namespace):
    all_view = driver_scores.loc[driver_scores["feature_view"] == FEATURE_VIEW].copy()
    stability = driver_scores.copy()
    stability["opposite_profile_bool"] = stability["opposite_profile"].astype(str).str.lower().eq("true")
    stability = (
        stability.groupby(["slide_id", "case_submitter_id"], as_index=False)
        .agg(
            views_available=("feature_view", "nunique"),
            opposite_profile_views=("opposite_profile_bool", "sum"),
        )
    )
    stability["expected_profile_views"] = stability["views_available"] - stability["opposite_profile_views"]
    all_view = all_view.merge(
        stability,
        on=["slide_id", "case_submitter_id"],
        how="left",
        validate="one_to_one",
    )
    low = (
        all_view.loc[
            (all_view["clinical_her2_group"] == "HER2-low")
            & (all_view["expected_profile_views"] >= 3)
        ]
        .sort_values("zero_like_score", ascending=True)
        .head(args.case_count_per_category)
        .copy()
    )
    low["review_category"] = "label_consistent_her2_low"
    low["selection_reason"] = "Lowest HER2-zero-like score among stable HER2-low cases"

    zero = (
        all_view.loc[
            (all_view["clinical_her2_group"] == "HER2-zero")
            & (all_view["expected_profile_views"] >= 3)
        ]
        .sort_values("zero_like_score", ascending=False)
        .head(args.case_count_per_category)
        .copy()
    )
    zero["review_category"] = "label_consistent_her2_zero"
    zero["selection_reason"] = "Highest HER2-zero-like score among stable HER2-zero cases"

    review = classifier_review.head(args.manual_review_count).copy()
    review["review_category"] = "opposite_profile_manual_review"
    review["selection_reason"] = "Opposite profile and classifier-error priority case"
    if "zero_like_score" not in review.columns and "all_sampled_tissue_zero_like_score" in review.columns:
        review["zero_like_score"] = review["all_sampled_tissue_zero_like_score"]
    review = review.merge(
        all_view[
            [
                "slide_id",
                "case_submitter_id",
                "expected_profile_score",
                "her2_detail_subgroup",
            ]
        ],
        on=["slide_id", "case_submitter_id"],
        how="left",
        suffixes=("", "_driver"),
    )

    selected = pd.concat([low, zero, review], ignore_index=True)
    selected = selected.drop_duplicates("case_submitter_id", keep="first")
    selected = selected.merge(
        high_trust[
            [
                "slide_id",
                "case_submitter_id",
                "slide_local_path",
                "her2_ihc_score",
                "her2_ish_status",
                "er_status",
                "pr_status",
                "erbb2_tpm",
                "histological_type",
                "pathologic_stage",
                "processed_tissue_qc",
            ]
        ],
        on=["slide_id", "case_submitter_id"],
        how="left",
        suffixes=("", "_cohort"),
        validate="one_to_one",
    )
    for column in [
        "slide_local_path",
        "her2_ihc_score",
        "her2_ish_status",
        "er_status",
        "pr_status",
        "erbb2_tpm",
        "histological_type",
        "pathologic_stage",
        "processed_tissue_qc",
    ]:
        cohort_column = f"{column}_cohort"
        if cohort_column not in selected.columns:
            continue
        if column in selected.columns:
            selected[column] = selected[column].fillna(selected[cohort_column])
        else:
            selected[column] = selected[cohort_column]
    return selected


def signal_channel_rows(signal_channels):
    rows = signal_channels.loc[signal_channels["feature_view"] == FEATURE_VIEW].copy()
    rows = rows.sort_values("mannwhitney_q_value_bh_within_view")
    return rows


def add_tile_driver_scores(tile_rows, signals):
    tile_rows = tile_rows.copy()
    score_columns = []
    for _, signal in signals.iterrows():
        channel = signal["channel"]
        column = f"mean_{channel}"
        if column not in tile_rows.columns:
            continue
        values = tile_rows[column].astype(float)
        center = float(values.mean())
        scale = float(values.std(ddof=0))
        if scale == 0 or scale != scale:
            continue
        score_column = f"zero_like_tile_z_{slug(channel)}"
        tile_rows[score_column] = ((values - center) / scale) * float(signal["zero_higher_direction"])
        score_columns.append(score_column)
    if not score_columns:
        tile_rows["tile_zero_like_score"] = 0.0
    else:
        tile_rows["tile_zero_like_score"] = tile_rows[score_columns].mean(axis=1)
    return tile_rows


def select_tiles_for_case(tile_qc, case_row, signals, n_tiles: int):
    slide_tiles = tile_qc.loc[tile_qc["slide_id"] == case_row["slide_id"]].copy()
    if slide_tiles.empty:
        raise ValueError(f"No tile QC rows found for {case_row['slide_id']}")
    slide_tiles = add_tile_driver_scores(slide_tiles, signals)
    slide_tiles = slide_tiles.sort_values(["tissue_fraction", "tile_zero_like_score"], ascending=[False, False])
    try:
        case_zero_like_score = float(case_row["zero_like_score"])
    except (TypeError, ValueError):
        case_zero_like_score = float("nan")
    if case_zero_like_score == case_zero_like_score:
        select_zero_like_tiles = case_zero_like_score >= 0
    else:
        select_zero_like_tiles = case_row["clinical_her2_group"] == "HER2-zero"

    if select_zero_like_tiles:
        selected = slide_tiles.sort_values(["tile_zero_like_score", "tissue_fraction"], ascending=[False, False]).head(
            n_tiles
        )
    else:
        selected = slide_tiles.sort_values(["tile_zero_like_score", "tissue_fraction"], ascending=[True, False]).head(
            n_tiles
        )
    return selected.copy(), slide_tiles


def run_batch(torch, model, batch, meta):
    tensor = torch.stack(batch, dim=0)
    maps = torch.sigmoid(model(tensor)).detach().cpu().numpy()
    return [{**item, "maps": maps[idx]} for idx, item in enumerate(meta)]


def infer_tile_maps(torch, model, openslide, slide_path: Path, tile_rows, tile_size: int, batch_size: int, device):
    slide = openslide.OpenSlide(str(slide_path))
    records: list[dict[str, Any]] = []
    batch = []
    meta = []
    with torch.no_grad():
        for _idx, row in tile_rows.iterrows():
            x = int(row["x"])
            y = int(row["y"])
            region = slide.read_region((x, y), 0, (tile_size, tile_size)).convert("RGB")
            rgb = np.asarray(region)
            batch.append(preprocess_tile(torch, rgb, device))
            meta.append(
                {
                    "x": x,
                    "y": y,
                    "he_rgb": rgb,
                    "tile_zero_like_score": float(row["tile_zero_like_score"]),
                    "tissue_fraction": float(row["tissue_fraction"]),
                    "mean_CK": float(row.get("mean_CK", float("nan"))),
                    "mean_CD68": float(row.get("mean_CD68", float("nan"))),
                    "mean_PD-L1": float(row.get("mean_PD-L1", float("nan"))),
                    "mean_CD11c": float(row.get("mean_CD11c", float("nan"))),
                }
            )
            if len(batch) == batch_size:
                records.extend(run_batch(torch, model, batch, meta))
                batch = []
                meta = []
        if batch:
            records.extend(run_batch(torch, model, batch, meta))
    slide.close()
    return records


def panel_limits(records: list[dict[str, Any]], panel: list[tuple[str, str]]) -> dict[str, tuple[float, float]]:
    limits = {}
    for marker, _color in panel:
        channel_idx = GIGATIME_CHANNELS.index(marker)
        values = np.concatenate([record["maps"][channel_idx].ravel() for record in records])
        lower_quantile = 0.55 if marker == "DAPI" else 0.80
        upper_quantile = 0.995
        low = float(np.quantile(values, lower_quantile))
        high = float(np.quantile(values, upper_quantile))
        if high <= low:
            high = low + 1e-6
        limits[marker] = (low, high)
    return limits


def make_composite(channel_maps: np.ndarray, panel: list[tuple[str, str]], limits: dict[str, tuple[float, float]]) -> np.ndarray:
    height, width = channel_maps.shape[1:]
    composite = np.zeros((height, width, 3), dtype=np.float32)
    for marker, color in panel:
        channel_idx = GIGATIME_CHANNELS.index(marker)
        low, high = limits[marker]
        signal = np.clip((channel_maps[channel_idx] - low) / (high - low), 0, 1)
        signal = np.power(signal, 0.75)
        composite += signal[..., None] * hex_to_rgb(color)
    return np.clip(composite, 0, 1)


def save_case_panel(records, case_row, out_path: Path, plt, mpatches) -> None:
    limits_by_panel = {name: panel_limits(records, panel) for name, panel in PANELS.items()}
    nrows = len(records)
    fig, axes = plt.subplots(nrows=nrows, ncols=3, figsize=(9.6, max(7.2, nrows * 2.55)), facecolor="black")
    axes = np.atleast_2d(axes)
    for row_idx, record in enumerate(records):
        axes[row_idx, 0].imshow(record["he_rgb"])
        axes[row_idx, 0].set_title(
            "H&E\n"
            f"tile score {record['tile_zero_like_score']:.2f} | tissue {record['tissue_fraction']:.2f}\n"
            f"CK {record['mean_CK']:.3f} CD68 {record['mean_CD68']:.3f} PD-L1 {record['mean_PD-L1']:.3f}",
            color="white",
            fontsize=8,
        )
        for col_idx, (panel_name, panel) in enumerate(PANELS.items(), start=1):
            axes[row_idx, col_idx].imshow(make_composite(record["maps"], panel, limits_by_panel[panel_name]))
            axes[row_idx, col_idx].set_title(panel_name.replace("_", " "), color="white", fontsize=10)
        for axis in axes[row_idx]:
            axis.set_xticks([])
            axis.set_yticks([])
            axis.set_facecolor("black")
            for spine in axis.spines.values():
                spine.set_edgecolor("#333333")

    handles = []
    seen = set()
    for panel in PANELS.values():
        for marker, color in panel:
            if marker in seen:
                continue
            handles.append(mpatches.Patch(color=color, label=marker))
            seen.add(marker)
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 8), facecolor="black", labelcolor="white")
    title = (
        f"{case_row['case_submitter_id']} | {case_row['clinical_her2_group']} | "
        f"{case_row['review_category'].replace('_', ' ')}"
    )
    subtitle = (
        f"case zero-like score {fmt(case_row.get('zero_like_score'))} | "
        f"expected views {fmt(case_row.get('expected_profile_views'), 0)} | "
        f"opposite views {fmt(case_row.get('opposite_profile_views'), 0)}"
    )
    fig.suptitle(title + "\n" + subtitle, color="white", fontsize=14, y=0.992)
    fig.text(
        0.5,
        0.035,
        "Virtual mIF composites are GigaTIME predictions from H&E tiles, not experimental mIF.",
        color="white",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    fig.savefig(out_path, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


def scale_rect(x: float, y: float, tile_size: int, scale: float) -> tuple[int, int, int, int]:
    return (
        int(round(x * scale)),
        int(round(y * scale)),
        int(round((x + tile_size) * scale)),
        int(round((y + tile_size) * scale)),
    )


def save_tile_overlay(
    openslide,
    all_tiles,
    selected_tiles,
    slide_path: Path,
    case_row,
    tile_size: int,
    thumbnail_width: int,
    out_path: Path,
    plt,
) -> None:
    slide = openslide.OpenSlide(str(slide_path))
    slide_width, slide_height = slide.dimensions
    thumbnail_height = max(1, int(round(thumbnail_width * slide_height / slide_width)))
    thumbnail = slide.get_thumbnail((thumbnail_width, thumbnail_height)).convert("RGB")
    slide.close()

    scale = thumbnail.width / slide_width
    values = all_tiles["tile_zero_like_score"].astype(float)
    norm = plt.Normalize(
        vmin=float(values.quantile(0.05)),
        vmax=float(values.quantile(0.95)) if values.quantile(0.95) > values.quantile(0.05) else float(values.max()) + 1e-9,
    )
    cmap = plt.get_cmap("coolwarm")

    fig, ax = plt.subplots(figsize=(10.5, max(3.5, 10.5 * thumbnail.height / thumbnail.width)))
    ax.imshow(thumbnail)
    for _, tile in all_tiles.iterrows():
        x0, y0, x1, y1 = scale_rect(tile["x"], tile["y"], tile_size, scale)
        color = cmap(norm(float(tile["tile_zero_like_score"])))
        rect = plt.Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=color, edgecolor="none", alpha=0.28)
        ax.add_patch(rect)
    for _, tile in selected_tiles.iterrows():
        x0, y0, x1, y1 = scale_rect(tile["x"], tile["y"], tile_size, scale)
        rect = plt.Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="none", edgecolor="white", linewidth=1.7)
        ax.add_patch(rect)
    ax.set_title(
        f"{case_row['case_submitter_id']} | {case_row['clinical_her2_group']} | all sampled tile driver scores"
    )
    ax.axis("off")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label("HER2-zero-like tile score")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def write_markdown(path: Path, asset_dir: Path, manifest, selected_cases, signal_channels) -> None:
    overview_rows = []
    for _, row in selected_cases.iterrows():
        overview_rows.append(
            [
                row["case_submitter_id"],
                row["clinical_her2_group"],
                row["review_category"],
                row.get("her2_detail_subgroup", ""),
                fmt(row.get("zero_like_score")),
                fmt(row.get("expected_profile_views"), 0),
                fmt(row.get("opposite_profile_views"), 0),
                row.get("selection_reason", ""),
            ]
        )
    signal_rows = [
        [
            row["channel"],
            fmt(row["delta_low_minus_zero"], 4),
            fmt(row["mannwhitney_q_value_bh_within_view"], 4),
        ]
        for _, row in signal_channels.iterrows()
    ]
    lines = [
        "# Case-Driver Visual QC Panels",
        "",
        "This report renders a small visual QC set from the strict high-trust HER2-low versus HER2-zero case-driver analysis.",
        "",
        "The goal is not to validate the whole cohort visually. The goal is to inspect representative label-consistent cases and opposite-profile/manual-review cases, then decide whether the GigaTIME signal appears tissue-plausible or artifact-prone.",
        "",
        "Important: the fluorescence-style panels are GigaTIME virtual predictions from H&E tiles, not real multiplex immunofluorescence.",
        "",
        "## Cases Rendered",
        "",
        markdown_table(
            [
                "Case",
                "Group",
                "Review category",
                "HER2 detail",
                "Zero-like score",
                "Expected views",
                "Opposite views",
                "Selection reason",
            ],
            overview_rows,
        ),
        "",
        "## Signal Channels Used For Tile Selection",
        "",
        markdown_table(["Channel", "Low-zero delta", "BH q"], signal_rows),
        "",
        "## Selected Tile Summary",
        "",
        "This summary is the key visual-QC caveat. Low-like selected tiles can be tissue-containing but nearly blank for virtual CK, CD68, PD-L1, and CD11c, so they may represent stromal/collagen-rich tissue context rather than tumor-rich biology.",
        "",
        markdown_table(
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
            selected_tile_summary_rows(manifest.attrs["selected_tiles"]),
        ),
        "",
        "## Visual Panels",
        "",
    ]
    for _, row in manifest.iterrows():
        if row["image_type"] != "he_vs_virtual_mif_qc":
            continue
        overlay = manifest.loc[
            (manifest["case_submitter_id"] == row["case_submitter_id"]) & (manifest["image_type"] == "tile_overlay")
        ]
        lines.extend(
            [
                f"### {row['case_submitter_id']} | {row['clinical_her2_group']} | {row['review_category']}",
                "",
                f"![{row['case_submitter_id']} H&E and virtual mIF panel]({asset_link(asset_dir, Path(row['path']).name)})",
                "",
            ]
        )
        if not overlay.empty:
            lines.extend(
                [
                    f"![{row['case_submitter_id']} tile score overlay]({asset_link(asset_dir, Path(overlay.iloc[0]['path']).name)})",
                    "",
                ]
            )
    lines.extend(
        [
            "## How To Read This",
            "",
            "- Label-consistent HER2-low cases should generally show low HER2-zero-like tile/channel scores.",
            "- Label-consistent HER2-zero cases should generally show higher HER2-zero-like tile/channel scores.",
            "- Opposite-profile cases are the most important QC cases. They may reflect label noise, slide artifact, non-tumor tissue sampling, or real biological exceptions.",
            "- If the selected H&E tiles are blank, folded, necrotic, mostly stroma, or not tumor-rich, the case should be flagged before making biological claims.",
            "",
            "## Output Files",
            "",
            f"- `{path}`",
            f"- `{path.parent.parent / 'results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/case_driver_visual_qc/case_driver_visual_qc_manifest.csv'}`",
            f"- `{path.parent.parent / 'results/gigatime_tcga_brca_clinical_her2_high_trust_tile128/case_driver_visual_qc/case_driver_visual_qc_selected_tiles.csv'}`",
            f"- `{asset_dir}/`",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def selected_tile_summary_rows(selected_tiles) -> list[list[str]]:
    order = [
        ("label_consistent_her2_low", "HER2-low"),
        ("label_consistent_her2_zero", "HER2-zero"),
        ("opposite_profile_manual_review", "HER2-low"),
        ("opposite_profile_manual_review", "HER2-zero"),
    ]
    rows = []
    for review_category, clinical_group in order:
        subset = selected_tiles.loc[
            (selected_tiles["review_category"] == review_category)
            & (selected_tiles["clinical_her2_group"] == clinical_group)
        ]
        rows.append(
            [
                review_category.replace("_", " "),
                clinical_group,
                str(len(subset)),
                fmt(subset["tissue_fraction"].mean()),
                fmt(subset["tile_zero_like_score"].mean()),
                fmt(subset["mean_CK"].mean(), 4),
                fmt(subset["mean_CD68"].mean(), 4),
                fmt(subset["mean_PD-L1"].mean(), 4),
                fmt(subset["mean_CD11c"].mean(), 4),
            ]
        )
    return rows


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    asset_dir = Path(args.asset_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)

    pd, plt, mpatches = require_runtime(out_dir / ".matplotlib")
    driver_scores, signal_channels, classifier_review, tile_qc, high_trust = load_inputs(pd, args)
    selected_cases = selected_case_rows(pd, driver_scores, classifier_review, high_trust, args)
    signals = signal_channel_rows(signal_channels)

    torch, gigatime_class, snapshot_download, _image, openslide = import_runtime(Path(args.gigatime_repo))
    device = resolve_device(torch, args.device)
    model = load_model(torch, gigatime_class, snapshot_download, device)

    manifest_rows = []
    selected_tile_rows = []
    for _, case_row in selected_cases.iterrows():
        slide_path = Path(case_row["slide_local_path"])
        selected_tiles, all_tiles = select_tiles_for_case(tile_qc, case_row, signals, args.tiles_per_case)
        records = infer_tile_maps(torch, model, openslide, slide_path, selected_tiles, args.tile_size, args.batch_size, device)

        prefix = f"{slug(case_row['review_category'])}_{case_row['case_submitter_id']}"
        panel_path = asset_dir / f"{prefix}_he_vs_virtual_mif_qc.png"
        overlay_path = asset_dir / f"{prefix}_tile_driver_overlay.png"
        save_case_panel(records, case_row, panel_path, plt, mpatches)
        save_tile_overlay(
            openslide,
            all_tiles,
            selected_tiles,
            slide_path,
            case_row,
            args.tile_size,
            args.thumbnail_width,
            overlay_path,
            plt,
        )

        for image_type, path in [("he_vs_virtual_mif_qc", panel_path), ("tile_overlay", overlay_path)]:
            manifest_rows.append(
                {
                    "case_submitter_id": case_row["case_submitter_id"],
                    "slide_id": case_row["slide_id"],
                    "clinical_her2_group": case_row["clinical_her2_group"],
                    "review_category": case_row["review_category"],
                    "image_type": image_type,
                    "path": str(path),
                }
            )
        for _, tile in selected_tiles.iterrows():
            row = {
                "case_submitter_id": case_row["case_submitter_id"],
                "slide_id": case_row["slide_id"],
                "clinical_her2_group": case_row["clinical_her2_group"],
                "review_category": case_row["review_category"],
                "x": int(tile["x"]),
                "y": int(tile["y"]),
                "tissue_fraction": float(tile["tissue_fraction"]),
                "tile_zero_like_score": float(tile["tile_zero_like_score"]),
                "mean_CK": float(tile.get("mean_CK", float("nan"))),
                "mean_CD68": float(tile.get("mean_CD68", float("nan"))),
                "mean_PD-L1": float(tile.get("mean_PD-L1", float("nan"))),
                "mean_CD11c": float(tile.get("mean_CD11c", float("nan"))),
            }
            selected_tile_rows.append(row)

    manifest = pd.DataFrame(manifest_rows)
    selected_tiles = pd.DataFrame(selected_tile_rows)
    manifest.attrs["selected_tiles"] = selected_tiles
    manifest.to_csv(out_dir / "case_driver_visual_qc_manifest.csv", index=False)
    selected_cases.to_csv(out_dir / "case_driver_visual_qc_selected_cases.csv", index=False)
    selected_tiles.to_csv(out_dir / "case_driver_visual_qc_selected_tiles.csv", index=False)
    write_markdown(Path(args.out_markdown), asset_dir, manifest, selected_cases, signals)

    print(f"Wrote case-driver visual QC images to {asset_dir}")
    print(f"Wrote case-driver visual QC tables to {out_dir}")
    print(f"Wrote case-driver visual QC markdown to {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
