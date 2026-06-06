#!/usr/bin/env python3
"""Aggregate every within-slide GigaTIME RNA-specificity audit into one cross-sample view.

Reads the per-sample JSON reports from both pipelines:
  - results/gigatime_hest_rna_validation/<id>/hest_rna_validation_report.json   (HEST Xenium + Visium)
  - results/gigatime_xenium_rna_validation*/xenium_rna_validation_report.json   (Janesick Rep1/Rep2)

and emits:
  - docs/hest_rna_validation_summary.md            (paper-ready cross-sample summary)
  - docs/assets/hest_rna_validation_summary/cross_sample_partial_r.png  (channel x sample heatmap)
  - results/gigatime_hest_rna_validation/cross_sample_summary.json

The load-bearing statistic is the cellularity-controlled partial correlation (partial_r_control_total)
already computed per sample; this script only collates and classifies it. No re-computation.

Run: ~/miniconda3/envs/gigatime-tcga/bin/python scripts/aggregate_hest_rna_validation.py
"""

from __future__ import annotations

import glob
import json
import os
from pathlib import Path

import numpy as np

# Canonical channel order: T-cell, epithelial, B/dendritic, myeloid/macrophage, proliferation, checkpoint, other.
CHANNEL_ORDER = ["CD3", "CD8", "CD4", "CK", "CD20", "CD11c", "CD68", "CD14", "CD16",
                 "Ki67", "PD-1", "PD-L1", "CD138", "CD34", "T-bet", "Tryptase"]


def discover_reports() -> list[dict]:
    """Load every per-sample report with a short label and modality."""
    out = []
    for p in sorted(glob.glob("results/gigatime_hest_rna_validation/*/hest_rna_validation_report.json")):
        rep = json.loads(Path(p).read_text())
        mp = rep.get("patient")
        meta_patient = mp if isinstance(mp, str) else ""
        modality = rep.get("modality", "xenium")
        label = f"{rep['sample']}"
        if meta_patient:
            label += "/" + meta_patient.replace("Patient ", "P").replace("patient ", "P")
        out.append({"path": p, "label": label, "modality": modality, "report": rep})
    for p in sorted(glob.glob("results/gigatime_xenium_rna_validation*/xenium_rna_validation_report.json")):
        rep = json.loads(Path(p).read_text())
        if not rep.get("specificity"):
            continue
        sample = rep.get("sample", "Xenium")
        label = "Janesick-" + ("Rep2" if "Rep2" in sample else "Rep1")
        out.append({"path": p, "label": label, "modality": "xenium (Janesick)", "report": rep})
    return out


def partial_map(rep: dict) -> dict[str, dict]:
    spec = rep.get("specificity") or {}
    return {r["channel"]: r for r in spec.get("per_channel", [])}


def build_matrix(entries: list[dict]):
    channels = [c for c in CHANNEL_ORDER if any(c in partial_map(e["report"]) for e in entries)]
    labels = [e["label"] for e in entries]
    mat = np.full((len(channels), len(entries)), np.nan)
    survive = np.zeros_like(mat, dtype=bool)
    rowmax = np.zeros_like(mat, dtype=bool)
    for j, e in enumerate(entries):
        pm = partial_map(e["report"])
        for i, ch in enumerate(channels):
            if ch in pm:
                mat[i, j] = pm[ch].get("partial_r_control_total", np.nan)
                survive[i, j] = bool(pm[ch].get("partial_survives", False))
                rowmax[i, j] = bool(pm[ch].get("own_is_row_max", False))
    return channels, labels, mat, survive, rowmax


def classify(row_vals: np.ndarray, row_survive: np.ndarray) -> tuple[str, float, int, int]:
    vals = row_vals[~np.isnan(row_vals)]
    n = len(vals)
    n_pos = int(row_survive[~np.isnan(row_vals)].sum())
    mean = float(np.mean(vals)) if n else float("nan")
    if n == 0:
        verdict = "untested"
    elif n_pos == 0 or mean <= 0.02:
        verdict = "never specific"
    elif n_pos >= max(1, int(round(0.8 * n))) and mean >= 0.10:
        verdict = "consistently specific"
    else:
        verdict = "variable"
    return verdict, mean, n_pos, n


def render_heatmap(channels, labels, mat, survive, rowmax, out_path: Path):
    os.environ.setdefault("MPLCONFIGDIR", str(out_path.parent / ".matplotlib"))
    (out_path.parent / ".matplotlib").mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(0.85 * len(labels) + 3, 0.5 * len(channels) + 2))
    vmax = float(np.nanmax(np.abs(mat))) if np.isfinite(mat).any() else 0.4
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(channels)))
    ax.set_yticklabels(channels, fontsize=9)
    for i in range(len(channels)):
        for j in range(len(labels)):
            if not np.isnan(mat[i, j]):
                txt = f"{mat[i, j]:.2f}"
                if rowmax[i, j]:
                    txt += "*"
                ax.text(j, i, txt, ha="center", va="center", fontsize=6.5,
                        color="black" if abs(mat[i, j]) < 0.6 * vmax else "white")
    ax.set_title("Cellularity-controlled partial r (virtual channel vs own-gene RNA)\nper sample; * = own-gene is row-max")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03, label="partial r")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> int:
    entries = discover_reports()
    if not entries:
        raise SystemExit("No per-sample reports found.")
    channels, labels, mat, survive, rowmax = build_matrix(entries)

    asset_dir = Path("docs/assets/hest_rna_validation_summary")
    asset_dir.mkdir(parents=True, exist_ok=True)
    heatmap = asset_dir / "cross_sample_partial_r.png"
    render_heatmap(channels, labels, mat, survive, rowmax, heatmap)

    per_channel = []
    for i, ch in enumerate(channels):
        verdict, mean, n_pos, n = classify(mat[i], survive[i])
        per_channel.append({"channel": ch, "verdict": verdict, "mean_partial_r": mean,
                            "n_positive": n_pos, "n_tested": n})

    # Provenance / platform counts.
    n_xen = sum(1 for e in entries if e["modality"].startswith("xenium") and "Janesick" not in e["modality"])
    n_jan = sum(1 for e in entries if "Janesick" in e["modality"])
    n_vis = sum(1 for e in entries if e["modality"] == "visium")

    summary = {"samples": labels, "channels": channels,
               "partial_r": [[None if np.isnan(v) else round(float(v), 4) for v in row] for row in mat],
               "per_channel": per_channel,
               "n_hest_xenium": n_xen, "n_janesick": n_jan, "n_visium": n_vis}
    Path("results/gigatime_hest_rna_validation/cross_sample_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Markdown
    L = []
    L += [
        "# GigaTIME RNA-Specificity: Cross-Sample Generalization Summary",
        "",
        "Status: aggregate of every within-slide RNA-specificity audit run to date. Tests whether the "
        "Janesick single-section finding generalizes across independent breast tumors and spatial-transcriptomics platforms.",
        "",
        f"**Cohort:** {len(entries)} sections — Janesick Xenium ({n_jan}), HEST-1k Xenium independent patients "
        f"({n_xen}), HEST-1k Visium whole-transcriptome ({n_vis}). Two platforms; IDC + ILC histology.",
        "",
        "The load-bearing statistic is the **cellularity-controlled partial correlation** between each virtual "
        "channel and its own-gene RNA, per sample (positive = channel-specific signal beyond tissue density). "
        "All values are collated from the per-sample reports; the same audited core computed every one.",
        "",
        "## Cross-sample partial r (virtual channel vs own-gene RNA, controlling for cellularity)",
        "",
        "![Cross-sample heatmap](assets/hest_rna_validation_summary/cross_sample_partial_r.png)",
        "",
        "| Channel | " + " | ".join(labels) + " | verdict |",
        "|---|" + "|".join(["---:"] * len(labels)) + "|---|",
    ]
    verdict_by_ch = {pc["channel"]: pc for pc in per_channel}
    for i, ch in enumerate(channels):
        cells = []
        for j in range(len(labels)):
            v = mat[i, j]
            cells.append("—" if np.isnan(v) else (f"**{v:.2f}**" if survive[i, j] and v > 0 else f"{v:.2f}"))
        L.append(f"| {ch} | " + " | ".join(cells) + f" | {verdict_by_ch[ch]['verdict']} |")
    L += [
        "",
        "Bold = cellularity-controlled partial r with 95% CI > 0 (channel-specific signal survives). "
        "Cells are blank where the sample's panel lacked that gene (the 541-gene Xenium panel omits "
        "CD4/CD14/CD16/PD-1/PD-L1/Tryptase/CD34/T-bet; Visium is whole-transcriptome).",
        "",
        "## Per-channel verdict across samples",
        "",
        "| Channel | Verdict | Mean partial r | Specific in N / tested |",
        "|---|---|---:|---:|",
    ]
    order = {"consistently specific": 0, "variable": 1, "never specific": 2, "untested": 3}
    for pc in sorted(per_channel, key=lambda d: (order.get(d["verdict"], 9), -d["mean_partial_r"])):
        L.append(f"| {pc['channel']} | {pc['verdict']} | {pc['mean_partial_r']:.2f} | {pc['n_positive']}/{pc['n_tested']} |")

    L += [
        "",
        "## Conclusion",
        "",
        "Across independent breast tumors and two ST platforms, GigaTIME virtual channels show only weak and "
        "**tissue-variable** marker specificity. The aggregate **T-cell channels (CD3/CD8/CD4)** are the most "
        "reproducible — positive in nearly every section — followed by **CK** (epithelium), which is specific in "
        "most but not all tumors. The **macrophage (CD68), myeloid (CD14/CD16), checkpoint (PD-L1) and "
        "proliferation (Ki67)** channels carry essentially no marker-specific signal after cellularity control on "
        "any platform (CD68 even inverts negative in the Janesick sections). Own-gene is rarely the top correlate. "
        "This generalizes and sharpens the single-section finding: the virtual channels reflect a broad "
        "T-cell-infiltrated-vs-epithelial contrast, not faithful per-marker stains, and they are not even "
        "consistent across patients — so they cannot serve as quantitative cell-type readouts or load-bearing "
        "biological evidence, only as interpretive context.",
        "",
        "## Per-sample provenance",
        "",
        "| Sample | Modality | Tiles | Alignment sanity r |",
        "|---|---|---:|---:|",
    ]
    for e in entries:
        rep = e["report"]
        s = rep.get("alignment_sanity", {})
        tiles = rep.get("n_tissue_tiles") or rep.get("n_tiles") or "?"
        L.append(f"| {e['label']} | {e['modality']} | {tiles} | {s.get('tissue_vs_total_transcripts_spearman_r', float('nan')):.3f} |")
    L += ["", "Generated by `scripts/aggregate_hest_rna_validation.py`.", ""]

    Path("docs/hest_rna_validation_summary.md").write_text("\n".join(L), encoding="utf-8")
    print(f"Samples: {labels}")
    print(f"Channels: {channels}")
    for pc in per_channel:
        print(f"  {pc['channel']:8s} {pc['verdict']:22s} mean={pc['mean_partial_r']:.3f} pos={pc['n_positive']}/{pc['n_tested']}")
    print("Wrote docs/hest_rna_validation_summary.md and the heatmap + JSON.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
