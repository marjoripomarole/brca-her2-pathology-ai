#!/usr/bin/env python3
"""Figure 1: study-design schematic for the cautionary-methods manuscript.

The apparent HER2-low-vs-zero H&E biomarker, the four falsification tests applied to it, and the
conclusion. Pure matplotlib (no data dependency); regenerate any time.

Run: ~/miniconda3/envs/gigatime-tcga/bin/python scripts/make_figure1_study_design.py
"""

from __future__ import annotations

import os
from pathlib import Path

OUT = Path("docs/assets/manuscript/figure1_study_design.png")


def main() -> int:
    os.environ.setdefault("MPLCONFIGDIR", str(OUT.parent / ".matplotlib"))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    (OUT.parent / ".matplotlib").mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    fig, ax = plt.subplots(figsize=(9.5, 11))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12.4)
    ax.axis("off")

    def box(cx, cy, w, h, text, face, edge, fontsize=9, weight="normal", textcolor="#111827"):
        ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                    boxstyle="round,pad=0.08,rounding_size=0.12",
                                    facecolor=face, edgecolor=edge, linewidth=1.6, zorder=2))
        ax.text(cx, cy, text, ha="center", va="center", fontsize=fontsize, color=textcolor,
                weight=weight, zorder=3, wrap=True)

    def arrow(x0, y0, x1, y1, color="#6b7280"):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=16,
                                     lw=1.6, color=color, zorder=1))

    # Header: the apparent biomarker.
    box(5, 11.4, 9.2, 1.5,
        "APPARENT BIOMARKER\nGigaTIME virtual-mIF channels separate HER2-low vs HER2-zero\n"
        "on TCGA-BRCA H&E (balanced accuracy 0.71; HER2-positive at chance)",
        face="#fef3c7", edge="#d97706", fontsize=10.5, weight="bold", textcolor="#7c2d12")

    ax.text(5, 9.95, "Four falsification tests", ha="center", va="center",
            fontsize=11, style="italic", color="#374151")

    tests = [
        (2.55, 8.35, "1 · Acquisition-confound baselines\n+ leave-source-site-out",
         "Slide size alone separates better (0.88) and stays\nportable under site holdout (0.88);\nGigaTIME collapses 0.71 → 0.62"),
        (7.45, 8.35, "2 · Generic foundation-model embeddings\n(H-Optimus-0, Virchow2)",
         "Reproduce the separation (0.73 / 0.69, p=0.005)\nand collapse identically under site holdout\n→ virtual-immune framing not required"),
        (2.55, 5.55, "3 · External single-scanner cohort\n(BCNB; grade/ER/PR/Ki67)",
         "Signal survives but is modest (AUC ~0.60–0.65),\ncomparable to & partly explained by\nclinical covariates"),
        (7.45, 5.55, "4 · Independent-modality specificity audit\n(Xenium/Visium RNA, 9 sections) + 2nd model",
         "Only CK + T-cell channels are marker-specific;\nCD68/myeloid/checkpoint fail. GigaTIME vs ROSIE\ndisagree on which channels work (r=0.12)"),
    ]
    for cx, cy, title, body in tests:
        ax.add_patch(FancyBboxPatch((cx - 2.25, cy - 1.05), 4.5, 2.1,
                                    boxstyle="round,pad=0.08,rounding_size=0.12",
                                    facecolor="#eff6ff", edgecolor="#2563eb", linewidth=1.6, zorder=2))
        ax.text(cx, cy + 0.62, title, ha="center", va="center", fontsize=9.0, weight="bold",
                color="#1e3a8a", zorder=3)
        ax.text(cx, cy - 0.28, body, ha="center", va="center", fontsize=8.2, color="#1f2937", zorder=3)

    # Arrows header -> tests, tests -> conclusion.
    for cx, cy, *_ in tests:
        arrow(5, 10.62, cx, cy + 1.12)
        arrow(cx, cy - 1.12, 5, 2.95)

    # Conclusion.
    box(5, 1.85, 9.2, 1.9,
        "CONCLUSION\nThe apparent biomarker = scanner/acquisition batch + tumor grade + immune context,\n"
        "NOT tumor-cell HER2 protein. H&E→virtual-mIF channels behave as a broad tissue-compartment\n"
        "contrast — interpretive context, not quantitative cell-type readouts.",
        face="#fee2e2", edge="#b91c1c", fontsize=10, weight="bold", textcolor="#7f1d1d")

    fig.tight_layout()
    fig.savefig(OUT, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
