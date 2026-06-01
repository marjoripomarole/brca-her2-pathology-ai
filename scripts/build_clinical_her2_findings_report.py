#!/usr/bin/env python3
"""Build a simple notebook and HTML report for the clinical HER2 findings."""

from __future__ import annotations

import argparse
import csv
import html
import json
import shutil
from pathlib import Path


ASSET_COPIES = [
    (
        "results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_channel_boxplots.png",
        "docs/assets/clinical_her2_findings/clinical_her2_channel_boxplots.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_group_mean_heatmap.png",
        "docs/assets/clinical_her2_findings/clinical_her2_group_mean_heatmap.png",
    ),
    (
        "results/gigatime_tcga_brca_clinical_her2/clinical_summary/erbb2_tpm_by_clinical_her2_group.png",
        "docs/assets/clinical_her2_findings/erbb2_tpm_by_clinical_her2_group.png",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-notebook", default="notebooks/clinical_her2_findings_simple.ipynb")
    parser.add_argument("--out-html", default="notebooks/clinical_her2_findings_simple.html")
    parser.add_argument(
        "--channel-summary",
        default="results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_channel_summary.csv",
    )
    parser.add_argument(
        "--pairwise-tests",
        default="results/gigatime_tcga_brca_clinical_her2/clinical_summary/clinical_her2_pairwise_tests.csv",
    )
    parser.add_argument(
        "--rna-correlations",
        default="results/gigatime_tcga_brca_clinical_her2/rna_validation/gigatime_rna_signature_correlations.csv",
    )
    parser.add_argument(
        "--visual-qc-cases",
        default="docs/assets/clinical_her2_visual_qc/clinical_her2_visual_qc_selected_cases.csv",
    )
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def fmt(value: float, digits: int = 3) -> str:
    if value != value:
        return ""
    if abs(value) < 0.001 and value != 0:
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def html_table(headers: list[str], rows: list[list[str]]) -> str:
    header_html = "".join(f"<th>{html.escape(str(header))}</th>" for header in headers)
    body_rows = []
    for row in rows:
        body_rows.append("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>")
    return f"<table><thead><tr>{header_html}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def copy_assets() -> None:
    for source, destination in ASSET_COPIES:
        src = Path(source)
        dst = Path(destination)
        if not src.exists():
            raise FileNotFoundError(src)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def build_content(args: argparse.Namespace):
    copy_assets()
    channels = read_rows(Path(args.channel_summary))
    pairwise = read_rows(Path(args.pairwise_tests))
    rna = read_rows(Path(args.rna_correlations))
    qc_cases = read_rows(Path(args.visual_qc_cases))

    top_channels = sorted(channels, key=lambda row: as_float(row, "kruskal_p_value"))[:8]
    top_pairs = sorted(pairwise, key=lambda row: as_float(row, "mannwhitney_p_value"))[:6]
    rna_sorted = sorted(rna, key=lambda row: as_float(row, "spearman_rho"), reverse=True)

    channel_rows = [
        [
            row["channel"],
            fmt(as_float(row, "kruskal_p_value"), 4),
            fmt(as_float(row, "kruskal_q_value_bh"), 4),
            row["highest_mean_group"],
            row["lowest_mean_group"],
            fmt(as_float(row, "max_minus_min_mean"), 4),
        ]
        for row in top_channels
    ]
    pair_rows = [
        [
            row["channel"],
            f"{row['group_a']} vs {row['group_b']}",
            fmt(as_float(row, "delta_mean_a_minus_b"), 4),
            fmt(as_float(row, "mannwhitney_p_value"), 4),
            fmt(as_float(row, "mannwhitney_q_value_bh"), 4),
        ]
        for row in top_pairs
    ]
    rna_rows = [
        [
            row["channel"],
            fmt(as_float(row, "spearman_rho"), 3),
            fmt(as_float(row, "spearman_p_value"), 4),
            fmt(as_float(row, "spearman_q_value_bh"), 4),
        ]
        for row in rna_sorted
    ]
    qc_rows = [
        [
            row["clinical_her2_group"],
            row["case_submitter_id"],
            fmt(as_float(row, "qc_signal"), 3),
            fmt(as_float(row, "mean_CD68"), 3),
            fmt(as_float(row, "mean_PD-L1"), 3),
            fmt(as_float(row, "mean_CD11c"), 3),
        ]
        for row in qc_cases
    ]
    return {
        "channel_rows": channel_rows,
        "pair_rows": pair_rows,
        "rna_rows": rna_rows,
        "qc_rows": qc_rows,
    }


def notebook_cell(source: str) -> dict[str, object]:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": [line + "\n" for line in source.strip().splitlines()],
    }


def build_notebook(args: argparse.Namespace, content: dict[str, list[list[str]]]) -> None:
    nb_path = Path(args.out_notebook)
    nb_path.parent.mkdir(parents=True, exist_ok=True)
    cells = [
        notebook_cell(
            """
# Clinical HER2 GigaTIME Findings

**Simple display notebook**  
Updated clinical HER2 pilot summary for TCGA-BRCA.

**Core message:** we completed a balanced 30-slide pilot across HER2-positive, HER2-low, and HER2-zero cases. GigaTIME predicted higher immune/checkpoint-like signal in HER2-zero than HER2-low for several channels, but RNA validation was weak. The result is interesting and useful for a proposal, but it is not validated biology yet.
            """
        ),
        notebook_cell(
            """
## Study Design

- Dataset: TCGA-BRCA breast cancer.
- Images: diagnostic H&E whole-slide images.
- Model: released GigaTIME model.
- Groups: 10 HER2-positive, 10 HER2-low, and 10 HER2-zero cases.
- Sampling: 64 random tissue tiles per slide.
- Main outputs: virtual mIF channel scores, RNA validation, and visual QC panels.

**Plain-language translation:** we asked whether an AI model can see immune-like tissue patterns on ordinary H&E slides that differ across clinical HER2 groups.
            """
        ),
        notebook_cell(
            """
## Main Finding

The strongest pilot pattern was:

> HER2-zero had higher GigaTIME-predicted immune/checkpoint signals than HER2-low, especially CD68, PD-L1, and CD11c.

HER2-positive was usually between HER2-zero and HER2-low for these channels.

This should be described as a **hypothesis-generating pilot signal**, not a final biological claim.
            """
        ),
        notebook_cell(
            "## Top Three-Group Differences\n\n"
            + markdown_table(
                ["Channel", "Kruskal p", "BH q", "Highest group", "Lowest group", "Max-min mean"],
                content["channel_rows"],
            )
        ),
        notebook_cell(
            "## Group-Level Virtual mIF Summary\n\n"
            "This heatmap shows mean GigaTIME virtual-channel activation by clinical HER2 group.\n\n"
            "![Clinical HER2 group mean heatmap](../docs/assets/clinical_her2_findings/clinical_her2_group_mean_heatmap.png)"
        ),
        notebook_cell(
            "## Channel Distributions\n\n"
            "Each dot is one TCGA case. This is useful for seeing how noisy the pilot still is.\n\n"
            "![Clinical HER2 channel boxplots](../docs/assets/clinical_her2_findings/clinical_her2_channel_boxplots.png)"
        ),
        notebook_cell(
            "## Top Pairwise Tests\n\n"
            + markdown_table(
                ["Channel", "Comparison", "Delta mean", "Mann-Whitney p", "BH q"],
                content["pair_rows"],
            )
            + "\n\nThe strongest pairwise tests were mostly HER2-low versus HER2-zero, but they did not survive multiple-testing correction."
        ),
        notebook_cell(
            "## RNA Validation Check\n\n"
            "We compared GigaTIME virtual channels with matched RNA-seq marker signatures from the same 30 cases.\n\n"
            + markdown_table(["Channel", "Spearman rho", "p", "BH q"], content["rna_rows"])
            + "\n\n**Interpretation:** RNA validation did not strongly confirm the virtual immune-channel signal. Ki67 had the strongest positive trend, but no channel was FDR-significant."
        ),
        notebook_cell(
            "## RNA Correlation Heatmap\n\n"
            "![GigaTIME RNA correlation heatmap](../docs/assets/clinical_her2_rna_validation/gigatime_rna_correlation_heatmap.png)"
        ),
        notebook_cell(
            "## Visual QC Check\n\n"
            "We selected the top case from each HER2 group by combined CD68 + PD-L1 + CD11c virtual signal.\n\n"
            + markdown_table(["Group", "Case", "Combined", "CD68", "PD-L1", "CD11c"], content["qc_rows"])
            + "\n\nVisual QC showed that high-scoring tiles were tissue-containing and cellular, not obvious blank background. This supports continuing the analysis, but it does not validate the virtual markers."
        ),
        notebook_cell(
            "## Example H&E vs Virtual mIF Panels\n\n"
            "### HER2-zero top case\n\n"
            "![HER2-zero visual QC](../docs/assets/clinical_her2_visual_qc/her2_zero_TCGA-A2-A0T2_he_vs_virtual_mif_qc.png)\n\n"
            "### HER2-low top case\n\n"
            "![HER2-low visual QC](../docs/assets/clinical_her2_visual_qc/her2_low_TCGA-A2-A04Q_he_vs_virtual_mif_qc.png)\n\n"
            "### HER2-positive top case\n\n"
            "![HER2-positive visual QC](../docs/assets/clinical_her2_visual_qc/her2_positive_TCGA-A2-A0EQ_he_vs_virtual_mif_qc.png)"
        ),
        notebook_cell(
            """
## What We Can Say

- The full 30-slide clinical HER2 pilot is complete.
- The most interesting signal is HER2-zero > HER2-low for virtual CD68, PD-L1, and CD11c.
- Visual QC suggests the signal is not simply blank background.
- RNA validation did not strongly confirm the virtual immune-channel signal.

## What We Should Not Say Yet

- Do not say GigaTIME validated real mIF in TCGA.
- Do not say HER2-zero definitively has more immune infiltration.
- Do not say the model can classify HER2 status.
- Do not overinterpret p values from 10 cases per group.
            """
        ),
        notebook_cell(
            """
## Next Step

The clean next step is robustness:

1. Rerun the 30 selected slides with more tiles per slide, ideally 256 or 512.
2. Repeat the clinical HER2 summary.
3. Repeat RNA validation.
4. Repeat visual QC on the new top-driving cases.
5. Ask an advisor/pathologist whether the high-signal H&E regions look biologically plausible.

**One-sentence proposal framing:** GigaTIME produced a plausible but unvalidated virtual immune/checkpoint signal separating HER2-zero from HER2-low in a balanced TCGA-BRCA pilot, motivating deeper tile sampling and orthogonal validation.
            """
        ),
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


def section(title: str, body: str) -> str:
    return f"<section><h2>{html.escape(title)}</h2>{body}</section>"


def build_html(args: argparse.Namespace, content: dict[str, list[list[str]]]) -> None:
    html_path = Path(args.out_html)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    css = """
body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #1f2933; background: #f6f7f9; }
main { max-width: 1080px; margin: 0 auto; padding: 42px 28px 70px; }
.hero { background: #111827; color: white; padding: 34px; border-radius: 8px; }
.hero h1 { margin: 0 0 10px; font-size: 34px; letter-spacing: 0; }
.hero p { font-size: 18px; line-height: 1.5; max-width: 880px; }
section { background: white; margin-top: 18px; padding: 26px; border-radius: 8px; border: 1px solid #e5e7eb; }
h2 { margin: 0 0 14px; font-size: 24px; }
h3 { margin-top: 24px; }
p, li { font-size: 16px; line-height: 1.55; }
.callout { background: #eef6ff; border-left: 5px solid #2563eb; padding: 14px 16px; margin: 16px 0; }
.warning { background: #fff7ed; border-left: 5px solid #f97316; padding: 14px 16px; margin: 16px 0; }
table { width: 100%; border-collapse: collapse; margin: 14px 0; font-size: 14px; }
th, td { text-align: left; border-bottom: 1px solid #e5e7eb; padding: 8px 9px; vertical-align: top; }
th { background: #f3f4f6; font-weight: 650; }
img { display: block; width: 100%; max-width: 980px; height: auto; margin: 14px auto; border: 1px solid #e5e7eb; border-radius: 6px; background: white; }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 16px; }
.small { color: #52606d; font-size: 14px; }
"""
    parts = [
        "<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>Clinical HER2 GigaTIME Findings</title>",
        f"<style>{css}</style></head><body><main>",
        """
<div class="hero">
  <h1>Clinical HER2 GigaTIME Findings</h1>
  <p>Simple summary of the TCGA-BRCA pilot so far: 30 slides, three clinical HER2 groups, GigaTIME virtual mIF outputs, RNA validation, and visual QC.</p>
</div>
""",
        section(
            "Bottom Line",
            """
<div class="callout">
  <p><strong>Main signal:</strong> HER2-zero had higher GigaTIME-predicted immune/checkpoint-like signal than HER2-low, especially CD68, PD-L1, and CD11c.</p>
</div>
<div class="warning">
  <p><strong>Important caution:</strong> RNA validation was weak, and the pairwise tests did not survive multiple-testing correction. This is a proposal-ready hypothesis, not a validated biological claim.</p>
</div>
""",
        ),
        section(
            "Study Design",
            """
<ul>
  <li>Dataset: TCGA-BRCA breast cancer.</li>
  <li>Groups: 10 HER2-positive, 10 HER2-low, and 10 HER2-zero cases.</li>
  <li>Input: diagnostic H&E whole-slide images.</li>
  <li>Model: released GigaTIME virtual mIF model.</li>
  <li>Sampling: 64 random tissue tiles per slide.</li>
</ul>
""",
        ),
        section(
            "Top Three-Group Differences",
            html_table(
                ["Channel", "Kruskal p", "BH q", "Highest group", "Lowest group", "Max-min mean"],
                content["channel_rows"],
            ),
        ),
        section(
            "Group-Level Virtual mIF Summary",
            "<p>This heatmap summarizes mean virtual-channel activation by clinical HER2 group.</p>"
            "<img src='../docs/assets/clinical_her2_findings/clinical_her2_group_mean_heatmap.png' alt='Clinical HER2 group mean heatmap'>",
        ),
        section(
            "Channel Distributions",
            "<p>Each dot is one TCGA case. The spread shows why this remains a pilot.</p>"
            "<img src='../docs/assets/clinical_her2_findings/clinical_her2_channel_boxplots.png' alt='Clinical HER2 channel boxplots'>",
        ),
        section(
            "Top Pairwise Tests",
            html_table(["Channel", "Comparison", "Delta mean", "Mann-Whitney p", "BH q"], content["pair_rows"])
            + "<p class='small'>The strongest pairwise tests were mostly HER2-low versus HER2-zero, but none were FDR-significant.</p>",
        ),
        section(
            "RNA Validation",
            "<p>We compared virtual channels with matched RNA-seq marker signatures. No channel was FDR-significant.</p>"
            + html_table(["Channel", "Spearman rho", "p", "BH q"], content["rna_rows"])
            + "<img src='../docs/assets/clinical_her2_rna_validation/gigatime_rna_correlation_heatmap.png' alt='GigaTIME RNA correlation heatmap'>",
        ),
        section(
            "Visual QC",
            "<p>Top cases were selected by combined CD68 + PD-L1 + CD11c virtual signal.</p>"
            + html_table(["Group", "Case", "Combined", "CD68", "PD-L1", "CD11c"], content["qc_rows"])
            + "<p>The high-scoring tiles contain tissue and cells rather than obvious blank background. This supports follow-up, but it is not biological validation.</p>",
        ),
        section(
            "Example QC Panels",
            """
<div class="grid">
  <div><h3>HER2-zero</h3><img src="../docs/assets/clinical_her2_visual_qc/her2_zero_TCGA-A2-A0T2_he_vs_virtual_mif_qc.png" alt="HER2-zero visual QC"></div>
  <div><h3>HER2-low</h3><img src="../docs/assets/clinical_her2_visual_qc/her2_low_TCGA-A2-A04Q_he_vs_virtual_mif_qc.png" alt="HER2-low visual QC"></div>
  <div><h3>HER2-positive</h3><img src="../docs/assets/clinical_her2_visual_qc/her2_positive_TCGA-A2-A0EQ_he_vs_virtual_mif_qc.png" alt="HER2-positive visual QC"></div>
</div>
""",
        ),
        section(
            "What To Say",
            """
<ul>
  <li>We completed the balanced 30-case clinical HER2 pilot.</li>
  <li>The leading signal is HER2-zero greater than HER2-low for virtual CD68, PD-L1, and CD11c.</li>
  <li>Visual QC makes the signal look plausible enough to follow.</li>
  <li>RNA validation did not confirm the signal strongly, so the claim must stay cautious.</li>
</ul>
""",
        ),
        section(
            "Next Step",
            """
<p>Rerun the same 30 slides with denser tile sampling, such as 256 or 512 tiles per slide. Then repeat the clinical HER2 summary, RNA validation, and visual QC.</p>
<p><strong>Proposal framing:</strong> GigaTIME produced a plausible but unvalidated virtual immune/checkpoint signal separating HER2-zero from HER2-low in a balanced TCGA-BRCA pilot, motivating deeper sampling and orthogonal validation.</p>
""",
        ),
        "</main></body></html>",
    ]
    html_path.write_text("\n".join(parts), encoding="utf-8")


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
