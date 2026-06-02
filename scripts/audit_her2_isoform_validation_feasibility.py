#!/usr/bin/env python3
"""Audit whether the current local TCGA/GigaTIME data can validate HER2 isoform biology."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


BASE_RESULT_DIR = Path("results/gigatime_tcga_brca_clinical_her2_high_trust_tile128")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--high-trust-slides",
        default="docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv",
    )
    parser.add_argument("--erbb2-expression", default="data/tcga_brca/erbb2_expression.csv")
    parser.add_argument("--star-files", default="data/tcga_brca/tcga_brca_star_counts_files.csv")
    parser.add_argument("--expression-dir", default="data/tcga_brca/expression_files")
    parser.add_argument(
        "--out-dir",
        default=str(BASE_RESULT_DIR / "her2_isoform_validation_feasibility"),
    )
    parser.add_argument(
        "--out-markdown",
        default="docs/her2_isoform_validation_feasibility.md",
    )
    return parser.parse_args()


def read_rows(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def case_group_counts(rows: list[dict[str, str]], cases: set[str]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for row in rows:
        case_id = row.get("case_submitter_id", "")
        if case_id in cases:
            counts[row.get("clinical_her2_group", "unknown")] += 1
    return dict(counts)


def inspect_expression_file(path: Path) -> dict[str, object]:
    if not path.exists():
        return {
            "exists": False,
            "gene_model": "",
            "columns": "",
            "has_gene_id": False,
            "has_gene_name": False,
            "has_tpm": False,
            "has_transcript_id": False,
            "has_junction_columns": False,
        }
    gene_model = ""
    columns: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("# gene-model:"):
                gene_model = line.replace("# gene-model:", "").strip()
                continue
            if line and not line.startswith("#"):
                columns = line.split("\t")
                break
    lower_columns = {column.lower() for column in columns}
    return {
        "exists": True,
        "gene_model": gene_model,
        "columns": ", ".join(columns),
        "has_gene_id": "gene_id" in lower_columns,
        "has_gene_name": "gene_name" in lower_columns,
        "has_tpm": "tpm_unstranded" in lower_columns,
        "has_transcript_id": "transcript_id" in lower_columns,
        "has_junction_columns": any("junction" in column for column in lower_columns),
    }


def build_audit(args: argparse.Namespace) -> tuple[dict[str, object], list[dict[str, object]]]:
    high_trust = read_rows(args.high_trust_slides)
    erbb2_rows = read_rows(args.erbb2_expression)
    star_rows = read_rows(args.star_files)
    expression_dir = Path(args.expression_dir)

    high_trust_cases = {row["case_submitter_id"] for row in high_trust}
    erbb2_cases = {row["case_submitter_id"] for row in erbb2_rows}
    star_cases = {row["case_submitter_id"] for row in star_rows}
    local_expression_files = sorted(expression_dir.glob("*/*.rna_seq.augmented_star_gene_counts.tsv"))
    local_expression_cases = {path.parent.name for path in local_expression_files}

    missing_patterns = {
        "bam": list(Path("data/tcga_brca").glob("**/*.bam")),
        "fastq": list(Path("data/tcga_brca").glob("**/*.fastq*")),
        "junction": list(Path("data/tcga_brca").glob("**/*junction*")),
        "isoform": list(Path("data/tcga_brca").glob("**/*isoform*")),
        "transcript": list(Path("data/tcga_brca").glob("**/*transcript*")),
    }
    example_expression = local_expression_files[0] if local_expression_files else Path()
    expression_schema = inspect_expression_file(example_expression)

    high_trust_expression_cases = high_trust_cases & local_expression_cases
    high_trust_erbb2_cases = high_trust_cases & erbb2_cases
    low_zero_cases = {
        row["case_submitter_id"]
        for row in high_trust
        if row.get("clinical_her2_group") in {"HER2-low", "HER2-zero"}
    }
    low_zero_expression_cases = low_zero_cases & local_expression_cases

    audit = {
        "n_high_trust_slides": len(high_trust),
        "n_high_trust_cases": len(high_trust_cases),
        "n_erbb2_expression_cases": len(erbb2_cases),
        "n_star_manifest_cases": len(star_cases),
        "n_local_star_gene_count_files": len(local_expression_files),
        "n_local_star_gene_count_cases": len(local_expression_cases),
        "n_high_trust_cases_with_local_star_gene_counts": len(high_trust_expression_cases),
        "n_high_trust_cases_with_erbb2_expression_table": len(high_trust_erbb2_cases),
        "n_low_zero_high_trust_cases": len(low_zero_cases),
        "n_low_zero_high_trust_cases_with_local_star_gene_counts": len(low_zero_expression_cases),
        "high_trust_expression_group_counts": case_group_counts(high_trust, high_trust_expression_cases),
        "low_zero_expression_group_counts": case_group_counts(high_trust, low_zero_expression_cases),
        "local_bam_files": len(missing_patterns["bam"]),
        "local_fastq_files": len(missing_patterns["fastq"]),
        "local_junction_files": len(missing_patterns["junction"]),
        "local_isoform_files": len(missing_patterns["isoform"]),
        "local_transcript_files": len(missing_patterns["transcript"]),
        "example_expression_file": str(example_expression),
        "expression_schema": expression_schema,
    }

    feasibility_rows = [
        {
            "analysis": "Gene-level ERBB2 expression check",
            "status": "available_now",
            "local_evidence": f"{len(erbb2_cases)} ERBB2 expression cases; {len(high_trust_erbb2_cases)} overlap strict high-trust cases",
            "what_it_can_support": "Sanity check that clinical HER2 labels broadly track ERBB2 RNA expression.",
            "what_it_cannot_support": "Cannot identify HER2 isoforms or antibody-binding-domain loss.",
            "next_action": "Use only as context, not as isoform validation.",
        },
        {
            "analysis": "Bulk RNA immune/stromal program validation",
            "status": "available_now_limited",
            "local_evidence": f"{len(local_expression_cases)} local STAR gene-count cases; {len(low_zero_expression_cases)} strict high-trust low/zero overlap cases",
            "what_it_can_support": "Indirect validation of immune, epithelial, stromal, endothelial, and proliferation programs.",
            "what_it_cannot_support": "Cannot prove GigaTIME virtual mIF channels are true spatial protein measurements.",
            "next_action": "Continue treating as indirect support only; adjust for confounders.",
        },
        {
            "analysis": "HER2 isoform quantification like Guardia et al.",
            "status": "not_available_from_current_local_files",
            "local_evidence": "Local files are GDC STAR augmented gene-count TSVs with gene_id/gene_name/TPM columns and no transcript_id column.",
            "what_it_can_support": "Nothing directly; gene-level STAR counts cannot estimate ERBB2 isoform proportions.",
            "what_it_cannot_support": "Cannot compute kallisto transcript TPM, SUPPA2 PSI, p95/Delta16 isoform states, or antibody-binding-domain loss.",
            "next_action": "Requires short-read RNA-seq reads or an existing transcript-level isoform matrix from the paper/authors.",
        },
        {
            "analysis": "rMATS junction confirmation like Guardia et al.",
            "status": "not_available_from_current_local_files",
            "local_evidence": f"{audit['local_bam_files']} local BAM, {audit['local_fastq_files']} local FASTQ, {audit['local_junction_files']} local junction files under data/tcga_brca",
            "what_it_can_support": "Nothing directly right now.",
            "what_it_cannot_support": "Cannot confirm exon inclusion/exclusion events from junction counts.",
            "next_action": "Would require controlled-access BAM/FASTQ or author-provided junction/count matrices.",
        },
        {
            "analysis": "GigaTIME association with HER2 isoform/state labels",
            "status": "blocked_until_isoform_labels",
            "local_evidence": "GigaTIME slide features and clinical HER2 labels exist; isoform labels do not.",
            "what_it_can_support": "Once labels exist, can test whether image features associate with isoform/state groups.",
            "what_it_cannot_support": "Cannot claim image AI detects isoforms without molecular labels.",
            "next_action": "Ask Galante/Guardia group for TCGA sample-level HER2 isoform group/PSI table or reproduce their RNA-seq pipeline on appropriate reads.",
        },
    ]
    return audit, feasibility_rows


def write_markdown(path: Path, out_dir: Path, audit: dict[str, object], feasibility_rows: list[dict[str, object]]) -> None:
    schema = audit["expression_schema"]
    status_table = [
        [
            row["analysis"],
            row["status"],
            row["what_it_can_support"],
            row["next_action"],
        ]
        for row in feasibility_rows
    ]
    lines = [
        "# HER2 Isoform Validation Feasibility Audit",
        "",
        "Status: current feasibility audit after correcting the target paper to Guardia et al., Genome Research 2025, PMID 40664477.",
        "",
        "## Bottom Line",
        "",
        "The current local project can support gene-level ERBB2 and bulk RNA program checks, but it cannot directly validate HER2 isoforms. The local RNA files are GDC STAR augmented gene-count TSVs. They contain gene-level counts/TPM, not transcript-level isoform quantification, raw reads, BAMs, or junction-count outputs.",
        "",
        "Therefore, we can say GigaTIME features associate with HER2-low/HER2-zero state hypotheses, but we cannot say GigaTIME detects HER2 isoforms unless we obtain or generate isoform labels from appropriate RNA-seq data.",
        "",
        "## Correct Paper Requirement",
        "",
        "Guardia et al. used unprocessed TCGA short-read RNA-seq data, quantified transcripts with kallisto against an expanded HER2 transcriptome, computed isoform PSI with SUPPA2, and confirmed splicing events with rMATS junction counts. That is a different data level from our current local GDC STAR gene-count TSV files.",
        "",
        "## Local Data Audit",
        "",
        markdown_table(
            ["Item", "Count / value"],
            [
                ["Strict high-trust slides", audit["n_high_trust_slides"]],
                ["Strict high-trust cases", audit["n_high_trust_cases"]],
                ["ERBB2 expression table cases", audit["n_erbb2_expression_cases"]],
                ["STAR manifest cases", audit["n_star_manifest_cases"]],
                ["Local STAR gene-count files", audit["n_local_star_gene_count_files"]],
                ["Local STAR gene-count cases", audit["n_local_star_gene_count_cases"]],
                ["High-trust cases with local STAR gene counts", audit["n_high_trust_cases_with_local_star_gene_counts"]],
                ["Low/zero high-trust cases with local STAR gene counts", audit["n_low_zero_high_trust_cases_with_local_star_gene_counts"]],
                ["Local BAM files under data/tcga_brca", audit["local_bam_files"]],
                ["Local FASTQ files under data/tcga_brca", audit["local_fastq_files"]],
                ["Local junction files under data/tcga_brca", audit["local_junction_files"]],
                ["Local isoform files under data/tcga_brca", audit["local_isoform_files"]],
            ],
        ),
        "",
        "## Example Local Expression File Schema",
        "",
        markdown_table(
            ["Field", "Value"],
            [
                ["Example file", audit["example_expression_file"]],
                ["Gene model", schema["gene_model"]],
                ["Has gene_id", schema["has_gene_id"]],
                ["Has gene_name", schema["has_gene_name"]],
                ["Has tpm_unstranded", schema["has_tpm"]],
                ["Has transcript_id", schema["has_transcript_id"]],
                ["Has junction columns", schema["has_junction_columns"]],
            ],
        ),
        "",
        "## Feasibility Table",
        "",
        markdown_table(
            ["Analysis", "Status", "What it can support", "Next action"],
            status_table,
        ),
        "",
        "## Practical Recommendation",
        "",
        "For the paper proposal, the strongest honest next move is not to claim isoform detection from H&E. Instead:",
        "",
        "1. Present the GigaTIME result as a HER2-low/HER2-zero tissue-context association with strong TCGA confounding caveats.",
        "2. Use gene-level ERBB2 and RNA programs only as indirect context.",
        "3. Request or reproduce sample-level HER2 isoform labels from the Guardia/Galante workflow if we want to test isoform biology.",
        "4. Once isoform labels exist, test whether GigaTIME slide features, tumor-rich tile features, or embeddings predict isoform/state groups under source-site/slide-size controls.",
        "5. Keep the language: image AI predicts or associates with HER2 isoform/state hypotheses; it does not detect HER2 isoforms from the current data.",
        "",
        "## Output Files",
        "",
        f"- `{path}`",
        f"- `{out_dir / 'her2_isoform_validation_feasibility_summary.json'}`",
        f"- `{out_dir / 'her2_isoform_validation_feasibility_table.csv'}`",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    audit, feasibility_rows = build_audit(args)
    summary_path = out_dir / "her2_isoform_validation_feasibility_summary.json"
    table_path = out_dir / "her2_isoform_validation_feasibility_table.csv"
    summary_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    write_csv(
        table_path,
        feasibility_rows,
        [
            "analysis",
            "status",
            "local_evidence",
            "what_it_can_support",
            "what_it_cannot_support",
            "next_action",
        ],
    )
    write_markdown(Path(args.out_markdown), out_dir, audit, feasibility_rows)
    print(f"Wrote {summary_path}")
    print(f"Wrote {table_path}")
    print(f"Wrote {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
