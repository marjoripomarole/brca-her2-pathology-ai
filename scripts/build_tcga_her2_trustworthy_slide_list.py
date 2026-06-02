#!/usr/bin/env python3
"""Build a trustworthiness/QC list for TCGA-BRCA clinical HER2 slides."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


HER2_GROUPS = ["HER2-positive", "HER2-low", "HER2-zero"]
MISSING_MARKERS = {"", "[Not Available]", "[Not Evaluated]", "[Unknown]", "Unknown", "Not Reported"}


def require_pandas():
    try:
        import pandas as pd
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "Missing Python package: pandas. Use `conda activate gigatime-tcga` "
            "or `conda run -n gigatime-tcga ...`."
        ) from exc
    return pd


def clean(value) -> str:
    if value is None:
        return ""
    return str(value).strip()


def is_missing(value) -> bool:
    return clean(value) in MISSING_MARKERS


def bool_from_value(value) -> bool:
    if isinstance(value, bool):
        return value
    return clean(value).lower() in {"true", "1", "yes"}


def her2_detail_subgroup(row) -> str:
    group = clean(row.get("clinical_her2_group"))
    ihc = clean(row.get("her2_ihc_score"))
    ish = clean(row.get("her2_ish_status"))

    if group == "HER2-positive":
        if ihc == "3+" and ish == "Negative":
            return "HER2-positive_IHC3_ISH-negative_discordant"
        if ihc == "3+":
            return "HER2-positive_IHC3"
        if ish == "Positive" and ihc == "2+":
            return "HER2-positive_IHC2_ISH-positive"
        if ish == "Positive" and ihc == "1+":
            return "HER2-positive_IHC1_ISH-positive_discordant"
        if ish == "Positive" and is_missing(ihc):
            return "HER2-positive_ISH-positive_IHC-missing"
        if clean(row.get("clinical_her2_group_confidence")) == "inferred":
            return "HER2-positive_inferred_receptor-positive"
        return "HER2-positive_other"

    if group == "HER2-low":
        if ihc == "1+" and ish == "Negative":
            return "HER2-low_IHC1_ISH-negative"
        if ihc == "1+" and is_missing(ish):
            return "HER2-low_IHC1_ISH-not-evaluated"
        if ihc == "2+" and ish == "Negative":
            return "HER2-low_IHC2_ISH-negative"
        return "HER2-low_other"

    if group == "HER2-zero":
        if ihc == "0" and ish == "Negative":
            return "HER2-zero_IHC0_ISH-negative"
        if ihc == "0" and is_missing(ish):
            return "HER2-zero_IHC0_ISH-not-evaluated"
        return "HER2-zero_other"

    return "HER2-unknown_or_other"


def discordance_flags(row) -> list[str]:
    flags: list[str] = []
    group = clean(row.get("clinical_her2_group"))
    confidence = clean(row.get("clinical_her2_group_confidence"))
    ihc_status = clean(row.get("her2_ihc_receptor_status"))
    ihc = clean(row.get("her2_ihc_score"))
    ish = clean(row.get("her2_ish_status"))

    if confidence == "inferred":
        flags.append("inferred_label_missing_detailed_ihc_ish")
    if ihc in {"0", "1+"} and ish == "Positive":
        flags.append("low_or_zero_ihc_with_positive_ish")
    if ihc == "3+" and ish == "Negative":
        flags.append("ihc3_with_negative_ish")
    if ish in {"Equivocal", "Indeterminate"}:
        flags.append(f"ish_{ish.lower()}")
    if ihc == "2+" and ish not in {"Negative", "Positive"}:
        flags.append("ihc2_without_definitive_ish")
    if group == "HER2-low" and ihc_status == "Positive":
        flags.append("low_group_with_positive_receptor_status")
    if group == "HER2-zero" and ihc_status == "Positive":
        flags.append("zero_group_with_positive_receptor_status")
    return flags


def load_optional_csv(pd, path: Path):
    if not path.exists():
        return None
    return pd.read_csv(path)


def load_patient_clinical_metadata(pd, path: Path):
    if not path.exists():
        return None
    clinical = pd.read_csv(path, sep="\t", header=1, skiprows=[2], low_memory=False)
    columns = [
        "bcr_patient_barcode",
        "gender",
        "histological_type",
        "pathologic_stage",
        "history_of_neoadjuvant_treatment",
    ]
    available = [column for column in columns if column in clinical.columns]
    if "bcr_patient_barcode" not in available:
        return None
    clinical = clinical[available].rename(
        columns={
            "bcr_patient_barcode": "case_submitter_id",
            "gender": "patient_gender",
            "history_of_neoadjuvant_treatment": "history_neoadjuvant_treatment",
        }
    )
    return clinical.drop_duplicates("case_submitter_id")


def add_openslide_checks(rows, check_openslide: bool):
    if not check_openslide:
        rows["openslide_checked"] = False
        rows["openslide_ok"] = ""
        rows["slide_width"] = ""
        rows["slide_height"] = ""
        rows["openslide_error"] = ""
        return rows

    try:
        import openslide
    except ModuleNotFoundError:
        rows["openslide_checked"] = False
        rows["openslide_ok"] = ""
        rows["slide_width"] = ""
        rows["slide_height"] = ""
        rows["openslide_error"] = "openslide Python package not available"
        return rows

    openslide_ok = []
    widths = []
    heights = []
    errors = []
    for path_value in rows["slide_local_path"]:
        path = Path(str(path_value))
        try:
            slide = openslide.OpenSlide(str(path))
            width, height = slide.dimensions
            slide.close()
            openslide_ok.append(True)
            widths.append(int(width))
            heights.append(int(height))
            errors.append("")
        except Exception as exc:  # noqa: BLE001 - record QC failure details.
            openslide_ok.append(False)
            widths.append("")
            heights.append("")
            errors.append(str(exc))
    rows["openslide_checked"] = True
    rows["openslide_ok"] = openslide_ok
    rows["slide_width"] = widths
    rows["slide_height"] = heights
    rows["openslide_error"] = errors
    return rows


def summarize_group_counts(pd, rows, column: str) -> str:
    if column not in rows.columns:
        return "Not available."
    counts = rows.groupby(["clinical_her2_group", column]).size().reset_index(name="n")
    if counts.empty:
        return "No rows."
    lines = ["| HER2 group | Value | Slides |", "|---|---|---:|"]
    for _, row in counts.sort_values(["clinical_her2_group", column]).iterrows():
        lines.append(f"| {row['clinical_her2_group']} | {row[column]} | {int(row['n'])} |")
    return "\n".join(lines)


def fmt_median(series) -> str:
    values = series.dropna()
    if values.empty:
        return "NA"
    return f"{float(values.median()):.3g}"


def write_markdown(path: Path, rows, summary: dict, sources: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    high = rows[rows["label_slide_trust"] == "high_label_and_slide_trust"].copy()
    review = rows[rows["label_slide_trust"] == "review_before_primary_analysis"].copy()
    exclude = rows[rows["label_slide_trust"] == "exclude_from_primary_analysis"].copy()

    lines: list[str] = [
        "# Trustworthy TCGA-BRCA HER2 Slide List",
        "",
        "Status: Generated QC/trustworthiness list for the 61/61/61 laptop-balanced TCGA-BRCA HER2 cohort.",
        "",
        "## Bottom Line",
        "",
        f"- Slides checked: {summary['n_slides']}",
        f"- High label+slide trust: {summary['n_high_label_slide_trust']}",
        f"- Review before primary analysis: {summary['n_review_before_primary_analysis']}",
        f"- Exclude from primary analysis: {summary['n_exclude_from_primary_analysis']}",
        f"- Output CSV: `{summary['output_csv']}`",
        f"- Tracked CSV for GitHub: `{summary['output_tracked_csv']}`",
        f"- High-trust-only CSV: `{summary['output_high_trust_csv']}`",
        f"- Tracked high-trust-only CSV for GitHub: `{summary['output_tracked_high_trust_csv']}`",
        "",
        "High label+slide trust means: known clinical HER2 group, direct HER2 label rule, no flagged IHC/ISH discordance, female patient, primary-tumor diagnostic SVS metadata, local file present, exact file size match, and readable by OpenSlide when that check is available.",
        "",
        "This does not mean the slide has pathologist-confirmed tumor-rich sampled tiles. For that, we still need GigaTIME/tile-level QC or human review.",
        "",
        "## What We Borrowed From Guardia et al.",
        "",
        "The paper the advisor pointed to is Guardia et al., Genome Research 2025, PMID 40664477: `Alternative splicing generates HER2 isoform diversity underlying antibody-drug conjugate resistance in breast cancer`.",
        "",
        "Guardia et al. profiled HER2 isoforms in 561 TCGA-BRCA primary breast tumors from female patients and cell-line models. For the TCGA part, they selected and stratified samples using technical considerations such as library preparation and mapped reads, excluded male samples, and grouped tumors by hormone-receptor status plus HER2-high, HER2-low, or HER2-zero status from IHC and/or FISH.",
        "",
        "We use this paper as a TCGA cohort-cleanup and HER2 biology reference, not as validation that GigaTIME virtual mIF channels are real measured mIF.",
        "",
        "For our H&E/GigaTIME project, the practical translation is: use only traceable direct HER2 labels for primary analysis, exclude or review IHC/ISH discordance, keep female primary-tumor diagnostic slides, stratify or adjust by ER and PR, separate HER2-low and HER2-zero detail subgroups, check slide integrity and tissue quality, and treat RNA/isoform evidence as validation rather than as a substitute for clinical IHC/ISH.",
        "",
        "## Group Counts",
        "",
        summarize_group_counts(None, rows, "label_slide_trust"),
        "",
        "## HER2 Detail Subgroups",
        "",
        summarize_group_counts(None, rows, "her2_detail_subgroup"),
        "",
        "## Selection Bias Snapshot",
        "",
        "| HER2 group | Slides | Female | ER positive | PR positive | ERBB2 TPM median | Slide size median MB |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]

    for group in HER2_GROUPS:
        group_rows = rows[rows["clinical_her2_group"] == group]
        female = int((group_rows["patient_gender"] == "FEMALE").sum()) if "patient_gender" in group_rows else 0
        er_pos = int((group_rows["er_status"] == "Positive").sum()) if "er_status" in group_rows else 0
        pr_pos = int((group_rows["pr_status"] == "Positive").sum()) if "pr_status" in group_rows else 0
        lines.append(
            f"| {group} | {len(group_rows)} | {female} | {er_pos} | {pr_pos} | "
            f"{fmt_median(group_rows['erbb2_tpm'])} | {fmt_median(group_rows['slide_file_size_mb'])} |"
        )

    lines.extend(
        [
            "",
            "These ER/PR and ERBB2 RNA differences are important because a classifier could accidentally learn hormone-receptor or molecular-subtype context instead of HER2-low versus HER2-zero biology.",
            "",
            "## Clinical Context Checks",
            "",
            "Patient sex/gender:",
            "",
            summarize_group_counts(None, rows, "patient_gender"),
            "",
            "Histology:",
            "",
            summarize_group_counts(None, rows, "histological_type"),
            "",
            "Pathologic stage:",
            "",
            summarize_group_counts(None, rows, "pathologic_stage"),
            "",
            "## Available GigaTIME Tissue QC",
            "",
        ]
    )

    if summary["n_processed_with_tissue_qc"]:
        lines.extend(
            [
                f"- Slides with existing GigaTIME tissue QC: {summary['n_processed_with_tissue_qc']}",
                f"- Processed slides passing strict tissue QC: {summary['n_processed_strict_tissue_qc']}",
                "",
                "Strict processed tissue QC currently means at least 200 sampled tiles, mean tissue fraction at least 0.70, and QC-cellular retained fraction at least 0.50.",
            ]
        )
    else:
        lines.append("- No existing GigaTIME tissue QC was found for this cohort.")

    lines.extend(
        [
            "",
            "## High-Trust Slides By Group",
            "",
            "The full machine-readable list is in the CSV above. The table below lists the high-trust slide IDs selected for primary analysis.",
            "",
            "| HER2 group | Case | Detail subgroup | Slide ID | ER | PR | ERBB2 TPM | Processed tissue QC |",
            "|---|---|---|---|---|---|---:|---|",
        ]
    )

    for _, row in high.sort_values(["clinical_her2_group", "case_submitter_id"]).iterrows():
        lines.append(
            f"| {row['clinical_her2_group']} | {row['case_submitter_id']} | "
            f"{row['her2_detail_subgroup']} | {row['slide_id']} | {row.get('er_status', '')} | "
            f"{row.get('pr_status', '')} | {row.get('erbb2_tpm', '')} | {row['processed_tissue_qc']} |"
        )

    if not review.empty:
        lines.extend(
            [
                "",
                "## Review Before Primary Analysis",
                "",
                "| HER2 group | Case | Slide ID | Reason |",
                "|---|---|---|---|",
            ]
        )
        for _, row in review.sort_values(["clinical_her2_group", "case_submitter_id"]).iterrows():
            lines.append(
                f"| {row['clinical_her2_group']} | {row['case_submitter_id']} | {row['slide_id']} | "
                f"{row['trust_reasons']} |"
            )

    if not exclude.empty:
        lines.extend(
            [
                "",
                "## Excluded From Primary Analysis",
                "",
                "| HER2 group | Case | Slide ID | Reason |",
                "|---|---|---|---|",
            ]
        )
        for _, row in exclude.sort_values(["clinical_her2_group", "case_submitter_id"]).iterrows():
            lines.append(
                f"| {row['clinical_her2_group']} | {row['case_submitter_id']} | {row['slide_id']} | "
                f"{row['trust_reasons']} |"
            )

    lines.extend(
        [
            "",
            "## Sources Checked",
            "",
            *[f"- {source}" for source in sources],
            "",
            "## Recommended Use",
            "",
            "- Primary analysis: use `label_slide_trust == high_label_and_slide_trust`.",
            "- Strict sensitivity analysis: exclude `processed_tissue_qc != strict_pass` after GigaTIME is run on the full 183-slide cohort.",
            "- HER2-low sensitivity analysis: run IHC `1+` and IHC `2+`/ISH-negative separately.",
            "- HER2-zero sensitivity analysis: run IHC `0`/ISH-negative and IHC `0`/ISH-not-evaluated separately.",
            "- Discordance sensitivity analysis: exclude all cases with non-empty `discordance_flags`.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", default="data/tcga_brca/clinical_her2_labels.csv")
    parser.add_argument("--cases", default="data/tcga_brca/clinical_her2_laptop_balanced61_cases.csv")
    parser.add_argument("--slides", default="data/tcga_brca/clinical_her2_laptop_balanced61_slides_files.csv")
    parser.add_argument("--patient-clinical", default="data/tcga_brca/clinical/nationwidechildrens.org_clinical_patient_brca.txt")
    parser.add_argument(
        "--slide-scores",
        default="results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/slide_scores.csv",
    )
    parser.add_argument(
        "--retention-summary",
        default="results/gigatime_tcga_brca_clinical_her2_expanded20_tile256/gigatime_cleanup/filter_retention_summary.csv",
    )
    parser.add_argument("--out-csv", default="data/tcga_brca/clinical_her2_laptop_balanced61_trustworthy_slides.csv")
    parser.add_argument("--out-summary", default="data/tcga_brca/clinical_her2_laptop_balanced61_trustworthy_summary.json")
    parser.add_argument(
        "--out-tracked-csv",
        default="docs/assets/clinical_her2_trustworthy_slide_list/trustworthy_slides.csv",
    )
    parser.add_argument(
        "--out-tracked-summary",
        default="docs/assets/clinical_her2_trustworthy_slide_list/trustworthy_summary.json",
    )
    parser.add_argument(
        "--out-high-trust-csv",
        default="data/tcga_brca/clinical_her2_laptop_balanced61_high_trust_slides.csv",
    )
    parser.add_argument(
        "--out-tracked-high-trust-csv",
        default="docs/assets/clinical_her2_trustworthy_slide_list/high_trust_slides.csv",
    )
    parser.add_argument("--out-markdown", default="docs/clinical_her2_trustworthy_slide_list.md")
    parser.add_argument("--check-openslide", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    pd = require_pandas()

    labels = pd.read_csv(args.labels)
    cases = pd.read_csv(args.cases)
    slides = pd.read_csv(args.slides)
    patient_clinical = load_patient_clinical_metadata(pd, Path(args.patient_clinical))

    slide_columns = set(slides.columns)
    case_keep = [
        column
        for column in cases.columns
        if column in {"case_submitter_id", "cohort_group", "selection_rank"} or column not in slide_columns
    ]
    rows = slides.merge(cases[case_keep], on=["case_submitter_id", "cohort_group", "selection_rank"], how="left")
    if "clinical_her2_group" not in rows.columns:
        rows = rows.merge(labels, on="case_submitter_id", how="left")
    if patient_clinical is not None:
        rows = rows.merge(patient_clinical, on="case_submitter_id", how="left")
    else:
        rows["patient_gender"] = ""
        rows["histological_type"] = ""
        rows["pathologic_stage"] = ""
        rows["history_neoadjuvant_treatment"] = ""

    rows["slide_id"] = rows["slide_file_name"].map(lambda value: Path(str(value)).stem)
    rows["slide_path_exists_now"] = rows["slide_local_path"].map(lambda value: Path(str(value)).exists())
    rows["slide_actual_size"] = rows["slide_local_path"].map(
        lambda value: Path(str(value)).stat().st_size if Path(str(value)).exists() else None
    )
    rows["slide_file_size"] = pd.to_numeric(rows["slide_file_size"], errors="coerce")
    rows["slide_file_size_mb"] = rows["slide_file_size"] / (1024 * 1024)
    rows["slide_size_matches"] = rows["slide_actual_size"] == rows["slide_file_size"]
    rows["slide_local_exists"] = rows["slide_local_exists"].map(bool_from_value)
    rows["is_primary_tumor_slide"] = rows["slide_sample_type"] == "Primary Tumor"
    rows["is_diagnostic_svs"] = (
        (rows["slide_data_type"] == "Slide Image")
        & (rows["slide_data_format"] == "SVS")
        & (rows["slide_experimental_strategy"].isin({"Tissue Slide", "Diagnostic Slide"}))
    )
    rows["sample_barcode_primary_tumor_code"] = rows["slide_sample_submitter_id"].map(
        lambda value: str(value).split("-")[3][:2] if len(str(value).split("-")) >= 4 else ""
    )
    rows["sample_barcode_is_01"] = rows["sample_barcode_primary_tumor_code"] == "01"
    rows["known_her2_group"] = rows["clinical_her2_group"].isin(HER2_GROUPS)
    rows["direct_label"] = rows["clinical_her2_group_confidence"] == "direct"
    rows["female_patient"] = rows["patient_gender"] == "FEMALE"
    rows["her2_detail_subgroup"] = rows.apply(her2_detail_subgroup, axis=1)
    rows["discordance_flags"] = rows.apply(lambda row: ";".join(discordance_flags(row)), axis=1)
    rows["has_discordance_flag"] = rows["discordance_flags"] != ""

    rows = add_openslide_checks(rows, args.check_openslide)
    if "openslide_ok" in rows.columns and rows["openslide_ok"].astype(str).isin({"True", "False"}).any():
        openslide_ok_or_unchecked = rows["openslide_ok"].astype(str) == "True"
    else:
        openslide_ok_or_unchecked = True

    mandatory_slide_ok = (
        rows["slide_path_exists_now"]
        & rows["slide_size_matches"]
        & rows["is_primary_tumor_slide"]
        & rows["sample_barcode_is_01"]
        & rows["is_diagnostic_svs"]
        & rows["female_patient"]
        & openslide_ok_or_unchecked
    )
    rows["mandatory_slide_checks_pass"] = mandatory_slide_ok

    trust_labels = []
    trust_reasons = []
    for _, row in rows.iterrows():
        reasons = []
        if not bool(row["known_her2_group"]):
            reasons.append("unknown HER2 group")
        if not bool(row["direct_label"]):
            reasons.append("label is inferred rather than direct")
        if bool(row["has_discordance_flag"]):
            reasons.append(f"HER2 field flag: {row['discordance_flags']}")
        if not bool(row["slide_path_exists_now"]):
            reasons.append("local slide file missing")
        if not bool(row["slide_size_matches"]):
            reasons.append("local slide size mismatch")
        if not bool(row["is_primary_tumor_slide"]):
            reasons.append("slide metadata is not Primary Tumor")
        if not bool(row["sample_barcode_is_01"]):
            reasons.append("sample barcode is not TCGA primary tumor code 01")
        if not bool(row["is_diagnostic_svs"]):
            reasons.append("slide metadata is not a diagnostic/primary SVS slide")
        if not bool(row["female_patient"]):
            reasons.append("patient gender is not FEMALE")
        if clean(row.get("openslide_ok")) == "False":
            reasons.append("OpenSlide could not read slide")

        if not bool(row["mandatory_slide_checks_pass"]) or not bool(row["known_her2_group"]):
            trust = "exclude_from_primary_analysis"
        elif bool(row["direct_label"]) and not bool(row["has_discordance_flag"]):
            trust = "high_label_and_slide_trust"
        else:
            trust = "review_before_primary_analysis"

        trust_labels.append(trust)
        trust_reasons.append("; ".join(reasons) if reasons else "passes label and slide metadata checks")

    rows["label_slide_trust"] = trust_labels
    rows["trust_reasons"] = trust_reasons

    slide_scores = load_optional_csv(pd, Path(args.slide_scores))
    if slide_scores is not None:
        keep = ["case_submitter_id", "slide_id", "n_tiles", "mean_tissue_fraction"]
        rows = rows.merge(slide_scores[keep], on=["case_submitter_id", "slide_id"], how="left")
    else:
        rows["n_tiles"] = float("nan")
        rows["mean_tissue_fraction"] = float("nan")

    retention = load_optional_csv(pd, Path(args.retention_summary))
    if retention is not None:
        qc = retention[retention["feature_view"] == "qc_cellular_tissue"].copy()
        keep = [
            "case_submitter_id",
            "slide_id",
            "n_tiles_total",
            "n_tiles_retained",
            "retained_fraction",
            "mean_DAPI",
            "mean_CK",
        ]
        rows = rows.merge(qc[keep], on=["case_submitter_id", "slide_id"], how="left")
    else:
        rows["retained_fraction"] = float("nan")

    rows["processed_tissue_qc"] = "not_processed_yet"
    processed = rows["n_tiles"].notna()
    strict = processed & (rows["n_tiles"] >= 200) & (rows["mean_tissue_fraction"] >= 0.70) & (rows["retained_fraction"] >= 0.50)
    moderate = processed & ~strict & (rows["n_tiles"] >= 128) & (rows["mean_tissue_fraction"] >= 0.60) & (rows["retained_fraction"] >= 0.35)
    weak = processed & ~(strict | moderate)
    rows.loc[strict, "processed_tissue_qc"] = "strict_pass"
    rows.loc[moderate, "processed_tissue_qc"] = "moderate_pass"
    rows.loc[weak, "processed_tissue_qc"] = "weak_or_review"

    output_columns = [
        "label_slide_trust",
        "processed_tissue_qc",
        "trust_reasons",
        "clinical_her2_group",
        "her2_detail_subgroup",
        "case_submitter_id",
        "slide_id",
        "slide_file_id",
        "slide_file_name",
        "slide_local_path",
        "clinical_her2_group_rule",
        "clinical_her2_group_confidence",
        "discordance_flags",
        "her2_ihc_receptor_status",
        "her2_ihc_score",
        "her2_ish_status",
        "her2_cep17_ratio",
        "er_status",
        "pr_status",
        "patient_gender",
        "histological_type",
        "pathologic_stage",
        "history_neoadjuvant_treatment",
        "erbb2_tpm",
        "slide_sample_submitter_id",
        "slide_sample_type",
        "slide_file_size_mb",
        "slide_path_exists_now",
        "slide_size_matches",
        "openslide_checked",
        "openslide_ok",
        "slide_width",
        "slide_height",
        "n_tiles",
        "mean_tissue_fraction",
        "retained_fraction",
        "mean_DAPI",
        "mean_CK",
    ]

    for column in output_columns:
        if column not in rows.columns:
            rows[column] = ""

    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    rows[output_columns].to_csv(args.out_csv, index=False)
    Path(args.out_tracked_csv).parent.mkdir(parents=True, exist_ok=True)
    rows[output_columns].to_csv(args.out_tracked_csv, index=False)
    high_rows = rows[rows["label_slide_trust"] == "high_label_and_slide_trust"].copy()
    Path(args.out_high_trust_csv).parent.mkdir(parents=True, exist_ok=True)
    high_rows[output_columns].to_csv(args.out_high_trust_csv, index=False)
    Path(args.out_tracked_high_trust_csv).parent.mkdir(parents=True, exist_ok=True)
    high_rows[output_columns].to_csv(args.out_tracked_high_trust_csv, index=False)

    summary = {
        "n_slides": int(len(rows)),
        "n_high_label_slide_trust": int((rows["label_slide_trust"] == "high_label_and_slide_trust").sum()),
        "n_review_before_primary_analysis": int((rows["label_slide_trust"] == "review_before_primary_analysis").sum()),
        "n_exclude_from_primary_analysis": int((rows["label_slide_trust"] == "exclude_from_primary_analysis").sum()),
        "n_processed_with_tissue_qc": int(rows["n_tiles"].notna().sum()),
        "n_processed_strict_tissue_qc": int((rows["processed_tissue_qc"] == "strict_pass").sum()),
        "n_female_patients": int((rows["patient_gender"] == "FEMALE").sum()),
        "label_slide_trust_counts": rows["label_slide_trust"].value_counts().to_dict(),
        "processed_tissue_qc_counts": rows["processed_tissue_qc"].value_counts().to_dict(),
        "her2_detail_subgroup_counts": rows["her2_detail_subgroup"].value_counts().to_dict(),
        "patient_gender_counts": rows["patient_gender"].value_counts(dropna=False).to_dict(),
        "output_csv": args.out_csv,
        "output_tracked_csv": args.out_tracked_csv,
        "output_high_trust_csv": args.out_high_trust_csv,
        "output_tracked_high_trust_csv": args.out_tracked_high_trust_csv,
        "output_markdown": args.out_markdown,
    }
    Path(args.out_summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    Path(args.out_tracked_summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")

    sources = [
        "Guardia et al., PMID 40664477: https://pubmed.ncbi.nlm.nih.gov/40664477/",
        "Guardia et al., Genome Research abstract: https://genome.cshlp.org/content/35/9/1942.abstract",
        "Guardia et al., supplemental material: https://genome.cshlp.org/content/35/9/1942/suppl/DC1",
        "Guardia et al., medRxiv full-text preprint used to verify methods language: https://www.medrxiv.org/content/10.1101/2024.11.25.24317569v1.full-text",
        "GDC Clinical Data overview: https://docs.gdc.cancer.gov/Encyclopedia/pages/Clinical_Data/",
        "GDC Data Dictionary overview: https://docs.gdc.cancer.gov/Data_Dictionary/",
        "TCGA-BRCA Breast Enrollment Form: https://gdc.cancer.gov/system/files/public/file/Breast%20Enrollment%20Form.pdf",
        "CAP/ASCO HER2 testing guideline page: https://www.cap.org/protocols-and-guidelines/cap-guidelines/current-cap-guidelines/recommendations-for-human-epidermal-growth-factor-2-testing-in-breast-cancer",
    ]
    write_markdown(Path(args.out_markdown), rows, summary, sources)

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
