#!/usr/bin/env python3
"""Compare GigaTIME virtual programs with broader RNA immune programs."""

from __future__ import annotations

import argparse
import gzip
import math
import os
from pathlib import Path


RNA_PROGRAMS: dict[str, dict[str, object]] = {
    "t_cell_cytotoxic": {
        "label": "T cell / cytotoxic",
        "genes": ["CD3D", "CD3E", "CD3G", "TRAC", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7", "GNLY"],
    },
    "checkpoint_ifng": {
        "label": "Checkpoint / IFNG",
        "genes": ["CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "IFNG", "CXCL9", "CXCL10", "IDO1", "STAT1"],
    },
    "myeloid_macrophage": {
        "label": "Myeloid / macrophage",
        "genes": ["CD68", "CD163", "MRC1", "CSF1R", "LST1", "AIF1", "TYROBP", "C1QA", "C1QB", "C1QC"],
    },
    "dendritic_apc": {
        "label": "Dendritic / APC",
        "genes": ["ITGAX", "IRF8", "BATF3", "CLEC9A", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74"],
    },
    "b_cell": {
        "label": "B cell",
        "genes": ["MS4A1", "CD79A", "CD79B", "CD19", "CD22", "BANK1", "CD74"],
    },
    "proliferation": {
        "label": "Proliferation",
        "genes": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "CCNB1", "BIRC5"],
    },
    "epithelial": {
        "label": "Epithelial / tumor",
        "genes": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1"],
    },
    "stromal_fibroblast": {
        "label": "Stromal / fibroblast",
        "genes": ["COL1A1", "COL1A2", "ACTA2", "VIM", "DCN", "LUM", "TAGLN"],
    },
    "endothelial": {
        "label": "Endothelial",
        "genes": ["PECAM1", "VWF", "CD34", "KDR", "ENG"],
    },
}

VIRTUAL_PROGRAMS: dict[str, dict[str, object]] = {
    "myeloid_checkpoint": {
        "label": "Virtual myeloid/checkpoint",
        "channels": ["CD68", "CD11c", "PD-L1"],
    },
    "t_cell_checkpoint": {
        "label": "Virtual T cell/checkpoint",
        "channels": ["CD3", "CD4", "CD8", "PD-1"],
    },
    "immune_checkpoint": {
        "label": "Virtual all immune/checkpoint",
        "channels": ["CD3", "CD4", "CD8", "CD20", "CD68", "CD11c", "PD-1", "PD-L1"],
    },
    "proliferation": {
        "label": "Virtual proliferation",
        "channels": ["Ki67"],
    },
    "epithelial": {
        "label": "Virtual epithelial",
        "channels": ["CK"],
    },
}

GROUP_ORDER = ["HER2-positive", "HER2-low", "HER2-zero"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--joined",
        default=(
            "results/gigatime_tcga_brca_clinical_her2_tile256/clinical_summary/"
            "joined_slide_clinical_her2_gigatime.csv"
        ),
        help="Joined clinical HER2 + GigaTIME slide score table.",
    )
    parser.add_argument("--expression-dir", default="data/tcga_brca/expression_files")
    parser.add_argument("--out-dir", default="results/gigatime_tcga_brca_clinical_her2_tile256/rna_program_validation")
    return parser.parse_args()


def require_analysis_libs(mpl_config_dir: Path):
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
        from scipy import stats
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"Missing Python package: {exc.name}. Use `conda activate gigatime-tcga` "
            "or `conda run -n gigatime-tcga ...`."
        ) from exc
    return pd, plt, sns, stats


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def find_expression_file(expression_dir: Path, case_id: str) -> Path:
    case_dir = expression_dir / case_id
    candidates = sorted(case_dir.glob("*.rna_seq.augmented_star_gene_counts.tsv*"))
    if not candidates:
        candidates = sorted(case_dir.glob("*.tsv*"))
    if not candidates:
        raise FileNotFoundError(f"No STAR-count expression file found for {case_id} under {case_dir}")
    return candidates[0]


def read_gene_tpms(path: Path, target_genes: set[str]) -> dict[str, float]:
    header: list[str] | None = None
    values: dict[str, float] = {}
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            columns = line.split("\t")
            if columns[0] == "gene_id":
                header = columns
                continue
            if header is None:
                continue
            row = dict(zip(header, columns))
            gene_name = row.get("gene_name", "")
            if gene_name not in target_genes:
                continue
            try:
                values[gene_name] = float(row.get("tpm_unstranded", "nan"))
            except ValueError:
                values[gene_name] = float("nan")
    return values


def log2p1(value: float) -> float:
    if value is None or math.isnan(value):
        return float("nan")
    return math.log2(value + 1.0)


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    indexed = [(idx, p) for idx, p in enumerate(p_values) if not math.isnan(p)]
    if not indexed:
        return [float("nan")] * len(p_values)
    ranked = sorted(indexed, key=lambda item: item[1])
    m = len(ranked)
    adjusted = [float("nan")] * len(p_values)
    prev = 1.0
    for rank, (idx, p_value) in reversed(list(enumerate(ranked, start=1))):
        q_value = min(prev, p_value * m / rank)
        adjusted[idx] = q_value
        prev = q_value
    return adjusted


def describe_tile_sampling(joined) -> str:
    if "n_tiles" not in joined.columns:
        return "Tile sampling information was not available."
    counts = joined.drop_duplicates("case_submitter_id")["n_tiles"].dropna()
    if counts.empty:
        return "Tile sampling information was not available."
    unique_counts = sorted({int(value) for value in counts})
    if len(unique_counts) == 1:
        return f"{unique_counts[0]} tissue tiles per slide"
    return f"median {int(counts.median())} tissue tiles per slide (range {int(counts.min())}-{int(counts.max())})"


def build_rna_program_table(pd, joined, expression_dir: Path):
    target_genes = {gene for spec in RNA_PROGRAMS.values() for gene in spec["genes"]}
    rows = []
    gene_rows = []
    for _, joined_row in joined.drop_duplicates("case_submitter_id").iterrows():
        case_id = joined_row["case_submitter_id"]
        expression_file = find_expression_file(expression_dir, case_id)
        gene_tpms = read_gene_tpms(expression_file, target_genes)
        row = {
            "case_submitter_id": case_id,
            "clinical_her2_group": joined_row["clinical_her2_group"],
            "expression_file": str(expression_file),
        }
        for program, spec in RNA_PROGRAMS.items():
            genes = list(spec["genes"])
            available = [gene for gene in genes if gene in gene_tpms and not math.isnan(gene_tpms[gene])]
            log_values = [log2p1(gene_tpms[gene]) for gene in available]
            row[f"rna_{program}"] = sum(log_values) / len(log_values) if log_values else float("nan")
            row[f"rna_{program}_genes_found"] = ";".join(available)
            for gene in genes:
                value = gene_tpms.get(gene, float("nan"))
                gene_rows.append(
                    {
                        "case_submitter_id": case_id,
                        "clinical_her2_group": joined_row["clinical_her2_group"],
                        "rna_program": program,
                        "rna_program_label": spec["label"],
                        "gene": gene,
                        "tpm": value,
                        "log2_tpm_plus_1": log2p1(value),
                    }
                )
        rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(gene_rows)


def build_virtual_program_table(pd, joined):
    rows = []
    for _, joined_row in joined.drop_duplicates("case_submitter_id").iterrows():
        row = {
            "case_submitter_id": joined_row["case_submitter_id"],
            "clinical_her2_group": joined_row["clinical_her2_group"],
            "n_tiles": joined_row.get("n_tiles", float("nan")),
        }
        for program, spec in VIRTUAL_PROGRAMS.items():
            values = []
            channels_found = []
            for channel in spec["channels"]:
                col = f"mean_{channel}"
                if col in joined_row.index and not math.isnan(float(joined_row[col])):
                    values.append(float(joined_row[col]))
                    channels_found.append(channel)
            row[f"virtual_{program}"] = sum(values) / len(values) if values else float("nan")
            row[f"virtual_{program}_channels_found"] = ";".join(channels_found)
        rows.append(row)
    return pd.DataFrame(rows)


def build_correlation_table(pd, stats, merged):
    rows = []
    for virtual_program, virtual_spec in VIRTUAL_PROGRAMS.items():
        virtual_col = f"virtual_{virtual_program}"
        for rna_program, rna_spec in RNA_PROGRAMS.items():
            rna_col = f"rna_{rna_program}"
            pair_df = merged[["case_submitter_id", "clinical_her2_group", virtual_col, rna_col]].dropna()
            if len(pair_df) >= 3 and pair_df[virtual_col].nunique() > 1 and pair_df[rna_col].nunique() > 1:
                spearman = stats.spearmanr(pair_df[virtual_col], pair_df[rna_col])
                rho = float(spearman.statistic)
                p_value = float(spearman.pvalue)
            else:
                rho = float("nan")
                p_value = float("nan")
            rows.append(
                {
                    "virtual_program": virtual_program,
                    "virtual_program_label": virtual_spec["label"],
                    "rna_program": rna_program,
                    "rna_program_label": rna_spec["label"],
                    "n_cases": int(len(pair_df)),
                    "spearman_rho": rho,
                    "spearman_p_value": p_value,
                }
            )
    q_values = benjamini_hochberg([row["spearman_p_value"] for row in rows])
    for row, q_value in zip(rows, q_values):
        row["spearman_q_value_bh"] = q_value
    return pd.DataFrame(rows)


def build_group_test_table(pd, stats, merged, prefix: str, programs: dict[str, dict[str, object]]):
    rows = []
    for program, spec in programs.items():
        col = f"{prefix}_{program}"
        groups = []
        group_means = {}
        group_medians = {}
        group_counts = {}
        for group in GROUP_ORDER:
            values = merged.loc[merged["clinical_her2_group"] == group, col].dropna()
            groups.append(values)
            group_counts[group] = int(values.shape[0])
            group_means[group] = float(values.mean()) if len(values) else float("nan")
            group_medians[group] = float(values.median()) if len(values) else float("nan")
        if all(len(values) >= 2 for values in groups):
            result = stats.kruskal(*groups)
            p_value = float(result.pvalue)
        else:
            p_value = float("nan")
        valid_means = {group: value for group, value in group_means.items() if not math.isnan(value)}
        highest_group = max(valid_means, key=valid_means.get) if valid_means else ""
        lowest_group = min(valid_means, key=valid_means.get) if valid_means else ""
        rows.append(
            {
                "program": program,
                "program_label": spec["label"],
                "n_cases": int(sum(group_counts.values())),
                "kruskal_p_value": p_value,
                "highest_mean_group": highest_group,
                "lowest_mean_group": lowest_group,
                "max_minus_min_mean": (
                    max(valid_means.values()) - min(valid_means.values()) if valid_means else float("nan")
                ),
                "her2_positive_mean": group_means["HER2-positive"],
                "her2_low_mean": group_means["HER2-low"],
                "her2_zero_mean": group_means["HER2-zero"],
                "her2_positive_median": group_medians["HER2-positive"],
                "her2_low_median": group_medians["HER2-low"],
                "her2_zero_median": group_medians["HER2-zero"],
            }
        )
    q_values = benjamini_hochberg([row["kruskal_p_value"] for row in rows])
    for row, q_value in zip(rows, q_values):
        row["kruskal_q_value_bh"] = q_value
    return pd.DataFrame(rows)


def plot_correlation_heatmap(plt, sns, correlation, out_dir: Path) -> None:
    matrix = correlation.pivot(index="virtual_program_label", columns="rna_program_label", values="spearman_rho")
    matrix = matrix.reindex([spec["label"] for spec in VIRTUAL_PROGRAMS.values()])
    matrix = matrix[[spec["label"] for spec in RNA_PROGRAMS.values()]]
    plt.figure(figsize=(12.8, 5.8))
    sns.heatmap(matrix, cmap="vlag", center=0, annot=True, fmt=".2f", vmin=-1, vmax=1, linewidths=0.35)
    plt.title("GigaTIME virtual programs vs RNA programs")
    plt.xlabel("RNA program")
    plt.ylabel("GigaTIME virtual program")
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plt.savefig(out_dir / "virtual_rna_program_correlation_heatmap.png", dpi=180)
    plt.close()


def plot_group_programs(plt, sns, pd, merged, summary, prefix: str, programs: dict[str, dict[str, object]], out_path: Path) -> None:
    ordered_programs = summary.sort_values("kruskal_p_value")["program"].tolist()
    label_lookup = {program: str(programs[program]["label"]) for program in programs}
    records = []
    for program in ordered_programs:
        col = f"{prefix}_{program}"
        for _, row in merged[["case_submitter_id", "clinical_her2_group", col]].dropna().iterrows():
            records.append(
                {
                    "case_submitter_id": row["case_submitter_id"],
                    "clinical_her2_group": row["clinical_her2_group"],
                    "program": label_lookup[program],
                    "score": row[col],
                }
            )
    if not records:
        return
    plot_df = pd.DataFrame(records)
    plot_df["clinical_her2_group"] = pd.Categorical(plot_df["clinical_her2_group"], GROUP_ORDER, ordered=True)
    grid = sns.catplot(
        data=plot_df,
        x="clinical_her2_group",
        y="score",
        col="program",
        col_wrap=3,
        kind="box",
        height=3.0,
        aspect=1.05,
        sharey=False,
        color="#b9d8c2",
    )
    for ax in grid.axes.flat:
        program = ax.get_title().replace("program = ", "")
        subset = plot_df.loc[plot_df["program"] == program]
        sns.stripplot(
            data=subset,
            x="clinical_her2_group",
            y="score",
            order=GROUP_ORDER,
            ax=ax,
            size=4,
            color="#2f4858",
            alpha=0.75,
        )
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=25)
    grid.set_titles("{col_name}")
    grid.set_axis_labels("", "Score")
    grid.fig.suptitle("Program scores by clinical HER2 group", y=1.02)
    grid.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(grid.fig)


def format_float(value: float, digits: int = 3) -> str:
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
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def write_markdown(path: Path, joined, correlation, rna_group_summary, virtual_group_summary) -> None:
    case_count = joined["case_submitter_id"].nunique()
    tile_sampling = describe_tile_sampling(joined)
    top_corr = correlation.assign(abs_rho=correlation["spearman_rho"].abs()).sort_values(
        ["spearman_q_value_bh", "spearman_p_value", "abs_rho"], ascending=[True, True, False]
    ).head(10)
    rna_top = rna_group_summary.sort_values(["kruskal_q_value_bh", "kruskal_p_value"]).head(8)
    virtual_top = virtual_group_summary.sort_values(["kruskal_q_value_bh", "kruskal_p_value"]).head(5)

    lines = [
        "# RNA Program Validation",
        "",
        f"- Cases with paired GigaTIME and RNA-seq data: {case_count}",
        f"- Tile sampling: {tile_sampling}.",
        f"- RNA programs tested: {len(RNA_PROGRAMS)}",
        f"- GigaTIME virtual programs tested: {len(VIRTUAL_PROGRAMS)}",
        "- RNA program scores are means of log2(TPM + 1) across marker genes found in each case.",
        "- This is still indirect validation, not real mIF ground truth.",
        "",
        "## Strongest Virtual-vs-RNA Program Associations",
        "",
        markdown_table(
            ["Virtual program", "RNA program", "Spearman rho", "p", "BH q"],
            [
                [
                    row["virtual_program_label"],
                    row["rna_program_label"],
                    format_float(row["spearman_rho"], 3),
                    format_float(row["spearman_p_value"], 4),
                    format_float(row["spearman_q_value_bh"], 4),
                ]
                for _, row in top_corr.iterrows()
            ],
        ),
        "",
        "## RNA Program Differences Across Clinical HER2 Groups",
        "",
        markdown_table(
            ["RNA program", "Kruskal p", "BH q", "Highest group", "Lowest group", "Max-min mean"],
            [
                [
                    row["program_label"],
                    format_float(row["kruskal_p_value"], 4),
                    format_float(row["kruskal_q_value_bh"], 4),
                    row["highest_mean_group"],
                    row["lowest_mean_group"],
                    format_float(row["max_minus_min_mean"], 3),
                ]
                for _, row in rna_top.iterrows()
            ],
        ),
        "",
        "## Virtual Program Differences Across Clinical HER2 Groups",
        "",
        markdown_table(
            ["Virtual program", "Kruskal p", "BH q", "Highest group", "Lowest group", "Max-min mean"],
            [
                [
                    row["program_label"],
                    format_float(row["kruskal_p_value"], 4),
                    format_float(row["kruskal_q_value_bh"], 4),
                    row["highest_mean_group"],
                    row["lowest_mean_group"],
                    format_float(row["max_minus_min_mean"], 4),
                ]
                for _, row in virtual_top.iterrows()
            ],
        ),
        "",
        "## Interpretation Guardrails",
        "",
        "- A strong positive virtual-vs-RNA correlation would support biological plausibility of a GigaTIME program.",
        "- A weak or negative correlation does not automatically prove the virtual signal is wrong, because TCGA bulk RNA-seq and H&E tile sampling measure different tissue material.",
        "- RNA program differences across HER2 groups help distinguish whether the RNA data itself supports the same HER2-zero versus HER2-low immune pattern.",
        "- These results should guide validation planning and advisor review, not final biological claims.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, stats = require_analysis_libs(out_dir / ".matplotlib")

    joined = pd.read_csv(args.joined)
    rna_programs, gene_expression = build_rna_program_table(pd, joined, Path(args.expression_dir))
    virtual_programs = build_virtual_program_table(pd, joined)
    merged = virtual_programs.merge(rna_programs, on=["case_submitter_id", "clinical_her2_group"], how="inner")
    correlation = build_correlation_table(pd, stats, merged)
    rna_group_summary = build_group_test_table(pd, stats, merged, "rna", RNA_PROGRAMS)
    virtual_group_summary = build_group_test_table(pd, stats, merged, "virtual", VIRTUAL_PROGRAMS)

    rna_programs.to_csv(out_dir / "case_rna_programs.csv", index=False)
    gene_expression.to_csv(out_dir / "case_rna_program_gene_expression_long.csv", index=False)
    virtual_programs.to_csv(out_dir / "case_virtual_programs.csv", index=False)
    merged.to_csv(out_dir / "joined_virtual_rna_programs.csv", index=False)
    correlation.to_csv(out_dir / "virtual_rna_program_correlations.csv", index=False)
    rna_group_summary.to_csv(out_dir / "rna_program_group_summary.csv", index=False)
    virtual_group_summary.to_csv(out_dir / "virtual_program_group_summary.csv", index=False)

    plot_correlation_heatmap(plt, sns, correlation, out_dir)
    plot_group_programs(
        plt,
        sns,
        pd,
        merged,
        rna_group_summary,
        "rna",
        RNA_PROGRAMS,
        out_dir / "rna_programs_by_her2_group.png",
    )
    plot_group_programs(
        plt,
        sns,
        pd,
        merged,
        virtual_group_summary,
        "virtual",
        VIRTUAL_PROGRAMS,
        out_dir / "virtual_programs_by_her2_group.png",
    )
    write_markdown(out_dir / "rna_program_validation_summary.md", joined, correlation, rna_group_summary, virtual_group_summary)
    print(f"Wrote RNA program validation outputs to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
