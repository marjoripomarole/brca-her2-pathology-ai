#!/usr/bin/env python3
"""Compare clinical HER2 GigaTIME virtual channels with RNA marker signatures."""

from __future__ import annotations

import argparse
import gzip
import math
import os
from pathlib import Path

CHANNEL_SIGNATURES: dict[str, list[str]] = {
    "CD3": ["CD3D", "CD3E", "CD3G", "TRAC"],
    "CD8": ["CD8A", "CD8B"],
    "CD4": ["CD4"],
    "CD20": ["MS4A1", "CD79A", "CD79B"],
    "CD68": ["CD68", "CD163", "MRC1"],
    "CD11c": ["ITGAX", "IRF8", "BATF3"],
    "PD-1": ["PDCD1"],
    "PD-L1": ["CD274"],
    "CK": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Ki67": ["MKI67", "TOP2A"],
}

GROUP_ORDER = ["HER2-positive", "HER2-low", "HER2-zero"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--joined",
        default="results/gigatime_tcga_brca_clinical_her2/clinical_summary/joined_slide_clinical_her2_gigatime.csv",
        help="Joined clinical HER2 + GigaTIME slide score table.",
    )
    parser.add_argument("--expression-dir", default="data/tcga_brca/expression_files")
    parser.add_argument("--out-dir", default="results/gigatime_tcga_brca_clinical_her2/rna_validation")
    parser.add_argument("--channels", default=",".join(CHANNEL_SIGNATURES))
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


def build_case_signature_table(pd, joined, expression_dir: Path, channels: list[str]):
    target_genes = {gene for channel in channels for gene in CHANNEL_SIGNATURES[channel]}
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
        for channel in channels:
            genes = CHANNEL_SIGNATURES[channel]
            available = [gene for gene in genes if gene in gene_tpms and not math.isnan(gene_tpms[gene])]
            log_values = [log2p1(gene_tpms[gene]) for gene in available]
            row[f"{channel}_rna_signature_log2_tpm_mean"] = (
                sum(log_values) / len(log_values) if log_values else float("nan")
            )
            row[f"{channel}_rna_signature_genes_found"] = ";".join(available)
            for gene in genes:
                value = gene_tpms.get(gene, float("nan"))
                gene_rows.append(
                    {
                        "case_submitter_id": case_id,
                        "clinical_her2_group": joined_row["clinical_her2_group"],
                        "channel": channel,
                        "gene": gene,
                        "tpm": value,
                        "log2_tpm_plus_1": log2p1(value),
                    }
                )
        rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(gene_rows)


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


def build_correlation_table(pd, stats, joined, case_signatures, channels: list[str]):
    merged = joined.merge(case_signatures, on=["case_submitter_id", "clinical_her2_group"], how="inner")
    rows = []
    for channel in channels:
        gigatime_col = f"mean_{channel}"
        rna_col = f"{channel}_rna_signature_log2_tpm_mean"
        if gigatime_col not in merged.columns or rna_col not in merged.columns:
            continue
        pair_df = merged[["case_submitter_id", "clinical_her2_group", gigatime_col, rna_col]].dropna()
        if len(pair_df) >= 3 and pair_df[gigatime_col].nunique() > 1 and pair_df[rna_col].nunique() > 1:
            spearman = stats.spearmanr(pair_df[gigatime_col], pair_df[rna_col])
            rho = float(spearman.statistic)
            p_value = float(spearman.pvalue)
        else:
            rho = float("nan")
            p_value = float("nan")
        rows.append(
            {
                "channel": channel,
                "n_cases": int(len(pair_df)),
                "spearman_rho": rho,
                "spearman_p_value": p_value,
                "rna_signature_genes": ";".join(CHANNEL_SIGNATURES[channel]),
                "median_gigatime_mean": float(pair_df[gigatime_col].median()) if len(pair_df) else float("nan"),
                "median_rna_signature_log2_tpm_mean": float(pair_df[rna_col].median()) if len(pair_df) else float("nan"),
            }
        )
    q_values = benjamini_hochberg([row["spearman_p_value"] for row in rows])
    for row, q_value in zip(rows, q_values):
        row["spearman_q_value_bh"] = q_value
    return pd.DataFrame(rows), merged


def build_group_summary(pd, merged, channels: list[str]):
    rows = []
    for channel in channels:
        gigatime_col = f"mean_{channel}"
        rna_col = f"{channel}_rna_signature_log2_tpm_mean"
        if gigatime_col not in merged.columns or rna_col not in merged.columns:
            continue
        for group in GROUP_ORDER:
            group_df = merged.loc[merged["clinical_her2_group"] == group]
            rows.append(
                {
                    "channel": channel,
                    "clinical_her2_group": group,
                    "n_cases": int(group_df["case_submitter_id"].nunique()),
                    "gigatime_mean": float(group_df[gigatime_col].mean()),
                    "gigatime_median": float(group_df[gigatime_col].median()),
                    "rna_signature_mean": float(group_df[rna_col].mean()),
                    "rna_signature_median": float(group_df[rna_col].median()),
                }
            )
    return pd.DataFrame(rows)


def plot_correlation_heatmap(plt, sns, correlation, out_dir: Path) -> None:
    if correlation.empty:
        return
    plot_df = correlation.set_index("channel")[["spearman_rho"]].sort_values("spearman_rho", ascending=False)
    plt.figure(figsize=(4.5, max(3.5, 0.42 * len(plot_df))))
    sns.heatmap(plot_df, cmap="vlag", center=0, annot=True, fmt=".2f", vmin=-1, vmax=1, linewidths=0.3)
    plt.title("GigaTIME vs RNA Signature Spearman Correlation")
    plt.ylabel("Virtual channel")
    plt.xlabel("")
    plt.tight_layout()
    plt.savefig(out_dir / "gigatime_rna_correlation_heatmap.png", dpi=180)
    plt.close()


def plot_top_scatter(plt, sns, merged, correlation, out_dir: Path) -> None:
    if correlation.empty:
        return
    candidates = correlation.dropna(subset=["spearman_rho"]).copy()
    if candidates.empty:
        return
    candidates["abs_rho"] = candidates["spearman_rho"].abs()
    top_channels = candidates.sort_values(["abs_rho", "channel"], ascending=[False, True]).head(6)["channel"].tolist()
    records = []
    for channel in top_channels:
        gigatime_col = f"mean_{channel}"
        rna_col = f"{channel}_rna_signature_log2_tpm_mean"
        for _, row in merged[["case_submitter_id", "clinical_her2_group", gigatime_col, rna_col]].dropna().iterrows():
            records.append(
                {
                    "case_submitter_id": row["case_submitter_id"],
                    "clinical_her2_group": row["clinical_her2_group"],
                    "channel": channel,
                    "gigatime_mean": row[gigatime_col],
                    "rna_signature": row[rna_col],
                }
            )
    if not records:
        return
    import pandas as pd

    plot_df = pd.DataFrame(records)
    grid = sns.lmplot(
        data=plot_df,
        x="rna_signature",
        y="gigatime_mean",
        hue="clinical_her2_group",
        col="channel",
        col_wrap=3,
        height=3.1,
        scatter_kws={"s": 32, "alpha": 0.8},
        line_kws={"alpha": 0.35},
        facet_kws={"sharex": False, "sharey": False},
    )
    grid.set_axis_labels("RNA signature mean log2(TPM + 1)", "GigaTIME mean activation")
    grid.fig.suptitle("Top GigaTIME-RNA Signature Associations", y=1.03)
    grid.savefig(out_dir / "top_gigatime_rna_signature_scatter.png", dpi=180, bbox_inches="tight")
    plt.close(grid.fig)


def describe_tile_sampling(merged) -> str:
    if "n_tiles" not in merged.columns:
        return "Tile sampling information was not available in the joined table."
    tile_counts = merged.drop_duplicates("case_submitter_id")["n_tiles"].dropna()
    if tile_counts.empty:
        return "Tile sampling information was not available in the joined table."
    unique_counts = sorted({int(value) for value in tile_counts})
    if len(unique_counts) == 1:
        return f"Tile sampling: {unique_counts[0]} tissue tiles per slide."
    return (
        f"Tile sampling: median {int(tile_counts.median())} tissue tiles per slide "
        f"(range {int(tile_counts.min())}-{int(tile_counts.max())})."
    )


def write_markdown(path: Path, correlation, group_summary, merged) -> None:
    case_count = merged["case_submitter_id"].nunique()
    tile_sampling = describe_tile_sampling(merged)
    top_positive = correlation.sort_values("spearman_rho", ascending=False).head(5)
    top_abs = correlation.assign(abs_rho=correlation["spearman_rho"].abs()).sort_values("abs_rho", ascending=False).head(5)
    lines = [
        "# GigaTIME RNA Signature Validation",
        "",
        f"- Cases with paired GigaTIME and RNA-seq data: {case_count}",
        f"- Channels tested: {len(correlation)}",
        f"- {tile_sampling}",
        "- RNA values are simple marker-signature means of log2(TPM + 1).",
        "- This is an indirect validation check, not real mIF ground truth.",
        "",
        "## Strongest Positive Correlations",
        "",
        "| Channel | Spearman rho | p | BH q | RNA genes |",
        "|---|---:|---:|---:|---|",
    ]
    for _, row in top_positive.iterrows():
        lines.append(
            f"| {row['channel']} | {row['spearman_rho']:.3f} | {row['spearman_p_value']:.4g} | "
            f"{row['spearman_q_value_bh']:.4g} | {row['rna_signature_genes']} |"
        )
    lines.extend(
        [
            "",
            "## Strongest Absolute Correlations",
            "",
            "| Channel | Spearman rho | p | BH q | RNA genes |",
            "|---|---:|---:|---:|---|",
        ]
    )
    for _, row in top_abs.iterrows():
        lines.append(
            f"| {row['channel']} | {row['spearman_rho']:.3f} | {row['spearman_p_value']:.4g} | "
            f"{row['spearman_q_value_bh']:.4g} | {row['rna_signature_genes']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Guardrails",
            "",
            "- A positive correlation means cases with higher RNA marker expression also tended to have higher GigaTIME virtual-channel activation.",
            "- A weak or negative correlation does not automatically mean GigaTIME is wrong; bulk RNA-seq and H&E-tile virtual mIF measure different biological layers.",
            f"- The current analysis has only {case_count} cases. {tile_sampling}",
            "- These results should guide follow-up validation, not final claims.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd, plt, sns, stats = require_analysis_libs(out_dir / ".matplotlib")
    joined = pd.read_csv(args.joined)
    channels = [channel.strip() for channel in args.channels.split(",") if channel.strip()]
    channels = [channel for channel in channels if channel in CHANNEL_SIGNATURES]
    if not channels:
        raise SystemExit("No valid channels requested.")
    case_signatures, gene_expression = build_case_signature_table(
        pd, joined, Path(args.expression_dir), channels
    )
    correlation, merged = build_correlation_table(pd, stats, joined, case_signatures, channels)
    group_summary = build_group_summary(pd, merged, channels)

    case_signatures.to_csv(out_dir / "case_rna_signatures.csv", index=False)
    gene_expression.to_csv(out_dir / "case_gene_expression_long.csv", index=False)
    merged.to_csv(out_dir / "joined_gigatime_rna_signatures.csv", index=False)
    correlation.to_csv(out_dir / "gigatime_rna_signature_correlations.csv", index=False)
    group_summary.to_csv(out_dir / "gigatime_rna_group_summary.csv", index=False)
    plot_correlation_heatmap(plt, sns, correlation, out_dir)
    plot_top_scatter(plt, sns, merged, correlation, out_dir)
    write_markdown(out_dir / "rna_validation_summary.md", correlation, group_summary, merged)
    print(f"Wrote RNA validation outputs to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
