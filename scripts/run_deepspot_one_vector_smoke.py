#!/usr/bin/env python3
"""Run a minimal DeepSpot smoke test from one feature vector."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import torch
import yaml


MARKER_GENES = [
    "ERBB2",
    "ESR1",
    "PGR",
    "MKI67",
    "EPCAM",
    "KRT8",
    "KRT18",
    "PTPRC",
    "CD3D",
    "CD8A",
    "MS4A1",
    "COL1A1",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deepspot-repo", default="external/DeepSpotTutorial")
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--genes-csv", required=True)
    parser.add_argument("--config", default="")
    parser.add_argument(
        "--embedding-csv",
        default="results/phikon_one_tile_smoke/tile_embedding.csv",
        help="Optional one-vector CSV. Uses it only when its dimension matches the checkpoint input size.",
    )
    parser.add_argument("--tile-summary", default="results/phikon_one_tile_smoke/phikon_one_tile_summary.json")
    parser.add_argument("--input-source", choices=["auto", "embedding", "synthetic"], default="auto")
    parser.add_argument("--synthetic-seed", type=int, default=2026)
    parser.add_argument("--out-dir", default="results/deepspot_one_vector_smoke")
    parser.add_argument("--top-n", type=int, default=25)
    parser.add_argument("--device", choices=["cpu", "mps", "cuda"], default="cpu")
    return parser.parse_args()


def load_yaml(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def load_tile_summary(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def load_gene_names(path: Path, expected_size: int) -> list[str]:
    genes: list[str] = []
    predicted: list[str] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            gene_name = row.get("gene_name") or row.get("gene") or ""
            if not gene_name:
                continue
            genes.append(gene_name)
            if str(row.get("isPredicted", "")).lower() == "true":
                predicted.append(gene_name)
    if len(predicted) == expected_size:
        return predicted
    if len(genes) == expected_size:
        return genes
    return [f"gene_{idx:04d}" for idx in range(expected_size)]


def load_embedding_csv(path: Path) -> tuple[np.ndarray | None, dict]:
    if not path.exists():
        return None, {"path": str(path), "status": "missing"}
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        return None, {"path": str(path), "status": "empty"}

    if {"feature_index", "value"}.issubset(rows[0]):
        ordered = sorted(rows, key=lambda row: int(row["feature_index"]))
        values = np.array([float(row["value"]) for row in ordered], dtype=np.float32)
        return values, {"path": str(path), "status": "loaded", "format": "feature_index_value", "dimension": int(values.size)}

    embedding_keys = sorted(key for key in rows[0] if key.startswith("embedding_"))
    if embedding_keys:
        values = np.array([float(rows[0][key]) for key in embedding_keys], dtype=np.float32)
        return values, {"path": str(path), "status": "loaded", "format": "embedding_columns_first_row", "dimension": int(values.size)}

    return None, {"path": str(path), "status": "unsupported_format", "columns": list(rows[0])}


def resolve_input_vector(args: argparse.Namespace, input_size: int) -> tuple[np.ndarray, dict]:
    embedding, metadata = load_embedding_csv(Path(args.embedding_csv))
    if args.input_source in {"auto", "embedding"} and embedding is not None:
        if embedding.size == input_size:
            metadata["used"] = True
            return embedding.astype(np.float32), metadata
        metadata["used"] = False
        metadata["reason"] = f"embedding dimension {embedding.size} does not match checkpoint input size {input_size}"
        if args.input_source == "embedding":
            raise ValueError(metadata["reason"])

    rng = np.random.default_rng(args.synthetic_seed)
    vector = rng.standard_normal(input_size).astype(np.float32)
    synthetic_metadata = {
        "status": "synthetic",
        "used": True,
        "seed": args.synthetic_seed,
        "dimension": int(input_size),
    }
    if metadata:
        synthetic_metadata["embedding_probe"] = metadata
    return vector, synthetic_metadata


def context_from_vector(vector: np.ndarray, spot_context: str, device: torch.device):
    tensor = torch.from_numpy(vector).to(device=device, dtype=torch.float32).reshape(1, 1, -1)
    if spot_context == "spot":
        return tensor
    if spot_context == "spot_subspot":
        return [tensor, tensor]
    if spot_context == "spot_neighbors":
        return [tensor, tensor]
    if spot_context == "spot_subspot_neighbors":
        return [tensor, tensor, tensor]
    raise ValueError(f"Unsupported DeepSpot spot_context: {spot_context}")


def stats(values: np.ndarray) -> dict[str, float]:
    values = values.astype(np.float64, copy=False)
    return {
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
    }


def prediction_rows(genes: list[str], raw: np.ndarray, inverse: np.ndarray | None) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for idx, gene in enumerate(genes):
        row: dict[str, object] = {
            "rank_input_order": idx,
            "gene_name": gene,
            "raw_prediction": float(raw[idx]),
        }
        if inverse is not None:
            row["inverse_prediction"] = float(inverse[idx])
        rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    deepspot_repo = Path(args.deepspot_repo)
    if str(deepspot_repo) not in sys.path:
        sys.path.insert(0, str(deepspot_repo))

    device = torch.device(args.device)
    with warnings.catch_warnings(record=True) as caught_warnings:
        warnings.simplefilter("always")
        model = torch.load(args.checkpoint, map_location=device, weights_only=False)
    model.to(device)
    model.eval()

    input_size = int(model.hparams.input_size)
    output_size = int(model.hparams.output_size)
    spot_context = str(model.hparams.spot_context)

    vector, input_metadata = resolve_input_vector(args, input_size)
    x = context_from_vector(vector, spot_context, device)
    with torch.inference_mode():
        raw_prediction = model(x).detach().cpu().float().numpy()[0]

    inverse_prediction = None
    inverse_error = ""
    try:
        inverse_prediction = np.asarray(model.inverse_transform(raw_prediction.reshape(1, -1))[0], dtype=np.float32)
    except Exception as exc:  # noqa: BLE001 - capture scaler/version issues in the smoke summary.
        inverse_error = repr(exc)

    genes = load_gene_names(Path(args.genes_csv), output_size)
    rows = prediction_rows(genes, raw_prediction, inverse_prediction)
    score_key = "inverse_prediction" if inverse_prediction is not None else "raw_prediction"
    top_rows = sorted(rows, key=lambda row: float(row[score_key]), reverse=True)[: args.top_n]
    bottom_rows = sorted(rows, key=lambda row: float(row[score_key]))[: args.top_n]
    marker_set = set(MARKER_GENES)
    marker_rows = [row for row in rows if str(row["gene_name"]) in marker_set]

    out_dir = Path(args.out_dir)
    write_csv(out_dir / "all_gene_predictions.csv", rows)
    write_csv(out_dir / "top_predicted_genes.csv", top_rows)
    write_csv(out_dir / "bottom_predicted_genes.csv", bottom_rows)
    write_csv(out_dir / "selected_marker_predictions.csv", marker_rows)

    summary = {
        "model": "DeepSpot one-vector smoke",
        "checkpoint": str(args.checkpoint),
        "config": load_yaml(Path(args.config)) if args.config else {},
        "checkpoint_hparams": {
            "input_size": input_size,
            "output_size": output_size,
            "spot_context": spot_context,
            "loss_func": str(model.hparams.loss_func),
        },
        "input": input_metadata,
        "tile_context": load_tile_summary(Path(args.tile_summary)),
        "device": str(device),
        "n_genes": len(rows),
        "raw_prediction_stats": stats(raw_prediction),
        "inverse_prediction_stats": stats(inverse_prediction) if inverse_prediction is not None else None,
        "inverse_transform_error": inverse_error,
        "warnings": [str(warning.message) for warning in caught_warnings],
        "top_gene_names": [row["gene_name"] for row in top_rows],
        "marker_gene_count": len(marker_rows),
        "outputs": {
            "all_gene_predictions": str(out_dir / "all_gene_predictions.csv"),
            "top_predicted_genes": str(out_dir / "top_predicted_genes.csv"),
            "bottom_predicted_genes": str(out_dir / "bottom_predicted_genes.csv"),
            "selected_marker_predictions": str(out_dir / "selected_marker_predictions.csv"),
            "summary_json": str(out_dir / "deepspot_one_vector_summary.json"),
        },
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "deepspot_one_vector_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
