#!/usr/bin/env python3
"""Extract patient-level BCNB embeddings from precomputed patch zip members."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import zipfile
from collections import Counter, defaultdict
from io import BytesIO
from pathlib import Path

import numpy as np

from run_hoptimus_tcga_brca import (
    autocast_context,
    features_from_output,
    import_runtime,
    load_model,
    make_transform as make_hoptimus_transform,
    read_existing_csv,
    resolve_device,
    selected_embedding_mode,
    tissue_fraction,
    write_csv,
)
from run_virchow2_one_slide_smoke import (
    load_virchow2,
    make_transform as make_virchow2_transform,
    virchow2_embedding,
)


DEFAULT_MANIFEST = Path("data/bcnb/bcnb_patch_manifest_capped10.csv")
DEFAULT_PATCH_ZIP = Path("data/bcnb/paper_patches.zip")

MODEL_CHOICES = {"hoptimus0", "h0-mini", "virchow2"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patch-manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--patch-zip", type=Path, default=DEFAULT_PATCH_ZIP)
    parser.add_argument("--out-dir", type=Path, default=Path("results/bcnb_patch_embeddings_hoptimus0_capped10"))
    parser.add_argument("--model", choices=sorted(MODEL_CHOICES), default="hoptimus0")
    parser.add_argument("--model-id", default=None, help="Override default Hugging Face/timm model id.")
    parser.add_argument("--embedding-mode", choices=["auto", "direct", "cls", "concat_cls_patch_mean", "patch_mean"], default="auto")
    parser.add_argument("--input-size", type=int, default=224)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--device", choices=["auto", "cuda", "mps", "cpu"], default="auto")
    parser.add_argument("--precision", choices=["auto", "float32", "float16"], default="auto")
    parser.add_argument("--tissue-threshold", type=float, default=0.0, help="Skip patches below this tissue fraction.")
    parser.add_argument(
        "--groups",
        default="HER2-zero,HER2-low,HER2-positive",
        help="Comma-separated clinical HER2 groups to include.",
    )
    parser.add_argument("--max-patients", type=int, default=0, help="Maximum patients to process after filtering. Use 0 for all.")
    parser.add_argument(
        "--max-patients-per-group",
        type=int,
        default=0,
        help="Optional deterministic cap per HER2 group, useful for balanced smokes.",
    )
    parser.add_argument("--max-patches-per-patient", type=int, default=0, help="Optional cap applied after manifest order.")
    parser.add_argument("--save-patch-csv", action="store_true", help="Write per-patch embeddings. Use only for small smokes.")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="Read manifest/zip and report planned work without loading a model.")
    return parser.parse_args()


def parse_groups(raw: str) -> set[str]:
    return {item.strip() for item in raw.split(",") if item.strip()}


def read_manifest(path: Path, groups: set[str], max_patches_per_patient: int) -> dict[str, list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Patch manifest not found: {path}")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"patient_id", "clinical_her2_group", "patch_zip_member"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: {', '.join(sorted(required))}")
        for row in reader:
            if row["clinical_her2_group"] not in groups:
                continue
            patient_id = row["patient_id"]
            if max_patches_per_patient and len(grouped[patient_id]) >= max_patches_per_patient:
                continue
            grouped[patient_id].append(row)
    return dict(sorted(grouped.items(), key=lambda item: int(item[0])))


def first_patient_metadata(rows: list[dict[str, str]]) -> dict[str, object]:
    first = rows[0]
    return {
        "patient_id": first["patient_id"],
        "clinical_her2_group": first.get("clinical_her2_group", ""),
        "her2_status": first.get("her2_status", ""),
        "her2_ihc": first.get("her2_ihc", ""),
        "grade": first.get("grade", ""),
        "ER": first.get("ER", ""),
        "PR": first.get("PR", ""),
        "ki67": first.get("ki67", ""),
        "molecular_subtype": first.get("molecular_subtype", ""),
        "aln_status": first.get("aln_status", ""),
    }


def load_patch_rgb(archive: zipfile.ZipFile, member: str, Image) -> np.ndarray:
    with archive.open(member) as handle:
        payload = handle.read()
    image = Image.open(BytesIO(payload)).convert("RGB")
    return np.asarray(image)


def make_model(args: argparse.Namespace, torch, timm, Image, resolve_data_config, create_transform, device):
    if args.model == "virchow2":
        model_id = args.model_id or "hf-hub:paige-ai/Virchow2"
        model = load_virchow2(torch, timm, model_id, device)
        transform = make_virchow2_transform(model, args.input_size, Image, resolve_data_config, create_transform)
        return model, transform, model_id, "concat(class_token, mean(patch_tokens_after_4_register_tokens))"

    model, model_id = load_model(torch, timm, args.model, args.model_id, device)
    transform = make_hoptimus_transform(args.model, model, args.input_size, Image, resolve_data_config, create_transform)
    embedding_mode = selected_embedding_mode(args.model, args.embedding_mode)
    return model, transform, model_id, embedding_mode


def infer_batch(torch, model, model_name: str, transform, rgb_batch, meta_batch, device, precision: str, embedding_mode: str):
    if model_name == "virchow2":
        tensor = torch.stack([transform(rgb) for rgb in rgb_batch], dim=0).to(device)
        with autocast_context(torch, device, precision):
            output = model(tensor)
            features = virchow2_embedding(torch, output)
    else:
        tensor = torch.stack([torch.from_numpy(transform(rgb)) for rgb in rgb_batch], dim=0).to(device)
        with autocast_context(torch, device, precision):
            output = model(tensor)
            features = features_from_output(torch, output, embedding_mode)

    features_np = features.detach().cpu().float().numpy()
    rows: list[dict[str, object]] = []
    for idx, meta in enumerate(meta_batch):
        row = dict(meta)
        for feature_idx, value in enumerate(features_np[idx]):
            row[f"embedding_{feature_idx:04d}"] = float(value)
        rows.append(row)
    return rows


def summarize_patient(patient_id: str, manifest_rows: list[dict[str, str]], patch_rows: list[dict[str, object]]) -> dict[str, object]:
    row = first_patient_metadata(manifest_rows)
    row["n_manifest_patches"] = len(manifest_rows)
    row["n_used_patches"] = len(patch_rows)
    row["n_skipped_low_tissue"] = len(manifest_rows) - len(patch_rows)
    if not patch_rows:
        return row
    row["mean_tissue_fraction"] = float(np.mean([float(patch["tissue_fraction"]) for patch in patch_rows]))
    row["min_tissue_fraction"] = float(np.min([float(patch["tissue_fraction"]) for patch in patch_rows]))
    row["max_tissue_fraction"] = float(np.max([float(patch["tissue_fraction"]) for patch in patch_rows]))
    embedding_cols = [key for key in patch_rows[0] if key.startswith("embedding_")]
    for key in embedding_cols:
        row[key] = float(np.mean([float(patch[key]) for patch in patch_rows]))
    row["patient_id"] = patient_id
    return row


def process_patient(
    *,
    torch,
    model,
    model_name: str,
    transform,
    archive: zipfile.ZipFile,
    Image,
    patient_id: str,
    manifest_rows: list[dict[str, str]],
    args: argparse.Namespace,
    device,
    embedding_mode: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    rgb_batch = []
    meta_batch = []
    patch_rows: list[dict[str, object]] = []
    with torch.inference_mode():
        for manifest_row in manifest_rows:
            rgb = load_patch_rgb(archive, manifest_row["patch_zip_member"], Image)
            frac = tissue_fraction(rgb)
            if frac < args.tissue_threshold:
                continue
            rgb_batch.append(rgb)
            meta_batch.append(
                {
                    "patient_id": patient_id,
                    "clinical_her2_group": manifest_row.get("clinical_her2_group", ""),
                    "patch_zip_member": manifest_row["patch_zip_member"],
                    "patch_filename": manifest_row.get("patch_filename", ""),
                    "tissue_fraction": frac,
                }
            )
            if len(rgb_batch) == args.batch_size:
                patch_rows.extend(
                    infer_batch(torch, model, model_name, transform, rgb_batch, meta_batch, device, args.precision, embedding_mode)
                )
                rgb_batch = []
                meta_batch = []
        if rgb_batch:
            patch_rows.extend(
                infer_batch(torch, model, model_name, transform, rgb_batch, meta_batch, device, args.precision, embedding_mode)
            )
    return summarize_patient(patient_id, manifest_rows, patch_rows), patch_rows


def planned_summary(grouped: dict[str, list[dict[str, str]]]) -> dict[str, object]:
    group_counts: Counter[str] = Counter()
    patch_counts: Counter[str] = Counter()
    for rows in grouped.values():
        group = rows[0]["clinical_her2_group"]
        group_counts[group] += 1
        patch_counts[group] += len(rows)
    return {
        "n_patients": len(grouped),
        "n_manifest_patches": sum(len(rows) for rows in grouped.values()),
        "patients_by_group": dict(sorted(group_counts.items())),
        "patches_by_group": dict(sorted(patch_counts.items())),
    }


def limit_patients(
    grouped: dict[str, list[dict[str, str]]],
    max_patients: int,
    max_patients_per_group: int,
) -> dict[str, list[dict[str, str]]]:
    if max_patients_per_group:
        group_counts: Counter[str] = Counter()
        limited = {}
        for patient_id, rows in grouped.items():
            group = rows[0]["clinical_her2_group"]
            if group_counts[group] >= max_patients_per_group:
                continue
            limited[patient_id] = rows
            group_counts[group] += 1
        grouped = limited
    if max_patients:
        grouped = dict(list(grouped.items())[:max_patients])
    return grouped


def main() -> int:
    args = parse_args()
    grouped = read_manifest(args.patch_manifest, parse_groups(args.groups), args.max_patches_per_patient)
    grouped = limit_patients(grouped, args.max_patients, args.max_patients_per_group)
    if not grouped:
        raise SystemExit("No patients selected from patch manifest.")
    if not args.patch_zip.exists():
        raise FileNotFoundError(f"Patch zip not found: {args.patch_zip}")

    if args.dry_run:
        print(json.dumps(planned_summary(grouped), indent=2))
        return 0

    if args.model in {"hoptimus0", "h0-mini"} and not os.environ.get("HF_TOKEN") and not os.environ.get("HUGGING_FACE_HUB_TOKEN"):
        print("HF_TOKEN/HUGGING_FACE_HUB_TOKEN is not set; gated H-Optimus download may fail.", file=sys.stderr)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    patient_embeddings_path = args.out_dir / "patient_embeddings.csv"
    patch_embeddings_path = args.out_dir / "patch_embeddings.csv"
    patient_rows = read_existing_csv(patient_embeddings_path) if args.resume else []
    all_patch_rows = read_existing_csv(patch_embeddings_path) if args.resume and args.save_patch_csv else []
    processed = {str(row.get("patient_id", "")) for row in patient_rows if row.get("patient_id")}

    torch, timm, Image, _openslide, resolve_data_config, create_transform = import_runtime()
    device = resolve_device(torch, args.device)
    print(f"Using device: {device}", file=sys.stderr)
    model, transform, model_id, embedding_mode = make_model(args, torch, timm, Image, resolve_data_config, create_transform, device)

    with zipfile.ZipFile(args.patch_zip) as archive:
        for index, (patient_id, rows) in enumerate(grouped.items(), start=1):
            if args.resume and patient_id in processed:
                print(f"[{index}/{len(grouped)}] Skipping existing patient {patient_id}", file=sys.stderr)
                continue
            print(f"[{index}/{len(grouped)}] Processing patient {patient_id} ({len(rows)} patches)", file=sys.stderr)
            patient_row, patch_rows = process_patient(
                torch=torch,
                model=model,
                model_name=args.model,
                transform=transform,
                archive=archive,
                Image=Image,
                patient_id=patient_id,
                manifest_rows=rows,
                args=args,
                device=device,
                embedding_mode=embedding_mode,
            )
            patient_rows.append(patient_row)
            all_patch_rows.extend(patch_rows)
            write_csv(patient_embeddings_path, patient_rows)
            if args.save_patch_csv:
                write_csv(patch_embeddings_path, all_patch_rows)

    embedding_dim = len([key for key in patient_rows[-1] if key.startswith("embedding_")]) if patient_rows else 0
    summary = {
        "task": "bcnb_patch_patient_embeddings",
        "model": args.model,
        "model_id": model_id,
        "embedding_mode": embedding_mode,
        "input_size": args.input_size,
        "patch_manifest": str(args.patch_manifest),
        "patch_zip": str(args.patch_zip),
        "groups": sorted(parse_groups(args.groups)),
        "tissue_threshold": args.tissue_threshold,
        "max_patients": args.max_patients,
        "max_patients_per_group": args.max_patients_per_group,
        "max_patches_per_patient": args.max_patches_per_patient,
        "n_patients": len(patient_rows),
        "embedding_dimensions": embedding_dim,
        "planned": planned_summary(grouped),
        "outputs": {
            "patient_embeddings": str(patient_embeddings_path),
            "patch_embeddings": str(patch_embeddings_path) if args.save_patch_csv else "",
        },
    }
    (args.out_dir / "bcnb_patch_embedding_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Done. Wrote {patient_embeddings_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
