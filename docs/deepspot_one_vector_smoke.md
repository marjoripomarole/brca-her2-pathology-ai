# DeepSpot one-vector smoke test

This smoke test checks that a published DeepSpot pretrained checkpoint can be
loaded locally and can produce a 5,000-gene prediction vector.

There are now two successful smoke levels:

- A synthetic-vector machinery check, used when only the 768-dimensional Phikon
  TCGA tile probe was available.
- A real TCGA one-tile check using a 1536-dimensional `bioptimus/H-optimus-0`
  embedding from `TCGA-A7-A26J`.

This is still not a biological TCGA-BRCA DeepSpot result. It is a one-tile
compatibility test using a lung Visium pretrained DeepSpot checkpoint.

## Local setup

The current setup uses the newer DeepSpot tutorial package:

```bash
git clone https://github.com/ratschlab/DeepSpot.git external/DeepSpotTutorial
```

The pretrained checkpoint archive was downloaded from Zenodo:

```bash
curl -L --fail -C - \
  -o /private/tmp/DeepSpot_pretrained_model_weights.zip \
  https://zenodo.org/api/records/15322099/files/DeepSpot_pretrained_model_weights.zip/content
```

For the first smoke, extract only the lung Visium spot-level model:

```bash
mkdir -p /private/tmp/deepspot_weights_lung_visium
unzip -o -j /private/tmp/DeepSpot_pretrained_model_weights.zip \
  DeepSpot_pretrained_model_weights/Lung_LUSC_LUAD_Visium/final_model.pkl \
  DeepSpot_pretrained_model_weights/Lung_LUSC_LUAD_Visium/top_param_overall.yaml \
  DeepSpot_pretrained_model_weights/Lung_LUSC_LUAD_Visium/info_highly_variable_genes.csv \
  -d /private/tmp/deepspot_weights_lung_visium
```

## Synthetic fallback command

```bash
conda run -n gigatime-tcga python scripts/run_deepspot_one_vector_smoke.py \
  --checkpoint /private/tmp/deepspot_weights_lung_visium/final_model.pkl \
  --genes-csv /private/tmp/deepspot_weights_lung_visium/info_highly_variable_genes.csv \
  --config /private/tmp/deepspot_weights_lung_visium/top_param_overall.yaml \
  --out-dir results/deepspot_one_vector_smoke \
  --device cpu \
  --top-n 25
```

## Expected behavior

By default, the script first probes
`results/phikon_one_tile_smoke/tile_embedding.csv`. If an H-Optimus tile CSV
with 1536 embedding dimensions is passed with `--embedding-csv`, that vector
will be used. With the Phikon CSV, the dimension mismatch is recorded and the
script uses a deterministic synthetic vector instead.

## Initial smoke-test output

The initial run used:

- checkpoint: `Lung_LUSC_LUAD_Visium/final_model.pkl`
- image feature model expected by checkpoint: `hoptimus0`
- input source: deterministic synthetic vector, seed `2026`
- checkpoint input size: `1536`
- output size: `5000` predicted gene values
- Phikon tile probe: loaded from `results/phikon_one_tile_smoke/tile_embedding.csv`, but not used because it has `768` dimensions
- raw prediction mean: `-0.0903433936`
- inverse-transformed prediction mean: `0.1467310426`
- top inverse-transformed genes: `SFTPC`, `FTL`, `SFTPA1`, `SFTPA2`, `S100A6`
- ERBB2 inverse-transformed value: `0.3198132515`

Ignored result files are written under:

- `results/deepspot_one_vector_smoke/all_gene_predictions.csv`
- `results/deepspot_one_vector_smoke/top_predicted_genes.csv`
- `results/deepspot_one_vector_smoke/bottom_predicted_genes.csv`
- `results/deepspot_one_vector_smoke/selected_marker_predictions.csv`
- `results/deepspot_one_vector_smoke/deepspot_one_vector_summary.json`

## Real H-Optimus one-tile smoke

The next smoke generated one H-Optimus-0 tile embedding from the same TCGA tile
context, then passed that real 1536-dimensional vector into DeepSpot.

H-Optimus-0 command:

```bash
conda run -n gigatime-tcga python scripts/run_hoptimus_tcga_brca.py \
  --model-preset hoptimus0 \
  --max-slides 1 \
  --tile-limit 1 \
  --save-tile-csv \
  --out-dir results/hoptimus_one_tile_for_deepspot \
  --device auto \
  --batch-size 1
```

DeepSpot command:

```bash
conda run -n gigatime-tcga python scripts/run_deepspot_one_vector_smoke.py \
  --checkpoint /private/tmp/deepspot_weights_lung_visium/final_model.pkl \
  --genes-csv /private/tmp/deepspot_weights_lung_visium/info_highly_variable_genes.csv \
  --config /private/tmp/deepspot_weights_lung_visium/top_param_overall.yaml \
  --embedding-csv results/hoptimus_one_tile_for_deepspot/tile_embeddings.csv \
  --tile-summary results/hoptimus_one_tile_for_deepspot/hoptimus_embedding_summary.json \
  --input-source embedding \
  --out-dir results/deepspot_hoptimus_one_tile_smoke \
  --device cpu \
  --top-n 25
```

The real-vector run used:

- case: `TCGA-A7-A26J`
- slide: `TCGA-A7-A26J-01B-02-BS2.2BDFB544-F62A-402C-9D97-DE2B6766DEDC`
- tile coordinates: `x=4961`, `y=2255`
- tile read size: `451 x 451` level-0 pixels
- MPP source: OpenSlide, `0.2485 x 0.2485`
- tissue fraction: `0.7153419962`
- H-Optimus model: `hf-hub:bioptimus/H-optimus-0`
- H-Optimus embedding dimension: `1536`
- DeepSpot checkpoint input size: `1536`
- DeepSpot output size: `5000` predicted gene values
- raw prediction mean: `0.0305611282`
- inverse-transformed prediction mean: `0.2215564781`
- top inverse-transformed genes: `FTL`, `RPL41`, `RPLP1`, `FTH1`, `S100A6`
- ERBB2 inverse-transformed value: `0.3993973732`

Selected marker outputs from the real-vector run:

| Gene | Inverse prediction |
|---|---:|
| ERBB2 | 0.3993973732 |
| EPCAM | 0.7018679380 |
| KRT8 | 1.2504146099 |
| KRT18 | 1.1665616035 |
| PTPRC | 0.2462048829 |
| CD3D | 0.0989107341 |
| CD8A | 0.0332667157 |
| MS4A1 | 0.0650566295 |
| ESR1 | 0.0058887554 |
| MKI67 | 0.0606881008 |
| COL1A1 | 0.7971688509 |

Ignored result files are written under:

- `results/hoptimus_one_tile_for_deepspot/slide_embeddings.csv`
- `results/hoptimus_one_tile_for_deepspot/tile_embeddings.csv`
- `results/hoptimus_one_tile_for_deepspot/hoptimus_embedding_summary.json`
- `results/deepspot_hoptimus_one_tile_smoke/all_gene_predictions.csv`
- `results/deepspot_hoptimus_one_tile_smoke/top_predicted_genes.csv`
- `results/deepspot_hoptimus_one_tile_smoke/bottom_predicted_genes.csv`
- `results/deepspot_hoptimus_one_tile_smoke/selected_marker_predictions.csv`
- `results/deepspot_hoptimus_one_tile_smoke/deepspot_one_vector_summary.json`

## Next step

The next useful DeepSpot experiment is not a full cohort run. Run a tiny,
balanced multi-tile smoke first, for example 1-3 tiles from one HER2-low and one
HER2-zero case, then inspect whether the predicted marker profiles are stable
enough to justify a larger virtual transcriptomics experiment.
