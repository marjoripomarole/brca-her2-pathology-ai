# DeepSpot one-vector smoke test

This smoke test checks that a published DeepSpot pretrained checkpoint can be
loaded locally and can produce a 5,000-gene prediction vector.

It is not yet a true TCGA-BRCA tile prediction. The downloaded pretrained
DeepSpot spot checkpoints use H-Optimus embeddings with input size `1536`,
whereas the current available TCGA tile smoke vector is the open Phikon model
with input size `768`. Until `bioptimus/H-optimus-0` weights are downloaded and
a real H-Optimus tile embedding is produced, the script falls back to a
deterministic synthetic 1536-dimensional vector and records that fallback in the
summary JSON.

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

## Command

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

The script first probes `results/phikon_one_tile_smoke/tile_embedding.csv`. If a
future H-Optimus tile CSV with 1536 embedding dimensions is passed, that vector
will be used. With the current Phikon CSV, the dimension mismatch is recorded
and the script uses a deterministic synthetic vector instead.

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

## Next step

For a true raw-tile DeepSpot smoke on `TCGA-A7-A26J`, download or otherwise
provide `bioptimus/H-optimus-0`, run one tile through
`scripts/run_hoptimus_tcga_brca.py --model-preset hoptimus0 --save-tile-csv`,
then pass that 1536-dimensional tile CSV to this smoke script with
`--embedding-csv`.
