# metabolabpytools
[![DOI](https://zenodo.org/badge/644108117.svg)](https://zenodo.org/doi/10.5281/zenodo.13342463)
Tools used for metabolism related data analysis

## Documentation

Two self-contained HTML handouts live in [`docs/`](docs/). **View them rendered
in your browser** (GitHub shows the raw source otherwise):

- **[▶ Quickstart](https://htmlpreview.github.io/?https://github.com/ludwigc/metabolabpytools/blob/main/docs/quickstart.html)** — how to install and run the train/predict workflow.
- **[▶ How it works](https://htmlpreview.github.io/?https://github.com/ludwigc/metabolabpytools/blob/main/docs/how-it-works.html)** — what the pipeline does and why accuracy plateaus for larger metabolites.

Or open the files locally after cloning:

```bash
# macOS
open docs/quickstart.html docs/how-it-works.html
# Linux
xdg-open docs/quickstart.html
# Windows
start docs\quickstart.html
```

Both are styled for print, so **⌘P / Ctrl+P → Save as PDF** turns them into
shareable handouts.

## Isotopomer analysis: train & predict

Predict isotopomer distributions from HSQC + GC-MS data with a small neural
network, one model per HSQC vector. There are two ways to run the workflow.

### Notebooks (walkthrough)

In `metabolabpytools/jupyter/`, run top-to-bottom:

1. `01_train.ipynb` — simulate training data, train a model, save it.
2. `02_predict.ipynb` — load real data, predict a distribution, save results.

### Command line (one command each)

```bash
# Train a model for an HSQC vector (samples scale with carbon count automatically)
mlp-train --hsqc-vector 0 1 1

# Predict on measured data (auto-loads the current best model for the vector)
mlp-predict --hsqc-vector 0 1 1 \
    --hsqc-data hsqcData1.xlsx --gcms-data gcmsData1.xlsx --metabolite L-LacticAcid
```

Both write to the current directory: models to `saved_models/`, a summary to
`model_summaries/`, and predictions to `nn_analysis_results/`. Pass `--work-dir`
to run elsewhere (e.g. `--work-dir metabolabpytools/jupyter`).

### Reproducibility and keep-best

Training is **deterministic**: the dataset is simulated from a fixed seed, and
`--seed` (default 42) fixes the training run, so the same command always gives
the same model. Each run evaluates on a **held-out test set** and only
overwrites the saved model if the test MAE improves — so the saved model is
always the best one found, and `mlp-predict` always loads that best model.

```bash
# Search for a better model by sweeping seeds; the best is kept automatically
for s in 1 2 3 4 5; do mlp-train --hsqc-vector 0 1 1 --seed $s; done

mlp-train --hsqc-vector 0 1 1 --seed 7 --force   # overwrite even if not better
```

Training uses a fixed, fast architecture (`Dense(128) → Dense(64)`), so there is
no hyperparameter search to wait on. The optional `IsotopomerAnalysisNN.tune_model`
Bayesian search is still available if you want to explore architectures.

### Labelled-only models (optional, higher precision on the labelled pattern)

For larger metabolites the unlabelled isotopomer dominates the distribution
(often 60–90%), and fitting it alongside the small labelled species spends most
of the model's error budget on that one large term. Passing `--labeled-only`
trains on the labelled species only — the unlabelled component is dropped and
the remaining isotopomers are renormalised to 100%, giving the labelling
*pattern* the full dynamic range. In tests this improves precision on the
labelled isotopomers by ~10–20%. (It does not overcome the underdetermination
limit for 4+ carbons — see *How it works* — it just spends precision where it
matters.)

In a 5-seed test this improved the labelling *pattern* (ratios among labelled
species) by 11–27% and absolute labelled percentages by 4–12%, consistently
across 3–5 carbons.

```bash
# Train a labelled-only model (saved as ..._labeled.keras, alongside any baseline)
mlp-train --hsqc-vector 0 1 1 --labeled-only

# Predict with it: auto-loads the ..._labeled model; the unlabelled fraction
# is taken from the measured GC-MS M+0
mlp-predict --hsqc-vector 0 1 1 --labeled-only \
    --hsqc-data hsqcData1.xlsx --gcms-data gcmsData1.xlsx --metabolite L-LacticAcid
```

Labelled-only models are trained on demand (not committed to the repo); run the
command above to create one for your vector.
