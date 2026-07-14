# metabolabpytools
[![DOI](https://zenodo.org/badge/644108117.svg)](https://zenodo.org/doi/10.5281/zenodo.13342463)
Tools used for metabolism related data analysis

## Documentation

Two self-contained HTML handouts live in [`docs/`](docs/) — open either in a
browser (no build step, no dependencies):

- **[Quickstart](docs/quickstart.html)** — how to install and run the train/predict workflow.
- **[How it works](docs/how-it-works.html)** — what the pipeline does and why accuracy plateaus for larger metabolites.

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

# Predict on measured data
mlp-predict --model saved_models/model_hsqc_0_1_1.keras --hsqc-vector 0 1 1 \
    --hsqc-data hsqcData1.xlsx --gcms-data gcmsData1.xlsx --metabolite L-LacticAcid
```

Both write to the current directory: models to `saved_models/`, a summary to
`model_summaries/`, and predictions to `nn_analysis_results/`. Pass `--work-dir`
to run elsewhere (e.g. `--work-dir metabolabpytools/jupyter`).

Training uses a fixed, fast architecture (`Dense(128) → Dense(64)`), so there is
no hyperparameter search to wait on. The optional `IsotopomerAnalysisNN.tune_model`
Bayesian search is still available if you want to explore architectures.
