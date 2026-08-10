#!/usr/bin/env python
"""Train an isotopomer-distribution model for one HSQC vector.

End-to-end: simulate training data -> collate features/labels -> train the
fixed-architecture network -> save the model and a summary.

Example:
    python -m metabolabpytools.scripts.train --hsqc-vector 0 1 1
    mlp-train --hsqc-vector 0 1 1 --n-distributions 20000
"""
import argparse
import os

import numpy as np

from metabolabpytools import isotopomerAnalysis

# Fixed seed for data simulation so the training set (and the held-out test
# split) are identical across runs. --seed varies only the model training,
# which is what best-of-N model selection should explore.
DATA_SEED = 1234


def default_n_distributions(n_carbons):
    """More output classes (2**n_carbons) need more training samples."""
    if n_carbons <= 3:
        return 10000
    if n_carbons == 4:
        return 20000
    return 40000


def main():
    parser = argparse.ArgumentParser(description="Train an isotopomer NN for one HSQC vector.")
    parser.add_argument("--hsqc-vector", nargs="+", type=int, required=True,
                        help="HSQC vector, e.g. --hsqc-vector 0 1 1")
    parser.add_argument("--n-distributions", type=int, default=None,
                        help="Number of synthetic samples to simulate (default scales with carbon count).")
    parser.add_argument("--epochs", type=int, default=500, help="Max epochs (early stopping usually stops sooner).")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--seed", type=int, default=42,
                        help="Training seed (weight init/dropout). Same seed -> identical model. "
                             "Try different seeds and keep-best retains the lowest test-MAE one.")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite the saved model even if this run's test MAE is not better.")
    parser.add_argument("--out-dir", default="saved_models", help="Directory to save the trained model.")
    parser.add_argument("--summary-dir", default="model_summaries",
                        help="Directory for the model summary CSV (holds the best test MAE for keep-best).")
    parser.add_argument("--work-dir", default=".",
                        help="Directory to run in; sim_data/ and saved_models/ are created here.")
    parser.add_argument("--labeled-only", action="store_true",
                        help="Train on the labelled-only target: drop the unlabelled isotopomer and "
                             "renormalise the rest to 100%%. Improves precision of the labelled pattern "
                             "(most for 3-4 carbons). Predict with mlp-predict --labeled-only.")
    args = parser.parse_args()

    os.chdir(args.work_dir)

    hsqc_vector = args.hsqc_vector
    n_carbons = len(hsqc_vector)
    n_distributions = args.n_distributions or default_n_distributions(n_carbons)

    analysis = isotopomerAnalysis.IsotopomerAnalysisNN()

    # Deterministic simulation: fixed data seed -> identical dataset every run.
    np.random.seed(DATA_SEED)
    print(f"Simulating {n_distributions} distributions for HSQC vector {hsqc_vector} ({n_carbons} carbons)...")
    distributions = analysis.generate_isotopomer_distributions(n_distributions=n_distributions, n_carbons=n_carbons)
    isotopomer_data, hsqc_data, gcms_data = analysis.simulate_hsqc_gcms(distributions, hsqc_vector)

    # Keep a copy of the simulated data for inspection/reproducibility.
    analysis.save_simulation_data(isotopomer_data, hsqc_data, gcms_data, hsqc_vector)

    # Collate features (X) and labels (Y) directly from the simulated data.
    all_possible_hsqc_multiplets = analysis.generate_possible_hsqc_multiplets(hsqc_vector)
    Y = analysis.collate_y_labels(isotopomer_data, n_carbons, exclude_unlabeled=args.labeled_only)
    X = analysis.collate_x_labels_without_noise(hsqc_data, gcms_data, all_possible_hsqc_multiplets)

    if args.labeled_only:
        # Drop samples with no labelled content (their conditional target is
        # undefined and would be an all-zero, sum-to-100-violating row).
        keep = Y.sum(axis=1) > 0
        dropped = int((~keep).sum())
        if dropped:
            print(f"Labelled-only: dropping {dropped} all-unlabelled samples.")
        X, Y = X[keep], Y[keep]

    print(f"Training on X{X.shape} -> Y{Y.shape}  (seed={args.seed}) ...")
    _, _, info = analysis.train_and_keep_best(
        X, Y, hsqc_vector, labeled=args.labeled_only, seed=args.seed,
        epochs=args.epochs, batch_size=args.batch_size, force=args.force,
        save_dir=args.out_dir, summary_dir=args.summary_dir,
    )

    path = analysis.best_model_path(hsqc_vector, labeled=args.labeled_only, directory=args.out_dir)
    verb = "saved" if info["saved"] else "unchanged (existing model was better)"
    print(f"Done. Best model {verb}: {path}")


if __name__ == "__main__":
    main()
