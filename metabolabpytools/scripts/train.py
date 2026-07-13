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

from metabolabpytools import isotopomerAnalysis


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
    parser.add_argument("--out-dir", default="saved_models", help="Directory to save the trained model.")
    parser.add_argument("--work-dir", default=".",
                        help="Directory to run in; sim_data/ and saved_models/ are created here.")
    args = parser.parse_args()

    os.chdir(args.work_dir)

    hsqc_vector = args.hsqc_vector
    n_carbons = len(hsqc_vector)
    n_distributions = args.n_distributions or default_n_distributions(n_carbons)

    analysis = isotopomerAnalysis.IsotopomerAnalysisNN()

    print(f"Simulating {n_distributions} distributions for HSQC vector {hsqc_vector} ({n_carbons} carbons)...")
    distributions = analysis.generate_isotopomer_distributions(n_distributions=n_distributions, n_carbons=n_carbons)
    isotopomer_data, hsqc_data, gcms_data = analysis.simulate_hsqc_gcms(distributions, hsqc_vector)

    # Keep a copy of the simulated data for inspection/reproducibility.
    analysis.save_simulation_data(isotopomer_data, hsqc_data, gcms_data, hsqc_vector)

    # Collate features (X) and labels (Y) directly from the simulated data.
    all_possible_hsqc_multiplets = analysis.generate_possible_hsqc_multiplets(hsqc_vector)
    Y = analysis.collate_y_labels(isotopomer_data, n_carbons)
    X = analysis.collate_x_labels_without_noise(hsqc_data, gcms_data, all_possible_hsqc_multiplets)

    print(f"Training on X{X.shape} -> Y{Y.shape} ...")
    analysis.train_neural_network(
        X, Y, hsqc_vector=hsqc_vector, epochs=args.epochs, batch_size=args.batch_size,
        save=True, save_dir=args.out_dir,
    )

    model_name = analysis.generate_model_filename(hsqc_vector)
    print(f"Done. Model saved to {os.path.join(args.out_dir, model_name)}")


if __name__ == "__main__":
    main()
