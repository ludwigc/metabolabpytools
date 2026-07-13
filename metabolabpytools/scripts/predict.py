#!/usr/bin/env python
"""Predict an isotopomer distribution from real HSQC + GC-MS data.

End-to-end: load a trained model -> build feature vectors from the measured
data -> Monte-Carlo-dropout prediction -> back-simulate HSQC/GC-MS from the
prediction -> save a results summary.

Example:
    python -m metabolabpytools.scripts.predict \\
        --model saved_models/model_hsqc_0_1_1.keras --hsqc-vector 0 1 1 \\
        --hsqc-data hsqcData1.xlsx --gcms-data gcmsData1.xlsx --metabolite L-LacticAcid
"""
import argparse
import os

from metabolabpytools import isotopomerAnalysis


def main():
    parser = argparse.ArgumentParser(description="Predict an isotopomer distribution from measured data.")
    parser.add_argument("--model", required=True, help="Path to a trained .keras model.")
    parser.add_argument("--hsqc-vector", nargs="+", type=int, required=True,
                        help="HSQC vector the model was trained on, e.g. --hsqc-vector 0 1 1")
    parser.add_argument("--hsqc-data", required=True, help="HSQC multiplets spreadsheet.")
    parser.add_argument("--gcms-data", required=True, help="GC-MS data spreadsheet.")
    parser.add_argument("--metabolite", required=True, help="Metabolite name to analyse.")
    parser.add_argument("--n-iter", type=int, default=200, help="Monte-Carlo dropout iterations.")
    parser.add_argument("--work-dir", default=".", help="Directory to run in (results are written here).")
    args = parser.parse_args()

    os.chdir(args.work_dir)

    hsqc_vector = args.hsqc_vector
    n_carbons = len(hsqc_vector)

    ia = isotopomerAnalysis.IsotopomerAnalysisNN()
    ia.load_hsqc_and_gcms_data(args.hsqc_data, args.gcms_data)
    ia.inspect_metabolite_data(args.metabolite)

    X_real_data = ia.create_feature_vectors(args.metabolite, hsqc_vector)
    if X_real_data is None:
        raise SystemExit(f"No data found for metabolite '{args.metabolite}'.")

    mean_predictions, std_dev_predictions, predicted_distributions = ia.load_model_and_predict(
        args.model, X_real_data, n_carbons, n_iter=args.n_iter,
    )

    # Back-simulate HSQC/GC-MS from the predicted distribution for comparison.
    predicted_hsqc_data, predicted_gcms_data = ia.simulate_from_predictions(predicted_distributions, hsqc_vector)

    ia.save_results_summary(
        X_real_data, predicted_distributions, std_dev_predictions,
        predicted_hsqc_data, predicted_gcms_data, hsqc_vector,
    )
    print("Done.")


if __name__ == "__main__":
    main()
