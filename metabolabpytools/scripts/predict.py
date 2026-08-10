#!/usr/bin/env python
"""Predict an isotopomer distribution from real HSQC + GC-MS data.

End-to-end: load a trained model -> build feature vectors from the measured
data -> Monte-Carlo-dropout prediction -> back-simulate HSQC/GC-MS from the
prediction -> save a results summary.

By default it loads the current best model for the given HSQC vector; pass
--model to point at a specific file.

Example:
    python -m metabolabpytools.scripts.predict --hsqc-vector 0 1 1 \\
        --hsqc-data hsqcData1.xlsx --gcms-data gcmsData1.xlsx --metabolite L-LacticAcid
"""
import argparse
import os

from metabolabpytools import isotopomerAnalysis


def main():
    parser = argparse.ArgumentParser(description="Predict an isotopomer distribution from measured data.")
    parser.add_argument("--model", default=None,
                        help="Path to a trained .keras model. Default: the current best saved model "
                             "for --hsqc-vector (in --model-dir).")
    parser.add_argument("--model-dir", default="saved_models",
                        help="Where to look for the best model when --model is not given.")
    parser.add_argument("--hsqc-vector", nargs="+", type=int, required=True,
                        help="HSQC vector the model was trained on, e.g. --hsqc-vector 0 1 1")
    parser.add_argument("--hsqc-data", required=True, help="HSQC multiplets spreadsheet.")
    parser.add_argument("--gcms-data", required=True, help="GC-MS data spreadsheet.")
    parser.add_argument("--metabolite", required=True, help="Metabolite name to analyse.")
    parser.add_argument("--n-iter", type=int, default=200, help="Monte-Carlo dropout iterations.")
    parser.add_argument("--work-dir", default=".", help="Directory to run in (results are written here).")
    parser.add_argument("--labeled-only", action="store_true",
                        help="The model was trained with mlp-train --labeled-only. Its conditional "
                             "labelled pattern is recombined with the measured GC-MS unlabelled "
                             "fraction to report the absolute distribution.")
    args = parser.parse_args()

    os.chdir(args.work_dir)

    hsqc_vector = args.hsqc_vector
    n_carbons = len(hsqc_vector)

    ia = isotopomerAnalysis.IsotopomerAnalysisNN()

    # Resolve the model: explicit --model, else the current best for this vector.
    model_path = args.model or ia.best_model_path(hsqc_vector, labeled=args.labeled_only, directory=args.model_dir)
    if not os.path.exists(model_path):
        raise SystemExit(
            f"Model not found: {model_path}\n"
            f"Train one first, e.g.  mlp-train --hsqc-vector {' '.join(map(str, hsqc_vector))}"
            + (" --labeled-only" if args.labeled_only else ""))
    print(f"Using model: {model_path}")

    ia.load_hsqc_and_gcms_data(args.hsqc_data, args.gcms_data)
    ia.inspect_metabolite_data(args.metabolite)

    X_real_data = ia.create_feature_vectors(args.metabolite, hsqc_vector)
    if X_real_data is None:
        raise SystemExit(f"No data found for metabolite '{args.metabolite}'.")

    if args.labeled_only:
        # The measured unlabelled fraction is the GC-MS M+0 percentage per experiment.
        measured_unlabeled = [ia.exp_gcms[args.metabolite][e][0] for e in range(ia.n_exps)]
        mean_predictions, std_dev_predictions, predicted_distributions = ia.load_model_and_predict_labeled(
            model_path, X_real_data, n_carbons, measured_unlabeled, n_iter=args.n_iter,
        )
    else:
        mean_predictions, std_dev_predictions, predicted_distributions = ia.load_model_and_predict(
            model_path, X_real_data, n_carbons, n_iter=args.n_iter,
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
