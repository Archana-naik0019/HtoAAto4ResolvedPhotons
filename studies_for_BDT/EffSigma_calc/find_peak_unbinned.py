import os
import glob
import pandas as pd
import numpy as np
import ROOT

def fit_continuous_mode(mass_array, weight_array, initial_center, sigma_eff, tol=1e-4, max_iters=5):
    """
    Gaussian fit in RooFit over the mass_array events to extract the peak
    """
    # mass_gggg within Signal Region bounds
    mass_gggg = ROOT.RooRealVar("mass_gggg", "m_{#gamma#gamma#gamma#gamma}", 115.0, 135.0, "GeV")
    weight_var = ROOT.RooRealVar("weight", "event_weight", -100.0, 1000.0)

    current_center = initial_center

    for iteration in range(max_iters):
        # initial window : Center +/- 2 * sigma_eff
        fit_min = current_center - 2.0 * sigma_eff
        fit_max = current_center + 2.0 * sigma_eff

        # ROOT.RooFit.WeightVar(weight_var) to specify the weight properly
        dataset = ROOT.RooDataSet(
            "dataset", 
            "Signal SubSet", 
            ROOT.RooArgSet(mass_gggg, weight_var), 
            ROOT.RooFit.WeightVar(weight_var)
        )

        for m_val, w_val in zip(mass_array, weight_array):
            if fit_min <= m_val <= fit_max:
                mass_gggg.setVal(float(m_val))
                weight_var.setVal(float(w_val))
                dataset.add(ROOT.RooArgSet(mass_gggg, weight_var), float(w_val))

        # RooFit Gaussian PDF in the window defined
        mean  = ROOT.RooRealVar("mean", "Fitted Mode", current_center, fit_min, fit_max)
        sigma = ROOT.RooRealVar("sigma", "Resolution Width", sigma_eff, 0.2 * sigma_eff, 2.0 * sigma_eff)
        gauss_pdf = ROOT.RooGaussian("gauss_pdf", "Local Gaussian Fit", mass_gggg, mean, sigma)

        # unbinned maximum likelihood fit (with SumW2 error correction for event weights)
        fit_result = gauss_pdf.fitTo(
            dataset,
            ROOT.RooFit.Range(fit_min, fit_max),
            ROOT.RooFit.SumW2Error(True),
            ROOT.RooFit.Save(True),
            ROOT.RooFit.PrintLevel(-1)
        )

        fitted_mode = mean.getVal()
        fitted_mode_err = mean.getError()

        print(f"  [Iter {iteration + 1}] Window [{fit_min:.3f}, {fit_max:.3f}] GeV -> Fitted Peak = {fitted_mode:.4f} +/- {fitted_mode_err:.4f} GeV")

        # stoping the iterative fit if peak shift is smaller than 0.1 MeV
        if abs(fitted_mode - current_center) < tol:
            print(f"  ------- Converged on Mode: {fitted_mode:.4f} GeV\n")
            return fitted_mode, fitted_mode_err

        current_center = fitted_mode

    return current_center, mean.getError()


# --------------------------------------------------------------------------------
# these values are taken from teh output of sigma_calc_fromParquet.py
if __name__ == "__main__":
    table_data = {
        15: (124.4827, 2.4857),
        20: (124.5152, 2.5042),
        25: (124.4598, 2.5029),
        30: (124.3958, 2.5445),
        35: (124.4146, 2.5586),
        40: (124.5448, 2.5918),
        45: (124.5406, 2.6262),
        50: (124.4768, 2.6172),
        55: (124.4519, 2.6075),
        60: (124.3018, 2.6073),
    }

    base_dir = "/eos/home-a/arnaik/Sergi_oldPresel_DataBkgSig_toTrain/signal_processedOldPresel"

    fitted_results = {}

    for m_point, (center_seed, eff_sig) in table_data.items():
        print(f"=== Reading and Fitting Mass Point: {m_point} GeV ===")

        dir_path = os.path.join(base_dir, f"Signal_{m_point}_GeV", "nominal")
        parquet_files = glob.glob(os.path.join(dir_path, "*.parquet"))

        if not parquet_files:
            print(f"WARNING: No parquet files found in {dir_path}")
            continue

        df_list = [pd.read_parquet(f) for f in parquet_files]
        combined_df = pd.concat(df_list, ignore_index=True)

        mass_array = combined_df['mass_gggg'].to_numpy()
        weight_array = combined_df['weight'].to_numpy() if 'weight' in combined_df.columns else np.ones(len(mass_array))

        mode_val, mode_err = fit_continuous_mode(
            mass_array=mass_array,
            weight_array=weight_array,
            initial_center=center_seed,
            sigma_eff=eff_sig
        )

        fitted_results[m_point] = (mode_val, mode_err)

    print("\n---------------------------------------------------------------------------")
    print("      CONTINUOUS FITTED PEAKS SUMMARY             ")
    for m_point, (mode_val, mode_err) in fitted_results.items():
        print(f"Mass Point {m_point:2d} GeV : Peak = {mode_val:.4f} +/- {mode_err:.4f} GeV")
