"""
Usage: python3 sigma_calc_fromParquet.py --apply-scaling --filter-mass
"""

import argparse
import glob
import os
import awkward as ak
import numpy as np
import ROOT


def getEffSigma_unbinned(masses, weights=None):
    """Calculates unbinned effective sigma (effSigma, center, [low_edge, high_edge]). The center mentioned here is just (m_high + m_low)/2, obtained just to see how far it is from the calculated mode(peak) value. The 'center' isn't used anywhere in the analysis.
    """
    if weights is None:
        sorted_masses = np.sort(masses)
        N = len(sorted_masses)
        k = int(np.ceil(0.683 * N))

        low_edges = sorted_masses[: N - k + 1]
        high_edges = sorted_masses[k - 1 :]
        widths = high_edges - low_edges

        best_idx = np.argmin(widths)

        best_low = low_edges[best_idx]
        best_high = high_edges[best_idx]
    else:
        sort_idx = np.argsort(masses)
        sorted_masses = masses[sort_idx]
        sorted_weights = weights[sort_idx]

        cum_weights = np.cumsum(sorted_weights)
        total_weight = cum_weights[-1]
        target_weight = 0.683 * total_weight

        min_width = float("inf")
        best_low, best_high = 0.0, 0.0

        j = 0
        N = len(sorted_masses)

        #--------------------altering this slightly to make a check--------------
        
        for i in range(N):
            while (
                j < N
                and (cum_weights[j] - (cum_weights[i - 1] if i > 0 else 0))
                < target_weight
            ):
                j += 1
            if j >= N:
                break

            width = sorted_masses[j] - sorted_masses[i]
            if width < min_width:
                min_width = width
                best_low = sorted_masses[i]
                best_high = sorted_masses[j]
                
        #-------------------------------------------------------------------------
        '''
        for i in range(N):
            w_i_prev = cum_weights[i - 1] if i > 0 else 0.0

            while j < N and (cum_weights[j] - w_i_prev) < target_weight:
                j += 1
            if j >= N:
                break

            """
            Trims the overshoot from event j so the window hits exactly 68.3% target weight----using linear interpolation, for a check
            This would have been more useful if I were using a binned algo to get sigma eff.
            """
            w_prev_j = cum_weights[j - 1] - w_i_prev if j > 0 else 0.0
            fraction = (target_weight - w_prev_j) / sorted_weights[j]
            m_high_interp = sorted_masses[j - 1] + fraction * (
                sorted_masses[j] - sorted_masses[j - 1]
            )

            width = m_high_interp - sorted_masses[i]
            if width < min_width:
                min_width = width
                best_low = sorted_masses[i]
                best_high = m_high_interp
                '''
        #---------------------------------------------------------------------------------------

    # captured percentage calculation inside getEffSigma_unbinned
    if weights is None:
        captured_pct = (k / N) * 100.0
    else:
        captured_pct = (target_weight / total_weight) * 100.0

    effSigma = 0.5 * (best_high - best_low)
    center = 0.5 * (best_high + best_low)

    return effSigma, center, (best_low, best_high), captured_pct


def get_root_mean_and_mode(masses, weights=None, nbins=160):
    """Fills a ROOT TH1F histogram to calculate ROOT's internal mean and mode (peak bin center).
    """
    xmin, xmax = float(np.min(masses)), float(np.max(masses))

    ROOT.gROOT.SetBatch(True)

    h = ROOT.TH1F("h_mass_tmp", "Mass Distribution", nbins, xmin, xmax)
    h.Sumw2()

    if weights is not None:
        for x, w in zip(masses, weights):
            h.Fill(float(x), float(w))
    else:
        for x in masses:
            h.Fill(float(x))

    root_mean = h.GetMean()
    max_bin = h.GetMaximumBin()
    root_mode = h.GetXaxis().GetBinCenter(max_bin)

    return root_mean, root_mode


def main():
    parser = argparse.ArgumentParser(
        description="Run effective sigma and mode/mean scan over signal mass points."
    )

    parser.add_argument(
        "-base",
        "--base-dir",
        type=str,
        default="/eos/user/a/arnaik/Sergi_oldPresel_DataBkgSig_toTrain/signal_processedOldPresel",
        help="Base path containing Signal_{i}_GeV directories",
    )
    parser.add_argument(
        "-m",
        "--mass-col",
        type=str,
        default="mass_gggg",
        help="Name of the mass column in Parquet (default: mass_gggg)",
    )
    parser.add_argument(
        "-w",
        "--weight-col",
        type=str,
        default="weight",
        help="Name of the weight column (optional)",
    )
    parser.add_argument(
        "--nbins",
        type=int,
        default=200,
        help="Number of bins for ROOT histogram (default: 200)",
    )
    parser.add_argument(
        "--filter-mass",
        action="store_true",
        help="Filter events to mass_gggg between 115 and 135 GeV",
    )
    parser.add_argument(
        "--apply-scaling",
        action="store_true",
        help="Apply lumi * xs / sumw event weight scaling",
    )
    parser.add_argument(
        "--lumi",
        type=float,
        default=1.0,
        help="Luminosity in fb^-1 (default: 1.0)",
    )

    args = parser.parse_args()

    # Signal samples mapping to gensumw, these are values directly taken from the current analysis script: https://gitlab.cern.ch/takumar/h4g-tools/-/blob/decoupled_from_cmssw/h4g_tools/utils/scales.py?ref_type=heads#L62
    mass_points = list(range(15, 65, 5))  # [15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
    sumw_values = [
        1025691,
        971967,
        969906,
        981050,
        1047791,
        981828,
        974661,
        1022301,
        994105,
        955078,
    ]
    sumw_dict = dict(zip(mass_points, sumw_values))

    print("=" * 75)
    print(
        f"{'Mass Point':<12} | {'Events':<7} | {'Center (GeV)':<12} | {'Mode (GeV)':<12} | {'effSigma (GeV)':<14} | {'Captured (%)':<12}"
    )
    print("=" * 75)

    for m_val in mass_points:
        folder_path = os.path.join(
            args.base_dir, f"Signal_{m_val}_GeV", "nominal"
        )
        parquet_files = glob.glob(os.path.join(folder_path, "*.parquet"))

        if not parquet_files:
            print(
                f"Signal_{m_val:<4}_GeV | WARNING: No parquet files found in {folder_path}"
            )
            continue

       # require cols
        cols_to_load = [args.mass_col]
        if args.weight_col:
            cols_to_load.append(args.weight_col)

        # ------
        events_list = [
            ak.from_parquet(f, columns=cols_to_load) for f in parquet_files
        ]
        events = ak.concatenate(events_list)

        masses = ak.to_numpy(events[args.mass_col])

        weights = None
        if args.weight_col and args.weight_col in events.fields:
            weights = ak.to_numpy(events[args.weight_col])

        # 1. Mass Range Filtering (taking only the SR for all these calculations)
        if args.filter_mass:
            mask = (masses >= 115.0) & (masses <= 135.0)
            masses = masses[mask]
            if weights is not None:
                weights = weights[mask]

        if len(masses) == 0:
            print(
                f"Signal_{m_val:<4}_GeV | 0       | No events remaining after filtering"
            )
            continue

        # 2. Lumi / gensumw Weight Scaling
        if args.apply_scaling:
            sumw = sumw_dict[m_val]
            scale_factor = args.lumi * 1.0 * 1.0 / sumw
            if weights is None:
                weights = np.full_like(masses, scale_factor, dtype=float)
            else:
                weights = weights * scale_factor

        eff_sigma, center, (low_edge, high_edge), captured_pct = (
        getEffSigma_unbinned(masses, weights=weights)
        )

        # Mean and Mode
        root_mean, root_mode = get_root_mean_and_mode(
            masses, weights=weights, nbins=args.nbins
        )

        print(
            f"{m_val:<5} GeV    | {len(masses):<7} | {center:<12.4f} | {root_mode:<12.4f} | {eff_sigma:<14.4f} | {captured_pct:<12.4f}%"
        )

    print("=" * 75)


if __name__ == "__main__":
    main()
