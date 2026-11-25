'''
usage: python3 AMS_crossCheck_windows.py -y 2024 -xgb ArcTrialBDT_12Oct2024 --fixed_cuts --nCats 1

'''

from ROOT import TGraph, TFile, gROOT, TTree, kRed, kOrange, kMagenta, kBlack, kBlue, kGreen, TCanvas, TLegend, gStyle, TH1D
import os
import json
import pandas as pd
from typing import List
import array
import numpy as np
from h4g_tools.bdt.categories import defineAMS, sumSignificance
from h4g_tools.bdt.reweight import reweight1Dim
from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist
from h4g_tools.utils.plotting import plot_smooth, plot_smooth_all, plot_reweighting, plotFancy, plot_smoothed_hists
from h4g_tools.utils.scales import Scales as ScalesCls, SideBand, SignalRegion
from h4g_tools.utils.runner_utils import get_parser_cats, handle_errors
gROOT.LoadMacro("/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_for_reweighting/h4g-tools/h4g_tools/bdt/smoother.cpp")
from ROOT import SmoothHist #added to check crash reason
import ROOT
_SmoothHist_cpp = SmoothHist

def SmoothHist(name, hist, smoothness, span):
    """Return (TH1D, TGraph) whether or not cppyy can unpack the std::tuple."""
    res = _SmoothHist_cpp(name, hist, smoothness, span)
    try:
        h, g = res                      
    except TypeError:
        h = ROOT.std.get[0](res)    
        g = ROOT.std.get[1](res)      
    return h, g
#from ROOT import SmoothHist, SplitSignalHistograms, MergeHistograms, MergeGraphs
from operator import itemgetter
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# No popups with matplotlib
matplotlib.use('Agg')

# ---------------------------------------------------------------------------
# Fixed (from teh current version of AN as of 16th Nov 2026) BDT cuts.
# Only used when running with --fixed_cuts. Taken from AN (Table 5), to check the AMS values.
# ---------------------------------------------------------------------------
FIXED_CUTS = {
    "15_GeV": 0.960,
    "20_GeV": 0.960,
    "25_GeV": 0.960,
    "30_GeV": 0.970,
    "35_GeV": 0.965,
    "40_GeV": 0.970,
    "45_GeV": 0.975,
    "50_GeV": 0.980,
    "55_GeV": 0.980,
    "60_GeV": 0.980,
}

# ---------------------------------------------------------------------------
# +- 2 sigma windows (calculated narrowest mass windows that contain 95.4% info ), to be added directly form the 2-point slider check, this would mostly take care of the m_4g dist shape asymmetry as well.
# ---------------------------------------------------------------------------
# full SR
'''
WINDOWS = {
    "15_GeV": [115, 135],
    "20_GeV": [115, 135],
    "25_GeV": [115, 135],
    "30_GeV": [115, 135],
    "35_GeV": [115, 135],
    "40_GeV": [115, 135],
    "45_GeV": [115, 135],
    "50_GeV": [115, 135],
    "55_GeV": [115, 135],
    "60_GeV": [115, 135],
}

'''

# newly calculated mass windows for new presel samples (these were obtained by fittinf using sum of 3 gaussians)

WINDOWS = {
    "15_GeV": [119.5306, 129.4734],
    "20_GeV": [119.4736, 129.4904],
    "25_GeV": [119.5782, 129.5898],
    "30_GeV": [119.4550, 129.6330],
    "35_GeV": [119.4808, 129.7152],
    "40_GeV": [119.3354, 129.7026],
    "45_GeV": [119.2786, 129.7834],
    "50_GeV": [119.4276, 129.8964],
    "55_GeV": [119.3450, 129.7750],
    "60_GeV": [119.2644, 129.6936],
}
'''
# these are from the original script used in analysis, calculate_categories.py. so, assuming its for old-preselections

WINDOWS = {
    "15_GeV": [119.4845, 129.4045],
    "20_GeV": [119.4737, 129.3137],
    "25_GeV": [119.6140, 129.4140],
    "30_GeV": [119.6707, 129.3107],
    "35_GeV": [119.6374, 129.3974],
    "40_GeV": [119.6517, 129.4917],
    "45_GeV": [119.3271, 129.1271],
    "50_GeV": [119.5782, 129.4582],
    "55_GeV": [119.4636, 129.3036],
    "60_GeV": [119.2573, 129.2573],
}
'''

'''
# these are obtained by be using an older approach that assumes ~symmetruc m_4g distribution--------------------------------

#for old preselections
WINDOWS = {
    "15_GeV": [119.4151, 129.3751],
    "20_GeV": [119.4051, 129.4051],
    "25_GeV": [119.4228, 129.4228],
    "30_GeV": [119.3723, 129.5323],
    "35_GeV": [119.3257, 129.5657],
    "40_GeV": [119.2991, 129.6591],
    "45_GeV": [119.2095, 129.7295],
    "50_GeV": [119.2365, 129.7165],
    "55_GeV": [119.1905, 129.6305],
    "60_GeV": [119.1025, 129.5425],
}
#for new preselections
WINDOWS = {
    "15_GeV": [119.4001, 129.3601],
    "20_GeV": [119.3903, 129.3503],
    "25_GeV": [119.4120, 129.4120],
    "30_GeV": [119.3313, 129.5713],
    "35_GeV": [119.3505, 129.5505],
    "40_GeV": [119.3293, 129.6493],
    "45_GeV": [119.1904, 129.6704],
    "50_GeV": [119.2579, 129.6979],
    "55_GeV": [119.2085, 129.6085],
    "60_GeV": [119.1085, 129.5085],
}
#----------------------------------------------------------------------------------------------------------------------------
'''
'''
# these are obtained by be using the 2-point slider for narrowest 95.4% info, approach that accounts for asymmetry in m_4g distribution shape-----------

# for old preselections
WINDOWS = {
    "15_GeV": [117.5704, 130.1368],
    "20_GeV": [117.8387, 130.4313],
    "25_GeV": [117.7149, 130.4936],
    "30_GeV": [117.4867, 130.6240],
    "35_GeV": [117.7412, 130.9233],
    "40_GeV": [117.5065, 131.0060],
    "45_GeV": [117.9179, 131.2926],
    "50_GeV": [117.3418, 130.6060],
    "55_GeV": [117.6578, 130.6263],
    "60_GeV": [117.5187, 130.4288],
}

# for new preselections
WINDOWS = {
    "15_GeV": [117.4746, 129.9795],
    "20_GeV": [118.0030, 130.5530],
    "25_GeV": [117.7324, 130.4429],
    "30_GeV": [117.5077, 130.6693],
    "35_GeV": [117.6522, 130.7018],
    "40_GeV": [117.4884, 131.0266],
    "45_GeV": [117.9190, 131.3227],
    "50_GeV": [117.3418, 130.5263],
    "55_GeV": [117.6943, 130.6305],
    "60_GeV": [117.5924, 130.4379],
}
# this is ~2GeV wider than the window width obtained using the older approach (that assumes symmetry)
#------------------------------------------------------------------------------------------------------------------------------------------------------
'''

def plot_ams_score(
    inputs: list,
    mass: str,
    translator_hist: TH1D,
    best_cut: float,
    year: str,
    xgb_name: str,
    cwd: str
) -> None:
    """
    Plots AMS vs BDT score for each mass point.
    """

    max_ams = 0.0
    max_ams_noRes = 0.0
    ams = []
    cut = []
    invalid = []
    invalid_cut = []
    s_sr = []
    b_sr = []
    d_sb = []
    for tup in inputs:
        if tup[0] > max_ams_noRes:
            max_ams_noRes = tup[0]
        if tup[1] >= 10:
        #if tup[1] >= 0: # this line is temporarily changed for a check (soon to be reverted back to 10)
            ams.append(tup[0])
            cut.append(translator_hist.GetXaxis().GetBinLowEdge(tup[2]))
            if tup[0] > max_ams:
                data_sb = tup[1]
                max_ams = tup[0]
        else:
            invalid.append(tup[0])
            invalid_cut.append(translator_hist.GetXaxis().GetBinLowEdge(tup[2]))
        s_sr.append(tup[3])
        b_sr.append(tup[4])
        d_sb.append(tup[1])

    fig, ax1 = plt.subplots(figsize = (8,8))
    title = f"AMS vs BDT Cut for m_{{a}} = {mass.replace('_',' ')}"
    ax1.scatter(cut, ams, s=5, c="b", label="AMS")
    ax1.scatter(invalid_cut, invalid, s=5, c="r", label="AMS")
    ax1.axvline(translator_hist.GetXaxis().GetBinLowEdge(best_cut))
    ax1.set_xlabel("Cut on BDT Score", fontsize=14)
    ax1.set_ylabel("AMS", fontsize=14)
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
    if year == "2022":
        ax1.set_ylim(top=23.0)
    elif year == "2024":
        ax1.set_ylim(top=39.0)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig(f"{cwd}/plots/BDT/{year}/cats/{xgb_name}/ams_vs_bdt_{mass}.pdf")

    # Save to json for future plotting/comparison
    to_json = {mass: {"ams": ams, "cut": cut, "invalid_cut": invalid_cut, "invalid": invalid, "max_ams": max_ams, "max_ams_noRestrictions": max_ams_noRes, "axline": translator_hist.GetXaxis().GetBinLowEdge(best_cut), "data_sb": d_sb, "signal_sr": s_sr, "bkg_sr": b_sr}}
    rw = "w"
    if os.path.exists(f"{cwd}/cats/{year}/{xgb_name}/ams_vs_bdt.json"):
        rw = "r+"
    if rw == "w":
        with open(f"{cwd}/cats/{year}/{xgb_name}/ams_vs_bdt.json", rw) as f:
            json.dump(to_json, f, indent=4)
    else:
        with open(f"{cwd}/cats/{year}/{xgb_name}/ams_vs_bdt.json", rw) as f:
            in_json = json.load(f)
            rw = "w"
        with open(f"{cwd}/cats/{year}/{xgb_name}/ams_vs_bdt.json", rw) as f:
            in_json.update(to_json)
            to_json = in_json
            json.dump(to_json, f, indent=4)


if __name__ == "__main__":

    parser = get_parser_cats()
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    parser.add_argument("-xgb", "--xgb_name", type=str, required=True, help="XGBoost model name.")
    parser.add_argument("--fixed_cuts", action="store_true",
                        help="Use FIXED_CUTS (single category [cut, 1.0]) instead of the optimized partition.")
    args = parser.parse_args()

    pathSetup(year=args.year, model=args.xgb_name)

    handle_errors(
        ["Span must be between 0.0 and 1.0", (args.span < 0.0 or args.span > 1.0), (args.span != -1.0)]
    )

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    Scales = ScalesCls(args.year if "2024" not in args.year else "2024")

    print("\nLoading Signal samples:")
    print(f"{cwd}/outputs/BDT/{args.year}/{args.xgb_name}/signal/")
    signal = loadPerMass(f"{cwd}/outputs/BDT/{args.year}/{args.xgb_name}/signal/", ["BDT_score", "mass_gggg", "weight"], bdt=True)
    print("\nLoading EM Background samples:")
    print(f"{cwd}/outputs/BDT/{args.year}/{args.xgb_name}/background/")
    background = loadPerMass(f"{cwd}/outputs/BDT/{args.year}/{args.xgb_name}/background/", ["BDT_score", "mass_gggg"], bdt=True)
    print("\nLoading Data samples:")
    print(f"{cwd}/outputs/BDT/{args.year}/{args.xgb_name}/data/")
    data = loadPerMass(f"{cwd}/outputs/BDT/{args.year}/{args.xgb_name}/data/", ["BDT_score", "mass_gggg"], bdt=True, nominalMassOnly=True)

    print("\nChecking if inputs are transformed:")
    transformed = False
    check_transformed = loadSamples(f"{cwd}/outputs/BDT/{args.year}/{args.xgb_name}/data/")
    if "BDT_score_orig" in check_transformed:
        transformed = True
    print(f"Transformed: {transformed}")

    # Cut off l_edge for smoothing
    l_edge = 0.05 if not transformed else 0.0
    r_edge = 1.0
    bins = 200.0
    retries = 10

    # Signal-region windows for calculating S and B: taken directly from WINDOWS
    windows = {key: list(val) for key, val in WINDOWS.items()}
    for _m in range(15, 65, 5):
        _k = f"{_m}_GeV"
        assert _k in windows, f"No window defined for mass point {_k}"
        assert len(windows[_k]) == 2 and windows[_k][0] < windows[_k][1], f"Bad window for {_k}: {windows[_k]}"
    print("\nUsing user-defined signal-region windows:")
    for _k, _w in windows.items():
        print(f"  {_k}: [{_w[0]:.4f}, {_w[1]:.4f}]")

    sig_samples_df = {}
    bkg_samples_df = {}
    data_samples_df = {}

    #-----------------ADDING THIS SECTION TO GET SOME INTERMEDIATE INFO-*********************
    # Create directory to store intermediate ROOT objects
    debug_dir = f"{cwd}/plots/BDT/{args.year}/cats/{args.xgb_name}/intermediates"
    os.makedirs(debug_dir, exist_ok=True)

    # Create a master ROOT file to store all steps, histograms, and smooth curves
    debug_root_file = TFile(f"{debug_dir}/all_intermediate_steps.root", "RECREATE")
    #---------------------********************************************************************

    smoothed_hists_mass = {}
    check_store = {}   # sanity-check info per mass point
    print(f"\nNumber of categories set to {args.nCats}.")
    print(f"Number of bins set to {bins} (minus those left of {l_edge}).")
    if args.fixed_cuts:
        print("\n*** --fixed_cuts: using user-defined BDT cuts (single category per mass point) ***")

    # ******%%%%---- 80/20 split: test = global row index % 5 == 0 (must match split_dataset in training) ----

    TEST_MOD = 2 # set this to 5 for 20/80 split
    sig_offset, _o = {}, 0
    for _m in range(15, 65, 5):
        _k = f"{_m}_GeV"
        sig_offset[_k] = _o   # signal was stacked 15, 20, ..., 60 in training, so the index is global
        _o += len(signal[_k])

    # ******%%%%---------------------------------------------------------------------------------------------

    for idx, key in enumerate([f"{m}_GeV" for m in range(15,65,5)]):
        print(f"\n\n ~~~~~~~  Starting smoothing procedure for mass point: {key.replace('_',' ')}  ~~~~~~~")

        #************%%%%%%%%%%%%%%%%%%%%%%%%%%%

        print(f"Keeping test events only (global index % {TEST_MOD} == 0) for signal and background.")
        sig_sample_df = signal[key][((signal[key].index + sig_offset[key]) % TEST_MOD == 0) & (signal[key]["weight"] < 10.0)]
        bkg_sample_df = background[key][background[key].index % TEST_MOD == 0]
        print(f"[check] test events: signal {len(sig_sample_df)}, background {len(bkg_sample_df)}")

        #*********%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sig_sample_df.reset_index(inplace=True)
        bkg_sample_df.reset_index(inplace=True)
        data_sample_df = data[key]

        print("Reweighting background to data...")
        bkg_sample_df = reweight1Dim(bkg_sample_df, data_sample_df, "BDT_score")

        bkg_sample_df[["BDT_score", "1dimreweight"]].to_csv(f"reweighted_weights_{key}.csv", index=False) # to plot these values distribution

        # Add bkg/data samples to dictionaries
        sig_samples_df.update({key: sig_sample_df})
        bkg_samples_df.update({key: bkg_sample_df})
        data_samples_df.update({key: data_sample_df})

        ########################################################################
        ## Reweight BDT scores to make future cuts accurate (apply category cuts to reweighted BDT scores)
        # Fill signal hists
        print("Filling signal hist...")
        # Signal in signal region
        # WITH WINDOWS
        sig_hist_sig = fillHist(f"BDT_score", SignalRegion(sig_sample_df, window=windows[key]), [l_edge, r_edge], name=f"BDT_score_signal_SR_{key}", binsScale=bins, normalize=False, preventOverFlow=True) # ***this was originally commented out
        # NO WINDOWS:
        #sig_hist_sig = fillHist(f"BDT_score", SignalRegion(sig_sample_df), [l_edge, r_edge], name=f"BDT_score_signal_SR_{key}", binsScale=bins, normalize=False, preventOverFlow=True) #***this was originally used in the code
        # Signal in sideband region
        sig_sideband = fillHist(f"BDT_score", SideBand(sig_sample_df), [l_edge, r_edge], name=f"BDT_score_signal_SB_{key}", binsScale=bins, normalize=False, preventOverFlow=True)

        # Fill EM background hist
        print("Filling EM background hist...")
        # EM background in signal region
        # WITH WINDOWS:
        bkg_hist_sig = fillHist(f"BDT_score", SignalRegion(bkg_sample_df, window=windows[key]), [l_edge, 1.0], name=f"BDT_score_bkg_SR_{key}", binsScale=bins, normalize=False, customWeights="1dimreweight", preventOverFlow=True) #***this was originally commented out
        # Full-SR background (no window), used ONLY to get the normalization to data SR
        bkg_hist_full = fillHist(f"BDT_score", SignalRegion(bkg_sample_df), [l_edge, 1.0], name=f"BDT_score_bkg_SRfull_{key}", binsScale=bins, normalize=False, customWeights="1dimreweight", preventOverFlow=True)
        # NO WINDOWS:
        #bkg_hist_sig = fillHist(f"BDT_score", SignalRegion(bkg_sample_df), [l_edge, 1.0], name=f"BDT_score_bkg_SR_{key}", binsScale=bins, normalize=False, customWeights="1dimreweight", preventOverFlow=True) #***this was originally used
        # EM background in sideband region
        bkg_sideband = fillHist(f"BDT_score", SideBand(bkg_sample_df), [l_edge, 1.0], name=f"BDT_score_bkg_SB_{key}", binsScale=bins, normalize=False, customWeights="1dimreweight", preventOverFlow=True)

        # Fill data hist
        print("Filling data hist...")
        print(len(data_sample_df), len(SignalRegion(data_sample_df)), len(SideBand(data_sample_df)))
        # Data in signal region for normalization
        data_hist_sig = fillHist(f"BDT_score", SignalRegion(data_sample_df), [l_edge, 1.0], name=f"BDT_score_data_SR_{key}", binsScale=bins, normalize=False, preventOverFlow=True)
        # Data in sideband region
        data_sideband = fillHist(f"BDT_score", SideBand(data_sample_df), [l_edge, 1.0], name=f"BDT_score_data_SB_{key}", binsScale=bins, normalize=False, preventOverFlow=True)

        #********************************-----------------------------------------------------------------
        # 1: Raw & Reweighted Binned Histograms (Unnormalized)
        # Create a sub-folder inside ROOT file for the current mass point
        mass_dir = debug_root_file.mkdir(f"mass_{key}")
        mass_dir.cd()
        # Raw / Reweighted filled ROOT histograms preserving exact binning
        sig_hist_sig.Clone(f"hist_sig_SR_raw").Write()
        sig_sideband.Clone(f"hist_sig_SB_raw").Write()
        bkg_hist_sig.Clone(f"hist_bkg_SR_reweighted").Write()   # Bkg SR after 1D reweighting
        bkg_sideband.Clone(f"hist_bkg_SB_reweighted").Write()   # Bkg SB after 1D reweighting
        data_hist_sig.Clone(f"hist_data_SR_raw").Write()
        data_sideband.Clone(f"hist_data_SB_raw").Write()

        # Raw bkg hists (BEFORE 1D reweight): same selection/binning, no customWeights
        bkg_hist_sig_raw = fillHist(f"BDT_score", SignalRegion(bkg_sample_df, window=windows[key]), [l_edge, 1.0], name=f"BDT_score_bkg_SR_raw_{key}", binsScale=bins, normalize=False, preventOverFlow=True)
        bkg_sideband_raw = fillHist(f"BDT_score", SideBand(bkg_sample_df), [l_edge, 1.0], name=f"BDT_score_bkg_SB_raw_{key}", binsScale=bins, normalize=False, preventOverFlow=True)
        bkg_hist_sig_raw.Clone("hist_bkg_SR_raw").Write()
        bkg_sideband_raw.Clone("hist_bkg_SB_raw").Write()
        # ==================================================================-------------------------

        # Normalize background to data
        #bkg_hist_sig.Scale(data_hist_sig.Integral() / bkg_hist_sig.Integral())
        bkg_hist_sig.Scale(data_hist_sig.Integral() / bkg_hist_full.Integral())   # scale from full SR, applied to the (windowed) hist
        bkg_sideband.Scale(data_sideband.Integral() / bkg_sideband.Integral())

        # lumi * xs(1fb) * BR(100%) / sum(# of generated events)
        sig_hist_sig.Scale(Scales.lumi_fb * 1.0 * 1.0 / Scales.sumw[args.year if "2024" not in args.year else "2024"][key])
        sig_sideband.Scale(Scales.lumi_fb * 1.0 * 1.0 / Scales.sumw[args.year if "2024" not in args.year else "2024"][key])

        check_store[key] = {
            "bkg": bkg_hist_sig.Clone(f"chk_bkg_{key}"),
            "sig": sig_hist_sig.Clone(f"chk_sig_{key}"),
            "data_sr_total": data_hist_sig.Integral(),
        }
        check_store[key]["bkg"].SetDirectory(0)
        check_store[key]["sig"].SetDirectory(0)

        # ************************------------------------------------------------------------------
        # 2: Normalized / Rescaled Histograms
        mass_dir.cd()
        sig_hist_sig.Clone(f"hist_sig_SR_normalized").Write()
        bkg_hist_sig.Clone(f"hist_bkg_SR_normalized").Write()
        #-------------------***********************************************


        ########################################################################
        ### Smooth signal region histograms
        #bkg_data_smoothness = 8
        bkg_data_smoothness = 7.5
        bkg_data_span = -1.0
        print("\nCalculating signficance with BDT distributions:")
        print(f"Signal smoothness set to {args.smoothness}")
        print(f"Signal span set to {args.span}")
        print(f"Background and data smoothness set to {bkg_data_smoothness}")
        print(f"Background and data span set to {bkg_data_span}")

        ## Smooth signal: signal region
        print("\nSmoothing signal hist in SIGNAL REGION...")
        sig_graph_hist, sig_graph = SmoothHist(f"sig_smooth_sig_{key}", sig_hist_sig, args.smoothness, args.span)
        count_retries = 0
        while sig_graph.GetN() <= 1 and count_retries < retries:
            print("Retrying to smooth signal.")
            sig_graph_hist, sig_graph = SmoothHist(f"sig_smooth_sig_{key}", sig_hist_sig, args.smoothness, args.span)
            count_retries += 1
        assert sig_graph.GetN() > 1, f"Failed to smooth background. N: {sig_graph.GetN()}"

        ## Smooth signal: sideband region
        print("Smoothing signal hist in SIDEBAND REGION...")
        sig_sb_graph_hist, sig_sb_graph = SmoothHist(f"sig_smooth_sb_{key}", sig_sideband, args.smoothness, args.span)
        count_retries = 0
        while sig_sb_graph.GetN() <= 1 and count_retries < retries:
            print("Retrying to smooth background.")
            sig_sb_graph_hist, sig_sb_graph = SmoothHist(f"sig_smooth_sb_{key}", sig_sideband, args.smoothness, args.span)
            count_retries += 1
        assert sig_sb_graph.GetN() > 1, f"Failed to smooth background. N: {sig_sb_graph.GetN()}"

        ## Smooth EM background: signal region
        print("Smoothing EM background hist in SIGNAL REGION...")
        bkg_graph_hist, bkg_graph = SmoothHist(f"bkg_smooth_{key}", bkg_hist_sig, bkg_data_smoothness, bkg_data_span)
        count_retries = 0
        while bkg_graph.GetN() <= 1 and count_retries < retries:
            print("Retrying to smooth background.")
            bkg_graph_hist, bkg_graph = SmoothHist(f"bkg_smooth_{key}", bkg_hist_sig, bkg_data_smoothness, bkg_data_span)
            count_retries += 1
        assert bkg_graph.GetN() > 1, f"Failed to smooth background. N: {bkg_graph.GetN()}"

        ## Smooth EM background: sideband region
        print("Smoothing EM background hist in SIDEBAND REGION...")
        bkg_sb_graph_hist, bkg_sb_graph = SmoothHist(f"bkg_smooth_sb_{key}", bkg_sideband, bkg_data_smoothness, bkg_data_span)
        count_retries = 0
        while bkg_sb_graph.GetN() <= 1 and count_retries < retries:
            print("Retrying to smooth background.")
            bkg_sb_graph_hist, bkg_sb_graph = SmoothHist(f"bkg_smooth_sb_{key}", bkg_sideband, bkg_data_smoothness, bkg_data_span)
            count_retries += 1
        assert bkg_sb_graph.GetN() > 1, f"Failed to smooth background. N: {bkg_sb_graph.GetN()}"

        ## Smooth data: sideband region
        print("Smoothing data hist in SIDEBAND REGION...")
        # Sideband region
        data_sb_graph_hist, data_sb_graph = SmoothHist(f"data_smooth_{key}", data_sideband, bkg_data_smoothness, bkg_data_span)
        count_retries = 0
        while data_sb_graph.GetN() <= 1 and count_retries < retries:
            print("Retrying to smooth data.")
            data_sb_graph_hist, data_sb_graph = SmoothHist(f"data_smooth_{key}", data_sideband, bkg_data_smoothness, bkg_data_span)
            count_retries += 1
        assert data_sb_graph.GetN() == bkg_sb_graph.GetN(), f"Failed to smooth data. Data N: {data_sb_graph.GetN()}\t Background N: {bkg_sb_graph.GetN()}"


        # *************************************---------------------------------------------
        #3: All Smoothed Histograms & Continuous TGraphs
        mass_dir.cd()

        # 1. Smoothed Histograms (TH1D)
        sig_graph_hist.Clone("hist_sig_SR_smoothed").Write()
        sig_sb_graph_hist.Clone("hist_sig_SB_smoothed").Write()
        bkg_graph_hist.Clone("hist_bkg_SR_smoothed").Write()
        bkg_sb_graph_hist.Clone("hist_bkg_SB_smoothed").Write()
        data_sb_graph_hist.Clone("hist_data_SB_smoothed").Write()

        # 2. Smoothed Continuous Curves (TGraph)
        sig_graph.Clone("graph_sig_SR_smooth_curve").Write()
        sig_sb_graph.Clone("graph_sig_SB_smooth_curve").Write()
        bkg_graph.Clone("graph_bkg_SR_smooth_curve").Write()
        bkg_sb_graph.Clone("graph_bkg_SB_smooth_curve").Write()
        data_sb_graph.Clone("graph_data_SB_smooth_curve").Write()

        debug_root_file.cd()
        #-----------------------------------------*****************************

        # Plot BDT distributions used in AMS calculation
        # Signal in SR, other in SB - input hist + smoothed graph
        plot_smooth_all(
            f"{cwd}/plots/BDT/{args.year}/cats/{args.xgb_name}/sideband_sigSR_plots_inputs_{key}.pdf",
            key.replace("_", " "),
            [sig_hist_sig, sig_graph, kRed, kOrange+3, "Signal SR", "l", "hist same"],
            [bkg_hist_sig, bkg_graph, kBlue, kGreen+3, "EM SR", "l", "hist same"],
            [data_sideband, None, kBlack, kMagenta, "Data SB", "ple", "hist P same"],
            Scales = Scales,
        )

        # ALL in SR - input hist + smoothed graph
        # KEEP THIS BLINDED UNTIL READY!
        """
        plot_smooth_all(
            f"{cwd}/plots/BDT/{args.year}/cats/{args.xgb_name}/allSR_plots_inputs_{key}.pdf",
            key.replace("_", " "),
            [sig_hist_sig, sig_graph, kRed, kOrange+3, "Signal SR", "l", "hist same"],
            [bkg_hist_sig, bkg_graph, kBlue, kGreen+3, "Bkg SR", "l", "hist same"],
            [data_hist_sig, None, kBlack, kMagenta, "Data SR", "ple", "hist P same"],
            Scales = Scales,
        )
        """

        # All in SB - input hist + smoothed graph
        plot_smooth_all(
            f"{cwd}/plots/BDT/{args.year}/cats/{args.xgb_name}/sideband_plots_inputs_{key}.pdf",
            key.replace("_", " "),
            [sig_sideband, sig_sb_graph, kRed, kOrange+3, "Signal SB", "l", "hist same"],
            [bkg_sideband, bkg_sb_graph, kBlue, kGreen+3, "EM SB", "l", "hist same"],
            [data_sideband, None, kBlack, kMagenta, "Data SB", "ple", "hist P same"],
            Scales = Scales,
        )

        # All in SB - smoothed graph + smoothed hist
        plot_smooth_all(
            f"{cwd}/plots/BDT/{args.year}/cats/{args.xgb_name}/sideband_plots_inputs_wErrors_{key}.pdf",
            key.replace("_", " "),
            [sig_sideband, sig_sb_graph, kRed, kOrange+3, "Signal SB", "l", "hist same"],
            [bkg_sideband, bkg_sb_graph, kBlue, kGreen+3, "EM SB", "l", "hist same"],
            [data_sideband, data_sb_graph, kBlack, kMagenta, "Data SB", "ple", "P E1 same"],
            Scales = Scales,
        )

        plot_smooth_all(
            f"{cwd}/plots/BDT/{args.year}/cats/{args.xgb_name}/sideband_plots_outputs_{key}.pdf",
            key.replace("_", " "),
            [sig_graph_hist, sig_graph, kRed, kOrange+3, "Signal SB", "l", "hist same"],
            [bkg_graph_hist, bkg_graph, kBlue, kGreen+3, "EM SB", "l", "hist same"],
            [data_sb_graph_hist, data_sb_graph, kBlack, kMagenta, "Data SB", "ple", "hist P same"],
            Scales = Scales,
        )

        sig_graph.Delete()
        sig_sb_graph.Delete()
        bkg_graph.Delete()
        bkg_sb_graph.Delete()
        data_sb_graph.Delete()

        ###########################################################################
        ## BDT Spectrum Method ##
        # Get histograms from smoothed TGraphs
        smoothed_hists = {}
        smoothed_hists["signal_SR"] = sig_graph_hist
        smoothed_hists["signal_SB"] = sig_sb_graph_hist
        smoothed_hists["bkg_SR"] = bkg_graph_hist
        smoothed_hists["bkg_SB"] = bkg_sb_graph_hist
        smoothed_hists["data_SB"] = data_sb_graph_hist

        plot_smoothed_hists(smoothed_hists, f"{cwd}/plots/BDT/{args.year}/cats/{args.xgb_name}/smoothed_hists_{key}.pdf", key, idx, Scales)
        smoothed_hists_mass.update({key: smoothed_hists})


    # Loop through number of categories all at once to save compute time on smoothing

    # ================= SANITY CHECK of Signal SR / Bkg SR at FIXED_CUTS =================
    yr = args.year if "2024" not in args.year else "2024"
    inr = lambda d: d[d.BDT_score >= l_edge]
    print("\n~~~~  SANITY CHECK (all yields above the fixed cut)  ~~~~")
    for key in [f"{m}_GeV" for m in range(15,65,5)]:
        hs, chk = smoothed_hists_mass[key], check_store[key]
        ax = hs["signal_SR"].GetXaxis()
        ib = ax.FindBin(FIXED_CUTS[key] + 1e-9); nb = hs["signal_SR"].GetNbinsX() + 1
        edge = ax.GetBinLowEdge(ib)
        tail = lambda d: d[d.BDT_score >= edge]
        bdf, ddf, sdf = bkg_samples_df[key], data_samples_df[key], sig_samples_df[key]
        win = windows[key]

        # ---- background: histograms (what the table uses) ----
        B_sm, B_un = hs["bkg_SR"].Integral(ib, nb), chk["bkg"].Integral(ib, nb)
        # ---- background: independent recomputation from the dataframe ----
        w = "1dimreweight"
        b_sr, b_wn = inr(SignalRegion(bdf)), inr(SignalRegion(bdf, window=win))
        d_sr_tot = chk["data_sr_total"]
        d_wn_tot = len(inr(SignalRegion(ddf, window=win)))
        B_sr       = tail(b_sr)[w].sum() * d_sr_tot / b_sr[w].sum()   # SR, normalized to data SR
        B_wn_pipe  = tail(b_wn)[w].sum() * d_sr_tot / b_sr[w].sum()   # window, scale taken from full SR (pipeline)
        B_wn_match = tail(b_wn)[w].sum() * d_wn_tot / b_wn[w].sum()   # window, normalized to data in window
        # ---- background: data-driven expectation from the sidebands ----
        d_sr, d_sb = inr(SignalRegion(ddf)), inr(SideBand(ddf))
        r, D_tail = len(d_sr) / len(d_sb), len(tail(d_sb))
        B_dd, B_dd_err = r * D_tail, r * np.sqrt(max(D_tail, 1))

        # ---- signal ----
        scale = Scales.lumi_fb / Scales.sumw[yr][key]
        S_sm, S_un = hs["signal_SR"].Integral(ib, nb), chk["sig"].Integral(ib, nb)
        s_sr, s_wn = inr(SignalRegion(sdf)), inr(SignalRegion(sdf, window=win))
        S_sr, S_wn = tail(s_sr).weight.sum() * scale, tail(s_wn).weight.sum() * scale
        eff = tail(s_sr).weight.sum() / s_sr.weight.sum()
        sw_ratio = signal[key]["weight"].sum() / Scales.sumw[yr][key]
        S_all = tail(inr(SignalRegion(signal[key]))).weight.sum() * scale   # all events, no split

        print(f"\n{key}  cut {FIXED_CUTS[key]:.3f} (bin edge {edge:.4f})")
        print(f"  BKG  smoothed {B_sm:7.3f} | unsmoothed {B_un:7.3f} | df: SR {B_sr:7.3f}  win(pipeline norm) {B_wn_pipe:7.3f}  win(matched norm) {B_wn_match:7.3f}")
        print(f"       data-driven (SB tail x SR/SB ratio) {B_dd:6.2f} +- {B_dd_err:5.2f}   |  totals smoothed/unsmoothed/data SR: {hs['bkg_SR'].Integral():.1f} / {chk['bkg'].Integral():.1f} / {d_sr_tot:.1f}")
        print(f"  SIG  smoothed {S_sm:7.3f} | unsmoothed {S_un:7.3f} | df: SR {S_sr:7.3f}  win {S_wn:7.3f}  (win<=SR: {'OK' if S_wn <= S_sr else 'PROBLEM'})  eff above cut {eff:.3f}")
        print(f"       sum(weight, all events)/sumw = {sw_ratio:.3f}   |   S(test)/S(all) = {S_sr / S_all:.3f}  (~1 if scale is split-corrected, ~{1/TEST_MOD:.2f} if not)")
    # =====================================================================================
    print("~~~~    Begin categorization procedure!   ~~~~")
    # Need to loop through them so I can do AMS plots together
    # With --fixed_cuts only a single category [cut, 1.0] per mass point is evaluated
    n_cats_list = [1] if args.fixed_cuts else range(args.nCats, 0, -1)
    for nCats in n_cats_list:
        significance_mass = {}
        partition_mass = {}

        for idx, key in enumerate([f"{m}_GeV" for m in range(15,65,5)]):
            print(f"~  Processing: {key.replace('_',' ')}  ~")
            bkg_sample_df = bkg_samples_df[key]
            data_sample_df = data_samples_df[key]
            # Loop through points from right to left
            nBins = smoothed_hists_mass[key]["bkg_SR"].GetNbinsX()
            significance_final = -999.
            partition_final = []
            plotting = []

            # Start first category at [0.05, 0.70] for nCats > 1 since there's no way the optimal cut is that far out
            # The first category would then be [bin 1, bin i]. For nCats > 1, bin i starts at BDT score 0.70 (see next line)
            nBins_start = int(nBins * 0.70)
            if nCats == 1:
                for i in range(1, nBins+1):
                    partition = [[i, nBins+1]]
                    significance, significance_all, S_out, B_out, Dsb_out = sumSignificance(
                        partition,
                        smoothed_hists_mass[key]["signal_SR"],
                        smoothed_hists_mass[key]["bkg_SR"],
                        #smoothed_hists_mass[key]["bkg_SB"],
                        #smoothed_hists_mass[key]["data_SB"]
                        SideBand(bkg_sample_df),
                        SideBand(data_sample_df),
                        return_all=True
                    )

                    plotting.append((significance_all, Dsb_out, partition[0][0], S_out, B_out))  # Keep boundary from i to 1 (BDT score)

                    if significance > significance_final:
                        #print(significance, significance_all, Dsb_out)
                        significance_final = significance
                        partition_final = partition

            elif nCats == 2:
                for i in range(nBins_start, nBins+1):
                    for j in range(i+1, nBins+1):
                        # Ensure that the "empty region" aligns with the start of the categories
                        if abs(j-i) == 1:
                            partition = [[1, i], [j, nBins+1]]
                            significance = sumSignificance(
                                partition,
                                smoothed_hists_mass[key]["signal_SR"],
                                smoothed_hists_mass[key]["bkg_SR"],
                                #smoothed_hists_mass[key]["bkg_SB"],
                                #smoothed_hists_mass[key]["data_SB"]
                                SideBand(bkg_sample_df),
                                SideBand(data_sample_df),
                            )

                            if significance > significance_final:
                                significance_final = significance
                                partition_final = partition

            elif nCats == 3:
                for i in range(nBins_start, nBins+1):
                    for j in range(i+1, nBins+1):
                        for k in range(j+1, nBins+1):
                            # Ensure that the "empty region" aligns with the start of the categories
                            if abs(j-i) == 1:
                                partition = [[1, i], [j, k-1], [k, nBins+1]]
                                significance = sumSignificance(
                                    partition,
                                    smoothed_hists_mass[key]["signal_SR"],
                                    smoothed_hists_mass[key]["bkg_SR"],
                                    #smoothed_hists_mass[key]["bkg_SB"],
                                    #smoothed_hists_mass[key]["data_SB"]
                                    SideBand(bkg_sample_df),
                                    SideBand(data_sample_df),
                                )

                                if significance > significance_final:
                                    significance_final = significance
                                    partition_final = partition

            elif nCats == 4:
                for i in range(nBins_start, nBins+1):
                    for j in range(i+1, nBins+1):
                        for k in range(j+1, nBins+1):
                            for l in range(k+1, nBins+1):
                                # Ensure that the "empty region" aligns with the start of the categories
                                if abs(j-i) == 1:
                                    partition = [[1, i], [j, k-1], [k, l-1], [l, nBins+1]]
                                    significance = sumSignificance(
                                        partition,
                                        smoothed_hists_mass[key]["signal_SR"],
                                        smoothed_hists_mass[key]["bkg_SR"],
                                        #smoothed_hists_mass[key]["bkg_SB"],
                                        #smoothed_hists_mass[key]["data_SB"]
                                        SideBand(bkg_sample_df),
                                        SideBand(data_sample_df),
                                    )

                                    if significance > significance_final:
                                        significance_final = significance
                                        partition_final = partition

            elif nCats == 5:
                for i in range(nBins_start, nBins+1):
                    for j in range(i+1, nBins+1):
                        for k in range(j+1, nBins+1):
                            for l in range(k+1, nBins+1):
                                for m in range(l+1, nBins+1):
                                    # Ensure that the "empty region" aligns with the start of the categories
                                    if abs(j-i) == 1:
                                        partition = [[1, i], [j, k-1], [k, l-1], [l, m-1], [m, nBins+1]]
                                        significance = sumSignificance(
                                            partition,
                                            smoothed_hists_mass[key]["signal_SR"],
                                            smoothed_hists_mass[key]["bkg_SR"],
                                            #smoothed_hists_mass[key]["bkg_SB"],
                                            #smoothed_hists_mass[key]["data_SB"]
                                            SideBand(bkg_sample_df),
                                            SideBand(data_sample_df),
                                        )

                                        if significance > significance_final:
                                            significance_final = significance
                                            partition_final = partition

            else:
                assert nCats <= 5, nCats

            # ---- Override the optimized partition with the user-defined fixed cut ----
            if args.fixed_cuts:
                h = smoothed_hists_mass[key]["signal_SR"]
                # +1e-9 avoids floating-point trouble when the cut sits exactly on a bin edge
                ib = h.GetXaxis().FindBin(FIXED_CUTS[key] + 1e-9)
                partition_final = [[ib, nBins+1]]
                significance_final = sumSignificance(
                    partition_final,
                    smoothed_hists_mass[key]["signal_SR"],
                    smoothed_hists_mass[key]["bkg_SR"],
                    SideBand(bkg_sample_df),
                    SideBand(data_sample_df),
                )
                print(f"[fixed cut] {key}: requested {FIXED_CUTS[key]:.4f}, "
                      f"bin {ib}, actual low edge {h.GetXaxis().GetBinLowEdge(ib):.4f}")

            significance_mass.update({key: significance_final})
            partition_mass.update({key: partition_final})

            # Plotting of AMS vs BDT score goes here
            if nCats == 1:
                plot_ams_score(plotting, key, smoothed_hists_mass[key]["signal_SR"], partition_mass[key][0][0], args.year, args.xgb_name, cwd)


        suffix = "_fixed" if args.fixed_cuts else ""
        cats_path = f"{cwd}/cats/{args.year}/{args.xgb_name}/cats{nCats}{suffix}.txt"
        with open(cats_path, "w") as f:
            pass
        for key in significance_mass.keys():
            sig_sample_df = sig_samples_df[key]
            bkg_sample_df = bkg_samples_df[key]
            data_sample_df = data_samples_df[key]

            significance = sumSignificance(
                partition_mass[key],
                smoothed_hists_mass[key]["signal_SR"],
                smoothed_hists_mass[key]["bkg_SR"],
                #smoothed_hists_mass[key]["bkg_SB"],
                #smoothed_hists_mass[key]["data_SB"],
                SideBand(bkg_sample_df),
                SideBand(data_sample_df),
            )

            with open(cats_path, "a") as f:
                print(f"\nMass Point: {key.replace('_',' ')}  # of cats: {nCats}  AMS: {significance:4.04f}  Signal Total: {smoothed_hists_mass[key]['signal_SR'].Integral(): .2f}")
                f.write(f"\nMass Point: {key.replace('_',' ')}  # of cats: {nCats}  AMS: {significance:4.04f}  Signal Total: {smoothed_hists_mass[key]['signal_SR'].Integral(): .2f}\n")

            for idx, pair in enumerate(partition_mass[key]):
                pair_translated = [smoothed_hists_mass[key]["signal_SR"].GetXaxis().GetBinLowEdge(pair[0]), smoothed_hists_mass[key]["signal_SR"].GetXaxis().GetBinLowEdge(pair[1])]
                S = smoothed_hists_mass[key]["signal_SR"].Integral(pair[0], pair[1])
                B = smoothed_hists_mass[key]["bkg_SR"].Integral(pair[0], pair[1])
                #B_sb = smoothed_hists_mass[key]["bkg_SB"].Integral(pair[0], pair[1])
                #D_sb = smoothed_hists_mass[key]["data_SB"].Integral(pair[0], pair[1])
                B_sb = sum(SideBand(bkg_sample_df)[(bkg_sample_df.BDT_score > pair_translated[0]) & (bkg_sample_df.BDT_score <= pair_translated[1])]["1dimreweight"])
                D_sb = len(SideBand(data_sample_df)[(data_sample_df.BDT_score > pair_translated[0]) & (data_sample_df.BDT_score <= pair_translated[1])])

                with open(cats_path, "a") as f:
                    print(f"Signal SR: {S:6.6}  Bkg SR: {float(B):4.3}  Bkg SB: {float(B_sb):4.3}  Data SB: {float(D_sb):4.3}  Sig: {sum(sig_sample_df[(sig_sample_df.BDT_score > pair_translated[0]) & (sig_sample_df.BDT_score <= pair_translated[1])].weight): 4.2f}")
                    f.write(f"Signal SR: {S:6.6}  Bkg SR: {float(B):4.3}  Bkg SB: {float(B_sb):4.3}  Data SB: {float(D_sb):4.3}  Sig: {sum(sig_sample_df[(sig_sample_df.BDT_score > pair_translated[0]) & (sig_sample_df.BDT_score <= pair_translated[1])].weight): 4.2f}\n")

            with open(cats_path, "a") as f:
                print(f"Cuts: ", end="")
                f.write(f"Cuts: ")
                for idx, partition in enumerate(partition_mass[key]):
                    print(f"{smoothed_hists_mass[key]['signal_SR'].GetXaxis().GetBinLowEdge(partition[0]):<4.04f}, {smoothed_hists_mass[key]['signal_SR'].GetXaxis().GetBinLowEdge(partition[1]):<4.04f}", end="")
                    f.write(f"{smoothed_hists_mass[key]['signal_SR'].GetXaxis().GetBinLowEdge(partition[0]):<4.04f}, {smoothed_hists_mass[key]['signal_SR'].GetXaxis().GetBinLowEdge(partition[1]):<4.04f}")
                    if idx+1 == len(partition_mass[key]):
                        print()
                        f.write("\n")
                    else:
                        print(", ", end="")
                        f.write(", ")


    print("\n\n  ~~~  Finished categorization procedure!  ~~~\n")

    #--------------------------------------------------------------------------------a check if ams evaluation on the test sample alone is the culprit------------------------------
    # =====================================================================================
    # EXTRA TABLE: same quantities, but computed on the WHOLE sample (no train/test split)
    # Uses FIXED_CUTS, the same WINDOWS and the same smoothing/normalization as above.
    # Nothing above is modified.
    # =====================================================================================
    print("\n~~~~  Building FULL-SAMPLE table (all signal and background events)  ~~~~")

    APPLY_WEIGHT_CUT_FULL = True   # same weight < 10 requirement as the main table; set False to keep all events

    def _smooth_retry(name, hist, smoothness, span):
        h_, g_ = SmoothHist(name, hist, smoothness, span)
        n_try = 0
        while g_.GetN() <= 1 and n_try < retries:
            h_, g_ = SmoothHist(name, hist, smoothness, span)
            n_try += 1
        assert g_.GetN() > 1, f"Failed to smooth {name}. N: {g_.GetN()}"
        return h_, g_

    full_rows = []
    for key in [f"{m}_GeV" for m in range(15, 65, 5)]:
        print(f"\n[full sample] {key.replace('_', ' ')}")

        # ---- all events, no split ----
        sig_full = signal[key].copy()
        if APPLY_WEIGHT_CUT_FULL:
            sig_full = sig_full[sig_full["weight"] < 10.0]
        sig_full = sig_full.reset_index(drop=True)
        bkg_full_df = background[key].copy().reset_index(drop=True)
        data_full_df = data[key]
        print(f"  events: signal {len(sig_full)}, background {len(bkg_full_df)}, data {len(data_full_df)}")

        # ---- reweight bkg to data (in memory only, nothing written to disk) ----
        bkg_full_df = reweight1Dim(bkg_full_df, data_full_df, "BDT_score")

        # ---- histograms (unique names so nothing existing is overwritten) ----
        h_sig = fillHist("BDT_score", SignalRegion(sig_full, window=windows[key]), [l_edge, r_edge],
                         name=f"BDT_score_signal_SR_FULLSAMPLE_{key}", binsScale=bins, normalize=False, preventOverFlow=True)
        h_bkg_win = fillHist("BDT_score", SignalRegion(bkg_full_df, window=windows[key]), [l_edge, 1.0],
                             name=f"BDT_score_bkg_SR_FULLSAMPLE_{key}", binsScale=bins, normalize=False,
                             customWeights="1dimreweight", preventOverFlow=True)
        h_bkg_sr = fillHist("BDT_score", SignalRegion(bkg_full_df), [l_edge, 1.0],
                            name=f"BDT_score_bkg_SRfull_FULLSAMPLE_{key}", binsScale=bins, normalize=False,
                            customWeights="1dimreweight", preventOverFlow=True)
        h_data_sr = fillHist("BDT_score", SignalRegion(data_full_df), [l_edge, 1.0],
                             name=f"BDT_score_data_SR_FULLSAMPLE_{key}", binsScale=bins, normalize=False, preventOverFlow=True)

        # ---- normalization (same as main pipeline) ----
        h_bkg_win.Scale(h_data_sr.Integral() / h_bkg_sr.Integral())
        h_sig.Scale(Scales.lumi_fb * 1.0 * 1.0 / Scales.sumw[yr][key])   # all events kept -> no split correction needed

        # ---- smooth (same settings as main pipeline) ----
        sig_h_sm, sig_g = _smooth_retry(f"sig_smooth_sig_FULLSAMPLE_{key}", h_sig, args.smoothness, args.span)
        bkg_h_sm, bkg_g = _smooth_retry(f"bkg_smooth_FULLSAMPLE_{key}", h_bkg_win, bkg_data_smoothness, bkg_data_span)

        # ---- fixed cut -> partition ----
        nBins_f = bkg_h_sm.GetNbinsX()
        if args.fixed_cuts:
            cut_f = FIXED_CUTS[key]
            cut_src = "fixed"
        else:
            # partition_mass still holds the nCats=1 result (last iteration of the loop above)
            cut_f = smoothed_hists_mass[key]["signal_SR"].GetXaxis().GetBinLowEdge(partition_mass[key][0][0])
            cut_src = "optimized"
        ib_f = sig_h_sm.GetXaxis().FindBin(cut_f + 1e-9)   # same binning as the main histograms
        part_f = [[ib_f, nBins_f + 1]]
        lo_edge = sig_h_sm.GetXaxis().GetBinLowEdge(ib_f)
        hi_edge = sig_h_sm.GetXaxis().GetBinLowEdge(nBins_f + 1)
        print(f"  [{cut_src} cut] {key}: {cut_f:.4f} -> bin {ib_f}, low edge {lo_edge:.4f}")

        ams_f = sumSignificance(part_f, sig_h_sm, bkg_h_sm, SideBand(bkg_full_df), SideBand(data_full_df))
        S_f = sig_h_sm.Integral(ib_f, nBins_f + 1)
        B_f = bkg_h_sm.Integral(ib_f, nBins_f + 1)
        D_sb_f = len(SideBand(data_full_df)[(data_full_df.BDT_score > lo_edge) & (data_full_df.BDT_score <= hi_edge)])

        # independent check: unweighted-by-hist signal efficiency above the cut in the SR window
        s_win = SignalRegion(sig_full, window=windows[key])
        eff_f = s_win[s_win.BDT_score >= lo_edge].weight.sum() / s_win.weight.sum()

        full_rows.append((key, lo_edge, hi_edge, S_f, B_f, D_sb_f, ams_f, eff_f, significance_mass.get(key, float("nan"))))

        sig_g.Delete()
        bkg_g.Delete()

    # ---- table ----
    header = f"{'m_a [GeV]':<10}{'Category':<16}{'Signal (SR)':>12}{'Background (SR)':>17}{'Data (SB)':>11}{'Significance':>14}{'sig eff>cut':>13}{'AMS (test-only)':>17}"
    lines = ["FULL-SAMPLE TABLE (all events, no train/test split)", header]
    lines[0] += f"  [BDT cuts: {'FIXED_CUTS' if args.fixed_cuts else 'optimized (main pipeline)'}]"
    for (key, lo, hi, S_f, B_f, D_f, A_f, e_f, A_test) in full_rows:
        lines.append(f"{key.replace('_GeV',' GeV'):<10}{f'{lo:.3f}, {hi:.0f}':<16}{S_f:>12.2f}{float(B_f):>17.2f}{D_f:>11d}{A_f:>14.2f}{e_f:>13.3f}{A_test:>17.2f}")
    table_txt = "\n".join(lines)
    print("\n" + table_txt)

    full_path = f"{cwd}/cats/{args.year}/{args.xgb_name}/table_fullsample{'_fixed' if args.fixed_cuts else ''}.txt"
    with open(full_path, "w") as f:
        f.write(table_txt + "\n")
    print(f"\nSaved full-sample table to {full_path}")
    # ----------------------------------------------------------------------------------------------------------------------lets see if it works-----------------------------
