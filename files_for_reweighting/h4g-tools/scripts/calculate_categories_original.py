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
gROOT.LoadMacro("../h4g_tools/bdt/smoother.cpp")
from ROOT import SmoothHist, SplitSignalHistograms, MergeHistograms, MergeGraphs
from operator import itemgetter
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# No popups with matplotlib
matplotlib.use('Agg')


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

    # Run pathSetup and get parser
    parser = get_parser_cats()
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    parser.add_argument("-xgb", "--xgb_name", type=str, required=True, help="XGBoost model name.")
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
    print(f"{cwd}/../../scripts/outputs/BDT/{args.year}/{args.xgb_name}/signal/")
    signal = loadPerMass(f"{cwd}/../../scripts/outputs/BDT/{args.year}/{args.xgb_name}/signal/", ["BDT_score", "mass_gggg", "weight"], bdt=True)
    print("\nLoading EM Background samples:")
    print(f"{cwd}/../../scripts/outputs/BDT/{args.year}/{args.xgb_name}/background/")
    background = loadPerMass(f"{cwd}/../../scripts/outputs/BDT/{args.year}/{args.xgb_name}/background/", ["BDT_score", "mass_gggg"], bdt=True)
    print("\nLoading Data samples:")
    print(f"{cwd}/../../scripts/outputs/BDT/{args.year}/{args.xgb_name}/data/")
    data = loadPerMass(f"{cwd}/../../scripts/outputs/BDT/{args.year}/{args.xgb_name}/data/", ["BDT_score", "mass_gggg"], bdt=True, nominalMassOnly=True)

    print("\nChecking if inputs are transformed:")
    transformed = False
    check_transformed = loadSamples(f"{cwd}/../../scripts/outputs/BDT/{args.year}/{args.xgb_name}/data/")
    if "BDT_score_orig" in check_transformed:
        transformed = True
    print(f"Transformed: {transformed}")

    # Cut off l_edge for smoothing
    l_edge = 0.05 if not transformed else 0.0
    r_edge = 1.0
    bins = 200.0
    retries = 10

    # 2sigma windows for calculating S and B
    mean = {
        "15_GeV": 124.4445,
        "20_GeV": 124.3937,
        "25_GeV": 124.5140,
        "30_GeV": 124.4907,
        "35_GeV": 124.5174,
        "40_GeV": 124.5717,
        "45_GeV": 124.2271,
        "50_GeV": 124.5182,
        "55_GeV": 124.3836,
        "60_GeV": 124.2573
    }
    res = {
        "15_GeV": 2.48,
        "20_GeV": 2.46,
        "25_GeV": 2.45,
        "30_GeV": 2.41,
        "35_GeV": 2.44,
        "40_GeV": 2.46,
        "45_GeV": 2.45,
        "50_GeV": 2.47,
        "55_GeV": 2.46,
        "60_GeV": 2.50
    }
    windows = {}
    for key in mean.keys():
        windows.update({key: [mean[key] - 2*res[key], mean[key] + 2*res[key]]})

    sig_samples_df = {}
    bkg_samples_df = {}
    data_samples_df = {}
    smoothed_hists_mass = {}
    print(f"\nNumber of categories set to {args.nCats}.")
    print(f"Number of bins set to {bins} (minus those left of {l_edge}).")
    for idx, key in enumerate([f"{m}_GeV" for m in range(15,65,5)]):
        print(f"\n\n ~~~~~~~  Starting smoothing procedure for mass point: {key.replace('_',' ')}  ~~~~~~~")
        print("Keeping test (even) events only for signal and background.")
        #sig_sample_df = signal[key]
        #bkg_sample_df = background[key]
        sig_sample_df = signal[key].iloc[lambda x: x.index % 2 == 0]  # Even rows only (test set)
        bkg_sample_df = background[key].iloc[lambda x: x.index % 2 == 0]  # Even rows only (test set)
        sig_sample_df.reset_index(inplace=True)
        bkg_sample_df.reset_index(inplace=True)
        data_sample_df = data[key]

        print("Reweighting background to data...")
        bkg_sample_df = reweight1Dim(bkg_sample_df, data_sample_df, "BDT_score")

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
        #sig_hist_sig = fillHist(f"BDT_score", SignalRegion(sig_sample_df, window=windows[key]), [l_edge, r_edge], name=f"BDT_score_signal_SR_{key}", binsScale=bins, normalize=False, preventOverFlow=True)
        # NO WINDOWS:
        sig_hist_sig = fillHist(f"BDT_score", SignalRegion(sig_sample_df), [l_edge, r_edge], name=f"BDT_score_signal_SR_{key}", binsScale=bins, normalize=False, preventOverFlow=True)
        # Signal in sideband region
        sig_sideband = fillHist(f"BDT_score", SideBand(sig_sample_df), [l_edge, r_edge], name=f"BDT_score_signal_SB_{key}", binsScale=bins, normalize=False, preventOverFlow=True)

        # Fill EM background hist
        print("Filling EM background hist...")
        # EM background in signal region
        # WITH WINDOWS:
        #bkg_hist_sig = fillHist(f"BDT_score", SignalRegion(bkg_sample_df, window=windows[key]), [l_edge, 1.0], name=f"BDT_score_bkg_SR_{key}", binsScale=bins, normalize=False, customWeights="1dimreweight", preventOverFlow=True)
        # NO WINDOWS:
        bkg_hist_sig = fillHist(f"BDT_score", SignalRegion(bkg_sample_df), [l_edge, 1.0], name=f"BDT_score_bkg_SR_{key}", binsScale=bins, normalize=False, customWeights="1dimreweight", preventOverFlow=True)
        # EM background in sideband region
        bkg_sideband = fillHist(f"BDT_score", SideBand(bkg_sample_df), [l_edge, 1.0], name=f"BDT_score_bkg_SB_{key}", binsScale=bins, normalize=False, customWeights="1dimreweight", preventOverFlow=True)

        # Fill data hist
        print("Filling data hist...")
        print(len(data_sample_df), len(SignalRegion(data_sample_df)), len(SideBand(data_sample_df)))
        # Data in signal region for normalization
        data_hist_sig = fillHist(f"BDT_score", SignalRegion(data_sample_df), [l_edge, 1.0], name=f"BDT_score_data_SR_{key}", binsScale=bins, normalize=False, preventOverFlow=True)
        # Data in sideband region
        data_sideband = fillHist(f"BDT_score", SideBand(data_sample_df), [l_edge, 1.0], name=f"BDT_score_data_SB_{key}", binsScale=bins, normalize=False, preventOverFlow=True)
        
        # Normalize background to data
        bkg_hist_sig.Scale(data_hist_sig.Integral() / bkg_hist_sig.Integral())
        bkg_sideband.Scale(data_sideband.Integral() / bkg_sideband.Integral())

        # lumi * xs(1fb) * BR(100%) / sum(# of generated events)
        sig_hist_sig.Scale(Scales.lumi_fb * 1.0 * 1.0 / Scales.sumw[args.year if "2024" not in args.year else "2024"][key])
        sig_sideband.Scale(Scales.lumi_fb * 1.0 * 1.0 / Scales.sumw[args.year if "2024" not in args.year else "2024"][key])


        ########################################################################
        ### Smooth signal region histograms
        bkg_data_smoothness = 8
        #bkg_data_smoothness = 7.5
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
        """
        split = 0.85
        sig_sideband1, sig_sideband2 = SplitSignalHistograms(sig_sideband, split, f"{cwd}/plots/BDT/{args.year}/cats/{args.xgb_name}/merged_hists_{key}.pdf")
        sig_sb_graph_hist1, sig_sb_graph1 = SmoothHist(f"sig_smooth_sb1_{key}", sig_sideband1, args.smoothness, args.span)
        sig_sb_graph_hist2, sig_sb_graph2 = SmoothHist(f"sig_smooth_sb2_{key}", sig_sideband2, args.smoothness, 0.02)
        sig_sb_graph_hist = MergeHistograms(sig_sb_graph_hist1, sig_sb_graph_hist2)
        sig_sb_graph = MergeGraphs(sig_sb_graph1, sig_sb_graph2, split, f"{cwd}/plots/BDT/{args.year}/cats/{args.xgb_name}/merged_graphs_{key}.pdf")
        plot_smooth_all(
            f"{cwd}/plots/BDT/{args.year}/cats/{args.xgb_name}/sig_hists_{key}.pdf",
            key.replace("_", " "),
            [sig_sideband1, sig_sb_graph1, kRed, kOrange, "hist1", "l", "hist same"],
            [sig_sideband2, sig_sb_graph2, kBlue, kGreen+3, "hist2", "l", "hist same"],
            [sig_sideband, sig_sb_graph, kBlack, kMagenta, "combined hist", "p", "hist P same"],
            Scales = Scales,
        )
        #"""
        #"""
        sig_sb_graph_hist, sig_sb_graph = SmoothHist(f"sig_smooth_sb_{key}", sig_sideband, args.smoothness, args.span)
        count_retries = 0
        while sig_sb_graph.GetN() <= 1 and count_retries < retries:
            print("Retrying to smooth background.")
            sig_sb_graph_hist, sig_sb_graph = SmoothHist(f"sig_smooth_sb_{key}", sig_sideband, args.smoothness, args.span)
            count_retries += 1
        #"""
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

        """
        print("Plotting smoothed BDT distributions...")
        plot_smooth(bkg_sb_graph, bkg_sideband, sig_sideband, data_sideband, key, f"cat_bkg_{key.replace('Signal_','')}_smooth.pdf", [[0]], l_edge, Scales=Scales)
        """

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
    print("~~~~    Begin categorization procedure!   ~~~~")
    # Need to loop through them so I can do AMS plots together
    for nCats in range(args.nCats, 0, -1):
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

            significance_mass.update({key: significance_final})
            partition_mass.update({key: partition_final})

            # Plotting of AMS vs BDT score goes here
            if nCats == 1:
                plot_ams_score(plotting, key, smoothed_hists_mass[key]["signal_SR"], partition_mass[key][0][0], args.year, args.xgb_name, cwd)


        cats_path = f"{cwd}/cats/{args.year}/{args.xgb_name}/cats{nCats}.txt"
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
