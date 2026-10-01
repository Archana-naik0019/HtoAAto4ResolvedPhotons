import gc
import numpy as np
import copy
import argparse
import os
from ROOT import TMath
from h4g_tools.utils.scales import Scales as ScalesCls
from h4g_tools.bdt.event_selection_bdt import redefine_dataframe
from h4g_tools.utils.runner_utils import get_parser_plotting, handle_errors
from h4g_tools.utils.standard_plots import plotPsMasses, plot4Object, plotPsKinematics, plotPsMassCombined, plotHMassCombined, plotMHyp, plotdRCombined
from h4g_tools.utils.plotting import plot_data_sig_bkg_comparison, plot_data_sig_bkg_comparison_combined, compare_four_masses_data_sig_bkg, loadPlottingParameters, compare_eras_data_bkg_only, compare_four_masses_data_sig_bkg_reweight_overlay, plot_corr_matrix, plotFourMassesCompared, make_scatter_mva
from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist
from h4g_tools.bdt.reweight import loadNDimWeights


"""
Plots comparison plots for 
python3 plot_BDT_cuts_compared.py -y [2022, 2022preEE, 2022postEE]
"""



def plotBDTCutCompared(
    args: argparse.ArgumentParser
) -> None:

    handle_errors(
        ["Must provide a year for processing.", args.year is None]
    )

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    Scales = ScalesCls(args.year) 
    masses = [x for x in range(15, 65, 5)]

    # Handle path setup after checking year
    pathSetup(args.year)

    print("Loading signal...")
    signal = loadPerMass(os.path.join(*(cwd, "outputs/BDT", args.year, "signal")), skipInterMass=True)
    print()
    print("Loading background...")
    bkg = loadPerMass(os.path.join(*(cwd, "outputs/BDT", args.year, "background")), skipInterMass=True)
    print()
    print("Loading data...")
    data = loadPerMass(os.path.join(*(cwd, "outputs/BDT", args.year, "data")), skipInterMass=True)
    print()
    
    # Cut on BDT score
    bdt_cut = 0.5
    assert set(signal.keys()) == set(bkg.keys()) == set(data.keys())
    for key in signal.keys():
        if "BDT_score" not in signal[key].columns or "BDT_score" not in bkg[key].columns or "BDT_score" not in data[key].columns:
            print(f"No column for BDT score found in {key}. Continuing without cuts.")
            return

        signal[key] = signal[key][signal[key]["BDT_score"] > bdt_cut]
        bkg[key] = bkg[key][bkg[key]["BDT_score"] > bdt_cut]
        data[key] = data[key][data[key]["BDT_score"] > bdt_cut]

    # Make scatter plot of photon 3/4 MVA IDs
    mva_samples = {}
    for key in signal.keys():
        if "Signal" in key:
            continue
        mva_samples.update({f"Signal_{key}": signal[key]})
    for key in bkg.keys():
        mva_samples.update({f"background_{key}": bkg[key]})
    for key in data.keys():
        mva_samples.update({f"data_{key}": data[key]})
    make_scatter_mva(mva_samples, args.year, bdt_cut=True, cut=bdt_cut)

    # Setup for lists
    branches, boundsList, binsScaleList, names, norms, titles = loadPlottingParameters()
     
    # Sample setup for per mass comparisons
    max_order = 6
    weights = bkg["15_GeV"]["Ndimreweight"]
    for canvas_num, mass in enumerate(signal.keys()):
        paths = [
            f"./plots/BDT/{args.year}/cuts/{mass}/pho1_mvaID_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/pho2_mvaID_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/pho3_mvaID_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/pho4_mvaID_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/LeadPs_pt_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/SubleadPs_pt_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/dR_aa_mass_gggg_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/LeadPs_interMass_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/SubleadPs_interMass_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/LeadPs_interMass_abs_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/SubleadPs_interMass_abs_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/Ps_massDiff_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/cos_ag_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/pT1_ma1_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/pT2_ma1_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/pT1_ma2_compare.png",
            f"./plots/BDT/{args.year}/cuts/{mass}/pT2_ma2_compare.png",
        ]

        assert len(branches) == len(boundsList) == len(binsScaleList) == len(names) == len(norms) == len(titles) == len(paths), f"{len(branches)} {len(boundsList)} {len(binsScaleList)} {len(names)} {len(norms)} {len(titles)} {len(paths)}"

        # Fill histograms
        for branch, bounds, binsScale, name, norm, title, path in zip(branches, boundsList, binsScaleList, names, norms, titles, paths):
            hists_signal = {}
            hists_bkg = {}
            hists_bkg_rw = {}
            hists_data = {}

            hists_signal.update({mass: fillHist(branch, signal[mass], bounds, binsScale=binsScale, name=f"{name}_{mass}_signal", normalize=False)})
            hists_bkg.update({mass: fillHist(branch, bkg[mass], bounds, binsScale=binsScale, name=f"{name}_{mass}_bkg", normalize=False)})
            hists_bkg_rw.update({mass: fillHist(branch, bkg[mass], bounds, binsScale=binsScale, name=f"{name}_{mass}_bkg", normalize=False, customWeights="Ndimreweight")})
            hists_data.update({mass: fillHist(branch, data[mass], bounds, binsScale=binsScale, name=f"{name}_{mass}_data", normalize=False)})
            
            # Plot a comparison for each mHyp (corresponding to pseudoscalar mass points)
            # No bkg RW
            hists = {"background": hists_bkg[mass], "data": hists_data[mass], "signal": hists_signal[mass]}
            plot_data_sig_bkg_comparison(copy.deepcopy(hists), title, path.replace(".png", "_pre_reweight.png"), Scales.sumw[args.year[:4]][mass[:6]], log=True, maximum=10**max_order - 5*10**(max_order-1), mass=mass[:6].replace("_"," "), canvas_num=str(canvas_num)+branch, year=args.year, Scales=Scales)

            # With bkg RW
            hists = {"background": hists_bkg_rw[mass], "data": hists_data[mass], "signal": hists_signal[mass]}
            plot_data_sig_bkg_comparison(copy.deepcopy(hists), title, path.replace(".png", "_reweight.png"), Scales.sumw[args.year[:4]][mass[:6]], log=True, maximum=10**max_order - 5*10**(max_order-1), mass=mass[:6].replace("_"," "), canvas_num=str(canvas_num)+branch, year=args.year, Scales=Scales)
    
    # Load flat mHyp bkg and data samples
    samples_general = signal.copy()
    for key, inputs in zip(["background", "data"], [os.path.join(*(cwd, f"outputs/BDT/{args.year}/generic/background")), os.path.join(*(cwd, f"outputs/BDT/{args.year}/generic/data"))]):
        print(f"Loading {key}...")
        assert os.path.exists(inputs), inputs
        samples_general.update({key: loadSamples(inputs)})

    paths = [
        f"./plots/BDT/{args.year}/cuts/compare/pho1_mvaID_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/pho2_mvaID_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/pho3_mvaID_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/pho4_mvaID_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/LeadPs_pt_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/SubleadPs_pt_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/dR_aa_mass_gggg_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/LeadPs_interMass_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/SubleadPs_interMass_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/LeadPs_interMass_abs_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/SubleadPs_interMass_abs_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/Ps_massDiff_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/cos_ag_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/pT1_ma1_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/pT2_ma1_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/pT1_ma2_compare.png",
        f"./plots/BDT/{args.year}/cuts/compare/pT2_ma2_compare.png",
    ]

    # Plot reweighted versions in entire range
    # Fill hists
    # Comparing only four maases
    for branch, bounds, binsScale, name, norm, title, path in zip(branches, boundsList, binsScaleList, names, norms, titles, paths):
        hists = {}
        hists_rw = {}
        for key, sample in samples_general.items():
            if "background" in key or "data" in key:
                continue
            tmp = fillHist(branch, sample, bounds, binsScale=binsScale, name=f"{name}_{key}", normalize=norm)
            hists.update({key: tmp})
        
        # Normalizes histogram by lumi * xs(1fb) * BR(1.0) * efficiency (after-sel/pre-sel genWeight) since histograms are already scaled to after-sel weight
        # NOTE: weight column is based on genWeight with adjustments from corrections and/or systematics
        gw = [Scales.sumw[args.year[:4]][f"{mass}_GeV"] for mass in [15, 30, 40, 60]]

        # Plot reweighted versions in entire range
        # Fill hists
        bkg_hist = fillHist(branch, samples_general["background"], bounds, binsScale=binsScale, name=f"{name}_background", normalize=norm)
        data_hist = fillHist(branch, samples_general["data"], bounds, binsScale=binsScale, name=f"{name}_data", normalize=norm)
        if not args.noRW:
            bkg_rw = fillHist(branch, samples_general["background"], bounds, binsScale=binsScale, name=f"{name}_background_rw", normalize=norm, customWeights="Ndimreweight")

        bkg_hist.Sumw2()
        data_hist.Sumw2()
        if not args.noRW:
            bkg_rw.Sumw2()

        # Reorder hists dictionary to play nice with the colors
        hists_pre_rw = {
            "background": bkg_hist,
            "data": data_hist,
            "Signal_15_GeV": hists["15_GeV"],
            "Signal_30_GeV": hists["30_GeV"],
            "Signal_40_GeV": hists["40_GeV"],
            "Signal_60_GeV": hists["60_GeV"]
        }

        if not args.noRW:
            hists_rw = {
                "background": bkg_rw,
                "data": data_hist,
                "Signal_15_GeV": hists["15_GeV"],
                "Signal_30_GeV": hists["30_GeV"],
                "Signal_40_GeV": hists["40_GeV"],
                "Signal_60_GeV": hists["60_GeV"]
            }

        compare_four_masses_data_sig_bkg(copy.deepcopy(hists_pre_rw), title.replace(' [GeV]'," [GeV]"), path[:-4]+"_pre_reweight.png", gw, log=True, maximum=10**max_order - 5*10**(max_order-1), canvas_num=branch+"entire_region_prerw", ignore_scales=False, year=args.year, Scales=Scales)

        new_title = title + " with Reweighting"
        if "[GeV]" in new_title:
            new_title = new_title.replace("[GeV]","") + " [GeV]"
        if "[GeV^{-1}]" in new_title:
            new_title = new_title.replace("[GeV^{-1}]","") + " [GeV^{-1}]"
        if not args.noRW:
            compare_four_masses_data_sig_bkg(copy.deepcopy(hists_rw), new_title, path[:-4]+"_reweight.png", gw, log=True, maximum=10**max_order - 5*10**(max_order-1), canvas_num=branch+"entire_region_rw", ignore_scales=False, year=args.year, Scales=Scales)

    """
    # Make scatter plot of photon 3/4 MVA IDs
    mva_samples = {}
    for key in signal.keys():
        if "Signal" in key:
            continue
        mva_samples.update({f"Signal_{key}": signal[key]})
    mva_samples.update({f"background": samples_general["background"]})
    mva_samples.update({f"data": samples_general["data"]})
    make_scatter_mva(mva_samples, args.year, bdt_cut=True, cut=bdt_cut)
    """



if __name__ == "__main__":
    parser = get_parser_plotting()
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    args = parser.parse_args()

    plotBDTCutCompared(args)
