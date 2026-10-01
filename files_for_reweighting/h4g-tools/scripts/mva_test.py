import numpy as np
import copy
import argparse
import os
from h4g_tools.utils.plotting import compare_four_masses_data_sig_bkg, loadPlottingParameters
from h4g_tools.utils.scales import Scales as ScalesCls
from h4g_tools.bdt.event_selection_bdt import redefine_dataframe
from h4g_tools.utils.runner_utils import get_parser_plotting, handle_errors
from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist

_inputs = "/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal_method1_SaS_Pileup_NFs_doublePreselDiphoton_11Jun2025/"
_bkg = "/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_bkg_method1_doublePreselDiphoton_11Jun2025/"
_data = "/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_data_method1_doublePreselDiphoton_11Jun2025/"


def plot(inputs, background, data, year = "2018"):

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    Scales = ScalesCls(year) 
    masses = [x for x in range(15, 65, 5)]
    if year == "2022preEE":
        eras={
            "data": ["Run2022C","Run2022D"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_preEE" for m in masses]
        }
    elif year == "2022postEE":
        eras={
            "data": ["Run2022E","Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_postEE" for m in masses]
        }
    elif year == "2022":
        eras={
            "data": ["Run2022C","Run2022D", "Run2022E","Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D", "RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preEE", "postEE"]]
        }
    else:
        assert False, "Must specify year in Scales class!"

    # Handle path setup after checking year
    pathSetup(year)

    # Setup for lists
    branches, boundsList, binsScaleList, names, norms, titles = loadPlottingParameters()
    max_order = 6

    print("Loading samples...")
    samples = loadPerMass(inputs, eras = eras, year = year)

    # Load samples
    samples_mHyp = {}
    for key, inputs in zip(["background", "data"], [background, data]):
        print(f"Loading {key}...")
        assert os.path.exists(inputs)
        samples_mHyp.update({key: loadSamples(inputs, eras = eras)})

    # Sample setup for general comparisons
    samples_general = samples_mHyp
    for canvas_num, mass in enumerate(samples.keys()):
        # Replace signal with new masses
        for m in ["Signal_15_GeV", "Signal_30_GeV", "Signal_40_GeV", "Signal_60_GeV"]:
            if m in mass:
                samples_general.update({m: samples[mass]})

    paths = [
        f"./plots/BDT/{year}/compare/pho1_mvaID_compare.png",
        f"./plots/BDT/{year}/compare/pho2_mvaID_compare.png",
        f"./plots/BDT/{year}/compare/pho3_mvaID_compare.png",
        f"./plots/BDT/{year}/compare/pho4_mvaID_compare.png",
        f"./plots/BDT/{year}/compare/LeadPs_pt_compare.png",
        f"./plots/BDT/{year}/compare/SubleadPs_pt_compare.png",
        f"./plots/BDT/{year}/compare/dR_aa_mass_gggg_compare.png",
        f"./plots/BDT/{year}/compare/LeadPs_interMass_compare.png",
        f"./plots/BDT/{year}/compare/SubleadPs_interMass_compare.png",
        f"./plots/BDT/{year}/compare/LeadPs_interMass_abs_compare.png",
        f"./plots/BDT/{year}/compare/SubleadPs_interMass_abs_compare.png",
        f"./plots/BDT/{year}/compare/Ps_massDiff_compare.png",
        f"./plots/BDT/{year}/compare/cos_ag_compare.png",
    ]

    branches = branches[2:4]
    boundsList = boundsList[2:4]
    binsScaleList = binsScaleList[2:4]
    names = names[2:4]
    norms = norms[2:4]
    titles = titles[2:4]
    paths = paths[2:4]
    assert len(branches) == len(boundsList) == len(binsScaleList) == len(names) == len(norms) == len(titles) == len(paths), f"{len(branches)} {len(boundsList)} {len(binsScaleList)} {len(names)} {len(norms)} {len(titles)} {len(paths)}"

    # Plot reweighted versions in entire range
    # Fill hists
    # Comparing only four maases
    for branch, bounds, binsScale, name, norm, title, path in zip(branches, boundsList, binsScaleList, names, norms, titles, paths):
        hists = {}
        hists_rw = {}
        for key, sample in samples_general.items():
            sample_df = sample
            print(f"key: {key}\t # events in SB: {len(sample[((sample.mass_gggg >= 110.0) & (sample.mass_gggg <= 115.0)) | ((sample.mass_gggg >= 135.0) & (sample.mass_gggg <= 180.0))])}\t # events in SR: {len(sample[(sample.mass_gggg > 115.0) & (sample.mass_gggg < 135.0)])}")
            if "background" in key or "data" in key:
                continue
            tmp = fillHist(branch, sample_df, bounds, binsScale=binsScale, name=f"{name}_{key}", normalize=norm)
            hists.update({key: tmp})
        
        # Normalizes histogram by lumi * xs(1fb) * BR(1.0) * efficiency (after-sel/pre-sel genWeight) since histograms are already scaled to after-sel weight
        # NOTE: weight column is based on genWeight with adjustments from corrections and/or systematics
        gw = [Scales.sumw[year[:4]][f"{mass}_GeV"] for mass in [15, 30, 40, 60]]

        # Fill hists
        bkg_hist = fillHist(branch, samples_general["background"], bounds, binsScale=binsScale, name=f"{name}_background", normalize=norm)
        data_hist = fillHist(branch, samples_general["data"], bounds, binsScale=binsScale, name=f"{name}_data", normalize=norm)
        bkg_hist.Sumw2()
        data_hist.Sumw2()

        # Reorder hists dictionary to play nice with the colors
        hists_pre_rw = {
            "background": bkg_hist,
            "data": data_hist,
            "Signal_15_GeV": hists["Signal_15_GeV"],
            "Signal_30_GeV": hists["Signal_30_GeV"],
            "Signal_40_GeV": hists["Signal_40_GeV"],
            "Signal_60_GeV": hists["Signal_60_GeV"]
        }

        # Reweighting normalizes data to bkg for some reason...

        # Remove first bin to see what that does to normalization
        hists_pre_rw["background"].SetBinContent(1, 0)
        hists_pre_rw["data"].SetBinContent(1, 0)

        print(title, "before norm", hists_pre_rw["background"].Integral(), hists_pre_rw["data"].Integral())
        print("bin 0", "bkg", hists_pre_rw["background"].Integral(0, 1), "data", hists_pre_rw["data"].Integral(0, 1))
        print("bin 1", "bkg", hists_pre_rw["background"].Integral(1, 2), "data", hists_pre_rw["data"].Integral(1, 2))
        print("bin 2-30", "bkg", hists_pre_rw["background"].Integral(2, 30), "data", hists_pre_rw["data"].Integral(2, 30))
        print("bin 31", "bkg", hists_pre_rw["background"].Integral(30, 31), "data", hists_pre_rw["data"].Integral(30, 31))

        hists = compare_four_masses_data_sig_bkg(copy.deepcopy(hists_pre_rw), title.replace(' [GeV]'," [GeV]"), path[:-4]+"_pre_reweight.png", gw, log=True, maximum=10**max_order - 5*10**(max_order-1), canvas_num=branch+"entire_region_prerw", ignore_scales=False, year=year, Scales=Scales)

        print(title, "after norm", hists["background"].Integral(), hists["data"].Integral())
        print("bin 0", "bkg", hists["background"].Integral(0, 1), "data", hists["data"].Integral(0, 1))
        print("bin 1", "bkg", hists["background"].Integral(1, 2), "data", hists["data"].Integral(1, 2))
        print("bin 2-30", "bkg", hists["background"].Integral(2, 30), "data", hists["data"].Integral(2, 30))
        print("bin 31", "bkg", hists["background"].Integral(30, 31), "data", hists["data"].Integral(30, 31))


if __name__ == "__main__":
    parser = get_parser_plotting()
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    args = parser.parse_args()

    plot(args.inputs, args.background, args.data, args.year)
