import gc
import os
import numpy as np
import copy
import argparse
from h4g_tools.bdt.event_selection_bdt import redefine_dataframe
from h4g_tools.utils.runner_utils import get_parser_plotting, handle_errors
from h4g_tools.utils.plotting import compare_eras_data_bkg_only, loadPlottingParameters
from h4g_tools.utils.loading import loadSamples, fillHist
from h4g_tools.utils.scales import Scales as ScalesCls

"""
python3 plot_preEE_postEE_comparison.py -bkg ~/Private/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_bkg2022_noCorr_noSyst_noMass55_Nov2/ -d ~/Private/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_data2022_noCorr_noSyst_noMass55_Nov2/
"""


def compare(
    args: argparse.ArgumentParser
) -> None:

    ## SET YEAR FOR PLOTS!
    year = "2022"

    # Scales instance setup
    Scales = ScalesCls(year) 
    max_order = 6

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))
    if not os.path.exists(os.path.join(cwd, "plots/BDT/2022/", "compare_prepost")):
        os.mkdir(os.path.join(cwd, "plots/BDT/2022/", "compare_prepost"))

    if args.background is None or args.data is None:
        return

    # Load samples
    samples = {}
    masses = [x for x in range(15, 65, 5)]
    for key, inputs in zip(["bkg_preEE", "data_preEE"], [args.background, args.data]):
        print(f"Loading {key}...")
        eras_preEE = {
            "data": ["Run2022C","Run2022D"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_preEE" for m in masses]
        }
        samples.update({key: loadSamples(inputs, eras=eras_preEE)})
    for key, inputs in zip(["bkg_postEE", "data_postEE"], [args.background, args.data]):
        eras_postEE = {
            "data": ["Run2022E","Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_preEE" for m in masses]
        }
        print(f"Loading {key}...")
        samples.update({key: loadSamples(inputs, eras=eras_postEE)})

    # Setup for lists
    branches, boundsList, binsScaleList, names, norms, titles = loadPlottingParameters()

    # Collect previous objects and disable circular garbage collector (ref counter still works)
    gc.collect()
    gc.disable()

    paths = [
        f"./plots/BDT/2022/compare_prepost/pho1_mvaID_compare.png",
        f"./plots/BDT/2022/compare_prepost/pho2_mvaID_compare.png",
        f"./plots/BDT/2022/compare_prepost/pho3_mvaID_compare.png",
        f"./plots/BDT/2022/compare_prepost/pho4_mvaID_compare.png",
        f"./plots/BDT/2022/compare_prepost/LeadPs_pt_compare.png",
        f"./plots/BDT/2022/compare_prepost/SubleadPs_pt_compare.png",
        f"./plots/BDT/2022/compare_prepost/dR_aa_mass_gggg_compare.png",
        f"./plots/BDT/2022/compare_prepost/LeadPs_interMass_compare.png",
        f"./plots/BDT/2022/compare_prepost/SubleadPs_interMass_compare.png",
        f"./plots/BDT/2022/compare_prepost/LeadPs_interMass_abs_compare.png",
        f"./plots/BDT/2022/compare_prepost/SubleadPs_interMass_abs_compare.png",
        f"./plots/BDT/2022/compare_prepost/Ps_massDiff_compare.png",
        f"./plots/BDT/2022/compare_prepost/cos_ag_compare.png",
    ]

    # Fill histograms
    for branch, bounds, binsScale, name, norm, title, path in zip(branches, boundsList, binsScaleList, names, norms, titles, paths):
        hists = {}
        hypMass_preEE = np.random.choice([float(x) for x in range(15, 65, 5)], len(samples["data_preEE"]) if len(samples["data_preEE"]) > len(samples["bkg_preEE"]) else len(samples["bkg_preEE"]))
        hypMass_postEE = np.random.choice([float(x) for x in range(15, 65, 5)], len(samples["data_postEE"]) if len(samples["data_postEE"]) > len(samples["bkg_postEE"]) else len(samples["bkg_postEE"]))
        for key, sample in samples.items():
            if "preEE" in key:
                #sample_df = redefine_dataframe(sample, hypMass_preEE, rand=True)
                sample_df = redefine_dataframe(sample, hypMass_postEE[:len(hypMass_preEE)], rand=True)
            elif "postEE" in key:
                sample_df = redefine_dataframe(sample, hypMass_postEE, rand=True)
            hists.update({key: fillHist(branch, sample_df, bounds, binsScale=binsScale, name=f"{name}_{key}", normalize=norm)})

        compare_eras_data_bkg_only(copy.deepcopy(hists), title, path[:-4]+"_preEE_postEE.png", log=True, maximum=10**max_order, canvas_num=branch+"entire_region", year=year, Scales=Scales)

if __name__ == "__main__":
    parser = get_parser_plotting()
    args = parser.parse_args()

    compare(args)
