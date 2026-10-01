import gc
import numpy as np
import copy
import argparse
import os
from ROOT import TMath
from h4g_tools.utils.scales import Scales as ScalesCls
from h4g_tools.bdt.event_selection_bdt import redefine_dataframe
from h4g_tools.utils.runner_utils import handle_errors
from h4g_tools.utils.plotting import plotCorrectedComparison, loadPlottingParameters
from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist


"""
For comparing per-event corrections for the BDT training variables + plus others:
python3 compare_perEvent_corrections.py --sigUC <sig-uncorrected> --sigC <sig-corrected> --bkgUC <bkg-uncorrected> --bkgC <bkg-corrected> --dataUC <data-uncorrected> --dataC <data-corrected> --year [2022, 2022preEE, 2022postEE]

Signal:
--sigUC /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/merged_outputs_signal_noCorr_noSyst_PsuedosArray_Jan30/ --sigC /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/merged_outputs_signal_MCSmear_DataEtDepScale_wSyst_PsuedosArray_Jan30/

Bkg:
--bkgUC /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/merged_outputs_bkg_noCorr_noSyst_PsuedosArray_Jan31/ --bkgC /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/merged_outputs_bkg_MCSmear_DataEtDepScale_wSyst_PsuedosArray_Jan31/

Data:
--dataUC /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/merged_outputs_data_noCorr_noSyst_PsuedosArray_Jan31/ --dataC /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/merged_outputs_data_MCSmear_DataEtDepScale_wSyst_PsuedosArray_Jan31/


Do not perform reweighting here! Want to show the difference in corrected vs uncorrected samples.
"""

def compareCorrections(
    args: argparse.ArgumentParser
) -> None:

    handle_errors(
        ["Need both corrected and uncorrected signal samples.", (args.sigUC is not None and args.sigC is None) or (args.sigUC is None and args.sigC is not None)],
        ["Need both corrected and uncorrected bkg samples.", (args.bkgUC is not None and args.bkgC is None) or (args.bkgUC is None and args.bkgC is not None)],
        ["Need both corrected and uncorrected data samples.", (args.dataUC is not None and args.dataC is None) or (args.dataUC is None and args.dataC is not None)],
        ["Need at least one set of samples to compare.", args.sigUC is None, args.sigC is None, args.bkgUC is None, args.bkgC is None, args.dataUC is None, args.dataC is None],
        ["Must provide a year for processing.", args.year is None],
    )

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    Scales = ScalesCls(args.year) 
    masses = [x for x in range(15, 65, 5)]
    if args.year == "2022preEE":
        eras={
            "data": ["Run2022C","Run2022D"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_preEE" for m in masses]
        }
    elif args.year == "2022postEE":
        eras={
            "data": ["Run2022E","Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_postEE" for m in masses]
        }
    elif args.year == "2022":
        eras={
            "data": ["Run2022C","Run2022D", "Run2022E","Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D", "RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preEE", "postEE"]]
        }
    else:
        assert False, "Must specify year in Scales class!"

    # Handle path setup after checking year
    pathSetup(args.year)

    samples = {}

    if args.sigUC is not None and args.sigC is not None:
        assert os.path.exists(os.path.abspath(args.sigUC))
        assert os.path.exists(os.path.abspath(args.sigC))

        print("Loading uncorrected signal samples...")
        samples.update({"sig-uncorrected": loadPerMass(args.sigUC, eras = eras, year = args.year, era_type_arg="signal")})

        print("Loading corrected signal samples...")
        samples.update({"sig-corrected": loadPerMass(args.sigC, eras = eras, year = args.year, era_type_arg="signal")})
        
    # Load bkg/data samples
    if args.bkgUC is not None and args.bkgC is not None:
        assert os.path.exists(os.path.abspath(args.bkgUC))
        assert os.path.exists(os.path.abspath(args.bkgC))

        for key, inputs in zip(["bkg-uncorrected", "bkg-corrected"], [args.bkgUC, args.bkgC]):
            print(f"Loading {key}...")
            samples.update({key: loadSamples(inputs, eras = eras, era_type_arg="bkg")})

    if args.dataUC is not None and args.dataC is not None:
        assert os.path.exists(os.path.abspath(args.dataUC))
        assert os.path.exists(os.path.abspath(args.dataC))

        for key, inputs in zip(["data-uncorrected", "data-corrected"], [args.dataUC, args.dataC]):
            print(f"Loading {key}...")
            samples.update({key: loadSamples(inputs, eras = eras, era_type_arg="data")})

    # Setup for branches to print
    branches, boundsList, binsScaleList, names, norms, titles = loadPlottingParameters()
    branches.extend(["pho1_pt", "pho2_pt", "pho3_pt", "pho4_pt"])

    outputs = {"sigUC": dict().copy(), "sigC": dict().copy(), "bkgUC": list().copy(), "bkgC": list().copy(), "dataUC": list().copy(), "dataC": list().copy()}

    # Setup for signal comparison
    if args.sigUC is not None and args.sigC is not None:
        # Updates to include masses for smearing checks
        branches.extend(["mass_gggg", "LeadPs_mass", "SubleadPs_mass"])

        for canvas_num, sig in enumerate([samples["sig-uncorrected"], samples["sig-corrected"]]):
            for original_key in sig.keys():
                mass = original_key.replace("Signal_","")
                mass = mass[:6]

                # Add entries to signal histogram dict for mass points
                sig_type = "sigUC" if canvas_num == 0 else "sigC"
                outputs[sig_type].update({mass: list().copy()})

                # Loop through branches to print events
                for branch in branches:
                    # Store some per-event values here for sigUC and sigC!
                    outputs[sig_type][mass].append(set(sig[original_key][branch].iloc[0:5].to_list()))

        print("\n ~~~   Signal   ~~~")
        for mass in [f"{m}_GeV" for m in range(15,65,5)]:
            print(f"\n Mass: {mass.replace('_',' ')}")
            for suc, sc, branch in zip(outputs["sigUC"][mass], outputs["sigC"][mass], branches):
                # Print some per-event values here in a fancy way
                print(f"\t Branch: {branch}")
                for evt, (uc, c) in enumerate(zip(suc, sc)):
                    print(f"\t\t Event {evt}   UC: {uc:7.5f}  C: {c:7.5f}")

        branches = branches[:-4]

    # Setup for background comparison
    if args.bkgUC is not None and args.bkgC is not None:
        hypMass_bkg = np.random.choice([float(x) for x in range(15, 65, 5)], len(samples[f"bkg-corrected"]) if len(samples[f"bkg-corrected"]) > len(samples[f"bkg-uncorrected"]) else len(samples[f"bkg-uncorrected"]))
        for canvas_num, df in enumerate([samples["bkg-uncorrected"], samples["bkg-corrected"]]):
            bd_type = f"bkgUC" if canvas_num == 0 else f"bkgC"

            # Set hypMass to same random distribution to make proper comparison between samples
            df = redefine_dataframe(df, hypMass_bkg, rand=True)

            # Store per-event values here
            for branch in branches:
                outputs[bd_type].append(set(df[branch].iloc[0:5].to_list()))

        print("\n ~~~   Background   ~~~")
        for buc, bc, branch in zip(outputs["bkgUC"], outputs["bkgC"], branches):
            # Print some per-event values here in a fancy way
            print(f"\t Branch: {branch}")
            for evt, (uc, c) in enumerate(zip(buc, bc)):
                print(f"\t\t Event {evt}   UC: {uc:7.5f}  C: {c:7.5f}")

    # Setup for data comparison
    if args.dataUC is not None and args.dataC is not None:
        hypMass_data = np.random.choice([float(x) for x in range(15, 65, 5)], len(samples[f"data-corrected"]) if len(samples[f"data-corrected"]) > len(samples[f"data-uncorrected"]) else len(samples[f"data-uncorrected"]))
        for canvas_num, df in enumerate([samples["data-uncorrected"], samples["data-corrected"]]):
            bd_type = f"dataUC" if canvas_num == 0 else f"dataC"

            # Set hypMass to same random distribution to make proper comparison between samples
            df = redefine_dataframe(df, hypMass_data, rand=True)

            # Store per-event values here
            for branch in branches:
                outputs[bd_type].append(set(df[branch].iloc[0:5].to_list()))

        print("\n ~~~   Data   ~~~")
        for duc, dc, branch in zip(outputs["dataUC"], outputs["dataC"], branches):
            # Print some per-event values here in a fancy way
            print(f"\t Branch: {branch}")
            for evt, (uc, c) in enumerate(zip(duc, dc)):
                print(f"\t\t Event {evt}   UC: {uc:7.5f}  C: {c:7.5f}")


if __name__ == "__main__":
    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser", conflict_handler="resolve")

    # Input arguments
    parser.add_argument("-su", "--sigUC", help="Directory to import uncorrected signal samples from.", type=str, default=None, required=False)
    parser.add_argument("-suc", "--sigC", help="Directory to import corrected signal samples from.", type=str, default=None, required=False)
    parser.add_argument("-bu", "--bkgUC", help="Directory to import uncorrected bkg samples from.", type=str, default=None, required=False)
    parser.add_argument("-buc", "--bkgC", help="Directory to import corrected bkg samples from.", type=str, default=None, required=False)
    parser.add_argument("-du", "--dataUC", help="Directory to import uncorrected data samples from.", type=str, default=None, required=False)
    parser.add_argument("-duc", "--dataC", help="Directory to import corrected data samples from.", type=str, default=None, required=False)
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    args = parser.parse_args()

    compareCorrections(args)
