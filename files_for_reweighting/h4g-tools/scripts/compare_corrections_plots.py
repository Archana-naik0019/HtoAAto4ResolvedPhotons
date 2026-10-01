import gc
import json
import numpy as np
import copy
import argparse
import os
from ROOT import TMath, TH2D
from h4g_tools.utils.scales import Scales as ScalesCls
from h4g_tools.bdt.event_selection_bdt import redefine_dataframe
from h4g_tools.utils.runner_utils import handle_errors
from h4g_tools.utils.plotting import plotCorrectedComparison, loadPlottingParameters
from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist
from h4g_tools.utils.remove_smearing import recalc_nosmear


"""
For generating plots to compare corrections on samples:
python3 compare_corrections_plots.py --sigUC <sig-uncorrected> --sigC <sig-corrected> --bkgUC <bkg-uncorrected> --bkgC <bkg-corrected> --dataUC <data-uncorrected> --dataC <data-corrected> --year [2022, 2022preEE, 2022postEE]

Signal:
--sigUC /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/merged_outputs_signal_noCorr_noSyst_PsuedosArray_Jan30/ --sigC /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/merged_outputs_signal_MCSmear_DataEtDepScale_wSyst_PsuedosArray_Jan30/

Bkg:
--bkgUC /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/merged_outputs_bkg_noCorr_noSyst_PsuedosArray_Jan31/ --bkgC /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/merged_outputs_bkg_MCSmear_DataEtDepScale_wSyst_PsuedosArray_Jan31/

Data:
--dataUC /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/merged_outputs_data_noCorr_noSyst_PsuedosArray_Jan31/ --dataC /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/merged_outputs_data_MCSmear_DataEtDepScale_wSyst_PsuedosArray_Jan31/


Do not perform reweighting here! Want to show the difference in corrected vs uncorrected samples.
"""


scale_order = 4
scale = 10**scale_order


def compareCorrections(
    args: argparse.ArgumentParser
) -> None:

    handle_errors(
        #["Need both corrected and uncorrected signal samples.", (args.sigUC is not None and args.sigC is None) or (args.sigUC is None and args.sigC is not None)],
        ["Need both corrected and uncorrected bkg samples.", (args.bkgUC is not None and args.bkgC is None) or (args.bkgUC is None and args.bkgC is not None)],
        ["Need both corrected and uncorrected data samples.", (args.dataUC is not None and args.dataC is None) or (args.dataUC is None and args.dataC is not None)],
        ["Need at least one set of samples to compare.", args.sigUC is None, args.sigC is None, args.bkgUC is None, args.bkgC is None, args.dataUC is None, args.dataC is None],
        ["Must provide a year for processing.", args.year is None],
        ["Can only check mass_0 for signal MC!", (args.bkgUC is not None and args.bkgC is not None) or (args.dataUC is not None and args.dataC is not None), args.mass_0]
    )

    if args.sigUC is None and args.sigC is None and args.bdtJSON is not None:
        print("~~ WARNING: BDT cuts are only applied to signal MC for this test!")

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
    elif args.year == "2024":
        eras={
            "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
            "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
            "signal": [f"Signal_{m}_GeV" for m in masses]
        }
    else:
        assert False, "Must specify year in Scales class!"

    # Handle path setup after checking year
    pathSetup(args.year)

    samples = {}

    if (args.sigUC is not None and args.sigC is not None and not args.recalcNoSmearing) or (args.sigC is not None and args.recalcNoSmearing):
        if args.sigUC is not None:
            assert os.path.exists(os.path.abspath(args.sigUC))
        assert os.path.exists(os.path.abspath(args.sigC))
        
        print("Loading corrected signal samples...")
        samples.update({"sig-corrected": loadPerMass(args.sigC, eras = eras, year = args.year, era_type_arg="signal", branches=["mass_0", "lead_a_mass", "sublead_a_mass"] if args.mass_0 else None)})

        if args.recalcNoSmearing:
            print("Recalculating masses without smearing applied from corrected samples.")
            samples.update({"sig-uncorrected": {}})
            for k in samples["sig-corrected"].keys():
                print(f"Recalculating {k}...")
                samples["sig-uncorrected"].update({k: recalc_nosmear(samples["sig-corrected"][k])})
                assert len(samples["sig-uncorrected"][k]) == len(samples["sig-corrected"][k]), (len(samples["sig-uncorrected"][k]), len(samples["sig-corrected"][k]))
        else:
            print("Loading uncorrected signal samples...")
            samples.update({"sig-uncorrected": loadPerMass(args.sigUC, eras = eras, year = args.year, era_type_arg="signal", branches=["mass_0", "lead_a_mass", "sublead_a_mass"] if args.mass_0 else None)})

        if args.bdtJSON is not None:
            # Load and apply BDT cuts to signal MC samples
            # Used to check if this explains the weird signal fitting
            with open(args.bdtJSON, "r") as f:
                bdt_eff = json.load(f)
            
            print("\n~~ Applying BDT cut to corrected samples:")
            # Do this twice in case the keys are different!
            corr_keys = sorted(samples["sig-corrected"].keys())
            for k in corr_keys:
                key = k.replace("Signal_","")
                key = key[:key.find("GeV")+3]
                orig_events = len(samples["sig-corrected"][k])
                samples["sig-corrected"][k] = samples["sig-corrected"][k][samples["sig-corrected"][k].BDT_score > bdt_eff[key]["cut0"][0]]
                print(f"Reduced {k} events from {orig_events} to {len(samples['sig-corrected'][k])}")

            print("\n~~ Applying BDT cut to uncorrected samples:\n")
            uncorr_keys = sorted(samples["sig-uncorrected"].keys())
            for k in uncorr_keys:
                key = k.replace("Signal_","")
                key = key[:key.find("GeV")+3]
                orig_events = len(samples["sig-uncorrected"][k])
                samples["sig-uncorrected"][k] = samples["sig-uncorrected"][k][samples["sig-uncorrected"][k].BDT_score > bdt_eff[key]["cut0"][0]]
                print(f"Reduced {k} events from {orig_events} to {len(samples['sig-uncorrected'][k])}")
         
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

    print()

    # Setup for plot parameters
    max_order = 6
    branches, boundsList, binsScaleList, names, norms, titles = loadPlottingParameters(bdt=True)
    branches.extend(["pho1_pt", "pho2_pt", "pho3_pt", "pho4_pt"])
    boundsList.extend([[20.0, 100.0], [10.0, 80.0], [0.0, 70.0], [0.0, 70.0]])
    binsScaleList.extend([5.0/4.0, 10.0/7.0, 10.0/7.0, 10.0/7.0])
    names.extend(["pho1_pt", "pho2_pt", "pho3_pt", "pho4_pt"])
    norms.extend([False] * 4)
    titles.extend(["#gamma_{1} p_{T}", "#gamma_{2} p_{T}", "#gamma_{3} p_{T}", "#gamma_{4} p_{T}"])

    # Collect previous objects and disable circular garbage collector (ref counter still works)
    gc.collect()
    gc.disable()

    # Place to store histograms - use dicts for signal to account for mass points
    hists = {"sigUC": dict().copy(), "sigC": dict().copy(), "bkgUC": list().copy(), "bkgC": list().copy(), "dataUC": list().copy(), "dataC": list().copy()}

    if (args.sigUC is not None and args.sigC is not None and not args.recalcNoSmearing) or (args.sigC is not None and args.recalcNoSmearing):
        # Updates to include masses for smearing checks
        branches.extend(["mass_gggg", "LeadPs_mass", "SubleadPs_mass"])
        boundsList.extend([[110.0, 135.0], [0.0, 85.0], [0.0, 85.0]])
        binsScaleList.extend([4.0, 200.0/85.0, 200.0/85.0])
        names.extend(["mass_gggg", "LeadPs_mass", "SubleadPs_mass"])
        norms.extend([False, False, False])
        titles.extend(["m_{4#gamma}", "m_{a1}", "m_{a2}"])
        assert len(branches) == len(boundsList) == len(binsScaleList) == len(names) == len(norms) == len(titles), f"{len(branches)} {len(boundsList)} {len(binsScaleList)} {len(names)} {len(norms)} {len(titles)}"

        for canvas_num, sig in enumerate([samples["sig-uncorrected"], samples["sig-corrected"]]):
            for mass in sig.keys():
                mass = mass.replace("Signal_","")
                mass = mass[:6]

                # Add entries to signal histogram dict for mass points
                sig_type = "sigUC" if canvas_num == 0 else "sigC"
                hists[sig_type].update({mass: list().copy()})

                # Test case for args.mass_0 with Tanay's version of HiggsDNA
                if args.mass_0:
                    branches = ["mass_0", "lead_a_mass", "sublead_a_mass"]
                    boundsList = [[110.0, 135.0], [0.0, 85.0], [0.0, 85.0]]
                    binsScaleList = [4.0, 200.0/85.0, 200.0/85.0]
                    names = ["mass_0", "lead_a_mass", "sublead_a_mass"]
                    norms = [True] * 3
                    titles = ["mass_{4#gamma}", "mass_{a1}", "mass_{a2}"]

                # Fill histograms
                for branch, bounds, binsScale, name, norm, title in zip(branches, boundsList, binsScaleList, names, norms, titles):
                    binsScale_edit = binsScale
                    bounds_edit = bounds
                    if "LeadPs_mass" in branch or "SubleadPs_mass" in branch:
                        true_mass = float(mass[:2])
                        mass_shift = 5.0 if true_mass <= 30 else 10.0
                        bounds_edit = [true_mass-mass_shift, true_mass+mass_shift]
                        binsScale_edit = 10  # 4 = 0.25 GeV bins, 10 = 0.1 GeV bins. Used 4 for mass_gggg and 10 for mass_a1/a2 linear
                        # Rebin to 0.2 GeV bins after plotting linear histograms (i.e., use Rebin(2))
                    if (args.recalcNoSmearing or args.massOnly) and branch not in ["mass_gggg", "LeadPs_mass", "SubleadPs_mass"]:
                        continue
                    hists[sig_type][mass].append(fillHist(branch, sig[f"Signal_{mass}" + ("_{args.year[4:]}" if "2022" in args.year else "")], bounds_edit, binsScale=binsScale_edit, name=f"{name}_{sig_type}", normalize=norm, preventOverFlow=True))
                    hists[sig_type][mass][-1].Scale(scale / hists[sig_type][mass][-1].Integral())

        # Plotting histograms for signal comparisons
        hists_sig = {}
        for mass in hists["sigUC"].keys():
            paths_sig = [
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/pho1_mvaID_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/pho2_mvaID_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/pho3_mvaID_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/pho4_mvaID_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/LeadPs_pt_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/SubleadPs_pt_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/dR_aa_mass_gggg_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/LeadPs_interMass_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/SubleadPs_interMass_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/Ps_massDiff_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/cos_ag_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/pho1_pT_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/pho2_pT_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/pho3_pT_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/pho4_pT_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/mass_gggg_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/mass_a1_correction_compare.png",
                f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/mass_a2_correction_compare.png",
            ]

            if args.mass_0:
                paths_sig = [
                    f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/mass_0_correction_compare.png",
                    f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/mass_a1_correction_compare.png",
                    f"./plots/BDT/{args.year}/compare_corrected/{mass.replace('Signal_','')[:6]}/mass_a2_correction_compare.png",
                ]

            if args.recalcNoSmearing or args.massOnly:
                paths_sig = paths_sig[-3:]
                titles = titles[-3:]
                branches = branches[-3:]

            for path, suc, sc, title, branch in zip(paths_sig, hists["sigUC"][mass], hists["sigC"][mass], titles, branches):
                print()
                hists_sig = {
                    f"sigUC_{mass}": suc,
                    f"sigC_{mass}": sc
                }

                if "pT" in path:
                    for pho, mx in zip(["pho1", "pho2", "pho3", "pho4"], [1000, 1400, 1800, 3000]):
                        if pho in path:
                            plotCorrectedComparison(hists_sig, title, path, Scales.sumw[args.year[:4]][mass.replace('Signal_','')[:6]], log=False, minimum=-0.9, maximum=mx/1.5, canvas_num=str(canvas_num)+branch, year=args.year, Scales=Scales)
                        else:
                            pass
                elif "pT" not in path:
                    # Also plot linear y-scale for mass plot
                    if "mass_gggg" in path or "mass_a1" in path or "mass_a2" in path and "dR_aa" not in path:
                        # Linear plot
                        print(f"Original Linear Binning: sigUC {hists_sig[list(hists_sig.keys())[0]].GetNbinsX()}  sigC {hists_sig[list(hists_sig.keys())[1]].GetNbinsX()}")
                        print(f"Original Linear Integral: sigUC {hists_sig[list(hists_sig.keys())[0]].Integral()}  sigC {hists_sig[list(hists_sig.keys())[1]].Integral()}")
                        print(f"Original Linear Entries: sigUC {hists_sig[list(hists_sig.keys())[0]].GetEntries()}  sigC {hists_sig[list(hists_sig.keys())[1]].GetEntries()}")
                        print(f"Original Linear Low Edge: sigUC {hists_sig[list(hists_sig.keys())[0]].GetBinLowEdge(1)}  sigC {hists_sig[list(hists_sig.keys())[1]].GetBinLowEdge(1)}")
                        plotCorrectedComparison(hists_sig, title, path.replace(".png", "_linear.png"), log=False if not args.mass_0 else False, minimum=0.0, canvas_num=str(canvas_num)+branch, year=args.year, Scales=Scales, sim=True, legend_text={"sigUC": "All corrections without smearing", "sigC": "All corrections with smearing"})

                        # Also make 2D plots
                        plotCorrectedComparison(hists_sig, title, path, log=True if not args.mass_0 else False, minimum=0.0 if args.mass_0 else 10**-1 + 10**0, maximum=10**5 - 5*10**(5-1) if not args.mass_0 else 0.2, canvas_num=str(canvas_num)+branch, year=args.year, Scales=Scales, sim=True, legend_text={"sigUC": "All corrections without smearing", "sigC": "All corrections with smearing"}, df_corr=samples["sig-corrected"][f"Signal_{mass}"], df_uncorr=samples["sig-uncorrected"][f"Signal_{mass}"], branch=branch, unsmeared = args.recalcNoSmearing)

                    # Plot log version after rebinning
                    if "mass_a1" in path or "mass_a2" in path:
                        for hk,h in hists_sig.items():
                            h.Rebin(2)  # From 0.1 GeV bins to 0.2 GeV bins
                    plotCorrectedComparison(hists_sig, title, path, log=True if not args.mass_0 else False, minimum=0.0 if args.mass_0 else 10**-1 + 10**0, maximum=10**5 - 5*10**(5-1) if not args.mass_0 else 0.2, canvas_num=str(canvas_num)+branch, year=args.year, Scales=Scales, sim=True, legend_text={"sigUC": "All corrections without smearing", "sigC": "All corrections with smearing"}, branch=branch)

        # Remove masses from plotting parameters after validating signal (DO NOT LOOK AT MASS IN DATA OR BKG!!!)
        if not args.recalcNoSmearing and not args.massOnly:
            branches, boundsList, binsScaleList, names, norms, titles = branches[:-3], boundsList[:-3], binsScaleList[:-3], names[:-3], norms[:-3], titles[:-3]
            assert len(branches) == len(boundsList) == len(binsScaleList) == len(names) == len(norms) == len(titles), f"{len(branches)} {len(boundsList)} {len(binsScaleList)} {len(names)} {len(norms)} {len(titles)}"

    # Plotting histograms for background/data comparisons
    if args.bkgUC is not None and args.bkgC is not None:
        hypMass_bkg = np.random.choice([float(x) for x in range(15, 65, 5)], len(samples[f"bkg-corrected"]) if len(samples[f"bkg-corrected"]) > len(samples[f"bkg-uncorrected"]) else len(samples[f"bkg-uncorrected"]))
        for canvas_num, df in enumerate([samples["bkg-uncorrected"], samples["bkg-corrected"]]):
            bd_type = f"bkgUC" if canvas_num == 0 else f"bkgC"

            # Set hypMass to same random distribution to make proper comparison between samples
            if not args.mass_0:
                df = redefine_dataframe(df, hypMass_bkg, rand=True)

            # Fill histograms
            for branch, bounds, binsScale, name, norm, title in zip(branches, boundsList, binsScaleList, names, norms, titles):
                binsScale_edit = binsScale
                bounds_edit = bounds
                if "LeadPs_mass" in branch or "SubleadPs_mass" in branch:
                    bounds_edit = [15.0, 65.0]
                    binsScale_edit = 4.0  # 0.25 GeV bins
                hists[bd_type].append(fillHist(branch, df, bounds_edit, binsScale=binsScale_edit, name=f"{name}_{bd_type}", normalize=norm))

        # Setup for bkg comparison plots
        paths_bkg = [
            f"./plots/BDT/{args.year}/compare_corrected/bkg/pho1_mvaID_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/pho2_mvaID_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/pho3_mvaID_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/pho4_mvaID_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/LeadPs_pt_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/SubleadPs_pt_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/dR_aa_mass_gggg_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/LeadPs_interMass_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/SubleadPs_interMass_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/Ps_massDiff_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/cos_ag_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/pho1_pT_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/pho2_pT_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/pho3_pT_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/bkg/pho4_pT_correction_compare.png",
        ]

        hists_bkg = {}
        for path, buc, bc, title in zip(paths_bkg, hists["bkgUC"], hists["bkgC"], titles):
            hists_bkg = {
                "bkgUC": buc,
                "bkgC": bc
            }

            if "pT" in path:
                plotCorrectedComparison(hists_bkg, title, path, log=False, minimum=-0.9, maximum=10000, canvas_num=str(canvas_num)+branch, year=args.year, Scales=Scales, fix_res_scale=True)
            elif "pT" not in path:
                # Also plot linear y-scale for mass plot
                if "mass_gggg" in path:
                    plotCorrectedComparison(hists_bkg, title, path.replace(".png", "_linear.png"), log=False, canvas_num=str(canvas_num)+branch, year=args.year, Scales=Scales, fix_res_scale=True)
                plotCorrectedComparison(hists_bkg, title, path, log=True, maximum=10**max_order - 5*10**(max_order-1), canvas_num=str(canvas_num)+branch, year=args.year, Scales=Scales, fix_res_scale=True)

    # Setup for data comparison plots
    if args.dataUC is not None and args.dataC is not None:
        hypMass_data = np.random.choice([float(x) for x in range(15, 65, 5)], len(samples[f"data-corrected"]) if len(samples[f"data-corrected"]) > len(samples[f"data-uncorrected"]) else len(samples[f"data-uncorrected"]))
        for canvas_num, df in enumerate([samples["data-uncorrected"], samples["data-corrected"]]):
            bd_type = f"dataUC" if canvas_num == 0 else f"dataC"

            # Set hypMass to same random distribution to make proper comparison between samples
            if not args.mass_0:
                df = redefine_dataframe(df, hypMass_data, rand=True)

            # Fill histograms
            for branch, bounds, binsScale, name, norm, title in zip(branches, boundsList, binsScaleList, names, norms, titles):
                hists[bd_type].append(fillHist(branch, df, bounds, binsScale=binsScale, name=f"{name}_{bd_type}", normalize=norm))

        paths_data = [
            f"./plots/BDT/{args.year}/compare_corrected/data/pho1_mvaID_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/pho2_mvaID_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/pho3_mvaID_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/pho4_mvaID_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/LeadPs_pt_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/SubleadPs_pt_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/dR_aa_mass_gggg_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/LeadPs_interMass_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/SubleadPs_interMass_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/Ps_massDiff_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/cos_ag_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/pho1_pT_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/pho2_pT_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/pho3_pT_correction_compare.png",
            f"./plots/BDT/{args.year}/compare_corrected/data/pho4_pT_correction_compare.png",
        ]

        hists_data = {}
        for path, duc, dc, title in zip(paths_data, hists["dataUC"], hists["dataC"], titles):
            hists_data = {
                "dataUC": duc,
                "dataC": dc
            }

            if "pT" in path:
                plotCorrectedComparison(hists_data, title, path, log=False, minimum=-0.9, maximum=10000, canvas_num=str(canvas_num)+branch, year=args.year, Scales=Scales, fix_res_scale=True)
            elif "pT" not in path:
                # Also plot linear y-scale for mass plot
                if "mass_gggg" in path:
                    plotCorrectedComparison(hists_data, title, path.replace(".png", "_linear.png"), log=False, canvas_num=str(canvas_num)+branch, year=args.year, Scales=Scales, fix_res_scale=True)
                plotCorrectedComparison(hists_data, title, path, log=True, maximum=10**max_order - 5*10**(max_order-1), canvas_num=str(canvas_num)+branch, year=args.year, Scales=Scales, fix_res_scale=True)
    

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
    parser.add_argument("-m0", "--mass_0", action="store_true", required=False, help="Load mass_0 from Tanay's samples.")
    parser.add_argument("-rc", "--recalcNoSmearing", action="store_true", required=False, help="Recalculate masses without smearing and only generate those plots.")
    parser.add_argument("-bdt", "--bdtJSON", type=str, required=False, help="BDT cuts json to apply. Used for signal ONLY!")
    parser.add_argument("-m", "--massOnly", action="store_true", required=False, help="Plot masses only.")
    args = parser.parse_args()

    compareCorrections(args)
