import argparse
import pandas as pd
import numpy as np
import xgboost
import os
import time
import copy
import random
import matplotlib
import matplotlib.pyplot as plt
from typing import Dict, List
from h4g_tools.utils.runner_utils import handle_errors
from h4g_tools.bdt.event_selection_bdt import define_variables, redefine_dataframe, reduce_dataframe, train_ES_bdt, predict_BDT
from h4g_tools.bdt.reweight import reweight1Dim, loadNDimWeights
from h4g_tools.utils.loading import fillHist, loadSamples
from h4g_tools.utils.plotting import plot_all, plot_data_sig_bkg_comparison, plot_even_odd_BDT, plot_variations
from h4g_tools.utils.saving import dump_pandas
from h4g_tools.utils.scales import all4EB, SideBand, SignalRegion
from h4g_tools.utils.save_plots_to_root import savePlots



def gen_BDT(
    year: str,
    sig: str,
    scale: float,
    debug: bool,
    gen: bool,
    train_no_params: bool,
    xgb_name: str,
    eras: list = None,
    Scales: str = None,
    noRW: bool = False,
    mva_cut: float = -1.0,
    EB_only: bool = False,
    bad: bool = False,
    training_json: str = None
) -> None:
    """
    Trains BDT over various masses for pseudoscalars.
    """

    args = argparse.Namespace(
        train = True,
        model = f"../models/{xgb_name}.xgb",
        sigInput = sig,
        bkgInput = None,
        data = None,
        debug = debug,
        save = True,
        gen = ("BDT_input_params.json" if training_json is None else training_json) if gen else None,
        train_no_params = train_no_params,
        year = year,
        mva_cut = mva_cut,
        EB_only = EB_only,
        bad = bad
    )

    # Reweighting loading that contains all bkg events + weights
    print("\nRunning N-dim reweighting procedure for BDT inputs...")
    bkg_weights = loadNDimWeights(year)
    if noRW:
        print(f"No reweight requested! Replacing weights with 1.", end=" ")
        bkg_weights["Ndimreweight"] = pd.Series([1.0] * len(bkg_weights))
        print(f"Sum: {sum(bkg_weights['Ndimreweight'])} Len: {len(bkg_weights)}")

    # BDT training
    print(f"\nRunning BDT setup...")
    run_bdt(args, background=bkg_weights, eras = eras, Scales=Scales)
    print("\n")


def run_samples_BDT(
    masses: List[int],
    year: int,
    d: str,
    xgb_name: str,
    sig: bool = False,
    eras: list = None,
    Scales: str = None,
    nominalOnly: bool = True,
    pick_variation: str = None,
    is_data: bool = False,
    sig_mix: bool = False
) -> None:
    """
    Runs BDT over signal, data, or event mixed background.
    """

    # Need to load full dataset and reduced dataset. Process with reduced and concat the BDT score to the full dataset.
    samples = {}
    if not sig:
        samples.update({"nominal": {}})
        print("Loading full dataset:")
        samples["nominal"].update({"data": define_variables(d, reduce_df=False, perMass=False, eras=eras)})

        print("\nLoading reduced dataset:")
        # Need to reduce dataframe later!!! Keeping it unreduced for recalculating mHyp variables.
        samples["nominal"].update({"data_reduced": define_variables(d, reduce_df=True, perMass=False, eras=eras)})
        print()

    elif sig:
        # Need to sort through all variations applied to the signal MC
        if nominalOnly:
            variations = ["nominal"]
        else:
            variations = os.listdir(os.path.join(*(d, os.listdir(d)[0])))
            if pick_variation is not None:
                variations = [v for v in variations if v == pick_variation]

        for variation in variations:
            print(f"~~~ Loading variation: {variation} ~~~")
            samples.update({variation: {}})
            print("Loading full dataset:")
            samples[variation].update({"data": define_variables(d, reduce_df=False, perMass=True, eras=eras, variation=variation if not sig_mix else None)})

            to_remove = []
            for mass in list(samples[variation]["data"].keys()).copy():
                samples[variation]["data"][int(mass[mass.find("GeV")-3 : mass.find("GeV")-1])] = samples[variation]["data"][mass]
                to_remove.append(mass)
            for mass in to_remove:
                del samples[variation]["data"][mass]

            print("\nLoading reduced dataset:")
            samples[variation].update({"data_reduced": define_variables(d, reduce_df=True, perMass=True, eras=eras, variation=variation, applyReduce=True)})
            print()

    # Generic model processing
    args = argparse.Namespace(
        train = False,
        model = f"../models/{xgb_name}.xgb",
        sigInput = None,
        bkgInput = None,
        data = d,
        hypMass = None,
        save = True,
        year = year,
        generic = "generic",
    )

    if not sig:
        print("~~~ Running BDT over samples for flat mHyp distribution ~~~")
        samples["nominal"]["data_reduced"] = reduce_dataframe(samples["nominal"]["data_reduced"])
        # Removing existing BDT scores (i.e. always recalculate)
        samples["nominal"]["data"].drop("BDT_score", axis=1, inplace=True)
        processSamples(args, samples["nominal"], "nominal", Scales=Scales)
        print()

    # Specific mass processing
    for mass in masses:
        args = argparse.Namespace(
            train = False,
            model = f"{xgb_name}.xgb",
            sigInput = None,
            bkgInput = None,
            data = d,
            hypMass = mass,
            save = True,
            year = year,
            generic = None,
        )

        print(f"~~~ Running BDT over samples for {mass} GeV mass point ~~~")
        for variation in samples.keys():
            print(f"\nProcessing variation: {variation}")

            # Handles different m_Hyp values for background and data. Variations are just nominal
            if not sig:
                samples[variation]["data"] = redefine_dataframe(samples[variation]["data"], mass)
                samples[variation]["data_reduced"] = reduce_dataframe(redefine_dataframe(samples[variation]["data"], mass))
                #print(mass)
                #print(samples[variation]["data"].LeadPs_interMass)

            # Removing existing BDT scores (i.e. always recalculate)
            samples[variation]["data"][mass].drop("BDT_score", axis=1, inplace=True, errors="ignores") if sig else samples[variation]["data"].drop("BDT_score", axis=1, inplace=True, errors="ignore")
            processSamples(args, samples[variation], variation, mass = None if not sig else mass, Scales=Scales)
        print("\n")


def plot_results(
    masses: List[int],
    year: str,
    xgb_name: str,
    sig: str,
    bkg: str,
    d: str,
    Scales: str = None,
    extra_ratio: bool = False,
    EB_only: bool = False,
    nominalOnly: bool = True,
    original: bool = False,
    blind: bool = True
) -> None:
    """
    Generates plots of BDT score distribution for a given pseudoscalar mass, event mixing (background), and data.
    """

    print()
    # Not generic - i.e. mass specific
    for mass in masses:
        args = argparse.Namespace(
            hypMass = float(mass),
            sigInput = f"{sig.rstrip('/')}/{year}/{xgb_name}/signal/{mass}_GeV/",
            bkgInput = f"{bkg.rstrip('/')}/{year}/{xgb_name}/background/{mass}_GeV/",
            data = f"{d.rstrip('/')}/{year}/{xgb_name}/data/{mass}_GeV/",
            title = "BDT Score",
            year = year,
            model = xgb_name,
            load_1dimreweight = True,
            blind = blind
        )

        print(f"Generating plot for {mass} GeV pseudoscalar:")
        generateESPlots(args, Scales=Scales, extra_ratio=extra_ratio, EB_only=EB_only, original=original)
        print("\n")

    # Generic - i.e. flat mHyp distribution
    args = argparse.Namespace(
        hypMass = 0.0,
        sigInput = f"{sig.rstrip('/')}/{year}/{xgb_name}/signal/30_GeV/",
        bkgInput = f"{bkg.rstrip('/')}/{year}/{xgb_name}/generic/background/",
        data = f"{d.rstrip('/')}/{year}/{xgb_name}/generic/data/",
        title = "BDT Score",
        year = year,
        model = xgb_name,
        load_1dimreweight = True,
        blind = blind
    )

    print("Generating plot for generic BDT model (unmodified mHyp):")
    generateESPlots(args, Scales=Scales, generic=True, extra_ratio=extra_ratio, EB_only=EB_only, original=original)
    print("\n")

    if not nominalOnly:
        print("Generating plots for systematic variations:")
        plot_stats(args.model, args.year, f"{sig.rstrip('/')}/{year}/{xgb_name}/signal/", Scales=Scales, original=original)


def run_bdt(
    args: argparse.Namespace,
    background: Dict = None,
    data: Dict = None,
    data_reduced: Dict = None,
    bkg_weights: pd.DataFrame = None,
    eras: list = None,
    Scales: str = None,
) -> None:
    """
    Run BDT for a variety of inputs.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Check for parser issues
    handle_errors(
        ["Need to specify samples in order to make predictions.", not args.train, args.data is None, data is None, data_reduced is None],
        ["Need to provide a model to load from file.", not args.train, args.model is None],
        ["Cannot specify training samples when not training.", not args.train, args.sigInput is not None or args.bkgInput is not None],
        ["Need to provide a signal samples to train with.", args.train, args.sigInput is None],
        ["Need to provide a background samples to train with.", args.train, args.bkgInput is None, background is None],
        ["Need to specify path/filename to save model.", args.train, args.model is None, args.save],
        ["Cannot specify model to load while training.", args.train, args.model is not None, not args.save],
    )

    # Load samples
    print("Loading provided samples...")
    samples = {}
    for idx, samplePath in enumerate([args.sigInput, args.bkgInput]):
        if samplePath is not None:
            if samplePath is args.sigInput:
                print("Loading signal:")
                key = "signal"
                extra = ["weight", "pho1_isScEtaEB", "pho2_isScEtaEB", "pho3_isScEtaEB", "pho4_isScEtaEB"]
            elif samplePath is args.bkgInput and background is None:
                print("Loading background:")
                key = "background"
                extra = None
            else:
                assert False

            samples.update({key:
                define_variables(
                    samplePath,
                    reduce_df = True,
                    extra = extra,
                    eras = eras,
                    mass = True if args.train else False
                )
            })

        elif samplePath is None and idx == 1 and background is not None:
            key = "background"
            samples.update({key: reduce_dataframe(background, extra = ["mass_gggg", "Ndimreweight", "pho1_isScEtaEB", "pho2_isScEtaEB", "pho3_isScEtaEB", "pho4_isScEtaEB"])})
            print("Reformatted background.")

    if args.EB_only:
        print("Training on all 4 photons in EB.")
        for key, sample in samples.items():
            samples[key] = all4EB(sample)

    # Train BDT from dataset
    if args.sigInput is not None and (args.bkgInput is not None or background is not None) and args.train:
        if args.mva_cut > -1.0:
            print(f"Applying photon MVA cut > {args.mva_cut} to all samples!")
            print(f"Before cut - Signal: {len(samples['signal'])}  Bkg: {len(samples['background'])}")
            samples["signal"] = samples["signal"][(samples["signal"]["pho3_mvaID"] > args.mva_cut) & (samples["signal"]["pho4_mvaID"] > args.mva_cut)]
            samples["background"] = samples["background"][(samples["background"]["pho3_mvaID"] > args.mva_cut) & (samples["background"]["pho4_mvaID"] > args.mva_cut)]
            print(f"After cut - Signal: {len(samples['signal'])}  Bkg: {len(samples['background'])}")

        assert (args.model is None and not args.save) or (args.model is not None and args.save)

        """"
        ### TEST ###
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("sig:", len(samples["signal"]), "bkg:", len(samples["background"]))
        samples["signal"] = samples["signal"][samples["signal"]["pho3_mvaID"] > -0.93]
        samples["signal"] = samples["signal"][samples["signal"]["pho4_mvaID"] > -0.93]
        samples["background"] = samples["background"][samples["background"]["pho3_mvaID"] > -0.93]
        samples["background"] = samples["background"][samples["background"]["pho4_mvaID"] > -0.93]
        print("sig:", len(samples["signal"]), "bkg:", len(samples["background"]))
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        ### TEST ###
        """

        # Keep m_4gamma for correlation check
        mass_gggg = {"signal": samples["signal"].mass_gggg, "background": samples["background"].mass_gggg}
        samples["signal"].drop(columns=["mass_gggg"], inplace=True)
        samples["background"].drop(columns=["mass_gggg"], inplace=True)

        bdt = train_ES_bdt(
            samples["signal"],
            samples["background"],
            args.gen,
            args.year,
            modelPath = args.model,
            save = args.save,
            debug = args.debug,
            train_no_params = args.train_no_params,
            Scales = Scales,
            mass_gggg = mass_gggg
        )


def processSamples(
    args: argparse.Namespace,
    samples: dict = None,
    variation: str = None,
    mass: float = None,
    Scales: str = None
) -> None:

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Set required directory for model loading/saving
    if args.model is not None and f"{cwd}/../../models/" not in args.model:
        args.model = f"{cwd}/../../models/{args.model}"
        assert os.path.exists(args.model), f"Requested model {args.model} does not exist."

    # Check for parser issues
    handle_errors(
        ["Need to provide a model to load from file.", args.model is None],
        ["Need to provide samples to process.", samples is None],
    )
    
    # Load BDT from file
    bdt = xgboost.XGBClassifier()
    bdt.load_model(args.model)
    print(f"Loaded model: {args.model.split('/')[-1]}")

    # Run BDT over provided samples
    print(f"Number of events: {len(samples['data_reduced'] if mass is None else samples['data_reduced'][mass])}")

    start_time = time.time()
    results = {"data": predict_BDT(bdt, samples["data_reduced"] if mass is None else samples["data_reduced"][mass])}

    # Validation that the new scores are applied done May 5, 2026. Uncomment below to confirm that it works
    #print("NEW PREDICITONS")
    #print(results["data"])
    results["data"] = pd.concat([samples["data"].reset_index() if mass is None else samples["data"][mass].reset_index(), pd.Series(results["data"], name="BDT_score")], axis=1)
    results["data"] = results["data"].loc[:, ~results["data"].columns.duplicated()].copy()
    #print("CONCATTED SCORES")
    #print(results["data"].BDT_score)
    assert type(results["data"]) is pd.DataFrame, "Prediction set is not a DataFrame."

    # Apply transformation to BDT scores
    results["data"] = Scales.TransformBDT(results["data"], "BDT_score")

    elapsed_time = time.time() - start_time
    print(f"Elapsed time to run BDT: {elapsed_time:6.4}s")

    # Save outputs
    if args.save:
        # Note that the output is much smaller than the total input directory size due to a larger overhead
        # with parquet files. To run a check, simply dump the merged input to a parqet file and compare the
        # sizes of the file with merged inputs and the file with result and BDT score.
        output = f"{args.data.strip('/').split('/')[-1]}"
        assert os.path.exists(f"{cwd}/outputs/BDT/{args.year}"), f"{cwd}/outputs/BDT/{args.year}"

        total = len(results["data"])
        if "signal" in output.lower():
            output = "outputs_signal_BDT.parquet"
            sub = "signal"
        elif "data" in output and ("EM" not in output or "bkg" not in output or "background" not in output):
            output = "outputs_data_BDT.parquet"
            sub = "data"
        else:
            output = "outputs_background_BDT.parquet"
            sub = "background"

        path = ""
        if args.generic == None:
            path = f"{cwd}/outputs/BDT/{args.year}/{args.model.split('/')[-1].replace('.xgb','')}/{sub}/{int(args.hypMass)}_GeV/"
        elif args.generic == "generic":
            path = f"{cwd}/outputs/BDT/{args.year}/{args.model.split('/')[-1].replace('.xgb','')}/{args.generic}/{sub}/"
        
        if not os.path.exists(f"{path}/{variation}"):
            os.mkdir(f"{path}/{variation}")
        assert os.path.exists(f"{path}/{variation}"), f"{path}/{variation}"

        assert len(results["data"]) > 0, len(results["data"])
        for t in range(3):
            dump_pandas(results["data"], output, os.path.join(*(path, variation)), pyarrow = True if t > 0 else False)
            if os.stat(os.path.join(*(path, variation, output))).st_size != 0:
                break
            print(f"File size is zero! Retrying after removing file and directory... [{t+1}/3]")
            if os.path.exists(os.path.join(*(path, variation, output))):
                os.remove(os.path.join(*(path, variation, output)))
                os.rmdir(os.path.join(*(path, variation)))
            if t == 0:
                print(results["data"])
                print("~~~ WARNING: Check storage space! ~~~") 

        print(f"Saved {len(results['data'])} events to output file: {os.path.join(*(path, variation, output))}")
        assert len(os.listdir(os.path.join(*(path, variation)))) == 1, os.listdir(os.path.join(*(path, variation)))
        assert os.stat(os.path.join(*(path, variation, output))).st_size != 0, os.stat(os.path.join(*(path, variation, output))).st_size
            

def generateESPlots(
    args: argparse.Namespace,
    Scales: str = None,
    generic: bool = None,
    extra_ratio: bool = False,
    EB_only: bool = False,
    original: bool = False,
) -> None:
    """
    Setup function for running plot_full_ES_BDT.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Check for parser issues
    handle_errors(
        ["Need to specify signal samples.", args.sigInput is None],
        ["Need to specify background samples.", args.bkgInput is None],
        ["Need to specify data samples.", args.data is None],
    )

    # Load samples
    bdt_score = "BDT_score_orig" if original else "BDT_score"
    print("Loading provided samples...")
    samples = {}
    for samplePath in [args.bkgInput, args.data, args.sigInput]:
        # NOTE: This automatically only looks at nominal variations, as intended!
        assert len(os.listdir(samplePath)) == 1 or "nominal" in os.listdir(samplePath), os.listdir(samplePath)
        if samplePath is args.sigInput:
            key = "signal"
        elif samplePath is args.bkgInput:
            key = "background"
        elif samplePath is args.data:
            key = "data"

        if key in ["background", "data"]:
            samples.update({key:
                (loadSamples(samplePath, [bdt_score, "mass_gggg", "pho1_isScEtaEB", "pho2_isScEtaEB", "pho3_isScEtaEB", "pho4_isScEtaEB"])) if EB_only else loadSamples(samplePath, [bdt_score, "mass_gggg"])
            })
        elif key in ["signal"]:
            samples.update({key if not generic else "Signal_30_GeV":
                (loadSamples(samplePath, [bdt_score, "mass_gggg", "weight", "genWeight", "pho1_isScEtaEB", "pho2_isScEtaEB", "pho3_isScEtaEB", "pho4_isScEtaEB"])) if EB_only else loadSamples(samplePath, [bdt_score, "mass_gggg", "weight", "genWeight"])
            })
    
    assert len(samples) == 3 or len(samples) == 12, len(samples)

    # Fill hists
    hists = {}
    bounds = [0.0, 1.0]
    bins = 30.0  # Consistent with binning from BDT training plots
    #bins = 60.0
    for key, sample in samples.items():
        if key not in ["background", "data"]:
            hists.update({
                f"{key}": fillHist(
                    bdt_score,
                    all4EB(sample) if EB_only else sample,
                    bounds = bounds,
                    binsScale = bins / (bounds[1] - bounds[0]),
                    normalize = False,
                    name = f"BDT_score_{key}_{int(args.hypMass)}_GeV" if not generic else f"BDT_score_{key}_30_GeV"
                )
            })

            hists.update({
                f"{key}_SB": fillHist(
                    bdt_score,
                    all4EB(sample) if EB_only else sample,
                    bounds = bounds,
                    binsScale = bins / (bounds[1] - bounds[0]),
                    normalize = False,
                    sb = True,
                    name = f"BDT_score_{key}_SB_{int(args.hypMass)}_GeV" if not generic else f"BDT_score_SB_{key}_30_GeV"
                )
            })

            hists.update({
                f"{key}_SR": fillHist(
                    bdt_score,
                    all4EB(sample[(sample.mass_gggg >= 115.0) & (sample.mass_gggg <= 135.0)]) if EB_only else sample[(sample.mass_gggg >= 115.0) & (sample.mass_gggg <= 135.0)],
                    bounds = bounds,
                    binsScale = bins / (bounds[1] - bounds[0]),
                    normalize = False,
                    name = f"BDT_score_{key}_SR_{int(args.hypMass)}_GeV" if not generic else f"BDT_score_SR_{key}_30_GeV"
                )
            })

    # Reweight bkg to data
    samples["background_rw"] = reweight1Dim(samples["background"], samples["data"], col=bdt_score)

    # Generates hists for pre-reweighting
    hists["background_SB"] = fillHist(
        bdt_score,
        all4EB(samples["background"]) if EB_only else samples["background"],
        bounds = bounds,
        binsScale = bins / (bounds[1] - bounds[0]),
        normalize = False,
        sb = True,
        name = f"BDT_score_background_sideband_{int(args.hypMass)}_GeV"
    )

    hists["background_SR"] = fillHist(
        bdt_score,
        all4EB(SignalRegion(samples["background"])) if EB_only else SignalRegion(samples["background"]),
        bounds = bounds,
        binsScale = bins / (bounds[1] - bounds[0]),
        normalize = False,
        name = f"BDT_score_background_signal_region_{int(args.hypMass)}_GeV"
    )

    hists["background_all"] = fillHist(
        bdt_score,
        all4EB(samples["background"]) if EB_only else samples["background"],
        bounds = bounds,
        binsScale = bins / (bounds[1] - bounds[0]),
        normalize = False,
        name = f"BDT_score_background_{int(args.hypMass)}_GeV"
    )

    hists["data_SB"] = fillHist(
        bdt_score,
        all4EB(samples["data"]) if EB_only else samples["data"],
        bounds = bounds,
        binsScale = bins / (bounds[1] - bounds[0]),
        normalize = False,
        sb = True,
        name = f"BDT_score_data_sideband_{int(args.hypMass)}_GeV"
    )

    hists["data_SR"] = fillHist(
        bdt_score,
        all4EB(SignalRegion(samples["data"])) if EB_only else SignalRegion(samples["data"]),
        bounds = bounds,
        binsScale = bins / (bounds[1] - bounds[0]),
        normalize = False,
        name = f"BDT_score_data_signal_region_{int(args.hypMass)}_GeV"
    )

    hists["data_all"] = fillHist(
        bdt_score,
        all4EB(samples["data"]) if EB_only else samples["data"],
        bounds = bounds,
        binsScale = bins / (bounds[1] - bounds[0]),
        normalize = False,
        name = f"BDT_score_data_{int(args.hypMass)}_GeV"
    )

    # Generates hists for post-reweighting
    hists["background_reweight_SB"] = fillHist(
        bdt_score,
        all4EB(samples["background_rw"]) if EB_only else samples["background_rw"],
        bounds = bounds,
        binsScale = bins / (bounds[1] - bounds[0]),
        normalize = False,
        name = f"BDT_score_background_sideband_{int(args.hypMass)}_GeV",
        sb = True,
        customWeights="1dimreweight"
    )

    hists["background_reweight_SR"] = fillHist(
        bdt_score,
        samples["background_rw"][(samples["background_rw"].mass_gggg >= 115.0) & (samples["background_rw"].mass_gggg <= 135.0)],
        bounds = bounds,
        binsScale = bins / (bounds[1] - bounds[0]),
        normalize = False,
        name = f"BDT_score_background_signal_region_{int(args.hypMass)}_GeV",
        customWeights="1dimreweight"
    )

    hists["background_reweight_all"] = fillHist(
        bdt_score,
        samples["background_rw"],
        bounds = bounds,
        binsScale = bins / (bounds[1] - bounds[0]),
        normalize = False,
        name = f"BDT_score_background_{int(args.hypMass)}_GeV",
        customWeights="1dimreweight"
    )

    print(f"Chi2/ndf for {args.hypMass}: {hists['data_SB'].Chi2Test(hists['background_SB'], 'UU CHI2/NDF')}")

    # Reorder hists dictionary to play nice with the color lists
    # Data/bkg are from sidebands
    hists_reordered_SB = {
        "background_reweight": hists["background_reweight_SB"],
        "background": hists["background_SB"],
        "data": hists["data_SB"],
        "signal": hists["signal" if not generic else "Signal_30_GeV"]  # Show signal from full region
    }

    hists_reordered_all_dataBlind = {
        "background_reweight": hists["background_reweight_all"],
        "background": hists["background_all"],
        "data": hists["data_SB"],
        "signal": hists["signal" if not generic else "Signal_30_GeV"]
    }

    #"""
    hists_reordered_all = {
        "background_reweight": hists["background_reweight_all"],
        "background": hists["background_all"],
        "data": hists["data_all"],
        "signal": hists["signal" if not generic else "Signal_30_GeV"]
    }
    #"""

    # Scale signal histograms to reweighted bkg integral
    #hists_reordered_SB["signal"].Scale(hists_reordered_SB["data"].Integral() / hists_reordered_SB["signal"].Integral())
    #hists_reordered_SR["signal"].Scale(hists_reordered_SR["data"].Integral() / hists_reordered_SR["signal"].Integral())
    #hists_reordered_all["signal"].Scale(hists_reordered_all["data"].Integral() / hists_reordered_all["signal"].Integral())
        
    # Generate plots
    if not generic:
        path_SB = f"{cwd}/../../scripts/plots/BDT/{args.year}/{int(args.hypMass)}_GeV/{args.model}/BDT_SB_{int(args.hypMass)}_GeV.pdf"
        path_SR = f"{cwd}/../../scripts/plots/BDT/{args.year}/{int(args.hypMass)}_GeV/{args.model}/BDT_SR_{int(args.hypMass)}_GeV.pdf"
        path_all = f"{cwd}/../../scripts/plots/BDT/{args.year}/{int(args.hypMass)}_GeV/{args.model}/BDT_{int(args.hypMass)}_GeV.pdf"
    else:
        path_SB = f"{cwd}/../../scripts/plots/BDT/{args.year}/BDT_pred_generic/{args.model}/BDT_SB_pred.pdf"
        path_SR = f"{cwd}/../../scripts/plots/BDT/{args.year}/BDT_pred_generic/{args.model}/BDT_SR_pred.pdf"
        path_all = f"{cwd}/../../scripts/plots/BDT/{args.year}/BDT_pred_generic/{args.model}/BDT_pred.pdf"

    if original:
        path_SB = path_SB.replace(".pdf", "_orig.pdf")
        path_SR = path_SR.replace(".pdf", "_orig.pdf")
        path_all = path_all.replace(".pdf", "_orig.pdf")

    max_order = 6
    signal_str = "signal" if not generic else "Signal_30_GeV"
    if "2024" in args.year:
        yr = "2024"
    else:
        yr = args.year
    
    savePlots(
        "BDT",
        #plot_data_sig_bkg_comparison(hists_reordered_SB, args.title, path_SB, Scales.sumw[yr][f"{int(args.hypMass) if args.hypMass != 0 else 30}_GeV"], maximum=(10**max_order - 5*10**(max_order-1)), mass=str(int(args.hypMass))+" GeV" if not generic else "30 GeV", reweight=True, sub="Data/Bkg restricted to SB", Scales = Scales, extra_ratio=extra_ratio),
        plot_data_sig_bkg_comparison(hists_reordered_all_dataBlind if args.blind else hists_reordered_all, args.title, path_all, Scales.sumw[yr][f"{int(args.hypMass) if args.hypMass != 0 else 30}_GeV"], maximum=(10**max_order - 5*10**(max_order-1)), mass=str(int(args.hypMass))+" GeV" if not generic else "30 GeV", reweight=True, sub="Data blinded" if args.blind else "", Scales = Scales, extra_ratio=extra_ratio),
        [f"{int(args.hypMass)}_GeV_{k}" if args.hypMass != 0 else f"generic_{k}" for k in list(hists_reordered_all_dataBlind.keys() if args.blind else hists_reordered_all.keys())],
        "paper_plots/histograms_dataSB.root" if args.blind else "paper_plots/histograms.root"
    )

    # Plot even/odd events comparison
    if not generic:
        even_sample = samples["signal"].iloc[lambda x: x.index % 2 == 0]
        evens = fillHist(
            bdt_score,
            even_sample,
            bounds = bounds,
            binsScale = (bins / 2) / (bounds[1] - bounds[0]),
            normalize = False,
            name = f"BDT_score_signal_{int(args.hypMass)}_even_GeV"
        )

        odd_sample = samples["signal"].iloc[lambda x: x.index % 2 == 1]
        odds = fillHist(
            bdt_score,
            odd_sample,
            bounds = bounds,
            binsScale = (bins / 2) / (bounds[1] - bounds[0]),
            normalize = False,
            name = f"BDT_score_signal_{int(args.hypMass)}_odd_GeV"
        )

        path = f"{cwd}/../../scripts/plots/BDT/{args.year}/{int(args.hypMass)}_GeV/{args.model}/even_odd_events{str(int(args.hypMass)) if not original else '_orig'}.pdf"
        plot_even_odd_BDT(path, evens, odds, even_sample["weight"].sum(), odd_sample["weight"].sum(), int(args.hypMass), args.year, Scales=Scales)


def plot_stats(
    model: str,
    year: int,
    d: str,
    eras: list = None,
    Scales: str = None,
    original: bool = False,
) -> None:

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Need to sort through all variations applied to the signal MC
    print(d, os.listdir(d)[0])
    variations = os.listdir(os.path.join(d, os.listdir(d)[0]))
    main_variations = set(["_".join(list(var.split("_")[:-1])) if "nominal" not in var else var for var in variations])
    print("===", ", ".join(main_variations), "===\n")

    bdt_score = "BDT_score_orig" if original else "BDT_score"
    samples = {}
    for variation in variations:
        print(f"~~~ Loading variation: {variation} ~~~")
        samples.update({variation: {}})
        print("Loading full dataset:")
        samples[variation].update(define_variables(d, reduce_df=True, applyReduce=True, perMass=True, variation=variation, extra=[bdt_score, "weight"]))

        to_del = []
        for key in range(15,65,5):
            samples[variation][f"{key}_GeV"] = samples[variation][key]
            to_del.append(key)
        for key in to_del:
            del samples[variation][key]

        print()

    # Specific mass processing
    for mass in [f"{m}_GeV" for m in range(15,65,5)]:
        print(f"~~~ Plotting samples for {mass} GeV mass point ~~~")
        for variation in main_variations:
            if "nominal" in variation:
                continue

            print(f"\nProcessing variation: {variation}")

            up = samples[f"{variation}_up"][mass]
            down = samples[f"{variation}_down"][mass]

            path = f"{cwd}/../../scripts/plots/BDT/{year}/{mass}/{model}/{variation}/"
            if not os.path.exists(path):
                os.mkdir(path)

            plot_variations(samples["nominal"][mass], up, down, path, Scales) 

        print()

