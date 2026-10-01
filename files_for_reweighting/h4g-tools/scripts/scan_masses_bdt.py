import time
import os
from h4g_tools.bdt.run_bdt import gen_BDT, run_samples_BDT, plot_results
from h4g_tools.utils.runner_utils import get_parser_scan, handle_errors
from h4g_tools.utils.loading import pathSetup
from h4g_tools.utils.scales import Scales as ScalesCls

"""
~~~ Must run the reweighting script first!!! ~~~

Change year and xgb model as necessary.

For generating BDT from parameters JSON:
python3 scan_masses_bdt.py --gen -fp -sig <sig-sample> -y 2024 -xgb ES_BDT_2024 

Without parameters:
python3 scan_masses_bdt.py --gen  -tnp -sig <sig-sample> -xgb ES_BDT_2024 -y 2024

For running BDT over samples:
python3 scan_masses_bdt.py --run -sig <sig-sample> -y 2024 -xgb ES_BDT_2024

For plotting:
python3 scan_masses_bdt.py --plot -y 2024 -xgb ES_BDT_2024

"""


if __name__ == "__main__":
    parser = get_parser_scan()
    parser.add_argument("-noRW", "--noRW", default=False, action="store_true", help="Skip reweighting. Useful for testing efficacy in BDT training.")
    parser.add_argument("-mvaCut", "--mvaCut", default=-1.0, type=float, required=False, help="Cut to apply on photon MVA ID score.")
    parser.add_argument("-extra_ratio", "--extra_ratio", default=False, action="store_true", help="Add ratio of unreweighted bkg scores to ratio plots.")
    parser.add_argument("-EB", "--EB_only", default=False, action="store_true", help="Use only events with all 4 photons in EB.")
    parser.add_argument("-bad", "--bad", default=False, action="store_true", help="Create artificially 'bad' BDT model.")
    parser.add_argument("-transform", "--transform", default=False, action="store_true", help="Transform BDT distributions for smoothing.")
    parser.add_argument("-tj", "--training_json", type=str, required=False, help="Provide a training json for input parameters.")
    parser.add_argument("--sigMix", default=False, action="store_true", help="Process mixed signal for checks.")
    parser.add_argument("-unblind", "--unblind", default=False, action="store_true", help="Unblind data BDT spectrum.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-includeStats", "--includeStats", default=False, action="store_true", help="Evaluate the statistical variations too.")
    group.add_argument("-var", "--variation", type=str, required=False, help="Only process a specific variation.")
    args = parser.parse_args()
    
    # Catch parser issues
    handle_errors(
        ["Need training samples.", args.gen, args.sigInput is None],
        ["Only need to provide signal samples and run reweighting before this.", args.gen, args.bkgInput is not None, args.data is not None],
        ["Need samples to run over.", args.run, args.data is None, args.sigInput is None, args.bkgInput is None],
        ["Cannot debug when running or plotting BDT results.", args.debug, args.plot or args.run],
        ["Must provide a year for processing.", args.year is None]
    )

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    pathSetup(args.year, args.xgb_name)

    # Scan over masses 15-60 GeV in steps of 5 GeV
    masses = [x for x in range(15, 65, 5)]
    
    # Scales instance setup
    Scales = ScalesCls(args.year) 

    if args.year == "2022preEE":
        eras={
            "data": ["Run2022C", "Run2022D"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_preEE" for m in masses]
        }
    elif args.year == "2022postEE":
        eras={
            "data": ["Run2022E", "Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_postEE" for m in masses]
        }
    elif args.year == "2022":
        eras={
            "data": ["Run2022C", "Run2022D", "Run2022E", "Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D", "RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preEE", "postEE"]]
        }
    elif args.year == "2023preBPix":
        eras={
            "data": ["Run2023C"],
            "bkg": ["RunEvtMix2023C"],
            "signal": [f"Signal_{m}_GeV_preBPix" for m in masses]
        }
    elif args.year == "2022postBPix":
        eras={
            "data": ["Run2022D"],
            "bkg": ["RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_postBPix" for m in masses]
        }
    elif args.year == "2023":
        eras={
            "data": ["Run2023C", "Run2023D"],
            "bkg": ["RunEvtMix2023C", "RunEvtMix2023D"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preBPix", "postBPix"]]
        }
    elif "2024" in args.year:
        eras={
            "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
            "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
            "signal": [f"Signal_{m}_GeV" for m in masses]
        }
    else:
        assert False, "Must specify year in Scales class!"
    # ADD NEW YEAR HERE!

    start_time = time.time()
    if args.gen:
        # Don't pass eras for generating. Samples are already split b/c of reweighting reasons. Still need year for naming
        gen_BDT(args.year, args.sigInput.replace("year", args.year), 1.0, args.debug, args.from_params, args.train_no_params, args.xgb_name, eras = eras, Scales=Scales, noRW=args.noRW, mva_cut=args.mvaCut, EB_only=args.EB_only, bad=args.bad, training_json=args.training_json)

    elif args.run:
        masses = [x for x in range(15,63, 5)]
        masses_data = [x for x in range(15, 63, 1)]
        print()

        if args.transform:
            print("Deriving transformation...")
            Scales.transform = args.transform
            Scales.PrepareTransformation(args.sigInput.replace("year", args.year), eras, args.year, args.xgb_name)
            print()

        print("Running BDT model over signal samples...")
        run_samples_BDT(masses, args.year, args.sigInput.replace("year", args.year), args.xgb_name, sig=True, eras = eras, Scales=Scales, nominalOnly = not args.includeStats and args.variation is None, pick_variation=args.variation, sig_mix=args.sigMix)

        if not args.sigMix:
            print("Running BDT model over background samples...")
            run_samples_BDT(masses, args.year, "{cwd}/../../scripts/bdtIO/BDT_year_input_weights.parquet".replace("year", args.year), args.xgb_name, eras = eras, Scales=Scales)

            print("Running BDT model over data samples...")
            run_samples_BDT(masses_data, args.year, "{cwd}/../../scripts/bdtIO/BDT_year_input_weights_data.parquet".replace("year", args.year), args.xgb_name, eras = eras, Scales=Scales)

    elif args.plot:
        masses = [x for x in range(15,63, 5)]
        for original in [args.transform, not args.transform] if args.transform else [args.transform]:
            plot_results(masses, args.year, args.xgb_name, f"{cwd}/../../scripts/outputs/BDT/", f"{cwd}/../../scripts/outputs/BDT/", f"{cwd}/../../scripts/outputs/BDT/", Scales = Scales, extra_ratio=args.extra_ratio, EB_only=args.EB_only, nominalOnly = not args.includeStats, original = original, blind = not args.unblind)
 
    elapsed_time = time.time() - start_time
    print(f"Total elapsed time: {elapsed_time}s")
