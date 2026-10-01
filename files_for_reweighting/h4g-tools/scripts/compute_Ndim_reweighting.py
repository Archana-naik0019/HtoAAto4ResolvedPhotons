import argparse
import gc
from h4g_tools.utils.runner_utils import handle_errors
from h4g_tools.utils.loading import pathSetup
from h4g_tools.bdt.reweight import performNDimReweighting


"""
Directories must already be split into eras, BUT you must specify year for naming reasons.
MUST MERGE SAMPLES USING convert_parquet_to_root_Haa_oneFile.py FIRST!

To run:
python3 compute_Ndim_reweighting.py -bkg /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/<merged-root-dir> -d /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/<merged-root-dir> -y <year> [-dw to include default weight in filename, -dwf to specify default weight, -nb to specify number of bins in original 4-dim variables]
"""


if __name__ == "__main__":
    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser", conflict_handler="resolve")
    parser.add_argument("-bkg", "--bkgInput", help="Directory to import background samples from.", type=str, default=None, required=False)
    parser.add_argument("-d", "--data", help="Directory to import data samples from.", type=str, default=None, required=False)
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    parser.add_argument("-nb", "--nBins", type=int, required=False, help="Number of bins to use (excluding under/overflow).")
    parser.add_argument("-mvaCut", "--mvaCut", type=float, required=False, help="Apply cut to photons 3/4 MVA ID >= mvaCut.")
    parser.add_argument("-EB", "--EB_only", default=False, action="store_true", help="Use only events with all 4 photons in EB.")
    parser.add_argument("-dw", "--default_weight", default=False, action="store_true", help="Include default weight in filename.")
    parser.add_argument("-dwf", "--default_weight_float", type=float, default=0.1, required=False, help="Default weight to use.")

    args = parser.parse_args()
    
    # Catch parser issues
    handle_errors(
        ["Need samples to run over.", args.data is None and args.bkgInput is None],
        #["Must provide a year for processing.", args.year is None],
    )

    pathSetup(args.year)

    # b/c of merging, no need to specify eras other than for naming
    masses = [x for x in range(15, 65, 5)]
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
            "data": ["Run2024C"], #, "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
            "bkg": ["RunEvtMix2024C"], #, "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
            "signal": [f"Signal_{m}_GeV" for m in masses]
        }
    else:
        assert False, "Must specify year!"

    # Account for year in sample name. This must be consistent if the replace should work
    args.bkgInput = args.bkgInput.replace("year", args.year)
    args.data = args.data.replace("year", args.year)

    # Collect previous objects and disable circular garbage collector (ref counter still works)
    gc.collect()
    gc.disable()

    performNDimReweighting(args.bkgInput, args.data, args.year, eras, nBins=args.nBins, mvaCut=args.mvaCut, EB_only=args.EB_only, save_default=args.default_weight, default_weight=args.default_weight_float)
