from ROOT import TChain
import pandas as pd
import os
import gc
from typing import List, Dict, Optional, Union
import pyarrow
from pyarrow.lib import ArrowInvalid
import warnings
import argparse
warnings.filterwarnings('ignore')


def loadROOTSamples(
    directory: Union[str, list],
    eras: dict = None,
    tree_name: Optional[str] = "Events",
):
    """
    Loads in samples from ROOT files. Provide directory as list to load individual files directly.
    """

    do_eras = False
    era_type = ""
    if eras is not None and type(eras) is dict:
        do_eras = True
        if "bkg" in directory or "background" in directory:
            print(f"Loading eras: {', '.join(eras['bkg'])}")
            era_type = "bkg"
        elif "data" in directory:
            print(f"Loading eras: {', '.join(eras['data'])}")
            era_type = "data"
        elif "signal" in directory:
            print(f"Loading eras: {', '.join(eras['signal'])}")
            era_type = "signal"
        else:
            assert era_type not in ["data", "bkg", "signal"]

    if type(directory) is str and "root://" not in directory:
        assert os.path.exists(directory), f"Directory ({directory}) does not exist."

    files_set = set()

    # Load from directory
    if ".root" not in directory:
        proceed = True
        # Scan directory and fill samples_files with contents of parquet files
        for (root, dirs, files) in os.walk(directory):
            if do_eras:
                proceed = False
                for era in eras[era_type]:
                    if era in root:
                        proceed = True
                        break
            if proceed:            
                for name in files:
                    if os.path.isfile(os.path.join(root, name)) and ".root" in name:
                        files_set.add(os.path.join(root, name))

    # Collect previous objects and disable circular garbage collector (ref counter still works)
    gc.collect()
    gc.disable()

    tree = TChain(tree_name)
    for file in files_set:
        tree.Add(file)
    print(f"Loaded {len(files_set)} ROOT files with {tree.GetEntries()} events.")

    assert type(tree) is TChain
    assert tree.GetEntries() != 0, "Number of events = 0"

    return tree


def loadSamples(
    directory: str,
    eras: dict = None,
    era_type_arg = None,
) -> pd.DataFrame:
    """
    Scans directory (or a single file) and reads in parquet files. Creates list of dataframes with parquet file contents and merges them.
    """

    assert os.path.exists(directory), f"Directory ({directory}) does not exist."

    samples_files = []
    do_eras = False
    era_type = "" if era_type_arg is None else era_type_arg
    if eras is not None and type(eras) is dict:
        do_eras = True
        if "data" in directory and "signal" not in directory and "bkg" not in directory and "background" not in directory:
            print(f"Loading eras: {', '.join(eras['data'])}")
            era_type = "data"
        elif "bkg" in directory or "background" in directory:
            print(f"Loading eras: {', '.join(eras['bkg'])}")
            era_type = "bkg"
        elif "signal" in directory:
            print(f"Loading eras: {', '.join(eras['signal'])}")
            era_type = "signal"
        else:
            assert era_type not in ["data", "bkg", "signal"]

    # Load from directory
    proceed = True
    walk = os.walk(directory)
    # Scan directory and fill samples_files with contents of parquet files
    for (root, dirs, files) in walk:
        if do_eras:
            proceed = False
            for era in eras[era_type]:
                if era in root:
                    proceed = True
                    break
        if proceed:            
            for name in files:
                if os.path.isfile(os.path.join(root, name)) and "nominal" in os.path.join(root, name) and ".parquet" in name:
                    try:
                        samples_files.append(pd.read_parquet(os.path.join(root, name)))
                    except ArrowInvalid as err:
                        print(f"Could not load {name}")
    

    # Fill df with all file contents
    assert len(samples_files) != 0, f"Samples size is 0!\t {directory}"
    df = pd.concat(samples_files, axis=0, ignore_index=True)

    # Now remove any known NaN columns
    for column in df.columns:
        if "veto" in column or "superclusterEta" in column:
            df = df.drop(column, axis=1)

    assert df.shape[0] != 0, "Number of events = 0"

    schema = "pd.DataFrame"
    print(f"Loaded {len(samples_files)} parquet file(s) with {len(df)} events using {schema} schema.")

    return df


def loadingHelper(args):

    # Eras definitions
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
    elif args.year == "2023postBPix":
        eras={
            "data": ["Run2023D"],
            "bkg": ["RunEvtMix2023D"],
            "signal": [f"Signal_{m}_GeV_postBPix" for m in masses]
        }
    elif args.year == "2023":
        eras={
            "data": ["Run2023C", "Run2023D"],
            "bkg": ["RunEvtMix2023C", "RunEvtMix2023D"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preBPix", "postBPix"]]
        }
    elif args.year == "2024":
        eras={
            "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
            "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
            "signal": [f"Signal_{m}_GeV" for m in masses]
        }
    else:
        assert False, "Must specify year!"

    return eras

    if args.dataType == "root":
        samples = loadROOTSamples(args.inputs, eras)

    elif args.dataType == "parquet":
        samples = loadSamples(args.inputs, eras)

    return samples


if __name__ == "__main__":

    ## Example use cases ##
    # python3 loadingScript.py root <samples-directory> -y 2022C
    # python3 loadingScript.py parquet <samples-directory> -y 2022postEE

    ####  Use this if you want command-line arguments  ####
    parser = argparse.ArgumentParser(description="Command line options parser")
    parser.add_argument("dataType", choices=["root", "parquet"], help="Specify how to load files.")
    parser.add_argument("-i", "--inputs", default=None, help="Input directory for samples.")
    parser.add_argument("-y", "--year", type=str, required=True, help="Year/eras of provided samples.")
    args = parser.parse_args()

    ####   ALTERNATIVE METHOD BELOW   ####

    ####  Use this if you want to specify manually via some script  ####
    """
    args = argparse.Namespace()
    args.dataType = "root"  # or "parquet"
    args.inputs = "<samples-directory>"
    args.year = "2022postEE"  # or any other era
    """

    ####  loadingHelper will return either a TTree or a Pandas DataFrame depending on what you ask it for.  ####
    ## You can call it as many times as you want to get all eras for bkg and data separately. Remember to use an updated "args" for each call!! ##
    samples = loadingHelper(args)
    for key in samples.keys():
        print(key, sum(samples[key]["genWeight"]))
