from ROOT import TH1D, gStyle, RDataFrame, RDataFrame, TTree, TChain, SetOwnership
from ROOT.std import vector
import numpy as np
import pandas as pd
import copy
import os
import json
import gc
from typing import List, Dict, Optional, Union
import pyarrow
import pyarrow.parquet as pq
from pyarrow.lib import ArrowInvalid
import fnmatch
import warnings
from array import array
import h4g_tools.utils.scales as scales
warnings.filterwarnings('ignore')


def pathSetup(year: str, model: str = None) -> None:
    """
    Makes directories as needed and prepares for plots.
    """

    # Get directory of file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Check main plotting directory
    if not os.path.exists(f"{cwd}/../../scripts/plots/"):
        os.mkdir(f"{cwd}/../../scripts/plots/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/combine_fit/"):
        os.mkdir(f"{cwd}/../../scripts/plots/combine_fit/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/gen/"):
        os.mkdir(f"{cwd}/../../scripts/plots/gen/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/gen/{year}/"):
        os.mkdir(f"{cwd}/../../scripts/plots/gen/{year}/")
    
    # Check cats files directory
    if not os.path.exists(f"{cwd}/../../scripts/cats/"):
        os.mkdir(f"{cwd}/../../scripts/cats/")
    if not os.path.exists(f"{cwd}/../../scripts/cats/{year}"):
        os.mkdir(f"{cwd}/../../scripts/cats/{year}")

    # Check signal modelling directory
    if not os.path.exists(f"{cwd}/../../scripts/plots/gauss/"):
        os.mkdir(f"{cwd}/../../scripts/plots/gauss/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/gauss/{year}/"):
        os.mkdir(f"{cwd}/../../scripts/plots/gauss/{year}/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/dcb/"):
        os.mkdir(f"{cwd}/../../scripts/plots/dcb/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/dcb/{year}/"):
        os.mkdir(f"{cwd}/../../scripts/plots/dcb/{year}/")

    # Check BDT plots directory
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/sig_mix/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/sig_mix/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/m4g_bdt_checks/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/m4g_bdt_checks/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/compare/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/compare/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/compare_all/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/compare_all/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/bdt_inputs/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/bdt_inputs/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/chi2/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/chi2/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/Ndimreweights/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/Ndimreweights/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/BDT_pred_generic/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/BDT_pred_generic/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/compare_mva_lt_m0p85/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/compare_mva_lt_m0p85/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/compare_mva_gt_m0p85/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/compare_mva_gt_m0p85/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/compare_mva_lt_m0p85_linear/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/compare_mva_lt_m0p85_linear/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/compare_mva_gt_m0p85_linear/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/compare_mva_gt_m0p85_linear/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/compare_reweight/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/compare_reweight/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/cats/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/cats/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/compare_corrected/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/compare_corrected/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/compare_corrected/bkg/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/compare_corrected/bkg/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/compare_corrected/data/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/compare_corrected/data/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/mva_correlation/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/mva_correlation/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/cuts/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/cuts/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/cuts/compare/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/cuts/compare/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/cuts/mva_correlation/"):
        os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/cuts/mva_correlation/")

    # Check BDT models directory
    if not os.path.exists(f"{cwd}/../../models/"):
        os.mkdir(f"{cwd}/../../models/")

    # Check JSON directory for BDT
    if not os.path.exists(f"{cwd}/../../scripts/bdtIO/"):
        os.mkdir(f"{cwd}/../../scripts/bdtIO/")

    # Check outputs directory
    if not os.path.exists(f"{cwd}/../../scripts/outputs/"):
        os.mkdir(f"{cwd}/../../scripts/outputs/")
    if not os.path.exists(f"{cwd}/../../scripts/outputs/BDT/"):
        os.mkdir(f"{cwd}/../../scripts/outputs/BDT/")
    if not os.path.exists(f"{cwd}/../../scripts/outputs/cats/"):
        os.mkdir(f"{cwd}/../../scripts/outputs/cats/")
    if not os.path.exists(f"{cwd}/../../scripts/outputs/BDT/{year}/"):
        os.mkdir(f"{cwd}/../../scripts/outputs/BDT/{year}/")
    if not os.path.exists(f"{cwd}/../../scripts/outputs/cats/{year}/"):
        os.mkdir(f"{cwd}/../../scripts/outputs/cats/{year}/")
    
    # Check standard plots directory
    if not os.path.exists(f"{cwd}/../../scripts/plots/standard/"):
        os.mkdir(f"{cwd}/../../scripts/plots/standard/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/standard/{year}/"):
        os.mkdir(f"{cwd}/../../scripts/plots/standard/{year}/")
    if not os.path.exists(f"{cwd}/../../scripts/plots/standard/{year}/compare"):
        os.mkdir(f"{cwd}/../../scripts/plots/standard/{year}/compare/")
    
    # Check individual mass plots directory
    for mass in [f"{x}_GeV" for x in range(15,63,5)]:
        if not os.path.exists(f"{cwd}/../../scripts/plots/standard/{year}/{mass}/"):
            os.mkdir(f"{cwd}/../../scripts/plots/standard/{year}/{mass}/")
        if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/{mass}/"):
            os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/{mass}/")
        if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/cuts/{mass}/"):
            os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/cuts/{mass}/")
        if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/compare_corrected/{mass}/"):
            os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/compare_corrected/{mass}/")
        if not os.path.exists(f"{cwd}/../../scripts/plots/gen/{year}/{mass}/"):
            os.mkdir(f"{cwd}/../../scripts/plots/gen/{year}/{mass}/")
        if model is not None:
            if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/BDT_pred_generic/{model}"):
                os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/BDT_pred_generic/{model}")
            if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/{mass}/{model}"):
                os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/{mass}/{model}")
    
    if model is not None:
        if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/{model}/"):
            os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/{model}/")
        if not os.path.exists(f"{cwd}/../../scripts/outputs/BDT/{year}/{model}/"):
            os.mkdir(f"{cwd}/../../scripts/outputs/BDT/{year}/{model}/")
        if not os.path.exists(f"{cwd}/../../scripts/outputs/BDT/{year}/{model}/generic/"):
            os.mkdir(f"{cwd}/../../scripts/outputs/BDT/{year}/{model}/generic/")
        if not os.path.exists(f"{cwd}/../../scripts/cats/{year}/{model}"):
            os.mkdir(f"{cwd}/../../scripts/cats/{year}/{model}")
        if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/cats/{model}"):
            os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/cats/{model}")

        for sub in ["signal", "background", "data"]:
            if not os.path.exists(f"{cwd}/../../scripts/outputs/BDT/{year}/{model}/{sub}/"):
                os.mkdir(f"{cwd}/../../scripts/outputs/BDT/{year}/{model}/{sub}/")
            if not os.path.exists(f"{cwd}/../../scripts/outputs/BDT/{year}/{model}/generic/{sub}/") and sub != "signal":
                os.mkdir(f"{cwd}/../../scripts/outputs/BDT/{year}/{model}/generic/{sub}/")

            if sub in ["signal", "background"]:
                for mass in [f"{x}_GeV" for x in range(15,63,5)]:
                    if not os.path.exists(f"{cwd}/../../scripts/outputs/BDT/{year}/{model}/{sub}/{mass}/"):
                        os.mkdir(f"{cwd}/../../scripts/outputs/BDT/{year}/{model}/{sub}/{mass}/")
            elif sub == "data":
                for mass in [f"{x}_GeV" for x in range(15,63,1)]:
                    if not os.path.exists(f"{cwd}/../../scripts/outputs/BDT/{year}/{model}/{sub}/{mass}/"):
                        os.mkdir(f"{cwd}/../../scripts/outputs/BDT/{year}/{model}/{sub}/{mass}/")

    for sub in ["signal", "background", "data"]:
        if not os.path.exists(f"{cwd}/../../scripts/outputs/cats/{year}/{sub}/"):
            os.mkdir(f"{cwd}/../../scripts/outputs/cats/{year}/{sub}/")
        if not os.path.exists(f"{cwd}/../../scripts/plots/BDT/{year}/cuts/mva_correlation/{sub}/"):
            os.mkdir(f"{cwd}/../../scripts/plots/BDT/{year}/cuts/mva_correlation/{sub}/")
        for mass in [f"{x}_GeV" for x in range(15,63,5)]:
            if not os.path.exists(f"{cwd}/../../scripts/outputs/cats/{year}/{sub}/{mass}/"):
                os.mkdir(f"{cwd}/../../scripts/outputs/cats/{year}/{sub}/{mass}/")


def loadROOTSamples(
    directory: Union[str, list],
    eras: dict = None,
    tree_name: Optional[str] = "Events",
    use_uproot: Optional[bool] = False,
    mvaCut: float = None
#) -> dak.Array:
):
    """
    Loads in samples from ROOT files. Provide directory as list to load individual files directly.
    """

    do_eras = False
    era_type = ""
    if eras is not None and type(eras) is dict:
        do_eras = True
        if "bkg" in directory.split("/")[-2] or "background" in directory.split("/")[-2]:
            print(f"Loading eras: {', '.join(eras['bkg'])}")
            era_type = "bkg"
        elif "data" in directory.split("/")[-2]:
            print(f"Loading eras: {', '.join(eras['data'])}")
            era_type = "data"
        elif "signal" in directory.split("/")[-2]:
            print(f"Loading eras: {', '.join(eras['signal'])}")
            era_type = "signal"
        else:
            assert era_type not in ["data", "bkg", "signal"]

    if type(directory) is str and "root:/" not in directory:
        assert os.path.exists(directory), f"Directory ({directory}) does not exist."

    files_set = set()

    if type(directory) is list:
        for f in directory:
            files_set.add(f)

    else:
        # Load from directory
        use_nominal = False
        if ".root" not in directory:
            proceed = True
            # Scan directory and fill samples_files with contents of parquet files
            for (root, dirs, files) in os.walk(directory):
                if do_eras:
                    proceed = False
                    for era in eras[era_type]:
                        if era in root:
                            if "signal" in era.lower():
                                use_nominal = True
                            proceed = True
                            break
                if proceed:            
                    for name in files:
                        if os.path.isfile(os.path.join(root, name)) and ".root" in name:
                            if use_nominal:
                                if "nominal" in root:
                                    files_set.add(os.path.join(root, name))
                            else:
                                files_set.add(os.path.join(root, name))


    # Collect previous objects and disable circular garbage collector (ref counter still works)
    gc.collect()
    gc.disable()

    if not use_uproot:
        tree = TChain(tree_name)
        for file in files_set:
            tree.Add(file)
        print(f"Loaded {len(files_set)} ROOT files with {tree.GetEntries()} events.")
    """
    else:
        # Check this!
        df_list = []
        for file in files_set:
            df_list.append(uproot.open(file)["Events"].pandas.df()
        out_df = pd.concatenate(df_list, axis=0)
        print(f"Loaded {len(files_set)} ROOT files with {len(out_df)} events.")
    """

    assert type(tree) is TChain
    assert tree.GetEntries() != 0, "Number of events = 0"

    return tree



def loadROOTPerMass(
    inputs: str,
    tree: Optional[str] = "Events"
#) -> Dict[str, dak.Array]:
):
    """
    Loads ROOT samples in a per-mass basis.
    """

    assert os.path.exists(inputs), "Input does not exist!"
    
    samples = {}

    if ".json" in inputs:
        with open(inputs) as f:
            js = json.load(f)
            js_masses = [m.replace("Signal_","") for m in js.keys()]

        for mass in js_masses:
            print(f"Loading {mass}:")
            samples.update({mass: loadROOTSamples(js[f"Signal_{mass}"], tree)})

    else:
        # Load samples with mass branch
        samples = {}
        for (root, dirs, files) in os.walk(inputs):
            for dname in dirs:
                full_path = os.path.join(root, dname)
                if os.path.isdir(full_path) and "GeV" in dname:
                    mass = dname.strip("/").replace("Signal_","")
                    print(f"Loading {mass}:")
                    samples.update({mass: loadROOTSamples(full_path, tree)})
    
    print()
    return samples


def loadSamples(
    directory: str,
    branches: Optional[List[str]] = None,
    bdt: bool = False,
    eras: dict = None,
    era_type_arg = None,
    meta: bool = False,
    skipInterMass: bool = False,
    variation: str = "nominal"
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
    if ".parquet" not in directory:
        # Scan directory and fill samples_files with contents of parquet files
        proceed = True
        walk = os.walk(directory)
        for (root, dirs, files) in walk:
            if do_eras:
                proceed = False
                for era in eras[era_type]:
                    if era in root:
                        proceed = True
                        break
            if proceed:            
                for name in files:
                    if os.path.isfile(os.path.join(root, name)) and (variation in os.path.join(root, name) if variation is not None else True) and ".parquet" in name:
                        try:
                            if not meta:
                                samples_files.append(pd.read_parquet(os.path.join(root, name), columns=branches))
                            else:
                                samples_files.append(pq.read_table(os.path.join(root, name), columns=branches))
                        except ArrowInvalid as err:
                            print(f"Could not load {name}")
    
    # Load from file
    elif ".parquet" in directory:
        try:
            if not meta:
                samples_files.append(pd.read_parquet(directory, columns=branches))
            else:
                samples_files.append(pq.read_table(directory, columns=branches))
        except ArrowInvalid:
            print(f"Could not load {name}")

    # Fill df with all file contents
    assert len(samples_files) != 0, f"Samples size is 0!\t {directory}"
    if not meta:
        df = pd.concat(samples_files, axis=0, ignore_index=True)
    else:
        # Add metadata together
        meta_combined = dict.fromkeys([k for k in samples_files[0].schema.metadata.keys() if "Nphotons" not in k.decode("utf-8")], 0)
        if b"ak:parameters" in meta_combined.keys():
            del meta_combined[b"ak:parameters"]
        if samples_files[0].schema.metadata[b'sum_genw_presel'] == b'Data':
            del meta_combined[b'sum_genw_presel']
        for sample in samples_files:
            for key in sample.schema.metadata.keys():
                if key not in meta_combined.keys():
                    continue
                meta_combined[key] += float(sample.schema.metadata[key])
                #if key == b'initial_events':
                #    print(meta_combined[key])

        for key in meta_combined.keys():
            meta_combined[key] = str(meta_combined[key]).encode("utf-8")

        # NOTE: The branches in each file are slightly different so the concat will break. Load only one branch if you only need metadata! (e.g. pho1_pt)
        # Compare with plotEfficiency.py in HiggsDNA
        table = pyarrow.concat_tables(samples_files)
        table = table.replace_schema_metadata(meta_combined)

    # Now remove any known NaN columns
    for column in df.columns if not meta else table.columns:
        if "veto" in column or "superclusterEta" in column:
            if not meta:
                df = df.drop(column, axis=1)
            else:
                table = table.drop_columns(column)

    # Add interMass_abs fields
    if not meta and not skipInterMass:
        if not bdt and branches is None:
            df = addInterMassAbs(df)
        if branches is not None:
            if "LeadPs_interMass" in branches and "SubleadPs_interMass" in branches:
                df = addInterMassAbs(df)

    if not meta:
        assert df.shape[0] != 0, "Number of events = 0"
    else:
        assert table.num_rows != 0, "Number of events = 0"

    schema = "pd.DataFrame" if not meta else "pyarrow.table"
    print(f"Loaded {len(samples_files)} parquet file(s) with {len(df if not meta else table)} events using {schema} schema.")

    return df if not meta else table


def loadPerMass(
    inputs: str,
    branches: Optional[List[str]] = None,
    bdt: bool = False,
    eras: dict = None,
    year: str = None,
    era_type_arg = "",
    meta: bool = False,
    skipInterMass: bool = False,
    variation: str = "nominal",
    bdt_outputs: bool = False,
    nominalMassOnly: bool = True
) -> Dict[str, pd.DataFrame]:
    """
    Loads samples in a per-mass basis.
    """

    assert type(inputs) is str
    assert os.path.exists(inputs), f"Input does not exist! {inputs}"

    do_eras = False
    era_type = "" if era_type_arg is None else era_type_arg
    if eras is not None and type(eras) is dict:
        do_eras = True
        if "data" in inputs and "signal" not in inputs and "bkg" not in inputs and "background" not in inputs:
            print(f"Loading eras: {', '.join(eras['data'])}")
            era_type = "data"
        elif "bkg" in inputs or "background" in inputs:
            print(f"Loading eras: {', '.join(eras['bkg'])}")
            era_type = "bkg"
        elif "signal" in inputs:
            if bdt_outputs:
                eras_out = [f"Signal_{e.split('/')[-1]}" + (f"_{e.split('/')[0][4:]}" if len(f"{e.split('/')[0][4:]}") > 0 else "") for e in eras["signal"]]
                print(f"Loading eras: {', '.join(eras_out)}")
            else:
                print(f"Loading eras: {', '.join(eras['signal'])}")
            era_type = "signal"
        else:
            assert era_type in ["data", "bkg", "signal"], "era_type must be either signal, bkg, or data!"

    # Load samples with mass branch
    samples = {}
    for (root, dirs, files) in os.walk(inputs):
        for dname in dirs:
            full_path = os.path.join(root, dname)
            # Need to avoid grabbing postprocessed "merged" and "root" directories, but also can't skip over my merged samples for plotting
            if os.path.isdir(full_path) and "GeV" in dname and "/merged/" not in full_path and "/root/" not in full_path:
                proceed = True
                if do_eras:
                    proceed = False
                    for era in eras[era_type]:
                        if era in full_path:
                            proceed = True
                            break
                if proceed:
                    mass = dname.strip("/")
                    if nominalMassOnly and mass[mass.find("GeV")-3:mass.find("GeV")+3].replace("Signal_","") not in [f"{m}_GeV" for m in range(15,65,5)]:
                        continue
                    print(f"Loading {mass}:")
                    samples.update({mass: loadSamples(full_path, branches, bdt, meta=meta, skipInterMass=skipInterMass, variation=variation)})
    
    samples_merged = {}
    if do_eras:
        if (year == "2022" or year == "2023" or any(["postEE" in k for k in samples.keys()]) or any(["postBPix" in k for k in samples.keys()])):
            for m in (range(15,63,5) if nominalMassOnly else range(15,63,1)):
                filtered_keys = [keys for keys in samples.keys() if f"{m}_GeV" in keys]
                if not meta:
                    samples_merged.update({
                        f"Signal_{m}_GeV" if "Signal_" in filtered_keys[0] else f"{m}_GeV": pd.concat([copy.deepcopy(samples[fk]) for fk in filtered_keys], axis=0)
                    })
                else:
                    samples_merged.update({
                        f"Signal_{m}_GeV" if "Signal_" in filtered_keys[0] else f"{m}_GeV": pyarrow.concat_tables([copy.deepcopy(samples[fk]) for fk in filtered_keys])
                    })
                    
                    # Sum the metadata too
                    generic_key = f"Signal_{m}_GeV" if "Signal_" in filtered_keys[0] else f"{m}_GeV"
                    meta_combined = dict.fromkeys(list(samples_merged[generic_key].schema.metadata.keys()), 0)
                    for fk in filtered_keys:
                        for key in samples[fk].schema.metadata.keys():
                            if key not in meta_combined.keys():
                                continue
                            meta_combined[key] += float(samples[fk].schema.metadata[key])
                            #if key == b'initial_events':
                            #    print(meta_combined[key])

                    for key in meta_combined.keys():
                        meta_combined[key] = str(meta_combined[key]).encode("utf-8")
                    
                    # Assign summed metadata to combined samples
                    samples_merged[generic_key] = samples_merged[generic_key].replace_schema_metadata(meta_combined)
                
        else:
            samples_merged = copy.deepcopy(samples)
    else:
        samples_merged = copy.deepcopy(samples)
    
    if not meta:
        assert np.all([df.shape[0] != 0 for df in samples_merged.values()]), "Number of events = 0"
    else:
        assert np.all([table.num_rows != 0 for table in samples_merged.values()]), "Number of events = 0"

    return samples_merged


def fillHist(
    branch: str,
    sample: pd.DataFrame,
    bounds: List[float],
    name: Optional[str] = None,
    binsScale: Optional[float] = 1.0,
    normalize: Optional[bool] = True,
    customWeights: Optional[pd.DataFrame] = None,
    includeStats: Optional[bool] = False,
    noWeights: Optional[bool] = False,
    sb: Optional[bool] = False,
    preventOverFlow: bool = False,
    bins: list = None
) -> TH1D:
    """
    Fills histogram with data from samples dataframe. Binning and edges are based on bounds input.
    """

    if sb:
        sample = scales.SideBand(sample)

    # Set name to branch if not provided
    if name is None:
        name = branch

    # Histogram with custom number of equal bins, bounded by (bounds[0], bounds[1])
    if bins is None:
        hist = TH1D(name, name, int((bounds[1]-bounds[0])*binsScale), bounds[0], bounds[1])
    else:
        hist = TH1D(name, name, int((bounds[1]-bounds[0])*binsScale), array("d", bins))

    if type(sample) is pd.DataFrame or type(sample) is pd.Series:
        branch_df = sample[branch].tolist()
    elif type(sample) is list:
        branch_df = sample
    else:
        print(type(sample))

    try:
        if type(sample) is pd.DataFrame:
            weights = sample["weight"].tolist()
        else:
            raise KeyError
    except KeyError as err:
        weights = np.ones_like(branch_df)
        if not noWeights:
            print("~~~ WARNING ~~~ Could not find weights! Using 1.0 for weights.")

    if customWeights is not None:
        if type(customWeights) is not str:
            weights = customWeights.tolist() if type(customWeights) is pd.DataFrame else customWeights
        else:
            weights = sample[customWeights]
        print("Using additional custom weights.")
    
    # Fill histogram
    for idx, (data, weight) in enumerate(zip(branch_df, weights)):
        if type(data) is np.ndarray:
            data2 = data.tolist()
            weights2 = [weight] * len(data2)
            for d2, w2 in zip(data2, weights2):
                hist.Fill(d2, w2)
        elif noWeights:
            hist.Fill(data)
        else:
            hist.Fill(data, weight)

    # Normalize
    if normalize:
        if hist.Integral("width") > 0.0:
            hist.Scale(1.0 / hist.Integral())
        else:
            hist.Scale(1.0)
    
    # Turn of stat box by default
    if not includeStats:
        gStyle.SetOptStat(0)

    # Fill first bin with underflow and set underflow to zero
    if not preventOverFlow:
        hist.SetBinContent(1, hist.GetBinContent(0) + hist.GetBinContent(1))
        hist.SetBinContent(0, 0)

    # Fill last bin with overflow and set overflow to zero
    if not preventOverFlow:
        hist.SetBinContent(hist.GetNbinsX(), hist.GetBinContent(hist.GetNbinsX()) + hist.GetBinContent(hist.GetNbinsX()+1))
        hist.SetBinContent(hist.GetNbinsX()+1, 0)

    assert hist.GetEntries() != 0, f"No entries in {branch} histogram within bounds! {sample[branch]}"
    note = f"Filled {branch} histogram with {int(hist.GetEntries())} events ({hist.Integral()} with weights and bounds)."
    print(f"{note:<50}\t Bounds: [{bounds[0]},{bounds[1]}]\t NBins: {int((bounds[1]-bounds[0])*binsScale)}")

    SetOwnership(hist, False)

    return hist

def addInterMassAbs(
    df: pd.DataFrame
) -> pd.DataFrame:
    """
    Adds interMass_abs field for leading and subleading pseudoscalars.
    """

    df["LeadPs_interMass_abs"] = abs(df["LeadPs_interMass"])
    df["SubleadPs_interMass_abs"] = abs(df["SubleadPs_interMass"])

    return df
