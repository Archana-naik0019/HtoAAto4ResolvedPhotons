import os
from array import array
from numba import jit
import uproot
import numpy as np
import pandas as pd
import itertools
from typing import List, Dict, Union, Tuple
from ROOT import TH1D, TMath, TTree, TFile, gROOT, TCanvas, kRed, kBlue, kGreen, TLegend, kBlack, TH2D, gStyle, SetOwnership
from h4g_tools.bdt.event_selection_bdt import define_variables
from h4g_tools.utils.loading import loadROOTSamples
from h4g_tools.utils.plotting import plot_Nreweighting, plot_Ndim_weights, loadPlottingParameters, plotFancy
from h4g_tools.utils.scales import Scales


def reweight1Dim(
    bkg: pd.DataFrame,
    data: pd.DataFrame,
    col: str = None,
) -> TH1D:
    """
    Reweights background events to data in sidebands to improve BDT score agreement. These weights are thrown out after use here since we only need signal and data after this.
    """

    nBins = 50
    BDT_score_hist = TH1D("BDT_score_bkg", "BDT_score_bkg", nBins, 0.0, 1.0)

    perm = []
    # Bins + 2 to account for under/overflow bins
    for bn in range(0,nBins+2):
        perm.append(bn)

    weight_data = {}
    weight_bkg = {}
    weight_ratio = {}
    for subset in perm:
        weight_data[subset] = 0.0
        weight_bkg[subset] = 0.0
        weight_ratio[subset] = 1.0

    # Use the ratio of the cut region that is being reweighted!
    sb_bkg_entries = len(bkg[(bkg["mass_gggg"] <= 115.0) | (bkg["mass_gggg"] >= 135.0)])
    sb_data_entries = len(data[(data["mass_gggg"] <= 115.0) | (data["mass_gggg"] >= 135.0)])

    sb_scaling_weight = 1.0
    if sb_bkg_entries > 0.0:
        sb_scaling_weight = sb_data_entries / sb_bkg_entries

    #@jit
    def sparseBins1D(
        hist: TH1D,
        df: pd.DataFrame,
        weight_dict: dict,
        sidebands: bool = True,
        weight_scale: float = 1.0,
    )-> dict:
        entries = len(df)
        sbcut = (df["mass_gggg"] <= 115.0) | (df["mass_gggg"] >= 135.0)

        for idx in range(0, entries):
            if (sbcut.iloc[idx] if sidebands else True):
                ibin = hist.FindBin(df[col].iloc[idx])
                weight_dict[ibin] += 1.0 * weight_scale

        return weight_dict

    print("Filling sparse bins for data...")
    weight_data = sparseBins1D(BDT_score_hist, data, weight_data, sidebands=True)

    print("Filling sparse bins for bkg...")
    weight_bkg = sparseBins1D(BDT_score_hist, bkg, weight_bkg, sidebands=True, weight_scale=sb_scaling_weight)

    print("Calculating ratio of data/bkg...")

    for subset in weight_ratio.keys():
        if weight_bkg[subset] > 0.0 and weight_data[subset] > 0.0:
            weight_ratio[subset] = weight_data[subset] / weight_bkg[subset]
        elif weight_bkg[subset] > 0.0 and weight_data[subset] == 0.0:
            weight_ratio[subset] = 1.0

    weights = []
    for idx in range(0, len(bkg)):
        ibin = BDT_score_hist.FindBin(bkg[col].iloc[idx])
        weights.append(weight_ratio[ibin])

    assert len(weights) == len(bkg), "{len(weights)} {len(bkg)}"
    reweight_df = pd.concat([bkg, pd.Series(weights, name="1dimreweight")], axis=1)

    return reweight_df


# REMOVED: findVal, getEntries, sparseBins (per-event ROOT loops, no longer used)


def loadNDimWeights(
    year: str
) -> pd.Series:
    """
    Loads N-dim reweight values for a given year. The reweighting process must be completed and the file must exist for this to work.
    """

    cwd = os.path.dirname(os.path.abspath(__file__))

    weights_path = f"{cwd}/../../scripts/bdtIO/BDT_{year}_input_weights.parquet"
    weights_path_data = f"{cwd}/../../scripts/bdtIO/BDT_{year}_input_weights_data.parquet"
    assert os.path.exists(weights_path)
    assert os.path.exists(weights_path_data)

    weights_df = pd.read_parquet(weights_path)
    weights_data_df = pd.read_parquet(weights_path_data)
    weights_df["Ndimreweight"] = (len(weights_data_df) / len(weights_df)) * weights_df["Ndimreweight"]
    print("Read weights from file (normalized to data).")

    return weights_df


# NEW: small helpers
def _flat_bin_index(df, edges_list, names):
    """ROOT-like bin number per variable (0=underflow, n+1=overflow), combined into one integer per event."""
    flat = np.zeros(len(df), dtype=np.int64)
    for name, edges in zip(names, edges_list):
        b = np.digitize(df[name].to_numpy(dtype=np.float64), edges)
        flat = flat * (len(edges) + 1) + b
    return flat


def _sideband(df):
    m = df["mass_gggg"].to_numpy()
    return (m <= 115.0) | (m >= 135.0)


def performNDimReweighting(
    b: str,
    d: str,
    year: str,
    eras: list,
    nBins: int = None,
    mvaCut: float = None,
    EB_only: bool = False,
    save_default: bool = False,
    default_weight: float = 0.1
) -> None:
    """
    Performs N-dim reweighting for event mixing from data/event mixing in sideband region (vectorized with numpy).
    """

    # NEW: read straight into DataFrames with uproot (replaces loadROOTSamples + CopyTree)
    def read_events(dirpath):
        files = [n for n in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath, n))]
        assert len(files) == 1, (dirpath, files)
        with uproot.open(os.path.join(dirpath, files[0])) as f:
            df = f["Events"].arrays(library="pd")
        if mvaCut is not None:
            assert type(mvaCut) is float, type(mvaCut)
            df = df[(df.pho3_mvaID >= mvaCut) & (df.pho4_mvaID >= mvaCut)]
        if EB_only:
            df = df[df.pho1_isScEtaEB & df.pho2_isScEtaEB & df.pho3_isScEtaEB & df.pho4_isScEtaEB]
        return df.reset_index(drop=True)

    print("Loading background samples:")
    evtmix_df = read_events(b)
    print("Loading data samples:")
    data_df = read_events(d)

    cwd = os.path.dirname(os.path.abspath(__file__))
    weights_path = f"{cwd}/../../scripts/bdtIO/BDT_{year}_input_weights.root"   # only used to build output names now

    nBins_default_2022 = 8
    nBins_default_2024 = 10
    if nBins is None:
        if "2022" in year:
            nBins = nBins_default_2022
        elif "2024" in year:
            nBins = nBins_default_2024
    print(f"Using {nBins}+2 bins for TH1s.")

    # Same variable selection as before, but store bin EDGES instead of building TH1Ds
    names = []
    edges_list = []
    branches, boundsList, _, _, _, titles = loadPlottingParameters()
    dims = ["LeadPs_pt", "SubleadPs_pt", "Ps_massDiff", "dR_aa_mass_gggg", "pho1_mvaID", "pho2_mvaID", "pho3_mvaID", "pho4_mvaID"]

    for idx, branch in enumerate(branches):
        if branch in dims:
            if "mvaID" in branch:
                if "pho3" in branch or "pho4" in branch:
                    edges = np.array([-1, -0.93, -0.45, 0.2, 1.0])
                else:
                    edges = np.array([-1.0, 1.0])
            else:
                edges = np.linspace(boundsList[idx][0], boundsList[idx][1], nBins + 1)
            names.append(branch)
            edges_list.append(edges)

    # NEW: one integer cell index per event (replaces findVal / FindBin / tuple(subset))
    idx_data = _flat_bin_index(data_df, edges_list, names)
    idx_bkg = _flat_bin_index(evtmix_df, edges_list, names)
    sb_data = _sideband(data_df)
    sb_bkg = _sideband(evtmix_df)

    # Scale s = N_data(SB) / N_evtmix(SB)   (replaces the getEntries loops)
    sb_scaling_weight = 1.0
    if Scales.useScaledCounts:
        print("Scaling counts to data/bkg ratio in SB region.")
        evtmix_entries, data_entries = len(evtmix_df), len(data_df)
        sb_evtmix_entries, sb_data_entries = int(sb_bkg.sum()), int(sb_data.sum())
        sr_evtmix_entries = evtmix_entries - sb_evtmix_entries
        sr_data_entries = data_entries - sb_data_entries

        print(f"\t\t {'Event Mixing': <12}   {'Data': <12}")
        print(f"{'Full region: ':<18} {evtmix_entries:<9} {data_entries:<9}")
        print(f"{'Signal region: ':<18} {sr_evtmix_entries:<9} {sr_data_entries:<9}")
        print(f"{'Sideband region: ':<18} {sb_evtmix_entries:<9} {sb_data_entries:<9}")

        sb_scaling_weight = sb_data_entries / sb_evtmix_entries

    # NEW: counts per occupied cell, sidebands only (replaces sparseBins + THnSparse)
    uniq_d, cnt_d = np.unique(idx_data[sb_data], return_counts=True)
    uniq_b, cnt_b = np.unique(idx_bkg[sb_bkg], return_counts=True)
    cnt_b = cnt_b * sb_scaling_weight

    # NEW: ratio per evtmix-occupied cell
    print("Calculating ratio of data/bkg...")
    ratio_b = np.full(len(uniq_b), default_weight)          # evtmix>0, data==0
    has_data = np.zeros(len(uniq_b), dtype=bool)
    if len(uniq_d):
        pos = np.clip(np.searchsorted(uniq_d, uniq_b), 0, len(uniq_d) - 1)
        has_data = uniq_d[pos] == uniq_b
        ratio_b[has_data] = cnt_d[pos[has_data]] / cnt_b[has_data]
    hist_counts = {"zeros": int(has_data.sum()), "no_data_yes_bkg": int((~has_data).sum())}
    print(hist_counts)

    # NEW: weight for every evtmix event (SR included); cells never seen in sideband evtmix -> 1.0
    weights = np.ones(len(evtmix_df), dtype=np.float32)
    if len(uniq_b):
        pos = np.clip(np.searchsorted(uniq_b, idx_bkg), 0, len(uniq_b) - 1)
        found = uniq_b[pos] == idx_bkg
        weights[found] = ratio_b[pos[found]]

    # REMOVED: THnSparse, fillObj, evtmix_tree, intermediate ROOT file, makeWeightPlots block,
    #          and the second re-read of the data/evtmix files with uproot (already in memory)

    # Save (paths and layout identical to before)
    if not save_default:
        data_path = weights_path.replace(".root", "_data.parquet")
        data_df.to_parquet(data_path)
    else:
        data_path = weights_path.replace(".root", f"_{str(default_weight).replace('.','p')}.root").replace("weights_", "weights_data_").replace(".root", ".parquet").replace("bdtIO",f"chi2_checks/{year}/")
        if ("2022" in year and nBins != nBins_default_2022) or (year == "2024" and nBins != nBins_default_2024):
            data_path = data_path.replace(".parquet", f"_nBins{nBins}.parquet")
        data_df.to_parquet(data_path)
    print(f"Saved data to {data_path}")

    # NEW: weights column first, then reset_index() exactly as in the original concat
    evtmix_plus_weights = evtmix_df.reset_index()
    evtmix_plus_weights.insert(0, "Ndimreweight", weights)
    if not save_default:
        evtmix_path = weights_path.replace(".root", ".parquet")
        evtmix_plus_weights.to_parquet(evtmix_path)
    else:
        evtmix_path = weights_path.replace(".root", f"_{str(default_weight).replace('.','p')}.root").replace(".root", ".parquet").replace("bdtIO",f"chi2_checks/{year}/")
        if ("2022" in year and nBins != nBins_default_2022) or (year == "2024" and nBins != nBins_default_2024):
            evtmix_path = evtmix_path.replace(".parquet", f"_nBins{nBins}.parquet")
        evtmix_plus_weights.to_parquet(evtmix_path)
    print(f"Saved data to {evtmix_path}")

    print(f"Total sum of weights: {sum(evtmix_plus_weights.Ndimreweight):.4f}")
    print("~~~~    Finished reweighting procedure!    ~~~~")
