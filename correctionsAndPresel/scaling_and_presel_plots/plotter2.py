#!/usr/bin/env python3

import os
import glob

import awkward as ak
import numpy as np
import pandas as pd
import uproot

import matplotlib.pyplot as plt

##########################################################
# Inputs
##########################################################

INPUT_ROOT = "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/skimmedData2024Full"

INPUT_PARQUET = "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/correctionsAndPresel/bulk_trial1"

TREE = "Events"

##########################################################

os.makedirs("plots/normalized", exist_ok=True)
os.makedirs("plots/unnormalized", exist_ok=True)

variables = [
    "pho1_pt",
    "pho2_pt",
    "pho3_pt",
    "pho4_pt",
]

xmax = {
    "pho1_pt": 400,
    "pho2_pt": 200,
    "pho3_pt": 100,
    "pho4_pt": 50,
}
##########################################################
# Find common files
##########################################################

print("Finding common files...")

root_map = {}

for f in glob.glob(os.path.join(INPUT_ROOT, "*", "*.root")):
    key = os.path.splitext(os.path.basename(f))[0]
    root_map[key] = f

parquet_map = {}

for f in glob.glob(os.path.join(INPUT_PARQUET, "*", "*_photons.parquet")):
    key = os.path.basename(f).replace("_photons.parquet", "")
    parquet_map[key] = f

common_keys = sorted(set(root_map.keys()) & set(parquet_map.keys()))

print(f"ROOT files    : {len(root_map)}")
print(f"Parquet files : {len(parquet_map)}")
print(f"Common files  : {len(common_keys)}")

if len(common_keys) == 0:
    raise RuntimeError("No common files found.")

root_files = [root_map[k] for k in common_keys]
parquet_files = [parquet_map[k] for k in common_keys]

##########################################################
# Read ROOT
##########################################################

print("Reading ROOT files...")

root_arrays = {v: [] for v in variables}

for f in root_files:

    try:

        tree = uproot.open(f)[TREE]

        pt = tree["Photon_pt"].array(library="ak")

        # Sort photons by descending pT
        order = ak.argsort(pt, axis=1, ascending=False)
        pt = pt[order]

        # Keep only events with at least 4 photons
        mask = ak.num(pt) >= 4
        pt = pt[mask]

        root_arrays["pho1_pt"].append(ak.to_numpy(pt[:, 0]))
        root_arrays["pho2_pt"].append(ak.to_numpy(pt[:, 1]))
        root_arrays["pho3_pt"].append(ak.to_numpy(pt[:, 2]))
        root_arrays["pho4_pt"].append(ak.to_numpy(pt[:, 3]))

    except Exception as e:

        print(f"Skipping {f}")
        print(e)

root_data = {}

for v in variables:

    if len(root_arrays[v]) == 0:
        raise RuntimeError(f"No data found for {v}")

    root_data[v] = pd.Series(
        np.concatenate(root_arrays[v]),
        name=v,
    )

##########################################################
# Read parquet
##########################################################

print("Reading parquet files...")

parquet_arrays = {v: [] for v in variables}

for f in parquet_files:

    try:

        df = pd.read_parquet(f)

        for v in variables:
            parquet_arrays[v].append(df[v].dropna())

    except Exception as e:

        print(f"Skipping {f}")
        print(e)

parquet_data = {}

for v in variables:

    if len(parquet_arrays[v]) == 0:
        raise RuntimeError(f"No parquet data found for {v}")

    parquet_data[v] = pd.concat(
        parquet_arrays[v],
        ignore_index=True,
    )

##########################################################
# Plotting
##########################################################

print("Making plots...")

for var in variables:

    before = root_data[var].dropna()
    after = parquet_data[var].dropna()

    xmin = 0
    xmax_plot = xmax[var]

    #bins = np.linspace(xmin, xmax_plot, 60)
    bin_width = 3.0
    bins = np.arange(xmin, xmax_plot + bin_width, bin_width)

    ###################################################
    # Statistics
    ###################################################

    n_before = len(before)
    mean_before = before.mean()
    median_before = before.median()
    overflow_before = (before > xmax_plot).sum()

    n_after = len(after)
    mean_after = after.mean()
    median_after = after.median()
    overflow_after = (after > xmax_plot).sum()

    #######################################################
    # Unnormalized
    #######################################################

    plt.figure(figsize=(7, 6))

    plt.hist(
        before,
        bins=bins,
        histtype="step",
        linewidth=2,
        label="HLT & ≥4 photons",
    )

    plt.hist(
        after,
        bins=bins,
        histtype="step",
        linewidth=2,
        label=f"After Scale & Presel",
    )

    plt.xlim(0, xmax_plot)

    stats = (
        f"ROOT\n"
        f"Entries : {n_before:,}\n"
        f"Mean    : {mean_before:.2f}\n"
        f"Median  : {median_before:.2f}\n"
        f"Overflow: {overflow_before:,}\n\n"
        f"Parquet\n"
        f"Entries : {n_after:,}\n"
        f"Mean    : {mean_after:.2f}\n"
        f"Median  : {median_after:.2f}\n"
        f"Overflow: {overflow_after:,}"
    )
    
    plt.text(
        0.97,
        0.97,
        stats,
        transform=plt.gca().transAxes,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(facecolor="white", alpha=0.8),
    )

    #plt.xlabel(var)
    xlabel = {
        "pho1_pt": r"Photon 1 $p_T$ (GeV)",
        "pho2_pt": r"Photon 2 $p_T$ (GeV)",
        "pho3_pt": r"Photon 3 $p_T$ (GeV)",
        "pho4_pt": r"Photon 4 $p_T$ (GeV)",
    }
    
    plt.xlabel(xlabel[var])

    
    plt.ylabel("Events")
    plt.title(var)
    plt.legend()

    plt.tight_layout()

    plt.savefig(
        f"plots/unnormalized/{var}.png",
        dpi=200,
    )

    plt.close()

    #######################################################
    # Normalized
    #######################################################

    plt.figure(figsize=(7, 6))

    plt.hist(
        before,
        bins=bins,
        density=True,
        histtype="step",
        linewidth=2,
        label="HLT & ≥4 photons",
    )

    plt.hist(
        after,
        bins=bins,
        density=True,
        histtype="step",
        linewidth=2,
        label="After Scale & Presel",
    )

    plt.xlim(0, xmax_plot)

    stats = (
        f"ROOT\n"
        f"Entries : {n_before:,}\n"
        f"Mean    : {mean_before:.2f}\n"
        f"Median  : {median_before:.2f}\n"
        f"Overflow: {overflow_before:,}\n\n"
        f"Parquet\n"
        f"Entries : {n_after:,}\n"
        f"Mean    : {mean_after:.2f}\n"
        f"Median  : {median_after:.2f}\n"
        f"Overflow: {overflow_after:,}"
    )
    
    plt.text(
        0.97,
        0.97,
        stats,
        transform=plt.gca().transAxes,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(facecolor="white", alpha=0.8),
    )

    #plt.xlabel(var)
    xlabel = {
        "pho1_pt": r"Photon 1 $p_T$ (GeV)",
        "pho2_pt": r"Photon 2 $p_T$ (GeV)",
        "pho3_pt": r"Photon 3 $p_T$ (GeV)",
        "pho4_pt": r"Photon 4 $p_T$ (GeV)",
    }

    plt.xlabel(xlabel[var])

    
    plt.ylabel("Normalized")
    plt.title(var)
    plt.legend()

    plt.tight_layout()

    plt.savefig(
        f"plots/normalized/{var}.png",
        dpi=200,
    )

    plt.close()

print()
print("Done.")
print(f"Compared {len(common_keys)} common files.")
print("Plots written to:")
print("  plots/unnormalized/")
print("  plots/normalized/")