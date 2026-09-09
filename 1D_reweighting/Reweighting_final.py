#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

numerator_file = (
    "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/files_to_analyze/2024C_noMassCut/skimmed_files/nominal/960b74ae-9d64-11f1-9c9c-4f69b8bcbeef_Events_0-24234.parquet"
)

denominator_file = (
    "/eos/user/a/arnaik/HiggsDNA_LCG/higgs-dna-2024/higgs_dna/h4g_mixedBkg_2022EE/nominal/mixed_events_output_postmix.parquet"
)

output_file = (
    "FallbackPt0001_pTgggg_1Dreweighted_v2.parquet"
)


xvar = "pT_gggg_vec"

weight_column = "pT_gggg_1d_weight"


# ============================================================
# DERIVED COLUMN: pT_gggg_vec
# This quantity is not present in the input files a priori, so it needs to be computed and added before it can be used as the reweighting variable.
# ============================================================

def compute_derived_columns(df):
    """Reproduce the MET / vector-pT calculation from the original ROOT-making script."""
    px_total = (df["pho1_pt"] * np.cos(df["pho1_phi"]) + df["pho2_pt"] * np.cos(df["pho2_phi"]) +
                df["pho3_pt"] * np.cos(df["pho3_phi"]) + df["pho4_pt"] * np.cos(df["pho4_phi"]))
    py_total = (df["pho1_pt"] * np.sin(df["pho1_phi"]) + df["pho2_pt"] * np.sin(df["pho2_phi"]) +
                df["pho3_pt"] * np.sin(df["pho3_phi"]) + df["pho4_pt"] * np.sin(df["pho4_phi"]))
    pz_total = (df["pho1_pt"] * np.sinh(df["pho1_eta"]) + df["pho2_pt"] * np.sinh(df["pho2_eta"]) +
                df["pho3_pt"] * np.sinh(df["pho3_eta"]) + df["pho4_pt"] * np.sinh(df["pho4_eta"]))
    df["pT_gggg_vec"] = np.sqrt(px_total**2 + py_total**2)
    return df


xedges = np.array([
    0, 10, 20, 30, 40, 50, 65, 80, 100, 125, 160, 250, 850
], dtype=float)

nbins_x = len(xedges) - 1


# -------------------------------------------------
print("Reading numerator...")

num = pd.read_parquet(
    numerator_file
)

print("Numerator events:", len(num))


print("Reading denominator...")

den = pd.read_parquet(
    denominator_file
)

print("Denominator events:", len(den))


# ---------------------------------------------------
# add the derived column to both dataframes
print()
print("Computing derived column:", xvar)

num = compute_derived_columns(num)
den = compute_derived_columns(den)

# ---------------------------------------------------

# 1D-histogram

num_hist, _ = np.histogram(
    num[xvar],
    bins=xedges
)

den_hist, _ = np.histogram(
    den[xvar],
    bins=xedges
)


# weights calculation

weights = np.zeros_like(num_hist, dtype=float)

mask = den_hist > 0

weights[mask] = num_hist[mask] / den_hist[mask]

# Fallback for any bin that ends up with a weight of exactly 0 -- this covers both num=0 & den>0 (a real ratio of 0) and den=0.
weights[weights == 0] = 0.0001


# checker to denominator 0 cases -----------------------

uncovered_bins = (num_hist > 0) & (den_hist == 0)

print()
print("============================================")
print("UNCOVERED NUMERATOR BINS")
print("============================================")

print(
    "Number of bins with numerator > 0 "
    "but denominator = 0:",
    np.sum(uncovered_bins)
)

print(
    "Numerator events in those bins:",
    num_hist[uncovered_bins].sum()
)

# ------------------------------------------------------
# print summary
print()
print("============================================")
print("1D REWEIGHTING")
print("============================================")

print("Number of X bins:", nbins_x)
print("Bin edges:", xedges)

print("Numerator events:", num_hist.sum())
print("Denominator events:", den_hist.sum())

print("Non-empty numerator bins:", np.count_nonzero(num_hist))
print("Non-empty denominator bins:", np.count_nonzero(den_hist))

print("Bins with denominator = 0:", np.sum(den_hist == 0))

print("Bins with numerator = 0:", np.sum(num_hist == 0))

print()
print("Total numerator / denominator:")
print(num_hist.sum() / den_hist.sum())

print()
print("Net sum of per-bin weights (sum of the weight-map array itself):")
print(weights.sum())

# 1D VALIDATION PLOTS
# ------------------------------------------------------------


def plot_1d(hist, title, ylabel, outname, log=False):
    fig, ax = plt.subplots(figsize=(6, 5))
    centers = 0.5 * (xedges[:-1] + xedges[1:])
    widths = np.diff(xedges)
    ax.bar(centers, hist, width=widths, align="center",
           edgecolor="black", linewidth=0.5)
    if log:
        ax.set_yscale("log")
        ax.set_ylim(bottom=1e-6)
    ax.set_xlabel(xvar)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outname, dpi=150)
    plt.close(fig)
    print(f"  -> wrote {outname}")


print()
print("============================================")
print("1D VALIDATION PLOTS")
print("============================================")

plot_1d(num_hist, "Numerator: Data", "Events",
        "numerator_1Dhist.png", log=True)
plot_1d(den_hist, "Denominator: Event Mixing (nc20)", "Events",
        "denominator_1Dhist.png", log=True)
plot_1d(weights, "1D Reweighting Map (Numerator / Denominator)", "Weight",
        "weight_map_1Dhist.png", log=False)

# applying weights to the mixed events-----

ix = np.digitize(den[xvar], xedges) - 1

# Handle values exactly on upper boundary
ix = np.clip(ix, 0, nbins_x - 1)

event_weights = weights[ix]


# adding a additional column to the mixed events file with the calculated weights

den[weight_column] = event_weights


# ============================================================
# summary
# ============================================================

print()
print("============================================")
print("WEIGHT STATISTICS")
print("============================================")

print("Minimum weight:", event_weights.min())
print("Maximum weight:", event_weights.max())
print("Mean weight:", event_weights.mean())
print("Sum of weights:", event_weights.sum())

print()
print("Target numerator events:", len(num))
print("Weighted denominator events:", event_weights.sum())

print()
print("Number of zero-weight events:", np.sum(event_weights == 0))
print("Number of non-zero-weight events:", np.sum(event_weights > 0))


print()
print("Writing output:")
print(output_file)

den.to_parquet(
    output_file,
    index=False
)

print()
print("DONE")
