#!/usr/bin/env python3

import os
import glob

import numpy as np
import pandas as pd
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

##########################################################
# Inputs
##########################################################

INPUT_ROOT = "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/skimmedData2024Full"

INPUT_PARQUET = "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/correctionsAndPresel/bulk_trial1"

TREE = "Events"

# >>> CHECK / ADJUST THESE <<<
# Name of the event-id branch in the ROOT tree, and the column
# name in the parquet file. Change if they don't match your files.
EVENT_ID_BRANCH_ROOT = "event"
EVENT_ID_COL_PARQUET = "event"

VAR = "pho1_pt"

# Range/binning of the scale-factor histogram
SF_MIN = 0.5
SF_MAX = 1.5
SF_BINS = 100

##########################################################

os.makedirs("plots/scale_factor", exist_ok=True)

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

##########################################################
# Read ROOT: event id + leading photon pt
##########################################################

print("Reading ROOT files...")

ROOT.gInterpreter.Declare(
    """
    std::vector<float> sortedPt(const ROOT::VecOps::RVec<float>& pt) {
        std::vector<float> v(pt.begin(), pt.end());
        std::sort(v.begin(), v.end(), std::greater<float>());
        return v;
    }
    """
)

root_data = {}  # key -> DataFrame(event, pho1_pt_root)

for key in common_keys:
    f = root_map[key]
    try:
        df = ROOT.RDataFrame(TREE, f)
        df = df.Filter("Photon_pt.size() >= 4", ">=4 photons")
        df = df.Define("pho_pt_sorted", "sortedPt(Photon_pt)")
        df = df.Define("pho1_pt", "pho_pt_sorted[0]")

        cols = df.AsNumpy([EVENT_ID_BRANCH_ROOT, "pho1_pt"])

        root_data[key] = pd.DataFrame(
            {
                "event": np.asarray(cols[EVENT_ID_BRANCH_ROOT]).astype(np.int64),
                "pho1_pt_root": np.asarray(cols["pho1_pt"], dtype=np.float64),
            }
        )

    except Exception as e:
        print(f"Skipping {f}")
        print(e)

##########################################################
# Read parquet: event id + pho1_pt
##########################################################

print("Reading parquet files...")

parquet_data = {}  # key -> DataFrame(event, pho1_pt_parquet)

for key in common_keys:
    f = parquet_map[key]
    try:
        df = pd.read_parquet(f, columns=[EVENT_ID_COL_PARQUET, VAR])
        df = df.dropna()

        parquet_data[key] = pd.DataFrame(
            {
                "event": df[EVENT_ID_COL_PARQUET].to_numpy().astype(np.int64),
                "pho1_pt_parquet": df[VAR].to_numpy(dtype=np.float64),
            }
        )

    except Exception as e:
        print(f"Skipping {f}")
        print(e)

##########################################################
# Match by event id (per file, to avoid cross-file id collisions)
# and compute the scale factor
##########################################################

print("Matching events and computing scale factors...")

sf_list = []
n_matched_total = 0
n_root_total = 0

for key in common_keys:
    if key not in root_data or key not in parquet_data:
        continue

    merged = pd.merge(
        root_data[key],
        parquet_data[key],
        on="event",
        how="inner",
        validate="one_to_one",
    )

    n_root_total += len(root_data[key])
    n_matched_total += len(merged)

    # guard against division by zero / bad values
    merged = merged[merged["pho1_pt_root"] > 0]

    sf = merged["pho1_pt_parquet"] / merged["pho1_pt_root"]
    sf_list.append(sf.to_numpy(dtype=np.float64))

if len(sf_list) == 0:
    raise RuntimeError("No matched events found — check EVENT_ID_BRANCH_ROOT / EVENT_ID_COL_PARQUET.")

scale_factors = np.concatenate(sf_list)

print(f"ROOT events (>=4 photons) : {n_root_total:,}")
print(f"Matched events            : {n_matched_total:,}")
print(f"Unmatched / dropped       : {n_root_total - n_matched_total:,}")

##########################################################
# Plot
##########################################################

print("Making plot...")

h = ROOT.TH1F("h_sf_pho1_pt", "pho1_pt Scale Factor", SF_BINS, SF_MIN, SF_MAX)

arr = np.ascontiguousarray(scale_factors, dtype=np.float64)
w = np.ones_like(arr)
h.FillN(len(arr), arr, w)

n_entries = len(scale_factors)
mean_sf = scale_factors.mean()
median_sf = np.median(scale_factors)
std_sf = scale_factors.std()
underflow = int((scale_factors < SF_MIN).sum())
overflow = int((scale_factors > SF_MAX).sum())

c = ROOT.TCanvas("c_sf_pho1_pt", "pho1_pt scale factor", 750, 650)

h.SetLineColor(ROOT.kBlue + 1)
h.SetLineWidth(2)
h.GetXaxis().SetTitle("pho1_pt scale factor  (parquet / ROOT)")
h.GetYaxis().SetTitle("Events")
h.GetXaxis().SetTitleSize(0.045)
h.GetYaxis().SetTitleSize(0.045)
h.Draw("HIST E")

pave = ROOT.TPaveText(0.62, 0.68, 0.90, 0.88, "NDC")
pave.SetFillColor(ROOT.kWhite)
pave.SetFillStyle(1001)
pave.SetBorderSize(1)
pave.SetTextAlign(12)
pave.SetTextSize(0.032)
pave.AddText(f"Entries : {n_entries:,}")
pave.AddText(f"Mean    : {mean_sf:.4f}")
pave.AddText(f"Median  : {median_sf:.4f}")
pave.AddText(f"Std Dev : {std_sf:.4f}")
pave.AddText(f"Under/Overflow : {underflow:,} / {overflow:,}")
pave.Draw()

c.SaveAs("plots/scale_factor/pho1_pt_scale_factor.png")

print()
print("Done.")
print(f"Compared {len(common_keys)} common files.")
print("Plot written to: plots/scale_factor/pho1_pt_scale_factor.png")
