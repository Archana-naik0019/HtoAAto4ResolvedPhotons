#!/usr/bin/env python3

import os
import glob

import numpy as np
import pandas as pd
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(1)

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

xlabel = {
    "pho1_pt": "Photon 1 p_{T} (GeV)",
    "pho2_pt": "Photon 2 p_{T} (GeV)",
    "pho3_pt": "Photon 3 p_{T} (GeV)",
    "pho4_pt": "Photon 4 p_{T} (GeV)",
}

bin_width = 1.0

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
# Read ROOT files with RDataFrame
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

root_arrays = {v: [] for v in variables}

for f in root_files:
    try:
        df = ROOT.RDataFrame(TREE, f)
        df = df.Filter("Photon_pt.size() >= 4", ">=4 photons")
        df = df.Define("pho_pt_sorted", "sortedPt(Photon_pt)")

        for i, v in enumerate(variables):
            df = df.Define(v, f"pho_pt_sorted[{i}]")

        cols = df.AsNumpy(variables)

        for v in variables:
            root_arrays[v].append(np.asarray(cols[v], dtype=np.float64))

    except Exception as e:
        print(f"Skipping {f}")
        print(e)

root_data = {}
for v in variables:
    if len(root_arrays[v]) == 0:
        raise RuntimeError(f"No data found for {v}")
    root_data[v] = np.concatenate(root_arrays[v])

##########################################################
# Read parquet
##########################################################

print("Reading parquet files...")

parquet_arrays = {v: [] for v in variables}

for f in parquet_files:
    try:
        df = pd.read_parquet(f)
        for v in variables:
            parquet_arrays[v].append(df[v].dropna().to_numpy(dtype=np.float64))
    except Exception as e:
        print(f"Skipping {f}")
        print(e)

parquet_data = {}
for v in variables:
    if len(parquet_arrays[v]) == 0:
        raise RuntimeError(f"No parquet data found for {v}")
    parquet_data[v] = np.concatenate(parquet_arrays[v])


##########################################################
# Helpers
##########################################################

def make_bins(xmax_plot):
    """Same edges as np.arange(0, xmax_plot + bin_width, bin_width)."""
    edges = np.arange(0, xmax_plot + bin_width, bin_width)
    return edges


def fill_hist(hist, arr):
    arr = np.ascontiguousarray(arr, dtype=np.float64)
    w = np.ones_like(arr)
    hist.FillN(len(arr), arr, w)


def make_stats_text(n, mean, median, overflow, label):
    return (
        f"#splitline{{{label}}}{{"
        f"Entries: {n:,}  "
        f"Mean: {mean:.2f}  "
        f"Median: {median:.2f}  "
        f"Overflow: {overflow:,}}}"
    )


def draw_canvas(var, h_before, h_after, stats_before, stats_after, ylabel, outpath):

    c = ROOT.TCanvas(f"c_{var}_{outpath}", var, 750, 750)

    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.13)
    pad1.Draw()

    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.3)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.13)
    pad2.SetGridy()
    pad2.Draw()

    # ---------------- Upper pad ----------------
    pad1.cd()

    h_before.SetLineColor(ROOT.kBlue + 1)
    h_before.SetLineWidth(2)
    h_after.SetLineColor(ROOT.kRed + 1)
    h_after.SetLineWidth(2)

    ymax = max(h_before.GetMaximum(), h_after.GetMaximum())
    h_before.SetMaximum(ymax * 1.4)
    h_before.SetMinimum(0)

    h_before.GetYaxis().SetTitle(ylabel)
    h_before.GetYaxis().SetTitleSize(0.05)
    h_before.GetYaxis().SetLabelSize(0.04)
    h_before.GetXaxis().SetLabelSize(0)
    h_before.SetTitle(var)

    h_before.Draw("HIST E")
    h_after.Draw("HIST E SAME")

    legend = ROOT.TLegend(0.15, 0.75, 0.45, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(h_before, "HLT & #geq4 photons", "l")
    legend.AddEntry(h_after, "After Scale & Presel", "l")
    legend.Draw()


    # Put stat box in the right margin
    pave = ROOT.TPaveText(0.60, 0.55, 0.99, 0.95, "NDC")
    pave.SetFillColor(ROOT.kWhite)
    pave.SetFillStyle(1001)
    pave.SetBorderSize(1)
    pave.SetTextAlign(12)
    pave.SetTextSize(0.03)
    pave.AddText(stats_before)
    pave.AddText(stats_after)
    pave.Draw()

    # ---------------- Lower pad: ratio ----------------
    pad2.cd()

    h_ratio = h_before.Clone(f"ratio_{var}_{outpath}")
    h_ratio.Divide(h_after)

    h_ratio.SetLineColor(ROOT.kBlack)
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetMarkerSize(0.7)
    h_ratio.SetTitle("")

    h_ratio.GetYaxis().SetTitle("Skimmed/Scale&Presel")
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetTitleSize(0.11)
    h_ratio.GetYaxis().SetTitleOffset(0.5)
    h_ratio.GetYaxis().SetLabelSize(0.09)
    h_ratio.SetMinimum(0.0)
    h_ratio.SetMaximum(2.0)

    h_ratio.GetXaxis().SetTitle(xlabel[var])
    h_ratio.GetXaxis().SetTitleSize(0.12)
    h_ratio.GetXaxis().SetLabelSize(0.10)
    h_ratio.GetXaxis().SetTitleOffset(1.15)

    h_ratio.Draw("EP")

    line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1.0, h_ratio.GetXaxis().GetXmax(), 1.0)
    line.SetLineColor(ROOT.kGray + 2)
    line.SetLineStyle(2)
    line.Draw()

    c.SaveAs(outpath)


##########################################################
# Plotting
##########################################################

print("Making plots...")

# Keep Python references alive for the duration of the loop iteration
_keepalive = []

for var in variables:

    before = root_data[var]
    after = parquet_data[var]

    xmax_plot = xmax[var]
    edges = make_bins(xmax_plot)
    nbins = len(edges) - 1

    n_before = len(before)
    mean_before = before.mean()
    median_before = np.median(before)
    overflow_before = int((before > xmax_plot).sum())

    n_after = len(after)
    mean_after = after.mean()
    median_after = np.median(after)
    overflow_after = int((after > xmax_plot).sum())

    stats_before = make_stats_text(n_before, mean_before, median_before, overflow_before, "Skimmed")
    stats_after = make_stats_text(n_after, mean_after, median_after, overflow_after, "Scale & Presel")

    #######################################################
    # Unnormalized
    #######################################################

    h_before_raw = ROOT.TH1F(f"h_before_raw_{var}", var, nbins, edges)
    h_after_raw = ROOT.TH1F(f"h_after_raw_{var}", var, nbins, edges)

    fill_hist(h_before_raw, before)
    fill_hist(h_after_raw, after)

    draw_canvas(
        var,
        h_before_raw,
        h_after_raw,
        stats_before,
        stats_after,
        "Events",
        f"plots/unnormalized/{var}.png",
    )

    #######################################################
    # Normalized (matches matplotlib density=True)
    #######################################################

    h_before_norm = h_before_raw.Clone(f"h_before_norm_{var}")
    h_after_norm = h_after_raw.Clone(f"h_after_norm_{var}")

    int_before = h_before_norm.Integral()
    int_after = h_after_norm.Integral()

    if int_before > 0:
        h_before_norm.Scale(1.0 / (int_before * bin_width))
    if int_after > 0:
        h_after_norm.Scale(1.0 / (int_after * bin_width))

    draw_canvas(
        var,
        h_before_norm,
        h_after_norm,
        stats_before,
        stats_after,
        "Normalized",
        f"plots/normalized/{var}.png",
    )

    _keepalive.extend(
        [h_before_raw, h_after_raw, h_before_norm, h_after_norm]
    )

print()
print("Done.")
print(f"Compared {len(common_keys)} common files.")
print("Plots written to:")
print("  plots/unnormalized/")
print("  plots/normalized/")