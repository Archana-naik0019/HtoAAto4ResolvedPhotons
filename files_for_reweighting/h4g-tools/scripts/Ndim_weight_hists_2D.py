import os
import argparse
import pandas as pd
from ROOT import gPad, TCanvas, TLine, kBlack, kRed, gStyle, TH2D
from h4g_tools.utils.loading import pathSetup, loadSamples, fillHist
from h4g_tools.utils.plotting import plotFancy

"""
python3 Ndim_weight_hists_2D.py -i1 /cms/cephfs/data/store/user/castells/bdtIO/BDT_2022_input_weights.parquet -i2 /cms/cephfs/data/store/user/castells/bdtIO/BDT_2024_input_weights.parquet -y1 2022 -y2 2024
"""


def makeHist(
    path1: str,
    path2: str,
    outpath: str,
    title: str,
    year1: str,
    year2: str,
) -> None:
    """
    Makes histogram of Ndimreweight column for a given file.
    """

    sample_y1 = loadSamples(path1, branches=["Ndimreweight"])
    sample_y2 = loadSamples(path2, branches=["Ndimreweight"])

    bounds1 = [0.0, 23500]
    bounds2 = [0.0, 23500]
    binsScale = 1/250
    hist_y1 = fillHist("Ndimreweight", sample_y1, [0.0, 3.5], binsScale=25.0, normalize=False, noWeights=True)
    hist_y2 = fillHist("Ndimreweight", sample_y2, [0.0, 3.5], binsScale=25.0, normalize=False, noWeights=True)
    hist2D = TH2D("hist2D", "hist2D", int((bounds2[1]-bounds2[0])*binsScale), bounds2[0], bounds2[1], int((bounds1[1]-bounds1[0])*binsScale), bounds1[0], bounds1[1])

    for b in range(hist_y1.GetNbinsX()):
        hist2D.Fill(hist_y1.GetBinContent(b), hist_y2.GetBinContent(b))

    canvas = TCanvas("canvas", "cavnas", 1000, 1000)
    gStyle.SetOptStat(1)
    
    hist2D.SetStats(0)
    hist2D.SetMarkerStyle(8)
    hist2D.SetMarkerSize(0.5)
    hist2D.SetTitle("")
    hist2D.GetXaxis().SetTitle(f"{year1} N-dim Reweight (Data / Bkg)")
    hist2D.GetYaxis().SetTitle(f"{year2} N-dim Reweight (Data / Bkg)")
    hist2D.GetYaxis().SetTitleFont(42)
    hist2D.GetXaxis().SetTitleFont(42)
    hist2D.GetXaxis().SetTitleOffset(1.5)
    hist2D.GetYaxis().SetMaxDigits(3)
    hist2D.GetXaxis().SetMaxDigits(3)
    hist2D.Draw("COLZ")

    plotFancy(canvas, title, prelim=True, inPlot=True, colz=True)
    canvas.SaveAs(outpath)


if __name__ == "__main__":

    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser", conflict_handler="resolve")
    parser.add_argument("-i1", "--input1", required=True, help="Load first file or directory to produce N-dim weights histograms.", type=str)
    parser.add_argument("-i2", "--input2", required=True, help="Load second file or directory to produce N-dim weights histograms.", type=str)
    parser.add_argument("-y1", "--year1", required=True, help="Year of first sample.", type=str)
    parser.add_argument("-y2", "--year2", required=True, help="Year of second sample.", type=str)

    args = parser.parse_args()

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    samples_y1 = []
    samples_y2 = []
    if os.path.isfile(args.input1):
        assert "_data" not in args.input1 and ".parquet" in args.input1
        samples_y1.append(args.input1)
    if os.path.isfile(args.input2):
        assert "_data" not in args.input2 and ".parquet" in args.input2
        samples_y2.append(args.input2)
    if os.path.isdir(args.input1):
        for path in os.listdir(os.path.join(cwd, args.input1)):
            if "_data" not in path and ".parquet" in path and os.path.isfile(os.path.join(cwd, args.input1, path)):
                samples_y1.append(os.path.join(cwd, args.input1, path))
    if os.path.isdir(args.input2):
        for path in os.listdir(os.path.join(cwd, args.input2)):
            if "_data" not in path and ".parquet" in path and os.path.isfile(os.path.join(cwd, args.input2, path)):
                samples_y2.append(os.path.join(cwd, args.input2, path))

    for path1, path2, in zip(samples_y1, samples_y2):
        outfile = os.path.join(cwd, f"plots/BDT/compare_{args.year1}_{args.year2}_2D.png")
        makeHist(path1, path2, outfile, f"{args.year1} vs {args.year2} 8-dim RWing", args.year1, args.year2)
