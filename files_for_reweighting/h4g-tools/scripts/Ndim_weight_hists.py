import os
import argparse
import pandas as pd
from ROOT import gPad, TCanvas, TLine, kBlack, kRed, gStyle
from h4g_tools.utils.loading import pathSetup, loadSamples, fillHist
from h4g_tools.utils.plotting import plotFancy

def makeHist(
    inpath: str,
    outpath: str,
    title: str,
    era: str,
) -> None:
    """
    Makes histogram of Ndimreweight column for a given file.
    """

    sample = loadSamples(inpath, branches=["Ndimreweight"])
    hist = fillHist("Ndimreweight", sample, [0.0, 10.0], binsScale=25, normalize=False, noWeights=True)
    hist_small = fillHist("Ndimreweight", sample, [0.0, 3.5], binsScale=25, normalize=False, noWeights=True)

    canvas = TCanvas("canvas", "cavnas", 1000, 1000)
    gStyle.SetOptStat(1)
    
    for t, h in zip(["all", "small"], [hist, hist_small]):
        h.SetTitle("")
        h.GetYaxis().SetTitle("Events")
        h.GetXaxis().SetTitle("N-dim Reweight (Data / Bkg)")
        h.SetLineColor(kBlack)
        h.SetLineWidth(2)
        h.GetYaxis().SetRangeUser(0.0, h.GetMaximum()*1.1)
        h.Draw("hist")

        if t == "all":
            print("All events: ", h.Integral())
        elif t == "small":
            print("Subset: ", h.Integral(1, h.GetNbinsX()-2))
        plotFancy(canvas, title, prelim=True, inPlot=True)
        canvas.SaveAs(outpath if t == "all" else outpath.replace(".png", "_subset.png"))


if __name__ == "__main__":

    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser", conflict_handler="resolve")
    parser.add_argument("-i", "--inputs", required=True, help="Load a file or directory to produce N-dim weights histograms.", type=str)
    parser.add_argument("-dw", "--dw", required=True, help="Default weight used.", type=str)
    parser.add_argument("-y", "--year", required=True, help="Year of sample weights.", type=str)

    args = parser.parse_args()

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    preEE = []
    postEE = []
    y2022 = []
    y2024 = []
    if os.path.isfile(args.inputs):
        assert "_data" not in args.inputs and ".parquet" in args.inputs
        if "preEE" in args.year:
            preEE.append(args.inputs)
        elif "postEE" in args.year:
            postEE.append(args.inputs)
        elif "2022" in args.year:
            y2022.append(args.inputs)
        elif "2024" in args.year:
            y2024.append(args.inputs)
    elif os.path.isdir(args.inputs):
        for path in os.listdir(os.path.join(cwd, args.inputs)):
            if "_data" not in path and ".parquet" in path and os.path.isfile(os.path.join(cwd, args.inputs, path)):
                if "preEE" in args.year and "preEE" in path:
                    preEE.append(os.path.join(cwd, args.inputs, path))
                elif "postEE" in args.year and "postEE" in path:
                    postEE.append(os.path.join(cwd, args.inputs, path))
                elif "2022" in args.year:
                    y2022.append(os.path.join(cwd, args.inputs, path))
                elif "2024" in args.year:
                    y2024.append(os.path.join(cwd, args.inputs, path))

    if len(preEE) > 0:
        pathSetup("2022preEE")
    if len(postEE) > 0:
        pathSetup("2022postEE")
    if len(y2022) > 0:
        pathSetup("2022")
    if len(y2024) > 0:
        pathSetup("2024")

    for pathList, era in zip([preEE, postEE, y2022, y2024], ["2022preEE", "2022postEE", "2022", "2024"]):
        outfile_base = os.path.relpath(os.path.join(cwd, "plots/BDT", era, "Ndimreweights/"), cwd)
        for path in pathList:
            outfile = os.path.join(outfile_base, path.split("_")[-1].replace(".parquet", ""))
            makeHist(path, outfile + ".png", f"8-dim RWing with {args.dw} default weight", era)
