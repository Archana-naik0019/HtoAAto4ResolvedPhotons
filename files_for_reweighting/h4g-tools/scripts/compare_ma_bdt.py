import pandas as pd
import xgboost
import os
import argparse
from h4g_tools.bdt.event_selection_bdt import reduce_dataframe, predict_BDT
from h4g_tools.utils.scales import Scales as ScalesCls
from h4g_tools.utils.loading import loadPerMass, loadSamples, fillHist
from h4g_tools.utils.plotting import plotFancy
from ROOT import TCanvas, TLegend, gStyle, kGreen, kMagenta, kBlue, kRed, kOrange


if __name__ == "__main__":

    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser")
    parser.add_argument("-i", "--inputs", type=str, required=True, help="Samples to load.")
    args = parser.parse_args()

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    year = args.inputs.split("/")[-3]
    Scales = ScalesCls(year) 

    stype = ""
    if "signal" in args.inputs:
        stype = "signal"
    elif "background" in args.inputs:
        stype = "background"
    elif "data" in args.inputs:
        stype = "data"

    # Load signal MC for each mass point
    print("Loading samples...")
    samples = loadPerMass(args.inputs)
    if stype == "background":
        data_samples = loadPerMass(os.path.join(*(args.inputs.split("/")[:-2] + ["data"])))

    # Histogram setup
    bins = 30
    hists = {}
    data_hists = {}
    for mass in samples.keys():
        if mass[:mass.find("GeV")+3] in ["15_GeV", "25_GeV", "35_GeV", "45_GeV", "60_GeV"]:
            hists.update({mass: fillHist("BDT_score", samples[mass][(samples[mass].mass_gggg < 115.0) | (samples[mass].mass_gggg > 135.0)], [0, 1], binsScale=bins, normalize=False)})
            if stype == "background":
                data_hists.update({mass: fillHist("BDT_score", data_samples[mass][(data_samples[mass].mass_gggg < 115.0) | (data_samples[mass].mass_gggg > 135.0)], [0, 1], binsScale=bins, normalize=False)})
                hists[mass].Scale(data_hists[mass].Integral() / hists[mass].Integral())

    # Canvas setup
    canvas = TCanvas("canvas", "canvas", 1000, 1000)
    canvas.SetLogy()

    # Legend setup
    legend = TLegend(0.65,0.60,0.87,0.8)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # Main plotting
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    colors = [kGreen+3, kMagenta, kOrange, kRed, kBlue]
    for (key, hist), color in zip(hists.items(), colors):
        hist.SetTitle("")
        hist.GetYaxis().SetTitleFont(42)
        hist.GetYaxis().SetTitle("Events")
        hist.GetXaxis().SetTitle("BDT Score")
        hist.GetXaxis().SetTitleOffset(1.1)
        hist.Sumw2()
        hist.SetLineWidth(2)
        hist.SetLineColor(color)
        
        if stype == "background":
            hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, 10**6 - 5*10**5)
        elif stype == "signal":
            hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, 1.01)

        legend.AddEntry(hist, key.replace("_"," "), "l")
        hist.Draw("hist same")

    legend.Draw()
    plotFancy(canvas, "BDT Score Comparison", lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, sub="Sideband Only", Scales=Scales, shiftLeft=0.05)

    canvas.Update()
    canvas.SaveAs(os.path.join(*(cwd, "plots/BDT", year, "BDT_pred_generic", f"{stype}_ma_comparisons.png")))
