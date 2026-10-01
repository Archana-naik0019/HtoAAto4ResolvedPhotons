import os
import json
import xgboost
import argparse
import pandas as pd
from h4g_tools.bdt.event_selection_bdt import reduce_dataframe, predict_BDT
from h4g_tools.utils.scales import Scales as ScalesCls, SideBand
from h4g_tools.utils.loading import loadPerMass, loadSamples, fillHist
from h4g_tools.utils.plotting import plotFancy, residuals_error_prop
from h4g_tools.bdt.reweight import reweight1Dim
from ROOT import TCanvas, TLegend, gStyle, kGreen, kMagenta, kBlue, kAzure, kRed, kOrange, kBlack, TPad, gPad, TLine, TPaveText



def plot(masses, sig_hists, xTitle, title, path, bdt_cuts): 
    # Canvas setup
    canvas = TCanvas("canvas", "canvas", 1000, 1000)
    
    canvas.SetLogy()

    # Legend setup
    # NOTE: Make more legend to split signal, bkg, and data into something easier to read...
    legend_sig = TLegend(0.2,0.7,0.43,0.87)
    legend_sig.SetTextFont(42)
    legend_sig.SetBorderSize(0)
    legend_sig.SetFillStyle(0)

    # Main plotting
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    colors = [kGreen+3, kMagenta, kRed, kOrange, kBlue]
    shapes = [23, 22, 20, 47, 21]
    for mass, color, shape in zip([masses[i] for i in range(len(masses)-1, -1, -1)], colors, shapes):
        sig_hist = sig_hists[mass]
        sig_hist.SetTitle("")
        sig_hist.GetYaxis().SetTitleFont(42)
        sig_hist.GetYaxis().SetTitle("Events")
        sig_hist.GetXaxis().SetTitleOffset(1.1)
        sig_hist.Sumw2()
        sig_hist.SetLineWidth(2)
        sig_hist.SetLineColor(color)
        sig_hist.GetXaxis().SetTitle(xTitle)
        
        min_order = -2
        max_order = 8
        sig_hist.GetYaxis().SetRangeUser(10**min_order + 10**(min_order-1), 10**max_order - 5*10**(max_order-1))

    # Force draw order
    for mass, color in zip([masses[i] for i in range(len(masses)-1, -1, -1)], colors):
        sig_hist = sig_hists[mass]
        sig_hist.Draw("hist same")

    # Reverse order of legend
    for mass in [masses[i] for i in range(len(masses)-1, -1, -1)]:
        legend_sig.AddEntry(sig_hists[mass], "Signal " + mass.replace("_"," "), "l")

    legend_sig.Draw("same")

    plotFancy(canvas, title, lumiTxt=Scales.lumiTxt, pad=False, prelim=True, inPlot=False, Scales=Scales, shiftLeft=0.05)

    def printTxt(lowX, lowY, text, align=12):
        subtxt = TPaveText(lowX, lowY+0.05, lowX+0.15, lowY+0.035, "NDC")
        subtxt.SetTextFont(42)
        subtxt.SetTextSize(0.028)
        subtxt.SetTextColor(1)
        subtxt.SetTextAlign(align)
        subtxt.SetFillStyle(0)
        subtxt.SetBorderSize(0)

        subtxt.AddText(text)
        subtxt.DrawClone("same")

    lowX, lowY = [0.7, 0.8]
    printTxt(lowX-0.2, lowY, f"Before BDT Cut", align=22)
    printTxt(lowX, lowY, f"After BDT Cut", align=22)
    for idx, m in enumerate([masses[i] for i in range(len(masses)-1, -1, -1)]):
        bdt_cut = bdt_cuts[m]["cut0"][0]
        printTxt(lowX, lowY-0.036*(idx+1), f"{m.replace('_',' ')}: {sig_hists[m].Integral(sig_hists[m].GetXaxis().FindBin(bdt_cut), sig_hists[m].GetXaxis().FindBin(1.0)): .0f}")
        printTxt(lowX-0.2, lowY-0.036*(idx+1), f"{m.replace('_',' ')}: {sig_hists[m].Integral(0, sig_hists[m].GetNbinsX()): .0f}")
    canvas.SaveAs(path)
    assert os.path.exists(path), "Saved plot does not exist."


if __name__ == "__main__":

    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser")
    parser.add_argument("-y", "--year", type=str, required=True, help="Year to load.")
    parser.add_argument("-xgb", "--xgb_name", type=str, required=True, help="XGBoost model name.")
    parser.add_argument("-bdt", "--bdt_path", required=True, type=str, help="Path to file containing efficiency of BDT cuts. Assumes 1 category.")
    args = parser.parse_args()

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    Scales = ScalesCls(args.year) 

    sig_path = f"{cwd}/../../scripts/outputs/BDT/{args.year}/{args.xgb_name}/signal/"
    bkg_path = f"{cwd}/../../scripts/outputs/BDT/{args.year}/{args.xgb_name}/background/"

    # Load signal MC for each mass point
    print("=== Loading samples (nominal only) ===")
    print("~~~ Signal samples ~~~")
    sig_samples = loadPerMass(sig_path, branches=["BDT_score"])
    print()
    print("~~~ Background samples ~~~")
    bkg_samples = loadPerMass(bkg_path, branches=["BDT_score"])
    print()

    print("Keeping all events only for signal BDT cut check.")
    print()

    with open(args.bdt_path, "r") as f:
        bdt_eff = json.load(f)

    # Histogram setup
    bins = 200  # Enough bins to cut on any optimized BDT value
    masses1 = ["55_GeV", "45_GeV", "35_GeV", "25_GeV", "15_GeV"]  # Do first half of masses together
    masses2 = ["60_GeV", "50_GeV", "40_GeV", "30_GeV", "20_GeV"]  # Do second half of masses together
    for idx, masses in enumerate([masses1, masses2]):
        sig_bdts = {}
        bkg_bdts = {}
        for mass in masses:
            if mass[:mass.find("GeV")+3] in masses:
                sig_bdts.update({mass: fillHist("BDT_score", sig_samples[mass], [0, 1], binsScale=bins, normalize=False)})
                bkg_bdts.update({mass: fillHist("BDT_score", bkg_samples[mass], [0, 1], binsScale=bins, normalize=False)})

        plot(masses, sig_bdts, "BDT Score", "BDT Cut Checks", os.path.join(*(cwd, "plots/BDT", args.year, "BDT_pred_generic", args.xgb_name, f"bdt_cut_checks{idx}.png")), bdt_cuts=bdt_eff)
        plot(masses, bkg_bdts, "BDT Score", "BDT Cut Checks", os.path.join(*(cwd, "plots/BDT", args.year, "BDT_pred_generic", args.xgb_name, f"bdt_cut_bkg_checks{idx}.png")), bdt_cuts=bdt_eff)


