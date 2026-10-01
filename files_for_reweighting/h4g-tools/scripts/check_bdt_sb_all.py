import os
from ROOT import TCanvas, gStyle, TLegend, TLine, kDashed, kAzure, kBlack, kGreen, kMagenta, kYellow, kSpring, kOrange, kCyan, kRed, kViolet, kOrange, kBlue, kGray, SetOwnership
from h4g_tools.utils.loading import pathSetup, loadSamples, fillHist
from h4g_tools.utils.plotting import plotFancy
from h4g_tools.utils.scales import Scales as ScalesCls, SideBand

if __name__ == "__main__":

    pathSetup("2024")
    Scales = ScalesCls("2024")

    # Look at 30, 35, and 40 GeV data in SB
    df15 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/data/15_GeV/", branches=["BDT_score", "mass_gggg"])
    df20 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/data/20_GeV/", branches=["BDT_score", "mass_gggg"])
    df25 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/data/25_GeV/", branches=["BDT_score", "mass_gggg"])
    df30 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/data/30_GeV/", branches=["BDT_score", "mass_gggg"])
    df35 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/data/35_GeV/", branches=["BDT_score", "mass_gggg"])
    df40 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/data/40_GeV/", branches=["BDT_score", "mass_gggg"])
    df45 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/data/45_GeV/", branches=["BDT_score", "mass_gggg"])
    df50 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/data/50_GeV/", branches=["BDT_score", "mass_gggg"])
    df55 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/data/55_GeV/", branches=["BDT_score", "mass_gggg"])
    df60 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/data/60_GeV/", branches=["BDT_score", "mass_gggg"])

    binsScale = 30  # 200 for zoomed version, 100 for full range
    bounds = [0.0, 1.0]
    bins = None
    h15 = fillHist("BDT_score", df15, bounds, binsScale=binsScale, bins=bins, normalize=False, sb = True, preventOverFlow = True)
    h20 = fillHist("BDT_score", df20, bounds, binsScale=binsScale, bins=bins, normalize=False, sb = True, preventOverFlow = True)
    h25 = fillHist("BDT_score", df25, bounds, binsScale=binsScale, bins=bins, normalize=False, sb = True, preventOverFlow = True)
    h30 = fillHist("BDT_score", df30, bounds, binsScale=binsScale, bins=bins, normalize=False, sb = True, preventOverFlow = True)
    h35 = fillHist("BDT_score", df35, bounds, binsScale=binsScale, bins=bins, normalize=False, sb = True, preventOverFlow = True)
    h40 = fillHist("BDT_score", df40, bounds, binsScale=binsScale, bins=bins, normalize=False, sb = True, preventOverFlow = True)
    h45 = fillHist("BDT_score", df45, bounds, binsScale=binsScale, bins=bins, normalize=False, sb = True, preventOverFlow = True)
    h50 = fillHist("BDT_score", df50, bounds, binsScale=binsScale, bins=bins, normalize=False, sb = True, preventOverFlow = True)
    h55 = fillHist("BDT_score", df55, bounds, binsScale=binsScale, bins=bins, normalize=False, sb = True, preventOverFlow = True)
    h60 = fillHist("BDT_score", df60, bounds, binsScale=binsScale, bins=bins, normalize=False, sb = True, preventOverFlow = True)
 
    # Canvas setup
    canvas = TCanvas("canvas", "canvas", 1000, 1000)
    canvas.SetLogy()

    # Legend setup
    legend = TLegend(0.55,0.50,0.87,0.83)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    max_order = 6

    colors = [kBlack,kRed, kOrange, kMagenta, kViolet-5, kBlue, kAzure-4, kCyan+1, kGreen+3, kSpring-3]
    for h, c, m in zip([h15, h20, h25, h30, h35, h40, h45, h50, h55, h60], colors, [m for m in range(15,65,5)]):
        h.SetTitle("")
        h.SetMarkerStyle(8)
        #h.SetMarkerSize(2)
        h.SetMarkerColor(c)
        h.SetLineColor(c)
        h.GetYaxis().SetTitle("Events")
        h.GetXaxis().SetTitle("BDT Score")
        #h.GetYaxis().SetRangeUser(0.0, 7.5)
        #h.GetXaxis().SetRangeUser(0.9, 1.0)

        legend.AddEntry(h, f"Data (SB) for m_{{a}} = {m} GeV", "p")

        h.Draw("P E1 same")

    legend.Draw("same")

    plotFancy(canvas, "BDT Distribution Comparison", lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, sub="Sidebands only", Scales=Scales)

    # Save figure
    path = "plots/BDT/2024/compare_BDT_SB_all_masses.png"
    canvas.SaveAs(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path} to file.")

