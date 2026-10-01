import os
from ROOT import TCanvas, gStyle, TLegend, kBlack, kRed, kBlue, TLine, kDashed
from h4g_tools.utils.loading import pathSetup, loadSamples, fillHist
from h4g_tools.utils.plotting import plotFancy
from h4g_tools.utils.scales import Scales as ScalesCls, SideBand

if __name__ == "__main__":

    pathSetup("2024")
    Scales = ScalesCls("2024")

    # Look at 30, 35, and 40 GeV data in SB
    df30 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/background/30_GeV/", branches=["BDT_score", "mass_gggg"])
    df35 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/background/35_GeV/", branches=["BDT_score", "mass_gggg"])
    df40 = loadSamples("{cwd}/../../scripts/outputs/BDT/2024/ES_BDT_2024/background/40_GeV/", branches=["BDT_score", "mass_gggg"])

    binsScale = 200  # 200 for zoomed version, 100 for full range
    #bounds = [0.8976, 1.0]
    #bounds = [0.9000, 1.0]  # Integrals match, but the bins are shifted left by 0.005
    #bins = [bounds[0] + 0.005*i for i in range(20+1)]
    #print(bins, 1/binsScale)
    bounds = [0.0, 1.0]
    bins = None
    h30 = fillHist("BDT_score", df30, bounds, binsScale=binsScale, bins=bins, normalize=True, sb = True, preventOverFlow = True)
    h35 = fillHist("BDT_score", df35, bounds, binsScale=binsScale, bins=bins, normalize=True, sb = True, preventOverFlow = True)
    h40 = fillHist("BDT_score", df40, bounds, binsScale=binsScale, bins=bins, normalize=True, sb = True, preventOverFlow = True)

    h30.SetTitle("")
    h35.SetTitle("")
    h40.SetTitle("")
    
    # Canvas setup
    canvas = TCanvas("canvas", "canvas", 1000, 1000)
    canvas.SetLogy()

    # Legend setup
    legend = TLegend(0.55-0.35,0.65,0.87-0.35,0.83)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    max_order = 6

    # 30 GeV
    h30.SetLineWidth(2)
    h30.SetMarkerStyle(2)
    h30.SetMarkerSize(2)
    h30.SetMarkerColor(kBlack)
    h30.SetLineColor(kBlack)
    h30.GetYaxis().SetTitle("Events")
    h30.GetXaxis().SetTitle("BDT Score")
    h30.GetYaxis().SetRangeUser(10**-4, 10**-2)
    h30.GetXaxis().SetRangeUser(0.95, 1.0)

    # 35 GeV
    h35.SetLineWidth(2)
    h35.SetMarkerStyle(24)
    h35.SetMarkerSize(2)
    h35.SetMarkerColor(kRed)
    h35.SetLineColor(kRed)
    h35.GetYaxis().SetTitle("Events")
    h35.GetXaxis().SetTitle("BDT Score")
    h35.GetYaxis().SetRangeUser(10**-4, 10**-2)
    h35.GetXaxis().SetRangeUser(0.95, 1.0)

    # 40 GeV
    h40.SetLineWidth(2)
    h40.SetMarkerStyle(5)
    h40.SetMarkerSize(2)
    h40.SetMarkerColor(kBlue)
    h40.SetLineColor(kBlue)
    h40.GetYaxis().SetTitle("Events")
    h40.GetXaxis().SetTitle("BDT Score")
    h40.GetYaxis().SetRangeUser(10**-4, 10**-2)
    h40.GetXaxis().SetRangeUser(0.95, 1.0)

    legend.AddEntry(h30, "Event Mix (SB) for m_{a} = 30 GeV", "l")
    legend.AddEntry(h35, "Event Mix (SB) for m_{a} = 35 GeV", "l")
    legend.AddEntry(h40, "Event Mix (SB) for m_{a} = 40 GeV", "l")

    h30.Draw("hist")
    h40.Draw("hist same")
    h35.Draw("hist same")
    legend.Draw("same")

    plotFancy(canvas, "BDT Distribution Comparison", lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, sub="Sidebands only", Scales=Scales)

    # Save figure
    path = "plots/BDT/2024/compare_BDT_30_35_40_bkg.png"
    canvas.SaveAs(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path} to file.")

