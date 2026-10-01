import os
from random import random
import argparse
from h4g_tools.utils.scales import Scales as ScalesCls
from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist
from h4g_tools.utils.plotting import plotFancy, residuals_error_prop
from ROOT import TCanvas, TPad, TLegend, kAzure, kBlack, kGreen, kMagenta, kYellow, kSpring, kOrange, kCyan, kRed, kViolet, kOrange, kBlue, kGray, SetOwnership, gPad
import gc


def plotInter(
    hists: dict,
    title: str,
    outpath: str,
) -> None:
 
    """
    pad1.SetBottomMargin(0.0001)
    pad1.SetBorderMode(0)
    pad1.SetLogy()
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(0.3)
    pad2.SetBorderMode(0)

    pad1.Draw()
    pad2.Draw()
    """
    canvas.SetLogy()

    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    colors = [kAzure-4, kRed]
    # Histogram settings
    hists[0].SetLineWidth(3)
    hists[1].SetLineWidth(3)
    hists[0].SetLineColor(colors[0])
    hists[1].SetLineColor(colors[1])
    hists[0].SetFillColor(colors[0])
    hists[0].SetTitle("")
    hists[0].GetXaxis().SetTitle(title)
    hists[0].GetYaxis().SetTitle("Events")
    hists[0].GetYaxis().SetTitleOffset(1.5)
    hists[0].GetXaxis().SetTitleOffset(1.2)
    hists[0].GetYaxis().SetRangeUser(10**-3 + 10**-4, ((10**6) - (5*10**5))*10)
    hists[1].SetTitle("")
    hists[1].GetXaxis().SetTitle(title)
    hists[1].GetYaxis().SetTitle("Events")
    hists[1].GetYaxis().SetTitleOffset(1.5)
    hists[1].GetXaxis().SetTitleOffset(1.2)
    hists[1].GetYaxis().SetRangeUser(10**-3 + 10**-4, ((10**6) - (5*10**5))*10)

    # Reset canvas if first in branch
    hists[0].Draw("hist")
    hists[1].Draw("same hist")

    # Clone background to include statistical errors
    background_error = hists[0].Clone("background_error")
    background_error.SetLineColor(kGray+3)
    background_error.SetFillColor(kGray+3)
    background_error.SetFillStyle(3008)
    background_error.Draw("E2 same")

    legend.AddEntry(hists[0], "Signal MC Combined", "l")
    legend.AddEntry(hists[1], "Event Mixing", "f")
    legend.AddEntry(background_error, "Statistical Uncertainty", "f")    
    legend.Draw()

    """
    # Residuals histogram setup
    residuals = residuals_error_prop(hists[1], hists[0], title, autoScale=False)

    # Plotting residuals on pad2
    pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
    pad2.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    gPad.SetFrameBorderMode(0)
    gPad.SetTickx(1)
    gPad.SetTicky(1)
    residuals.SetTitle("")
    residuals.GetXaxis().SetTitle(title)
    residuals.GetYaxis().SetTitle("Signal / Bkg")
    residuals.GetXaxis().SetTitleSize(9.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
    residuals.GetYaxis().SetTitleSize(8.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
    residuals.GetXaxis().SetTitleOffset(1.1)
    residuals.GetYaxis().SetTitleOffset(0.4)
    residuals.GetYaxis().SetNdivisions(6)
    residuals.GetYaxis().CenterTitle(True)
    residuals.GetXaxis().SetLabelSize(2.5*residuals.GetXaxis().GetLabelSize())
    residuals.GetYaxis().SetLabelSize(2.5*residuals.GetYaxis().GetLabelSize())
    residuals.SetMarkerSize(1)
    residuals.SetMarkerStyle(20)
    residuals.SetMarkerColor(kBlack)
    residuals.SetLineColor(kBlack)
    residuals.Draw("P")

    # Line at y=1
    max_edge = residuals.GetBinCenter(residuals.GetNbinsX())
    min_edge = residuals.GetBinCenter(1)
    line = TLine(min_edge, 1.0, max_edge, 1.0)
    line.SetLineStyle(7)
    line.SetLineColor(kBlack)
    line.Draw("same")

    pad1.Modified()
    pad2.Modified()
    """
    canvas.Update()
    plotFancy(
        #pad1,
        canvas,
        title,
        lumiTxt=Scales.defaultLumiTxt,
        prelim=True,
        inPlot=False,
        Scales=Scales,
        #pad=True,
        sub="Sidebands only",
    )

    # Save canvas
    canvas.SaveAs(f"{cwd}/plots/BDT/{Scales.yr}/sig_mix/{outpath}")
    assert os.path.exists(f"{cwd}/plots/BDT/{Scales.yr}/sig_mix/{outpath}")

    # Clear legend for next branch
    legend.Clear()
    canvas.Clear()


def plotMass(
    hist: dict,
    title: str,
    outpath: str,
    log: bool
) -> None:

    canvasX = 1000
    canvasY = 1000
    canvas = TCanvas(str(random()), "canvas", canvasX, canvasY)
    legend = TLegend(0.175,0.575,0.55,0.875)

    if log:
        canvas.SetLogy(True)
    else:
        canvas.SetLogy(False)

    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    colors = [kAzure-4, kBlack, kGreen+3, kMagenta, kYellow, kSpring-3, kOrange, kCyan+1, kRed, kViolet-5, kOrange+4, kBlue]
    for idx, (mass, color) in enumerate(zip(hist.keys(), colors)):
        # Histogram settings
        hist[mass].SetLineWidth(3)
        hist[mass].SetLineColor(color)
        hist[mass].SetTitle("")
        hist[mass].GetXaxis().SetTitle(title)
        hist[mass].GetYaxis().SetTitle("Events / 1 GeV")
        hist[mass].GetYaxis().SetTitleOffset(1.5 if log else 2.0)
        hist[mass].GetXaxis().SetTitleOffset(1.2)

        # Reset canvas if first in branch
        if idx == 0:
            hist[mass].Draw("hist")
        else:
            hist[mass].Draw("same hist")

        legend.AddEntry(hist[mass], f"m(a) = {mass[7:9]} GeV", "l")

        if not log:
            hist[mass].GetYaxis().SetRangeUser(0.0, 0.04)
        else:
            hist[mass].GetYaxis().SetRangeUser(10**-4, 1)
        
    legend.Draw("same")
    plotFancy(
        canvas,
        title,
        lumiTxt=Scales.defaultLumiTxt,
        prelim=True,
        inPlot=False,
        Scales=Scales,
        simulation=True
    )

    # Save canvas
    canvas.SaveAs(f"{cwd}/plots/standard/{Scales.yr}/{outpath}")
    assert os.path.exists(f"{cwd}/plots/standard/{Scales.yr}/{outpath}")

    # Clear legend for next branch
    legend.Clear()
    canvas.Clear()

if __name__ == "__main__":

    gc.disable()

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", type=str, required=True, help="Signal samples to plot per mass.")
    parser.add_argument("-bkg", "--bkg", type=str, required=True, help="Background samples to plot.")
    parser.add_argument("-d", "--data", type=str, required=True, help="Data samples for normalization.")
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    args = parser.parse_args()
    
    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    Scales = ScalesCls(args.year) 
    masses = [x for x in range(15, 65, 5)]
    if args.year == "2022preEE":
        eras={
            "data": ["Run2022C", "Run2022D"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_preEE" for m in masses]
        }
    elif args.year == "2022postEE":
        eras={
            "data": ["Run2022E", "Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_postEE" for m in masses]
        }
    elif args.year == "2022":
        eras={
            "data": ["Run2022C", "Run2022D", "Run2022E", "Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D", "RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preEE", "postEE"]]
        }
    elif args.year == "2023preBPix":
        eras={
            "data": ["Run2023C"],
            "bkg": ["RunEvtMix2023C"],
            "signal": [f"Signal_{m}_GeV_preBPix" for m in masses]
        }
    elif args.year == "2023postBPix":
        eras={
            "data": ["Run2023D"],
            "bkg": ["RunEvtMix2023D"],
            "signal": [f"Signal_{m}_GeV_postBPix" for m in masses]
        }
    elif args.year == "2023":
        eras={
            "data": ["Run2023C", "Run2023D"],
            "bkg": ["RunEvtMix2023C", "RunEvtMix2023D"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preBPix", "postBPix"]]
        }
    elif args.year == "2024":
        eras={
            "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
            "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
            "signal": [f"Signal_{m}_GeV" for m in masses]
        }
    else:
        assert False, "Must specify year in Scales class!"

    # Handle path setup after checking year
    pathSetup(args.year)

    # Load samples
    print("Loading signal samples...")
    samples = loadPerMass(args.inputs, eras = eras, year = args.year)

    """
    print("Loading background samples...")
    bkg_samples = loadSamples(args.bkg)

    print("Loading data samples...")
    data_samples = loadSamples(args.data)
    """

    # Fill histograms
    hists_m4g = {}
    hists_ma1 = {}
    hists_ma2 = {}
    sorted_samples = sorted(samples.keys())
    for idx, mass in enumerate(sorted_samples):
        sample = samples[mass]
        #hists_m4g[mass] = fillHist("mass_gggg", sample, [110.0, 180.0], binsScale=1, name=f"mass_gggg_{mass[7:13]}", normalize=False)
        #hists_ma1[mass] = fillHist("LeadPs_mass", sample, [10.0, 90.0], binsScale=1, name=f"Ps1mass_{mass[7:13]}", normalize=False)
        #hists_ma2[mass] = fillHist("SubleadPs_mass", sample, [10.0, 90.0], binsScale=1, name=f"Ps2mass_{mass[7:13]}", normalize=False)

        hists_m4g[mass] = fillHist("mass_gggg", sample, [110.0, 180.0], binsScale=1, name=f"mass_gggg_{mass[7:13]}", normalize=True, preventOverFlow=True)
        hists_ma1[mass] = fillHist("LeadPs_mass", sample, [0.0, 100.0], binsScale=1, name=f"Ps1mass_{mass[7:13]}", normalize=True, preventOverFlow=True)
        hists_ma2[mass] = fillHist("SubleadPs_mass", sample, [0.0, 100.0], binsScale=1, name=f"Ps2mass_{mass[7:13]}", normalize=True, preventOverFlow=True)

        #scale = Scales.lumi_fb * 1.0 / Scales.sumw[Scales.yr][mass.replace("Signal_","")]
        #scale = Scales.lumi_fb * 1.0 / sum(sample.weight)  # stupid original scaling
        #hists_m4g[mass].Scale(scale)
        #hists_ma1[mass].Scale(scale)
        #hists_ma2[mass].Scale(scale)


        """
        if idx < 1:
            sig_a1_inter = fillHist("LeadPs_interMass", sample, [-0.6, 0.6], binsScale=1.0/0.04, name=f"sig_LeadPs_interMass", normalize=False)
            sig_a2_inter = fillHist("SubleadPs_interMass", sample, [-0.6, 0.6], binsScale=1.0/0.04, name=f"sig_SubleadPs_interMass", normalize=False)
        else:
            sig_a1_inter.Add(fillHist("LeadPs_interMass", sample, [-0.6, 0.6], binsScale=1.0/0.04, name=f"{idx}sig_LeadPs_interMass", normalize=False))
            sig_a2_inter.Add(fillHist("SubleadPs_interMass", sample, [-0.6, 0.6], binsScale=1.0/0.04, name=f"{idx}sig_SubleadPs_interMass", normalize=False))
        """
    
    """
    bkg_hists = {}
    bkg_hists_rw = {}
    data_hists = {}
    for branch in ["LeadPs_interMass", "SubleadPs_interMass"]:
        bkg_hists[branch] = fillHist(branch, bkg_samples, [-0.6, 0.6], binsScale=1.0/0.04, name=f"bkg_{branch}", normalize=False, sb = True)
        bkg_hists_rw[branch] = fillHist(branch, bkg_samples, [-0.6, 0.6], binsScale=1.0/0.04, name=f"bkg_{branch}", normalize=False, sb = True, customWeights="Ndimreweight")
        data_hists[branch] = fillHist(branch, data_samples, [-0.6, 0.6], binsScale=1.0/0.04, name=f"data_{branch}", normalize=False, sb = True)
        bkg_hists[branch].Scale(data_hists[branch].Integral() / bkg_hists[branch].Integral())
        bkg_hists_rw[branch].Scale(data_hists[branch].Integral() / bkg_hists_rw[branch].Integral())

    sig_a1_inter.Scale(data_hists["LeadPs_interMass"].Integral() / sig_a1_inter.Integral())
    sig_a2_inter.Scale(data_hists["SubleadPs_interMass"].Integral() / sig_a2_inter.Integral())
    """

    plotMass(hists_m4g, "m_{4#gamma}", "m4g_all.pdf", log=True)
    plotMass(hists_m4g, "m_{4#gamma}", "m4g_all_linear.pdf", log=False)

    plotMass(hists_ma1, "m_{a1}", "ma1_all.pdf", log=True)
    plotMass(hists_ma1, "m_{a1}", "ma1_all_linear.pdf", log=False)

    plotMass(hists_ma2, "m_{a2}", "ma2_all.pdf", log=True)
    plotMass(hists_ma2, "m_{a2}", "ma2_all_linear.pdf", log=False)

    """
    canvas = TCanvas(f"canvas", "canvas", canvasX, canvasY)
    legend = TLegend(0.62,0.68, 0.87,0.83)

    plotInter([bkg_hists["SubleadPs_interMass"], sig_a2_inter], "(m_{a2} - m_{Hyp}) / m_{#gamma#gamma#gamma#gamma}", "SubleadPs_interMass_compare_pre_reweight.pdf")
    plotInter([bkg_hists["LeadPs_interMass"], sig_a1_inter], "(m_{a1} - m_{Hyp}) / m_{#gamma#gamma#gamma#gamma}", "LeadPs_interMass_compare_pre_reweight.pdf")
    plotInter([bkg_hists_rw["LeadPs_interMass"], sig_a1_inter], "(m_{a1} - m_{Hyp}) / m_{#gamma#gamma#gamma#gamma}", "LeadPs_interMass_compare_reweight.pdf")
    plotInter([bkg_hists_rw["SubleadPs_interMass"], sig_a2_inter], "(m_{a2} - m_{Hyp}) / m_{#gamma#gamma#gamma#gamma}", "SubleadPs_interMass_compare_reweight.pdf")
    """

