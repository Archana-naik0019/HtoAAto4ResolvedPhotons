import os
from ROOT import TCanvas, gStyle, TLegend, kAzure, kGray, kMagenta, kOrange, kRed, kBlue
from h4g_tools.bdt.reweight import reweight1Dim
from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist
from h4g_tools.utils.plotting import plotFancy
from h4g_tools.utils.scales import Scales as ScalesCls

if __name__ == "__main__":

    pathSetup("2024")
    Scales = ScalesCls("2024")

    # Look at 30, 35, and 40 GeV data in SB
    samples = loadPerMass("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_mixedCheck_02Mar2026", branches=["BDT_score"])
    samples_unmixed = loadPerMass("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_wBDT_genMatching_06Feb2026", branches=["BDT_score"])
    bkg = loadSamples("/cms/cephfs/data/store/user/castells/outputs/BDT/2024/ES_BDT_2024/generic/background/", branches=["BDT_score", "mass_gggg"])
    data = loadSamples("/cms/cephfs/data/store/user/castells/outputs/BDT/2024/ES_BDT_2024/generic/data/", branches=["BDT_score", "mass_gggg"])

    bkg = reweight1Dim(bkg, data, "BDT_score")

    hists = {}
    hists_unmixed = {}
    bkg_hist = fillHist("BDT_score", bkg, [0.0, 1.0], binsScale=30, normalize=False, sb = True, preventOverFlow = True, customWeights="1dimreweight")
    data_hist = fillHist("BDT_score", data, [0.0, 1.0], binsScale=30, normalize=False, sb = True, preventOverFlow = True)
    for m, s in samples.items():
        # Can look at full region since it's signal
        hists.update({m: fillHist("BDT_score", s, [0.0, 1.0], binsScale=30, normalize=False, sb = False, preventOverFlow = True)})
        hists[m].SetTitle("")
        hists[m].SetLineWidth(2)
    for m, s in samples_unmixed.items():
        # Can look at full region since it's signal
        hists_unmixed.update({m: fillHist("BDT_score", s, [0.0, 1.0], binsScale=30, normalize=False, sb = False, preventOverFlow = True)})
        hists_unmixed[m].SetTitle("")
        hists_unmixed[m].SetLineWidth(2)
    bkg_hist.Scale(data_hist.Integral() / bkg_hist.Integral())
 
    # Canvas setup
    canvas = TCanvas("canvas", "canvas", 1000, 1000)
    canvas.SetLogy()

    # Legend setup
    legend = TLegend(0.6,0.58,0.87,0.83)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    max_order = 6
    
    bkg_hist.SetTitle("")
    bkg_hist.SetFillColor(kAzure-4)
    bkg_hist.GetYaxis().SetTitle("Events")
    bkg_hist.GetXaxis().SetTitle("BDT Score")
    bkg_hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, 10 * (10**max_order - 5*10**(max_order-1)))
    legend.AddEntry(bkg_hist, f"Background", "f")

    background_error = bkg_hist.Clone("background_error")
    background_error.SetTitle("")
    background_error.SetLineColor(kGray+3)
    background_error.SetFillColor(kGray+3)
    background_error.SetFillStyle(3008)
    legend.AddEntry(background_error, "Statistical Uncertainty", "f")

    new_hists = {}
    for k in sorted(hists.keys()):
        if k in [f"Signal_{skip}_GeV" for skip in [20, 25, 35, 45, 50, 55]]:
            continue

        new_hists.update({k: hists[k]})
        
    bkg_hist.Draw("hist")
    background_error.Draw("E2 same")
    for (k, h), c in zip(new_hists.items(), [kMagenta, kOrange, kRed, kBlue]):
        m = int(k.replace("Signal_","").replace("_GeV",""))
        h.SetTitle("")
        h.SetLineWidth(2)
        h.SetLineColor(c)
        h.GetYaxis().SetTitle("Events")
        h.GetXaxis().SetTitle("BDT Score")
        h.GetYaxis().SetRangeUser(10**-3 + 10**-4, 10 * (10**max_order - 5*10**(max_order-1)))
        h.Scale(Scales.lumi_fb * 1.0 / Scales.sumw["2024"][f"{m}_GeV"])
        legend.AddEntry(h, f"Mixed Signal MC {m} GeV", "l")

        h.Draw("hist same") 

    legend.Draw("same")
    plotFancy(canvas, "BDT Score (Mixed)", lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, sub="Bkg restricted to SB", Scales=Scales)

    # Save figure
    path = "plots/BDT/2024/sig_mix/bdt_spectra.pdf"
    canvas.SaveAs(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path} to file.")

    ### Plot unmixed version
    # Canvas setup
    canvas = TCanvas("canvas2", "canvas2", 1000, 1000)
    canvas.SetLogy()

    # Legend setup
    legend = TLegend(0.6,0.58,0.87,0.83)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    legend.AddEntry(bkg_hist, f"Background", "f")
    legend.AddEntry(background_error, "Statistical Uncertainty", "f")

    new_hists_unmixed = {}
    for k in sorted(hists.keys()):
        if k in [f"Signal_{skip}_GeV" for skip in [20, 25, 35, 45, 50, 55]]:
            continue

        new_hists_unmixed.update({k: hists_unmixed[k]})
        
    bkg_hist.Draw("hist")
    background_error.Draw("E2 same")
    for (k, h), c in zip(new_hists_unmixed.items(), [kMagenta, kOrange, kRed, kBlue]):
        m = int(k.replace("Signal_","").replace("_GeV",""))
        h.SetTitle("")
        h.SetLineWidth(2)
        h.SetLineColor(c)
        h.GetYaxis().SetTitle("Events")
        h.GetXaxis().SetTitle("BDT Score")
        h.GetYaxis().SetRangeUser(10**-3 + 10**-4, 10 * (10**max_order - 5*10**(max_order-1)))
        h.Scale(Scales.lumi_fb * 1.0 / Scales.sumw["2024"][f"{m}_GeV"])
        legend.AddEntry(h, f"Signal MC {m} GeV      ", "l")

        h.Draw("hist same") 

    legend.Draw("same")
    plotFancy(canvas, "BDT Score (Standard)", lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, sub="Bkg restricted to SB", Scales=Scales)

    # Save figure
    path = "plots/BDT/2024/sig_mix/bdt_spectra_unmixed.pdf"
    canvas.SaveAs(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path} to file.")
