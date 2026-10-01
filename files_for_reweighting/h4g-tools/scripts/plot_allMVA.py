from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist
from h4g_tools.utils.scales import Scales as ScalesCls
from ROOT import TCanvas, TLegend, kAzure, kBlack, kGreen, kMagenta, kYellow, kSpring, kOrange, kCyan, kRed, kViolet, kOrange, kBlue
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import json
import os

if __name__ == "__main__":

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    year = "2024"
    masses = [x for x in range(15, 65, 5)]
    eras={
        "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
        "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
        "signal": [f"Signal_{m}_GeV" for m in masses]
    }

    # Handle path setup after checking year
    Scales = ScalesCls(year) 

    # Load samples
    print("Loading signal samples...")
    samples = loadPerMass("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSFpostBPix_TriggerEff_updatedPresel_wBDT_19May2026", branches=["weight", "pho1_mvaID", "pho2_mvaID", "pho3_mvaID", "pho4_mvaID", "BDT_score"])

    bdt = {
        "Signal_15_GeV": 0.9600,
        "Signal_20_GeV": 0.9600,
        "Signal_25_GeV": 0.9600,
        "Signal_30_GeV": 0.9700,
        "Signal_35_GeV": 0.9650,
        "Signal_40_GeV": 0.9700,
        "Signal_45_GeV": 0.9750,
        "Signal_50_GeV": 0.9800,
        "Signal_55_GeV": 0.9800,
        "Signal_60_GeV": 0.9800,
    }

    canvas = TCanvas("canvas", "canvas", 1000, 1000)
    legend = TLegend(0.55,0.6,0.88,0.9)
    legend.SetNColumns(2)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    canvas.SetLogy()
    colors = [kAzure-4, kBlack, kGreen+3, kMagenta, kYellow, kSpring-3, kOrange, kCyan+1, kRed, kViolet-5, kOrange+4, kBlue]
    for N in range(1,5):
        hists = {}
        for idx, (mass, color) in enumerate(zip(sorted(samples.keys()), colors)):
            # Make plots of photon pT comparing the two samples
            hists.update({mass: fillHist(f"pho{N}_mvaID", samples[mass], [-1.0, 1.0], binsScale=30.0, normalize=False)})

            # Scale plots
            #hists[mass].Scale(Scales.lumi_fb * 1.0 / Scales.sumw["2024"][mass.replace("Signal_","")])

            hists[mass].SetLineColor(color)
            hists[mass].SetLineWidth(2)
            hists[mass].GetXaxis().SetTitle(f"#gamma_{{{N}}} MVA ID")
            hists[mass].GetYaxis().SetTitle("Events")
            hists[mass].SetTitle("")
            #hists[mass].GetYaxis().SetRangeUser(10**-2 - 5*10**-3, 10**5 - 5*10**4)
            hists[mass].GetYaxis().SetRangeUser(0.5, 10**8 - 5*10**7)
            legend.AddEntry(hists[mass], f"m_{{a}} = {mass.replace('Signal_','').replace('_',' ')}", "l")
            hists[mass].Draw("hist same" if idx > 0 else "hist")

        legend.Draw("same")
        canvas.SaveAs(f"plots/BDT/2024/compare_pho{N}_mvaID_wBDTCuts.pdf")
        legend.Clear()


