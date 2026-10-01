from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist
from h4g_tools.utils.scales import Scales as ScalesCls
from ROOT import TCanvas, kBlue, kRed, TLegend
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
    ## Only trigger on/off ##
    """
    samples = loadPerMass("/users/scastel2/HiggsDNA/higgs_dna/scripts/outputs_signal_2024_TriggerEff_25May2026", branches=["weight", "pho1_pt", "pho2_pt", "pho3_pt", "pho4_pt", "BDT_score"])
    samples_noTrigEff = loadPerMass("/users/scastel2/HiggsDNA/higgs_dna/scripts/outputs_signal_2024_noCorr_25May2026", branches=["weight", "pho1_pt", "pho2_pt", "pho3_pt", "pho4_pt", "BDT_score"])
    """
    samples = loadPerMass("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSFpostBPix_TriggerEff_updatedPresel_wBDT_19May2026", branches=["weight", "pho1_pt", "pho2_pt", "pho3_pt", "pho4_pt", "BDT_score"])
    samples_noTrigEff = loadPerMass("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSFpostBPix_updatedPresel_wBDT_04May2026", branches=["weight", "pho1_pt", "pho2_pt", "pho3_pt", "pho4_pt", "BDT_score"])

    # below: wrong samples that I sent Tanay a while back. Got confused
    #samples_noTrigEff = loadPerMass("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSFpostBPix_wBDT_wHLT_28Apr2026/", branches=["weight", "pho1_pt", "pho2_pt", "pho3_pt", "pho4_pt", "BDT_score"])

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
    for mass in samples.keys():
        # Apply bdt cut
        samples[mass] = samples[mass][samples[mass].BDT_score > bdt[mass]]
        samples_noTrigEff[mass] = samples_noTrigEff[mass][samples_noTrigEff[mass].BDT_score > bdt[mass]]
        print(mass, f"With Trigger Eff:{sum(samples[mass].weight):.2f}", f"No Trigger Eff {sum(samples_noTrigEff[mass].weight):.2f}")

        # Plot weights distribution
        x = samples[mass].weight.to_list()
        print(f"Average: {sum(x)/len(x)}")
        fig, ax1 = plt.subplots(figsize = (6,6))
        title = f"Weight value per event (TriggerSF only)"
        bins = [0.70 + 0.01 * i for i in range(0, 31)]
        plt.hist(x, bins)
        plt.xlabel("Event Weight (Trigger SF)", fontsize=14)
        plt.ylabel("Events", fontsize=14)
        plt.grid()

        ax = plt.gca()

        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        plt.savefig(f"{cwd}/plots/BDT/2024/{mass.replace('Signal_','')}/triggerSF_weights_dist.pdf")
        print(f"Saved plot to {cwd}/plots/BDT/2024/{mass.replace('Signal_','')}/triggerSF_weights_dist.pdf")

        # Make plots of photon pT comparing the two samples
        pho1_trig = fillHist("pho1_pt", samples[mass], [0.0, 300.0], binsScale=0.2,  normalize=False)
        pho2_trig = fillHist("pho4_pt", samples[mass], [0.0, 200.0], binsScale=0.2, normalize=False)
        pho3_trig = fillHist("pho4_pt", samples[mass], [0.0, 100.0], binsScale=0.2, normalize=False)
        pho4_trig = fillHist("pho4_pt", samples[mass], [0.0, 80.0], binsScale=0.2, normalize=False)

        pho1_noTrig = fillHist("pho1_pt", samples_noTrigEff[mass], [0.0, 300.0], binsScale=0.2,  normalize=False)
        pho2_noTrig = fillHist("pho2_pt", samples_noTrigEff[mass], [0.0, 200.0], binsScale=0.2, normalize=False)
        pho3_noTrig = fillHist("pho3_pt", samples_noTrigEff[mass], [0.0, 100.0], binsScale=0.2, normalize=False)
        pho4_noTrig = fillHist("pho4_pt", samples_noTrigEff[mass], [0.0, 80.0], binsScale=0.2, normalize=False)

        for N, (h1, h2) in enumerate(zip([pho1_trig, pho2_trig, pho3_trig, pho4_trig], [pho1_noTrig, pho2_noTrig, pho3_noTrig, pho4_noTrig])):
            canvas = TCanvas("canvas", "canvas", 1000, 1000)
            legend = TLegend(0.65,0.65,0.87,0.85)
            legend.SetTextFont(42)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            canvas.SetLogy()
            h1.SetLineColor(kBlue)
            h1.SetLineWidth(2)
            h2.SetLineColor(kRed)
            h2.SetLineWidth(2)
            h2.GetXaxis().SetTitle("p_{T}")
            h2.GetYaxis().SetTitle("Events / 5 GeV")
            h1.SetTitle("")
            h2.SetTitle(f"#gamma_{{{N+1}}} p_{{T}}")
            h2.GetYaxis().SetRangeUser(0.5, 10**5 - 5*10**4)
            legend.AddEntry(h2, "No Trigger Eff", "l")
            legend.AddEntry(h1, "With Trigger Eff", "l")
            #print(f"With Trigger Eff:{h1.Integral():.2f}", f"No Trigger Eff {h2.Integral():.2f}")

            h2.Draw("hist")
            h1.Draw("hist same")
            legend.Draw("same")
            canvas.SaveAs(f"plots/BDT/2024/{mass.replace('Signal_','')}/trigger_compare_pho{N+1}.pdf")


