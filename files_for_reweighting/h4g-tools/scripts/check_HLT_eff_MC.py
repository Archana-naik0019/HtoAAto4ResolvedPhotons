from h4g_tools.utils.scales import Scales as ScalesCls
from h4g_tools.utils.loading import pathSetup, loadPerMass
from h4g_tools.utils.plotting import plotFancy
from ROOT import TLegend, kBlack, kBlue, kRed, kGreen, TCanvas, TH1D
from array import array
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import json
import os

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", type=str, required=True, help="Signal samples to plot per mass.")
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    parser.add_argument("-bdt", "--bdt_path", type=str, required=True, help="Path to file containing efficiency of BDT cuts. Assumes 1 category.")
    args = parser.parse_args()

    # Load bdt cuts
    with open(args.bdt_path, "r") as f:
        bdt_cuts = json.load(f)

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
    elif "2024" in args.year:
        eras={
            "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
            "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
            "signal": [f"Signal_{m}_GeV" for m in masses]
        }
    else:
        assert False, "Must specify year in Scales class!"
    # WILL NEED TO ADD YEAR HERE!

    # Handle path setup after checking year
    pathSetup(args.year)

    # Load samples
    print("Loading signal samples...")
    samples = loadPerMass(args.inputs, eras = eras, year = args.year, branches=["HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId", "BDT_score", "pho1_pt", "pho2_pt", "pho3_pt", "pho4_pt", "weight"])
    samples_merged = None
    samples_merged_bdt = None

    eff = {}
    eff_bdt = {}
    pho1_pt = []
    pho2_pt = []
    pho3_pt = []
    pho4_pt = []
    sorted_samples = sorted(samples.keys())
    for idx, mass in enumerate(sorted_samples):
        sample = samples[mass]
        if samples_merged is None:
            samples_merged = sample
        else:
            samples_merged = pd.concat([samples_merged, sample], axis=0)

        eff.update({mass: sum(sample[sample["HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId"] == 1].weight) / sum(sample.weight)})
        print(f"eff     {mass}: {sum(sample[sample['HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId'] == 1].weight):.2f} {sum(sample.weight):.2f} {eff[mass]*100:.2f}%")
        if args.bdt_path is not None:
            sample = sample[sample.BDT_score > bdt_cuts[mass.replace("Signal_","")]["cut0"][0]]
            eff_bdt.update({mass: sum(sample[sample["HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId"] == 1].weight) / sum(sample.weight)})
            #print(f"eff_bdt {mass}: {sum(sample[sample['HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId'] == 1].weight)} {sum(sample.weight)} {eff_bdt[mass]}")

        if samples_merged_bdt is None:
            samples_merged_bdt = sample
        else:
            samples_merged_bdt = pd.concat([samples_merged_bdt, sample], axis=0)

    # Make pt-bins
    pt_bins = [b for b in range(0, 310, 5)]
    pt_bins = array('d', pt_bins)
    pho1_tmp1 = TH1D("pho1_pt_tmp1", "pho1_pt", int(300/5), pt_bins)
    pho2_tmp1 = TH1D("pho2_pt_tmp1", "pho2_pt", int(300/5), pt_bins)
    pho3_tmp1 = TH1D("pho3_pt_tmp1", "pho3_pt", int(300/5), pt_bins)
    pho4_tmp1 = TH1D("pho4_pt_tmp1", "pho4_pt", int(300/5), pt_bins)
    pho1_tmp2 = TH1D("pho1_pt_tmp2", "pho1_pt", int(300/5), pt_bins)
    pho2_tmp2 = TH1D("pho2_pt_tmp2", "pho2_pt", int(300/5), pt_bins)
    pho3_tmp2 = TH1D("pho3_pt_tmp2", "pho3_pt", int(300/5), pt_bins)
    pho4_tmp2 = TH1D("pho4_pt_tmp2", "pho4_pt", int(300/5), pt_bins)
    pho1_tmp1_bdt = TH1D("pho1_pt_tmp1_bdt", "pho1_pt", int(300/5), pt_bins)
    pho2_tmp1_bdt = TH1D("pho2_pt_tmp1_bdt", "pho2_pt", int(300/5), pt_bins)
    pho3_tmp1_bdt = TH1D("pho3_pt_tmp1_bdt", "pho3_pt", int(300/5), pt_bins)
    pho4_tmp1_bdt = TH1D("pho4_pt_tmp1_bdt", "pho4_pt", int(300/5), pt_bins)
    pho1_tmp2_bdt = TH1D("pho1_pt_tmp2_bdt", "pho1_pt", int(300/5), pt_bins)
    pho2_tmp2_bdt = TH1D("pho2_pt_tmp2_bdt", "pho2_pt", int(300/5), pt_bins)
    pho3_tmp2_bdt = TH1D("pho3_pt_tmp2_bdt", "pho3_pt", int(300/5), pt_bins)
    pho4_tmp2_bdt = TH1D("pho4_pt_tmp2_bdt", "pho4_pt", int(300/5), pt_bins)
    pho1_pt_hist = TH1D("pho1_pt", "pho1_pt", int(300/5), pt_bins)
    pho2_pt_hist = TH1D("pho2_pt", "pho2_pt", int(300/5), pt_bins)
    pho3_pt_hist = TH1D("pho3_pt", "pho3_pt", int(300/5), pt_bins)
    pho4_pt_hist = TH1D("pho4_pt", "pho4_pt", int(300/5), pt_bins)
    pho1_pt_hist_bdt = TH1D("pho1_pt_bdt", "pho1_pt_bdt", int(300/5), pt_bins)
    pho2_pt_hist_bdt = TH1D("pho2_pt_bdt", "pho2_pt_bdt", int(300/5), pt_bins)
    pho3_pt_hist_bdt = TH1D("pho3_pt_bdt", "pho3_pt_bdt", int(300/5), pt_bins)
    pho4_pt_hist_bdt = TH1D("pho4_pt_bdt", "pho4_pt_bdt", int(300/5), pt_bins)
    for idx, pt in enumerate(pt_bins[1:]):
        for i, (n, d) in enumerate(zip([pho1_tmp1, pho2_tmp1, pho3_tmp1, pho4_tmp1], [pho1_tmp2, pho2_tmp2, pho3_tmp2, pho4_tmp2])):
            v = len(samples_merged[(samples_merged[f"pho{i+1}_pt"] >= pt_bins[idx]) & (samples_merged[f"pho{i+1}_pt"] < pt) & (samples_merged["HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId"] == 1)])
            t = len(samples_merged[(samples_merged[f"pho{i+1}_pt"] >= pt_bins[idx]) & (samples_merged[f"pho{i+1}_pt"] < pt)])
            n.SetBinContent(idx+1, v)
            d.SetBinContent(idx+1, t)
        for i, (n, d) in enumerate(zip([pho1_tmp1_bdt, pho2_tmp1_bdt, pho3_tmp1_bdt, pho4_tmp1_bdt], [pho1_tmp2_bdt, pho2_tmp2_bdt, pho3_tmp2_bdt, pho4_tmp2_bdt])):
            v = len(samples_merged_bdt[(samples_merged_bdt[f"pho{i+1}_pt"] >= pt_bins[idx]) & (samples_merged_bdt[f"pho{i+1}_pt"] < pt) & (samples_merged_bdt["HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId"] == 1)])
            t = len(samples_merged_bdt[(samples_merged_bdt[f"pho{i+1}_pt"] >= pt_bins[idx]) & (samples_merged_bdt[f"pho{i+1}_pt"] < pt)])
            n.SetBinContent(idx+1, v)
            d.SetBinContent(idx+1, t)

    pho1_pt_hist.Divide(pho1_tmp1, pho1_tmp2, 1.0, 1.0, "B")
    pho2_pt_hist.Divide(pho2_tmp1, pho2_tmp2, 1.0, 1.0, "B")
    pho3_pt_hist.Divide(pho3_tmp1, pho3_tmp2, 1.0, 1.0, "B")
    pho4_pt_hist.Divide(pho4_tmp1, pho4_tmp2, 1.0, 1.0, "B")
    pho1_pt_hist_bdt.Divide(pho1_tmp1_bdt, pho1_tmp2_bdt, 1.0, 1.0, "B")
    pho2_pt_hist_bdt.Divide(pho2_tmp1_bdt, pho2_tmp2_bdt, 1.0, 1.0, "B")
    pho3_pt_hist_bdt.Divide(pho3_tmp1_bdt, pho3_tmp2_bdt, 1.0, 1.0, "B")
    pho4_pt_hist_bdt.Divide(pho4_tmp1_bdt, pho4_tmp2_bdt, 1.0, 1.0, "B")

    pho1_pt_hist.Scale(100.0)
    pho2_pt_hist.Scale(100.0)
    pho3_pt_hist.Scale(100.0)
    pho4_pt_hist.Scale(100.0)
    pho1_pt_hist_bdt.Scale(100.0)
    pho2_pt_hist_bdt.Scale(100.0)
    pho3_pt_hist_bdt.Scale(100.0)
    pho4_pt_hist_bdt.Scale(100.0)

    pho1_pt_hist.GetYaxis().SetRangeUser(40.0, 105.0)
    pho2_pt_hist.GetYaxis().SetRangeUser(40.0, 105.0)
    pho3_pt_hist.GetYaxis().SetRangeUser(40.0, 105.0)
    pho4_pt_hist.GetYaxis().SetRangeUser(40.0, 105.0)
    pho1_pt_hist_bdt.GetYaxis().SetRangeUser(40.0, 105.0)
    pho2_pt_hist_bdt.GetYaxis().SetRangeUser(40.0, 105.0)
    pho3_pt_hist_bdt.GetYaxis().SetRangeUser(40.0, 105.0)
    pho4_pt_hist_bdt.GetYaxis().SetRangeUser(40.0, 105.0)

    for path, (pho1, pho2, pho3, pho4) in zip(["", "_bdt"], [[pho1_pt_hist, pho2_pt_hist, pho3_pt_hist, pho4_pt_hist], [pho1_pt_hist_bdt, pho2_pt_hist_bdt, pho3_pt_hist_bdt, pho4_pt_hist_bdt]]):
        wbdt = ""
        if "bdt" in path:
            wbdt = " (with BDT cuts)"

        # Make legend
        legend = TLegend(0.8,0.6-0.1,0.92,0.75-0.1)
        legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)

        # Set line colors and width
        pho1.SetLineColor(kBlack)
        pho2.SetLineColor(kRed)
        pho3.SetLineColor(kBlue)
        pho4.SetLineColor(kGreen+3)
        pho1.SetLineWidth(2)
        pho2.SetLineWidth(2)
        pho3.SetLineWidth(2)
        pho4.SetLineWidth(2)
        pho1.SetStats(0)
        pho2.SetStats(0)
        pho3.SetStats(0)
        pho4.SetStats(0)
        pho1.SetTitle("")
        pho2.SetTitle("")
        pho3.SetTitle("")
        pho4.SetTitle("")
        pho1.GetXaxis().SetTitle("p_{T}")
        pho1.GetYaxis().SetTitle("HLT Efficiency (Events / 5 GeV)")

        canvas = TCanvas("canvas", "canvas", 1000, 1000)
        canvas.SetLeftMargin(0.14)
        canvas.SetRightMargin(0.05)
        canvas.SetBottomMargin(0.14)
        pho1.Draw("hist same")
        pho2.Draw("hist same")
        pho3.Draw("hist same")
        pho4.Draw("hist same")
        legend.AddEntry(pho1, "#gamma_{1} p_{T}", "l")
        legend.AddEntry(pho2, "#gamma_{2} p_{T}", "l")
        legend.AddEntry(pho3, "#gamma_{3} p_{T}", "l")
        legend.AddEntry(pho4, "#gamma_{4} p_{T}", "l")
        legend.Draw("same")

        plotFancy(canvas, f"HLT Eff vs #gamma pT{wbdt}", prelim=True, simulation=True, inPlot=False, Scales=Scales)
        canvas.SaveAs(f"{cwd}/plots/standard/{args.year}/HLT_bit_eff_pho_pt{path}.pdf")

    # Plot HLT bit efficiency from signal MC
    fig, ax1 = plt.subplots(figsize = (6,6))
    title = f"HLT Bit Efficiency"
    x = list(range(15,65,5))
    y = [e*100 for e in list(eff.values())]
    plt.plot(x, y, "-o")
    if args.bdt_path is not None:
        y_bdt = [e*100 for e in list(eff_bdt.values())]
    #plt.plot(x, y_bdt, "-o")
    plt.xlabel("m(a) [GeV]", fontsize=14)
    plt.ylabel("HLT Bit Efficiency (%)", fontsize=14)
    plt.grid()
    #plt.legend(labels=["Without BDT Cut", "With BDT cut"], frameon=True, loc="upper right", fontsize=14)

    ax = plt.gca()
    ymin = 90.0
    ymax = 100.0
    ax.set_ylim(ymin, ymax)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig(f"{cwd}/plots/standard/{args.year}/HLT_bit_eff.pdf")
    print(f"Saved plot to {cwd}/plots/standard/{args.year}/HLT_bit_eff.pdf")

