from h4g_tools.utils.loading import loadPerMass, fillHist
from h4g_tools.bdt.event_selection_bdt import reduce_dataframe, predict_BDT
from h4g_tools.utils.scales import Scales as ScalesCls
from pathlib import Path
import xgboost
import os
import time
from copy import deepcopy
from ROOT import TCanvas, TLegend, kBlue, kBlack, kRed

if __name__ == "__main__":

    variation_path = "PATH-TO-VARIATIONS"

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
    # Hard-coded for 2024 --> will need to update for new years

    # Handle path setup after checking year
    Scales = ScalesCls(year) 

    # NEED TO UPDATE WITH NEW VALUES FOR NEW YEARS!
    bdt_cuts = {
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
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.14)

    legend = TLegend(0.55,0.7,0.93,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    canvas.SetLogy()

    bdt_out = {}
    variations = {}
    variations.update({"nominal": loadPerMass(original_path, branches=["weight", "pho1_mvaID", "pho2_mvaID", "pho3_mvaID", "pho4_mvaID", "BDT_score"])})
    for var in ["mvaID_Up", "mvaID_Down"]:
        bdt_out.update({var: {}})
        for N in range(1,5):
            bdt_out[var][f"{N}"] = {}
            key = f"pho{N}_mvaID_{var.replace('mvaID_','var')}"
            variations.update({key: loadPerMass(f"{variation_path}/{key}/", branches=["weight", "pho1_mvaID", "pho2_mvaID", "pho3_mvaID", "pho4_mvaID", "BDT_score"], variation=var)})
            
            for mass in sorted(variations[key].keys()):
                bdt_out[var][f"{N}"][mass] = {}
                legend.Clear()
                mva_nominal = fillHist(f"pho{N}_mvaID", variations["nominal"][mass], [-1.0, 1.0], binsScale=15.0, normalize=False)
                mva_original = fillHist(f"pho{N}_mvaID", variations[key][mass], [-1.0, 1.0], binsScale=15.0, normalize=False)
                mva_cut = fillHist(f"pho{N}_mvaID", variations[key][mass][variations[key][mass].BDT_score > bdt_cuts[mass]], [-1.0, 1.0], binsScale=15.0, normalize=False)

                mva_nominal.SetLineWidth(2)
                mva_nominal.SetLineColor(kBlack)
                mva_original.SetLineWidth(2)
                mva_original.SetLineColor(kBlue)
                mva_cut.SetLineWidth(2)
                mva_cut.SetLineColor(kRed)

                mva_nominal.Scale(Scales.lumi_fb * 1.0 / Scales.sumw["2024"][mass.replace("Signal_","")])
                mva_original.Scale(Scales.lumi_fb * 1.0 / Scales.sumw["2024"][mass.replace("Signal_","")])
                mva_cut.Scale(Scales.lumi_fb * 1.0 / Scales.sumw["2024"][mass.replace("Signal_","")])
                mva_original.GetYaxis().SetRangeUser(10**-3 - 5*10**-4, 10**3 - 5*10**2)
                mva_cut.GetYaxis().SetRangeUser(10**-3 - 5*10**-4, 10**3 - 5*10**2)
                mva_original.SetTitle("")
                mva_original.GetYaxis().SetTitle("Events")
                mva_original.GetXaxis().SetTitle(f"#gamma_{{{N}}} MVA ID")

                legend.AddEntry(mva_nominal, f"Nominal ({mva_original.Integral():.3f})", "l")
                legend.AddEntry(mva_original, f"{var} ({mva_original.Integral():.3f})", "l")
                legend.AddEntry(mva_cut, f"With BDT Cut ({mva_cut.Integral():.3f})", "l")

                mva_original.Draw("hist")
                mva_cut.Draw("hist same")
                mva_nominal.Draw("hist same")
                legend.Draw()
                canvas.SaveAs(f"{cwd}/plots/BDT/2024_updatedPresel/{mass.replace('Signal_','')}/pho{N}_{var}.pdf")

                # Plot BDT for each variation to see change as each photon is varied
                legend.Clear()
                bdt_nom = fillHist(f"BDT_score", variations["nominal"][mass], [0.0, 1.0], binsScale=30.0, normalize=False)
                bdt_original = fillHist(f"BDT_score", variations[key][mass], [0.0, 1.0], binsScale=30.0, normalize=False)
                bdt_cut = fillHist(f"BDT_score", variations[key][mass][variations[key][mass].BDT_score > bdt_cuts[mass]], [0.0, 1.0], binsScale=30.0, normalize=False)
                bdt_nom_cut = fillHist(f"BDT_score", variations["nominal"][mass][variations["nominal"][mass].BDT_score > bdt_cuts[mass]], [0.0, 1.0], binsScale=30.0, normalize=False)

                bdt_nom.SetTitle("")
                bdt_nom.GetYaxis().SetTitle("Events")
                bdt_nom.GetXaxis().SetTitle("BDT Score")

                bdt_original.SetLineWidth(2)
                bdt_original.SetLineColor(kBlue)
                bdt_nom.SetLineWidth(2)
                bdt_nom.SetLineColor(kRed)

                bdt_nom.Scale(Scales.lumi_fb * 1.0 / Scales.sumw["2024"][mass.replace("Signal_","")])
                bdt_original.Scale(Scales.lumi_fb * 1.0 / Scales.sumw["2024"][mass.replace("Signal_","")])
                bdt_cut.Scale(Scales.lumi_fb * 1.0 / Scales.sumw["2024"][mass.replace("Signal_","")])
                bdt_nom_cut.Scale(Scales.lumi_fb * 1.0 / Scales.sumw["2024"][mass.replace("Signal_","")])
                bdt_nom.GetYaxis().SetRangeUser(10**-3 - 5*10**-4, 10**3 - 5*10**2)
                bdt_original.GetYaxis().SetRangeUser(10**-3 - 5*10**-4, 10**3 - 5*10**2)
                bdt_cut.GetYaxis().SetRangeUser(10**-3 - 5*10**-4, 10**3 - 5*10**2)

                legend.AddEntry(bdt_nom, f"Nominal ({bdt_nom.Integral():.3f})", "l")
                legend.AddEntry(bdt_original, f"{var} ({bdt_original.Integral():.3f})", "l")
                legend.AddEntry(bdt_nom_cut, f"Cut on nominal: ({bdt_nom_cut.Integral():.3f})", "")
                legend.AddEntry(bdt_cut, f"Cut on {var}: ({bdt_cut.Integral():.3f})", "")

                bdt_out[var][f"{N}"][mass]["bdt_nom_cut"] = f"{bdt_nom_cut.Integral():.3f}"
                bdt_out[var][f"{N}"][mass][f"bdt_{var}_cut"] = f"{bdt_cut.Integral():.3f}"

                bdt_nom.Draw("hist")
                bdt_original.Draw("hist same")
                legend.Draw()
                canvas.SaveAs(f"{cwd}/plots/BDT/2024_updatedPresel/{mass.replace('Signal_','')}/bdt_pho{N}_{var}.pdf")



    for idx, var in enumerate(["mvaID_Up", "mvaID_Down"]):
        print(f"\nProcessing {var} variation:")
        for N in range(1,5):
            print(f"~ Photon {N}")
            for mass in sorted(variations[key].keys()):
                print(f"~~~ Mass: {mass.replace('Signal_','').replace('_',' ')}", end='\t')
                print(bdt_out[var][f"{N}"][mass]["bdt_nom_cut"], bdt_out[var][f"{N}"][mass][f"bdt_{var}_cut"], (bdt_out[var][f"{N}"][mass]["bdt_nom_cut"] <= bdt_out[var][f"{N}"][mass][f"bdt_{var}_cut"]) if var == "mvaID_Up" else (bdt_out[var][f"{N}"][mass]["bdt_nom_cut"] >= bdt_out[var][f"{N}"][mass][f"bdt_{var}_cut"]))

