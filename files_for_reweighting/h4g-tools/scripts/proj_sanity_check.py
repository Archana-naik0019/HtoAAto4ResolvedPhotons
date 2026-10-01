import os
import argparse
import pandas as pd
from ROOT import TH1D, TH2D, TCanvas, TPad, kRed, kBlack, TLine
from h4g_tools.utils.remove_smearing import recalc_nosmear
from h4g_tools.utils.scales import Scales as ScalesCls
from h4g_tools.utils.loading import fillHist
from h4g_tools.utils.plotting import plotCorrectedComparison, plotFancy


scale_order = 4
scale = 10**scale_order


if __name__ == "__main__":

    # Sample path -- change this as desired
    year = "2024"
    sp = "/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_wBDT_04Feb2026/Signal_{0}_GeV/nominal/"

    if not os.path.exists("smearing_sanity_checks/"):
        os.mkdir("smearing_sanity_checks/")

    canvas = TCanvas("c", "c", 1000, 1000)
    for mass, shift in zip([m for m in range(15,65,5)], [5, 5, 5, 5, 10, 10, 10, 10, 10, 10]):
        canvas.Clear()
        print(f"Loading only {mass} GeV corrected sample...")
        sigC = pd.read_parquet(sp.format(mass))

        print("\nRecalculating masses for unsmeared events...")
        sigUC = recalc_nosmear(sigC)

        # Look at m4g only  --> bins = (135-110)*4 = 100
        m4g_corr = fillHist("mass_gggg", sigC, [110.0, 135.0], binsScale=4, name=f"sig{mass}_corr", normalize=False, preventOverFlow=True)
        m4g_uncorr = fillHist("mass_gggg", sigUC, [110.0, 135.0], binsScale=4, name=f"sig{mass}_uncorr", normalize=False, preventOverFlow=True)
        ma1_corr = fillHist("LeadPs_mass", sigC, [mass-shift, mass+shift], binsScale=10, name=f"sig{mass}_corr", normalize=False, preventOverFlow=True)
        ma1_uncorr = fillHist("LeadPs_mass", sigUC, [mass-shift, mass+shift], binsScale=10, name=f"sig{mass}_uncorr", normalize=False, preventOverFlow=True)
        ma2_corr = fillHist("SubleadPs_mass", sigC, [mass-shift, mass+shift], binsScale=10, name=f"sig{mass}_corr", normalize=False, preventOverFlow=True)
        ma2_uncorr = fillHist("SubleadPs_mass", sigUC, [mass-shift, mass+shift], binsScale=10, name=f"sig{mass}_uncorr", normalize=False, preventOverFlow=True)

        m4g_corr.Scale(scale / m4g_corr.Integral())
        m4g_uncorr.Scale(scale / m4g_uncorr.Integral())
        ma1_corr.Scale(scale / ma1_corr.Integral())
        ma1_uncorr.Scale(scale / ma1_uncorr.Integral())
        ma2_corr.Scale(scale / ma2_corr.Integral())
        ma2_uncorr.Scale(scale / ma2_uncorr.Integral())
    

        ### Plot 1D
        plotCorrectedComparison({"sigUC": m4g_uncorr, "sigC": m4g_corr}, "m_{4#gamma}", f"smearing_sanity_checks/{mass}_GeV_m4g_1D.png", log=False, minimum=0.0, year=year, Scales=ScalesCls(year), sim=True, legend_text={"sigUC": "All corrections without smearing", "sigC": "All corrections with smearing"})
        plotCorrectedComparison({"sigUC": ma1_uncorr, "sigC": ma1_corr}, "m_{a1}", f"smearing_sanity_checks/{mass}_GeV_ma1_1D.png", log=False, minimum=0.0, year=year, Scales=ScalesCls(year), sim=True, legend_text={"sigUC": "All corrections without smearing", "sigC": "All corrections with smearing"})
        plotCorrectedComparison({"sigUC": ma2_uncorr, "sigC": ma2_corr}, "m_{a2}", f"smearing_sanity_checks/{mass}_GeV_ma2_1D.png", log=False, minimum=0.0, year=year, Scales=ScalesCls(year), sim=True, legend_text={"sigUC": "All corrections without smearing", "sigC": "All corrections with smearing"})


        ### Manually plot 2D
        hist2D = TH2D(f"hist2D_{mass}", f"hist2D_{mass}", 100, 110.0, 135.0, 100, 110.0, 135.0)
        hist2Da1 = TH2D(f"hist2Da1_{mass}", f"hist2Da1_{mass}", 2*shift*10, mass-shift, mass+shift, 2*shift*10, mass-shift, mass+shift)
        hist2Da2 = TH2D(f"hist2Da2_{mass}", f"hist2Da2_{mass}", 2*shift*10, mass-shift, mass+shift, 2*shift*10, mass-shift, mass+shift)
        for idx in range(len(sigC)):
            hist2D.Fill(sigUC["mass_gggg"].iloc[idx], sigC["mass_gggg"].iloc[idx], sigC["weight"].iloc[idx])
            hist2Da1.Fill(sigUC["LeadPs_mass"].iloc[idx], sigC["LeadPs_mass"].iloc[idx], sigC["weight"].iloc[idx])
            hist2Da2.Fill(sigUC["SubleadPs_mass"].iloc[idx], sigC["SubleadPs_mass"].iloc[idx], sigC["weight"].iloc[idx])

        hist2D.SetStats(0)
        hist2D.SetMarkerStyle(8)
        hist2D.SetMarkerSize(0.5)
        hist2D.SetTitle("")
        hist2D.GetYaxis().SetTitleFont(42)
        hist2D.GetXaxis().SetTitleFont(42)
        hist2D.GetXaxis().SetTitleOffset(1.5)
        hist2D.GetXaxis().CenterTitle(True)
        hist2D.GetYaxis().CenterTitle(True)
        hist2D.GetXaxis().SetTitle("Recalculated Masses m_{4#gamma}")
        hist2D.GetYaxis().SetTitle("Corrected m_{4#gamma}")
        
        hist2Da1.SetStats(0)
        hist2Da1.SetMarkerStyle(8)
        hist2Da1.SetMarkerSize(0.5)
        hist2Da1.SetTitle("")
        hist2Da1.GetYaxis().SetTitleFont(42)
        hist2Da1.GetXaxis().SetTitleFont(42)
        hist2Da1.GetXaxis().SetTitleOffset(1.5)
        hist2Da1.GetXaxis().CenterTitle(True)
        hist2Da1.GetYaxis().CenterTitle(True)
        hist2Da1.GetXaxis().SetTitle("Recalculated Masses m_{a1}")
        hist2Da1.GetYaxis().SetTitle("Corrected m_{a1}")
        
        hist2Da2.SetStats(0)
        hist2Da2.SetMarkerStyle(8)
        hist2Da2.SetMarkerSize(0.5)
        hist2Da2.SetTitle("")
        hist2Da2.GetYaxis().SetTitleFont(42)
        hist2Da2.GetXaxis().SetTitleFont(42)
        hist2Da2.GetXaxis().SetTitleOffset(1.5)
        hist2Da2.GetXaxis().CenterTitle(True)
        hist2Da2.GetYaxis().CenterTitle(True)
        hist2Da2.GetXaxis().SetTitle("Recalculated Masses m_{a2}")
        hist2Da2.GetYaxis().SetTitle("Corrected m_{a2}")

        line = TLine(110.0, 110.0, 135.0, 135.0)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.SetLineColor(kBlack)

        hist2D.Draw("COLZ")
        line.Draw("same")
        plotFancy(canvas, "m_{4#gamma}", prelim=True, simulation=True, inPlot=False, colz=True)
        canvas.SaveAs(f"smearing_sanity_checks/{mass}_GeV_m4g_2D.png")

        line_ma = TLine(mass-shift, mass+shift, mass-shift, mass+shift)
        line_ma.SetLineStyle(2)
        line_ma.SetLineWidth(2)
        line_ma.SetLineColor(kBlack)

        hist2Da1.Draw("COLZ")
        line_ma.Draw("same")
        plotFancy(canvas, "m_{a1}", prelim=True, simulation=True, inPlot=False, colz=True)
        canvas.SaveAs(f"smearing_sanity_checks/{mass}_GeV_ma1_2D.png")

        hist2Da2.Draw("COLZ")
        line_ma.Draw("same")
        plotFancy(canvas, "m_{a2}", prelim=True, simulation=True, inPlot=False, colz=True)
        canvas.SaveAs(f"smearing_sanity_checks/{mass}_GeV_ma2_2D.png")


        ### Plot 2D projections
        px = hist2D.ProjectionX()
        py = hist2D.ProjectionY()
        px.Scale(scale / px.Integral())
        py.Scale(scale / py.Integral()) 
        plotCorrectedComparison({"sigUC": px, "sigC": py}, "m_{4#gamma}", f"smearing_sanity_checks/{mass}_GeV_m4g_2D_proj.png", log=False, minimum=0.0, year=year, Scales=ScalesCls(year), sim=True, legend_text={"sigUC": "All corrections without smearing", "sigC": "All corrections with smearing"})

        px = hist2Da1.ProjectionX()
        py = hist2Da1.ProjectionY()
        px.Scale(scale / px.Integral())
        py.Scale(scale / py.Integral()) 
        plotCorrectedComparison({"sigUC": px, "sigC": py}, "m_{a1}", f"smearing_sanity_checks/{mass}_GeV_ma1_2D_proj.png", log=False, minimum=0.0, year=year, Scales=ScalesCls(year), sim=True, legend_text={"sigUC": "All corrections without smearing", "sigC": "All corrections with smearing"})

        px = hist2Da2.ProjectionX()
        py = hist2Da2.ProjectionY()
        px.Scale(scale / px.Integral())
        py.Scale(scale / py.Integral()) 
        plotCorrectedComparison({"sigUC": px, "sigC": py}, "m_{a2}", f"smearing_sanity_checks/{mass}_GeV_ma2_2D_proj.png", log=False, minimum=0.0, year=year, Scales=ScalesCls(year), sim=True, legend_text={"sigUC": "All corrections without smearing", "sigC": "All corrections with smearing"})
