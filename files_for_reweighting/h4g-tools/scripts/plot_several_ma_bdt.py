import pandas as pd
import xgboost
import os
import argparse
from h4g_tools.bdt.event_selection_bdt import reduce_dataframe, predict_BDT
from h4g_tools.utils.scales import Scales as ScalesCls, SideBand
from h4g_tools.utils.loading import loadPerMass, loadSamples, fillHist
from h4g_tools.utils.plotting import plotFancy, residuals_error_prop
from h4g_tools.bdt.reweight import reweight1Dim
from ROOT import TCanvas, TLegend, gStyle, kGreen, kMagenta, kBlue, kAzure, kRed, kOrange, kBlack, TPad, gPad, TLine


if __name__ == "__main__":

    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser")
    parser.add_argument("-y", "--year", type=str, required=True, help="Year to load.")
    parser.add_argument("-xgb", "--xgb_name", type=str, required=True, help="XGBoost model name.")
    args = parser.parse_args()

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    Scales = ScalesCls(args.year) 

    sig_path = f"/cms/cephfs/data/store/user/castells/outputs/BDT/{args.year}/{args.xgb_name}/signal/"
    bkg_path = f"/cms/cephfs/data/store/user/castells/outputs/BDT/{args.year}/{args.xgb_name}/background/"
    data_path = f"/cms/cephfs/data/store/user/castells/outputs/BDT/{args.year}/{args.xgb_name}/data/"

    # Load signal MC for each mass point
    print("=== Loading samples (nominal only) ===")
    print("~~~ Signal samples ~~~")
    sig_samples = loadPerMass(sig_path, branches=["BDT_score", "LeadPs_interMass", "mass_gggg"])
    print("\n~~~ Background samples ~~~")
    bkg_samples = loadPerMass(bkg_path, branches=["BDT_score", "LeadPs_interMass", "mass_gggg"])
    print("\n~~~ Data samples ~~~")
    data_samples = loadPerMass(data_path, branches=["BDT_score", "LeadPs_interMass", "mass_gggg"])
    print()

    doRW = False
    sig_sb = False

    print("Signal from sideband only!" if sig_sb else "Signal from full region!")
    print("Keeping test (even) events only for signal and background.")
    for key in bkg_samples.keys():
        sig_samples[key].iloc[lambda x: x.index % 2 == 0]  # Even rows only (test set)
        bkg_samples[key].iloc[lambda x: x.index % 2 == 0]  # Even rows only (test set)
        if doRW:
            print(f"Reweighting background samples for {key.replace('_',' ')}:")
            bkg_samples[key] = reweight1Dim(bkg_samples[key], data_samples[key], "BDT_score")
    print()

    # Histogram setup
    bins = 30
    sig_bdts = {}
    bkg_bdts = {}
    data_bdts = {}
    sig_interMass1 = {}
    bkg_interMass1 = {}
    data_interMass1 = {}
    sig_interMass2 = {}
    bkg_interMass2 = {}
    data_interMass2 = {}
    masses = ["60_GeV", "45_GeV", "35_GeV", "25_GeV", "15_GeV"]
    #masses = ["15_GeV", "25_GeV", "35_GeV", "45_GeV", "60_GeV"]
    for mass in masses:
        if mass[:mass.find("GeV")+3] in masses:
            if sig_sb:
                sig_bdts.update({mass: fillHist("BDT_score", SideBand(sig_samples[mass]), [0, 1], binsScale=bins, normalize=False)})
            else:
                sig_bdts.update({mass: fillHist("BDT_score", sig_samples[mass], [0, 1], binsScale=bins, normalize=False)})
            bkg_bdts.update({mass: fillHist("BDT_score", SideBand(bkg_samples[mass]), [0, 1], binsScale=bins, normalize=False, customWeights="1dimreweight" if doRW else None)})
            data_bdts.update({mass: fillHist("BDT_score", SideBand(data_samples[mass]), [0, 1], binsScale=bins, normalize=False)})

            #sig_interMass1.update({mass: fillHist("LeadPs_interMass", SideBand(sig_samples[mass]), [-0.6, 0.6], binsScale=bins, normalize=False)})
            #bkg_interMass1.update({mass: fillHist("LeadPs_interMass", SideBand(bkg_samples[mass]), [-0.6, 0.6], binsScale=bins, normalize=False, customWeights="Ndimreweight" if "noRW" not in args.xgb_name else None)})
            #data_interMass1.update({mass: fillHist("LeadPs_interMass", SideBand(data_samples[mass]), [-0.6, 0.6], binsScale=bins, normalize=False)})

            #sig_interMass2.update({mass: fillHist("LeadPs_interMass", SideBand(sig_samples[mass]), [-0.6, 0.6], binsScale=bins, normalize=False)})
            #bkg_interMass2.update({mass: fillHist("LeadPs_interMass", SideBand(bkg_samples[mass]), [-0.6, 0.6], binsScale=bins, normalize=False, customWeights="Ndimreweight" if "noRW" not in args.xgb_name else None)})
            #data_interMass2.update({mass: fillHist("LeadPs_interMass", SideBand(data_samples[mass]), [-0.6, 0.6], binsScale=bins, normalize=False)})

            # Scale sig/bkg like normal to 1fb and to data, respectively
            # This is a slight nonsensical scaling to 1fb as is usually plotted since in the BDT training they need to be literally 1:1
            sig_bdts[mass].Scale(Scales.lumi_fb * 1.0 / Scales.sumw[args.year[:4]][mass])
            bkg_bdts[mass].Scale(data_bdts[mass].Integral() / bkg_bdts[mass].Integral())

            #sig_interMass1[mass].Scale(Scales.lumi_fb * 1.0 / Scales.sumw[args.year[:4]][mass])
            #sig_interMass2[mass].Scale(Scales.lumi_fb * 1.0 / Scales.sumw[args.year[:4]][mass])
            #bkg_interMass1[mass].Scale(data_interMass1[mass].Integral() / bkg_interMass1[mass].Integral())
            #bkg_interMass2[mass].Scale(data_interMass2[mass].Integral() / bkg_interMass2[mass].Integral())


    def plot(sig_hists, bkg_hists, data_hists, xTitle, title, path): 
        # Canvas setup
        canvas = TCanvas("canvas", "canvas", 1000, 1000)
        
        canvas.SetLogy()
        """
        pad1 = TPad("pad1" , "pad1", 0, 0.25, 1, 1)
        pad2 = TPad("pad2" , "pad2", 0, 0, 1, 0.23)

        pad1.SetLogy()
        pad1.SetBottomMargin(0.0001)
        pad1.SetBorderMode(0)
        pad2.SetTopMargin(0.01)
        pad2.SetBottomMargin(0.3)
        pad2.SetBorderMode(0)

        pad1.Draw()
        pad2.Draw()
        pad1.cd()
        """

        # Legend setup
        # NOTE: Make more legend to split signal, bkg, and data into something easier to read...
        legend_sig = TLegend(0.2,0.7,0.34,0.87)
        legend_sig.SetTextFont(42)
        legend_sig.SetBorderSize(0)
        legend_sig.SetFillStyle(0)
        legend_bkg = TLegend(0.35,0.7,0.49,0.87)
        legend_bkg.SetTextFont(42)
        legend_bkg.SetBorderSize(0)
        legend_bkg.SetFillStyle(0)
        #legend_data = TLegend(0.5,0.7,0.64,0.87)
        legend_data = TLegend(0.35,0.7,0.49,0.87)
        legend_data.SetTextFont(42)
        legend_data.SetBorderSize(0)
        legend_data.SetFillStyle(0)

        # Main plotting
        gStyle.SetTitleFont(42)
        gStyle.SetLabelFont(42)
        colors = [kGreen+3, kMagenta, kRed, kOrange, kBlue]
        shapes = [23, 22, 20, 47, 21]
        for (key, sig_hist), bkg_hist, data_hist, color, shape in zip(sig_hists.items(), bkg_hists.values(), data_hists.values(), colors, shapes):
            sig_hist.SetTitle("")
            sig_hist.GetYaxis().SetTitleFont(42)
            sig_hist.GetYaxis().SetTitle("Events")
            sig_hist.GetXaxis().SetTitleOffset(1.1)
            sig_hist.Sumw2()
            sig_hist.SetLineWidth(2)
            sig_hist.SetLineColor(color)
            
            bkg_hist.SetTitle("")
            bkg_hist.GetYaxis().SetTitleFont(42)
            bkg_hist.GetYaxis().SetTitle("Events")
            bkg_hist.GetXaxis().SetTitleOffset(1.1)
            bkg_hist.Sumw2()
            bkg_hist.SetLineWidth(2)
            bkg_hist.SetLineColor(color)
            #bkg_hist.SetFillColor(color)
            
            data_hist.SetTitle("")
            data_hist.GetYaxis().SetTitleFont(42)
            data_hist.GetYaxis().SetTitle("Events")
            data_hist.GetXaxis().SetTitleOffset(1.1)
            data_hist.Sumw2()
            data_hist.SetMarkerSize(1)
            data_hist.SetMarkerStyle(shape)
            data_hist.SetMarkerColor(color)
            data_hist.SetLineColor(color)

            sig_hist.GetXaxis().SetTitle(xTitle)
            bkg_hist.GetXaxis().SetTitle(xTitle)
            data_hist.GetXaxis().SetTitle(xTitle)
            
            min_order = -5
            max_order = 5
            sig_hist.GetYaxis().SetRangeUser(10**min_order + 10**(min_order-1), 10**max_order - 5*10**(max_order-1))
            bkg_hist.GetYaxis().SetRangeUser(10**min_order + 10**(min_order-1), 10**max_order - 5*10**(max_order-1))
            data_hist.GetYaxis().SetRangeUser(10**min_order + 10**(min_order-1), 10**max_order - 5*10**(max_order-1))

            bkg_hist.Draw("hist same")

        # Force draw order
        for (key, sig_hist), bkg_hist, data_hist, color in zip(sig_hists.items(), bkg_hists.values(), data_hists.values(), colors):
            sig_hist.Draw("hist same")

        # Force draw order
        """
        for (key, sig_hist), bkg_hist, data_hist, color in zip(sig_hists.items(), bkg_hists.values(), data_hists.values(), colors):
            if "35_GeV" in key and "interMass" not in path:
                data_hist.Draw("E1 P same")
            else:
                data_hist.Draw("E1 P same")
        """

        # Reverse order of legend
        for mass in [masses[i] for i in range(len(masses)-1, -1, -1)]:
            legend_sig.AddEntry(sig_hists[mass], "Signal " + mass.replace("_"," "), "l")
            legend_bkg.AddEntry(bkg_hists[mass], "Bkg " + mass.replace("_"," "), "f")
            """
            if "35_GeV" in mass and "interMass" not in path:
                legend_data.AddEntry(data_hists[mass], "Data " + mass.replace("_"," "), "lpe")
            else:
                legend_data.AddEntry(data_hists[mass], "Data " + mass.replace("_"," "), "lpe")
            """

        legend_sig.Draw("same")
        legend_bkg.Draw("same")
        #legend_data.Draw("same")

        """
        pad2.cd()

        # Residuals histogram setup
        residuals = []
        for mass in masses:
            residuals.append(residuals_error_prop(bkg_hists[mass], data_hists[mass], mass))

        # Plotting residuals on pad2
        pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
        pad2.cd()
        gStyle.SetTitleFont(42)
        gStyle.SetLabelFont(42)
        gPad.SetFrameBorderMode(0)
        gPad.SetTickx(1)
        gPad.SetTicky(1)

        for idx, (mass, rs, color, shape) in enumerate(zip(masses, residuals, colors, shapes)):
            rs.SetTitle("")
            rs.GetXaxis().SetTitle(xTitle)
            rs.GetYaxis().SetTitle("Bkg/Data")
            rs.GetXaxis().SetTitleSize(9.0*sig_hists["15_GeV"].GetYaxis().GetTitleSize() * pad1Scale)
            rs.GetYaxis().SetTitleSize(8.0*sig_hists["15_GeV"].GetYaxis().GetTitleSize() * pad1Scale)
            rs.GetXaxis().SetTitleOffset(1.1)
            rs.GetYaxis().SetTitleOffset(0.4)
            rs.GetYaxis().SetNdivisions(6)
            rs.GetYaxis().CenterTitle(True)
            #rs.GetYaxis().SetRangeUser(0.0, 2.0)
            rs.GetYaxis().SetRangeUser(rs.GetBinContent(rs.GetMinimumBin()) * (0.9 if rs.GetBinContent(rs.GetMinimumBin()) != 0 else rs.GetBinContent(rs.GetMinimumBin())), rs.GetBinContent(rs.GetMaximumBin()) * 1.105)
            rs.GetXaxis().SetLabelSize(2.5*rs.GetXaxis().GetLabelSize())
            rs.GetYaxis().SetLabelSize(2.5*rs.GetYaxis().GetLabelSize())
            rs.SetMarkerSize(1)
            rs.SetMarkerStyle(shape)
            rs.SetLineWidth(1)
            rs.SetLineColor(color)
            rs.SetMarkerColor(color)
            if mass == "35_GeV" and "interMass" not in path:
                rs.Draw("P same")
            else:
                rs.Draw("P same")

        # Line at y=1
        max_edge = residuals[0].GetBinCenter(residuals[0].GetNbinsX())
        min_edge = residuals[0].GetBinCenter(1)
        line = TLine(min_edge, 1.0, max_edge, 1.0)
        line.SetLineStyle(7)
        line.SetLineColor(kBlack)
        line.Draw("same")

        pad1.Modified()
        pad2.Modified()
        canvas.Update()
        """

        plotFancy(canvas, title, lumiTxt=Scales.lumiTxt, pad=False, prelim=True, inPlot=False, sub="Sideband Only" if sig_sb else "Bkg from SB", Scales=Scales, shiftLeft=0.05)
        #plotFancy(pad1, title, lumiTxt=Scales.lumiTxt, pad=True, prelim=True, inPlot=False, sub="Sideband Only" if sig_sb else "Bkg from SB", Scales=Scales, shiftLeft=0.05)
        canvas.SaveAs(path)
        assert os.path.exists(path), "Saved plot does not exist."


    plot(sig_bdts, bkg_bdts, data_bdts, "BDT Score (no RW)" if not doRW else "BDT Score (RW)", "BDT Score Comparison", os.path.join(*(cwd, "plots/BDT", args.year, "BDT_pred_generic", args.xgb_name, f"many_ma_comparisons.png")))
    #plot(sig_interMass1, bkg_interMass1, data_interMass1, "(m_{a1} - m_{Hyp}) / m_{#gamma#gamma#gamma#gamma}", "(m_{a1} - m_{Hyp}) / m_{#gamma#gamma#gamma#gamma}", os.path.join(*(cwd, "plots/BDT", args.year, "BDT_pred_generic", args.xgb_name, f"LeadPs_interMass_compare_multiple_ma.png")))
    #plot(sig_interMass2, bkg_interMass2, data_interMass2, "(m_{a2} - m_{Hyp}) / m_{#gamma#gamma#gamma#gamma}", "(m_{a2} - m_{Hyp}) / m_{#gamma#gamma#gamma#gamma}", os.path.join(*(cwd, "plots/BDT", args.year, "BDT_pred_generic", args.xgb_name, f"SubleadPs_interMass_compare_multiple_ma.png")))



