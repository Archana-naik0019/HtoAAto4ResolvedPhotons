import pandas as pd
import xgboost
import os
import argparse
from h4g_tools.bdt.event_selection_bdt import reduce_dataframe, predict_BDT
from h4g_tools.utils.scales import Scales as ScalesCls, SideBand
from h4g_tools.utils.loading import loadPerMass, loadSamples, fillHist, pathSetup
from h4g_tools.utils.plotting import plotFancy, residuals_error_prop
from ROOT import TCanvas, TLegend, gStyle, kGreen, kMagenta, kBlue, kAzure, kRed, kOrange, kBlack, TPad, gPad, TLine, MakeNullPointer, TH1D, TPaveText


if __name__ == "__main__":

    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser")
    parser.add_argument("-y", "--year", type=str, required=True, help="Year to load.")
    parser.add_argument("-xgb1", "--xgb_name1", type=str, default=None, help="First XGBoost model name.")
    parser.add_argument("-xgb2", "--xgb_name2", type=str, default=None, help="Second XGBoost model name.")
    args = parser.parse_args()

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    Scales = ScalesCls(args.year) 

    pathSetup(args.year)

    model1 = {}
    model2 = {}
    for idx, (model, samples) in enumerate(zip([args.xgb_name1, args.xgb_name2], [model1, model2])):
        sig_path = f"{cwd}/outputs/BDT/{args.year}/{model}/signal/"
        bkg_path = f"{cwd}/outputs/BDT/{args.year}/{model}/background/"
        data_path = f"{cwd}/outputs/BDT/{args.year}/{model}/data/"

        # Load signal MC for each mass point
        print(f"=== Loading samples for Model {idx+1}: {model} (nominal only) ===")
        print("~~~ Signal samples ~~~")
        samples.update({"signal": loadPerMass(sig_path, branches=["weight", "BDT_score", "mass_gggg"])})
        print("\n~~~ Background samples ~~~")
        samples.update({"background": loadPerMass(bkg_path, branches=["BDT_score", "mass_gggg"])})
        print("\n~~~ Data samples ~~~")
        samples.update({"data": loadPerMass(data_path, branches=["BDT_score", "mass_gggg"])})
        print()

    # Histogram setup
    bins = 30
    sig_bdts = {"model1": {}, "model2": {}}
    bkg_bdts = {"model1": {}, "model2": {}}
    data_bdts = {"model1": {}, "model2": {}}
    masses = [f"{m}_GeV" for m in range(15,65,5)]
    for model, samples in zip(["model1", "model2"], [model1, model2]):
        for mass in masses:
            sig_bdts[model].update({mass: fillHist("BDT_score", SideBand(samples["signal"][mass]), [0, 1], binsScale=bins, normalize=False)})
            bkg_bdts[model].update({mass: fillHist("BDT_score", SideBand(samples["background"][mass]), [0, 1], binsScale=bins, normalize=False)})
            data_bdts[model].update({mass: fillHist("BDT_score", SideBand(samples["data"][mass]), [0, 1], binsScale=bins, normalize=False)})

            # Scale sig/bkg like normal to 1fb and to data, respectively
            # This is a slight nonsensical scaling of signal to bkg then to 1fb as is usually plotted, but in the BDT training they need to be literally 1:1
            sig_bdts[model][mass].Scale(bkg_bdts[model][mass].Integral() / sig_bdts[model][mass].Integral())
            sig_bdts[model][mass].Scale(Scales.lumi_fb * 1.0 / Scales.sumw[args.year[:4]][mass])
            bkg_bdts[model][mass].Scale(data_bdts[model][mass].Integral() / bkg_bdts[model][mass].Integral())

    print()

    def plot(sig_model1, sig_model2, bkg_model1, bkg_model2, data_model1, data_model2, model1, model2, mass, xTitle, title, path): 
        # Canvas setup
        canvas = TCanvas("canvas", "canvas", 1000, 1000)
        
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

        # Legend setup
        legend_model1 = TLegend(0.2,0.7,0.40,0.82)
        legend_model1.SetTextFont(42)
        legend_model1.SetBorderSize(0)
        legend_model1.SetFillStyle(0)
        legend_model2 = TLegend(0.45,0.7,0.65,0.82)
        legend_model2.SetTextFont(42)
        legend_model2.SetBorderSize(0)
        legend_model2.SetFillStyle(0)

        # Main plotting
        gStyle.SetTitleFont(42)
        gStyle.SetLabelFont(42)
        shapes = [24, 21]
        for sig_hist, bkg_hist, data_hist, shape in zip([sig_model1, sig_model2], [bkg_model1, bkg_model2], [data_model1, data_model2], shapes):
            sig_hist.SetTitle("")
            sig_hist.GetYaxis().SetTitleFont(42)
            sig_hist.GetYaxis().SetTitle("Events")
            sig_hist.GetXaxis().SetTitleOffset(1.1)
            sig_hist.Sumw2()
            sig_hist.SetLineWidth(2)
            
            bkg_hist.SetTitle("")
            bkg_hist.GetYaxis().SetTitleFont(42)
            bkg_hist.GetYaxis().SetTitle("Events")
            bkg_hist.GetXaxis().SetTitleOffset(1.1)
            bkg_hist.Sumw2()
            bkg_hist.SetLineWidth(2)
            
            data_hist.SetTitle("")
            data_hist.GetYaxis().SetTitleFont(42)
            data_hist.GetYaxis().SetTitle("Events")
            data_hist.GetXaxis().SetTitleOffset(1.1)
            data_hist.Sumw2()
            data_hist.SetMarkerSize(1)
            data_hist.SetMarkerStyle(shape)
            #data_hist.SetLineWidth(2)
            
            sig_hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, 10**7 - 5*10**6)
            bkg_hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, 10**7 - 5*10**6)
            data_hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, 10**7 - 5*10**6) 

        # Manually specify colors for plots
        colors = [kAzure-4, kRed]
        sig_model1.SetLineColor(kBlue)
        sig_model2.SetLineColor(kGreen+3)
        #bkg_model1.SetFillColor(colors[0])
        bkg_model1.SetLineColor(colors[0])
        #bkg_model2.SetFillColor(colors[1])
        bkg_model2.SetLineColor(colors[1])
        #data_model1.SetMarkerColor(colors[0])
        #data_model1.SetLineColor(colors[0])
        #data_model2.SetMarkerColor(colors[1])
        #data_model2.SetLineColor(colors[1])
        data_model1.SetMarkerColor(kBlack)
        data_model1.SetLineColor(kBlack)
        data_model2.SetMarkerColor(kBlack)
        data_model2.SetLineColor(kBlack)
        
        subtxt_model1 = TPaveText(0.17, 0.815, 0.42, 0.865, "NDC")
        subtxt_model1.SetTextFont(42)
        subtxt_model1.SetTextSize(0.025)
        subtxt_model1.SetTextColor(1)
        subtxt_model1.SetTextAlign(12)
        subtxt_model1.SetFillStyle(0)
        subtxt_model1.SetBorderSize(0)
        subtxt_model1.AddText("postEE Model w/ RWing")

        subtxt_model2 = TPaveText(0.42, 0.815, 0.67, 0.865, "NDC")
        subtxt_model2.SetTextFont(42)
        subtxt_model2.SetTextSize(0.025)
        subtxt_model2.SetTextColor(1)
        subtxt_model2.SetTextAlign(12)
        subtxt_model2.SetFillStyle(0)
        subtxt_model2.SetBorderSize(0)
        subtxt_model2.AddText("postEE Model w/out RWing")

        # Add objects to model1 legend
        legend_model1.AddEntry(sig_model1, "Signal " + mass.replace("_"," "), "l")
        legend_model1.AddEntry(bkg_model1, "Bkg " + mass.replace("_"," "), "f")
        legend_model1.AddEntry(data_model1, "Data " + mass.replace("_"," "), "lpe")
        
        # Add objects to model2 legend
        legend_model2.AddEntry(sig_model2, "Signal " + mass.replace("_"," "), "l")
        legend_model2.AddEntry(bkg_model2, "Bkg " + mass.replace("_"," "), "f")
        legend_model2.AddEntry(data_model2, "Data " + mass.replace("_"," "), "lpe")
        
        # Draw everything
        bkg_model2.Draw("hist")
        bkg_model1.Draw("hist same")
        data_model2.Draw("E1 P same")
        data_model1.Draw("E1 P same")
        sig_model2.Draw("hist same")
        sig_model1.Draw("hist same")
        legend_model1.Draw("same")
        legend_model2.Draw("same")
        subtxt_model1.Draw("same")
        subtxt_model2.Draw("same")

        pad2.cd()

        # Residuals histogram setup
        res_model1 = residuals_error_prop(bkg_model1, bkg_model2, mass)
        res_model2 = residuals_error_prop(data_model1, data_model2, mass)

        # Plotting residuals on pad2
        pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
        pad2.cd()
        gStyle.SetTitleFont(42)
        gStyle.SetLabelFont(42)
        gPad.SetFrameBorderMode(0)
        gPad.SetTickx(1)
        gPad.SetTicky(1)

        mx = -999
        mn = 999
        colors = [kRed, kBlack]
        shapes = [20, 24]
        for idx, (mass, rs, color, shape) in enumerate(zip(masses, [res_model2, res_model1], [colors[1], colors[0]], [shapes[1], shapes[0]])):
            rs.SetTitle("")
            rs.GetXaxis().SetTitle(xTitle)
            rs.GetYaxis().SetTitle("RWing / No RWing")
            rs.GetXaxis().SetTitleSize(9.0*sig_model1.GetYaxis().GetTitleSize() * pad1Scale)
            rs.GetYaxis().SetTitleSize(8.0*sig_model1.GetYaxis().GetTitleSize() * pad1Scale)
            rs.GetXaxis().SetTitleOffset(1.1)
            rs.GetYaxis().SetTitleOffset(0.4)
            rs.GetYaxis().SetNdivisions(6)
            rs.GetYaxis().CenterTitle(True)
            #rs.GetYaxis().SetRangeUser(0.0, 2.0)
            mx = rs.GetBinContent(rs.GetMaximumBin()) if rs.GetBinContent(rs.GetMaximumBin()) > mx else mx
            mn = rs.GetBinContent(rs.GetMinimumBin()) if rs.GetBinContent(rs.GetMinimumBin()) < mn else mn
            rs.GetYaxis().SetRangeUser(mn * 0.99, mx * 1.105)
            rs.GetXaxis().SetLabelSize(2.5*rs.GetXaxis().GetLabelSize())
            rs.GetYaxis().SetLabelSize(2.5*rs.GetYaxis().GetLabelSize())
            rs.SetMarkerSize(1)
            rs.SetMarkerStyle(shape)
            rs.SetLineWidth(1)
            rs.SetLineColor(color)
            rs.SetMarkerColor(color)
            rs.Draw("P same")

        # Line at y=1
        max_edge = res_model1.GetBinCenter(res_model1.GetNbinsX())
        min_edge = res_model2.GetBinCenter(1)
        line = TLine(min_edge, 1.0, max_edge, 1.0)
        line.SetLineStyle(7)
        line.SetLineColor(kBlack)
        line.Draw("same")

        pad1.Modified()
        pad2.Modified()
        canvas.Update()

        plotFancy(pad1, title, lumiTxt=Scales.lumiTxt, pad=True, prelim=True, inPlot=False, sub="Sideband Only", Scales=Scales, shiftLeft=0.05)
        canvas.SaveAs(path)
        assert os.path.exists(path), "Saved plot does not exist."


    for mass in masses:
        plot(
            sig_bdts["model1"][mass],
            sig_bdts["model2"][mass],
            bkg_bdts["model1"][mass],
            bkg_bdts["model2"][mass],
            data_bdts["model1"][mass],
            data_bdts["model2"][mass],
            args.xgb_name1,
            args.xgb_name2,
            mass,
            "BDT Score",
            "BDT Model Comparison",
            os.path.join(*(cwd, "plots/BDT", args.year, mass, f"compare_{args.xgb_name1}_{args.xgb_name2}.png"))
        )



