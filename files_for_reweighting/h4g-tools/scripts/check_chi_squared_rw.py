import os
import gc
import json
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
from ROOT import TH1D, TCanvas, gROOT, kBlack, gPad, TPad, gStyle, kGray, TLine, kRed, TLegend, TGaxis, kGreen, kAzure
from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist
from h4g_tools.utils.plotting import loadPlottingParameters, plotFancy


# Get current directory
cwd = os.path.dirname(os.path.abspath(__file__))

bdt_inputs = [
    "pho1_mvaID",
    "pho2_mvaID",
    "pho3_mvaID",
    "pho4_mvaID",
    "LeadPs_pt",
    "SubleadPs_pt",
    "dR_aa_mass_gggg",
    "LeadPs_interMass",
    "SubleadPs_interMass",
    "Ps_massDiff",
    "cos_ag",
]


def doPlot(
    #canvas: TCanvas,
    #pad1: TPad,
    #pad2: TPad,
    scatter: TH1D,
    residuals: TH1D,
    rw0p1: TH1D,
    min1: TH1D,
    min2: TH1D,
    legend: TLegend,
    title: str,
    outdir: str,
) -> None:
 
    # Make canvas for plotting
    canvas = TCanvas("canvas", "canvas", 1000, 1000)
    canvas.SetLeftMargin(0.18)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.14)

    # Pad setup
    pad1 = TPad(f"pad1", "pad1", 0, 0.25, 1, 1)
    pad2 = TPad(f"pad2", "pad2", 0, 0, 1, 0.23)

    #pad1.SetBottomMargin(0.0001)
    top_bottom = 0.035
    pad1.SetBottomMargin(top_bottom)
    pad1.SetBorderMode(0)
    #pad2.SetTopMargin(0.01)
    pad2.SetTopMargin(top_bottom)
    pad2.SetBottomMargin(0.3)
    pad2.SetBorderMode(0)
    pad1.Draw()
    pad2.Draw()

    # Main plotting on pad1
    pad1.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)

    scatter.GetYaxis().SetTitle("#chi^{2}")
    scatter.GetYaxis().CenterTitle(True)
    scatter.SetTitle("")
    scatter.SetMarkerColor(kBlack)
    scatter.SetMarkerStyle(20)
    scatter.SetMarkerSize(1)
    rw0p1.SetTitle("")
    rw0p1.GetYaxis().SetTitle("")
    rw0p1.SetMarkerColor(kGreen)
    rw0p1.SetMarkerStyle(21)
    rw0p1.SetMarkerSize(1.5)
    min1.SetTitle("")
    min1.GetYaxis().SetTitle("")
    min1.SetMarkerColor(kRed)
    min1.SetMarkerStyle(21)
    min1.SetMarkerSize(1.5)

    scatter.GetYaxis().SetRangeUser(scatter.GetBinContent(scatter.GetMinimumBin()) - (scatter.GetBinContent(scatter.GetMaximumBin()) - scatter.GetBinContent(scatter.GetMinimumBin()))*0.1, scatter.GetBinContent(scatter.GetMaximumBin()) + (scatter.GetBinContent(scatter.GetMaximumBin()) - scatter.GetBinContent(scatter.GetMinimumBin()))*0.1)
    rw0p1.GetYaxis().SetRangeUser(scatter.GetBinContent(scatter.GetMinimumBin()) - (scatter.GetBinContent(scatter.GetMaximumBin()) - scatter.GetBinContent(scatter.GetMinimumBin()))*0.1, scatter.GetBinContent(scatter.GetMaximumBin()) + (scatter.GetBinContent(scatter.GetMaximumBin()) - scatter.GetBinContent(scatter.GetMinimumBin()))*0.1)
    min1.GetYaxis().SetRangeUser(scatter.GetBinContent(scatter.GetMinimumBin()) - (scatter.GetBinContent(scatter.GetMaximumBin()) - scatter.GetBinContent(scatter.GetMinimumBin()))*0.1, scatter.GetBinContent(scatter.GetMaximumBin()) + (scatter.GetBinContent(scatter.GetMaximumBin()) - scatter.GetBinContent(scatter.GetMinimumBin()))*0.1)

    scatter.Draw("hist p")
    #rw0p1.Draw("same hist p")
    min1.Draw("same hist p")
    legend.Draw()

    # Plotting residuals on pad2
    pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
    pad2.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    gPad.SetFrameBorderMode(0)
    gPad.SetTickx(1)
    gPad.SetTicky(1)

    # Make subplot - not actually residuals
    residuals.SetTitle("")
    residuals.GetYaxis().SetTitle("#chi^{2} / Min #chi^{2}")
    residuals.GetXaxis().SetTitleSize(9.0*scatter.GetYaxis().GetTitleSize() * pad1Scale)
    residuals.GetYaxis().SetTitleSize(8.0*scatter.GetYaxis().GetTitleSize() * pad1Scale)
    residuals.GetYaxis().SetTitleOffset(0.4)
    residuals.GetYaxis().CenterTitle(True)
    residuals.GetXaxis().SetLabelSize(2.5*residuals.GetXaxis().GetLabelSize())
    residuals.GetYaxis().SetLabelSize(2.5*residuals.GetYaxis().GetLabelSize())
    residuals.GetXaxis().SetLabelOffset(1.2)
    residuals.SetMarkerSize(1)
    residuals.SetMarkerStyle(20)
    residuals.SetMarkerColor(kBlack)
    #residuals.GetYaxis().SetRangeUser(-abs(residuals.GetBinContent(residuals.GetMaximumBin()))*0.10, residuals.GetBinContent(residuals.GetMaximumBin())*1.05)
    #residuals.GetYaxis().SetRangeUser(0, residuals.GetBinContent(residuals.GetMaximumBin())*1.05)
    residuals.Draw("P hist")
    min2.SetTitle("")
    min2.SetMarkerSize(1.5)
    min2.SetMarkerStyle(21)
    min2.SetMarkerColor(kRed)
    #min2.GetYaxis().SetRangeUser(-abs(residuals.GetBinContent(residuals.GetMaximumBin()))*0.10, residuals.GetBinContent(residuals.GetMaximumBin())*1.05)
    #min2.GetYaxis().SetRangeUser(0, residuals.GetBinContent(residuals.GetMaximumBin())*1.05)
    min2.Draw("same P hist")

    # Line at y=1
    max_edge = residuals.GetBinCenter(residuals.GetNbinsX()) + residuals.GetBinWidth(residuals.GetNbinsX())/2
    min_edge = residuals.GetBinLowEdge(1)
    line = TLine(min_edge, 1.0, max_edge, 1.0)
    line.SetLineStyle(7)
    line.SetLineColor(kBlack)
    line.Draw("same")

    pad1.Modified()
    pad2.Modified()
    canvas.Update()

    plotFancy(pad1, title.replace("[GeV]","").replace("[GeV^{-1}]","").strip(), prelim=True, pad=True, inPlot=False)#, shiftRight_all=0.05)
    canvas.SaveAs(f"{cwd}/plots/BDT/{outdir}.png")

    return 1


def draw(
    chi2: dict,
    era: str,
    reweightStyles: list,
) -> None:

    # Legend setup
    legend = TLegend(0.40,0.78,0.55,0.85)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    
    for varIndex, var in enumerate(bdt_inputs):

        # Handle filling here...
        legend.Clear()
        scatter = TH1D(f"scatter_{varIndex}", f"scatter_{varIndex}", len(reweightStyles), -1.0, len(reweightStyles))
        residuals = TH1D(f"residuals_{varIndex}", f"residuals_{varIndex}", len(reweightStyles), -1.0, len(reweightStyles))
        rw0p1 = TH1D(f"rw0p1_{varIndex}", f"rw0p1_{varIndex}", len(reweightStyles), -1.0, len(reweightStyles))
        rw0p1.SetBinContent(reweightStyles.index("0p1")+1, rw0p1.GetBinContent(reweightStyles.index("0p1")+1))
        min1 = TH1D(f"min1_{varIndex}", f"min1_{varIndex}", len(reweightStyles), -1.0, len(reweightStyles))
        min2 = TH1D(f"min2_{varIndex}", f"min2_{varIndex}", len(reweightStyles), -1.0, len(reweightStyles))
        for rs_orig in chi2.keys():
            if era not in chi2[rs_orig].keys():
                continue

            #rs = f"Default Bin = {float(rs_orig[rs_orig.find('W')+1:].replace('p','.'))}"
            rs = rs_orig
            scatter.Fill(rs, chi2[rs_orig][era][var][0])
            residuals.Fill(rs, chi2[rs_orig][era][var][0])
            rw0p1.Fill(rs, -999)
            min1.Fill(rs, -999)
            min2.Fill(rs, -999)

        for b in range(residuals.GetNbinsX()):
            #residuals.SetBinContent(b+1, residuals.GetBinContent(b+1) - scatter.GetBinContent(scatter.GetMinimumBin()))
            residuals.SetBinContent(b+1, residuals.GetBinContent(b+1) / scatter.GetBinContent(scatter.GetMinimumBin()))

        rw0p1.SetBinContent(reweightStyles.index("0p1")+1, scatter.GetBinContent(reweightStyles.index("0p1")+1))
        min1.SetBinContent(scatter.GetMinimumBin(), scatter.GetBinContent(scatter.GetMinimumBin()))
        min2.SetBinContent(scatter.GetMinimumBin(), 1.0)
        #legend.AddEntry(rw0p1, "8-dim w/ 0.1", "p")
        legend.AddEntry(min1, "Min #chi^{2}", "p")

        if scatter.GetEntries() == 0:
            continue

        # Do common plotting
        outdir = f"{era}/chi2/{var}"
        title = chi2[rs_orig][era][var][1]
        doPlot(
            scatter,
            residuals,
            rw0p1,
            min1,
            min2,
            legend,
            title,
            outdir
        )


def drawTotals(
    sums: dict,
    era: str,
    reweightStyles: list,
) -> None:

    if era not in sums.keys():
        return 0

    # Legend setup
    legend = TLegend(0.55,0.70,0.70,0.85)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    
    # Handle filling here...
    legend.Clear()
    scatter = TH1D(f"scatter", f"scatter", len(reweightStyles), -1.0, len(reweightStyles))
    residuals = TH1D(f"residuals", f"residuals", len(reweightStyles), -1.0, len(reweightStyles))
    rw0p1 = TH1D(f"rw0p1", f"rw0p1", len(reweightStyles), -1.0, len(reweightStyles))
    min1 = TH1D(f"min1", f"min1", len(reweightStyles), -1.0, len(reweightStyles))
    min2 = TH1D(f"min2", f"min2", len(reweightStyles), -1.0, len(reweightStyles))
    for rs in sums[era].keys():
        scatter.Fill(rs, sums[era][rs])
        residuals.Fill(rs, sums[era][rs])
        rw0p1.Fill(rs, -999)
        min1.Fill(rs, -999)
        min2.Fill(rs, -999)

    for b in range(residuals.GetNbinsX()):
        #residuals.SetBinContent(b+1, residuals.GetBinContent(b+1) - scatter.GetBinContent(scatter.GetMinimumBin()))
        residuals.SetBinContent(b+1, residuals.GetBinContent(b+1) / scatter.GetBinContent(scatter.GetMinimumBin()))

    rw0p1.SetBinContent(reweightStyles.index("0p1")+1, scatter.GetBinContent(reweightStyles.index("0p1")+1))
    min1.SetBinContent(scatter.GetMinimumBin(), scatter.GetBinContent(scatter.GetMinimumBin()))
    min2.SetBinContent(scatter.GetMinimumBin(), 1.0)
    #legend.AddEntry(rw0p1, "8-dim w/ 0.1", "p")
    legend.AddEntry(min1, "Min #chi^{2}", "p")

    # Do common plotting
    outdir = f"{era}/chi2/totals"
    title = f"Sum of #chi^{{2}} for {era}"
    doPlot(
        scatter,
        residuals,
        rw0p1,
        min1,
        min2,
        legend,
        title,
        outdir
    )


def do_chi_squared(
    args: argparse.ArgumentParser,
    reweightStyle: str,
    year: str,
) -> dict:

    # Load samples from bdtIO
    path = os.path.join("{cwd}/../../scripts/chi2_checks", args.year, f"BDT_{args.year}_input_weights_{reweightStyle}.parquet")
    if year == "2024":
        path = path.replace(".parquet", "_nBins10.parquet")
    bkg = pd.read_parquet(path)
    print(f"bkg path: {path}")
    data = pd.read_parquet(path.replace("weights_", "weights_data_"))
    print(f"data path: {path.replace('weights_', 'weights_data_')}")

    chi2_out = {}

    # Keep only the branches in "bdt_inputs"
    all_branches, all_boundsList, all_binsScaleList, all_names, all_norms, all_titles = loadPlottingParameters()
    branches, boundsList, binsScaleList, names, norms, titles = [], [], [], [], [], []
    for idx in range(len(all_branches)):
        if all_branches[idx] in bdt_inputs:
            branches.append(all_branches[idx])
            boundsList.append(all_boundsList[idx])
            binsScaleList.append(all_binsScaleList[idx])
            names.append(all_names[idx])
            norms.append(all_norms[idx])
            titles.append(all_titles[idx])

    # Make canvas for plotting
    canvas = TCanvas(f"canvas_{reweightStyle}", f"canvas_{reweightStyle}", 1000, 1000)
    canvas.SetLogy()

    era = args.year
    # Check Chi2 of the data/bkg for the BDT input variables
    print(f"~~~~    {era} {reweightStyle}    ~~~~")
    chi2_out.update({era: {}})
    for idx, var in enumerate(bdt_inputs):
        bkg_hist = fillHist(var, bkg, boundsList[idx], binsScale=binsScaleList[idx], name=f"{names[idx]}_bkg_rw_{era}_bkg", normalize=norms[idx], customWeights="Ndimreweight", sb=True)
        data_hist = fillHist(var, data, boundsList[idx], binsScale=binsScaleList[idx], name=f"{names[idx]}_data_{era}_data", normalize=norms[idx], sb=True)

        # Normalize bkg to data
        bkg_hist.Scale(data_hist.Integral() / bkg_hist.Integral())
        
        # Legend setup
        legend = TLegend(0.65,0.60,0.87,0.80)
        legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)

        # Formatting
        bkg_hist.SetFillColor(kAzure-4)
        legend.AddEntry(bkg_hist, "Event Mixing", "f")

        # Clone background to include statistical errors
        background_error = bkg_hist.Clone("background_error")
        background_error.SetLineColor(kGray+3)
        background_error.SetFillColor(kGray+3)
        background_error.SetFillStyle(3008)
        legend.AddEntry(background_error, "Statistical Uncertainty", "f")

        data_hist.SetMarkerSize(1)
        data_hist.SetMarkerStyle(20)
        data_hist.SetMarkerColor(kBlack)
        data_hist.SetFillStyle(0)
        legend.AddEntry(data_hist, "Data", "lep")

        max_order = 6
        maximum = 10**max_order - 5*10**(max_order-1)
        bkg_hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum)
        data_hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum)

        bkg_hist.Draw("hist")
        background_error.Draw("E2 same")
        data_hist.Draw("P E1 same")
        legend.Draw("same")

        # Plotting of histograms for checks
        path_inputs_dir = f"{cwd}/plots/BDT/{year}/chi2/{reweightStyle}/"
        path_inputs = f"{var}.png"
        if not os.path.exists(path_inputs_dir):
            Path(path_inputs_dir).mkdir(parents=True, exist_ok=True)
        canvas.SaveAs(os.path.join(path_inputs_dir, path_inputs))

        # Compute Chi2 using built-in ROOT function
        # NOTE: Look into normalized residuals from Chi2Test!
        chi2 = data_hist.Chi2Test(bkg_hist, "UW CHI2/NDF P")
        chi2_out[era].update({var: (chi2, titles[idx])})

        print()

    return chi2_out 




if __name__ == "__main__":
    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser", conflict_handler="resolve")
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    args = parser.parse_args()

    if args.year == "2022preEE":
        pass
    elif args.year == "2022postEE":
        pass
    elif args.year == "2022":
        pass
    elif args.year == "2023preBPix":
        pass
    elif args.year == "2023postBPix":
        pass
    elif args.year == "2023":
        pass
    elif args.year == "2024":
        pass
    else:
        assert False, "Must specify valid year!"
    # ADD A NEW YEAR HERE!

    pathSetup(args.year)

    # Turn off garbage collector
    gc.disable()

    chi2 = {}
    # Do 8-dim only!!!
    reweightStyles = ["1p0", "0p75", "0p50", "0p25", "0p1", "0p05", "0p025", "0p015"]
    for rs in reweightStyles:
        chi2.update({rs: do_chi_squared(args, rs, args.year)})

    # Print everything
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    sums = {args.year: dict.fromkeys(reweightStyles, 0.0).copy()}

    for rs in reweightStyles:
        for era in chi2[rs].keys():
            print(f"\t~~~  Chi2 Table for {era} using {rs} ~~~")
            for var in chi2[rs][era].keys():
                sums[era][rs] += chi2[rs][era][var][0]
                print(f"{var:<20}  {chi2[rs][era][var][0]:6.4f}")

        print()
    

    # Print total chi2 for each reweight scheme (is this a sensible number at all??)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    for era in sums.keys():
        print(f"\n\t~~~  Total Chi2 Table for {era} ~~~")
        for rs in sums[era].keys():
            print(f"{rs:<10}  {sums[era][rs]:6.4f}")

        print()

    # Save outputs to JSON
    print()
    with open("chi2_checks.json", "w") as f:
        print("Saved Chi2 test results to JSON!")
        json.dump(chi2, f, indent=4)


    # Plot chi2
    print(f"\nPlotting {args.year}...")
    draw(chi2, args.year, reweightStyles)
    drawTotals(sums, args.year, reweightStyles)


