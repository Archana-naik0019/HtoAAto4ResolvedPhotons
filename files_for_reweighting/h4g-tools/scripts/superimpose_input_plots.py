import argparse
import os
import pandas as pd
from pathlib import Path
from ROOT import TH1D, TCanvas, TLegend, kBlack, kGreen, kRed, kBlue, TPad, gStyle, gPad, TLine
from h4g_tools.utils.scales import Scales as ScalesCls, SideBand, SignalRegion
from h4g_tools.bdt.event_selection_bdt import redefine_dataframe
from h4g_tools.utils.plotting import plotFancy, loadPlottingParameters, residuals_error_prop
from h4g_tools.utils.loading import loadPerMass, loadSamples, fillHist
from h4g_tools.bdt.reweight import loadNDimWeights

"""
python3 superimpose_input_plots.py -sig22 <signal-path-2022> -sig24 <signal-path-2024> -bkg22 <bkg-path-2022> -bkg24 <bkg-path-20242>

python3 superimpose_input_plots.py -sig22 /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2022_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSF_TriggerSF_oldXqcutQcut_20Dec2025/ -sig24 /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_14Nov2025/ -bkg22 /cms/cephfs/data/store/user/castells/bdtIO/BDT_2024_input_weights.parquet -bkg24 /cms/cephfs/data/store/user/castells/bdtIO/BDT_2022_input_weights.parquet
"""


# Get current directory
cwd = os.path.dirname(os.path.abspath(__file__))
scale_order = 4
scale = 10**scale_order


def includeMasses(
    branches: list,
    boundsList: list,
    binsScaleList: list,
    names: list,
    norms: list,
    titles: list,
) -> None:
    """
    Added masses to plotting variables.
    """
    
    branches.extend(["mass_gggg", "LeadPs_mass", "SubleadPs_mass"])
    boundsList.extend([[110.0, 180.0], [0.0, 65.0], [0.0, 65.0]])
    binsScaleList.extend([1.0, 1.0, 1.0])
    names.extend(["mass_gggg", "LeadPs_mass", "SubleadPs_mass"])
    norms.extend([False]*3)
    titles.extend(["m_{#gamma#gamma#gamma#gamma}", "m_{a1}", "m_{a2}"])



def plot_ratio(
    pad1: TPad,
    pad2: TPad,
    legend: TLegend,
    hist22: TH1D,
    hist24: TH1D,
    res_title: str,
    mass: int = None,
    res: TH1D = None,
    line: TLine = None,
    rng: tuple = None
) -> None:

    # Pad formatting
    pad1.SetBottomMargin(0.0001)
    pad1.SetBorderMode(0)
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(0.3)
    pad2.SetBorderMode(0)
    pad1.SetLogy()

    pad1.Draw()
    pad2.Draw()

    legend.Clear()
    legend.AddEntry(hist22, "Signal MC 2022" if mass is not None else "Background 2022")
    legend.AddEntry(hist24, "Signal MC 2024" if mass is not None else "Background 2024")

    pad1.cd()
    hist22.Draw("hist")
    hist24.Draw("hist same")
    legend.Draw()

    max_order = 6
    maximum = 10**max_order - 5*10**(max_order-1)
    if rng is None:
        hist22.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum*10)
        hist24.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum*10)
    else:
        hist22.GetYaxis().SetRangeUser(rng[0], rng[1])
        hist24.GetYaxis().SetRangeUser(rng[0], rng[1])

    # Plotting residuals on pad2
    pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
    pad2.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    gPad.SetFrameBorderMode(0)
    gPad.SetTickx(1)
    gPad.SetTicky(1)
    res.SetTitle("")
    res.GetXaxis().SetTitle(res_title)
    res.GetYaxis().SetTitle("2024/2022")
    res.GetXaxis().SetTitleSize(9.0*hist22.GetYaxis().GetTitleSize() * pad1Scale)
    res.GetYaxis().SetTitleSize(8.0*hist22.GetYaxis().GetTitleSize() * pad1Scale)
    res.GetXaxis().SetTitleOffset(1.1)
    res.GetYaxis().SetTitleOffset(0.4)
    res.GetYaxis().SetNdivisions(6)
    res.GetYaxis().CenterTitle(True)
    #res.GetYaxis().SetRangeUser(0.89, 1.11)
    res.GetYaxis().SetRangeUser(res.GetBinContent(res.GetMinimumBin()) * (0.9 if res.GetBinContent(res.GetMinimumBin()) != 0 else res.GetBinContent(res.GetMinimumBin())), res.GetBinContent(res.GetMaximumBin()) * 1.105)
    res.GetXaxis().SetLabelSize(2.5*res.GetXaxis().GetLabelSize())
    res.GetYaxis().SetLabelSize(2.5*res.GetYaxis().GetLabelSize())
    res.SetMarkerSize(1)
    res.SetMarkerStyle(20)
    res.SetMarkerColor(kBlack)
    res.SetLineColor(kBlack)
    res.Draw("P")

    # Line at y=1
    line.SetLineStyle(7)
    line.SetLineColor(kBlack)
    line.Draw("same")


def plot(
    canvas: TCanvas,
    legend: TLegend,
    hist22: TH1D,
    hist24: TH1D,
    mass: int = None,
) -> None:
    
    legend.Clear()
    legend.AddEntry(hist22, "Signal MC 2022" if mass is not None else "Background 2022")
    legend.AddEntry(hist24, "Signal MC 2024" if mass is not None else "Background 2024")

    hist22.Draw("hist")
    hist24.Draw("hist same")
    legend.Draw()

    max_order = 6
    maximum = 10**max_order - 5*10**(max_order-1)
    hist22.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum*10)
    hist24.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum*10)


def superimpose(
    args: argparse.ArgumentParser,
) -> None:

    # Path setup
    if not os.path.exists(os.path.join(cwd, "superimposed_2022_2024_bdt_inputs")):
        Path(os.path.join(cwd, "superimposed_2022_2024_bdt_inputs")).mkdir(parents=True, exist_ok=True)
    if not os.path.exists(os.path.join(cwd, "superimposed_2022_2024_bdt_inputs", "signal")):
        Path(os.path.join(cwd, "superimposed_2022_2024_bdt_inputs", "signal")).mkdir(parents=True, exist_ok=True)
    if not os.path.exists(os.path.join(cwd, "superimposed_2022_2024_bdt_inputs", "signal_ratio")):
        Path(os.path.join(cwd, "superimposed_2022_2024_bdt_inputs", "signal_ratio")).mkdir(parents=True, exist_ok=True)
    if not os.path.exists(os.path.join(cwd, "superimposed_2022_2024_bdt_inputs", "background")):
        Path(os.path.join(cwd, "superimposed_2022_2024_bdt_inputs", "background")).mkdir(parents=True, exist_ok=True)
    if not os.path.exists(os.path.join(cwd, "superimposed_2022_2024_bdt_inputs", "background_ratio")):
        Path(os.path.join(cwd, "superimposed_2022_2024_bdt_inputs", "background_ratio")).mkdir(parents=True, exist_ok=True)

    # Define BDT inputs
    branches, boundsList, binsScaleList, names, norms, titles = loadPlottingParameters(bdt=True)
    includeMasses(branches, boundsList, binsScaleList, names, norms, titles)
    paths = [os.path.join(cwd, "superimposed_2022_2024_bdt_inputs", "sample_type", f"{branch}.png") for branch in branches]

    # Canvas setup
    canvas = TCanvas("canvas", "canvas", 1000, 1000)
    canvas.SetLogy()

    # Legend setup
    legend = TLegend(0.65,0.75,0.87,0.85)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    if args.sig22 is not None and args.sig24 is not None:
        # Load signal samples
        print("Loading 2022 signal samples:")
        signal2022 = loadPerMass(args.sig22, branches = branches, eras = {"signal": [f"Signal_{m}_GeV_{prepost}" for m in range(15,65,5) for prepost in ["preEE", "postEE"]]})
        print("\nLoading 2024 signal samples:")
        signal2024 = loadPerMass(args.sig24, branches = branches)
        print()

        # Fill histograms for signal samples
        signal2022_hists = {}
        signal2024_hists = {}

        for branch, bounds, binsScale, name, norm, title, path in zip(branches, boundsList, binsScaleList, names, norms, titles, paths):
            for mass in range(15,65,5):
                # Set colors for 2022 to red and colors for 2024 to blue
                if "mass_gggg" == branch:
                    bounds = [110.0, 140.0]
                signal2022_hists.update({f"{branch}_{mass}_GeV": fillHist(branch, signal2022[f"Signal_{mass}_GeV"], bounds, binsScale=binsScale, name=f"{name}_{mass}_sig2022", normalize=norm, sb=False, preventOverFlow=True)})
                signal2024_hists.update({f"{branch}_{mass}_GeV": fillHist(branch, signal2024[f"Signal_{mass}_GeV"], bounds, binsScale=binsScale, name=f"{name}_{mass}_sig2024", normalize=norm, sb=False, preventOverFlow=True)})

                # Formatting 2022
                signal2022_hists[f"{branch}_{mass}_GeV"].SetTitle("")
                signal2022_hists[f"{branch}_{mass}_GeV"].GetYaxis().SetTitleFont(42)
                signal2022_hists[f"{branch}_{mass}_GeV"].GetYaxis().SetTitle("Events")
                signal2022_hists[f"{branch}_{mass}_GeV"].GetXaxis().SetTitle(title)
                signal2022_hists[f"{branch}_{mass}_GeV"].SetLineColor(kRed)
                signal2022_hists[f"{branch}_{mass}_GeV"].SetLineWidth(2)
                signal2022_hists[f"{branch}_{mass}_GeV"].Sumw2()

                # Formatting 2024
                signal2024_hists[f"{branch}_{mass}_GeV"].SetTitle("")
                signal2024_hists[f"{branch}_{mass}_GeV"].GetYaxis().SetTitleFont(42)
                signal2024_hists[f"{branch}_{mass}_GeV"].GetYaxis().SetTitle("Events")
                signal2024_hists[f"{branch}_{mass}_GeV"].GetXaxis().SetTitle(title)
                signal2024_hists[f"{branch}_{mass}_GeV"].SetLineColor(kBlue)
                signal2024_hists[f"{branch}_{mass}_GeV"].SetLineWidth(2)
                signal2024_hists[f"{branch}_{mass}_GeV"].Sumw2()

                # Set scales identically for both years
                signal2022_hists[f"{branch}_{mass}_GeV"].Scale(scale / signal2022_hists[f"{branch}_{mass}_GeV"].Integral())
                signal2024_hists[f"{branch}_{mass}_GeV"].Scale(scale / signal2024_hists[f"{branch}_{mass}_GeV"].Integral())

                # Plot and save canvas
                canvas.Clear()
                plot(canvas, legend, signal2022_hists[f"{branch}_{mass}_GeV"], signal2024_hists[f"{branch}_{mass}_GeV"], mass)
                plotFancy(canvas, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=f"Normalized to 10^{{{scale_order}}}", prelim=True, inPlot=False, simulation=True)
                canvas.SaveAs(path.replace("sample_type/", f"signal/{mass}_GeV_"))

                # Pad setup
                pad1 = TPad(f"pad1", "pad1", 0, 0.25, 1, 1)
                pad2 = TPad(f"pad2", "pad2", 0, 0, 1, 0.23)
                canvas.Clear()
                pad1.Draw()
                pad2.Draw()

                # Residuals histogram setup
                residuals = residuals_error_prop(signal2024_hists[f"{branch}_{mass}_GeV"], signal2022_hists[f"{branch}_{mass}_GeV"])
                max_edge = residuals.GetBinLowEdge(residuals.GetNbinsX()) + residuals.GetBinWidth(1)
                min_edge = residuals.GetBinLowEdge(1)
                line = TLine(min_edge, 1.0, max_edge, 1.0)

                # Plot with ratio and save
                plot_ratio(pad1, pad2, legend, signal2022_hists[f"{branch}_{mass}_GeV"], signal2024_hists[f"{branch}_{mass}_GeV"], title, mass, res=residuals, line=line, rng=None if "mass_gggg" not in path else (1.0, 10**5 - (5*10**4)))
                pad1.Modified()
                pad2.Modified()
                canvas.Update()
                plotFancy(pad1, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=f"Normalized to 10^{{{scale_order}}}", prelim=True, inPlot=False, simulation=True, pad=True)
                canvas.SaveAs(path.replace("sample_type/", f"signal_ratio/{mass}_GeV_"))
                print()

    if args.bkg22 is not None and args.bkg24 is not None:
        # Load bkg samples
        print("\nLoading 2022 background samples:")
        bkg2022 = loadSamples(args.bkg22, branches = branches+["Ndimreweight", "mass_gggg"])
        print("\nLoading 2024 background samples:")
        bkg2024 = loadSamples(args.bkg24, branches = branches+["Ndimreweight", "mass_gggg"])
        print()

        # Fill histograms for background samples
        bkg2022_hists = {}
        bkg2024_hists = {}

        for branch, bounds, binsScale, name, norm, title, path in zip(branches, boundsList, binsScaleList, names, norms, titles, paths):
            # Set colors for 2022 to red and colors for 2024 to blue
            bkg2022_hists.update({branch: fillHist(branch, bkg2022, bounds, binsScale=binsScale, name=f"{name}_bkg2022", normalize=norm, sb=True, customWeights="Ndimreweight")})
            bkg2024_hists.update({branch: fillHist(branch, bkg2024, bounds, binsScale=binsScale, name=f"{name}_bkg2024", normalize=norm, sb=True, customWeights="Ndimreweight")})

            # Formatting 2022
            bkg2022_hists[branch].SetTitle("")
            bkg2022_hists[branch].GetYaxis().SetTitleFont(42)
            bkg2022_hists[branch].GetYaxis().SetTitle("Events")
            bkg2022_hists[branch].GetXaxis().SetTitle(title)
            bkg2022_hists[branch].SetLineColor(kRed)
            bkg2022_hists[branch].SetLineWidth(2)
            bkg2022_hists[branch].Sumw2()

            # Formatting 2024
            bkg2024_hists[branch].SetTitle("")
            bkg2024_hists[branch].GetYaxis().SetTitleFont(42)
            bkg2024_hists[branch].GetYaxis().SetTitle("Events")
            bkg2024_hists[branch].GetXaxis().SetTitle(title)
            bkg2024_hists[branch].SetLineColor(kBlue)
            bkg2024_hists[branch].SetLineWidth(2)
            bkg2024_hists[branch].Sumw2()

            # Set scales identically for both years
            bkg2022_hists[branch].Scale(scale / bkg2022_hists[branch].Integral())
            bkg2024_hists[branch].Scale(scale / bkg2024_hists[branch].Integral())

            # Plot and save canvas
            canvas.Clear()
            plot(canvas, legend, bkg2022_hists[branch], bkg2024_hists[branch])
            plotFancy(canvas, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=f"Normalized to 10^{{{scale_order}}}", prelim=True, inPlot=False)
            canvas.SaveAs(path.replace("sample_type/", f"background/"))

            # Pad setup
            pad1 = TPad(f"pad1", "pad1", 0, 0.25, 1, 1)
            pad2 = TPad(f"pad2", "pad2", 0, 0, 1, 0.23)
            canvas.Clear()
            
            # Residuals histogram setup
            residuals = residuals_error_prop(bkg2024_hists[branch], bkg2022_hists[branch])
            max_edge = residuals.GetBinLowEdge(residuals.GetNbinsX()) + residuals.GetBinWidth(1)
            min_edge = residuals.GetBinLowEdge(1)
            line = TLine(min_edge, 1.0, max_edge, 1.0)

            # Plot with ratio and save
            plot_ratio(pad1, pad2, legend, bkg2022_hists[branch], bkg2024_hists[branch], title, res=residuals, line=line)
            pad1.Modified()
            pad2.Modified()
            canvas.Update()
            plotFancy(pad1, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=f"Normalized to 10^{{{scale_order}}}", prelim=True, inPlot=False, pad=True)
            canvas.SaveAs(path.replace("sample_type/", f"background_ratio/"))
            print()


if __name__ == "__main__":
    # Command-line arguments setup
    parser = argparse.ArgumentParser(description="Command line options parser", conflict_handler="resolve")

    # Input arguments
    parser.add_argument("-sig22", "--sig22", help="Directory to import signal samples for 2022.", type=str, default=None, required=False)
    parser.add_argument("-sig24", "--sig24", help="Directory to import signal samples for 2024.", type=str, default=None, required=False)
    parser.add_argument("-bkg22", "--bkg22", help="Directory to import background samples for 2022.", type=str, default=None, required=False)
    parser.add_argument("-bkg24", "--bkg24", help="Directory to import background samples for 2024.", type=str, default=None, required=False)
    args = parser.parse_args()

    assert (args.sig22 is not None and args.sig24 is not None) or (args.bkg22 is not None and args.bkg24 is not None)

    superimpose(args)
