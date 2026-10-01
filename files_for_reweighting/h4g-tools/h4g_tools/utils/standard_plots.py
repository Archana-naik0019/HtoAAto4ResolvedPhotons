from ROOT import gROOT, TH1D, TCanvas, TLegend, kBlue, kRed, kOrange, kGreen, kYellow, kAzure, kBlack, gStyle, kRainBow
import os
import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from h4g_tools.utils.loading import fillHist, addInterMassAbs
from h4g_tools.utils.plotting import plotFancy

# Set batch so no plots pop up
gROOT.SetBatch(True)


def plotHMassCombined(
    samples: Dict[str, pd.DataFrame],
    branch: str,
    Scales: str = None
) -> Dict[str, Dict[str, TH1D]]:
    """
    Plots of Higgs mass distribution for 15, 30, 40, 60 GeV to compare to analysis note plots.
    """

    # Get directory of file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Fill histograms
    hist = {}
    valid = ["Signal_15_GeV", "Signal_30_GeV", "Signal_40_GeV", "Signal_60_GeV"]
    maximum = 0.0

    # Fill histograms
    hist = {}
    valid = ["Signal_15_GeV", "Signal_30_GeV", "Signal_40_GeV", "Signal_60_GeV"]
    for mass, sample in samples.items():
        # Only load desired masses into histograms
        fail = True
        for v in valid:
            if v in mass:
                fail = False
                adjusted_mass = v
        if fail:
            continue

        hist[adjusted_mass] = fillHist(branch, sample, [60.0, 200.0], binsScale=1.0/2.5, name=f"{branch}_{mass[7:13]}", normalize=True)

        # Get maximums of histograms
        if hist[adjusted_mass].GetMaximum() > maximum:
            maximum = hist[adjusted_mass].GetMaximum()
    
    canvasX = 1000
    canvasY = 1000
    canvas = TCanvas("canvas", "canvas", canvasX, canvasY)

    legend = TLegend(0.6,0.6,0.9,0.85)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    colors = [kBlue, kRed, kOrange, kGreen]
    for idx, (mass, color) in enumerate(zip(sorted(valid, reverse=True), colors)):
        # Histogram settings
        hist[mass].SetLineWidth(3)
        hist[mass].SetLineColor(color)

        hist[mass].SetTitle("")
        hist[mass].GetXaxis().SetTitle("m_{#gamma#gamma#gamma#gamma}")
        hist[mass].GetYaxis().SetTitle("Events / 2.5 GeV")
        hist[mass].GetYaxis().SetTitleOffset(1.6)
        hist[mass].GetXaxis().SetTitleOffset(1.2)
        hist[mass].GetYaxis().SetLabelOffset(0.01)
        hist[mass].GetXaxis().SetLabelOffset(0.01)
        hist[mass].SetAxisRange(0, maximum*1.1, "Y")
        
        hist[mass].Draw("hist same")

    for idx, mass in enumerate(valid):
        legend.AddEntry(hist[mass], f"m(a) = {mass.replace('Signal_',' ').replace('_',' ')}", "l")

    plotFancy(canvas, "", lumiTxt=Scales.defaultLumiTxt, prelim=True, inPlot=True, Scales=Scales, simulation=True)
    legend.Draw("same")

    # Save canvas
    canvas.SaveAs(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/combined_{branch}.pdf")
    assert os.path.exists(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/combined_{branch}.pdf")

    canvas.Destructor()
    return hist


def plotPsMassCombined(
    samples: Dict[str, pd.DataFrame],
    branches: List[str],
    Scales: str = None
) -> Dict[str, Dict[str, TH1D]]:
    """
    Plots of pseudoscalar masses for 15, 30, 40, 60 GeV to compare to analysis note plots.
    """

    # Get directory of file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Fill histograms
    hist = {}
    valid = ["Signal_15_GeV", "Signal_30_GeV", "Signal_40_GeV", "Signal_60_GeV"]
    for mass, sample in samples.items():
        # Only load desired masses into histograms
        fail = True
        for v in valid:
            if v in mass:
                fail = False
                adjusted_mass = v
        if fail:
            continue
        else:
            hist.update({adjusted_mass: {}})
        for branch in branches:
            hist[adjusted_mass].update({branch: fillHist(branch, sample, [0.0, 100.0], binsScale=0.5, name=f"{branch}_{mass[7:13]}", normalize=True, preventOverFlow=True)})
    
    canvasX = 1000
    canvasY = 1000
    canvas = TCanvas("canvas", "canvas", canvasX, canvasY)

    legend = TLegend(0.6,0.6,0.9,0.85)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    colors = [kGreen, kOrange, kRed, kBlue]
    for ps_num, branch in enumerate(branches):
        for idx, (mass, color) in enumerate(zip(valid, colors)):
            # Histogram settings
            hist[mass][branch].SetLineWidth(3)
            hist[mass][branch].SetLineColor(color)
            hist[mass][branch].SetTitle("")
            hist[mass][branch].GetXaxis().SetTitle("m_{a" + str(ps_num+1) + "}")
            hist[mass][branch].GetYaxis().SetTitle("Events / 2 GeV")
            hist[mass][branch].GetYaxis().SetTitleOffset(1.5)
            hist[mass][branch].GetXaxis().SetTitleOffset(1.2)
            hist[mass][branch].GetYaxis().SetRangeUser(0, 0.95)
            #hist[mass][branch].GetYaxis().SetRangeUser(0, 0.22)

            # Reset canvas if first in branch
            if idx == 0:
                hist[mass][branch].Draw("hist")
            else:
                hist[mass][branch].Draw("same hist")

            if ps_num == 0:
                header = "Leading Pseudoscalar"
            elif ps_num == 1:
                header = "Subleading Pseudoscalar"

            #legend.SetHeader(header)
            start = 7
            end = 9
            if "2018" in mass:
                start = 12
                end = 14
            legend.AddEntry(hist[mass][branch], f"m(a) = {mass[7:9]} GeV", "l")
        
        plotFancy(
            canvas,
            #"a_{1} Masses" if ps_num == 0 else "a_{2} Masses",
            "",
            lumiTxt=Scales.defaultLumiTxt,
            prelim=True,
            inPlot=True,
            Scales=Scales,
            simulation=True
        )
        legend.Draw("same")

        # Save canvas
        canvas.SaveAs(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/combined_{branch}.pdf")
        assert os.path.exists(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/combined_{branch}.pdf")

        # Clear legend for next branch
        legend.Clear()

    # Deal with wrapper saving canvas to heap (i.e. using new keyword in C++)
    canvas.Destructor()
    return hist


def massChecks(
    samples: Dict[str, pd.DataFrame],
    Scales: str = None
) -> Dict[str, Dict[str, TH1D]]:
    """
    Plots for psuedoscalar mass distributions.
    """

    # Get directory of file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Fill histograms
    hist = {}
    reco_hist = {}
    maximum = 0.0
    maximum_reco = 0.0
    gStyle.SetOptStat(0)
    for mass, sample in sorted(samples.items()):
        hist.update({mass: {}})
        reco_hist.update({mass: {}})
        genpart = sample.GenPart_mass.to_numpy()
        pdgID = sample.GenPart_pdgId.to_numpy()
        hist[mass] = TH1D(mass[7:13], mass[7:13], int((65.0 - 0.0)*1.0), 0.0, 65.0)
        for row_id, row_mass in zip(pdgID, genpart):
            assert len(row_id) == len(row_mass), (len(row_id), len(row_mass))
            for item_id, item_mass in zip(row_id, row_mass):
                if item_id == 35:
                    hist[mass].Fill(item_mass)
                
        hist[mass].Scale(1.0 / hist[mass].Integral())

        for branch in ["LeadPs_mass", "SubleadPs_mass"]:
            reco_hist[mass].update({branch: fillHist(branch, sample, [0.0, 65.0], name=f"{branch}_{mass[7:13]}", normalize=True)})

        # Get maximum of histograms
        if hist[mass].GetMaximum() > maximum:
            maximum = hist[mass].GetMaximum()
    
    canvasX = 1000
    canvasY = 1000
    canvas = TCanvas("canvas", "canvas", canvasX, canvasY)

    legend = TLegend(0.70,0.75,0.9,0.85)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    colors = [kGreen, kBlue]
    maximum = 0.0
    for idx, mass in enumerate(hist.keys()):
        # Plot setup
        hist[mass].SetLineWidth(3)
        hist[mass].SetTitle("")
        hist[mass].GetYaxis().SetTitle("Events (Normalized) / 1 GeV")
        hist[mass].GetXaxis().SetTitle("m_{a}")
        hist[mass].SetAxisRange(0, maximum*1.1, "Y")

        hist[mass].SetLineColor(colors[1])
        hist[mass].Draw("hist")

        canvas = plotFancy(canvas, f"m_{{a}} = {mass[7:13]} (gen)", lumiTxt=Scales.defaultLumiTxt, prelim=True, inPlot=False, Scales=Scales, simulation=True)

        # Save canvas to file
        canvas.SaveAs(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/{mass[7:13]}/{mass[7:13]}_combined_ps_genmass.pdf")
        assert os.path.exists(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/{mass[7:13]}/{mass[7:13]}_combined_ps_genmass.pdf")

        # Plot setup
        for idx, (branch, color) in enumerate(zip(["LeadPs_mass", "SubleadPs_mass"], colors)):
            reco_hist[mass][branch].SetLineWidth(3)
            reco_hist[mass][branch].SetTitle("")
            reco_hist[mass][branch].GetYaxis().SetTitle("Events (Normalized) / 1 GeV")
            reco_hist[mass][branch].GetXaxis().SetTitle("m_{a}")
            reco_hist[mass][branch].SetAxisRange(0, maximum_reco*1.1, "Y")

            reco_hist[mass][branch].SetLineColor(color)
            reco_hist[mass][branch].Draw("hist" if idx == 0 else "hist same")

            ltext = ""
            if "LeadPs" in branch:
                ltext = "m_{a1}"
            elif "SubleadPs" in branch:
                ltext = "m_{a2}"
            legend.AddEntry(reco_hist[mass][branch], ltext, "l")

        legend.Draw()
        canvas = plotFancy(canvas, f"m_{{a}} = {mass[7:13]} (reco)", lumiTxt=Scales.defaultLumiTxt, prelim=True, inPlot=False, Scales=Scales, simulation=True)

        # Save canvas to file
        canvas.SaveAs(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/{mass[7:13]}/{mass[7:13]}_combined_ps_recomass.pdf")
        assert os.path.exists(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/{mass[7:13]}/{mass[7:13]}_combined_ps_recomass.pdf")

        # Clear canvas for next mass
        canvas.Clear()
        legend.Clear()

    # Deal with wrapper saving canvas to heap (i.e. using new keyword in C++)
    canvas.Destructor()
    return hist


def plotPsMasses(
    samples: Dict[str, pd.DataFrame],
    branches: List[str],
    Scales: str = None
) -> Dict[str, Dict[str, TH1D]]:
    """
    Plots for psuedoscalar mass distributions.
    """

    # Get directory of file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Fill histograms
    hist = {}
    maximum = 0.0
    for mass, sample in sorted(samples.items()):
        hist.update({mass: {}})
        for branch in branches:
            hist[mass].update({branch: fillHist(branch, sample, [0.0, 100.0], name=f"{branch}_{mass[7:13]}", normalize=True)})

            # Get maximum of histograms
            if hist[mass][branch].GetMaximum() > maximum:
                maximum = hist[mass][branch].GetMaximum()
    
    canvasX = 1000
    canvasY = 1000
    canvas_even = TCanvas("canvas_even", "canvas_even", canvasX, canvasY)
    canvas_odd = TCanvas("canvas_odd", "canvas_odd", canvasX, canvasY)
    individual = TCanvas("individual", "individual", canvasX, canvasY)

    # Deal with legend
    legend_even = TLegend(0.7,0.62,0.97,0.87)
    legend_odd = TLegend(0.7,0.62,0.97,0.87)
    legend_even.SetBorderSize(0)
    legend_even.SetFillStyle(0)
    legend_odd.SetBorderSize(0)
    legend_odd.SetFillStyle(0)

    colors_raw = [kGreen, kYellow, kOrange, kRed, kBlue]
    colors = [color for color in colors_raw for _ in (0,1)] # 5x2 colors -> len == 10
    maximum = 0.0
    for branch in hist[list(hist.keys())[-1]].keys():
        # Set header of legend
        if "Lead" in branch:
            header = "Leading Pseudoscalar"
        elif "Sub" in branch:
            header = "Subleading Pseudoscalar"
        
        for idx, (mass, color) in enumerate(zip(hist.keys(), colors)):
            # Plot setup
            hist[mass][branch].SetLineWidth(3)
            hist[mass][branch].SetTitle("")
            hist[mass][branch].GetYaxis().SetTitle("Events (Normalized) / 1 GeV")
            hist[mass][branch].GetXaxis().SetTitle("m_{a1}  [GeV]" if "Lead" in branch else "m_{a2}  [GeV]")
            hist[mass][branch].SetAxisRange(0, maximum*1.1, "Y")
            hist[mass][branch].SetLineColor(color)

            # Separate even/odd masses
            start = 7
            end = 9
            if "2018" in mass:
                start = 12
                end = 14
            if int(mass[start:end]) % 2 == 0:
                canvas_even.cd()
                legend_even.AddEntry(hist[mass][branch], f"m(a) = {mass.replace('Signal_','').replace('_',' ')}", "l")
            elif int(mass[start:end]) % 2 != 0:
                canvas_odd.cd()
                legend_odd.AddEntry(hist[mass][branch], f"m(a) = {mass.replace('Signal_','').replace('_',' ')}", "l")
            
            # Reset canvas for each branch
            if idx == 0:
                hist[mass][branch].Draw("hist")
            else:
                hist[mass][branch].Draw("hist same")

            # Save individual plot for mass
            individual.cd()
            thist = hist[mass][branch].Clone()
            thist.SetLineColor(kBlue)
            thist.Draw("hist")
            plotFancy(individual, mass.replace("_"," "), lumiTxt=Scales.defaultLumiTxt, prelim=True, inPlot=False, Scales=Scales, simulation=True)
            individual.SaveAs(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/{mass[7:13]}/{branch}.pdf")

        canvas_even.cd()
        legend_even.SetHeader(header)
        legend_even.Draw("same")
        canvas_even = plotFancy(canvas_even, "Even Pseudoscalar Masses", lumiTxt=Scales.defaultLumiTxt, prelim=True, inPlot=False, Scales=Scales, simulation=True)

        canvas_odd.cd()
        legend_odd.SetHeader(header)
        legend_odd.Draw("same")
        canvas_odd = plotFancy(canvas_odd, "Odd Pseudoscalar Masses", lumiTxt=Scales.defaultLumiTxt, prelim=True, inPlot=False, Scales=Scales, simulation=True)

        # Save canvas to file
        canvas_even.SaveAs(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/pseudo_{branch}_even.pdf")
        canvas_odd.SaveAs(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/pseudo_{branch}_odd.pdf")
        assert os.path.exists(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/pseudo_{branch}_even.pdf"), "Even masses plot does not exist."
        assert os.path.exists(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/pseudo_{branch}_odd.pdf"), "Odd masses plot does not exist."

        # Clear legend for next branch
        legend_even.Clear()
        legend_odd.Clear()
        canvas_even.Clear()
        canvas_odd.Clear()

    # Deal with wrapper saving canvas to heap (i.e. using new keyword in C++)
    individual.Destructor()
    canvas_even.Destructor()
    canvas_odd.Destructor()
    return hist


def plot4Object(
    samples: Dict[str, pd.DataFrame],
    branches: List[str],
    Scales: str = None
) -> Dict[str, Dict[str, TH1D]]:
    """
    Plotting script for 4 photon object distributions.
    """

    # Get directory of file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Combine samples
    sample = pd.concat(samples.values())

    # Set variables
    hist = {}
    xtitle = {}
    bounds = {}
    bins = {}
    units = dict.fromkeys(branches, "GeV")
    for branch in branches:
        if branch == "mass_gggg":
            bins[branch] = 30.0/70.0
            bounds[branch] = [110.0, 180.0]
            xtitle[branch] = "m_{#gamma#gamma#gamma#gamma}"
        elif branch == "sumPt_gggg":
            bins[branch] = 1.0
            bounds[branch] = [50.0, 300.0]
            xtitle[branch] = "#sum{p_{T #gamma}} "
        elif branch == "dM":
            bins[branch] = 5.0
            bounds[branch] = [0.0, 25.0]
            xtitle[branch] = "#Delta M "
        elif branch == "dR_aa_mass_gggg":
            bins[branch] = 50.0
            bounds[branch] = [0.0, 1.0]
            units[branch] = "GeV^{-1}"
            xtitle[branch] = "#frac{#Delta R_{a1a2}}{m_{#gamma#gamma#gamma#gamma}} "
        elif branch == "cos_ag":
            bins[branch] = 15.0
            bounds[branch] = [-1.0, 1.0]
            units[branch] = ""
            xtitle[branch] = "cos(#theta_{#gamma a})"  
        elif branch == "LeadPs_interMass":
            bins[branch] = 20.0
            bounds[branch] = [-2.0, 2.0]
            units[branch] = ""
            xtitle[branch] = "#frac{m_{a1} - m_{Hyp}}{m_{#gamma#gamma#gamma#gamma}}"
        elif branch == "SubleadPs_interMass":
            bins[branch] = 20.0
            bounds[branch] = [-2.0, 2.0]
            units[branch] = ""
            xtitle[branch] = "#frac{m_{a2} - m_{Hyp}}{m_{#gamma#gamma#gamma#gamma}}"
        elif branch == "LeadPs_interMass_abs":
            bins[branch] = 20.0
            bounds[branch] = [-2.0, 2.0]
            units[branch] = ""
            xtitle[branch] = "#frac{|m_{a1} - m_{Hyp}|}{m_{#gamma#gamma#gamma#gamma}}"
        elif branch == "SubleadPs_interMass_abs":
            bins[branch] = 20.0
            bounds[branch] = [-2.0, 2.0]
            units[branch] = ""
            xtitle[branch] = "#frac{|m_{a2} - m_{Hyp}|}{m_{#gamma#gamma#gamma#gamma}}"
        elif branch == "Ps_massDiff":
            bins[branch] = 5.0
            bounds[branch] = [-10.0, 10.0]
            xtitle[branch] = "m_{a1} - m_{a2} "
        elif branch == "pT1_ma1":
            bins[branch] = 10.0
            bounds[branch] = [0.0, 10.0]
            units[branch] = ""
            xtitle[branch] = "p_{T}^{#gamma 1} / m_{a1}"
        elif branch == "pT2_ma1":
            bins[branch] = 10.0
            bounds[branch] = [0.0, 10.0]
            units[branch] = ""
            xtitle[branch] = "p_{T}^{#gamma 2} / m_{a1}"
        elif branch == "pT1_ma2":
            bins[branch] = 10.0
            bounds[branch] = [0.0, 10.0]
            units[branch] = ""
            xtitle[branch] = "p_{T}^{#gamma 1} / m_{a2}"
        elif branch == "pT2_ma2":
            bins[branch] = 10.0
            bounds[branch] = [0.0, 10.0]
            units[branch] = ""
            xtitle[branch] = "p_{T}^{#gamma 2} / m_{a2}"
        else:
            bins[branch] = 1.0
            bounds[branch] = [0.0, 100.0]
            units[branch] = ""
        
    # Fill histograms
    for mass, sample in samples.items():
        hist.update({mass: {}})
        for branch in branches:
            if branch == "LeadPs_interMass_abs" or branch == "LeadPs_interMass_abs":
                hist[mass].update({branch: fillHist(branch.replace("_abs",""), sample, bounds[branch], name=f"{branch}_{mass[7:13]}", binsScale=bins[branch], normalize=True)})
            hist[mass].update({branch: fillHist(branch, sample, bounds[branch], name=f"{branch}_{mass[7:13]}", binsScale=bins[branch], normalize=True)})

    # Canvas setup
    canvasX = 1000
    canvasY = 1000
    canvas = TCanvas("canvas", "canvas", canvasX, canvasY)

    for mass in samples.keys():
        for branch in branches:
            if branch in ["cos_ag", "cos_ag_no_boost"]:
                canvas.SetLogy()  # Log scale
            else:
                canvas.SetLogy(0)  # Linear scale

            # Plot setup
            hist[mass][branch].SetLineWidth(3)
            hist[mass][branch].SetTitle("")
            hist[mass][branch].GetYaxis().SetTitle(f"Events (Normalized) / {1.0 / bins[branch]:0.2f} {units[branch].replace('[','').replace(']','')}")
            hist[mass][branch].GetXaxis().SetTitle(f"{xtitle[branch]} [{units[branch]}]" if units[branch] != "" else f"{xtitle[branch]} {units[branch]}")
            hist[mass][branch].GetXaxis().SetTitleOffset(1.5)
            hist[mass][branch].Draw("hist")
            plotFancy(canvas, branch, sub="m_{a} = "+mass.replace("_"," "), lumiTxt=Scales.defaultLumiTxt, prelim=True, inPlot=False, Scales=Scales, simulation=True)

            # Save canvas
            canvas.SaveAs(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/{mass[7:13]}/{branch.replace('/','_')}.pdf")
            assert os.path.exists(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/{mass[7:13]}/{branch.replace('/','_')}.pdf"), f"{cwd}/../../scripts/plots/standard/{Scales.yr}/{mass[7:]}/{branch.replace('/','_')}.pdf"
    
    # Deal with wrapper saving canvas to heap (i.e. using new keyword in C++)
    canvas.Destructor()
    return hist


def plotPsKinematics(
    samples: Dict[str, pd.DataFrame],
    branches: List[str],
    Scales: str = None,
    skipPs: bool = False
) -> Dict[str, Dict[str, TH1D]]:
    """
    Plotting script for psuedoscalar kinematics distributions.
    """

    # Get directory of file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Set kinematic names
    names = []
    invalid = ["mvaID", "hoe", "sieie", "r9", "pfPhoIso03", "pfPhoIso03_rhoCorrected", "trkSumPtHollowConeDR03", "pixelSeed", "pfRelIso03_chg_quadratic", "pfRelIso_photon_pt"]
    diphoton = ["pT1_m_gg", "pT2_m_gg"]
    ps = ["leading_pho_pt", "subleading_pho_pt"]

    for branch in branches:
        for prefix in ["LeadPs", "SubleadPs", "pho1", "pho2", "pho3", "pho4"] if not skipPs else ["pho1", "pho2", "pho3", "pho4"]:
            if branch in invalid and (prefix == "LeadPs" or prefix == "SubleadPs"):
                continue
            elif branch in ps and ("pho" in prefix):
                continue
            elif branch == "dR_gg" and "pho" in prefix:
                continue
            elif branch in diphoton:
                continue

            names.append(f"{prefix}_{branch}")

        if diphoton[0] in branch or diphoton[1] in branch:
            names.append(branch)
            names.append(f"{branch}_EB")
            names.append(f"{branch}_EE")
            names.append(f"{branch}_2EB")
            names.append(f"{branch}_2EE")
            names.append(f"{branch}_1EB1EE")


    # Fill histogram
    hist = {}
    xtitle =  {}
    bounds = {}
    bins = {}
    units = dict.fromkeys(names, "[GeV]")
    for mass, sample in samples.items():
        hist.update({mass: {}})
        for name in names:
            if "pt" in name or "pT" in name:
                if "pT1_m_gg" in name:
                    bins[name] = 100.0/7.0
                    bounds[name] = [0.0, 7.0]
                    units[name] = ""
                    xtitle[name] = "p_{T}^{#gamma 1} / m_{#gamma #gamma}"

                    if "_EB" in name:
                        xtitle[name] = "p_{T}^{#gamma 1} / m_{#gamma #gamma} (#gamma_{1} EB)"
                    elif "_EE" in name:
                        xtitle[name] = "p_{T}^{#gamma 1} / m_{#gamma #gamma} (#gamma_{1} EE)"
                    elif "_2EB" in name:
                        xtitle[name] = "p_{T}^{#gamma 1} / m_{#gamma #gamma} (2EB)"
                    elif "_2EE" in name:
                        xtitle[name] = "p_{T}^{#gamma 1} / m_{#gamma #gamma} (2EE)"
                    elif "_1EB1EE" in name:
                        xtitle[name] = "p_{T}^{#gamma 1} / m_{#gamma #gamma} (1EB1EE)"
                elif "pT2_m_gg" in name:
                    bins[name] = 100.0/7.0
                    bounds[name] = [0.0, 7.0]
                    units[name] = ""
                    xtitle[name] = "p_{T}^{#gamma 2} / m_{#gamma #gamma}"

                    if "_EB" in name:
                        xtitle[name] = "p_{T}^{#gamma 2} / m_{#gamma #gamma} (#gamma_{2} EB)"
                    elif "_EE" in name:
                        xtitle[name] = "p_{T}^{#gamma 2} / m_{#gamma #gamma} (#gamma_{2} EE)"
                    elif "_2EB" in name:
                        xtitle[name] = "p_{T}^{#gamma 2} / m_{#gamma #gamma} (2EB)"
                    elif "_2EE" in name:
                        xtitle[name] = "p_{T}^{#gamma 2} / m_{#gamma #gamma} (2EE)"
                    elif "_1EB1EE" in name:
                        xtitle[name] = "p_{T}^{#gamma 2} / m_{#gamma #gamma} (1EB1EE)"
                elif "pfRelIso_photon_pt" in name:
                    bins[name] = 50.0/0.5
                    bounds[name] = [0.0, 0.5]
                    units[name] = ""
                    xtitle[name] = "pfRelIso03_chg_quadratic * p_{T}"
                elif "leading_pho_pt" in name:
                    bins[name] = 1.0
                    bounds[name] = [0.0, 200.0]
                    units[name] = ""
                    xtitle[name] = "p_{T}^{#gamma 1}"
                elif "subleading_pho_pt" in name:
                    bins[name] = 1.0
                    bounds[name] = [0.0, 200.0]
                    units[name] = ""
                    xtitle[name] = "p_{T}^{#gamma 1}"
                else:
                    bins[name] = 1.0
                    bounds[name] = [0.0, 200.0]
                    units[name] = "GeV"
                    xtitle[name] = "p_{T}"
            elif "eta" in name:
                bins[name] = 10.0
                bounds[name] = [-3.0, 3.0]
                units[name] = ""
                xtitle[name] = "#eta"
            elif "phi" in name:
                bins[name] = 2.0
                bounds[name] = [-5.0, 5.0]
                units[name] = ""
                xtitle[name] = "#phi"
            elif "mass" in name:
                bins[name] = 1.0
                bounds[name] = [0.0, 80.0]
                units[name] = ""
                xtitle[name] = "m"
            elif "mvaID" in name:
                bins[name] = 15.0
                bounds[name] = [-1.0, 1.0]
                units[name] = ""
                xtitle[name] = "#gamma_{" + name[3] + "}  mvaID"
            elif "dR_gg" in name:
                bins[name] = 2.0
                bounds[name] = [0.0, 25.0]
                units[name] = ""
                xtitle[name] = "#Delta R_{#gamma#gamma}"
            elif "hoe" in name:
                bins[name] = 500.0
                bounds[name] = [0.0, 0.1]
                units[name] = ""
                xtitle[name] = "H/E"
            elif "sieie" in name:
                bins[name] = 500.0
                bounds[name] = [0.0, 0.1]
                units[name] = ""
                xtitle[name] = "#sigma_{i#eta i #eta}"
            elif "r9" in name:
                bins[name] = 25.0
                bounds[name] = [0.0, 2.0]
                units[name] = ""
                xtitle[name] = "R_{9}"
            elif "pfPhoIso03_rhoCorrected" in name:
                bins[name] = 50.0/4.5
                bounds[name] = [0.0, 4.5]
                units[name] = ""
                xtitle[name] = "pfPhoIso03 #rho Corrected"
            elif "pfPhoIso03" in name:
                bins[name] = 50.0/4.5
                bounds[name] = [0.0, 4.5]
                units[name] = ""
                xtitle[name] = "pfPhoIso03"
            elif "trkSumPtHollowConeDR03" in name:
                bins[name] = 50.0/6.5
                bounds[name] = [0.0, 6.5]
                units[name] = ""
                xtitle[name] = "trkSumPtHollowConeDR03"
            elif "pixelSeed" in name:
                bins[name] = 1.0
                bounds[name] = [-0.5, 1.5]
                units[name] = ""
                xtitle[name] = "pixelSeed"
            elif "pfRelIso03_chg_quadratic" in name:
                bins[name] = 200.0/21.0
                bounds[name] = [0.0, 21.0]
                units[name] = ""
                xtitle[name] = "pfRelIso03_chg_quadratic"
            else:
                bins[name] = 1.0
                bounds[name] = [0.0, 100.0]
                units[name] = ""
            
            if "LeadPs" in name and "dR_gg" not in name:
                xtitle[name] += "_{a1}"
            elif "SubleadPs" in name and "dR_gg" not in name:
                xtitle[name] += "_{a2}"
            elif "pho" in name and not "mvaID" in name:
                xtitle[name] += "_{#gamma" + name[3] + "}"

            if "dR_gg" in name:
                if "LeadPs" in name:
                    if "LeadPs_gg_inEB" in sample.keys():
                        sample = sample.drop(columns=["LeadPs_gg_inEB"])
                    hist[mass][name] = fillHist(name, sample, bounds[name], name=f"{mass[7:13]}_{name}", binsScale=bins[name], normalize=True)
                elif "SubleadPs" in name:
                    if "SubleadPs_gg_inEB" in sample.keys():
                        sample = sample.drop(columns=["SubleadPs_gg_inEB"])
                    hist[mass][name] = fillHist(name, sample, bounds[name], name=f"{mass[7:13]}_{name}", binsScale=bins[name], normalize=True)
            elif "pT1_m_gg" in name or "pT2_m_gg" in name:
                # Split into EB/EE for photons
                if "_EB" in name:
                    if "pT1" in name:
                        hist[mass][name] = fillHist(name[:name.find("m_gg")+4], sample[sample.lead_isScEtaEB], bounds[name], name=f"{mass[7:13]}_{name}", binsScale=bins[name], normalize=False)
                    elif "pT2" in name:
                        hist[mass][name] = fillHist(name[:name.find("m_gg")+4], sample[sample.sublead_isScEtaEB], bounds[name], name=f"{mass[7:13]}_{name}", binsScale=bins[name], normalize=False)
                elif "_EE" in name:
                    if "pT1" in name:
                        hist[mass][name] = fillHist(name[:name.find("m_gg")+4], sample[sample.lead_isScEtaEE], bounds[name], name=f"{mass[7:13]}_{name}", binsScale=bins[name], normalize=False)
                    elif "pT2" in name:
                        hist[mass][name] = fillHist(name[:name.find("m_gg")+4], sample[sample.sublead_isScEtaEE], bounds[name], name=f"{mass[7:13]}_{name}", binsScale=bins[name], normalize=False)
                elif "_2EB" in name:
                    hist[mass][name] = fillHist(name[:name.find("m_gg")+4], sample[sample.lead_isScEtaEB & sample.sublead_isScEtaEB], bounds[name], name=f"{mass[7:13]}_{name}", binsScale=bins[name], normalize=False)
                elif "_2EE" in name:
                    hist[mass][name] = fillHist(name[:name.find("m_gg")+4], sample[sample.lead_isScEtaEE & sample.sublead_isScEtaEE], bounds[name], name=f"{mass[7:13]}_{name}", binsScale=bins[name], normalize=False)
                elif "_1EB1EE" in name:
                    hist[mass][name] = fillHist(name[:name.find("m_gg")+4], sample[(sample.lead_isScEtaEB & sample.sublead_isScEtaEE) | (sample.lead_isScEtaEE & sample.sublead_isScEtaEB)], bounds[name], name=f"{mass[7:13]}_{name}", binsScale=bins[name], normalize=False)
                else:
                    hist[mass][name] = fillHist(name, sample, bounds[name], name=f"{mass[7:13]}_{name}", binsScale=bins[name], normalize=False)
            elif "ma1" in name or "ma2" in name:
                hist[mass][name] = fillHist(name, sample, bounds[name], name=f"{mass[7:13]}_{name}", binsScale=bins[name], normalize=False)
            else:
                hist[mass][name] = fillHist(name, sample, bounds[name], name=f"{mass[7:13]}_{name}", binsScale=bins[name], normalize=True)
    
    canvasX = 1000
    canvasY = 1000
    canvas = TCanvas("canvas", "canvas", canvasX, canvasY)

    for mass in samples.keys():
        for name in names:
            # Set Logy status
            if "dR_gg" in name:
                canvas.SetLogy()
            else:
                canvas.SetLogy(0)

            # Plot setup
            hist[mass][name].SetLineWidth(3)
            hist[mass][name].SetTitle("")
            hist[mass][name].GetYaxis().SetTitle(f"Events" + ("" if ("pT1_m_gg" in name or "pT2_m_gg" in name or "ma1" in name or "ma2" in name) else " (Normalized)") + f" / {1.0 / bins[name]} {units[name].replace('[','').replace(']','')}")
            hist[mass][name].GetXaxis().SetTitle(f"{xtitle[name]} {units[name]}")
            hist[mass][name].GetXaxis().SetTitleOffset(1.2)
            hist[mass][name].Draw("hist")

            pT_mgg_extra = f" {hist[mass][name].Integral(0, hist[mass][name].GetNbinsX()+1)}" if ("pT1_m_gg" in name or "pT2_m_gg" in name) else ""
            plotFancy(canvas, xtitle[name], sub="m_{a} = "+mass.replace("Signal_","").replace("_"," ") + pT_mgg_extra, lumiTxt=Scales.defaultLumiTxt, prelim=True, inPlot=False, Scales=Scales, simulation=True)

            #"""
            if "mva" in name.lower():
                canvas.SetLogy()
                hist[mass][name].GetYaxis().SetRangeUser(10**-4 - 10**-5, hist[mass][name].GetMaximum() * 10**10)
            #"""

            # Save canvas
            canvas.SaveAs(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/{mass[7:13]}/{name}.pdf")
            assert os.path.exists(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/{mass[7:13]}/{name}.pdf")
    
    # Deal with wrapper saving canvas to heap (i.e. using new keyword in C++)
    canvas.Destructor()
    return hist


def plotMHyp(
    signal: pd.DataFrame,
    samples: Dict[str, pd.DataFrame], # For data, ps masses, and EM background
    branch: str,
    name: str,
    Scales: str = None
) -> Dict[str, TH1D]:
    """
    Plot of m_hyp for samples.
    """

    # Get directory of file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Fill histograms
    hists = {}
    signal_hists = {}
    for key, sample in samples.items():
        hists.update({key: fillHist(branch, sample, [10.0, 66.0], binsScale=1.0, name=f"m_{{hyp}}_{key}", normalize=False)})
    for key, sample in signal.items():
        #if key in ["Signal_15_GeV", "Signal_30_GeV", "Signal_40_GeV", "Signal_60_GeV"]:
        signal_hists.update({key: fillHist(branch, sample, [10.0, 66.0], binsScale=1.0, name=f"m_{{hyp}}_{key}", normalize=False)})
    
    scale = {"background": 0, "data": 0}
    for key, sample in signal.items():
        # Signal
        scale.update({key: Scales.lumi_fb * 1.0 / Scales.sumw[Scales.yr][key.replace("Signal_","")]})
    for key in samples.keys():
        # Background / Data
        if "background" in key:
            scale[key] = hists["data"].Integral() / hists["background"].Integral()
        elif "data" in key:
            scale["data"] = 1.0

    canvasX = 1000
    canvasY = 1000
    canvas = TCanvas("canvas", "canvas", canvasX, canvasY)
    canvas.SetLogy()

    # Legend setup
    legend_left = TLegend(0.52,0.71,0.7,0.79)
    legend_right = TLegend(0.7,0.55,0.9,0.8)
    for legend in [legend_left, legend_right]:
        legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
    
    # Histogram settings
    for key, hist in hists.items():
        hist.SetFillColor(kAzure-4)
        hist.SetTitle("")
        hist.GetXaxis().SetTitle("m_{Hyp}  [GeV]")
        hist.GetYaxis().SetTitle("Events / 1 GeV")
        hist.GetYaxis().SetTitleOffset(1.5)
        hist.GetXaxis().SetTitleOffset(1.2)
        hist.GetYaxis().SetRangeUser(10**-2, 10**10)
        hist.Print()
        hist.Scale(scale[key])

        if key == "background":
            hist.Draw("hist same")
            legend_left.AddEntry(hist, key.title(), "f")
        elif key == "data":
            hist.SetMarkerSize(1)
            hist.SetMarkerStyle(20)
            hist.SetMarkerColor(kBlack)
            hist.SetLineColor(kBlack)
            hist.Draw("P E1 same")
            legend_left.AddEntry(hist, key.title(), "lep")

    gStyle.SetPalette(kRainBow)
    for idx, (key, hist) in enumerate(signal_hists.items()):
    #for (key, hist), color in zip(signal_hists.items(), [kGreen+3, kOrange, kRed, kBlue]):
        #hist.SetLineColor(color)
        hists.update({key: hist})
        hist.SetLineWidth(3)
        hist.Scale(scale[key])
        if 2+idx < 6:
            legend_left.AddEntry(hist, key.replace("Signal_","").replace("_"," "), "l")
        if 2+idx > 6:
            legend_right.AddEntry(hist, key.replace("Signal_","").replace("_"," "), "l")
        hist.Draw("hist same")

    legend.Draw("same")
    plotFancy(canvas, "m_{Hyp} Distribution", lumiTxt=Scales.defaultLumiTxt, prelim=True, inPlot=True, Scales=Scales, simulation=True)

    # Save canvas
    canvas.SaveAs(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/{name}.pdf")
    assert os.path.exists(f"{cwd}/../../scripts/plots/standard/{Scales.yr}/{name}.pdf")

    # Deal with wrapper saving canvas to heap (i.e. using new keyword in C++)
    canvas.Destructor()
    return hists



def plotdRCombined(
    samples: Dict[str, pd.DataFrame],
    branch: str,
    pseudoscalar: int,
    Scales: str = None
) -> Dict[str, Dict[str, TH1D]]:
    """
    Plots of dR_gg distribution for 15, 30, 40, 60 GeV to compare to analysis note plots.
    """

    # Get directory of file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Fill histograms
    hist = {}
    valid = ["Signal_15_GeV", "Signal_30_GeV", "Signal_40_GeV", "Signal_60_GeV"]
    maximum = 0.0
    for mass, sample in samples.items():
        # Only load desired masses into histograms
        fail = True
        for v in valid:
            if v in mass:
                fail = False
                adjusted_mass = v
        if fail:
            continue
        else:
            hist.update({adjusted_mass: fillHist(branch, sample, [0.0, 5.0], binsScale=12.0, name=f"{branch}_{mass[7:13]}", normalize=True)})
            print(branch, mass, len(sample), len(sample[sample[branch] < 0.14]), len(sample[sample[branch] < 0.14]) / len(sample) * 100)

            # Get maximums of histograms
            if hist[adjusted_mass].GetMaximum() > maximum:
                maximum = hist[adjusted_mass].GetMaximum()
    
    canvasX = 1000
    canvasY = 1000
    canvas = TCanvas("canvas", "canvas", canvasX, canvasY)

    legend = TLegend(0.7,0.65,0.9,0.85)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    colors = [kBlue, kRed, kOrange, kGreen]
    for idx, (mass, color) in enumerate(zip(sorted(valid, reverse=True), colors)):
        # Histogram settings
        hist[mass].SetLineWidth(3)
        hist[mass].SetLineColor(color)

        hist[mass].SetTitle("")
        hist[mass].GetXaxis().SetTitle(f"a_{{{pseudoscalar}}} #Delta R_{{#gamma#gamma}}")
        hist[mass].GetYaxis().SetTitle("Events (Normalized)")
        hist[mass].GetYaxis().SetTitleOffset(1.6)
        hist[mass].GetXaxis().SetTitleOffset(1.2)
        hist[mass].GetYaxis().SetLabelOffset(0.01)
        hist[mass].GetXaxis().SetLabelOffset(0.01)
        hist[mass].SetAxisRange(0, maximum*1.1, "Y")
        
        hist[mass].Draw("hist same")
        legend.AddEntry(hist[mass], f"m(a) = {mass.replace('Signal_',' ').replace('_',' ')}", "l")

    plotFancy(canvas, f"a_{{{pseudoscalar}}} #Delta R_{{#gamma#gamma}}", lumiTxt=Scales.defaultLumiTxt, prelim=True, inPlot=False, Scales=Scales, simulation=True)
    legend.Draw("same")

    # Save canvas
    canvas.SaveAs(f"{cwd}/../../scripts/plots/BDT/{Scales.yr}/compare/a{pseudoscalar}_dR_gg_compare.pdf")
    assert os.path.exists(f"{cwd}/../../scripts/plots/BDT/{Scales.yr}/compare/a{pseudoscalar}_dR_gg_compare.pdf")

    canvas.Destructor()
    return hist
