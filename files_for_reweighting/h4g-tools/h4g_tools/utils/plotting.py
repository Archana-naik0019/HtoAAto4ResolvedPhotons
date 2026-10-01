from ROOT import gROOT, TGaxis, TFile, TCanvas, gPad, TPad, TLatex, TPaveText, TArrow, kBlack, TLegend, MakeNullPointer, TObject, TMath, TH1D, kAzure, kRed, kOrange, kGreen, kBlue, kGray, TGraph, gStyle, TLine, kMagenta, kYellow, kSpring, kCyan, kViolet, TGraphAsymmErrors, kPink, TH2D, TProfile, BindObject
import os
import json
from array import array
import subprocess
import shutil
from math import fsum
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
import seaborn as sns
from typing import List, Optional, Dict, Tuple
import h4g_tools.utils.loading as loading


# Set batch so no plots pop up
gROOT.SetBatch(True)

# No popups with matplotlib
matplotlib.use('Agg')



def loadPlottingParameters(
    bdt: bool = False
) -> Tuple:
    """
    Returns branches with corresponding parameters for plotting.
    """

    branches = [
        "pho1_mvaID",
        "pho2_mvaID",
        "pho3_mvaID",
        "pho4_mvaID",
        "LeadPs_pt",
        "SubleadPs_pt",
        "dR_aa_mass_gggg",
        "LeadPs_interMass",
        "SubleadPs_interMass",
        "LeadPs_interMass_abs",
        "SubleadPs_interMass_abs",
        "Ps_massDiff",
        "cos_ag",
        "m_hyp",
        "pT1_ma1",
        "pT1_ma2",
        "pT2_ma1",
        "pT2_ma2",
        "pho1_s4",
        "pho2_s4",
        "pho3_s4",
        "pho4_s4",
        "pho1_sieie",
        "pho2_sieie",
        "pho3_sieie",
        "pho4_sieie",
        "pho1_r9",
        "pho2_r9",
        "pho3_r9",
        "pho4_r9",
        "pho1_esEnergyOverRawE",
        "pho2_esEnergyOverRawE",
        "pho3_esEnergyOverRawE",
        "pho4_esEnergyOverRawE",
        "pho1_esEffSigmaRR",
        "pho2_esEffSigmaRR",
        "pho3_esEffSigmaRR",
        "pho4_esEffSigmaRR",
        "pho1_pfPhoIso03",
        "pho2_pfPhoIso03",
        "pho3_pfPhoIso03",
        "pho4_pfPhoIso03",
        "pho1_pfPhoIso03_rhoCorrected",
        "pho2_pfPhoIso03_rhoCorrected",
        "pho3_pfPhoIso03_rhoCorrected",
        "pho4_pfPhoIso03_rhoCorrected",
        "pho1_ecalPFClusterIso",
        "pho2_ecalPFClusterIso",
        "pho3_ecalPFClusterIso",
        "pho4_ecalPFClusterIso",
        "pho1_hcalPFClusterIso",
        "pho2_hcalPFClusterIso",
        "pho3_hcalPFClusterIso",
        "pho4_hcalPFClusterIso",
        "pho1_pfChargedIso",
        "pho2_pfChargedIso",
        "pho3_pfChargedIso",
        "pho4_pfChargedIso",
        "pho1_pfChargedIsoWorstVtx",
        "pho2_pfChargedIsoWorstVtx",
        "pho3_pfChargedIsoWorstVtx",
        "pho4_pfChargedIsoWorstVtx",
        "pho1_eta",
        "pho2_eta",
        "pho3_eta",
        "pho4_eta",
        "pho1_trkSumPtHollowConeDR03",
        "pho2_trkSumPtHollowConeDR03",
        "pho3_trkSumPtHollowConeDR03",
        "pho4_trkSumPtHollowConeDR03",
        "pho1_trkSumPtSolidConeDR04",
        "pho2_trkSumPtSolidConeDR04",
        "pho3_trkSumPtSolidConeDR04",
        "pho4_trkSumPtSolidConeDR04",
        "pho1_etaWidth",
        "pho2_etaWidth",
        "pho3_etaWidth",
        "pho4_etaWidth",
        "pho1_phiWidth",
        "pho2_phiWidth",
        "pho3_phiWidth",
        "pho4_phiWidth",
    ]

    boundsList = [
        [-1.0, 1.0],
        [-1.0, 1.0],
        [-1.0, 1.0],
        [-1.0, 1.0],
        [0.0, 300.0],
        [0.0, 200.0],
        [0.0, 60 * 10**-3],
        [-0.6, 0.6],
        [-0.6, 0.6],
        [0.0, 0.6],
        [0.0, 0.6],
        [-100.0, 100.0],
        [-1.0, 1.0],
        [10.0, 66.0],
        [0.0, 15.0],
        [0.0, 15.0],
        [0.0, 15.0],
        [0.0, 15.0],
        [0.0, 1.0],
        [0.0, 1.0],
        [0.0, 1.0],
        [0.0, 1.0],
        [0.0, 0.06],
        [0.0, 0.06],
        [0.0, 0.06],
        [0.0, 0.06],
        [0.0, 1.2],
        [0.0, 1.2],
        [0.0, 1.2],
        [0.0, 1.2],
        [0.0, 0.2],
        [0.0, 0.2],
        [0.0, 0.2],
        [0.0, 0.2],
        [0.0, 15.0],
        [0.0, 15.0],
        [0.0, 15.0],
        [0.0, 15.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [-3.0, 3.0],
        [-3.0, 3.0],
        [-3.0, 3.0],
        [-3.0, 3.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 10.0],
        [0.0, 0.06],
        [0.0, 0.06],
        [0.0, 0.06],
        [0.0, 0.06],
        [0.0, 0.1],
        [0.0, 0.1],
        [0.0, 0.1],
        [0.0, 0.1],
    ]

    # Want 30 bins each usually
    binsScaleList = [
        15.0,
        15.0,
        15.0,
        15.0,
        0.1,
        0.15,
        0.5 / 10**-3,
        1.0/0.04,
        1.0/0.04,
        1.0/0.02,
        1.0/0.02,
        0.15,
        15.0,
        1.0,
        2.0,
        2.0,
        2.0,
        2.0,
        100.0,
        100.0,
        100.0,
        100.0,
        100.0/0.06,
        100.0/0.06,
        100.0/0.06,
        100.0/0.06,
        100.0/1.2,
        100.0/1.2,
        100.0/1.2,
        100.0/1.2,
        100.0/0.2,
        100.0/0.2,
        100.0/0.2,
        100.0/0.2,
        100.0/15.0,
        100.0/15.0,
        100.0/15.0,
        100.0/15.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        100.0/6.0,
        100.0/6.0,
        100.0/6.0,
        100.0/6.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        100.0/0.06,
        100.0/0.06,
        100.0/0.06,
        100.0/0.06,
        1000.0,
        1000.0,
        1000.0,
        1000.0,
    ]

    names = [
        "pho1_mvaID",
        "pho2_mvaID",
        "pho3_mvaID",
        "pho4_mvaID",
        "LeadPs_pt",
        "SubleadPs_pt",
        "dR_aa_mass_gggg",
        "LeadPs_interMass",
        "SubleadPs_interMass",
        "LeadPs_interMass_abs",
        "SubleadPs_interMass_abs",
        "Ps_massDiff",
        "cos_ag",
        "m_Hyp",
        "pT1_ma1",
        "pT1_ma2",
        "pT2_ma1",
        "pT2_ma2",
        "pho1_s4",
        "pho2_s4",
        "pho3_s4",
        "pho4_s4",
        "pho1_sieie",
        "pho2_sieie",
        "pho3_sieie",
        "pho4_sieie",
        "pho1_r9",
        "pho2_r9",
        "pho3_r9",
        "pho4_r9",
        "pho1_esEnergyOverRawE",
        "pho2_esEnergyOverRawE",
        "pho3_esEnergyOverRawE",
        "pho4_esEnergyOverRawE",
        "pho1_esEffSigmaRR",
        "pho2_esEffSigmaRR",
        "pho3_esEffSigmaRR",
        "pho4_esEffSigmaRR",
        "pho1_pfPhoIso03",
        "pho2_pfPhoIso03",
        "pho3_pfPhoIso03",
        "pho4_pfPhoIso03",
        "pho1_pfPhoIso03_rhoCorrected",
        "pho2_pfPhoIso03_rhoCorrected",
        "pho3_pfPhoIso03_rhoCorrected",
        "pho4_pfPhoIso03_rhoCorrected",
        "pho1_ecalPFClusterIso",
        "pho2_ecalPFClusterIso",
        "pho3_ecalPFClusterIso",
        "pho4_ecalPFClusterIso",
        "pho1_hcalPFClusterIso",
        "pho2_hcalPFClusterIso",
        "pho3_hcalPFClusterIso",
        "pho4_hcalPFClusterIso",
        "pho1_pfChargedIso",
        "pho2_pfChargedIso",
        "pho3_pfChargedIso",
        "pho4_pfChargedIso",
        "pho1_pfChargedIsoWorstVtx",
        "pho2_pfChargedIsoWorstVtx",
        "pho3_pfChargedIsoWorstVtx",
        "pho4_pfChargedIsoWorstVtx",
        "pho1_eta",
        "pho2_eta",
        "pho3_eta",
        "pho4_eta",
        "pho1_trkSumPtHollowConeDR03",
        "pho2_trkSumPtHollowConeDR03",
        "pho3_trkSumPtHollowConeDR03",
        "pho4_trkSumPtHollowConeDR03",
        "pho1_trkSumPtSolidConeDR04",
        "pho2_trkSumPtSolidConeDR04",
        "pho3_trkSumPtSolidConeDR04",
        "pho4_trkSumPtSolidConeDR04",
        "pho1_etaWidth",
        "pho2_etaWidth",
        "pho3_etaWidth",
        "pho4_etaWidth",
        "pho1_phiWidth",
        "pho2_phiWidth",
        "pho3_phiWidth",
        "pho4_phiWidth",
    ]

    norms = [False] * len(branches)

    titles = [
        "#gamma_{1}  mvaID",
        "#gamma_{2}  mvaID",
        "#gamma_{3}  mvaID",
        "#gamma_{4}  mvaID",
        "p_{T}  a_{1}  [GeV]",
        "p_{T}  a_{2}  [GeV]",
        "#Delta R_{a1a2} / m_{#gamma#gamma#gamma#gamma}  [GeV^{-1}]",
        "(m_{a1} - m_{Hyp}) / m_{#gamma#gamma#gamma#gamma}",
        "(m_{a2} - m_{Hyp}) / m_{#gamma#gamma#gamma#gamma}",
        "|(m_{a1} - m_{Hyp})| / m_{#gamma#gamma#gamma#gamma}",
        "|(m_{a2} - m_{Hyp})| / m_{#gamma#gamma#gamma#gamma}",
        "m_{a1} - m_{a2}  [GeV]",
        "cos(#theta_{#gammaa})",
        "m_{Hyp}",
        "p_{T}^{#gamma_{1}} / m_{a1}",
        "p_{T}^{#gamma_{2}} / m_{a1}",
        "p_{T}^{#gamma_{1}} / m_{a2}",
        "p_{T}^{#gamma_{2}} / m_{a2}",
        "#gamma_{1}  E_{2x2} / E_{5x5}",
        "#gamma_{2}  E_{2x2} / E_{5x5}",
        "#gamma_{3}  E_{2x2} / E_{5x5}",
        "#gamma_{4}  E_{2x2} / E_{5x5}",
        "#gamma_{1}  #sigma_{i#eta i#eta}",
        "#gamma_{2}  #sigma_{i#eta i#eta}",
        "#gamma_{3}  #sigma_{i#eta i#eta}",
        "#gamma_{4}  #sigma_{i#eta i#eta}",
        "#gamma_{1}  R_{9}",
        "#gamma_{2}  R_{9}",
        "#gamma_{3}  R_{9}",
        "#gamma_{4}  R_{9}",
        "#gamma_{1}  E_{ES} / Raw E_{SC}",
        "#gamma_{2}  E_{ES} / Raw E_{SC}",
        "#gamma_{3}  E_{ES} / Raw E_{SC}",
        "#gamma_{4}  E_{ES} / Raw E_{SC}",
        "#gamma_{1}  ES #sigma_{eff RR}",
        "#gamma_{2}  ES #sigma_{eff RR}",
        "#gamma_{3}  ES #sigma_{eff RR}",
        "#gamma_{4}  ES #sigma_{eff RR}",
        "#gamma_{1}  PF Pho Iso (Cone 03)",
        "#gamma_{2}  PF Pho Iso (Cone 03)",
        "#gamma_{3}  PF Pho Iso (Cone 03)",
        "#gamma_{4}  PF Pho Iso (Cone 03)",
        "#gamma_{1}  PF Pho Iso (Cone 03) #rho Corr",
        "#gamma_{2}  PF Pho Iso (Cone 03) #rho Corr",
        "#gamma_{3}  PF Pho Iso (Cone 03) #rho Corr",
        "#gamma_{4}  PF Pho Iso (Cone 03) #rho Corr",
        "#gamma_{1}  PF ECAL Cluster Iso",
        "#gamma_{2}  PF ECAL Cluster Iso",
        "#gamma_{3}  PF ECAL Cluster Iso",
        "#gamma_{4}  PF ECAL Cluster Iso",
        "#gamma_{1}  PF HCAL Cluster Iso",
        "#gamma_{2}  PF HCAL Cluster Iso",
        "#gamma_{3}  PF HCAL Cluster Iso",
        "#gamma_{4}  PF HCAL Cluster Iso",
        "#gamma_{1}  PF Charged Iso (Chosen Vertex)",
        "#gamma_{2}  PF Charged Iso (Chosen Vertex)",
        "#gamma_{3}  PF Charged Iso (Chosen Vertex)",
        "#gamma_{4}  PF Charged Iso (Chosen Vertex)",
        "#gamma_{1}  PF Charged Iso (Worst Vertex)",
        "#gamma_{2}  PF Charged Iso (Worst Vertex)",
        "#gamma_{3}  PF Charged Iso (Worst Vertex)",
        "#gamma_{4}  PF Charged Iso (Worst Vertex)",
        "#gamma_{1}  #eta_{SC}",
        "#gamma_{2}  #eta_{SC}",
        "#gamma_{3}  #eta_{SC}",
        "#gamma_{4}  #eta_{SC}",
        "#gamma_{1}  PF Track Iso (Hollow Cone 03)",
        "#gamma_{2}  PF Track Iso (Hollow Cone 03)",
        "#gamma_{3}  PF Track Iso (Hollow Cone 03)",
        "#gamma_{4}  PF Track Iso (Hollow Cone 03)",
        "#gamma_{1}  PF Track Iso (Solid Cone 04)",
        "#gamma_{2}  PF Track Iso (Solid Cone 04)",
        "#gamma_{3}  PF Track Iso (Solid Cone 04)",
        "#gamma_{4}  PF Track Iso (Solid Cone 04)",
        "#gamma_{1}  SC #eta Width (#sigma_{#eta})",
        "#gamma_{2}  SC #eta Width (#sigma_{#eta})",
        "#gamma_{3}  SC #eta Width (#sigma_{#eta})",
        "#gamma_{4}  SC #eta Width (#sigma_{#eta})",
        "#gamma_{1}  SC #phi Width (#sigma_{#phi})",
        "#gamma_{2}  SC #phi Width (#sigma_{#phi})",
        "#gamma_{3}  SC #phi Width (#sigma_{#phi})",
        "#gamma_{4}  SC #phi Width (#sigma_{#phi})",
    ]

    if bdt:
        branches = branches[:13]
        boundsList = boundsList[:13]
        binsScaleList = binsScaleList[:13]
        names = names[:13]
        norms = norms[:13]
        titles = titles[:13]

        for remove_var in ["LeadPs_interMass_abs", "SubleadPs_interMass_abs"]:
            idx = branches.index(remove_var)
            del branches[idx]
            del boundsList[idx]
            del binsScaleList[idx]
            del names[idx]
            del norms[idx]
            del titles[idx]

    return branches, boundsList, binsScaleList, names, norms, titles


def plotFancy(
    canvas: TCanvas,
    title: str,
    prelim: bool,
    lumiTxt: Optional[str] = "",# = f"     {Scales.defaultLumiTxt}",
    sub: Optional[str] = None,
    inPlot: Optional[bool] = True,
    pad: Optional[bool] = False,
    adjust_sub_fit: Optional[bool] = False,
    Scales: str = None,
    colz: bool = False,
    shiftLeft = 0.0,
    simulation: bool = False,
    shiftRight_all = 0.0,
    overlay: bool = False,
) -> None:
    """
    Sets up CMS style canvas. Tested with 1000x1000 canvas.
    """

    # Select canvas
    canvas.cd()

    # General canvas settings
    gPad.SetFillColor(0)
    gPad.SetBorderMode(0)
    gPad.SetBorderSize(10)
    gPad.SetTickx(1)
    gPad.SetTicky(1)
    gPad.SetFrameFillStyle(0)
    gPad.SetFrameLineStyle(0)
    # gPad.SetFrameLineWidth(3)
    gPad.SetFrameBorderMode(0)
    gPad.SetFrameBorderSize(10)

    if not pad:
        canvas.SetLeftMargin(0.14)
        canvas.SetRightMargin(0.05)
        canvas.SetBottomMargin(0.14)
    if colz:
        canvas.SetRightMargin(0.15)

    # Latex settings
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextFont(43)
    latex.SetTextSize(20)
    latex.SetTextAlign(31)
    latex.SetTextAlign(11)

    low = {"lumi": [], "cms": [], "title": [], "prelim": []}
    textSize = {"lumi": 0, "cms": 0, "title": 0, "prelim": 0}
    if inPlot:
        low["lumi"] = [0.72 if not overlay else 0.58, 0.818]
        low["title"] = [0.18, 0.70]
        low["sub"] = [0.68 if not overlay else 0.68-0.1, 0.75]
        low["cms"] = [0.182, 0.74]
        low["prelim"] = [0.27, 0.745]
        textSize["lumi"] = 0.03
        textSize["cms"] = 0.04
        textSize["title"] = 0.025
        textSize["prelim"] = 0.031
        textSize["sub"] = 0.05*0.5
    else:
        low["title"] = [0.36+shiftRight_all, 0.805]
        low["sub"] = [0.68, 0.755]
        if not pad:
            low["prelim"] = [0.22, 0.815]
            low["cms"] = [0.13, 0.807]
            low["lumi"] = [0.68-shiftLeft, 0.818]
        else:
            low["prelim"] = [0.19+shiftRight_all, 0.815]
            low["cms"] = [0.10+shiftRight_all, 0.807]
            low["lumi"] = [0.71, 0.818]

        textSize["lumi"] = 0.03
        textSize["cms"] = 0.038
        textSize["title"] = 0.025
        textSize["prelim"] = 0.03
        textSize["sub"] = 0.05*0.5

    if adjust_sub_fit:
        low["sub"][0] += 0.1

    # Write lumi text to canvas
    lowX = low["lumi"][0]
    lowY = low["lumi"][1]
    lumi = TPaveText(lowX,lowY, lowX+0.3, lowY+0.2, "NDC")
    lumi.SetTextFont(42)
    lumi.SetTextSize(textSize["lumi"])
    lumi.SetTextColor(1)
    lumi.SetTextAlign(12)
    lumi.SetFillStyle(0)
    lumi.SetBorderSize(0)
    lumi.AddText(lumiTxt)
    lumi.DrawClone("same")

    # Write CMS text to canvas
    lowX = low["cms"][0]
    lowY = low["cms"][1]
    cmstxt = TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    cmstxt.SetTextFont(61)
    cmstxt.SetTextSize(textSize["cms"])
    cmstxt.SetTextColor(1)
    cmstxt.SetTextAlign(12)
    cmstxt.SetFillStyle(0)
    cmstxt.SetBorderSize(0)
    cmstxt.AddText("CMS")
    cmstxt.DrawClone("same")

    # Write title text to canvas
    lowX = low["title"][0]
    lowY = low["title"][1]
    samplestxt = TPaveText(lowX, lowY+0.06, lowX+0.3, lowY+0.16, "NDC")
    samplestxt.SetTextFont(42)
    samplestxt.SetTextSize(textSize["title"])
    samplestxt.SetTextColor(1)
    samplestxt.SetTextAlign(12)
    samplestxt.SetFillStyle(0)
    samplestxt.SetBorderSize(0)
    samplestxt.AddText(title)
    samplestxt.DrawClone("same")

    #prelim_str = "Preliminary"
    prelim_str = "Private Work"
    if simulation:
        prelim_str = "Simulation"
    if prelim and not simulation:
        # Write prelim text to canvas
        lowX = low["prelim"][0]
        lowY = low["prelim"][1]
        pretxt = TPaveText(lowX, lowY+0.05, lowX+0.15, lowY+0.15, "NDC")
        pretxt.SetTextFont(52)
        pretxt.SetTextSize(textSize["prelim"])
        pretxt.SetTextColor(1)
        pretxt.SetTextAlign(12)
        pretxt.SetFillStyle(0)
        pretxt.SetBorderSize(0)
        pretxt.AddText(prelim_str)
        pretxt.DrawClone("same")
    elif simulation:
        # Write prelim text to canvas
        lowX = low["prelim"][0]
        lowY = low["prelim"][1]
        pretxt = TPaveText(lowX, lowY+0.05, lowX+0.15, lowY+0.15, "NDC")
        pretxt.SetTextFont(52)
        pretxt.SetTextSize(textSize["prelim"])
        pretxt.SetTextColor(1)
        pretxt.SetTextAlign(12)
        pretxt.SetFillStyle(0)
        pretxt.SetBorderSize(0)
        pretxt.AddText(prelim_str)
        pretxt.DrawClone("same")

    if sub is not None:
        # Write sub text to canvas
        lowX = low["sub"][0]
        lowY = low["sub"][1]
        subtxt = TPaveText(lowX, lowY+0.05, lowX+0.15, lowY+0.15, "NDC")
        subtxt.SetTextFont(42)
        subtxt.SetTextSize(textSize["sub"])
        subtxt.SetTextColor(1)
        subtxt.SetTextAlign(12)
        subtxt.SetFillStyle(0)
        subtxt.SetBorderSize(0)
        subtxt.AddText(sub)
        subtxt.DrawClone("same")
    
    # Redraw axes ticks
    gPad.RedrawAxis()

    return canvas


def residuals_error_prop(
    numer: TH1D,
    denom: TH1D,
    var: str = None,
    corrections: bool = False,
    autoScale: bool = True,
    returnNewScale: bool = False,
) -> TH1D:
    """
    Calculates statistical error on residuals (ratio) plots based on input histograms using error propagation.
    """

    # Calculate Sumw2 as required: https://root.cern.ch/doc/master/classTH1.html#ac782a09c31b4f7de40f8fd4f77efa090
    numer.Sumw2()
    denom.Sumw2()

    # Clone numerator into residuals histogram
    residuals = numer.Clone()

    # Divide by denominator in place using binomial error
    # See https://root.cern.ch/doc/v632/TH1_8cxx_source.html#l03019
    # and https://root-forum.cern.ch/t/how-to-calculate-binomial-efficiency-error-with-weights/3650/5
    # for details on why Divide(..., "B") works!
    # Produces symmetric errors. Asymmetric errors require TGraphAsymErrors: https://root.cern.ch/doc/master/classTGraphAsymmErrors.html#a37a202762b286cf4c7f5d34046be8c0b
    # Also here: https://root.cern.ch/doc/master/classTGraphAsymmErrors.html#ac9a2403d1297546c603f5cf1511a5ca5
    residuals.Divide(residuals, denom, 1.0, 1.0, "B")

    binom = True
    # Use error propagation formula for errors instead of binomial errors if binom is set to False
    # Want error_res = |binContent_res| * sqrt( [var(sum of squared weights)_numer / binContent_numer^2 + var(sum of squared weights)_denom / binContent_denom^2] )
    if not binom:
        for b in range(residuals.GetNbinsX()):
            # Account for divide by zero issues
            if numer.GetBinContent(b) < 0.001 or denom.GetBinContent(b) < 0.001:
                continue

            residuals.SetBinError(
                b,
                TMath.Abs(residuals.GetBinContent(b)) * TMath.Sqrt(numer.GetBinError(b) / TMath.Sq(numer.GetBinContent(b)) + denom.GetBinError(b) / TMath.Sq(denom.GetBinContent(b)))
            )

    ylow = 0.89
    yhigh = 1.11
    if var is not None:
        if "#gamma_{1}  mvaID" in var:
            ylow, yhigh = (0.59, 1.41) if not corrections else (0.79, 1.11)
        elif "#gamma_{2}  mvaID" in var:
            ylow, yhigh = (0.59, 1.41) if not corrections else (0.79, 1.11)
        elif "#gamma_{3}  mvaID" in var:
            ylow, yhigh = (0.59, 6.11) if not corrections else (0.79, 1.11)
        elif "#gamma_{4}  mvaID" in var:
            ylow, yhigh = (0.59, 6.11) if not corrections else (0.79, 1.11)
        elif "p_{T}  a_{1}" in var:
            ylow, yhigh = (0.09, 5.11) if not corrections else (0.79, 1.11)
        elif "p_{T}  a_{2}" in var:
            ylow, yhigh = (0.29, 5.11) if not corrections else (0.79, 1.11)
        elif "#Delta R_{a1a2} / m_{#gamma#gamma#gamma#gamma}" in var:
            ylow, yhigh = (0.09, 2.01) if not corrections else (0.79, 1.11)
        elif "(m_{a1} - m_{Hyp}) / m_{#gamma#gamma#gamma#gamma}" in var:
            ylow, yhigh = (0.01, 1.51)
        elif "(m_{a2} - m_{Hyp}) / m_{#gamma#gamma#gamma#gamma}" in var:
            ylow, yhigh = (0.01, 1.51)
        elif "|(m_{a1} - m_{Hyp})| / m_{#gamma#gamma#gamma#gamma}" in var:
            ylow, yhigh = (0.39, 1.51)
        elif "|(m_{a2} - m_{Hyp})| / m_{#gamma#gamma#gamma#gamma}" in var:
            ylow, yhigh = (0.49, 2.31)
        elif "m_{a1} - m_{a2}" in var:
            ylow, yhigh = (0.09, 1.41) if not corrections else (0.79, 1.21)
        elif "cos" in var:
            ylow, yhigh = (0.49, 2.01)

    #mx = -999
    #mn = 999
    #mx = residuals.GetBinContent(residuals.GetMaximumBin()) if residuals.GetBinContent(residuals.GetMaximumBin()) > mx else mx
    #mn = residuals.GetBinContent(residuals.GetMinimumBin()) if residuals.GetBinContent(residuals.GetMinimumBin()) < mn else mn
    #residuals.GetYaxis().SetRangeUser(mn * 0.99, mx * 1.105)
    if autoScale:
        residuals.GetYaxis().SetRangeUser(ylow,yhigh)

    if not returnNewScale:
        return residuals
    else:
        return residuals, ylow, yhigh


def plot_test_train(
    results: Dict[str, pd.DataFrame],
    suppressSave: bool = False,
    year: str = None,
    xgb_name: str = None,
) -> plt.figure:
    """
    Plot signal and background predictions.
    """

    plt.close()
    assert type(results["train"]) is pd.DataFrame, f"Training set is not a DataFrame: ({type(results['train'])})"
    assert type(results["test"]) is pd.DataFrame, f"Testing set is not a DataFrame: ({type(results['test'])})"

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # List to store train plot
    store_train = []

    for idx, flav in enumerate(results.keys()):
        # Set save paths for plots
        paths = [
            f"{cwd}/../../scripts/plots/BDT/{year}/{xgb_name}/events_BDT_{flav}.pdf",
            f"{cwd}/../../scripts/plots/BDT/{year}/{xgb_name}/sig_bkg_compare_BDT_{flav}.pdf",
        ]

        sig_scale = (len(results[flav][results[flav]["category"] == 0]) / len(results[flav][results[flav]["category"] == 1])) * np.ones_like(results[flav][results[flav]["category"] == 1].BDT_score)
        assert len(sig_scale) == len(results[flav][results[flav]["category"] == 1].BDT_score), f"{len(sig_scale)}\t {len(results[flav][results[flav]['category'] == 1].BDT_score)}"

        # Plot all predictions
        if not suppressSave:
            # sig_scale is already handled in plot_all
            plot_all(results, flav, paths[0])

        # Plot signal and background separately
        fig = plt.figure(0+idx if not suppressSave else 9+idx)
        bins = 70
        plt.hist(
            results[flav][results[flav]["category"] == 1].BDT_score,
            bins = np.linspace(0, 1, bins),
            histtype = "step",
            color = "midnightblue" if idx == 0 else "green",
            label = f"Signal {flav.title()}",
            weights = sig_scale
        )
        plt.hist(
            results[flav][results[flav]["category"] == 0].BDT_score,
            bins = np.linspace(0, 1, bins),
            histtype = "step",
            color = "firebrick" if idx == 0 else "black",
            label = f"Background {flav.title()}"
        )

        plt.xlabel("Prediction from BDT", fontsize=12)
        plt.ylabel("Events", fontsize=12)
        plt.title(year)
        plt.yscale("log")
        plt.legend(frameon=False)

        ax = plt.gca()
        ax.set_ylim(10**-2, ax.get_ylim()[1])

        # Save figures
        if not suppressSave:
            # Save log scale plot
            fig.savefig(paths[1])
            assert os.path.isfile(paths[1]), "Saved plot does not exist."
            print(f"Saved {paths[1].split('/')[-1]} to file.")
 
            # Close figures
            plt.close(fig)

        store_train.append(fig)
        
        # Plot signal and background LINEAR
        fig2 = plt.figure(2+idx if not suppressSave else 91+idx)
        bins = 40
        plt.hist(
            results[flav][results[flav]["category"] == 1].BDT_score,
            bins = np.linspace(0, 1, bins),
            histtype = "step",
            color = "midnightblue",
            label = f"Signal {flav.title()}",
            weights = sig_scale
        )
        plt.hist(
            results[flav][results[flav]["category"] == 0].BDT_score,
            bins = np.linspace(0, 1, bins),
            histtype = "step",
            color = "firebrick",
            label = f"Background {flav.title()}"
        )

        plt.xlabel("Prediction from BDT", fontsize=12)
        plt.ylabel("Events", fontsize=12)
        plt.yscale("linear")
        plt.legend(frameon=False)

        ax = plt.gca()
        ax.set_ylim(0, ax.get_ylim()[1])

        # Save figures
        if not suppressSave:
            # Save linear scale plot
            plt.yscale("linear")
            fig2.savefig(paths[1].replace(".pdf", "_linear.pdf"))
            assert os.path.isfile(paths[1].replace(".pdf", "_linear.pdf")), "Saved plot does not exist."
            print(f"Saved {paths[1].replace('.pdf', '_linear.pdf').split('/')[-1]} to file.")
            
            # Close figures
            plt.close(fig2)

    sig_scale1 = (len(results[list(results.keys())[0]][results[list(results.keys())[0]]["category"] == 0]) / len(results[list(results.keys())[0]][results[list(results.keys())[0]]["category"] == 1])) * np.ones_like(results[list(results.keys())[0]][results[list(results.keys())[0]]["category"] == 1].BDT_score)
    sig_scale2 = (len(results[list(results.keys())[1]][results[list(results.keys())[1]]["category"] == 0]) / len(results[list(results.keys())[1]][results[list(results.keys())[1]]["category"] == 1])) * np.ones_like(results[list(results.keys())[1]][results[list(results.keys())[1]]["category"] == 1].BDT_score)

    fig3 = plt.figure(3+idx if not suppressSave else 991+idx)
    bins = 40
    plt.hist(
        results[list(results.keys())[0]][results[list(results.keys())[0]]["category"] == 0].BDT_score,
        bins = np.linspace(0, 1, bins),
        histtype = "step",
        color = "orange",
        label = f"Background {list(results.keys())[0]}"
    )
    plt.hist(
        results[list(results.keys())[1]][results[list(results.keys())[1]]["category"] == 0].BDT_score,
        bins = np.linspace(0, 1, bins),
        histtype = "step",
        color = "firebrick",
        label = f"Background {list(results.keys())[1]}"
    )
    plt.hist(
        results[list(results.keys())[0]][results[list(results.keys())[0]]["category"] == 1].BDT_score,
        bins = np.linspace(0, 1, bins),
        histtype = "step",
        color = "violet",
        label = f"Signal {list(results.keys())[0]}",
        weights = sig_scale1
    )
    plt.hist(
        results[list(results.keys())[1]][results[list(results.keys())[1]]["category"] == 1].BDT_score,
        bins = np.linspace(0, 1, bins),
        histtype = "step",
        color = "midnightblue",
        label = f"Signal {list(results.keys())[1]}",
        weights = sig_scale2
    )

    plt.xlabel("Prediction from BDT", fontsize=12)
    plt.ylabel("Events", fontsize=12)
    plt.yscale("linear")
    plt.legend(frameon=False)

    ax = plt.gca()
    ax.set_ylim(0, ax.get_ylim()[1])

    # Save figures
    if not suppressSave:
        # Save linear scale plot
        path_superimposed = f"{cwd}/../../scripts/plots/BDT/{year}/{xgb_name}/sig_bkg_compare_BDT_superimposed.pdf"
        plt.yscale("linear")
        fig3.savefig(paths[1].replace(".pdf", "_linear.pdf"))
        assert os.path.isfile(paths[1].replace(".pdf", "_linear.pdf")), "Saved plot does not exist."
        print(f"Saved {paths[1].replace('.pdf', '_linear.pdf').split('/')[-1]} to file.")
        
        # Close figures
        plt.close(fig3)
        
    return store_train[0]


def plot_test_train_superimposed(
    results: Dict[str, pd.DataFrame],
    year: str,
    xgb_name: str = None,
    suppressSave: bool = False
) -> plt.figure:
    """
    Plot signal and background predictions for training and testing superimposed.
    """

    assert type(results["train"]) is pd.DataFrame, f"Training set is not a DataFrame: ({type(results['train'])})"
    assert type(results["test"]) is pd.DataFrame, f"Testing set is not a DataFrame: ({type(results['test'])})"

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Set save paths for plots
    path = f"{cwd}/../../scripts/plots/BDT/{year}/{xgb_name}/train_test_superimposed_BDT.pdf"

    sig_scale_train = (len(results["train"][results["train"]["category"] == 0]) / len(results["train"][results["train"]["category"] == 1])) * np.ones_like(results["train"][results["train"]["category"] == 1].BDT_score)
    #sig_scale_test = (len(results["test"][results["test"]["category"] == 0]) / len(results["test"][results["test"]["category"] == 1])) * np.ones_like(results["test"][results["test"]["category"] == 1].BDT_score)
    
    # Normalize signal test to bkg train and bkg test to bkg train so everything is normalized to bkg train.
    sig_scale_test = (len(results["train"][results["train"]["category"] == 0]) / len(results["test"][results["test"]["category"] == 1])) * np.ones_like(results["test"][results["test"]["category"] == 1].BDT_score)
    bkg_scale_test = (len(results["train"][results["train"]["category"] == 0]) / len(results["test"][results["test"]["category"] == 0])) * np.ones_like(results["test"][results["test"]["category"] == 0].BDT_score)

    for sig_scale, flav in zip([sig_scale_train, sig_scale_test], ["train", "test"]):
        assert len(sig_scale) == len(results[flav][results[flav]["category"] == 1].BDT_score), f"{len(sig_scale)}\t {len(results[flav][results[flav]['category'] == 1].BDT_score)}"

    # Plot signal and background separately
    import random
    fig = plt.figure(random.randint(100,100000))
    bins = 70
    x_axis = np.linspace(0, 1, bins)

    plt.hist(
        results["train"][results["train"]["category"] == 1].BDT_score,
        bins = x_axis,
        histtype = "step",
        edgecolor = "violet",
        color = "violet",
        label = f"Signal_train",
        weights = sig_scale_train
    )
    plt.hist(
        results["train"][results["train"]["category"] == 0].BDT_score,
        bins = x_axis,
        histtype = "step",
        edgecolor = "orange",
        color = "orange",
        label = f"Background_train",
    )
    plt.hist(
        results["test"][results["test"]["category"] == 1].BDT_score,
        bins = x_axis,
        histtype = "step",
        edgecolor = "midnightblue",
        color = "midnightblue",
        label = f"Signal_test",
        weights = sig_scale_test
    )
    plt.hist(
        results["test"][results["test"]["category"] == 0].BDT_score,
        bins = x_axis,
        histtype = "step",
        edgecolor = "firebrick",
        color = "firebrick",
        label = f"Background_test",
        weights = bkg_scale_test
    )

    plt.xlabel("Prediction from BDT", fontsize=12)
    plt.ylabel("Events", fontsize=12)
    plt.yscale("log")
    plt.legend(frameon=False)

    from matplotlib.ticker import AutoMinorLocator
    ax = plt.gca()
    ax.set_ylim(10**-2, ax.get_ylim()[1])
    ax.xaxis.set_minor_locator(AutoMinorLocator()) 
    ax.tick_params(which='minor', length=4, color='r')

    # Save figure
    if not suppressSave:
        fig.savefig(path)
        assert os.path.isfile(path), "Saved plot does not exist."
        print(f"Saved {path.split('/')[-1]} to file.")
        plt.close(fig)
    
    return fig


def plot_all(
    results: Dict[str, pd.DataFrame],
    flav: str,
    path: str
) -> None:
    """
    Plot combined plot to see BDT results distribution.
    """

    sig_scale = (len(results[flav][results[flav]["category"] == 0]) / len(results[flav][results[flav]["category"] == 1])) * np.ones_like(results[flav].BDT_score)
    assert len(sig_scale) == len(results[flav].BDT_score), f"{len(sig_scale)}\t {len(results[flav].BDT_score)}"

    plt.figure()
    bins = 70
    plt.hist(
        results[flav].BDT_score,
        bins = np.linspace(0, 1, bins),
        histtype = "step",
        color = "darkgreen",
        label = "All events",
        weights = sig_scale
    )
    plt.xlabel("Prediction from BDT", fontsize=12)
    plt.ylabel("Events", fontsize=12)
    plt.yscale("log")
    plt.legend(frameon=False)

    ax = plt.gca()
    ax.set_ylim(10**-2, ax.get_ylim()[1])

    # Save figure
    plt.savefig(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path.split('/')[-1]} to file.")
    plt.close()


def plot_roc(
    results: Dict[str, pd.DataFrame],
    path: str,
    sig: Optional[float] = None,
    bkg: Optional[float] = None,
) -> None:
    """
    Plot efficiency 
    """

    plt.figure()

    stats = {}
    auc = []
    colors = ["red", "blue"]
    for c, (flav, df) in zip(colors, results.items()):
        # FPR = false positive rate 
        # TPR = true positive rate
        # Thresholds = cuts on BDT score
        # AUC = area under ROC curve
        # ROC = receiver operating characteristic curve
        fpr, tpr, thresholds = roc_curve(df["category"], df["BDT_score"])
        auc_tmp = roc_auc_score(df["category"], df["BDT_score"])
        auc.append(auc_tmp)
        stats.update({flav: {"fpr": fpr.tolist(), "tpr": tpr.tolist(), "thresholds": thresholds.tolist(), "auc": auc_tmp, "events": {"signal": len(df[df.category == 1]), "bkg": len(df[df.category == 0])}}})

        plt.plot(
            1-fpr,
            tpr,
            color = c,
            label = flav
        )

    plt.text(0.0, 0.2, f"{list(results.keys())[0].title()} AUC: {auc[0]:.4}")
    plt.text(0.0, 0.15, f"{list(results.keys())[1].title()} AUC: {auc[1]:.4}")

    plt.xlabel("Signal Efficiency", fontsize=12)
    plt.ylabel("Background Rejection", fontsize=12)
    plt.legend(frameon=False)

    # Save figure
    plt.savefig(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path.split('/')[-1]} to file.")

    # Make a second copy of the plot zoomed to 0.85 and 1.0
    ax = plt.gca()
    ax.set_xlim(0.85, 1.0)
    ax.set_ylim(0.85, 1.0)

    # Save zoomed figure
    path_zoom = path.replace(".pdf", "_zoomed.pdf")
    plt.savefig(path_zoom)
    assert os.path.isfile(path_zoom), "Saved plot does not exist."
    print(f"Saved {path_zoom.split('/')[-1]} to file.")
    plt.close()

    # Save JSON of scatter plot weights in case I need to remake this!
    with open(path.replace(".pdf", ".json"), "w") as f:
        json.dump(
            stats,
            f,
            indent = 4
        )
    assert os.path.isfile(path.replace(".pdf", ".json")), "Saved JSON does not exist."
    print(f"Saved {path.replace('.pdf', '.json')} to file.")


def plot_corr_matrix(
    results: Dict[str, pd.DataFrame],
    path: str,
    full: bool = False,
    sig: pd.DataFrame = None,
    bkg: pd.DataFrame = None,
) -> None:
    """
    Plots correlation matrix for input DataFrame.
    """

    # Combine training/testing signal/background (respectively) to look at full samples instead of half.
    def plotHeatMap(pruned_df, key, cat):
        heatmap = sns.heatmap(pruned_df.corr(), vmin=-1, vmax=1, annot=True)
        heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=60)
        heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=60)

        # Save figure
        plt.savefig(path.replace(".pdf", f"_{key}_{cat}.pdf"))
        assert os.path.isfile(path.replace(".pdf", f"_{key}_{cat}.pdf")), "Saved plot does not exist."
        print(f"Saved {path.replace('.pdf', '_' + key + '_' + cat + '.pdf').split('/')[-1]} to file.")
        plt.close()


    # Plot training/testing sets. Separates signal and background sample plots.
    if not full:
        for i, (key, df) in enumerate(results.items()):
            # Do signal/bkg separately
            for j, cat in enumerate(["sig", "bkg"]):
                plt.figure(i+j, figsize=(14,14))
                cp_df = df.copy()
                cp_df = cp_df[cp_df.category == (0 if cat == "bkg" else 1)]
                assert len(df) > 0

                pruned_df = cp_df
                prune = []
                for remove in ["BDT_score", "category", "m_hyp", "mass_gggg", "weight"]:
                    if remove in cp_df.columns:
                        prune.append(remove)
                pruned_df = cp_df.drop(columns=prune)

                plotHeatMap(pruned_df, key, cat)

    if not full:
        combined_sig = pd.concat([results["train"][results["train"].category == 1], results["test"][results["test"].category == 1]], axis=0)
        combined_bkg = pd.concat([results["train"][results["train"].category == 0], results["test"][results["test"].category == 0]], axis=0)
    else:
        combined_sig = sig
        combined_bkg = bkg
    
    for i, (key, df) in enumerate(zip(["sig", "bkg"], [combined_sig, combined_bkg])):
        prune = []
        for remove in ["BDT_score", "category", "m_hyp", "mass_gggg", "weight"]:
            if remove in df.columns:
                prune.append(remove)

        # Remove extra columns and plot
        df = df.drop(columns=prune) 
        plt.figure(5+i, figsize=(14,14))
        plotHeatMap(df, "fullSample", key)
        

def plot_data_sig_bkg_comparison(
    hists: Dict[str, TH1D],
    title: str,
    path: str,
    genWeight: float,
    log: Optional[bool] = True,
    maximum: Optional[float] = None,
    mass: Optional[str] = None,
    canvas_num: Optional[int] = None,
    reweight: Optional[bool] = False,
    sub: Optional[str] = None,
    scaleLumi: bool = False,
    year: str = "",
    Scales: str = None,
    extra_ratio: bool = False,
    omit_data: bool = False,
    ignore_scale: bool = False,
) -> Dict[str, TH1D]:
    """
    Plots distributions for a given pseudoscalar mass, event mixing (background), and data.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Canvas setup
    canvas = TCanvas(f"canvas_{canvas_num}" if canvas_num is not None else "canvas", "canvas", 1000, 1000)
    
    # Pad setup
    if not omit_data:
        pad1 = TPad(f"pad1_{canvas_num}" if canvas_num is not None else "pad1" , "pad1", 0, 0.25, 1, 1)
        pad2 = TPad(f"pad2_{canvas_num}" if canvas_num is not None else "pad2" , "pad2", 0, 0, 1, 0.23)
        
        pad1.SetBottomMargin(0.0001)
        pad1.SetBorderMode(0)
        pad2.SetTopMargin(0.01)
        pad2.SetBottomMargin(0.3)
        pad2.SetBorderMode(0)
        if log:
            pad1.SetLogy()

        pad1.Draw()
        pad2.Draw()
    else:
        if log:
            canvas.SetLogy()

    # Legend setup
    if not reweight:
        lsize = [0.65, 0.65, 0.85, 0.80]
    else:
        lsize = [0.5, 0.55, 0.87, 0.80]
    legend = TLegend(lsize[0], lsize[1], lsize[2], lsize[3])
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # Handle scaling/max
    scale = {"signal": 1.0, "background": 1.0, "data": 1.0, "background_reweight": 1.0}
    for key, hist in hists.items():
        hist.Sumw2()
        # Get maximum of histograms
        if maximum is None:
            maximum = 0.0
            if maximum < hist.GetBinContent(hist.GetMaximumBin()):
                maximum = hist.GetBinContent(hist.GetMaximumBin())
        
        if not ignore_scale:
            # Signal
            if "signal" in key:
                if mass is None:
                    if "15" in key:
                        #scale["signal"] = Scales.lumi * 1.0 / Scales.sumw[year]["15_GeV"]  # Signal scaled to xs = 1 fb (https://indico.cern.ch/event/1057512/contributions/4449049/attachments/2281868/3879051/H4G_Approval_15Jul2021.pdf)
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)  # Signal scaled to xs = 1 fb
                    elif "20" in key:
                        #scale["signal"] = Scales.lumi * 1.0 / Scales.sumw[year]["20_GeV"]
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "25" in key:
                        #scale["signal"] = Scales.lumi * 1.0 / Scales.sumw[year]["25_GeV"]
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "30" in key:
                        #scale["signal"] = Scales.lumi * 1.0 / Scales.sumw[year]["30_GeV"]
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "35" in key:
                        #scale["signal"] = Scales.lumi * 1.0 / Scales.sumw[year]["35_GeV"]
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "40" in key:
                        #scale["signal"] = Scales.lumi * 1.0 / Scales.sumw[year]["40_GeV"]
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "45" in key:
                        #scale["signal"] = Scales.lumi * 1.0 / Scales.sumw[year]["45_GeV"]
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "50" in key:
                        #scale["signal"] = Scales.lumi * 1.0 / Scales.sumw[year]["50_GeV"]
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "55" in key:
                        #scale["signal"] = Scales.lumi * 1.0 / Scales.sumw[year]["55_GeV"]
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "60" in key:
                        #scale["signal"] = Scales.lumi * 1.0 / Scales.sumw[year]["60_GeV"]
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "signal" in key:
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                elif mass is not None:
                    if "15" in mass:
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)  # Signal scaled to xs = 1 fb
                    elif "20" in mass:
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "25" in mass:
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "30" in mass:
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "35" in mass:
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "40" in mass:
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "45" in mass:
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "50" in mass:
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "55" in mass:
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
                    elif "60" in mass:
                        scale["signal"] = Scales.lumi_fb * 1.0 / (genWeight)
            # Background / Data
            if "background" in key:
                scale[key] = Scales.lumi if scaleLumi else hists["data"].Integral() / hists["background"].Integral()
            elif "data" in key:
                scale["data"] = Scales.lumi if scaleLumi else 1.0


    # Main plotting on pad1
    if not omit_data:
        pad1.cd()
    colors = [kAzure-4, kBlack, kRed]
    if reweight:
        colors.insert(1, kMagenta)
    if omit_data:
        colors = [colors[0]] + list(colors[2:])
        assert len(hists) == len(colors), (len(hists), len(colors))

    import math
    for (key, hist), color in zip(hists.items(), colors):
        hist.SetTitle("")
        hist.GetYaxis().SetTitleFont(42)
        hist.GetYaxis().SetTitle("Events")
        hist.SetLineColor(color)
        hist.Scale(scale[key])
        hist.Sumw2()
        for bi in range(hist.GetNbinsX()):
            b = bi+1
            #assert hist.GetBinError(b) - math.sqrt(hist.GetBinContent(b)) < 0.001, (b, hist.GetBinError(b), math.sqrt(hist.GetBinContent(b)))

        if "signal" in key:
            hist.SetLineWidth(3)
            if mass == None:
                legend.AddEntry(hist, "Signal MC Merged", "l")
            else:
                legend.AddEntry(hist, "m_{a} = "+mass, "l")
        
        if "background" in key and "error" not in key:
            if not reweight:
                hist.SetFillColor(color)
                legend.AddEntry(hist, "Event Mixing", "f")
            elif reweight and "reweight" not in key:
                hist.SetLineWidth(3)
                legend.AddEntry(hist, "Event Mixing (w/out weights)", "l")
            elif reweight and "reweight" in key:
                hist.SetFillColor(color)
                legend.AddEntry(hist, "Event Mixing (w/ weights)", "f")
            
        elif "data" in key:
            hist.SetMarkerSize(1)
            hist.SetMarkerStyle(20)
            hist.SetMarkerColor(kBlack)
            hist.SetFillStyle(0)
            legend.AddEntry(hist, "Data", "lep")

        # Draw histograms
        if "background_reweight" in key:
            background_error = hist.Clone("background_error")
            background_error.SetLineColor(kGray+3)
            background_error.SetFillColor(kGray+3)
            background_error.SetFillStyle(3008)
            legend.AddEntry(background_error, "Statistical Uncertainty", "f")
            hist.Draw("hist same")
            background_error.Draw("E2 same")
        elif "data" in key:
            hist.Draw("P E1 same")
        else:
            hist.Draw("hist same")
        if not omit_data:
            hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum if log else maximum*1.5)
        else:
            hist.GetYaxis().SetRangeUser(0.1, maximum*10 if log else maximum*1.5)
        #hist.GetYaxis().SetRangeUser(0.9, maximum if log else maximum*1.5)

    # Draw legend
    legend.Draw()
    
    if not omit_data:
        # Residuals histogram setup
        residuals = residuals_error_prop(hists["background" if not reweight else "background_reweight"], hists["data"], title)
        if extra_ratio:
            residuals_noRW = residuals_error_prop(hists["background"], hists["data"], title)

        # Plotting residuals on pad2
        pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
        pad2.cd()
        gStyle.SetTitleFont(42)
        gStyle.SetLabelFont(42)
        gPad.SetFrameBorderMode(0)
        gPad.SetTickx(1)
        gPad.SetTicky(1)

        mn = 999
        mx = -999
        for rs, color in zip([residuals] if not extra_ratio else [residuals, residuals_noRW], [kBlack] if not extra_ratio else [kBlack, kMagenta]):
            if "Delta R" in title:
                rs.GetXaxis().SetMaxDigits(1)
            rs.SetTitle("")
            rs.GetXaxis().SetTitle(title)
            rs.GetYaxis().SetTitle("Bkg/Data")
            rs.GetXaxis().SetTitleSize(9.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
            rs.GetYaxis().SetTitleSize(8.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
            rs.GetXaxis().SetTitleOffset(1.1)
            rs.GetYaxis().SetTitleOffset(0.4)
            rs.GetYaxis().SetNdivisions(6)
            rs.GetYaxis().CenterTitle(True)
            if mn > rs.GetBinContent(rs.GetMinimumBin()):
                mn = rs.GetBinContent(rs.GetMinimumBin())
            if mx < rs.GetBinContent(rs.GetMaximumBin()):
                mx = rs.GetBinContent(rs.GetMaximumBin())
            rs.GetXaxis().SetLabelSize(2.5*rs.GetXaxis().GetLabelSize())
            rs.GetYaxis().SetLabelSize(2.5*rs.GetYaxis().GetLabelSize())
            rs.SetMarkerStyle(20)
            rs.SetMarkerSize(1)
            rs.SetLineWidth(1)
            rs.SetMarkerColor(color)
            rs.SetLineColor(color)

        mn_scale = 0.5
        mx_scale = 1.5
        residuals.GetYaxis().SetRangeUser(mn * (mn_scale if residuals.GetBinContent(residuals.GetMinimumBin()) != 0 else mn), mx * mx_scale)
        if extra_ratio:
            residuals_noRW.GetYaxis().SetRangeUser(mn * (mn_scale if residuals_noRW.GetBinContent(residuals_noRW.GetMinimumBin()) != 0 else mn), mx * mx_scale)
        residuals.Draw("P same")
        if extra_ratio:
            residuals_noRW.Draw("P same")

        # Line at y=1
        max_edge = residuals.GetBinCenter(residuals.GetNbinsX())
        min_edge = residuals.GetBinCenter(1)
        line = TLine(min_edge, 1.0, max_edge, 1.0)
        line.SetLineStyle(7)
        line.SetLineColor(kBlack)
        line.Draw("same")

        pad1.Modified()
        pad2.Modified()
        canvas.Update()
        plotFancy(pad1, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, pad=True, sub=sub, Scales=Scales, shiftLeft=0.02)
    else:
        plotFancy(canvas, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, sub=sub, Scales=Scales)

    # Save figure
    canvas.SaveAs(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path.split('/')[-1]} to file.")

    return hists



def plot_data_sig_bkg_comparison_combined(
    hists: Dict[str, TH1D],
    title: str,
    path: str,
    genWeight: Dict[str, float],
    maximum: Optional[float] = None,
    canvas_num: Optional[int] = None,
    sub: Optional[str] = None,
    scaleLumi: bool = False,
    year: str = "",
    Scales: str = None
) -> Dict[str, TH1D]:
    """
    Plots distributions for a given pseudoscalar mass, event mixing (background), and data.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Canvas setup
    canvas = TCanvas(f"canvas_{canvas_num}" if canvas_num is not None else "canvas", "canvas", 1000, 1000)
    
    # Pad setup
    pad1 = TPad(f"pad1_{canvas_num}" if canvas_num is not None else "pad1", "pad1", 0, 0.25, 1, 1)
    pad2 = TPad(f"pad2_{canvas_num}" if canvas_num is not None else "pad2", "pad2", 0, 0, 1, 0.23)
    
    pad1.SetBottomMargin(0.0001)
    pad1.SetBorderMode(0)
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(0.3)
    pad2.SetBorderMode(0)
    pad1.SetLogy()

    pad1.Draw()
    pad2.Draw()

    # Legend setup
    legend_left = TLegend(0.52,0.71,0.7,0.79)
    legend_right = TLegend(0.7,0.55,0.9,0.80)
    for legend in [legend_left, legend_right]:
        legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)

    # Handle scaling/max
    scale = {"background": 0.0, "data": 0.0}
    for key, hist in hists.items():
        hist.Sumw2()
        # Get maximum of histograms
        if maximum is None:
            maximum = 0.0
            if maximum < hist.GetBinContent(hist.GetMaximumBin()):
                maximum = hist.GetBinContent(hist.GetMaximumBin())
        
        # Signal
        if "Signal" in key:
            #scale.update({key: 1.0 / Scales.sumw[year][key.replace("Signal_","")]})
            scale.update({key: 1.0 / (genWeight)})

        # Background / Data
        if "background" in key or "data" in key:
            scale.update({key: Scales.lumi if scaleLumi else 1.0})


    # Residuals histogram setup
    residuals = residuals_error_prop(hists["background"], hists["data"], title)
    """
    residuals = hists["background"] / hists["data"]
    for b in range(1, int(residuals.GetNbinsX())+1):
        # residuals.SetBinError(b, hists["data"].GetBinError(b))
        residuals.SetBinError(b, 0)
    """

    # Main plotting on pad1
    pad1.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    colors = [kAzure-4, kBlack, kGreen+3, kMagenta, kYellow, kSpring-3, kOrange, kCyan+1, kRed, kViolet-5, kOrange+4, kBlue]
    assert len(colors) == len(hists)
    for idx, ((key, hist), color) in enumerate(zip(hists.items(), colors)):
        hist.SetTitle("")
        hist.GetYaxis().SetTitleFont(42)
        hist.GetYaxis().SetTitle("Events")
        hist.SetLineColor(color)
        hist.Scale(scale[key])
        if "GeV" in key:
            hist.SetLineWidth(3)
            legend_right.AddEntry(hist, "m_{a} = "+key.replace("Signal_","").replace("_", " "), "l")
        if key == "background":
            hist.SetFillColor(color)
            legend_left.AddEntry(hist, "Event Mixing", "f")

            # Clone background to include statistical errors
            background_error = hists["background"].Clone("background_error")
            background_error.SetLineColor(kGray+3)
            background_error.SetFillColor(kGray+3)
            background_error.SetFillStyle(3008)
            legend_left.AddEntry(background_error, "Statistical Uncertainty", "f")

        elif key == "data":
            hist.SetMarkerSize(1)
            hist.SetMarkerStyle(20)
            hist.SetMarkerColor(kBlack)
            hist.SetFillStyle(0)
            legend_left.AddEntry(hist, "Data", "lep")
            
        hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum*10)
        #hist.GetYaxis().SetRangeUser(0.9, maximum*10)

    # Force draw order so Signal_15_GeV is on top
    hists["background"].Draw("hist same")
    background_error.Draw("E2 same")
    hists["data"].Draw("P E1 same")
    for mass in range(60, 10, -5):
        hists[f"Signal_{mass}_GeV"].Draw("hist same")
    legend_left.Draw()
    legend_right.Draw()
    
    # Plotting residuals on pad2
    pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
    pad2.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    gPad.SetFrameBorderMode(0)
    gPad.SetTickx(1)
    gPad.SetTicky(1)
    residuals.SetTitle("")
    residuals.GetXaxis().SetTitle(title)
    residuals.GetYaxis().SetTitle("Bkg/Data")
    residuals.GetXaxis().SetTitleSize(9.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
    residuals.GetYaxis().SetTitleSize(8.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
    residuals.GetXaxis().SetTitleOffset(1.1)
    residuals.GetYaxis().SetTitleOffset(0.4)
    residuals.GetYaxis().SetNdivisions(6)
    residuals.GetYaxis().CenterTitle(True)
    residuals.GetXaxis().SetLabelSize(2.5*residuals.GetXaxis().GetLabelSize())
    residuals.GetYaxis().SetLabelSize(2.5*residuals.GetYaxis().GetLabelSize())
    residuals.SetMarkerSize(1)
    residuals.SetMarkerStyle(20)
    residuals.SetMarkerColor(kBlack)
    residuals.SetLineColor(kBlack)
    residuals.SetFillColor(kGray+2)
    residuals.Draw("P E2")
    #residuals.Draw("P")

    # Line at y=1
    max_edge = residuals.GetBinCenter(residuals.GetNbinsX())
    min_edge = residuals.GetBinCenter(1)
    line = TLine(min_edge, 1.0, max_edge, 1.0)
    line.SetLineStyle(7)
    line.SetLineColor(kBlack)
    line.Draw("same")

    pad1.Modified()
    pad2.Modified()
    canvas.Update()
    plotFancy(pad1, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, pad=True, sub=sub, Scales=Scales)

    # Save figure
    canvas.SaveAs(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path.split('/')[-1]} to file.")

    return hists



def compare_eras_data_bkg_only(
    hists: Dict[str, TH1D],
    title: str,
    path: str,
    log: Optional[bool] = True,
    maximum: Optional[float] = None,
    canvas_num: Optional[int] = None,
    sub: str = None,
    year: str = "",
    Scales: str = None
) -> Dict[str, TH1D]:
    """
    Plots distributions for event mixing (background) and data comparing different eras. No residuals.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Canvas setup
    canvas = TCanvas(f"canvas_{canvas_num}" if canvas_num is not None else "canvas", "canvas", 1000, 1000)

    # Pad setup
    pad1 = TPad(f"pad1_{canvas_num}" if canvas_num is not None else "pad1", "pad1", 0, 0.25, 1, 1)
    pad2 = TPad(f"pad2_{canvas_num}" if canvas_num is not None else "pad2", "pad2", 0, 0, 1, 0.23)
    
    pad1.SetBottomMargin(0.0001)
    pad1.SetBorderMode(0)
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(0.3)
    pad2.SetBorderMode(0)
    if log:
        pad1.SetLogy()

    pad1.Draw()
    pad2.Draw()

    # Legend setup
    legend = TLegend(0.65,0.65,0.87,0.85)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # Handle scaling/max
    for key, hist in hists.items():
        hist.Sumw2()
        # Get maximum of histograms
        if maximum is None:
            maximum = 0.0
            if maximum < hist.GetBinContent(hist.GetMaximumBin()):
                maximum = hist.GetBinContent(hist.GetMaximumBin())

    residuals_preEE = residuals_error_prop(hists["bkg_preEE"], hists["data_preEE"], title)
    residuals_postEE = residuals_error_prop(hists["bkg_postEE"], hists["data_postEE"], title)
    """
    residuals_preEE = hists["bkg_preEE"] / hists["data_preEE"]
    residuals_postEE = hists["bkg_postEE"] / hists["data_postEE"]
    assert int(residuals_preEE.GetNbinsX()) == int(residuals_postEE.GetNbinsX())
    for b in range(1, int(residuals_preEE.GetNbinsX())+1):
        residuals_preEE.SetBinError(b, 0.0)
        residuals_postEE.SetBinError(b, 0.0)
        # errors are more complicated... ask about how tf to do this. it's clearly not sqrt(N). Erorr prop for bkg/data maybe is correct?
    """

    # Main plotting on pad1
    pad1.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    colors = [kAzure-4, kPink, kGreen+3, kBlack]
    for (key, hist), color in zip(hists.items(), colors):
        hist.SetTitle("")
        hist.GetYaxis().SetTitleFont(42)
        hist.GetYaxis().SetTitle("Events")
        hist.SetLineColor(color)
        if "background" in key or "bkg" in key:
            hist.SetFillColor(color)
            stat_unc = "Stat. Uncertainty"
            if "preEE" in key:
                legend.AddEntry(hist, "Event Mixing preEE", "f")
                # Clone background to include statistical errors
                background_error1 = hist.Clone("background_error1")
                background_error1.SetLineColor(kGray+3)
                background_error1.SetFillColor(kGray+3)
                background_error1.SetFillStyle(3004)
                legend.AddEntry(background_error1, stat_unc+" preEE", "f")
            elif "postEE" in key:
                legend.AddEntry(hist, "Event Mixing postEE", "f")
                # Clone background to include statistical errors
                background_error2 = hist.Clone("background_error2")
                background_error2.SetLineColor(kGray+3)
                background_error2.SetFillColor(kGray+3)
                background_error2.SetFillStyle(3005)
                legend.AddEntry(background_error2, stat_unc+" postEE", "f")

        elif "data" in key:
            hist.SetMarkerSize(1)
            hist.SetMarkerStyle(20)
            hist.SetMarkerColor(color)
            hist.SetFillStyle(0)
            if "preEE" in key:
                legend.AddEntry(hist, "Data preEE", "lep")
            elif "postEE" in key:
                legend.AddEntry(hist, "Data postEE", "lep")
            
        hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum*10 if log else maximum*1.5)
        #hist.GetYaxis().SetRangeUser(0.9, maximum*10 if log else maximum*1.5)

    # Force draw order so bigger histogram is on top
    hists["bkg_postEE"].Draw("hist same")
    hists["bkg_preEE"].Draw("hist same")
    background_error1.Draw("E2 same")
    background_error2.Draw("E2 same")
    hists["data_postEE"].Draw("P E1 same")
    hists["data_preEE"].Draw("P E1 same")
    legend.Draw()
    
    # Plotting residuals on pad2
    pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
    pad2.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    gPad.SetFrameBorderMode(0)
    gPad.SetTickx(1)
    gPad.SetTicky(1)

    for residuals, color in zip([residuals_preEE, residuals_postEE], [kAzure-4, kGreen+3]):
        if "Delta R" in title:
            residuals.GetXaxis().SetMaxDigits(1)
        residuals.SetTitle("")
        residuals.GetXaxis().SetTitle(title)
        residuals.GetYaxis().SetTitle("Bkg/Data")
        residuals.GetXaxis().SetTitleSize(9.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
        residuals.GetYaxis().SetTitleSize(8.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
        residuals.GetXaxis().SetTitleOffset(1.1)
        residuals.GetYaxis().SetTitleOffset(0.4)
        residuals.GetYaxis().SetNdivisions(6)
        residuals.GetYaxis().CenterTitle(True)
        residuals.GetXaxis().SetLabelSize(2.5*residuals.GetXaxis().GetLabelSize())
        residuals.GetYaxis().SetLabelSize(2.5*residuals.GetYaxis().GetLabelSize())
        residuals.SetMarkerSize(1)
        residuals.SetMarkerStyle(20)
        residuals.SetMarkerColor(color)
        residuals.SetLineColor(color)
        #residuals.SetFillColor(kGray)
        residuals.Draw("P E2 same")
        #residuals.Draw("P same")

    # Line at y=1
    max_edge = residuals_preEE.GetBinCenter(residuals.GetNbinsX())
    min_edge = residuals_preEE.GetBinCenter(1)
    line = TLine(min_edge, 1.0, max_edge, 1.0)
    line.SetLineStyle(7)
    line.SetLineColor(kBlack)
    line.Draw("same")

    pad1.Modified()
    pad2.Modified()
    canvas.Update()
    plotFancy(pad1, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, pad=True, sub=sub, Scales=Scales)

    # Save figure
    canvas.SaveAs(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path.split('/')[-1]} to file.")

    return hists


def plotFourMassesCompared(
    hists: Dict[str, TH1D],
    title: str,
    path: str,
    genWeight: List[float],
    log: Optional[bool] = True,
    maximum: Optional[float] = None,
    canvas_num: Optional[int] = None,
    ignore_scales: Optional[bool] = False,
    sub: str = None,
    scaleLumi: bool = False,
    year: str = "",
    Scales: str = None
) -> Dict[str, TH1D]:
    """
    Plots distributions for pseudoscalar masses (15, 30, 40, and 60 GeV) only. Very similar to other compare plotters.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Canvas setup
    canvas = TCanvas(f"canvas_{canvas_num}" if canvas_num is not None else "canvas", "canvas", 1000, 1000)
    
    if log:
        canvas.SetLogy()

    # Legend setup
    legend = TLegend(0.65,0.82-(0.2*4.0/6.0),0.87,0.82)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # Handle scaling/max
    scale = {}
    for key, hist in hists.items():
        hist.Sumw2()
        # Get maximum of histograms
        if maximum is None:
            maximum = 0.0
            if maximum < hist.GetBinContent(hist.GetMaximumBin()):
                maximum = hist.GetBinContent(hist.GetMaximumBin())
        
        # Set scale for signal
        print(genWeight)
        if "15_GeV" in key:
            #scale["Signal_15_GeV"] = Scales.lumi * 1.0 / Scales.sumw[year]["15_GeV"]  # Signal scaled to xs = 1 fb (https://indico.cern.ch/event/1057512/contributions/4449049/attachments/2281868/3879051/H4G_Approval_15Jul2021.pdf)
            scale["Signal_15_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[0])  # Signal scaled to xs = 1 fb
        elif "30_GeV" in key:
            scale["Signal_30_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[1])
        elif "40_GeV" in key:
            scale["Signal_40_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[2])
        elif "60_GeV" in key:
            scale["Signal_60_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[3])

    # Main plotting on pad1
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    colors = [kMagenta, kOrange, kRed, kBlue]
    for (key, hist), color in zip(hists.items(), colors):
        hist.SetTitle("")
        hist.GetYaxis().SetTitleFont(42)
        hist.GetYaxis().SetTitle("Events")
        hist.SetLineColor(color)
        hist.GetXaxis().SetTitleOffset(1.1)
        hist.GetXaxis().SetTitle(title)
        if not ignore_scales:
            hist.Scale(scale[key])
        if "GeV" in key:
            hist.SetLineWidth(3)
            legend.AddEntry(hist, "m_{a} = "+key.replace("Signal_","").replace("_", " "), "l")
            
        hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum*10 if log else maximum*1.5)
            
    # Force draw order so 15_GeV is on top
    new_keys = {}
    for key in hists.keys():
        for m in range(15,63,5):
            if f"{m}_GeV" in key:
                new_keys.update({f"Signal_{m}_GeV": key})
                break

    if "mva" in new_keys["Signal_60_GeV"].lower():
        canvas.SetLogy()
        hists[new_keys["Signal_60_GeV"]].GetYaxis().SetRangeUser(10**-4 - 10**-5, hists[new_keys["Signal_60_GeV"]].GetMaximum() * 10**2)
        hists[new_keys["Signal_40_GeV"]].GetYaxis().SetRangeUser(10**-4 - 10**-5, hists[new_keys["Signal_40_GeV"]].GetMaximum() * 10**2)
        hists[new_keys["Signal_30_GeV"]].GetYaxis().SetRangeUser(10**-4 - 10**-5, hists[new_keys["Signal_30_GeV"]].GetMaximum() * 10**2)
        hists[new_keys["Signal_15_GeV"]].GetYaxis().SetRangeUser(10**-4 - 10**-5, hists[new_keys["Signal_15_GeV"]].GetMaximum() * 10**2)

    hists[new_keys["Signal_60_GeV"]].Draw("hist same")
    hists[new_keys["Signal_40_GeV"]].Draw("hist same")
    hists[new_keys["Signal_30_GeV"]].Draw("hist same")
    hists[new_keys["Signal_15_GeV"]].Draw("hist same")
    legend.Draw()
    
    plotFancy(canvas, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, sub=sub, Scales=Scales, shiftLeft=0.07, simulation=True)

    # Save figure
    canvas.SaveAs(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path.split('/')[-1]} to file.")

    return hists


def compare_four_masses_mva(
    hists: Dict[str, TH1D],
    title: str,
    path: str,
    genWeight: List[float],
    log: Optional[bool] = True,
    maximum: Optional[float] = None,
    canvas_num: Optional[int] = None,
    ignore_scales: Optional[bool] = False,
    sub: str = None,
    scaleLumi: bool = False,
    year: str = "",
    Scales: str = None,
) -> Dict[str, TH1D]:
    """
    Plots distributions for pseudoscalar masses (15, 30, 40, and 60 GeV) and event mixing (background) for MVA ID training variables individually. Includes error bars on signal MC for one mass point per plot.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Canvas setup
    canvas = TCanvas(f"canvas_{canvas_num}" if canvas_num is not None else "canvas", "canvas", 1000, 1000)
    
    canvas.SetFrameLineWidth(0)
    if log:
        canvas.SetLogy()

    # Legend setup
    legend = TLegend(0.65,0.70,0.87,0.80)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # Handle scaling/max
    scale = {"signal": 1.0, "background": 1.0}
    for key, hist in hists.items():
        hist.Sumw2()
        
        # Set scale for signal
        if "15_GeV" in key:
            scale["Signal_15_GeV"] = hists["background"].Integral() / hist.Integral()
        elif "30_GeV" in key:
            scale["Signal_30_GeV"] = hists["background"].Integral() / hist.Integral()
        elif "40_GeV" in key:
            scale["Signal_40_GeV"] = hists["background"].Integral() / hist.Integral()
        elif "60_GeV" in key:
            scale["Signal_60_GeV"] = hists["background"].Integral() / hist.Integral()
        elif "background" in key:
            # Do this since we reweight using scaled weights (i.e., no normalization in the reweighting itself)
            scale["background"] = 1.0

    # Main plotting
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    colors = [kAzure-4, kMagenta, kOrange, kRed, kBlue]
    for (key, hist), color in zip(hists.items(), colors):
        hist.SetTitle("")
        hist.GetYaxis().SetTitleFont(42)
        hist.GetYaxis().SetTitle("Events")
        hist.GetXaxis().SetTitle(title)
        hist.GetXaxis().SetTitleOffset(1.1)
        hist.Sumw2()

        if not ignore_scales:
            hist.Scale(scale[key])
        if "GeV" in key:
            hist.SetLineWidth(2)
        if key == "background":
            hist.SetFillColor(color)

            # Clone background to include statistical errors
            background_error = hists["background"].Clone("background_error")
            background_error.SetLineColor(kGray+3)
            background_error.SetFillColor(kGray+3)
            background_error.SetFillStyle(3008)
 
    new_keys = {}
    for key in hists.keys():
        for m in range(15,63,5):
            if f"{m}_GeV" in key:
                new_keys.update({f"Signal_{m}_GeV": key})
                break

    for mass in ["60_GeV", "40_GeV", "30_GeV", "15_GeV"]:
        legend.Clear()

        if maximum is None:
            maximum_bkg = hists["background"].GetBinContent(hists["background"].GetMaximumBin())
            maximum_sig = hists[new_keys[f"Signal_{mass}"]].GetBinContent(hists[new_keys[f"Signal_{mass}"]].GetMaximumBin())
            maximum = maximum_bkg if maximum_bkg > maximum_sig else maximum_sig
            if "esEffSigmaRR" in path or "esEnergyOverRawE" in path and not log:
                maximum = maximum / 300  # Scales by 1/200 = 1.5/300
            elif "PFCluster" in path or "pfPhoIso03_compare" in path and not log:
                maximum = maximum * 0.5/1.5
            elif "trkSum" in path or "pfChargedIso_" in path and not log:
                maximum = maximum * 0.25/1.5 

        # Set min/max for both background and signal MC histograms
        hists["background"].GetYaxis().SetRangeUser(0.9 if log else -0.1, maximum*10 if log else maximum*1.5)
        hists[new_keys[f"Signal_{mass}"]].GetYaxis().SetRangeUser(0.9 if log else -0.1, maximum*10 if log else maximum*1.5)

        # Force draw order so 15_GeV is on top
        hists["background"].Draw("hist")
        legend.AddEntry(hists["background"], "Event Mixing", "f")
        background_error.Draw("E2 same")
        legend.AddEntry(background_error, "Statistical Uncertainty", "f")

        hists[new_keys[f"Signal_{mass}"]].SetLineColor(kBlack)
        hists[new_keys[f"Signal_{mass}"]].SetMarkerColor(kBlack)
        hists[new_keys[f"Signal_{mass}"]].Draw("P E1 same")
        legend.AddEntry(hists[new_keys[f"Signal_{mass}"]], f"m_{{a}} = {mass.replace('_',' ')}", "lpe")
        legend.Draw()
    
        canvas.Update()
        plotFancy(canvas, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, sub=sub, Scales=Scales, shiftLeft=0.05)

        # Save figure
        new_path = path.replace("/mass/", f"/{mass}/")
        if not os.path.exists(os.path.join(*new_path.split("/")[:-1])):
            os.mkdir(os.path.join(*new_path.split("/")[:-1]))
        canvas.SaveAs(new_path)
        assert os.path.isfile(new_path), "Saved plot does not exist."
        print(f"Saved {new_path.split('/')[-1]} to file.")

    return hists


def compare_sigMixing_bkg(
    hists: Dict[str, TH1D],
    data: TH1D,
    title: str,
    path: str,
    genWeight: List[float],
    log: Optional[bool] = True,
    maximum: Optional[float] = None,
    canvas_num: Optional[int] = None,
    ignore_scales: Optional[bool] = False,
    sub: str = None,
    scaleLumi: bool = False,
    year: str = "",
    Scales: str = None,
    skipSigScale: bool = False
) -> Dict[str, TH1D]:
    """
    Plots distributions for pseudoscalar masses (15, 30, 40, and 60 GeV), event mixing (background), and data.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Canvas setup
    canvas = TCanvas(f"canvas_{canvas_num}" if canvas_num is not None else "canvas", "canvas", 1000, 1000)
    
    # Pad setup
    pad1 = TPad(f"pad1_{canvas_num}" if canvas_num is not None else "pad1", "pad1", 0, 0.25, 1, 1)
    pad2 = TPad(f"pad2_{canvas_num}" if canvas_num is not None else "pad2", "pad2", 0, 0, 1, 0.23)
    
    pad1.SetBottomMargin(0.0001)
    pad1.SetBorderMode(0)
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(0.3)
    pad2.SetBorderMode(0)
    if log:
        pad1.SetLogy()

    pad1.Draw()
    pad2.Draw()

    # Legend setup
    legend = TLegend(0.62,0.58,0.87,0.83)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    del hists["data"]

    # Handle scaling/max
    scale = {"signal": 1.0, "background": 1.0}
    for key, hist in hists.items():
        hist.Sumw2()
        
        # Set scale for signal
        if "15_GeV" in key:
            if not skipSigScale:
                scale["Signal_15_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[0])  # Signal scaled to xs = 1 fb
            else:
                scale[key] = 1.0
        elif "30_GeV" in key:
            if not skipSigScale:
                scale["Signal_30_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[1])
            else:
                scale[key] = 1.0
        elif "40_GeV" in key:
            if not skipSigScale:
                scale["Signal_40_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[2])
            else:
                scale[key] = 1.0
        elif "60_GeV" in key:
            if not skipSigScale:
                scale["Signal_60_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[3])
            else:
                scale[key] = 1.0
        elif "background" in key:
            # Do this since we reweight using scaled weights (i.e., no normalization in the reweighting itself)
            scale["background"] = data.Integral() / hists["background"].Integral()

    if "interMass" in path:
        hists["Signal_15_GeV"].Add(hists["Signal_30_GeV"])
        hists["Signal_15_GeV"].Add(hists["Signal_40_GeV"])
        hists["Signal_15_GeV"].Add(hists["Signal_60_GeV"])
        scale["Signal_15_GeV"] = data.Integral() / hists["Signal_15_GeV"].Integral()

    # Main plotting on pad1
    pad1.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    colors = [kAzure-4, kMagenta, kOrange, kRed, kBlue]
    if "interMass" in path:
        # Force only the bkg and combined signal MC to be plotted
        colors = [kAzure-4, kRed]
    for idx, ((key, hist), color) in enumerate(zip(hists.items(), colors)):
        hist.SetTitle("")
        hist.GetYaxis().SetTitleFont(42)
        hist.GetYaxis().SetTitle("Events")
        hist.SetLineColor(color)
        hist.Sumw2()

        if not ignore_scales:
            hist.Scale(scale[key])
        if "GeV" in key:
            hist.SetLineWidth(3)
            if "interMass" not in path:
                legend.AddEntry(hist, "m_{a} = "+key.replace("Signal_","").replace("_", " "), "l")
            else:
                legend.AddEntry(hist, "m_{a} = 15, 30, 40, 60", "l")
        if key == "background":
            hist.SetFillColor(color)
            legend.AddEntry(hist, "Event Mixing", "f")

            # Clone background to include statistical errors
            background_error = hists["background"].Clone("background_error")
            background_error.SetLineColor(kGray+3)
            background_error.SetFillColor(kGray+3)
            background_error.SetFillStyle(3008)
            legend.AddEntry(background_error, "Statistical Uncertainty", "f")

        # Get maximum of histograms
        if maximum is None:
            maximum = 0.0
            if maximum < hist.GetBinContent(hist.GetMaximumBin()):
                maximum = hist.GetBinContent(hist.GetMaximumBin())

    for (key, hist), color in zip(hists.items(), colors):
        hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum*10 if log else maximum*1.5)

    # Force draw order so 15_GeV is on top
    hists["background"].Draw("hist same")
    background_error.Draw("E2 same")

    new_keys = {}
    for key in hists.keys():
        for m in range(15,63,5):
            if f"{m}_GeV" in key:
                new_keys.update({f"Signal_{m}_GeV": key})
                break

    if "interMass" not in path:
        for key in [f"Signal_{m}_GeV" for m in range(15,65,5)]:
            if key in new_keys.keys():
                hists[new_keys[key]].Draw("hist same")
    else:
        hists["Signal_15_GeV"].Draw("hist same")
    legend.Draw()
    
    # Residuals histogram setup
    residuals = residuals_error_prop(hists["Signal_30_GeV" if "interMass" not in path else "Signal_15_GeV"], hists["background"], title, autoScale=False)

    # Plotting residuals on pad2
    pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
    pad2.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    gPad.SetFrameBorderMode(0)
    gPad.SetTickx(1)
    gPad.SetTicky(1)
    if "Delta R" in title:
        residuals.GetXaxis().SetMaxDigits(1)
    residuals.SetTitle("")
    residuals.GetXaxis().SetTitle(title)
    residuals.GetYaxis().SetTitle("Signal / Bkg")
    residuals.GetXaxis().SetTitleSize(9.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
    residuals.GetYaxis().SetTitleSize(8.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
    residuals.GetXaxis().SetTitleOffset(1.1)
    residuals.GetYaxis().SetTitleOffset(0.4)
    residuals.GetYaxis().SetNdivisions(6)
    residuals.GetYaxis().CenterTitle(True)
    residuals.GetXaxis().SetLabelSize(2.5*residuals.GetXaxis().GetLabelSize())
    residuals.GetYaxis().SetLabelSize(2.5*residuals.GetYaxis().GetLabelSize())
    residuals.SetMarkerSize(1)
    residuals.SetMarkerStyle(20)
    residuals.SetMarkerColor(kBlack)
    residuals.SetLineColor(kBlack)
    #residuals.SetFillColor(kGray)
    #residuals.Draw("P E2")
    residuals.Draw("P")

    # Line at y=1
    max_edge = residuals.GetBinCenter(residuals.GetNbinsX())
    min_edge = residuals.GetBinCenter(1)
    line = TLine(min_edge, 1.0, max_edge, 1.0)
    line.SetLineStyle(7)
    line.SetLineColor(kBlack)
    line.Draw("same")

    pad1.Modified()
    pad2.Modified()
    canvas.Update()
    plotFancy(pad1, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, pad=True, sub=sub, Scales=Scales)

    # Save figure
    canvas.SaveAs(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path.split('/')[-1]} to file.")
    
    return hists


def compare_four_masses_data_sig_bkg(
    hists: Dict[str, TH1D],
    title: str,
    path: str,
    genWeight: List[float],
    log: Optional[bool] = True,
    maximum: Optional[float] = None,
    canvas_num: Optional[int] = None,
    ignore_scales: Optional[bool] = False,
    sub: str = None,
    scaleLumi: bool = False,
    year: str = "",
    Scales: str = None,
    use_five: bool = False,
    skipSigScale: bool = False,
    skipBkgScale: bool = False,
    noRatioScaling: bool = False,
) -> Dict[str, TH1D]:
    """
    Plots distributions for pseudoscalar masses (15, 30, 40, and 60 GeV), event mixing (background), and data.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Canvas setup
    canvas = TCanvas(f"canvas_{canvas_num}" if canvas_num is not None else "canvas", "canvas", 1000, 1000)
    
    # Pad setup
    pad1 = TPad(f"pad1_{canvas_num}" if canvas_num is not None else "pad1", "pad1", 0, 0.25, 1, 1)
    pad2 = TPad(f"pad2_{canvas_num}" if canvas_num is not None else "pad2", "pad2", 0, 0, 1, 0.23)
    
    pad1.SetBottomMargin(0.0001)
    pad1.SetBorderMode(0)
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(0.3)
    pad2.SetBorderMode(0)
    if log:
        pad1.SetLogy()

    pad1.Draw()
    pad2.Draw()

    # Legend setup
    if "mHyp" not in path:
        if len(hists) < 12:
            legend = TLegend(0.62,0.58,0.87,0.83)
        elif len(hists) == 12:
            legend = TLegend(0.5,0.57,0.9,0.87)
            legend.SetNColumns(2)
        legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
    else:
        legend_left = TLegend(0.45,0.55,0.68,0.83)
        legend_right = TLegend(0.68,0.55,0.9,0.83)
        for legend in [legend_left, legend_right]:
            legend.SetTextFont(42)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)

    # Handle scaling/max
    scale = {"signal": 1.0, "background": 1.0, "data": 1.0}
    for key, hist in hists.items():
        hist.Sumw2()
        
        # Set scale for signal
        if "15_GeV" in key:
            #scale["Signal_15_GeV"] = Scales.lumi * 1.0 / Scales.sumw[year]["15_GeV"]  # Signal scaled to xs = 1 fb (https://indico.cern.ch/event/1057512/contributions/4449049/attachments/2281868/3879051/H4G_Approval_15Jul2021.pdf)
            if not skipSigScale:
                scale["Signal_15_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[0])  # Signal scaled to xs = 1 fb
            else:
                scale[key] = 1.0
        elif "20_GeV" in key:
            if not skipSigScale:
                scale["Signal_20_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[1])
            else:
                scale[key] = 1.0
        elif "25_GeV" in key:
            if not skipSigScale:
                scale["Signal_25_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[2 if not use_five else 1])
            else:
                scale[key] = 1.0
        elif "30_GeV" in key:
            if not skipSigScale:
                scale["Signal_30_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[3 if "mHyp" in path else 1])
            else:
                scale[key] = 1.0
        elif "35_GeV" in key:
            if not skipSigScale:
                scale["Signal_35_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[4 if not use_five else 2])
            else:
                scale[key] = 1.0
        elif "40_GeV" in key:
            if not skipSigScale:
                scale["Signal_40_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[5 if "mHyp" in path else 2])
            else:
                scale[key] = 1.0
        elif "45_GeV" in key:
            if not skipSigScale:
                scale["Signal_45_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[6 if not use_five else 3])
            else:
                scale[key] = 1.0
        elif "50_GeV" in key:
            if not skipSigScale:
                scale["Signal_50_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[7])
            else:
                scale[key] = 1.0
        elif "55_GeV" in key:
            if not skipSigScale:
                scale["Signal_55_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[8])
            else:
                scale[key] = 1.0
        elif "60_GeV" in key:
            if not skipSigScale:
                scale["Signal_60_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[-1])
            else:
                scale[key] = 1.0
        elif "background" in key:
            # Do this since we reweight using scaled weights (i.e., no normalization in the reweighting itself)
            if not skipBkgScale:
                scale["background"] = hists["data"].Integral() / hists["background"].Integral()
            else:
                scale["background"] = 1.0
        elif "data" in key:
            scale["data"] = 1.0

    # Main plotting on pad1
    pad1.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    colors = [kAzure-4, kBlack, kMagenta, kOrange, kRed, kBlue]
    if len(hists) == 12:
        colors = [kAzure-4, kBlack, kMagenta, kOrange, kRed, kBlue, kGreen, kYellow, kViolet, kCyan, kPink, kGray]
    if "mHyp" in path:
        colors.extend([kGreen+3, kYellow, kViolet, kCyan, kPink, kGray])
        print("mHyp", hists.keys())
    if use_five:
        colors = [kAzure-4, kBlack, kBlue, kOrange, kRed, kMagenta, kGreen+3]
    if "mHyp" not in path:
        legend.AddEntry(hists["data"], "Data", "lep")
        legend.AddEntry(BindObject(0, "TH1D"), "", "")
    for idx, ((key, hist), color) in enumerate(zip(hists.items(), colors)):
        hist.SetTitle("")
        hist.GetYaxis().SetTitleFont(42)
        hist.GetYaxis().SetTitle("Events")
        hist.SetLineColor(color)
        hist.Sumw2()

        if not ignore_scales:
            hist.Scale(scale[key])
        if "GeV" in key:
            hist.SetLineWidth(3)
            if "mHyp" in path:
                if idx < 5:
                    legend_left.AddEntry(hist, "m_{a} = "+key.replace("Signal_","").replace("_", " "), "l")
                else:
                    legend_right.AddEntry(hist, "m_{a} = "+key.replace("Signal_","").replace("_", " "), "l")
            else:
                legend.AddEntry(hist, "m_{a} = "+key.replace("Signal_","").replace("_", " "), "l")
        if key == "background":
            hist.SetFillColor(color)
            if "mHyp" in path:
                legend_left.AddEntry(hist, "Event Mixing", "f")
            else:
                legend.AddEntry(hist, "Event Mixing", "f")

            # Clone background to include statistical errors
            background_error = hists["background"].Clone("background_error")
            background_error.SetLineColor(kGray+3)
            background_error.SetFillColor(kGray+3)
            background_error.SetFillStyle(3008)
            if "mHyp" in path:
                legend_left.AddEntry(background_error, "Stat. Uncertainty", "f")
            else:
                legend.AddEntry(background_error, "Statistical Uncertainty", "f")

        elif key == "data":
            hist.SetMarkerSize(1)
            hist.SetMarkerStyle(20)
            hist.SetMarkerColor(kBlack)
            hist.SetFillStyle(0)
            if "mHyp" in path:
                legend_left.AddEntry(hist, "Data", "lep")
            else:
                pass
            
        # Get maximum of histograms
        if maximum is None:
            maximum = 0.0
            if maximum < hist.GetBinContent(hist.GetMaximumBin()):
                maximum = hist.GetBinContent(hist.GetMaximumBin())

    for (key, hist), color in zip(hists.items(), colors):
        hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum*10 if log else maximum*1.5)
        if "mHyp" in path:
            hist.GetYaxis().SetRangeUser(0.75, maximum*10 if log else maximum*1.5)
        #print(f"Maximum: {maximum}")
        #print("integral", key, hist.Integral())
            

    # Force draw order so 15_GeV is on top
    hists["background"].Draw("hist same")
    background_error.Draw("E2 same")
    hists["data"].Draw("P E1 same")

    new_keys = {}
    for key in hists.keys():
        for m in range(15,63,5):
            if f"{m}_GeV" in key:
                new_keys.update({f"Signal_{m}_GeV": key})
                break

    for key in [f"Signal_{m}_GeV" for m in range(15,65,5)]:
        if key in new_keys.keys():
            hists[new_keys[key]].Draw("hist same")
            #print(new_keys[key], hists[new_keys[key]].Integral(), scale[key])
    if "mHyp" in path:
        legend_left.Draw()
        legend_right.Draw()
    else:
        legend.Draw()
    
    # Residuals histogram setup
    residuals = residuals_error_prop(hists["background"], hists["data"], title, autoScale=False)
    """
    residuals = hists["background"] / hists["data"]
    for b in range(1, int(residuals.GetNbinsX())+1):
        residuals.SetBinError(b, 0.0)
        # errors are more complicated... ask about how tf to do this. it's clearly not sqrt(N). Erorr prop for bkg/data maybe is correct?
    """

    # Plotting residuals on pad2
    pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
    pad2.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    gPad.SetFrameBorderMode(0)
    gPad.SetTickx(1)
    gPad.SetTicky(1)
    if "Delta R" in title:
        residuals.GetXaxis().SetMaxDigits(1)
    residuals.SetTitle("")
    residuals.GetXaxis().SetTitle(title)
    residuals.GetYaxis().SetTitle("Bkg/Data")
    residuals.GetXaxis().SetTitleSize(9.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
    residuals.GetYaxis().SetTitleSize(8.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
    residuals.GetXaxis().SetTitleOffset(1.1)
    residuals.GetYaxis().SetTitleOffset(0.4)
    residuals.GetYaxis().SetNdivisions(6)
    residuals.GetYaxis().CenterTitle(True)
    residuals.GetXaxis().SetLabelSize(2.5*residuals.GetXaxis().GetLabelSize())
    residuals.GetYaxis().SetLabelSize(2.5*residuals.GetYaxis().GetLabelSize())
    residuals.SetMarkerSize(1)
    residuals.SetMarkerStyle(20)
    residuals.SetMarkerColor(kBlack)
    residuals.SetLineColor(kBlack)
    #residuals.SetFillColor(kGray)
    #residuals.Draw("P E2")
    residuals.Draw("P")

    # Line at y=1
    max_edge = residuals.GetBinCenter(residuals.GetNbinsX())
    min_edge = residuals.GetBinCenter(1)
    line = TLine(min_edge, 1.0, max_edge, 1.0)
    line.SetLineStyle(7)
    line.SetLineColor(kBlack)
    line.Draw("same")

    pad1.Modified()
    pad2.Modified()
    canvas.Update()
    plotFancy(pad1, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, pad=True, sub=sub, Scales=Scales)

    # Save figure
    canvas.SaveAs(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path.split('/')[-1]} to file.")

    if not noRatioScaling:
        _, mn, mx = residuals_error_prop(hists["background"], hists["data"], title, returnNewScale=True)
        pad2.cd()
        residuals.GetYaxis().SetRangeUser(mn, mx)
        residuals.Draw("P")

        # Line at y=1
        max_edge = residuals.GetBinCenter(residuals.GetNbinsX())
        min_edge = residuals.GetBinCenter(1)
        line = TLine(min_edge, 1.0, max_edge, 1.0)
        line.SetLineStyle(7)
        line.SetLineColor(kBlack)
        line.Draw("same")

        pad1.Modified()
        pad2.Modified()
        canvas.Update()

        path_scaled = path.replace(".pdf", "_ratioScaled.pdf")
        canvas.SaveAs(path_scaled)
        assert os.path.isfile(path_scaled), "Saved plot does not exist."
        print(f"Saved {path_scaled.split('/')[-1]} to file.")
    
    return hists


def compare_four_masses_data_sig_bkg_reweight_overlay(
    hists: Dict[str, TH1D],
    title: str,
    path: str,
    genWeight: List[float],
    log: Optional[bool] = True,
    maximum: Optional[float] = None,
    canvas_num: Optional[int] = None,
    ignore_scales: Optional[bool] = False,
    sub: str = None,
    scaleLumi: bool = False,
    year: str = "",
    Scales: str = None,
    suppressNonBkg: bool = False
) -> Dict[str, TH1D]:
    """
    Plots distributions for pseudoscalar masses (15, 30, 40, and 60 GeV), event mixing (background), and data.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Canvas setup
    canvas = TCanvas(f"canvas_{canvas_num}" if canvas_num is not None else "canvas", "canvas", 1000, 1000)
    
    # Pad setup
    pad1 = TPad(f"pad1_{canvas_num}" if canvas_num is not None else "pad1", "pad1", 0, 0.25, 1, 1)
    pad2 = TPad(f"pad2_{canvas_num}" if canvas_num is not None else "pad2", "pad2", 0, 0, 1, 0.23)
    
    pad1.SetBottomMargin(0.0001)
    pad1.SetBorderMode(0)
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(0.3)
    pad2.SetBorderMode(0)
    if log:
        pad1.SetLogy()

    pad1.Draw()
    pad2.Draw()

    # Legend setup
    legend = TLegend(0.65,0.65,0.87,0.85)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # Handle scaling/max
    scale = {"signal": 1.0, "background": 1.0, "background_rw": 1.0,  "data": 1.0}
    for key, hist in hists.items():
        hist.Sumw2()
        # Get maximum of histograms
        if maximum is None:
            maximum = 0.0
            if maximum < hist.GetBinContent(hist.GetMaximumBin()):
                maximum = hist.GetBinContent(hist.GetMaximumBin())
        
        # Set scale for signal
        if "15_GeV" in key:
            scale["Signal_15_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[0])  # Signal scaled to xs = 1 fb (https://indico.cern.ch/event/1057512/contributions/4449049/attachments/2281868/3879051/H4G_Approval_15Jul2021.pdf)
        elif "30_GeV" in key:
            scale["Signal_30_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[1])
        elif "40_GeV" in key:
            scale["Signal_40_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[2])
        elif "60_GeV" in key:
            scale["Signal_60_GeV"] = Scales.lumi_fb * 1.0 / (genWeight[3])
        elif "background" in key:
            # Account for background_rw - make it a generic key
            # Do this since we reweight using scaled weights (i.e., no normalization in the reweighting itself)
            scale[key] = Scales.lumi if scaleLumi else hists["data"].Integral() / hists[key].Integral()
        elif "data" in key:
            scale["data"] = Scales.lumi if scaleLumi else 1.0


    # Residuals histogram setup
    residuals_prerw = residuals_error_prop(hists["background"], hists["data"], title)
    residuals_rw = residuals_error_prop(hists["background_rw"], hists["data"], title)
    """
    residuals_prerw = hists["background"] / hists["data"]
    residuals_rw = hists["background_rw"] / hists["data"]
    assert int(residuals_prerw.GetNbinsX()) == int(residuals_rw.GetNbinsX())
    for b in range(1, int(residuals_prerw.GetNbinsX())+1):
        residuals_prerw.SetBinError(b, 0.0)
        residuals_rw.SetBinError(b, 0.0)
        # errors are more complicated... ask about how tf to do this. it's clearly not sqrt(N). Erorr prop for bkg/data maybe is correct?
    """

    # Main plotting on pad1
    pad1.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    colors = [kAzure-4, kGreen+3, kBlack, kMagenta, kOrange, kRed, kBlue]
    for (key, hist), color in zip(hists.items(), colors):
        hist.SetTitle("")
        hist.GetYaxis().SetTitleFont(42)
        hist.GetYaxis().SetTitle("Events")
        hist.SetLineColor(color)
        if not ignore_scales:
            hist.Scale(scale[key])
        if "GeV" in key:
            hist.SetLineWidth(3)
            if not suppressNonBkg:
                legend.AddEntry(hist, "m_{a} = "+key.replace("Signal_","").replace("_", " "), "l")

        if key == "background":
            hist.SetFillColor(color)
            legend.AddEntry(hist, "Event Mixing", "f")

        if key == "background_rw":
            hist.SetMarkerSize(1)
            hist.SetMarkerStyle(20)
            hist.SetMarkerColor(color)
            legend.AddEntry(hist, "Event Mixing Reweighted", "lp")

        elif key == "data":
            hist.SetMarkerSize(1)
            hist.SetMarkerStyle(20)
            hist.SetMarkerColor(color)
            hist.SetFillStyle(0)
            if not suppressNonBkg:
                legend.AddEntry(hist, "Data", "lep")
            
        hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum*10 if log else maximum*1.5)

    # Force draw order so 15_GeV is on top
    hists["background"].Draw("hist same")
    hists["background_rw"].Draw("hist P same")
    if not suppressNonBkg:
        hists["data"].Draw("P E1 same")

    if not suppressNonBkg:
        new_keys = {}
        for key in hists.keys():
            for m in range(15,63,5):
                if f"{m}_GeV" in key:
                    new_keys.update({f"Signal_{m}_GeV": key})
                    break

        hists[new_keys["Signal_60_GeV"]].Draw("hist same")
        hists[new_keys["Signal_40_GeV"]].Draw("hist same")
        hists[new_keys["Signal_30_GeV"]].Draw("hist same")
        hists[new_keys["Signal_15_GeV"]].Draw("hist same")

    # Draw legend
    legend.Draw()
    
    # Plotting residuals on pad2
    pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
    pad2.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    gPad.SetFrameBorderMode(0)
    gPad.SetTickx(1)
    gPad.SetTicky(1)

    for residuals, color in zip([residuals_prerw, residuals_rw], [kAzure-4, kGreen+3]):
        if "Delta R" in title:
            residuals.GetXaxis().SetMaxDigits(1)
        residuals.SetTitle("")
        residuals.GetXaxis().SetTitle(title)
        residuals.GetYaxis().SetTitle("Bkg/Data")
        residuals.GetXaxis().SetTitleSize(9.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
        residuals.GetYaxis().SetTitleSize(8.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
        residuals.GetXaxis().SetTitleOffset(1.1)
        residuals.GetYaxis().SetTitleOffset(0.4)
        residuals.GetYaxis().SetNdivisions(6)
        residuals.GetYaxis().CenterTitle(True)
        residuals.GetXaxis().SetLabelSize(2.5*residuals.GetXaxis().GetLabelSize())
        residuals.GetYaxis().SetLabelSize(2.5*residuals.GetYaxis().GetLabelSize())
        residuals.SetMarkerSize(1)
        residuals.SetMarkerStyle(20)
        residuals.SetMarkerColor(color)
        residuals.SetLineColor(color)
        #residuals.SetFillColor(kGray)
        #residuals.Draw("P E2 same")
        residuals.Draw("P same")

    # Line at y=1
    max_edge = residuals_prerw.GetBinCenter(residuals.GetNbinsX())
    min_edge = residuals_prerw.GetBinCenter(1)
    line = TLine(min_edge, 1.0, max_edge, 1.0)
    line.SetLineStyle(7)
    line.SetLineColor(kBlack)
    line.Draw("same")

    pad1.Modified()
    pad2.Modified()
    canvas.Update()
    plotFancy(pad1, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, pad=True, sub=sub, Scales=Scales)

    # Save figure
    canvas.SaveAs(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path.split('/')[-1]} to file.")

    return hists


def plotCorrectedComparison(
    hists: Dict[str, TH1D],
    title: str,
    path: str,
    genWeight: float = None,
    log: Optional[bool] = True,
    minimum: Optional[float] = 10**-3 + 10**-4,
    maximum: Optional[float] = None,
    canvas_num: Optional[int] = None,
    year: str = "",
    Scales: str = None,
    fix_res_scale: bool = False,
    legend_text: dict = None,
    sim: bool = False,
    df_corr: pd.DataFrame = None,
    df_uncorr: pd.DataFrame = None,
    branch: str = None,
    unsmeared: bool = False
) -> Dict[str, TH1D]:
    """
    Plot comparison of corrected and uncorrected samples for BDT training variables.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Canvas setup
    canvas = TCanvas(f"canvas_{canvas_num}" if canvas_num is not None else "canvas", "canvas", 1000, 1000)
    
    # Pad setup
    pad1 = TPad(f"pad1_{canvas_num}" if canvas_num is not None else "pad1", "pad1", 0, 0.25, 1, 1)
    pad2 = TPad(f"pad2_{canvas_num}" if canvas_num is not None else "pad2", "pad2", 0, 0, 1, 0.23)
    
    pad1.SetBottomMargin(0.0001)
    pad1.SetBorderMode(0)
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(0.3)
    pad2.SetBorderMode(0)
    if log:
        pad1.SetLogy()

    pad1.Draw()
    pad2.Draw()

    # Legend setup
    legend_box = (0.55,0.65,0.87,0.85)
    if not log and "mass_a1" not in path and "mass_a2" not in path and "ma1" not in path and "ma2" not in path:
        legend_box = (0.15,0.65,0.47,0.85)
    legend = TLegend(legend_box[0], legend_box[1], legend_box[2], legend_box[3])
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # Handle scaling/max
    scale = {"background": 1.0, "data": 1.0}
    for key, hist in hists.items():
        hist.Sumw2()
        # Get maximum of histograms
        if maximum is None:
            maximum = 0.0
            if maximum < hist.GetBinContent(hist.GetMaximumBin()):
                maximum = hist.GetBinContent(hist.GetMaximumBin())
        
        # Set scale for bkg/data
        if "bkg" in key:
            scale["background"] = 1.0
        elif "data" in key:
            scale["data"] = 1.0
        elif "sig" in key:
            if genWeight is not None:
                scale["sig"] = Scales.lumi_fb * 1.0 / (genWeight)
            else:
                scale["sig"] = Scales.lumi_fb

    # Main plotting on pad1
    pad1.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    colors = [kBlack, kRed]  # kGreen+3 is for corrected sample (now kRed b/c Nancy asked for it)
    """
    # Reverse draw order if desired
    colors2 = [kRed, kBlack]
    hists2 = {list(hists.keys())[1]: hists[list(hists.keys())[1]], list(hists.keys())[0]: hists[list(hists.keys())[0]]}
    #"""
    for (key, hist), color in zip(hists.items(), colors):
        hist.SetTitle("")
        hist.GetYaxis().SetTitleFont(42)
        hist.GetYaxis().SetTitle("Events")
        hist.SetLineColor(color)
        sc = 1.0
        if "bkg" in key:
            sc = scale["background"]
        elif "data" in key:
            sc = scale["data"]
        elif "sig" in key:
            sc = scale["sig"]

        # Scale to 1.0 b/c it's easier to see
        hist.Scale(1.0)

        hist.SetMarkerSize(1)
        hist.SetMarkerStyle(20)
        hist.SetMarkerColor(color)
        hist.SetFillStyle(0)
        hist.SetLineWidth(2)

        if legend_text is None:
            legt = {
                "bkgUC": "Event Mixing (Uncorrected)",
                "bkgC": "Event Mixing (Corrected)",
                "sigUC": f"Signal {key[-6:].replace('_',' ')} (Uncorrected)",
                "sigC": f"Signal {key[-6:].replace('_',' ')} (Corrected)",
                "dataUC": "Data (Uncorrected)",
                "dataC": "Data (Corrected)",
            }
        else:
            legt = legend_text

        if "bkgUC" in key:
            legend.AddEntry(hist, legt["bkgUC"], "lep")
        elif "bkgC" in key:
            legend.AddEntry(hist, legt["bkgC"], "lep")
        if "dataUC" in key:
            legend.AddEntry(hist, legt["dataUC"], "lep")
        elif "dataC" in key:
            legend.AddEntry(hist, legt["dataC"], "lep")
        if "sigUC" in key:
            legend.AddEntry(hist, legt["sigUC"], "lep")
        elif "sigC" in key:
            legend.AddEntry(hist, legt["sigC"], "lep")
            
        hist.GetYaxis().SetRangeUser(minimum, maximum*10 if log else maximum*1.2)
        if "mass_gggg" in path and "dR_aa" not in path and not log:
            hist.GetYaxis().SetRangeUser(minimum, 650.0)
        #if not log:
        #    TGaxis.SetMaxDigits(3)

        hist.Draw("P E1 same")        
        if "pT" in path:
            hist.Draw("hist same")        
    legend.Draw()
    
    # Residuals histogram setup
    uc, c = tuple(hists.keys())
    residuals = residuals_error_prop(
        hists[c],
        hists[uc],
        title,
        corrections=fix_res_scale,
        autoScale=False
    )

    # Plotting residuals on pad2
    pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
    pad2.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    gPad.SetFrameBorderMode(0)
    gPad.SetTickx(1)
    gPad.SetTicky(1)
    if "Delta R" in title:
        residuals.GetXaxis().SetMaxDigits(1)
    residuals.SetTitle("")
    residuals.GetXaxis().SetTitle(title)
    residuals.GetYaxis().SetTitle("Corrected/Uncorrected")
    residuals.GetXaxis().SetTitleSize(9.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
    residuals.GetYaxis().SetTitleSize(8.0*hist.GetYaxis().GetTitleSize() * pad1Scale)
    residuals.GetXaxis().SetTitleOffset(1.1)
    residuals.GetYaxis().SetTitleOffset(0.4)
    residuals.GetYaxis().SetNdivisions(6)
    residuals.GetYaxis().CenterTitle(True)
    if "pT" in path:
        residuals.GetYaxis().SetRangeUser(0.94, 1.06)
    if "mass_a" in path or "ma1" in path or "ma2" in path:
        residuals.GetYaxis().SetRangeUser(0.0, 2.0)
    #"""
    if "mass_gggg" in path or "m4g" in path:
        residuals.GetYaxis().SetRangeUser(0.0, 2.0)
    #"""
    residuals.GetXaxis().SetLabelSize(2.5*residuals.GetXaxis().GetLabelSize())
    residuals.GetYaxis().SetLabelSize(2.5*residuals.GetYaxis().GetLabelSize())
    residuals.SetMarkerSize(1)
    residuals.SetMarkerStyle(20)
    residuals.SetMarkerColor(kBlack)
    residuals.SetLineColor(kBlack)
    residuals.Draw("P")

    # Line at y=1
    max_edge = residuals.GetBinCenter(residuals.GetNbinsX())
    min_edge = residuals.GetBinCenter(1)
    line = TLine(min_edge, 1.0, max_edge, 1.0)
    line.SetLineStyle(7)
    line.SetLineColor(kBlack)
    line.Draw("same")

    pad1.Modified()
    pad2.Modified()
    canvas.Update()
    plotFancy(pad1, title.replace("[GeV]","").replace("[GeV^{-1}]",""), lumiTxt=Scales.lumiTxt, prelim=True if not sim else False, simulation=sim, inPlot=False, pad=True, Scales=Scales)

    # Save figure
    canvas.SaveAs(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path.split('/')[-1]} to file.")

    if df_corr is not None and df_uncorr is not None:
        # Make TH2D with events on x/y-axes and mass on z-axis. Can use TH1s for this.
        # Do it again with mass on x/y-axes and events on z-axis. Have to use the sample dataframe for this.
        # Doesn't matter that "hist" is an old variable since the binning is consistent
        hist_zEvents = TH2D("hist_zEvents", "hist_zEvents", hist.GetNbinsX(), 0, 650, hist.GetNbinsX(), 0, 650)
        for b in range(1,hists[list(hists.keys())[0]].GetNbinsX()+1):
            hist_zEvents.Fill(hists[list(hists.keys())[0]].GetBinContent(b), hists[list(hists.keys())[1]].GetBinContent(b))

        # Use a similar granularity at Nancy's request. Needs finer binning for ma1/ma2 to look comparatively similar granularity
        binsScale = 10 if "mass_gggg" in path else 20
        hist_zMass = TH2D("hist_zMass", "hist_zMass", int((hist.GetBinLowEdge(hist.GetNbinsX()) - hist.GetBinLowEdge(1)) * binsScale), hist.GetBinLowEdge(1), hist.GetBinLowEdge(hist.GetNbinsX()), int((hist.GetBinLowEdge(hist.GetNbinsX()) - hist.GetBinLowEdge(1)) * binsScale), hist.GetBinLowEdge(1), hist.GetBinLowEdge(hist.GetNbinsX()))
        hist_proj = TH2D("hist_zMass", "hist_zMass", hist.GetNbinsX(), hist.GetBinLowEdge(1), hist.GetBinLowEdge(hist.GetNbinsX()), hist.GetNbinsX(), hist.GetBinLowEdge(1), hist.GetBinLowEdge(hist.GetNbinsX()))
        rng = 0
        if len(df_uncorr) > len(df_corr):
            rng = len(df_corr)
            print(f"Using smaller of the two number of events for plots! df_uncorr: {len(df_uncorr)} events")
        elif len(df_corr) > len(df_uncorr):
            rng = len(df_uncorr)
            print(f"Using smaller of the two number of events for plots! df_corr: {len(df_corr)} events")
        else:
            rng = len(df_uncorr)
        for idx in range(rng):
            hist_zMass.Fill(df_uncorr[branch].iloc[idx], df_corr[branch].iloc[idx], df_uncorr["weight"].iloc[idx])
            hist_proj.Fill(df_uncorr[branch].iloc[idx], df_corr[branch].iloc[idx], df_uncorr["weight"].iloc[idx])

        over = under = 0
        for i in range(hist_proj.GetNbinsX()):
            for j in range(hist_proj.GetNbinsY()):
                if hist_proj.IsBinUnderflow(hist_proj.GetBin(i,j)):
                    under += 1
                elif hist_proj.IsBinOverflow(hist_proj.GetBin(i,j)):
                    over += 1
        print(f"2D fill #: {rng}  2D entries: {hist_proj.GetEntries()}  underflow {under}  overflow {over}")

        canvas.Clear()
        hist_zEvents.SetStats(0)
        hist_zEvents.SetMarkerStyle(8)
        hist_zEvents.SetMarkerSize(0.5)
        hist_zEvents.SetTitle("")
        hist_zEvents.GetYaxis().SetTitleFont(42)
        hist_zEvents.GetXaxis().SetTitleFont(42)
        hist_zEvents.GetXaxis().SetTitleOffset(1.5)
        hist_zEvents.GetXaxis().CenterTitle(True)
        hist_zEvents.GetYaxis().CenterTitle(True)
        hist_zEvents.GetXaxis().SetTitle("Recalculated Masses" if unsmeared else "No smearing")
        hist_zEvents.GetYaxis().SetTitle("Corrected Events")
        hist_zEvents.Draw("COLZ")

        plotFancy(canvas, "Frequency of Events per Bin", prelim=True, simulation=True, inPlot=False, colz=True)
        canvas.SaveAs(path.replace(".pdf","_2D_events.pdf"))

        canvas.Clear()
        hist_zMass.SetStats(0)
        hist_zMass.SetMarkerStyle(8)
        hist_zMass.SetMarkerSize(0.5)
        hist_zMass.SetTitle("")
        hist_zMass.GetYaxis().SetTitleFont(42)
        hist_zMass.GetXaxis().SetTitleFont(42)
        hist_zMass.GetXaxis().SetTitleOffset(1.5)
        hist_zMass.GetXaxis().CenterTitle(True)
        hist_zMass.GetYaxis().CenterTitle(True)
        hist_zMass.GetXaxis().SetTitle(f"Recalculated Masses {title}" if unsmeared else f"No smearing {title}")
        hist_zMass.GetYaxis().SetTitle(f"Corrected {title}")
        hist_zMass.Draw("COLZ")

        line = TLine(hist.GetBinLowEdge(1), hist.GetBinLowEdge(1), hist.GetBinLowEdge(hist.GetNbinsX()), hist.GetBinLowEdge(hist.GetNbinsX()))
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.SetLineColor(kBlack)
        line.Draw("same")

        plotFancy(canvas, title, prelim=True, simulation=True, inPlot=False, colz=True)
        canvas.SaveAs(path.replace(".pdf","_2D_mass.pdf"))

        # Don't do projection plots here! Use proj_sanity_check.py instead! Somehow, things are mismatched here...

    return hists


def plot_smooth_all(
    path: str,
    mass: str,
    sig: tuple,
    bkg: tuple,
    data: tuple,
    Scales: str = None,
) -> None:
    """
    Plots histograms and TGraphSmooth objects for BDT score distributions.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Canvas setup
    canvas = TCanvas("canvas", "canvas", 1000, 1000)
    canvas.SetLogy()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    maxY = 10**5 - 10**4
    minY = 10**-3 + 10**-4

    # Legend setup
    lsize = [0.65, 0.68, 0.92, 0.83]
    legend = TLegend(lsize[0], lsize[1], lsize[2], lsize[3])
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    for hist, graph, hcolor, gcolor, ltext, lopt, dopt in [bkg, data, sig]:
        # Hist and graph setup
        if graph is not None:
            graph.SetTitle("")
            graph.SetLineColor(gcolor)
            graph.SetLineWidth(2)
        hist.SetTitle("")
        hist.GetYaxis().SetTitleFont(42)
        hist.GetYaxis().SetTitle("Events")
        hist.GetXaxis().SetTitle("BDT Score")
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(0.5)
        hist.SetMarkerColor(hcolor)
        hist.SetLineColor(hcolor)
        hist.SetLineWidth(2)
        hist.GetYaxis().SetRangeUser(minY, maxY)

        if graph is not None:
            legend.AddEntry(graph, f"Smoothed {ltext}", "l")
        legend.AddEntry(hist, ltext, lopt)

        # Draw hist, graph, and legend
        hist.Draw(dopt)
        if graph is not None:
            graph.DrawClone("LP")

    legend.Draw()
    plotFancy(canvas, "Smoothed BDT Scores", lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, sub=f"m_{{a}} = {mass.replace('_',' ')}", Scales=Scales, shiftLeft=0.05)

    # Save canvas to file
    canvas.SaveAs(path)

    # Cleanup
    canvas.Destructor()


def plot_smooth(
    bkg_smooth: TGraph,
    bkg_hist: TH1D,
    sig_hist: TH1D,
    data_hist: TH1D,
    mass: str,
    path: str,
    cats: List[float],
    l_edge: float,
    Scales: str = None
) -> None:
    """
    Plots TGraphSmooth for BDT scores.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Canvas setup
    canvas = TCanvas("canvas", "canvas", 1000, 1000)
    canvas.SetLogy()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)

    # Hist and graph setup
    maxY = 10**5 - 10**4
    minY = 10**-3 + 10**-4
    bkg_smooth.SetTitle("")
    bkg_smooth.SetLineColor(kBlue)
    bkg_smooth.SetLineWidth(2)
    bkg_hist.SetTitle("")
    bkg_hist.GetYaxis().SetTitleFont(42)
    bkg_hist.GetYaxis().SetTitle("Events")
    bkg_hist.GetXaxis().SetTitle("BDT Score")
    bkg_hist.SetLineColor(kGreen+3)
    bkg_hist.SetLineWidth(2)
    bkg_hist.GetYaxis().SetRangeUser(minY, maxY)
    sig_hist.SetTitle("")
    sig_hist.GetYaxis().SetTitleFont(42)
    sig_hist.GetYaxis().SetTitle("Events")
    sig_hist.GetXaxis().SetTitle("BDT Score")
    sig_hist.SetLineColor(kRed)
    sig_hist.SetFillColor(kRed)
    sig_hist.SetLineWidth(3)
    sig_hist.GetYaxis().SetRangeUser(minY, maxY)
    data_hist.SetTitle("")
    data_hist.GetYaxis().SetTitleFont(42)
    data_hist.GetYaxis().SetTitle("Events")
    data_hist.GetXaxis().SetTitle("BDT Score")
    data_hist.SetMarkerColor(kBlack)
    data_hist.SetLineColor(kBlack)
    data_hist.SetMarkerStyle(20)
    data_hist.SetMarkerSize(1)
    data_hist.GetYaxis().SetRangeUser(minY, maxY)

    # Legend setup
    lsize = [0.65, 0.68, 0.92, 0.83]
    legend = TLegend(lsize[0], lsize[1], lsize[2], lsize[3])
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(bkg_smooth, "Smoothed EM Background", "l")
    legend.AddEntry(bkg_hist, "EM Background", "l")
    legend.AddEntry(sig_hist, "m_{a} = "+f"{mass.replace('_',' ')}", "l")
    legend.AddEntry(data_hist, "Data", "ple")

    # Adjust x-range if necessary
    sig_hist.GetXaxis().SetRangeUser(l_edge, 1.0)
    bkg_hist.GetXaxis().SetRangeUser(l_edge, 1.0)
    data_hist.GetXaxis().SetRangeUser(l_edge, 1.0)

    # Draw hist, graph, and legend
    sig_hist.Draw("hist same")
    bkg_hist.Draw("hist same")
    data_hist.Draw("P E1 same")
    bkg_smooth.DrawClone("LP")
    legend.Draw()

    # Draw category lines
    for cat in cats[-1]:
        if cat > l_edge:
            line = TLine(cat, minY, cat, maxY)
            line.SetLineColor(kAzure-4)
            line.SetLineWidth(2)
            line.Draw("same")

    plotFancy(canvas, "Smoothed BDT Score", lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, sub="Sideband Only", Scales=Scales)

    # Save canvas to file
    canvas.SaveAs(f"{cwd}/../../scripts/plots/BDT/{Scales.yr}/cats/{path}")

    # Cleanup
    canvas.Destructor()


def plot_reweighting(
    rw_hist: TH1D,
    pre_rw_hist: TH1D,
    target_hist: TH1D,
    rw_title: str,
    pre_rw_title: str,
    target_title: str,
    path: str,
    Scales: str = None
) -> None:
    """
    Plot test of reweighting procedure.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Save plots of ratio histograms
    canvas = TCanvas("canvas", "canvas", 1000, 1000)

    # Pad setup
    pad1 = TPad("pad1", "pad1", 0, 0.25, 1, 1)
    pad2 = TPad("pad2", "pad2", 0, 0, 1, 0.23)
    
    pad1.SetBottomMargin(0.0001)
    pad1.SetBorderMode(0)
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(0.3)
    pad2.SetBorderMode(0)
    pad1.SetLogy()

    pad1.Draw()
    pad2.Draw()

    maximum = 0.0
    # Get maximum of histograms
    for hist in [rw_hist, target_hist, pre_rw_hist]:
        if maximum < hist.GetBinContent(hist.GetMaximumBin()):
            maximum = hist.GetBinContent(hist.GetMaximumBin())

    # Background histogram setup
    rw_hist.GetYaxis().SetRangeUser(10**-3 + 10**-4, maximum + 10**5)
    rw_hist.SetTitle("")
    rw_hist.GetYaxis().SetTitleFont(42)
    rw_hist.GetYaxis().SetTitle("Events")
    rw_hist.GetXaxis().SetTitle("BDT Score")
    rw_hist.SetLineColor(kBlue)
    rw_hist.SetLineWidth(2)

    # Data histogram setup
    target_hist.SetTitle("")
    target_hist.GetYaxis().SetTitleFont(42)
    target_hist.SetMarkerColor(kBlack)
    target_hist.SetLineColor(kBlack)
    target_hist.SetMarkerStyle(20)
    target_hist.SetMarkerSize(1)
    target_hist.SetLineWidth(2)

    # Pre-reweight background histogram setup
    pre_rw_hist.SetTitle("")
    pre_rw_hist.GetYaxis().SetTitleFont(42)
    pre_rw_hist.SetLineColor(kGreen+3)
    pre_rw_hist.SetLineWidth(2)

    # Legend setup
    lsize = [0.65, 0.68, 0.87, 0.83]
    legend = TLegend(lsize[0], lsize[1], lsize[2], lsize[3])
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(rw_hist, rw_title, "l")
    legend.AddEntry(pre_rw_hist, pre_rw_title, "l")
    legend.AddEntry(target_hist, target_title, "ple")

    # Draw everything
    pad1.cd()
    rw_hist.Draw("hist")
    target_hist.Draw("P E1 same")
    pre_rw_hist.Draw("hist same")
    legend.Draw()

    # Residuals histogram setup
    rw_res = residuals_error_prop(rw_hist, target_hist)
    pre_rw_res = residuals_error_prop(pre_rw_hist, target_hist)
    """
    rw_res = rw_hist / target_hist
    pre_rw_res = pre_rw_hist / target_hist

    # Formatting trick for residuals
    for b in range(int(rw_hist.GetNbinsX())):
        # rw_res.SetBinError(b, rw_hist.GetBinError(b))
        #rw_res.SetBinError(b, 0.001)
        rw_res.SetBinError(b, 0.0)
    
    for b in range(int(pre_rw_hist.GetNbinsX())):
        # pre_rw_res.SetBinError(b, pre_rw_hist.GetBinError(b))
        #pre_rw_res.SetBinError(b, 0.001)
        pre_rw_res.SetBinError(b, 0.0)
    """

    # Plotting rw_res/pre_rw_res on pad2
    pad1Scale = (pad2.GetWNDC() * pad2.GetHNDC()) / (pad1.GetWNDC() * pad1.GetHNDC())
    pad2.cd()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    gPad.SetFrameBorderMode(0)
    gPad.SetTickx(1)
    gPad.SetTicky(1)
    rw_res.SetTitle("")
    rw_res.GetYaxis().SetTitle("Bkg/Data")
    rw_res.GetXaxis().SetTitle("BDT Score")
    rw_res.GetXaxis().SetTitleSize(2.5*rw_res.GetYaxis().GetLabelSize())
    #rw_res.GetYaxis().SetTitleSize(2.5*rw_res.GetYaxis().GetLabelSize())
    rw_res.GetXaxis().SetTitleSize(9.0*rw_hist.GetYaxis().GetTitleSize() * pad1Scale)
    rw_res.GetYaxis().SetTitleSize(8.0*rw_hist.GetYaxis().GetTitleSize() * pad1Scale)
    rw_res.GetXaxis().SetTitleOffset(1.1)
    rw_res.GetYaxis().SetTitleOffset(0.4)
    rw_res.GetYaxis().CenterTitle(True)
    rw_res.GetYaxis().SetNdivisions(6)
    rw_res.GetXaxis().SetLabelSize(2.5*rw_res.GetXaxis().GetLabelSize())
    rw_res.GetYaxis().SetLabelSize(2.5*rw_res.GetYaxis().GetLabelSize())
    rw_res.SetLineColor(kGreen+3)
    rw_res.SetMarkerColor(kGreen+3)
    rw_res.SetMarkerStyle(20)
    rw_res.SetMarkerSize(0.5)

    pre_rw_res.GetYaxis().SetTitle("Bkg/Data")
    pre_rw_res.GetXaxis().SetTitle("BDT Score")
    pre_rw_res.GetXaxis().SetTitleSize(2.5*pre_rw_res.GetYaxis().GetLabelSize())
    pre_rw_res.GetYaxis().SetNdivisions(6)
    pre_rw_res.GetXaxis().SetLabelSize(2.5*pre_rw_res.GetXaxis().GetLabelSize())
    pre_rw_res.GetYaxis().SetLabelSize(2.5*pre_rw_res.GetYaxis().GetLabelSize())
    pre_rw_res.SetLineColor(kBlue)
    pre_rw_hist.SetLineWidth(2)

    pre_rw_res.Draw("E same")
    rw_res.Draw("PE same")

    pad1.Modified()
    pad2.Modified()
    canvas.Update()

    # Plot fancy
    plotFancy(pad1, f"BDT Score for m_{{a}} = {path[path.find('GeV')-3: path.find('GeV')+3].replace('_',' ')}", lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, pad=True, sub="110 < m_{#gamma#gamma#gamma#gamma} < 180 GeV", Scales=Scales)

    # Save canvas to file
    canvas.SaveAs(f"{cwd}/../../scripts/plots/BDT/{Scales.yr}/cats/{path}")

    # Cleanup
    canvas.Destructor()


def plot_Nreweighting(
    ratios: List[TH1D],
    titles: List[str],
    ytitle: str,
    Scales: str = None
) -> None:
    """
    Plot ratio plots for N dim reweighting.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Save plots of ratio histograms
    canvas = TCanvas("canvas", "canvas", 1000, 1000)
    for title, ratio in zip(titles, ratios):
        ratio.SetTitle("")
        ratio.GetYaxis().SetTitleFont(42)
        ratio.GetYaxis().SetTitle(ytitle)
        ratio.GetXaxis().SetTitle(title)
        ratio.SetLineColor(kBlack)

        # Draw and make fancy
        ratio.Draw("hist")
        plotFancy(canvas, title, lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, Scales=Scales)

        # Save plot to file
        canvas.SaveAs(f"{cwd}/../../scripts/plots/BDT/{Scales.yr}/ratios/{ratio.GetName().replace('/','_').replace('data_','')}.pdf")
    
    # Cleanup
    canvas.Destructor()


def plot_Ndim_weights(
    weights: TH1D,
    Scales = None
) -> None:
    """
    Plot weight multipliers from N dim reweighting.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Save plots of ratio histograms
    canvas = TCanvas("canvas_Ndim_weights", "canvas_Ndim_weights", 1000, 1000)


    title = f"N-dim Reweighting Factors {Scales.yr}"
    weights.SetTitle("")
    weights.GetYaxis().SetTitleFont(42)
    weights.GetYaxis().SetTitle("Multiplicative Weight Factor")
    weights.GetXaxis().SetTitle(f"Bins ({int(weights.GetNbinsX())})")
    weights.SetLineColor(kBlack)
    weights.GetYaxis().SetRangeUser(-0.3, 2.0)

    # Draw and make fancy
    weights.Draw("P E1")
    plotFancy(canvas, title, lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, Scales=Scales)

    # Save plot to file
    canvas.SaveAs(f"{cwd}/../../scripts/plots/BDT/{Scales.yr}/ratios/weights_Ndim.pdf")
    
    # Cleanup
    canvas.Destructor()


def plot_AMS(
    AMS: List,
    path: str,
    year: str
) -> None:
    """
    Plots AMS versus number of categories.
    """

    # Get directory of file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Plotting setup
    nCats = len(AMS)
    x = [0.0]
    y = [0.0]

    for idx, ams in enumerate(AMS):
        x.append(idx+1)
        y.append(fsum(ams))
    
    # Plot AMS vs number of categories
    plt.figure()
    plt.plot(x, y, "-o")

    plt.xlabel("Number of categories", fontsize=12)
    plt.ylabel("AMS", fontsize=12)
    plt.grid()

    plt.locator_params(axis="x", integer=True)

    # Save figure
    savePath = f"{cwd}/../../scripts/plots/BDT/{year}/cats/{path}.pdf"
    plt.savefig(savePath)
    assert os.path.isfile(savePath), "Saved plot does not exist."

    plt.close()


def plot_even_odd_BDT(
    path: str,
    evens: TH1D,
    odds: TH1D,
    even_scale: int,
    odd_scale: int,
    mass: int,
    year: str,
    Scales: str = None,
) -> None:
    """
    Plots BDT score of even and odd events.
    """

    # Get directory of file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Canvas setup
    canvas = TCanvas("canvas", "canvas", 1000, 1000)
    canvas.SetLogy()

    # Legend setup
    legend = TLegend(0.2,0.75,0.45,0.86)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # Scaling
    evens.Scale(1.0 / even_scale)
    odds.Scale(1.0 / odd_scale)

    # Histogram setup
    maxY = 10**1 - 10**0
    minY = 10**-3 + 10**-4
    evens.SetTitle("")
    evens.GetYaxis().SetTitleFont(42)
    evens.GetYaxis().SetTitle("Events")
    evens.GetXaxis().SetTitle("BDT Score")
    evens.GetYaxis().SetRangeUser(minY, maxY)
    evens.SetLineWidth(4)
    evens.SetLineColor(kOrange)
    legend.AddEntry(evens, "Even Events", "l")

    odds.SetTitle("")
    odds.GetYaxis().SetTitleFont(42)
    odds.GetYaxis().SetTitle("Events")
    odds.GetXaxis().SetTitle("BDT Score")
    odds.GetYaxis().SetRangeUser(minY, maxY)
    odds.SetLineWidth(4)
    odds.SetLineColor(kBlue)
    legend.AddEntry(odds, "Odd Events", "l")

    odds.Draw("hist")
    evens.Draw("hist same")
    legend.Draw()

    plotFancy(canvas, f"BDT Predictions for m_{{a}} = {mass}", prelim=True, inPlot=False, sub="Predictions on Signal MC", Scales=Scales)#"110 < m_{#gamma#gamma#gamma#gamma} < 180 GeV")

    # Save canvas to file
    canvas.SaveAs(path)

    # Cleanup
    canvas.Destructor()


def plotEfficiency(
    samples: Dict[str, pd.DataFrame],
    year: str
) -> None:
    """
    Plots efficiency of cuts on signal samples for all mass points.
    """

    # Get directory of file
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Fill scatter arrays
    Nevents = {
        "initial_events": [],
        "HLT Flag": [],
        "Pre-selections": [],
        "Pre-selections + At least 4 photons": [],
        "Pre-selections + At least 4 photons + Photon selections": []
    }
    masses = []
    efficiencies = {}

    sorted_samples = sorted(samples.items())
    for mass, sample in sorted_samples:
        # Issue here when doing 2022 merged for some reason... Probably duplicate columns
        masses.append(int(mass[7:9]))
        mass_in = int(mass[7:9])
        efficiencies.update({mass_in: Nevents.copy()})
        for i, key in enumerate(Nevents.keys()):
            # Sum of genWeight should be the same for full samples!
            N = np.sum(sample.Nevents.dropna(how="any")[0::5])
            Ni = np.sum(sample.Nevents.dropna(how="any")[i::5])
            efficiencies[mass_in][key] = Ni / N * 100.0

    plt.figure()
    label_cut = 1
    for eff in list(Nevents.keys())[label_cut:2] + list(Nevents.keys())[3:]:
        x = [m for m in range(15,65,5)]
        y = [efficiencies[m][eff] for m in range(15,65,5)]
        plt.plot(x, y, "-o")
    plt.xlabel("m(a) [GeV]", fontsize=12)
    plt.ylabel("Selection Efficiency (%)", fontsize=12)
    plt.legend(labels=list(Nevents.keys())[label_cut:2] + list(Nevents.keys())[3:], frameon=True, loc="upper left")
    # plt.xticks(masses)
    plt.grid()

    ax = plt.gca()
    ax.set_ylim(0.0, 75.0)

    # Save figure
    path = f"{cwd}/../../scripts/plots/standard/{year}/signal_efficiency_{year}.pdf"
    plt.savefig(path)
    assert os.path.isfile(path), "Saved plot does not exist."
    print(f"Saved {path.split('/')[-1]} to file.")
    plt.close()


def plotLimits(
    limits: dict,
    bounds: Optional[List[int]] = [15,60],
    binsScale: Optional[int] = 1.0,
    canvasX: Optional[int] = 1000,
    canvasY: Optional[int] = 1000,
    year: int = "2018",
    Scales: str = None
) -> None:
    """
    Plots limits for all nominal mass points.
    """

    cwd = os.path.dirname(os.path.abspath(__file__))

    N = len(limits)
    x = [m for m in range(15,65,5)]
    y = [[] for _ in range(6)]

    # Fill graph lists
    #print(f"Scale: {scale}")
    for mass in range(15,65,5):
        scale = Scales.ggH_xs * Scales.BR[year][f"{mass}_GeV"]
        ml = limits[mass]
        y[0].append(ml["obs"] * scale)
        y[1].append(ml["exp0"] * scale)
        y[2].append(ml["exp-1"] * scale)
        y[3].append(ml["exp+1"] * scale)
        y[4].append(ml["exp-2"] * scale)
        y[5].append(ml["exp+2"] * scale)
        print(f"~~~~~~~~  {mass} GeV  ~~~~~~~~")
        print(f"Scale: {scale}")
        print(f"Obs: {y[0][-1]:.5f}", f"Exp: {y[1][-1]:.5f}", f"Exp -1s: {y[2][-1]:.5f}", f"Exp +1s: {y[3][-1]:.5f}", f"Exp -2s: {y[4][-1]:.5f}", f"Exp +2s: {y[5][-1]:.5f}")

    # Create graphs
    axisHist = TH1D("axisHist", "axisHist", int((bounds[1]-bounds[0])*binsScale), bounds[0], bounds[1]) 
    obs = TGraph(N, array("d", x), array("d", y[0]))
    exp0 = TGraph(N, array("d", x), array("d", y[1]))
    # Note: flipped y[2] & y[4] since it plots correctly when they're swapped... wtf?
    #exp1 = TGraphAsymmErrors(N, array("d", x), array("d", y[1]), array("d", [0]*N), array("d", [0]*N), array("d", y[2]), array("d", y[3]))
    #exp2 = TGraphAsymmErrors(N, array("d", x), array("d", y[1]), array("d", [0]*N), array("d", [0]*N), array("d", y[4]), array("d", y[5]))
    exp1 = TGraphAsymmErrors(N, array("d", x), array("d", y[1]), array("d", [0]*N), array("d", [0]*N), array("d", y[4]), array("d", y[3]))
    exp2 = TGraphAsymmErrors(N, array("d", x), array("d", y[1]), array("d", [0]*N), array("d", [0]*N), array("d", y[2]), array("d", y[5]))

    canvas = TCanvas("limits_canvas", "limits_canvas", canvasX, canvasY)
    
    # Set styles
    axisHist.SetTitle("")
    axisHist.SetStats(0)
    axisHist.GetYaxis().SetRangeUser(-1.0, 15)
    axisHist.GetYaxis().SetTitle("#sigma(pp#rightarrow H) #times BR(H#rightarrow aa#rightarrow #gamma#gamma#gamma#gamma)  (fb)")
    axisHist.GetYaxis().SetTitleOffset(1.5)
    axisHist.GetYaxis().SetLabelOffset(0.01)
    axisHist.GetYaxis().CenterTitle(True)
    axisHist.GetXaxis().SetTitle("m_{a} [GeV]")
    axisHist.GetXaxis().SetTitleOffset(1.2)
    axisHist.GetXaxis().SetLabelOffset(0.01)
    obs.SetLineColor(kBlack)
    obs.SetLineWidth(2)
    exp0.SetLineWidth(2)
    exp0.SetLineColor(kBlack)
    exp0.SetLineStyle(2)
    exp1.SetFillColor(kGreen+1)
    exp2.SetFillColor(kOrange)

    # Legend setup
    legend = TLegend(0.68,0.70,0.90,0.82)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(obs, "Observed", "pl")
    legend.AddEntry(exp0, "Median Expected", "l")
    legend.AddEntry(exp1, "#pm#sigma Expected", "f")
    legend.AddEntry(exp2, "#pm2#sigma Expected", "f")

    # Plotting
    axisHist.Draw("AXIS")
    exp2.Draw("E3 same")
    exp1.Draw("E3 same")
    exp0.Draw("L same")
    obs.Draw("L same")
    legend.Draw("same")

    plotFancy(canvas, "", True, lumiTxt=Scales.lumiTxt, inPlot=False, sub="95% CL Upper Limits", Scales=Scales)
    canvas.SaveAs(f"{cwd}/../../scripts/plots/limits.pdf")


def make_scatter_mva(
    samples: dict,
    year: str,
    bdt_cut: bool = False,
    cut: int = 0
) -> None:
    """
    Make scatter plot of MVA ID for photon 3 and 4.
    """

    cwd = os.path.dirname(os.path.abspath(__file__))

    def do_scatter(x, y, xtitle, ytitle, title, xlim, new_path, figNum):
        from matplotlib.colors import LinearSegmentedColormap

        # "Viridis-like" colormap with white background
        white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
            (0, '#ffffff'),
            (1e-20, '#440053'),
            (0.2, '#404388'),
            (0.4, '#2a788e'),
            (0.6, '#21a784'),
            (0.8, '#78d151'),
            (1, '#fde624'),
        ], N=256)

        import mpl_scatter_density
        fig = plt.figure(figNum, figsize=(14,14))
        ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
        density = ax.scatter_density(x, y, cmap=white_viridis, dpi=25)
        fig.colorbar(density, label='Number of points per dot')

        plt.ylabel(ytitle, fontsize=20)
        plt.xlabel(xtitle, fontsize=26)
        plt.title(f"{year}  (N = {len(y)})", fontsize=26)
        ax.set_xlim([-1, xlim])
        ax.set_ylim([-1.1, 1.1])

        # Save figure
        fig.savefig(new_path)
        assert os.path.isfile(new_path), "Saved plot does not exist."
        print(f"Saved {new_path.split('/')[-1]} to file.")
        plt.close(fig)


    for idx, key in enumerate(samples.keys()):
        sample_type = ""
        if "Signal" in key:
            sample_type = "signal"
        elif "data" in key:
            sample_type = "data"
        elif "background" in key:
            sample_type = "background"

        path = f"{cwd}/../../scripts/plots/BDT/{year}/mva_correlation/{sample_type}_mva_scatter.pdf"
        if bdt_cut:
            path = f"{cwd}/../../scripts/plots/BDT/{year}/cuts/mva_correlation/{sample_type}/{key.replace('Signal_','')}_mva_scatter.pdf"

        #samples[key] = samples[key][(samples[key].pho3_mvaID >= -0.85) | (samples[key].pho4_mvaID >= -0.85)]
        mva1 = samples[key].pho1_mvaID
        mva2 = samples[key].pho2_mvaID
        mva3 = samples[key].pho3_mvaID
        mva4 = samples[key].pho4_mvaID
        pt_a1 = samples[key].LeadPs_pt
        pt_a2 = samples[key].SubleadPs_pt
        assert len(mva3) == len(mva4) == len(pt_a1) == len(pt_a2)

        fig = plt.figure(idx, figsize=(14,14))
        plt.scatter(mva3, mva4, s=1/2)

        plt.xlabel("Photon 3 MVA ID", fontsize=20)
        plt.ylabel("Photon 4 MVA ID", fontsize=20)
        plt.title(f"{year}  (N = {len(mva3)}{'  BDT Score > '+str(cut) if bdt_cut else ''})", fontsize=20)

        # Save figure
        fig.savefig(path)
        assert os.path.isfile(path), "Saved plot does not exist."
        print(f"Saved {path.split('/')[-1]} to file.")
        plt.close(fig)


        # Do MVA correlation with cut on > -0.9
        """
        mva_cut = -0.9
        samples[key] = samples[key][(samples[key].pho3_mvaID > mva_cut) & (samples[key].pho4_mvaID > mva_cut)]
        mva3 = samples[key].pho3_mvaID
        mva4 = samples[key].pho4_mvaID
        assert len(mva3) == len(mva4)

        fig = plt.figure(idx+100, figsize=(14,14))
        plt.scatter(mva3, mva4, s=1/2)

        plt.xlabel("Photon 3 MVA ID", fontsize=20)
        plt.ylabel("Photon 4 MVA ID", fontsize=20)
        plt.title(f"{year} Events with Photon 3/4 MVA ID > {mva_cut}  (N = {len(mva3)})", fontsize=20)

        # Save figure
        path = path.replace(".pdf", "_m0p9cut.pdf")
        fig.savefig(path)
        assert os.path.isfile(path), "Saved plot does not exist."
        print(f"Saved {path.split('/')[-1]} to file.")
        """

        # Do correlation for pT a1/a2 and MVA ID
        # pt a1 vs gamma 1 MVA ID
        do_scatter(pt_a1, mva1, "Lead Pseudoscalar pT", "Photon 1 MVA ID", f"{year}(N = {len(mva1)}", 300, path.replace(".pdf", "_pta1_pho1.pdf"), idx+198)

        # pt a1 vs gamma 2 MVA ID
        do_scatter(pt_a1, mva2, "Lead Pseudoscalar pT", "Photon 2 MVA ID", f"{year}(N = {len(mva2)}", 300, path.replace(".pdf", "_pta1_pho2.pdf"), idx+199)

        # pt a1 vs gamma 3 MVA ID
        do_scatter(pt_a1, mva3, "Lead Pseudoscalar pT", "Photon 3 MVA ID", f"{year}(N = {len(mva3)}", 300, path.replace(".pdf", "_pta1_pho3.pdf"), idx+200)

        # pt a1 vs gamma 4 MVA ID
        do_scatter(pt_a1, mva4, "Lead Pseudoscalar pT", "Photon 4 MVA ID", f"{year}(N = {len(mva4)}", 300, path.replace(".pdf", "_pta1_pho4.pdf"), idx+201)
        
        # pt a2 vs gamma 1 MVA ID
        do_scatter(pt_a2, mva1, "Sublead Pseudoscalar pT", "Photon 1 MVA ID", f"{year}(N = {len(mva1)}", 200, path.replace(".pdf", "_pta2_pho1.pdf"), idx+195)

        # pt a2 vs gamma 2 MVA ID
        do_scatter(pt_a2, mva2, "Sublead Pseudoscalar pT", "Photon 2 MVA ID", f"{year}(N = {len(mva2)}", 200, path.replace(".pdf", "_pta2_pho2.pdf"), idx+196)

        # pt a2 vs gamma 3 MVA ID
        do_scatter(pt_a2, mva3, "Sublead Pseudoscalar pT", "Photon 3 MVA ID", f"{year}(N = {len(mva3)}", 200, path.replace(".pdf", "_pta2_pho3.pdf"), idx+202)

        # pt a1 vs gamma 4 MVA ID
        do_scatter(pt_a2, mva4, "Sublead Pseudoscalar pT", "Photon 4 MVA ID", f"{year}(N = {len(mva4)}", 200, path.replace(".pdf", "_pta2_pho4.pdf"), idx+203)


def plot_variations(
    nominal: pd.DataFrame,
    up: pd.DataFrame,
    down: pd.DataFrame,
    path: str,
    Scales: str,
) -> None:
    """
    Plot statistical variations relative to nominal for inputs and BDT scores.
    """

    branches, boundsList, binsScaleList, names, norms, titles = loadPlottingParameters()
    keep1 = 9
    keep2 = 11
    keep3 = 13
    branches = branches[:keep1] + branches[keep2:keep3]
    boundsList = boundsList[:keep1] + boundsList[keep2:keep3]
    binsScaleList = binsScaleList[:keep1] + binsScaleList[keep2:keep3]
    names = names[:keep1] + names[keep2:keep3]
    norms = norms[:keep1] + norms[keep2:keep3]
    titles = titles[:keep1] + titles[keep2:keep3]

    # Fill histograms for all BDT input variables and BDT score
    nhists = {}
    uhists = {}
    dhists = {}
    for branch, bounds, binsScale, name, norm, title in zip(branches, boundsList, binsScaleList, names, norms, titles):
        nhists.update({branch: loading.fillHist(branch, nominal, bounds, binsScale=binsScale, name=f"{branch}_nom", normalize=norm)})
        uhists.update({branch: loading.fillHist(branch, up, bounds, binsScale=binsScale, name=f"{branch}_up", normalize=norm)})
        dhists.update({branch: loading.fillHist(branch, down, bounds, binsScale=binsScale, name=f"{branch}_down", normalize=norm)})

    nhists.update({"BDT_score": loading.fillHist("BDT_score", nominal, [0.0, 1.0], binsScale=30.0, name="BDT_score_nom", normalize=False)})
    uhists.update({"BDT_score": loading.fillHist("BDT_score", up, [0.0, 1.0], binsScale=30.0, name="BDT_score_up", normalize=False)})
    dhists.update({"BDT_score": loading.fillHist("BDT_score", down, [0.0, 1.0], binsScale=30.0, name="BDT_score_down", normalize=False)})
    
    # Plot the nominal, up, and down cases together for each variable and the BDT score (3 hists per plot!)
    def plot(nhist, uhist, dhist, branch, title, path, Scales): 
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
        legend = TLegend(0.2,0.7,0.40,0.82)
        legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)

        # Main plotting
        gStyle.SetTitleFont(42)
        gStyle.SetLabelFont(42)

        nhist.SetTitle("")
        nhist.GetYaxis().SetTitleFont(42)
        nhist.GetYaxis().SetTitle("Events")
        nhist.GetXaxis().SetTitleOffset(1.1)
        nhist.Sumw2()
        nhist.SetLineWidth(2)
        
        uhist.SetTitle("")
        uhist.GetYaxis().SetTitleFont(42)
        uhist.GetYaxis().SetTitle("Events")
        uhist.GetXaxis().SetTitleOffset(1.1)
        uhist.Sumw2()
        uhist.SetLineWidth(2)
        
        dhist.SetTitle("")
        dhist.GetYaxis().SetTitleFont(42)
        dhist.GetYaxis().SetTitle("Events")
        dhist.GetXaxis().SetTitleOffset(1.1)
        dhist.Sumw2()
        dhist.SetLineWidth(2)
        
        nhist.GetYaxis().SetRangeUser(10**-3 + 10**-4, 10**7 - 5*10**6)
        uhist.GetYaxis().SetRangeUser(10**-3 + 10**-4, 10**7 - 5*10**6)
        dhist.GetYaxis().SetRangeUser(10**-3 + 10**-4, 10**7 - 5*10**6) 

        # Manually specify colors for plots
        nhist.SetLineColor(kGreen+3)
        uhist.SetLineColor(kRed)
        dhist.SetLineColor(kBlue)
        
        # Add objects to model1 legend
        #legend_model1.AddEntry(MakeNullPointer(TH1D), model1, "")
        legend.AddEntry(nhist, "Nominal", "l")
        legend.AddEntry(uhist, "Up Variation", "l")
        legend.AddEntry(dhist, "Down Variation", "l")

        nhist.Draw("hist same")
        uhist.Draw("hist same")
        dhist.Draw("hist same")
        legend.Draw()
        
        pad2.cd()

        # Residuals histogram setup
        ures = residuals_error_prop(uhist, nhist)
        dres = residuals_error_prop(dhist, nhist)

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
        for rs, color, shape in zip([ures, dres], [kRed, kBlue], [20, 21]):
            rs.SetTitle("")
            rs.GetXaxis().SetTitle(title)
            rs.GetYaxis().SetTitle("Var/Nominal")
            rs.GetXaxis().SetTitleSize(9.0*nhist.GetYaxis().GetTitleSize() * pad1Scale)
            rs.GetYaxis().SetTitleSize(8.0*nhist.GetYaxis().GetTitleSize() * pad1Scale)
            rs.GetXaxis().SetTitleOffset(1.1)
            rs.GetYaxis().SetTitleOffset(0.4)
            rs.GetYaxis().SetNdivisions(6)
            rs.GetYaxis().CenterTitle(True)
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
        max_edge = ures.GetBinCenter(ures.GetNbinsX())
        min_edge = ures.GetBinCenter(1)
        line = TLine(min_edge, 1.0, max_edge, 1.0)
        line.SetLineStyle(7)
        line.SetLineColor(kBlack)
        line.Draw("same")

        pad1.Modified()
        pad2.Modified()
        canvas.Update()

        plotFancy(pad1, title, lumiTxt=Scales.lumiTxt, pad=True, prelim=True, inPlot=False, Scales=Scales, shiftLeft=0.05)
        canvas.SaveAs(path)
        assert os.path.exists(path), "Saved plot does not exist."

    # Update the path with a filename for each plot - stored in each mass folder + variation folder
    branches.append("BDT_score")
    titles.append("BDT Score")
    for branch, title in zip(branches, titles):
        new_path = path + branch + ".pdf"
        plot(nhists[branch], uhists[branch], dhists[branch], branch, title, new_path, Scales) 


def plot_smoothed_hists(
    smoothed_hists: Dict,
    path: str,
    key: str,
    n: int,
    Scales: str
) -> None:
    """
    Plots smoothed histograms from categorization procedure.
    """
    
    # Canvas setup
    canvas = TCanvas("canvas_cats{n}", "canvas_cats{n}", 1000, 1000)
    canvas.SetLogy()
    gStyle.SetTitleFont(42)
    gStyle.SetLabelFont(42)
    maxY = 10**5 - 10**4
    minY = 10**-3 + 10**-4

    lsize = [0.65, 0.68, 0.92, 0.80]
    legend = TLegend(lsize[0], lsize[1], lsize[2], lsize[3])
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    for hk, hist in smoothed_hists.items():
        if "data_SR" in hk or "bkg_SB" in hk or "signal_SB" in hk:
            continue

        if "signal" in hk:
            hist.SetLineColor(kRed)
        elif "bkg" in hk:
            hist.SetLineColor(kAzure+3)
        elif "data" in hk:
            continue

        hist.SetTitle("")
        hist.GetYaxis().SetTitleFont(42)
        hist.GetYaxis().SetTitle("Events")
        hist.GetXaxis().SetTitle("BDT Score")
        hist.SetLineWidth(2)
        hist.GetYaxis().SetRangeUser(minY, maxY)
        legend.AddEntry(hist, hk, "l")
        hist.Draw("hist same")
    plotFancy(canvas, f"Smoothed BDT Histograms", lumiTxt=Scales.lumiTxt, prelim=True, inPlot=False, sub=f"Sideband only", Scales=Scales, shiftLeft=0.05)
    canvas.SaveAs(path)

    canvas.Destructor()

