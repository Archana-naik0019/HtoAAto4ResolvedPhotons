import os
import argparse
import pandas as pd
import numpy as np
import pyarrow
import time
from pathlib import Path
import pyarrow.parquet as pq
from pyarrow.lib import ArrowInvalid
from typing import List, Dict, Optional
from ROOT import TH1D, TFile, TCanvas, TLegend, kBlue, kRed, kBlack, kGreen, gPad, TLatex, TPaveText, gROOT
import h4g_tools.utils.plotting as plotting
import h4g_tools.utils.loading as loading

"""
To change sample over which to run, go to end of file and change the mc_path variable for uncorrected and/or corrected samples. 
To change the pileup variable to check, change [0] or [1] in pileup_var
"""


def loadDataPileup(
    f: TFile,
    key: str = "pileup"
) -> TH1D:
    """
    Load data pileup histogram.
    """

    # Use cd to ensure that we're looking at the correct file. ROOT assume that we're looking at the last opened file otherwise and can cause problems
    f.cd()

    # Can't use with-statement in this version of ROOT since TFile hasn't been made more Pythonic yet - update this if you want
    hist = TH1D(f.Get(key))
    assert type(hist) is TH1D, f"type(hist): {type(hist)}"

    return hist


def prepareMCHist(
    path: str,
    hist: TH1D,
    var: str,
    era: str
) -> dict[TH1D]:
    """
    Prepares MC pileup histograms (per mass point) using the data histogram as a template to keep binning and range consistent.
    """

    # Load samples per point
    mc_sample = loading.loadPerMass(path, branches=[var, "weight"])

    # Make empty/clean copy of data histogram
    clean = hist.Clone()
    clean.Reset()

    # Helper function to fill the histograms with the MC pileup information
    def fillHist(hist_to_fill: TH1D, sample: pd.DataFrame) -> TH1D:
        assert len(sample[var]) == len(sample["weight"])
        for val, weight in zip(sample[var].to_list(), sample["weight"].to_list()):
            hist_to_fill.Fill(val, weight)

        return hist_to_fill

    # Fill hists dictionary using same keys as mc_sample
    print(era)
    if era == "2022preEE":
        keys = [k for k in mc_sample.keys() if "preEE" in k]
    elif era == "2022postEE":
        keys = [k for k in mc_sample.keys() if "postEE" in k]
    elif era == "2024":
        keys = mc_sample.keys()
    else:
        print("Keys not specified!")
        sys.exit()
    hists = dict.fromkeys(keys)
    for key in keys:
        # Clone clean histogram each time to avoid redefining histograms
        """
        hists[mass] = clean.Clone()
        for key in mc_sample.keys():
            if mass in key:
                print(mass, key)
                hists[mass].Add(fillHist(hists[mass], mc_sample[key]))
        """
        hists[key] = clean.Clone()
        hists[key] = fillHist(hists[key], mc_sample[key])

        # Scale MC histogram to data for easier comparison of pileup
        hists[key].Scale(hist.Integral() / hists[key].Integral())

    return hists
    

def plotPileup(
    mass: str,
    mc_hist: TH1D,
    data_hist: TH1D,
    path: str,
    filename: str,
    corrected: bool,
    xs: str
) -> None:
    """
    Plot pileup of MC vs data using centrally produced histogram.
    """

    # Canvas setup
    canvas = TCanvas(f"{mass}_canvas", "{mass}_canvas", 1000, 1000)

    # Legend setup
    # Double the height when we go from 2 items in the legend to 4 to keep font consistently sized
    lsize = [0.60, 0.72 if type(data_hist) is not list else 0.64, 0.92, 0.8]
    legend = TLegend(lsize[0], lsize[1], lsize[2], lsize[3])
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    # Turn off stats box b/c it's not needed (I think)
    mc_hist.SetStats(0)
    if type(data_hist) is list:
        for dh in data_hist:
            dh.SetStats(0)
    elif type(data_hist) is TH1D:
        data_hist.SetStats(0)

    # Formatting for MC/data
    mc_hist.SetMarkerStyle(20)
    mc_hist.SetMarkerSize(0.6)
    mc_hist.SetLineWidth(2)
    mc_hist.SetTitle(f"")
    if type(data_hist) is list:
        for dh in data_hist:
            dh.SetMarkerStyle(20)
            dh.SetMarkerSize(0.6)
            dh.SetLineWidth(2)
            dh.SetTitle(f"")
    elif type(data_hist) is TH1D:
        data_hist.SetMarkerStyle(20)
        data_hist.SetMarkerSize(0.6)
        data_hist.SetLineWidth(2)
        data_hist.SetTitle(f"")

    # Y-axis formatting
    ytitle = "Events / 1.0"
    mc_hist.GetYaxis().SetTitle(ytitle)
    mc_hist.GetYaxis().SetTitleFont(40)
    if type(data_hist) is list:
        dh_max = 0
        for dh in data_hist:
            dh.GetYaxis().SetTitle(ytitle)
            dh.GetYaxis().SetTitleFont(40)
            if dh.GetBinContent(dh.GetMaximumBin()) > dh_max:
                dh_max = dh.GetBinContent(dh.GetMaximumBin())
        for dh in data_hist:
            dh.GetYaxis().SetRangeUser(0.0, dh.GetBinContent(dh.GetMaximumBin()) * 1.15)
    elif type(data_hist) is TH1D:
        data_hist.GetYaxis().SetTitle(ytitle)
        data_hist.GetYaxis().SetTitleFont(40)
        data_hist.GetYaxis().SetRangeUser(0.0, data_hist.GetBinContent(data_hist.GetMaximumBin()) * 1.15)

    # X-axis formatting
    xtitle = "Pileup" if not corrected else "Pileup (MC corrected)"
    mc_hist.GetXaxis().SetTitle(xtitle)
    mc_hist.GetXaxis().SetTitleFont(40)
    mc_hist.GetXaxis().SetTitleOffset(1.5)
    if type(data_hist) is list:
        for dh in data_hist:
            dh.GetXaxis().SetTitle(xtitle)
            dh.GetXaxis().SetTitleFont(40)
            dh.GetXaxis().SetTitleOffset(1.5)
    elif type(data_hist) is TH1D:
        data_hist.GetXaxis().SetTitle(xtitle)
        data_hist.GetXaxis().SetTitleFont(40)
        data_hist.GetXaxis().SetTitleOffset(1.5)

    # Use distinct colors for MC/data histograms for clarity
    # Note: for overlay, the order is central, up, down.
    mc_hist.SetMarkerColor(kGreen+3)
    mc_hist.SetLineColor(kGreen+3)
    if type(data_hist) is list:
        data_hist[0].SetMarkerColor(kBlack)
        data_hist[0].SetLineColor(kBlack)
        data_hist[1].SetMarkerColor(kBlue)
        data_hist[1].SetLineColor(kBlue)
        data_hist[2].SetMarkerColor(kRed)
        data_hist[2].SetLineColor(kRed)
    elif type(data_hist) is TH1D:
        data_hist.SetMarkerColor(kBlack)
        data_hist.SetLineColor(kBlack)

    # Add MC and data to legend
    # Note: for variations, the order is central, up, down.
    legend.AddEntry(mc_hist, f"MC (m_{{a}} = {mass.replace('Signal_','').replace('_',' ')})", "lep")
    if type(data_hist) is list:
        legend.AddEntry(data_hist[0], "Data", "lep")
        legend.AddEntry(data_hist[1], "Data (#sigma(pp) + 4.6%)", "lep")
        legend.AddEntry(data_hist[2], "Data (#sigma(pp) - 4.6%)", "lep")
    elif type(data_hist) is TH1D:
        legend.AddEntry(data_hist, "Data", "lep")

    # Draw with points and error bars b/c statistical error is important here!
    # Draw data first so it's underneath MC
    if type(data_hist) is list:
        # Draw variations first so they're underneath the central xs pileup
        data_hist[2].Draw("P E1 same")
        data_hist[1].Draw("P E1 same")
        data_hist[0].Draw("P E1 same")
    elif type(data_hist) is TH1D:
        data_hist.Draw("P E1")
    mc_hist.Draw("P E1 same")
    legend.Draw("same")

    sub = "         MC normalized to Data" if type(data_hist) is not list else "MC normalized to Central Data"
    lumiTxt = f" #sigma(pp) = {xs} mb^{{-1}}  (13.6 TeV)".strip() if type(data_hist) is not list else "(13.6 TeV)"
    plotting.plotFancy(canvas, "", lumiTxt=lumiTxt, prelim=True, inPlot=True, sub=sub, overlay=True)

    # Update filename in case of corrected weights
    if corrected:
        filename = filename.replace(".pdf", "_corrected.pdf")

    # Ensure directories exist, save plot, and ensure saving worked
    if not os.path.exists(path):
        Path(path).mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(f"{path}/{filename}")
    assert os.path.isfile(f"{path}/{filename}"), "Saved plot does not exist."


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--year", required=True, help="Specify year/era of pileup Ntuples to use.")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("-up", "--up", default=False, required=False, action="store_true", help="Use +4% variation from central inelastic pp xs recommended by LUM.")
    group.add_argument("-down", "--down", default=False, required=False, action="store_true", help="Use -4% variation from central inelastic pp xs recommended by LUM.")
    group.add_argument("-overlay", "--overlay", default=False, required=False, action="store_true", help="Overlay central and variation pp xs histograms for data.")
    args = parser.parse_args()

    if "2022" not in args.year and "2024" not in args.year:
        print("Must provide valid year - based on which pileup files I have.")
        sys.exit()

    # Set batch so no plots pop up
    gROOT.SetBatch(True)

    # Get directory of this file
    cwd = os.path.dirname(os.path.abspath(__file__))

    eras = []
    if args.year == "2022":
        eras = ["2022postEE", "2022preEE"]
    elif args.year == "2024":
        eras = ["2024"]
    # ADD DETAILS FOR NEW YEAR HERE

    for era in eras:
        # Central (recommended by LUM) inelastic pp xs file
        if "2022" in args.year:
            if era == "2022preEE":
                pileup_file_central = "pileupHistogram-Cert_Collisions2022_355100_357900_eraBCD_GoldenJson-13p6TeV-69200ub-99bins.root"
                # pp xs +4%
                pileup_file_up = "pileupHistogram-Cert_Collisions2022_355100_357900_eraBCD_GoldenJson-13p6TeV-72400ub-99bins.root"
                # pp xs -4%
                pileup_file_down = "pileupHistogram-Cert_Collisions2022_355100_357900_eraBCD_GoldenJson-13p6TeV-66000ub-99bins.root"
            elif era == "2022postEE":
                pileup_file_central = "pileupHistogram-Cert_Collisions2022_359022_362760_eraEFG_GoldenJson-13p6TeV-69200ub-99bins.root"
                # pp xs +4%
                pileup_file_up = "pileupHistogram-Cert_Collisions2022_359022_362760_eraEFG_GoldenJson-13p6TeV-72400ub-99bins.root"
                # pp xs -4%
                pileup_file_down = "pileupHistogram-Cert_Collisions2022_359022_362760_eraEFG_GoldenJson-13p6TeV-66000ub-99bins.root"
        elif "2024" in args.year:
            pileup_file_central = "dataPileupHistogram-2024BCDEFGHI-69200ub.root"
            # pp xs +4%
            pileup_file_up = "dataPileupHistogram-2024BCDEFGHI-72400ub.root"
            # pp xs -4%
            pileup_file_down = "dataPileupHistogram-2024BCDEFGHI-66000ub.root"

        # Note: if I don't have the pileup histograms, look here for the relevant CollisionXX/PileUp/ directory:
        # https://cms-service-dqmdc.web.cern.ch/CAF/certification/
        # Additional information on these central histograms here:
        # https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData#Centrally_produced_ROOT_histogra

        # Choose pileup file to use for data
        xs = None
        pileup_file = None
        xs_split = 2 if "2022" in args.year else 1
        if not args.up and not args.down and not args.overlay:
            pileup_file = os.path.join("pileup", pileup_file_central)
            xs = pileup_file.split("-")[-xs_split].replace("ub","").replace(".root","")
            xs = str(float(xs) / 1000)
            print(f"\nUsing central {xs} ub pp xs file!\n")
        elif args.up:
            pileup_file = os.path.join("pileup", pileup_file_up)
            xs = pileup_file.split("-")[-1*xs_split].replace("ub","").replace(".root","")
            print("\n~~~~~~~~~~~~~~~~~")
            print(f"WARNING: Using +4% variation, {xs} ub pp xs file!")
            print("~~~~~~~~~~~~~~~~~\n")
            time.sleep(2)
        elif args.down:
            pileup_file = os.path.join("pileup", pileup_file_down)
            xs = pileup_file.split("-")[-xs_split].replace("ub","").replace(".root","")
            print("\n~~~~~~~~~~~~~~~~~")
            print(f"WARNING: Using -4% variation, {xs} ub pp xs file!")
            print("~~~~~~~~~~~~~~~~~\n")
            time.sleep(2)
        elif args.overlay:
            pileup_file = {"central": os.path.join("pileup", pileup_file_central), "up": os.path.join("pileup", pileup_file_up), "down": os.path.join("pileup", pileup_file_down)}

        assert pileup_file is not None

        # Do a very un-Pythonic thing and define TFile outside of function so the histogram stays in scope
        # Also, only define this once since it's a common object between uncorr/corr MC samples.
        f = TFile()
        fup = TFile()
        fdown = TFile()
        if not args.overlay:
            f = TFile(os.path.join(cwd, pileup_file))
            dataPU = loadDataPileup(f)

            # Close fup/fdown since we don't need them here
            fup.Close()
            fdown.Close()
        elif args.overlay:
            f = TFile(os.path.join(cwd, pileup_file["central"]))
            fup = TFile(os.path.join(cwd, pileup_file["up"]))
            fdown = TFile(os.path.join(cwd, pileup_file["down"]))
            dataPU_central = loadDataPileup(f)
            dataPU_up = loadDataPileup(fup)
            dataPU_down = loadDataPileup(fdown)

            # Rename the variations for debugging purposes
            dataPU_central.SetName("pileup_central")
            dataPU_up.SetName("pileup_up")
            dataPU_down.SetName("pileup_down")


        # NOTE: EDIT ME TO CHANGE INPUT SAMPLES!
        if "2022" in args.year:
            uncorr_path = "/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2022_noCorr_noSel_26Jan2026"
            corr_path = "/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2022_PUonly_noSel_26Jan2026"
        elif "2024" in args.year:
            uncorr_path = "/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_noCorr_noSel_26Jan2026"
            corr_path = "/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_PUonly_noSel_26Jan2026"
        # ADD SAMPLES FOR NEW YEAR HERE!
        for corr, mc_path in {"uncorr": uncorr_path, "corr": corr_path}.items():
            print()
            # Specify pileup variable to check - nTrueInt is the one modified with pileup reweighting
            pileup_var = ["nPU", "nTrueInt"]
            # An assumption is made that the binning is the same for central and variation histograms here!
            mcPU = prepareMCHist(mc_path, dataPU if not args.overlay else dataPU_central, pileup_var[1], era)
            assert type(mcPU) is dict

            # Adjust naming of output directories so outputs are easier to read
            uncorrected = "uncorrected"
            corrected = "corrected"
            if args.up:
                uncorrected += "_up"
                corrected += "_up"
            elif args.down:
                uncorrected += "_down"
                corrected += "_down"
            elif args.overlay:
                uncorrected += "_overlay"
                corrected += "_overlay"

            # Plot per mass so we can see how/if it changes per MC sample
            for mass in mcPU.keys():
                plotPileup(
                    mass,
                    mcPU[mass],
                    dataPU if not args.overlay else [dataPU_central, dataPU_up, dataPU_down],
                    str(os.path.join(cwd, f"pileup_checks_{args.year}", uncorrected if corr == "uncorr" else corrected)),
                    f"{mass}_pileup.pdf",
                    False if corr == "uncorr" else True,
                    xs
                )

                fhist_out = TFile(os.path.join(cwd, "pileupHists.root"), "RECREATE")
                fhist_out.cd()
                d = fhist_out.mkdir("pileup")
                d.cd()
                mcPU[mass].SetName("pileup_MC_Signal_30_GeV")
                mcPU[mass].Write()
                dataPU.SetName("pileup_data")
                dataPU.Write()
                fhist_out.Close()
                f.cd()

        # Close TFile to conclude scope problem from before
        f.Close()
        if args.overlay:
            fup.Close()
            fdown.Close()
