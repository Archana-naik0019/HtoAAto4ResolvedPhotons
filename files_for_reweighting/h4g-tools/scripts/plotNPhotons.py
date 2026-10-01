import os
import pandas as pd
import numpy as np
from typing import List, Dict, Optional
import argparse
from h4g_tools.utils.scales import Scales as ScalesCls
from ROOT import kBlack, kRed, kGreen, kOrange, kBlue, TH1D, TLegend, TCanvas
from h4g_tools.utils.loading import fillHist, loadPerMass
from h4g_tools.utils.plotting import plotFancy


def plotNPhotons(
    samples: Dict[str, pd.DataFrame],
    filename: str,
    era: str,
    path: str,
    scales: ScalesCls,
    region: str,
    log: bool = False
) -> None:
    """
    Plots number of photons at different cut stages.
    """

    # Region setup
    region_str = ""
    if "EB" in region and "One" not in region:
        region_str = "EB"
    elif "EE" in region:
        region_str = "EE"
    if "OneEB" in region:
        region_str = "OneEB"

    # Fill scatter arrays
    if region_str == "":
        Nphotons = {
            "Nphotons Initial Events": 0,
            "Nphotons Pre-selections": 0,
            "Nphotons Pre-selections + At least 4 photons": 0,
            "Nphotons Pre-selections + At least 4 photons + Pseudoscalar selections": 0
        }
    elif region_str == "EB":
        Nphotons = {
            "Nphotons_EB Initial Events": 0,
            "Nphotons_EB Pre-selections": 0,
            "Nphotons_EB Pre-selections + At least 4 photons": 0,
            "Nphotons_EB Pre-selections + At least 4 photons + Pseudoscalar selections": 0
        }
    elif region_str == "EE":
        Nphotons = {
            "Nphotons_EE Initial Events": 0,
            "Nphotons_EE Pre-selections": 0,
            "Nphotons_EE Pre-selections + At least 4 photons": 0,
            "Nphotons_EE Pre-selections + At least 4 photons + Pseudoscalar selections": 0
        }
    elif region_str == "OneEB":
        Nphotons = {
            "Nphotons_OneEB Initial Events": 0,
            "Nphotons_OneEB Pre-selections": 0,
            "Nphotons_OneEB Pre-selections + At least 4 photons": 0,
            "Nphotons_OneEB Pre-selections + At least 4 photons + Pseudoscalar selections": 0
        }

    masses = []
    photons = {}

    if region_str != "":
        print(f"\n~~ {era} {region_str} ~~")
    else:
        print(f"\n~~ {era} ~~")

    sorted_samples = sorted(samples.items())
    for idx, (mass, sample) in enumerate(sorted_samples):
        mass_in = int(mass[7:9])
        photons.update({mass_in: Nphotons.copy()})
        print(f"Mass: {mass_in} GeV")

        maximum = 0.0
        for i, key in enumerate(Nphotons.keys()):
            Nphotons[key] = {}
            # Need to translate tuple of (N photons, counts) into lists to fill hists
            # The metadata keys don't match Nphotons keys - also translate from binary string to integers
            for metadata_key in sample.schema.metadata:
                if key == metadata_key.decode("utf-8")[:-2].strip():
                    Nphotons[key].update({int(metadata_key.decode("utf-8")[-2:]): int(float(sample.schema.metadata[metadata_key].decode("utf-8")))})

            bins = list(Nphotons[key].keys())
            weights = list(Nphotons[key].values())
            #print(bins)
            #print(weights)
            photons[mass_in][key] = fillHist(f"Nphotons_{region}_{i}", bins, customWeights=weights, bounds=[-0.5,9.5], normalize=False)

            # Set bin errors to zero because it's hard to read otherwise. Also don't know how to turn them off...
            for b in range(photons[mass_in][key].GetNbinsX()):
                photons[mass_in][key].SetBinError(b, 0.0)

            if maximum < photons[mass_in][key].GetBinContent(photons[mass_in][key].GetMaximumBin()):
                maximum = photons[mass_in][key].GetBinContent(photons[mass_in][key].GetMaximumBin())

        # Canvas setup
        canvas = TCanvas("canvas", "canvas", 1000, 1000)

        # Set log scale if required
        if log:
            canvas.SetLogy()

        # Legend setup
        lsize = [0.57, 0.65, 0.92, 0.85]
        legend = TLegend(lsize[0], lsize[1], lsize[2], lsize[3])
        legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)

        for idx, draw in enumerate(["hist", "p"]):
            #for key, color, marker in zip(list(Nphotons.keys()), [kBlack, kRed, kBlue, kOrange+1, kGreen+3], [5, 25, 46, 24, 26, 32]):
            for key, color, marker in zip(list(Nphotons.keys()), [kBlack, kBlue, kOrange+1, kGreen+3], [5, 46, 24, 26, 32]):
                if idx == 0:
                    legend_key = " ".join(key.split(" ")[1:])
                    if "Pseudo" in legend_key:
                        legend.AddEntry(photons[mass_in][key], f"#splitline{{{legend_key[:legend_key.find('+ Ps')]}}}{{{legend_key[legend_key.find('+ Ps'):]}}}", "lp")
                    else:
                        legend.AddEntry(photons[mass_in][key], legend_key, "lp")
                    
                photons[mass_in][key].SetMarkerStyle(marker)
                photons[mass_in][key].SetMarkerColor(color)
                photons[mass_in][key].SetMarkerSize(2)
                photons[mass_in][key].SetLineWidth(2)
                photons[mass_in][key].SetTitle(f"")
                photons[mass_in][key].GetYaxis().SetTitleFont(42)
                photons[mass_in][key].GetYaxis().SetTitle("Events")
                photons[mass_in][key].GetXaxis().SetTitleFont(42)
                photons[mass_in][key].GetXaxis().SetTitle("# of Photons")
                photons[mass_in][key].SetLineColor(color)
                photons[mass_in][key].GetYaxis().SetRangeUser(0.0 if not log else 0.01, maximum * (1.2 if not log else 10**2))
                photons[mass_in][key].Draw(f"same {draw}")

        legend.Draw("same")
        plotFancy(canvas, f"m_{{a}} = {mass_in} GeV" if region_str == "" else f"m_{{a}} = {mass_in} GeV in {region_str}", lumiTxt=Scales.lumiTxt.strip(), prelim=True, inPlot=True)

        new_filename = filename.replace("_nPhotons", f"{mass_in}GeV_{era}_nPhotons{region_str}")
        if not os.path.exists(f"{path}/{mass_in}_GeV/"):
            os.mkdir(f"{path}/{mass_in}_GeV/")
        canvas.SaveAs(f"{path}/{mass_in}_GeV/{new_filename}")
        assert os.path.isfile(f"{path}/{mass_in}_GeV/{new_filename}"), "Saved plot does not exist."
        print(f"Saved {path}/{mass_in}_GeV/{new_filename.split('/')[-1]} to file.")



def plotNPhotons_individualPlots(
    samples: Dict[str, pd.DataFrame],
    filename: str,
    era: str,
    path: str,
    scales: ScalesCls,
    region: str,
    log: bool = False
) -> None:
    """
    Plots number of photons at different cut stages as individual plots.
    """
    
    # Region setup
    region_str = ""
    if "EB" in region and "One" not in region:
        region_str = "EB"
    elif "EE" in region:
        region_str = "EE"
    if "OneEB" in region:
        region_str = "OneEB"

    # Fill scatter arrays
    if region_str == "":
        Nphotons = {
            "Nphotons Initial Events": 0,
            "Nphotons Pre-selections": 0,
            "Nphotons Pre-selections + At least 4 photons": 0,
            "Nphotons Pre-selections + At least 4 photons + Pseudoscalar selections": 0
        }
    elif region_str == "EB":
        Nphotons = {
            "Nphotons_EB Initial Events": 0,
            "Nphotons_EB Pre-selections": 0,
            "Nphotons_EB Pre-selections + At least 4 photons": 0,
            "Nphotons_EB Pre-selections + At least 4 photons + Pseudoscalar selections": 0
        }
    elif region_str == "EE":
        Nphotons = {
            "Nphotons_EE Initial Events": 0,
            "Nphotons_EE Pre-selections": 0,
            "Nphotons_EE Pre-selections + At least 4 photons": 0,
            "Nphotons_EE Pre-selections + At least 4 photons + Pseudoscalar selections": 0
        }
    elif region_str == "OneEB":
        Nphotons = {
            "Nphotons_OneEB Initial Events": 0,
            "Nphotons_OneEB Pre-selections": 0,
            "Nphotons_OneEB Pre-selections + At least 4 photons": 0,
            "Nphotons_OneEB Pre-selections + At least 4 photons + Pseudoscalar selections": 0
        }

    masses = []
    photons = {}

    if region_str != "":
        print(f"\n~~ {era} {region_str} ~~")
    else:
        print(f"\n~~ {era} ~~")

    sorted_samples = sorted(samples.items())
    for idx, (mass, sample) in enumerate(sorted_samples):
        mass_in = int(mass[7:9])
        photons.update({mass_in: Nphotons.copy()})
        print(f"Mass: {mass_in} GeV")

        if not os.path.exists(f"{path}/{mass_in}_GeV/"):
            os.mkdir(f"{path}/{mass_in}_GeV/")

        maximum = 0.0
        for i, key in enumerate(Nphotons.keys()):
            Nphotons[key] = {}
            # Need to translate tuple of (N photons, counts) into lists to fill hists
            # The metadata keys don't match Nphotons keys - also translate from binary string to integers
            for metadata_key in sample.schema.metadata:
                if key == metadata_key.decode("utf-8")[:-2].strip():
                    Nphotons[key].update({int(metadata_key.decode("utf-8")[-2:]): int(float(sample.schema.metadata[metadata_key].decode("utf-8")))})

            bins = list(Nphotons[key].keys())
            weights = list(Nphotons[key].values())
            photons[mass_in][key] = fillHist(region, bins, customWeights=weights, bounds=[-0.5,9.5], normalize=False)

            # Set bin errors to zero because it's hard to read otherwise. Also don't know how to turn them off...
            for b in range(photons[mass_in][key].GetNbinsX()):
                photons[mass_in][key].SetBinError(b, 0.0)

            if maximum < photons[mass_in][key].GetBinContent(photons[mass_in][key].GetMaximumBin()):
                maximum = photons[mass_in][key].GetBinContent(photons[mass_in][key].GetMaximumBin())

        #for sel_index, (key, color, marker) in enumerate(zip(list(Nphotons.keys()), [kBlack, kRed, kBlue, kOrange+1, kGreen+3], [5, 25, 46, 24, 26, 32])):
        for sel_index, (key, color, marker) in enumerate(zip(list(Nphotons.keys()), [kBlack, kBlue, kOrange+1, kGreen+3], [5, 46, 24, 26, 32])):
            # Canvas setup
            canvas = TCanvas("canvas", "canvas", 1000, 1000)

            # Set log scale if required
            if log:
                canvas.SetLogy()

            # Legend setup
            lsize = [0.57, 0.65, 0.92, 0.85]
            legend = TLegend(lsize[0], lsize[1], lsize[2], lsize[3])
            legend.SetTextFont(42)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)

            for idx, draw in enumerate(["hist", "p"]):
                if idx == 0:
                    legend_key = " ".join(key.split(" ")[1:])
                    if "Pseudo" in legend_key:
                        legend.AddEntry(photons[mass_in][key], f"#splitline{{{legend_key[:legend_key.find('+ Ps')]}}}{{{legend_key[legend_key.find('+ Ps'):]}}}", "lp")
                    else:
                        legend.AddEntry(photons[mass_in][key], legend_key, "lp")
                    
                photons[mass_in][key].SetMarkerStyle(marker)
                photons[mass_in][key].SetMarkerColor(color)
                photons[mass_in][key].SetMarkerSize(2)
                photons[mass_in][key].SetLineWidth(2)
                photons[mass_in][key].SetTitle(f"")
                photons[mass_in][key].GetYaxis().SetTitleFont(42)
                photons[mass_in][key].GetYaxis().SetTitle("Events")
                photons[mass_in][key].GetXaxis().SetTitleFont(42)
                photons[mass_in][key].GetXaxis().SetTitle("# of Photons")
                photons[mass_in][key].SetLineColor(color)
                photons[mass_in][key].GetYaxis().SetRangeUser(0.0 if not log else 0.01, maximum * (1.2 if not log else 10.0))
                photons[mass_in][key].Draw(f"same {draw}")

            legend.Draw("same")
            plotFancy(canvas, f"m_{{a}} = {mass_in} GeV" if region_str == "" else f"m_{{a}} = {mass_in} GeV in {region_str}", lumiTxt=Scales.lumiTxt.strip(), prelim=True, inPlot=True, sub=f"Integral: {photons[mass_in][key].Integral()}")

            new_filename = filename.replace("_nPhotons", f"_{mass_in}GeV_{era}_nPhotons{sel_index}{region_str}")
            if not os.path.exists(f"{path}/{mass_in}_GeV/split_{region}/"):
                os.mkdir(f"{path}/{mass_in}_GeV/split_{region}")
            canvas.SaveAs(f"{path}/{mass_in}_GeV/split_{region}/{new_filename}")
            assert os.path.isfile(f"{path}/{mass_in}_GeV/split_{region}/{new_filename}"), "Saved plot does not exist."
            print(f"Saved {path}/{mass_in}_GeV/split_{region}/{new_filename.split('/')[-1]} to file.")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", required=False, type=str, help="Only run photon plots for this subdirectory.")
    parser.add_argument("-l", "--log", required=False, action="store_true", help="Plot with log scale on y-axis.")
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    args = parser.parse_args()
 
    # Do year/eras setup
    Scales = ScalesCls(args.year) 
    masses = [x for x in range(15, 65, 5)]
    if args.year == "2022preEE":
        eras={
            "data": ["Run2022C","Run2022D"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_preEE" for m in masses]
        }
    elif args.year == "2022postEE":
        eras={
            "data": ["Run2022E","Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_postEE" for m in masses]
        }
    elif args.year == "2022":
        eras={
            "data": ["Run2022C","Run2022D", "Run2022E","Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D", "RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preEE", "postEE"]]
        }
    else:
        assert False, "Must specify year in Scales class!"


    # Get specified directory in signal_tests
    d = "/afs/crc.nd.edu/user/s/scastel2/Private/higgs-dna-4-gamma-tanays-copy/scripts/signal_tests/"
    lst_dir = os.listdir(d)
    assert args.directory in lst_dir if args.directory is not None else True, f"Provided directory not found in {lst_dir}!"

    for directory in lst_dir:
        if args.directory is not None:
            if directory != args.directory.strip("/"):
                continue

        if directory != "backup" and directory != "photon_plots" and directory != "photon_plots_log":
            if not os.path.exists(os.path.join(d, "photon_plots" if not args.log else "photon_plots_log")):
                os.mkdir(os.path.join(d, "photon_plots" if not args.log else "photon_plots_log"))

            samples = loadPerMass(f"{d}/{directory}", eras=eras, meta=True)
            for region in ["Nphotons", "Nphotons_EB", "Nphotons_EE", "Nphotons_OneEB"]:
                plotNPhotons(samples, f"{directory}_nPhotons.png", args.year, os.path.join(d, "photon_plots" if not args.log else "photon_plots_log"), Scales, region, log=args.log)
                plotNPhotons_individualPlots(samples, f"{directory}_nPhotons.png", args.year, os.path.join(d, "photon_plots" if not args.log else "photon_plots_log"), Scales, region, log=args.log)
