import os
from pathlib import Path
import numpy as np
import pandas as pd
from typing import List, Dict, Optional
import argparse
from h4g_tools.utils.scales import Scales as ScalesCls
from ROOT import kBlack, kRed, kGreen, kOrange, kBlue, TCanvas, TH2D, TLegend
from h4g_tools.utils.loading import fillHist, loadPerMass
from h4g_tools.utils.plotting import plotFancy


def plot2D(
    samples_uncorr: Dict[str, pd.DataFrame],
    samples_corr: Dict[str, pd.DataFrame],
    filename: str,
    era: str,
    path: str,
    scales: ScalesCls,
    var: str
) -> None:
    """
    Plot 2D histogram of var comparing corr/uncorr samples. Also plot ratio of corr/uncorr in 1D and corr values for any zero uncorr values to see movement.
    """ 

    sorted_samples_uncorr = sorted(samples_uncorr.items())
    sorted_samples_corr = sorted(samples_corr.items())
    for idx, ((mass_uncorr, sample_uncorr), (mass_corr, sample_corr)) in enumerate(zip(sorted_samples_uncorr, sorted_samples_corr)):
        assert mass_uncorr == mass_corr
        mass = mass_uncorr
        mass_in = int(mass[7:9])
        print(f"Mass: {mass_in} GeV")

        # Canvas setup
        canvas = TCanvas("canvas", "canvas", 1000, 1000)

        # Histogram with custom number of equal bins, bounded by (bounds[0], bounds[1])
        normalize = False

        if var == "r9":
            names = [var] * 4
            branches = [f"pho{N+1}_{var}" for N in range(4)]
            # Binning will be added to titles later
            ytitles = [f"R_{{9}} #gamma_{{{N+1}}} (Uncorrected)" for N in range(4)]
            xtitles = [f"R_{{9}} #gamma_{{{N+1}}} (Corrected)" for N in range(4)]
            ratio_titles = [f"R_{{9}} #gamma_{{{N+1}}} (Corr/Uncorr)" for N in range(4)]
            boundsList = [[0.0, 2.0, 0.0, 2.0, 0.0, 2.0]] * 4
            binsScales = [[10.0, 10.0, 10.0]] * 4
        elif var == "hoe":
            names = [var] * 4
            branches = [f"pho{N+1}_{var}" for N in range(4)]
            # Binning will be added to titles later
            ytitles = [f"H/E #gamma_{{{N+1}}} (Uncorrected)" for N in range(4)]
            xtitles = [f"H/E #gamma_{{{N+1}}} (Corrected)" for N in range(4)]
            ratio_titles = [f"H/E #gamma_{{{N+1}}} (Corr/Uncorr)" for N in range(4)]
            boundsList = [[0.0, 1.0, 0.0, 1.0, 0.0, 3.0]] * 4
            binsScales = [[100.0, 100.0, 100.0]] * 4
        elif var == "sieie":
            names = [var] * 4
            branches = [f"pho{N+1}_{var}" for N in range(4)]
            # Binning will be added to titles later
            ytitles = [f"#sigma_{{i#etai#eta}} #gamma_{{{N+1}}} (Uncorrected)" for N in range(4)]
            xtitles = [f"#sigma_{{i#etai#eta}} #gamma_{{{N+1}}} (Corrected)" for N in range(4)]
            ratio_titles = [f"#sigma_{{i#etai#eta}} #gamma_{{{N+1}}} (Corr/Uncorr)" for N in range(4)]
            boundsList = [[0.0, 0.1, 0.0, 0.1, 0.0, 3.0]] * 4
            binsScales = [[1000.0, 1000.0, 1000.0]] * 4
        elif var == "mvaID":
            names = [var] * 4
            branches = [f"pho{N+1}_{var}" for N in range(4)]
            # Binning will be added to titles later
            ytitles = [f"MVA ID #gamma_{{{N+1}}} (Uncorrected)" for N in range(4)]
            xtitles = [f"MVA ID #gamma_{{{N+1}}} (Corrected)" for N in range(4)]
            ratio_titles = [f"MVA ID #gamma_{{{N+1}}} (Corr/Uncorr)" for N in range(4)]
            boundsList = [[-1.0, -1.0, -1.0, 1.0, -1.0, 1.0]] * 4
            binsScales = [[100.0, 100.0, 100.0]] * 4

        assert len(names) == len(branches) == len(xtitles) == len(ytitles) == len(ratio_titles) == len(boundsList) == len(binsScales), f"{len(names)} {len(branches)} {len(xtitles)} {len(ytitles)} {len(ratio_titles)} {len(boundsList)} {len(binsScales)}"
        for N, (name, branch, xtitle, ytitle, ratio_title, bounds, binsScale) in enumerate(zip(names, branches, xtitles, ytitles, ratio_titles, boundsList, binsScales)):
            ### Do 2D hist comparison
            canvas.Clear()

            # Fill 2D hist with corr/uncorr samples
            hist = TH2D(name, name, int((bounds[1]-bounds[0])*binsScale[0]), bounds[0], bounds[1], int((bounds[3]-bounds[2])*binsScale[1]), bounds[2], bounds[3])
            ybranch_df = sample_uncorr[branch].tolist()
            xbranch_df = sample_corr[branch].tolist()

            # Fill histogram
            for x, y in zip(xbranch_df, ybranch_df):
                # Omit weights since they should not be changed by NFs (i.e. the only correction applied)
                hist.Fill(x,y)

            # Normalize
            if normalize:
                if hist.Integral() > 0.0:
                    hist.Scale(1.0 / hist.Integral())
                else:
                    hist.Scale(1.0)
                    print("Cannot normalize empty histogram.")

            hist.SetStats(True)  # Want to see stats box for entries, mean x/y, std dev x/y, and integral
            hist.SetMarkerStyle(8)
            hist.SetMarkerSize(0.5)
            hist.SetTitle(f"")
            hist.GetYaxis().SetTitleFont(42)
            hist.GetYaxis().SetTitle(ytitle + f"  /  {1.0 / binsScale[1]: 1.3f}")
            hist.GetXaxis().SetTitleFont(42)
            hist.GetXaxis().SetTitle(xtitle + f"  /  {1.0 / binsScale[1]: 1.3f}")
            hist.GetXaxis().SetTitleOffset(1.5)
            hist.Draw("COLZ")

            gPad.Update()
            st = TPaveStats(hist.FindObject("stats"))
            #st.SetX1NDC(newx1)
            st.SetX2NDC(0.5)

            plotFancy(canvas, "", lumiTxt=Scales.lumiTxt.strip(), prelim=True, inPlot=False, colz=True)

            new_filename = filename.replace("X", str(N+1))
            directory = os.path.join(path, f"{mass_in}_GeV", "2D_compare")
            if not os.path.exists(directory):
                Path(directory).mkdir(parents=True, exist_ok=True)

            full_filepath = os.path.join(directory, new_filename)
            canvas.SaveAs(full_filepath)
            assert os.path.isfile(full_filepath), "Saved plot does not exist."

            ### Organize samples for zero values of uncorr
            zero_check = 0.00001
            zeros_uncorr = sample_uncorr[sample_uncorr[f"pho{N+1}_{var}"] <= zero_check]
            zeros_corr = sample_corr[sample_uncorr[f"pho{N+1}_{var}"] <= zero_check]
            nonzero_uncorr = sample_uncorr[sample_uncorr[f"pho{N+1}_{var}"] > zero_check]
            nonzero_corr = sample_corr[sample_uncorr[f"pho{N+1}_{var}"] > zero_check]

            ### Now do 1D ratio of corr/uncorr samples
            canvas.Clear()

            # Fill 1D hist ratio of corrected / uncorrected samples
            hist_ratio = fillHist(branch, nonzero_corr / nonzero_uncorr, bounds=[bounds[4], bounds[5]], binsScale=binsScale[2], normalize=False)

            #hist_ratio.SetStats(0)  # Keep on stats box to compare with 2D plots
            hist_ratio.SetMarkerStyle(8)
            hist_ratio.SetMarkerSize(0.5)
            hist_ratio.SetTitle(f"")
            hist_ratio.GetYaxis().SetTitleFont(42)
            hist_ratio.GetYaxis().SetTitle(f"Events / {1.0 / binsScale[2]: 1.3f}")
            hist_ratio.GetXaxis().SetTitleFont(42)
            hist_ratio.GetXaxis().SetTitle(ratio_title)
            hist_ratio.GetXaxis().SetTitleOffset(1.5)
            hist_ratio.Draw("hist")

            plotFancy(canvas, "", lumiTxt=Scales.lumiTxt.strip(), prelim=True, inPlot=False)

            new_filename = filename.replace("X", str(N+1))
            directory = os.path.join(path, f"{mass_in}_GeV", "ratio")
            if not os.path.exists(directory):
                Path(directory).mkdir(parents=True, exist_ok=True)

            full_filepath = os.path.join(directory, new_filename)
            canvas.SaveAs(full_filepath)
            assert os.path.isfile(full_filepath), "Saved plot does not exist."

            ### Now do corrected value for which there is a corresponding zero in uncorrected samples
            print(f"# of zeros in uncorr: {len(zeros_uncorr)}")
            if len(zeros_uncorr) == 0:
                print("Skipping zeros plot...")
            elif len(zeros_uncorr) > 0:
                canvas.Clear()

                # Fill 1D hist ratio of corrected / uncorrected samples
                hist_zeros = fillHist(branch, zeros_corr, bounds=[bounds[2], bounds[3]], binsScale=binsScale[1], normalize=False)

                #hist_zeros.SetStats(0)  # Keep on stats box to compare with 2D plots
                hist_zeros.SetMarkerStyle(8)
                hist_zeros.SetMarkerSize(0.5)
                hist_zeros.SetTitle(f"")
                hist_zeros.GetYaxis().SetTitleFont(42)
                hist_zeros.GetYaxis().SetTitle(f"Events / {1.0 / binsScale[1]: 1.3f}")
                hist_zeros.GetXaxis().SetTitleFont(42)
                hist_zeros.GetXaxis().SetTitle(f"{xtitle} (Uncorr = 0)")
                hist_zeros.GetXaxis().SetTitleOffset(1.5)
                hist_zeros.Draw("hist")

                plotFancy(canvas, "", lumiTxt=Scales.lumiTxt.strip(), prelim=True, inPlot=False)

                new_filename = filename.replace("X", str(N+1))
                directory = os.path.join(path, f"{mass_in}_GeV", "zeros")
                if not os.path.exists(directory):
                    Path(directory).mkdir(parents=True, exist_ok=True)

                full_filepath = os.path.join(directory, new_filename)
                canvas.SaveAs(full_filepath)
                assert os.path.isfile(full_filepath), "Saved plot does not exist."


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    parser.add_argument("-su", "--sigUC", help="Directory to import uncorrected signal samples from.", type=str, default=None, required=False)
    parser.add_argument("-suc", "--sigC", help="Directory to import corrected signal samples from.", type=str, default=None, required=False)
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


    # Get specified directory
    d = "/afs/crc.nd.edu/user/s/scastel2/Private/higgs-dna-4-gamma-tanays-copy/scripts"
    assert os.path.exists(os.path.join(d, args.sigUC))
    assert os.path.exists(os.path.join(d, args.sigC))

    variables = [
        "r9",
        "hoe",
        "sieie",
        "mvaID"
    ]

    # Get uncorrected samples
    samples_uncorr = loadPerMass(os.path.join(d, args.sigUC), branches=[f"pho{N+1}_{var}" for N in range(4) for var in variables], eras=eras)

    # Get corrected samples
    samples_corr = loadPerMass(os.path.join(d, args.sigUC), branches=[f"pho{N+1}_{var}" for N in range(4) for var in variables], eras=eras)

    # Plot corr/uncorr samples for these variables
    for var in variables:
        plot2D(samples_uncorr, samples_corr, f"2D_phoX_{var}.png", args.year, os.path.join(d, "signal_tests", "noSel_compare2D"), Scales, var)
