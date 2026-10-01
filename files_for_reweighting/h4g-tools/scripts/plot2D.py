import os
import numpy as np
import pandas as pd
from typing import List, Dict, Optional
import argparse
from h4g_tools.utils.scales import Scales as ScalesCls
from ROOT import kBlack, kRed, kGreen, kOrange, kBlue, TCanvas, TH2D, TLegend
from h4g_tools.utils.loading import fillHist, loadPerMass
from h4g_tools.utils.plotting import plotFancy


def plot2D(
    samples: Dict[str, pd.DataFrame],
    filename: str,
    era: str,
    path: str,
    scales: ScalesCls,
    region: str,
    ybranch: str
) -> None:
    """
    Plot 2D histogram of pT/mgg vs mgg and pT vs mgg in both EB/EE.
    """ 

    sorted_samples = sorted(samples.items())
    for idx, (mass, sample) in enumerate(sorted_samples):
        # Run EB/EE configuration before plotting
        if region == "combined":
            cut_p1 = [True] * len(sample)
            cut_p2 = [True] * len(sample)
        if ybranch == "dipho_mass":
            if region == "EB":
                cut_p1 = sample.lead_isScEtaEB
                cut_p2 = sample.sublead_isScEtaEB
            elif region == "EE":
                cut_p1 = sample.lead_isScEtaEE
                cut_p2 = sample.sublead_isScEtaEE
            elif region == "2EB":
                cut_p1 = sample.lead_isScEtaEB & sample.sublead_isScEtaEB
                cut_p2 = sample.lead_isScEtaEB & sample.sublead_isScEtaEB
            elif region == "2EE":
                cut_p1 = sample.lead_isScEtaEE & sample.sublead_isScEtaEE
                cut_p2 = sample.lead_isScEtaEE & sample.sublead_isScEtaEE
            elif region == "1EB1EE":
                cut_p1 = (sample.lead_isScEtaEB & sample.sublead_isScEtaEE) | (sample.lead_isScEtaEE & sample.sublead_isScEtaEB)
                cut_p2 = (sample.lead_isScEtaEB & sample.sublead_isScEtaEE) | (sample.lead_isScEtaEE & sample.sublead_isScEtaEB)

        mass_in = int(mass[7:9])
        print(f"Mass: {mass_in} GeV")

        # Canvas setup
        canvas = TCanvas("canvas", "canvas", 1000, 1000)

        # Histogram with custom number of equal bins, bounded by (bounds[0], bounds[1])
        normalize = False

        if ybranch == "dipho_mass":
            names = [f"pT1_mgg_vs_mgg_{region}", f"pT2_mgg_vs_mgg_{region}", f"pT1_vs_mgg_{region}", f"pT2_vs_mgg_{region}"]
            branches = ["pT1_m_gg", "pT2_m_gg", "lead_pt", "sublead_pt"]
            xtitles = ["p_{T}^{#gamma_{1}} / m_{#gamma #gamma}", "p_{T}^{#gamma_{2}} / m_{#gamma #gamma}", "p_{T}^{#gamma_{1}}", "p_{T}^{#gamma_{2}}"]
            boundsList = [[0.0, 7.0, 0.0, 125.0], [0.0, 5.0, 0.0, 125.0], [0.0, 200.0, 0.0, 125.0], [0.0, 200.0, 0.0, 125.0]]
            binsScales = [[100.0/7.0, 1.0], [20.0, 1.0], [1.0, 1.0], [1.0, 1.0]]
            cut_list = [cut_p1, cut_p2, cut_p1, cut_p2]
        elif ybranch == "mass":
            names = [f"pT1_ma1_vs_ma1_{region}", f"pT2_ma1_vs_ma1_{region}", f"pT1_ma2_vs_ma2_{region}", f"pT2_ma2_vs_ma2_{region}", f"pT1_vs_ma1_{region}", f"pT2_vs_ma1_{region}", f"pT1_vs_ma2_{region}", f"pT2_vs_ma2_{region}"]
            branches = ["pT1_ma1", "pT2_ma1", "pT1_ma2", "pT2_ma2", "LeadPs_leading_pho_pt", "LeadPs_subleading_pho_pt", "SubleadPs_leading_pho_pt", "SubleadPs_subleading_pho_pt"]
            xtitles = ["p_{T}^{#gamma_{1}} / m_{a1}", "p_{T}^{#gamma_{2}} / m_{a1}", "p_{T}^{#gamma_{1}} / m_{a2}", "p_{T}^{#gamma_{2}} / m_{a2}", "p_{T}^{#gamma_{1}}", "p_{T}^{#gamma_{2}}", "p_{T}^{#gamma_{1}}", "p_{T}^{#gamma_{2}}"]
            boundsList = [[0.0, 10.0, 10.0, 65.0], [0.0, 10.0, 10.0, 65.0], [0.0, 10.0, 10.0, 65.0], [0.0, 10.0, 10.0, 65.0], [0.0, 160.0, 10.0, 65.0], [0.0, 160.0, 10.0, 65.0], [0.0, 160.0, 10.0, 65.0], [0.0, 160.0, 10.0, 65.0]]
            binsScales = [[10.0, 5.0], [10.0, 5.0], [10.0, 5.0], [10.0, 5.0], [1.0, 5.0], [1.0, 5.0], [1.0, 5.0], [1.0, 5.0]]
            cut_list = [[True] * len(sample) for _ in range(len(names))]

        assert len(names) == len(branches) == len(xtitles) == len(boundsList) == len(binsScales) == len(cut_list)
        for name, branch, xtitle, bounds, binsScale, cut in zip(names, branches, xtitles, boundsList, binsScales, cut_list):
            if ybranch == "mass":
                if "ma1" in name:
                    ybranch_loop = "LeadPs_mass"
                elif "ma2" in name:
                    ybranch_loop = "SubleadPs_mass"
            else:
                ybranch_loop = ybranch

            hist = TH2D(name, name, int((bounds[1]-bounds[0])*binsScale[0]), bounds[0], bounds[1], int((bounds[3]-bounds[2])*binsScale[1]), bounds[2], bounds[3])
            xbranch_df = sample[cut][branch].tolist()
            ybranch_df = sample[cut][ybranch_loop].tolist()

            # Fill histogram
            for x, y in zip(xbranch_df, ybranch_df):
                hist.Fill(x,y)

            # Normalize
            if normalize:
                if hist.Integral() > 0.0:
                    hist.Scale(1.0 / hist.Integral())
                else:
                    hist.Scale(1.0)
                    print("Cannot normalize empty histogram.")

            hist.SetStats(0)
            hist.SetMarkerStyle(8)
            hist.SetMarkerSize(0.5)
            hist.SetTitle(f"")
            hist.GetYaxis().SetTitleFont(42)
            if ybranch == "dipho_mass":
                hist.GetYaxis().SetTitle("m_{#gamma #gamma}" + f"  /  {1.0 / binsScale[1]: 1.3f}")
            elif ybranch_loop == "LeadPs_mass" or ybranch_loop == "SubleadPs_mass":
                if "ma1" in name:
                    hist.GetYaxis().SetTitle("m_{a1}" + f"  /  {1.0 / binsScale[1]: 1.3f}")
                elif "ma2" in name:
                    hist.GetYaxis().SetTitle("m_{a2}" + f"  /  {1.0 / binsScale[1]: 1.3f}")
            hist.GetXaxis().SetTitleFont(42)
            hist.GetXaxis().SetTitle(xtitle + f"  /  {1.0 / binsScale[0]: 1.3f}")
            hist.GetXaxis().SetTitleOffset(1.5)
            hist.Draw("COLZ")

            plotFancy(canvas, name, lumiTxt=Scales.lumiTxt.strip(), prelim=True, inPlot=False, colz=True)

            #new_filename = filename.replace("_replace", f"_{name}_{mass_in}GeV_{era}")
            new_filename = filename.replace("_replace", f"{name}_{mass_in}GeV_{era}")

            if not os.path.exists(f"{path}/{mass_in}_GeV/"):
                os.mkdir(f"{path}/{mass_in}_GeV/")
            canvas.SaveAs(f"{path}/{mass_in}_GeV/{new_filename}")
            assert os.path.isfile(f"{path}/{mass_in}_GeV/{new_filename}"), "Saved plot does not exist."

            for plot in [f"pT1_mgg_vs_mgg_{region}", f"pT2_mgg_vs_mgg_{region}", f"pT1_ma1_vs_ma1_{region}", f"pT2_ma1_vs_ma1_{region}", f"pT1_ma2_vs_ma2_{region}", f"pT2_ma2_vs_ma2_{region}"]:
                if name == plot:
                    canvas.Clear()

                    # X Projection
                    hist1D = hist.ProjectionX()
                    hist1D.SetStats(0)
                    hist1D.SetMarkerStyle(8)
                    hist1D.SetMarkerSize(0.5)
                    hist1D.SetTitle(f"")
                    hist1D.GetYaxis().SetTitleFont(42)
                    hist1D.GetYaxis().SetTitle(f"Events / {1.0 / binsScale[0]: 1.3f}")
                    hist1D.GetXaxis().SetTitleFont(42)
                    hist1D.GetXaxis().SetTitle(xtitle)
                    hist1D.Draw("hist")

                    plotFancy(canvas, f"{name}", sub="projX", lumiTxt=Scales.lumiTxt.strip(), prelim=True, inPlot=False, colz=True)

                    #new_filename = filename.replace("_replace", f"_{name}_{mass_in}GeV_{era}")
                    new_filename = filename.replace("_replace", f"{name}_projX_{mass_in}GeV_{era}")

                    if not os.path.exists(f"{path}/{mass_in}_GeV/projX/"):
                        os.mkdir(f"{path}/{mass_in}_GeV/projX/")
                    canvas.SaveAs(f"{path}/{mass_in}_GeV/projX/{new_filename}")
                    assert os.path.isfile(f"{path}/{mass_in}_GeV/projX/{new_filename}"), "Saved plot does not exist."
                    canvas.Clear()


                    # Legend setup
                    lsize = [0.65, 0.7, 0.9, 0.8]
                    legend = TLegend(lsize[0], lsize[1], lsize[2], lsize[3])
                    legend.SetTextFont(42)
                    legend.SetBorderSize(0)
                    legend.SetFillStyle(0)

                    # X Projection
                    hist1D = hist.ProjectionX()
                    hist1D.SetStats(0)
                    hist1D.SetTitle(f"")
                    hist1D.SetLineWidth(2)
                    hist1D.GetYaxis().SetTitleFont(42)
                    hist1D.GetYaxis().SetTitle(f"Events / {1.0 / binsScale[0]: 1.3f}")
                    hist1D.GetXaxis().SetTitleFont(42)
                    hist1D.GetXaxis().SetTitle(xtitle)
                    hist1D.SetLineColor(kRed)
                    legend.AddEntry(hist1D, f"Projection-X: {hist1D.Integral()}", "l")
                    hist1D.Draw("hist")

                    # To remove weird vertical line: apply plotFancy on one canvas then add additional histograms
                    plotFancy(canvas, f"{name}", lumiTxt=Scales.lumiTxt.strip(), prelim=True, inPlot=False)

                    # Original pT/mgg plot
                    original_hist = fillHist(branch, sample[cut], bounds=bounds[:2], binsScale=binsScale[0], normalize=False)
                    original_hist.SetStats(0)
                    original_hist.SetTitle(f"")
                    original_hist.SetLineWidth(2)
                    legend.AddEntry(original_hist, f"Original Plot {original_hist.Integral()}", "l")
                    original_hist.Draw("same hist")
                    legend.Draw("same")

                    #new_filename = filename.replace("_replace", f"_{name}_{mass_in}GeV_{era}")
                    new_filename = filename.replace("_replace", f"{name}_overlay_{mass_in}GeV_{era}")

                    if not os.path.exists(f"{path}/{mass_in}_GeV/overlay/"):
                        os.mkdir(f"{path}/{mass_in}_GeV/overlay/")
                    canvas.SaveAs(f"{path}/{mass_in}_GeV/overlay/{new_filename}")
                    assert os.path.isfile(f"{path}/{mass_in}_GeV/overlay/{new_filename}"), "Saved plot does not exist."
                


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", required=False, type=str, help="Only run photon plots for this subdirectory.")
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

        if directory != "backup" and directory != "2d_plots":
            if not os.path.exists(os.path.join(d, "2d_plots")):
                os.mkdir(os.path.join(d, "2d_plots"))

            #ybranch = "mass"
            for ybranch in ["mass", "dipho_mass"]:
                if ybranch == "dipho_mass":
                    samples = loadPerMass(f"{d}/{directory}", branches=["dipho_mass", "pT1_m_gg", "pT2_m_gg", "lead_isScEtaEB", "lead_isScEtaEE", "sublead_isScEtaEB", "sublead_isScEtaEE", "lead_pt", "sublead_pt"], eras=eras)
                elif ybranch == "mass":
                    samples = loadPerMass(f"{d}/{directory}", branches=["LeadPs_mass", "SubleadPs_mass", "pT1_ma1", "pT2_ma1", "pT1_ma2", "pT2_ma2", "LeadPs_leading_pho_pt", "LeadPs_subleading_pho_pt", "SubleadPs_leading_pho_pt", "SubleadPs_subleading_pho_pt"], eras=eras)
                for region in ["combined", "EB", "EE", "2EB", "2EE", "1EB1EE"] if ybranch == "dipho_mass" else ["combined"]:
                    plot2D(samples, f"_replace.png", args.year, os.path.join(d, "2d_plots"), Scales, region, ybranch=ybranch)
