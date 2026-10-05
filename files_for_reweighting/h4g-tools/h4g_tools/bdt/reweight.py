import os
from array import array
from numba import jit
import uproot
import numpy as np
import pandas as pd
import itertools
from typing import List, Dict, Union, Tuple
from ROOT import TH1D, TMath, TTree, TFile, gROOT, TCanvas, kRed, kBlue, kGreen, TLegend, kBlack, TH2D, gStyle, SetOwnership, THnSparseD
from h4g_tools.bdt.event_selection_bdt import define_variables
from h4g_tools.utils.loading import loadROOTSamples
from h4g_tools.utils.plotting import plot_Nreweighting, plot_Ndim_weights, loadPlottingParameters, plotFancy
from h4g_tools.utils.scales import Scales
import glob
from array import array
import sys

def reweight1Dim(
    bkg: pd.DataFrame,
    data: pd.DataFrame,
    col: str = None,
) -> TH1D:
    """
    Reweighting background events wrt data in sidebands to improve BDT score agreement
    """

    nBins = 10
    BDT_score_hist = TH1D("BDT_score_bkg", "BDT_score_bkg", nBins, 0.0, 1.0)

    perm = []
    # + 2 to account for under/overflow bins
    for bn in range(0,nBins+2):
        perm.append(bn)

    weight_data = {}
    weight_bkg = {}
    weight_ratio = {}
    for subset in perm:
        # Data/evtmix are set to 0.0 so the true values from Events tree are used
        weight_data[subset] = 0.0
        weight_bkg[subset] = 0.0
        # Ratio is set to 1.0 
        weight_ratio[subset] = 1.0

    # Using the ratio of the cut region that is being reweighted!
    sb_bkg_entries = len(bkg[(bkg["mass_gggg"] <= 115.0) | (bkg["mass_gggg"] >= 135.0)]) # used for scaling
    sb_data_entries = len(data[(data["mass_gggg"] <= 115.0) | (data["mass_gggg"] >= 135.0)])

    # handles division-by-zero issue for testing MVA cuts
    sb_scaling_weight = 1.0
    if sb_bkg_entries > 0.0:
        sb_scaling_weight = sb_data_entries / sb_bkg_entries

    #@jit
    def sparseBins1D(
        hist: TH1D,
        df: pd.DataFrame,
        weight_dict: dict,
        sidebands: bool = True,
        weight_scale: float = 1.0,
    )-> dict:
        """
        Populate sparse bins with data from a TTree in the sideband region.
        """

        entries = len(df)
        sbcut = (df["mass_gggg"] <= 115.0) | (df["mass_gggg"] >= 135.0)

        for idx in range(0, entries):
            # Reweight from sidebands only
            if (sbcut.iloc[idx] if sidebands else True):
                ibin = hist.FindBin(df[col].iloc[idx])
                weight_dict[ibin] += 1.0 * weight_scale

        return weight_dict


    # Fill sparse bins for data
    print("Filling sparse bins for data...")
    weight_data = sparseBins1D(BDT_score_hist, data, weight_data, sidebands=True)

    # Fill sparse bins for bkg
    print("Filling sparse bins for bkg...")
    weight_bkg = sparseBins1D(BDT_score_hist, bkg, weight_bkg, sidebands=True, weight_scale=sb_scaling_weight)

    # Get ratio for reweights
    print("Calculating ratio of data/bkg...")

    # Fill 2D hist even if the weights don't hit bkg events (i.e. w = 0.2 and 0.3)
    for subset in weight_ratio.keys():
        if weight_bkg[subset] > 0.0 and weight_data[subset] > 0.0:
            # Regular ratio (excludes zero data to not remove background bins)
            weight_ratio[subset] = weight_data[subset] / weight_bkg[subset]
        elif weight_bkg[subset] > 0.0 and weight_data[subset] == 0.0:
            # Number of bins with zero data and non-zero event mixing. Don't remove event mixed background if data is nonzero.
            weight_ratio[subset] = 1.0

    # apply weights to each event
    weights = []
    for idx in range(0, len(bkg)):
        ibin = BDT_score_hist.FindBin(bkg[col].iloc[idx])
        weights.append(weight_ratio[ibin])

    assert len(weights) == len(bkg), "{len(weights)} {len(bkg)}"
    reweight_df = pd.concat([bkg, pd.Series(weights, name="1dimreweight")], axis=1)

    return reweight_df


def findVal(
    tree: TTree,
    var: str
) -> float:
    """
    Get value of variable from a TTree from a pre-set entry.
    """

    if var == "pho1_mvaID":
        return tree.pho1_mvaID
    elif var == "pho2_mvaID":
        return tree.pho2_mvaID
    elif var == "pho3_mvaID":
        return tree.pho3_mvaID
    elif var == "pho4_mvaID":
        return tree.pho4_mvaID
    elif var == "LeadPs_pt":
        return tree.LeadPs_pt
    elif var == "SubleadPs_pt":
        return tree.SubleadPs_pt
    elif var == "dR_aa_mass_gggg":
        return tree.dR_aa_mass_gggg
    elif var == "dR_aa":
        return tree.dR_aa
    elif var == "LeadPs_interMass":
        return tree.LeadPs_interMass
    elif var == "SubleadPs_interMass":
        return tree.SubleadPs_interMass
    elif var == "Ps_massDiff":
        return tree.Ps_massDiff
    elif var == "cos_ag":
        return tree.cos_ag
    else:
        assert False


#@jit
def getEntries(
    tree: TTree,
    region: str,
    mva_cut: float = -1
)-> dict:
    """
    Get entries in sideband to properly scale the new weights for bkg.
    """

    counter = 0.0
    entries = tree.GetEntries()

    for idx in range(0, entries):
        tree.GetEntry(idx)

        # Edges of 110.0 and 180.0 GeV are handled by HiggsDNA!
        if region.lower() == "sb":
            # Get entries from sidebands only. Generally what we want. Others are for testing.
            cut = tree.mass_gggg <= 115.0 or tree.mass_gggg >= 135.0
        elif region.lower() == "sr":
            # Get entries from signal region only. For testing.
            cut = tree.mass_gggg > 115.0 and tree.mass_gggg < 135.0
        elif region.lower() == "full":
            cut = True
        elif region.lower() == "sb_mva":
            # Get entries from signal region only. For testing.
            cut = tree.mass_gggg > 115.0 and tree.mass_gggg < 135.0 and tree.pho3_mvaID < mva_cut and tree.pho4_mvaID < mva_cut

        if cut:
            counter += 1.0

    return counter


#@jit
def sparseBins(
    tree: TTree,
    plots: list,
    hist: THnSparseD,
    sidebands: bool = True,
    weight_scale: float = 1.0,
):
    entries = tree.GetEntries()
    values = array("d", [0.0] * len(plots))

    for idx in range(entries):
        tree.GetEntry(idx)

        sbcut = tree.mass_gggg <= 115.0 or tree.mass_gggg >= 135.0

        if sbcut if sidebands else True:
            for ivar, plot in enumerate(plots):
                values[ivar] = findVal(tree, plot[0])

            hist.Fill(values, weight_scale)


def loadNDimWeights(
    year: str
) -> pd.Series:
    """
    Loads N-dim reweight values for a given year. The reweighting process must be completed and the file must exist for this to work.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Weights are already converted to parquet in the previous step
    weights_path = f"{cwd}/../../scripts/bdtIO/BDT_{year}_input_weights.parquet"
    weights_path_data = f"{cwd}/../../scripts/bdtIO/BDT_{year}_input_weights_data.parquet"
    assert os.path.exists(weights_path)
    assert os.path.exists(weights_path_data)

    # Read in both bkg and data to multiply bkg event weights by data/bkg
    weights_df = pd.read_parquet(weights_path)
    weights_data_df = pd.read_parquet(weights_path_data)
    weights_df["Ndimreweight"] = (len(weights_data_df) / len(weights_df)) * weights_df["Ndimreweight"]
    print("Read weights from file (normalized to data).")

    return weights_df


def performNDimReweighting(
    b: str,
    d: str,
    year: str,
    eras: list,
    nBins: int = None,
    mvaCut: float = None,
    EB_only: bool = False,
    save_default: bool = False,
    default_weight: float = 0.1
) -> None:
    """
    Performs N-dim reweighting for event mixing from data/event mixing in sideband region. (Used to take ~30 minutes before jit).
    """

    makeWeightPlots = False   # need to alter this part; for now kept False

    # Load samples from ROOT conversions
    # Eras must already be split in directories!
    print("Loading background samples:")
    bkg = loadROOTSamples(b)
    if mvaCut is not None:
        assert type(mvaCut) is float, type(mvaCut)
        bkg = bkg.CopyTree(f"pho3_mvaID >= {mvaCut} & pho4_mvaID >= {mvaCut}")
    if EB_only:
        bkg = bkg.CopyTree("pho1_isScEtaEB && pho2_isScEtaEB && pho3_isScEtaEB && pho4_isScEtaEB")
    SetOwnership(bkg, False)

    print("Loading data samples:")
    data = loadROOTSamples(d)
    if mvaCut is not None:
        assert type(mvaCut) is float, type(mvaCut)
        data = data.CopyTree(f"pho3_mvaID >= {mvaCut} & pho4_mvaID >= {mvaCut}")
    if EB_only:
        data = data.CopyTree("pho1_isScEtaEB && pho2_isScEtaEB && pho3_isScEtaEB && pho4_isScEtaEB")
    SetOwnership(data, False)

    # Get current directory and set weights path
    cwd = os.path.dirname(os.path.abspath(__file__))
    weights_path = f"{cwd}/../../scripts/bdtIO/BDT_{year}_input_weights.root"

    nBins_default_2022 = 8
    nBins_default_2024 = 10
    if nBins is None:
        if "2022" in year:
            nBins = nBins_default_2022  # The total # of bins is nBins+2 to account for under/overflow bins in the TH1s
        elif "2024" in year:
            nBins = nBins_default_2024
    print(f"Using {nBins}+2 bins for TH1s.")

    plots = []
    mva_plot_indices = []
    branches, boundsList, _, _, _, titles = loadPlottingParameters()
    dims = ["LeadPs_pt", "SubleadPs_pt", "Ps_massDiff", "dR_aa_mass_gggg", "pho1_mvaID", "pho2_mvaID", "pho3_mvaID", "pho4_mvaID"]

    for idx, branch in enumerate(branches):
        if branch in dims:
            plots.append([branches[idx], titles[idx], nBins, boundsList[idx][0], boundsList[idx][1]])
            if "mvaID" in branch:
                mva_plot_indices.append(idx)

    histos = []
    # variable-bin edges 
    edges_list = []

    for idx, plot in enumerate(plots):
        if idx not in mva_plot_indices:
            histos.append(TH1D(plot[0], plot[1], plot[2], plot[3], plot[4]))
            edges_list.append(None)
        else:
            if "pho3" in plot[0] or "pho4" in plot[0]:
                mva_bins = array("d", [-1, -0.93, -0.45, 0.2, 1])
            else:
                mva_bins = array("d", [-1, 1])

            histos.append(TH1D(plot[0], plot[1], len(mva_bins) - 1, mva_bins))
            edges_list.append(mva_bins)

    ndim = len(histos)

    # renamed nbins to axis_nbins
    axis_nbins = array("i", [hist.GetNbinsX() for hist in histos])
    xmin = array("d", [hist.GetXaxis().GetXmin() for hist in histos])
    xmax = array("d", [hist.GetXaxis().GetXmax() for hist in histos])

    weight_data = THnSparseD("weight_data", "Data", ndim, axis_nbins, xmin, xmax)
    weight_evtmix = THnSparseD("weight_evtmix", "Event Mixing", ndim, axis_nbins, xmin, xmax)
    weight_ratio = THnSparseD("weight_ratio", "Data / Event Mixing", ndim, axis_nbins, xmin, xmax)

    # SetBinEdges
    for idx, edges in enumerate(edges_list):
        if edges is not None:
            for h in (weight_data, weight_evtmix, weight_ratio):
                h.SetBinEdges(idx, edges)

    print("Filling sparse bins for data...")
    sparseBins(data, plots, weight_data)

    print("Filling sparse bins for bkg...")
    sb_scaling_weight = None
    
    # Event-mixed sb scaled to data sb before reweighting
    if Scales.useScaledCounts:
        print("Scaling counts to data/bkg ratio in SB region.")

        evtmix_entries = getEntries(bkg, region="full")
        data_entries = getEntries(data, region="full")
        sr_evtmix_entries = getEntries(bkg, region="sr")
        sr_data_entries = getEntries(data, region="sr")
        sb_evtmix_entries = getEntries(bkg, region="sb")
        sb_data_entries = getEntries(data, region="sb")

        print(f"\t\t {'Event Mixing': <12}   {'Data': <12}")
        print(f"{'Full region: ':<18} {evtmix_entries:<9} {data_entries:<9}")
        print(f"{'Signal region: ':<18} {sr_evtmix_entries:<9} {sr_data_entries:<9}")
        print(f"{'Sideband region: ':<18} {sb_evtmix_entries:<9} {sb_data_entries:<9}")

        sb_scaling_weight = sb_data_entries / sb_evtmix_entries

        sparseBins(bkg, plots, weight_evtmix, weight_scale=sb_scaling_weight)
    else:
        sparseBins(bkg, plots, weight_evtmix)

    print("Calculating ratio of data/bkg...")

    hist_counts = {
        "zeros": 0,
        "neither": 0,
        "no_data_yes_bkg": 0,
        "yes_data_no_bkg": 0,
    }

    if makeWeightPlots:
        yrange = (0, 2)
        ybins = 200

        hist2D = TH2D(
            "hist2D",
            f"Weights (bins {nBins+2})",
            weight_evtmix.GetNbins(),
            0,
            weight_evtmix.GetNbins(),
            ybins,
            yrange[0],
            yrange[1],
        )

    coord = array("i", [0] * ndim)

    for ibin in range(weight_evtmix.GetNbins()):
        evtmix_value = weight_evtmix.GetBinContent(ibin, coord)
        data_value = weight_data.GetBinContent(coord)

        if evtmix_value > 0.0 and data_value > 0.0:
            ratio = data_value / evtmix_value
            hist_counts["zeros"] += 1

        elif evtmix_value > 0.0 and data_value == 0.0:
            ratio = default_weight
            hist_counts["no_data_yes_bkg"] += 1

        else:
            continue

        weight_ratio.SetBinContent(coord, ratio)

        if makeWeightPlots:
            hist2D.Fill(ibin, ratio, evtmix_value)

    # Writing weights as new branch in identical event mixing TTree - saved as parquet for sorting with other DataFrames
    print(f"Converting event mixing dataset and weights to DataFrame for {len(dims)} dimensions...")
    evtmix_tree = TTree("Events", "Events")
    SetOwnership(evtmix_tree, False)
    mix_weight = array("f", [0])
    bname = "Ndimreweight"
    _mix_weight = evtmix_tree.Branch(bname, mix_weight, f"{bname}/F")
    
    print("check 1")

    def fillObj(obj_to_fill, obj_to_read):


        entries = obj_to_read.GetEntries()

        values = array("d", [0.0] * len(plots))
        coord_ev = array("i", [0] * len(plots))     # fixes the segmentation fault
        for idx in range(entries):
            obj_to_read.GetEntry(idx)

            for ivar, plot in enumerate(plots):
                values[ivar] = findVal(obj_to_read, plot[0])
                coord_ev[ivar] = histos[ivar].FindBin(values[ivar])

            w = weight_ratio.GetBinContent(coord_ev)   # returns 0 for unfilled bins, no memory allocation
            mix_weight[0] = w if w > 0.0 else 1.0

            obj_to_fill.Fill()
            
    fillObj(evtmix_tree, bkg)

    print("check 2")
    
    assert evtmix_tree.GetEntries() == bkg.GetEntries(), (
        evtmix_tree.GetEntries(),
        bkg.GetEntries(),
    )
    print("check 3")
    outfile = TFile(weights_path, "RECREATE")
    outfile.cd()
    evtmix_tree.Write()
    print("check 4")
    outfile.Close()
    print("Wrote weights to event mixed tree!")

    # This block still uses the old dict (len(weight_ratio), .items()) and will NOT work with THnSparse. for now not touched as makeWeightPlots = False.------
    if makeWeightPlots:
        zeros = TH1D("zeros", f"Weights (bins {nBins}+2)", len(weight_ratio), 0, len(weight_ratio))
        neither = TH1D("neither", f"Weights (bins {nBins}+2)", len(weight_ratio), 0, len(weight_ratio))
        no_data_yes_bkg = TH1D("no_data_yes_bkg", f"Weights (bins {nBins}+2)", len(weight_ratio), 0, len(weight_ratio))
        yes_data_no_bkg = TH1D("yes_data_no_bkg", f"Weights (bins {nBins}+2)", len(weight_ratio), 0, len(weight_ratio))
        for bn,w in weight_ratio.items():
            if w == 0.2:
                zeros.Fill(f"{bn}", -1)
                neither.Fill(f"{bn}", 1.0)
                no_data_yes_bkg.Fill(f"{bn}", -1)
                yes_data_no_bkg.Fill(f"{bn}", -1)
            elif w == 0.3:
                zeros.Fill(f"{bn}", -1)
                neither.Fill(f"{bn}", -1)
                no_data_yes_bkg.Fill(f"{bn}", -1)
                yes_data_no_bkg.Fill(f"{bn}", 1.0)
            elif w == 1.0:
                zeros.Fill(f"{bn}", -1)
                neither.Fill(f"{bn}", -1)
                no_data_yes_bkg.Fill(f"{bn}", w)
                yes_data_no_bkg.Fill(f"{bn}", -1)
            else:
                zeros.Fill(f"{bn}", w)
                neither.Fill(f"{bn}", -1)
                no_data_yes_bkg.Fill(f"{bn}", -1)
                yes_data_no_bkg.Fill(f"{bn}", -1)

        assert sum(hist_counts.values()) == (nBins+2)**len(dims), (sum(hist_counts.values()), f"Total {nBins+2}D bins: {(nBins+2)**len(dims)}")

        gROOT.SetBatch(True)
        canvas = TCanvas(f"canvas_{year}", f"canvas_{year}", 1000, 1000)
        canvas.SetLeftMargin(0.2)

        lsize = [0.55, 0.8, 0.99, 0.925]
        legend = TLegend(lsize[0], lsize[1], lsize[2], lsize[3])

        for hist in [zeros, neither, no_data_yes_bkg, yes_data_no_bkg, hist2D]:
            hist.GetYaxis().SetRangeUser(-0.18, 2.5)
            hist.SetTitleSize(10)
            hist.GetYaxis().SetTitle("Data/Bkg")
            hist.GetXaxis().SetTitle("4D Bin")
            hist.GetXaxis().SetTitleOffset(1.2)
            hist.GetXaxis().LabelsOption("u")
            hist.SetStats(0)

            for bn in range(hist.GetNbinsX()+1):
                hist.SetBinError(bn, 0.0)
                if bn not in [hist.GetNbinsX()/(nBins+2) * n for n in range(nBins+3)]:
                    if bn == 1:
                        pass
                    else:
                        hist.GetXaxis().SetBinLabel(bn, "")

        neither.SetMarkerStyle(8)
        neither.SetMarkerSize(0.5)
        neither.SetMarkerColor(kRed)
        no_data_yes_bkg.SetMarkerStyle(8)
        no_data_yes_bkg.SetMarkerSize(0.5)
        no_data_yes_bkg.SetMarkerColor(kBlue)
        yes_data_no_bkg.SetMarkerStyle(8)
        yes_data_no_bkg.SetMarkerSize(0.5)
        yes_data_no_bkg.SetMarkerColor(kGreen+1)
        zeros.SetMarkerStyle(8)
        zeros.SetMarkerSize(0.5)
        zeros.SetMarkerColor(kBlack)

        legend.AddEntry(zeros, f"Ratio (N = {hist_counts['zeros']})", "p")
        legend.AddEntry(yes_data_no_bkg, f"Zero Bkg, Nonzero Data (N = {hist_counts['yes_data_no_bkg']})", "p")
        legend.AddEntry(no_data_yes_bkg, f"Nonzero Bkg, Zero Data (N = {hist_counts['no_data_yes_bkg']})", "p")
        legend.AddEntry(neither, f"Zero Bkg, Zero Data (N = {hist_counts['neither']})", "p")

        neither.Draw("P same")
        no_data_yes_bkg.Draw("P same")
        yes_data_no_bkg.Draw("P same")
        zeros.Draw("P same")
        legend.Draw("same")

        neither.GetXaxis().SetTickLength(0.0)
        no_data_yes_bkg.GetXaxis().SetTickLength(0.0)
        yes_data_no_bkg.GetXaxis().SetTickLength(0.0)
        zeros.GetXaxis().SetTickLength(0.0)

        canvas.Update()
        canvas.SaveAs(f"{cwd}/../../scripts/plots/BDT/{year}/reweight_zeros_{nBins+2}bins.png")

        canvas.Clear()

        hist2D.SetStats(0)
        hist2D.SetMarkerStyle(8)
        hist2D.SetMarkerSize(0.5)
        hist2D.SetTitle("")
        hist2D.GetYaxis().SetTitleFont(42)
        hist2D.GetXaxis().SetTitleFont(42)
        hist2D.GetXaxis().SetTitleOffset(1.5)
        hist2D.GetXaxis().SetTickLength(0.0)
        hist2D.Draw("COLZ")

        plotFancy(canvas, f"Bkg events per weight (bins {nBins+2})", lumiTxt="", prelim=True, inPlot=False, colz=True)
        canvas.SaveAs(f"{cwd}/../../scripts/plots/BDT/{year}/reweight_2D_{nBins+2}bins.png")

    # Convert data tree back to parquet so no changes need to be made to BDT loading--------------
    dlist = glob.glob(os.path.join(d, "Run2024*", "nominal", "*.root"))

    data_df = uproot.concatenate(
        [f"{data_file}:Events" for data_file in dlist],
        library="pd"
    )

    if mvaCut is not None:
        data_df = data_df[(data_df.pho3_mvaID >= mvaCut) & (data_df.pho4_mvaID >= mvaCut)]

    if EB_only:
        data_df = data_df[data_df.pho1_isScEtaEB & data_df.pho2_isScEtaEB & data_df.pho3_isScEtaEB & data_df.pho4_isScEtaEB]

    if not save_default:
        data_path = weights_path.replace(".root", "_data.parquet")
        data_df.to_parquet(data_path)
    else:
        data_path = weights_path.replace(".root", f"_{str(default_weight).replace('.','p')}.root").replace("weights_", "weights_data_").replace(".root", ".parquet").replace("bdtIO",f"chi2_checks/{year}/")
        if ("2022" in year and nBins != nBins_default_2022) or (year == "2024" and nBins != nBins_default_2024):
            data_path = data_path.replace(".parquet", f"_nBins{nBins}.parquet")
        data_df.to_parquet(data_path)

    print(f"Saved data to {data_path}")

    
    # Do the same for the event mixing tree
    '''
    blist = [name for name in os.listdir(b) if os.path.isfile(os.path.join(b,name))]
    assert len(blist) == 1, blist
    bkg_file = os.path.join(b, blist[0])

    with uproot.open(bkg_file) as f:
        evtmix_df = f["Events"].arrays(library="pd")
        if mvaCut is not None:
            evtmix_df = evtmix_df[(evtmix_df.pho3_mvaID >= mvaCut) & (evtmix_df.pho4_mvaID >= mvaCut)]
        if EB_only:
            evtmix_df = evtmix_df[evtmix_df.pho1_isScEtaEB & evtmix_df.pho2_isScEtaEB & evtmix_df.pho3_isScEtaEB & evtmix_df.pho4_isScEtaEB]
        assert len(evtmix_df) == bkg.GetEntries(), (len(evtmix_df), bkg.GetEntries())
    with uproot.open(weights_path) as f:
        weights_df = f["Events"].arrays(library="pd")

    assert len(weights_df) == len(evtmix_df), (len(weights_df), len(evtmix_df))
    evtmix_plus_weights = pd.concat([weights_df["Ndimreweight"], evtmix_df.reset_index()], axis=1)
    if not save_default:
        evtmix_path = weights_path.replace(".root", ".parquet")
        evtmix_plus_weights.to_parquet(evtmix_path)
    else:
        evtmix_path = weights_path.replace(".root", f"_{str(default_weight).replace('.','p')}.root").replace(".root", ".parquet").replace("bdtIO",f"chi2_checks/{year}/")
        if ("2022" in year and nBins != nBins_default_2022) or (year == "2024" and nBins != nBins_default_2024):
            evtmix_path = evtmix_path.replace(".parquet", f"_nBins{nBins}.parquet")
        evtmix_plus_weights.to_parquet(evtmix_path)
    print(f"Saved data to {evtmix_path}")
    '''

    blist = glob.glob(os.path.join(b, "RunEvtMix2024*", "nominal", "*.root"))

    evtmix_df = uproot.concatenate(
        [bkg_file + ":Events" for bkg_file in blist],
        library="pd"
    )

    if mvaCut is not None:
        evtmix_df = evtmix_df[(evtmix_df.pho3_mvaID >= mvaCut) & (evtmix_df.pho4_mvaID >= mvaCut)]

    if EB_only:
        evtmix_df = evtmix_df[evtmix_df.pho1_isScEtaEB & evtmix_df.pho2_isScEtaEB & evtmix_df.pho3_isScEtaEB & evtmix_df.pho4_isScEtaEB]

    assert len(evtmix_df) == bkg.GetEntries(), (len(evtmix_df), bkg.GetEntries())

    with uproot.open(weights_path) as f:
        weights_df = f["Events"].arrays(library="pd")

    assert len(weights_df) == len(evtmix_df), (len(weights_df), len(evtmix_df))
    evtmix_plus_weights = pd.concat([weights_df["Ndimreweight"], evtmix_df.reset_index()], axis=1)

    if not save_default:
        evtmix_path = weights_path.replace(".root", ".parquet")
        evtmix_plus_weights.to_parquet(evtmix_path)
    else:
        evtmix_path = weights_path.replace(".root", f"_{str(default_weight).replace('.','p')}.root").replace(".root", ".parquet").replace("bdtIO",f"chi2_checks/{year}/")
        if ("2022" in year and nBins != nBins_default_2022) or (year == "2024" and nBins != nBins_default_2024):
            evtmix_path = evtmix_path.replace(".parquet", f"_nBins{nBins}.parquet")
        evtmix_plus_weights.to_parquet(evtmix_path)

    print(f"Saved data to {evtmix_path}")
    print(f"Total sum of weights: {sum(evtmix_plus_weights.Ndimreweight):.4f}")
    print("~~~~    Finished reweighting procedure!    ~~~~")
