import ROOT #adding this as I  also want to save a root file here

import os
import pandas as pd
import h4g_tools.utils.plotting as plotting
import h4g_tools.utils.loading as loading
import h4g_tools.utils.scales as scales

import sys

def plotInputs(
    signal_df: pd.DataFrame,
    bkg_df: pd.DataFrame,
    year: str,
) -> None:
    """
    Plots inputs directly as they are given to the event selection BDT. All signal MC is merged in this case.
    """

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    Scales = scales.Scales(year)
    masses = [x for x in range(15, 65, 5)]
    if year == "2022preEE":
        eras={
            "data": ["Run2022C", "Run2022D"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_preEE" for m in masses]
        }
    elif year == "2022postEE":
        eras={
            "data": ["Run2022E", "Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_postEE" for m in masses]
        }
    elif year == "2022":
        eras={
            "data": ["Run2022C", "Run2022D", "Run2022E", "Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D", "RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preEE", "postEE"]]
        }
    elif year == "2023preBPix":
        eras={
            "data": ["Run2023C"],
            "bkg": ["RunEvtMix2023C"],
            "signal": [f"Signal_{m}_GeV_preBPix" for m in masses]
        }
    elif year == "2022postBPix":
        eras={
            "data": ["Run2022D"],
            "bkg": ["RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_postBPix" for m in masses]
        }
    elif year == "2023":
        eras={
            "data": ["Run2023C", "Run2023D"],
            "bkg": ["RunEvtMix2023C", "RunEvtMix2023D"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preBPix", "postBPix"]]
        }
    elif "2024" in year:
        eras={
            "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
            "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
            "signal": [f"Signal_{m}_GeV" for m in masses]
        }
    else:
        assert False, "Must specify year in Scales class!"


    # Setup for lists
    branches, boundsList, binsScaleList, names, norms, titles = plotting.loadPlottingParameters()

    paths = [
        f"./plots/BDT/{year}/bdt_inputs/pho1_mvaID.png",
        f"./plots/BDT/{year}/bdt_inputs/pho2_mvaID.png",
        f"./plots/BDT/{year}/bdt_inputs/pho3_mvaID.png",
        f"./plots/BDT/{year}/bdt_inputs/pho4_mvaID.png",
        f"./plots/BDT/{year}/bdt_inputs/LeadPs_pt.png",
        f"./plots/BDT/{year}/bdt_inputs/SubleadPs_pt.png",
        f"./plots/BDT/{year}/bdt_inputs/dR_aa_mass_gggg.png",
        f"./plots/BDT/{year}/bdt_inputs/LeadPs_interMass.png",
        f"./plots/BDT/{year}/bdt_inputs/SubleadPs_interMass.png",
        f"./plots/BDT/{year}/bdt_inputs/LeadPs_interMass_abs.png",
        f"./plots/BDT/{year}/bdt_inputs/SubleadPs_interMass_abs.png",
        f"./plots/BDT/{year}/bdt_inputs/Ps_massDiff.png",
        f"./plots/BDT/{year}/bdt_inputs/cos_ag.png",
    ]

    branches = branches[:len(paths)]
    boundsList = boundsList[:len(paths)]
    binsScaleList = binsScaleList[:len(paths)]
    names = names[:len(paths)]
    norms = norms[:len(paths)]
    titles = titles[:len(paths)]

    assert len(branches) == len(boundsList) == len(binsScaleList) == len(names) == len(norms) == len(titles) == len(paths), f"{len(branches)} {len(boundsList)} {len(binsScaleList)} {len(names)} {len(norms)} {len(titles)} {len(paths)}"

    # Load data to compare with inputs (demonstrate RWing - or lack thereof - effects)
    data_path = os.path.join(*["/cms/cephfs/data/store/user/castells/bdtIO/", f"BDT_{year}_input_weights_data.parquet"])
    data_df = loading.loadSamples(data_path, eras=eras)

    # Make histograms for each variable
    hists = {}
    for branch, bounds, binsScale, name, norm, title, path in zip(branches, boundsList, binsScaleList, names, norms, titles, paths):
        if "interMass_abs" in path:
            continue

        hists.update({f"{branch}_sig_SB": loading.fillHist(branch, scales.SideBand(signal_df), bounds, binsScale=binsScale, name=f"{name}_signal_combined_SB", normalize=norm, customWeights="weight")})
        hists.update({f"{branch}_bkg_SB": loading.fillHist(branch, scales.SideBand(bkg_df), bounds, binsScale=binsScale, name=f"{name}_bkg_SB", normalize=norm, customWeights="weight")})
        hists.update({f"{branch}_data_SB": loading.fillHist(branch, scales.SideBand(data_df), bounds, binsScale=binsScale, name=f"{name}_data_SB", normalize=norm)})
        hists.update({f"{branch}_data": loading.fillHist(branch, scales.SideBand(data_df), bounds, binsScale=binsScale, name=f"{name}_data", normalize=norm)})
        hists.update({f"{branch}_sig": loading.fillHist(branch, signal_df, bounds, binsScale=binsScale, name=f"{name}_signal_combined", normalize=norm, customWeights="weight")})
        hists.update({f"{branch}_bkg": loading.fillHist(branch, bkg_df, bounds, binsScale=binsScale, name=f"{name}_bkg", normalize=norm, customWeights="weight")})
        hists.update({f"{branch}_bkg_pre_rw": loading.fillHist(branch, scales.SideBand(bkg_df), bounds, binsScale=binsScale, name=f"{name}_bkg_pre_rw", normalize=norm, customWeights="weight")})

        # Plotting SB will scale background to data and signal to 1fb. Absolute scale does not matter
        # Scale signal to background for full spectrum
        hists[f"{branch}_sig_SB"].Scale(hists[f"{branch}_data_SB"].Integral() / hists[f"{branch}_sig_SB"].Integral())
        hists[f"{branch}_bkg_SB"].Scale(hists[f"{branch}_data_SB"].Integral() / hists[f"{branch}_bkg_SB"].Integral())
        hists[f"{branch}_bkg_pre_rw"].Scale(hists[f"{branch}_data_SB"].Integral() / hists[f"{branch}_bkg_pre_rw"].Integral())
        hists[f"{branch}_bkg"].Scale(hists[f"{branch}_data"].Integral() / hists[f"{branch}_bkg"].Integral())
        hists[f"{branch}_sig"].Scale(hists[f"{branch}_bkg"].Integral() / hists[f"{branch}_sig"].Integral())

    #*********************************************** input histograms to a ROOT file ----
    outdir_root = f"./plots/BDT/{year}/bdt_inputs"
    os.makedirs(outdir_root, exist_ok=True)
    fout = ROOT.TFile(f"{outdir_root}/bdt_inputs_hists_{year}.root", "RECREATE")
    d_merged = fout.mkdir("merged")
    d_merged.cd()
    for key, h in hists.items():
        hc = h.Clone(key)
        hc.SetDirectory(d_merged)
        hc.SetStats(1)
        hc.Write()
    fout.cd()
    #--------------------*****************************************************

    max_order = 6
    for idx, (branch, title, path) in enumerate(zip(branches, titles, paths)):
        if "interMass_abs" in path:
            continue

        hists_SB = {
            "background_reweight": hists[f"{branch}_bkg_pre_rw"],
            "background": hists[f"{branch}_bkg_SB"],
            "data": hists[f"{branch}_data_SB"],
            "signal": hists[f"{branch}_sig_SB"],
        }

        hists_full = {
            "background": hists[f"{branch}_bkg"],
            "signal": hists[f"{branch}_sig"],
        }

        # Plot with data in sidebands only
        plotting.plot_data_sig_bkg_comparison(hists_SB, title, path.replace(".png", "_wDataSB.png"), sum([Scales.sumw[year[:4]][f"{mass}_GeV"] for mass in masses]), maximum=(10**max_order - 5*10**(max_order-1)), canvas_num = idx, reweight = True, sub = "Sideband Only", Scales = Scales, ignore_scale = True)

        # Plot signal and background only in full region
        plotting.plot_data_sig_bkg_comparison(hists_full, title, path, sum([Scales.sumw[year[:4]][f"{mass}_GeV"] for mass in masses]), maximum=(10**max_order - 5*10**(max_order-1)), canvas_num = idx, reweight = False, Scales = Scales, ignore_scale = True, omit_data = True)

    # ---- per-mass-point overlays: single-mass signal vs. (unchanged) evtmix bkg, full region ----*************
    # Recover each signal event's mass point from interMass = (m_reco - m_hyp)/m_gggg
    sig_mpt = ((signal_df["LeadPs_mass"] - signal_df["LeadPs_interMass"] * signal_df["mass_gggg"]) / 5).round() * 5

    # ---- events per signal mass point entering the BDT training ----
    counts = signal_df.groupby(sig_mpt).agg(
        n_events=("weight", "size"),
        sum_weights=("weight", "sum"),
    )
    counts["frac_events_%"] = 100.0 * counts["n_events"] / counts["n_events"].sum()
    counts["frac_weight_%"] = 100.0 * counts["sum_weights"] / counts["sum_weights"].sum()
    print(f"\n[plotInputs] Signal events entering training ({year}):")
    print(counts.to_string(float_format=lambda x: f"{x:.3f}"))
    print(f"Total: {counts['n_events'].sum()} events\n")
    #----------------------------------------------------------------

    os.makedirs(f"./plots/BDT/{year}/bdt_inputs/per_mass", exist_ok=True)
    counts.to_csv(f"./plots/BDT/{year}/bdt_inputs/per_mass/signal_training_counts.csv")

    for k, mass in enumerate(masses):
        sig_m = signal_df[sig_mpt == mass]
        if len(sig_m) == 0:
            print(f"[plotInputs] No signal events found for {mass} GeV, skipping.")
            continue
        print(f"[plotInputs] {mass} GeV: {len(sig_m)} signal events")

        outdir = f"./plots/BDT/{year}/bdt_inputs/per_mass/{mass}GeV"
        os.makedirs(outdir, exist_ok=True)

        for idx, (branch, bounds, binsScale, name, norm, title, path) in enumerate(zip(branches, boundsList, binsScaleList, names, norms, titles, paths)):
            if "interMass_abs" in path:
                continue

            h_sig = loading.fillHist(branch, sig_m, bounds, binsScale=binsScale, name=f"{name}_signal_{mass}GeV", normalize=norm, customWeights="weight")
            h_sig.Scale(hists[f"{branch}_bkg"].Integral() / h_sig.Integral())

            #*******************************************----------------------------
            d_mass = fout.GetDirectory(f"per_mass/{mass}GeV") or fout.mkdir(f"per_mass/{mass}GeV")
            
            hc = h_sig.Clone(f"{branch}_sig_{mass}GeV")
            hc.SetDirectory(d_mass)
            hc.SetStats(1) #without this get no statbox in my root file, dont remove
            hc.Write()
            #******************************************------------------------------

            hists_mass = {
                "background": hists[f"{branch}_bkg"],   # same evtmix bkg as the merged plot
                "signal": h_sig,
            }
            out_path = os.path.join(outdir, os.path.basename(path))
            plotting.plot_data_sig_bkg_comparison(hists_mass, f"{title} ({mass} GeV)", out_path, Scales.sumw[year[:4]][f"{mass}_GeV"], maximum=(10**max_order - 5*10**(max_order-1)), canvas_num=idx + 100*(k+1), reweight=False, Scales=Scales, ignore_scale=True, omit_data=True)

    fout.Close()

            #-----------------------------------------------------------------------------*************************
