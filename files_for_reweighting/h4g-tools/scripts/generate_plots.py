import gc
import numpy as np
import copy
import argparse
import os
from ROOT import TMath
from h4g_tools.utils.save_plots_to_root import savePlots
from h4g_tools.utils.scales import Scales as ScalesCls, SideBand, SignalRegion
from h4g_tools.bdt.event_selection_bdt import redefine_dataframe
from h4g_tools.utils.runner_utils import get_parser_plotting, handle_errors
from h4g_tools.utils.standard_plots import plotPsMasses, plot4Object, plotPsKinematics, plotPsMassCombined, plotHMassCombined, plotMHyp, plotdRCombined, massChecks
from h4g_tools.utils.plotting import plot_data_sig_bkg_comparison, plot_data_sig_bkg_comparison_combined, compare_four_masses_data_sig_bkg, loadPlottingParameters, compare_eras_data_bkg_only, compare_four_masses_data_sig_bkg_reweight_overlay, plot_corr_matrix, plotFourMassesCompared, make_scatter_mva, compare_four_masses_mva, compare_sigMixing_bkg
from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist
from h4g_tools.bdt.reweight import loadNDimWeights


"""
For generating BDT inputs plus others (basicallly everything!), omit the -bdt option. Remember to run reweighting first!
python3 generate_plots.py -i /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2022_SaS_Pileup_NFs_eVetoSF_Material_FNUF_26Nov2025/ -bkg /cms/cephfs/data/store/user/castells/bdtIO/BDT_2022_input_weights.parquet -d /cms/cephfs/data/store/user/castells/bdtIO/BDT_2022_input_weights_data.parquet -y 2022 -bdt

python3 generate_plots.py -i /cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_wBDT_genMatching_06Feb2026/ -bkg /cms/cephfs/data/store/user/castells/bdtIO/BDT_2024_input_weights.parquet -d /cms/cephfs/data/store/user/castells/bdtIO/BDT_2024_input_weights_data.parquet -y 2024 -bdt

Add -sigOnly to run only signal plots.
"""



def generatePlots(
    args: argparse.ArgumentParser
) -> None:

    handle_errors(
        ["Need at least one argument for samples.", args.inputs is None, args.data is None, args.background is None],
        ["Signal samples are required.", args.inputs is None],
        #["Need background to generate plots.", args.inputs is not None, args.data is not None, args.background is None],
        #["Need data to generate plots.", args.inputs is not None, args.data is None, args.background is not None],
        ["Must provide a year for processing.", args.year is None]
    )

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    masses = [x for x in range(15, 65, 5)]
    if args.year == "2022preEE":
        eras={
            "data": ["Run2022C", "Run2022D"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D"],
            "signal": [f"Signal_{m}_GeV_preEE" for m in masses]
        }
    elif args.year == "2022postEE":
        eras={
            "data": ["Run2022E", "Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_postEE" for m in masses]
        }
    elif args.year == "2022":
        eras={
            "data": ["Run2022C", "Run2022D", "Run2022E", "Run2022F", "Run2022G"],
            "bkg": ["RunEvtMix2022C", "RunEvtMix2022D", "RunEvtMix2022E", "RunEvtMix2022F", "RunEvtMix2022G"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preEE", "postEE"]]
        }
    elif args.year == "2023preBPix":
        eras={
            "data": ["Run2023C"],
            "bkg": ["RunEvtMix2023C"],
            "signal": [f"Signal_{m}_GeV_preBPix" for m in masses]
        }
    elif args.year == "2023postBPix":
        eras={
            "data": ["Run2023D"],
            "bkg": ["RunEvtMix2023D"],
            "signal": [f"Signal_{m}_GeV_postBPix" for m in masses]
        }
    elif args.year == "2023":
        eras={
            "data": ["Run2023C", "Run2023D"],
            "bkg": ["RunEvtMix2023C", "RunEvtMix2023D"],
            "signal": [f"Signal_{m}_GeV_{prepost}" for m in masses for prepost in ["preBPix", "postBPix"]]
        }
    elif "2024" in args.year:
        eras={
            "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
            "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
            "signal": [f"Signal_{m}_GeV" for m in masses]
        }
        year_adjust = "2024"
    else:
        assert False, "Must specify year in Scales class!"

    # Handle path setup after checking year
    Scales = ScalesCls(year_adjust) 
    pathSetup(args.year)

    if args.inputs is None:
        return

    print("Loading signal samples...")
    samples = loadPerMass(args.inputs, eras = eras if "outputs/BDT" not in args.inputs else None, year = args.year)
    to_del = []
    for key in samples.keys():
        if "Signal_" not in key:
            to_del.append(key)

    for key in to_del:
        samples.update({f"Signal_{key}": samples[key]})
        del samples[key]
    
    #massChecks(samples, Scales=Scales)
    savePlots(
        "noRatio",
        plot4Object(samples, ["mass_gggg", "sumPt_gggg", "dM", "dR_aa_mass_gggg", "cos_ag", "LeadPs_interMass", "SubleadPs_interMass", "LeadPs_interMass_abs", "SubleadPs_interMass_abs", "Ps_massDiff", "pT1_ma1", "pT2_ma1", "pT1_ma2", "pT2_ma2"], Scales=Scales),
        ["mass_gggg", "sumPt_gggg", "dM", "dR_aa_mass_gggg", "cos_ag", "LeadPs_interMass", "SubleadPs_interMass", "LeadPs_interMass_abs", "SubleadPs_interMass_abs", "Ps_massDiff", "pT1_ma1", "pT2_ma1", "pT1_ma2", "pT2_ma2"],
        "paper_plots/histograms_dataSB.root" if not args.fullRegion else "paper_plots/histograms.root"
    )
    savePlots(
        "noRatio",
        plotPsMasses(samples, ["LeadPs_mass", "SubleadPs_mass"], Scales=Scales),
        ["LeadPs_mass", "SubleadPs_mass"],
        "paper_plots/histograms_dataSB.root" if not args.fullRegion else "paper_plots/histograms.root"
    )
    # Add cos_ag_star if it's in the samples!
    savePlots(
        "noRatio",
        plotPsKinematics(samples, ["pt", "eta", "phi", "mvaID", "dR_gg", "pT1_m_gg", "pT2_m_gg", "leading_pho_pt", "subleading_pho_pt", "hoe", "sieie", "r9", "pfPhoIso03", "pfPhoIso03_rhoCorrected", "trkSumPtHollowConeDR03", "pixelSeed", "pfRelIso_photon_pt", "pfRelIso03_chg_quadratic"], Scales=Scales),
        ["pt", "eta", "phi", "mvaID", "dR_gg", "pT1_m_gg", "pT2_m_gg", "leading_pho_pt", "subleading_pho_pt", "hoe", "sieie", "r9", "pfPhoIso03", "pfPhoIso03_rhoCorrected", "trkSumPtHollowConeDR03", "pixelSeed", "pfRelIso_photon_pt", "pfRelIso03_chg_quadratic"],
        "paper_plots/histograms_dataSB.root" if not args.fullRegion else "paper_plots/histograms.root"
    )
    savePlots(
        "combined",
        plotPsMassCombined(samples, ["LeadPs_mass", "SubleadPs_mass"], Scales=Scales),
        ["LeadPs_mass", "SubleadPs_mass"]*4,
        "paper_plots/histograms_dataSB.root" if not args.fullRegion else "paper_plots/histograms.root"
    )
    savePlots(
        "combined",
        plotHMassCombined(samples, "mass_gggg", Scales=Scales),
        [f"mass_gggg_{m}" for m in [15, 30, 40, 60]],
        "paper_plots/histograms_dataSB.root" if not args.fullRegion else "paper_plots/histograms.root"
    )
    for ps, br in enumerate(["LeadPs_dR_gg", "SubleadPs_dR_gg"]):
        plotdRCombined(samples, br, ps+1, Scales=Scales),

    # Setup for lists
    branches, boundsList, binsScaleList, names, norms, titles = loadPlottingParameters()
    
    # Comparing 4 signal histograms together without bkg or data
    four_sig_hists = {}
    max_order_sig = 4
    #gw = [sum(samples[f"Signal_{mass}_GeV"].weight) for mass in [15, 30, 40, 60]]
    gw = [Scales.sumw[args.year if "2024" not in args.year else year_adjust][f"{mass}_GeV"] for mass in [15, 30, 40, 60]]

    sig_paths = [
        f"./plots/standard/{args.year}/compare/pho1_mvaID.pdf",
        f"./plots/standard/{args.year}/compare/pho2_mvaID.pdf",
        f"./plots/standard/{args.year}/compare/pho3_mvaID.pdf",
        f"./plots/standard/{args.year}/compare/pho4_mvaID.pdf",
        f"./plots/standard/{args.year}/compare/LeadPs_pt.pdf",
        f"./plots/standard/{args.year}/compare/SubleadPs_pt.pdf",
        f"./plots/standard/{args.year}/compare/dR_aa_mass_gggg.pdf",
        f"./plots/standard/{args.year}/compare/LeadPs_interMass.pdf",
        f"./plots/standard/{args.year}/compare/SubleadPs_interMass.pdf",
        f"./plots/standard/{args.year}/compare/LeadPs_interMass_abs.pdf",
        f"./plots/standard/{args.year}/compare/SubleadPs_interMass_abs.pdf",
        f"./plots/standard/{args.year}/compare/Ps_massDiff.pdf",
        f"./plots/standard/{args.year}/compare/cos_ag.pdf",
    ]

    # Fill histograms
    gc.collect()
    samples = dict(sorted(samples.items()))
    compare_lists = {
        "branches": branches[:len(sig_paths)] + ["mass_gggg"],  # Remove extra branches from loop!
        "boundsList": boundsList[:len(sig_paths)] + [[110.0, 180.0]],
        "binsScaleList": binsScaleList[:len(sig_paths)] + [30/70],
        "names": names[:len(sig_paths)] + ["mass_gggg"],
        "norms": norms[:len(sig_paths)] + [False],
        "titles": titles[:len(sig_paths)] + ["m_{#gamma#gamma#gamma#gamma}"],
        "sig_paths": sig_paths + [f"./plots/standard/{args.year}/compare/mass_gggg.pdf"],
    }
    for branch, bounds, binsScale, name, norm, title, path in zip(compare_lists["branches"], compare_lists["boundsList"], compare_lists["binsScaleList"], compare_lists["names"], compare_lists["norms"], compare_lists["titles"], compare_lists["sig_paths"]):
        for key, sample in samples.items():
            for m in ["Signal_15_GeV", "Signal_30_GeV", "Signal_40_GeV", "Signal_60_GeV"]:
                if m in key:
                    four_sig_hists.update({m: fillHist(branch, sample, bounds, binsScale=binsScale, name=f"{name}_{key}", normalize=norm, sb = False)})
                    four_sig_hists[m].Sumw2()
                    
        # Requires four_sig_hists to be set properly first
        if branch == "mass_gggg":
            print("Integrals:")
            for m in four_sig_hists:
                print(m, len(SideBand(samples[m])), len(SignalRegion(samples[m])), len(samples[m]))

        savePlots(
            "fourMassesSigOnly",
            plotFourMassesCompared(copy.deepcopy(four_sig_hists), title.replace(' [GeV]'," [GeV]"), path, gw, log=True, maximum=10**max_order_sig - 5*10**(max_order_sig-1), canvas_num=branch+"full_search_region", ignore_scales=False, year=args.year, Scales=Scales, sub="Full search region"),
            [f"{branch}_{m}" for m in ["Signal_15_GeV", "Signal_30_GeV", "Signal_40_GeV", "Signal_60_GeV"]],
            "paper_plots/histograms_dataSB.root" if not args.fullRegion else "paper_plots/histograms.root"
        )

    # Stop if requested to only do signal plots
    if args.sigOnly:
        return

    # Load samples
    samples_mHyp = {}
    if args.background is not None and args.data is not None:
        for key, inputs in zip(["background", "data"], [args.background, args.data]):
            print(f"Loading {key}...")
            assert os.path.exists(inputs)
            samples_mHyp.update({key: loadSamples(inputs, eras = eras if "outputs/BDT" not in inputs else None)})
            
            if args.cutBDT is not None:
                samples_mHyp[key] = samples_mHyp[key][samples_mHyp[key]["BDT_score"] > args.cutBDT]
    else:
        assert args.background is not None and args.data is not None, "Need to specify background and data samples!"

    savePlots(
        "mHyp",
        plotMHyp(samples, samples_mHyp, "m_hyp", "mHyp_EM_background", Scales=Scales),
        ["mHyp_EM_background", "mHyp_data"] + [f"mHyp_{m}_GeV" for m in range(15,65,5)],
        "paper_plots/histograms_dataSB.root" if not args.fullRegion else "paper_plots/histograms.root"
    )

    # Collect previous objects and disable circular garbage collector (ref counter still works)
    gc.disable()

    # Sample setup for per mass comparisons
    max_order = 6
    samples_general = samples_mHyp
    for canvas_num, mass in enumerate(samples.keys()):
        # Replace signal with new mass
        samples_general.update({f"signal": samples[mass]})
        mass = mass.replace("Signal_","")

        paths = [
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_mvaID_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_mvaID_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_mvaID_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_mvaID_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/LeadPs_pt_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/SubleadPs_pt_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/dR_aa_mass_gggg_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/LeadPs_interMass_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/SubleadPs_interMass_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/LeadPs_interMass_abs_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/SubleadPs_interMass_abs_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/Ps_massDiff_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/cos_ag_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/m_Hyp_{mass[:6]}.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pT1_ma1.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pT2_ma1.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pT1_ma2.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pT2_ma2.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_s4_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_s4_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_s4_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_s4_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_sieie_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_sieie_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_sieie_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_sieie_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_r9_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_r9_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_r9_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_r9_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_esEnergyOverRawE_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_esEnergyOverRawE_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_esEnergyOverRawE_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_esEnergyOverRawE_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_esEffSigmaRR_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_esEffSigmaRR_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_esEffSigmaRR_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_esEffSigmaRR_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_pfPhoIso03_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_pfPhoIso03_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_pfPhoIso03_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_pfPhoIso03_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_pfPhoIso03_rhoCorrected_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_pfPhoIso03_rhoCorrected_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_pfPhoIso03_rhoCorrected_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_pfPhoIso03_rhoCorrected_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_ecalPFClusterIso_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_ecalPFClusterIso_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_ecalPFClusterIso_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_ecalPFClusterIso_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_hcalPFClusterIso_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_hcalPFClusterIso_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_hcalPFClusterIso_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_hcalPFClusterIso_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_pfChargedIso_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_pfChargedIso_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_pfChargedIso_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_pfChargedIso_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_pfChargedIsoWorstVtx_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_pfChargedIsoWorstVtx_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_pfChargedIsoWorstVtx_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_pfChargedIsoWorstVtx_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_eta_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_eta_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_eta_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_eta_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_trkSumPtHollowConeDR03_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_trkSumPtHollowConeDR03_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_trkSumPtHollowConeDR03_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_trkSumPtHollowConeDR03_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_trkSumPtSolidConeDR04_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_trkSumPtSolidConeDR04_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_trkSumPtSolidConeDR04_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_trkSumPtSolidConeDR04_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_etaWidth_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_etaWidth_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_etaWidth_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_etaWidth_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho1_phiWidth_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho2_phiWidth_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho3_phiWidth_compare.pdf",
            f"./plots/BDT/{args.year}/{mass[:6]}/pho4_phiWidth_compare.pdf",
        ]

        if args.bdtOnly:
            keep = 14
            branches = branches[:keep]
            boundsList = boundsList[:keep]
            binsScaleList = binsScaleList[:keep]
            names = names[:keep]
            norms = norms[:keep]
            titles = titles[:keep]
            paths = paths[:keep]

        assert len(branches) == len(boundsList) == len(binsScaleList) == len(names) == len(norms) == len(titles) == len(paths), f"{len(branches)} {len(boundsList)} {len(binsScaleList)} {len(names)} {len(norms)} {len(titles)} {len(paths)}"

    # Make scatter plot of photon 3/4 MVA IDs
    mva_samples = samples_mHyp.copy()
    mva_samples.pop("signal")
    for canvas_num, mass in enumerate(samples.keys()):
        # Replace signal with new masses
        if "30_GeV" in mass:
            mva_samples.update({m: samples[mass]})
    make_scatter_mva(mva_samples, args.year)

    # Sample setup for general comparisons
    samples_general = samples_mHyp
    samples_general.pop("signal")
    for canvas_num, mass in enumerate(samples.keys()):
        samples_general.update({mass: samples[mass]})

    paths = [
        f"./plots/BDT/{args.year}/compare/pho1_mvaID_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_mvaID_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_mvaID_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_mvaID_compare.pdf",
        f"./plots/BDT/{args.year}/compare/LeadPs_pt_compare.pdf",
        f"./plots/BDT/{args.year}/compare/SubleadPs_pt_compare.pdf",
        f"./plots/BDT/{args.year}/compare/dR_aa_mass_gggg_compare.pdf",
        f"./plots/BDT/{args.year}/compare/LeadPs_interMass_compare.pdf",
        f"./plots/BDT/{args.year}/compare/SubleadPs_interMass_compare.pdf",
        f"./plots/BDT/{args.year}/compare/LeadPs_interMass_abs_compare.pdf",
        f"./plots/BDT/{args.year}/compare/SubleadPs_interMass_abs_compare.pdf",
        f"./plots/BDT/{args.year}/compare/Ps_massDiff_compare.pdf",
        f"./plots/BDT/{args.year}/compare/cos_ag_compare.pdf",
        f"./plots/BDT/{args.year}/compare/mHyp.pdf",
        f"./plots/BDT/{args.year}/compare/mass_gggg.pdf",
        f"./plots/BDT/{args.year}/compare/pT1_ma1_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pT2_ma1_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pT1_ma2_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pT2_ma2_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_s4_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_s4_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_s4_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_s4_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_sieie_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_sieie_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_sieie_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_sieie_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_r9_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_r9_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_r9_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_r9_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_esEnergyOverRawE_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_esEnergyOverRawE_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_esEnergyOverRawE_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_esEnergyOverRawE_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_esEffSigmaRR_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_esEffSigmaRR_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_esEffSigmaRR_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_esEffSigmaRR_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_pfPhoIso03_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_pfPhoIso03_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_pfPhoIso03_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_pfPhoIso03_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_pfPhoIso03_rhoCorrected_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_pfPhoIso03_rhoCorrected_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_pfPhoIso03_rhoCorrected_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_pfPhoIso03_rhoCorrected_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_ecalPFClusterIso_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_ecalPFClusterIso_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_ecalPFClusterIso_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_ecalPFClusterIso_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_hcalPFClusterIso_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_hcalPFClusterIso_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_hcalPFClusterIso_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_hcalPFClusterIso_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_pfChargedIso_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_pfChargedIso_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_pfChargedIso_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_pfChargedIso_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_pfChargedIsoWorstVtx_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_pfChargedIsoWorstVtx_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_pfChargedIsoWorstVtx_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_pfChargedIsoWorstVtx_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_eta_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_eta_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_eta_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_eta_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_trkSumPtHollowConeDR03_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_trkSumPtHollowConeDR03_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_trkSumPtHollowConeDR03_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_trkSumPtHollowConeDR03_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_trkSumPtSolidConeDR04_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_trkSumPtSolidConeDR04_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_trkSumPtSolidConeDR04_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_trkSumPtSolidConeDR04_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_etaWidth_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_etaWidth_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_etaWidth_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_etaWidth_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho1_phiWidth_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho2_phiWidth_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho3_phiWidth_compare.pdf",
        f"./plots/BDT/{args.year}/compare/pho4_phiWidth_compare.pdf",
    ]

    if args.bdtOnly:
        keep = 14
        branches = branches[:keep] + ["mass_gggg"]
        boundsList = boundsList[:keep] + [[110.0, 180.0]]
        binsScaleList = binsScaleList[:keep] + [30/70]
        names = names[:keep] + ["mass_gggg"]
        norms = norms[:keep] + [False]
        titles = titles[:keep] + ["m_{#gamma#gamma#gamma#gamma}"]
        paths = paths[:keep] + [f"./plots/BDT/{args.year}/compare/mass_gggg_compare.pdf"]

    assert len(branches) == len(boundsList) == len(binsScaleList) == len(names) == len(norms) == len(titles) == len(paths), f"{len(branches)} {len(boundsList)} {len(binsScaleList)} {len(names)} {len(norms)} {len(titles)} {len(paths)}"

    # Comment me out to remove BDT cuts!
    apply_bdt_cuts = args.applyBDTCuts
    if apply_bdt_cuts:
        bdt_2024_updatedPresel = {
            "Signal_15_GeV": 0.9600,
            "Signal_20_GeV": 0.9600,
            "Signal_25_GeV": 0.9600,
            "Signal_30_GeV": 0.9700,
            "Signal_35_GeV": 0.9650,
            "Signal_40_GeV": 0.9700,
            "Signal_45_GeV": 0.9750,
            "Signal_50_GeV": 0.9800,
            "Signal_55_GeV": 0.9800,
            "Signal_60_GeV": 0.9800,
        }
        for key in samples_general.keys():
            if not "BDT_score" in samples_general[key].columns:
                print(f"No column for BDT score found in {key}. Continuing without cuts.")
                continue
            else:
                if "GeV" in key:
                    print(f"Applying BDT cuts to {key}!", f"before: {Scales.lumi_fb * 1.0 * sum(samples_general[key].weight) / Scales.sumw[args.year.replace('_updatedPresel','')][key.replace('Signal_','')]}", len(samples_general[key]), end='\t')
                    print(f"cut {bdt_2024_updatedPresel[key]}", end='\t')
                    samples_general[key] = samples_general[key][samples_general[key]["BDT_score"] > bdt_2024_updatedPresel[key]]
                    print(f"after: {Scales.lumi_fb * 1.0 * sum(samples_general[key].weight) / Scales.sumw[args.year.replace('_updatedPresel','')][key.replace('Signal_','')]}", len(samples_general[key]))
                else:
                    print(f"Applying BDT cuts to {key}!", f"before: {sum(samples_general[key].weight)}", len(samples_general[key]), end='\t')
                    samples_general[key] = samples_general[key][samples_general[key]["BDT_score"] > bdt_2024_updatedPresel["Signal_30_GeV"]]
                    print(f"after: {sum(samples_general[key].weight)}", len(samples_general[key]))

    # Plot reweighted versions in sideband
    # Fill hists
    # Comparing only four maases
    for branch, bounds, binsScale, name, norm, title, path in zip(branches, boundsList, binsScaleList, names, norms, titles, paths):
        hists = {}
        hists_rw = {}
        for key, sample in samples_general.items():
            #for m in ["Signal_15_GeV", "Signal_30_GeV", "Signal_40_GeV", "Signal_60_GeV"] if not args.useFive else ["Signal_15_GeV", "Signal_25_GeV", "Signal_35_GeV", "Signal_45_GeV", "Signal_60_GeV"]:
            for m in [f"Signal_{x}_GeV" for x in range(15,65,5)]:
                if "mHyp" not in path and m not in key:
                    continue
                elif "mHyp" in path:
                    pass
            sample_df = sample
            print(f"key: {key}\t # events in SB: {len(sample[((sample.mass_gggg >= 110.0) & (sample.mass_gggg <= 115.0)) | ((sample.mass_gggg >= 135.0) & (sample.mass_gggg <= 180.0))])}\t # events in SR: {len(sample[(sample.mass_gggg > 115.0) & (sample.mass_gggg < 135.0)])}")
            if "background" in key or "data" in key:
                continue

            if "mHyp" in path:
                print(key, "mHyp")
            # BDT CUT
            if apply_bdt_cuts:
                sample_df = sample_df[sample_df["BDT_score"] > bdt_2024_updatedPresel[key]]
            tmp = fillHist(branch, sample_df, bounds, binsScale=binsScale, name=f"{name}_{key}", normalize=norm, sb = False)
            hists.update({key: tmp})

        # Normalizes histogram by lumi * xs(1fb) * BR(1.0) * efficiency (after-sel/pre-sel genWeight) since histograms are already scaled to after-sel weight
        # NOTE: weight column is based on genWeight with adjustments from corrections and/or systematics
        #gw = [sum(samples_general[f"Signal_{mass}_GeV"].weight) for mass in [15, 30, 40, 60]]
        gw = [Scales.sumw[args.year if "2024" not in args.year else year_adjust][f"{mass}_GeV"] for mass in range(15,65,5)]
        if args.useFive:
            gw5 = [sum(samples_general[f"Signal_{mass}_GeV"].weight) for mass in [15, 25, 35, 45, 60]]

        # Plot reweighted versions in sideband
        # Fill hists
        # BDT CUT
        if apply_bdt_cuts:
            samples_general["background"] = samples_general["background"][samples_general["background"]["BDT_score"] > bdt_2024_updatedPresel["Signal_30_GeV"]]
            samples_general["data"] = samples_general["data"][samples_general["data"]["BDT_score"] > bdt_2024_updatedPresel["Signal_30_GeV"]]
        bkg_hist = fillHist(branch, samples_general["background"], bounds, binsScale=binsScale, name=f"{name}_background", normalize=norm, sb = False)
        data_hist = fillHist(branch, samples_general["data"], bounds, binsScale=binsScale, name=f"{name}_data", normalize=norm, sb = False if args.fullRegion else True)
        if not args.noRW:
            bkg_rw = fillHist(branch, samples_general["background"], bounds, binsScale=binsScale, name=f"{name}_background_rw", normalize=norm, sb = False if args.fullRegion else True, customWeights="Ndimreweight")

        bkg_hist.Sumw2()
        data_hist.Sumw2()
        if not args.noRW:
            bkg_rw.Sumw2()

        # Reorder hists dictionary to play nice with the colors
        if not args.useFive:
            hists_pre_rw = {
                "background": bkg_hist,
                "data": data_hist,
                "Signal_15_GeV": hists["Signal_15_GeV"],
                "Signal_20_GeV": hists["Signal_20_GeV"],
                "Signal_25_GeV": hists["Signal_25_GeV"],
                "Signal_30_GeV": hists["Signal_30_GeV"],
                "Signal_35_GeV": hists["Signal_35_GeV"],
                "Signal_40_GeV": hists["Signal_40_GeV"],
                "Signal_45_GeV": hists["Signal_45_GeV"],
                "Signal_50_GeV": hists["Signal_50_GeV"],
                "Signal_55_GeV": hists["Signal_55_GeV"],
                "Signal_60_GeV": hists["Signal_60_GeV"]
            }
            
            if "mHyp" in path:
                hists_pre_rw_fixed = {
                    "background": bkg_hist,
                    "data": data_hist
                }
                for m in range(15,65,5):
                    hists_pre_rw_fixed.update({f"Signal_{m}_GeV": hists[f"Signal_{m}_GeV"]})
                hists_pre_rw = hists_pre_rw_fixed

            if args.sigScaleBkg:
                for k in [k for k in hists_pre_rw.keys() if "Signal_" in k]:
                    hists_pre_rw[k].Scale(hists_pre_rw["data"].Integral() / hists_pre_rw[k].Integral())

        else:
            hists_pre_rw = {
                "background": bkg_hist,
                "data": data_hist,
                "Signal_15_GeV": hists["Signal_15_GeV"],
                "Signal_25_GeV": hists["Signal_25_GeV"],
                "Signal_35_GeV": hists["Signal_35_GeV"],
                "Signal_45_GeV": hists["Signal_45_GeV"],
                "Signal_60_GeV": hists["Signal_60_GeV"]
            }

        if not args.noRW:
            if not args.useFive:
                hists_rw = {
                    "background": bkg_rw,
                    "data": data_hist,
                    "Signal_15_GeV": hists["Signal_15_GeV"],
                    "Signal_20_GeV": hists["Signal_20_GeV"],
                    "Signal_25_GeV": hists["Signal_25_GeV"],
                    "Signal_30_GeV": hists["Signal_30_GeV"],
                    "Signal_35_GeV": hists["Signal_35_GeV"],
                    "Signal_40_GeV": hists["Signal_40_GeV"],
                    "Signal_45_GeV": hists["Signal_45_GeV"],
                    "Signal_50_GeV": hists["Signal_50_GeV"],
                    "Signal_55_GeV": hists["Signal_55_GeV"],
                    "Signal_60_GeV": hists["Signal_60_GeV"]
                }

                if "mHyp" in path:
                    hists_rw_fixed = {
                        "background": bkg_hist,
                        "data": data_hist
                    }
                    for m in range(15,65,5):
                        hists_rw_fixed.update({f"Signal_{m}_GeV": hists[f"Signal_{m}_GeV"]})
                    hists_rw = hists_rw_fixed

            else:
                hists_rw = {
                    "background": bkg_rw,
                    "data": data_hist,
                    "Signal_15_GeV": hists["Signal_15_GeV"],
                    "Signal_25_GeV": hists["Signal_25_GeV"],
                    "Signal_35_GeV": hists["Signal_35_GeV"],
                    "Signal_45_GeV": hists["Signal_45_GeV"],
                    "Signal_60_GeV": hists["Signal_60_GeV"]
                }

            if args.sigScaleBkg:
                for k in [k for k in hists_rw.keys() if "Signal_" in k]:
                    hists_rw[k].Scale(hists_rw["data"].Integral() / hists_rw[k].Integral())

        if "mHyp" in path:
            #gw = [sum(samples_general[f"Signal_{mass}_GeV"].weight) for mass in range(15,65,5)]
            gw = [Scales.sumw[args.year if "2024" not in args.year else year_adjust][f"{mass}_GeV"] for mass in range(15,65,5)]
        else:
            #gw = [sum(samples_general[f"Signal_{mass}_GeV"].weight) for mass in [15, 30, 40, 60]]
            gw = [Scales.sumw[args.year if "2024" not in args.year else year_adjust][f"{mass}_GeV"] for mass in range(15,65,5)]
            if args.useFive:
                gw5 = [sum(samples_general[f"Signal_{mass}_GeV"].weight) for mass in [15, 25, 35, 45, 60]]

        savePlots(
            "compare_pre_reweight",
            compare_four_masses_data_sig_bkg(copy.deepcopy(hists_pre_rw), title.replace(' [GeV]'," [GeV]"), path[:-4]+"_pre_reweight.pdf", gw if not args.useFive else gw5, log=True, maximum=(10**max_order - 5*10**(max_order-1)), canvas_num=branch+"entire_region_prerw", ignore_scales=False, year=args.year, Scales=Scales, use_five=args.useFive, sub="Data blinded" if not args.fullRegion else None, skipSigScale=args.sigScaleBkg),
            [f"{branch}_{k}" for k in list(hists_pre_rw.keys())],
            "paper_plots/histograms_dataSB.root" if not args.fullRegion else "paper_plots/histograms.root"
        )

        """
        if "mHyp" not in path and "abs" not in path:
            compare_sigMixing_bkg(copy.deepcopy(hists_pre_rw), copy.deepcopy(hists_pre_rw["data"]), title.replace(' [GeV]'," [GeV]"), (path[:-4]+"_pre_reweight.pdf").replace("/compare/", "/sig_mix/"), gw, log=True, maximum=(10**max_order - 5*10**(max_order-1)), canvas_num=branch+"entire_region_prerw", ignore_scales=False, year=args.year, Scales=Scales, sub="Data blinded" if not args.fullRegion else None, skipSigScale=args.sigScaleBkg)
        """

        new_title = title + " with Reweighting"
        if "[GeV]" in new_title:
            new_title = new_title.replace("[GeV]","") + " [GeV]"
        if "[GeV^{-1}]" in new_title:
            new_title = new_title.replace("[GeV^{-1}]","") + " [GeV^{-1}]"
        if not args.noRW:
            savePlots(
                "compare_reweight",
                compare_four_masses_data_sig_bkg(copy.deepcopy(hists_rw), new_title, path[:-4]+"_reweight.pdf", gw if not args.useFive else gw5, log=True, maximum=(10**max_order - 5*10**(max_order-1)) , canvas_num=branch+"sb_region_rw", ignore_scales=False, year=args.year, Scales=Scales, use_five=args.useFive, sub="Data blinded" if not args.fullRegion else None, skipSigScale=args.sigScaleBkg),
                [f"{branch}_{k}" for k in list(hists_rw.keys())],
                "paper_plots/histograms_dataSB.root" if not args.fullRegion else "paper_plots/histograms.root"
            )

            """
            if "mHyp" not in path and "abs" not in path:
                print("keys:", hists_rw.keys())
                compare_sigMixing_bkg(copy.deepcopy(hists_rw), copy.deepcopy(hists_rw["data"]), new_title, (path[:-4]+"_reweight.pdf").replace("/compare/", "/sig_mix/"), gw, log=True, maximum=(10**max_order - 5*10**(max_order-1)), canvas_num=branch+"sb_region_rw", ignore_scales=False, year=args.year, Scales=Scales, sub="Data blinded" if not args.fullRegion else None, skipSigScale=args.sigScaleBkg)
            """


if __name__ == "__main__":
    parser = get_parser_plotting()
    parser.add_argument("-y", "--year", type=str, required=True, help="Year of provided samples.")
    parser.add_argument("-sigOnly", "--sigOnly", default=False, action="store_true", help="Runs only signal plots.")
    parser.add_argument("-five", "--useFive", default=False, action="store_true", help="Plot five masses instead for BDT inputs.")
    parser.add_argument("-full", "--fullRegion", default=False, action="store_true", help="Plot full samples. Defaults to sideband only.")
    parser.add_argument("--sigScaleBkg", default=False, action="store_true", help="Scale signal hists to bkg for direct comparisons.")
    parser.add_argument("--applyBDTCuts", default=False, action="store_true", help="Apply BDT cuts before plotting. Checks values!")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-bdt", "--bdtOnly", default=False, action="store_true", help="Restricts plotting to BDT input variables only.")
    args = parser.parse_args()

    generatePlots(args)
