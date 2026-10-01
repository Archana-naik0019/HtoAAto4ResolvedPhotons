import os
import copy
from h4g_tools.utils.save_plots_to_root import savePlots
from h4g_tools.utils.scales import Scales as ScalesCls, SideBand, SignalRegion
from h4g_tools.utils.runner_utils import get_parser_plotting, handle_errors
from h4g_tools.utils.plotting import compare_four_masses_data_sig_bkg
from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist
from h4g_tools.bdt.reweight import loadNDimWeights

if __name__ == "__main__":
    signal = loadPerMass("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSFpostBPix_wBDT_wHLT_28Apr2026/", branches=["BDT_score", "mass_gggg"])
    bkg = loadSamples("/cms/cephfs/data/store/user/castells/outputs/BDT/2024/ES_BDT_2024/generic/background/", branches=["BDT_score", "mass_gggg", "Ndimreweight"])
    data = loadSamples("/cms/cephfs/data/store/user/castells/outputs/BDT/2024/ES_BDT_2024/generic/data/", branches=["BDT_score", "mass_gggg"])

    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Scales instance setup
    masses = [x for x in range(15, 65, 5)]
    eras={
        "data": ["Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"],
        "bkg": ["RunEvtMix2024C", "RunEvtMix2024D", "RunEvtMix2024E", "RunEvtMix2024F", "RunEvtMix2024G", "RunEvtMix2024H", "RunEvtMix2024I"],
        "signal": [f"Signal_{m}_GeV" for m in masses]
    }

    # Handle path setup after checking year
    Scales = ScalesCls("2024") 
    pathSetup("2024")

    bkg_scale = 1.0
    for BDT_cut in [0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98]:
        signal_copy = copy.deepcopy(signal)
        bkg_copy = copy.deepcopy(bkg)
        data_copy = copy.deepcopy(data)

        # Apply BDT cut
        for k in signal_copy.keys():
            signal_copy[k] = signal_copy[k][signal_copy[k].BDT_score > BDT_cut]
        bkg_copy = bkg_copy[bkg_copy.BDT_score > BDT_cut]
        data_copy = data_copy[data_copy.BDT_score > BDT_cut]

        signal_15_GeV = fillHist("mass_gggg", signal_copy["Signal_15_GeV"], [110.0, 180.0], binsScale=(30.0/70.0), name=f"mass_gggg_15_GeV_{BDT_cut}", normalize=False, sb = False)
        signal_30_GeV = fillHist("mass_gggg", signal_copy["Signal_30_GeV"], [110.0, 180.0], binsScale=(30.0/70.0), name=f"mass_gggg_30_GeV_{BDT_cut}", normalize=False, sb = False)
        signal_40_GeV = fillHist("mass_gggg", signal_copy["Signal_40_GeV"], [110.0, 180.0], binsScale=(30.0/70.0), name=f"mass_gggg_40_GeV_{BDT_cut}", normalize=False, sb = False)
        signal_60_GeV = fillHist("mass_gggg", signal_copy["Signal_60_GeV"], [110.0, 180.0], binsScale=(30.0/70.0), name=f"mass_gggg_60_GeV_{BDT_cut}", normalize=False, sb = False)
        bkg_unblind   = fillHist("mass_gggg", bkg_copy, [110.0, 180.0], binsScale=1.0, name=f"mass_gggg_bkg_{BDT_cut}", normalize=False, sb = False, customWeights="Ndimreweight")
        bkg_blind     = fillHist("mass_gggg", bkg_copy, [110.0, 180.0], binsScale=1.0, name=f"mass_gggg_bkg_{BDT_cut}", normalize=False, sb = True, customWeights="Ndimreweight")
        data_unblind  = fillHist("mass_gggg", data_copy, [110.0, 180.0], binsScale=1.0, name=f"mass_gggg_data_{BDT_cut}", normalize=False, sb = False)
        data_blind    = fillHist("mass_gggg", data_copy, [110.0, 180.0], binsScale=1.0, name=f"mass_gggg_data_{BDT_cut}", normalize=False, sb = True)

        # Scale bkg to blinded data
        sub = "Data blinded"
        if BDT_cut == 0.0:
            bkg_scale = data_unblind.Integral() / bkg_unblind.Integral()
            sub = ""
        print(f"BKG scaled by: {bkg_scale}")
        bkg_unblind.Scale(bkg_scale)
        print(f"BKG (unblind) integral: {bkg_unblind.Integral()}")
        print(f"DATA (unblind) integral: {data_unblind.Integral()}")
        print(f"DATA (blind) integral: {data_blind.Integral()}")
        
        hists_rw = {
            "background": bkg_unblind,
            "data": data_unblind if BDT_cut == 0.0 else data_blind,
            "Signal_15_GeV": signal_15_GeV,
            "Signal_30_GeV": signal_30_GeV,
            "Signal_40_GeV": signal_40_GeV,
            "Signal_60_GeV": signal_60_GeV
        }
        max_order = 6 
        gw = [Scales.sumw["2024"][f"{m}_GeV"] for m in [15, 30, 40, 60]]
        compare_four_masses_data_sig_bkg(hists_rw, "m_{#gamma#gamma#gamma#gamma}", f"./plots/BDT/2024/m4g_bdt_checks/mass_gggg_compare_reweight_{str(BDT_cut).replace('.','p')}.pdf", gw, log=True, maximum=(10**max_order - 5*10**(max_order-1)) , canvas_num=f"mass_gggg_sb_region_rw_{BDT_cut}", ignore_scales=False, year="2024", Scales=Scales, use_five=False, sub=sub, skipSigScale=False, skipBkgScale=True, noRatioScaling=True),

        compare_four_masses_data_sig_bkg(hists_rw, "m_{#gamma#gamma#gamma#gamma}", f"./plots/BDT/2024/m4g_bdt_checks/mass_gggg_compare_reweight_{str(BDT_cut).replace('.','p')}_linear.pdf", gw, log=False, canvas_num=f"mass_gggg_sb_region_rw_{BDT_cut}", ignore_scales=False, year="2024", Scales=Scales, use_five=False, sub=sub, skipSigScale=False, skipBkgScale=True, noRatioScaling=True),


