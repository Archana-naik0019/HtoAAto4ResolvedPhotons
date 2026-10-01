import os
from ROOT import TCanvas, gStyle, TLegend, kAzure, kGray, kMagenta, kOrange, kRed, kBlue
from h4g_tools.bdt.reweight import reweight1Dim
from h4g_tools.utils.loading import pathSetup, loadPerMass, loadSamples, fillHist
from h4g_tools.utils.plotting import plot_corr_matrix, loadPlottingParameters

if __name__ == "__main__":

    # 2022
    year = "2022"
    xgb_name = "ES_BDT_2022"
    pathSetup(year)

    old_branches, _, _, _, _, _ = loadPlottingParameters(bdt=True)
    branches = []
    for br in old_branches:
        if "interMass" not in br:
            branches.append(br)


    sig = loadSamples("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2022_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSF_TriggerSF_oldXqcutQcut_20Dec2025/", branches=branches)
    bkg = loadSamples("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_bkg2022_SaS_05Aug2025/", branches=branches)
    
    results = {}
    plot_corr_matrix(results, f"plots/BDT/{year}/{xgb_name}/corr_matrix_{year}.png", full=True, sig=sig, bkg=bkg)

    # 2024
    year = "2024"
    xgb_name = "ES_BDT_2024"
    pathSetup(year)

    sig = loadSamples("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_wBDT_genMatching_06Feb2026", branches=branches)
    bkg = loadSamples("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_bkg2024_SaS_01Dec2025/", branches=branches)
    
    results = {}
    plot_corr_matrix(results, f"plots/BDT/{year}/{xgb_name}/corr_matrix_{year}.png", full=True, sig=sig, bkg=bkg)
