import os, xgboost, array
from typing import List
import pandas as pd, numpy as np
from ROOT import TH1D, TGraph, SetOwnership, TCanvas, gPad, kRed, kBlue, TLegend, TList
import h4g_tools.bdt.event_selection_bdt as event_selection_bdt
import h4g_tools.utils.loading as loading

def SideBand(df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns DataFrame for only the sideband region.
    """
    return df[(df.mass_gggg >= 110.0) & (df.mass_gggg < 115.0) | (df.mass_gggg > 135.0) & (df.mass_gggg <= 180.0)]


def SignalRegion(df: pd.DataFrame, window=None) -> pd.DataFrame:
    """
    Returns DataFrame for only the signal region.
    """

    return df[(df.mass_gggg >= 115.0) & (df.mass_gggg <= 135.0)] if window is None else df[(df.mass_gggg >= window[0]) & (df.mass_gggg <= window[1])]


def all4EB(df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns DataFrame for only events with all 4 photons in the barrel.
    """
    print("4 photons in EB: ", len(df[df.pho1_isScEtaEB & df.pho2_isScEtaEB & df.pho3_isScEtaEB & df.pho4_isScEtaEB]))
    return df[df.pho1_isScEtaEB & df.pho2_isScEtaEB & df.pho3_isScEtaEB & df.pho4_isScEtaEB]


class Scales:
    __doc__ = "\n    Stores scales for cross-section and luminosity for H->aa->gggg analysis.\n    "

    def __init__(self, _yr: str):
        years = {
            '2018': {'lumiTxt': '    62.5 fb^{-1} (13 TeV)', 'lumi': 62471.907487772, 'ggH_xs': 48.61, 'defaultLumiTxt': '     13 TeV (2018)'},
            '2022': {'lumiTxt': '   35 fb^{-1} (13.6 TeV)', 'lumi': 34748.6, 'ggH_xs': 51.96, 'defaultLumiTxt': '   13.6 TeV (2022)'},
            '2022preEE': {'lumiTxt': '    7.98 fb^{-1} (13.6 TeV)', 'lumi': 7980.4, 'ggH_xs': 51.96, 'defaultLumiTxt': '   13.6 TeV (2022)'},
            '2022postEE': {'lumiTxt': '   26.67 fb^{-1} (13.6 TeV)', 'lumi': 26671.7, 'ggH_xs': 51.96, 'defaultLumiTxt': ' 13.6 TeV (2022EE)'},
            '2024': {'lumiTxt': '    110 fb^{-1} (13.6 TeV)', 'lumi': 109820.0, 'ggH_xs': 51.96, 'defaultLumiTxt': '   13.6 TeV (2024)'},
            '2024_mixed': {'lumiTxt': '    110 fb^{-1} (13.6 TeV)', 'lumi': 109820.0, 'ggH_xs': 51.96, 'defaultLumiTxt': '   13.6 TeV (2024)'},
            # ADD THE LUMI INFO LIKE ABOVE FOR EACH NEW YEAR
        }

        if "2024" in _yr:
            yr = "2024"
        else:
            yr = _yr
        self.defaultLumiTxt = years[yr]["defaultLumiTxt"]
        self.lumiTxt = years[yr]["lumiTxt"]
        self.lumi = years[yr]["lumi"]
        self.lumi_fb = years[yr]["lumi"] / 1000.0

        self.sumw = {"2018": (dict(zip([f"{m}_GeV" for m in range(15, 65, 5)], [494609, 495244, 495688, 495872, 495530, 494765, 494511, 495205, 494471, 495597])))}
        if yr == "2022preEE":
            self.sumw.update({"2022": (dict(zip([f"{m}_GeV" for m in range(15, 65, 5)], [608224, 604425, 591575, 605230, 594887, 596384, 595630, 597008, 592309, 601131])))})
        elif yr == "2022postEE":
            self.sumw.update({"2022": (dict(zip([f"{m}_GeV" for m in range(15, 65, 5)], [894716, 896493, 894327, 902834, 888381, 895125, 898545, 899078, 893452, 907428])))})
        elif yr == "2022":
            self.sumw.update({"2022": (dict(zip([f"{m}_GeV" for m in range(15, 65, 5)], [1502940, 1500918, 1485902, 1508064, 1483268, 1491509, 1494175, 1455585, 1485761, 1508559])))})
        elif "2024" in yr:
            self.sumw.update({"2024": (dict(zip([f"{m}_GeV" for m in range(15, 65, 5)], [1025691, 971967, 969906, 981050, 1047791, 981828, 974661, 1022301, 994105, 955078])))})

        # ADD THE SUM OF WEIGHTS INFO FOR EACH NEW YEAR HERE

        self.yr = yr

    useScaledCounts = True
    BR = {"2018": (dict(zip([f"{m}_GeV" for m in range(15, 65, 5)], [1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05])))}





