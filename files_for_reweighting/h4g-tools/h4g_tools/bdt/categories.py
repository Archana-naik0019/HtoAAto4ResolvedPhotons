from ROOT import TGraph, TH1D, TH1
import numpy as np
from math import sqrt, log, fsum
import pandas as pd

def defineAMS(
    S: float,
    B: float,
    B_sb: float,
    D_sb: float,
    return_all: bool
) -> float:
    """
    Defines AMS with float inputs (in signal region).
    """

    if B < 1e-4 or S < 1e-4:
        return -1
    
    log_e = log(1.0 + S/B)
    expr = 2.0 * ((S + B) * log_e - S)
    sq = sqrt(expr)

    min_evt_sb = 10
    #if sq > 0.0 and B_sb >= min_evt_sb and D_sb >= min_evt_sb and B > 0:
    if sq > 0.0 and D_sb >= min_evt_sb and B > 0 and not return_all:
        return sq
    if sq > 0.0 and B > 0 and return_all:
        return -999 if D_sb < min_evt_sb else sq, sq 
    else:
        return -999


def sumSignificance(
    partition: list,
    hist_sig_SR: TH1D,
    hist_bkg_SR: TH1D,
    #hist_bkg_SB: TH1D,
    #hist_data_SB: TH1D,
    B_sb_df: pd.DataFrame,
    D_sb_df: pd.DataFrame,
    return_all: bool = False,
) -> float:
    """
    Finds best significance for an arbitary number of categories.
    """

    sum_significance = 0.0
    all_significance = 0.0

    for pair in partition:
        S = hist_sig_SR.Integral(pair[0], pair[1])
        #B_sb = hist_bkg_SB.Integral(pair[0], pair[1])
        #D_sb = hist_data_SB.Integral(pair[0], pair[1])

        pair_translated = [hist_sig_SR.GetXaxis().GetBinLowEdge(pair[0]), hist_sig_SR.GetXaxis().GetBinLowEdge(pair[1])]
        if isinstance(hist_bkg_SR, TH1):
            B = hist_bkg_SR.Integral(pair[0], pair[1])
        elif type(hist_bkg_SR) is pd.DataFrame:
            B = sum(hist_bkg_SR[(hist_bkg_SR.BDT_score > pair_translated[0]) & (hist_bkg_SR.BDT_score <= pair_translated[1])]["1dimreweight"])
        B_sb = sum(B_sb_df[(B_sb_df.BDT_score > pair_translated[0]) & (B_sb_df.BDT_score <= pair_translated[1])]["1dimreweight"])
        D_sb = len(D_sb_df[(D_sb_df.BDT_score > pair_translated[0]) & (D_sb_df.BDT_score <= pair_translated[1])])

        significance = defineAMS(S, B, B_sb, D_sb, return_all)
        #print(f"BDT cut: {hist_sig_SR.GetXaxis().GetBinLowEdge(pair[0]):<5.4}", f"AMS: {significance:<7.4f}", f"S: {S:<4.2f}", f"B: {float(B):<4.2f}", f"D_sb: {float(D_sb):<4.2f}")
        #print(f"AMS: {significance:<7.4f}", f"S: {S:<4.2f}", f"B: {float(B_sb):<4.2f}", f"D_sb: {float(D_sb):<4.2f}", [hist_sig_SR.GetBinCenter(pair[0])-hist_sig_SR.GetBinWidth(pair[0])/2, hist_sig_SR.GetBinCenter(pair[1])-hist_sig_SR.GetBinWidth(pair[1])/2])

        if (significance if not return_all else significance[1]) > 0.0:
            if not return_all:
                sum_significance += significance * significance
            else:
                if significance[0] > 0.0:
                    sum_significance += significance[0] * significance[0]
                else:
                    sum_significance = -999
                all_significance += significance[1] * significance[1]
        elif not return_all:
            return -999

    if return_all:
        return sqrt(sum_significance) if sum_significance > 0.0 else -999, sqrt(all_significance), S, B, D_sb
    else:
        return sqrt(sum_significance)

