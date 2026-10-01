from h4g_tools.utils.loading import loadPerMass, fillHist
from ROOT import TCanvas, TLegend, kBlack, kBlue, kRed, kMagenta, kOrange


def setcolors(title, hists=[], colors=[], titles=None):
    for h,c in zip(hists, colors):
        h.SetTitle(title)
        h.SetLineColor(c)
        h.SetLineWidth(2)
        
        if titles is not None:
            h.GetXaxis().SetTitle(titles[0])
            h.GetYaxis().SetTitle(titles[1])
        

if __name__ == "__main__":
    #samples = loadPerMass("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSFpostBPix_updatedPresel_onlyMVA90SFs_wBDT_04Jun2024", branches=["weight", "weight_PhotonIdwp90SFUp", "weight_PhotonIdwp90SFDown", "pho1_pt", "pho2_pt", "pho3_pt", "pho4_pt", "pho1_ScEta", "pho2_ScEta", "pho3_ScEta", "pho4_ScEta"])
    samples = loadPerMass("/cms/cephfs/data/store/user/castells/_DO_NOT_DELETE_outputs_HiggsDNA/backup/outputs_signal2024_SaS_Pileup_NFs_Material_FNUF_eVetoSF_PreselSFpostBPix_TriggerEff_updatedPresel_wBDT_04Jun2026", branches=["pho1_pt", "pho2_pt", "pho3_pt", "pho4_pt", "pho1_ScEta", "pho2_ScEta", "pho3_ScEta", "pho4_ScEta", "weight"])
   
    for mass in samples.keys():
        # Check structures in regions of weights
        s = samples[mass]
        #s = samples[mass][(samples[mass].weight > 0.8) & (samples[mass].weight < 0.95)]
        #s = samples[mass][(samples[mass].weight > 1.0) & (samples[mass].weight < 1.1)]
        #s = samples[mass][(samples[mass].weight > 1.15) & (samples[mass].weight < 1.25)]

        #peak = None
        #peak = [0.8, 0.95]
        #peak = [1.0, 1.1]
        #peak = [1.15, 1.25]
        peak = [0.6, 0.85]
        if peak is not None:
            print("!!! Applying cut on peak !!!")
            s = samples[mass][(samples[mass].weight > peak[0]) & (samples[mass].weight < peak[1])]

        # Fill pT histograms (put all together in one plot)
        hist_pt1 = fillHist("pho1_pt", s, [0.0, 250.0], binsScale=0.24, normalize=False)
        hist_pt2 = fillHist("pho2_pt", s, [0.0, 250.0], binsScale=0.24, normalize=False)
        hist_pt3 = fillHist("pho3_pt", s, [0.0, 250.0], binsScale=0.24, normalize=False)
        hist_pt4 = fillHist("pho4_pt", s, [0.0, 250.0], binsScale=0.24, normalize=False)

        # Fill eta histograms (put all together in one plot)
        hist_eta1 = fillHist("pho1_ScEta", s, [-2.5, 2.5], binsScale=6, normalize=False)
        hist_eta2 = fillHist("pho2_ScEta", s, [-2.5, 2.5], binsScale=6, normalize=False)
        hist_eta3 = fillHist("pho3_ScEta", s, [-2.5, 2.5], binsScale=6, normalize=False)
        hist_eta4 = fillHist("pho4_ScEta", s, [-2.5, 2.5], binsScale=6, normalize=False)

        """
        hist_nom = fillHist("weight", s, [0.5, 2.0], binsScale=50, normalize=False, includeStats=True)
        hist_up = fillHist("weight_PhotonIdwp90SFUp", s, [0.5, 2.0], binsScale=50, normalize=False, includeStats=True)
        hist_down = fillHist("weight_PhotonIdwp90SFDown", s, [0.5, 2.0], binsScale=50, normalize=False, includeStats=True)

        setcolors(
            f"m_{{a}} = {mass.replace('Signal_','').replace('_',' ')}",
            [hist_nom, hist_up, hist_down],
            [kBlack, kBlue, kRed]
        )
        """

        c = TCanvas("c","c", 600,600)
        c.SetLeftMargin(0.14)
        c.SetRightMargin(0.05)
        c.SetBottomMargin(0.14)

        #legend = TLegend(0.35,0.6,0.75,0.85)
        legend = TLegend(0.65,0.7,0.85,0.85)
        legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)

        """
        legend.AddEntry(hist_up, f"Up (mean = {hist_up.GetMean(): .3f})", "l")
        legend.AddEntry(hist_nom, f"Nominal (mean = {hist_nom.GetMean(): .3f})", "l")
        legend.AddEntry(hist_down, f"Down (mean = {hist_down.GetMean(): .3f})", "l")

        hist_nom.Draw("hist")
        hist_up.Draw("hist same")
        hist_down.Draw("hist same")
        legend.Draw()

        c.SaveAs(f"plots/BDT/2024_updatedPresel/{mass.replace('Signal_','')}/weights_wp90_{mass.replace('Signal_','')}.pdf")
        """

        # Do pT plotting
        legend.Clear()
        c.SetLogy()

        setcolors(
            f"m_{{a}} = {mass.replace('Signal_','').replace('_',' ')}",
            [hist_pt1, hist_pt2, hist_pt3, hist_pt4],
            [kMagenta, kOrange, kRed, kBlue],
            ["p_{T}", "Events"]
        )

        legend.AddEntry(hist_pt1, "#gamma_{1} p_{T}", "l")
        legend.AddEntry(hist_pt2, "#gamma_{2} p_{T}", "l")
        legend.AddEntry(hist_pt3, "#gamma_{3} p_{T}", "l")
        legend.AddEntry(hist_pt4, "#gamma_{4} p_{T}", "l")

        #hist_pt1.GetYaxis().SetRangeUser(10**-1 - 5*10**-2, 10**5 - 5*10**4)
        #hist_pt1.GetYaxis().SetRangeUser(10**-1 - 5*10**-2, 10**5 - 5*10**4)
        hist_pt1.GetYaxis().SetRangeUser(10**-1 - 5*10**-2, 10**5 - 5*10**4)

        hist_pt1.Draw("hist")
        hist_pt2.Draw("hist same")
        hist_pt3.Draw("hist same")
        hist_pt4.Draw("hist same")
        legend.Draw()

        c.SaveAs(f"plots/BDT/2024_updatedPresel/{mass.replace('Signal_','')}/pt_wp90SubRegion_{mass.replace('Signal_','')}.pdf")

        # Do eta plotting
        legend.Clear()
        c.SetLogy(False)

        setcolors(
            f"m_{{a}} = {mass.replace('Signal_','').replace('_',' ')}",
            [hist_eta1, hist_eta2, hist_eta3, hist_eta4],
            [kMagenta, kOrange, kRed, kBlue],
            ["#eta", "Events"]
        )

        legend.AddEntry(hist_eta1, "#gamma_{1} p_{T}", "l")
        legend.AddEntry(hist_eta2, "#gamma_{2} p_{T}", "l")
        legend.AddEntry(hist_eta3, "#gamma_{3} p_{T}", "l")
        legend.AddEntry(hist_eta4, "#gamma_{4} p_{T}", "l")

        #hist_eta1.GetYaxis().SetRangeUser(0.0, 2000)
        #hist_eta1.GetYaxis().SetRangeUser(0.0, 1500)
        hist_eta1.GetYaxis().SetRangeUser(0.0, 800)
        #hist_eta1.GetYaxis().SetRangeUser(0.0, 2300)

        hist_eta1.Draw("hist")
        hist_eta2.Draw("hist same")
        hist_eta3.Draw("hist same")
        hist_eta4.Draw("hist same")
        legend.Draw()

        c.SaveAs(f"plots/BDT/2024_updatedPresel/{mass.replace('Signal_','')}/eta_wp90SubRegion_{mass.replace('Signal_','')}.pdf")

