from h4g_tools.utils.loading import loadPerMass, fillHist
from ROOT import TCanvas, TLegend, kBlack, kBlue, kRed

if __name__ == "__main__":
    path = "signal-sample-path"
    samples = loadPerMass(path)
   
    for mass in samples.keys():
        # Split based on EB/EE for various photons!!!
        print(len(samples[mass]), end='\t')
        eta = None
        if eta == "allEB":
            samples[mass] = samples[mass][(samples[mass].pho1_isScEtaEB == True) & (samples[mass].pho2_isScEtaEB == True) & (samples[mass].pho3_isScEtaEB == True) & (samples[mass].pho4_isScEtaEB == True)]
        elif eta == "allEE":
            samples[mass] = samples[mass][(samples[mass].pho1_isScEtaEE == True) & (samples[mass].pho2_isScEtaEE == True) & (samples[mass].pho3_isScEtaEE == True) & (samples[mass].pho4_isScEtaEE == True)]
        elif eta == "oneEB":
            samples[mass] = samples[mass][(samples[mass].pho1_isScEtaEB == True) | (samples[mass].pho2_isScEtaEB == True) | (samples[mass].pho3_isScEtaEB == True) | (samples[mass].pho4_isScEtaEB == True)]
        elif eta == "oneEE":
            samples[mass] = samples[mass][(samples[mass].pho1_isScEtaEE == True) | (samples[mass].pho2_isScEtaEE == True) | (samples[mass].pho3_isScEtaEE == True) | (samples[mass].pho4_isScEtaEE == True)]

        """
        #peak = [0.6, 0.85]
        peak = None
        if peak is not None:
            samples[mass] = samples[mass][(samples[mass].weight > peak[0]) & (samples[mass].weight < peak[1])]
        """
        print(len(samples[mass]))
        hist_nom = fillHist("weight", samples[mass], [0.0, 3.0], binsScale=50, normalize=False, includeStats=True)
        hist_nom.SetTitle(f"m_{{a}} = {mass.replace('Signal_','').replace('_',' ')}")

        hist_nom.SetLineColor(kBlack)
        hist_nom.SetLineWidth(2)

        c = TCanvas("c","c", 600,600)
        legend = TLegend(0.35,0.6,0.75,0.85)
        legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)

        legend.AddEntry(hist_nom, f"Nominal (mean = {hist_nom.GetMean(): .3f})", "l")

        hist_nom.Draw("hist")
        legend.Draw()

        c.SaveAs(f"plots/BDT/2024_updatedPresel/{mass.replace('Signal_','')}/weights_{mass.replace('Signal_','')}.pdf")

