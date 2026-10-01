import os
import pandas as pd
from ROOT import TH1D, TCanvas, gROOT, TLegend, kGreen, kBlue, THStack
from h4g_tools.utils.loading import fillHist


# Set batch so no plots pop up
gROOT.SetBatch(True)

if __name__ == "__main__":
    # Get current directory
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Load samples from bdtIO
    preEE = pd.read_parquet(f"/cms/cephfs/data/store/user/castells/bdtIO/BDT_2022preEE_input_weights.parquet")
    postEE = pd.read_parquet(f"/cms/cephfs/data/store/user/castells/bdtIO/BDT_2022postEE_input_weights.parquet")

    # Set bounds
    low = 0.0
    high = 0.5
    #low = 0.9
    #high = 1.4
    binsScale = 500  # Default of 1000 scale means 500 bins

    # Make a histogram with values of weights and plot it
    stack_preEE = THStack("stack_preEE", "Ndimreweight")
    stack_postEE = THStack("stack_postEE", "Ndimreweight")
    hist_preEE = fillHist("Ndimreweight", preEE, [low, high], binsScale=binsScale, name="preEE", normalize=False)
    hist_postEE = fillHist("Ndimreweight", postEE, [low, high], binsScale=binsScale, name="postEE", normalize=False)
    plotSB = True
    if plotSB:
        preEE = preEE[((preEE.mass_gggg >= 110.0) & (preEE.mass_gggg <= 115.0)) | ((preEE.mass_gggg >= 135.0) & (preEE.mass_gggg <= 180.0))]
        postEE = postEE[((postEE.mass_gggg >= 110.0) & (postEE.mass_gggg <= 115.0)) | ((postEE.mass_gggg >= 135.0) & (postEE.mass_gggg <= 180.0))]
    hist_preEE_SB = fillHist("Ndimreweight", preEE, [low, high], binsScale=binsScale, name="preEE", normalize=False)
    hist_postEE_SB = fillHist("Ndimreweight", postEE, [low, high], binsScale=binsScale, name="postEE", normalize=False)

    for year, (hist, hist_SB), stack, maxY in zip(["2022preEE", "2022postEE"], [(hist_preEE, hist_preEE_SB), (hist_postEE, hist_postEE_SB)], [stack_preEE, stack_postEE], [7000, 30000]):
        canvas = TCanvas(f"canvas_{year}", f"canvas_{year}", 1000, 1000)
        canvas.SetLeftMargin(0.2)

        hist.GetYaxis().SetRangeUser(0.0, maxY)  # Adjust this as needed
        hist_SB.GetYaxis().SetRangeUser(0.0, maxY)  # Adjust this as needed
        hist.GetYaxis().SetTitle("Counts")
        hist.GetYaxis().SetTitleOffset(2)
        hist.GetXaxis().SetTitle("Weight Values")
        hist.SetLineColor(kBlue)
        hist_SB.SetLineColor(kGreen+3)
        hist.SetFillColor(kBlue)
        hist_SB.SetFillColor(kGreen+3)
        
        print(f"Std Dev Full Spectrum: {hist.GetStdDev()}")

        legend = TLegend(0.5,0.65,0.87,0.85)
        legend.SetTextFont(42)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(hist, "Weights from full m_{4#gamma} spectrum", "l")
        legend.AddEntry(hist_SB, "Weights from m_{4#gamma} sidebands", "l")

        #hist.Draw("hist")
        #hist_SB.Draw("hist same")
        stack.Add(hist)
        stack.Add(hist_SB)
        stack.Draw("nostack")
        legend.Draw("same")
        canvas.SaveAs(f"{cwd}/plots/BDT/{year}_reweights.png")
