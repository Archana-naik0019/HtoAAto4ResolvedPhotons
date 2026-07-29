import ROOT
import os

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(11111)

infile = ROOT.TFile.Open("PhotonPlots_improv.root")

dirs = ["SR", "CR"]
mass_points = [15,20,25,30,35,40,45,50,55,60]

outdir = "InterMassPlots"
os.makedirs(outdir, exist_ok=True)

colors = [
    ROOT.kRed+1,
    ROOT.kBlue+1,
    ROOT.kGreen+2,
    ROOT.kMagenta+1,
    ROOT.kOrange+7,
    ROOT.kCyan+2,
    ROOT.kViolet+1,
    ROOT.kPink+7,
    ROOT.kAzure+2,
    ROOT.kSpring+5
]


def normalize(hist):
    if hist and hist.Integral() > 0:
        hist.Scale(1.0/hist.Integral())


###########################################################
# Loop over SR / CR
###########################################################

for dirname in dirs:

    directory = infile.Get(dirname)

    #######################################################
    # 1. LEADING : Total + all masses
    #######################################################

    c = ROOT.TCanvas("c1","",900,800)

    leg = ROOT.TLegend(0.58,0.52,0.88,0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    hTot = directory.Get("h_LeadPs_interMass")
    normalize(hTot)

    hTot.SetLineColor(ROOT.kBlack)
    hTot.SetLineWidth(4)
    hTot.SetTitle("Leading pseudoscalar")

    hTot.GetXaxis().SetTitle("(m_{Lead}-m_{hyp})/m_{4#gamma}")
    hTot.GetYaxis().SetTitle("Normalized Events")

    hTot.Draw("hist")
    leg.AddEntry(hTot,"Total (random m_{hyp})","l")

    for i,m in enumerate(mass_points):

        h = directory.Get(f"h_LeadPs_interMass_m{m}")
        normalize(h)

        h.SetLineColor(colors[i])
        h.SetLineWidth(2)

        h.Draw("hist sames")
        leg.AddEntry(h,f"m_{{hyp}} = {m} GeV","l")

    leg.Draw()

    c.SaveAs(f"{outdir}/{dirname}_Lead_TotalPlusMasses.png")


    #######################################################
    # 2. LEADING : only masses
    #######################################################

    c = ROOT.TCanvas("c2","",900,800)

    leg = ROOT.TLegend(0.58,0.52,0.88,0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    first = True

    for i,m in enumerate(mass_points):

        h = directory.Get(f"h_LeadPs_interMass_m{m}")
        normalize(h)

        h.SetLineColor(colors[i])
        h.SetLineWidth(3)

        if first:
            h.SetTitle("Leading pseudoscalar")
            h.GetXaxis().SetTitle("(m_{Lead}-m_{hyp})/m_{4#gamma}")
            h.GetYaxis().SetTitle("Normalized Events")
            h.Draw("hist")
            first = False
        else:
            h.Draw("hist sames")

        leg.AddEntry(h,f"m_{{hyp}} = {m} GeV","l")

    leg.Draw()

    c.SaveAs(f"{outdir}/{dirname}_Lead_OnlyMasses.png")


    #######################################################
    # 3. SUBLEADING : Total + all masses
    #######################################################

    c = ROOT.TCanvas("c3","",900,800)

    leg = ROOT.TLegend(0.58,0.52,0.88,0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    hTot = directory.Get("h_SubleadPs_interMass")
    normalize(hTot)

    hTot.SetLineColor(ROOT.kBlack)
    hTot.SetLineWidth(4)
    hTot.SetTitle("Subleading pseudoscalar")

    hTot.GetXaxis().SetTitle("(m_{Sublead}-m_{hyp})/m_{4#gamma}")
    hTot.GetYaxis().SetTitle("Normalized Events")

    hTot.Draw("hist")
    leg.AddEntry(hTot,"Total (random m_{hyp})","l")

    for i,m in enumerate(mass_points):

        h = directory.Get(f"h_SubleadPs_interMass_m{m}")
        normalize(h)

        h.SetLineColor(colors[i])
        h.SetLineWidth(2)

        h.Draw("hist sames")
        leg.AddEntry(h,f"m_{{hyp}} = {m} GeV","l")

    leg.Draw()

    c.SaveAs(f"{outdir}/{dirname}_Sublead_TotalPlusMasses.png")


    #######################################################
    # 4. SUBLEADING : only masses
    #######################################################

    c = ROOT.TCanvas("c4","",900,800)

    leg = ROOT.TLegend(0.58,0.52,0.88,0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    first = True

    for i,m in enumerate(mass_points):

        h = directory.Get(f"h_SubleadPs_interMass_m{m}")
        normalize(h)

        h.SetLineColor(colors[i])
        h.SetLineWidth(3)

        if first:
            h.SetTitle("Subleading pseudoscalar")
            h.GetXaxis().SetTitle("(m_{Sublead}-m_{hyp})/m_{4#gamma}")
            h.GetYaxis().SetTitle("Normalized Events")
            h.Draw("hist")
            first = False
        else:
            h.Draw("hist sames")

        leg.AddEntry(h,f"m_{{hyp}} = {m} GeV","l")

    leg.Draw()

    c.SaveAs(f"{outdir}/{dirname}_Sublead_OnlyMasses.png")


infile.Close()

print("\nDone!")
print(f"Plots saved in : {outdir}/")
