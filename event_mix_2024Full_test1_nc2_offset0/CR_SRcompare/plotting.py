import ROOT
import os

ROOT.gROOT.SetBatch(True)

ROOT.gStyle.SetOptStat(111111)
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetTitleBorderSize(0)

# -------------------------------------------------------
# Input / Output
# -------------------------------------------------------

infile = ROOT.TFile.Open("PhotonPlots.root")

outfile = ROOT.TFile("OverlayPlots.root","RECREATE")
outfile.mkdir("Normalized")
outfile.mkdir("Unnormalized")

os.makedirs("plots/normalized", exist_ok=True)
os.makedirs("plots/unnormalized", exist_ok=True)

# -------------------------------------------------------
# Histograms
# -------------------------------------------------------

histograms = [

    # photons
    "h_pt1","h_pt2","h_pt3","h_pt4",

    "h_eta1","h_eta2","h_eta3","h_eta4",

    "h_phi1","h_phi2","h_phi3","h_phi4",

    "h_mvaid1","h_mvaid2","h_mvaid3","h_mvaid4",

    # pseudoscalars
    "h_ma1",
    "h_ma2",

    "h_avgMass",
    "h_massDiff",

    "h_pt_a1",
    "h_pt_a2",

    "h_eta_a1",
    "h_eta_a2",

    "h_phi_a1",
    "h_phi_a2",

    "h_dr_a1",
    "h_dr_a2",

    "h_dr_a1a2",

    "h_drOverM4g",

    # Opening angles
    "h_opening_a1a2",
    "h_opening_a1",
    "h_opening_a2",

    # Helicity variable
    "h_cos_ag",

    # Intermediate mass variables
    "h_LeadPs_interMass",
    "h_SubleadPs_interMass"
]

# =======================================================
# Loop over histograms
# =======================================================

for hname in histograms:

    print("Drawing",hname)

    hSR = infile.Get("SR/"+hname)
    hCR = infile.Get("CR/"+hname)

    if (not hSR) or (not hCR):
        print("Missing",hname)
        continue

    # Clone originals
    hSR = hSR.Clone(hname+"_SR")
    hCR = hCR.Clone(hname+"_CR")

    # Clone normalized copies
    hSR_norm = hSR.Clone(hname+"_SRnorm")
    hCR_norm = hCR.Clone(hname+"_CRnorm")

    if hSR_norm.Integral()>0:
        hSR_norm.Scale(1./hSR_norm.Integral())

    if hCR_norm.Integral()>0:
        hCR_norm.Scale(1./hCR_norm.Integral())

    ####################################################################
    # Draw twice
    ####################################################################

    for mode in ["Normalized","Unnormalized"]:

        outfile.cd(mode)

        canvas = ROOT.TCanvas(mode+"_"+hname,hname,900,900)

        pad1 = ROOT.TPad("pad1","pad1",0,0.30,1,1)
        pad2 = ROOT.TPad("pad2","pad2",0,0.00,1,0.30)

        pad1.SetBottomMargin(0.02)
        pad2.SetTopMargin(0.05)
        pad2.SetBottomMargin(0.35)

        pad1.SetRightMargin(0.28)
        pad2.SetRightMargin(0.28)

        pad1.Draw()
        pad2.Draw()

        if mode=="Normalized":

            h1 = hSR_norm
            h2 = hCR_norm
            ytitle = "Normalized Entries"

        else:

            h1 = hSR
            h2 = hCR
            ytitle = "Events"

        # ---------------- Style ----------------

        h1.SetLineColor(ROOT.kBlack)
        h2.SetLineColor(ROOT.kRed)

        h1.SetLineWidth(2)
        h2.SetLineWidth(2)

        h1.SetStats(1)
        h2.SetStats(1)

        ymax = 1.25*max(h1.GetMaximum(),h2.GetMaximum())
        h1.SetMaximum(ymax)

        h1.GetYaxis().SetTitle(ytitle)

        # Fix x-axis titles for special variables
        if hname == "h_LeadPs_interMass":
            h1.GetXaxis().SetTitle("(m_{Lead}-m_{hyp})/m_{4#gamma}")
            h2.GetXaxis().SetTitle("(m_{Lead}-m_{hyp})/m_{4#gamma}")

        if hname == "h_SubleadPs_interMass":
            h1.GetXaxis().SetTitle("(m_{Sublead}-m_{hyp})/m_{4#gamma}")
            h2.GetXaxis().SetTitle("(m_{Sublead}-m_{hyp})/m_{4#gamma}")
        
        # Hide x-axis labels/titles on upper pad
        h1.GetXaxis().SetLabelSize(0)
        h1.GetXaxis().SetTitleSize(0)

        h2.GetXaxis().SetLabelSize(0)
        h2.GetXaxis().SetTitleSize(0)

        # Restrict x-axis for leading 4 photon pT
        if hname in ["h_pt1","h_pt2","h_pt3","h_pt4"]:
            h1.GetXaxis().SetRangeUser(0,100)
            h2.GetXaxis().SetRangeUser(0,100)

        # Use log-scale y-axis for photon MVA ID histograms
        if hname in ["h_mvaid1", "h_mvaid2", "h_mvaid3", "h_mvaid4"]:
            pad1.SetLogy()

        # Draw

        pad1.cd()

        h1.Draw("HIST")
        h2.Draw("HIST SAMES")

        canvas.Update()

        # ---------------- Stat boxes ----------------

        st1 = h1.GetListOfFunctions().FindObject("stats")
        st2 = h2.GetListOfFunctions().FindObject("stats")

        if st1:

            st1.SetTextColor(ROOT.kBlack)

            st1.SetX1NDC(0.73)
            st1.SetX2NDC(0.98)

            st1.SetY1NDC(0.78)
            st1.SetY2NDC(0.93)

        if st2:

            st2.SetTextColor(ROOT.kRed)

            st2.SetX1NDC(0.73)
            st2.SetX2NDC(0.98)

            st2.SetY1NDC(0.58)
            st2.SetY2NDC(0.73)

        canvas.Modified()
        canvas.Update()

        # ---------------- Legend ----------------

        leg = ROOT.TLegend(0.73,0.46,0.97,0.56)

        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.028)

        leg.AddEntry(h1,"SR","l")
        leg.AddEntry(h2,"CR","l")

        leg.Draw()

        # Ratio pad

        pad2.cd()

        ratio = h1.Clone(f"ratio_{hname}_{mode}")
        ratio.Sumw2()
        ratio.Divide(h2)

        if hname in ["h_pt1","h_pt2","h_pt3","h_pt4"]:
            ratio.GetXaxis().SetRangeUser(0,100)

        ratio.SetStats(0)
        ratio.SetTitle("")
        ratio.GetXaxis().SetTitle(h1.GetXaxis().GetTitle())

        ratio.GetYaxis().SetTitle("SR/CR")
        ratio.GetYaxis().CenterTitle()

        ratio.GetYaxis().SetNdivisions(503)

        ratio.GetYaxis().SetTitleSize(0.08)
        ratio.GetYaxis().SetTitleOffset(0.60)

        ratio.GetYaxis().SetLabelSize(0.07)

        ratio.GetXaxis().SetLabelSize(0.08)
        ratio.GetXaxis().SetTitleSize(0.09)
        ratio.GetXaxis().SetTitleOffset(1.15)

        ratio.SetMinimum(0.0)
        ratio.SetMaximum(2.0)

        ratio.SetLineColor(ROOT.kBlack)
        ratio.SetMarkerStyle(20)
        ratio.SetMarkerSize(0.8)

        ratio.Draw("E1")

        line = ROOT.TLine(
             ratio.GetXaxis().GetXmin(),
             1.0,
             ratio.GetXaxis().GetXmax(),
             1.0
         )

        line.SetLineStyle(2)
        line.Draw()

        canvas.cd()
        
        canvas.Modified()
        canvas.Update()

        canvas.Write()

        if mode=="Normalized":
            canvas.SaveAs("plots/normalized/"+hname+".png")
        else:
            canvas.SaveAs("plots/unnormalized/"+hname+".png")

        canvas.Close()

# =======================================================

outfile.Close()
infile.Close()

print("Done!")
