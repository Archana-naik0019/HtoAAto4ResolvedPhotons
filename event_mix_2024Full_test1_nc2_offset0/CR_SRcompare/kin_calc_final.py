import ROOT
import random

ROOT.gROOT.SetBatch(True)

# best diphoton pairing, this is to get the a1, a2 candidates------------

def findBestPairing(g1, g2, g3, g4):

    pairings = [

        ((g1+g2), (g3+g4), (0,1,2,3)),

        ((g1+g3), (g2+g4), (0,2,1,3)),

        ((g1+g4), (g2+g3), (0,3,1,2))

    ]

    best = None
    bestDiff = 1e9

    for a1, a2, indices in pairings:

        diff = abs(a1.M() - a2.M())

        if diff < bestDiff:

            bestDiff = diff
            best = (a1, a2, indices)

    return best
#-------------------------------------------------------------------

files = {
    "SR" : "m4g_115to135.root",
    "CR" : "m4g_110to115_and_115to135.root"
}

mass_points = [15,20,25,30,35,40,45,50,55,60] # added since I want the "inter-mass" distributions separated for the different mass points

outfile = ROOT.TFile("PhotonPlots_improv.root","RECREATE")



for dirname, filename in files.items():

    print("Processing", filename)

    infile = ROOT.TFile.Open(filename)
    tree = infile.Get("Events")

    outfile.mkdir(dirname)
    outfile.cd(dirname)


    h_pt  = []
    h_eta = []
    h_phi = []
    h_mva = []

    h_ma1 = None
    h_ma2 = None

    h_avgMass = None
    h_massDiff = None

    h_pt_a1 = None
    h_pt_a2 = None

    h_eta_a1 = None
    h_eta_a2 = None

    h_phi_a1 = None
    h_phi_a2 = None

    h_dr_a1 = None
    h_dr_a2 = None

    h_dr_a1a2 = None
    h_drOverM4g = None

    for i in range(4):

        h_pt.append(
            ROOT.TH1F(
                f"h_pt{i+1}",
                f"Photon {i+1} pT;Photon p_{{T}} [GeV];Events",
                300,0,300
            )
        )

        h_eta.append(
            ROOT.TH1F(
                f"h_eta{i+1}",
                f"Photon {i+1} eta;#eta;Events",
                60,-3,3
            )
        )

        h_phi.append(
            ROOT.TH1F(
                f"h_phi{i+1}",
                f"Photon {i+1} phi;#phi;Events",
                64,-3.2,3.2
            )
        )

        h_mva.append(
            ROOT.TH1F(
                f"h_mvaid{i+1}",
                f"Photon {i+1} mvaID;mvaID;Events",
                100,-1,1
            )
        )

        #------------------------------------------------------------
    h_ma1 = ROOT.TH1F("h_ma1","a1 mass;m_{a1} (GeV);Events",100,0,100)

    h_ma2 = ROOT.TH1F("h_ma2","a2 mass;m_{a2} (GeV);Events",100,0,100)

    h_avgMass = ROOT.TH1F("h_avgMass",
                      "(m_{a1}+m_{a2})/2;Average mass (GeV);Events",
                      100,0,100)

    h_massDiff = ROOT.TH1F("h_massDiff",
                       "|m_{a1}-m_{a2}|;Mass difference (GeV);Events",
                       100,0,50)

    h_pt_a1 = ROOT.TH1F("h_pt_a1","a1 p_{T};p_{T} (GeV);Events",100,0,200)
    h_pt_a2 = ROOT.TH1F("h_pt_a2","a2 p_{T};p_{T} (GeV);Events",100,0,200)
       
    h_eta_a1 = ROOT.TH1F("h_eta_a1","a1 #eta;#eta;Events",60,-3,3)
    h_eta_a2 = ROOT.TH1F("h_eta_a2","a2 #eta;#eta;Events",60,-3,3)

    h_phi_a1 = ROOT.TH1F("h_phi_a1","a1 #phi;#phi;Events",64,-3.2,3.2)
    h_phi_a2 = ROOT.TH1F("h_phi_a2","a2 #phi;#phi;Events",64,-3.2,3.2)

    h_dr_a1 = ROOT.TH1F("h_dr_a1",
                    "#DeltaR(photons in a1);#DeltaR;Events",
                    80,0,8)
    h_dr_a2 = ROOT.TH1F("h_dr_a2",
                    "#DeltaR(photons in a2);#DeltaR;Events",
                    80,0,8)

    h_dr_a1a2 = ROOT.TH1F("h_dr_a1a2",
                      "#DeltaR(a1,a2);#DeltaR;Events",
                      80,0,8)
    h_drOverM4g = ROOT.TH1F("h_drOverM4g",
                        "#DeltaR(a1,a2)/m_{4#gamma};#DeltaR/m_{4#gamma};Events",
                        100,0,0.05)
    #opening angles
    h_opening_a1a2 = ROOT.TH1F("h_opening_a1a2",
                         "Opening angle(a1,a2);Opening angle (rad);Events",
                         64,0,3.2)

    h_opening_a1 = ROOT.TH1F("h_opening_a1",
                         "Opening angle(#gamma,#gamma in a1);Opening angle (rad);Events",
                         64,0,3.2)

    h_opening_a2 = ROOT.TH1F("h_opening_a2",
                         "Opening angle(#gamma,#gamma in a2);Opening angle (rad);Events",
                         64,0,3.2)

    h_cos_ag = ROOT.TH1F("h_cos_ag",
                         "cos#theta*;cos#theta*;Events",
                         100,-1,1)

    h_LeadPs_interMass = ROOT.TH1F("h_LeadPs_interMass",
                          "(m_{Lead}-m_{hyp})/m_{4#gamma};Variable;Events",
                          100,-0.5,0.5)

    h_SubleadPs_interMass = ROOT.TH1F("h_SubleadPs_interMass",
                          "(m_{Sublead}-m_{hyp})/m_{4#gamma};Variable;Events",
                          100,-0.5,0.5)

    h_m4g = ROOT.TH1F("h_m4g",
                          "Four-photon invariant mass;m_{4#gamma} (GeV);Events",
                          700,110,180)
    # 2D histograms

    h2_ma1_vs_opening_a1 = ROOT.TH2F("h2_ma1_vs_opening_a1",
                          "m_{a1} vs opening angle(#gamma,#gamma in a1);m_{a1} (GeV);Opening angle (rad)",
                         100,0,100,
                        64,0,3.2)

    h2_ma2_vs_opening_a2 = ROOT.TH2F("h2_ma2_vs_opening_a2",
                          "m_{a2} vs opening angle(#gamma,#gamma in a2);m_{a2} (GeV);Opening angle (rad)",
                           100,0,100,
                           64,0,3.2)

    # histograms for "inter-mass" of separate 'a' mass points
    h_LeadPs_interMass_each = {}
    h_SubleadPs_interMass_each = {}

    for m in mass_points:

        h_LeadPs_interMass_each[m] = ROOT.TH1F(
            f"h_LeadPs_interMass_m{m}",
            f"(m_{{Lead}}-{m})/m_{{4#gamma}} (m_{{hyp}}={m});Variable;Events",
            100,-0.5,0.5)

        h_SubleadPs_interMass_each[m] = ROOT.TH1F(
            f"h_SubleadPs_interMass_m{m}",
            f"(m_{{Sublead}}-{m})/m_{{4#gamma}} (m_{{hyp}}={m});Variable;Events",
            100,-0.5,0.5)


    nentries = tree.GetEntries()

    for iev in range(nentries):

        tree.GetEntry(iev)

        m_hyp = random.choice([15,20,25,30,35,40,45,50,55,60])

        if iev%50000==0:
            print(iev,"/",nentries)

        g1 = ROOT.TLorentzVector()
        g2 = ROOT.TLorentzVector()
        g3 = ROOT.TLorentzVector()
        g4 = ROOT.TLorentzVector()

        g1.SetPtEtaPhiM(tree.Photon_pt[0],
                    tree.Photon_eta[0],
                    tree.Photon_phi[0],
                    0.)

        g2.SetPtEtaPhiM(tree.Photon_pt[1],
                    tree.Photon_eta[1],
                    tree.Photon_phi[1],
                    0.)

        g3.SetPtEtaPhiM(tree.Photon_pt[2],
                    tree.Photon_eta[2],
                    tree.Photon_phi[2],
                    0.)

        g4.SetPtEtaPhiM(tree.Photon_pt[3],
                    tree.Photon_eta[3],
                    tree.Photon_phi[3],
                    0.)

        # Compute m4g
        m4g = (g1+g2+g3+g4).M()

        # Best pairing
        a1, a2, idx = findBestPairing(g1, g2, g3, g4)
        i1, i2, i3, i4 = idx

        # Define a1 as the leading-pT pseudoscalar
        if a2.Pt() > a1.Pt():

            a1, a2 = a2, a1

            i1, i2, i3, i4 = i3, i4, i1, i2

        # Compute pseudoscalar variables
        ma1 = a1.M()
        ma2 = a2.M()

        avgMass = 0.5*(ma1 + ma2)
        massDiff = abs(ma1 - ma2)

        #defining the "inter variables"
        LeadPs_interMass = (a1.M() - m_hyp)/m4g
        SubleadPs_interMass = (a2.M() - m_hyp)/m4g

        pt_a1 = a1.Pt()
        pt_a2 = a2.Pt()

        eta_a1 = a1.Eta()
        eta_a2 = a2.Eta()

        phi_a1 = a1.Phi()
        phi_a2 = a2.Phi()

        # DeltaR of photons inside a1 and a2
        photons = [g1, g2, g3, g4]

        dr_a1 = photons[i1].DeltaR(photons[i2])
        dr_a2 = photons[i3].DeltaR(photons[i4])

        # DeltaR between pseudoscalars
        dr_a1a2 = a1.DeltaR(a2)

        # DeltaR(a1,a2)/m4g
        drOverM4g = dr_a1a2 / m4g

        # Opening angles (3D)
        opening_a1a2 = a1.Vect().Angle(a2.Vect())
        opening_a1 = photons[i1].Vect().Angle(photons[i2].Vect())
        opening_a2 = photons[i3].Vect().Angle(photons[i4].Vect())

        # Helicity angle
        leadPho = photons[i1]
        pho_rest = ROOT.TLorentzVector(leadPho)
        pho_rest.Boost(-a1.BoostVector())
        cos_ag = ROOT.TMath.Cos(pho_rest.Theta())

        # Existing photon histograms
        for i in range(4):

            h_pt[i].Fill(tree.Photon_pt[i])
            h_eta[i].Fill(tree.Photon_eta[i])
            h_phi[i].Fill(tree.Photon_phi[i])
            h_mva[i].Fill(tree.Photon_mvaID[i])

        # Fill pseudoscalar histograms
        h_ma1.Fill(ma1)
        h_ma2.Fill(ma2)

        h_avgMass.Fill(avgMass)
        h_massDiff.Fill(massDiff)

        h_pt_a1.Fill(pt_a1)
        h_pt_a2.Fill(pt_a2)

        h_eta_a1.Fill(eta_a1)
        h_eta_a2.Fill(eta_a2)

        h_phi_a1.Fill(phi_a1)
        h_phi_a2.Fill(phi_a2)

        h_dr_a1.Fill(dr_a1)
        h_dr_a2.Fill(dr_a2)

        h_dr_a1a2.Fill(dr_a1a2)

        h_drOverM4g.Fill(drOverM4g)

        h_opening_a1a2.Fill(opening_a1a2)
        h_opening_a1.Fill(opening_a1)
        h_opening_a2.Fill(opening_a2)

        h2_ma1_vs_opening_a1.Fill(ma1, opening_a1)
        h2_ma2_vs_opening_a2.Fill(ma2, opening_a2)

        h_cos_ag.Fill(cos_ag)
        h_LeadPs_interMass.Fill(LeadPs_interMass)
        h_SubleadPs_interMass.Fill(SubleadPs_interMass)

        for m in mass_points:
            h_LeadPs_interMass_each[m].Fill((a1.M() - m)/m4g)
            h_SubleadPs_interMass_each[m].Fill((a2.M() - m)/m4g)

        h_m4g.Fill(m4g)

    for h in h_pt:
        h.Write()

    for h in h_eta:
        h.Write()

    for h in h_phi:
        h.Write()

    for h in h_mva:
        h.Write()

    # Write pseudoscalar histograms
    h_ma1.Write()
    h_ma2.Write()

    h_avgMass.Write()
    h_massDiff.Write()

    h_pt_a1.Write()
    h_pt_a2.Write()

    h_eta_a1.Write()
    h_eta_a2.Write()

    h_phi_a1.Write()
    h_phi_a2.Write()

    h_dr_a1.Write()
    h_dr_a2.Write()

    h_dr_a1a2.Write()
    h_drOverM4g.Write()

    h_opening_a1a2.Write()
    h_opening_a1.Write()
    h_opening_a2.Write()

    h2_ma1_vs_opening_a1.Write()
    h2_ma2_vs_opening_a2.Write()

    h_cos_ag.Write()
    h_LeadPs_interMass.Write()
    h_SubleadPs_interMass.Write()

    for m in mass_points:
        h_LeadPs_interMass_each[m].Write()
        h_SubleadPs_interMass_each[m].Write()

    h_m4g.Write()

    infile.Close()

outfile.Close()

print("Done!")
