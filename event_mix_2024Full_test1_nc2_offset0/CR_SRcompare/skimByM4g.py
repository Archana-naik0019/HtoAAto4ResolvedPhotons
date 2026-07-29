import ROOT
import os
from glob import glob
from array import array

ROOT.gROOT.SetBatch(True)

input_dir = "/eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/event_mix_2024Full_test1_nc2_offset0"

# Build TChain
chain = ROOT.TChain("Events")

nfiles = 0

for root, dirs, files in os.walk(input_dir):
    for f in files:
        if f.endswith(".root"):
            chain.Add(os.path.join(root, f))
            nfiles += 1

print(f"Added {nfiles} ROOT files")
print(f"Total events = {chain.GetEntries()}")

# Output files
f_sr = ROOT.TFile("m4g_115to135.root", "RECREATE")
t_sr = chain.CloneTree(0)

f_notsr = ROOT.TFile("m4g_not115to135.root", "RECREATE")
t_notsr = chain.CloneTree(0)

f_srandcr = ROOT.TFile("m4g_110to135.root", "RECREATE")
t_srandcr = chain.CloneTree(0)

f_cr = ROOT.TFile("m4g_110to115_and_115to135.root", "RECREATE")
t_cr = chain.CloneTree(0)

# m4g branches
m4g_sr = array('f', [0.])
m4g_notsr = array('f', [0.])
m4g_srandcr = array('f', [0.])
m4g_cr = array('f', [0.])

t_sr.Branch("m4g", m4g_sr, "m4g/F")
t_notsr.Branch("m4g", m4g_notsr, "m4g/F")
t_srandcr.Branch("m4g", m4g_srandcr, "m4g/F")
t_cr.Branch("m4g", m4g_cr, "m4g/F")

# Lorentz vectors
g1 = ROOT.TLorentzVector()
g2 = ROOT.TLorentzVector()
g3 = ROOT.TLorentzVector()
g4 = ROOT.TLorentzVector()

# Event loop
N = chain.GetEntries()

for i in range(N):

    chain.GetEntry(i)

    if i % 100000 == 0:
        print(f"{i}/{N}")

    g1.SetPtEtaPhiM(chain.Photon_pt[0],
                    chain.Photon_eta[0],
                    chain.Photon_phi[0],
                    0.0)

    g2.SetPtEtaPhiM(chain.Photon_pt[1],
                    chain.Photon_eta[1],
                    chain.Photon_phi[1],
                    0.0)

    g3.SetPtEtaPhiM(chain.Photon_pt[2],
                    chain.Photon_eta[2],
                    chain.Photon_phi[2],
                    0.0)

    g4.SetPtEtaPhiM(chain.Photon_pt[3],
                    chain.Photon_eta[3],
                    chain.Photon_phi[3],
                    0.0)

    m4g = (g1 + g2 + g3 + g4).M()


    # SR
    if 115. < m4g < 135.:
        m4g_sr[0] = m4g
        t_sr.Fill()

    # NOT SR
    if m4g <= 115. or m4g >= 135.:
        m4g_notsr[0] = m4g
        t_notsr.Fill()

    # 110-135 (SR+CR)
    if 110. < m4g < 135.:
        m4g_srandcr[0] = m4g
        t_srandcr.Fill()

    # sidebands (CR)
    if (110. < m4g < 115.) or (115. < m4g < 135.):
        m4g_cr[0] = m4g
        t_cr.Fill()


# output written
f_sr.cd()
t_sr.Write()
f_sr.Close()

f_notsr.cd()
t_notsr.Write()
f_notsr.Close()

f_srandcr.cd()
t_srandcr.Write()
f_srandcr.Close()

f_cr.cd()
t_cr.Write()
f_cr.Close()

print("\nDone!")
