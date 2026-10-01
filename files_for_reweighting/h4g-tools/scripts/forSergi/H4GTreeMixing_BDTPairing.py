#!/usr/bin/python
import numpy as n
import ROOT
# from ROOT import *
import sys, getopt
from array import array
from H4GSkimTools import *
from optparse import OptionParser
import operator

if __name__ == '__main__':


  parser = OptionParser()
  parser.add_option(   "-i", "--inputFile",   dest="inputFile",    default="input.root",       type="string",  help="Input file" )
  parser.add_option(   "-o", "--outputFile",  dest="outputFile",   default="output.root",      type="string",  help="Output file")
  parser.add_option(   "-t", "--tree",  dest="tree",   default="tree",      type="string",  help="tree")
  parser.add_option(   "-w", "--weight",  dest="weight",   default="",      type="string",  help="weight")
  parser.add_option(   "-m", "--mass",  dest="mass",   default="60",      type="string",  help="mass")

  #/afs/cern.ch/user/a/amassiro/public/H4G/signal_m_50.root

  (options, args) = parser.parse_args()

  itree = ROOT.TChain("h4gCandidateDumper/trees/"+str(options.tree))
  itree.AddFile(options.inputFile)

  print "Total number of events to be analyzed:", itree.GetEntries()

  outRoot = ROOT.TFile(options.outputFile, "RECREATE")
  outtree = itree.CloneTree(0)

  weight = options.weight
  mass = options.mass

  def getCosThetaStar_CS():
      h_lor = a1 + a2
      h = ROOT.TLorentzVector(0,0,0,0)
      h.SetPxPyPzE(h_lor.Px(),h_lor.Py(),h_lor.Pz(),h_lor.E())
      a1_lor = a1
      a_1 = ROOT.TLorentzVector(0,0,0,0)
      a_1.SetPxPyPzE(a1_lor.Px(),a1_lor.Py(),a1_lor.Pz(),a1_lor.E())
      a_1.Boost(-h.BoostVector())
      return a_1.CosTheta()

  def helicityCosTheta( Booster,  Boosted):
      BoostVector = Booster.BoostVector()
      Boosted.Boost( -BoostVector.x(), -BoostVector.y(), -BoostVector.z() );
      return Boosted.CosTheta()

  def CosThetaAngles():
      helicityThetas = []
      Boosted_a1 = ROOT.TLorentzVector(0,0,0,0)
      Boosted_a1.SetPxPyPzE(a1.Px(),a1.Py(),a1.Pz(),a1.E())
      BoostedLeadingPhoton_a1 = ROOT.TLorentzVector(0,0,0,0)
      BoostedLeadingPhoton_a1.SetPxPyPzE(a1_pho1.Px(),a1_pho1.Py(),a1_pho1.Pz(),a1_pho1.E())
      helicityThetas.append( helicityCosTheta(Boosted_a1, BoostedLeadingPhoton_a1))

      Boosted_a2 = ROOT.TLorentzVector(0,0,0,0)
      Boosted_a2.SetPxPyPzE(a2.Px(),a2.Py(),a2.Pz(),a2.E())
      BoostedLeadingPhoton_a2 = ROOT.TLorentzVector(0,0,0,0)
      BoostedLeadingPhoton_a2.SetPxPyPzE(a2_pho1.Px(),a2_pho1.Py(),a2_pho1.Pz(),a2_pho1.E())
      helicityThetas.append( helicityCosTheta(Boosted_a2, BoostedLeadingPhoton_a2))
      return helicityThetas

  pho1_pt_mix = array('f', [0])
  pho2_pt_mix = array('f', [0])
  pho3_pt_mix = array('f', [0])
  pho4_pt_mix = array('f', [0])
  pho1_eta_mix = array('f', [0])
  pho2_eta_mix = array('f', [0])
  pho3_eta_mix = array('f', [0])
  pho4_eta_mix = array('f', [0])
  pho1_MVA_mix = array('f', [0])
  pho2_MVA_mix = array('f', [0])
  pho3_MVA_mix = array('f', [0])
  pho4_MVA_mix = array('f', [0])
  pho_minMVA_mix = array('f', [0])
  pho_maxMVA_mix = array('f', [0])
  pho12_DR_mix = array('f', [0])
  pho13_DR_mix = array('f', [0])
  pho14_DR_mix = array('f', [0])
  pho23_DR_mix = array('f', [0])
  pho24_DR_mix = array('f', [0])
  pho34_DR_mix = array('f', [0])
  pho12_M_mix = array('f', [0])
  pho13_M_mix = array('f', [0])
  pho14_M_mix = array('f', [0])
  pho23_M_mix = array('f', [0])
  pho24_M_mix = array('f', [0])
  pho34_M_mix = array('f', [0])
  pairMVAscore_mix = array('f', [0])
  a1_mass_mix = array('f',[0])
  a2_mass_mix = array('f',[0])
  avg_a_mass_mix = array('f',[0])
  a1_Pt_mix = array('f',[0])
  a2_Pt_mix = array('f',[0])
  a1_Eta_mix = array('f',[0])
  a2_Eta_mix = array('f',[0])
  a1_Phi_mix = array('f',[0])
  a2_Phi_mix = array('f',[0])
  a1_Energy_mix = array('f',[0])
  a2_Energy_mix = array('f',[0])
  a1_DR_mix = array('f',[0])
  a2_DR_mix = array('f',[0])
  dM1_gen1_mix = array('f',[0])
  dM2_gen1_mix = array('f',[0])
  a1_a2_DR_mix = array('f',[0])
  a1_iPho1_mix = array('f',[0])
  a1_iPho2_mix = array('f',[0])
  a2_iPho1_mix = array('f',[0])
  a2_iPho2_mix = array('f',[0])
  a1_PtOverMass_mix = array('f',[0])
  a2_PtOverMass_mix = array('f',[0])
  a1_Pho1PtOvera1Mass_mix = array('f',[0])
  a1_Pho2PtOvera1Mass_mix = array('f',[0])
  a2_Pho1PtOvera2Mass_mix = array('f',[0])
  a2_Pho2PtOvera2Mass_mix = array('f',[0])
  CTStarCS_mix = array('f',[0])
  CT_a1Pho1_mix = array('f',[0])
  CT_a2Pho1_mix = array('f',[0])
  tp_mass_mix  = array('f',[0])
  tp_pt_mix  = array('f',[0])
  tp_eta_mix = array('f',[0])
  tp_phi_mix  = array('f',[0])

  _pho1_pt_mix = outtree.Branch('pho1_pt_mix',pho1_pt_mix,'pho1_pt_mix/F')
  _pho2_pt_mix = outtree.Branch('pho2_pt_mix',pho2_pt_mix,'pho2_pt_mix/F')
  _pho3_pt_mix = outtree.Branch('pho3_pt_mix',pho3_pt_mix,'pho3_pt_mix/F')
  _pho4_pt_mix = outtree.Branch('pho4_pt_mix',pho4_pt_mix,'pho4_pt_mix/F')
  _pho1_eta_mix = outtree.Branch('pho1_eta_mix',pho1_eta_mix,'pho1_eta_mix/F')
  _pho2_eta_mix = outtree.Branch('pho2_eta_mix',pho2_eta_mix,'pho2_eta_mix/F')
  _pho3_eta_mix = outtree.Branch('pho3_eta_mix',pho3_eta_mix,'pho3_eta_mix/F')
  _pho4_eta_mix = outtree.Branch('pho4_eta_mix',pho4_eta_mix,'pho4_eta_mix/F')
  _pho1_MVA_mix = outtree.Branch('pho1_MVA_mix',pho1_MVA_mix,'pho1_MVA_mix/F')
  _pho2_MVA_mix = outtree.Branch('pho2_MVA_mix',pho2_MVA_mix,'pho2_MVA_mix/F')
  _pho3_MVA_mix = outtree.Branch('pho3_MVA_mix',pho3_MVA_mix,'pho3_MVA_mix/F')
  _pho4_MVA_mix = outtree.Branch('pho4_MVA_mix',pho4_MVA_mix,'pho4_MVA_mix/F')
  _pho_minMVA_mix = outtree.Branch('pho_minMVA_mix',pho_minMVA_mix,'pho_minMVA_mix/F')
  _pho_maxMVA_mix = outtree.Branch('pho_maxMVA_mix',pho_maxMVA_mix,'pho_maxMVA_mix/F')
  _pho12_DR_mix = outtree.Branch('pho12_DR_mix',pho12_DR_mix,'pho12_DR_mix/F')
  _pho13_DR_mix = outtree.Branch('pho13_DR_mix',pho13_DR_mix,'pho13_DR_mix/F')
  _pho14_DR_mix = outtree.Branch('pho14_DR_mix',pho14_DR_mix,'pho14_DR_mix/F')
  _pho23_DR_mix = outtree.Branch('pho23_DR_mix',pho23_DR_mix,'pho23_DR_mix/F')
  _pho24_DR_mix = outtree.Branch('pho24_DR_mix',pho24_DR_mix,'pho24_DR_mix/F')
  _pho34_DR_mix = outtree.Branch('pho34_DR_mix',pho34_DR_mix,'pho34_DR_mix/F')
  _pho12_M_mix = outtree.Branch('pho12_M_mix',pho12_M_mix,'pho12_M_mix/F')
  _pho13_M_mix = outtree.Branch('pho13_M_mix',pho13_M_mix,'pho13_M_mix/F')
  _pho14_M_mix = outtree.Branch('pho14_M_mix',pho14_M_mix,'pho14_M_mix/F')
  _pho23_M_mix = outtree.Branch('pho23_M_mix',pho23_M_mix,'pho23_M_mix/F')
  _pho24_M_mix = outtree.Branch('pho24_M_mix',pho24_M_mix,'pho24_M_mix/F')
  _pho34_M_mix = outtree.Branch('pho34_M_mix',pho34_M_mix,'pho34_M_mix/F')
  _pairMVAscore_mix = outtree.Branch('pairMVAscore_mix',pairMVAscore_mix,'pairMVAscore_mix/F')
  _a1_mass_mix = outtree.Branch('a1_mass_mix',a1_mass_mix,'a1_mass_mix/F')
  _a2_mass_mix = outtree.Branch('a2_mass_mix',a2_mass_mix,'a2_mass_mix/F')
  _avg_a_mass_mix = outtree.Branch('avg_a_mass_mix',avg_a_mass_mix,'avg_a_mass_mix/F')
  _a1_Pt_mix = outtree.Branch('a1_Pt_mix',a1_Pt_mix,'a1_Pt_mix/F')
  _a2_Pt_mix = outtree.Branch('a2_Pt_mix',a2_Pt_mix,'a2_Pt_mix/F')
  _a1_Eta_mix = outtree.Branch('a1_Eta_mix',a1_Eta_mix,'a1_Eta_mix/F')
  _a2_Eta_mix = outtree.Branch('a2_Eta_mix',a2_Eta_mix,'a2_Eta_mix/F')
  _a1_Phi_mix = outtree.Branch('a1_Phi_mix',a1_Phi_mix,'a1_Phi_mix/F')
  _a2_Phi_mix = outtree.Branch('a2_Phi_mix',a2_Phi_mix,'a2_Phi_mix/F')
  _a1_Energy_mix = outtree.Branch('a1_Energy_mix',a1_Energy_mix,'a1_Energy_mix/F')
  _a2_Energy_mix = outtree.Branch('a2_Energy_mix',a2_Energy_mix,'a2_Energy_mix/F')
  _dM1_gen1_mix = outtree.Branch('dM1_gen1_mix',dM1_gen1_mix,'dM1_gen1_mix/F')
  _dM2_gen1_mix = outtree.Branch('dM2_gen1_mix',dM2_gen1_mix,'dM2_gen1_mix/F')
  _a1_DR_mix = outtree.Branch('a1_DR_mix',a1_DR_mix,'a1_DR_mix/F')
  _a2_DR_mix = outtree.Branch('a2_DR_mix',a2_DR_mix,'a2_DR_mix/F')
  _a1_a2_DR_mix = outtree.Branch('a1_a2_DR_mix',a1_a2_DR_mix,'a1_a2_DR_mix/F')
  _a1_iPho1_mix = outtree.Branch('a1_iPho1_mix',a1_iPho1_mix,'a1_iPho1_mix/F')
  _a1_iPho2_mix = outtree.Branch('a1_iPho2_mix',a1_iPho2_mix,'a1_iPho2_mix/F')
  _a2_iPho1_mix = outtree.Branch('a2_iPho1_mix',a2_iPho1_mix,'a2_iPho1_mix/F')
  _a2_iPho2_mix = outtree.Branch('a2_iPho2_mix',a2_iPho2_mix,'a2_iPho2_mix/F')
  _a1_PtOverMass_mix = outtree.Branch('a1_PtOverMass_mix',a1_PtOverMass_mix,'a1_PtOverMass_mix/F')
  _a2_PtOverMass_mix = outtree.Branch('a2_PtOverMass_mix',a2_PtOverMass_mix,'a2_PtOverMass_mix/F')
  _CTStarCS_mix = outtree.Branch('CTStarCS_mix',CTStarCS_mix,'CTStarCS_mix/F')
  _CT_a1Pho1_mix = outtree.Branch('CT_a1Pho1_mix',CT_a1Pho1_mix,'CT_a1Pho1_mix/F')
  _CT_a2Pho1_mix = outtree.Branch('CT_a2Pho1_mix',CT_a2Pho1_mix,'CT_a2Pho1_mix/F')
  _a1_Pho1PtOvera1Mass_mix = outtree.Branch('a1_Pho1PtOvera1Mass_mix',a1_Pho1PtOvera1Mass_mix,'a1_Pho1PtOvera1Mass_mix/F')
  _a1_Pho2PtOvera1Mass_mix = outtree.Branch('a1_Pho2PtOvera1Mass_mix',a1_Pho2PtOvera1Mass_mix,'a1_Pho2PtOvera1Mass_mix/F')
  _a2_Pho1PtOvera2Mass_mix = outtree.Branch('a2_Pho1PtOvera2Mass_mix',a2_Pho1PtOvera2Mass_mix,'a2_Pho1PtOvera2Mass_mix/F')
  _a2_Pho2PtOvera2Mass_mix = outtree.Branch('a2_Pho2PtOvera2Mass_mix',a2_Pho2PtOvera2Mass_mix,'a2_Pho2PtOvera2Mass_mix/F')
  _tp_mass_mix = outtree.Branch('tp_mass_mix',tp_mass_mix,'tp_mass_mix/F')
  _tp_pt_mix = outtree.Branch('tp_pt_mix',tp_pt_mix,'tp_pt_mix/F')
  _tp_eta_mix = outtree.Branch('tp_eta_mix',tp_eta_mix,'tp_eta_mix/F')
  _tp_phi_mix = outtree.Branch('tp_phi_mix',tp_phi_mix,'tp_phi_mix/F')

  reader = ROOT.TMVA.Reader()
  reader.AddVariable('dipho1_energy',a1_Energy_mix)
  reader.AddVariable('dipho1_pt',a1_Pt_mix)
  reader.AddVariable('dipho1_eta',a1_Eta_mix)
  reader.AddVariable('dipho1_dR',a1_DR_mix)
  reader.AddVariable('deltaM1_gen1',dM1_gen1_mix)
  reader.AddVariable('dipho2_energy',a2_Energy_mix)
  reader.AddVariable('dipho2_pt',a2_Pt_mix)
  reader.AddVariable('dipho2_eta',a2_Eta_mix)
  reader.AddVariable('dipho2_dR',a2_DR_mix)
  reader.AddVariable('deltaM2_gen1',dM2_gen1_mix)
  reader.AddVariable('dipair_dR',a1_a2_DR_mix)

  reader.BookMVA("BDT",str(weight))

  treeSkimmer = SkimmedTreeTools()
  otree = treeSkimmer.MakeSkimmedTree()


  for b in itree.GetListOfBranches():
    b.SetStatus(0)

  otree = itree.CloneTree(0)
  # BUT keep all branches "active" in the old tree
  itree.SetBranchStatus('*'  ,1)


  # prepare for the loop ...
  eventsToRun = itree.GetEntries()

  # print " itree.GetListOfBranches = ", itree.GetListOfBranches

  # and loop
  for evt in range(0, eventsToRun):
  # for evt in range(0,100):

    if evt%1 == 0: print "## Analyzing event ", evt
    tempPhos = []
    itree.GetEntry(evt)

    dp1_p1i = itree.dp1_p1i
    dp1_p2i = itree.dp1_p2i
    dp2_p1i = itree.dp2_p1i
    dp2_p2i = itree.dp2_p2i

    # treeSkimmer.FillEvent(itree)

    if evt == eventsToRun :
      itree.GetEntry(0)
    else :
      itree.GetEntry(evt)
      # print "evt", evt

    ObjList = [key.GetName() for key in  itree.GetListOfBranches()]
    for branch in ObjList:
      nameToSearch1 = "pho" + str(int(dp1_p1i)+1) + "_"
      if branch.startswith(nameToSearch1) :
        getattr(treeSkimmer, branch)[0] = getattr(itree, branch)


    if evt == eventsToRun :
      itree.GetEntry(0)
    else :
      itree.GetEntry(evt+1)

    ObjList = [key.GetName() for key in  itree.GetListOfBranches()]
    for branch in ObjList:
      nameToSearch1 = "pho" + str(int(dp1_p2i)+1) + "_"
      if branch.startswith(nameToSearch1) :
        getattr(treeSkimmer, branch)[0] = getattr(itree, branch)


    if evt == eventsToRun :
      itree.GetEntry(1)
    elif evt == (eventsToRun-1) :
      itree.GetEntry(0)
    else :
      itree.GetEntry(evt+2)

    ObjList = [key.GetName() for key in  itree.GetListOfBranches()]
    for branch in ObjList:
      nameToSearch1 = "pho" + str(int(dp2_p1i)+1) + "_"
      if branch.startswith(nameToSearch1) :
        getattr(treeSkimmer, branch)[0] = getattr(itree, branch)


    if evt == eventsToRun :
      itree.GetEntry(2)
    elif evt == (eventsToRun-1) :
      itree.GetEntry(1)
    elif evt == (eventsToRun-2) :
      itree.GetEntry(0)
    else :
      itree.GetEntry(evt+3)

    ObjList = [key.GetName() for key in  itree.GetListOfBranches()]
    for branch in ObjList:
      nameToSearch1 = "pho" + str(int(dp2_p2i)+1) + "_"
      if branch.startswith(nameToSearch1) :
        getattr(treeSkimmer, branch)[0] = getattr(itree, branch)


    order_photons = {
          1: treeSkimmer.pho1_pt,
          2: treeSkimmer.pho2_pt,
          3: treeSkimmer.pho3_pt,
          4: treeSkimmer.pho4_pt
          }
    # print "treeSkimmer.pho1_pt ", treeSkimmer.pho1_pt
    # print "treeSkimmer.pho2_pt ", treeSkimmer.pho2_pt
    # print "treeSkimmer.pho3_pt ", treeSkimmer.pho3_pt
    # print "treeSkimmer.pho4_pt ", treeSkimmer.pho4_pt
    sorted_order_photons = sorted(order_photons.items(), key=operator.itemgetter(1))
    sorted_order_photons.reverse()   # from high to low


    sPhos = []
    sPhos_mva = []
    sPhos_full5x5_r9 = []
    sPhos_chHadIso = []
    sPhos_HoE = []
    sPhos_SCEta = []
    sPhos_PSV = []
    pho1_vec = []
    pho2_vec = []
    pho3_vec = []
    pho4_vec = []
    #for iphoton in range(0, 4):  # 4 photons
    for iphoton in sorted_order_photons:
      p4 = ROOT.TLorentzVector(0,0,0,0)
      name_pt  = "pho" + str( int(iphoton [0] ) ) + "_pt"
      name_eta = "pho" + str( int(iphoton [0] ) ) + "_eta"
      name_phi = "pho" + str( int(iphoton [0] ) ) + "_phi"
      name_mva = "pho" + str( int(iphoton [0] ) ) + "_MVA"
      name_full5x5_r9 = "pho" + str(int(iphoton [0])) + "_full5x5_r9"
      name_chHadIso = "pho" + str(int(iphoton [0])) + "_ChHadIso"
      name_HoE = "pho" + str(int(iphoton [0])) + "_HoE"
      name_SCEta = "pho" + str(int(iphoton [0])) + "_SC_eta"
      name_PSV = "pho" + str(int(iphoton [0])) + "_pixelseed"

      p4.SetPtEtaPhiM( getattr(treeSkimmer, name_pt),  getattr(treeSkimmer, name_eta),    getattr(treeSkimmer, name_phi) , 0 )
      sPhos.append(p4)
      sPhos_mva.append(getattr(treeSkimmer, name_mva))
      sPhos_full5x5_r9.append(getattr(treeSkimmer, name_full5x5_r9))
      sPhos_chHadIso.append(getattr(treeSkimmer, name_chHadIso))
      sPhos_HoE.append(getattr(treeSkimmer, name_HoE))
      sPhos_PSV.append(getattr(treeSkiimer, name_PSV))

    pho1_vec.append(sPhos[0])
    pho1_vec.append(sPhos_full5x5_r9[0])
    pho1_vec.append(sPhos_chHadIso[0])
    pho1_vec.append(sPhos_HoE[0])
    pho1_vec.append(sPhos_SCEta[0])
    pho1_vec.append(sPhos_PSV[0])

    pho2_vec.append(sPhos[1])
    pho2_vec.append(sPhos_full5x5_r9[1])
    pho2_vec.append(sPhos_chHadIso[1])
    pho2_vec.append(sPhos_HoE[1])
    pho2_vec.append(sPhos_SCEta[1])
    pho2_vec.append(sPhos_PSV[1])

    pho3_vec.append(sPhos[2])
    pho3_vec.append(sPhos_full5x5_r9[2])
    pho3_vec.append(sPhos_chHadIso[2])
    pho3_vec.append(sPhos_HoE[2])
    pho3_vec.append(sPhos_SCEta[2])
    pho3_vec.append(sPhos_PSV[2])

    pho4_vec.append(sPhos[3])
    pho4_vec.append(sPhos_full5x5_r9[3])
    pho4_vec.append(sPhos_chHadIso[3])
    pho4_vec.append(sPhos_HoE[3])
    pho4_vec.append(sPhos_SCEta[3])
    pho4_vec.append(sPhos_PSV[3])


    # print "pt ", sPhos[0].Pt()
    # print sPhos[1].Pt()
    # print sPhos[2].Pt()
    # print sPhos[3].Pt()
    # print sPhos_mva[0], "  ", sPhos_mva[1], "  ", sPhos_mva[2], "  ", sPhos_mva[3]

    pho1_pt_mix[0] = sPhos[0].Pt()
    pho2_pt_mix[0] = sPhos[1].Pt()
    pho3_pt_mix[0] = sPhos[2].Pt()
    pho4_pt_mix[0] = sPhos[3].Pt()
    pho1_eta_mix[0] = sPhos[0].Eta()
    pho2_eta_mix[0] = sPhos[1].Eta()
    pho3_eta_mix[0] = sPhos[2].Eta()
    pho4_eta_mix[0] = sPhos[3].Eta()
    pho1_MVA_mix[0] = sPhos_mva[0]
    pho2_MVA_mix[0] = sPhos_mva[1]
    pho3_MVA_mix[0] = sPhos_mva[2]
    pho4_MVA_mix[0] = sPhos_mva[3]
    pho_minMVA_mix[0] = min(sPhos_mva)
    pho_maxMVA_mix[0] = max(sPhos_mva)
    pho12_DR_mix[0] = sPhos[0].DeltaR(sPhos[1])
    pho13_DR_mix[0] = sPhos[0].DeltaR(sPhos[2])
    pho14_DR_mix[0] = sPhos[0].DeltaR(sPhos[3])
    pho23_DR_mix[0] = sPhos[1].DeltaR(sPhos[2])
    pho24_DR_mix[0] = sPhos[1].DeltaR(sPhos[3])
    pho34_DR_mix[0] = sPhos[2].DeltaR(sPhos[3])
    pho12_M_mix[0] = (sPhos[0]+sPhos[1]).M()
    pho13_M_mix[0] = (sPhos[0]+sPhos[2]).M()
    pho14_M_mix[0] = (sPhos[0]+sPhos[3]).M()
    pho23_M_mix[0] = (sPhos[1]+sPhos[2]).M()
    pho24_M_mix[0] = (sPhos[1]+sPhos[3]).M()
    pho34_M_mix[0] = (sPhos[2]+sPhos[3]).M()

    dipair_index_map = [[],[]]
    mvaVal = []
    for i1 in range(0,len(sPhos)):
        pho1 = sPhos[i1]
        for i2 in range(0,len(sPhos)):
            if (i2<=i1): continue
            pho2 = sPhos[i2]
            for i3 in range(0,len(sPhos)):
                if (i3<=i2): continue
                if (i3==i1): continue
                pho3 = sPhos[i3]
                for i4 in range(0,len(sPhos)):
                    if (i4<=i3): continue
                    if (i4 == i1 or i4 == i2): continue
                    pho4 = sPhos[i4]
                    # print i1, i2, i3, i4

                    dipho1 = pho1 + pho2
                    dipho1_energy = dipho1.E()
                    dipho1_pt = dipho1.Pt()
                    dipho1_eta = dipho1.Eta()
                    dipho1_Phi = dipho1.Phi()
                    dipho1_dR = pho1.DeltaR(pho2)
                    dipho1_mass = dipho1.M()
                    deltaM1_gen1 = abs(dipho1.M()-float(mass))
                    deltaM1_gen2 = abs(dipho1.M()-float(mass))

                    dipho2 = pho3 + pho4
                    dipho2_energy = dipho2.E()
                    dipho2_pt = dipho2.Pt()
                    dipho2_eta = dipho2.Eta()
                    dipho2_Phi = dipho2.Phi()
                    dipho2_dR = pho3.DeltaR(pho4)
                    dipho2_mass = dipho2.M()
                    deltaM2_gen1 = abs(dipho2.M()-float(mass))
                    deltaM2_gen2 = abs(dipho2.M()-float(mass))

                    dipair_dR = dipho1.DeltaR(dipho2)

                    a1_Energy_mix[0] = dipho1_energy
                    a1_Pt_mix[0] = dipho1_pt
                    a1_Eta_mix[0] = dipho1_eta
                    a1_DR_mix[0] = dipho1_dR
                    dM1_gen1_mix[0] = deltaM1_gen1
                    a2_Energy_mix[0] = dipho2_energy
                    a2_Pt_mix[0] = dipho2_pt
                    a2_Eta_mix[0] = dipho2_eta
                    a2_DR_mix[0] = dipho2_dR
                    dM2_gen1_mix[0] = deltaM2_gen1
                    a1_a2_DR_mix[0] = dipair_dR

                    mva_value = reader.EvaluateMVA( "BDT" )
                    mvaVal.append(mva_value)
                    dipair_index_map[0].append([i1,i2,i3,i4])
                    dipair_index_map[1].append(mva_value)

                    dipho1 = pho1 + pho3
                    dipho1_energy = dipho1.E()
                    dipho1_pt = dipho1.Pt()
                    dipho1_eta = dipho1.Eta()
                    dipho1_Phi = dipho1.Phi()
                    dipho1_dR = pho1.DeltaR(pho3)
                    dipho1_mass = dipho1.M()
                    deltaM1_gen1 = abs(dipho1.M()-float(mass))
                    deltaM1_gen2 = abs(dipho1.M()-float(mass))

                    dipho2 = pho2 + pho4
                    dipho2_energy = dipho2.E()
                    dipho2_pt = dipho2.Pt()
                    dipho2_eta = dipho2.Eta()
                    dipho2_Phi = dipho2.Phi()
                    dipho2_dR = pho2.DeltaR(pho4)
                    dipho2_mass = dipho2.M()
                    deltaM2_gen1 = abs(dipho2.M()-float(mass))
                    deltaM2_gen2 = abs(dipho2.M()-float(mass))

                    dipair_dR = dipho1.DeltaR(dipho2)

                    a1_Energy_mix[0] = dipho1_energy
                    a1_Pt_mix[0] = dipho1_pt
                    a1_Eta_mix[0] = dipho1_eta
                    a1_DR_mix[0] = dipho1_dR
                    dM1_gen1_mix[0] = deltaM1_gen1
                    a2_Energy_mix[0] = dipho2_energy
                    a2_Pt_mix[0] = dipho2_pt
                    a2_Eta_mix[0] = dipho2_eta
                    a2_DR_mix[0] = dipho2_dR
                    dM2_gen1_mix[0] = deltaM2_gen1
                    a1_a2_DR_mix[0] = dipair_dR
                    mva_value = reader.EvaluateMVA( "BDT" )
                    # print "mva 2 ", mva_value
                    mvaVal.append(mva_value)
                    dipair_index_map[0].append([i1,i3,i2,i4])
                    dipair_index_map[1].append(mva_value)

                    dipho1 = pho1 + pho4
                    dipho1_energy = dipho1.E()
                    dipho1_pt = dipho1.Pt()
                    dipho1_eta = dipho1.Eta()
                    dipho1_Phi = dipho1.Phi()
                    dipho1_dR = pho1.DeltaR(pho4)
                    dipho1_mass = dipho1.M()
                    deltaM1_gen1 = abs(dipho1.M()-float(mass))
                    deltaM1_gen2 = abs(dipho1.M()-float(mass))

                    dipho2 = pho2 + pho3
                    dipho2_energy = dipho2.E()
                    dipho2_pt = dipho2.Pt()
                    dipho2_eta = dipho2.Eta()
                    dipho2_Phi = dipho2.Phi()
                    dipho2_dR = pho2.DeltaR(pho3)
                    dipho2_mass = dipho2.M()
                    deltaM2_gen1 = abs(dipho2.M()-float(mass))
                    deltaM2_gen2 = abs(dipho2.M()-float(mass))

                    dipair_dR = dipho1.DeltaR(dipho2)

                    a1_Energy_mix[0] = dipho1_energy
                    a1_Pt_mix[0] = dipho1_pt
                    a1_Eta_mix[0] = dipho1_eta
                    a1_DR_mix[0] = dipho1_dR
                    dM1_gen1_mix[0] = deltaM1_gen1
                    a2_Energy_mix[0] = dipho2_energy
                    a2_Pt_mix[0] = dipho2_pt
                    a2_Eta_mix[0] = dipho2_eta
                    a2_DR_mix[0] = dipho2_dR
                    dM2_gen1_mix[0] = deltaM2_gen1
                    a1_a2_DR_mix[0] = dipair_dR
                    mva_value = reader.EvaluateMVA( "BDT" )
                    # print "mva 3 ", mva_value
                    mvaVal.append(mva_value)
                    dipair_index_map[0].append([i1,i4,i2,i3])
                    dipair_index_map[1].append(mva_value)

    indexlist = dipair_index_map[0]
    mvalist = dipair_index_map[1]
    # print mvalist

    pairMVAscore_mix[0] = max(mvalist)

    Pho1 = sPhos[indexlist[max(xrange(len(mvalist)), key=mvalist.__getitem__)][0]]
    Pho2 = sPhos[indexlist[max(xrange(len(mvalist)), key=mvalist.__getitem__)][1]]
    Pho3 = sPhos[indexlist[max(xrange(len(mvalist)), key=mvalist.__getitem__)][2]]
    Pho4 = sPhos[indexlist[max(xrange(len(mvalist)), key=mvalist.__getitem__)][3]]

    if (Pho1.Pt() > Pho2.Pt()):
        a1_pho1 = Pho1
        a1_pho2 = Pho2
        a1_ipho1 = indexlist[max(xrange(len(mvalist)), key=mvalist.__getitem__)][0]
        a1_ipho2 = indexlist[max(xrange(len(mvalist)), key=mvalist.__getitem__)][1]
    else:
        a1_pho1 = Pho2
        a1_pho2 = Pho1
        a1_ipho1 = indexlist[max(xrange(len(mvalist)), key=mvalist.__getitem__)][1]
        a1_ipho2 = indexlist[max(xrange(len(mvalist)), key=mvalist.__getitem__)][0]
    if (Pho3.Pt() > Pho4.Pt()):
        a2_pho1 = Pho3
        a2_pho2 = Pho4
        a2_ipho1 = indexlist[max(xrange(len(mvalist)), key=mvalist.__getitem__)][2]
        a2_ipho2 = indexlist[max(xrange(len(mvalist)), key=mvalist.__getitem__)][3]
    else:
        a2_pho1 = Pho4
        a2_pho2 = Pho3
        a2_ipho1 = indexlist[max(xrange(len(mvalist)), key=mvalist.__getitem__)][3]
        a2_ipho2 = indexlist[max(xrange(len(mvalist)), key=mvalist.__getitem__)][2]

    a1_temp = Pho1 + Pho2
    a2_temp = Pho3 + Pho4
    if (a1_temp.Pt() > a2_temp.Pt()):
          a1 = a1_temp
          a2 = a2_temp
    elif (a1_temp.Pt() < a2_temp.Pt()):
          a1 = a2_temp
          a2 = a1_temp

    a1_mass_mix[0] = a1.M()
    a2_mass_mix[0] = a2.M()
    a1_Pt_mix[0] = a1.Pt()
    a2_Pt_mix[0] = a2.Pt()
    a1_Eta_mix[0] = a1.Eta()
    a2_Eta_mix[0] = a2.Eta()
    a1_Phi_mix[0] = a1.Phi()
    a2_Phi_mix[0] = a2.Phi()
    a1_Energy_mix[0] = a1.Energy()
    a2_Energy_mix[0] = a2.Energy()
    a1_DR_mix[0] = a1_pho1.DeltaR(a1_pho2)
    a2_DR_mix[0] = a2_pho1.DeltaR(a2_pho2)
    a1_a2_DR_mix[0] = a1.DeltaR(a2)
    a1_iPho1_mix[0] = float(a1_ipho1)
    a1_iPho2_mix[0] = float(a1_ipho2)
    a2_iPho1_mix[0] = float(a2_ipho1)
    a2_iPho2_mix[0] = float(a2_ipho2)
    a1_PtOverMass_mix[0] = a1.Pt()/a1.M()
    a2_PtOverMass_mix[0] = a2.Pt()/a2.M()
    a1_Pho1PtOvera1Mass_mix[0] = a1_pho1.Pt()/a1.M()
    a1_Pho2PtOvera1Mass_mix[0] = a1_pho2.Pt()/a1.M()
    a2_Pho1PtOvera2Mass_mix[0] = a2_pho1.Pt()/a2.M()
    a2_Pho2PtOvera2Mass_mix[0] = a2_pho2.Pt()/a2.M()
    CTStarCS_mix[0] =  abs(getCosThetaStar_CS())
    CT_a1Pho1_mix[0] = abs(CosThetaAngles()[0])
    CT_a2Pho1_mix[0] = abs(CosThetaAngles()[1])

    Pgggg = sPhos[0] + sPhos[1] + sPhos[2] + sPhos[3]
    tp_mass_mix[0] = Pgggg.M()
    tp_pt_mix[0] = Pgggg.Pt()
    tp_eta_mix[0] = Pgggg.Eta()
    tp_phi_mix[0] = Pgggg.Phi()


    outtree.Fill()

  outRoot.cd()
  outtree.Write()
  outRoot.Close()
