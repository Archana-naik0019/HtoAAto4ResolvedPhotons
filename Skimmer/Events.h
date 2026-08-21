//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 20 13:47:55 2026 by ROOT version 6.38.00
// from TTree Events/
// found on file: root://eoshome-a.cern.ch//eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/NanoAOD_skimming/Skimmer/updated_test.root
//////////////////////////////////////////////////////////

#ifndef Events_h
#define Events_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Events {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   Bool_t          HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId;
   Bool_t          HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55;
   Float_t         Rho_fixedGridRhoAll;
   Int_t           nJet;
   Float_t         Jet_chEmEF[10];   //[nJet]
   Float_t         Jet_eta[10];   //[nJet]
   Float_t         Jet_neEmEF[10];   //[nJet]
   Float_t         Jet_phi[10];   //[nJet]
   Float_t         Jet_pt[10];   //[nJet]
   Int_t           nPhoton;
   UChar_t         Photon_cutBased[6];   //[nPhoton]
   Bool_t          Photon_electronVeto[6];   //[nPhoton]
   Bool_t          Photon_hasConversionTracks[6];   //[nPhoton]
   Bool_t          Photon_isScEtaEB[6];   //[nPhoton]
   Bool_t          Photon_isScEtaEE[6];   //[nPhoton]
   Bool_t          Photon_mvaID_WP80[6];   //[nPhoton]
   Bool_t          Photon_mvaID_WP90[6];   //[nPhoton]
   Bool_t          Photon_pixelSeed[6];   //[nPhoton]
   UChar_t         Photon_seedGain[6];   //[nPhoton]
   Short_t         Photon_electronIdx[6];   //[nPhoton]
   Short_t         Photon_jetIdx[6];   //[nPhoton]
   Short_t         Photon_seediEtaOriX[6];   //[nPhoton]
   Short_t         Photon_seediPhiOriY[6];   //[nPhoton]
   Int_t           Photon_vidNestedWPBitmap[6];   //[nPhoton]
   Float_t         Photon_ecalPFClusterIso[6];   //[nPhoton]
   Float_t         Photon_energyErr[6];   //[nPhoton]
   Float_t         Photon_energyRaw[6];   //[nPhoton]
   Float_t         Photon_esEffSigmaRR[6];   //[nPhoton]
   Float_t         Photon_esEnergyOverRawE[6];   //[nPhoton]
   Float_t         Photon_eta[6];   //[nPhoton]
   Float_t         Photon_etaWidth[6];   //[nPhoton]
   Float_t         Photon_haloTaggerMVAVal[6];   //[nPhoton]
   Float_t         Photon_hcalPFClusterIso[6];   //[nPhoton]
   Float_t         Photon_hoe[6];   //[nPhoton]
   Float_t         Photon_hoe_PUcorr[6];   //[nPhoton]
   Float_t         Photon_hoe_Tower[6];   //[nPhoton]
   Float_t         Photon_mvaID[6];   //[nPhoton]
   Float_t         Photon_pfChargedIso[6];   //[nPhoton]
   Float_t         Photon_pfChargedIsoPFPV[6];   //[nPhoton]
   Float_t         Photon_pfChargedIsoWorstVtx[6];   //[nPhoton]
   Float_t         Photon_pfPhoIso03[6];   //[nPhoton]
   Float_t         Photon_pfRelIso03_all_quadratic[6];   //[nPhoton]
   Float_t         Photon_pfRelIso03_chg_quadratic[6];   //[nPhoton]
   Float_t         Photon_phi[6];   //[nPhoton]
   Float_t         Photon_phiWidth[6];   //[nPhoton]
   Float_t         Photon_pt[6];   //[nPhoton]
   Float_t         Photon_r9[6];   //[nPhoton]
   Float_t         Photon_s4[6];   //[nPhoton]
   Float_t         Photon_sieie[6];   //[nPhoton]
   Float_t         Photon_sieip[6];   //[nPhoton]
   Float_t         Photon_sipip[6];   //[nPhoton]
   Float_t         Photon_superclusterEta[6];   //[nPhoton]
   Float_t         Photon_trkSumPtHollowConeDR03[6];   //[nPhoton]
   Float_t         Photon_trkSumPtSolidConeDR04[6];   //[nPhoton]
   Float_t         Photon_x_calo[6];   //[nPhoton]
   Float_t         Photon_y_calo[6];   //[nPhoton]
   Float_t         Photon_z_calo[6];   //[nPhoton]
   Int_t           nElectron;
   Bool_t          Electron_convVeto[3];   //[nElectron]
   UChar_t         Electron_cutBased[3];   //[nElectron]
   Bool_t          Electron_cutBased_HEEP[3];   //[nElectron]
   Bool_t          Electron_isEB[3];   //[nElectron]
   Bool_t          Electron_isEcalDriven[3];   //[nElectron]
   Bool_t          Electron_isPFcand[3];   //[nElectron]
   UChar_t         Electron_jetNDauCharged[3];   //[nElectron]
   UChar_t         Electron_lostHits[3];   //[nElectron]
   Bool_t          Electron_mvaIso_WP80[3];   //[nElectron]
   Bool_t          Electron_mvaIso_WP90[3];   //[nElectron]
   Bool_t          Electron_mvaIso_WPHZZ[3];   //[nElectron]
   Bool_t          Electron_mvaNoIso_WP80[3];   //[nElectron]
   Bool_t          Electron_mvaNoIso_WP90[3];   //[nElectron]
   UChar_t         Electron_seedGain[3];   //[nElectron]
   UChar_t         Electron_tightCharge[3];   //[nElectron]
   Short_t         Electron_jetIdx[3];   //[nElectron]
   Short_t         Electron_photonIdx[3];   //[nElectron]
   Short_t         Electron_seediEtaOriX[3];   //[nElectron]
   Short_t         Electron_seediPhiOriY[3];   //[nElectron]
   Short_t         Electron_svIdx[3];   //[nElectron]
   Short_t         Electron_fsrPhotonIdx[3];   //[nElectron]
   Int_t           Electron_charge[3];   //[nElectron]
   Int_t           Electron_pdgId[3];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmap[3];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmapHEEP[3];   //[nElectron]
   Float_t         Electron_PreshowerEnergy[3];   //[nElectron]
   Float_t         Electron_deltaEtaSC[3];   //[nElectron]
   Float_t         Electron_dr03EcalRecHitSumEt[3];   //[nElectron]
   Float_t         Electron_dr03HcalDepth1TowerSumEt[3];   //[nElectron]
   Float_t         Electron_dr03TkSumPt[3];   //[nElectron]
   Float_t         Electron_dr03TkSumPtHEEP[3];   //[nElectron]
   Float_t         Electron_dxy[3];   //[nElectron]
   Float_t         Electron_dxyErr[3];   //[nElectron]
   Float_t         Electron_dz[3];   //[nElectron]
   Float_t         Electron_dzErr[3];   //[nElectron]
   Float_t         Electron_eInvMinusPInv[3];   //[nElectron]
   Float_t         Electron_ecalEnergy[3];   //[nElectron]
   Float_t         Electron_ecalEnergyError[3];   //[nElectron]
   Float_t         Electron_energyErr[3];   //[nElectron]
   Float_t         Electron_eta[3];   //[nElectron]
   Float_t         Electron_fbrem[3];   //[nElectron]
   Float_t         Electron_gsfTrketaMode[3];   //[nElectron]
   Float_t         Electron_gsfTrkpMode[3];   //[nElectron]
   Float_t         Electron_gsfTrkpModeErr[3];   //[nElectron]
   Float_t         Electron_gsfTrkphiMode[3];   //[nElectron]
   Float_t         Electron_hoe[3];   //[nElectron]
   Float_t         Electron_ip3d[3];   //[nElectron]
   Float_t         Electron_jetDF[3];   //[nElectron]
   Float_t         Electron_jetPtRelv2[3];   //[nElectron]
   Float_t         Electron_jetRelIso[3];   //[nElectron]
   Float_t         Electron_mass[3];   //[nElectron]
   Float_t         Electron_miniPFRelIso_all[3];   //[nElectron]
   Float_t         Electron_miniPFRelIso_chg[3];   //[nElectron]
   Float_t         Electron_mvaHZZIso[3];   //[nElectron]
   Float_t         Electron_mvaIso[3];   //[nElectron]
   Float_t         Electron_mvaNoIso[3];   //[nElectron]
   Float_t         Electron_pfRelIso03_all[3];   //[nElectron]
   Float_t         Electron_pfRelIso03_chg[3];   //[nElectron]
   Float_t         Electron_pfRelIso04_all[3];   //[nElectron]
   Float_t         Electron_phi[3];   //[nElectron]
   Float_t         Electron_pt[3];   //[nElectron]
   Float_t         Electron_r9[3];   //[nElectron]
   Float_t         Electron_rawEnergy[3];   //[nElectron]
   Float_t         Electron_scEtOverPt[3];   //[nElectron]
   Float_t         Electron_sieie[3];   //[nElectron]
   Float_t         Electron_sip3d[3];   //[nElectron]
   Float_t         Electron_superclusterEta[3];   //[nElectron]
   Float_t         Electron_promptMVA[3];   //[nElectron]
   Float_t         Electron_IPx[3];   //[nElectron]
   Float_t         Electron_IPy[3];   //[nElectron]
   Float_t         Electron_IPz[3];   //[nElectron]
   Float_t         Electron_ipLengthSig[3];   //[nElectron]
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_pt;
   UChar_t         PV_npvs;
   UChar_t         PV_npvsGood;
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Float_t         PV_sumpt2;
   Float_t         PV_sumpx;
   Float_t         PV_sumpy;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId;   //!
   TBranch        *b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55;   //!
   TBranch        *b_Rho_fixedGridRhoAll;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_neEmEF;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_nPhoton;   //!
   TBranch        *b_Photon_cutBased;   //!
   TBranch        *b_Photon_electronVeto;   //!
   TBranch        *b_Photon_hasConversionTracks;   //!
   TBranch        *b_Photon_isScEtaEB;   //!
   TBranch        *b_Photon_isScEtaEE;   //!
   TBranch        *b_Photon_mvaID_WP80;   //!
   TBranch        *b_Photon_mvaID_WP90;   //!
   TBranch        *b_Photon_pixelSeed;   //!
   TBranch        *b_Photon_seedGain;   //!
   TBranch        *b_Photon_electronIdx;   //!
   TBranch        *b_Photon_jetIdx;   //!
   TBranch        *b_Photon_seediEtaOriX;   //!
   TBranch        *b_Photon_seediPhiOriY;   //!
   TBranch        *b_Photon_vidNestedWPBitmap;   //!
   TBranch        *b_Photon_ecalPFClusterIso;   //!
   TBranch        *b_Photon_energyErr;   //!
   TBranch        *b_Photon_energyRaw;   //!
   TBranch        *b_Photon_esEffSigmaRR;   //!
   TBranch        *b_Photon_esEnergyOverRawE;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_etaWidth;   //!
   TBranch        *b_Photon_haloTaggerMVAVal;   //!
   TBranch        *b_Photon_hcalPFClusterIso;   //!
   TBranch        *b_Photon_hoe;   //!
   TBranch        *b_Photon_hoe_PUcorr;   //!
   TBranch        *b_Photon_hoe_Tower;   //!
   TBranch        *b_Photon_mvaID;   //!
   TBranch        *b_Photon_pfChargedIso;   //!
   TBranch        *b_Photon_pfChargedIsoPFPV;   //!
   TBranch        *b_Photon_pfChargedIsoWorstVtx;   //!
   TBranch        *b_Photon_pfPhoIso03;   //!
   TBranch        *b_Photon_pfRelIso03_all_quadratic;   //!
   TBranch        *b_Photon_pfRelIso03_chg_quadratic;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_phiWidth;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_s4;   //!
   TBranch        *b_Photon_sieie;   //!
   TBranch        *b_Photon_sieip;   //!
   TBranch        *b_Photon_sipip;   //!
   TBranch        *b_Photon_superclusterEta;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR03;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR04;   //!
   TBranch        *b_Photon_x_calo;   //!
   TBranch        *b_Photon_y_calo;   //!
   TBranch        *b_Photon_z_calo;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_Electron_convVeto;   //!
   TBranch        *b_Electron_cutBased;   //!
   TBranch        *b_Electron_cutBased_HEEP;   //!
   TBranch        *b_Electron_isEB;   //!
   TBranch        *b_Electron_isEcalDriven;   //!
   TBranch        *b_Electron_isPFcand;   //!
   TBranch        *b_Electron_jetNDauCharged;   //!
   TBranch        *b_Electron_lostHits;   //!
   TBranch        *b_Electron_mvaIso_WP80;   //!
   TBranch        *b_Electron_mvaIso_WP90;   //!
   TBranch        *b_Electron_mvaIso_WPHZZ;   //!
   TBranch        *b_Electron_mvaNoIso_WP80;   //!
   TBranch        *b_Electron_mvaNoIso_WP90;   //!
   TBranch        *b_Electron_seedGain;   //!
   TBranch        *b_Electron_tightCharge;   //!
   TBranch        *b_Electron_jetIdx;   //!
   TBranch        *b_Electron_photonIdx;   //!
   TBranch        *b_Electron_seediEtaOriX;   //!
   TBranch        *b_Electron_seediPhiOriY;   //!
   TBranch        *b_Electron_svIdx;   //!
   TBranch        *b_Electron_fsrPhotonIdx;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_pdgId;   //!
   TBranch        *b_Electron_vidNestedWPBitmap;   //!
   TBranch        *b_Electron_vidNestedWPBitmapHEEP;   //!
   TBranch        *b_Electron_PreshowerEnergy;   //!
   TBranch        *b_Electron_deltaEtaSC;   //!
   TBranch        *b_Electron_dr03EcalRecHitSumEt;   //!
   TBranch        *b_Electron_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_Electron_dr03TkSumPt;   //!
   TBranch        *b_Electron_dr03TkSumPtHEEP;   //!
   TBranch        *b_Electron_dxy;   //!
   TBranch        *b_Electron_dxyErr;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Electron_dzErr;   //!
   TBranch        *b_Electron_eInvMinusPInv;   //!
   TBranch        *b_Electron_ecalEnergy;   //!
   TBranch        *b_Electron_ecalEnergyError;   //!
   TBranch        *b_Electron_energyErr;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_fbrem;   //!
   TBranch        *b_Electron_gsfTrketaMode;   //!
   TBranch        *b_Electron_gsfTrkpMode;   //!
   TBranch        *b_Electron_gsfTrkpModeErr;   //!
   TBranch        *b_Electron_gsfTrkphiMode;   //!
   TBranch        *b_Electron_hoe;   //!
   TBranch        *b_Electron_ip3d;   //!
   TBranch        *b_Electron_jetDF;   //!
   TBranch        *b_Electron_jetPtRelv2;   //!
   TBranch        *b_Electron_jetRelIso;   //!
   TBranch        *b_Electron_mass;   //!
   TBranch        *b_Electron_miniPFRelIso_all;   //!
   TBranch        *b_Electron_miniPFRelIso_chg;   //!
   TBranch        *b_Electron_mvaHZZIso;   //!
   TBranch        *b_Electron_mvaIso;   //!
   TBranch        *b_Electron_mvaNoIso;   //!
   TBranch        *b_Electron_pfRelIso03_all;   //!
   TBranch        *b_Electron_pfRelIso03_chg;   //!
   TBranch        *b_Electron_pfRelIso04_all;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_r9;   //!
   TBranch        *b_Electron_rawEnergy;   //!
   TBranch        *b_Electron_scEtOverPt;   //!
   TBranch        *b_Electron_sieie;   //!
   TBranch        *b_Electron_sip3d;   //!
   TBranch        *b_Electron_superclusterEta;   //!
   TBranch        *b_Electron_promptMVA;   //!
   TBranch        *b_Electron_IPx;   //!
   TBranch        *b_Electron_IPy;   //!
   TBranch        *b_Electron_IPz;   //!
   TBranch        *b_Electron_ipLengthSig;   //!
   TBranch        *b_PuppiMET_phi;   //!
   TBranch        *b_PuppiMET_pt;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_sumpt2;   //!
   TBranch        *b_PV_sumpx;   //!
   TBranch        *b_PV_sumpy;   //!

   Events(TTree *tree=0);
   virtual ~Events();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Events_cxx
Events::Events(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoshome-a.cern.ch//eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/NanoAOD_skimming/Skimmer/updated_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoshome-a.cern.ch//eos/user/a/arnaik/Higgs_AA_3Photons_analysis/Data_2024/resolved_4photons/NanoAOD_skimming/Skimmer/updated_test.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
}

Events::~Events()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Events::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Events::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Events::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId", &HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId, &b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId);
   fChain->SetBranchAddress("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55", &HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55, &b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55);
   fChain->SetBranchAddress("Rho_fixedGridRhoAll", &Rho_fixedGridRhoAll, &b_Rho_fixedGridRhoAll);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   fChain->SetBranchAddress("Photon_cutBased", Photon_cutBased, &b_Photon_cutBased);
   fChain->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);
   fChain->SetBranchAddress("Photon_hasConversionTracks", Photon_hasConversionTracks, &b_Photon_hasConversionTracks);
   fChain->SetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB, &b_Photon_isScEtaEB);
   fChain->SetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE, &b_Photon_isScEtaEE);
   fChain->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80);
   fChain->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90);
   fChain->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);
   fChain->SetBranchAddress("Photon_seedGain", Photon_seedGain, &b_Photon_seedGain);
   fChain->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);
   fChain->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);
   fChain->SetBranchAddress("Photon_seediEtaOriX", Photon_seediEtaOriX, &b_Photon_seediEtaOriX);
   fChain->SetBranchAddress("Photon_seediPhiOriY", Photon_seediPhiOriY, &b_Photon_seediPhiOriY);
   fChain->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);
   fChain->SetBranchAddress("Photon_ecalPFClusterIso", Photon_ecalPFClusterIso, &b_Photon_ecalPFClusterIso);
   fChain->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr);
   fChain->SetBranchAddress("Photon_energyRaw", Photon_energyRaw, &b_Photon_energyRaw);
   fChain->SetBranchAddress("Photon_esEffSigmaRR", Photon_esEffSigmaRR, &b_Photon_esEffSigmaRR);
   fChain->SetBranchAddress("Photon_esEnergyOverRawE", Photon_esEnergyOverRawE, &b_Photon_esEnergyOverRawE);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_etaWidth", Photon_etaWidth, &b_Photon_etaWidth);
   fChain->SetBranchAddress("Photon_haloTaggerMVAVal", Photon_haloTaggerMVAVal, &b_Photon_haloTaggerMVAVal);
   fChain->SetBranchAddress("Photon_hcalPFClusterIso", Photon_hcalPFClusterIso, &b_Photon_hcalPFClusterIso);
   fChain->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe);
   fChain->SetBranchAddress("Photon_hoe_PUcorr", Photon_hoe_PUcorr, &b_Photon_hoe_PUcorr);
   fChain->SetBranchAddress("Photon_hoe_Tower", Photon_hoe_Tower, &b_Photon_hoe_Tower);
   fChain->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID);
   fChain->SetBranchAddress("Photon_pfChargedIso", Photon_pfChargedIso, &b_Photon_pfChargedIso);
   fChain->SetBranchAddress("Photon_pfChargedIsoPFPV", Photon_pfChargedIsoPFPV, &b_Photon_pfChargedIsoPFPV);
   fChain->SetBranchAddress("Photon_pfChargedIsoWorstVtx", Photon_pfChargedIsoWorstVtx, &b_Photon_pfChargedIsoWorstVtx);
   fChain->SetBranchAddress("Photon_pfPhoIso03", Photon_pfPhoIso03, &b_Photon_pfPhoIso03);
   fChain->SetBranchAddress("Photon_pfRelIso03_all_quadratic", Photon_pfRelIso03_all_quadratic, &b_Photon_pfRelIso03_all_quadratic);
   fChain->SetBranchAddress("Photon_pfRelIso03_chg_quadratic", Photon_pfRelIso03_chg_quadratic, &b_Photon_pfRelIso03_chg_quadratic);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_phiWidth", Photon_phiWidth, &b_Photon_phiWidth);
   fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon_s4", Photon_s4, &b_Photon_s4);
   fChain->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie);
   fChain->SetBranchAddress("Photon_sieip", Photon_sieip, &b_Photon_sieip);
   fChain->SetBranchAddress("Photon_sipip", Photon_sipip, &b_Photon_sipip);
   fChain->SetBranchAddress("Photon_superclusterEta", Photon_superclusterEta, &b_Photon_superclusterEta);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", Photon_trkSumPtHollowConeDR03, &b_Photon_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", Photon_trkSumPtSolidConeDR04, &b_Photon_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("Photon_x_calo", Photon_x_calo, &b_Photon_x_calo);
   fChain->SetBranchAddress("Photon_y_calo", Photon_y_calo, &b_Photon_y_calo);
   fChain->SetBranchAddress("Photon_z_calo", Photon_z_calo, &b_Photon_z_calo);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
   fChain->SetBranchAddress("Electron_cutBased_HEEP", Electron_cutBased_HEEP, &b_Electron_cutBased_HEEP);
   fChain->SetBranchAddress("Electron_isEB", Electron_isEB, &b_Electron_isEB);
   fChain->SetBranchAddress("Electron_isEcalDriven", Electron_isEcalDriven, &b_Electron_isEcalDriven);
   fChain->SetBranchAddress("Electron_isPFcand", Electron_isPFcand, &b_Electron_isPFcand);
   fChain->SetBranchAddress("Electron_jetNDauCharged", Electron_jetNDauCharged, &b_Electron_jetNDauCharged);
   fChain->SetBranchAddress("Electron_lostHits", Electron_lostHits, &b_Electron_lostHits);
   fChain->SetBranchAddress("Electron_mvaIso_WP80", Electron_mvaIso_WP80, &b_Electron_mvaIso_WP80);
   fChain->SetBranchAddress("Electron_mvaIso_WP90", Electron_mvaIso_WP90, &b_Electron_mvaIso_WP90);
   fChain->SetBranchAddress("Electron_mvaIso_WPHZZ", Electron_mvaIso_WPHZZ, &b_Electron_mvaIso_WPHZZ);
   fChain->SetBranchAddress("Electron_mvaNoIso_WP80", Electron_mvaNoIso_WP80, &b_Electron_mvaNoIso_WP80);
   fChain->SetBranchAddress("Electron_mvaNoIso_WP90", Electron_mvaNoIso_WP90, &b_Electron_mvaNoIso_WP90);
   fChain->SetBranchAddress("Electron_seedGain", Electron_seedGain, &b_Electron_seedGain);
   fChain->SetBranchAddress("Electron_tightCharge", Electron_tightCharge, &b_Electron_tightCharge);
   fChain->SetBranchAddress("Electron_jetIdx", Electron_jetIdx, &b_Electron_jetIdx);
   fChain->SetBranchAddress("Electron_photonIdx", Electron_photonIdx, &b_Electron_photonIdx);
   fChain->SetBranchAddress("Electron_seediEtaOriX", Electron_seediEtaOriX, &b_Electron_seediEtaOriX);
   fChain->SetBranchAddress("Electron_seediPhiOriY", Electron_seediPhiOriY, &b_Electron_seediPhiOriY);
   fChain->SetBranchAddress("Electron_svIdx", Electron_svIdx, &b_Electron_svIdx);
   fChain->SetBranchAddress("Electron_fsrPhotonIdx", Electron_fsrPhotonIdx, &b_Electron_fsrPhotonIdx);
   fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmap", Electron_vidNestedWPBitmap, &b_Electron_vidNestedWPBitmap);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmapHEEP", Electron_vidNestedWPBitmapHEEP, &b_Electron_vidNestedWPBitmapHEEP);
   fChain->SetBranchAddress("Electron_PreshowerEnergy", Electron_PreshowerEnergy, &b_Electron_PreshowerEnergy);
   fChain->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
   fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electron_dr03TkSumPt", Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
   fChain->SetBranchAddress("Electron_dr03TkSumPtHEEP", Electron_dr03TkSumPtHEEP, &b_Electron_dr03TkSumPtHEEP);
   fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_eInvMinusPInv", Electron_eInvMinusPInv, &b_Electron_eInvMinusPInv);
   fChain->SetBranchAddress("Electron_ecalEnergy", Electron_ecalEnergy, &b_Electron_ecalEnergy);
   fChain->SetBranchAddress("Electron_ecalEnergyError", Electron_ecalEnergyError, &b_Electron_ecalEnergyError);
   fChain->SetBranchAddress("Electron_energyErr", Electron_energyErr, &b_Electron_energyErr);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_fbrem", Electron_fbrem, &b_Electron_fbrem);
   fChain->SetBranchAddress("Electron_gsfTrketaMode", Electron_gsfTrketaMode, &b_Electron_gsfTrketaMode);
   fChain->SetBranchAddress("Electron_gsfTrkpMode", Electron_gsfTrkpMode, &b_Electron_gsfTrkpMode);
   fChain->SetBranchAddress("Electron_gsfTrkpModeErr", Electron_gsfTrkpModeErr, &b_Electron_gsfTrkpModeErr);
   fChain->SetBranchAddress("Electron_gsfTrkphiMode", Electron_gsfTrkphiMode, &b_Electron_gsfTrkphiMode);
   fChain->SetBranchAddress("Electron_hoe", Electron_hoe, &b_Electron_hoe);
   fChain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_jetDF", Electron_jetDF, &b_Electron_jetDF);
   fChain->SetBranchAddress("Electron_jetPtRelv2", Electron_jetPtRelv2, &b_Electron_jetPtRelv2);
   fChain->SetBranchAddress("Electron_jetRelIso", Electron_jetRelIso, &b_Electron_jetRelIso);
   fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
   fChain->SetBranchAddress("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
   fChain->SetBranchAddress("Electron_mvaHZZIso", Electron_mvaHZZIso, &b_Electron_mvaHZZIso);
   fChain->SetBranchAddress("Electron_mvaIso", Electron_mvaIso, &b_Electron_mvaIso);
   fChain->SetBranchAddress("Electron_mvaNoIso", Electron_mvaNoIso, &b_Electron_mvaNoIso);
   fChain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
   fChain->SetBranchAddress("Electron_pfRelIso04_all", Electron_pfRelIso04_all, &b_Electron_pfRelIso04_all);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_r9", Electron_r9, &b_Electron_r9);
   fChain->SetBranchAddress("Electron_rawEnergy", Electron_rawEnergy, &b_Electron_rawEnergy);
   fChain->SetBranchAddress("Electron_scEtOverPt", Electron_scEtOverPt, &b_Electron_scEtOverPt);
   fChain->SetBranchAddress("Electron_sieie", Electron_sieie, &b_Electron_sieie);
   fChain->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_superclusterEta", Electron_superclusterEta, &b_Electron_superclusterEta);
   fChain->SetBranchAddress("Electron_promptMVA", Electron_promptMVA, &b_Electron_promptMVA);
   fChain->SetBranchAddress("Electron_IPx", Electron_IPx, &b_Electron_IPx);
   fChain->SetBranchAddress("Electron_IPy", Electron_IPy, &b_Electron_IPy);
   fChain->SetBranchAddress("Electron_IPz", Electron_IPz, &b_Electron_IPz);
   fChain->SetBranchAddress("Electron_ipLengthSig", Electron_ipLengthSig, &b_Electron_ipLengthSig);
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_sumpt2", &PV_sumpt2, &b_PV_sumpt2);
   fChain->SetBranchAddress("PV_sumpx", &PV_sumpx, &b_PV_sumpx);
   fChain->SetBranchAddress("PV_sumpy", &PV_sumpy, &b_PV_sumpy);
   Notify();
}

bool Events::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be for a new TTree in a TChain. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void Events::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Events::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Events_cxx
