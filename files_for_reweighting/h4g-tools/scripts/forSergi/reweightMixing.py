import ROOT
import itertools
import argparse
import os
import math
from array import array

def findVal(trre,var):
   val = -999.
   if(var=="cosThetaStarCS_dM"): val = trre.cosThetaStarCS_dM
   elif(var=="cosTheta_a1_dM"): val = trre.cosTheta_a1_dM
   elif(var=="cosTheta_a2_dM"): val = trre.cosTheta_a2_dM 
   elif(var=="a1_mass_dM"): val = trre.a1_mass_dM
   elif(var=="a2_mass_dM"): val = trre.a2_mass_dM
   elif(var=="a1_mass_dM_M_a2_mass_dM"): val = trre.a1_mass_dM-trre.a2_mass_dM
   elif(var=="a1_pt_dM"): val = trre.a1_pt_dM
   elif(var=="a2_pt_dM"): val = trre.a2_pt_dM
   elif(var=="a1_eta_dM"): val = trre.a1_eta_dM
   elif(var=="a2_eta_dM"): val = trre.a2_eta_dM
   elif(var=="a1_dR_dM"): val = trre.a1_dR_dM
   elif(var=="a2_dR_dM"): val = trre.a2_dR_dM
   elif(var=="a1_a2_dR_dM"): val = trre.a1_a2_dR_dM
   elif(var=="tp_pt"): val = trre.tp_pt
   elif(var=="tp_eta"): val = trre.tp_eta
   elif(var=="pho1_MVA"): 
        if trre.pho1_MVA<-1.: val = -1.1
        else: val = trre.pho1_MVA
   elif(var=="pho2_MVA"): 
        if trre.pho2_MVA<-1.: val = -1.1
        else: val = trre.pho2_MVA
   elif(var=="pho3_MVA"): 
        if trre.pho3_MVA<-1.: val = -1.1
        else: val = trre.pho3_MVA
   elif(var=="pho4_MVA"): 
        if trre.pho4_MVA<-1.: val = -1.1
        else: val = trre.pho4_MVA
   elif(var=="pho1_pt"): val = trre.pho1_pt 
   elif(var=="pho2_pt"): val = trre.pho2_pt 
   elif(var=="pho3_pt"): val = trre.pho3_pt 
   elif(var=="pho4_pt"): val = trre.pho4_pt 
   elif(var=="pho1_eta"): val = trre.pho1_eta  
   elif(var=="pho2_eta"): val = trre.pho2_eta 
   elif(var=="pho3_eta"): val = trre.pho3_eta 
   elif(var=="pho4_eta"): val = trre.pho4_eta 
   else: print "findVal ---> WARNING MISSING VAR: ",var

   return val

if __name__ == '__main__':

 ROOT.gROOT.SetBatch(ROOT.kTRUE)

 parser =  argparse.ArgumentParser(description='cat MVA')
 parser.add_argument('-d', '--data', dest='data', required=True, type=str)
 parser.add_argument('-m', '--mix',  dest='mix',  required=True, type=str)
 args = parser.parse_args()

 data = args.data+"/tagsDumper/trees/Data_13TeV_H4GTag_0"
 #data = '/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2016/hadd/data_2016.root'

 datamix = args.mix+"/Data_13TeV_H4GTag_0"
 #datamix = '/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2016/hadd/data_mix.root'
 #datamix = '/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2016/hadd/data_mix_2016_highstat_skim.root'

 print "Data:   ",data
 print "Datamix:",datamix

 nbin = 10
 plots = []
 #plots.append(["cosThetaStarCS_dM","cosThetaStarCS_dM","Cos #theta*",nbin,-1.,1.])
 #plots.append(["cosTheta_a1_dM","cosTheta_a1_dM","Cos #theta_{#gamma a_{1}}",nbin,-1.,1.])
 #plots.append(["cosTheta_a2_dM","cosTheta_a2_dM","Cos #theta_{#gamma a_{2}}",nbin,-1.,1.])
 plots.append(["a1_a2_dR_dM","a1_a2_dR_dM","#Delta R (a1,a2)",nbin,0,7.])
 #plots.append(["a1_dR_dM","a1_dR_dM","#Delta R (#gamma1, #gamma2)",nbin,0,9.])
 #plots.append(["a2_dR_dM","a2_dR_dM","#Delta R (#gamma3, #gamma4)",nbin,0,9.])
 plots.append(["a1_pt_dM","a1_pt_dM","a1 (pT) [GeV]",nbin,0,200])
 plots.append(["a2_pt_dM","a2_pt_dM","a2 (pT) [GeV]",nbin,0,200])
 #plots.append(["a1_eta_dM","a1_eta_dM","a1 (#eta)",nbin,-6.,6.])
 #plots.append(["a2_eta_dM","a2_eta_dM","a2 (#eta)",nbin,-6.,6.]) 
 #plots.append(["a1_mass_dM","a1_mass_dM","a1 (Mass) [GeV]",nbin,0,150.])
 #plots.append(["a2_mass_dM","a2_mass_dM","a2 (Mass) [GeV]",nbin,0,150.])
 plots.append(["a1_mass_dM_M_a2_mass_dM","a1_mass_dM-a2_mass_dM","",nbin,-100, 100]) 
 #plots.append(["a1_energy_dM","a1_energy_dM","a1 (Energy) [GeV]",nbin,0,800.])
 #plots.append(["a2_energy_dM","a2_energy_dM","a2 (Energy) [GeV]",nbin,0,600.])
 #plots.append(["tp_pt","tp_pt","Higgs pT [GeV]",nbin,0,500.])
 #plots.append(["tp_eta","tp_eta","Higgs #eta",nbin,-5.,5.])
 #plots.append(["tp_energy","tp_energy","Higgs Energy [GeV]",nbin,-5.,5.])
 #plots.append(["pho1_MVA","pho1_MVA","#gamma1 MVA",nbin,-1.1,1])
 #plots.append(["pho2_MVA","pho2_MVA","#gamma2 MVA",nbin,-1.1,1])
 #plots.append(["pho3_MVA","pho3_MVA","#gamma3 MVA",nbin,-1.1,1])
 #plots.append(["pho4_MVA","pho4_MVA","#gamma4 MVA",nbin,-1.1,1])
 #plots.append(["pho1_pt","pho1_pt","#gamma1 pT [GeV]",nbin,30,400])
 #plots.append(["pho2_pt","pho2_pt","#gamma2 pT [GeV]",nbin,18,150])
 #plots.append(["pho3_pt","pho3_pt","#gamma3 pT [GeV]",nbin,15,100])
 #plots.append(["pho4_pt","pho4_pt","#gamma4 pT [GeV]",nbin,15,70])
 #plots.append(["pho1_energy","pho1_energy","#gamma1 Energy [GeV]",nbin,0.,2000.])
 #plots.append(["pho2_energy","pho2_energy","#gamma2 Energy [GeV]",nbin,0.,1500])
 #plots.append(["pho3_energy","pho3_energy","#gamma3 Energy [GeV]",nbin,0.,700])
 #plots.append(["pho4_energy","pho4_energy","#gamma4 Energy [GeV]",nbin,0.,500.])
 #plots.append(["pho1_eta","pho1_eta","#gamma1 #eta",nbin,-2.5,2.5])
 #plots.append(["pho2_eta","pho2_eta","#gamma2 #eta",nbin,-2.5,2.5])
 #plots.append(["pho3_eta","pho3_eta","#gamma3 #eta",nbin,-2.5,2.5])
 #plots.append(["pho4_eta","pho4_eta","#gamma4 #eta",nbin,-2.5,2.5])

 
 tree_datamix = ROOT.TChain() 
 tree_datamix.AddFile(datamix) 
 
 tree_data = ROOT.TChain() 
 tree_data.AddFile(data)   

 #Cut = '1>0' 
 Cut = ' pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && fabs(pho1_eta) < 2.5 && fabs(pho2_eta) < 2.5 && fabs(pho3_eta) < 2.5 && fabs(pho4_eta) < 2.5 && (fabs(pho1_eta) < 1.4442 || fabs(pho1_eta) > 1.566) && (fabs(pho2_eta) < 1.4442 || fabs(pho2_eta) > 1.566) && (fabs(pho3_eta) < 1.4442 || fabs(pho3_eta) > 1.566) && (fabs(pho4_eta) < 1.4442 || fabs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 100 && tp_mass < 180 && !(tp_mass > 115 && tp_mass < 135) '
 
 Histos = []
 weight_data = {}
 weight_datamix = {}
 weight_ratio = {}

 for plot in plots:
    Histos.append(ROOT.TH1F("histo_"+plot[0], plot[2], int(plot[3]), float(plot[4]), float(plot[5])))

 bins = []
 for bin in range(0,nbin+2):
    bins.append(bin) 
 
 perm = itertools.product(bins,repeat=len(Histos)) 
 for subset in list(perm):  
    weight_datamix[tuple(subset)] = 0.
    weight_data[tuple(subset)] = 0.
    weight_ratio[tuple(subset)] = 1.       
  
 #histo_scale_data = ROOT.TH1F("histo_scale_data","",100000,-1.1,1.)
 #histo_scale_datamix = ROOT.TH1F("histo_scale_datamix","",100000,-1.1,1.)
 #tree_datamix.Draw('pho1_MVA<-1.?-1.1:pho1_MVA >> histo_scale_datamix',ROOT.TCut(Cut))   
 #tree_data.Draw('pho1_MVA<-1.?-1.1:pho1_MVA >> histo_scale_data',ROOT.TCut(Cut))  
 #scale = float(histo_scale_data.Integral())/float(histo_scale_datamix.Integral())
  
 print "Computing occupancies for datamix..." 
 nentries = tree_datamix.GetEntries()
 for i in range(0, nentries):
    #if i%1000 == 0: print i
    tree_datamix.GetEntry(i)
    #if 1>0: 
    if tree_datamix.pho1_pt > 30 and tree_datamix.pho2_pt > 18 and tree_datamix.pho3_pt > 15 and tree_datamix.pho4_pt > 15 and  abs(tree_datamix.pho1_eta) < 2.5 and abs(tree_datamix.pho2_eta) < 2.5 and abs(tree_datamix.pho3_eta) < 2.5 and abs(tree_datamix.pho4_eta) < 2.5 and (abs(tree_datamix.pho1_eta) < 1.4442 or abs(tree_datamix.pho1_eta) > 1.566) and (abs(tree_datamix.pho2_eta) < 1.4442 or abs(tree_datamix.pho2_eta) > 1.566) and (abs(tree_datamix.pho3_eta) < 1.4442 or abs(tree_datamix.pho3_eta) > 1.566) and (abs(tree_datamix.pho4_eta) < 1.4442 or abs(tree_datamix.pho4_eta) > 1.566) and tree_datamix.pho1_electronveto==1 and tree_datamix.pho2_electronveto==1 and tree_datamix.pho3_electronveto==1 and tree_datamix.pho4_electronveto==1 and tree_datamix.tp_mass > 100 and tree_datamix.tp_mass < 180 and (tree_datamix.tp_mass <= 115 or tree_datamix.tp_mass >= 135): 
       subset = []
       for iVar,plot in enumerate(plots):       
          ibin = Histos[iVar].FindBin(findVal(tree_datamix,plot[0]))    
          subset.append(ibin)
          iVar+=1
       weight_datamix[tuple(subset)]+=1.
       
 print "Computing occupancies for data..." 
 nentries = tree_data.GetEntries()
 for i in range(0, nentries):
    #if i%1000 == 0: print i
    tree_data.GetEntry(i)
    #if 1>0: 
    if tree_data.pho1_pt > 30 and tree_data.pho2_pt > 18 and tree_data.pho3_pt > 15 and tree_data.pho4_pt > 15 and  abs(tree_data.pho1_eta) < 2.5 and abs(tree_data.pho2_eta) < 2.5 and abs(tree_data.pho3_eta) < 2.5 and abs(tree_data.pho4_eta) < 2.5 and (abs(tree_data.pho1_eta) < 1.4442 or abs(tree_data.pho1_eta) > 1.566) and (abs(tree_data.pho2_eta) < 1.4442 or abs(tree_data.pho2_eta) > 1.566) and (abs(tree_data.pho3_eta) < 1.4442 or abs(tree_data.pho3_eta) > 1.566) and (abs(tree_data.pho4_eta) < 1.4442 or abs(tree_data.pho4_eta) > 1.566) and tree_data.pho1_electronveto==1 and tree_data.pho2_electronveto==1 and tree_data.pho3_electronveto==1 and tree_data.pho4_electronveto==1 and tree_data.tp_mass > 100 and tree_data.tp_mass < 180 and (tree_data.tp_mass <= 115 or tree_data.tp_mass >= 135):  
       subset = []
       for iVar,plot in enumerate(plots):        
          ibin = Histos[iVar].FindBin(findVal(tree_data,plot[0]))    
          subset.append(ibin)
       weight_data[tuple(subset)]+=1.

 for subset in weight_ratio:
    if weight_datamix[subset]>0.: weight_ratio[subset]=float(weight_data[subset]/weight_datamix[subset])  

 print "Filling datamix weights..." 
 name = args.mix.replace('.root','_kinWeight.root')
 outfile = ROOT.TFile(name, "RECREATE")
 outtree = tree_datamix.CloneTree(0)
 mix_weight = array('f', [0])
 _mix_weight = outtree.Branch('weight', mix_weight, 'weight/F')     
 nentries = tree_datamix.GetEntries()
 for i in range(0, nentries):
    #if i%1000 == 0: print i
    tree_datamix.GetEntry(i)
    weight=1.
    if 1>0: 
    #if tree_datamix.pho1_pt > 30 and tree_datamix.pho2_pt > 18 and tree_datamix.pho3_pt > 15 and tree_datamix.pho4_pt > 15 and  abs(tree_datamix.pho1_eta) < 2.5 and abs(tree_datamix.pho2_eta) < 2.5 and abs(tree_datamix.pho3_eta) < 2.5 and abs(tree_datamix.pho4_eta) < 2.5 and (abs(tree_datamix.pho1_eta) < 1.4442 or abs(tree_datamix.pho1_eta) > 1.566) and (abs(tree_datamix.pho2_eta) < 1.4442 or abs(tree_datamix.pho2_eta) > 1.566) and (abs(tree_datamix.pho3_eta) < 1.4442 or abs(tree_datamix.pho3_eta) > 1.566) and (abs(tree_datamix.pho4_eta) < 1.4442 or abs(tree_datamix.pho4_eta) > 1.566) and tree_datamix.pho1_electronveto==1 and tree_datamix.pho2_electronveto==1 and tree_datamix.pho3_electronveto==1 and tree_datamix.pho4_electronveto==1 and tree_datamix.tp_mass > 100 and tree_datamix.tp_mass < 180 and (tree_datamix.tp_mass <= 115 or tree_datamix.tp_mass >= 135): 
       subset = []
       for iVar,plot in enumerate(plots):    
          ibin = Histos[iVar].FindBin(findVal(tree_datamix,plot[0]))         
          subset.append(ibin)
       weight=weight_ratio[tuple(subset)]
    mix_weight[0] = weight
    outtree.Fill() 
 outfile.cd()
 outtree.Write()
 outfile.Close()

      

