import ROOT
import argparse
import os

ROOT.gROOT.SetBatch(ROOT.kTRUE)

parser =  argparse.ArgumentParser(description='cat MVA')
parser.add_argument('-d', '--outdir', dest='outdir', required=True, type=str)

args = parser.parse_args()
outdir = args.outdir
os.system('mkdir '+outdir)

Cut_Signal = '(pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && fabs(pho1_eta) < 2.5 && fabs(pho2_eta) < 2.5 && fabs(pho3_eta) < 2.5 && fabs(pho4_eta) < 2.5 && (fabs(pho1_eta) < 1.4442 || fabs(pho1_eta) > 1.566) && (fabs(pho2_eta) < 1.4442 || fabs(pho2_eta) > 1.566) && (fabs(pho3_eta) < 1.4442 || fabs(pho3_eta) > 1.566) && (fabs(pho4_eta) < 1.4442 || fabs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 100 && tp_mass < 180)'

Cut_Background = '(pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && fabs(pho1_eta) < 2.5 && fabs(pho2_eta) < 2.5 && fabs(pho3_eta) < 2.5 && fabs(pho4_eta) < 2.5 && (fabs(pho1_eta) < 1.4442 || fabs(pho1_eta) > 1.566) && (fabs(pho2_eta) < 1.4442 || fabs(pho2_eta) > 1.566) && (fabs(pho3_eta) < 1.4442 || fabs(pho3_eta) > 1.566) && (fabs(pho4_eta) < 1.4442 || fabs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 100 && tp_mass < 180)'

Cut_Data = '(pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && fabs(pho1_eta) < 2.5 && fabs(pho2_eta) < 2.5 && fabs(pho3_eta) < 2.5 && fabs(pho4_eta) < 2.5 && (fabs(pho1_eta) < 1.4442 || fabs(pho1_eta) > 1.566) && (fabs(pho2_eta) < 1.4442 || fabs(pho2_eta) > 1.566) && (fabs(pho3_eta) < 1.4442 || fabs(pho3_eta) > 1.566) && (fabs(pho4_eta) < 1.4442 || fabs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 100 && tp_mass < 180 && !(tp_mass > 115 && tp_mass < 135))'

#Tail cut: && a1_pt_dM<140. && a2_pt_dM<80. && fabs(a1_mass_dM-a2_mass_dM)<55. && a1_a2_dR_dM<6.2 && a1_pt_dM/tp_mass<0.8 && a2_pt_dM/tp_mass<0.8

histo_scale_data = ROOT.TH1F("histo_scale_data","",100000,-1.1,1.)
histo_scale_datamix = ROOT.TH1F("histo_scale_datamix","",100000,-1.1,1.)
histo_scale_signal = ROOT.TH1F("histo_scale_signal","",100000,-1.1,1.)

sig_list = [15,20,25,30,40,45,50,55,60]
sig_nums_2016 = []
sig_nums_2017 = []
sig_nums_2018 = []
sig_file_2016 = []
sig_file_2017 = []
sig_file_2018 = []

for mass in sig_list:
   sig_file_2016_tmp = ROOT.TChain()
   tree_name_2016 = '/eos/user/t/twamorka/H4G_Signal_Campaign_18Jan2021/2016/hadd/signal_m_'+str(mass)+'_odd_genMass.root/HAHMHToAA_AToGG_MA_'+str(mass)+'GeV_TuneCUETP8M1_PSweights_13TeV_madgraph_pythia8_13TeV_H4GTag_0'
   sig_file_2016_tmp.AddFile(tree_name_2016)
   histo_scale_signal.Reset() 
   sig_file_2016_tmp.Draw("pho1_MVA<-1.?-1.1:pho1_MVA>>histo_scale_signal","weight*"+Cut_Signal)
   sig_file_2016.append(sig_file_2016_tmp)      
   sig_nums_2016.append(float(histo_scale_signal.Integral())) 
   
   sig_file_2017_tmp = ROOT.TChain()
   tree_name_2017 = '/eos/user/t/twamorka/H4G_Signal_Campaign_18Jan2021/2017/hadd/signal_m_'+str(mass)+'_odd_genMass.root/HAHMHToAA_AToGG_MA_'+str(mass)+'GeV_TuneCP5_PSweights_13TeV_madgraph_pythia8_13TeV_H4GTag_0'
   sig_file_2017_tmp.AddFile(tree_name_2017)
   histo_scale_signal.Reset() 
   sig_file_2017_tmp.Draw("pho1_MVA<-1.?-1.1:pho1_MVA>>histo_scale_signal","weight*"+Cut_Signal) 
   sig_file_2017.append(sig_file_2017_tmp)      
   sig_nums_2017.append(float(histo_scale_signal.Integral())) 

   sig_file_2018_tmp = ROOT.TChain()
   tree_name_2018 = '/eos/user/t/twamorka/H4G_Signal_Campaign_18Jan2021/2018/hadd/signal_m_'+str(mass)+'_odd_genMass.root/HAHMHToAA_AToGG_MA_'+str(mass)+'GeV_TuneCP5_PSweights_13TeV_madgraph_pythia8_13TeV_H4GTag_0'
   sig_file_2018_tmp.AddFile(tree_name_2018)
   histo_scale_signal.Reset() 
   sig_file_2018_tmp.Draw("pho1_MVA<-1.?-1.1:pho1_MVA>>histo_scale_signal","weight*"+Cut_Signal)
   sig_file_2018.append(sig_file_2018_tmp)      
   sig_nums_2018.append(float(histo_scale_signal.Integral())) 

lumi_2016 = 35.9
data_file_2016 = ROOT.TChain()
#data_file_2016.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2016/hadd/data_2016.root/tagsDumper/trees/Data_13TeV_H4GTag_0')
data_file_2016.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2016/hadd/data_2016_genMass.root/Data_13TeV_H4GTag_0')
histo_scale_data.Reset() 
data_file_2016.Draw("pho1_MVA<-1.?-1.1:pho1_MVA>>histo_scale_data",Cut_Data)
data_2016 = float(histo_scale_data.Integral())
datamix_file_2016 = ROOT.TChain()
datamix_file_2016.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/DataMix_HighStat/hadd/data_mix_2016_kinWeight_odd_genMass.root/Data_13TeV_H4GTag_0')
histo_scale_datamix.Reset() 
datamix_file_2016.Draw("pho1_MVA<-1.?-1.1:pho1_MVA>>histo_scale_datamix","weight*"+Cut_Data)
datamix_2016_SB = float(histo_scale_datamix.Integral())
histo_scale_datamix.Reset() 
datamix_file_2016.Draw("pho1_MVA<-1.?-1.1:pho1_MVA>>histo_scale_datamix","weight*"+Cut_Background)
datamix_2016 = float(histo_scale_datamix.Integral())

lumi_2017 = 41.5
data_file_2017 = ROOT.TChain()
#data_file_2017.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2017/hadd/data_2017.root/tagsDumper/trees/Data_13TeV_H4GTag_0')
data_file_2017.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2017/hadd/data_2017_genMass.root/Data_13TeV_H4GTag_0')
histo_scale_data.Reset() 
data_file_2017.Draw("pho1_MVA<-1.?-1.1:pho1_MVA>>histo_scale_data",Cut_Data)
data_2017 = float(histo_scale_data.Integral())
datamix_file_2017 = ROOT.TChain()
datamix_file_2017.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/DataMix_HighStat/hadd/data_mix_2017_kinWeight_odd_genMass.root/Data_13TeV_H4GTag_0')
histo_scale_datamix.Reset() 
datamix_file_2017.Draw("pho1_MVA<-1.?-1.1:pho1_MVA>>histo_scale_datamix","weight*"+Cut_Data)
datamix_2017_SB = float(histo_scale_datamix.Integral())
histo_scale_datamix.Reset() 
datamix_file_2017.Draw("pho1_MVA<-1.?-1.1:pho1_MVA>>histo_scale_datamix","weight*"+Cut_Background)
datamix_2017 = float(histo_scale_datamix.Integral())

lumi_2018 = 54.38
data_file_2018 = ROOT.TChain()
#data_file_2018.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2018/hadd/data_2018.root/tagsDumper/trees/Data_13TeV_H4GTag_0')
data_file_2018.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2018/hadd/data_2018_genMass.root/Data_13TeV_H4GTag_0')
histo_scale_data.Reset() 
data_file_2018.Draw("pho1_MVA<-1.?-1.1:pho1_MVA>>histo_scale_data",Cut_Data)
data_2018 = float(histo_scale_data.Integral())
datamix_file_2018 = ROOT.TChain()
datamix_file_2018.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/DataMix_HighStat/hadd/data_mix_2018_kinWeight_odd_genMass.root/Data_13TeV_H4GTag_0')
histo_scale_datamix.Reset() 
datamix_file_2018.Draw("pho1_MVA<-1.?-1.1:pho1_MVA>>histo_scale_datamix","weight*"+Cut_Data)
datamix_2018_SB = float(histo_scale_datamix.Integral())
histo_scale_datamix.Reset() 
datamix_file_2018.Draw("pho1_MVA<-1.?-1.1:pho1_MVA>>histo_scale_datamix","weight*"+Cut_Background)
datamix_2018 = float(histo_scale_datamix.Integral())

print "datamix_2016_SB:",datamix_2016_SB
print "datamix_2017_SB:",datamix_2017_SB
print "datamix_2018_SB:",datamix_2018_SB

f_out = ROOT.TFile(outdir+'/categorization_training.root','RECREATE')
ROOT.TMVA.Tools.Instance()
factory = ROOT.TMVA.Factory("TMVAClassification", f_out,"AnalysisType=Classification")

dataloader = ROOT.TMVA.DataLoader(outdir)

#mvaVars = ['pho1_MVA<-1.? -1.1: pho1_MVA', 'pho2_MVA<-1.? -1.1: pho2_MVA', 'pho3_MVA<-1.? -1.1: pho3_MVA', 'pho4_MVA<-1.? -1.1: pho4_MVA', 'a1_mass_dM-a2_mass_dM', 'a1_mass_dM-genMass', 'a2_mass_dM-genMass', 'cosTheta_a1_dM', 'a1_pt_dM', 'a2_pt_dM', 'a1_a2_dR_dM', 'a1_pt_dM/tp_mass', 'a2_pt_dM/tp_mass']

mvaVars = ['pho1_MVA<-1.? -1.1: pho1_MVA', 'pho2_MVA<-1.? -1.1: pho2_MVA', 'pho3_MVA<-1.? -1.1: pho3_MVA', 'pho4_MVA<-1.? -1.1: pho4_MVA', 'a1_mass_dM-a2_mass_dM', '(a1_mass_dM-genMass)/tp_mass', '(a2_mass_dM-genMass)/tp_mass', 'cosTheta_a1_dM', 'a1_pt_dM', 'a2_pt_dM', 'a1_a2_dR_dM/tp_mass']

for x in mvaVars:
   #factory.AddVariable(x,"F")
   dataloader.AddVariable(x,"F")

for i,mass in enumerate(sig_list):
   print "---> Signal normalization 2016: mass",mass,"-",float(sig_nums_2016[i])
   #dataloader.AddSignalTree(sig_file_2016[i],lumi_2016/float(sig_nums_2016[i]))
   dataloader.AddSignalTree(sig_file_2016[i],lumi_2016) 
   print "---> Signal normalization 2017: mass",mass,"-",float(sig_nums_2017[i])
   #dataloader.AddSignalTree(sig_file_2017[i],lumi_2017/float(sig_nums_2017[i])) 
   dataloader.AddSignalTree(sig_file_2017[i],lumi_2017)   
   print "---> Signal normalization 2018: mass",mass,"-",float(sig_nums_2018[i])
   #dataloader.AddSignalTree(sig_file_2018[i],lumi_2018/float(sig_nums_2018[i]))
   dataloader.AddSignalTree(sig_file_2018[i],lumi_2018) 

scale_2016 = data_2016/datamix_2016_SB
#scale_2016 = lumi_2016/datamix_2016
scale_2017 = data_2017/datamix_2017_SB
#scale_2017 = lumi_2017/datamix_2017
scale_2018 = data_2018/datamix_2018_SB
#scale_2018 = lumi_2018/datamix_2018

print "---> Background normalization 2016: ",scale_2016
dataloader.AddBackgroundTree(datamix_file_2016,scale_2016)
print "---> Background normalization 2017: ",scale_2017
dataloader.AddBackgroundTree(datamix_file_2017,scale_2017)
print "---> Background normalization 2018: ",scale_2018
dataloader.AddBackgroundTree(datamix_file_2018,scale_2018)

sigCut = ROOT.TCut(Cut_Signal)
bkgCut = ROOT.TCut(Cut_Background)

print "S Cut: ", sigCut
print "B Cut: ", bkgCut

dataloader.SetWeightExpression("weight","Signal")
dataloader.SetWeightExpression("weight","Background")

dataloader.PrepareTrainingAndTestTree(sigCut,bkgCut,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V")
method = factory.BookMethod( dataloader,ROOT.TMVA.Types.kBDT, "BDTG", "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" )

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

