from ROOT import *

# Aprire il file
mc_file = TFile('/eos/cms/store/group/phys_tracking/None/crab_AOD_Run2022_deltaR_0p3_ptrel_2p0_deltaEta_0p3_tagdeltaz_1p0_triggerchain_fakerate_alltracks_allsubera_true_twomethods_or_isZmass/muon/Z/Run2022/AOD/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/total.root')
data_file = TFile('/afs/cern.ch/user/b/bdanzi/spark_tnp_updated/spark_tnp/pileup/data/Run2022.root')
# Create file with normalized pileup distribution
theFile = TFile('mc_pileup_2022.root',"RECREATE")
# Prendere il tree
mc_tree = mc_file.Get("muon/StandAloneEvents")
histo_data_pileup = data_file.Get("pileup")
# Definire l'istogramma
histo_pileup_mc = TH1D('pileup','pileup', histo_data_pileup.GetNbinsX(), (histo_data_pileup.GetXaxis()).GetBinLowEdge(1), (histo_data_pileup.GetXaxis()).GetBinUpEdge(histo_data_pileup.GetXaxis().GetNbins()))
#3. looppare sul tree
for entry in mc_tree :
    histo_pileup_mc.Fill(entry.nTrueInteractions)
    if entry.nTrueInteractions > -999.:  histo_pileup_mc.Fill(entry.nTrueInteractions)
gStyle.SetOptTitle(0)
#######################################
#4. fuori dal loop for ################
histo_pileup_mc.Scale(1./histo_pileup_mc.Integral(0,-1))
#######################################
theFile.Write()
theFile.Close()
#######################################
mc_file.Close()
data_file.Close()


