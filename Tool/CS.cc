#include <iostream>
#include <algorithm>
#include <cstring>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/baselineDef.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "SusyAnaTools/Tools/ISRCorrector.h"
#include "TStopwatch.h"
#include "TString.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"
#include "CS.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraphAsymmErrors.h"

#include "PrepSystematics.h"
#include "CommonShared.h"

#include "SusyAnaTools/Tools/PileupWeights.h"

using namespace std;

const bool doISR = true;
const bool dobSF = true;
const bool dolepSF = true;
const bool doPU = false;

const bool dospecialFix = false;
const bool useTTbarHTsamples = true;

TFile * bTagEffFile =0;

TH2D * mu_mediumID_SF = 0, * mu_miniISO_SF = 0;
TGraphAsymmErrors * mu_trkptGT10_SF = 0, * mu_trkptLT10_SF = 0;

TH2D * ele_VetoID_SF = 0, * ele_miniISO_SF = 0;
TH2D * ele_trkpt_SF = 0;

std::shared_ptr<TopTagger> ttPtr(nullptr);

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
  {
    std::cerr <<"Please give 5 arguments "<<"SubsampleName"<<" Input Template" <<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
    std::cerr <<" Valid configurations are " << std::endl;
    std::cerr <<" ./CS TTbarInc 1000 0 1" << std::endl;
    return -1;
  }
  const char *subsamplename = argv[1];
  const char *Maxevent = argv[2];
  const  char *Stratfile = argv[3];
  const  char *Filerun = argv[4];
  const  int startfile = std::atoi(Stratfile);
  const int filerun = std::atoi(Filerun);
  bool isData = false;
  TChain *fChain = 0;
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(subsamplename, startfile);
  const string condorSpec = argc==6 ? argv[5]: "";

  if(!FillChain(fChain, subsamplename, condorSpec, startfile, filerun))
  {
    std::cerr << "Cannot get the tree " << std::endl;
  }
  const int maxevent = std::atoi(Maxevent);  
  
  TString sampleString(subsamplename);
  if(sampleString.Contains("Data")){Lumiscale = 1.0; isData = true;}
  //Searchbin                                                                                                                                                                                    
  // SearchBins SB("SB_v1_2017");
  SearchBins SB("SB_Aggregate_2017");

//  const int nTotBins = SB.nSearchBins();
  std::vector<int> cached_MT2_binIdx_mu_3DVec, cached_HT_binIdx_mu_3DVec;
  std::vector<TH1D*> cached_MT2_hist_mu_3DVec, cached_HT_hist_mu_3DVec;
  std::vector<int> cached_MT2_binIdx_el_3DVec, cached_HT_binIdx_el_3DVec;
  std::vector<TH1D*> cached_MT2_hist_el_3DVec, cached_HT_hist_el_3DVec;

  std::vector<int> cached_MT2_binIdx_mu_2DVec, cached_HT_binIdx_mu_2DVec;
  std::vector<TH1D*> cached_MT2_hist_mu_2DVec, cached_HT_hist_mu_2DVec, cached_met_hist_mu_2DVec, cached_nJets_hist_mu_2DVec, cached_recoTopPt_hist_mu_2DVec;
  std::vector<int> cached_MT2_binIdx_el_2DVec, cached_HT_binIdx_el_2DVec;
  std::vector<TH1D*> cached_MT2_hist_el_2DVec, cached_HT_hist_el_2DVec, cached_met_hist_el_2DVec, cached_nJets_hist_el_2DVec, cached_recoTopPt_hist_el_2DVec;

  TH1D * h1_genHT = new TH1D("genHT", "genHT", 60, 0, 3000);

  std::string spec = "CS";

  AnaFunctions::prepareForNtupleReader();
  //AnaFunctions::prepareTopTagger();
  NTupleReader *tr =0;
  if( isData ) tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames_DataOnly);
  else tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames);
  CSBaseline = new BaselineVessel(*tr, spec);
  //CSBaseline->SetupTopTagger(true,"TopTagger.cfg");
  CSBaseline->SetupTopTagger(true,"TopTagger_Simplified.cfg");
  tr->registerFunction((*CSBaseline));

  ttPtr = CSBaseline->GetTopTaggerPtr();

  if( !isData ){
    // Lepton SF
    TFile * allINone_leptonSF_file = new TFile("allINone_leptonSF_Moriond17.root");
    if( !allINone_leptonSF_file->IsZombie() ){
      mu_mediumID_SF = (TH2D*) allINone_leptonSF_file->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0");
      mu_miniISO_SF = (TH2D*) allINone_leptonSF_file->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass");
// SF for muon tracks
      mu_trkptGT10_SF = (TGraphAsymmErrors*) allINone_leptonSF_file->Get("ratio_eff_eta3_dr030e030_corr");
      mu_trkptLT10_SF = (TGraphAsymmErrors*) allINone_leptonSF_file->Get("ratio_eff_eta3_tk0_dr030e030_corr");
  
      ele_VetoID_SF = (TH2D*) allINone_leptonSF_file->Get("GsfElectronToVeto");
      ele_miniISO_SF = (TH2D*) allINone_leptonSF_file->Get("MVAVLooseElectronToMini");
      ele_trkpt_SF = (TH2D*) allINone_leptonSF_file->Get("EGamma_SF2D");
    }

    //BTag SF
    BTagCorrector btagcorr = BTagCorrector("allINone_bTagEff.root");
    if( sampleString == "TTbarSingleLepT" || sampleString == "TTbarSingleLepT_deepTrimmed" ) btagcorr.resetEffs("TTbarSingleLepT");
    if( sampleString == "TTbarSingleLepTbar" || sampleString == "TTbarSingleLepTbar_deepTrimmed" ) btagcorr.resetEffs("TTbarSingleLepTbar");
    if( sampleString == "TTbarDiLep" || sampleString == "TTbarDiLep_deepTrimmed" ) btagcorr.resetEffs("TTbarDiLep");
    if( sampleString == "TTbar_HT-600to800" || sampleString == "TTbar_HT-600to800_deepTrimmed" ) btagcorr.resetEffs("TTbar_HT-600to800");
    if( sampleString == "TTbar_HT-800to1200" || sampleString == "TTbar_HT-800to1200_deepTrimmed" ) btagcorr.resetEffs("TTbar_HT-800to1200");
    if( sampleString == "TTbar_HT-1200to2500" || sampleString == "TTbar_HT-1200to2500_deepTrimmed" ) btagcorr.resetEffs("TTbar_HT-1200to2500");
    if( sampleString == "TTbar_HT-2500toInf" || sampleString == "TTbar_HT-2500toInf_deepTrimmed" ) btagcorr.resetEffs("TTbar_HT-2500toInf");
    if( sampleString == "WJetsToLNu_HT_70to100" || sampleString == "WJetsToLNu_HT_70to100_deepTrimmed" ) btagcorr.resetEffs("WJetsToLNu_HT_70to100");
    if( sampleString == "WJetsToLNu_HT_100to200" || sampleString == "WJetsToLNu_HT_100to200_deepTrimmed" ) btagcorr.resetEffs("WJetsToLNu_HT_100to200");
    if( sampleString == "WJetsToLNu_HT_200to400" || sampleString == "WJetsToLNu_HT_200to400_deepTrimmed" ) btagcorr.resetEffs("WJetsToLNu_HT_200to400");
    if( sampleString == "WJetsToLNu_HT_400to600" || sampleString == "WJetsToLNu_HT_400to600_deepTrimmed" ) btagcorr.resetEffs("WJetsToLNu_HT_400to600");
    if( sampleString == "WJetsToLNu_HT_600to800" || sampleString == "WJetsToLNu_HT_600to800_deepTrimmed" ) btagcorr.resetEffs("WJetsToLNu_HT_600to800");
    if( sampleString == "WJetsToLNu_HT_800to1200" || sampleString == "WJetsToLNu_HT_800to1200_deepTrimmed" ) btagcorr.resetEffs("WJetsToLNu_HT_800to1200");
    if( sampleString == "WJetsToLNu_HT_1200to2500" || sampleString == "WJetsToLNu_HT_1200to2500_deepTrimmed" ) btagcorr.resetEffs("WJetsToLNu_HT_1200to2500");
    if( sampleString == "WJetsToLNu_HT_2500toInf" || sampleString == "WJetsToLNu_HT_2500toInf_deepTrimmed" ) btagcorr.resetEffs("WJetsToLNu_HT_2500toInf");
    if( sampleString == "tW_top_incl" || sampleString == "tW_top_incl_deepTrimmed" ) btagcorr.resetEffs("tW_top_incl");
    if( sampleString == "tW_antitop_incl" || sampleString == "tW_antitop_incl_deepTrimmed" ) btagcorr.resetEffs("tW_antitop_incl");
    tr->registerFunction(btagcorr);
  
    //ISR Correcotr
    ISRCorrector *isrcorr = new ISRCorrector("allINone_ISRJets.root", "", "");
    if( sampleString == "TTbarSingleLepT" || sampleString == "TTbarSingleLepT_deepTrimmed" ) isrcorr->resetSample("TTbarSingleLepT");
    if( sampleString == "TTbarSingleLepTbar" || sampleString == "TTbarSingleLepTbar_deepTrimmed" ) isrcorr->resetSample("TTbarSingleLepTbar");
    if( sampleString == "TTbarDiLep" || sampleString == "TTbarDiLep_deepTrimmed" ) isrcorr->resetSample("TTbarDiLep");
    if( sampleString == "TTbar_HT-600to800" || sampleString == "TTbar_HT-600to800_deepTrimmed" ) isrcorr->resetSample("TTbar_HT-600to800");
    if( sampleString == "TTbar_HT-800to1200" || sampleString == "TTbar_HT-800to1200_deepTrimmed" ) isrcorr->resetSample("TTbar_HT-800to1200");
    if( sampleString == "TTbar_HT-1200to2500" || sampleString == "TTbar_HT-1200to2500_deepTrimmed" ) isrcorr->resetSample("TTbar_HT-1200to2500");
    if( sampleString == "TTbar_HT-2500toInf" || sampleString == "TTbar_HT-2500toInf_deepTrimmed" ) isrcorr->resetSample("TTbar_HT-2500toInf");
    if( sampleString == "WJetsToLNu_HT_70to100" || sampleString == "WJetsToLNu_HT_70to100_deepTrimmed" ) isrcorr->resetSample("WJetsToLNu_HT_70to100");
    if( sampleString == "WJetsToLNu_HT_100to200" || sampleString == "WJetsToLNu_HT_100to200_deepTrimmed" ) isrcorr->resetSample("WJetsToLNu_HT_100to200");
    if( sampleString == "WJetsToLNu_HT_200to400" || sampleString == "WJetsToLNu_HT_200to400_deepTrimmed" ) isrcorr->resetSample("WJetsToLNu_HT_200to400");
    if( sampleString == "WJetsToLNu_HT_400to600" || sampleString == "WJetsToLNu_HT_400to600_deepTrimmed" ) isrcorr->resetSample("WJetsToLNu_HT_400to600");
    if( sampleString == "WJetsToLNu_HT_600to800" || sampleString == "WJetsToLNu_HT_600to800_deepTrimmed" ) isrcorr->resetSample("WJetsToLNu_HT_600to800");
    if( sampleString == "WJetsToLNu_HT_800to1200" || sampleString == "WJetsToLNu_HT_800to1200_deepTrimmed" ) isrcorr->resetSample("WJetsToLNu_HT_800to1200");
    if( sampleString == "WJetsToLNu_HT_1200to2500" || sampleString == "WJetsToLNu_HT_1200to2500_deepTrimmed" ) isrcorr->resetSample("WJetsToLNu_HT_1200to2500");
    if( sampleString == "WJetsToLNu_HT_2500toInf" || sampleString == "WJetsToLNu_HT_2500toInf_deepTrimmed" ) isrcorr->resetSample("WJetsToLNu_HT_2500toInf");
    if( sampleString == "tW_top_incl" || sampleString == "tW_top_incl_deepTrimmed" ) isrcorr->resetSample("tW_top_incl");
    if( sampleString == "tW_antitop_incl" || sampleString == "tW_antitop_incl_deepTrimmed" ) isrcorr->resetSample("tW_antitop_incl");
    tr->registerFunction((*isrcorr));

    systBaseline::pdfScale = new PDFUncertainty();

    systBaseline::SRblv_jecUp = new BaselineVessel(*tr, systBaseline::spec_jecUp);
    systBaseline::SRblv_jecDn = new BaselineVessel(*tr, systBaseline::spec_jecDn);
    systBaseline::SRblv_metMagUp = new BaselineVessel(*tr, systBaseline::spec_metMagUp);
    systBaseline::SRblv_metMagDn = new BaselineVessel(*tr, systBaseline::spec_metMagDn);
    systBaseline::SRblv_metPhiUp = new BaselineVessel(*tr, systBaseline::spec_metPhiUp);
    systBaseline::SRblv_metPhiDn = new BaselineVessel(*tr, systBaseline::spec_metPhiDn);

    tr->registerFunction((*systBaseline::pdfScale));
 
    PrepSystematics sysPrep;

    tr->registerFunction(sysPrep);

    tr->registerFunction((*systBaseline::SRblv_jecUp));
    tr->registerFunction((*systBaseline::SRblv_jecDn));
    tr->registerFunction((*systBaseline::SRblv_metMagUp));
    tr->registerFunction((*systBaseline::SRblv_metMagDn));
    tr->registerFunction((*systBaseline::SRblv_metPhiUp));
    tr->registerFunction((*systBaseline::SRblv_metPhiDn));

    Pileup_Sys pileup("PileupHistograms_0121_69p2mb_pm4p6.root");
    tr->registerFunction(pileup);
  }

  //Add cleanJets function
  //stopFunctions::cleanJets.setMuonIso("mini");
  //stopFunctions::cleanJets.setElecIso("mini");
  //stopFunctions::cleanJets.setRemove(false);
  //tr->registerFunction(&stopFunctions::cleanJets);
 
  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr->getNEntries();
  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<Lumiscale<<std::endl; 
  cout<<"maxevent: "<<maxevent<<endl;
  int k =0;
  // Loop over the events (tree entries)
  while(tr->getNextEvent())
  {
    k++;
    //std::cout<<"evt: "<<k<<std::endl;
    if(maxevent>=0 && tr->getEvtNum() > maxevent ) break;
    // Add print out of the progress of looping
    if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;

    const double isrWght = (doISR && !isData)? tr->getVar<double>("isr_Unc_Cent") : 1.0;
    const double bSF = (dobSF && !isData)? tr->getVar<double>("bTagSF_EventWeightSimple_Central") : 1.0;
    const double bSF_up = (dobSF && !isData)? tr->getVar<double>("bTagSF_EventWeightSimple_Up") : 1.0;
    const double bSF_down = (dobSF && !isData)? tr->getVar<double>("bTagSF_EventWeightSimple_Down") : 1.0;

    const double isr_up = (doISR && !isData)? tr->getVar<double>("isr_Unc_Up"): 1.0;
    const double isr_down = (doISR && !isData)? tr->getVar<double>("isr_Unc_Down"):1.0;

    const double NNPDF_From_Median_Central = isData? 1.0 : tr->getVar<double>("NNPDF_From_Median_Central");
    const double NNPDF_From_Median_Up = isData? 1.0 : tr->getVar<double>("NNPDF_From_Median_Up");
    const double NNPDF_From_Median_Down = isData? 1.0 : tr->getVar<double>("NNPDF_From_Median_Down");

    const double Scaled_Variations_Up = isData? 1.0 : tr->getVar<double>("Scaled_Variations_Up");
    const double Scaled_Variations_Down = isData? 1.0 : tr->getVar<double>("Scaled_Variations_Down");

    const double puWght = (doPU && !isData)? tr->getVar<double>("_PUweightFactor") : 1.0;
    const double puSys_up = (doPU && !isData)? tr->getVar<double>("_PUSysUp") : 1.0;
    const double puSys_down = (doPU && !isData)? tr->getVar<double>("_PUSysDown") : 1.0;

    const double genHT = tr->hasVar("genHT") ? tr->getVar<double>("genHT") : -999;

    const vector<TLorentzVector> &muonsLVec = tr->getVec<TLorentzVector>("muonsLVec");
    const vector<double> &muonsRelIso = tr->getVec<double>("muonsRelIso");
    const vector<double> &muonsMiniIso = tr->getVec<double>("muonsMiniIso");
    const vector<double> &muonsMtw = tr->getVec<double>("muonsMtw");    
    const vector<TLorentzVector> &elesLVec = tr->getVec<TLorentzVector>("elesLVec");
    const vector<double> &elesRelIso = tr->getVec<double>("elesRelIso");
    const vector<double> &elesMiniIso = tr->getVec<double>("elesMiniIso");
    const vector<double> &elesMtw = tr->getVec<double>("elesMtw");
    const vector<unsigned int> &elesisEB = tr->getVec<unsigned int>("elesisEB");    
    const vector<TLorentzVector> &jetsLVec = tr->getVec<TLorentzVector>("jetsLVec");
    const vector<double> &recoJetsBtag_0 = tr->getVec<double>("recoJetsBtag_0");
    const vector<int> &looseisoTrksMatchedJetIdx = tr->getVec<int>("looseisoTrksMatchedJetIdx");
    const vector<TLorentzVector> &loose_isoTrksLVec = tr->getVec<TLorentzVector>("loose_isoTrksLVec");
    const vector<double> &loose_isoTrks_iso = tr->getVec<double>("loose_isoTrks_iso");
    const vector<double> &loose_isoTrks_mtw = tr->getVec<double>("loose_isoTrks_mtw");
    const vector<int> &loose_isoTrks_pdgId = tr->getVec<int>("loose_isoTrks_pdgId");
    const std::vector<double> muonspfActivity = tr->getVec<double>("muonspfActivity");
    const std::vector<int> & muonsFlagIDVec = tr->getVec<int>("muonsFlagMedium");
    const std::vector<int> & elesFlagIDVec = tr->getVec<int>("elesFlagVeto");
    const std::vector<double>& recoJetschargedHadronEnergyFraction = tr->getVec<double>("recoJetschargedHadronEnergyFraction");
    const std::vector<double>& recoJetschargedEmEnergyFraction = tr->getVec<double>("recoJetschargedEmEnergyFraction");
    const std::vector<std::string> & TriggerNames = tr->getVec<std::string>("TriggerNames");
    const std::vector<int> & PassTrigger = tr->getVec<int>("PassTrigger");

    int run = tr->getVar<int>("run");
    int lumi = tr->getVar<int>("lumi");
    int event = tr->getVar<int>("event");
    int vtxSize = tr->getVar<int>("vtxSize");

    bool passNoiseEventFilter = tr->getVar<bool>("passNoiseEventFilter"+spec);
    const double EvtWt = tr->getVar<double>("evtWeight");
    //change event weight for MC sample
    EventWeight = EvtWt;
    Lumiscale = Lumiscale * EventWeight;
    //Event Filter                                                                                                                           
    //if(!passNoiseEventFilter) continue;
    bool passBaseline = tr->getVar<bool>("passBaseline"+spec);
    bool passMuonVeto = tr->getVar<bool>("passMuonVeto"+spec);
    bool passEleVeto = tr->getVar<bool>("passEleVeto"+spec);
    bool passIsoTrkVeto = tr->getVar<bool>("passIsoTrkVeto"+spec);
    bool passLeptVeto = tr->getVar<bool>("passLeptVeto"+spec);
    bool passnJets = tr->getVar<bool>("passnJets"+spec);
    bool passdPhis = tr->getVar<bool>("passdPhis"+spec);
    bool passMET =  tr->getVar<bool>("passMET"+spec);
    bool passBJets = tr->getVar<bool>("passBJets"+spec);
    bool passTagger = tr->getVar<bool>("passTagger"+spec);
    bool passHT =  tr->getVar<bool>("passHT"+spec);
    bool passMT2 = tr->getVar<bool>("passMT2" + spec);
    const int nElectrons = tr->getVar<int>("nElectrons_CUT"+spec);
    const int nMuons = tr->getVar<int>("nMuons_CUT"+spec);
    const int nJets = tr->getVar<int>("cntNJetsPt30Eta24"+spec);
    const int nbJets = tr->getVar<int>("cntCSVS"+spec);
    const int nTops = tr->getVar<int>("nTopCandSortedCnt"+spec);
    const double MT2 = tr->getVar<double>("best_had_brJet_MT2"+spec);
    const vector<double> dPhiVec = tr->getVec<double>("dPhiVec"+spec);
    const double HT = tr->getVar<double>("HT"+spec);
    double met=tr->getVar<double>("met");
    double metphi=tr->getVar<double>("metphi");
    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);

    const double met_metMagUp = isData? met : tr->getVar<double>("met_metMagUp");
    const double metphi_metMagUp = isData? metphi : tr->getVar<double>("metphi");
    const double MT2_metMagUp = isData? MT2 : tr->getVar<double>("best_had_brJet_MT2" + systBaseline::spec_metMagUp);
    const int nbJets_metMagUp = isData? nbJets : tr->getVar<int>("cntCSVS" + systBaseline::spec_metMagUp);
    const int nTops_metMagUp = isData? nTops : tr->getVar<int>("nTopCandSortedCnt" + systBaseline::spec_metMagUp);
    const double HT_metMagUp = isData? HT : tr->getVar<double>("HT" + systBaseline::spec_metMagUp);
    const bool passBaselineNoLepVeto_metMagUp = isData? tr->getVar<bool>("passBaselineNoLepVeto"+spec) : tr->getVar<bool>("passBaselineNoLepVeto" + systBaseline::spec_metMagUp);

    const double met_metMagDn = isData? met : tr->getVar<double>("met_metMagDn");
    const double metphi_metMagDn = isData? metphi : tr->getVar<double>("metphi");
    const double MT2_metMagDn = isData? MT2 : tr->getVar<double>("best_had_brJet_MT2" + systBaseline::spec_metMagDn);
    const int nbJets_metMagDn = isData? nbJets : tr->getVar<int>("cntCSVS" + systBaseline::spec_metMagDn);
    const int nTops_metMagDn = isData? nTops : tr->getVar<int>("nTopCandSortedCnt" + systBaseline::spec_metMagDn);
    const double HT_metMagDn = isData? HT : tr->getVar<double>("HT" + systBaseline::spec_metMagDn);
    const bool passBaselineNoLepVeto_metMagDn = isData? tr->getVar<bool>("passBaselineNoLepVeto"+spec) : tr->getVar<bool>("passBaselineNoLepVeto" + systBaseline::spec_metMagDn);

    const double met_metPhiUp = isData? met : tr->getVar<double>("met");
    const double metphi_metPhiUp = isData? metphi : tr->getVar<double>("metphi_metPhiUp");
    const double MT2_metPhiUp = isData? MT2 : tr->getVar<double>("best_had_brJet_MT2" + systBaseline::spec_metPhiUp);
    const int nbJets_metPhiUp = isData? nbJets : tr->getVar<int>("cntCSVS" + systBaseline::spec_metPhiUp);
    const int nTops_metPhiUp = isData? nTops : tr->getVar<int>("nTopCandSortedCnt" + systBaseline::spec_metPhiUp);
    const double HT_metPhiUp = isData? HT : tr->getVar<double>("HT" + systBaseline::spec_metPhiUp);
    const bool passBaselineNoLepVeto_metPhiUp = isData? tr->getVar<bool>("passBaselineNoLepVeto"+spec) : tr->getVar<bool>("passBaselineNoLepVeto" + systBaseline::spec_metPhiUp);

    const double met_metPhiDn = isData? met : tr->getVar<double>("met");
    const double metphi_metPhiDn = isData? metphi : tr->getVar<double>("metphi_metPhiDn");
    const double MT2_metPhiDn = isData? MT2 : tr->getVar<double>("best_had_brJet_MT2" + systBaseline::spec_metPhiDn);
    const int nbJets_metPhiDn = isData? nbJets : tr->getVar<int>("cntCSVS" + systBaseline::spec_metPhiDn);
    const int nTops_metPhiDn = isData? nTops : tr->getVar<int>("nTopCandSortedCnt" + systBaseline::spec_metPhiDn);
    const double HT_metPhiDn = isData? HT : tr->getVar<double>("HT" + systBaseline::spec_metPhiDn);
    const bool passBaselineNoLepVeto_metPhiDn = isData? tr->getVar<bool>("passBaselineNoLepVeto"+spec) : tr->getVar<bool>("passBaselineNoLepVeto" + systBaseline::spec_metPhiDn);

    const double met_jecUp = isData? met : tr->getVar<double>("met");
    const double metphi_jecUp = isData? metphi : tr->getVar<double>("metphi");
    const double MT2_jecUp = isData? MT2 : tr->getVar<double>("best_had_brJet_MT2" + systBaseline::spec_jecUp);
    const int nbJets_jecUp = isData? nbJets : tr->getVar<int>("cntCSVS" + systBaseline::spec_jecUp);
    const int nTops_jecUp = isData? nTops : tr->getVar<int>("nTopCandSortedCnt" + systBaseline::spec_jecUp);
    const double HT_jecUp = isData? HT : tr->getVar<double>("HT" + systBaseline::spec_jecUp);
    const bool passBaselineNoLepVeto_jecUp = isData? tr->getVar<bool>("passBaselineNoLepVeto"+spec) : tr->getVar<bool>("passBaselineNoLepVeto" + systBaseline::spec_jecUp);

    const double met_jecDn = isData? met : tr->getVar<double>("met");
    const double metphi_jecDn = isData? metphi : tr->getVar<double>("metphi");
    const double MT2_jecDn = isData? MT2 : tr->getVar<double>("best_had_brJet_MT2" + systBaseline::spec_jecDn);
    const int nbJets_jecDn = isData? nbJets : tr->getVar<int>("cntCSVS" + systBaseline::spec_jecDn);
    const int nTops_jecDn = isData? nTops : tr->getVar<int>("nTopCandSortedCnt" + systBaseline::spec_jecDn);
    const double HT_jecDn = isData? HT : tr->getVar<double>("HT" + systBaseline::spec_jecDn);
    const bool passBaselineNoLepVeto_jecDn = isData? tr->getVar<bool>("passBaselineNoLepVeto"+spec) : tr->getVar<bool>("passBaselineNoLepVeto" + systBaseline::spec_jecDn);

    //mht calculation
    TLorentzVector MhtLVec;	      
    for(unsigned int ij=0; ij<jetsLVec.size(); ij++)
    {
      if( !AnaFunctions::jetPassCuts(jetsLVec[ij], AnaConsts::pt30Arr) ) continue;
      MhtLVec -= jetsLVec[ij];
    }
    const double Mht = MhtLVec.Pt();

    std::vector<TopObject*> topObjVec;
    if( ttPtr )
    {
      const TopTaggerResults& ttr = ttPtr->getResults();
      topObjVec = ttr.getTops();
    }

// Do genHT split ONLY for TTbar
    if( sampleString.Contains("TTbar") && useTTbarHTsamples )
    {
       if( !( (sampleString.Contains("HT") && genHT >=600) || (sampleString.Contains("Lep") && genHT < 600 ) ) ) continue;
    }

    h1_genHT->Fill(genHT, Lumiscale);

    bool passBaselineCS = passnJets && passdPhis && passMET && passBJets && passTagger && passHT && passMT2;
    bool passTrigger = true;        
    if(isData)
    {
      bool foundTrigger = false;
      for(unsigned it=0; it<TriggerNames.size(); it++)
      {
        if( sampleString.Contains("MET") )
        {
          if( TriggerNames[it].find("HLT_PFMET170_NoiseCleaned_v")!= string::npos || TriggerNames[it].find("HLT_PFMET170_JetIdCleaned_v") != string::npos || TriggerNames[it].find("HLT_PFMET170_HBHECleaned_v") != string::npos || TriggerNames[it].find("HLT_PFMET100_PFMHT100_IDTight_v") != string::npos || TriggerNames[it].find("HLT_PFMET110_PFMHT110_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMET120_PFMHT120_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMET130_PFMHT130_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMET140_PFMHT140_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMET150_PFMHT150_IDTight_v")!= string::npos
              || TriggerNames[it].find("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v") != string::npos
              || TriggerNames[it].find("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v") != string::npos
              || TriggerNames[it].find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") != string::npos
           )
          {
            if( PassTrigger[it] ) foundTrigger = true;
	  }
        }
        if( sampleString.Contains("SingleMuon") )
        {
          if( TriggerNames[it].find("HLT_Mu45_eta2p1_v") != string::npos || TriggerNames[it].find("HLT_Mu50_eta2p1_v") != string::npos || TriggerNames[it].find("HLT_Mu50_v") != string::npos || TriggerNames[it].find("HLT_Mu55_v") != string::npos )
          {
            if( PassTrigger[it] ) foundTrigger = true;
          }
        }
      }
      if( !foundTrigger ) passTrigger = false;
    }
    if(!passTrigger) continue;
    if(dospecialFix)
    {
      const std::vector<int> & specialFixtype = tr->getVec<int>("specialFixtype");
      const std::vector<TLorentzVector> & specialFixMuonsLVec = tr->getVec<TLorentzVector>("specialFixMuonsLVec");
      const std::vector<double> & specialFixMuonsCharge = tr->getVec<double>("specialFixMuonsCharge");
      if( !specialFixtype.empty() ) continue;
    }
    // Muon Control sample
    // The kinematic properties of the well-reconstructed, isolated muon                                                         
    vector<TLorentzVector> isomuonsLVec;
    vector<int> isomuonsIdxVec;
    for(unsigned int im=0; im<muonsLVec.size(); im++)
    {
      if( AnaFunctions::passMuon(muonsLVec.at(im), muonsMiniIso.at(im), muonsMtw.at(im), muonsFlagIDVec.at(im), AnaConsts::muonsMiniIsoArr))
      {
        isomuonsLVec.push_back(muonsLVec.at(im)); isomuonsIdxVec.push_back(im);
      }
    }
    // Require one and only one muon                                                                                         
    // Veto events with additional electrons (same veto criteria as baseline for electrons)                           
    if( nMuons == 1 && nElectrons == AnaConsts::nElectronsSel )
    {
      if( nMuons != isomuonsLVec.size() )
      {
        std::cout<<"ERROR ... mis-matching between veto muon and selected muon! Skipping..."<<std::endl; continue; 
      }

      //mtW cut
      const TLorentzVector muLVec = isomuonsLVec.at(0);
      const double mtw = calcMT(muLVec, metLVec);
      const bool pass_mtw = mtw<100 ? true : false;

      const double eta = muLVec.Eta(), pt = muLVec.Pt();
      const double abseta = std::abs(eta);
      double mu_id_SF = 1.0, mu_iso_SF = 1.0, mu_trk_SF = 1.0;
      double mu_id_SF_err = 0.0, mu_iso_SF_err = 0.0, mu_trk_SF_err = 0.0;
      if( dolepSF && !isData )
      {
        if( mu_mediumID_SF ){ mu_id_SF = mu_mediumID_SF->GetBinContent(mu_mediumID_SF->FindBin(pt, abseta)); if( mu_id_SF == 0 ) mu_id_SF = 1.0; } // very simple way dealing with out of range issue of the TH2D
        if( mu_miniISO_SF ){ mu_iso_SF = mu_miniISO_SF->GetBinContent(mu_miniISO_SF->FindBin(pt, abseta)); if( mu_iso_SF == 0 ) mu_iso_SF = 1.0; }
        if( pt < 10 && mu_trkptLT10_SF ){ mu_trk_SF = mu_trkptLT10_SF->Eval(eta); if( mu_trk_SF == 0 ) mu_trk_SF = 1.0; }
        if( pt >= 10 && mu_trkptGT10_SF ){ mu_trk_SF = mu_trkptGT10_SF->Eval(eta); if( mu_trk_SF == 0 ) mu_trk_SF = 1.0; }
        mu_id_SF_err = mu_id_SF * 0.03;
      }
      const double mu_SF = mu_id_SF * mu_iso_SF * mu_trk_SF;
      const double rel_mu_SF_err = sqrt(mu_id_SF_err*mu_id_SF_err/mu_id_SF/mu_id_SF + mu_iso_SF_err*mu_iso_SF_err/mu_iso_SF/mu_iso_SF + mu_trk_SF_err*mu_trk_SF_err/mu_trk_SF/mu_trk_SF);
      const double mu_SF_up = mu_SF*(1 + rel_mu_SF_err);
      const double mu_SF_dn = mu_SF*(1 - rel_mu_SF_err);

      const double corr_SF = bSF * isrWght * mu_SF;

      const double corr_bSF_up = bSF_up * isrWght * mu_SF;
      const double corr_bSF_down = bSF_down * isrWght * mu_SF;
  
      const double corr_isr_up = isr_up * bSF  * mu_SF;
      const double corr_isr_down = isr_down * bSF * mu_SF;

      const double corr_SF_up = bSF * isrWght * mu_SF_up;
      const double corr_SF_dn = bSF * isrWght * mu_SF_dn;

      FillDouble(myBaseHistgram.hvtxSize_mu_noCuts, vtxSize, Lumiscale);
      FillDouble(myBaseHistgram.hvtxSize_mu_noCuts_aft_puWght, vtxSize, Lumiscale*puWght);

      FillDouble(myBaseHistgram.hMET_mu_noCuts, met, Lumiscale*puWght*corr_SF);
      FillDouble(myBaseHistgram.hMET_mu_noCuts_no_corr_SF, met, Lumiscale*puWght);
      FillDouble(myBaseHistgram.hMHT_mu_noCuts, Mht, Lumiscale*puWght*corr_SF);
      FillDouble(myBaseHistgram.hMHT_mu_noCuts_no_corr_SF, Mht, Lumiscale*puWght);

      FillDouble(myBaseHistgram.hNbJets_mu_noCuts, nbJets, Lumiscale*puWght*corr_SF);
      FillDouble(myBaseHistgram.hNbJets_mu_noCuts_no_corr_SF, nbJets, Lumiscale*puWght);

      FillDouble(myBaseHistgram.hNJets_mu_noCuts, nJets, Lumiscale*puWght*corr_SF);
      FillDouble(myBaseHistgram.hNJets_mu_noCuts_no_corr_SF, nJets, Lumiscale*puWght);

      FillDouble(myBaseHistgram.hNTops_mu_noCuts, nTops, Lumiscale*puWght*corr_SF);
      FillDouble(myBaseHistgram.hNTops_mu_noCuts_no_corr_SF, nTops, Lumiscale*puWght);

      if( passMET ){
         FillDouble(myBaseHistgram.hMET_mu_passMET, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_mu_passMET_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_mu_passMET, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_mu_passMET_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_mu_passMET, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_mu_passMET_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_mu_passMET, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_mu_passMET_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passMET ){
         FillDouble(myBaseHistgram.hMET_mu_passnJets, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_mu_passnJets_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_mu_passnJets, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_mu_passnJets_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_mu_passnJets, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_mu_passnJets_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_mu_passnJets, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_mu_passnJets_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passdPhis && passMET ){
         FillDouble(myBaseHistgram.hMET_mu_passdPhis, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_mu_passdPhis_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_mu_passdPhis, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_mu_passdPhis_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_mu_passdPhis, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_mu_passdPhis_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_mu_passdPhis, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_mu_passdPhis_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passdPhis && passMET && passBJets ){
         FillDouble(myBaseHistgram.hMET_mu_passBJets, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_mu_passBJets_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_mu_passBJets, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_mu_passBJets_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_mu_passBJets, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_mu_passBJets_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_mu_passBJets, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_mu_passBJets_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passdPhis && passMET && passBJets && passTagger ){
         FillDouble(myBaseHistgram.hMET_mu_passTagger, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_mu_passTagger_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_mu_passTagger, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_mu_passTagger_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_mu_passTagger, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_mu_passTagger_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_mu_passTagger, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_mu_passTagger_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passdPhis && passMET && passBJets && passTagger && passHT ){
         FillDouble(myBaseHistgram.hMET_mu_passHT, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_mu_passHT_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_mu_passHT, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_mu_passHT_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_mu_passHT, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_mu_passHT_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_mu_passHT, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_mu_passHT_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passdPhis && passMET && passBJets && passTagger && passHT && passMT2 ){
         FillDouble(myBaseHistgram.hMET_mu_passMT2, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_mu_passMT2_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_mu_passMT2, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_mu_passMT2_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_mu_passMT2, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_mu_passMT2_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_mu_passMT2, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_mu_passMT2_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passdPhis && passMET && passBJets && passTagger && passHT && passMT2 && pass_mtw ){
         FillDouble(myBaseHistgram.hMET_mu_pass_mtw, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_mu_pass_mtw_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_mu_pass_mtw, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_mu_pass_mtw_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_mu_pass_mtw, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_mu_pass_mtw_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_mu_pass_mtw, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_mu_pass_mtw_no_corr_SF, nTops, Lumiscale*puWght);
      }

      //Dist.
      if(passBaselineCS && passNoiseEventFilter && pass_mtw)
      {
        int jSR = SB.find_Binning_Index(nbJets, nTops, MT2, met, HT);
        if( jSR!= -1 )
        {
          FillDouble(myBaseHistgram.hvtxSize_mu, vtxSize, Lumiscale);
          FillDouble(myBaseHistgram.hvtxSize_mu_aft_puWght, vtxSize, Lumiscale*puWght);

// Do NOT apply puWght on yields since for some bins with small MC stat, the PU weights can cause 
// large un-expected flucturation of the yields.
          myBaseHistgram.hYields_mu_no_corr_SF->Fill(jSR, Lumiscale);
          myBaseHistgram.hYields_mu_bSF->Fill(jSR, Lumiscale*bSF);
          myBaseHistgram.hYields_mu_isrWght->Fill(jSR, Lumiscale*isrWght);
          myBaseHistgram.hYields_mu_mu_SF->Fill(jSR, Lumiscale*mu_SF);
          myBaseHistgram.hYields_mu_bSF_isrWght->Fill(jSR, Lumiscale*bSF*isrWght);
          myBaseHistgram.hYields_mu_bSF_mu_SF->Fill(jSR, Lumiscale*bSF*mu_SF);
          myBaseHistgram.hYields_mu_isrWght_mu_SF->Fill(jSR, Lumiscale*isrWght*mu_SF);

          myBaseHistgram.hYields_mu->Fill(jSR, Lumiscale*corr_SF);
          myBaseHistgram.hYields_mu_aft_puWght->Fill(jSR, Lumiscale*puWght*corr_SF);
          myBaseHistgram.hYields_mu_aft_puSysup->Fill(jSR, Lumiscale*puSys_up*corr_SF);
          myBaseHistgram.hYields_mu_aft_puSysdown->Fill(jSR, Lumiscale*puSys_down*corr_SF);
	  //bSF systematics
          myBaseHistgram.hYields_mu_bSFup->Fill(jSR, Lumiscale*corr_bSF_up);
	  myBaseHistgram.hYields_mu_bSFdown->Fill(jSR, Lumiscale*corr_bSF_down);
	  //ISR Systematics
	  myBaseHistgram.hYields_mu_isrup->Fill(jSR, Lumiscale*corr_isr_up);
          myBaseHistgram.hYields_mu_isrdown->Fill(jSR, Lumiscale*corr_isr_down);
          //lepton SF systematics
          myBaseHistgram.hYields_mu_SFup->Fill(jSR, Lumiscale*corr_SF_up);
          myBaseHistgram.hYields_mu_SFdn->Fill(jSR, Lumiscale*corr_SF_dn);
          //Scale
          myBaseHistgram.hYields_mu_scaleUncup->Fill(jSR, Lumiscale*corr_SF*Scaled_Variations_Up);
          myBaseHistgram.hYields_mu_scaleUncdn->Fill(jSR, Lumiscale*corr_SF*Scaled_Variations_Down);
          //PDF
          myBaseHistgram.hYields_mu_pdfUncup->Fill(jSR, Lumiscale*corr_SF*NNPDF_From_Median_Up*NNPDF_From_Median_Central);
          myBaseHistgram.hYields_mu_pdfUnccen->Fill(jSR, Lumiscale*corr_SF*NNPDF_From_Median_Central);
          myBaseHistgram.hYields_mu_pdfUncdn->Fill(jSR, Lumiscale*corr_SF*NNPDF_From_Median_Down*NNPDF_From_Median_Central);
        }
  	  
        FillDouble(myBaseHistgram.hMET_mu, met, Lumiscale*puWght*corr_SF);
        FillDouble(myBaseHistgram.hMHT_mu, Mht, Lumiscale*puWght*corr_SF);
        FillDouble(myBaseHistgram.hMT2_mu, MT2, Lumiscale*puWght*corr_SF);
        FillInt(myBaseHistgram.hNbJets_mu, nbJets, Lumiscale*puWght*corr_SF);
        FillInt(myBaseHistgram.hNTops_mu, nTops, Lumiscale*puWght*corr_SF);	
  	  
        FillInt(myBaseHistgram.hNJets_mu, nJets, Lumiscale*puWght*corr_SF);
        FillDouble(myBaseHistgram.hHT_mu, HT, Lumiscale*puWght*corr_SF);
        FillDouble(myBaseHistgram.hdPhi0_mu, dPhiVec[0], Lumiscale*puWght*corr_SF);
        FillDouble(myBaseHistgram.hdPhi1_mu, dPhiVec[1], Lumiscale*puWght*corr_SF);
        FillDouble(myBaseHistgram.hdPhi2_mu, dPhiVec[2], Lumiscale*puWght*corr_SF);
  
        FillDouble(myBaseHistgram.hMET_mu_no_corr_SF, met, Lumiscale*puWght);
        FillDouble(myBaseHistgram.hMHT_mu_no_corr_SF, Mht, Lumiscale*puWght);
        FillDouble(myBaseHistgram.hMT2_mu_no_corr_SF, MT2, Lumiscale*puWght);
        FillInt(myBaseHistgram.hNbJets_mu_no_corr_SF, nbJets, Lumiscale*puWght);
        FillInt(myBaseHistgram.hNTops_mu_no_corr_SF, nTops, Lumiscale*puWght);	
  	  
        FillInt(myBaseHistgram.hNJets_mu_no_corr_SF, nJets, Lumiscale*puWght);
        FillDouble(myBaseHistgram.hHT_mu_no_corr_SF, HT, Lumiscale*puWght);
        FillDouble(myBaseHistgram.hdPhi0_mu_no_corr_SF, dPhiVec[0], Lumiscale*puWght);
        FillDouble(myBaseHistgram.hdPhi1_mu_no_corr_SF, dPhiVec[1], Lumiscale*puWght);
        FillDouble(myBaseHistgram.hdPhi2_mu_no_corr_SF, dPhiVec[2], Lumiscale*puWght);
  
        if( nbJets <=2 && nTops<=2 )
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, 250, met, HT); // use the lowest MT2 bin in (nb, ntop, met) to collapse the MT2 bins
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          auto it = std::find(cached_MT2_binIdx_mu_3DVec.begin(), cached_MT2_binIdx_mu_3DVec.end(), pseudo_SR);
          if( it != cached_MT2_binIdx_mu_3DVec.end() )
          {
            auto index = std::distance(cached_MT2_binIdx_mu_3DVec.begin(), it);
            FillDouble(cached_MT2_hist_mu_3DVec[index], MT2, Lumiscale*corr_SF);
          }else
          {
            cached_MT2_binIdx_mu_3DVec.push_back(pseudo_SR);
            char tmp_str[200];
            if( pseudo_binDef.met_hi_ != -1 )
            {
              sprintf(tmp_str, "muCS_MT2_3D_nb%d_nt%d_%3.0fmet%3.0f", nbJets, nTops, pseudo_binDef.met_lo_, pseudo_binDef.met_hi_);
            }else
            {
              sprintf(tmp_str, "muCS_MT2_3D_nb%d_nt%d_%3.0fmetInf", nbJets, nTops, pseudo_binDef.met_lo_);
            }
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 28, 200, 900);
            cached_MT2_hist_mu_3DVec.push_back(h1_tmp);
            FillDouble(cached_MT2_hist_mu_3DVec.back(), MT2, Lumiscale*corr_SF);
          }
        }else
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, MT2, met, 350); // use the lowest HT bin in (nb, ntop, met) to collapse the HT bins
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          const int nbJetsCopy = nbJets >=3 ? 3 : nbJets; const int nTopsCopy = nTops >=3 ? 3 : nTops;
          auto it = std::find(cached_HT_binIdx_mu_3DVec.begin(), cached_HT_binIdx_mu_3DVec.end(), pseudo_SR);
          if( it != cached_HT_binIdx_mu_3DVec.end() )
          {
            auto index = std::distance(cached_HT_binIdx_mu_3DVec.begin(), it);
            FillDouble(cached_HT_hist_mu_3DVec[index], HT, Lumiscale*corr_SF);
          }else
          {
            cached_HT_binIdx_mu_3DVec.push_back(pseudo_SR);
            char tmp_str[200];
            if( pseudo_binDef.met_hi_ != -1 )
            {
              sprintf(tmp_str, "muCS_HT_3D_nb%d_nt%d_%3.0fmet%3.0f", nbJetsCopy, nTopsCopy, pseudo_binDef.met_lo_, pseudo_binDef.met_hi_);
            }else
            {
              sprintf(tmp_str, "muCS_HT_3D_nb%d_nt%d_%3.0fmetInf", nbJetsCopy, nTopsCopy, pseudo_binDef.met_lo_);
            }
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 34, 300, 2000);
            cached_HT_hist_mu_3DVec.push_back(h1_tmp);
            FillDouble(cached_HT_hist_mu_3DVec.back(), HT, Lumiscale*corr_SF);
          }
        }
  
  // 2D in (nb, nt)
        if( nbJets <=2 && nTops<=2 )
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, 250, 300, HT); // use the lowest (MT2, met) bin in (nb, ntop) to collapse the MT2 and met
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          auto it = std::find(cached_MT2_binIdx_mu_2DVec.begin(), cached_MT2_binIdx_mu_2DVec.end(), pseudo_SR);
          if( it != cached_MT2_binIdx_mu_2DVec.end() )
          {
            auto index = std::distance(cached_MT2_binIdx_mu_2DVec.begin(), it);
            FillDouble(cached_MT2_hist_mu_2DVec[index], MT2, Lumiscale*corr_SF);
            FillDouble(cached_met_hist_mu_2DVec[index], met, Lumiscale*corr_SF);
            FillDouble(cached_nJets_hist_mu_2DVec[index], nJets, Lumiscale*corr_SF);
            for(auto top : topObjVec)
            {
              FillDouble(cached_recoTopPt_hist_mu_2DVec.back(), top->P().Pt(), Lumiscale*corr_SF);
            }
          }else
          {
            cached_MT2_binIdx_mu_2DVec.push_back(pseudo_SR);
            char tmp_str[200];
            sprintf(tmp_str, "muCS_MT2_2D_nb%d_nt%d", nbJets, nTops);
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 28, 200, 900);
            cached_MT2_hist_mu_2DVec.push_back(h1_tmp);
            FillDouble(cached_MT2_hist_mu_2DVec.back(), MT2, Lumiscale*corr_SF);
  
            sprintf(tmp_str, "muCS_met_2D_nb%d_nt%d", nbJets, nTops);
            TH1D * h1_tmp2 = new TH1D(tmp_str, tmp_str, 24, 250, 850);
            cached_met_hist_mu_2DVec.push_back(h1_tmp2);
            FillDouble(cached_met_hist_mu_2DVec.back(), met, Lumiscale*corr_SF);

            sprintf(tmp_str, "muCS_nJets_2D_nb%d_nt%d", nbJets, nTops);
            TH1D * h1_tmp3 = new TH1D(tmp_str, tmp_str, 6, 4, 10);
            cached_nJets_hist_mu_2DVec.push_back(h1_tmp3);
            FillDouble(cached_nJets_hist_mu_2DVec.back(), nJets, Lumiscale*corr_SF);

            sprintf(tmp_str, "muCS_recoTopPt_2D_nb%d_nt%d", nbJets, nTops);
            TH1D * h1_tmp4 = new TH1D(tmp_str, tmp_str, 40, 0, 1000);
            cached_recoTopPt_hist_mu_2DVec.push_back(h1_tmp4);
            for(auto top : topObjVec)
            {
              FillDouble(cached_recoTopPt_hist_mu_2DVec.back(), top->P().Pt(), Lumiscale*corr_SF);
            }
          }
        }else
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, MT2, 300, 350); // use the lowest (HT, met) bin in (nb, ntop) to collapse the HT and met
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          const int nbJetsCopy = nbJets >=3 ? 3 : nbJets; const int nTopsCopy = nTops >=3 ? 3 : nTops;
          auto it = std::find(cached_HT_binIdx_mu_2DVec.begin(), cached_HT_binIdx_mu_2DVec.end(), pseudo_SR);
          if( it != cached_HT_binIdx_mu_2DVec.end() )
          {
            auto index = std::distance(cached_HT_binIdx_mu_2DVec.begin(), it);
            FillDouble(cached_HT_hist_mu_2DVec[index], HT, Lumiscale*corr_SF);
          }else
          {
            cached_HT_binIdx_mu_2DVec.push_back(pseudo_SR);
            char tmp_str[200];
            sprintf(tmp_str, "muCS_HT_2D_nb%d_nt%d", nbJetsCopy, nTopsCopy);
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 34, 300, 2000);
            cached_HT_hist_mu_2DVec.push_back(h1_tmp);
            FillDouble(cached_HT_hist_mu_2DVec.back(), HT, Lumiscale*corr_SF);
          }
        }
      }

      if( !isData )
      {
// metMagUp
        TLorentzVector metLVec_metMagUp; metLVec_metMagUp.SetPtEtaPhiM(met_metMagUp, 0, metphi_metMagUp, 0);
        const double mtw_metMagUp = calcMT(muLVec, metLVec_metMagUp);
        const bool pass_mtw_metMagUp = mtw_metMagUp<100 ? true : false;
        if( passBaselineNoLepVeto_metMagUp && passNoiseEventFilter && pass_mtw_metMagUp)
        {
           int jSR = SB.find_Binning_Index(nbJets_metMagUp, nTops_metMagUp, MT2_metMagUp, met_metMagUp, HT_metMagUp);
           if( jSR!= -1 ) myBaseHistgram.hYields_mu_metMagUp->Fill(jSR, Lumiscale*corr_SF);
        }
      
// metMagDn
        TLorentzVector metLVec_metMagDn; metLVec_metMagDn.SetPtEtaPhiM(met_metMagDn, 0, metphi_metMagDn, 0);
        const double mtw_metMagDn = calcMT(muLVec, metLVec_metMagDn);
        const bool pass_mtw_metMagDn = mtw_metMagDn<100 ? true : false;
        if( passBaselineNoLepVeto_metMagDn && passNoiseEventFilter && pass_mtw_metMagDn)
        {
           int jSR = SB.find_Binning_Index(nbJets_metMagDn, nTops_metMagDn, MT2_metMagDn, met_metMagDn, HT_metMagDn);
           if( jSR!= -1 ) myBaseHistgram.hYields_mu_metMagDn->Fill(jSR, Lumiscale*corr_SF);
        }
      
// metPhiUp
        TLorentzVector metLVec_metPhiUp; metLVec_metPhiUp.SetPtEtaPhiM(met_metPhiUp, 0, metphi_metPhiUp, 0);
        const double mtw_metPhiUp = calcMT(muLVec, metLVec_metPhiUp);
        const bool pass_mtw_metPhiUp = mtw_metPhiUp<100 ? true : false;
        if( passBaselineNoLepVeto_metPhiUp && passNoiseEventFilter && pass_mtw_metPhiUp)
        {
           int jSR = SB.find_Binning_Index(nbJets_metPhiUp, nTops_metPhiUp, MT2_metPhiUp, met_metPhiUp, HT_metPhiUp);
           if( jSR!= -1 ) myBaseHistgram.hYields_mu_metPhiUp->Fill(jSR, Lumiscale*corr_SF);
        }
     
// metPhiDn
        TLorentzVector metLVec_metPhiDn; metLVec_metPhiDn.SetPtEtaPhiM(met_metPhiDn, 0, metphi_metPhiDn, 0);
        const double mtw_metPhiDn = calcMT(muLVec, metLVec_metPhiDn);
        const bool pass_mtw_metPhiDn = mtw_metPhiDn<100 ? true : false;
        if( passBaselineNoLepVeto_metPhiDn && passNoiseEventFilter && pass_mtw_metPhiDn)
        {
           int jSR = SB.find_Binning_Index(nbJets_metPhiDn, nTops_metPhiDn, MT2_metPhiDn, met_metPhiDn, HT_metPhiDn);
           if( jSR!= -1 ) myBaseHistgram.hYields_mu_metPhiDn->Fill(jSR, Lumiscale*corr_SF);
        }
    
// jecUp
        TLorentzVector metLVec_jecUp; metLVec_jecUp.SetPtEtaPhiM(met_jecUp, 0, metphi_jecUp, 0);
        const double mtw_jecUp = calcMT(muLVec, metLVec_jecUp);
        const bool pass_mtw_jecUp = mtw_jecUp<100 ? true : false;
        if( passBaselineNoLepVeto_jecUp && passNoiseEventFilter && pass_mtw_jecUp)
        {
           int jSR = SB.find_Binning_Index(nbJets_jecUp, nTops_jecUp, MT2_jecUp, met_jecUp, HT_jecUp);
           if( jSR!= -1 ) myBaseHistgram.hYields_mu_jecUp->Fill(jSR, Lumiscale*corr_SF);
        }
   
// jecDn
        TLorentzVector metLVec_jecDn; metLVec_jecDn.SetPtEtaPhiM(met_jecDn, 0, metphi_jecDn, 0);
        const double mtw_jecDn = calcMT(muLVec, metLVec_jecDn);
        const bool pass_mtw_jecDn = mtw_jecDn<100 ? true : false;
        if( passBaselineNoLepVeto_jecDn && passNoiseEventFilter && pass_mtw_jecDn)
        {
           int jSR = SB.find_Binning_Index(nbJets_jecDn, nTops_jecDn, MT2_jecDn, met_jecDn, HT_jecDn);
           if( jSR!= -1 ) myBaseHistgram.hYields_mu_jecDn->Fill(jSR, Lumiscale*corr_SF);
        }
      }
    }//end of muon CS

    // Electron Control sample
    //The kinematic properties of the well-reconstructed, isolated muon                                                         
    vector<TLorentzVector> isoelesLVec;
    vector<int> isoelesIdxVec;
    for(unsigned int im=0; im<elesLVec.size(); im++)
    {
      if( AnaFunctions::passElectron(elesLVec.at(im), elesMiniIso.at(im), elesMtw.at(im), elesisEB.at(im), elesFlagIDVec.at(im), AnaConsts::elesMiniIsoArr))
      {
        isoelesLVec.push_back(elesLVec.at(im)); isoelesIdxVec.push_back(im); 
      }
    }
    // Require one and only one electron                                                                                         
    // Veto events with additional muons (same veto criteria as baseline for muons)                           
    if( nElectrons == 1 && nMuons == AnaConsts::nMuonsSel )
    {
      if( nElectrons != isoelesLVec.size() )
      { 
        std::cout<<"ERROR ... mis-matching between veto ele and selected ele! Skipping..."<<std::endl; continue; 
      }
	
      //mtW cut
      const TLorentzVector eleLVec = isoelesLVec.at(0);
      const double elemtw = calcMT(eleLVec, metLVec);
      const bool pass_mtwele = elemtw<100 ? true : false;
	
      const double eta = eleLVec.Eta(), pt = eleLVec.Pt();
      const double abseta = std::abs(eta);
      double ele_id_SF = 1.0, ele_iso_SF = 1.0, ele_trk_SF = 1.0;
      double ele_id_SF_err = 0.0, ele_iso_SF_err = 0.0, ele_trk_SF_err = 0.0;
      if( dolepSF && !isData )
      {
        if( ele_VetoID_SF )
        {
          ele_id_SF = ele_VetoID_SF->GetBinContent(ele_VetoID_SF->FindBin(pt, abseta));
          ele_id_SF_err = ele_VetoID_SF->GetBinError(ele_VetoID_SF->FindBin(pt, abseta));
          if( ele_id_SF == 0 ){ ele_id_SF = 1.0; ele_id_SF_err = 0.0; }
        } // very simple way dealing with out of range issue of the TH2D
        if( ele_miniISO_SF )
        {
          ele_iso_SF = ele_miniISO_SF->GetBinContent(ele_miniISO_SF->FindBin(pt, abseta));
          ele_iso_SF_err = ele_miniISO_SF->GetBinError(ele_miniISO_SF->FindBin(pt, abseta));
          if( ele_iso_SF == 0 ){ ele_iso_SF = 1.0; ele_iso_SF_err = 0.0; }
        }
        if( ele_trkpt_SF )
        {
          ele_trk_SF = ele_trkpt_SF->GetBinContent(ele_trkpt_SF->FindBin(eta, pt));
          ele_trk_SF_err = pt<20 || pt>80? 0.01: 0.00;
          if( ele_trk_SF == 0 ){ ele_trk_SF = 1.0; ele_trk_SF_err = 0.0; }
        }
      }
      const double ele_SF = ele_id_SF * ele_iso_SF * ele_trk_SF;
      const double rel_ele_SF_err = sqrt(ele_id_SF_err*ele_id_SF_err/ele_id_SF/ele_id_SF + ele_iso_SF_err*ele_iso_SF_err/ele_iso_SF/ele_iso_SF + ele_trk_SF_err*ele_trk_SF_err/ele_trk_SF/ele_trk_SF);
      const double ele_SF_up = ele_SF*(1 + rel_ele_SF_err);
      const double ele_SF_dn = ele_SF*(1 - rel_ele_SF_err);
	
      const double corr_SF = bSF * isrWght * ele_SF;

      const double corr_bSF_up = bSF_up * isrWght * ele_SF;
      const double corr_bSF_down = bSF_down * isrWght * ele_SF;
	
      const double corr_isr_up = bSF * isr_up * ele_SF;
      const double corr_isr_down = bSF * isr_down * ele_SF;

      const double corr_SF_up = bSF * isrWght * ele_SF_up;
      const double corr_SF_dn = bSF * isrWght * ele_SF_dn;

      FillDouble(myBaseHistgram.hvtxSize_el_noCuts, vtxSize, Lumiscale);
      FillDouble(myBaseHistgram.hvtxSize_el_noCuts_aft_puWght, vtxSize, Lumiscale*puWght);
	
      FillDouble(myBaseHistgram.hMET_el_noCuts, met, Lumiscale*puWght*corr_SF);
      FillDouble(myBaseHistgram.hMET_el_noCuts_no_corr_SF, met, Lumiscale*puWght);
      FillDouble(myBaseHistgram.hMHT_el_noCuts, Mht, Lumiscale*puWght*corr_SF);
      FillDouble(myBaseHistgram.hMHT_el_noCuts_no_corr_SF, Mht, Lumiscale*puWght);

      FillDouble(myBaseHistgram.hNbJets_el_noCuts, nbJets, Lumiscale*puWght*corr_SF);
      FillDouble(myBaseHistgram.hNbJets_el_noCuts_no_corr_SF, nbJets, Lumiscale*puWght);

      FillDouble(myBaseHistgram.hNJets_el_noCuts, nJets, Lumiscale*puWght*corr_SF);
      FillDouble(myBaseHistgram.hNJets_el_noCuts_no_corr_SF, nJets, Lumiscale*puWght);

      FillDouble(myBaseHistgram.hNTops_el_noCuts, nTops, Lumiscale*puWght*corr_SF);
      FillDouble(myBaseHistgram.hNTops_el_noCuts_no_corr_SF, nTops, Lumiscale*puWght);

      if( passMET ){
         FillDouble(myBaseHistgram.hMET_el_passMET, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_el_passMET_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_el_passMET, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_el_passMET_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_el_passMET, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_el_passMET_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_el_passMET, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_el_passMET_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passMET ){
         FillDouble(myBaseHistgram.hMET_el_passnJets, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_el_passnJets_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_el_passnJets, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_el_passnJets_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_el_passnJets, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_el_passnJets_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_el_passnJets, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_el_passnJets_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passdPhis && passMET ){
         FillDouble(myBaseHistgram.hMET_el_passdPhis, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_el_passdPhis_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_el_passdPhis, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_el_passdPhis_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_el_passdPhis, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_el_passdPhis_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_el_passdPhis, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_el_passdPhis_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passdPhis && passMET && passBJets ){
         FillDouble(myBaseHistgram.hMET_el_passBJets, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_el_passBJets_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_el_passBJets, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_el_passBJets_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_el_passBJets, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_el_passBJets_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_el_passBJets, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_el_passBJets_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passdPhis && passMET && passBJets && passTagger ){
         FillDouble(myBaseHistgram.hMET_el_passTagger, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_el_passTagger_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_el_passTagger, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_el_passTagger_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_el_passTagger, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_el_passTagger_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_el_passTagger, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_el_passTagger_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passdPhis && passMET && passBJets && passTagger && passHT ){
         FillDouble(myBaseHistgram.hMET_el_passHT, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_el_passHT_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_el_passHT, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_el_passHT_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_el_passHT, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_el_passHT_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_el_passHT, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_el_passHT_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passdPhis && passMET && passBJets && passTagger && passHT && passMT2 ){
         FillDouble(myBaseHistgram.hMET_el_passMT2, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_el_passMT2_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_el_passMT2, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_el_passMT2_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_el_passMT2, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_el_passMT2_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_el_passMT2, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_el_passMT2_no_corr_SF, nTops, Lumiscale*puWght);
      }
      if( passnJets && passdPhis && passMET && passBJets && passTagger && passHT && passMT2 && pass_mtwele ){
         FillDouble(myBaseHistgram.hMET_el_pass_mtw, met, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hMET_el_pass_mtw_no_corr_SF, met, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNbJets_el_pass_mtw, nbJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNbJets_el_pass_mtw_no_corr_SF, nbJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNJets_el_pass_mtw, nJets, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNJets_el_pass_mtw_no_corr_SF, nJets, Lumiscale*puWght);

         FillDouble(myBaseHistgram.hNTops_el_pass_mtw, nTops, Lumiscale*puWght*corr_SF);
         FillDouble(myBaseHistgram.hNTops_el_pass_mtw_no_corr_SF, nTops, Lumiscale*puWght);
      }

      //Dist.
      if(passBaselineCS && passNoiseEventFilter && pass_mtwele)
      {
        int kSR = SB.find_Binning_Index(nbJets, nTops, MT2, met, HT);
        if( kSR!= -1 )
        {
          FillDouble(myBaseHistgram.hvtxSize_el, vtxSize, Lumiscale);
          FillDouble(myBaseHistgram.hvtxSize_el_aft_puWght, vtxSize, Lumiscale*puWght);

// Do NOT apply puWght on yields since for some bins with small MC stat, the PU weights can cause 
// large un-expected flucturation of the yields.
          myBaseHistgram.hYields_el_no_corr_SF->Fill(kSR, Lumiscale);
          myBaseHistgram.hYields_el_bSF->Fill(kSR, Lumiscale*bSF);
          myBaseHistgram.hYields_el_isrWght->Fill(kSR, Lumiscale*isrWght);
          myBaseHistgram.hYields_el_ele_SF->Fill(kSR, Lumiscale*ele_SF);
          myBaseHistgram.hYields_el_bSF_isrWght->Fill(kSR, Lumiscale*bSF*isrWght);
          myBaseHistgram.hYields_el_bSF_ele_SF->Fill(kSR, Lumiscale*bSF*ele_SF);
          myBaseHistgram.hYields_el_isrWght_ele_SF->Fill(kSR, Lumiscale*isrWght*ele_SF);

          myBaseHistgram.hYields_el->Fill(kSR, Lumiscale*corr_SF);
          myBaseHistgram.hYields_el_aft_puWght->Fill(kSR, Lumiscale*puWght*corr_SF);
          myBaseHistgram.hYields_el_aft_puSysup->Fill(kSR, Lumiscale*puSys_up*corr_SF);
          myBaseHistgram.hYields_el_aft_puSysdown->Fill(kSR, Lumiscale*puSys_down*corr_SF);
	  //bSF systematics
          myBaseHistgram.hYields_el_bSFup->Fill(kSR, Lumiscale*corr_bSF_up);
	  myBaseHistgram.hYields_el_bSFdown->Fill(kSR, Lumiscale*corr_bSF_down);

	  myBaseHistgram.hYields_el_isrup->Fill(kSR, Lumiscale*corr_isr_up);
          myBaseHistgram.hYields_el_isrdown->Fill(kSR, Lumiscale*corr_isr_down);

	  myBaseHistgram.hYields_el_SFup->Fill(kSR, Lumiscale*corr_SF_up);
          myBaseHistgram.hYields_el_SFdn->Fill(kSR, Lumiscale*corr_SF_dn);
          //Scale
          myBaseHistgram.hYields_el_scaleUncup->Fill(kSR, Lumiscale*corr_SF*Scaled_Variations_Up);
          myBaseHistgram.hYields_el_scaleUncdn->Fill(kSR, Lumiscale*corr_SF*Scaled_Variations_Down);
          //PDF
          myBaseHistgram.hYields_el_pdfUncup->Fill(kSR, Lumiscale*corr_SF*NNPDF_From_Median_Up);
          myBaseHistgram.hYields_el_pdfUnccen->Fill(kSR, Lumiscale*corr_SF*NNPDF_From_Median_Central);
          myBaseHistgram.hYields_el_pdfUncdn->Fill(kSR, Lumiscale*corr_SF*NNPDF_From_Median_Down);
        }
  	  
        FillDouble(myBaseHistgram.hMET_el, met, Lumiscale*puWght*corr_SF);
        FillDouble(myBaseHistgram.hMHT_el, Mht, Lumiscale*puWght*corr_SF);
        FillDouble(myBaseHistgram.hMT2_el, MT2, Lumiscale*puWght*corr_SF);
        FillInt(myBaseHistgram.hNbJets_el, nbJets, Lumiscale*puWght*corr_SF);
        FillInt(myBaseHistgram.hNTops_el, nTops, Lumiscale*puWght*corr_SF);	
  	  
        FillInt(myBaseHistgram.hNJets_el, nJets, Lumiscale*puWght*corr_SF);
        FillDouble(myBaseHistgram.hHT_el, HT, Lumiscale*puWght*corr_SF);
        FillDouble(myBaseHistgram.hdPhi0_el, dPhiVec[0], Lumiscale*puWght*corr_SF);
        FillDouble(myBaseHistgram.hdPhi1_el, dPhiVec[1], Lumiscale*puWght*corr_SF);
        FillDouble(myBaseHistgram.hdPhi2_el, dPhiVec[2], Lumiscale*puWght*corr_SF);

        FillDouble(myBaseHistgram.hMET_el_no_corr_SF, met, Lumiscale*puWght);
        FillDouble(myBaseHistgram.hMHT_el_no_corr_SF, Mht, Lumiscale*puWght);
        FillDouble(myBaseHistgram.hMT2_el_no_corr_SF, MT2, Lumiscale*puWght);
        FillInt(myBaseHistgram.hNbJets_el_no_corr_SF, nbJets, Lumiscale*puWght);
        FillInt(myBaseHistgram.hNTops_el_no_corr_SF, nTops, Lumiscale*puWght);	
  	  
        FillInt(myBaseHistgram.hNJets_el_no_corr_SF, nJets, Lumiscale*puWght);
        FillDouble(myBaseHistgram.hHT_el_no_corr_SF, HT, Lumiscale*puWght);
        FillDouble(myBaseHistgram.hdPhi0_el_no_corr_SF, dPhiVec[0], Lumiscale*puWght);
        FillDouble(myBaseHistgram.hdPhi1_el_no_corr_SF, dPhiVec[1], Lumiscale*puWght);
        FillDouble(myBaseHistgram.hdPhi2_el_no_corr_SF, dPhiVec[2], Lumiscale*puWght);
  
        if( nbJets <=2 && nTops<=2 )
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, 250, met, HT); // use the lowest MT2 bin in (nb, ntop, met) to collapse the MT2 bins
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          auto it = std::find(cached_MT2_binIdx_el_3DVec.begin(), cached_MT2_binIdx_el_3DVec.end(), pseudo_SR);
          if( it != cached_MT2_binIdx_el_3DVec.end() )
          {
            auto index = std::distance(cached_MT2_binIdx_el_3DVec.begin(), it);
            FillDouble(cached_MT2_hist_el_3DVec[index], MT2, Lumiscale*corr_SF);
          }else
          {
            cached_MT2_binIdx_el_3DVec.push_back(pseudo_SR);
            char tmp_str[200];
            if( pseudo_binDef.met_hi_ != -1 )
            {
              sprintf(tmp_str, "eleCS_MT2_3D_nb%d_nt%d_%3.0fmet%3.0f", nbJets, nTops, pseudo_binDef.met_lo_, pseudo_binDef.met_hi_);
            }else
            {
              sprintf(tmp_str, "eleCS_MT2_3D_nb%d_nt%d_%3.0fmetInf", nbJets, nTops, pseudo_binDef.met_lo_);
            }
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 28, 200, 900);
            cached_MT2_hist_el_3DVec.push_back(h1_tmp);
            FillDouble(cached_MT2_hist_el_3DVec.back(), MT2, Lumiscale*corr_SF);
          }
        }else
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, MT2, met, 350); // use the lowest HT bin in (nb, ntop, met) to collapse the HT bins
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          const int nbJetsCopy = nbJets >=3 ? 3 : nbJets; const int nTopsCopy = nTops >=3 ? 3 : nTops;
          auto it = std::find(cached_HT_binIdx_el_3DVec.begin(), cached_HT_binIdx_el_3DVec.end(), pseudo_SR);
          if( it != cached_HT_binIdx_el_3DVec.end() )
          {
            auto index = std::distance(cached_HT_binIdx_el_3DVec.begin(), it);
            FillDouble(cached_HT_hist_el_3DVec[index], HT, Lumiscale*corr_SF);
          }else
          {
            cached_HT_binIdx_el_3DVec.push_back(pseudo_SR);
            char tmp_str[200];
            if( pseudo_binDef.met_hi_ != -1 )
            {
              sprintf(tmp_str, "eleCS_HT_3D_nb%d_nt%d_%3.0fmet%3.0f", nbJetsCopy, nTopsCopy, pseudo_binDef.met_lo_, pseudo_binDef.met_hi_);
            }else
            {
              sprintf(tmp_str, "eleCS_HT_3D_nb%d_nt%d_%3.0fmetInf", nbJetsCopy, nTopsCopy, pseudo_binDef.met_lo_);
            }
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 34, 300, 2000);
            cached_HT_hist_el_3DVec.push_back(h1_tmp);
            FillDouble(cached_HT_hist_el_3DVec.back(), HT, Lumiscale*corr_SF);
          }
        }
  
  // 2D in (nb, nt)
        if( nbJets <=2 && nTops<=2 )
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, 250, 300, HT); // use the lowest (MT2, met) bin in (nb, ntop) to collapse the MT2 and met
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          auto it = std::find(cached_MT2_binIdx_el_2DVec.begin(), cached_MT2_binIdx_el_2DVec.end(), pseudo_SR);
          if( it != cached_MT2_binIdx_el_2DVec.end() )
          {
            auto index = std::distance(cached_MT2_binIdx_el_2DVec.begin(), it);
            FillDouble(cached_MT2_hist_el_2DVec[index], MT2, Lumiscale*corr_SF);
            FillDouble(cached_met_hist_el_2DVec[index], met, Lumiscale*corr_SF);
            FillDouble(cached_nJets_hist_el_2DVec[index], nJets, Lumiscale*corr_SF);
            for(auto top : topObjVec)
            {
              FillDouble(cached_recoTopPt_hist_el_2DVec.back(), top->P().Pt(), Lumiscale*corr_SF);
            }
          }else
          {
            cached_MT2_binIdx_el_2DVec.push_back(pseudo_SR);
            char tmp_str[200];
            sprintf(tmp_str, "eleCS_MT2_2D_nb%d_nt%d", nbJets, nTops);
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 28, 200, 900);
            cached_MT2_hist_el_2DVec.push_back(h1_tmp);
            FillDouble(cached_MT2_hist_el_2DVec.back(), MT2, Lumiscale*corr_SF);
  
            sprintf(tmp_str, "eleCS_met_2D_nb%d_nt%d", nbJets, nTops);
            TH1D * h1_tmp2 = new TH1D(tmp_str, tmp_str, 24, 250, 850);
            cached_met_hist_el_2DVec.push_back(h1_tmp2);
            FillDouble(cached_met_hist_el_2DVec.back(), met, Lumiscale*corr_SF);

            sprintf(tmp_str, "eleCS_nJets_2D_nb%d_nt%d", nbJets, nTops);
            TH1D * h1_tmp3 = new TH1D(tmp_str, tmp_str, 6, 4, 10);
            cached_nJets_hist_el_2DVec.push_back(h1_tmp3);
            FillDouble(cached_nJets_hist_el_2DVec.back(), nJets, Lumiscale*corr_SF);

            sprintf(tmp_str, "eleCS_recoTopPt_2D_nb%d_nt%d", nbJets, nTops);
            TH1D * h1_tmp4 = new TH1D(tmp_str, tmp_str, 40, 0, 1000);
            cached_recoTopPt_hist_el_2DVec.push_back(h1_tmp4);
            for(auto top : topObjVec)
            {
              FillDouble(cached_recoTopPt_hist_el_2DVec.back(), top->P().Pt(), Lumiscale*corr_SF);
            }
          }
        }else
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, MT2, 300, 350); // use the lowest (HT, met) bin in (nb, ntop) to collapse the HT and met
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          const int nbJetsCopy = nbJets >=3 ? 3 : nbJets; const int nTopsCopy = nTops >=3 ? 3 : nTops;
          auto it = std::find(cached_HT_binIdx_el_2DVec.begin(), cached_HT_binIdx_el_2DVec.end(), pseudo_SR);
          if( it != cached_HT_binIdx_el_2DVec.end() )
          {
            auto index = std::distance(cached_HT_binIdx_el_2DVec.begin(), it);
            FillDouble(cached_HT_hist_el_2DVec[index], HT, Lumiscale*corr_SF);
          }else
          {
            cached_HT_binIdx_el_2DVec.push_back(pseudo_SR);
            char tmp_str[200];
            sprintf(tmp_str, "eleCS_HT_2D_nb%d_nt%d", nbJetsCopy, nTopsCopy);
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 34, 300, 2000);
            cached_HT_hist_el_2DVec.push_back(h1_tmp);
            FillDouble(cached_HT_hist_el_2DVec.back(), HT, Lumiscale*corr_SF);
          }
        }
      }

      if( !isData )
      {
// metMagUp
        TLorentzVector metLVec_metMagUp; metLVec_metMagUp.SetPtEtaPhiM(met_metMagUp, 0, metphi_metMagUp, 0);
        const double mtw_metMagUp = calcMT(eleLVec, metLVec_metMagUp);
        const bool pass_mtw_metMagUp = mtw_metMagUp<100 ? true : false;
        if( passBaselineNoLepVeto_metMagUp && passNoiseEventFilter && pass_mtw_metMagUp)
        {
           int jSR = SB.find_Binning_Index(nbJets_metMagUp, nTops_metMagUp, MT2_metMagUp, met_metMagUp, HT_metMagUp);
           if( jSR!= -1 ) myBaseHistgram.hYields_el_metMagUp->Fill(jSR, Lumiscale*corr_SF);
        }
      
// metMagDn
        TLorentzVector metLVec_metMagDn; metLVec_metMagDn.SetPtEtaPhiM(met_metMagDn, 0, metphi_metMagDn, 0);
        const double mtw_metMagDn = calcMT(eleLVec, metLVec_metMagDn);
        const bool pass_mtw_metMagDn = mtw_metMagDn<100 ? true : false;
        if( passBaselineNoLepVeto_metMagDn && passNoiseEventFilter && pass_mtw_metMagDn)
        {
           int jSR = SB.find_Binning_Index(nbJets_metMagDn, nTops_metMagDn, MT2_metMagDn, met_metMagDn, HT_metMagDn);
           if( jSR!= -1 ) myBaseHistgram.hYields_el_metMagDn->Fill(jSR, Lumiscale*corr_SF);
        }
     
// metPhiUp
        TLorentzVector metLVec_metPhiUp; metLVec_metPhiUp.SetPtEtaPhiM(met_metPhiUp, 0, metphi_metPhiUp, 0);
        const double mtw_metPhiUp = calcMT(eleLVec, metLVec_metPhiUp);
        const bool pass_mtw_metPhiUp = mtw_metPhiUp<100 ? true : false;
        if( passBaselineNoLepVeto_metPhiUp && passNoiseEventFilter && pass_mtw_metPhiUp)
        {
           int jSR = SB.find_Binning_Index(nbJets_metPhiUp, nTops_metPhiUp, MT2_metPhiUp, met_metPhiUp, HT_metPhiUp);
           if( jSR!= -1 ) myBaseHistgram.hYields_el_metPhiUp->Fill(jSR, Lumiscale*corr_SF);
        }
    
// metPhiDn
        TLorentzVector metLVec_metPhiDn; metLVec_metPhiDn.SetPtEtaPhiM(met_metPhiDn, 0, metphi_metPhiDn, 0);
        const double mtw_metPhiDn = calcMT(eleLVec, metLVec_metPhiDn);
        const bool pass_mtw_metPhiDn = mtw_metPhiDn<100 ? true : false;
        if( passBaselineNoLepVeto_metPhiDn && passNoiseEventFilter && pass_mtw_metPhiDn)
        {
           int jSR = SB.find_Binning_Index(nbJets_metPhiDn, nTops_metPhiDn, MT2_metPhiDn, met_metPhiDn, HT_metPhiDn);
           if( jSR!= -1 ) myBaseHistgram.hYields_el_metPhiDn->Fill(jSR, Lumiscale*corr_SF);
        }
   
// jecUp
        TLorentzVector metLVec_jecUp; metLVec_jecUp.SetPtEtaPhiM(met_jecUp, 0, metphi_jecUp, 0);
        const double mtw_jecUp = calcMT(eleLVec, metLVec_jecUp);
        const bool pass_mtw_jecUp = mtw_jecUp<100 ? true : false;
        if( passBaselineNoLepVeto_jecUp && passNoiseEventFilter && pass_mtw_jecUp)
        {
           int jSR = SB.find_Binning_Index(nbJets_jecUp, nTops_jecUp, MT2_jecUp, met_jecUp, HT_jecUp);
           if( jSR!= -1 ) myBaseHistgram.hYields_el_jecUp->Fill(jSR, Lumiscale*corr_SF);
        }
  
// jecDn
        TLorentzVector metLVec_jecDn; metLVec_jecDn.SetPtEtaPhiM(met_jecDn, 0, metphi_jecDn, 0);
        const double mtw_jecDn = calcMT(eleLVec, metLVec_jecDn);
        const bool pass_mtw_jecDn = mtw_jecDn<100 ? true : false;
        if( passBaselineNoLepVeto_jecDn && passNoiseEventFilter && pass_mtw_jecDn)
        {
           int jSR = SB.find_Binning_Index(nbJets_jecDn, nTops_jecDn, MT2_jecDn, met_jecDn, HT_jecDn);
           if( jSR!= -1 ) myBaseHistgram.hYields_el_jecDn->Fill(jSR, Lumiscale*corr_SF);
        }
      }

    }//end of electron CS
  }//event loop
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->cd();
  if( !isData )
  {
    h1_genHT->Write();

    for(auto hist: cached_MT2_hist_mu_3DVec) hist->Write();
    for(auto hist: cached_HT_hist_mu_3DVec) hist->Write();
    for(auto hist: cached_MT2_hist_el_3DVec) hist->Write();
    for(auto hist: cached_HT_hist_el_3DVec) hist->Write();

    for(auto hist: cached_MT2_hist_mu_2DVec) hist->Write();
    for(auto hist: cached_HT_hist_mu_2DVec) hist->Write();
    for(auto hist: cached_met_hist_mu_2DVec) hist->Write();
    for(auto hist: cached_nJets_hist_mu_2DVec) hist->Write();
    for(auto hist: cached_recoTopPt_hist_mu_2DVec) hist->Write();

    for(auto hist: cached_MT2_hist_el_2DVec) hist->Write();
    for(auto hist: cached_HT_hist_el_2DVec) hist->Write();
    for(auto hist: cached_met_hist_el_2DVec) hist->Write();
    for(auto hist: cached_nJets_hist_el_2DVec) hist->Write();
    for(auto hist: cached_recoTopPt_hist_el_2DVec) hist->Write();
  }
  (myBaseHistgram.oFile)->Write();
  fChain->Reset();  
  return 0;
}
