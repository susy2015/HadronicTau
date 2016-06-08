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
#include "TStopwatch.h"
#include "TString.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "DataMC.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TChain.h"
#include "TauResponse.h"
#include "utils.h"
#include "Efficiency.h"
#include "TRandom3.h"

using namespace std;

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 6)
    {
      std::cerr <<"Please give 5 arguments "<<"SubsampleName"<<" Input Template" <<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Closure TTbarInc TTbarInc_Template.root 1000 0 1" << std::endl;
      return -1;
    }
  const char *subsamplename = argv[1];
  const char *respTempl = argv[2];
  const char *Maxevent = argv[3];
  const  char *Stratfile = argv[4];
  const  char *Filerun = argv[5];
  const  int startfile = std::atoi(Stratfile);
  const int filerun = std::atoi(Filerun);
  bool isData = false;
  TChain *fChain = 0;
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(subsamplename, startfile);
  if(argc==7){
    const char *samplename = argv[6];
     if(!FillChain(fChain, samplename, subsamplename, startfile, filerun))
       {
	 std::cerr << "Cannot get the tree " << std::endl;
       }
  }
  else {

    if(!FillChain(fChain, "null", subsamplename, startfile, filerun))
       {
	 std::cerr << "Cannot get the tree " << std::endl;
       }
  }
  const int maxevent = std::atoi(Maxevent);  
  double Lumiscale = allSamples[subsamplename].getWeight(); 
  
  TString sampleString(subsamplename);
  if(sampleString.Contains("Data")){Lumiscale = 1.0; isData = true;}

  //Use BaselineVessel class for baseline variables and selections
  std::string spec = "ControlSample";
  std::string filterevent = "SingleMuon_csc2015.txt";
  ExpBaselineVessel = new BaselineVessel(spec, filterevent);
  AnaFunctions::prepareForNtupleReader();
  AnaFunctions::prepareTopTagger();
  NTupleReader *tr =0;
  
  if( isData ) tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames_DataOnly);
  else tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames);
  // tr = new NTupleReader(fChain);
    
  tr->registerFunction(&passBaselineFuncExp);
// Add cleanJets function
  stopFunctions::cjh.setMuonIso("mini");
  stopFunctions::cjh.setElecIso("mini");
  stopFunctions::cjh.setRemove(false);
  tr->registerFunction(&stopFunctions::cleanJets);
 
  TauResponse tauResp(respTempl);

  TRandom3 * rndm = new TRandom3(12345);

  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;

  int entries = tr->getNEntries();
  
  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<Lumiscale<<std::endl; 
  cout<<"maxevent: "<<maxevent<<endl;
  // Loop over the events (tree entries)
  int k = 0;
  while(tr->getNextEvent()){
    k++;
    if(maxevent>=0 && tr->getEvtNum() > maxevent ) break;
    // Add print out of the progress of looping
    if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;

    const vector<TLorentzVector> &muonsLVec = tr->getVec<TLorentzVector>("muonsLVec");
    const vector<double> &muonsRelIso = tr->getVec<double>("muonsRelIso");
    const vector<double> &muonsMiniIso = tr->getVec<double>("muonsMiniIso");
    const vector<double> &muonsMtw = tr->getVec<double>("muonsMtw");    
    const vector<TLorentzVector> &jetsLVec = tr->getVec<TLorentzVector>("jetsLVec");
    const vector<double> &recoJetsBtag_0 = tr->getVec<double>("recoJetsBtag_0");
    const vector<int> &looseisoTrksMatchedJetIdx = tr->getVec<int>("looseisoTrksMatchedJetIdx");
    const vector<TLorentzVector> &loose_isoTrksLVec = tr->getVec<TLorentzVector>("loose_isoTrksLVec");
    const vector<double> &loose_isoTrks_iso = tr->getVec<double>("loose_isoTrks_iso");
    const vector<double> &loose_isoTrks_mtw = tr->getVec<double>("loose_isoTrks_mtw");
    const vector<int> &loose_isoTrks_pdgId = tr->getVec<int>("loose_isoTrks_pdgId");
    const std::vector<int> & muonsFlagIDVec = tr->getVec<int>("muonsFlagMedium");
    const std::vector<double>& recoJetschargedHadronEnergyFraction = tr->getVec<double>("recoJetschargedHadronEnergyFraction");
    const std::vector<double>& recoJetschargedEmEnergyFraction = tr->getVec<double>("recoJetschargedEmEnergyFraction");
    const std::vector<std::string> & TriggerNames = tr->getVec<std::string>("TriggerNames");
    const std::vector<int> & PassTrigger = tr->getVec<int>("PassTrigger");

    double met=tr->getVar<double>("met");
    double metphi=tr->getVar<double>("metphi");
    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);
    int run = tr->getVar<int>("run");
    int lumi = tr->getVar<int>("lumi");
    int event = tr->getVar<int>("event");
    int nElectrons = tr->getVar<int>("nElectrons_CUT"+spec);
    int nMuons = tr->getVar<int>("nMuons_CUT"+spec);
    
    bool passBaseline_tru = tr->getVar<bool>("passBaseline"+spec);
    bool passMuonVeto_tru = tr->getVar<bool>("passMuonVeto"+spec);
    bool passEleVeto_tru = tr->getVar<bool>("passEleVeto"+spec);
    bool passIsoTrkVeto_tru = tr->getVar<bool>("passIsoTrkVeto"+spec);
    bool passLeptVeto_tru = tr->getVar<bool>("passLeptVeto"+spec);
    bool passnJets_tru = tr->getVar<bool>("passnJets"+spec);
    bool passdPhis_tru = tr->getVar<bool>("passdPhis"+spec);
    bool passMET_tru =  tr->getVar<bool>("passMET"+spec);
    bool passBJets_tru = tr->getVar<bool>("passBJets"+spec);
    bool passTagger_tru = tr->getVar<bool>("passTagger"+spec);
    bool passHT_tru =  tr->getVar<bool>("passHT"+spec);
    bool passMT2_tru = tr->getVar<bool>("passMT2" + spec);
    bool passNoiseEventFilter = tr->getVar<bool>("passNoiseEventFilter" + spec);

    bool passBaselineClosure = passMuonVeto_tru && passEleVeto_tru && passIsoTrkVeto_tru && passnJets_tru && passdPhis_tru && passMET_tru && passBJets_tru && passTagger_tru && passHT_tru;
      //Event Filter
    if(!passNoiseEventFilter) continue;  

    if(isData){
      bool foundTrigger = false;
      for(unsigned it=0; it<TriggerNames.size(); it++){
	if( sampleString.Contains("SingleMuon") ){
	  if( TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT350_v") != string::npos){
	    if( PassTrigger[it] ) foundTrigger = true;
	  }
	}
      }
      if( !foundTrigger ) continue;
    }
    //Control sample
    // The kinematic properties of the well-reconstructed, isolated muon                                                                    
    vector<TLorentzVector> isomuonsLVec;
    vector<int> isomuonsIdxVec;
    for(unsigned int im=0; im<muonsLVec.size(); im++){
      if( AnaFunctions::passMuon(muonsLVec.at(im), muonsMiniIso.at(im), muonsMtw.at(im), muonsFlagIDVec.at(im), AnaConsts::muonsMiniIsoArr) ){ isomuonsLVec.push_back(muonsLVec.at(im)); isomuonsIdxVec.push_back(im); }
    }
    
    // Require one and only one muon
    // Veto events with additional electrons (same veto criteria as baseline for electrons)
    if( nMuons == 1 && nElectrons == AnaConsts::nElectronsSel ) {
      if( nMuons != isomuonsLVec.size() ){ std::cout<<"ERROR ... mis-matching between veto muon and selected muon! Skipping..."<<std::endl; continue; }
      const TLorentzVector muLVec = isomuonsLVec.at(0);
      // Use only events where the muon is inside acceptance                                                                             
      if( muLVec.Pt() < TauResponse::ptMin() ) continue;
      if( fabs(muLVec.Eta()) > TauResponse::etaMax() ) continue;
      //mtW correction
      const double mtw = calcMT(muLVec, metLVec);
      bool pass_mtw = false;
      if(mtw<100)pass_mtw = true;
      
      //control sample distribution
      const double HT = tr->getVar<double>("HT"+spec);
      const int nJets30 = tr->getVar<int>("cntNJetsPt30Eta24"+spec);
      const int nJets50 = tr->getVar<int>("cntNJetsPt50Eta24"+spec);
      const double MT2 = tr->getVar<double>("best_had_brJet_MT2"+spec);
      if(HT < AnaConsts::defaultHTcut) continue;
      if(nJets30 < AnaConsts::nJetsSelPt30Eta24 - 1) continue;
      if(nJets50 < AnaConsts::nJetsSelPt50Eta24 - 1) continue;
      //      if(met < AnaConsts::defaultMETcut) continue;
      //            if(met < 100) continue;
      if(!isData){
	FillDouble(myBaseHistgram.hMCMuonPt, muLVec.Pt(), Lumiscale);
	FillDouble(myBaseHistgram.hMCMuonEta, muLVec.Eta(), Lumiscale);
	FillDouble(myBaseHistgram.hMCmet, met, Lumiscale);
	FillDouble(myBaseHistgram.hMCmetphi, metphi, Lumiscale);
	if(met >= 300) FillDouble(myBaseHistgram.hMCmetphi_met, metphi, Lumiscale);
	FillDouble(myBaseHistgram.hMCht, HT, Lumiscale);
	FillDouble(myBaseHistgram.hMCmt, mtw, Lumiscale);
	FillInt(myBaseHistgram.hMCnJet, nJets30, Lumiscale);
      }
      else{
	FillDouble(myBaseHistgram.hDataMuonPt, muLVec.Pt(), Lumiscale);
	FillDouble(myBaseHistgram.hDataMuonEta, muLVec.Eta(), Lumiscale);
        FillDouble(myBaseHistgram.hDatamet, met, Lumiscale);
        FillDouble(myBaseHistgram.hDatametphi, metphi, Lumiscale);
	if(met >= 300) FillDouble(myBaseHistgram.hDatametphi_met, metphi, Lumiscale);
        FillDouble(myBaseHistgram.hDataht, HT, Lumiscale);
	FillDouble(myBaseHistgram.hDatamt, mtw, Lumiscale);
	FillInt(myBaseHistgram.hDatanJet, nJets30, Lumiscale);
      }
      int jetsort = 0;
      for(unsigned int j = 0; j<jetsLVec.size(); j++){
	if(!AnaFunctions::jetPassCuts(jetsLVec[j], AnaConsts::pt30Eta24Arr)) continue;
	jetsort++;
	if(!isData){
	  FillDouble(myBaseHistgram.hMCJetPt, jetsLVec[j].Pt(), Lumiscale);
	  FillDouble(myBaseHistgram.hMCJetEta, jetsLVec[j].Eta(), Lumiscale);
	  FillDouble(myBaseHistgram.hMCJetPhi, jetsLVec[j].Phi(), Lumiscale);
	  FillDouble(myBaseHistgram.hMCJetMass, jetsLVec[j].M(), Lumiscale);
	  if(jetsort==1)FillDouble(myBaseHistgram.hMC1stJetMass, jetsLVec[j].M(), Lumiscale);
	  if(jetsort==2)FillDouble(myBaseHistgram.hMC2ndJetMass, jetsLVec[j].M(), Lumiscale);
	  if(jetsort==3)FillDouble(myBaseHistgram.hMC3rdJetMass, jetsLVec[j].M(), Lumiscale);
	}	
	else{
	  FillDouble(myBaseHistgram.hDataJetPt, jetsLVec[j].Pt(), Lumiscale);
	  FillDouble(myBaseHistgram.hDataJetEta, jetsLVec[j].Eta(), Lumiscale);
	  FillDouble(myBaseHistgram.hDataJetPhi, jetsLVec[j].Phi(), Lumiscale);
	  FillDouble(myBaseHistgram.hDataJetMass, jetsLVec[j].M(), Lumiscale);
	  if(jetsort==1)FillDouble(myBaseHistgram.hData1stJetMass, jetsLVec[j].M(), Lumiscale);
          if(jetsort==2)FillDouble(myBaseHistgram.hData2ndJetMass, jetsLVec[j].M(), Lumiscale);
          if(jetsort==3)FillDouble(myBaseHistgram.hData3rdJetMass, jetsLVec[j].M(), Lumiscale);
	}
      }
      
      // "Cross cleaning": find the jet that corresponds to the muon                                                                    
      const std::vector<TLorentzVector>& cleanJetVec      = tr->getVec<TLorentzVector>("cleanJetVec");
      const std::vector<double>& cleanJetBtag             = tr->getVec<double>("cleanJetBTag");
      // Get the cleaned jet indice (pointing to the jetsLVec) for the corresponding muons                                              
      const std::vector<int>& rejectJetIdx_formuVec = tr->getVec<int>("rejectJetIdx_formuVec");
      
    }//control sample loop
    
  }	//event loop
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();

  return 0;  
}


