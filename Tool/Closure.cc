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
#include "Closure.h"
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
  const string condorSpec = argc==7 ? argv[6]: "";

  if(!FillChain(fChain, subsamplename, condorSpec, startfile, filerun))
       {
	 std::cerr << "Cannot get the tree " << std::endl;
       }
  const int maxevent = std::atoi(Maxevent);  
  
  TString sampleString(subsamplename);
  if(sampleString.Contains("Data")){Lumiscale = 1.0; isData = true;}

  //Use BaselineVessel class for baseline variables and selections
  std::string spec = "ClosureExp";
  std::string filterevent = "SingleMuon_csc2015.txt";
  ExpBaselineVessel = new BaselineVessel(spec, filterevent);
  AnaFunctions::prepareForNtupleReader();
  AnaFunctions::prepareTopTagger();
  NTupleReader *tr =0;
  
  if( isData ) tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames_DataOnly);
  else tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames);
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
    double prntwgt = -1;
    double prntcor = -1;
    double wgtSB[45] = {0};
    for(int b=0; b<45; b++){
      wgtSB[b]= -1;
    }
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
    const std::vector<double> muonspfActivity = tr->getVec<double>("muonspfActivity");
    const std::vector<int> & muonsFlagIDVec = tr->getVec<int>("muonsFlagMedium");
    const std::vector<double>& recoJetschargedHadronEnergyFraction = tr->getVec<double>("recoJetschargedHadronEnergyFraction");
    const std::vector<double>& recoJetschargedEmEnergyFraction = tr->getVec<double>("recoJetschargedEmEnergyFraction");
    const std::vector<std::string> & TriggerNames = tr->getVec<std::string>("TriggerNames");
    const std::vector<int> & PassTrigger = tr->getVec<int>("PassTrigger");

    double met=tr->getVar<double>("met");
    double metphi=tr->getVar<double>("metphi");
    //const std::vector<double> &metMagUp = tr->getVec<double>("metMagUp");
    //const std::vector<double> &metMagDown = tr->getVec<double>("metMagDown");

    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);
    /*double metjecUp = metMagUp[1];//met Uncertainty
    double metjecLow = metMagDown[1];//met Uncertainty
    double metjerUp = metMagUp[0];//met Uncertainty
    double metjerLow = metMagDown[0];//met Uncertainty
    TLorentzVector metLVecjecUp;metLVecjecUp.SetPtEtaPhiM(metjecUp, 0, metphi, 0);
    TLorentzVector metLVecjecLow;metLVecjecLow.SetPtEtaPhiM(metjecLow, 0, metphi, 0);
    TLorentzVector metLVecjerUp;metLVecjerUp.SetPtEtaPhiM(metjerUp, 0, metphi, 0);
    TLorentzVector metLVecjerLow;metLVecjerLow.SetPtEtaPhiM(metjerLow, 0, metphi, 0);
    */
    int run = tr->getVar<int>("run");
    int lumi = tr->getVar<int>("lumi");
    int event = tr->getVar<int>("event");
    int nElectrons = tr->getVar<int>("nElectrons_CUT"+spec);
    int nMuons = tr->getVar<int>("nMuons_CUT"+spec);
    bool passNoiseEventFilter = tr->getVar<bool>("passNoiseEventFilter"+spec);
    const double EvtWt = tr->getVar<double>("evtWeight");
    //change event weight for MC sample
    EventWeight = EvtWt;
    // if(EventWeight<1) cout<<"EventWeight: "<<EventWeight<<endl;
    //Lumiscale = Lumiscale * EventWeight;
    //Event Filter                                                                                                                           
    if(!passNoiseEventFilter) continue;
    
    //Expectation part -- do it before prediction & do NOT skipping events

    if(!isData){      
      const vector<TLorentzVector> &genDecayLVec = tr->getVec<TLorentzVector>("genDecayLVec");
      const vector<int> &genDecayIdxVec =  tr->getVec<int>("genDecayIdxVec");
      const vector<int> &genDecayPdgIdVec = tr->getVec<int>("genDecayPdgIdVec");
      const vector<int> &W_emuVec =  tr->getVec<int>("W_emuVec");
      const vector<int> &W_tau_emuVec =  tr->getVec<int>("W_tau_emuVec");
      const vector<int> &W_tau_prongsVec =  tr->getVec<int>("W_tau_prongsVec");
      
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
      const int nJets_tru = tr->getVar<int>("cntNJetsPt30Eta24"+spec);
      const int nbJets_tru = tr->getVar<int>("cntCSVS"+spec);
      const int nTops_tru = tr->getVar<int>("nTopCandSortedCnt"+spec);
      const double MT2_tru = tr->getVar<double>("best_had_brJet_MT2"+spec);
      const vector<double> dPhiVec_tru = tr->getVec<double>("dPhiVec"+spec);
      const double ht_tru = tr->getVar<double>("HT"+spec);

      //mht calculation
      TLorentzVector Mht_truLVec;	      
	for(unsigned int ij=0; ij<jetsLVec.size(); ij++){
	  if( !AnaFunctions::jetPassCuts(jetsLVec[ij], AnaConsts::pt30Arr) ) continue;
	  Mht_truLVec -= jetsLVec[ij];
	}
	const double Mht_tru = Mht_truLVec.Pt();

      bool passBaselineClosure = passMuonVeto_tru && passEleVeto_tru && passIsoTrkVeto_tru && passnJets_tru && passdPhis_tru && passMET_tru && passBJets_tru && passTagger_tru && passHT_tru && passMT2_tru;
        
      //Select only events where the W decayed into a hadronically decaying tau      
      if(W_tau_prongsVec.size()!=0){
      //	myBaseHistgram.hTruecutFlow->Fill("original", Lumiscale);
	myBaseHistgram.hTruecutFlow->Fill("Filter", passNoiseEventFilter * Lumiscale);
	myBaseHistgram.hTruecutFlow->Fill("nJets", passNoiseEventFilter * passnJets_tru * Lumiscale);
	myBaseHistgram.hTruecutFlow->Fill("MuonVeto", passNoiseEventFilter * passnJets_tru * passMuonVeto_tru * Lumiscale);
	myBaseHistgram.hTruecutFlow->Fill("EleVeto", passNoiseEventFilter * passnJets_tru * passMuonVeto_tru * passEleVeto_tru * Lumiscale);
	myBaseHistgram.hTruecutFlow->Fill("IskVeto",passNoiseEventFilter * passnJets_tru * passMuonVeto_tru * passEleVeto_tru * passIsoTrkVeto_tru * Lumiscale);	
	myBaseHistgram.hTruecutFlow->Fill("dPhis", passNoiseEventFilter * passnJets_tru * passMuonVeto_tru * passEleVeto_tru * passIsoTrkVeto_tru *  passdPhis_tru * Lumiscale);
	myBaseHistgram.hTruecutFlow->Fill("BJets", passNoiseEventFilter * passnJets_tru * passMuonVeto_tru * passEleVeto_tru * passIsoTrkVeto_tru *  passdPhis_tru * passBJets_tru * Lumiscale);	
	myBaseHistgram.hTruecutFlow->Fill("MET", passNoiseEventFilter * passnJets_tru * passMuonVeto_tru * passEleVeto_tru * passIsoTrkVeto_tru *  passdPhis_tru * passBJets_tru * passMET_tru * Lumiscale);
	myBaseHistgram.hTruecutFlow->Fill("Tagger", passNoiseEventFilter * passnJets_tru * passMuonVeto_tru * passEleVeto_tru * passIsoTrkVeto_tru *  passdPhis_tru * passBJets_tru * passMET_tru * passTagger_tru * Lumiscale);
	myBaseHistgram.hTruecutFlow->Fill("MT2", passNoiseEventFilter * passnJets_tru * passMuonVeto_tru * passEleVeto_tru * passIsoTrkVeto_tru *  passdPhis_tru * passBJets_tru * passMET_tru * passTagger_tru * passMT2_tru* Lumiscale);
	myBaseHistgram.hTruecutFlow->Fill("HT", passNoiseEventFilter * passnJets_tru * passMuonVeto_tru * passEleVeto_tru * passIsoTrkVeto_tru *  passdPhis_tru * passBJets_tru * passMET_tru * passTagger_tru * passMT2_tru*  passHT_tru * Lumiscale);	

	//Dihad fraction calculation to calculate lostlepton contribution
	if(sampleString.Contains("DiLep")){
	  int gentau(0);
	  for(unsigned ig=0; ig<genDecayLVec.size(); ig++){
	    int pdgId = genDecayPdgIdVec.at(ig);
	    if(abs(pdgId)==15){ 
	      gentau++;
	    }
	  }
	  bool isdihadtau = false;
	  if(gentau==2 && W_tau_emuVec.size()==0) isdihadtau = true;
	  if(passBaselineClosure){
	    int kSR = find_Binning_Index(nbJets_tru, nTops_tru, MT2_tru, met);
	    myBaseHistgram.hnodihadtau->Fill(kSR, Lumiscale);
	    if(isdihadtau)myBaseHistgram.hdihadtau->Fill(kSR, Lumiscale);
	  }
	}
      }
      //Exp Dist.
      if(W_tau_prongsVec.size() !=0 && passBaselineClosure && passNoiseEventFilter){
	int jSR = find_Binning_Index(nbJets_tru, nTops_tru, MT2_tru, met);
	if( jSR!= -1 ) {
	  myBaseHistgram.hTrueYields->Fill(jSR, Lumiscale);
	}
	FillDouble(myBaseHistgram.hTrueHt, ht_tru, Lumiscale);
	FillDouble(myBaseHistgram.hTruemet, met, Lumiscale);
	FillDouble(myBaseHistgram.hTruemht, Mht_tru, Lumiscale);
	FillDouble(myBaseHistgram.hTrueMT2, MT2_tru, Lumiscale);
	FillInt(myBaseHistgram.hTrueNJets, nJets_tru, Lumiscale);
	FillInt(myBaseHistgram.hTrueNbJets, nbJets_tru, Lumiscale);
	FillInt(myBaseHistgram.hTrueNTops, nTops_tru, Lumiscale);	
	FillDouble(myBaseHistgram.hTruedPhi0, dPhiVec_tru[0], Lumiscale);
	FillDouble(myBaseHistgram.hTruedPhi1, dPhiVec_tru[1], Lumiscale);
	FillDouble(myBaseHistgram.hTruedPhi2, dPhiVec_tru[2], Lumiscale);
      }
    }//end of expectation

    //Prediction part	
    if(isData){
      bool foundTrigger = false;
      for(unsigned it=0; it<TriggerNames.size(); it++){
	if( sampleString.Contains("SingleMuon") ){
	  //	  if( TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT350_v") || TriggerNames[it].find("HLT_PFMET170_NoiseCleaned_v") || TriggerNames[it].find("HLT_PFMET170_JetIdCleaned_v") || TriggerNames[it].find("HLT_Mu50_v")){
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
     
      myBaseHistgram.hPredcutFlow->Fill("CS", Lumiscale);
      
      //mtW correction
      const double mtw = calcMT(muLVec, metLVec);
      bool pass_mtw = false;
      if(mtw<100)pass_mtw = true;

      /*const double mtwjecUp = calcMT(muLVec, metLVecjecUp);
      const double mtwjecLow = calcMT(muLVec, metLVecjecLow);
      const double mtwjerUp = calcMT(muLVec, metLVecjerUp);
      const double mtwjerLow = calcMT(muLVec, metLVecjerLow);
      */
      bool istaumu = false, istaumu_genRecoMatch = false;
      if(!isData){
	const vector<TLorentzVector> &genDecayLVec1 = tr->getVec<TLorentzVector>("genDecayLVec");
	const vector<int> &genDecayPdgIdVec1 = tr->getVec<int>("genDecayPdgIdVec");
	const vector<int> &W_tau_emuVec1 =  tr->getVec<int>("W_tau_emuVec");
	// Find events that contain W->tau->mu
	// Note that any code using gen info should be checked if they work for data or not!
	for(unsigned int ig=0; ig<W_tau_emuVec1.size(); ig++){
	  int genIdx = W_tau_emuVec1.at(ig);
	  if( std::abs(genDecayPdgIdVec1.at(genIdx)) == 13 ){
	    istaumu = true;
	    const TLorentzVector & genLVec = genDecayLVec1.at(genIdx);
	    if( muLVec.DeltaR(genLVec) < 0.2 ) istaumu_genRecoMatch = true;
	  }
	} 
      }
      
      // "Cross cleaning": find the jet that corresponds to the muon                                                                    
      const std::vector<TLorentzVector>& cleanJetVec      = tr->getVec<TLorentzVector>("cleanJetVec");
      const std::vector<double>& cleanJetBtag             = tr->getVec<double>("cleanJetBTag");
      // Get the cleaned jet indice (pointing to the jetsLVec) for the corresponding muons                                              
      const std::vector<int>& rejectJetIdx_formuVec = tr->getVec<int>("rejectJetIdx_formuVec");
      
      //Implement IsoTrackVeto
      if(looseisoTrksMatchedJetIdx.size()!=loose_isoTrksLVec.size())cout<<"Error: isotrack vetor size mismatch"<<endl;
      int nisoTotal = 0, nisomu = 0;	    
      for(unsigned int it=0; it< looseisoTrksMatchedJetIdx.size();it++){
	if(!passIsoTrks1(loose_isoTrksLVec[it], loose_isoTrks_iso[it], loose_isoTrks_mtw[it], loose_isoTrks_pdgId[it])) continue;
	nisoTotal++;
	// Do the matching
	const double dr = muLVec.DeltaR(loose_isoTrksLVec[it]);
	if(dr<0.1) nisomu++; // muon tagged as isotrack
      }//finish isotrack loop
      if((nisoTotal-nisomu)!=0)continue;//isotrackveto on the remaing part (exclude mu)  
      
      // Force the mass to be 0 for met and mht
      TLorentzVector selmetLVec; selmetLVec.SetVectM( (metLVec+ muLVec).Vect(), 0 );
      
      int selNJetPt30Eta24 = AnaFunctions::countJets(cleanJetVec, AnaConsts::pt30Eta24Arr);
      int selNJetPt50Eta24 = AnaFunctions::countJets(cleanJetVec, AnaConsts::pt50Eta24Arr);
      
      // rejecting events with nJets less than requirements even adding one more tau jet
      if( selNJetPt30Eta24 < AnaConsts::nJetsSelPt30Eta24 - 1 ) continue;
      if( selNJetPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 - 1 ) continue;
      myBaseHistgram.hPredcutFlow->Fill("CS njet-1 cut", Lumiscale);

      // Get random number from tau-response template
      // The template is chosen according to the muon pt 
      TH1F* temp = (TH1F*)tauResp.Resp(muLVec.Pt());
      //Loop over template bin  
      for(int ib = 1; ib<=50; ib++){
	//const double scale = tauResp.getRandom(muLVec.Pt());
	const double scale = temp->GetBinCenter(ib);
	const double weight = temp->GetBinContent(ib) * temp->GetBinWidth(ib);

	// Scale muon pt and energy with tau response --> simulate tau jet pt and energy
	const double simTauJetPt = scale * muLVec.Pt();
	const double simTauJetE = scale * muLVec.E();
	const double simTauJetEta = muLVec.Eta();
	const double simTauJetPhi = muLVec.Phi();
	
	TLorentzVector ori_tauJetLVec; ori_tauJetLVec.SetPtEtaPhiE(simTauJetPt, simTauJetEta, simTauJetPhi, simTauJetE);
	TLorentzVector tauJetLVec = ori_tauJetLVec;
	int muJetIdx = rejectJetIdx_formuVec.at(isomuonsIdxVec.at(0));
	
	// Decide the CSV value for the tau jet -> use the CSV of the associated muon-jet as the tau jet CSV
	// Default set to be 0 (low enough to be NOT a b jet)
	double oriJetCSVS = 0;
	if( muJetIdx >= 0) oriJetCSVS = recoJetsBtag_0[muJetIdx];
	double mistag = Efficiency::mistag(Efficiency::Ptbin1(simTauJetPt));
	double rno = rndm->Rndm();
	if( rno < mistag) oriJetCSVS = 1.0;
	//Adjustment of tau jet to the remaining part of mu cleaned jet    
	if( muJetIdx >=0 ) tauJetLVec += cleanJetVec[muJetIdx];
	vector<TLorentzVector> combNJetVec;
	vector<double> combJetsBtag;
	bool includeTauJet = false;
	for(unsigned int ij=0; ij<cleanJetVec.size(); ij++){
	  if( ij == muJetIdx ) continue;
	  if( tauJetLVec.Pt() > cleanJetVec.at(ij).Pt() && !includeTauJet ){
	    combNJetVec.push_back(tauJetLVec); combJetsBtag.push_back(oriJetCSVS);
	    includeTauJet = true;
	  }
	  combNJetVec.push_back(cleanJetVec.at(ij)); combJetsBtag.push_back(cleanJetBtag.at(ij));
	}
	// it's possible that the tau jet is the least energetic jet so that it's not added into the combNJetVec during the loop
	if( !includeTauJet ){ combNJetVec.push_back(tauJetLVec); combJetsBtag.push_back(oriJetCSVS); }
	
	// Taking into account the simulated tau jet, recompute HT, MHT, and N(jets)
	// If simulted tau-jet meets same criteria as as HT jets, recompute HT and MHT
	double simHt = AnaFunctions::calcHT(combNJetVec, AnaConsts::pt30Eta24Arr);
	
	bool passhtpred = true;
	if( simHt < AnaConsts::defaultHTcut )passhtpred = false;
	
	TLorentzVector simMhtLVec;	      
	for(unsigned int ij=0; ij<combNJetVec.size(); ij++){
	  if( !AnaFunctions::jetPassCuts(combNJetVec[ij], AnaConsts::pt30Arr) ) continue;
	  simMhtLVec -= combNJetVec[ij];
	}
	
	//recompute met                                                                                                                          
	TLorentzVector simmetLVec; simmetLVec.SetVectM( (selmetLVec - ori_tauJetLVec).Vect(), 0);
	const double simMht = simMhtLVec.Pt();
	const double simmet = simmetLVec.Pt();
	const double simmetPhi = simmetLVec.Phi();
	
	bool passmetPred = true;
	if(simmet<AnaConsts::defaultMETcut)passmetPred = false;
	//recompute jetVec
	int combNJetPt30Eta24 = AnaFunctions::countJets(combNJetVec, AnaConsts::pt30Eta24Arr);
	int combNJetPt50Eta24 = AnaFunctions::countJets(combNJetVec, AnaConsts::pt50Eta24Arr);
	      
	bool passNJetPred = true;
	if(combNJetPt50Eta24<AnaConsts::nJetsSelPt50Eta24)passNJetPred = false;
	if(combNJetPt30Eta24<AnaConsts::nJetsSelPt30Eta24)passNJetPred = false;
	
	//recompute deltaphi
	std::vector<double> deltaPhiVec = AnaFunctions::calcDPhi(combNJetVec, simmetPhi, 3, AnaConsts::dphiArr);
	double dPhi0_pred = deltaPhiVec.at(0);
	double dPhi1_pred = deltaPhiVec.at(1);
	double dPhi2_pred = deltaPhiVec.at(2);
	
	bool passdeltaPhi = true;
	if( deltaPhiVec.at(0) < AnaConsts::dPhi0_CUT || deltaPhiVec.at(1) < AnaConsts::dPhi1_CUT || deltaPhiVec.at(2) < AnaConsts::dPhi2_CUT){
	  passdeltaPhi = false;
	}
	
	//recompute bjet
	
	int cnt1CSVS = AnaFunctions::countCSVS(combNJetVec, combJetsBtag, AnaConsts::cutCSVS, AnaConsts::bTagArr);
	bool passbJets = true;
	if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cnt1CSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cnt1CSVS < AnaConsts::high_nJetsSelBtagged ) ) ){
	  passbJets = false;
	}
	
	//top tagger input
	int comb30_pre = AnaFunctions::countJets(combNJetVec, AnaConsts::pt30Arr);
	std::vector<TLorentzVector> jetsLVec_forTagger_pre; std::vector<double> recoJetsBtag_forTagger_pre;
	AnaFunctions::prepareJetsForTagger(combNJetVec, combJetsBtag, jetsLVec_forTagger_pre, recoJetsBtag_forTagger_pre);
	int nTopCandSortedCnt_pre = -1;
	double MT2_pre = -1;
	double mTcomb_pre = -1;
	
	myBaseHistgram.hPredcutFlow->Fill("ht cut", passhtpred * Lumiscale * weight);
	bool passNjet6 = true;
	bool passNjet9 = true;
	if(combNJetPt30Eta24 < 6) passNjet6 = false;
	if(combNJetPt30Eta24 < 9) passNjet9 = false;

	//Apply baseline cut
	if(!passNJetPred)continue;
	myBaseHistgram.hPredcutFlow->Fill("njet cut", passhtpred * Lumiscale * weight);
	if(!passdeltaPhi) continue;
	myBaseHistgram.hPredcutFlow->Fill("dPhi cut",passhtpred * Lumiscale * weight);
	FillDouble(myBaseHistgram.hCheckmet_dPhi, simmet, passhtpred * Lumiscale * weight);
	if(!passmetPred)continue;
	myBaseHistgram.hPredcutFlow->Fill("met cut",passhtpred * Lumiscale * weight);
	FillDouble(myBaseHistgram.hCheckmet_met, simmet, passhtpred * Lumiscale * weight);
	FillInt(myBaseHistgram.hChecknjet_met, combNJetPt30Eta24, passhtpred * Lumiscale * weight);
	for(int ij = 0; ij < combNJetVec.size(); ij++){
	  if(!AnaFunctions::jetPassCuts(combNJetVec[ij], AnaConsts::pt30Eta24Arr)) continue;
	  FillDouble(myBaseHistgram.hCheckjetpt_met, combNJetVec[ij].Pt(), passhtpred * Lumiscale * weight);
	  FillDouble(myBaseHistgram.hCheckjetmass_met, combNJetVec[ij].M(), passhtpred * Lumiscale * weight);
	  if(combJetsBtag[ij] > AnaConsts::cutCSVS) FillDouble(myBaseHistgram.hCheckBjetpt_met, combNJetVec[ij].Pt(), passhtpred * Lumiscale * weight);
	}
	if(!passbJets) continue;
	myBaseHistgram.hPredcutFlow->Fill("bjet cut", passhtpred * Lumiscale * weight);
	FillDouble(myBaseHistgram.hCheckmet_bjet, simmet, passhtpred * Lumiscale * weight);
	FillInt(myBaseHistgram.hChecknjet_bjet, combNJetPt30Eta24, passhtpred * Lumiscale * weight);
        for(int ij = 0; ij < combNJetVec.size(); ij++){
          if(!AnaFunctions::jetPassCuts(combNJetVec[ij], AnaConsts::pt30Eta24Arr)) continue;
          FillDouble(myBaseHistgram.hCheckjetpt_bjet, combNJetVec[ij].Pt(), passhtpred * Lumiscale * weight);
          FillDouble(myBaseHistgram.hCheckjetmass_bjet, combNJetVec[ij].M(), passhtpred * Lumiscale * weight);
          if(combJetsBtag[ij] > AnaConsts::cutCSVS) FillDouble(myBaseHistgram.hCheckBjetpt_bjet, combNJetVec[ij].Pt(), passhtpred * Lumiscale * weight);
        }
	myBaseHistgram.hPredcutFlow->Fill("6 njet cut", passNjet6 * passhtpred * Lumiscale * weight);
	myBaseHistgram.hPredcutFlow->Fill("9 njet cut", passNjet9 * passNjet6 * passhtpred * Lumiscale * weight);

	//Apply TopTagger        
	if(comb30_pre >= AnaConsts::nJetsSel ){
	  type3Ptr->processEvent(jetsLVec_forTagger_pre, recoJetsBtag_forTagger_pre, simmetLVec);
	  nTopCandSortedCnt_pre = type3Ptr->nTopCandSortedCnt;
	  MT2_pre = type3Ptr->best_had_brJet_MT2;
	  mTcomb_pre = type3Ptr->best_had_brJet_mTcomb;
	}
	bool passMT2pred = true;
	if(MT2_pre < AnaConsts::defaultMT2cut)passMT2pred = false;
	bool passTopTagger = type3Ptr->passNewTaggerReq() && nTopCandSortedCnt_pre >= AnaConsts::low_nTopCandSortedSel;
	if(!passTopTagger) continue;

	myBaseHistgram.hPredcutFlow->Fill("ntop cut", passhtpred * Lumiscale * weight);	
        FillDouble(myBaseHistgram.hCheckmet_ntop, simmet, passhtpred * Lumiscale * weight);
        FillInt(myBaseHistgram.hChecknjet_ntop, combNJetPt30Eta24, passhtpred * Lumiscale * weight);
	for(int ij = 0; ij < combNJetVec.size(); ij++){
          if(!AnaFunctions::jetPassCuts(combNJetVec[ij], AnaConsts::pt30Eta24Arr)) continue;
          FillDouble(myBaseHistgram.hCheckjetpt_ntop, combNJetVec[ij].Pt(), passhtpred * Lumiscale * weight);
          FillDouble(myBaseHistgram.hCheckjetmass_ntop, combNJetVec[ij].M(), passhtpred * Lumiscale * weight);
          if(combJetsBtag[ij] > AnaConsts::cutCSVS) FillDouble(myBaseHistgram.hCheckBjetpt_ntop, combNJetVec[ij].Pt(), passhtpred * Lumiscale * weight);
        }

	//if(!passMT2pred) continue;  	
	//Activity variable calculation:
	//double muact = AnaFunctions::getMuonActivity(muLVec, jetsLVec, recoJetschargedHadronEnergyFraction, recoJetschargedEmEnergyFraction, AnaConsts::muonsAct);
	const double muact = muonspfActivity[isomuonsIdxVec.at(0)];
	// iSR: this should be determined by search region requirement
	const int kSR = find_Binning_Index(cnt1CSVS, nTopCandSortedCnt_pre, MT2_pre, simmet);
	
	//correction factor:                                                                                                                      
	const double corrBRWToTauHad = 0.65;  // Correction for the BR of hadronic tau decays                                                   
	const double corrBRTauToMu = 1-Efficiency::taumucorMix(Efficiency::Njetbin(combNJetPt30Eta24), Efficiency::metbin(simmet));//correction from tauonic mu contamination   
	const double corrMuRecoEff = 1./Efficiency::reco(Efficiency::Ptbin(muLVec.Pt()), Efficiency::Actbin(muact)); // Correction for muon reconstruction efficiency                                        
	const double corrMuIsoEff = 1./Efficiency::iso(Efficiency::Ptbin(muLVec.Pt()), Efficiency::Actbin(muact)); // Correction for muon isolation efficiency
	const double corrMuAcc = 1./Efficiency::SBaccMix(kSR); // Correction for muon acceptance
	const double corrmtWEff =  1./Efficiency::mtwMix(Efficiency::Njetbin(combNJetPt30Eta24), Efficiency::metbin(simmet)); //correction for mtW cut
	const double corrisotrkEff = 1- Efficiency::isotrkeffMix_NjetNbjet(Efficiency::Njetbin(combNJetPt30Eta24), Efficiency::NBjetbin(cnt1CSVS));//correction for isotrackveto eff.
	const double corrllcont = 1-llcont;//correction for lost lepton overlap
	const double corrTrgeff = trgeff;//correction for trigger efficiency
	//const double corrllcont = 1;
	//const double corrTrgeff = 1;

	//The overall correction factor                                                                                                       
	const double corr = corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff * corrmtWEff * corrisotrkEff * corrBRTauToMu * corrllcont * corrTrgeff * weight;
	const double Evt_weight = Lumiscale * weight;
	const double Evt_corr = corr * Lumiscale;

	myBaseHistgram.hPredcutFlow->Fill("mtw cut", passhtpred * pass_mtw * Lumiscale  * weight);
	myBaseHistgram.hPredcutFlow->Fill("With corr", passhtpred * pass_mtw * Lumiscale  * corr);		
	myBaseHistgram.hPredcutFlow->Fill("With corr and MT2", passhtpred * pass_mtw *  passMT2pred * Lumiscale  * corr);		

	//Correction factor plot
	bool frstadd = true;
	bool frstadd1 = true;
	bool frstadd2 = true;
	if(pass_mtw && passhtpred && passMT2pred) {
	  myBaseHistgram.hcorrection->Fill(corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff * corrmtWEff * corrisotrkEff * corrBRTauToMu);
	  myBaseHistgram.hSB->Fill(kSR, Evt_corr);
	  myBaseHistgram.hweight_SB->Fill(kSR, Evt_corr);
	  if(prntwgt==-1){
	    prntwgt = weight;
	    frstadd = false;
	  }
          if(prntcor==-1){
            prntcor = corr;
            frstadd1 = false;
          }
	  // for(int b=0;b<45;b++){
	    if(wgtSB[kSR]==-1){
	      wgtSB[kSR]=corr;
	      frstadd2=false;
	    }
	    if(frstadd2)wgtSB[kSR] = wgtSB[kSR] + corr;
	    //}
	  if(frstadd) prntwgt = prntwgt + weight;	
	  if(frstadd1)prntcor = prntcor + corr;
	}
	// Fill the prediction dist
	if( pass_mtw && passhtpred && passMT2pred){
	  FillDouble(myBaseHistgram.hPredHt, simHt, Evt_corr);
	  FillDouble(myBaseHistgram.hPredmet_wt, simmet, Evt_corr);
	  FillDouble(myBaseHistgram.hPredmht_wt, simMht, Evt_corr);
	  FillDouble(myBaseHistgram.hPredMT2_wt, MT2_pre, Evt_corr);
	  FillDouble(myBaseHistgram.hPredmTcomb, mTcomb_pre, Evt_corr);
	  FillInt(myBaseHistgram.hPredNJets, combNJetPt30Eta24, Evt_corr);
	  FillInt(myBaseHistgram.hPredNbJets_wt, cnt1CSVS, Evt_corr);
	  FillInt(myBaseHistgram.hPredNTops_wt, nTopCandSortedCnt_pre, Evt_corr);
	  FillDouble(myBaseHistgram.hPreddPhi0_wt, dPhi0_pred, Evt_corr);
	  FillDouble(myBaseHistgram.hPreddPhi1_wt, dPhi1_pred, Evt_corr);
	  FillDouble(myBaseHistgram.hPreddPhi2_wt, dPhi2_pred, Evt_corr);
	  //corelation hist.
	  FillDouble(myBaseHistgram.h1DMET_wt, simmet, Evt_corr);
	  FillDouble(myBaseHistgram.h1DMT2_wt, simmet, Evt_corr);
	  FillInt(myBaseHistgram.h1DNb_wt, cnt1CSVS, Evt_corr);
	  FillInt(myBaseHistgram.h1DNt_wt, nTopCandSortedCnt_pre, Evt_corr);
	}
	
	//Fill search bin prediction
	if( pass_mtw && passhtpred && passMT2pred){
	  if( kSR!=-1) {
	    //correction in each SB  
	    const double SBcorr = corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff * corrmtWEff * corrisotrkEff * corrBRTauToMu * corrllcont * corrTrgeff * weight;
	    const double Evt_SBcorr = SBcorr * Lumiscale;
	    myBaseHistgram.hPredYields_wt->Fill(kSR, Evt_SBcorr);
	    myBaseHistgram.h1DSB_wt->Fill(kSR, Evt_SBcorr);//for corelation

	  }
	}

	//mtW correction
	if(!isData){ 
	  if( !istaumu_genRecoMatch && passhtpred && passMT2pred){
	    FillInt(myBaseHistgram.hnomtW_Njet_wt,combNJetPt30Eta24,Evt_weight);
	    FillDouble(myBaseHistgram.hnomtW_Ht_wt,simHt,Evt_weight);
	    FillDouble(myBaseHistgram.hnomtW_met_wt,simmet,Evt_weight);
	    FillDouble(myBaseHistgram.hnomtW_MT2_wt,MT2_pre,Evt_weight);
	    FillInt(myBaseHistgram.hnomtW_Nbjet_wt,cnt1CSVS,Evt_weight);
	    FillInt(myBaseHistgram.hnomtW_Ntop_wt,nTopCandSortedCnt_pre,Evt_weight);
	    Fill2D(myBaseHistgram.hnomtW_Njetmet_wt, combNJetPt30Eta24, simmet, Evt_weight);
	    Fill2D(myBaseHistgram.hnomtWSys_Njetmet, combNJetPt30Eta24, simmet, Evt_weight);
	    if(pass_mtw){
	      FillInt(myBaseHistgram.hmtW_Njet_wt,combNJetPt30Eta24,Evt_weight);
	      FillDouble(myBaseHistgram.hmtW_Ht_wt,simHt,Evt_weight);
	      FillDouble(myBaseHistgram.hmtW_met_wt,simmet,Evt_weight);
	      FillDouble(myBaseHistgram.hmtW_MT2_wt,MT2_pre,Evt_weight);
	      FillInt(myBaseHistgram.hmtW_Nbjet_wt,cnt1CSVS,Evt_weight);
	      FillInt(myBaseHistgram.hmtW_Ntop_wt,nTopCandSortedCnt_pre,Evt_weight);
	      Fill2D(myBaseHistgram.hmtW_Njetmet_wt, combNJetPt30Eta24, simmet, Evt_weight);
	    }
	    /*    if(mtwjecUp<100)Fill2D(myBaseHistgram.hmtWSysjecUp_Njetmet, combNJetPt30Eta24, simmet, Evt_weight);
	    if(mtwjecLow<100)Fill2D(myBaseHistgram.hmtWSysjecLow_Njetmet, combNJetPt30Eta24, simmet, Evt_weight);
	    if(mtwjerUp<100)Fill2D(myBaseHistgram.hmtWSysjerUp_Njetmet, combNJetPt30Eta24, simmet, Evt_weight);
	    if(mtwjerLow<100)Fill2D(myBaseHistgram.hmtWSysjerLow_Njetmet, combNJetPt30Eta24, simmet, Evt_weight);
	    */
	    if(kSR!=-1){
	      myBaseHistgram.hnomtW_wt->Fill(kSR, Evt_weight);
	      if(pass_mtw) myBaseHistgram.hmtW_wt->Fill(kSR, Evt_weight);
	    }
	  }
	  const unsigned int BjetBin = Efficiency::NBjetbin(cnt1CSVS);
	  const unsigned int HTBin = Efficiency::Htbin(simHt);
	  //tau mu contamination
	  if( passhtpred && passMT2pred){
            FillInt(myBaseHistgram.hnotaumu_Njet_wt,combNJetPt30Eta24,Evt_weight);
            FillDouble(myBaseHistgram.hnotaumu_Ht_wt,simHt,Evt_weight);
            FillDouble(myBaseHistgram.hnotaumu_met_wt,simmet,Evt_weight);
            FillDouble(myBaseHistgram.hnotaumu_MT2_wt,MT2_pre,Evt_weight);
            FillInt(myBaseHistgram.hnotaumu_Nbjet_wt,cnt1CSVS,Evt_weight);
            FillInt(myBaseHistgram.hnotaumu_Ntop_wt,nTopCandSortedCnt_pre,Evt_weight);
	    Fill2D(myBaseHistgram.hnotaumu_Njetmet_wt, combNJetPt30Eta24, simmet, Evt_weight);
	    Fill2D(myBaseHistgram.hnotaumu_Njetmet_Bjet.at(BjetBin), combNJetPt30Eta24, simmet, Evt_weight);
	    Fill2D(myBaseHistgram.hnotaumu_Njetmet_Ht.at(HTBin), combNJetPt30Eta24, simmet, Evt_weight);
	    if(istaumu_genRecoMatch) {
	      FillInt(myBaseHistgram.htaumu_Njet_wt,combNJetPt30Eta24,Evt_weight);
	      FillDouble(myBaseHistgram.htaumu_Ht_wt,simHt,Evt_weight);
	      FillDouble(myBaseHistgram.htaumu_met_wt,simmet,Evt_weight);
	      FillDouble(myBaseHistgram.htaumu_MT2_wt,MT2_pre,Evt_weight);
	      FillInt(myBaseHistgram.htaumu_Nbjet_wt,cnt1CSVS,Evt_weight);
	      FillInt(myBaseHistgram.htaumu_Ntop_wt,nTopCandSortedCnt_pre,Evt_weight);
	      Fill2D(myBaseHistgram.htaumu_Njetmet_wt, combNJetPt30Eta24, simmet, Evt_weight);
	      Fill2D(myBaseHistgram.htaumu_Njetmet_Bjet.at(BjetBin), combNJetPt30Eta24, simmet, Evt_weight);
	      Fill2D(myBaseHistgram.htaumu_Njetmet_Ht.at(HTBin), combNJetPt30Eta24, simmet, Evt_weight);
	    }	    
	    if(kSR!=-1){
	      myBaseHistgram.hnotaumu_wt->Fill(kSR, Evt_weight);
	      if(istaumu_genRecoMatch) myBaseHistgram.htaumu_wt->Fill(kSR, Evt_weight * corrBRTauToMu);
	    }
	  }
	}
	
      }//template bin loop
    }//control sample loop

    //Fill corelation histo
    TauResponse::HistfillCorr(myBaseHistgram.h1DSB_wt, myBaseHistgram.h2DSB_corr);
    TauResponse::HistfillCorr(myBaseHistgram.h1DMET_wt, myBaseHistgram.h2DMET_corr);
    TauResponse::HistfillCorr(myBaseHistgram.h1DMT2_wt, myBaseHistgram.h2DMT2_corr);
    TauResponse::HistfillCorr(myBaseHistgram.h1DNb_wt, myBaseHistgram.h2DNb_corr);
    TauResponse::HistfillCorr(myBaseHistgram.h1DNt_wt, myBaseHistgram.h2DNt_corr);

    //correct the uncertainties in pred histo
    TauResponse::Histfill(myBaseHistgram.hPredYields_wt, myBaseHistgram.hPredYields);
    TauResponse::Histfill(myBaseHistgram.hnomtW_wt, myBaseHistgram.hnomtW);    
    TauResponse::Histfill(myBaseHistgram.hmtW_wt, myBaseHistgram.hmtW);    
    TauResponse::Histfill(myBaseHistgram.hnotaumu_wt, myBaseHistgram.hnotaumu);
    TauResponse::Histfill(myBaseHistgram.htaumu_wt, myBaseHistgram.htaumu);
    TauResponse::Histfill(myBaseHistgram.hPredmet_wt, myBaseHistgram.hPredmet);
    TauResponse::Histfill(myBaseHistgram.hPredmht_wt, myBaseHistgram.hPredmht);
    TauResponse::Histfill(myBaseHistgram.hPredMT2_wt, myBaseHistgram.hPredMT2);
    TauResponse::Histfill(myBaseHistgram.hPredNbJets_wt, myBaseHistgram.hPredNbJets);
    TauResponse::Histfill(myBaseHistgram.hPredNTops_wt, myBaseHistgram.hPredNTops);
    TauResponse::Histfill(myBaseHistgram.hPreddPhi0_wt, myBaseHistgram.hPreddPhi0);  
    TauResponse::Histfill(myBaseHistgram.hPreddPhi1_wt, myBaseHistgram.hPreddPhi1);  
    TauResponse::Histfill(myBaseHistgram.hPreddPhi2_wt, myBaseHistgram.hPreddPhi2);  
    TauResponse::Histfill(myBaseHistgram.hnomtW_Njet_wt, myBaseHistgram.hnomtW_Njet);
    TauResponse::Histfill(myBaseHistgram.hnomtW_Ht_wt, myBaseHistgram.hnomtW_Ht);
    TauResponse::Histfill(myBaseHistgram.hnomtW_met_wt, myBaseHistgram.hnomtW_met);
    TauResponse::Histfill(myBaseHistgram.hnomtW_MT2_wt, myBaseHistgram.hnomtW_MT2);
    TauResponse::Histfill(myBaseHistgram.hnomtW_Nbjet_wt, myBaseHistgram.hnomtW_Nbjet);
    TauResponse::Histfill(myBaseHistgram.hnomtW_Ntop_wt, myBaseHistgram.hnomtW_Ntop);
    TauResponse::Histfill(myBaseHistgram.hmtW_Njet_wt, myBaseHistgram.hmtW_Njet);
    TauResponse::Histfill(myBaseHistgram.hmtW_Ht_wt, myBaseHistgram.hmtW_Ht);
    TauResponse::Histfill(myBaseHistgram.hmtW_met_wt, myBaseHistgram.hmtW_met);
    TauResponse::Histfill(myBaseHistgram.hmtW_MT2_wt, myBaseHistgram.hmtW_MT2);
    TauResponse::Histfill(myBaseHistgram.hmtW_Nbjet_wt, myBaseHistgram.hmtW_Nbjet);
    TauResponse::Histfill(myBaseHistgram.hmtW_Ntop_wt, myBaseHistgram.hmtW_Ntop);
    TauResponse::Histfill(myBaseHistgram.hnotaumu_Njet_wt, myBaseHistgram.hnotaumu_Njet);
    TauResponse::Histfill(myBaseHistgram.hnotaumu_Ht_wt, myBaseHistgram.hnotaumu_Ht);
    TauResponse::Histfill(myBaseHistgram.hnotaumu_met_wt, myBaseHistgram.hnotaumu_met);
    TauResponse::Histfill(myBaseHistgram.hnotaumu_MT2_wt, myBaseHistgram.hnotaumu_MT2);
    TauResponse::Histfill(myBaseHistgram.hnotaumu_Nbjet_wt, myBaseHistgram.hnotaumu_Nbjet);
    TauResponse::Histfill(myBaseHistgram.hnotaumu_Ntop_wt, myBaseHistgram.hnotaumu_Ntop);
    TauResponse::Histfill(myBaseHistgram.htaumu_Njet_wt, myBaseHistgram.htaumu_Njet);
    TauResponse::Histfill(myBaseHistgram.htaumu_Ht_wt, myBaseHistgram.htaumu_Ht);
    TauResponse::Histfill(myBaseHistgram.htaumu_met_wt, myBaseHistgram.htaumu_met);
    TauResponse::Histfill(myBaseHistgram.htaumu_MT2_wt, myBaseHistgram.htaumu_MT2);
    TauResponse::Histfill(myBaseHistgram.htaumu_Nbjet_wt, myBaseHistgram.htaumu_Nbjet);
    TauResponse::Histfill(myBaseHistgram.htaumu_Ntop_wt, myBaseHistgram.htaumu_Ntop);
    /*    TauResponse::Histfill2D(myBaseHistgram.hnomtW_Njetmet_wt, myBaseHistgram.hnomtW_Njetmet);
    TauResponse::Histfill2D(myBaseHistgram.hmtW_Njetmet_wt, myBaseHistgram.hmtW_Njetmet);
    TauResponse::Histfill2D(myBaseHistgram.hnotaumu_Njetmet_wt, myBaseHistgram.hnotaumu_Njetmet);
    TauResponse::Histfill2D(myBaseHistgram.htaumu_Njetmet_wt, myBaseHistgram.htaumu_Njetmet);
    */
    TauResponse::Histfill1Dto2D(myBaseHistgram.hSB, myBaseHistgram.hweight_SBCpy);
    myBaseHistgram.htempweight->Fill(prntwgt, Lumiscale);
    myBaseHistgram.htotweight->Fill(prntcor, Lumiscale);
    for(int b=0; b<45; b++){
      myBaseHistgram.htotweight_SB->Fill(b, wgtSB[b], Lumiscale);
    }
    myBaseHistgram.hevtWt->Fill(EventWeight, Lumiscale);
  }	//event loop
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();
  
  // This print out can be used to extract the corrBRTauToMu ratio
  
  return 0;
  
}
bool passIsoTrks1(const TLorentzVector isoTrksLVec, const double isoTrksIso, const double isoTrksMtw, const int isoTrkspdgId){
  bool passIsoTrks = false;
  if( std::abs(isoTrkspdgId) == 11 || std::abs(isoTrkspdgId) == 13 ){
    if( AnaFunctions::passIsoTrk(isoTrksLVec, isoTrksIso, isoTrksMtw, AnaConsts::isoLepTrksArr ) ) passIsoTrks = true;
  }
  if( std::abs(isoTrkspdgId) == 211 ){
    if(AnaFunctions::passIsoTrk(isoTrksLVec, isoTrksIso, isoTrksMtw, AnaConsts::isoHadTrksArr ) ) passIsoTrks = true;
  }
  return passIsoTrks;
}

