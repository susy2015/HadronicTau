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
  //Searchbin                                                                                                                                                                              
  SearchBins SB("SB_59_2016");
  //Use BaselineVessel class for baseline variables and selections
  std::string spec = "ClosureExp";
  ExpBaselineVessel = new BaselineVessel(spec);
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

    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);
    int run = tr->getVar<int>("run");
    int lumi = tr->getVar<int>("lumi");
    int event = tr->getVar<int>("event");
    int nElectrons = tr->getVar<int>("nElectrons_CUT"+spec);
    int nMuons = tr->getVar<int>("nMuons_CUT"+spec);
    bool passNoiseEventFilter = tr->getVar<bool>("passNoiseEventFilter"+spec);
    const double EvtWt = tr->getVar<double>("evtWeight");
    //change event weight for MC sample
    EventWeight = EvtWt;
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
	    int kSR = SB.find_Binning_Index(nbJets_tru, nTops_tru, MT2_tru, met);
	    myBaseHistgram.hnodihadtau->Fill(kSR, Lumiscale);
	    if(isdihadtau)myBaseHistgram.hdihadtau->Fill(kSR, Lumiscale);
	  }
	}
      }
      //Exp Dist.
      if(W_tau_prongsVec.size() !=0 && passBaselineClosure && passNoiseEventFilter){
	int jSR = SB.find_Binning_Index(nbJets_tru, nTops_tru, MT2_tru, met);
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
	  if( TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT350_v") != string::npos || TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT400_v")!= string::npos || TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT600_v")!= string::npos || TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v")!= string::npos ||TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT400_PFMET50_v")!= string::npos){
	    //if( TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT350_v") != string::npos){
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
      //jec, jer corr. on mtW
      /*      const std::vector<double> &metMagUp = tr->getVec<double>("metMagUp");
      const std::vector<double> &metMagDown = tr->getVec<double>("metMagDown");
      double metjecUp = metMagUp[1];//met Uncertainty
      double metjecLow = metMagDown[1];//met Uncertainty
      double metjerUp = metMagUp[0];//met Uncertainty
      double metjerLow = metMagDown[0];//met Uncertainty
      TLorentzVector metLVecjecUp;metLVecjecUp.SetPtEtaPhiM(metjecUp, 0, metphi, 0);
      TLorentzVector metLVecjecLow;metLVecjecLow.SetPtEtaPhiM(metjecLow, 0, metphi, 0);
      TLorentzVector metLVecjerUp;metLVecjerUp.SetPtEtaPhiM(metjerUp, 0, metphi, 0);
      TLorentzVector metLVecjerLow;metLVecjerLow.SetPtEtaPhiM(metjerLow, 0, metphi, 0);
      const double mtwjecUp = calcMT(muLVec, metLVecjecUp);
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
	
	
	//Apply baseline cut
	if(!passNJetPred)continue;
	if(!passdeltaPhi) continue;
	if(!passmetPred)continue;
	if(!passbJets) continue;
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

	//if(!passMT2pred) continue;  	
	//Activity variable calculation:
	//double muact = AnaFunctions::getMuonActivity(muLVec, jetsLVec, recoJetschargedHadronEnergyFraction, recoJetschargedEmEnergyFraction, AnaConsts::muonsAct);
	const double muact = muonspfActivity[isomuonsIdxVec.at(0)];
	// iSR: this should be determined by search region requirement
	const int kSR = SB.find_Binning_Index(cnt1CSVS, nTopCandSortedCnt_pre, MT2_pre, simmet);
	
	//correction factor:                                                                                                                      
	const double corrBRWToTauHad = 0.65;  // Correction for the BR of hadronic tau decays                                                   
	//	const double corrBRTauToMu = 1-Efficiency::taumucorMix(Efficiency::Njetbin(combNJetPt30Eta24), Efficiency::metbin(simmet));//correction from tauonic mu contamination   
	const double corrBRTauToMu = 1-Efficiency::SBtaumucorMix(kSR);//correction from tauonic mu contamination   
	const double corrMuRecoEff = 1./Efficiency::reco(Efficiency::Ptbin(muLVec.Pt()), Efficiency::Actbin(muact)); // Correction for muon reconstruction efficiency                
	const double corrMuIsoEff = 1./Efficiency::iso(Efficiency::Ptbin(muLVec.Pt()), Efficiency::Actbin(muact)); // Correction for muon isolation efficiency
	const double corrMuAcc = 1./Efficiency::SBaccMix(kSR); // Correction for muon acceptance
	//	const double corrmtWEff =  1./Efficiency::mtwMix(Efficiency::Njetbin(combNJetPt30Eta24), Efficiency::metbin(simmet)); //correction for mtW cut
	const double corrmtWEff =  1./Efficiency::SBmtwMix(kSR); //correction for mtW cut
	//	const double corrisotrkEff = 1- Efficiency::isotrkeffMix_NjetNbjet(Efficiency::Njetbin(combNJetPt30Eta24), Efficiency::NBjetbin(cnt1CSVS));//correction for isotrackveto eff.
	const double corrisotrkEff = 1- Efficiency::SBisotrkeffMix(kSR);//correction for isotrackveto eff.

	const double corrllcont = 1-llcont;//correction for lost lepton overlap
	const double corrCSTrgeff = trgeff;//correction for CS trigger efficiency
	const double corrSBTrgeff = Efficiency::HTMHT_trgEff(simHt, Efficiency::Trgmetbin(simmet));//correction for search trigger efficiency
	const double corrMuTrkSF = Efficiency::MuonTrkSF(Efficiency::etabin(muLVec.Eta()));
	const double corrMuIdSF = Efficiency::MuonIDSF(Efficiency::PtbinIDSF(muLVec.Pt()), Efficiency::etabinIDSF(fabs(muLVec.Eta())));
	cout<<"Pt: "<<muLVec.Pt()<<"\tEta: "<<muLVec.Eta()<<"\tcorrMuTrkSF: "<<corrMuTrkSF<<"\tcorrMuIdSF: "<<corrMuIdSF<<endl;
	//The overall correction factor      
	//const double corr = corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff * corrmtWEff * corrisotrkEff * corrBRTauToMu  * weight;
	const double corr = corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff * corrmtWEff * corrisotrkEff * corrBRTauToMu * corrllcont * corrCSTrgeff * corrSBTrgeff * corrMuTrkSF * corrMuIdSF * weight;
	const double Evt_weight = Lumiscale * weight;
	const double Evt_corr = corr * Lumiscale;

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
	}
	
	//Fill search bin prediction
	if( pass_mtw && passhtpred && passMT2pred){
	  if( kSR!=-1) {
	    myBaseHistgram.hPredYields_wt->Fill(kSR, Evt_corr);
	  }
	}

	//mtW correction
	if(!isData){ 
	  if( !istaumu_genRecoMatch && passhtpred && passMT2pred){
	    Fill2D(myBaseHistgram.hnomtW_Njetmet_wt, combNJetPt30Eta24, simmet, Evt_weight);
	    FillDouble(myBaseHistgram.hnomtW_met, simmet, Evt_weight);
	    FillDouble(myBaseHistgram.hnomtW_MT2, MT2_pre, Evt_weight);
	    FillInt(myBaseHistgram.hnomtW_Njet, combNJetPt30Eta24, Evt_weight);
	    FillInt(myBaseHistgram.hnomtW_Nbjet, cnt1CSVS, Evt_weight);
	    FillInt(myBaseHistgram.hnomtW_Ntop, nTopCandSortedCnt_pre, Evt_weight);
	    if(pass_mtw){
	      Fill2D(myBaseHistgram.hmtW_Njetmet_wt, combNJetPt30Eta24, simmet, Evt_weight);
	      FillDouble(myBaseHistgram.hmtW_met, simmet, Evt_weight);
	      FillDouble(myBaseHistgram.hmtW_MT2, MT2_pre, Evt_weight);
	      FillInt(myBaseHistgram.hmtW_Njet, combNJetPt30Eta24, Evt_weight);
	      FillInt(myBaseHistgram.hmtW_Nbjet, cnt1CSVS, Evt_weight);
	      FillInt(myBaseHistgram.hmtW_Ntop, nTopCandSortedCnt_pre, Evt_weight);
	      }
	   
	    if(kSR!=-1){
	      myBaseHistgram.hnomtW_wt->Fill(kSR, Evt_weight);
	      // myBaseHistgram.hnomtWSys->Fill(kSR, Evt_weight);
	      if(pass_mtw) myBaseHistgram.hmtW_wt->Fill(kSR, Evt_weight);
	      /*if(mtwjecUp<100)myBaseHistgram.hmtWSysjecUp->Fill(kSR, Evt_weight);
	      if(mtwjecLow<100)myBaseHistgram.hmtWSysjecLow->Fill(kSR, Evt_weight);
	      if(mtwjerUp<100)myBaseHistgram.hmtWSysjerUp->Fill(kSR, Evt_weight);
	      if(mtwjerLow<100)myBaseHistgram.hmtWSysjerLow->Fill(kSR, Evt_weight);
	      */
	    }
	  }
	  //tau mu contamination
	  if( passhtpred && passMT2pred){
	    Fill2D(myBaseHistgram.hnotaumu_Njetmet_wt, combNJetPt30Eta24, simmet, Evt_weight);
	    FillDouble(myBaseHistgram.hnotaumu_met, simmet, Evt_weight);
            FillDouble(myBaseHistgram.hnotaumu_MT2, MT2_pre, Evt_weight);
            FillInt(myBaseHistgram.hnotaumu_Njet, combNJetPt30Eta24, Evt_weight);
            FillInt(myBaseHistgram.hnotaumu_Nbjet, cnt1CSVS, Evt_weight);
            FillInt(myBaseHistgram.hnotaumu_Ntop, nTopCandSortedCnt_pre, Evt_weight);
	    if(istaumu_genRecoMatch) {
	      Fill2D(myBaseHistgram.htaumu_Njetmet_wt, combNJetPt30Eta24, simmet, Evt_weight);
	      FillDouble(myBaseHistgram.htaumu_met, simmet, Evt_weight);
	      FillDouble(myBaseHistgram.htaumu_MT2, MT2_pre, Evt_weight);
	      FillInt(myBaseHistgram.htaumu_Njet, combNJetPt30Eta24, Evt_weight);
	      FillInt(myBaseHistgram.htaumu_Nbjet, cnt1CSVS, Evt_weight);
	      FillInt(myBaseHistgram.htaumu_Ntop, nTopCandSortedCnt_pre, Evt_weight);
	    }	    
	    if(kSR!=-1){
	      myBaseHistgram.hnotaumu_wt->Fill(kSR, Evt_weight);
	      if(istaumu_genRecoMatch) myBaseHistgram.htaumu_wt->Fill(kSR, Evt_weight * corrBRTauToMu);
	    }
	  }
	}
	
      }//template bin loop
    }//control sample loop

    //correct the uncertainties in pred histo
    TauResponse::Histfill(myBaseHistgram.hPredYields_wt, myBaseHistgram.hPredYields);
    TauResponse::Histfill(myBaseHistgram.hPredmet_wt, myBaseHistgram.hPredmet);
    TauResponse::Histfill(myBaseHistgram.hPredmht_wt, myBaseHistgram.hPredmht);
    TauResponse::Histfill(myBaseHistgram.hPredMT2_wt, myBaseHistgram.hPredMT2);
    TauResponse::Histfill(myBaseHistgram.hPredNbJets_wt, myBaseHistgram.hPredNbJets);
    TauResponse::Histfill(myBaseHistgram.hPredNTops_wt, myBaseHistgram.hPredNTops);
    TauResponse::Histfill(myBaseHistgram.hPreddPhi0_wt, myBaseHistgram.hPreddPhi0);  
    TauResponse::Histfill(myBaseHistgram.hPreddPhi1_wt, myBaseHistgram.hPreddPhi1);  
    TauResponse::Histfill(myBaseHistgram.hPreddPhi2_wt, myBaseHistgram.hPreddPhi2);  
    TauResponse::Histfill(myBaseHistgram.hnomtW_wt, myBaseHistgram.hnomtW);    
    TauResponse::Histfill(myBaseHistgram.hmtW_wt, myBaseHistgram.hmtW);    
    TauResponse::Histfill(myBaseHistgram.hnotaumu_wt, myBaseHistgram.hnotaumu);
    TauResponse::Histfill(myBaseHistgram.htaumu_wt, myBaseHistgram.htaumu);
    /*    TauResponse::Histfill2D(myBaseHistgram.hnomtW_Njetmet_wt, myBaseHistgram.hnomtW_Njetmet);
    TauResponse::Histfill2D(myBaseHistgram.hmtW_Njetmet_wt, myBaseHistgram.hmtW_Njetmet);
    TauResponse::Histfill2D(myBaseHistgram.hnotaumu_Njetmet_wt, myBaseHistgram.hnotaumu_Njetmet);
    TauResponse::Histfill2D(myBaseHistgram.htaumu_Njetmet_wt, myBaseHistgram.htaumu_Njetmet);
    */
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

