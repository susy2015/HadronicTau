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
#include "SignalCont.h"
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
      std::cerr <<" ./SignalCont Signal_T2tt_mStop500_mLSP325 TTbarInc_Template.root 1000 0 1" << std::endl;
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
  //Searchbin                                                                                                                                                                      
  SearchBins SB("SB_69_2016");
  //Use BaselineVessel class for baseline variables and selections
  //std::string spec = "SignalCont";
  std::string spec = "usegenmet";
  ExpBaselineVessel = new BaselineVessel(spec , "fastsim");
  AnaFunctions::prepareForNtupleReader();
  AnaFunctions::prepareTopTagger();
  NTupleReader *tr =0;
  tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames);
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
    
    double met=tr->getVar<double>("met");
    double metphi=tr->getVar<double>("metphi");
    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);
    int run = tr->getVar<int>("run");
    int lumi = tr->getVar<int>("lumi");
    int event = tr->getVar<int>("event");
    int nElectrons = tr->getVar<int>("nElectrons_CUT"+spec);
    int nMuons = tr->getVar<int>("nMuons_CUT"+spec);
    bool passNoiseEventFilter = tr->getVar<bool>("passNoiseEventFilter"+spec);

    double stopmass = tr->getVar<double>("SusyMotherMass");
    double lspmass = tr->getVar<double>("SusyLSPMass");

    //Event Filter                                                                                                                                                                                                                         
    if(!passNoiseEventFilter) continue;
    
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
    
    //Prediction part	
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
	const double muact = muonspfActivity[isomuonsIdxVec.at(0)];
	// iSR: this should be determined by search region requirement
	const int kSR = SB.find_Binning_Index(cnt1CSVS, nTopCandSortedCnt_pre, MT2_pre, simmet);
	//correction factor:                                                                                              
	const double corrBRWToTauHad = 0.65;  // Correction for the BR of hadronic tau decays                    
	const double corrBRTauToMu = 1-Efficiency::SBtaumucorMix(kSR);//correction from tauonic mu contamination 
	const double corrMuRecoEff = 1./Efficiency::reco(Efficiency::Ptbin(muLVec.Pt()), Efficiency::Actbin(muact)); // Correction for muon reconstruction efficiency                     	
	//const double corrBRTauToMu = 1-Efficiency::taumucorMix(Efficiency::Njetbin(combNJetPt30Eta24), Efficiency::metbin(simmet));//correction from tauonic mu contamination   
	const double corrMuIsoEff = 1./Efficiency::iso(Efficiency::Ptbin(muLVec.Pt()), Efficiency::Actbin(muact)); // Correction for muon isolation efficiency
	const double corrMuAcc = 1./Efficiency::SBaccMix(kSR); // Correction for muon acceptance
	//const double corrmtWEff =  1./Efficiency::mtwMix(Efficiency::Njetbin(combNJetPt30Eta24), Efficiency::metbin(simmet)); //correction for mtW cut
	const double corrmtWEff =  1./Efficiency::SBmtwMix(kSR); //correction for mtW cut
	//const double corrisotrkEff = 1- Efficiency::isotrkeffMix_NjetNbjet(Efficiency::Njetbin(combNJetPt30Eta24), Efficiency::NBjetbin(cnt1CSVS));//correction for isotrackveto eff.
	const double corrisotrkEff = 1- Efficiency::SBisotrkeffMix(kSR);//correction for isotrackveto eff.
	const double corrSBTrgeff = Efficiency::HTMHT_trgEff(simHt, Efficiency::Trgmetbin(simmet));//correction for search trigger efficiency
        const double corrMuTrkSF = Efficiency::MuonTrkSF(Efficiency::etabin(muLVec.Eta()));
        const double corrMuIdSF = Efficiency::MuonIDSF(Efficiency::PtbinIDSF(muLVec.Pt()), Efficiency::etabinIDSF(fabs(muLVec.Eta())));

	//The overall correction factor                                                                                                       
	const double corr = corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff * corrmtWEff * corrisotrkEff * corrBRTauToMu * corrSBTrgeff * corrMuTrkSF * corrMuIdSF * weight;
	const double Evt_corr = corr * Lumiscale;

	//Prediction yield
	int hstidx = myBaseHistgram.find_idx(stopmass, lspmass);
	if( pass_mtw && passhtpred && passMT2pred){
	  // cout<<"stopmass: "<<stopmass<<"\tlspmass: "<<lspmass<<endl;
	  if( kSR!=-1) {
	    //correction in each SB  
	    const double SBcorr = corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff * corrmtWEff * corrisotrkEff * corrBRTauToMu * weight;
	    const double Evt_SBcorr = SBcorr * Lumiscale;
	    myBaseHistgram.hPredYields_wt.at(hstidx)->Fill(kSR, Evt_SBcorr);
	  }
	}
		
      }//template bin loop
    }//control sample loop
    
    //correct the uncertainties in pred histo
    for(int i=0; i<myBaseHistgram.hPredYields_wt.size(); i++){
      TauResponse::Histfill(myBaseHistgram.hPredYields_wt.at(i), myBaseHistgram.hPredYields.at(i));
    }
  }	//event loop
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();
  
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

