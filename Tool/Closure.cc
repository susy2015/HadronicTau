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
#include "SusyAnaTools/Tools/samples.h"
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
  if (argc < 4)
    {
      std::cerr <<"Please give 4 arguments " << "inputList " << " " <<" "<<"input template"<<" "<< "outputFileName" << "Lumiscale factor" << std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Closure ListSpring15_ttbar.txt HadTau_TauResponseTemplates.root HadTau_Closure.root #lumiscale" << std::endl;
      return -1;
    }
  const char *inputFileList = argv[1];
  const char *respTempl = argv[2];  
  const char *outFileName = argv[3];
  const char *lumiscale = argv[4];

  const double Lumiscale = std::atof(lumiscale);
  
  TChain *fChain = new TChain("stopTreeMaker/AUX");
  if(!FillChain(fChain, inputFileList))
    {
      std::cerr << "Cannot get the tree " << std::endl;
    }

  //Use BaselineVessel class for baseline variables and selections
  std::string spec = "HadtauExp";
  ExpBaselineVessel = new BaselineVessel(spec);

  NTupleReader tr(fChain);
  AnaFunctions::prepareTopTagger();
  tr.registerFunction(&passBaselineFuncExp);
// Add cleanJets function
  stopFunctions::cjh.setMuonIso("mini");
  stopFunctions::cjh.setElecIso("mini");
  stopFunctions::cjh.setRemove(false);
  tr.registerFunction(&stopFunctions::cleanJets);

  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(outFileName);

  TauResponse tauResp(respTempl);

  TRandom3 * rndm = new TRandom3(12345);
  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;

  int entries = tr.getNEntries();
  std::cout<<"\nentries : "<<entries<<std::endl; 

  // Loop over the events (tree entries)
  int k = 0;

  while(tr.getNextEvent()){
    k++;
    
// Add print out of the progress of looping
      if( tr.getEvtNum()-1 == 0 || tr.getEvtNum() == entries || (tr.getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr.getEvtNum()-1<<"th event ..."<<std::endl;

    const vector<TLorentzVector> &muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");
    const vector<double> &muonsRelIso = tr.getVec<double>("muonsRelIso");
    const vector<double> &muonsMiniIso = tr.getVec<double>("muonsMiniIso");
    const vector<double> &muonsMtw = tr.getVec<double>("muonsMtw");
    const vector<TLorentzVector> &genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
    const vector<int> &genDecayIdxVec = tr.getVec<int>("genDecayIdxVec");
    const vector<int> &genDecayPdgIdVec = tr.getVec<int>("genDecayPdgIdVec");
    const vector<int> &W_emuVec = tr.getVec<int>("W_emuVec");
    const vector<int> &W_tau_emuVec = tr.getVec<int>("W_tau_emuVec");
    const vector<TLorentzVector> &jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
    const vector<double> &recoJetsBtag_0 = tr.getVec<double>("recoJetsBtag_0");
    const vector<int> &looseisoTrksMatchedJetIdx = tr.getVec<int>("looseisoTrksMatchedJetIdx");
    const vector<TLorentzVector> &loose_isoTrksLVec = tr.getVec<TLorentzVector>("loose_isoTrksLVec");
    const vector<double> &loose_isoTrks_iso = tr.getVec<double>("loose_isoTrks_iso");
    const vector<double> &loose_isoTrks_mtw = tr.getVec<double>("loose_isoTrks_mtw");
    const vector<int> &loose_isoTrks_pdgId = tr.getVec<int>("loose_isoTrks_pdgId");
    const vector<int> &W_tau_prongsVec = tr.getVec<int>("W_tau_prongsVec");    
    const std::vector<int> & muonsFlagIDVec = tr.getVec<int>("muonsFlagMedium");
    double met=tr.getVar<double>("met");
    double metphi=tr.getVar<double>("metphi");
    double ht=tr.getVar<double>("ht");
    int run = tr.getVar<int>("run");
    int lumi = tr.getVar<int>("lumi");
    int event = tr.getVar<int>("event");
    
    int nElectrons = tr.getVar<int>("nElectrons_CUT2"+spec);
    int nMuons = tr.getVar<int>("nMuons_CUT2"+spec);
    bool passNoiseEventFilter = tr.getVar<bool>("passNoiseEventFilter"+spec);


    //Expectation part -- do it before prediction & do NOT skipping events
    bool passBaseline_tru = tr.getVar<bool>("passBaseline"+spec);
    bool passLeptVeto_tru = tr.getVar<bool>("passLeptVeto"+spec);
    bool passIsoTrkVeto_tru = tr.getVar<bool>("passIsoTrkVeto"+spec);
    bool passnJets_tru = tr.getVar<bool>("passnJets"+spec);
    bool passdPhis_tru = tr.getVar<bool>("passdPhis"+spec);
    bool passMET_tru =  tr.getVar<bool>("passMET"+spec);
    bool passBJets_tru = tr.getVar<bool>("passBJets"+spec);
    bool passTagger_tru = tr.getVar<bool>("passTagger"+spec);
    const int nJets_tru = tr.getVar<int>("cntNJetsPt30Eta24"+spec);
    const int nbJets_tru = tr.getVar<int>("cntCSVS"+spec);
    const int nTops_tru = tr.getVar<int>("nTopCandSortedCnt"+spec);
    const double MT2_tru = tr.getVar<double>("MT2_new"+spec);
    const vector<double> dPhiVec_tru = tr.getVec<double>("dPhiVec"+spec);

    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);

    //Event Filter

    //Select only events where the W decayed into a hadronically decaying tau
    //Note that for ttbar this includes (for two W's)
    // ] W->tau->had, W->tau->had
    // ] W->tau->had, W->qq
    // ] W->tau->had, W->e/mu (e or mu is lost)

    /*if(W_tau_prongsVec.size()!=0){
     myBaseHistgram.hTruecutFlow->Fill("original", Lumiscale);
     myBaseHistgram.hTruecutFlow->Fill("leptVeto", passLeptVeto_tru * Lumiscale);
     myBaseHistgram.hTruecutFlow->Fill("isotrackVeto", passLeptVeto_tru * passIsoTrkVeto_tru * Lumiscale);
      myBaseHistgram.hTruecutFlow->Fill("nJets", passLeptVeto_tru * passnJets_tru * Lumiscale);
      myBaseHistgram.hTruecutFlow->Fill("dPhis", passLeptVeto_tru * passnJets_tru * passdPhis_tru * Lumiscale);
      myBaseHistgram.hTruecutFlow->Fill("met", passLeptVeto_tru * passnJets_tru * passdPhis_tru * passMET_tru * Lumiscale);
      myBaseHistgram.hTruecutFlow->Fill("bJets", passLeptVeto_tru * passnJets_tru * passdPhis_tru * passMET_tru * passBJets_tru * Lumiscale);
      myBaseHistgram.hTruecutFlow->Fill("tagger", passLeptVeto_tru * passnJets_tru * passdPhis_tru * passMET_tru * passBJets_tru * passTagger_tru * Lumiscale);
        }
      if(W_tau_prongsVec.size()!=0 && passnJets_tru &&  passMET_tru && passBJets_tru &&  passTagger_tru){
    
      myBaseHistgram.hTruedPhi0->Fill(dPhi_0_tru, Lumiscale);
      myBaseHistgram.hTruedPhi1->Fill(dPhi_1_tru, Lumiscale);
      myBaseHistgram.hTruedPhi2->Fill(dPhi_2_tru, Lumiscale);
      }*/
    
      // ] W->tau->had, W->tau->had                                                                                                                                                                                                         
      // ] W->tau->had, W->qq                                                                                                                                                                                                               
      if(W_tau_prongsVec.size() !=0 && passBaseline_tru){
	int jSR = find_Binning_Index(nbJets_tru, nTops_tru, MT2_tru, met);
	if( jSR!= -1 ) {
	  myBaseHistgram.hTrueYields->Fill(jSR, Lumiscale);
	}
    
	FillDouble(myBaseHistgram.hTrueHt, ht, Lumiscale);
	FillDouble(myBaseHistgram.hTruemet, met, Lumiscale);
	FillDouble(myBaseHistgram.hTrueMT2, MT2_tru, Lumiscale);
	FillInt(myBaseHistgram.hTrueNJets, nJets_tru, Lumiscale);
	FillInt(myBaseHistgram.hTrueNbJets, nbJets_tru, Lumiscale);
	FillInt(myBaseHistgram.hTrueNTops, nTops_tru, Lumiscale);	
      }
      
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

	    // Find events that contain W->tau->mu
	    // Note that any code using gen info should be checked if they work for data or not!
	    bool istaumu = false, istaumu_genRecoMatch = false;
	    for(unsigned int ig=0; ig<W_tau_emuVec.size(); ig++){
	      int genIdx = W_tau_emuVec.at(ig);
	      if( std::abs(genDecayPdgIdVec.at(genIdx)) == 13 ){
		istaumu = true;
		const TLorentzVector & genLVec = genDecayLVec.at(genIdx);
		if( muLVec.DeltaR(genLVec) < 0.2 ) istaumu_genRecoMatch = true;
	      }
	    }

            // "Cross cleaning": find the jet that corresponds to the muon                                                                    
            const std::vector<TLorentzVector>& cleanJetVec      = tr.getVec<TLorentzVector>("cleanJetVec");
            const std::vector<double>& cleanJetBtag             = tr.getVec<double>("cleanJetBTag");
            // Get the cleaned jet indice (pointing to the jetsLVec) for the corresponding muons                                              
            const std::vector<int>& rejectJetIdx_formuVec = tr.getVec<int>("rejectJetIdx_formuVec");

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
	      // If simulted tau-jet meets same criteria as as HT jets,
	      // recompute HT and MHT
	      
	      double simHt = AnaFunctions::calcHT(combNJetVec, AnaConsts::pt50Eta24Arr);
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
	      
	      std::vector<double> * deltaPhiVec = new std::vector<double>();
	      (*deltaPhiVec) = AnaFunctions::calcDPhi(combNJetVec, simmetPhi, 3, AnaConsts::dphiArr);
	      double dPhi0_pred = deltaPhiVec->at(0);
	      double dPhi1_pred = deltaPhiVec->at(1);
	      double dPhi2_pred = deltaPhiVec->at(2);
	      bool passdeltaPhi = true;
	      if( deltaPhiVec->at(0) < AnaConsts::dPhi0_CUT || deltaPhiVec->at(1) < AnaConsts::dPhi1_CUT || deltaPhiVec->at(2) < AnaConsts::dPhi2_CUT){
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
	      std::vector<TLorentzVector> *jetsLVec_forTagger_pre = new std::vector<TLorentzVector>(); std::vector<double> *recoJetsBtag_forTagger_pre = new std::vector<double>();
	      AnaFunctions::prepareJetsForTagger(combNJetVec, combJetsBtag, (*jetsLVec_forTagger_pre), (*recoJetsBtag_forTagger_pre));
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
		type3Ptr->processEvent((*jetsLVec_forTagger_pre), (*recoJetsBtag_forTagger_pre), simmetLVec);
		nTopCandSortedCnt_pre = type3Ptr->nTopCandSortedCnt;
		MT2_pre = type3Ptr->best_had_brJet_MT2;
		mTcomb_pre = type3Ptr->best_had_brJet_mTcomb;
	      }
	      bool passTopTagger = type3Ptr->passNewTaggerReq() && nTopCandSortedCnt_pre >= AnaConsts::low_nTopCandSortedSel;
	      
	      if(!passTopTagger) continue;
	      
	      
	      //Activity variable calculation:
	      double muact = AnaFunctions::getMuonActivity(muLVec, jetsLVec, tr.getVec<double>("recoJetschargedHadronEnergyFraction"), tr.getVec<double>("recoJetschargedEmEnergyFraction"),AnaConsts::muonsAct);
	      
	      // iSR: this should be determined by search region requirement
	      const int kSR = find_Binning_Index(cnt1CSVS, nTopCandSortedCnt_pre, MT2_pre, simmet);
	      
	      //correction factor:                                                                                                                      
	      const double corrBRWToTauHad = 0.65;  // Correction for the BR of hadronic tau decays                                                   
	      //      const double corrBRTauToMu = Efficiency::taumucor(Efficiency::Ptbin1(muLVec.Pt()));//correction from tauonic mu contamination   
	      const double corrMuRecoEff = 1./Efficiency::reco(Efficiency::Ptbin(muLVec.Pt()), Efficiency::Actbin(muact)); // Correction for muon reconstruction efficiency                                                                                                                         
	      const double corrMuIsoEff = 1./Efficiency::iso(Efficiency::Ptbin(muLVec.Pt()), Efficiency::Actbin(muact)); // Correction for muon isolation efficiency                                                                                                                                
	      const double corrMuAcc = 1./Efficiency::acc(Efficiency::Njetbin(combNJetPt30Eta24)); // Correction for muon acceptance
	      //const double corrMuAcc = kSR==-1 ? 1/0.733 : 1./Efficiency::SBaccTT(kSR); // Correction for muon acceptance
	      const double corrmtWEff = kSR==-1 ? 1./Efficiency::MT2Bin_mtwcorr() : 1./Efficiency::SBmtwTT(kSR); //correction for mtW cut
	      const double corrisotrkEff = kSR==-1 ? (1-Efficiency::MT2Bin_isotrkcorr()) :(1- Efficiency::SBisotrkeffTT(kSR));//correction for isotrackveto eff.
	      //const double corrisotrkEff = 1- Efficiency::isotrkeffTT(Efficiency::NJetbin(combNJetPt30Eta24));//correction for isotrackveto eff.
	      //The overall correction factor                                                                                                       
	      //      const double corr = corrBRTauToMu * corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff;                                       
	      // For MC, no need of applying the corrBRTauToMu as you know if this event is from W->tau->mu or not
	      // However, we need know the fraction so that we can apply it to data (can be got from the printout)                                          
	      const double corr = corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff * corrmtWEff * corrisotrkEff * weight;
	      const double Evt_weight = Lumiscale * weight;
	      const double Evt_corr = corr * Lumiscale;
	      
	      if( !istaumu_genRecoMatch && pass_mtw){
	      myBaseHistgram.hPreddPhi0_wt->Fill(dPhi0_pred,Evt_corr);
	      myBaseHistgram.hPreddPhi1_wt->Fill(dPhi1_pred,Evt_corr);
	      myBaseHistgram.hPreddPhi2_wt->Fill(dPhi2_pred,Evt_corr);
	      }
	      //mtW correction in each SB
	      if(!istaumu_genRecoMatch){                                                                                            
		myBaseHistgram.hnomtW_wt->Fill(65, Evt_weight);
		if(pass_mtw) myBaseHistgram.hmtW_wt->Fill(65, Evt_weight);
	      
		if(cnt1CSVS>=1 && nTopCandSortedCnt_pre>=1 && simmet>=200 &&  MT2_pre<200){
		  myBaseHistgram.hnomtW_wt->Fill(64, Evt_weight);
		  if(pass_mtw) myBaseHistgram.hmtW_wt->Fill(64, Evt_weight);   
		}
	      
		if(kSR!=-1){
		  myBaseHistgram.hnomtW_wt->Fill(kSR, Evt_weight);
		  if(pass_mtw) myBaseHistgram.hmtW_wt->Fill(kSR, Evt_weight);
		}
	      }
	      
	      
	      // Fill the prediction distribution                                                                                                                 
	      if( !istaumu_genRecoMatch && pass_mtw){
		FillDouble(myBaseHistgram.hPredHt, simHt, Evt_corr);
		FillDouble(myBaseHistgram.hPredmet_wt, simmet, Evt_corr);
		FillDouble(myBaseHistgram.hPredMT2_wt, MT2_pre, Evt_corr);
		FillDouble(myBaseHistgram.hPredmTcomb, mTcomb_pre, Evt_corr);
		FillInt(myBaseHistgram.hPredNJets, combNJetPt30Eta24, Evt_corr);
		FillInt(myBaseHistgram.hPredNbJets_wt, cnt1CSVS, Evt_corr);
		FillInt(myBaseHistgram.hPredNTops_wt, nTopCandSortedCnt_pre, Evt_corr);
		}
	      
	      
	      //Fill search bin prediction
	      if( !istaumu_genRecoMatch && pass_mtw){
		if( kSR!=-1) {
		  //correction in each SB  
       
		  //const double SBcorrmtWEff = 1./Efficiency::SBmtwTT(kSR);
		  //const double SBcorrisotrkEff = 1-Efficiency::SBisotrkeffTT(kSR);
		  const double SBcorr = corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff * corrmtWEff * corrisotrkEff * weight;
		  const double Evt_SBcorr = SBcorr * Lumiscale;
		  myBaseHistgram.hPredYields_wt->Fill(kSR,Evt_SBcorr);
		}
	      }
	      //tau mu contamination
	      myBaseHistgram.hnotaumu_wt->Fill(65, Evt_weight);
	      if(istaumu_genRecoMatch) myBaseHistgram.htaumu_wt->Fill(65, Evt_weight);
	      
	      if(cnt1CSVS>=1 && nTopCandSortedCnt_pre>=1 && simmet>=200 &&  MT2_pre<200){
		myBaseHistgram.hnotaumu_wt->Fill(64, Evt_weight);
		if(istaumu_genRecoMatch) myBaseHistgram.htaumu_wt->Fill(64, Evt_weight);
	      }
	      
	      if(kSR!=-1){
		myBaseHistgram.hnotaumu_wt->Fill(kSR, Evt_weight);
		if(istaumu_genRecoMatch) myBaseHistgram.htaumu_wt->Fill(kSR, Evt_weight);
	      }
	     
	    }//template bin loop
	  }//control sample loop
	  //      }//spacial event criteria
	
	//correct the uncertainties in pred histo
      TauResponse::Histfill(myBaseHistgram.hPredYields_wt, myBaseHistgram.hPredYields);
      TauResponse::Histfill(myBaseHistgram.hnomtW_wt, myBaseHistgram.hnomtW);    
      TauResponse::Histfill(myBaseHistgram.hmtW_wt, myBaseHistgram.hmtW);    
      TauResponse::Histfill(myBaseHistgram.hnotaumu_wt, myBaseHistgram.hnotaumu);
      TauResponse::Histfill(myBaseHistgram.htaumu_wt, myBaseHistgram.htaumu);
      TauResponse::Histfill(myBaseHistgram.hPredmet_wt, myBaseHistgram.hPredmet);
      TauResponse::Histfill(myBaseHistgram.hPredMT2_wt, myBaseHistgram.hPredMT2);
      TauResponse::Histfill(myBaseHistgram.hPredNbJets_wt, myBaseHistgram.hPredNbJets);
      TauResponse::Histfill(myBaseHistgram.hPredNTops_wt, myBaseHistgram.hPredNTops);
      TauResponse::Histfill(myBaseHistgram.hPreddPhi0_wt, myBaseHistgram.hPreddPhi0);  
      TauResponse::Histfill(myBaseHistgram.hPreddPhi1_wt, myBaseHistgram.hPreddPhi1);  
      TauResponse::Histfill(myBaseHistgram.hPreddPhi2_wt, myBaseHistgram.hPreddPhi2);  

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

