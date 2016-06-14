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
#include <vector>
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/baselineDef.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "TStopwatch.h"
#include "TString.h"
#include "Acc.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "utils.h"
#include "TauResponse.h"
#include "Efficiency.h"
#include "TRandom3.h"

using namespace std;


// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 6)
    {
      std::cerr <<"Please give 5 arguments "<<"SubsampleName"<<" Input Template" <<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Acc TTbarSingleLep Template_2015.root 1000 0 1" << std::endl;
      return -1;
    }

  const char *subsamplename = argv[1];
  const char *respTempl = argv[2];
  const char *Maxevent = argv[3];
  const  char *Stratfile = argv[4];
  const  char *Filerun = argv[5];
  const  int startfile = std::atoi(Stratfile);
  const int filerun = std::atoi(Filerun);
  TChain *fChain = 0;
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(subsamplename, startfile);
  const string condorSpec = argc==7 ? argv[6]: "";
  if(!FillChain(fChain, subsamplename, condorSpec, startfile, filerun))
    {
      std::cerr << "Cannot get the tree " << std::endl;
    }
  const int maxevent = std::atoi(Maxevent);
  //Searchbin                                                                                                                                                                                     
  SearchBins SB("SB_69_2016");
  //Use BaselineVessel class for baseline variables and selections
  std::string spec = "Acc";
  AccBaselineVessel = new BaselineVessel(spec);
  AnaFunctions::prepareForNtupleReader();
  AnaFunctions::prepareTopTagger();
  NTupleReader *tr =0;
  tr = new NTupleReader(fChain);
  tr->registerFunction(&AccpassBaselineFunc);
  //Add cleanjet function and miniIsolation
  stopFunctions::cjh.setMuonIso("mini");
  stopFunctions::cjh.setElecIso("mini");
  stopFunctions::cjh.setRemove(false);
  tr->registerFunction(&stopFunctions::cleanJets);
  //Add PDF uncertainty
  tr->registerFunction(&mypdf);

  TauResponse tauResp(respTempl);

  TRandom3 * rndm = new TRandom3(12345);

  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr->getNEntries();
  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<Lumiscale<<std::endl;

  // Loop over the events (tree entries)
  int k = 0;
  while(tr->getNextEvent()){
    k++;

    if(maxevent>=0 && tr->getEvtNum() > maxevent ) break;
    if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;
    
    // Veto electrons -> to get acceptance in search regions as much as possible
    const double nElectrons = tr->getVar<double>("nElectrons_CUT"+spec);
    if( nElectrons != AnaConsts::nElectronsSel ) continue;

    const vector<TLorentzVector> &muonsLVec = tr->getVec<TLorentzVector>("muonsLVec");
    const vector<double> &muonsRelIso = tr->getVec<double>("muonsRelIso");
    const vector<double> &muonsMiniIso = tr->getVec<double>("muonsMiniIso");
    const vector<double> &muonsMtw = tr->getVec<double>("muonsMtw");
    const vector<TLorentzVector> &genDecayLVec = tr->getVec<TLorentzVector>("genDecayLVec");
    const vector<int> &genDecayIdxVec = tr->getVec<int>("genDecayIdxVec");
    const vector<int> &genDecayMomIdxVec = tr->getVec<int>("genDecayMomIdxVec");
    const vector<int> &genDecayPdgIdVec = tr->getVec<int>("genDecayPdgIdVec");
    const vector<string> &genDecayStrVec = tr->getVec<string>("genDecayStrVec");
    const vector<int> &W_emuVec = tr->getVec<int>("W_emuVec");
    const vector<int> &W_tau_emuVec = tr->getVec<int>("W_tau_emuVec");
    const vector<int> &W_tau_prongsVec = tr->getVec<int>("W_tau_prongsVec");
    const vector<int> &W_tau_nuVec = tr->getVec<int>("W_tau_nuVec");
    const vector<TLorentzVector> &jetsLVec = tr->getVec<TLorentzVector>("jetsLVec");
    const vector<double> &recoJetsBtag_0 = tr->getVec<double>("recoJetsBtag_0");
    const vector<int> & muonsFlagIDVec = tr->getVec<int>("muonsFlagMedium");
    const vector<int> &looseisoTrksMatchedJetIdx = tr->getVec<int>("looseisoTrksMatchedJetIdx");
    const vector<TLorentzVector> &loose_isoTrksLVec = tr->getVec<TLorentzVector>("loose_isoTrksLVec");
    const vector<double> &loose_isoTrks_iso = tr->getVec<double>("loose_isoTrks_iso");
    const vector<double> &loose_isoTrks_mtw = tr->getVec<double>("loose_isoTrks_mtw");
    const vector<int> &loose_isoTrks_pdgId = tr->getVec<int>("loose_isoTrks_pdgId");

    double met=tr->getVar<double>("met");
    double metphi=tr->getVar<double>("metphi");
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
    const int nJets_tru = tr->getVar<int>("cntNJetsPt30Eta24"+spec);
    const int nbJets_tru = tr->getVar<int>("cntCSVS"+spec);
    const int nTops_tru = tr->getVar<int>("nTopCandSortedCnt"+spec);
    const double MT2_tru = tr->getVar<double>("best_had_brJet_MT2"+spec);

    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);

    //tau acceptance
    bool passBaselineAcc = passMuonVeto_tru && passEleVeto_tru && passnJets_tru && passdPhis_tru && passMET_tru && passBJets_tru && passTagger_tru;
    if(W_tau_prongsVec.size()!=0){
      const unsigned int taumetbin = Efficiency::metbin(met);
      bool istauAcc = false;
      std::vector<TLorentzVector>genhadtauLVec;
      std::vector<TLorentzVector>genvisiblehadtauLVec;
      for(unsigned ig=0; ig<genDecayLVec.size(); ig++){
	int pdgId = genDecayPdgIdVec.at(ig);
	if(abs(pdgId)==15){
	  int flag=0;
	  if(W_tau_emuVec.size()!=0){
	    for(int k=0; k<W_tau_emuVec.size();k++){
	      int lIdx = W_tau_emuVec.at(k);
	      if( genDecayMomIdxVec.at(lIdx) == genDecayIdxVec.at(ig) )flag++;
	    }
	  }
	  if(!flag){
	    TLorentzVector gennu;
	    TLorentzVector genhadtau = genDecayLVec.at(ig);
	    for(int n=0; n<W_tau_nuVec.size();n++){
	      int nIdx = W_tau_nuVec.at(n);
	      if( genDecayMomIdxVec.at(nIdx) == genDecayIdxVec.at(ig)) gennu = genDecayLVec.at(nIdx);
	    }
	    genhadtauLVec.push_back(genhadtau);
	    genvisiblehadtauLVec.push_back(genhadtau-gennu);
	  }
	}
      }
      TLorentzVector usetauLVec = genhadtauLVec[0];
      if(passBaselineAcc && passMT2_tru && passHT_tru){
	//Isotrack on remaing part
	bool othrIsotrackVeto = false;
	if(looseisoTrksMatchedJetIdx.size()!=loose_isoTrksLVec.size())cout<<"Error: isotrack vetor size mismatch"<<endl;
	int nTotal = 0, nHadtau = 0;
	std::vector<int> IsotrkMatchedHadtau(genvisiblehadtauLVec.size(), 0);
	for(unsigned int it=0; it< looseisoTrksMatchedJetIdx.size();it++){
	  if(!passIsoTrks1(loose_isoTrksLVec[it], loose_isoTrks_iso[it], loose_isoTrks_mtw[it], loose_isoTrks_pdgId[it])) continue;
	  nTotal++;    
	  // Do the matching
	  const int isoJetIdx = looseisoTrksMatchedJetIdx[it]; // Will store the index of the jet matched to the isotrack
	  if(isoJetIdx == -1)cout<<"Event: "<<k<<" isoJetIdx: "<<isoJetIdx<<"  track pT: "<<loose_isoTrksLVec[it].Pt()<<"  track phi: "<<loose_isoTrksLVec[it].Phi()<<"  track Eta: "<<loose_isoTrksLVec[it].Eta()<<"  track id: "<<loose_isoTrks_pdgId[it]<<endl;
	  if(findTauMatchedisoJet(isoJetIdx,genvisiblehadtauLVec,jetsLVec, IsotrkMatchedHadtau)){
	    if(IsotrkMatchedHadtau[0]) nHadtau++; // considering only first had tau
	  }
	}//finish isotrack loop
	if((nTotal-nHadtau)==0) othrIsotrackVeto = true;//isotrackveto on the remaing part (exclude hadtau)
	if(usetauLVec.Pt()>TauResponse::ptMin() && fabs(usetauLVec.Eta())<TauResponse::etaMax()) istauAcc = true;
	int jSR = SB.find_Binning_Index(nbJets_tru,nTops_tru,MT2_tru,met);
	if(othrIsotrackVeto){
	  if(jSR!=-1)myBaseHistgram.hTaugen->Fill(jSR, Lumiscale);
	  if(istauAcc){
	    if(jSR!=-1)myBaseHistgram.hTauacc->Fill(jSR, Lumiscale);
	  }
	}
      }
    }
    //Muon Acceptance
    // Only use events with W->mu
    vector<TLorentzVector> genmuLVec;
    for(unsigned int ip=0; ip<W_emuVec.size(); ip++){
      int idx = W_emuVec.at(ip);
      int pdgId = genDecayPdgIdVec.at(idx);
      if( fabs(pdgId) == 13 ){
	genmuLVec.push_back(genDecayLVec.at(idx));
      }
    }
    
    std::vector<TLorentzVector> cleanJetVec;
    std::vector<double> cleanJetBtag;
    // ]W->mu W->qq                                                                   
    // ]W->mu W->tau->had 
    if(genmuLVec.size()==1){//control sample
      //pdf value
      const double pdfCentral = tr->getVar<double>("NNPDF_From_Median_Central");
      const double pdfUp = tr->getVar<double>("NNPDF_From_Median_Up");
      const double pdfDown = tr->getVar<double>("NNPDF_From_Median_Down");
      const double scaleUp = tr->getVar<double>("Scaled_Variations_Up");
      const double scaleDown = tr->getVar<double>("Scaled_Variations_Down");
      
      const TLorentzVector pergenmuLVec = genmuLVec.at(0);      
      bool useRecoMuon = false;
      int muJetIdx = -1;
      vector<TLorentzVector> isomuonsLVec;
      vector<int> isomuonsIdxVec;
      for(unsigned int im=0; im<muonsLVec.size(); im++){
	if( AnaFunctions::passMuon(muonsLVec.at(im), muonsMiniIso.at(im), muonsMtw.at(im), muonsFlagIDVec.at(im), AnaConsts::muonsMiniIsoArr) ){ isomuonsLVec.push_back(muonsLVec.at(im)); isomuonsIdxVec.push_back(im); }
      }
      // If there is a reco muon with pt and eta larger than requirement, use the reco muon
      // This way we get consistent results for wherever a reco muon is ready to be used
      if( nMuons == 1 ){
	if( isomuonsLVec.at(0).Pt() >= TauResponse::ptMin() && fabs(isomuonsLVec.at(0).Eta()) <= TauResponse::etaMax() ){
	  cleanJetVec = tr->getVec<TLorentzVector>("cleanJetVec");
	  cleanJetBtag = tr->getVec<double>("cleanJetBTag");
	  useRecoMuon = true;
	  const std::vector<int>& rejectJetIdx_formuVec = tr->getVec<int>("rejectJetIdx_formuVec");
	  muJetIdx = rejectJetIdx_formuVec.at(isomuonsIdxVec.at(0));
	}
      }
      
      // If no reco muon in the event, use the gen muon instead
      if( !useRecoMuon ){
	const float deltaRMax = 0.2;
	utils::findMatchedObject(muJetIdx, pergenmuLVec, jetsLVec, deltaRMax);
	for(int jetIdx = 0; jetIdx < jetsLVec.size(); ++jetIdx) { // Loop over reco jets
	  cleanJetVec.push_back(jetsLVec.at(jetIdx));
	  cleanJetBtag.push_back(recoJetsBtag_0.at(jetIdx));
	  // Substract muon part from this jet if its overlap with the muon   
	  if( jetIdx == muJetIdx ){
	    cleanJetVec[muJetIdx] =- pergenmuLVec;
	  }
	}
      }
      
      TLorentzVector usedmuLVec;
      if( useRecoMuon ) usedmuLVec = isomuonsLVec.at(0);
      else usedmuLVec = pergenmuLVec;
      
      // Determine if the gen muon passes the kinematic cuts or not
      bool passKinCuts = false;
      if(pergenmuLVec.Pt()>TauResponse::ptMin() && fabs(pergenmuLVec.Eta())<TauResponse::etaMax()) passKinCuts = true;
      
      TLorentzVector cleanmetLVec; cleanmetLVec.SetVectM( (metLVec + usedmuLVec).Vect(), 0);
      
      //Implement IsoTrackVeto
      if(looseisoTrksMatchedJetIdx.size()!=loose_isoTrksLVec.size())cout<<"Error: isotrack vetor size mismatch"<<endl;
      int nisoTotal = 0, nisomu = 0;    
      for(unsigned int it=0; it< looseisoTrksMatchedJetIdx.size();it++){
	if(!passIsoTrks1(loose_isoTrksLVec[it], loose_isoTrks_iso[it], loose_isoTrks_mtw[it], loose_isoTrks_pdgId[it])) continue;
	nisoTotal++;
	// Do the matching
	const double dr = usedmuLVec.DeltaR(loose_isoTrksLVec[it]);
	if(dr<0.1) nisomu++; // muon tagged as isotrack
      }//finish isotrack loop
      if((nisoTotal-nisomu)!=0)continue;//isotrackveto on the remaing part (exclude mu)  

      // Get random number from tau-response template
      //      const double scale = tauResp.getRandom(usedmuLVec.Pt());
      TH1F* temp = (TH1F*)tauResp.Resp(usedmuLVec.Pt());
      //Loop over template bin
      for(int ib = 1; ib<=50; ib++){
	const double scale = temp->GetBinCenter(ib);
	const double weight = temp->GetBinContent(ib) * temp->GetBinWidth(ib);
	
      // Scale muon pt with tau response --> simulate tau jet pt
	const double simTauJetPt = scale * usedmuLVec.Pt();
	const double simTauJetE = scale * usedmuLVec.E();

	TLorentzVector ori_tauJetLVec; ori_tauJetLVec.SetPtEtaPhiE(simTauJetPt, usedmuLVec.Eta(), usedmuLVec.Phi(), simTauJetE);
	TLorentzVector tauJetLVec = ori_tauJetLVec;

	// See comments for similar code in Closure.cc
	double oriJetCSVS = 0;
	if( muJetIdx != -1 ) oriJetCSVS = recoJetsBtag_0.at(muJetIdx);

	double mistag = Efficiency::mistag(Efficiency::Ptbin1(simTauJetPt));
	double rno = rndm->Rndm();
	if( rno < mistag) oriJetCSVS = 1.0;
	
	//Adjustment of tau jet to the remaining part of mu cleaned jet
	if( muJetIdx >=0 ) tauJetLVec += cleanJetVec[muJetIdx];

	vector<TLorentzVector> combJetVec;
	vector<double> combJetsBtag;

	// See comments for similar code in Closure.cc
	bool includeTauJet = false;
	for(unsigned int ij=0; ij<cleanJetVec.size(); ij++){
	  if( ij == muJetIdx ) continue;
	  if( tauJetLVec.Pt() > cleanJetVec.at(ij).Pt() && !includeTauJet ){
            combJetVec.push_back(tauJetLVec); combJetsBtag.push_back(oriJetCSVS);
            includeTauJet = true;
	  }
	  combJetVec.push_back(cleanJetVec.at(ij)); combJetsBtag.push_back(cleanJetBtag.at(ij));
	}
	// it's possible that the tau jet is the least energetic jet so that it's not added into the combNJetVec during the loop
	if( !includeTauJet ){ combJetVec.push_back(tauJetLVec); combJetsBtag.push_back(oriJetCSVS); }
	
	// recompute HT
	double Ht = AnaFunctions::calcHT(combJetVec, AnaConsts::pt30Eta24Arr);
	bool passht = true;
	if( Ht < AnaConsts::defaultHTcut )passht = false;

	//recompute met
	TLorentzVector combmetLVec;  combmetLVec.SetVectM( (cleanmetLVec - ori_tauJetLVec).Vect(), 0 );
	
	const double combmet = combmetLVec.Pt();
	const double combmetPhi = combmetLVec.Phi();
	
	bool passmet = true;
	if(combmet<AnaConsts::defaultMETcut) passmet = false; 

	//recompute jet
	int nJetPt30Eta24 = AnaFunctions::countJets(combJetVec, AnaConsts::pt30Eta24Arr);
	int nJetPt50Eta24 = AnaFunctions::countJets(combJetVec, AnaConsts::pt50Eta24Arr);
	
	bool passnjets = true;
	if(nJetPt30Eta24<AnaConsts::nJetsSelPt30Eta24) passnjets = false;
	if(nJetPt50Eta24<AnaConsts::nJetsSelPt50Eta24) passnjets = false;

	//recompute deltaphi
	std::vector<double> deltaPhiVec = AnaFunctions::calcDPhi(combJetVec, combmetPhi, 3, AnaConsts::dphiArr);
	bool passdeltaPhi = true;
	if( deltaPhiVec.at(0) < AnaConsts::dPhi0_CUT || deltaPhiVec.at(1) < AnaConsts::dPhi1_CUT || deltaPhiVec.at(2) < AnaConsts::dPhi2_CUT){
	  passdeltaPhi = false;
	}

      //recompute bjet
	int cnt1CSVS = AnaFunctions::countCSVS(combJetVec, combJetsBtag, AnaConsts::cutCSVS, AnaConsts::bTagArr);
	bool passbJets = true;
	if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cnt1CSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cnt1CSVS < AnaConsts::high_nJetsSelBtagged ) ) ){
	  passbJets = false;
	}

	if(!passnjets) continue;
	if(!passmet) continue;
	if(!passdeltaPhi) continue;
	if(!passbJets) continue;

	//top tagger input
	int comb30_acc = AnaFunctions::countJets(combJetVec, AnaConsts::pt30Arr);
	std::vector<TLorentzVector> jetsLVec_forTagger_acc; std::vector<double> recoJetsBtag_forTagger_acc;
	AnaFunctions::prepareJetsForTagger(combJetVec, combJetsBtag, (jetsLVec_forTagger_acc), (recoJetsBtag_forTagger_acc));
      
	int nTopCandSortedCnt_acc = -1;
	double MT2_acc = -1;
	
	//Apply Top tagger
	if(comb30_acc >= AnaConsts::nJetsSel ){
	  type3Ptr->processEvent((jetsLVec_forTagger_acc), (recoJetsBtag_forTagger_acc), combmetLVec);
	  nTopCandSortedCnt_acc = type3Ptr->nTopCandSortedCnt;
	  MT2_acc = type3Ptr->best_had_brJet_MT2;
	}
	bool passTopTagger = type3Ptr->passNewTaggerReq() && nTopCandSortedCnt_acc >= AnaConsts::low_nTopCandSortedSel;
	bool passMT2 = true;
	if(MT2_acc < AnaConsts::defaultMT2cut)passMT2 = false;

	if(!passTopTagger) continue;
	//500 GeV HT cut..
	if(!passht)continue;

	// iSR: this should be determined by search region requirement
	int iSR = SB.find_Binning_Index(cnt1CSVS, nTopCandSortedCnt_acc, MT2_acc, combmet);

	const double corrBRWToTauHad = 0.65;  // Correction for the BR of hadronic tau decays                          
	
	//const double corrMuAcc = 1./Efficiency::accMix_NjetMT2(Efficiency::Njetbin(nJetPt30Eta24), Efficiency::MT2bin(MT2_acc)); // Correction for muon acceptance
	const double corr = corrBRWToTauHad;
	const double Evt_weight = Lumiscale * weight * corr;

	//histogram
	if(passMT2){
	  if(iSR!=-1){
	    myBaseHistgram.hgen_wt->Fill(iSR, Evt_weight);
	    myBaseHistgram.hgen_pdfCentral_wt->Fill(iSR, Evt_weight*pdfCentral);
	    myBaseHistgram.hgen_pdfUp_wt->Fill(iSR, Evt_weight*pdfUp);
	    myBaseHistgram.hgen_pdfDown_wt->Fill(iSR, Evt_weight*pdfDown);
	    myBaseHistgram.hgen_scaleUp_wt->Fill(iSR, Evt_weight*scaleUp);
	    myBaseHistgram.hgen_scaleDown_wt->Fill(iSR, Evt_weight*scaleDown);
	  }
	  if(passKinCuts){
	    if(iSR!=-1){
	      myBaseHistgram.hacc_wt->Fill(iSR, Evt_weight);
	      myBaseHistgram.hacc_pdfCentral_wt->Fill(iSR, Evt_weight*pdfCentral);
	      myBaseHistgram.hacc_pdfUp_wt->Fill(iSR, Evt_weight*pdfUp);
	      myBaseHistgram.hacc_pdfDown_wt->Fill(iSR, Evt_weight*pdfDown);
	      myBaseHistgram.hacc_scaleUp_wt->Fill(iSR, Evt_weight*scaleUp);
	      myBaseHistgram.hacc_scaleDown_wt->Fill(iSR,Evt_weight*scaleDown);
	    }
	  }
	}
      }//template bin loop
    }//muon control sample loop

    //correct the uncertainties in pred histo
    TauResponse::Histfill(myBaseHistgram.hacc_wt, myBaseHistgram.hacc);
    TauResponse::Histfill(myBaseHistgram.hgen_wt, myBaseHistgram.hgen);
    /*    TauResponse::Histfill(myBaseHistgram.hacc_pdfCentral_wt, myBaseHistgram.hacc_pdfCentral);
    TauResponse::Histfill(myBaseHistgram.hacc_pdfUp_wt, myBaseHistgram.hacc_pdfUp);
    TauResponse::Histfill(myBaseHistgram.hacc_pdfDown_wt, myBaseHistgram.hacc_pdfDown);
    TauResponse::Histfill(myBaseHistgram.hacc_scaleUp_wt, myBaseHistgram.hacc_scaleUp);
    TauResponse::Histfill(myBaseHistgram.hacc_scaleDown_wt, myBaseHistgram.hacc_scaleDown);
    TauResponse::Histfill(myBaseHistgram.hgen_pdfCentral_wt, myBaseHistgram.hgen_pdfCentral);
    TauResponse::Histfill(myBaseHistgram.hgen_pdfUp_wt, myBaseHistgram.hgen_pdfUp);
    TauResponse::Histfill(myBaseHistgram.hgen_pdfDown_wt, myBaseHistgram.hgen_pdfDown);
    TauResponse::Histfill(myBaseHistgram.hgen_scaleUp_wt, myBaseHistgram.hgen_scaleUp);
    TauResponse::Histfill(myBaseHistgram.hgen_scaleDown_wt, myBaseHistgram.hgen_scaleDown);
    */
  }//event loop
  (myBaseHistgram.oFile)->Write();

  cout<<"ToTal Event: "<<k<<endl;
  return 0;
}
bool find_mother(int momIdx, int dauIdx, const std::vector<int> &genDecayIdxVec, const std::vector<int> &genDecayMomIdxVec){
  if( momIdx == -1 || dauIdx == -1 ) return false;
  if( dauIdx == momIdx ) return true;

  int thisIdx = dauIdx;
  while( thisIdx >=0 ){
    int momGenIdx = genDecayMomIdxVec[thisIdx];
    thisIdx = find_idx(momGenIdx, genDecayIdxVec);
    if( thisIdx == momIdx ) return true;
  }
  return false;
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
bool findTauMatchedisoJet(int matchedObjIdx, const std::vector<TLorentzVector> &genvisiblehadtauLVec, const std::vector<TLorentzVector> &jetsLVec, std::vector<int> &IsotrkMatchedHadtau){
  bool isisotrack = false;
  if(matchedObjIdx == -1){
    cout<<"Input matchedObjIdx is -1, please add extra function.."<<endl;
    return isisotrack;
  }
  TLorentzVector jetLVec = jetsLVec[matchedObjIdx];
  for(int is =0; is<genvisiblehadtauLVec.size(); is++){
    const double deltaRMax = deltaRmax(genvisiblehadtauLVec.at(is).Pt());  
    const double dr = jetLVec.DeltaR(genvisiblehadtauLVec.at(is));
    if( dr < deltaRMax ) {
      isisotrack = true;
      IsotrkMatchedHadtau[is]=1;
    }
  }
  return isisotrack;
}
double deltaRmax(double pt){
  double rmax = 0.25;
  if(pt>30) rmax = 0.2;
  if(pt>50) rmax = 0.15;
  if(pt>100) rmax = 0.1;
  return rmax;
}
