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
#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/baselineDef.h"
#include "TStopwatch.h"
#include "TString.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
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

using namespace std;


void passBaselineFunc(NTupleReader &tr)
{
  bool passBaseline_nolepveto = true;

  //Form TLorentzVector of MET
  TLorentzVector metLVec; metLVec.SetPtEtaPhiM(tr.getVar<double>("met"), 0, tr.getVar<double>("metphi"), 0);

  /*Calculate number of leptons
  int nMuons = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsRelIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsArr);
  int nElectrons = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesRelIso"), tr.getVec<double>("elesMtw"), AnaConsts::elesArr);
  int nIsoTrks = AnaFunctions::countIsoTrks(tr.getVec<TLorentzVector>("loose_isoTrksLVec"), tr.getVec<double>("loose_isoTrks_iso"), tr.getVec<double>("loose_isoTrks_mtw"), AnaConsts::isoTrksArr);
  */
  int nMuons = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsMiniIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsMiniIsoArr);
  int nElectrons = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesMiniIso"), tr.getVec<double>("elesMtw"), tr.getVec<unsigned int>("elesisEB"), AnaConsts::elesArr);

  //Calculate number of jets and b-tagged jets
  int cntCSVS = AnaFunctions::countCSVS(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), AnaConsts::cutCSVS, AnaConsts::bTagArr);
  int cntNJetsPt50Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt50Eta24Arr);
  int cntNJetsPt30Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Eta24Arr);
  int cntNJetsPt30 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Arr);

  //Calculate deltaPhi
  std::vector<double> * dPhiVec = new std::vector<double>();
  (*dPhiVec) = AnaFunctions::calcDPhi(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVar<double>("metphi"), 3, AnaConsts::dphiArr);

  //Pass number of jets?
  bool passnJets = true;
  if( cntNJetsPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 ){passBaseline_nolepveto = false; passnJets = false;}
  if( cntNJetsPt30Eta24 < AnaConsts::nJetsSelPt30Eta24 ){passBaseline_nolepveto = false; passnJets = false;}

  //Pass deltaPhi?
  bool passdPhis = true;
  if( dPhiVec->at(0) < AnaConsts::dPhi0_CUT || dPhiVec->at(1) < AnaConsts::dPhi1_CUT || dPhiVec->at(2) < AnaConsts::dPhi2_CUT ){passBaseline_nolepveto = false; passdPhis = false; }

  //Pass number of b-tagged jets?
  bool passBJets = true;
  if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cntCSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cntCSVS < AnaConsts::high_nJetsSelBtagged ) ) ){ passBaseline_nolepveto = false; passBJets = false; }

  //Pass the baseline MET requirement?
  bool passMET = true;
  if( tr.getVar<double>("met") < AnaConsts::defaultMETcut ){passBaseline_nolepveto = false; passMET = false; }

  //Register all the calculated variables
  tr.registerDerivedVar("nMuons_CUT2", nMuons);
  tr.registerDerivedVar("nElectrons_CUT2", nElectrons);
  tr.registerDerivedVar("cntNJetsPt30Eta24", cntNJetsPt30Eta24);
  tr.registerDerivedVar("passnJets", passnJets);
  tr.registerDerivedVar("passdPhis", passdPhis);
  tr.registerDerivedVar("passBJets", passBJets);
  tr.registerDerivedVar("passMET", passMET);
  tr.registerDerivedVar("passBaseline_nolepveto", passBaseline_nolepveto);

}

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 2)
    {
      std::cerr <<"Please give 2 arguments " << "inputList " <<" "<<"input template"<< std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Acc List1_ttbar.txt HadTau_TauResponseTemplates.root" << std::endl;
      return -1;
    }
  const char *inputFileList = argv[1];
  const char *respTempl = argv[2];

  TChain *fChain = new TChain("stopTreeMaker/AUX");
  if(!FillChain(fChain, inputFileList))
    {
      std::cerr << "Cannot get the tree " << std::endl;
    }
  NTupleReader tr(fChain);
  AnaFunctions::prepareTopTagger(); 
  tr.registerFunction(&passBaselineFunc);

  //Add cleanjet function and miniIsolatio
  stopFunctions::cjh.setMuonIso("mini");
  stopFunctions::cjh.setElecIso("mini");
  tr.registerFunction(&stopFunctions::cleanJets);

  TauResponse tauResp(respTempl);

  int genmu=0, accmu=0;

  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;

  int entries = tr.getNEntries();
  std::cout<<"\nentries : "<<entries<<std::endl;

  // Loop over the events (tree entries)
  int k = 0;
  while(tr.getNextEvent()){

    k++;

    if( tr.getEvtNum()-1 == 0 || tr.getEvtNum() == entries || (tr.getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr.getEvtNum()-1<<"th event ..."<<std::endl;

// Veto electrons -> to get acceptance in search regions as much as possible
    const double nElectrons = tr.getVar<double>("nElectrons_CUT2");
    if( nElectrons != AnaConsts::nElectronsSel ) continue;

    const vector<TLorentzVector> &muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");
    const vector<double> &muonsRelIso = tr.getVec<double>("muonsRelIso");
    const vector<double> &muonsMiniIso = tr.getVec<double>("muonsMiniIso");
    const vector<double> &muonsMtw = tr.getVec<double>("muonsMtw");
    const vector<TLorentzVector> &genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
    const vector<int> &genDecayIdxVec = tr.getVec<int>("genDecayIdxVec");
    const vector<int> &genDecayPdgIdVec = tr.getVec<int>("genDecayPdgIdVec");
    const vector<int> &W_emuVec = tr.getVec<int>("W_emuVec");
    const vector<TLorentzVector> &jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
    const vector<double> &recoJetsBtag_0 = tr.getVec<double>("recoJetsBtag_0");
    int nMuons = tr.getVar<int>("nMuons_CUT2");
    double met=tr.getVar<double>("met");
    double metphi=tr.getVar<double>("metphi");

    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);

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

    if(genmuLVec.size()==1){//control sample
      const TLorentzVector pergenmuLVec = genmuLVec.at(0);

      bool useRecoMuon = false;
      int muJetIdx = -1;

      vector<TLorentzVector> isomuonsLVec;
      vector<int> isomuonsIdxVec;
      for(unsigned int im=0; im<muonsLVec.size(); im++){
         if( AnaFunctions::passMuon(muonsLVec.at(im), muonsMiniIso.at(im), muonsMtw.at(im), AnaConsts::muonsMiniIsoArr) ){ isomuonsLVec.push_back(muonsLVec.at(im)); isomuonsIdxVec.push_back(im); }
      }
// If there is a reco muon with pt and eta larger than requirement, use the reco muon
// This way we get consistent results for wherever a reco muon is ready to be used
      if( nMuons == 1 ){
         if( isomuonsLVec.at(0).Pt() >= TauResponse::ptMin() && fabs(isomuonsLVec.at(0).Eta()) <= TauResponse::etaMax() ){
            cleanJetVec = tr.getVec<TLorentzVector>("cleanJetVec");
            cleanJetBtag = tr.getVec<double>("cleanJetBTag");
            useRecoMuon = true;
            const std::vector<int>& rejectJetIdx_formuVec = tr.getVec<int>("rejectJetIdx_formuVec");
            muJetIdx = rejectJetIdx_formuVec.at(isomuonsIdxVec.at(0));
         }
      }
// If no reco muon in the event, use the gen muon instead
      if( !useRecoMuon ){
         const float deltaRMax = 0.2;
         utils::findMatchedObject(muJetIdx, pergenmuLVec, jetsLVec, deltaRMax);
         for(int jetIdx = 0; jetIdx < jetsLVec.size(); ++jetIdx) { // Loop over reco jets
	// Skip this jet if it is the muon
  	   if( jetIdx == muJetIdx ) continue;
	// Calculate NJet
	   cleanJetVec.push_back(jetsLVec.at(jetIdx));
	   cleanJetBtag.push_back(recoJetsBtag_0.at(jetIdx));
         }
      }

      TLorentzVector usedmuLVec;
      if( useRecoMuon ) usedmuLVec = isomuonsLVec.at(0);
      else usedmuLVec = pergenmuLVec;

// Determine if the gen muon passes the kinematic cuts or not
      bool passKinCuts = false;
      if(pergenmuLVec.Pt()>TauResponse::ptMin() && fabs(pergenmuLVec.Eta())<TauResponse::etaMax()) passKinCuts = true;

      TLorentzVector cleanmetLVec; cleanmetLVec.SetVectM( (metLVec + usedmuLVec).Vect(), 0);

      // Get random number from tau-response template
      const double scale = tauResp.getRandom(usedmuLVec.Pt());
      // Scale muon pt with tau response --> simulate tau jet pt
      const double simTauJetPt = scale * usedmuLVec.Pt();
      const double simTauJetE = scale * usedmuLVec.E();

      TLorentzVector tauJetLVec; tauJetLVec.SetPtEtaPhiE(simTauJetPt, usedmuLVec.Eta(), usedmuLVec.Phi(), simTauJetE);
// See comments for similar code in Closure.cc
      double oriJetCSVS = 0;
      if( muJetIdx != -1 ) oriJetCSVS = recoJetsBtag_0.at(muJetIdx);

      vector<TLorentzVector> combJetVec;
      vector<double> combJetsBtag;

// See comments for similar code in Closure.cc
      bool includeTauJet = false;
      for(unsigned int ij=0; ij<cleanJetVec.size(); ij++){
         if( tauJetLVec.Pt() > cleanJetVec.at(ij).Pt() && !includeTauJet ){
            combJetVec.push_back(tauJetLVec); combJetsBtag.push_back(oriJetCSVS);
            includeTauJet = true;
         }
         combJetVec.push_back(cleanJetVec.at(ij)); combJetsBtag.push_back(cleanJetBtag.at(ij));
      }
      // it's possible that the tau jet is the least energetic jet so that it's not added into the combNJetVec during the loop
      if( !includeTauJet ){ combJetVec.push_back(tauJetLVec); combJetsBtag.push_back(oriJetCSVS); }

      TLorentzVector combmetLVec;  combmetLVec.SetVectM( (cleanmetLVec - tauJetLVec).Vect(), 0 );

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
      std::vector<double> * deltaPhiVec = new std::vector<double>();
      (*deltaPhiVec) = AnaFunctions::calcDPhi(combJetVec, combmetPhi, 3, AnaConsts::dphiArr);
      bool passdeltaPhi = true;
      if( deltaPhiVec->at(0) < AnaConsts::dPhi0_CUT || deltaPhiVec->at(1) < AnaConsts::dPhi1_CUT || deltaPhiVec->at(2) < AnaConsts::dPhi2_CUT){
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
      std::vector<TLorentzVector> *jetsLVec_forTagger_acc = new std::vector<TLorentzVector>(); std::vector<double> *recoJetsBtag_forTagger_acc = new std::vector<double>();
      AnaFunctions::prepareJetsForTagger(combJetVec, combJetsBtag, (*jetsLVec_forTagger_acc), (*recoJetsBtag_forTagger_acc));
      
    //Apply Top tagger
      int bestTopJetIdx_acc = -1;
      bool remainPassCSVS_acc = false;
      int pickedRemainingCombfatJetIdx_acc = -1;
      double bestTopJetMass_acc = -1;

      if(comb30_acc >= AnaConsts::nJetsSel ){
	type3Ptr->processEvent((*jetsLVec_forTagger_acc), (*recoJetsBtag_forTagger_acc), combmetLVec);
	bestTopJetIdx_acc = type3Ptr->bestTopJetIdx;
	remainPassCSVS_acc = type3Ptr->remainPassCSVS;
	pickedRemainingCombfatJetIdx_acc = type3Ptr->pickedRemainingCombfatJetIdx;
	if( bestTopJetIdx_acc != -1 ) bestTopJetMass_acc = type3Ptr->bestTopJetLVec.M();
      }
      bool passTopTagger = true;
      if( bestTopJetIdx_acc == -1 ){passTopTagger = false; }
      if( ! remainPassCSVS_acc ){passTopTagger = false; }
      if( pickedRemainingCombfatJetIdx_acc == -1 && jetsLVec_forTagger_acc->size()>=6 ){passTopTagger = false; }
      if( ! (bestTopJetMass_acc > AnaConsts::lowTopCut_ && bestTopJetMass_acc < AnaConsts::highTopCut_ ) ){ passTopTagger = false; }

      if(!passTopTagger) continue;

      genmu++;
      if( passKinCuts ) accmu++;
    }
  }
  cout<<"ToTal Event: "<<k<<endl;
  cout<<"Gen Mu: "<<genmu<<" "<<"Acc Mu: "<<accmu<<endl;
  cout<<"Acceptance: "<<(float)accmu/genmu<<endl;
  return 0;
}
