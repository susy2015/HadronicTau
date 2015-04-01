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
  //  if( cntNJetsPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 ){ passBaseline = false; passBaseline_nolepveto = false; passnJets = false;}
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
  TauResponse tauResp(respTempl);

  int genmu=0, accmu=0;

  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;

  // Loop over the events (tree entries)
  int k = 0;
  while(tr.getNextEvent()){
    k++;

    vector<TLorentzVector> genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
    vector<int> genDecayIdxVec = tr.getVec<int>("genDecayIdxVec");
    vector<int> genDecayPdgIdVec = tr.getVec<int>("genDecayPdgIdVec");
    vector<TLorentzVector> jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
    double met=tr.getVar<double>("met");
    double metphi=tr.getVar<double>("metphi");

    int nmu=0, nele=0;
    vector<double> genmueta, genmuphi, genmupt, genmum;
    genmueta.clear();
    genmuphi.clear();
    genmupt.clear();
    genmum.clear();

    for(unsigned ig=0; ig<genDecayLVec.size(); ig++){
      int pdgId = genDecayPdgIdVec.at(ig);

      if(fabs(pdgId)==13){
	nmu++;
	genmueta.push_back(genDecayLVec.at(ig).Eta());
	genmuphi.push_back(genDecayLVec.at(ig).Phi());
	genmupt.push_back(genDecayLVec.at(ig).Pt());
	genmum.push_back(genDecayLVec.at(ig).M());
      }
      if(fabs(pdgId)==11) nele++;
    }

    if(nmu ==1 && nele==0){//control sample
      
      const double muEta = genmueta.at(0);
      const double muPhi = genmuphi.at(0);
      const double muPt = genmupt.at(0);
      const double muM = genmum.at(0);

      vector<double> jetseta, jetsphi;
      jetseta.clear();
      jetsphi.clear();
      for(unsigned ij=0; ij<jetsLVec.size(); ij++){
	double eta1 = jetsLVec.at(ij).Eta();
	double phi1 = jetsLVec.at(ij).Phi();
	jetseta.push_back(eta1);
	jetsphi.push_back(phi1);
      }
      int muJetIdx = -1;
      const float deltaRMax = 0.2;
      const unsigned nObj = jetsLVec.size();
      if( !utils::findMatchedObject(muJetIdx,muEta,muPhi,jetseta,jetsphi,nObj,deltaRMax) ) continue;
      vector<TLorentzVector> cleanJetVec;
      cleanJetVec.clear();
      for(int jetIdx = 0; jetIdx < jetsLVec.size(); ++jetIdx) { // Loop over reco jets
	// Skip this jet if it is the muon
	if( jetIdx == muJetIdx ) continue;
	// Calculate NJet
	double Pt = jetsLVec.at(jetIdx).Pt();
	double Eta = jetsLVec.at(jetIdx).Eta();
	double Phi = jetsLVec.at(jetIdx).Phi();
	double M = jetsLVec.at(jetIdx).M();
	TLorentzVector cleanLVec; cleanLVec.SetPtEtaPhiM(Pt, Eta, Phi, M);
	cleanJetVec.push_back(cleanLVec);
      }

      // Get random number from tau-response template
      const double scale = tauResp.getRandom(muPt);
      // Scale muon pt with tau response --> simulate tau jet pt
      const double simJetPt = scale * muPt;

      //compute combined jet
      vector<TLorentzVector> simJetVec;
      simJetVec.clear();
      TLorentzVector simLVec; simLVec.SetPtEtaPhiM(simJetPt, muEta, muPhi, muM);
      simJetVec.push_back(simLVec);

      vector<TLorentzVector> combJetVec;
      combJetVec.clear();
      combJetVec = combjet(cleanJetVec, simJetVec);

      //recompute met
      double simmet = met + muPt - simJetPt;
      bool passmet = true;
      if(simmet<AnaConsts::defaultMETcut) passmet = false; 

      //recompute jet
      int nJet = AnaFunctions::countJets(combJetVec, AnaConsts::pt30Eta24Arr);
      bool passnjets = true;
      if(nJet<AnaConsts::nJetsSelPt30Eta24)passnjets = false;
      //      if(nJet<3)passnjets = false;

      //recompute deltaphi
      std::vector<double> * deltaPhiVec = new std::vector<double>();
      (*deltaPhiVec) = AnaFunctions::calcDPhi(combJetVec, metphi, 3, AnaConsts::dphiArr);
      bool passdeltaPhi = true;
      if( deltaPhiVec->at(0) < AnaConsts::dPhi0_CUT || deltaPhiVec->at(1) < AnaConsts::dPhi1_CUT || deltaPhiVec->at(2) < AnaConsts::dPhi2_CUT){
	passdeltaPhi = false;
      }
      //recompute bjet
      int cnt1CSVS = AnaFunctions::countCSVS(combJetVec, tr.getVec<double>("recoJetsBtag_0"), AnaConsts::cutCSVS, AnaConsts::bTagArr);
      bool passbJets = true;
      if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cnt1CSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cnt1CSVS < AnaConsts::high_nJetsSelBtagged ) ) ){
	passbJets = false;
      }

      //top tagger input
      TLorentzVector metLVec_acc; metLVec_acc.SetPtEtaPhiM(met, 0, metphi, 0);
      int comb30_acc = AnaFunctions::countJets(combJetVec, AnaConsts::pt30Arr);
      std::vector<TLorentzVector> *jetsLVec_forTagger_acc = new std::vector<TLorentzVector>(); std::vector<double> *recoJetsBtag_forTagger_acc = new std::vector<double>();
      AnaFunctions::prepareJetsForTagger(combJetVec, tr.getVec<double>("recoJetsBtag_0"), (*jetsLVec_forTagger_acc), (*recoJetsBtag_forTagger_acc));
      int bestTopJetIdx_acc = -1;
      bool remainPassCSVS_acc = false;
      int pickedRemainingCombfatJetIdx_acc = -1;
      double bestTopJetMass_acc = -1;

      if(!passnjets) continue;
      if(!passmet) continue;
      if(!passdeltaPhi) continue;
      if(!passbJets) continue;

      // if(!passBaseline_nolepveto) continue;

    //Apply Top tagger
      if(comb30_acc >= AnaConsts::nJetsSel ){
	type3Ptr->processEvent((*jetsLVec_forTagger_acc), (*recoJetsBtag_forTagger_acc), metLVec_acc);
	bestTopJetIdx_acc = type3Ptr->bestTopJetIdx;
	remainPassCSVS_acc = type3Ptr->remainPassCSVS;
	pickedRemainingCombfatJetIdx_acc = type3Ptr->pickedRemainingCombfatJetIdx;
	if( bestTopJetIdx_acc != -1 ) bestTopJetMass_acc = type3Ptr->bestTopJetLVec.M();
      }
      bool passTopTagger = true;
      //bestTopJetIdx_pre != -1 means at least 1 top candidate!
      if( bestTopJetIdx_acc == -1 ){passTopTagger = false; }
      if( ! remainPassCSVS_acc ){passTopTagger = false; }
      if( pickedRemainingCombfatJetIdx_acc == -1 && jetsLVec_forTagger_acc->size()>=6 ){passTopTagger = false; }
      if( ! (bestTopJetMass_acc > AnaConsts::lowTopCut_ && bestTopJetMass_acc < AnaConsts::highTopCut_ ) ){ passTopTagger = false; }
      
      if(!passTopTagger) continue;

      genmu++;
      if(simJetPt>20 && fabs(muEta)<2.4)accmu++;
    }

  }
  cout<<"ToTal Event: "<<k<<endl;
  cout<<"Gen Mu: "<<genmu<<" "<<"Acc Mu: "<<accmu<<endl;
  cout<<"Acceptance: "<<(float)accmu/genmu<<endl;
  return 0;
}

vector<TLorentzVector> combjet (const vector<TLorentzVector> &seljet, const vector<TLorentzVector> &simjet){
  vector<TLorentzVector> combNJet;
  combNJet.clear();
  unsigned int idx;
  for(unsigned int i=0; i<seljet.size(); i++){
    TLorentzVector comb;
    if(seljet.at(i).Pt()>simjet.at(0).Pt()){
      comb.SetPtEtaPhiM(seljet.at(i).Pt(), seljet.at(i).Eta(), seljet.at(i).Phi(), seljet.at(i).M());
      combNJet.push_back(comb);
    }
    else {
      comb.SetPtEtaPhiM(simjet.at(0).Pt(), simjet.at(0).Eta(), simjet.at(0).Phi(), simjet.at(0).M());
      combNJet.push_back(comb);
      idx=i;
      break;
    }
  }
  if(idx<seljet.size()){
    for(unsigned int j=idx; j<seljet.size(); j++){
      TLorentzVector comb;
      comb.SetPtEtaPhiM(seljet.at(j).Pt(), seljet.at(j).Eta(), seljet.at(j).Phi(), seljet.at(j).M());
      combNJet.push_back(comb);
    }
  }
  if(combNJet.size()==seljet.size()){
    TLorentzVector comb;
    comb.SetPtEtaPhiM(simjet.at(0).Pt(), simjet.at(0).Eta(), simjet.at(0).Phi(), simjet.at(0).M());
    combNJet.push_back(comb);
  }
  return combNJet;
}
