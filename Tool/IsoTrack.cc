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
#include "IsoTrack.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TChain.h"

using namespace std;

void passBaselineFunc1(NTupleReader &tr)
{
  bool passBaseline = true;
  //Form TLorentzVector of MET
  TLorentzVector metLVec; metLVec.SetPtEtaPhiM(tr.getVar<double>("met"), 0, tr.getVar<double>("metphi"), 0);
  //Calculate number of leptons
  int nMuons = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsMiniIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsMiniIsoArr);
  int nElectrons = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesMiniIso"), tr.getVec<double>("elesMtw"), tr.getVec<unsigned int>("elesisEB"), AnaConsts::elesArr);
  // int nIsoTrks = AnaFunctions::countIsoTrks(tr.getVec<TLorentzVector>("loose_isoTrksLVec"), tr.getVec<double>("loose_isoTrks_iso"), tr.getVec<double>("loose_isoTrks_mtw"), AnaConsts::isoTrksArr);
  //Calculate number of jets and b-tagged jets
  int cntCSVS = AnaFunctions::countCSVS(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), AnaConsts::cutCSVS, AnaConsts::bTagArr);
  int cntNJetsPt50Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt50Eta24Arr);
  int cntNJetsPt30Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Eta24Arr);
  int cntNJetsPt30 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Arr);
  //Calculate deltaPhi
  std::vector<double> * dPhiVec = new std::vector<double>();
  (*dPhiVec) = AnaFunctions::calcDPhi(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVar<double>("metphi"), 3, AnaConsts::dphiArr);
  std::vector<TLorentzVector> *jetsLVec_forTagger = new std::vector<TLorentzVector>(); std::vector<double> *recoJetsBtag_forTagger = new std::vector<double>();
  AnaFunctions::prepareJetsForTagger(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), (*jetsLVec_forTagger), (*recoJetsBtag_forTagger));
  //Pass lepton veto?
  bool passLeptVeto = true;
  if( nMuons != AnaConsts::nMuonsSel ){ passBaseline = false; passLeptVeto = false; }
  if( nElectrons != AnaConsts::nElectronsSel ){ passBaseline = false; passLeptVeto = false; }
  //Pass number of jets?
  bool passnJets = true;
  if( cntNJetsPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 ){ passBaseline = false; passnJets = false;}
  if( cntNJetsPt30Eta24 < AnaConsts::nJetsSelPt30Eta24 ){ passBaseline = false; passnJets = false;}
  //Pass deltaPhi?
  bool passdPhis = true;
  if( dPhiVec->at(0) < AnaConsts::dPhi0_CUT || dPhiVec->at(1) < AnaConsts::dPhi1_CUT || dPhiVec->at(2) < AnaConsts::dPhi2_CUT ){ passBaseline = false; passdPhis = false; }
  //Pass number of b-tagged jets?
  bool passBJets = true;
  if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cntCSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cntCSVS < AnaConsts::high_nJetsSelBtagged ) ) ){ passBaseline = false; passBJets = false; }
  //Pass the baseline MET requirement?
  bool passMET = true;
  if( tr.getVar<double>("met") < AnaConsts::defaultMETcut ){ passBaseline = false; passMET = false; }
  //Calculate top tagger related variables.
  //Note that to save speed, only do the calculation after previous base line requirements.
  int nTopCandSortedCnt = -1;
  double MT2 = -1;
  double mTcomb = -1;
  if( passnJets && cntNJetsPt30 >= AnaConsts::nJetsSel ){
    type3Ptr->processEvent((*jetsLVec_forTagger), (*recoJetsBtag_forTagger), metLVec);
    nTopCandSortedCnt = type3Ptr->nTopCandSortedCnt;
    MT2 = type3Ptr->best_had_brJet_MT2;
    mTcomb = type3Ptr->best_had_brJet_mTcomb;
  }
  //Pass top tagger requirement?
  bool passTagger = type3Ptr->passNewTaggerReq() && nTopCandSortedCnt >= AnaConsts::low_nTopCandSortedSel;
  //bestTopJetIdx != -1 means at least 1 top candidate!
  if( !passTagger ) passBaseline = false;
    //register new var
  tr.registerDerivedVar("cntNJetsPt30Eta24", cntNJetsPt30Eta24);
  tr.registerDerivedVar("passBaseline", passBaseline);
  tr.registerDerivedVar("cntCSVS", cntCSVS);
  tr.registerDerivedVar("nTopCandSortedCnt", nTopCandSortedCnt);
  tr.registerDerivedVar("MT2", MT2);
}

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 2)
    {
      std::cerr <<"Please give 2 arguments " << "inputList " << " " << "outputFileName" << std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./IsoTrack List2_ttbar.txt HadTau_IsoTrack.root" << std::endl;
      return -1;
    }
  const char *inputFileList = argv[1];
  const char *outFileName = argv[2];
  TChain *fChain = new TChain("stopTreeMaker/AUX");
  if(!FillChain(fChain, inputFileList))
    {
      std::cerr << "Cannot get the tree " << std::endl;
    }
  NTupleReader tr(fChain);
  AnaFunctions::prepareTopTagger();
  tr.registerFunction(&passBaselineFunc1);

  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(outFileName);

  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr.getNEntries();
  std::cout<<"\nentries : "<<entries<<std::endl;

  while(tr.getNextEvent()){

    //Add print out of the progress of looping
    if( tr.getEvtNum()-1 == 0 || tr.getEvtNum() == entries || (tr.getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr.getEvtNum()-1<<"th event ..."<<std::endl;

  const vector<TLorentzVector> &genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
  const vector<int> &genDecayIdxVec = tr.getVec<int>("genDecayIdxVec");
  const vector<int> &genDecayMomIdxVec = tr.getVec<int>("genDecayMomIdxVec");
  const vector<int> &genDecayPdgIdVec = tr.getVec<int>("genDecayPdgIdVec");
  const vector<int> &W_emuVec = tr.getVec<int>("W_emuVec");
  const vector<int> &W_tau_emuVec = tr.getVec<int>("W_tau_emuVec");
  const vector<int> &W_tau_nuVec = tr.getVec<int>("W_tau_nuVec");
  const vector<int> &W_tau_prongsVec = tr.getVec<int>("W_tau_prongsVec");
  const vector<TLorentzVector> &jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
  const vector<double> &recoJetsBtag_0 = tr.getVec<double>("recoJetsBtag_0");
  const vector<int> &looseisoTrksMatchedJetIdx = tr.getVec<int>("looseisoTrksMatchedJetIdx");
  const vector<TLorentzVector> &loose_isoTrksLVec = tr.getVec<TLorentzVector>("loose_isoTrksLVec");
  const vector<double> &loose_isoTrks_iso = tr.getVec<double>("loose_isoTrks_iso");
  const vector<double> &loose_isoTrks_mtw = tr.getVec<double>("loose_isoTrks_mtw");
  const vector<int> &loose_isoTrks_pdgId = tr.getVec<int>("loose_isoTrks_pdgId");

  const double met = tr.getVar<double>("met");
  const int njets = tr.getVar<int>("cntNJetsPt30Eta24");
  const int nbjets = tr.getVar<int>("cntCSVS");
  const int nTops = tr.getVar<int>("nTopCandSortedCnt");
  const double MT2 = tr.getVar<double>("MT2");
  bool passBaseline = tr.getVar<bool>("passBaseline");

  // Select only events where the W decayed into a hadronically decaying tau
  if(W_tau_prongsVec.size()==0 ) continue;
  if(!passBaseline) continue;

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
	TLorentzVector gennuLVec;
	TLorentzVector genhadtauLVec = genDecayLVec.at(ig);
	for(int n=0; n<W_tau_nuVec.size();n++){
	  int nIdx = W_tau_nuVec.at(n);
	  if( genDecayMomIdxVec.at(nIdx) == genDecayIdxVec.at(ig)) gennuLVec = genDecayLVec.at(nIdx);
	}
	genvisiblehadtauLVec.push_back(genhadtauLVec-gennuLVec);
      }
    }
  }

  bool isIsotrack = false;
  if(looseisoTrksMatchedJetIdx.size()!=loose_isoTrksLVec.size())cout<<"Error: isotrack vetor size mismatch"<<endl;
  int nTotal = 0, nHadtau = 0;
  std::vector<int> IsotrkMatchedHadtau(genvisiblehadtauLVec.size(), 0);

  for(unsigned int it=0; it< looseisoTrksMatchedJetIdx.size();it++){
    if(!passIsoTrks(loose_isoTrksLVec[it], loose_isoTrks_iso[it], loose_isoTrks_mtw[it], loose_isoTrks_pdgId[it])) continue;
    nTotal++;    
// Do the matching
      const int isoJetIdx = looseisoTrksMatchedJetIdx[it]; // Will store the index of the jet matched to the isotrack
      if(findTauMatchedisoJet(isoJetIdx,genvisiblehadtauLVec,jetsLVec, IsotrkMatchedHadtau)){
	if(IsotrkMatchedHadtau[0]) nHadtau++; // considering only first had tau
      }
    }//finish isotrack loop
  
  if((nTotal-nHadtau)!=0)continue;//isotrackveto on the remaing part (exclude hadtau)
  if(nHadtau)isIsotrack = true;

  myBaseHistgram.hSB_den->Fill(65);
    if( isIsotrack ){
      myBaseHistgram.hSB_num->Fill(65);
    }

    // iSR: this should be determined by search region requirement
    int iSR = find_Binning_Index(nbjets, nTops, MT2, met);
    if(iSR!=-1) {
      if( isIsotrack ){
	myBaseHistgram.hSB_num->Fill(iSR);
      }
      myBaseHistgram.hSB_den->Fill(iSR);
    }
    if(nbjets>=1 && nTops>=1 && met>=200 && MT2<200){
      myBaseHistgram.hSB_den->Fill(64);
      if( isIsotrack ){
	myBaseHistgram.hSB_num->Fill(64);
      }
    }

    //histogram
    myBaseHistgram.hmet_den->Fill(met);
    myBaseHistgram.hMT2_den->Fill(MT2);
    myBaseHistgram.hNbjet_den->Fill(nbjets);
    myBaseHistgram.hNtop_den->Fill(nTops);
    myBaseHistgram.hNjet_den->Fill(njets);
    myBaseHistgram.hNjetNbjet_den->Fill(njets, nbjets);
    if(isIsotrack){
      myBaseHistgram.hmet_num->Fill(met);
      myBaseHistgram.hMT2_num->Fill(MT2);
      myBaseHistgram.hNbjet_num->Fill(nbjets);
      myBaseHistgram.hNtop_num->Fill(nTops);
      myBaseHistgram.hNjet_num->Fill(njets);
      myBaseHistgram.hNjetNbjet_num->Fill(njets, nbjets);
    }

  
  
  }//finish event loop

  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();
  return 0;
}
double deltaRmax(double pt){
  double rmax = 0.25;
  if(pt>30) rmax = 0.2;
  if(pt>50) rmax = 0.15;
  if(pt>100) rmax = 0.1;
  return rmax;
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
bool passIsoTrks(const TLorentzVector isoTrksLVec, const double isoTrksIso, const double isoTrksMtw, const int isoTrkspdgId){
  bool passIsoTrks = false;
    if( std::abs(isoTrkspdgId) == 11 || std::abs(isoTrkspdgId) == 13 ){
      if( AnaFunctions::passIsoTrk(isoTrksLVec, isoTrksIso, isoTrksMtw, AnaConsts::isoLepTrksArr ) ) passIsoTrks = true;
    }
    if( std::abs(isoTrkspdgId) == 211 ){
      if(AnaFunctions::passIsoTrk(isoTrksLVec, isoTrksIso, isoTrksMtw, AnaConsts::isoHadTrksArr ) ) passIsoTrks = true;
    }
  return passIsoTrks;
}
