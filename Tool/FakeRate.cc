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
#include "TStopwatch.h"
#include "TString.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "FakeRate.h"
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

using namespace std;

void passBaselineFunc(NTupleReader &tr)
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
  /*
  //Pass deltaPhi?
  bool passdPhis = true;
  if( dPhiVec->at(0) < AnaConsts::dPhi0_CUT || dPhiVec->at(1) < AnaConsts::dPhi1_CUT || dPhiVec->at(2) < AnaConsts::dPhi2_CUT ){ passBaseline = false; passdPhis = false; }
  //Pass the baseline MET requirement?
  bool passMET = true;
  if( tr.getVar<double>("met") < AnaConsts::defaultMETcut ){ passBaseline = false; passMET = false; }
  
//Pass number of b-tagged jets?
  bool passBJets = true;
  if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cntCSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cntCSVS < AnaConsts::high_nJetsSelBtagged ) ) ){ passBaseline = false; passBJets = false; }
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
  bool passTagger = type3Ptr->passNewTaggerReq();
  //bestTopJetIdx != -1 means at least 1 top candidate!
  if( !passTagger ) passBaseline = false;
  */
    //register new var
  tr.registerDerivedVar("passBaseline", passBaseline);
  tr.registerDerivedVar("passLeptVeto", passLeptVeto);
  tr.registerDerivedVar("passnJets", passnJets);
}

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 3)
    {
      std::cerr <<"Please give 3 arguments " << "inputList " << " " << "outputFileName" << " " << "lumiscale" << std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./FakeRate List1_ttbar.txt HadTau_TaujetBjetFakerate.root scalevalue" << std::endl;
      return -1;
    }
  const char *inputFileList = argv[1];
  const char *outFileName = argv[2];
  const char *lumiscale = argv[3];
  
  const double scale = std::atof(lumiscale);

  TChain *fChain = new TChain("stopTreeMaker/AUX");
  if(!FillChain(fChain, inputFileList))
    {
      std::cerr << "Cannot get the tree " << std::endl;
    }
  NTupleReader tr(fChain);
  AnaFunctions::prepareTopTagger();
  tr.registerFunction(&passBaselineFunc);

  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(outFileName);

  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr.getNEntries();
  std::cout<<"\nentries : "<<entries<<std::endl;
  int k = 0;
  while(tr.getNextEvent()){
    k++;
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
  bool passBaseline = tr.getVar<bool>("passBaseline");
  bool passLeptVeto = tr.getVar<bool>("passLeptVeto");
  bool passnJets = tr.getVar<bool>("passnJets");

  // Select only events where the W decayed into a hadronically decaying tau
  if(W_tau_prongsVec.size()==0 ) continue;
  if(!passLeptVeto) continue;
  if(!passnJets) continue;

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

  for(unsigned int it=0; it< genvisiblehadtauLVec.size();it++){

    std::vector<TLorentzVector> cleanJetVec;
    std::vector<double> cleanJetBtag;
// Kinematic variables of generator-level tau
    TLorentzVector genTauLVec = genvisiblehadtauLVec.at(it);

  // Do the matching
    int tauJetIdx = -1; // Will store the index of the jet matched to the tau
    //const float deltaRMax = genTauLVec.Pt() < 50. ? 0.2 : 0.1; // Increase deltaRMax at low pt to maintain high-enought matching efficiency
    const float deltaRMax = deltaRmax(genTauLVec.Pt()); // Increase deltaRMax at low pt to maintain high-enought matching efficiency
    if( !utils::findTauMatchedJet(tauJetIdx,genTauLVec,jetsLVec,deltaRMax) ) continue;

    myBaseHistgram.htauBjetEta_den->Fill(jetsLVec.at(tauJetIdx).Eta(), scale);
    myBaseHistgram.htauBjetPt_den->Fill(jetsLVec.at(tauJetIdx).Pt(), scale);

      // select the had-tau jet
      cleanJetVec.push_back(jetsLVec.at(tauJetIdx));
      cleanJetBtag.push_back(recoJetsBtag_0.at(tauJetIdx));

    int cnt1CSVS = AnaFunctions::countCSVS(cleanJetVec, cleanJetBtag, AnaConsts::cutCSVS, AnaConsts::bTagArr);
    //if( (AnaConsts::low_nJetsSelBtagged == -1 || cnt1CSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cnt1CSVS < AnaConsts::high_nJetsSelBtagged ) ){
    if(cnt1CSVS>=1){
      myBaseHistgram.htauBjetEta_num->Fill(jetsLVec.at(tauJetIdx).Eta(), scale);
      myBaseHistgram.htauBjetPt_num->Fill(jetsLVec.at(tauJetIdx).Pt(), scale);
    }
  }//finish gen hadtau loop
  
  }//finish event loop
  std::cout<<"final event: "<<k<<std::endl;
  /*  myBaseHistgram.hEff_Pt = (TH1D*)myBaseHistgram.htauBjetPt_num->Clone("Ratio");
  myBaseHistgram.hEff_Pt->GetYaxis()->SetRangeUser(0. , 0.7);
  myBaseHistgram.hEff_Pt->Divide(myBaseHistgram.htauBjetPt_den);
  myBaseHistgram.hEff_Eta = (TH1D*)myBaseHistgram.htauBjetEta_num->Clone("Ratio");
  myBaseHistgram.hEff_Eta->Divide(myBaseHistgram.htauBjetEta_den);*/
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();

  /*  cout<<myBaseHistgram.hEff_Pt->GetBinContent(1)<<endl;
  cout<<myBaseHistgram.hEff_Pt->GetBinContent(2)<<endl;
  cout<<myBaseHistgram.hEff_Pt->GetBinContent(3)<<endl;
  cout<<myBaseHistgram.hEff_Pt->GetBinContent(4)<<endl;
  cout<<myBaseHistgram.hEff_Pt->GetBinContent(5)<<endl;
  cout<<myBaseHistgram.hEff_Pt->GetBinContent(6)<<endl;
  cout<<myBaseHistgram.hEff_Pt->GetBinContent(7)<<endl;
  cout<<myBaseHistgram.hEff_Pt->GetBinContent(8)<<endl;
  */

  return 0;
}
double deltaRmax(double pt){
  double rmax = 0.25;
  if(pt>30) rmax = 0.2;
  if(pt>50) rmax = 0.15;
  if(pt>100) rmax = 0.1;
  return rmax;
}
