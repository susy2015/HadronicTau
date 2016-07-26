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
#include "Efficiency.h"

using namespace std;


// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
    {
      std::cerr <<"Please give 4 arguments "<<"SubsampleName"<<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./IsoTrack TTbarInc 1000 0 1" << std::endl;
      return -1;
    }
  const char *subsamplename = argv[1];
  const char *Maxevent = argv[2];
  const  char *Stratfile = argv[3];
  const  char *Filerun = argv[4];
  const  int startfile = std::atoi(Stratfile);
  const int filerun = std::atoi(Filerun);
  TChain *fChain = 0;
  const string condorSpec = argc==6 ? argv[5]: "";

  if(!FillChain(fChain, subsamplename, condorSpec, startfile, filerun))
      {
	std::cerr << "Cannot get the tree " << std::endl;
      }

  const int maxevent = std::atoi(Maxevent);

  std::string spec = "IsoTrack";
  IsoTrackBaselineVessel = new BaselineVessel(spec);
  AnaFunctions::prepareForNtupleReader();
  AnaFunctions::prepareTopTagger();
  NTupleReader *tr =0;
  tr = new NTupleReader(fChain);
  tr->registerFunction(&passBaselineFuncIsoTrack);
  //Searchbin                                                                                                                                                                             
  SearchBins SB("SB_59_2016");
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(subsamplename, startfile);

  int nisotrk=0, nunmatch=0;

  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr->getNEntries();
  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<Lumiscale<<std::endl;
  cout<<"maxevent: "<<maxevent<<endl;
  while(tr->getNextEvent()){
    if(maxevent>=0 && tr->getEvtNum() > maxevent ) break;
    //Add print out of the progress of looping
    if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;

  const vector<TLorentzVector> &genDecayLVec = tr->getVec<TLorentzVector>("genDecayLVec");
  const vector<int> &genDecayIdxVec = tr->getVec<int>("genDecayIdxVec");
  const vector<int> &genDecayMomIdxVec = tr->getVec<int>("genDecayMomIdxVec");
  const vector<int> &genDecayPdgIdVec = tr->getVec<int>("genDecayPdgIdVec");
  const vector<int> &W_emuVec = tr->getVec<int>("W_emuVec");
  const vector<int> &W_tau_emuVec = tr->getVec<int>("W_tau_emuVec");
  const vector<int> &W_tau_nuVec = tr->getVec<int>("W_tau_nuVec");
  const vector<int> &W_tau_prongsVec = tr->getVec<int>("W_tau_prongsVec");
  const vector<TLorentzVector> &jetsLVec = tr->getVec<TLorentzVector>("jetsLVec");
  const vector<double> &recoJetsBtag_0 = tr->getVec<double>("recoJetsBtag_0");
  const vector<int> &looseisoTrksMatchedJetIdx = tr->getVec<int>("looseisoTrksMatchedJetIdx");
  const vector<TLorentzVector> &loose_isoTrksLVec = tr->getVec<TLorentzVector>("loose_isoTrksLVec");
  const vector<double> &loose_isoTrks_iso = tr->getVec<double>("loose_isoTrks_iso");
  const vector<double> &loose_isoTrks_mtw = tr->getVec<double>("loose_isoTrks_mtw");
  const vector<int> &loose_isoTrks_pdgId = tr->getVec<int>("loose_isoTrks_pdgId");

  const double met = tr->getVar<double>("met");
  const double ht = tr->getVar<double>("HT"+spec);
  const int njets = tr->getVar<int>("cntNJetsPt30Eta24"+spec);
  const int nbjets = tr->getVar<int>("cntCSVS"+spec);
  const int nTops = tr->getVar<int>("nTopCandSortedCnt"+spec);
  const double MT2 = tr->getVar<double>("best_had_brJet_MT2"+spec);
  bool passBaseline = tr->getVar<bool>("passBaseline"+spec);
  bool passLeptVeto = tr->getVar<bool>("passLeptVeto"+spec);
  bool passMuonVeto = tr->getVar<bool>("passMuonVeto"+spec);
  bool passEleVeto = tr->getVar<bool>("passEleVeto"+spec);
  bool passnJets = tr->getVar<bool>("passnJets"+spec);
  bool passdPhis = tr->getVar<bool>("passdPhis"+spec);
  bool passMET =  tr->getVar<bool>("passMET"+spec);
  bool passBJets = tr->getVar<bool>("passBJets"+spec);
  bool passTagger = tr->getVar<bool>("passTagger"+spec);
  bool passHT =  tr->getVar<bool>("passHT"+spec);
  bool passMT2 = tr->getVar<bool>("passMT2" + spec);
  bool passNoiseEventFilter = tr->getVar<bool>("passNoiseEventFilter"+spec);

  // Select only events where the W decayed into a hadronically decaying tau
  if(W_tau_prongsVec.size()==0 ) continue;
  bool IsotrkBaseline = passMuonVeto && passEleVeto && passnJets && passdPhis && passMET && passBJets && passTagger && passHT;
  if(!IsotrkBaseline) continue;
  //Event Filter                                                                                                                           
  if(!passNoiseEventFilter) continue;

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
  bool oneprng = false, threeprng = false;
  if(W_tau_prongsVec.size()%2==0 && W_tau_prongsVec.size()!=0) cout<<"genvisiblehadtauLVec: "<<genvisiblehadtauLVec.size()<<endl;
  if(W_tau_prongsVec.size()==1) oneprng =true;
  if(W_tau_prongsVec.size()==3) threeprng =true;
  bool isIsotrack = false;
  if(looseisoTrksMatchedJetIdx.size()!=loose_isoTrksLVec.size())cout<<"Error: isotrack vetor size mismatch"<<endl;
  int nTotal = 0, nHadtau = 0;
  std::vector<int> IsotrkMatchedHadtau(genvisiblehadtauLVec.size(), 0);
  for(unsigned int it=0; it< looseisoTrksMatchedJetIdx.size();it++){
    if(!passIsoTrks(loose_isoTrksLVec[it], loose_isoTrks_iso[it], loose_isoTrks_mtw[it], loose_isoTrks_pdgId[it])) continue;
    nTotal++; nisotrk++;    
    // Do the matching
    const int isoJetIdx = looseisoTrksMatchedJetIdx[it]; // Will store the index of the jet matched to the isotrack
    /*    if(isoJetIdx == -1){
      cout<<" isoJetIdx: "<<isoJetIdx<<"  track pT: "<<loose_isoTrksLVec[it].Pt()<<"  track phi: "<<loose_isoTrksLVec[it].Phi()<<"  track Eta: "<<loose_isoTrksLVec[it].Eta()<<"  track id: "<<loose_isoTrks_pdgId[it]<<endl;
      nunmatch++;
      }*/
    if(findTauMatchedisoJet(isoJetIdx,genvisiblehadtauLVec,jetsLVec, IsotrkMatchedHadtau)){
      if(IsotrkMatchedHadtau[0]) nHadtau++; // considering only first had tau
    }
  }//finish isotrack loop
  
  if((nTotal-nHadtau)!=0)continue;//isotrackveto on the remaing part (exclude hadtau)
  if(nHadtau)isIsotrack = true;
  
  //  const double isotrkEff = 1./Efficiency::isotrkeffMix_NjetNbjet(Efficiency::Njetbin(njets), Efficiency::NBjetbin(nbjets));//correction for isotrackveto eff.
  
  // iSR: this should be determined by search region requirement
  int iSR = SB.find_Binning_Index(nbjets, nTops, MT2, met);
  if(iSR!=-1) {
    if( isIsotrack ){
      myBaseHistgram.hSB_num->Fill(iSR, Lumiscale);
      if(oneprng)  myBaseHistgram.hSB_num1pr->Fill(iSR, Lumiscale);
      if(threeprng)  myBaseHistgram.hSB_num3pr->Fill(iSR, Lumiscale);
    }
    myBaseHistgram.hSB_den->Fill(iSR, Lumiscale);
      if(oneprng)  myBaseHistgram.hSB_den1pr->Fill(iSR, Lumiscale);
      if(threeprng)  myBaseHistgram.hSB_den3pr->Fill(iSR, Lumiscale);
  }
  
  //histogram
  if(passMT2){
    myBaseHistgram.hmet_den->Fill(met, Lumiscale);
    myBaseHistgram.hMT2_den->Fill(MT2, Lumiscale);
    myBaseHistgram.hNbjet_den->Fill(nbjets, Lumiscale);
    myBaseHistgram.hNtop_den->Fill(nTops, Lumiscale);
    FillInt(myBaseHistgram.hNjet_den, njets, Lumiscale);
    FillDouble(myBaseHistgram.hht_den, ht, Lumiscale);
    Fill2DIsoTrk(myBaseHistgram.hNjetNbjet_den, njets, nbjets, Lumiscale);
    if(isIsotrack){
      myBaseHistgram.hmet_num->Fill(met, Lumiscale);
      myBaseHistgram.hMT2_num->Fill(MT2, Lumiscale);
      myBaseHistgram.hNbjet_num->Fill(nbjets, Lumiscale);
      myBaseHistgram.hNtop_num->Fill(nTops, Lumiscale);
      FillInt(myBaseHistgram.hNjet_num, njets, Lumiscale);
      FillDouble(myBaseHistgram.hht_num, ht, Lumiscale);
      Fill2DIsoTrk(myBaseHistgram.hNjetNbjet_num, njets, nbjets, Lumiscale);
    }
  }
  
  
  
  }//finish event loop
  
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();

  //  cout<<"Event with Isotrk: "<<nisotrk<<"  Event with unmatched track: "<<nunmatch<<endl;
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
