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

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
    {
      std::cerr <<"Please give 4 arguments "<<"SubsampleName"<<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Closure TTbarInc 1000 0 1" << std::endl;
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

  //Use BaselineVessel class for baseline variables and selections                                                                                                                                                                           
  std::string spec = "FakeRate";
  FakeRateBaselineVessel = new BaselineVessel(spec);
  AnaFunctions::prepareForNtupleReader();
  AnaFunctions::prepareTopTagger();
  NTupleReader *tr =0;
  tr = new NTupleReader(fChain);
  tr->registerFunction(&passBaselineFuncFakeRate);

  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(subsamplename, startfile);

  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr->getNEntries();
  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<Lumiscale<<std::endl;
  cout<<"maxevent: "<<maxevent<<endl;
  int k = 0;
  while(tr->getNextEvent()){
    k++;
    if(maxevent>=0 && tr->getEvtNum() > maxevent ) break;
    // Add print out of the progress of looping
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
    bool passBaseline = tr->getVar<bool>("passBaseline"+spec);
    bool passLeptVeto = tr->getVar<bool>("passLeptVeto"+spec);
    bool passnJets = tr->getVar<bool>("passnJets"+spec);
    
    // Select only events where the W decayed into a hadronically decaying tau
    if(W_tau_prongsVec.size()==0 ) continue;
    if(!passLeptVeto) continue;
    if(!passnJets) continue;
    
    std::vector<TLorentzVector>genvisiblehadtauLVec;
    std::vector<TLorentzVector>genbLVec;
    for(unsigned ig=0; ig<genDecayLVec.size(); ig++){
      int pdgId = genDecayPdgIdVec.at(ig);
      if(abs(pdgId)==5){
	genbLVec.push_back(genDecayLVec.at(ig));
      }
    }    
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

      //dR betwn tau-b
      for(int i=0; i<genbLVec.size(); i++){
	float dR_taub = genTauLVec.DeltaR(genbLVec[i]);
	myBaseHistgram.htaub->Fill(dR_taub, Lumiscale);      
	myBaseHistgram.htaub_dRpt->Fill(genTauLVec.Pt(), dR_taub, Lumiscale);      
      }

      // Do the matching
      int tauJetIdx = -1; // Will store the index of the jet matched to the tau
      const double DeltaR = 0.4; //deltaR value to throw taujet overlapped with bjet
      const float deltaRMax = deltaRmax(genTauLVec.Pt()); // Increase deltaRMax at low pt to maintain high-enought matching efficiency
      if( !utils::findTauMatchedJet(tauJetIdx,genTauLVec,jetsLVec,deltaRMax) ) continue;
      
      //dR betwn taujet-b     

      vector<float>dR_taub(genbLVec.size(),0);

      for(int i=0; i<genbLVec.size(); i++){
        float dR_taujetb = jetsLVec.at(tauJetIdx).DeltaR(genbLVec[i]);
        myBaseHistgram.htaujetb->Fill(dR_taujetb, Lumiscale);
	dR_taub[i] = genTauLVec.DeltaR(genbLVec[i]);
      }
      vector<float>::iterator dR_taubmin = std::min_element(dR_taub.begin(),dR_taub.end());
      double taujetPt = jetsLVec.at(tauJetIdx).Pt();
      if(dR_taub.size()){
	if(taujetPt>=30 && taujetPt<50  ) FillDoubleFakeRate(myBaseHistgram.htaub1, *dR_taubmin, Lumiscale);
	if(taujetPt>=50 && taujetPt<100 ) FillDoubleFakeRate(myBaseHistgram.htaub2, *dR_taubmin, Lumiscale);
	if(taujetPt>=100 && taujetPt<200) FillDoubleFakeRate(myBaseHistgram.htaub3, *dR_taubmin, Lumiscale);
	if(taujetPt>=200 && taujetPt<400) FillDoubleFakeRate(myBaseHistgram.htaub4, *dR_taubmin, Lumiscale);
	if(taujetPt>=400                ) FillDoubleFakeRate(myBaseHistgram.htaub5, *dR_taubmin, Lumiscale);
      }
      // if(utils::findBMatchedTauJet(tauJetIdx, genbLVec,jetsLVec, DeltaR))continue;//exclude taujet-B overlapping
      myBaseHistgram.htauBjetEta_den->Fill(jetsLVec.at(tauJetIdx).Eta(), Lumiscale);
      myBaseHistgram.htauBjetPt_den->Fill(jetsLVec.at(tauJetIdx).Pt(), Lumiscale);
      
      // select the had-tau jet
      cleanJetVec.push_back(jetsLVec.at(tauJetIdx));
      cleanJetBtag.push_back(recoJetsBtag_0.at(tauJetIdx));
      
      int cnt1CSVS = AnaFunctions::countCSVS(cleanJetVec, cleanJetBtag, AnaConsts::cutCSVS, AnaConsts::bTagArr);
      //if( (AnaConsts::low_nJetsSelBtagged == -1 || cnt1CSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cnt1CSVS < AnaConsts::high_nJetsSelBtagged ) ){
      if(cnt1CSVS>=1){
	myBaseHistgram.htauBjetEta_num->Fill(jetsLVec.at(tauJetIdx).Eta(), Lumiscale);
	myBaseHistgram.htauBjetPt_num->Fill(jetsLVec.at(tauJetIdx).Pt(), Lumiscale);
      }
  }//finish gen hadtau loop
  
  }//finish event loop
  std::cout<<"final event: "<<k<<std::endl;
  myBaseHistgram.hEff_Pt = (TH1D*)myBaseHistgram.htauBjetPt_num->Clone("Ratio");
  myBaseHistgram.hEff_Pt->GetYaxis()->SetRangeUser(0. , 0.7);
  myBaseHistgram.hEff_Pt->Divide(myBaseHistgram.htauBjetPt_den);
  myBaseHistgram.hEff_Eta = (TH1D*)myBaseHistgram.htauBjetEta_num->Clone("Ratio");
  myBaseHistgram.hEff_Eta->Divide(myBaseHistgram.htauBjetEta_den);
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
