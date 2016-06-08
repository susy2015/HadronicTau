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
#include "TStopwatch.h"
#include "TString.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Template.h"

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




// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
    {
      std::cerr <<"Please give 4 arguments "<<"SubsampleName"<<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Template TTbarInc 1000 0 1" << std::endl;
      return -1;
    }
  const char *subsamplename = argv[1]; 
  const char *Maxevent = argv[2];
  const  char *Stratfile = argv[3];
  const  char *Filerun = argv[4];
  const  int startfile = std::atoi(Stratfile);
  const int filerun = std::atoi(Filerun);  
  const int maxevent = std::atoi(Maxevent);
  TChain *fChain = 0;
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(subsamplename, startfile);
  const string condorSpec = argc==6 ? argv[5]: "";  
  if(!FillChain(fChain, subsamplename, condorSpec, startfile, filerun))
    {
      std::cerr << "Cannot get the tree " << std::endl;
    }

  AnaFunctions::prepareForNtupleReader();
  NTupleReader *tr =0;

  tr = new NTupleReader(fChain);



  
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

    const vector<TLorentzVector> &genDecayLVec = tr->getVec<TLorentzVector>("genDecayLVec");
    const vector<int> &genDecayIdxVec = tr->getVec<int>("genDecayIdxVec");
    const vector<int> &genDecayPdgIdVec = tr->getVec<int>("genDecayPdgIdVec");
    const vector<int> &genDecayMomIdxVec = tr->getVec<int>("genDecayMomIdxVec");
    const vector<int> &W_emuVec = tr->getVec<int>("W_emuVec");
    const vector<int> &W_tau_emuVec = tr->getVec<int>("W_tau_emuVec");
    const vector<int> &W_tau_prongsVec = tr->getVec<int>("W_tau_prongsVec");
    const vector<int> &W_tau_nuVec = tr->getVec<int>("W_tau_nuVec");
    const vector<TLorentzVector> &jetsLVec = tr->getVec<TLorentzVector>("jetsLVec");
    const vector<double> &recoJetsJecUnc = tr->getVec<double>("recoJetsJecUnc");

    double HT = AnaFunctions::calcHT(jetsLVec, AnaConsts::pt30Eta24Arr);
    // Select only events where the W decayed into a hadronically decaying tau
    if(W_tau_prongsVec.size()==0) continue;
    if(HT<100) continue;

    vector<TLorentzVector> genhadtauLVec;
    vector<TLorentzVector> genvisiblehadtauLVec;
    vector<TLorentzVector> bLVec;
    
    //gen Tau information
    
    for(unsigned ig=0; ig<genDecayLVec.size(); ig++){
      int pdgId = genDecayPdgIdVec.at(ig);
      if(abs(pdgId)==15){
	int flag1=0;
	if(W_tau_emuVec.size()!=0){
	  for(int k=0; k<W_tau_emuVec.size();k++){
	    int lIdx = W_tau_emuVec.at(k);
	    if( genDecayMomIdxVec.at(lIdx) == genDecayIdxVec.at(ig) )flag1++;
	  }
	}
	if(!flag1){
	  TLorentzVector gennuVec;
	  TLorentzVector genhadtauVec = genDecayLVec.at(ig);
	  genhadtauLVec.push_back(genhadtauVec);
	  for(int n=0; n<W_tau_nuVec.size();n++){
	    int nIdx = W_tau_nuVec.at(n);
	    if( genDecayMomIdxVec.at(nIdx) == genDecayIdxVec.at(ig)) gennuVec = genDecayLVec.at(nIdx);
	  }
	  TLorentzVector genvisiblehadtauVec = genhadtauVec-gennuVec;
	  genvisiblehadtauLVec.push_back(genvisiblehadtauVec);
	}
      }
    }
    
    for(unsigned ib=0; ib<genDecayLVec.size(); ib++){
      int pdgId = genDecayPdgIdVec.at(ib);
      if(abs(pdgId)==5) bLVec.push_back(genDecayLVec.at(ib));
    }
    if(genhadtauLVec.size()==0)continue;//To skip rare events where same tau decays both into hadron and lepton!
    TLorentzVector genTauLVec = genhadtauLVec[0];
    TLorentzVector genvisibleLVec = genvisiblehadtauLVec[0];
    // Use only events where the tau is inside the muon acceptance (pt>20, eta>2.1)because lateron we will apply the response to muon+jet events
    if( genTauLVec.Pt() < TauResponse::ptMin() ) continue;
    if( fabs(genTauLVec.Eta()) > TauResponse::etaMax() ) continue;
    
    // Do the matching
    int tauJetIdx = -1;		// Will store the index of the jet matched to the tau
    const double DeltaR = 0.4; //deltaR value to throw taujet overlapped with bjet
    const float deltaRMax = deltaRmax(genTauLVec.Pt()); // Increase deltaRMax at low pt to maintain high-enought matching efficiency
    if( !utils::findTauMatchedJet(tauJetIdx, genvisibleLVec, jetsLVec, deltaRMax) ) continue;
    if(utils::findBMatchedTauJet(tauJetIdx, bLVec, jetsLVec, DeltaR))continue;
    // Fill histogram with relative visible energy of the tau
    // ("tau response template") for hadronically decaying taus
    for(unsigned jetIdx = 0; jetIdx < jetsLVec.size(); ++jetIdx) {	// Loop over reco jets
      // Select tau jet
      if( jetIdx == tauJetIdx ) {
	// Get the response pt bin for the tau
	const unsigned int ptBin = TauResponse::ptBin(genTauLVec.Pt());
	// Fill the corresponding response template
	const double tauJetPt = jetsLVec[jetIdx].Pt();
	myBaseHistgram.hTauResp.at(ptBin)->Fill( tauJetPt / genTauLVec.Pt(), Lumiscale );
	myBaseHistgram.hTauvisible.at(ptBin)->Fill(tauJetPt, Lumiscale);
	myBaseHistgram.hTaugenerated.at(ptBin)->Fill(genTauLVec.Pt(), Lumiscale);
	// tau energy scale corr for systematics.
	const double tauJetPtUp = jetsLVec[jetIdx].Pt() + (jetsLVec[jetIdx].Pt()*recoJetsJecUnc[jetIdx]);
	const double tauJetPtDown = jetsLVec[jetIdx].Pt() - (jetsLVec[jetIdx].Pt()*recoJetsJecUnc[jetIdx]);
	myBaseHistgram.hTauRespUp.at(ptBin)->Fill( tauJetPtUp / genTauLVec.Pt(), Lumiscale );
	myBaseHistgram.hTauRespDown.at(ptBin)->Fill( tauJetPtDown / genTauLVec.Pt(), Lumiscale );
	break;		// End the jet loop once the tau jet has been found
      }
    }	// End of loop over reco jets
  } // End of loop over tree entries
  
  /*  // Normalize the response distributions to get the probability density
      for(unsigned int i = 0; i <myBaseHistgram. hTauResp.size(); ++i) {
      if(myBaseHistgram. hTauResp.at(i)->Integral("width") > 0. ) {
      myBaseHistgram.hTauResp.at(i)->Scale(1./myBaseHistgram.hTauResp.at(i)->Integral("width"));
      }
      }*/
  
  
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
