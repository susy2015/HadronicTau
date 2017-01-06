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
#include "SusyAnaTools/Tools/BTagCorrector.h"
#include "HadTauLL.h"
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

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
    {
      std::cerr <<"Please give 5 arguments "<<"SubsampleName"<<" Input Template" <<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./HadTauMC TTbarInc 1000 0 1" << std::endl;
      return -1;
    }
  const char *subsamplename = argv[1];
  const char *Maxevent = argv[2];
  const  char *Stratfile = argv[3];
  const  char *Filerun = argv[4];
  const  int startfile = std::atoi(Stratfile);
  const int filerun = std::atoi(Filerun);
  TChain *fChain = 0;
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(subsamplename, startfile);
  const string condorSpec = argc==6 ? argv[5]: "";

  if(!FillChain(fChain, subsamplename, condorSpec, startfile, filerun))
       {
	 std::cerr << "Cannot get the tree " << std::endl;
       }
  const int maxevent = std::atoi(Maxevent);  
  TString sampleString(subsamplename);

  //Searchbin                          
  SearchBins SB("SB_v1_2017");
  
  //Use BaselineVessel class for baseline variables and selections
  std::string spec = "MCExp";
  AnaFunctions::prepareForNtupleReader();
  NTupleReader *tr =0;
  tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames);
  ExpBaseline = new BaselineVessel(*tr, spec);
  ExpBaseline->SetupTopTagger(true,"TopTagger.cfg");
  tr->registerFunction((*ExpBaseline));
 
  //BTag SF
  BTagCorrector btagcorr;
  tr->registerFunction(btagcorr);
  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr->getNEntries();
  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<Lumiscale<<std::endl; 
  cout<<"maxevent: "<<maxevent<<endl;
  // Loop over the events (tree entries)
  int k = 0;
  while(tr->getNextEvent()){
    k++;
    std::cout<<"Evt: "<<k<<std::endl;
    if(maxevent>=0 && tr->getEvtNum() > maxevent ) break;
    // Add print out of the progress of looping
    if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;

    const vector<TLorentzVector> &jetsLVec = tr->getVec<TLorentzVector>("jetsLVec");
    const std::vector<std::string> & TriggerNames = tr->getVec<std::string>("TriggerNames");
    const std::vector<int> & PassTrigger = tr->getVec<int>("PassTrigger");

    int run = tr->getVar<int>("run");
    int lumi = tr->getVar<int>("lumi");
    int event = tr->getVar<int>("event");
    bool passNoiseEventFilter = tr->getVar<bool>("passNoiseEventFilter"+spec);
    const double EvtWt = tr->getVar<double>("evtWeight");
    //change event weight for MC sample
    EventWeight = EvtWt;
    Lumiscale = Lumiscale * EventWeight;
    //Event Filter                                                                                                                           
    //if(!passNoiseEventFilter) continue;
    
    const vector<TLorentzVector> &genDecayLVec = tr->getVec<TLorentzVector>("genDecayLVec");
    const vector<int> &genDecayIdxVec =  tr->getVec<int>("genDecayIdxVec");
    const vector<int> &genDecayPdgIdVec = tr->getVec<int>("genDecayPdgIdVec");
    const vector<int> &W_emuVec =  tr->getVec<int>("W_emuVec");
    const vector<int> &W_tau_emuVec =  tr->getVec<int>("W_tau_emuVec");
    const vector<int> &W_tau_prongsVec =  tr->getVec<int>("W_tau_prongsVec");
    
    bool passBaseline = tr->getVar<bool>("passBaseline"+spec);
    bool passMuonVeto = tr->getVar<bool>("passMuonVeto"+spec);
    bool passEleVeto = tr->getVar<bool>("passEleVeto"+spec);
    bool passIsoTrkVeto = tr->getVar<bool>("passIsoTrkVeto"+spec);
    bool passLeptVeto = tr->getVar<bool>("passLeptVeto"+spec);
    bool passnJets = tr->getVar<bool>("passnJets"+spec);
    bool passdPhis = tr->getVar<bool>("passdPhis"+spec);
    bool passMET =  tr->getVar<bool>("passMET"+spec);
    bool passBJets = tr->getVar<bool>("passBJets"+spec);
    bool passTagger = tr->getVar<bool>("passTagger"+spec);
    bool passHT =  tr->getVar<bool>("passHT"+spec);
    bool passMT2 = tr->getVar<bool>("passMT2" + spec);
    const int nElectrons = tr->getVar<int>("nElectrons_CUT"+spec);
    const int nMuons = tr->getVar<int>("nMuons_CUT"+spec);
    const int nJets = tr->getVar<int>("cntNJetsPt30Eta24"+spec);
    const int nbJets = tr->getVar<int>("cntCSVS"+spec);
    const int nTops = tr->getVar<int>("nTopCandSortedCnt"+spec);
    const double MT2 = tr->getVar<double>("best_had_brJet_MT2"+spec);
    const vector<double> dPhiVec = tr->getVec<double>("dPhiVec"+spec);
    const double ht = tr->getVar<double>("HT"+spec);
    double met=tr->getVar<double>("met");
    double metphi=tr->getVar<double>("metphi");
    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);
    const double bSF = tr->getVar<double>("bTagSF_EventWeightSimple_Central");
    //mht calculation
    TLorentzVector Mht_LVec;	      
    for(unsigned int ij=0; ij<jetsLVec.size(); ij++){
      if( !AnaFunctions::jetPassCuts(jetsLVec[ij], AnaConsts::pt30Arr) ) continue;
      Mht_LVec -= jetsLVec[ij];
    }
    const double Mht = Mht_LVec.Pt();
    
    bool passBaselineFull = passMuonVeto && passEleVeto && passIsoTrkVeto && passnJets && passdPhis && passMET && passBJets && passTagger && passHT && passMT2;
    
    if(!(passBaselineFull && passNoiseEventFilter)) continue;
    //Exp LostLepton Dist.
    if(W_emuVec.size() !=0 || W_tau_emuVec.size() !=0){
      int kSR = SB.find_Binning_Index(nbJets, nTops, MT2, met, ht);
      if( kSR!= -1 ) {
	myBaseHistgram.hYields_LL->Fill(kSR, bSF*Lumiscale);
      }
      
      FillDouble(myBaseHistgram.hMET_LL, met, Lumiscale);
      FillDouble(myBaseHistgram.hMT2_LL, MT2, Lumiscale);
      FillInt(myBaseHistgram.hNbJets_LL, nbJets, bSF*Lumiscale);
      FillInt(myBaseHistgram.hNTops_LL, nTops, Lumiscale);	
      FillInt(myBaseHistgram.hNJets_LL, nJets, Lumiscale);
      FillDouble(myBaseHistgram.hHT_LL, ht, Lumiscale);
      FillDouble(myBaseHistgram.hdPhi0_LL, dPhiVec[0], Lumiscale);
      FillDouble(myBaseHistgram.hdPhi1_LL, dPhiVec[1], Lumiscale);
      FillDouble(myBaseHistgram.hdPhi2_LL, dPhiVec[2], Lumiscale);
    }        
    //Exp Hadtau Dist.
    else if(W_tau_prongsVec.size() !=0){
      int jSR = SB.find_Binning_Index(nbJets, nTops, MT2, met, ht);
      if( jSR!= -1 ) {
	myBaseHistgram.hYields_tau->Fill(jSR, bSF*Lumiscale);
      }
      
      FillDouble(myBaseHistgram.hMET_tau, met, Lumiscale);
      FillDouble(myBaseHistgram.hMT2_tau, MT2, Lumiscale);
      FillInt(myBaseHistgram.hNbJets_tau, nbJets, bSF*Lumiscale);
      FillInt(myBaseHistgram.hNTops_tau, nTops, Lumiscale);	
      FillInt(myBaseHistgram.hNJets_tau, nJets, Lumiscale);
      FillDouble(myBaseHistgram.hHT_tau, ht, Lumiscale);
      FillDouble(myBaseHistgram.hdPhi0_tau, dPhiVec[0], Lumiscale);
      FillDouble(myBaseHistgram.hdPhi1_tau, dPhiVec[1], Lumiscale);
      FillDouble(myBaseHistgram.hdPhi2_tau, dPhiVec[2], Lumiscale);
    }

  }	//event loop
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();
  
  return 0;
  
}



