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
  if (argc < 1)
    {
      std::cerr <<"Please give 1 arguments " << "inputList " << std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Acc List1_ttbar.txt " << std::endl;
      return -1;
    }
  const char *inputFileList = argv[1];
  TChain *fChain = new TChain("stopTreeMaker/AUX");
  if(!FillChain(fChain, inputFileList))
    {
      std::cerr << "Cannot get the tree " << std::endl;
    }
  NTupleReader tr(fChain);
  tr.registerFunction(&passBaselineFunc);

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
    double metphi=tr.getVar<double>("metphi");

    // bool passBaseline_nolepveto = tr.getVar<bool>("passBaseline_nolepveto");
    //bool passnJets = tr.getVar<bool>("passnJets");
    bool passMET = tr.getVar<bool>("passMET");
    //bool passdPhis = tr.getVar<bool>("passdPhis");
    //bool passBJets = tr.getVar<bool>("passBJets");



    int nmu=0, nele=0;
    vector<double> genmueta, genmuphi;
    genmueta.clear();
    genmuphi.clear();

    for(unsigned ig=0; ig<genDecayLVec.size(); ig++){
      int pdgId = genDecayPdgIdVec.at(ig);

      if(fabs(pdgId)==13){
	nmu++;
	genmueta.push_back(genDecayLVec.at(ig).Eta());
	genmuphi.push_back(genDecayLVec.at(ig).Phi());
      }
      if(fabs(pdgId)==11) nele++;
    }

    if(nmu ==1 && nele==0){//control sample
      
      const double muEta = genmueta.at(0);
      const double muPhi = genmuphi.at(0);
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

      //recompute jet
      int nJet = AnaFunctions::countJets(cleanJetVec, AnaConsts::pt30Eta24Arr);
      bool passnjets = true;
      if(nJet<AnaConsts::nJetsSelPt30Eta24)passnjets = false;
      //      if(nJet<3)passnjets = false;

      //recompute deltaphi
      std::vector<double> * deltaPhiVec = new std::vector<double>();
      (*deltaPhiVec) = AnaFunctions::calcDPhi(cleanJetVec, metphi, 3, AnaConsts::dphiArr);
      bool passdeltaPhi = true;
      if( deltaPhiVec->at(0) < AnaConsts::dPhi0_CUT || deltaPhiVec->at(1) < AnaConsts::dPhi1_CUT || deltaPhiVec->at(2) < AnaConsts::dPhi2_CUT){
	passdeltaPhi = false;
      }
      //recompute bjet
      int cnt1CSVS = AnaFunctions::countCSVS(cleanJetVec, tr.getVec<double>("recoJetsBtag_0"), AnaConsts::cutCSVS, AnaConsts::bTagArr);
      bool passbJets = true;
      if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cnt1CSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cnt1CSVS < AnaConsts::high_nJetsSelBtagged ) ) ){
	passbJets = false;
      }

      if(!passnjets) continue;
      if(!passMET) continue;
      if(!passdeltaPhi) continue;
      if(!passbJets) continue;

      // if(!passBaseline_nolepveto) continue;

      for(unsigned ig=0; ig<genDecayLVec.size(); ig++){
	int pdgId = genDecayPdgIdVec.at(ig);

	if(fabs(pdgId)==13){
        genmu++;
        double pt1= genDecayLVec.at(ig).Pt();
        double eta1 = genDecayLVec.at(ig).Eta();
        double phi1 = genDecayLVec.at(ig).Phi();
        if(pt1>20 && fabs(eta1)<2.4)accmu++;
	}
      }
    }

  }
  cout<<"ToTal Event: "<<k<<endl;
  cout<<"Gen Mu: "<<genmu<<" "<<"Acc Mu: "<<accmu<<endl;
  cout<<"Acceptance: "<<(float)accmu/genmu<<endl;
  return 0;
}
