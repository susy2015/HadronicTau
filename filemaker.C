#include <iostream>
#include<vector>
#include <cmath>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1F.h"
#include "TStyle.h"

void filemaker()
{
  TChain *t = new TChain("stopTreeMaker/AUX");
 for(int i=1; i<72; i++){
   t->Add(Form("/eos/uscms/store/user/lpcsusyhad/PHYS14_720_Mar14_2014_v2/benwu/PU20bx25_WJetsToLNu_HT-600toInf_madgraph-tauola/WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/150430_142404/0000/stopFlatNtuples_%d.root",i));
 }
 cout<<"merging done.."<<endl;
  
 
  //    TFile f("stopFlatNtuples_1.root");
  //TTree *t = (TTree*)f.Get("stopTreeMaker/AUX");

  std::vector<TLorentzVector> *genDecayLVec =0, *jetsLVec=0;
  std::vector<int> *genDecayIdxVec =0, *genDecayPdgIdVec =0, *genDecayMomIdxVec =0;
  std::vector<int> *W_emuVec=0, *W_tau_emuVec=0, *W_tau_prongsVec=0, *W_tau_nuVec=0;
  std::vector<string> *genDecayStrVec =0;

  t->SetBranchAddress("genDecayIdxVec", &genDecayIdxVec);
  t->SetBranchAddress("genDecayPdgIdVec", &genDecayPdgIdVec);
  t->SetBranchAddress("genDecayMomIdxVec", &genDecayMomIdxVec);
  t->SetBranchAddress("genDecayLVec", &genDecayLVec);
  t->SetBranchAddress("genDecayStrVec", &genDecayStrVec);
  t->SetBranchAddress("W_emuVec", &W_emuVec);
  t->SetBranchAddress("W_tau_emuVec", &W_tau_emuVec);
  t->SetBranchAddress("W_tau_prongsVec", &W_tau_prongsVec);
  t->SetBranchAddress("W_tau_nuVec", &W_tau_nuVec);
  t->SetBranchAddress("jetsLVec", &jetsLVec);

  std::vector<double> *jetspt, *jetseta, *jetsphi, *bjetpt, *bjeteta, *bjetphi;
  std::vector<double> *genhadtaupt, *genhadtaueta, *genhadtauphi, *genvisiblehadtaupt, *genvisiblehadtaueta, *genvisiblehadtauphi;
  int hadtauflag;
  TFile *f1 = new TFile("hadtau.root", "RECREATE");
  TTree *nt = new TTree("Hadtau","Hadtau info");
  TH1F *h1 = new TH1F("h1", "muon pT", 100, 0, 1000);
  TH1F *h2 = new TH1F("h2", "hadtau pT", 100, 0, 1000);
  nt->Branch("bjetpt", "std::vector<double>", &bjetpt, 3200, 0);
  nt->Branch("bjeteta", "std::vector<double>", &bjeteta, 3200, 0);
  nt->Branch("bjetphi", "std::vector<double>", &bjetphi, 3200, 0);
  nt->Branch("jetspt", "std::vector<double>", &jetspt, 3200, 0);
  nt->Branch("jetseta", "std::vector<double>", &jetseta, 3200, 0);
  nt->Branch("jetsphi", "std::vector<double>", &jetsphi, 3200, 0);
  nt->Branch("genhadtaupt", "std::vector<double>", &genhadtaupt, 3200, 0);
  nt->Branch("genhadtaueta", "std::vector<double>", &genhadtaueta, 3200, 0);
  nt->Branch("genhadtauphi", "std::vector<double>", &genhadtauphi, 3200, 0);
  nt->Branch("genvisiblehadtaupt", "std::vector<double>", &genvisiblehadtaupt, 3200, 0);
  nt->Branch("genvisiblehadtaueta", "std::vector<double>", &genvisiblehadtaueta, 3200, 0);
  nt->Branch("genvisiblehadtauphi", "std::vector<double>", &genvisiblehadtauphi, 3200, 0);
  nt->Branch("hadtauflag", &hadtauflag, "hadtauflag/I");


  Int_t nentries = (Int_t)t->GetEntries();
  std::cout<<nentries<<std::endl;
  for(Int_t i=0; i<nentries; i++) {
    t->GetEntry(i);

    bjetpt->clear();
    bjeteta->clear();
    bjetphi->clear();
    jetspt->clear();
    jetseta->clear();
    jetsphi->clear();
    genhadtaupt->clear();
    genhadtaueta->clear();
    genhadtauphi->clear();
    genvisiblehadtaupt->clear();
    genvisiblehadtaueta->clear();
    genvisiblehadtauphi->clear();

    for(unsigned ij=0; ij<jetsLVec->size(); ij++){
      double pt1 = jetsLVec->at(ij).Pt();
      double eta1 = jetsLVec->at(ij).Eta();
      double phi1 = jetsLVec->at(ij).Phi();
      jetspt->push_back(pt1);
      jetseta->push_back(eta1);
      jetsphi->push_back(phi1);
     }  

    int flag = 0;
    if(W_tau_prongsVec->size()!=0){
      flag = 1;
     }
    hadtauflag =flag;

    for(unsigned ig=0; ig<genDecayLVec->size(); ig++){
      int pdgId = genDecayPdgIdVec->at(ig);
      if(abs(pdgId)==5){
	TLorentzVector genbjetLVec = genDecayLVec->at(ig);
	bjetpt->push_back(genbjetLVec.Pt());
	bjeteta->push_back(genbjetLVec.Eta());
	bjetphi->push_back(genbjetLVec.Phi());
      }
    }
    
    if(W_tau_prongsVec->size()!=0){
	for(unsigned ig=0; ig<genDecayLVec->size(); ig++){
	  int pdgId = genDecayPdgIdVec->at(ig);
	  if(abs(pdgId)==15){
	    int flag1=0;
	    if(W_tau_emuVec->size()!=0){
	      for(int k=0; k<W_tau_emuVec->size();k++){
		int lIdx = W_tau_emuVec->at(k);
		if( genDecayMomIdxVec->at(lIdx) == genDecayIdxVec->at(ig) )flag1++;
	      }
	    }
	    if(!flag1){
	      TLorentzVector gennuLVec;
	      TLorentzVector genhadtauLVec = genDecayLVec->at(ig);
	      double pt2 = genhadtauLVec.Pt();
	      double eta2 = genhadtauLVec.Eta();
	      double phi2 = genhadtauLVec.Phi();
	      genhadtaupt->push_back(pt2);
	      genhadtaueta->push_back(eta2);
	      genhadtauphi->push_back(phi2);	
	      for(int n=0; n<W_tau_nuVec->size();n++){
		int nIdx = W_tau_nuVec->at(n);
		if( genDecayMomIdxVec->at(nIdx) == genDecayIdxVec->at(ig)) gennuLVec = genDecayLVec->at(nIdx);
	      }
	      TLorentzVector genvisiblehadtauLVec = genhadtauLVec-gennuLVec;    
	      double pt3 = genvisiblehadtauLVec.Pt();
	      double eta3 = genvisiblehadtauLVec.Eta();
	      double phi3 = genvisiblehadtauLVec.Phi();
	      genvisiblehadtaupt->push_back(pt3);
	      genvisiblehadtaueta->push_back(eta3);
	      genvisiblehadtauphi->push_back(phi3);
	      h2->Fill(pt2);
	    }
	  }
	}
    }
    if(W_emuVec->size()!=0 && W_tau_emuVec->size()==0){
      for(unsigned im=0; im<genDecayLVec->size(); im++){
	int pdgId = genDecayPdgIdVec->at(im);
	if(abs(pdgId)==13){
	  h1->Fill(genDecayLVec->at(im).Pt());
	}
      }
    }
    nt->Fill();

  }
  nt->Write();
  h1->Write();
  h2->Write();
  cout<< "total: "<<nentries<<endl;
}
