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
#include "Closure.h"
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
#include "Efficiency.h"

using namespace std;


void passBaselineFunc(NTupleReader &tr)
{
  bool passBaseline = true;
  //Form TLorentzVector of MET
  TLorentzVector metLVec; metLVec.SetPtEtaPhiM(tr.getVar<double>("met"), 0, tr.getVar<double>("metphi"), 0);
  //Calculate number of leptons
  int nMuons = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsRelIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsArr);
  int nElectrons = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesRelIso"), tr.getVec<double>("elesMtw"), AnaConsts::elesArr);
  //  int nIsoTrks = AnaFunctions::countIsoTrks(tr.getVec<TLorentzVector>("loose_isoTrksLVec"), tr.getVec<double>("loose_isoTrks_iso"), tr.getVec<double>("loose_isoTrks_mtw"), AnaConsts::isoTrksArr);

  //Calculate number of jets and b-tagged jets
  int cntCSVS = AnaFunctions::countCSVS(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), AnaConsts::cutCSVS, AnaConsts::bTagArr);
  //    int cntNJetsPt50Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt50Eta24Arr);
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
  //if( cntNJetsPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 ){ passBaseline = false; passnJets = false;}
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
  int bestTopJetIdx = -1;
  bool remainPassCSVS = false;
  int pickedRemainingCombfatJetIdx = -1;
  double bestTopJetMass = -1;
  int nTopCandSortedCnt = 0;
  double MT2 = -1;
  double mTcomb = -1;
  if( passBaseline && cntNJetsPt30 >= AnaConsts::nJetsSel )
    {
      type3Ptr->processEvent((*jetsLVec_forTagger), (*recoJetsBtag_forTagger), metLVec);
      bestTopJetIdx = type3Ptr->bestTopJetIdx;
      remainPassCSVS = type3Ptr->remainPassCSVS;
      pickedRemainingCombfatJetIdx = type3Ptr->pickedRemainingCombfatJetIdx;
      if( bestTopJetIdx != -1 ) bestTopJetMass = type3Ptr->bestTopJetLVec.M();
      nTopCandSortedCnt = type3Ptr->nTopCandSortedCnt;
      MT2 = type3Ptr->MT2;
      mTcomb = type3Ptr->mTbJet + 0.5*type3Ptr->mTbestTopJet;
    }

  //Pass top tagger requirement?
  bool passTagger = true;
  //bestTopJetIdx != -1 means at least 1 top candidate!
  if( bestTopJetIdx == -1 ){ passBaseline = false; passTagger = false; }
  if( ! remainPassCSVS ){ passBaseline = false; passTagger = false; }
  if( pickedRemainingCombfatJetIdx == -1 && jetsLVec_forTagger->size()>=6 ){ passBaseline = false; passTagger = false; }
  if( ! (bestTopJetMass > AnaConsts::lowTopCut_ && bestTopJetMass < AnaConsts::highTopCut_ ) ){ passBaseline = false; passTagger = false; }


  //register new var

  tr.registerDerivedVar("nMuons_CUT2", nMuons);
  tr.registerDerivedVar("nElectrons_CUT2", nElectrons);
  tr.registerDerivedVar("cntNJetsPt30Eta24", cntNJetsPt30Eta24);

  tr.registerDerivedVar("passBaseline", passBaseline);
  tr.registerDerivedVar("passLeptVeto", passLeptVeto);
  tr.registerDerivedVar("passMET", passMET);
  tr.registerDerivedVar("passnJets", passnJets);
  tr.registerDerivedVar("passdPhis", passdPhis);
  tr.registerDerivedVar("passBJets", passBJets);
  tr.registerDerivedVar("passTagger", passTagger);


}

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 3)
    {
      std::cerr <<"Please give 3 arguments " << "inputList " << " " <<" "<<"input template"<<" "<< "outputFileName" << std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Closure List1_ttbar.txt HadTau_TauResponseTemplates.root HadTau_Closure.root" << std::endl;
      return -1;
    }
  const char *inputFileList = argv[1];
  const char *respTempl = argv[2];  
  const char *outFileName = argv[3];
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

  TauResponse tauResp(respTempl);


  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  // Loop over the events (tree entries)
  int k = 0;
  while(tr.getNextEvent()){
    k++;
    vector<TLorentzVector> muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");
    vector<double> muonsRelIso = tr.getVec<double>("muonsRelIso");
    vector<TLorentzVector> jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
    int nElectrons = tr.getVar<int>("nElectrons_CUT2");
    int nMuons = tr.getVar<int>("nMuons_CUT2");
    double met=tr.getVar<double>("met");
    double metphi=tr.getVar<double>("metphi");
    double ht=tr.getVar<double>("ht");
    vector<int> W_tau_prongsVec = tr.getVec<int>("W_tau_prongsVec");


    //Prediction part
    //Control sample

      // The kinematic properties of the well-reconstructed, isolated muon                                                                    
      vector<TLorentzVector>isomuonsLVec;
      isomuonsLVec.clear();      

      for(unsigned im=0; im<muonsLVec.size(); im++){                                                                                   
	if(muonsRelIso.at(im)<0.2){
	  const double muPt1 = muonsLVec.at(im).Pt();
	  const double muEta1 = muonsLVec.at(im).Eta();
	  const double muPhi1 = muonsLVec.at(im).Phi();
	  const double muM1 = muonsLVec.at(im).M();
	  TLorentzVector isomu; isomu.SetPtEtaPhiM(muPt1, muEta1, muPhi1, muM1);
	  isomuonsLVec.push_back(isomu);
	}                                                                                                                                  
      }

    if( nMuons == 1 && nElectrons == 0 ) {

	  const double muPt = isomuonsLVec.at(0).Pt();
	  const double muEta = isomuonsLVec.at(0).Eta();
	  const double muPhi = isomuonsLVec.at(0).Phi();
	  const double muM = isomuonsLVec.at(0).M();
	  // Use only events where the muon is inside acceptance                                                                             
	  if( muPt < TauResponse::ptMin() ) continue;
	  if( fabs(muEta) > TauResponse::etaMax() ) continue;
      // "Cross cleaning": find the jet that corresponds to                                                                                  
      // the muon. Associate the jet that is closestin eta-phi space to the lepton                                                           

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
      // Calculate RA2 selection-variables from "cleaned" jets                                                                                
      vector<TLorentzVector> selNJetVec;
      selNJetVec.clear();
      double selHt   = 0.;
      double selMhtX = 0.;
      double selMhtY = 0.;
      double selmetX =0;
      double selmetY =0;

      for(int jetIdx = 0; jetIdx < jetsLVec.size(); ++jetIdx) {  // Loop over reco jets                                                        
        // Skip this jet if it is the muon                                                                                                    
        if( jetIdx == muJetIdx ) continue;

        // Calculate NJet                                                                                                             

	double Pt = jetsLVec.at(jetIdx).Pt();
	double Eta = jetsLVec.at(jetIdx).Eta();
	double Phi = jetsLVec.at(jetIdx).Phi();
	double M = jetsLVec.at(jetIdx).M();
	TLorentzVector selLVec; selLVec.SetPtEtaPhiM(Pt, Eta, Phi, M);
	selNJetVec.push_back(selLVec);

         // Calculate MHT and HT                                                                                                              
	if( jetsLVec.at(jetIdx).Pt() > htJetPtMin() && fabs(jetsLVec.at(jetIdx).Eta()) < htJetEtaMax()) {
          selHt += jetsLVec.at(jetIdx).Pt();
	}

        if( jetsLVec.at(jetIdx).Pt() > mhtJetPtMin() && fabs(jetsLVec.at(jetIdx).Eta()) < mhtJetEtaMax()) {
          selMhtX -= jetsLVec.at(jetIdx).Pt()*cos(jetsLVec.at(jetIdx).Phi());
          selMhtY -= jetsLVec.at(jetIdx).Pt()*sin(jetsLVec.at(jetIdx).Phi());
	}

      }

      //calculate met
      selmetX = met*cos(metphi)+ muPt*cos(muPhi);
      selmetY = met*sin(metphi)+ muPt*sin(muPhi);

      int selNJet = AnaFunctions::countJets(selNJetVec, AnaConsts::pt30Eta24Arr);
       //      if( selNJet < 3 ) continue;
       if( selNJetVec.size() < 3 ) continue;
      // Get random number from tau-response template                                                                                         
      // The template is chosen according to the muon pt                                                                                      
      const double scale = tauResp.getRandom(muPt);
      // Scale muon pt with tau response --> simulate tau jet pt                                                                              
      const double simTauJetPt = scale * muPt;
      const double simTauJetEta = muEta;
      const double simTauJetPhi = muPhi;

      // Taking into account the simulated tau jet, recompute                                                                                 
      // HT, MHT, and N(jets)                                                                                                                 
      vector<TLorentzVector> simNJetVec;
      simNJetVec.clear();
      double simHt = selHt;
      double simMhtX = selMhtX;
      double simMhtY = selMhtY;
      double simmetX = selmetX;
      double simmetY = selmetY;
      // recompute NJets                                                                                                               
      TLorentzVector simLVec; simLVec.SetPtEtaPhiM(simTauJetPt, simTauJetEta, simTauJetPhi, muM);
      simNJetVec.push_back(simLVec);

      // If simulted tau-jet meets same criteria as as HT jets,                                                                               
      // recompute HT and MH
     
      if( simTauJetPt > htJetPtMin() && fabs(muEta) < htJetEtaMax()) {
        simHt += simTauJetPt;
      }

      if( simTauJetPt > mhtJetPtMin() && fabs(simTauJetEta) < mhtJetEtaMax()) {
        simMhtX -= simTauJetPt*cos(simTauJetPhi);
        simMhtY -= simTauJetPt*sin(simTauJetPhi);
      }
      //recompute met                                                                                                                          
      simmetX -= simTauJetPt*cos(simTauJetPhi);
      simmetY -= simTauJetPt*sin(simTauJetPhi);


      const double simMht = sqrt( simMhtX*simMhtX + simMhtY*simMhtY );

      const double simmet = sqrt( simmetX*simmetX + simmetY*simmetY );
      const double simmetPhi = std::atan2(simmetY,simmetX);

      //recompute jetVec

      vector<TLorentzVector> combNJetVec;
      combNJetVec.clear();
      combNJetVec = combjet(selNJetVec, simNJetVec);
      int combNJet = AnaFunctions::countJets(combNJetVec, AnaConsts::pt30Eta24Arr);

      //recompute deltaphi

      std::vector<double> * deltaPhiVec = new std::vector<double>();
      (*deltaPhiVec) = AnaFunctions::calcDPhi(combNJetVec, simmetPhi, 3, AnaConsts::dphiArr);
      bool passdeltaPhi = true;
      if( deltaPhiVec->at(0) < AnaConsts::dPhi0_CUT || deltaPhiVec->at(1) < AnaConsts::dPhi1_CUT || deltaPhiVec->at(2) < AnaConsts::dPhi2_CUT){
	  passdeltaPhi = false;
	}

	//recompute bjet

	int cnt1CSVS = AnaFunctions::countCSVS(combNJetVec, tr.getVec<double>("recoJetsBtag_0"), AnaConsts::cutCSVS, AnaConsts::bTagArr);
	bool passbJets = true;
	if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cnt1CSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cnt1CSVS < AnaConsts::high_nJetsSelBtagged ) ) ){
	  passbJets = false;
	}
	
	//top tagger input
	TLorentzVector metLVec_pre; metLVec_pre.SetPtEtaPhiM(simmet, 0, simmetPhi, 0);
        int comb30_pre = AnaFunctions::countJets(combNJetVec, AnaConsts::pt30Arr);
	std::vector<TLorentzVector> *jetsLVec_forTagger_pre = new std::vector<TLorentzVector>(); std::vector<double> *recoJetsBtag_forTagger_pre = new std::vector<double>();
	AnaFunctions::prepareJetsForTagger(combNJetVec, tr.getVec<double>("recoJetsBtag_0"), (*jetsLVec_forTagger_pre), (*recoJetsBtag_forTagger_pre));

	int bestTopJetIdx_pre = -1;
	bool remainPassCSVS_pre = false;
	int pickedRemainingCombfatJetIdx_pre = -1;
	double bestTopJetMass_pre = -1;


      //Apply baseline cut

	if(combNJet<AnaConsts::nJetsSelPt30Eta24) continue;
      	if(simmet<AnaConsts::defaultMETcut) continue;
	if(!passdeltaPhi) continue;
	if(!passbJets) continue;
      	//Apply Top tagger
	if(comb30_pre >= AnaConsts::nJetsSel ){
	  type3Ptr->processEvent((*jetsLVec_forTagger_pre), (*recoJetsBtag_forTagger_pre), metLVec_pre);
	  bestTopJetIdx_pre = type3Ptr->bestTopJetIdx;
	  remainPassCSVS_pre = type3Ptr->remainPassCSVS;
	  pickedRemainingCombfatJetIdx_pre = type3Ptr->pickedRemainingCombfatJetIdx;
	  if( bestTopJetIdx_pre != -1 ) bestTopJetMass_pre = type3Ptr->bestTopJetLVec.M();
	}

	

		bool passTopTagger = true;
	//bestTopJetIdx_pre != -1 means at least 1 top candidate!
	if( bestTopJetIdx_pre == -1 ){passTopTagger = false; }
	if( ! remainPassCSVS_pre ){passTopTagger = false; }
	if( pickedRemainingCombfatJetIdx_pre == -1 && jetsLVec_forTagger_pre->size()>=6 ){passTopTagger = false; }
	if( ! (bestTopJetMass_pre > AnaConsts::lowTopCut_ && bestTopJetMass_pre < AnaConsts::highTopCut_ ) ){ passTopTagger = false; }

	if(!passTopTagger) continue;
      
      //correction factor:
      const double corrBRWToTauHad = 0.65;  // Correction for the BR of hadronic tau decays                          
      const double corrBRTauToMu = Efficiency::taumucor(Efficiency::Ptbin1(muPt));//correction from tauonic mu contamination
      const double corrMuAcc = 1./Efficiency::acc(); // Correction for muon acceptance                                                        
      const double corrMuRecoEff = 1./Efficiency::reco(Efficiency::Ptbin(muPt)); // Correction for muon reconstruction efficiency             
      const double corrMuIsoEff = 1./Efficiency::iso(Efficiency::Ptbin(muPt)); // Correction for muon isolation efficiency 

      //The overall correction factor                                                                                                          
      const double corr = corrBRTauToMu * corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff;

      // Fill the prediction
      myBaseHistgram.hPredHt->Fill(simHt,corr);
      myBaseHistgram.hPredmet->Fill(simmet,corr);
      myBaseHistgram.hPredNJets->Fill(combNJet,corr);

    }//control sample loop

      //Expectation part 
      // Select only events where the W decayed into a hadronically decaying tau

    bool passBaseline = tr.getVar<bool>("passBaseline");
    bool passLeptVeto = tr.getVar<bool>("passLeptVeto");
    bool passnJets = tr.getVar<bool>("passnJets");
    bool passMET = tr.getVar<bool>("passMET");
    bool passdPhis = tr.getVar<bool>("passdPhis");
    bool passBJets = tr.getVar<bool>("passBJets");
    bool passTagger = tr.getVar<bool>("passTagger");

    int nJets = tr.getVar<int>("cntNJetsPt30Eta24");

    if(W_tau_prongsVec.size()==0)continue;

    if(!passLeptVeto)continue;
    if(!passnJets)continue;
    if(!passMET)continue;
    if(!passdPhis)continue;
    if(!passBJets)continue;
    if(!passTagger) continue;

    myBaseHistgram.hTrueHt->Fill(ht);
    myBaseHistgram.hTruemet->Fill(met);
    myBaseHistgram.hTrueNJets->Fill(nJets);
  }
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();

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
