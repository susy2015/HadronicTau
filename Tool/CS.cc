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
#include "SusyAnaTools/Tools/ISRCorrector.h"
#include "TStopwatch.h"
#include "TString.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"
#include "CS.h"
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

const bool doISR = true;
const bool dobSF = true;
const bool dolepSF = true;

TFile * bTagEffFile =0;

TH2D * mu_mediumID_SF = 0, * mu_miniISO_SF = 0;
TH1D * mu_trkptGT10_SF = 0, * mu_trkptLT10_SF = 0;

TH2D * ele_VetoID_SF = 0, * ele_miniISO_SF = 0;
TH2D * ele_trkpt_SF = 0;

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
  {
    std::cerr <<"Please give 5 arguments "<<"SubsampleName"<<" Input Template" <<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
    std::cerr <<" Valid configurations are " << std::endl;
    std::cerr <<" ./CS TTbarInc 1000 0 1" << std::endl;
    return -1;
  }
  const char *subsamplename = argv[1];
  const char *Maxevent = argv[2];
  const  char *Stratfile = argv[3];
  const  char *Filerun = argv[4];
  const  int startfile = std::atoi(Stratfile);
  const int filerun = std::atoi(Filerun);
  bool isData = false;
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
  if(sampleString.Contains("Data")){Lumiscale = 1.0; isData = true;}
  //Searchbin                                                                                                                                                                                    
  SearchBins SB("SB_v1_2017");

//  const int nTotBins = SB.nSearchBins();
  std::vector<int> cachedMT2binIdx_mu_3DVec, cachedHTbinIdx_mu_3DVec;
  std::vector<TH1D*> cachedMT2hist_mu_3DVec, cachedHThist_mu_3DVec;
  std::vector<int> cachedMT2binIdx_el_3DVec, cachedHTbinIdx_el_3DVec;
  std::vector<TH1D*> cachedMT2hist_el_3DVec, cachedHThist_el_3DVec;

  std::vector<int> cachedMT2binIdx_mu_2DVec, cachedHTbinIdx_mu_2DVec;
  std::vector<TH1D*> cachedMT2hist_mu_2DVec, cachedHThist_mu_2DVec, cachedmethist_mu_2DVec;
  std::vector<int> cachedMT2binIdx_el_2DVec, cachedHTbinIdx_el_2DVec;
  std::vector<TH1D*> cachedMT2hist_el_2DVec, cachedHThist_el_2DVec, cachedmethist_el_2DVec;

  std::string spec = "CS";

  AnaFunctions::prepareForNtupleReader();
  //AnaFunctions::prepareTopTagger();
  NTupleReader *tr =0;
  if( isData ) tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames_DataOnly);
  else tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames);
  CSBaseline = new BaselineVessel(*tr, spec);
  CSBaseline->SetupTopTagger(true,"TopTagger.cfg");
  tr->registerFunction((*CSBaseline));

  if( !isData ){
    // Lepton SF
    TFile * allINone_leptonSF_file = new TFile("allINone_leptonSF.root");
    if( !allINone_leptonSF_file->IsZombie() ){
      mu_mediumID_SF = (TH2D*) allINone_leptonSF_file->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0");
      mu_miniISO_SF = (TH2D*) allINone_leptonSF_file->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass");
      mu_trkptGT10_SF = (TH1D*) allINone_leptonSF_file->Get("mutrksfptg10");
      mu_trkptLT10_SF = (TH1D*) allINone_leptonSF_file->Get("mutrksfptl10");
  
      ele_VetoID_SF = (TH2D*) allINone_leptonSF_file->Get("GsfElectronToVeto");
      ele_miniISO_SF = (TH2D*) allINone_leptonSF_file->Get("MVAVLooseElectronToMini");
      ele_trkpt_SF = (TH2D*) allINone_leptonSF_file->Get("EGamma_SF2D");
    }

    //BTag SF
    BTagCorrector btagcorr;
    if( bTagEffFile ) delete bTagEffFile;
    if( std::string(subsamplename).find("TTbar") != std::string::npos )
    {
       bTagEffFile = new TFile("TTbarNoHad_bTagEff.root");
    }else if( std::string(subsamplename).find("WJetsToLNu") != std::string::npos )
    {
       bTagEffFile = new TFile("WJetsToLNu_HT_bTagEff.root"); 
    }else if( std::string(subsamplename).find("tW") != std::string::npos )
    {
       bTagEffFile = new TFile("SingleTop_bTagEff.root");
    }else
    {
      std::cout<<"Error ... not supported subsamplename : "<<subsamplename<<std::endl;
      return 0; 
    }
    btagcorr.SetEffs(bTagEffFile); 
    tr->registerFunction(btagcorr);
  
    //ISR Correcotr
    ISRCorrector *isrcorr;
    if( std::string(subsamplename).find("TTbar") != std::string::npos )
    {
      isrcorr = new ISRCorrector("TTbarNoHad_NJetsISR.root", "", "");
    }else if( std::string(subsamplename).find("WJetsToLNu") != std::string::npos )
    {
      isrcorr = new ISRCorrector("WJetsToLNu_HT_NJetsISR.root", "", "");
    }else if( std::string(subsamplename).find("tW") != std::string::npos )
    {
      isrcorr = new ISRCorrector("SingleTop_NJetsISR.root", "", "");
    }else
    {
      std::cout<<"Error ... not supported subsamplename : "<<subsamplename<<std::endl;
      return 0; 
    }
    tr->registerFunction((*isrcorr));
  }

  //Add cleanJets function
  //stopFunctions::cleanJets.setMuonIso("mini");
  //stopFunctions::cleanJets.setElecIso("mini");
  //stopFunctions::cleanJets.setRemove(false);
  //tr->registerFunction(&stopFunctions::cleanJets);
 
  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr->getNEntries();
  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<Lumiscale<<std::endl; 
  cout<<"maxevent: "<<maxevent<<endl;
  int k =0;
  // Loop over the events (tree entries)
  while(tr->getNextEvent())
  {
    k++;
    //std::cout<<"evt: "<<k<<std::endl;
    if(maxevent>=0 && tr->getEvtNum() > maxevent ) break;
    // Add print out of the progress of looping
    if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;

    const double isrWght = (doISR && !isData)? tr->getVar<double>("isr_Unc_Cent") : 1.0;
    const double bSF = (dobSF && !isData)? tr->getVar<double>("bTagSF_EventWeightSimple_Central") : 1.0;
    const double bSF_up = (dobSF && !isData)? tr->getVar<double>("bTagSF_EventWeightSimple_Up") : 1.0;
    const double bSF_down = (dobSF && !isData)? tr->getVar<double>("bTagSF_EventWeightSimple_Down") : 1.0;

    const vector<TLorentzVector> &muonsLVec = tr->getVec<TLorentzVector>("muonsLVec");
    const vector<double> &muonsRelIso = tr->getVec<double>("muonsRelIso");
    const vector<double> &muonsMiniIso = tr->getVec<double>("muonsMiniIso");
    const vector<double> &muonsMtw = tr->getVec<double>("muonsMtw");    
    const vector<TLorentzVector> &elesLVec = tr->getVec<TLorentzVector>("elesLVec");
    const vector<double> &elesRelIso = tr->getVec<double>("elesRelIso");
    const vector<double> &elesMiniIso = tr->getVec<double>("elesMiniIso");
    const vector<double> &elesMtw = tr->getVec<double>("elesMtw");
    const vector<unsigned int> &elesisEB = tr->getVec<unsigned int>("elesisEB");    
    const vector<TLorentzVector> &jetsLVec = tr->getVec<TLorentzVector>("jetsLVec");
    const vector<double> &recoJetsBtag_0 = tr->getVec<double>("recoJetsBtag_0");
    const vector<int> &looseisoTrksMatchedJetIdx = tr->getVec<int>("looseisoTrksMatchedJetIdx");
    const vector<TLorentzVector> &loose_isoTrksLVec = tr->getVec<TLorentzVector>("loose_isoTrksLVec");
    const vector<double> &loose_isoTrks_iso = tr->getVec<double>("loose_isoTrks_iso");
    const vector<double> &loose_isoTrks_mtw = tr->getVec<double>("loose_isoTrks_mtw");
    const vector<int> &loose_isoTrks_pdgId = tr->getVec<int>("loose_isoTrks_pdgId");
    const std::vector<double> muonspfActivity = tr->getVec<double>("muonspfActivity");
    const std::vector<int> & muonsFlagIDVec = tr->getVec<int>("muonsFlagMedium");
    const std::vector<int> & elesFlagIDVec = tr->getVec<int>("elesFlagVeto");
    const std::vector<double>& recoJetschargedHadronEnergyFraction = tr->getVec<double>("recoJetschargedHadronEnergyFraction");
    const std::vector<double>& recoJetschargedEmEnergyFraction = tr->getVec<double>("recoJetschargedEmEnergyFraction");
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
    const double HT = tr->getVar<double>("HT"+spec);
    double met=tr->getVar<double>("met");
    double metphi=tr->getVar<double>("metphi");
    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);
    //mht calculation
    TLorentzVector MhtLVec;	      
    for(unsigned int ij=0; ij<jetsLVec.size(); ij++)
    {
      if( !AnaFunctions::jetPassCuts(jetsLVec[ij], AnaConsts::pt30Arr) ) continue;
      MhtLVec -= jetsLVec[ij];
    }
    const double Mht = MhtLVec.Pt();

    bool passBaselineCS = passnJets && passdPhis && passMET && passBJets && passTagger && passHT && passMT2;
    bool passTrigger = true;        
    if(isData)
    {
      bool foundTrigger = false;
      for(unsigned it=0; it<TriggerNames.size(); it++)
      {
        if( sampleString.Contains("MET") )
        {
          if( TriggerNames[it].find("HLT_PFMET170_NoiseCleaned_v")!= string::npos || TriggerNames[it].find("HLT_PFMET170_JetIdCleaned_v") != string::npos || TriggerNames[it].find("HLT_PFMET170_HBHECleaned_v") != string::npos || TriggerNames[it].find("HLT_PFMET100_PFMHT100_IDTight_v") != string::npos || TriggerNames[it].find("HLT_PFMET110_PFMHT110_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMET120_PFMHT120_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMET130_PFMHT130_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMET140_PFMHT140_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMET150_PFMHT150_IDTight_v")!= string::npos)
          {
            if( PassTrigger[it] ) foundTrigger = true;
	  }
        }
      }
      if( !foundTrigger ) passTrigger = false;
    }
    if(!passTrigger) continue;
    // Muon Control sample
    // The kinematic properties of the well-reconstructed, isolated muon                                                         
    vector<TLorentzVector> isomuonsLVec;
    vector<int> isomuonsIdxVec;
    for(unsigned int im=0; im<muonsLVec.size(); im++)
    {
      if( AnaFunctions::passMuon(muonsLVec.at(im), muonsMiniIso.at(im), muonsMtw.at(im), muonsFlagIDVec.at(im), AnaConsts::muonsMiniIsoArr))
      {
        isomuonsLVec.push_back(muonsLVec.at(im)); isomuonsIdxVec.push_back(im);
      }
    }
    // Require one and only one muon                                                                                         
    // Veto events with additional electrons (same veto criteria as baseline for electrons)                           
    if( nMuons == 1 && nElectrons == AnaConsts::nElectronsSel )
    {
      if( nMuons != isomuonsLVec.size() )
      {
        std::cout<<"ERROR ... mis-matching between veto muon and selected muon! Skipping..."<<std::endl; continue; 
      }

      //mtW cut
      const TLorentzVector muLVec = isomuonsLVec.at(0);
      const double mtw = calcMT(muLVec, metLVec);
      const bool pass_mtw = mtw<100 ? true : false;

      const double eta = muLVec.Eta(), pt = muLVec.Pt();
      const double abseta = std::abs(eta);
      double mu_id_SF = 1.0, mu_iso_SF = 1.0, mu_trk_SF = 1.0;
      if( dolepSF && !isData )
      {
        if( mu_mediumID_SF ){ mu_id_SF = mu_mediumID_SF->GetBinContent(mu_mediumID_SF->FindBin(pt, abseta)); if( mu_id_SF == 0 ) mu_id_SF = 1.0; } // very simple way dealing with out of range issue of the TH2D
        if( mu_miniISO_SF ){ mu_iso_SF = mu_miniISO_SF->GetBinContent(mu_miniISO_SF->FindBin(pt, abseta)); if( mu_iso_SF == 0 ) mu_iso_SF = 1.0; }
        if( pt < 10 && mu_trkptLT10_SF ){ mu_trk_SF = mu_trkptLT10_SF->GetBinContent(mu_trkptLT10_SF->FindBin(eta)); if( mu_trk_SF == 0 ) mu_trk_SF = 1.0; }
        if( pt >= 10 && mu_trkptGT10_SF ){ mu_trk_SF = mu_trkptGT10_SF->GetBinContent(mu_trkptGT10_SF->FindBin(eta)); if( mu_trk_SF == 0 ) mu_trk_SF = 1.0; }
      }
      const double mu_SF = mu_id_SF * mu_iso_SF * mu_trk_SF;

      const double corr_SF = bSF * isrWght * mu_SF;
      const double corr_bSF_up = bSF_up * isrWght * mu_SF;
      const double corr_bSF_down = bSF_down * isrWght * mu_SF;
  
      //Dist.
      if(passBaselineCS && passNoiseEventFilter && pass_mtw)
      {
        int jSR = SB.find_Binning_Index(nbJets, nTops, MT2, met, HT);
        if( jSR!= -1 )
        {
          myBaseHistgram.hYields_mu->Fill(jSR, Lumiscale*corr_SF);
	  //bSF systematics
          myBaseHistgram.hYields_mu_bSFup->Fill(jSR, Lumiscale*corr_bSF_up);
	  myBaseHistgram.hYields_mu_bSFdown->Fill(jSR, Lumiscale*corr_bSF_down);
        }
  	  
        FillDouble(myBaseHistgram.hMET_mu, met, Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hMT2_mu, MT2, Lumiscale*corr_SF);
        FillInt(myBaseHistgram.hNbJets_mu, nbJets, Lumiscale*corr_SF);
        FillInt(myBaseHistgram.hNTops_mu, nTops, Lumiscale*corr_SF);	
  	  
        FillInt(myBaseHistgram.hNJets_mu, nJets, Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hHT_mu, HT, Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hdPhi0_mu, dPhiVec[0], Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hdPhi1_mu, dPhiVec[1], Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hdPhi2_mu, dPhiVec[2], Lumiscale*corr_SF);
  
        if( nbJets <=2 && nTops<=2 )
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, 250, met, HT); // use the lowest MT2 bin in (nb, ntop, met) to collapse the MT2 bins
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          auto it = std::find(cachedMT2binIdx_mu_3DVec.begin(), cachedMT2binIdx_mu_3DVec.end(), pseudo_SR);
          if( it != cachedMT2binIdx_mu_3DVec.end() )
          {
            auto index = std::distance(cachedMT2binIdx_mu_3DVec.begin(), it);
            FillDouble(cachedMT2hist_mu_3DVec[index], MT2, Lumiscale*corr_SF);
          }else
          {
            cachedMT2binIdx_mu_3DVec.push_back(pseudo_SR);
            char tmp_str[200];
            if( pseudo_binDef.met_hi_ != -1 )
            {
              sprintf(tmp_str, "muCS_MT2_3D_nb%d_nt%d_%3.0fmet%3.0f", nbJets, nTops, pseudo_binDef.met_lo_, pseudo_binDef.met_hi_);
            }else
            {
              sprintf(tmp_str, "muCS_MT2_3D_nb%d_nt%d_%3.0fmetInf", nbJets, nTops, pseudo_binDef.met_lo_);
            }
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 28, 200, 900);
            cachedMT2hist_mu_3DVec.push_back(h1_tmp);
            FillDouble(cachedMT2hist_mu_3DVec.back(), MT2, Lumiscale*corr_SF);
          }
        }else
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, MT2, met, 350); // use the lowest HT bin in (nb, ntop, met) to collapse the HT bins
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          const int nbJetsCopy = nbJets >=3 ? 3 : nbJets; const int nTopsCopy = nTops >=3 ? 3 : nTops;
          auto it = std::find(cachedHTbinIdx_mu_3DVec.begin(), cachedHTbinIdx_mu_3DVec.end(), pseudo_SR);
          if( it != cachedHTbinIdx_mu_3DVec.end() )
          {
            auto index = std::distance(cachedHTbinIdx_mu_3DVec.begin(), it);
            FillDouble(cachedHThist_mu_3DVec[index], HT, Lumiscale*corr_SF);
          }else
          {
            cachedHTbinIdx_mu_3DVec.push_back(pseudo_SR);
            char tmp_str[200];
            if( pseudo_binDef.met_hi_ != -1 )
            {
              sprintf(tmp_str, "muCS_HT_3D_nb%d_nt%d_%3.0fmet%3.0f", nbJetsCopy, nTopsCopy, pseudo_binDef.met_lo_, pseudo_binDef.met_hi_);
            }else
            {
              sprintf(tmp_str, "muCS_HT_3D_nb%d_nt%d_%3.0fmetInf", nbJetsCopy, nTopsCopy, pseudo_binDef.met_lo_);
            }
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 34, 300, 2000);
            cachedHThist_mu_3DVec.push_back(h1_tmp);
            FillDouble(cachedHThist_mu_3DVec.back(), HT, Lumiscale*corr_SF);
          }
        }
  
  // 2D in (nb, nt)
        if( nbJets <=2 && nTops<=2 )
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, 250, 300, HT); // use the lowest (MT2, met) bin in (nb, ntop) to collapse the MT2 and met
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          auto it = std::find(cachedMT2binIdx_mu_2DVec.begin(), cachedMT2binIdx_mu_2DVec.end(), pseudo_SR);
          if( it != cachedMT2binIdx_mu_2DVec.end() )
          {
            auto index = std::distance(cachedMT2binIdx_mu_2DVec.begin(), it);
            FillDouble(cachedMT2hist_mu_2DVec[index], MT2, Lumiscale*corr_SF);
            FillDouble(cachedmethist_mu_2DVec[index], met, Lumiscale*corr_SF);
          }else
          {
            cachedMT2binIdx_mu_2DVec.push_back(pseudo_SR);
            char tmp_str[200];
            sprintf(tmp_str, "muCS_MT2_2D_nb%d_nt%d", nbJets, nTops);
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 28, 200, 900);
            cachedMT2hist_mu_2DVec.push_back(h1_tmp);
            FillDouble(cachedMT2hist_mu_2DVec.back(), MT2, Lumiscale*corr_SF);
  
            sprintf(tmp_str, "muCS_met_2D_nb%d_nt%d", nbJets, nTops);
            TH1D * h1_tmp2 = new TH1D(tmp_str, tmp_str, 24, 250, 850);
            cachedmethist_mu_2DVec.push_back(h1_tmp2);
            FillDouble(cachedmethist_mu_2DVec.back(), met, Lumiscale*corr_SF);
          }
        }else
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, MT2, 300, 350); // use the lowest (HT, met) bin in (nb, ntop) to collapse the HT and met
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          const int nbJetsCopy = nbJets >=3 ? 3 : nbJets; const int nTopsCopy = nTops >=3 ? 3 : nTops;
          auto it = std::find(cachedHTbinIdx_mu_2DVec.begin(), cachedHTbinIdx_mu_2DVec.end(), pseudo_SR);
          if( it != cachedHTbinIdx_mu_2DVec.end() )
          {
            auto index = std::distance(cachedHTbinIdx_mu_2DVec.begin(), it);
            FillDouble(cachedHThist_mu_2DVec[index], HT, Lumiscale*corr_SF);
          }else
          {
            cachedHTbinIdx_mu_2DVec.push_back(pseudo_SR);
            char tmp_str[200];
            sprintf(tmp_str, "muCS_HT_2D_nb%d_nt%d", nbJetsCopy, nTopsCopy);
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 34, 300, 2000);
            cachedHThist_mu_2DVec.push_back(h1_tmp);
            FillDouble(cachedHThist_mu_2DVec.back(), HT, Lumiscale*corr_SF);
          }
        }
      }
    }//end of muon CS

    // Electron Control sample
    //The kinematic properties of the well-reconstructed, isolated muon                                                         
    vector<TLorentzVector> isoelesLVec;
    vector<int> isoelesIdxVec;
    for(unsigned int im=0; im<elesLVec.size(); im++)
    {
      if( AnaFunctions::passElectron(elesLVec.at(im), elesMiniIso.at(im), elesMtw.at(im), elesisEB.at(im), elesFlagIDVec.at(im), AnaConsts::elesMiniIsoArr))
      {
        isoelesLVec.push_back(elesLVec.at(im)); isoelesIdxVec.push_back(im); 
      }
    }
    // Require one and only one electron                                                                                         
    // Veto events with additional muons (same veto criteria as baseline for muons)                           
    if( nElectrons == 1 && nMuons == AnaConsts::nMuonsSel )
    {
      if( nElectrons != isoelesLVec.size() )
      { 
        std::cout<<"ERROR ... mis-matching between veto ele and selected ele! Skipping..."<<std::endl; continue; 
      }
  
      //mtW cut
      const TLorentzVector eleLVec = isoelesLVec.at(0);
      const double elemtw = calcMT(eleLVec, metLVec);
      const bool pass_mtwele = elemtw<100 ? true : false;
  
      const double eta = eleLVec.Eta(), pt = eleLVec.Pt();
      const double abseta = std::abs(eta);
      double ele_id_SF = 1.0, ele_iso_SF = 1.0, ele_trk_SF = 1.0;
      if( dolepSF && !isData )
      {
        if( ele_VetoID_SF ){ ele_id_SF = ele_VetoID_SF->GetBinContent(ele_VetoID_SF->FindBin(pt, abseta)); if( ele_id_SF == 0 ) ele_id_SF = 1.0; } // very simple way dealing with out of range issue of the TH2D
        if( ele_miniISO_SF ){ ele_iso_SF = ele_miniISO_SF->GetBinContent(ele_miniISO_SF->FindBin(pt, abseta)); if( ele_iso_SF == 0 ) ele_iso_SF = 1.0; }
        if( ele_trkpt_SF ){ ele_trk_SF = ele_trkpt_SF->GetBinContent(ele_trkpt_SF->FindBin(eta, pt)); if( ele_trk_SF == 0 ) ele_trk_SF = 1.0; }
      }
      const double ele_SF = ele_id_SF * ele_iso_SF * ele_trk_SF;

      const double corr_SF = bSF * isrWght * ele_SF;
      const double corr_bSF_up = bSF_up * isrWght * ele_SF;
      const double corr_bSF_down = bSF_down * isrWght * ele_SF;
      //Dist.
      if(passBaselineCS && passNoiseEventFilter && pass_mtwele)
      {
        int kSR = SB.find_Binning_Index(nbJets, nTops, MT2, met, HT);
        if( kSR!= -1 )
        {
          myBaseHistgram.hYields_el->Fill(kSR, Lumiscale*corr_SF);
	  //bSF systematics
          myBaseHistgram.hYields_el_bSFup->Fill(kSR, Lumiscale*corr_bSF_up);
	  myBaseHistgram.hYields_el_bSFdown->Fill(kSR, Lumiscale*corr_bSF_down);
        }
  	  
        FillDouble(myBaseHistgram.hMET_el, met, Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hMT2_el, MT2, Lumiscale*corr_SF);
        FillInt(myBaseHistgram.hNbJets_el, nbJets, Lumiscale*corr_SF);
        FillInt(myBaseHistgram.hNTops_el, nTops, Lumiscale*corr_SF);	
  	  
        FillInt(myBaseHistgram.hNJets_el, nJets, Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hHT_el, HT, Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hdPhi0_el, dPhiVec[0], Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hdPhi1_el, dPhiVec[1], Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hdPhi2_el, dPhiVec[2], Lumiscale*corr_SF);

        if( nbJets <=2 && nTops<=2 )
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, 250, met, HT); // use the lowest MT2 bin in (nb, ntop, met) to collapse the MT2 bins
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          auto it = std::find(cachedMT2binIdx_el_3DVec.begin(), cachedMT2binIdx_el_3DVec.end(), pseudo_SR);
          if( it != cachedMT2binIdx_el_3DVec.end() )
          {
            auto index = std::distance(cachedMT2binIdx_el_3DVec.begin(), it);
            FillDouble(cachedMT2hist_el_3DVec[index], MT2, Lumiscale*corr_SF);
          }else
          {
            cachedMT2binIdx_el_3DVec.push_back(pseudo_SR);
            char tmp_str[200];
            if( pseudo_binDef.met_hi_ != -1 )
            {
              sprintf(tmp_str, "eleCS_MT2_3D_nb%d_nt%d_%3.0fmet%3.0f", nbJets, nTops, pseudo_binDef.met_lo_, pseudo_binDef.met_hi_);
            }else
            {
              sprintf(tmp_str, "eleCS_MT2_3D_nb%d_nt%d_%3.0fmetInf", nbJets, nTops, pseudo_binDef.met_lo_);
            }
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 28, 200, 900);
            cachedMT2hist_el_3DVec.push_back(h1_tmp);
            FillDouble(cachedMT2hist_el_3DVec.back(), MT2, Lumiscale*corr_SF);
          }
        }else
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, MT2, met, 350); // use the lowest HT bin in (nb, ntop, met) to collapse the HT bins
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          const int nbJetsCopy = nbJets >=3 ? 3 : nbJets; const int nTopsCopy = nTops >=3 ? 3 : nTops;
          auto it = std::find(cachedHTbinIdx_el_3DVec.begin(), cachedHTbinIdx_el_3DVec.end(), pseudo_SR);
          if( it != cachedHTbinIdx_el_3DVec.end() )
          {
            auto index = std::distance(cachedHTbinIdx_el_3DVec.begin(), it);
            FillDouble(cachedHThist_el_3DVec[index], HT, Lumiscale*corr_SF);
          }else
          {
            cachedHTbinIdx_el_3DVec.push_back(pseudo_SR);
            char tmp_str[200];
            if( pseudo_binDef.met_hi_ != -1 )
            {
              sprintf(tmp_str, "eleCS_HT_3D_nb%d_nt%d_%3.0fmet%3.0f", nbJetsCopy, nTopsCopy, pseudo_binDef.met_lo_, pseudo_binDef.met_hi_);
            }else
            {
              sprintf(tmp_str, "eleCS_HT_3D_nb%d_nt%d_%3.0fmetInf", nbJetsCopy, nTopsCopy, pseudo_binDef.met_lo_);
            }
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 34, 300, 2000);
            cachedHThist_el_3DVec.push_back(h1_tmp);
            FillDouble(cachedHThist_el_3DVec.back(), HT, Lumiscale*corr_SF);
          }
        }
  
  // 2D in (nb, nt)
        if( nbJets <=2 && nTops<=2 )
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, 250, 300, HT); // use the lowest (MT2, met) bin in (nb, ntop) to collapse the MT2 and met
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          auto it = std::find(cachedMT2binIdx_el_2DVec.begin(), cachedMT2binIdx_el_2DVec.end(), pseudo_SR);
          if( it != cachedMT2binIdx_el_2DVec.end() )
          {
            auto index = std::distance(cachedMT2binIdx_el_2DVec.begin(), it);
            FillDouble(cachedMT2hist_el_2DVec[index], MT2, Lumiscale*corr_SF);
            FillDouble(cachedmethist_el_2DVec[index], met, Lumiscale*corr_SF);
          }else
          {
            cachedMT2binIdx_el_2DVec.push_back(pseudo_SR);
            char tmp_str[200];
            sprintf(tmp_str, "eleCS_MT2_2D_nb%d_nt%d", nbJets, nTops);
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 28, 200, 900);
            cachedMT2hist_el_2DVec.push_back(h1_tmp);
            FillDouble(cachedMT2hist_el_2DVec.back(), MT2, Lumiscale*corr_SF);
  
            sprintf(tmp_str, "eleCS_met_2D_nb%d_nt%d", nbJets, nTops);
            TH1D * h1_tmp2 = new TH1D(tmp_str, tmp_str, 24, 250, 850);
            cachedmethist_el_2DVec.push_back(h1_tmp2);
            FillDouble(cachedmethist_el_2DVec.back(), met, Lumiscale*corr_SF);
          }
        }else
        {
          int pseudo_SR = SB.find_Binning_Index(nbJets, nTops, MT2, 300, 350); // use the lowest (HT, met) bin in (nb, ntop) to collapse the HT and met
          SearchBins::searchBinDef pseudo_binDef; SB.find_BinBoundaries(pseudo_SR, pseudo_binDef);
          const int nbJetsCopy = nbJets >=3 ? 3 : nbJets; const int nTopsCopy = nTops >=3 ? 3 : nTops;
          auto it = std::find(cachedHTbinIdx_el_2DVec.begin(), cachedHTbinIdx_el_2DVec.end(), pseudo_SR);
          if( it != cachedHTbinIdx_el_2DVec.end() )
          {
            auto index = std::distance(cachedHTbinIdx_el_2DVec.begin(), it);
            FillDouble(cachedHThist_el_2DVec[index], HT, Lumiscale*corr_SF);
          }else
          {
            cachedHTbinIdx_el_2DVec.push_back(pseudo_SR);
            char tmp_str[200];
            sprintf(tmp_str, "eleCS_HT_2D_nb%d_nt%d", nbJetsCopy, nTopsCopy);
            TH1D * h1_tmp = new TH1D(tmp_str, tmp_str, 34, 300, 2000);
            cachedHThist_el_2DVec.push_back(h1_tmp);
            FillDouble(cachedHThist_el_2DVec.back(), HT, Lumiscale*corr_SF);
          }
        }

      }
    }//end of electron CS
  }//event loop
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();
  fChain->Reset();  
  return 0;
}
