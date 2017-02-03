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

#include "SusyAnaTools/Tools/ISRCorrector.h"

#include "PrepSystematics.h"
#include "CommonShared.h"

using namespace std;

const bool doISR = true;
const bool dobSF = true;

TFile * bTagEffFile =0;

TH2D * mu_mediumID_SF = 0, * mu_miniISO_SF = 0;
TH1D * mu_trkptGT10_SF = 0, * mu_trkptLT10_SF = 0;

TH2D * ele_VetoID_SF = 0, * ele_miniISO_SF = 0;
TH2D * ele_trkpt_SF = 0;

// === Main Function ===================================================
int main(int argc, char* argv[])
{
  if (argc < 5)
  {
    std::cerr <<"Please give 5 arguments "<<"SubsampleName"<<" Input Template" <<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
    std::cerr <<" Valid configurations are " << std::endl;
    std::cerr <<" ./HadTauMC TTbarInc 1000 0 1" << std::endl;
    return -1;
  }
  const char *subsamplename = argv[1];
  const char *Maxevent = argv[2];
  const char *Stratfile = argv[3];
  const char *Filerun = argv[4];
  const int startfile = std::atoi(Stratfile);
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
  
  //Use BaselineVessel class for baseline variables and selections
  std::string spec = "MCExp";
  AnaFunctions::prepareForNtupleReader();
  NTupleReader *tr =0;
  tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames);
  ExpBaseline = new BaselineVessel(*tr, spec);
  ExpBaseline->SetupTopTagger(true,"TopTagger.cfg");
  tr->registerFunction((*ExpBaseline));
 
  //Assume ONLY run this on MC samples!
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

  systBaseline::pdfScale = new PDFUncertainty();
  
  systBaseline::SRblv_jecUp = new BaselineVessel(*tr, systBaseline::spec_jecUp);
  systBaseline::SRblv_jecDn = new BaselineVessel(*tr, systBaseline::spec_jecDn);
  systBaseline::SRblv_metMagUp = new BaselineVessel(*tr, systBaseline::spec_metMagUp);
  systBaseline::SRblv_metMagDn = new BaselineVessel(*tr, systBaseline::spec_metMagDn);
  systBaseline::SRblv_metPhiUp = new BaselineVessel(*tr, systBaseline::spec_metPhiUp);
  systBaseline::SRblv_metPhiDn = new BaselineVessel(*tr, systBaseline::spec_metPhiDn);
  
  tr->registerFunction((*systBaseline::pdfScale));
  
  PrepSystematics sysPrep;
    
  tr->registerFunction(sysPrep);
  
  tr->registerFunction((*systBaseline::SRblv_jecUp));
  tr->registerFunction((*systBaseline::SRblv_jecDn));
  tr->registerFunction((*systBaseline::SRblv_metMagUp));
  tr->registerFunction((*systBaseline::SRblv_metMagDn));
  tr->registerFunction((*systBaseline::SRblv_metPhiUp));
  tr->registerFunction((*systBaseline::SRblv_metPhiDn));

  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr->getNEntries();
  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<Lumiscale<<std::endl; 
  cout<<"maxevent: "<<maxevent<<endl;
  // Loop over the events (tree entries)
  int k = 0;
  while(tr->getNextEvent())
  {
    k++;
    if(maxevent>=0 && tr->getEvtNum() > maxevent ) break;
    // Add print out of the progress of looping
    if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;

    const double isrWght = doISR? tr->getVar<double>("isr_Unc_Cent") : 1.0;
    const double bSF = dobSF? tr->getVar<double>("bTagSF_EventWeightSimple_Central") : 1.0;
    const double bSF_up = dobSF? tr->getVar<double>("bTagSF_EventWeightSimple_Up") : 1.0;
    const double bSF_down = dobSF? tr->getVar<double>("bTagSF_EventWeightSimple_Down") : 1.0;

    const double isr_up = doISR? tr->getVar<double>("isr_Unc_Up") : 1.0;
    const double isr_down = doISR? tr->getVar<double>("isr_Unc_Down") : 1.0;

    const double NNPDF_From_Median_Central = tr->getVar<double>("NNPDF_From_Median_Central");
    const double NNPDF_From_Median_Up = tr->getVar<double>("NNPDF_From_Median_Up");
    const double NNPDF_From_Median_Down = tr->getVar<double>("NNPDF_From_Median_Down");

    const double Scaled_Variations_Up = tr->getVar<double>("Scaled_Variations_Up");
    const double Scaled_Variations_Down = tr->getVar<double>("Scaled_Variations_Down");

    const double genHT = tr->getVar<double>("genHT");

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

    const std::vector<TLorentzVector> & muonsLVec = tr->getVec<TLorentzVector>("muonsLVec");
    const std::vector<double> & muonsMiniIso = tr->getVec<double>("muonsMiniIso");
    const std::vector<double> & muonsMtw = tr->getVec<double>("muonsMtw");
    const std::vector<int> & muonsFlagMedium = tr->getVec<int>("muonsFlagMedium");

    const std::vector<TLorentzVector> & elesLVec = tr->getVec<TLorentzVector>("elesLVec");
    const std::vector<double> & elesMiniIso = tr->getVec<double>("elesMiniIso");
    const std::vector<double> & elesMtw = tr->getVec<double>("elesMtw");
    const std::vector<unsigned int> & elesisEB = tr->getVec<unsigned int>("elesisEB");
    const std::vector<int> & elesFlagVeto = tr->getVec<int>("elesFlagVeto");

    const std::vector<TLorentzVector> & isoTrksLVec = tr->getVec<TLorentzVector>("loose_isoTrksLVec");
    const std::vector<double> & isoTrksIso = tr->getVec<double>("loose_isoTrks_iso");
    const std::vector<double> & isoTrksMtw = tr->getVec<double>("loose_isoTrks_mtw");
    const std::vector<int> & isoTrkspdgId = tr->getVec<int>("loose_isoTrks_pdgId");
    
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

    const double met_metMagUp = tr->getVar<double>("met_metMagUp");
    const double metphi_metMagUp = tr->getVar<double>("metphi");
    const double MT2_metMagUp = tr->getVar<double>("best_had_brJet_MT2" + systBaseline::spec_metMagUp);
    const int nbJets_metMagUp = tr->getVar<int>("cntCSVS" + systBaseline::spec_metMagUp);
    const int nTops_metMagUp = tr->getVar<int>("nTopCandSortedCnt" + systBaseline::spec_metMagUp);
    const double HT_metMagUp = tr->getVar<double>("HT" + systBaseline::spec_metMagUp);
    const bool passBaseline_metMagUp = tr->getVar<bool>("passBaseline" + systBaseline::spec_metMagUp);

    const double met_metMagDn = tr->getVar<double>("met_metMagDn");
    const double metphi_metMagDn = tr->getVar<double>("metphi");
    const double MT2_metMagDn = tr->getVar<double>("best_had_brJet_MT2" + systBaseline::spec_metMagDn);
    const int nbJets_metMagDn = tr->getVar<int>("cntCSVS" + systBaseline::spec_metMagDn);
    const int nTops_metMagDn = tr->getVar<int>("nTopCandSortedCnt" + systBaseline::spec_metMagDn);
    const double HT_metMagDn = tr->getVar<double>("HT" + systBaseline::spec_metMagDn);
    const bool passBaseline_metMagDn = tr->getVar<bool>("passBaseline" + systBaseline::spec_metMagDn);

    const double met_metPhiUp = tr->getVar<double>("met");
    const double metphi_metPhiUp = tr->getVar<double>("metphi_metPhiUp");
    const double MT2_metPhiUp = tr->getVar<double>("best_had_brJet_MT2" + systBaseline::spec_metPhiUp);
    const int nbJets_metPhiUp = tr->getVar<int>("cntCSVS" + systBaseline::spec_metPhiUp);
    const int nTops_metPhiUp = tr->getVar<int>("nTopCandSortedCnt" + systBaseline::spec_metPhiUp);
    const double HT_metPhiUp = tr->getVar<double>("HT" + systBaseline::spec_metPhiUp);
    const bool passBaseline_metPhiUp = tr->getVar<bool>("passBaseline" + systBaseline::spec_metPhiUp);

    const double met_metPhiDn = tr->getVar<double>("met");
    const double metphi_metPhiDn = tr->getVar<double>("metphi_metPhiDn");
    const double MT2_metPhiDn = tr->getVar<double>("best_had_brJet_MT2" + systBaseline::spec_metPhiDn);
    const int nbJets_metPhiDn = tr->getVar<int>("cntCSVS" + systBaseline::spec_metPhiDn);
    const int nTops_metPhiDn = tr->getVar<int>("nTopCandSortedCnt" + systBaseline::spec_metPhiDn);
    const double HT_metPhiDn = tr->getVar<double>("HT" + systBaseline::spec_metPhiDn);
    const bool passBaseline_metPhiDn = tr->getVar<bool>("passBaseline" + systBaseline::spec_metPhiDn);

    const double met_jecUp = tr->getVar<double>("met");
    const double metphi_jecUp = tr->getVar<double>("metphi");
    const double MT2_jecUp = tr->getVar<double>("best_had_brJet_MT2" + systBaseline::spec_jecUp);
    const int nbJets_jecUp = tr->getVar<int>("cntCSVS" + systBaseline::spec_jecUp);
    const int nTops_jecUp = tr->getVar<int>("nTopCandSortedCnt" + systBaseline::spec_jecUp);
    const double HT_jecUp = tr->getVar<double>("HT" + systBaseline::spec_jecUp);
    const bool passBaseline_jecUp = tr->getVar<bool>("passBaseline" + systBaseline::spec_jecUp);

    const double met_jecDn = tr->getVar<double>("met");
    const double metphi_jecDn = tr->getVar<double>("metphi");
    const double MT2_jecDn = tr->getVar<double>("best_had_brJet_MT2" + systBaseline::spec_jecDn);
    const int nbJets_jecDn = tr->getVar<int>("cntCSVS" + systBaseline::spec_jecDn);
    const int nTops_jecDn = tr->getVar<int>("nTopCandSortedCnt" + systBaseline::spec_jecDn);
    const double HT_jecDn = tr->getVar<double>("HT" + systBaseline::spec_jecDn);
    const bool passBaseline_jecDn = tr->getVar<bool>("passBaseline" + systBaseline::spec_jecDn);

    //mht calculation
    TLorentzVector Mht_LVec;	      
    for(unsigned int ij=0; ij<jetsLVec.size(); ij++){
      if( !AnaFunctions::jetPassCuts(jetsLVec[ij], AnaConsts::pt30Arr) ) continue;
      Mht_LVec -= jetsLVec[ij];
    }
    const double Mht = Mht_LVec.Pt();

// Do genHT split ONLY for TTbar
    if( sampleString.Contains("TTbar") )
    {
       if( !( (sampleString.Contains("HT") && genHT >=600) || (sampleString.Contains("Lep") && genHT < 600 ) ) ) continue;
    }
    
    bool passBaselineFull = passMuonVeto && passEleVeto && passIsoTrkVeto && passnJets && passdPhis && passMET && passBJets && passTagger && passHT && passMT2;
    bool passBaselineNoLepVeto = passnJets && passdPhis && passMET && passBJets && passTagger && passHT && passMT2;

// Lepton scale factor cannot be applied per-event on the passing veto events since it makes sense only on veto'ed events where we have well-defined leptons (@ per-event level),
    if( passBaselineNoLepVeto )
    {
      if( !passMuonVeto || !passEleVeto || !passIsoTrkVeto )
      {
  // If multiple leptons selected, since any of them could trigger the veto we use the most probable lepton which is in general lowest in pt.
  // This is because the lepton pt spectrum drops quickly for high values.
  // Assume isotrk data/MC SF is 1.0 for now -> no consideration of isolated tracks for now
        int pickedMuonIdx = -1;
        for(unsigned int im=0; im<muonsLVec.size(); im++)
        {
          if( AnaFunctions::passMuon(muonsLVec[im], muonsMiniIso[im], muonsMtw[im], muonsFlagMedium[im], AnaConsts::muonsMiniIsoArr) )
          {
            if( pickedMuonIdx == -1 ) pickedMuonIdx = (int)im;
            else if( muonsLVec[pickedMuonIdx].Pt() > muonsLVec[im].Pt() ) pickedMuonIdx = (int)im;
          }
        }
  
        int pickedEleIdx = -1;
        for(unsigned int ie=0; ie<elesLVec.size(); ie++)
        {
          if( AnaFunctions::passElectron(elesLVec[ie], elesMiniIso[ie], elesMtw[ie], elesisEB[ie], elesFlagVeto[ie], AnaConsts::elesMiniIsoArr) )
          {
            if( pickedEleIdx == -1 ) pickedEleIdx = (int)ie;
            else if( elesLVec[pickedEleIdx].Pt() > elesLVec[ie].Pt() ) pickedEleIdx = (int)ie;
          }
        }

        int pickedIsoIdx = -1;
        for(unsigned int is=0; is<isoTrksLVec.size(); is++)
        {
           if( (std::abs(isoTrkspdgId[is]) == 11 || std::abs(isoTrkspdgId[is]) == 13) && AnaFunctions::passIsoTrk(isoTrksLVec[is], isoTrksIso[is], isoTrksMtw[is], AnaConsts::isoLepTrksArr) )
           {
              if( pickedIsoIdx == -1 ) pickedIsoIdx = (int)is;
              else if( isoTrksLVec[pickedIsoIdx].Pt() > isoTrksLVec[is].Pt() ) pickedIsoIdx = (int)is;
           } 
           if( std::abs(isoTrkspdgId[is]) == 211 && AnaFunctions::passIsoTrk(isoTrksLVec[is], isoTrksIso[is], isoTrksMtw[is], AnaConsts::isoHadTrksArr) )
           {
              if( pickedIsoIdx == -1 ) pickedIsoIdx = (int)is;
              else if( isoTrksLVec[pickedIsoIdx].Pt() > isoTrksLVec[is].Pt() ) pickedIsoIdx = (int)is;
           }
        }
  
        if( pickedMuonIdx == -1 && pickedEleIdx == -1 && pickedIsoIdx == -1 )
        {
           std::cout<<"Error ... mis-matching between passMuonVeto, passEleVeto and passIsoTrkVeto from baselineDef and local evaluation??"<<std::endl;
           return 0;
        }
  
        TLorentzVector pickedLepLVec;
        int pickedLepType = -1; // 0: muon   1 : electron  2 : isotrk
  
        if( pickedMuonIdx != -1 && pickedEleIdx != -1 )
        {
          const TLorentzVector pickedMuonLVec = muonsLVec[pickedMuonIdx];
          const TLorentzVector pickedEleLVec = elesLVec[pickedEleIdx];
          if( pickedMuonLVec.Pt() < pickedEleLVec.Pt() )
          {
             pickedLepLVec = pickedMuonLVec;
             pickedLepType = 0;
          }else
          {
             pickedLepLVec = pickedEleLVec;
             pickedLepType = 1;
          } 
        }else if( pickedMuonIdx != -1 || pickedEleIdx != -1)
        {
           pickedLepLVec = pickedMuonIdx == -1? elesLVec[pickedEleIdx] : muonsLVec[pickedMuonIdx];
           pickedLepType = pickedMuonIdx == -1? 1 : 0;
        }else{
           pickedLepLVec = isoTrksLVec[pickedIsoIdx];
           pickedLepType = 2;
        }
  
        const double eta = pickedLepLVec.Eta(), pt = pickedLepLVec.Pt();
        const double abseta = std::abs(eta);
  
        double lep_id_SF = 1.0, lep_iso_SF = 1.0, lep_trk_SF = 1.0;
        double lep_id_SF_err = 0.0, lep_iso_SF_err = 0.0, lep_trk_SF_err = 0.0;
        if( pickedLepType == 0 ) // muon
        {
          if( mu_mediumID_SF ){ lep_id_SF = mu_mediumID_SF->GetBinContent(mu_mediumID_SF->FindBin(pt, abseta)); if( lep_id_SF == 0 ) lep_id_SF = 1.0; } // very simple way dealing with out of range issue of the TH2D
          if( mu_miniISO_SF ){ lep_iso_SF = mu_miniISO_SF->GetBinContent(mu_miniISO_SF->FindBin(pt, abseta)); if( lep_iso_SF == 0 ) lep_iso_SF = 1.0; }
          if( pt < 10 && mu_trkptLT10_SF ){ lep_trk_SF = mu_trkptLT10_SF->GetBinContent(mu_trkptLT10_SF->FindBin(eta)); if( lep_trk_SF == 0 ) lep_trk_SF = 1.0; }
          if( pt >= 10 && mu_trkptGT10_SF ){ lep_trk_SF = mu_trkptGT10_SF->GetBinContent(mu_trkptGT10_SF->FindBin(eta)); if( lep_trk_SF == 0 ) lep_trk_SF = 1.0; }
          lep_id_SF_err = lep_id_SF * 0.03; // only use id_SF for error to assume a flat 3% for muon case
        }else if( pickedLepType ==1 )
        {
          if( ele_VetoID_SF )
          {
            lep_id_SF = ele_VetoID_SF->GetBinContent(ele_VetoID_SF->FindBin(pt, abseta));
            lep_id_SF_err = ele_VetoID_SF->GetBinError(ele_VetoID_SF->FindBin(pt, abseta));
            if( lep_id_SF == 0 ){ lep_id_SF = 1.0; lep_id_SF_err = 0.0; }
          }
          if( ele_miniISO_SF )
          {
            lep_iso_SF = ele_miniISO_SF->GetBinContent(ele_miniISO_SF->FindBin(pt, abseta));
            lep_iso_SF_err = ele_miniISO_SF->GetBinError(ele_miniISO_SF->FindBin(pt, abseta));
            if( lep_iso_SF == 0 ){ lep_iso_SF = 1.0; lep_iso_SF_err = 0.0; }
          }
          if( ele_trkpt_SF )
          {
            lep_trk_SF = ele_trkpt_SF->GetBinContent(ele_trkpt_SF->FindBin(eta, pt));
            lep_trk_SF_err = pt<20? 0.03: 0.00;
            if( lep_trk_SF == 0 ){ lep_trk_SF = 1.0; lep_trk_SF_err = 0.0; }
          }
        }else if( pickedLepType ==2 )
        {
           lep_id_SF = 1.0; lep_iso_SF = 1.0; lep_trk_SF = 1.0;
           lep_id_SF_err = lep_id_SF * 0.10; // assume 10% for isotrk
        }else
        {
           std::cout<<"NOT supported lepton type!"<<std::endl;
           return 0;
        }
        const double lep_SF = lep_id_SF * lep_iso_SF * lep_trk_SF;
        const double rel_lep_SF_err = sqrt(lep_id_SF_err*lep_id_SF_err/lep_id_SF/lep_id_SF + lep_iso_SF_err*lep_iso_SF_err/lep_iso_SF/lep_iso_SF + lep_trk_SF_err*lep_trk_SF_err/lep_trk_SF/lep_trk_SF);
        const double lep_SF_up = lep_SF * (1 + rel_lep_SF_err);
        const double lep_SF_dn = lep_SF * (1 - rel_lep_SF_err);

        const int kSR = SB.find_Binning_Index(nbJets, nTops, MT2, met, ht);
        if(W_emuVec.size() !=0 || W_tau_emuVec.size() !=0)
        {
          if( kSR!= -1 )
          {
            myBaseHistgram.hYields_Veto_LL->Fill(kSR, Lumiscale);
            myBaseHistgram.hYields_Veto_LL_SF->Fill(kSR, Lumiscale*lep_SF);
            myBaseHistgram.hYields_Veto_LL_SFup->Fill(kSR, Lumiscale*lep_SF_up);
            myBaseHistgram.hYields_Veto_LL_SFdn->Fill(kSR, Lumiscale*lep_SF_dn);
          }
        }
        else if(W_tau_prongsVec.size() !=0)
        {
          if( kSR!= -1 )
          {
            myBaseHistgram.hYields_Veto_tau->Fill(kSR, Lumiscale);
            myBaseHistgram.hYields_Veto_tau_SF->Fill(kSR, Lumiscale*lep_SF);
            myBaseHistgram.hYields_Veto_tau_SFup->Fill(kSR, Lumiscale*lep_SF_up);
            myBaseHistgram.hYields_Veto_tau_SFdn->Fill(kSR, Lumiscale*lep_SF_dn);
          }
        }
      }else
      {
        const int kSR = SB.find_Binning_Index(nbJets, nTops, MT2, met, ht);
        if(W_emuVec.size() !=0 || W_tau_emuVec.size() !=0)
        {
          if( kSR!= -1 )
          {
            myBaseHistgram.hYields_Pass_LL->Fill(kSR, Lumiscale);
          }
        }
        else if(W_tau_prongsVec.size() !=0)
        {
          if( kSR!= -1 )
          {
            myBaseHistgram.hYields_Pass_tau->Fill(kSR, Lumiscale);
          }
        }
      }
    }
    
//    if(!(passBaselineFull && passNoiseEventFilter)) continue;
    if(!passNoiseEventFilter) continue;

    const double corr_SF = bSF * isrWght;
    const double corr_bSF_up = bSF_up * isrWght;
    const double corr_bSF_down = bSF_down * isrWght;

    const double corr_isr_up = isr_up * bSF;
    const double corr_isr_down = isr_down * bSF;

    //Exp LostLepton Dist.
    if(W_emuVec.size() !=0 || W_tau_emuVec.size() !=0)
    {
      if( passBaselineFull )
      {
        const int kSR = SB.find_Binning_Index(nbJets, nTops, MT2, met, ht);
  
        if( kSR!= -1 )
        {
  	  myBaseHistgram.hYields_LL->Fill(kSR, Lumiscale*corr_SF);
  	  //bSF systematics                                                                                                                    
  	  myBaseHistgram.hYields_LL_bSFup->Fill(kSR, Lumiscale*corr_bSF_up);
  	  myBaseHistgram.hYields_LL_bSFdown->Fill(kSR, Lumiscale*corr_bSF_down);
  	  //ISR
  	  myBaseHistgram.hYields_LL_isrup->Fill(kSR, Lumiscale*corr_isr_up);
          myBaseHistgram.hYields_LL_isrdown->Fill(kSR, Lumiscale*corr_isr_down);
          //Scale
          myBaseHistgram.hYields_LL_scaleUncup->Fill(kSR, Lumiscale*corr_SF*Scaled_Variations_Up);
          myBaseHistgram.hYields_LL_scaleUncdn->Fill(kSR, Lumiscale*corr_SF*Scaled_Variations_Down);
          //PDF
          myBaseHistgram.hYields_LL_pdfUncup->Fill(kSR, Lumiscale*corr_SF*NNPDF_From_Median_Up);
          myBaseHistgram.hYields_LL_pdfUnccen->Fill(kSR, Lumiscale*corr_SF*NNPDF_From_Median_Central);
          myBaseHistgram.hYields_LL_pdfUncdn->Fill(kSR, Lumiscale*corr_SF*NNPDF_From_Median_Down);
        }
        
        FillDouble(myBaseHistgram.hMET_LL, met, Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hMT2_LL, MT2, Lumiscale*corr_SF);
        FillInt(myBaseHistgram.hNbJets_LL, nbJets, Lumiscale*corr_SF);
        FillInt(myBaseHistgram.hNTops_LL, nTops, Lumiscale*corr_SF);	
        FillInt(myBaseHistgram.hNJets_LL, nJets, Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hHT_LL, ht, Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hdPhi0_LL, dPhiVec[0], Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hdPhi1_LL, dPhiVec[1], Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hdPhi2_LL, dPhiVec[2], Lumiscale*corr_SF);
      }

      if( passBaseline_metMagUp )
      {
        const int kSR = SB.find_Binning_Index(nbJets_metMagUp, nTops_metMagUp, MT2_metMagUp, met_metMagUp, HT_metMagUp);
        if( kSR != -1 ) myBaseHistgram.hYields_LL_metMagUp->Fill(kSR, Lumiscale*corr_SF);
      }
      if( passBaseline_metMagDn )
      {
        const int kSR = SB.find_Binning_Index(nbJets_metMagDn, nTops_metMagDn, MT2_metMagDn, met_metMagDn, HT_metMagDn);
        if( kSR != -1 ) myBaseHistgram.hYields_LL_metMagDn->Fill(kSR, Lumiscale*corr_SF);
      }
      if( passBaseline_metPhiUp )
      {
        const int kSR = SB.find_Binning_Index(nbJets_metPhiUp, nTops_metPhiUp, MT2_metPhiUp, met_metPhiUp, HT_metPhiUp);
        if( kSR != -1 ) myBaseHistgram.hYields_LL_metPhiUp->Fill(kSR, Lumiscale*corr_SF);
      }
      if( passBaseline_metPhiDn )
      {
        const int kSR = SB.find_Binning_Index(nbJets_metPhiDn, nTops_metPhiDn, MT2_metPhiDn, met_metPhiDn, HT_metPhiDn);
        if( kSR != -1 ) myBaseHistgram.hYields_LL_metPhiDn->Fill(kSR, Lumiscale*corr_SF);
      }
      if( passBaseline_jecUp )
      {
        const int kSR = SB.find_Binning_Index(nbJets_jecUp, nTops_jecUp, MT2_jecUp, met_jecUp, HT_jecUp);
        if( kSR != -1 ) myBaseHistgram.hYields_LL_jecUp->Fill(kSR, Lumiscale*corr_SF);
      }
      if( passBaseline_jecDn )
      {
        const int kSR = SB.find_Binning_Index(nbJets_jecDn, nTops_jecDn, MT2_jecDn, met_jecDn, HT_jecDn);
        if( kSR != -1 ) myBaseHistgram.hYields_LL_jecDn->Fill(kSR, Lumiscale*corr_SF);
      }
    }
    //Exp Hadtau Dist.
    else if(W_tau_prongsVec.size() !=0)
    {
      if( passBaselineFull )
      {
        const int kSR = SB.find_Binning_Index(nbJets, nTops, MT2, met, ht);

        if( kSR!= -1 )
        {
  	  myBaseHistgram.hYields_tau->Fill(kSR, Lumiscale*corr_SF);
  	  //bSF systematics                                                                                                                    
  	  myBaseHistgram.hYields_tau_bSFup->Fill(kSR, Lumiscale*corr_bSF_up);
  	  myBaseHistgram.hYields_tau_bSFdown->Fill(kSR, Lumiscale*corr_bSF_down);
  	  //ISR
  	  myBaseHistgram.hYields_tau_isrup->Fill(kSR, Lumiscale*corr_isr_up);
          myBaseHistgram.hYields_tau_isrdown->Fill(kSR, Lumiscale*corr_isr_down);
          //Scale
          myBaseHistgram.hYields_tau_scaleUncup->Fill(kSR, Lumiscale*corr_SF*Scaled_Variations_Up);
          myBaseHistgram.hYields_tau_scaleUncdn->Fill(kSR, Lumiscale*corr_SF*Scaled_Variations_Down);
          //PDF
          myBaseHistgram.hYields_tau_pdfUncup->Fill(kSR, Lumiscale*corr_SF*NNPDF_From_Median_Up);
          myBaseHistgram.hYields_tau_pdfUnccen->Fill(kSR, Lumiscale*corr_SF*NNPDF_From_Median_Central);
          myBaseHistgram.hYields_tau_pdfUncdn->Fill(kSR, Lumiscale*corr_SF*NNPDF_From_Median_Down);
        }
        
        FillDouble(myBaseHistgram.hMET_tau, met, Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hMT2_tau, MT2, Lumiscale*corr_SF);
        FillInt(myBaseHistgram.hNbJets_tau, nbJets, Lumiscale*corr_SF);
        FillInt(myBaseHistgram.hNTops_tau, nTops, Lumiscale*corr_SF);	
        FillInt(myBaseHistgram.hNJets_tau, nJets, Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hHT_tau, ht, Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hdPhi0_tau, dPhiVec[0], Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hdPhi1_tau, dPhiVec[1], Lumiscale*corr_SF);
        FillDouble(myBaseHistgram.hdPhi2_tau, dPhiVec[2], Lumiscale*corr_SF);
      }
      if( passBaseline_metMagUp )
      {
        const int kSR = SB.find_Binning_Index(nbJets_metMagUp, nTops_metMagUp, MT2_metMagUp, met_metMagUp, HT_metMagUp);
        if( kSR != -1 ) myBaseHistgram.hYields_tau_metMagUp->Fill(kSR, Lumiscale*corr_SF);
      }
      if( passBaseline_metMagDn )
      {
        const int kSR = SB.find_Binning_Index(nbJets_metMagDn, nTops_metMagDn, MT2_metMagDn, met_metMagDn, HT_metMagDn);
        if( kSR != -1 ) myBaseHistgram.hYields_tau_metMagDn->Fill(kSR, Lumiscale*corr_SF);
      }
      if( passBaseline_metPhiUp )
      {
        const int kSR = SB.find_Binning_Index(nbJets_metPhiUp, nTops_metPhiUp, MT2_metPhiUp, met_metPhiUp, HT_metPhiUp);
        if( kSR != -1 ) myBaseHistgram.hYields_tau_metPhiUp->Fill(kSR, Lumiscale*corr_SF);
      }
      if( passBaseline_metPhiDn )
      {
        const int kSR = SB.find_Binning_Index(nbJets_metPhiDn, nTops_metPhiDn, MT2_metPhiDn, met_metPhiDn, HT_metPhiDn);
        if( kSR != -1 ) myBaseHistgram.hYields_tau_metPhiDn->Fill(kSR, Lumiscale*corr_SF);
      }
      if( passBaseline_jecUp )
      {
        const int kSR = SB.find_Binning_Index(nbJets_jecUp, nTops_jecUp, MT2_jecUp, met_jecUp, HT_jecUp);
        if( kSR != -1 ) myBaseHistgram.hYields_tau_jecUp->Fill(kSR, Lumiscale*corr_SF);
      }
      if( passBaseline_jecDn )
      {
        const int kSR = SB.find_Binning_Index(nbJets_jecDn, nTops_jecDn, MT2_jecDn, met_jecDn, HT_jecDn);
        if( kSR != -1 ) myBaseHistgram.hYields_tau_jecDn->Fill(kSR, Lumiscale*corr_SF);
      }

    }
  }	//event loop
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();
  
  return 0;
  
}
