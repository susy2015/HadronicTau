#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "TauResponse.h"
#include "SusyAnaTools/Tools/samples.h"

static const int nSB = 45; //We use nSB serach bins depending on Nbjet, Ntop, met and MT2 value.
static const int nTB = nSB + 2;// one extra bin for baseline and another bin for MT2 value less than 200 GeV 

static const double llcont = 0.024;
static const double trgeff = 100/95.1;

using namespace std;

static BaselineVessel *ExpBaselineVessel;
void passBaselineFuncExp(NTupleReader& tr)
{
  (*ExpBaselineVessel)(tr);
}
double Lumiscale = 1.0;
double EventWeight = 1.0;
class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;

  TH1D *hPredHt;
  TH1D *hPredmet;
  TH1D *hPredmet_wt;
  TH1D *hPredmht;
  TH1D *hPredmht_wt;
  TH1D *hPredNJets;
  TH1D *hPredNbJets;
  TH1D *hPredNbJets_wt;
  TH1D *hPredNTops;
  TH1D *hPredNTops_wt;
  TH1D *hPredMT2;
  TH1D *hPredMT2_wt;
  TH1D *hPredmTcomb;
  TH1D *hPreddPhi0;
  TH1D *hPreddPhi0_wt;
  TH1D *hPreddPhi1;
  TH1D *hPreddPhi1_wt;
  TH1D *hPreddPhi2;
  TH1D *hPreddPhi2_wt;
  
  TH1D *hPredYields;
  TH1D *hPredYields_wt;

  TH1D *hTrueHt;
  TH1D *hTruemet;
  TH1D *hTruemht;
  TH1D *hTrueNJets;
  TH1D *hTrueNbJets;
  TH1D *hTrueNTops;
  TH1D *hTrueMT2;
  TH1D *hTruemTcomb;
  TH1D *hTrueYields;
  TH1D *hTruedPhi0;
  TH1D *hTruedPhi1;
  TH1D *hTruedPhi2;

  TH1D * hTruecutFlow;
  TH1D * hTrueFilter;
  TH1D * hPredcutFlow;
  TH1D *hCheckmet_dPhi;
  TH1D *hCheckmet_met;
  TH1D *hCheckmet_bjet;
  TH1D *hCheckmet_ntop;
  TH1D *hCheckjetpt_met;
  TH1D *hCheckjetpt_bjet;
  TH1D *hCheckjetpt_ntop;
  TH1D *hCheckjetmass_met;
  TH1D *hCheckjetmass_bjet;
  TH1D *hCheckjetmass_ntop;
  TH1D *hCheckBjetpt_met;
  TH1D *hCheckBjetpt_bjet;
  TH1D *hCheckBjetpt_ntop;
  TH1D *hChecknjet_met;
  TH1D *hChecknjet_bjet;
  TH1D *hChecknjet_ntop;

  TH1D *htotweight;
  TH1D *hcorrection;
  TH1D *htempweight;
  TH1D *hSB;
  TH2D *hweight_SB;
  TH2D *hweight_SBCpy;
  TH2D *htotweight_SB;

  TH1D *h1DSB;
  TH1D *h1DSB_wt;
  TH1D *h1DMET;
  TH1D *h1DMET_wt;
  TH1D *h1DMT2;
  TH1D *h1DMT2_wt;
  TH1D *h1DNb;
  TH1D *h1DNb_wt;
  TH1D *h1DNt;
  TH1D *h1DNt_wt;
  TH2D *h2DSB_corr;
  TH2D *h2DMET_corr;
  TH2D *h2DMT2_corr;
  TH2D *h2DNb_corr;
  TH2D *h2DNt_corr;
  
  TH1D *hmtW;
  TH1D *hnomtW;
  TH1D *hmtW_wt;
  TH1D *hnomtW_wt;
  TH1D *hmtW_Njet_wt;
  TH1D *hmtW_Ht_wt;
  TH1D *hmtW_met_wt;
  TH1D *hmtW_MT2_wt;
  TH1D *hmtW_Nbjet_wt;
  TH1D *hmtW_Ntop_wt;
  TH1D *hmtW_Njet;
  TH1D *hmtW_Ht;
  TH1D *hmtW_met;
  TH1D *hmtW_MT2;
  TH1D *hmtW_Nbjet;
  TH1D *hmtW_Ntop;  
  TH1D *hnomtW_Njet_wt;
  TH1D *hnomtW_Ht_wt;
  TH1D *hnomtW_met_wt;
  TH1D *hnomtW_MT2_wt;
  TH1D *hnomtW_Nbjet_wt;
  TH1D *hnomtW_Ntop_wt;
  TH1D *hnomtW_Njet;
  TH1D *hnomtW_Ht;
  TH1D *hnomtW_met;
  TH1D *hnomtW_MT2;
  TH1D *hnomtW_Nbjet;
  TH1D *hnomtW_Ntop;
  TH1D *htaumu;
  TH1D *htaumu_wt;
  TH1D *hnotaumu;
  TH1D *hnotaumu_wt;
  TH1D *htaumu_Njet_wt;
  TH1D *htaumu_Ht_wt;
  TH1D *htaumu_met_wt;
  TH1D *htaumu_MT2_wt;
  TH1D *htaumu_Nbjet_wt;
  TH1D *htaumu_Ntop_wt;
  TH1D *hnotaumu_Njet_wt;
  TH1D *hnotaumu_Ht_wt;
  TH1D *hnotaumu_met_wt;
  TH1D *hnotaumu_MT2_wt;
  TH1D *hnotaumu_Nbjet_wt;
  TH1D *hnotaumu_Ntop_wt;
  TH1D *htaumu_Njet;
  TH1D *htaumu_Ht;
  TH1D *htaumu_met;
  TH1D *htaumu_MT2;
  TH1D *htaumu_Nbjet;
  TH1D *htaumu_Ntop;
  TH1D *hnotaumu_Njet;
  TH1D *hnotaumu_Ht;
  TH1D *hnotaumu_met;
  TH1D *hnotaumu_MT2;
  TH1D *hnotaumu_Nbjet;
  TH1D *hnotaumu_Ntop;

  TH2D *hmtW_Njetmet_wt;
  TH2D *hmtW_Njetmet;
  TH2D *hnomtW_Njetmet_wt;
  TH2D *hnomtW_Njetmet;
  TH2D *htaumu_Njetmet_wt;
  TH2D *hnotaumu_Njetmet_wt;
  TH2D *htaumu_Njetmet;
  TH2D *hnotaumu_Njetmet;
  std::vector<TH2*> htaumu_Njetmet_Bjet;
  std::vector<TH2*> hnotaumu_Njetmet_Bjet;
  std::vector<TH2*> htaumu_Njetmet_Ht;
  std::vector<TH2*> hnotaumu_Njetmet_Ht;
  TH2D *hnomtWSys_Njetmet;
  TH2D *hmtWSysjecUp_Njetmet;
  TH2D *hmtWSysjecLow_Njetmet;
  TH2D *hmtWSysjerUp_Njetmet;
  TH2D *hmtWSysjerLow_Njetmet;
  TH1D *hdihadtau;
  TH1D *hnodihadtau;
  TH1D *hevtWt;

  const double jetbins[7] = {4, 5, 6, 7, 8, 9, 10};
  const int njetbin = sizeof(jetbins)/sizeof(jetbins[0])-1;
  const double metbins[5] = {200, 275, 350, 450, 800};
  const int nmetbin = sizeof(metbins)/sizeof(metbins[0])-1;
  const int nbjetbin = 3;
  const int nhtbin = 4;

  const double METbin[5] ={200, 275, 350, 450, 800};
  const int nMETbin = sizeof(METbin)/sizeof(METbin[0])-1;
  const double MT2bin[4] ={200, 300, 400, 800};
  const int nMT2bin = sizeof(MT2bin)/sizeof(MT2bin[0])-1;
  const double Nbbin[3] ={1, 2, 3};
  const int nNbbin = sizeof(Nbbin)/sizeof(Nbbin[0])-1;
  const double Ntbin[3] ={1, 2, 3};
  const int nNtbin = sizeof(Ntbin)/sizeof(Ntbin[0])-1;

  const TString title = "Hadronic-Tau Closure Test";
  TString Title1(unsigned int i);
  TString Title2(unsigned int j);
  TString Title3(unsigned int i);
  TString Title4(unsigned int j);
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_Closure"+index+".root";
  oFile = new TFile(filename, "recreate");
  hPredHt = new TH1D("hPredHt",title+";H_{T} [GeV];Events",20,500.,1000.);
  hPredHt->Sumw2();
  hPredmet = new TH1D("hPredmet",title+";met [GeV];Events",24,200.,800.);
  hPredmet->Sumw2();
  hPredmet_wt = new TH1D("hPredmet_wt",title+";met [GeV];Events",24,200.,800.);
  hPredmet_wt->Sumw2();
  hPredmht = new TH1D("hPredmht",title+";mht [GeV];Events",24,200.,800.);
  hPredmht->Sumw2();
  hPredmht_wt = new TH1D("hPredmht_wt",title+";mht [GeV];Events",24,200.,800.);
  hPredmht_wt->Sumw2();
  hPredNJets = new TH1D("hPredNJets",title+";N_{jets};Events",6,4,10);
  hPredNJets->Sumw2();
  hPredNbJets = new TH1D("hPredNbJets",title+";N_{bjets};Events",4, 1, 5);
  hPredNbJets->Sumw2();
  hPredNbJets_wt = new TH1D("hPredNbJets_wt",title+";N_{bjets};Events",4, 1, 5);
  hPredNbJets_wt->Sumw2();
  hPredNTops = new TH1D("hPredNTops",title+";N_{tops};Events",4, 1, 5);
  hPredNTops->Sumw2();
  hPredNTops_wt = new TH1D("hPredNTops_wt",title+";N_{tops};Events",4, 1, 5);
  hPredNTops_wt->Sumw2();
  hPredMT2 = new TH1D("hPredMT2",title+";M_{T2}[GeV];Events",12,200,500);
  hPredMT2->Sumw2();
  hPredMT2_wt = new TH1D("hPredMT2_wt",title+";M_{T2}[GeV];Events",12,200,500);
  hPredMT2_wt->Sumw2();
  hPredmTcomb = new TH1D("hPredmTcomb",title+";M_{Tb}+0.5*M_{Tt}[GeV];Events",50,0,1000);
  hPredmTcomb->Sumw2();
  hPreddPhi0 = new TH1D("hPreddPhi0",title+";dPhi0;Events",16,0,3.2);
  hPreddPhi0->Sumw2();
  hPreddPhi0_wt = new TH1D("hPreddPhi0_wt",title+";dPhi0;Events",16,0,3.2);
  hPreddPhi0_wt->Sumw2();
  hPreddPhi1 = new TH1D("hPreddPhi1",title+";dPhi1;Events",16,0,3.2);
  hPreddPhi1->Sumw2();
  hPreddPhi1_wt = new TH1D("hPreddPhi1_wt",title+";dPhi1;Events",16,0,3.2);
  hPreddPhi1_wt->Sumw2();
  hPreddPhi2 = new TH1D("hPreddPhi2",title+";dPhi2;Events",16,0,3.2);
  hPreddPhi2->Sumw2();
  hPreddPhi2_wt = new TH1D("hPreddPhi2_wt",title+";dPhi2;Events",16,0,3.2);
  hPreddPhi2_wt->Sumw2();
  
  hPredYields_wt = new TH1D("hPredYields_wt", title+";search bin;Events",nSB,0,nSB);
  hPredYields_wt->Sumw2();
  hPredYields = new TH1D("hPredYields", title+";search bin;Events",nSB,0,nSB);
  hPredYields->Sumw2();

  hTrueHt = new TH1D("hTrueHt",title+";H_{T} [GeV];Events",20,500.,1000.);
  hTrueHt->Sumw2();
  hTruemet = new TH1D("hTruemet",title+";met [GeV];Events",24,200.,800.);
  hTruemet->Sumw2();
  hTruemht = new TH1D("hTruemht",title+";mht [GeV];Events",24,200.,800.);
  hTruemht->Sumw2();
  hTrueNJets = new TH1D("hTrueNJets",title+";N_{jets};Events",6 ,4,10);
  hTrueNJets->Sumw2();
  hTrueNbJets = new TH1D("hTrueNbJets",title+";N_{bjets};Events",4, 1, 5);
  hTrueNbJets->Sumw2();
  hTrueNTops = new TH1D("hTrueNTops",title+";N_{tops};Events",4, 1, 5);
  hTrueNTops->Sumw2();
  hTrueMT2 = new TH1D("hTrueMT2",title+";M_{T2}[GeV];Events",12,200,500);
  hTrueMT2->Sumw2();
  hTruemTcomb = new TH1D("hTruemTcomb",title+";M_{Tb}+0.5*M_{Tt}[GeV];Events",50,0,1000);
  hTruemTcomb->Sumw2();
  hTrueYields = new TH1D("hTrueYields", title+";search bin;Events",nSB,0,nSB);
  hTrueYields->Sumw2();

  hTruedPhi0 = new TH1D("hTruedPhi0", title+";dPhi0;Events", 16, 0, 3.2);
  hTruedPhi0->Sumw2();
  hTruedPhi1 = new TH1D("hTruedPhi1", title+";dPhi1;Events", 16, 0, 3.2);
  hTruedPhi1->Sumw2();
  hTruedPhi2 = new TH1D("hTruedPhi2", title+";dPhi2;Events", 16, 0, 3.2);
  hTruedPhi2->Sumw2();

  hTruecutFlow = new TH1D("hTruecutFlow", "cut flow table", 20, 0, 20);
  hTruecutFlow->SetBit(TH1::kCanRebin);                                                                                            
  hTruecutFlow->Sumw2();

  hPredcutFlow = new TH1D("hPredcutFlow", "cut flow table", 20, 0, 20);
  hPredcutFlow->SetBit(TH1::kCanRebin);
  hPredcutFlow->Sumw2();

  hTrueFilter = new TH1D("hTrueFilter", "Event filter table", 3, 0, 3);
  hTrueFilter->SetBit(TH1::kCanRebin);
  hTrueFilter->Sumw2();

  hCheckmet_dPhi = new TH1D("hCheckmet_dPhi", "met", 40, 0, 800);
  hCheckmet_dPhi->Sumw2();
  hCheckmet_met = new TH1D("hCheckmet_met", "met", 40, 0, 800);
  hCheckmet_met->Sumw2();
  hCheckmet_bjet = new TH1D("hCheckmet_bjet", "met", 40, 0, 800);
  hCheckmet_bjet->Sumw2();
  hCheckmet_ntop = new TH1D("hCheckmet_ntop", "met", 40, 0, 800);
  hCheckmet_ntop->Sumw2();
  hCheckjetpt_met = new TH1D("hCheckjetpt_met", "jet p_{T}", 50 , 0, 1000);
  hCheckjetpt_met->Sumw2();
  hCheckjetpt_bjet = new TH1D("hCheckjetpt_bjet", "jet p_{T}", 50 , 0, 1000);
  hCheckjetpt_bjet->Sumw2();
  hCheckjetpt_ntop = new TH1D("hCheckjetpt_ntop", "jet p_{T}", 50 , 0, 1000);
  hCheckjetpt_ntop->Sumw2();
  hCheckjetmass_met = new TH1D("hCheckjetmass_met", "jet mass", 50 , 0, 1000);
  hCheckjetmass_met->Sumw2();
  hCheckjetmass_bjet = new TH1D("hCheckjetmass_bjet", "jet mass", 50 , 0, 1000);
  hCheckjetmass_bjet->Sumw2();
  hCheckjetmass_ntop = new TH1D("hCheckjetmass_ntop", " jet mass", 50 , 0, 1000);
  hCheckjetmass_ntop->Sumw2();
  hCheckBjetpt_met = new TH1D("hCheckBjetpt_met", "Bjet p_{T}", 50 , 0, 1000);
  hCheckBjetpt_met->Sumw2();
  hCheckBjetpt_bjet = new TH1D("hCheckBjetpt_bjet", "Bjet p_{T}", 50 , 0, 1000);
  hCheckBjetpt_bjet->Sumw2();
  hCheckBjetpt_ntop = new TH1D("hCheckBjetpt_ntop", "Bjet p_{T}", 50 , 0, 1000);
  hCheckBjetpt_ntop->Sumw2();
  hChecknjet_met = new TH1D("hChecknjet_met", "N_{jet}", 6 , 4, 10);
  hChecknjet_met->Sumw2();
  hChecknjet_bjet = new TH1D("hChecknjet_bjet", "N_{jet}", 6 , 4, 10);
  hChecknjet_bjet->Sumw2();
  hChecknjet_ntop = new TH1D("hChecknjet_ntop", "N_{jet}", 6 , 4, 10);
  hChecknjet_ntop->Sumw2();

  htotweight = new TH1D("htotweight", "total weight", 30, 0, 1.5);
  htotweight->Sumw2();
  hcorrection = new TH1D("hcorrection", "correction", 30, 0, 1.5);
  hcorrection->Sumw2();
  htempweight = new TH1D("htempweight", "template weight", 30, 0, 1.5);
  htempweight->Sumw2();
  hweight_SB = new TH2D("hweight_SB", "weight vs SB", nSB, 0, nSB, 30, 0, 1.5);
  hweight_SB->Sumw2();
  hweight_SBCpy = new TH2D("hweight_SBCpy", "weight vs SB", nSB, 0, nSB, 30, 0, 1.5);
  hweight_SBCpy->Sumw2();
  hSB = new TH1D("hSB", "weight vs SB", nSB, 0, nSB);
  hSB->Sumw2();
  htotweight_SB = new TH2D("htotweight_SB", "weight vs SB", nSB, 0, nSB, 30, 0, 1.5);
  htotweight_SB->Sumw2();

  h1DSB = new TH1D("h1DSB", "Searchbin", nSB, 0, nSB);
  h1DSB_wt = new TH1D("h1DSB_wt", "Searchbin", nSB, 0, nSB);
  h2DSB_corr = new TH2D("h2DSB_corr", "SearchbinCorr", nSB, 0, nSB, nSB, 0, nSB);
  h1DMET = new TH1D("h1DMET", "SearchbinMET", nMETbin, METbin);
  h1DMET_wt = new TH1D("h1DMET_wt", "SearchbinMET", nMETbin, METbin);
  h2DMET_corr = new TH2D("h2DMET_corr", "SearchbinMETCorr", nMETbin, METbin, nMETbin, METbin);
  h1DMT2 = new TH1D("h1DMT2", "SearchbinMT2", nMT2bin, MT2bin);
  h1DMT2_wt = new TH1D("h1DMT2_wt", "SearchbinMET", nMT2bin, MT2bin);
  h2DMT2_corr = new TH2D("h2DMT2_corr", "SearchbinMT2Corr", nMT2bin, MT2bin, nMT2bin, MT2bin);
  h1DNb = new TH1D("h1DNb", "SearchbinNb", nNbbin, Nbbin);
  h1DNb_wt = new TH1D("h1DNb_wt", "SearchbinNb", nNbbin, Nbbin);
  h2DNb_corr = new TH2D("h2DNb_corr", "SearchbinNbCorr", nNbbin, Nbbin, nNbbin, Nbbin);
  h1DNt = new TH1D("h1DNt", "SearchbinNt", nNtbin, Ntbin);
  h1DNt_wt = new TH1D("h1DNt_wt", "SearchbinNt", nNtbin, Ntbin);
  h2DNt_corr = new TH2D("h2DNt_corr", "SearchbinNtCorr", nNtbin, Ntbin, nNtbin, Ntbin);

  hmtW = new TH1D("hmtW", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hmtW->Sumw2();
  hnomtW = new TH1D("hnomtW", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hnomtW->Sumw2();
  hmtW_wt = new TH1D("hmtW_wt", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hmtW_wt->Sumw2();
  hnomtW_wt = new TH1D("hnomtW_wt", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hnomtW_wt->Sumw2();
  
  hnomtW_Njet = new TH1D("hnomtW_Njet", "mtW correction;N_{jet};Events", 6, 4, 10);
  hnomtW_Njet->Sumw2();
  hnomtW_Njet_wt = new TH1D("hnomtW_Njet_wt", "mtW correction;N_{jet};Events", 6, 4, 10);
  hnomtW_Njet_wt->Sumw2();
  hmtW_Njet = new TH1D("hmtW_Njet", "mtW correction;N_{jet};Events", 6, 4, 10);
  hmtW_Njet->Sumw2();
  hmtW_Njet_wt = new TH1D("hmtW_Njet_wt", "mtW correction;N_{jet};Events", 6, 4, 10);
  hmtW_Njet_wt->Sumw2();

  hnomtW_Nbjet = new TH1D("hnomtW_Nbjet", "mtW correction;N_{bjet};Events", 4, 1, 5);
  hnomtW_Nbjet->Sumw2();
  hnomtW_Nbjet_wt = new TH1D("hnomtW_Nbjet_wt", "mtW correction;N_{bjet};Events", 4, 1, 5);
  hnomtW_Nbjet_wt->Sumw2();
  hmtW_Nbjet = new TH1D("hmtW_Nbjet", "mtW correction;N_{bjet};Events", 4, 1, 5);
  hmtW_Nbjet->Sumw2();
  hmtW_Nbjet_wt = new TH1D("hmtW_Nbjet_wt", "mtW correction;N_{bjet};Events", 4, 1, 5);
  hmtW_Nbjet_wt->Sumw2();

  hnomtW_Ntop = new TH1D("hnomtW_Ntop", "mtW correction;N_{top};Events", 4, 1, 5);
  hnomtW_Ntop->Sumw2();
  hnomtW_Ntop_wt = new TH1D("hnomtW_Ntop_wt", "mtW correction;N_{top};Events", 4, 1, 5);
  hnomtW_Ntop_wt->Sumw2();
  hmtW_Ntop = new TH1D("hmtW_Ntop", "mtW correction;N_{top};Events", 4, 1, 5);
  hmtW_Ntop->Sumw2();
  hmtW_Ntop_wt = new TH1D("hmtW_Ntop_wt", "mtW correction;N_{top};Events", 4, 1, 5);  
  hmtW_Ntop_wt->Sumw2();

  hmtW_Ht = new TH1D("hmtW_Ht", "mtW correction;HT;Events", 10, 500, 1000);
  hmtW_Ht->Sumw2();
  hmtW_Ht_wt = new TH1D("hmtW_Ht_wt", "mtW correction;HT;Events", 10, 500, 1000);
  hmtW_Ht_wt->Sumw2();
  hnomtW_Ht = new TH1D("hnomtW_Ht", "mtW correction;HT;Events", 10, 500, 1000);
  hnomtW_Ht->Sumw2();
  hnomtW_Ht_wt = new TH1D("hnomtW_Ht_wt", "mtW correction;HT;Events", 10, 500, 1000);
  hnomtW_Ht_wt->Sumw2();

  hmtW_met = new TH1D("hmtW_met", "mtW correction;met;Events", 24, 200, 800);
  hmtW_met->Sumw2();
  hmtW_met_wt = new TH1D("hmtW_met_wt", "mtW correction;met;Events", 24, 200, 800);
  hmtW_met_wt->Sumw2();
  hnomtW_met = new TH1D("hnomtW_met", "mtW correction;met;Events", 24, 200, 800);
  hnomtW_met->Sumw2();
  hnomtW_met_wt = new TH1D("hnomtW_met_wt", "mtW correction;met;Events", 24, 200, 800);
  hnomtW_met_wt->Sumw2();

  hmtW_MT2 = new TH1D("hmtW_MT2", "mtW correction;MT2;Events", 8, 100, 500);
  hmtW_MT2->Sumw2();
  hmtW_MT2_wt = new TH1D("hmtW_MT2_wt", "mtW correction;MT2;Events", 8, 100, 500);
  hmtW_MT2_wt->Sumw2();
  hnomtW_MT2 = new TH1D("hnomtW_MT2", "mtW correction;MT2;Events", 8, 100, 500);
  hnomtW_MT2->Sumw2();
  hnomtW_MT2_wt = new TH1D("hnomtW_MT2_wt", "mtW correction;MT2;Events", 8, 100, 500);
  hnomtW_MT2_wt->Sumw2();

  hmtW_Njetmet = new TH2D("hmtW_Njetmet", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hmtW_Njetmet->Sumw2();
  hmtW_Njetmet_wt = new TH2D("hmtW_Njetmet_wt", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hmtW_Njetmet_wt->Sumw2();
  hnomtW_Njetmet = new TH2D("hnomtW_Njetmet", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hnomtW_Njetmet->Sumw2();
  hnomtW_Njetmet_wt = new TH2D("hnomtW_Njetmet_wt", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hnomtW_Njetmet_wt->Sumw2();

  htaumu = new TH1D("htaumu", "taumu contamination;Search bin;Events", nSB, 0, nSB);
  htaumu->Sumw2();
  hnotaumu = new TH1D("hnotaumu", "taumu contamination;Search bin;Events", nSB, 0, nSB);
  hnotaumu->Sumw2();
  htaumu_wt = new TH1D("htaumu_wt", "taumu contamination;Search bin;Events", nSB, 0, nSB);
  htaumu_wt->Sumw2();
  hnotaumu_wt = new TH1D("hnotaumu_wt", "taumu contamination;Search bin;Events", nSB, 0, nSB);
  hnotaumu_wt->Sumw2();

  hnotaumu_Njet = new TH1D("hnotaumu_Njet", "taumu contamination;N_{jet};Events", 6, 4, 10);
  hnotaumu_Njet->Sumw2();
  hnotaumu_Njet_wt = new TH1D("hnotaumu_Njet_wt", "taumu contamination;N_{jet};Events", 6, 4, 10);
  hnotaumu_Njet_wt->Sumw2();
  htaumu_Njet = new TH1D("htaumu_Njet", "taumu contamination;N_{jet};Events", 6, 4, 10);
  htaumu_Njet->Sumw2();
  htaumu_Njet_wt = new TH1D("htaumu_Njet_wt", "taumu contamination;N_{jet};Events", 6, 4, 10);
  htaumu_Njet_wt->Sumw2();

  hnotaumu_Nbjet = new TH1D("hnotaumu_Nbjet", "taumu contamination;N_{bjet};Events", 4, 1, 5);
  hnotaumu_Nbjet->Sumw2();
  hnotaumu_Nbjet_wt = new TH1D("hnotaumu_Nbjet_wt", "taumu contamination;N_{bjet};Events", 4, 1, 5);
  hnotaumu_Nbjet_wt->Sumw2();
  htaumu_Nbjet = new TH1D("htaumu_Nbjet", "taumu contamination;N_{bjet};Events", 4, 1, 5);
  htaumu_Nbjet->Sumw2();
  htaumu_Nbjet_wt = new TH1D("htaumu_Nbjet_wt", "taumu contamination;N_{bjet};Events", 4, 1, 5);
  htaumu_Nbjet_wt->Sumw2();
  
  hnotaumu_Ntop = new TH1D("hnotaumu_Ntop", "taumu contamination;N_{top};Events", 4, 1, 5);
  hnotaumu_Ntop->Sumw2();
  hnotaumu_Ntop_wt = new TH1D("hnotaumu_Ntop_wt", "taumu contamination;N_{top};Events", 4, 1, 5);
  hnotaumu_Ntop_wt->Sumw2();
  htaumu_Ntop = new TH1D("htaumu_Ntop", "taumu contamination;N_{top};Events", 4, 1, 5);
  htaumu_Ntop->Sumw2();
  htaumu_Ntop_wt = new TH1D("htaumu_Ntop_wt", "taumu contamination;N_{top};Events", 4, 1, 5);
  htaumu_Ntop_wt->Sumw2();

  htaumu_Ht = new TH1D("htaumu_Ht", "taumu contamination;HT;Events", 10, 500, 1000);
  htaumu_Ht->Sumw2();
  htaumu_Ht_wt = new TH1D("htaumu_Ht_wt", "taumu contamination;HT;Events", 10, 500, 1000);
  htaumu_Ht_wt->Sumw2();
  hnotaumu_Ht = new TH1D("hnotaumu_Ht", "taumu contamination;HT;Events", 10, 500, 1000);
  hnotaumu_Ht->Sumw2();
  hnotaumu_Ht_wt = new TH1D("hnotaumu_Ht_wt", "taumu contamination;HT;Events", 10, 500, 1000);
  hnotaumu_Ht_wt->Sumw2();

  htaumu_met = new TH1D("htaumu_met", "taumu contamination;met;Events", 24, 200, 800);
  htaumu_met->Sumw2();
  htaumu_met_wt = new TH1D("htaumu_met_wt", "taumu contamination;met;Events", 24, 200, 800);
  htaumu_met_wt->Sumw2();
  hnotaumu_met = new TH1D("hnotaumu_met", "taumu contamination;met;Events", 24, 200, 800);
  hnotaumu_met->Sumw2();
  hnotaumu_met_wt = new TH1D("hnotaumu_met_wt", "taumu contamination;met;Events", 24, 200, 800);
  hnotaumu_met_wt->Sumw2();

  htaumu_MT2 = new TH1D("htaumu_MT2", "taumu contamination;MT2;Events", 8, 100, 500);
  htaumu_MT2->Sumw2();
  htaumu_MT2_wt = new TH1D("htaumu_MT2_wt", "taumu contamination;MT2;Events", 8, 100, 500);
  htaumu_MT2_wt->Sumw2();
  hnotaumu_MT2 = new TH1D("hnotaumu_MT2", "taumu contamination;MT2;Events", 8, 100, 500);
  hnotaumu_MT2->Sumw2();
  hnotaumu_MT2_wt = new TH1D("hnotaumu_MT2_wt", "taumu contamination;MT2;Events", 8, 100, 500);
  hnotaumu_MT2_wt->Sumw2();

  htaumu_Njetmet = new TH2D("htaumu_Njetmet", "taumu correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  htaumu_Njetmet->Sumw2();
  htaumu_Njetmet_wt = new TH2D("htaumu_Njetmet_wt", "taumu correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  htaumu_Njetmet_wt->Sumw2();
  hnotaumu_Njetmet = new TH2D("hnotaumu_Njetmet", "taumu correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hnotaumu_Njetmet->Sumw2();
  hnotaumu_Njetmet_wt = new TH2D("hnotaumu_Njetmet_wt", "taumu correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hnotaumu_Njetmet_wt->Sumw2();

  for(unsigned int i = 0; i < nbjetbin; ++i) {
    htaumu_Njetmet_Bjet.push_back(new TH2D(Title1(i),"taumu correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins));
    htaumu_Njetmet_Bjet.back()->Sumw2();
    hnotaumu_Njetmet_Bjet.push_back(new TH2D(Title2(i),"taumu correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins));
    hnotaumu_Njetmet_Bjet.back()->Sumw2();
  }
  for(unsigned int i = 0; i < nhtbin; ++i) {
    htaumu_Njetmet_Ht.push_back(new TH2D(Title3(i),"taumu correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins));
    htaumu_Njetmet_Ht.back()->Sumw2();
    hnotaumu_Njetmet_Ht.push_back(new TH2D(Title4(i),"taumu correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins));
    hnotaumu_Njetmet_Ht.back()->Sumw2();
  }   

  hnomtWSys_Njetmet = new TH2D("hnomtWSys_Njetmet", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hnomtWSys_Njetmet->Sumw2();
  hmtWSysjecUp_Njetmet = new TH2D("hmtWSysjecUp_Njetmet", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hmtWSysjecUp_Njetmet->Sumw2();
  hmtWSysjecLow_Njetmet = new TH2D("hmtWSysjecLow_Njetmet", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hmtWSysjecLow_Njetmet->Sumw2();
  hmtWSysjerUp_Njetmet = new TH2D("hmtWSysjerUp_Njetmet", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hmtWSysjerUp_Njetmet->Sumw2();
  hmtWSysjerLow_Njetmet = new TH2D("hmtWSysjerLow_Njetmet", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hmtWSysjerLow_Njetmet->Sumw2();
  hdihadtau = new TH1D("hdihadtau", "Di Hadronic tau fraction;Search bin;Events", nSB, 0, nSB);
  hdihadtau->Sumw2();
  hnodihadtau = new TH1D("hnodihadtau", "Di Hadronic tau fraction;Search bin;Events", nSB, 0, nSB);
  hnodihadtau->Sumw2();
  hevtWt = new TH1D("hevtWt", "Event weight form Theo. Xsec", 4,-2, 2);
  hevtWt->Sumw2();
}

TString BaseHistgram::Title1(unsigned int i){
  TString title = "htaumu_Njetmet_Bjet";
  title+=i;
  return title;
}
TString BaseHistgram::Title2(unsigned int j){
  TString title ="hnotaumu_Njetmet_Bjet";
  title+=j;
  return title;
}
TString BaseHistgram::Title3(unsigned int i){
  TString title ="htaumu_Njetmet_Ht";
  title+=i;
  return title;
}
TString BaseHistgram::Title4(unsigned int j){
  TString title ="hnotaumu_Njetmet_Ht";
  title+=j;
  return title;
}

bool FillChain(TChain* &chain, const char *subsample, const string condorSpec, const int& startfile, const int& filerun){
  
  AnaSamples::SampleSet        allSamples = condorSpec.empty()? AnaSamples::SampleSet():AnaSamples::SampleSet(condorSpec);
  AnaSamples::SampleCollection allCollections(allSamples);
  bool find = false;  
  TString subsamplename(subsample);
  
  chain = new TChain(allSamples[subsample].treePath.c_str());
  if(allSamples[subsample] != allSamples.null())
    {
      allSamples[subsample].addFilesToChain(chain, startfile, filerun);
      find = true;
      Lumiscale = allSamples[subsample].getWeight();
    }
    return find;
}

double htJetPtMin(){ return 50;}
double htJetEtaMax() {return 2.4;}
double mhtJetPtMin(){return 30;}
double mhtJetEtaMax() {return 5;}
double nJetPtMin(){return 30;}
double nJetEtaMax() {return 2.4;}

void drawOverFlowBin(TH1 *hist){

  int nbins = hist->GetXaxis()->GetNbins();

  double overflow = hist->GetBinContent(nbins+1);
  double lastCont = hist->GetBinContent(nbins);
  double ovrflweroor = hist->GetBinError(nbins+1);
  double lstbinerror = hist->GetBinError(nbins);
  hist->SetBinContent(nbins, overflow+lastCont);
  hist->SetBinError(nbins, TMath::Sqrt(ovrflweroor * ovrflweroor + lstbinerror * lstbinerror));
}
void FillDouble(TH1 *hist, const double &a, const double &w){
  int nbin = hist->GetNbinsX();
  double low = hist->GetBinLowEdge(nbin);
  double high = hist->GetBinLowEdge(nbin + 1);
  double copy = a;
  if(copy >= high) copy = low;
  hist->Fill(copy, w);
}
void FillInt(TH1 *hist, const int &a, const double &w){
  int nbin = hist->GetNbinsX();
  int low = (int)hist->GetBinLowEdge(nbin);
  int high = (int)hist->GetBinLowEdge(nbin + 1);
  int copy = a;
  if(copy >= high) copy = low;
  hist->Fill(copy, w);
}
void Fill2D(TH2 *hist, const int &a, const double &b, const double &w){
  int nbinx = hist->GetNbinsX();
  int nbiny = hist->GetNbinsY();
  int lowx = hist->GetXaxis()->GetBinLowEdge(nbinx);
  int highx = hist->GetXaxis()->GetBinLowEdge(nbinx + 1);
  double lowy = hist->GetYaxis()->GetBinLowEdge(nbiny);
  double highy = hist->GetYaxis()->GetBinLowEdge(nbiny + 1);
  int copyx = a;
  if(copyx >= highx) copyx = lowx;
  double copyy = b;
  if(copyy >= highy) copyy = lowy;
  hist->Fill(copyx, copyy, w);
}

double calcMT(const TLorentzVector &objLVec, const TLorentzVector &metLVec){

  const double objMass = objLVec.M(), objPt = objLVec.Pt(), objPx = objLVec.Px(), objPy = objLVec.Py();

  const double met = metLVec.Pt(), metphi = metLVec.Phi();

  double mt = sqrt( objMass*objMass + 2*( met*sqrt(objMass*objMass + objPt*objPt) -( met*cos(metphi)*objPx + met*sin(metphi)*objPy ) ) );

  return mt;

}
bool passIsoTrks1(const TLorentzVector isoTrksLVec, const double isoTrksIso, const double isoTrksMtw, const int isoTrkspdgId);
