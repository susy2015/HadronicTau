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

static const int nSB = 59; //We use nSB serach bins depending on Nbjet, Ntop, met and MT2 value.

static const double llcont = 0.035;
static const double trgeff = 100/91.4;

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
  TH1D *hTruedPhi0;
  TH1D *hTruedPhi1;
  TH1D *hTruedPhi2;
  TH1D *hTrueYields;

  TH1D *hdihadtau;
  TH1D *hnodihadtau;

  TH1D *hmtW;
  TH1D *hnomtW;
  TH1D *hmtW_wt;
  TH1D *hnomtW_wt;
  TH1D *htaumu;
  TH1D *hnotaumu;
  TH1D *htaumu_wt;
  TH1D *hnotaumu_wt;
  //jec, jer corr.
  TH1D *hnomtWSys;
  TH1D *hmtWSysjecUp;
  TH1D *hmtWSysjecLow;
  TH1D *hmtWSysjerUp;
  TH1D *hmtWSysjerLow;

  TH2D *hmtW_Njetmet_wt;
  TH2D *hmtW_Njetmet;
  TH2D *hnomtW_Njetmet_wt;
  TH2D *hnomtW_Njetmet;
  TH2D *htaumu_Njetmet_wt;
  TH2D *hnotaumu_Njetmet_wt;
  TH2D *htaumu_Njetmet;
  TH2D *hnotaumu_Njetmet;

  TH1D *hmtW_Njet;
  TH1D *hnomtW_Njet;
  TH1D *hmtW_Nbjet;
  TH1D *hnomtW_Nbjet;
  TH1D *hmtW_Ntop;
  TH1D *hnomtW_Ntop;
  TH1D *hmtW_met;
  TH1D *hnomtW_met;
  TH1D *hmtW_MT2;
  TH1D *hnomtW_MT2;
  TH1D *htaumu_Njet;
  TH1D *hnotaumu_Njet;
  TH1D *htaumu_Nbjet;
  TH1D *hnotaumu_Nbjet;
  TH1D *htaumu_Ntop;
  TH1D *hnotaumu_Ntop;
  TH1D *htaumu_met;
  TH1D *hnotaumu_met;
  TH1D *htaumu_MT2;
  TH1D *hnotaumu_MT2;

  const double jetbins[7] = {4, 5, 6, 7, 8, 9, 10};
  const int njetbin = sizeof(jetbins)/sizeof(jetbins[0])-1;
  const double metbins[5] = {200, 350, 500, 650, 800};
  const int nmetbin = sizeof(metbins)/sizeof(metbins[0])-1;
  const int nbjetbin = 3;
  const int nhtbin = 4;

  const TString title = "Hadronic-Tau Closure Test";
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

  hmtW = new TH1D("hmtW", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hmtW->Sumw2();
  hnomtW = new TH1D("hnomtW", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hnomtW->Sumw2();
  hmtW_wt = new TH1D("hmtW_wt", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hmtW_wt->Sumw2();
  hnomtW_wt = new TH1D("hnomtW_wt", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hnomtW_wt->Sumw2();

  hnomtWSys = new TH1D("hnomtWSys", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hnomtWSys->Sumw2();
  hmtWSysjecUp = new TH1D("hmtWSysjecUp", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hmtWSysjecUp->Sumw2();
  hmtWSysjecLow = new TH1D("hmtWSysjecLow", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hmtWSysjecLow->Sumw2();
  hmtWSysjerUp = new TH1D("hmtWSysjerUp", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hmtWSysjerUp->Sumw2();
  hmtWSysjerLow = new TH1D("hmtWSysjerLow", "mtW correction;Search bin;Events", nSB, 0, nSB);
  hmtWSysjerLow->Sumw2();

  hmtW_Njetmet = new TH2D("hmtW_Njetmet", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hmtW_Njetmet->Sumw2();
  hmtW_Njetmet_wt = new TH2D("hmtW_Njetmet_wt", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hmtW_Njetmet_wt->Sumw2();
  hnomtW_Njetmet = new TH2D("hnomtW_Njetmet", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hnomtW_Njetmet->Sumw2();
  hnomtW_Njetmet_wt = new TH2D("hnomtW_Njetmet_wt", "mtW correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hnomtW_Njetmet_wt->Sumw2();
  
  hmtW_Njet = new TH1D("hmtW_Njet", "taumu correction;N_{jet};Events", 6, 4, 10);
  hmtW_Njet->Sumw2();
  hnomtW_Njet = new TH1D("hnomtW_Njet", "taumu correction;N_{jet};Events", 6, 4, 10);
  hnomtW_Njet->Sumw2();
  hmtW_Nbjet = new TH1D("hmtW_Nbjet", "taumu correction;N_{bjet};Events", 4, 1, 5);
  hmtW_Nbjet->Sumw2();
  hnomtW_Nbjet = new TH1D("hnomtW_Nbjet", "taumu correction;N_{bjet};Events", 4, 1, 5);
  hnomtW_Nbjet->Sumw2();
  hmtW_Ntop = new TH1D("hmtW_Ntop", "taumu correction;N_{top};Events", 4, 1, 5);
  hmtW_Ntop->Sumw2();
  hnomtW_Ntop = new TH1D("hnomtW_Ntop", "taumu correction;N_{top};Events", 4, 1, 5);
  hnomtW_Ntop->Sumw2();
  hmtW_met = new TH1D("hmtW_met", "taumu correction;p_{T}^{miss};Events", 32, 200, 1000);
  hmtW_met->Sumw2();
  hnomtW_met = new TH1D("hnomtW_met", "taumu correction;p_{T}^{miss};Events", 32, 200, 1000);
  hnomtW_met->Sumw2();
  hmtW_MT2 = new TH1D("hmtW_MT2", "taumu correction;M_{T2};Events", 24, 200, 800);
  hmtW_MT2->Sumw2();
  hnomtW_MT2 = new TH1D("hnomtW_MT2", "taumu correction;M_{T2};Events", 24, 200, 800);
  hnomtW_MT2->Sumw2();

  htaumu = new TH1D("htaumu", "taumu contamination;Search bin;Events", nSB, 0, nSB);
  htaumu->Sumw2();
  hnotaumu = new TH1D("hnotaumu", "taumu contamination;Search bin;Events", nSB, 0, nSB);
  hnotaumu->Sumw2();
  htaumu_wt = new TH1D("htaumu_wt", "taumu contamination;Search bin;Events", nSB, 0, nSB);
  htaumu_wt->Sumw2();
  hnotaumu_wt = new TH1D("hnotaumu_wt", "taumu contamination;Search bin;Events", nSB, 0, nSB);
  hnotaumu_wt->Sumw2();

  htaumu_Njetmet = new TH2D("htaumu_Njetmet", "taumu correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  htaumu_Njetmet->Sumw2();
  htaumu_Njetmet_wt = new TH2D("htaumu_Njetmet_wt", "taumu correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  htaumu_Njetmet_wt->Sumw2();
  hnotaumu_Njetmet = new TH2D("hnotaumu_Njetmet", "taumu correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hnotaumu_Njetmet->Sumw2();
  hnotaumu_Njetmet_wt = new TH2D("hnotaumu_Njetmet_wt", "taumu correction;N_{jet};met", njetbin, jetbins, nmetbin, metbins);
  hnotaumu_Njetmet_wt->Sumw2();

  htaumu_Njet = new TH1D("htaumu_Njet", "taumu correction;N_{jet};Events", 6, 4, 10);
  htaumu_Njet->Sumw2();
  hnotaumu_Njet = new TH1D("hnotaumu_Njet", "taumu correction;N_{jet};Events", 6, 4, 10);
  hnotaumu_Njet->Sumw2();
  htaumu_Nbjet = new TH1D("htaumu_Nbjet", "taumu correction;N_{bjet};Events", 4, 1, 5);
  htaumu_Nbjet->Sumw2();
  hnotaumu_Nbjet = new TH1D("hnotaumu_Nbjet", "taumu correction;N_{bjet};Events", 4, 1, 5);
  hnotaumu_Nbjet->Sumw2();
  htaumu_Ntop = new TH1D("htaumu_Ntop", "taumu correction;N_{top};Events", 4, 1, 5);
  htaumu_Ntop->Sumw2();
  hnotaumu_Ntop = new TH1D("hnotaumu_Ntop", "taumu correction;N_{top};Events", 4, 1, 5);
  hnotaumu_Ntop->Sumw2();
  htaumu_met = new TH1D("htaumu_met", "taumu correction;p_{T}^{miss};Events", 32, 200, 1000);
  htaumu_met->Sumw2();
  hnotaumu_met = new TH1D("hnotaumu_met", "taumu correction;p_{T}^{miss};Events", 32, 200, 1000);
  hnotaumu_met->Sumw2();
  htaumu_MT2 = new TH1D("htaumu_MT2", "taumu correction;M_{T2};Events", 24, 200, 800);
  htaumu_MT2->Sumw2();
  hnotaumu_MT2 = new TH1D("hnotaumu_MT2", "taumu correction;M_{T2};Events", 24, 200, 800);
  hnotaumu_MT2->Sumw2();

  hdihadtau = new TH1D("hdihadtau", "Di Hadronic tau fraction;Search bin;Events", nSB, 0, nSB);
  hdihadtau->Sumw2();
  hnodihadtau = new TH1D("hnodihadtau", "Di Hadronic tau fraction;Search bin;Events", nSB, 0, nSB);
  hnodihadtau->Sumw2();
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
