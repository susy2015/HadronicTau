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
#include "SusyAnaTools/Tools/samples.h"

static const int nSB = 84; //We use nSB serach bins depending on Nbjet, Ntop, met and MT2 value.

using namespace std;
BaselineVessel *ExpBaseline = 0;
double Lumiscale = 1.0;
double EventWeight = 1.0;
class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;
  TH1D *hMET_tau;
  TH1D *hNbJets_tau;
  TH1D *hNTops_tau;
  TH1D *hMT2_tau;
  TH1D *hNJets_tau;
  TH1D *hHT_tau;
  TH1D *hYields_tau;
  TH1D *hdPhi0_tau;
  TH1D *hdPhi1_tau;
  TH1D *hdPhi2_tau;

  TH1D *hMET_LL;
  TH1D *hNbJets_LL;
  TH1D *hNTops_LL;
  TH1D *hMT2_LL;
  TH1D *hNJets_LL;
  TH1D *hHT_LL;
  TH1D *hYields_LL;
  TH1D *hdPhi0_LL;
  TH1D *hdPhi1_LL;
  TH1D *hdPhi2_LL;

  const TString title = "HadTau MC";
  const TString title_LL = "LL MC";
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_HadTauLL"+index+".root";
  oFile = new TFile(filename, "recreate");

  hMET_tau = new TH1D("hMET_tau",title+";met [GeV];Events",24,200.,800.);
  hMET_tau->Sumw2();
  hNbJets_tau = new TH1D("hNbJets_tau",title+";N_{bjets};Events",4, 1, 5);
  hNbJets_tau->Sumw2();
  hNTops_tau = new TH1D("hNTops_tau",title+";N_{tops};Events",4, 1, 5);
  hNTops_tau->Sumw2();
  hMT2_tau = new TH1D("hMT2_tau",title+";M_{T2}[GeV];Events",12,200,500);
  hMT2_tau->Sumw2();
  hYields_tau = new TH1D("hYields_tau", title+";search bin;Events",nSB,0,nSB);
  hYields_tau->Sumw2();
  hNJets_tau = new TH1D("hNJets_tau",title+";N_{jets};Events",6 ,4,10);
  hNJets_tau->Sumw2();
  hHT_tau = new TH1D("hHT_tau",title+";H_{T} [GeV];Events",20,500.,1000.);
  hHT_tau->Sumw2();  
  hdPhi0_tau = new TH1D("hdPhi0_tau", title+";dPhi0;Events", 16, 0, 3.2);
  hdPhi0_tau->Sumw2();
  hdPhi1_tau = new TH1D("hdPhi1_tau", title+";dPhi1;Events", 16, 0, 3.2);
  hdPhi1_tau->Sumw2();
  hdPhi2_tau = new TH1D("hdPhi2_tau", title+";dPhi2;Events", 16, 0, 3.2);
  hdPhi2_tau->Sumw2();
  
  hMET_LL = new TH1D("hMET_LL",title_LL+";met [GeV];Events",24,200.,800.);
  hMET_LL->Sumw2();
  hNbJets_LL = new TH1D("hNbJets_LL",title_LL+";N_{bjets};Events",4, 1, 5);
  hNbJets_LL->Sumw2();
  hNTops_LL = new TH1D("hNTops_LL",title_LL+";N_{tops};Events",4, 1, 5);
  hNTops_LL->Sumw2();
  hMT2_LL = new TH1D("hMT2_LL",title_LL+";M_{T2}[GeV];Events",12,200,500);
  hMT2_LL->Sumw2();
  hYields_LL = new TH1D("hYields_LL", title_LL+";search bin;Events",nSB,0,nSB);
  hYields_LL->Sumw2();
  hNJets_LL = new TH1D("hNJets_LL",title_LL+";N_{jets};Events",6 ,4,10);
  hNJets_LL->Sumw2();
  hHT_LL = new TH1D("hHT_LL",title_LL+";H_{T} [GeV];Events",20,500.,1000.);
  hHT_LL->Sumw2();  
  hdPhi0_LL = new TH1D("hdPhi0_LL", title_LL+";dPhi0;Events", 16, 0, 3.2);
  hdPhi0_LL->Sumw2();
  hdPhi1_LL = new TH1D("hdPhi1_LL", title_LL+";dPhi1;Events", 16, 0, 3.2);
  hdPhi1_LL->Sumw2();
  hdPhi2_LL = new TH1D("hdPhi2_LL", title_LL+";dPhi2;Events", 16, 0, 3.2);
  hdPhi2_LL->Sumw2();
   
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

