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
#include "TauResponse.h"

static const int nSB = 37; //We use 37 serach bins depending on Nbjet, Ntop, met and MT2 value.                                               

using namespace std;

static BaselineVessel *IsoTrackBaselineVessel;
void passBaselineFuncIsoTrack(NTupleReader& tr)
{
  (*IsoTrackBaselineVessel)(tr);
}

double Lumiscale = 1.0;

class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;
  TH1D *hSB_num;
  TH1D *hSB_den;
  TH1D *hmet_num;
  TH1D *hmet_den;
  TH1D *hht_num;
  TH1D *hht_den;
  TH1D *hMT2_num;
  TH1D *hMT2_den;
  TH1D *hNbjet_num;
  TH1D *hNbjet_den;
  TH1D *hNtop_num;
  TH1D *hNtop_den;
  TH1D *hNjet_num;
  TH1D *hNjet_den;
  TH2D *hNjetNbjet_num;
  TH2D *hNjetNbjet_den;
  TH1D *hSB_num1pr;
  TH1D *hSB_den1pr;
  TH1D *hSB_num3pr;
  TH1D *hSB_den3pr;
  const TString title1 = "With IsoTrackVeto";
  const TString title2 = "Without IsoTrackVeto";
 
 const double jetbins[7] = {4, 5, 6, 7, 8, 9, 10};
 const int njetbin = sizeof(jetbins)/sizeof(jetbins[0])-1;
 const double bjetbins[3] = {1, 2, 3};
 const int nbjetbin = sizeof(bjetbins)/sizeof(bjetbins[0])-1;
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_IsoTrack"+index+".root";
  oFile = new TFile(filename, "recreate");
  hSB_num = new TH1D("hSB_num", title1+";Search bin;Events", nSB, 0, 37);
  hSB_num->Sumw2();
  hSB_den = new TH1D("hSB_den", title2+";Search bin;Events", nSB, 0, 37);
  hSB_den->Sumw2();
  hmet_den = new TH1D("hmet_den", "met", 50, 200, 1200);
  hmet_den->Sumw2();
  hmet_num = new TH1D("hmet_num", "met", 50, 200, 1200);
  hmet_num->Sumw2();
  hht_den = new TH1D("hht_den", "ht", 10, 500, 1000);
  hht_den->Sumw2();
  hht_num = new TH1D("hht_num", "ht", 10, 500, 1000);
  hht_num->Sumw2();
  hMT2_den = new TH1D("hMT2_den", "MT2", 10, 100, 500);
  hMT2_den->Sumw2();
  hMT2_num = new TH1D("hMT2_num", "MT2", 10, 100, 500);
  hMT2_num->Sumw2();
  hNbjet_den = new TH1D("hNbjet_den", "Nbjet", 5, 1, 6);
  hNbjet_den->Sumw2();
  hNbjet_num = new TH1D("hNbjet_num", "Nbjet", 5, 1, 6);
  hNbjet_num->Sumw2();
  hNtop_den = new TH1D("hNtop_den", "Ntop", 5, 0, 5);
  hNtop_den->Sumw2();
  hNtop_num = new TH1D("hNtop_num", "Ntop", 5, 0, 5);
  hNtop_num->Sumw2();
  hNjet_den = new TH1D("hNjet_den", "Njet", 6, 4, 10);
  hNjet_den->Sumw2();
  hNjet_num = new TH1D("hNjet_num", "Njet", 6, 4, 10);
  hNjet_num->Sumw2();
  hNjetNbjet_num = new TH2D("hNjetNbjet_num", "Njet_Nbjet", njetbin, jetbins, nbjetbin, bjetbins);
  hNjetNbjet_num->Sumw2();
  hNjetNbjet_den = new TH2D("hNjetNbjet_den", "Njet_Nbjet", njetbin, jetbins, nbjetbin, bjetbins);
  hNjetNbjet_den->Sumw2();
  hSB_num1pr = new TH1D("hSB_num1pr", title1+";Search bin;Events", nSB, 0, 37);
  hSB_num1pr->Sumw2();
  hSB_den1pr = new TH1D("hSB_den1pr", title2+";Search bin;Events", nSB, 0, 37);
  hSB_den1pr->Sumw2();
  hSB_num3pr = new TH1D("hSB_num3pr", title1+";Search bin;Events", nSB, 0, 37);
  hSB_num3pr->Sumw2();
  hSB_den3pr = new TH1D("hSB_den3pr", title2+";Search bin;Events", nSB, 0, 37);
  hSB_den3pr->Sumw2();
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

double deltaRmax(double pt);
bool findTauMatchedisoJet(int matchedObjIdx, const std::vector<TLorentzVector> &genvisiblehadtauLVec, const std::vector<TLorentzVector> &jetsLVec, std::vector<int> &IsotrkMatchedHadtau);
bool passIsoTrks(const TLorentzVector isoTrksLVec, const double isoTrksIso, const double isoTrksMtw, const int isoTrkspdgId);
void Fill2DIsoTrk(TH2 *hist, const int &a, const int &b, const double &w){
  int nbinx = hist->GetNbinsX();
  int nbiny = hist->GetNbinsY();
  int lowx = hist->GetXaxis()->GetBinLowEdge(nbinx);
  int highx = hist->GetXaxis()->GetBinLowEdge(nbinx + 1);
  int lowy = hist->GetYaxis()->GetBinLowEdge(nbiny);
  int highy = hist->GetYaxis()->GetBinLowEdge(nbiny + 1);
  int copyx = a;
  if(copyx >= highx) copyx = lowx;
  int copyy = b;
  if(copyy >= highy) copyy = lowy;
  hist->Fill(copyx, copyy, w);
}
