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
#include "SusyAnaTools/Tools/PDFUncertainty.h"

using namespace std;

static const int nSB = 59; //We use 69 serach bins depending on Nbjet, Ntop, met and MT2 value.

static BaselineVessel *AccBaselineVessel;
void AccpassBaselineFunc(NTupleReader& tr)
{
  (*AccBaselineVessel)(tr);
}
//Add PDF uncertainty                                                                                                                         
static PDFUncertainty pdf;
void mypdf(NTupleReader& tr)
{
  pdf(tr);
}

double Lumiscale = 1.0;


class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;
  TH1D *hgen;
  TH1D *hacc;
  TH1D *hgen_wt;
  TH1D *hacc_wt;
  TH1D *hTauacc;
  TH1D *hTaugen;

  const TString title1 = "Generated muon ";
  const TString title2 = "Accepted muon ";

  const double jetbins[7] = {4, 5, 6, 7, 8, 9, 10};
  const int njetbin = sizeof(jetbins)/sizeof(jetbins[0])-1;
  const double metbins[5] = {200, 275, 350, 450, 800};
  const int nmetbin = sizeof(metbins)/sizeof(metbins[0])-1;

  TH1D *hgen_pdfCentral_wt;
  TH1D *hgen_pdfUp_wt;
  TH1D *hgen_pdfDown_wt;
  TH1D *hgen_scaleUp_wt;
  TH1D *hgen_scaleDown_wt;
  TH1D *hacc_pdfCentral_wt;
  TH1D *hacc_pdfUp_wt;
  TH1D *hacc_pdfDown_wt;
  TH1D *hacc_scaleUp_wt;
  TH1D *hacc_scaleDown_wt;

  //PDF & Scale check
  TH1D *hpdfUp;
  TH1D *hpdfDown;
  TH1D *hscaleUp;
  TH1D *hscaleDown;
  
};
void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_Acc"+index+".root";
  oFile = new TFile(filename, "recreate");
  hgen = new TH1D("hgen", title1+";Search bin;Events", nSB, 0, nSB);
  hgen->Sumw2();
  hacc = new TH1D("hacc", title2+";Search bin;Events", nSB, 0, nSB);
  hacc->Sumw2();
  hgen_wt = new TH1D("hgen_wt", title1+";Search bin;Events", nSB, 0, nSB);
  hgen_wt->Sumw2();
  hacc_wt = new TH1D("hacc_wt", title2+";Search bin;Events", nSB, 0, nSB);
  hacc_wt->Sumw2();

  hgen_pdfCentral_wt = new TH1D("hgen_pdfCentral_wt", title1+";Search bin;Events", nSB, 0, nSB);
  hgen_pdfCentral_wt->Sumw2();
  hgen_pdfUp_wt = new TH1D("hgen_pdfUp_wt", title1+";Search bin;Events", nSB, 0, nSB);
  hgen_pdfUp_wt->Sumw2();
  hgen_pdfDown_wt = new TH1D("hgen_pdfDown_wt", title1+";Search bin;Events", nSB, 0, nSB);
  hgen_pdfDown_wt->Sumw2();
  hgen_scaleUp_wt = new TH1D("hgen_scaleUp_wt", title1+";Search bin;Events", nSB, 0, nSB);
  hgen_scaleUp_wt->Sumw2();
  hgen_scaleDown_wt = new TH1D("hgen_scaleDown_wt", title1+";Search bin;Events", nSB, 0, nSB);
  hgen_scaleDown_wt->Sumw2();
  hacc_pdfCentral_wt = new TH1D("hacc_pdfCentral_wt", title1+";Search bin;Events", nSB, 0, nSB);
  hacc_pdfCentral_wt->Sumw2();
  hacc_pdfUp_wt = new TH1D("hacc_pdfUp_wt", title1+";Search bin;Events", nSB, 0, nSB);
  hacc_pdfUp_wt->Sumw2();
  hacc_pdfDown_wt = new TH1D("hacc_pdfDown_wt", title1+";Search bin;Events", nSB, 0, nSB);
  hacc_pdfDown_wt->Sumw2();
  hacc_scaleUp_wt = new TH1D("hacc_scaleUp_wt", title1+";Search bin;Events", nSB, 0, nSB);
  hacc_scaleUp_wt->Sumw2();
  hacc_scaleDown_wt = new TH1D("hacc_scaleDown_wt", title1+";Search bin;Events", nSB, 0, nSB);
  hacc_scaleDown_wt->Sumw2();

  
  hTaugen = new TH1D("hTaugen", "Gen Tau;Search bin;Events", nSB, 0, nSB);
  hTaugen->Sumw2();
  hTauacc = new TH1D("hTauacc", "Acc Tau;Search bin;Events", nSB, 0, nSB);
  hTauacc->Sumw2();

  hpdfUp = new TH1D("hpdfUp", "pdfUp;pdfUp;Event", 100, -100, 100);
  hpdfDown = new TH1D("hpdfDown", "pdfDown;pdfDown;Event", 100, -100, 100);
  hscaleUp = new TH1D("hscaleUp", "MCscaleUp;scaleUp;Event", 100, -100, 100);
  hscaleDown = new TH1D("hscaleDown", "MCscaleDown;scaleDown;Event", 100, -100, 100);

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

void FillDoubleAcc(TH1 *hist, const double &a, const double &w){
  int nbin = hist->GetNbinsX();
  double low = hist->GetBinLowEdge(nbin);
  double high = hist->GetBinLowEdge(nbin + 1);
  double copy = a;
  if(copy >= high) copy = low;
  hist->Fill(copy, w);
}
void FillIntAcc(TH1 *hist, const int &a, const double &w){
  int nbin = hist->GetNbinsX();
  int low = (int)hist->GetBinLowEdge(nbin);
  int high = (int)hist->GetBinLowEdge(nbin + 1);
  int copy = a;
  if(copy >= high) copy = low;
  hist->Fill(copy, w);
}
void Fill2DAcc(TH2 *hist, const int &a, const double &b, const double &w){
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

bool find_mother(int momIdx, int dauIdx, const std::vector<int> &genDecayIdxVec, const std::vector<int> &genDecayMomIdxVec);
int find_idx(int genIdx, const std::vector<int> &genDecayIdxVec){
  for(int ig=0; ig<(int)genDecayIdxVec.size(); ig++){
    if( genDecayIdxVec[ig] == genIdx ) return ig;
  }
  return -1;
}

bool passIsoTrks1(const TLorentzVector isoTrksLVec, const double isoTrksIso, const double isoTrksMtw, const int isoTrkspdgId);
bool findTauMatchedisoJet(int matchedObjIdx, const std::vector<TLorentzVector> &genvisiblehadtauLVec, const std::vector<TLorentzVector> &jetsLVec, std::vector<int> &IsotrkMatchedHadtau);
double deltaRmax(double pt);
