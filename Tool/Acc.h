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

static const int nSB = 37; //We use 45 serach bins depending on Nbjet, Ntop, met and MT2 value.

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
  TH1D *hNjet_gen;
  TH1D *hNjet_acc;
  TH1D *hNjet_gen_wt;
  TH1D *hNjet_acc_wt;
  TH1D *hmet_gen;
  TH1D *hmet_acc;
  TH1D *hMT2_gen;
  TH1D *hMT2_acc;
  TH1D *hNbjet_gen;
  TH1D *hNbjet_acc;
  TH1D *hNtop_gen;
  TH1D *hNtop_acc;
  TH1D *hht_gen;
  TH1D *hht_acc;
  TH1D *hTaugen;
  TH1D *hTauacc;
  TH1D *hTauNjet_gen;
  TH1D *hTauNjet_acc;
  TH1D *hTauNbjet_gen;
  TH1D *hTauNbjet_acc;
  TH1D *hTaumet_gen;
  TH1D *hTaumet_acc;
  TH1D *hTauMT2_gen;
  TH1D *hTauMT2_acc;
  TH1D *hTauNtop_gen;
  TH1D *hTauNtop_acc;
  TH2D *hTauNjetMT2_gen;
  TH2D *hTauNjetMT2_acc;
  TH2D *hNjetMT2_gen;
  TH2D *hNjetMT2_acc;
  TH2D *hNjetMT2_gen_wt;
  TH2D *hNjetMT2_acc_wt;
  TH2D *hNjetmet_gen;
  TH2D *hNjetmet_acc;
  TH2D *hNjetmet_gen_wt;
  TH2D *hNjetmet_acc_wt;
  TH2D *hTauNjetmet_gen;
  TH2D *hTauNjetmet_acc;

  const TString title1 = "Generated muon ";
  const TString title2 = "Accepted muon ";

  const double jetbins[7] = {4, 5, 6, 7, 8, 9, 10};
  const int njetbin = sizeof(jetbins)/sizeof(jetbins[0])-1;
  const double MT2bins[6] = {200, 250, 300, 350, 400, 800};
  const int nMT2bin = sizeof(MT2bins)/sizeof(MT2bins[0])-1;
  const double metbins[5] = {200, 275, 350, 450, 800};
  const int nmetbin = sizeof(metbins)/sizeof(metbins[0])-1;
  const double topbins[3] = {1, 2, 3};
  const int ntopbin = sizeof(topbins)/sizeof(topbins[0])-1;

  std::vector<TH2*>hNjetMT2_met_gen;
  std::vector<TH2*>hNjetMT2_met_acc;
  std::vector<TH2*>hNjetMT2_top_gen;
  std::vector<TH2*>hNjetMT2_top_acc;
  std::vector<TH2*>hTauNjetMT2_met_gen;
  std::vector<TH2*>hTauNjetMT2_met_acc;
  
  TH1D *hgen_pdfCentral;
  TH1D *hgen_pdfUp;
  TH1D *hgen_pdfDown;
  TH1D *hgen_scaleUp;
  TH1D *hgen_scaleDown;
  TH1D *hacc_pdfCentral;
  TH1D *hacc_pdfUp;
  TH1D *hacc_pdfDown;
  TH1D *hacc_scaleUp;
  TH1D *hacc_scaleDown;
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

  TString Title1(unsigned int i);
  TString Title2(unsigned int j);
  TString Title5(unsigned int j, TString spec);
  TString TauTitle1(unsigned int i);
  TString TauTitle2(unsigned int j);
};
void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_Acc"+index+".root";
  oFile = new TFile(filename, "recreate");
  hgen = new TH1D("hgen", title1+";Search bin;Events", nSB, 0, 37);
  hgen->Sumw2();
  hacc = new TH1D("hacc", title2+";Search bin;Events", nSB, 0, 37);
  hacc->Sumw2();
  hgen_wt = new TH1D("hgen_wt", title1+";Search bin;Events", nSB, 0, 37);
  hgen_wt->Sumw2();
  hacc_wt = new TH1D("hacc_wt", title2+";Search bin;Events", nSB, 0, 37);
  hacc_wt->Sumw2();

  hgen_pdfCentral = new TH1D("hgen_pdfCentral", title1+";Search bin;Events", nSB, 0, 37);
  hgen_pdfCentral->Sumw2();
  hgen_pdfUp = new TH1D("hgen_pdfUp", title1+";Search bin;Events", nSB, 0, 37);
  hgen_pdfUp->Sumw2();
  hgen_pdfDown = new TH1D("hgen_pdfDown", title1+";Search bin;Events", nSB, 0, 37);
  hgen_pdfDown->Sumw2();
  hgen_scaleUp = new TH1D("hgen_scaleUp", title1+";Search bin;Events", nSB, 0, 37);
  hgen_scaleUp->Sumw2();
  hgen_scaleDown = new TH1D("hgen_scaleDown", title1+";Search bin;Events", nSB, 0, 37);
  hgen_scaleDown->Sumw2();
  hacc_pdfCentral = new TH1D("hacc_pdfCentral", title1+";Search bin;Events", nSB, 0, 37);
  hacc_pdfCentral->Sumw2();
  hacc_pdfUp = new TH1D("hacc_pdfUp", title1+";Search bin;Events", nSB, 0, 37);
  hacc_pdfUp->Sumw2();
  hacc_pdfDown = new TH1D("hacc_pdfDown", title1+";Search bin;Events", nSB, 0, 37);
  hacc_pdfDown->Sumw2();
  hacc_scaleUp = new TH1D("hacc_scaleUp", title1+";Search bin;Events", nSB, 0, 37);
  hacc_scaleUp->Sumw2();
  hacc_scaleDown = new TH1D("hacc_scaleDown", title1+";Search bin;Events", nSB, 0, 37);
  hacc_scaleDown->Sumw2();
  hgen_pdfCentral_wt = new TH1D("hgen_pdfCentral_wt", title1+";Search bin;Events", nSB, 0, 37);
  hgen_pdfCentral_wt->Sumw2();
  hgen_pdfUp_wt = new TH1D("hgen_pdfUp_wt", title1+";Search bin;Events", nSB, 0, 37);
  hgen_pdfUp_wt->Sumw2();
  hgen_pdfDown_wt = new TH1D("hgen_pdfDown_wt", title1+";Search bin;Events", nSB, 0, 37);
  hgen_pdfDown_wt->Sumw2();
  hgen_scaleUp_wt = new TH1D("hgen_scaleUp_wt", title1+";Search bin;Events", nSB, 0, 37);
  hgen_scaleUp_wt->Sumw2();
  hgen_scaleDown_wt = new TH1D("hgen_scaleDown_wt", title1+";Search bin;Events", nSB, 0, 37);
  hgen_scaleDown_wt->Sumw2();
  hacc_pdfCentral_wt = new TH1D("hacc_pdfCentral_wt", title1+";Search bin;Events", nSB, 0, 37);
  hacc_pdfCentral_wt->Sumw2();
  hacc_pdfUp_wt = new TH1D("hacc_pdfUp_wt", title1+";Search bin;Events", nSB, 0, 37);
  hacc_pdfUp_wt->Sumw2();
  hacc_pdfDown_wt = new TH1D("hacc_pdfDown_wt", title1+";Search bin;Events", nSB, 0, 37);
  hacc_pdfDown_wt->Sumw2();
  hacc_scaleUp_wt = new TH1D("hacc_scaleUp_wt", title1+";Search bin;Events", nSB, 0, 37);
  hacc_scaleUp_wt->Sumw2();
  hacc_scaleDown_wt = new TH1D("hacc_scaleDown_wt", title1+";Search bin;Events", nSB, 0, 37);
  hacc_scaleDown_wt->Sumw2();

  hNjet_gen = new TH1D("hNjet_gen", "Njet", 6, 4, 10);                                                                                       
  hNjet_gen->Sumw2();                                                                                                                         
  hNjet_acc = new TH1D("hNjet_acc", "Njet", 6, 4, 10);                                                                                       
  hNjet_acc->Sumw2();
  hNjet_gen_wt = new TH1D("hNjet_gen_wt", "Njet", 6, 4, 10);
  hNjet_gen_wt->Sumw2();
  hNjet_acc_wt = new TH1D("hNjet_acc_wt", "Njet", 6, 4, 10);
  hNjet_acc_wt->Sumw2();

  hmet_gen = new TH1D("hmet_gen", "met", 50, 200, 1200);
  hmet_gen->Sumw2();
  hmet_acc = new TH1D("hmet_acc", "met", 50, 200, 1200);
  hmet_acc->Sumw2();
  hMT2_gen = new TH1D("hMT2_gen", "MT2", 20, 0, 500);
  hMT2_gen->Sumw2();
  hMT2_acc = new TH1D("hMT2_acc", "MT2", 20, 0, 500);
  hMT2_acc->Sumw2();
  hNbjet_gen = new TH1D("hNbjet_gen", "Nbjet", 4, 1, 5);
  hNbjet_gen->Sumw2();
  hNbjet_acc = new TH1D("hNbjet_acc", "Nbjet", 4, 1, 5);
  hNbjet_acc->Sumw2();
  hNtop_gen = new TH1D("hNtop_gen", "Ntop", 4, 1, 5);
  hNtop_gen->Sumw2();
  hNtop_acc = new TH1D("hNtop_acc", "Ntop", 4, 1, 5);
  hNtop_acc->Sumw2();
  hht_gen = new TH1D("hht_gen", "ht", 10, 500, 1000);
  hht_gen->Sumw2();
  hht_acc = new TH1D("hht_acc", "ht", 10, 500, 1000);
  hht_acc->Sumw2();
  
  hTaugen = new TH1D("hTaugen", "Gen Tau;Search bin;Events", nSB, 0, 37);
  hTaugen->Sumw2();
  hTauacc = new TH1D("hTauacc", "Acc Tau;Search bin;Events", nSB, 0, 37);
  hTauacc->Sumw2();
  hTauNjet_gen = new TH1D("hTauNjet_gen", "Njet", 6, 4, 10);
  hTauNjet_gen->Sumw2();
  hTauNjet_acc = new TH1D("hTauNjet_acc", "Njet", 6, 4, 10);
  hTauNjet_acc->Sumw2();
  hTaumet_gen = new TH1D("hTaumet_gen", "met", 50, 200, 1200);
  hTaumet_gen->Sumw2();
  hTaumet_acc = new TH1D("hTaumet_acc", "met", 50, 200, 1200);
  hTaumet_acc->Sumw2();
  hTauMT2_gen = new TH1D("hTauMT2_gen", "MT2", 20, 0, 500);
  hTauMT2_gen->Sumw2();
  hTauMT2_acc = new TH1D("hTauMT2_acc", "MT2", 20, 0, 500);
  hTauMT2_acc->Sumw2();
  hTauNbjet_gen = new TH1D("hTauNbjet_gen", "Nbjet", 4, 1, 5);
  hTauNbjet_gen->Sumw2();
  hTauNbjet_acc = new TH1D("hTauNbjet_acc", "Nbjet", 4, 1, 5);
  hTauNbjet_acc->Sumw2();
  hTauNtop_gen = new TH1D("hTauNtop_gen", "Ntop", 4, 1, 5);
  hTauNtop_gen->Sumw2();
  hTauNtop_acc = new TH1D("hTauNtop_acc", "Ntop", 4, 1, 5);
  hTauNtop_acc->Sumw2();
  hTauNjetMT2_gen = new TH2D("hTauNjetMT2_gen", "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins);
  hTauNjetMT2_gen->Sumw2();
  hTauNjetMT2_acc = new TH2D("hTauNjetMT2_acc", "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins);
  hTauNjetMT2_acc->Sumw2();
  hTauNjetmet_gen = new TH2D("hTauNjetmet_gen", "Njet vs met", njetbin, jetbins, nmetbin, metbins);
  hTauNjetmet_gen->Sumw2();
  hTauNjetmet_acc = new TH2D("hTauNjetmet_acc", "Njet vs met", njetbin, jetbins, nmetbin, metbins);
  hTauNjetmet_acc->Sumw2();
 
  hNjetMT2_gen = new TH2D("hNjetMT2_gen", "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins);
  hNjetMT2_gen->Sumw2();
  hNjetMT2_acc = new TH2D("hNjetMT2_acc", "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins);
  hNjetMT2_acc->Sumw2();
  hNjetMT2_gen_wt = new TH2D("hNjetMT2_gen_wt", "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins);
  hNjetMT2_gen_wt->Sumw2();
  hNjetMT2_acc_wt = new TH2D("hNjetMT2_acc_wt", "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins);
  hNjetMT2_acc_wt->Sumw2();
  hNjetmet_gen = new TH2D("hNjetmet_gen", "Njet vs met", njetbin, jetbins, nmetbin, metbins);
  hNjetmet_gen->Sumw2();
  hNjetmet_acc = new TH2D("hNjetmet_acc", "Njet vs met", njetbin, jetbins, nmetbin, metbins);
  hNjetmet_acc->Sumw2();
  hNjetmet_gen_wt = new TH2D("hNjetmet_gen_wt", "Njet vs met", njetbin, jetbins, nmetbin, metbins);
  hNjetmet_gen_wt->Sumw2();
  hNjetmet_acc_wt = new TH2D("hNjetmet_acc_wt", "Njet vs met", njetbin, jetbins, nmetbin, metbins);
  hNjetmet_acc_wt->Sumw2();
  for(unsigned int i = 0; i < nmetbin; ++i) {
    hNjetMT2_met_gen.push_back(new TH2D(Title1(i), "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins));
    hNjetMT2_met_gen.back()->Sumw2();

    hTauNjetMT2_met_gen.push_back(new TH2D(TauTitle1(i), "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins));
    hTauNjetMT2_met_gen.back()->Sumw2();


    hNjetMT2_met_acc.push_back(new TH2D(Title2(i), "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins));
    hNjetMT2_met_acc.back()->Sumw2();
    hTauNjetMT2_met_acc.push_back(new TH2D(TauTitle2(i), "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins));
    hTauNjetMT2_met_acc.back()->Sumw2();
  }
  for(unsigned int j = 0; j < ntopbin; ++j) {
    hNjetMT2_top_gen.push_back(new TH2D(Title5(j, "gen"), "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins));
    hNjetMT2_top_gen.back()->Sumw2();
    hNjetMT2_top_acc.push_back(new TH2D(Title5(j, "acc"), "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins));
    hNjetMT2_top_acc.back()->Sumw2();
  }
}
TString BaseHistgram::Title1(unsigned int i){
  TString title = "hNjetMT2_met_gen";
  title+=i;
  return title;
}
TString BaseHistgram::Title2(unsigned int j){
  TString title ="hNjetMT2_met_acc";
  title+=j;
  return title;
}
TString BaseHistgram::Title5(unsigned int j, TString spec){
  TString title ="hNjetMT2_top";
  title+="_"+spec;
  title+=j;
  return title;
}
TString BaseHistgram::TauTitle1(unsigned int i){
  TString title = "hTauNjetMT2_met_gen";
  title+=i;
  return title;
}
TString BaseHistgram::TauTitle2(unsigned int j){
  TString title ="hTauNjetMT2_met_acc";
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
