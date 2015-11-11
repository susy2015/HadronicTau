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
using namespace std;

static const int nSB = 45; //We use 45 serach bins depending on Nbjet, Ntop, met and MT2 value.
static const int nTB = nSB + 2;// one extra bin for baseline and another bin for MT2 value less than 200 GeV

static BaselineVessel *AccBaselineVessel;
void AccpassBaselineFunc(NTupleReader& tr)
{
  (*AccBaselineVessel)(tr);
}

AnaSamples::SampleSet        allSamples;
AnaSamples::SampleCollection allCollections(allSamples);

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

  const TString title1 = "Generated muon ";
  const TString title2 = "Accepted muon ";

  const double jetbins[7] = {4, 5, 6, 7, 8, 9, 10};
  const int njetbin = sizeof(jetbins)/sizeof(jetbins[0])-1;
  const double MT2bins[6] = {100, 150, 200, 300, 400, 800};
  const int nMT2bin = sizeof(MT2bins)/sizeof(MT2bins[0])-1;

};
void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_Acc"+index+".root";
  oFile = new TFile(filename, "recreate");
  hgen = new TH1D("hgen", title1+";Search bin;Events", nTB, -0.5, 46.5);
  hgen->Sumw2();
  hacc = new TH1D("hacc", title2+";Search bin;Events", nTB, -0.5, 46.5);
  hacc->Sumw2();
  hgen_wt = new TH1D("hgen_wt", title1+";Search bin;Events", nTB, -0.5, 46.5);
  hgen_wt->Sumw2();
  hacc_wt = new TH1D("hacc_wt", title2+";Search bin;Events", nTB, -0.5, 46.5);
  hacc_wt->Sumw2();

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
  hNjetMT2_gen = new TH2D("hNjetMT2_gen", "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins);
  hNjetMT2_gen->Sumw2();
  hNjetMT2_acc = new TH2D("hNjetMT2_acc", "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins);
  hNjetMT2_acc->Sumw2();
  hNjetMT2_gen_wt = new TH2D("hNjetMT2_gen_wt", "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins);
  hNjetMT2_gen_wt->Sumw2();
  hNjetMT2_acc_wt = new TH2D("hNjetMT2_acc_wt", "Njet vs MT2", njetbin, jetbins, nMT2bin, MT2bins);
  hNjetMT2_acc_wt->Sumw2();
}

/*bool FillChain(TChain *chain, const TString &inputFileList)
{
  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;
  if(!infile.is_open())
    {
      std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
      return false;
    }
  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1)
    {
      infile >> buffer;
      if(!infile.good()) break;
      //std::cout << "Adding tree from " << buffer.c_str() << std::endl;
      chain->Add(buffer.c_str());
    }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return true;
  }*/

bool FillChain(TChain* &chain, const char *sample, const char *subsample, const int& startfile, const int& filerun){
  bool find = false;
  TString samplename(sample), subsamplename(subsample);
  if(samplename == "null"){
    chain = new TChain(allSamples[subsample].treePath.c_str());
    if(allSamples[subsample] != allSamples.null())
      {
	allSamples[subsample].addFilesToChain(chain, startfile, filerun);
	find = true;
      }
  }
  else
    {
      for(const auto & filelist : allCollections){
	if(filelist.first!=samplename)continue;
	for(auto & file : filelist.second){
	  for(const auto & perST : allSamples ){
	    string perSubStr;
	    if(perST.second == file ) perSubStr = perST.first;
	    if(perSubStr!=subsamplename)continue;
	    find = true;
	    chain = new TChain(file.treePath.c_str());
	    file.addFilesToChain(chain, startfile, filerun);
	  }//file loop                                                                                                                   
	}//sample loop                                                                                                                     
      }//collection loop                                                                                                                     
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
