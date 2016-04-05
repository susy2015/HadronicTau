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

using namespace std;

static BaselineVessel *FakeRateBaselineVessel;
void passBaselineFuncFakeRate(NTupleReader& tr)
{
  (*FakeRateBaselineVessel)(tr);
}

double Lumiscale = 1.0;

class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;
  TH1D *htauBjetPt_den;
  TH1D *htauBjetPt_num;
  TH1D *htauBjetEta_den;
  TH1D *htauBjetEta_num;
  TH1D *hEff_Pt;
  TH1D *hEff_Eta;
  TH1D *htaub;
  TH1D *htaujetb;
  TH2D *htaub_dRpt;
  TH1D *htaub1;
  TH1D *htaub2;
  TH1D *htaub3;
  TH1D *htaub4;
  TH1D *htaub5;
  double bins[18] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 500, 1000};
  int binnum = 17;
};
void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_FakeRate"+index+".root";
  oFile = new TFile(filename, "recreate");
  //htauBjetPt_den = new TH1D("htauBjetPt_den", "Hadtau jet pT", 50, 20, 1020);
  htauBjetPt_den = new TH1D("htauBjetPt_den", "Hadtau jet pT",binnum, bins);
  htauBjetPt_den->Sumw2();
  htauBjetEta_den = new TH1D("htauBjetEta_den", "Hadtau jet eta", 50, -2.4, 2.4);
  htauBjetEta_den->Sumw2(); 
  //  htauBjetPt_num = new TH1D("htauBjetPt_num", "Hadtau Bjet pT", 50, 20, 1020);
  htauBjetPt_num = new TH1D("htauBjetPt_num", "Hadtau Bjet pT", binnum, bins);
  htauBjetPt_num->Sumw2();  
  htauBjetEta_num = new TH1D("htauBjetEta_num", "Hadtau Bjet eta", 50, -2.4, 2.4);
  htauBjetEta_num->Sumw2();
  
  htaub = new TH1D("htaub", "dR between hadtau and b", 10, 0, 1);
  htaub->Sumw2();
  htaujetb = new TH1D("htaujetb", "dR between taujet and b", 10, 0, 1);
  htaujetb->Sumw2();
  htaub_dRpt = new TH2D("htaub_dRpt", "dR vs hadtau pt", binnum, bins, 10, 0, 1);
  htaub_dRpt->Sumw2();
  htaub1 = new TH1D("htaub1", "dR(tau,b) for 30-50 Pt", 100, 0, 4);
  htaub1->Sumw2();
  htaub2 = new TH1D("htaub2", "dR(tau,b) for 50-100 Pt", 100, 0, 4);
  htaub2->Sumw2();
  htaub3 = new TH1D("htaub3", "dR(tau,b) for 100-200 Pt", 100, 0, 4);
  htaub3->Sumw2();
  htaub4 = new TH1D("htaub4", "dR(tau,b) for 200-400 Pt", 100, 0, 4);
  htaub4->Sumw2();
  htaub5 = new TH1D("htaub5", "dR(tau,b) for >400 Pt", 100, 0, 4);
  htaub5->Sumw2();
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

double deltaRmax(double pt);
void FillDoubleFakeRate(TH1 *hist, const double &a, const double &w){
  int nbin = hist->GetNbinsX();
  double low = hist->GetBinLowEdge(nbin);
  double high = hist->GetBinLowEdge(nbin + 1);
  double copy = a;
  if(copy >= high) copy = low;
  hist->Fill(copy, w);
}
