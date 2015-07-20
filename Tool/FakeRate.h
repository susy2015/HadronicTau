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
class BaseHistgram
{
 public:
  void BookHistgram(const char *);
  TFile *oFile;
  TH1D *htauBjetPt_den;
  TH1D *htauBjetPt_num;
  TH1D *htauBjetEta_den;
  TH1D *htauBjetEta_num;
  TH1D *hEff_Pt;
  TH1D *hEff_Eta;
  double bins[9] = {0, 50, 75, 100, 125, 150, 175, 200, 1000};
  int binnum = 8;
};
void BaseHistgram::BookHistgram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");
   htauBjetPt_den = new TH1D("htauBjetPt_den", "Hadtau jet pT", 50, 20, 1020);
  //  htauBjetPt_den = new TH1D("htauBjetPt_den", "Hadtau jet pT",binnum, bins);
  htauBjetPt_den->Sumw2();
  htauBjetEta_den = new TH1D("htauBjetEta_den", "Hadtau jet eta", 50, -2.4, 2.4);
  htauBjetEta_den->Sumw2(); 
    htauBjetPt_num = new TH1D("htauBjetPt_num", "Hadtau Bjet pT", 50, 20, 1020);
  //htauBjetPt_num = new TH1D("htauBjetPt_num", "Hadtau Bjet pT", binnum, bins);
  htauBjetPt_num->Sumw2();  
  htauBjetEta_num = new TH1D("htauBjetEta_num", "Hadtau Bjet eta", 50, -2.4, 2.4);
  htauBjetEta_num->Sumw2();
}
bool FillChain(TChain *chain, const TString &inputFileList)
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
}
double deltaRmax(double pt);
