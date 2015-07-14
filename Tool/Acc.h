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


class BaseHistgram
{
 public:
  void BookHistgram(const char *);
  TFile *oFile;
  TH1D *hgen;
  TH1D *hacc;
  TH1D *hmet_gen;
  TH1D *hmet_acc;
  TH1D *hMT2_gen;
  TH1D *hMT2_acc;
  TH1D *hNbjet_gen;
  TH1D *hNbjet_acc;
  TH1D *hNtop_gen;
  TH1D *hNtop_acc;

  const TString title1 = "Generated muon ";
  const TString title2 = "Accepted muon ";
};
void BaseHistgram::BookHistgram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");
  hgen = new TH1D("hgen", title1+";Search bin;Events", 65, -0.5, 64.5);
  hgen->Sumw2();
  hacc = new TH1D("hacc", title2+";Search bin;Events", 65, -0.5, 64.5);
  hacc->Sumw2();
  hmet_gen = new TH1D("hmet_gen", "met", 50, 200, 1200);
  hmet_gen->Sumw2();
  hmet_acc = new TH1D("hmet_acc", "met", 50, 200, 1200);
  hmet_acc->Sumw2();
  hMT2_gen = new TH1D("hMT2_gen", "MT2", 20, 0, 500);
  hMT2_gen->Sumw2();
  hMT2_acc = new TH1D("hMT2_acc", "MT2", 20, 0, 500);
  hMT2_acc->Sumw2();
  hNbjet_gen = new TH1D("hNbjet_gen", "Nbjet", 5, 0, 5);
  hNbjet_gen->Sumw2();
  hNbjet_acc = new TH1D("hNbjet_acc", "Nbjet", 5, 0, 5);
  hNbjet_acc->Sumw2();
  hNtop_gen = new TH1D("hNtop_gen", "Ntop", 5, 0, 5);
  hNtop_gen->Sumw2();
  hNtop_acc = new TH1D("hNtop_acc", "Ntop", 5, 0, 5);
  hNtop_acc->Sumw2();
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


