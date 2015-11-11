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

AnaSamples::SampleSet        allSamples;
AnaSamples::SampleCollection allCollections(allSamples);

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

double deltaRmax(double pt);
