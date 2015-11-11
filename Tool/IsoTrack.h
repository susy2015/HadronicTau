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

static const int nSB = 45; //We use 45 serach bins depending on Nbjet, Ntop, met and MT2 value.                                               
static const int nTB = nSB + 2;// one extra bin for baseline and another bin for MT2 value less than 200 GeV                                 

using namespace std;

static BaselineVessel *IsoTrackBaselineVessel;
void passBaselineFuncIsoTrack(NTupleReader& tr)
{
  (*IsoTrackBaselineVessel)(tr);
}

AnaSamples::SampleSet        allSamples;
AnaSamples::SampleCollection allCollections(allSamples);

class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;
  TH1D *hSB_num;
  TH1D *hSB_den;
  TH1D *hmet_num;
  TH1D *hmet_den;
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

  const TString title1 = "With IsoTrackVeto";
  const TString title2 = "Without IsoTrackVeto";
 
 const double jetbins[8] = {4, 5, 6, 7, 8, 9, 10, 100};
 const int njetbin = sizeof(jetbins)/sizeof(jetbins[0])-1;
 const double bjetbins[5] = {1, 2, 3, 4, 10};
 const int nbjetbin = sizeof(bjetbins)/sizeof(bjetbins[0])-1;
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_IsoTrack"+index+".root";
  oFile = new TFile(filename, "recreate");
  hSB_num = new TH1D("hSB_num", title1+";Search bin;Events", nTB, -0.5, 46.5);
  hSB_num->Sumw2();
  hSB_den = new TH1D("hSB_den", title2+";Search bin;Events", nTB, -0.5, 46.5);
  hSB_den->Sumw2();
  hmet_den = new TH1D("hmet_den", "met", 50, 200, 1200);
  hmet_den->Sumw2();
  hmet_num = new TH1D("hmet_num", "met", 50, 200, 1200);
  hmet_num->Sumw2();
  hMT2_den = new TH1D("hMT2_den", "MT2", 20, 0, 500);
  hMT2_den->Sumw2();
  hMT2_num = new TH1D("hMT2_num", "MT2", 20, 0, 500);
  hMT2_num->Sumw2();
  hNbjet_den = new TH1D("hNbjet_den", "Nbjet", 5, 1, 6);
  hNbjet_den->Sumw2();
  hNbjet_num = new TH1D("hNbjet_num", "Nbjet", 5, 1, 6);
  hNbjet_num->Sumw2();
  hNtop_den = new TH1D("hNtop_den", "Ntop", 5, 0, 5);
  hNtop_den->Sumw2();
  hNtop_num = new TH1D("hNtop_num", "Ntop", 5, 0, 5);
  hNtop_num->Sumw2();
  hNjet_den = new TH1D("hNjet_den", "Njet", 10, 4, 14);
  hNjet_den->Sumw2();
  hNjet_num = new TH1D("hNjet_num", "Njet", 10, 4, 14);
  hNjet_num->Sumw2();
  hNjetNbjet_num = new TH2D("hNjetNbjet_num", "Njet_Nbjet", njetbin, jetbins, nbjetbin, bjetbins);
  hNjetNbjet_num->Sumw2();
  hNjetNbjet_den = new TH2D("hNjetNbjet_den", "Njet_Nbjet", njetbin, jetbins, nbjetbin, bjetbins);
  hNjetNbjet_den->Sumw2();
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
bool findTauMatchedisoJet(int matchedObjIdx, const std::vector<TLorentzVector> &genvisiblehadtauLVec, const std::vector<TLorentzVector> &jetsLVec, std::vector<int> &IsotrkMatchedHadtau);
bool passIsoTrks(const TLorentzVector isoTrksLVec, const double isoTrksIso, const double isoTrksMtw, const int isoTrkspdgId);
