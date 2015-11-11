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

AnaSamples::SampleSet        allSamples;
AnaSamples::SampleCollection allCollections(allSamples);

class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;
  std::vector<TH1*> hTauResp;
  std::vector<TH1*> hTauvisible;
  std::vector<TH1*> hTaugenerated;
  TString Title1(int i);
  TString Title2(int j);
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_Template"+index+".root";
  oFile = new TFile(filename, "recreate");
  for(unsigned int i = 0; i < TauResponse::nBins(); ++i) {
    hTauResp.push_back(new TH1D(TauResponse::name(i),";p_{T}(visible) / p_{T}(generated);Count",50,0.,2.5));
    hTauResp.back()->Sumw2();

    hTauvisible.push_back(new TH1D(Title1(i),";p_{T}(visible)[GeV];Events",50,0.,500.));
    hTauvisible.back()->Sumw2();

    hTaugenerated.push_back(new TH1D(Title2(i),";p_{T}(generated)[GeV];Events",50,0.,500.));
    hTaugenerated.back()->Sumw2();

  }
}

TString BaseHistgram::Title1(int i){
  TString title = "htauvisible_";
  title += i;
  return title;
}
TString BaseHistgram::Title2(int j){
  TString title = "htaugenerated_";
  title += j;
  return title;
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
