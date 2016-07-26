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

double Lumiscale = 1.0;

class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;
  std::vector<TH1*> hTauResp;
  std::vector<TH1*> hTauvisible;
  std::vector<TH1*> hTaugenerated;
  std::vector<TH1*> hTauRespUp;
  std::vector<TH1*> hTauRespDown;
  TString Title1(int i);
  TString Title2(int j);
  TString Title3(int i);
  TString Title4(int j);
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

    hTauRespUp.push_back(new TH1D(Title3(i),";p_{T}(visible) / p_{T}(generated);Count",50,0.,2.5));
    hTauRespUp.back()->Sumw2();
    hTauRespDown.push_back(new TH1D(Title4(i),";p_{T}(visible) / p_{T}(generated);Count",50,0.,2.5));
    hTauRespDown.back()->Sumw2();
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

TString BaseHistgram::Title3(int j){
  TString title = "hTauRespUp_";
  title += j;
  return title;
}
TString BaseHistgram::Title4(int j){
  TString title = "hTauRespDown_";
  title += j;
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

double deltaRmax(double pt);
