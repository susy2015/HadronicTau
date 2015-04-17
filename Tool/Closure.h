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

  TH1D *hPredHt;
  TH1D *hPredmet;
  TH1D *hPredNJets;
  TH1D *hTrueHt;
  TH1D *hTruemet;
  TH1D *hTrueNJets;
  const TString title = "Hadronic-Tau Closure Test";
};

void BaseHistgram::BookHistgram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");
  hPredHt = new TH1D("hPredHt",title+";H_{T} [GeV];Events",25,200.,2000.);
  hPredHt->Sumw2();
  hPredmet = new TH1D("hPredmet",title+";met [GeV];Events",50,200.,1200.);
  hPredmet->Sumw2();
  hPredNJets = new TH1D("hPredNJets",title+";N_{jets};Events",10,4,14);
  hPredNJets->Sumw2();
  hTrueHt = new TH1D("hTrueHt",title+";H_{T} [GeV];Events",25,200.,2000.);
  hTrueHt->Sumw2();
  hTruemet = new TH1D("hTruemet",title+";met [GeV];Events",50,200.,1200.);
  hTruemet->Sumw2();
  hTrueNJets = new TH1D("hTrueNJets",title+";N_{jets};Events",10,4,14);
  hTrueNJets->Sumw2();
 
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
double htJetPtMin(){ return 50;}
double htJetEtaMax() {return 2.4;}
double mhtJetPtMin(){return 30;}
double mhtJetEtaMax() {return 5;}
double nJetPtMin(){return 30;}
double nJetEtaMax() {return 2.4;}
std::vector<TLorentzVector> combjet (const std::vector<TLorentzVector> &seljet, const std::vector<TLorentzVector> &simjet);
std::vector<double>combJetBtag(const std::vector<double> &selJetBtag, const std::vector<TLorentzVector> &seljet, const std::vector<TLorentzVector> &simjet);
