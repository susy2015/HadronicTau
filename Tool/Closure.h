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
  TH1D *hPredNbJets;
  TH1D *hPredNTops;
  TH1D *hPredMT2;
  TH1D *hPredmTcomb;

  TH1D *hTrueHt;
  TH1D *hTruemet;
  TH1D *hTrueNJets;
  TH1D *hTrueNbJets;
  TH1D *hTrueNTops;
  TH1D *hTrueMT2;
  TH1D *hTruemTcomb;
  const TString title = "Hadronic-Tau Closure Test";
};

void BaseHistgram::BookHistgram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");
  hPredHt = new TH1D("hPredHt",title+";H_{T} [GeV];Events",25,200.,2000.);
  hPredHt->Sumw2();
  hPredmet = new TH1D("hPredmet",title+";met [GeV];Events",24,200.,800.);
  hPredmet->Sumw2();
  hPredNJets = new TH1D("hPredNJets",title+";N_{jets};Events",10,4,14);
  hPredNJets->Sumw2();
  hPredNbJets = new TH1D("hPredNbJets",title+";N_{bjets};Events",5, 0, 5);
  hPredNbJets->Sumw2();
  hPredNTops = new TH1D("hPredNTops",title+";N_{tops};Events",5, 0, 5);
  hPredNTops->Sumw2();
  hPredMT2 = new TH1D("hPredMT2",title+";M_{T2}[GeV];Events",20,0,500);
  hPredMT2->Sumw2();
  hPredmTcomb = new TH1D("hPredmTcomb",title+";M_{Tb}+0.5*M_{Tt}[GeV];Events",50,0,1000);
  hPredmTcomb->Sumw2();

  hTrueHt = new TH1D("hTrueHt",title+";H_{T} [GeV];Events",25,200.,2000.);
  hTrueHt->Sumw2();
  hTruemet = new TH1D("hTruemet",title+";met [GeV];Events",24,200.,800.);
  hTruemet->Sumw2();
  hTrueNJets = new TH1D("hTrueNJets",title+";N_{jets};Events",10,4,14);
  hTrueNJets->Sumw2();
  hTrueNbJets = new TH1D("hTrueNbJets",title+";N_{bjets};Events",5, 0, 5);
  hTrueNbJets->Sumw2();
  hTrueNTops = new TH1D("hTrueNTops",title+";N_{tops};Events",5, 0, 5);
  hTrueNTops->Sumw2();
  hTrueMT2 = new TH1D("hTrueMT2",title+";M_{T2}[GeV];Events",20,0,500);
  hTrueMT2->Sumw2();
  hTruemTcomb = new TH1D("hTruemTcomb",title+";M_{Tb}+0.5*M_{Tt}[GeV];Events",50,0,1000);
  hTruemTcomb->Sumw2();
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

void drawOverFlowBin(TH1 *hist){

  int nbins = hist->GetXaxis()->GetNbins();

  double overflow = hist->GetBinContent(nbins+1);
  double lastCont = hist->GetBinContent(nbins);
  double ovrflweroor = hist->GetBinError(nbins+1);
  double lstbinerror = hist->GetBinError(nbins);
  hist->SetBinContent(nbins, overflow+lastCont);
  hist->SetBinError(nbins, TMath::Sqrt(ovrflweroor * ovrflweroor + lstbinerror * lstbinerror));
}
double calcMT(const TLorentzVector &objLVec, const TLorentzVector &metLVec){



  const double objMass = objLVec.M(), objPt = objLVec.Pt(), objPx = objLVec.Px(), objPy = objLVec.Py();

  const double met = metLVec.Pt(), metphi = metLVec.Phi();



  double mt = sqrt( objMass*objMass + 2*( met*sqrt(objMass*objMass + objPt*objPt) -( met*cos(metphi)*objPx + met*sin(metphi)*objPy ) ) );



  return mt;

}
