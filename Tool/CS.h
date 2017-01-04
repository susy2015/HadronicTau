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
#include "SusyAnaTools/Tools/baselineDef.h"

static const int nSB = 84; //We use nSB serach bins depending on Nbjet, Ntop, met and MT2 value.


using namespace std;

BaselineVessel *CSBaseline = 0;
double Lumiscale = 1.0;
double EventWeight = 1.0;
class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;
  TH1D *hMET_mu;
  TH1D *hNbJets_mu;
  TH1D *hNTops_mu;
  TH1D *hMT2_mu;
  TH1D *hNJets_mu;
  TH1D *hHT_mu;
  TH1D *hYields_mu;
  TH1D *hdPhi0_mu;
  TH1D *hdPhi1_mu;
  TH1D *hdPhi2_mu;

  TH1D *hMET_el;
  TH1D *hNbJets_el;
  TH1D *hNTops_el;
  TH1D *hMT2_el;
  TH1D *hNJets_el;
  TH1D *hHT_el;
  TH1D *hYields_el;
  TH1D *hdPhi0_el;
  TH1D *hdPhi1_el;
  TH1D *hdPhi2_el;


  const TString title = "Muon CS";
  const TString title_el = "Electron CS";

};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_CS"+index+".root";
  oFile = new TFile(filename, "recreate");
 
  hMET_mu = new TH1D("hMET_mu",title+";met [GeV];Events",24,200.,800.);
  hMET_mu->Sumw2();
  hNbJets_mu = new TH1D("hNbJets_mu",title+";N_{bjets};Events",4, 1, 5);
  hNbJets_mu->Sumw2();
  hNTops_mu = new TH1D("hNTops_mu",title+";N_{tops};Events",4, 1, 5);
  hNTops_mu->Sumw2();
  hMT2_mu = new TH1D("hMT2_mu",title+";M_{T2}[GeV];Events",12,200,500);
  hMT2_mu->Sumw2();
  hYields_mu = new TH1D("hYields_mu", title+";search bin;Events",nSB,0,nSB);
  hYields_mu->Sumw2();
  hNJets_mu = new TH1D("hNJets_mu",title+";N_{jets};Events",6 ,4,10);
  hNJets_mu->Sumw2();
  hHT_mu = new TH1D("hHT_mu",title+";H_{T} [GeV];Events",20,500.,1000.);
  hHT_mu->Sumw2();  
  hdPhi0_mu = new TH1D("hdPhi0_mu", title+";dPhi0;Events", 16, 0, 3.2);
  hdPhi0_mu->Sumw2();
  hdPhi1_mu = new TH1D("hdPhi1_mu", title+";dPhi1;Events", 16, 0, 3.2);
  hdPhi1_mu->Sumw2();
  hdPhi2_mu = new TH1D("hdPhi2_mu", title+";dPhi2;Events", 16, 0, 3.2);
  hdPhi2_mu->Sumw2();
  
  hMET_el = new TH1D("hMET_el",title_el+";met [GeV];Events",24,200.,800.);
  hMET_el->Sumw2();
  hNbJets_el = new TH1D("hNbJets_el",title_el+";N_{bjets};Events",4, 1, 5);
  hNbJets_el->Sumw2();
  hNTops_el = new TH1D("hNTops_el",title_el+";N_{tops};Events",4, 1, 5);
  hNTops_el->Sumw2();
  hMT2_el = new TH1D("hMT2_el",title_el+";M_{T2}[GeV];Events",12,200,500);
  hMT2_el->Sumw2();
  hYields_el = new TH1D("hYields_el", title_el+";search bin;Events",nSB,0,nSB);
  hYields_el->Sumw2();
  hNJets_el = new TH1D("hNJets_el",title_el+";N_{jets};Events",6 ,4,10);
  hNJets_el->Sumw2();
  hHT_el = new TH1D("hHT_el",title_el+";H_{T} [GeV];Events",20,500.,1000.);
  hHT_el->Sumw2();  
  hdPhi0_el = new TH1D("hdPhi0_el", title_el+";dPhi0;Events", 16, 0, 3.2);
  hdPhi0_el->Sumw2();
  hdPhi1_el = new TH1D("hdPhi1_el", title_el+";dPhi1;Events", 16, 0, 3.2);
  hdPhi1_el->Sumw2();
  hdPhi2_el = new TH1D("hdPhi2_el", title_el+";dPhi2;Events", 16, 0, 3.2);
  hdPhi2_el->Sumw2();
  
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

void FillDouble(TH1 *hist, const double &a, const double &w){
  int nbin = hist->GetNbinsX();
  double low = hist->GetBinLowEdge(nbin);
  double high = hist->GetBinLowEdge(nbin + 1);
  double copy = a;
  if(copy >= high) copy = low;
  hist->Fill(copy, w);
}
void FillInt(TH1 *hist, const int &a, const double &w){
  int nbin = hist->GetNbinsX();
  int low = (int)hist->GetBinLowEdge(nbin);
  int high = (int)hist->GetBinLowEdge(nbin + 1);
  int copy = a;
  if(copy >= high) copy = low;
  hist->Fill(copy, w);
}
void Fill2D(TH2 *hist, const int &a, const double &b, const double &w){
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

double calcMT(const TLorentzVector &objLVec, const TLorentzVector &metLVec){

  const double objMass = objLVec.M(), objPt = objLVec.Pt(), objPx = objLVec.Px(), objPy = objLVec.Py();

  const double met = metLVec.Pt(), metphi = metLVec.Phi();

  double mt = sqrt( objMass*objMass + 2*( met*sqrt(objMass*objMass + objPt*objPt) -( met*cos(metphi)*objPx + met*sin(metphi)*objPy ) ) );

  return mt;

}
