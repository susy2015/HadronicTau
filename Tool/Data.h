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
#include "SusyAnaTools/Tools/samples.h"

static const int nSB = 69; //We use 69 serach bins depending on Nbjet, Ntop, met and MT2 value.

using namespace std;

static BaselineVessel *ExpBaselineVessel;
void passBaselineFuncExp(NTupleReader& tr)
{
  (*ExpBaselineVessel)(tr);
}
double Lumiscale = 1.0;

class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;

  TH1D *hHT;
  TH1D *hHT_wt;
  TH1D *hMET;
  TH1D *hMET_wt;
  TH1D *hNJets30;
  TH1D *hNJets30_wt;
  TH1D *hNJets50;
  TH1D *hNJets50_wt;
  TH1D *hNbJets;
  TH1D *hNbJets_wt;
  TH1D *hNTops;
  TH1D *hNTops_wt;
  TH1D *hMT2;
  TH1D *hMT2_wt;
  TH1D *hSearchBins;
  TH1D *hSearchBins_wt;

};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_Data"+index+".root";
  oFile = new TFile(filename, "recreate");
  hHT = new TH1D("hHT","HT;H_{T} [GeV];Events",20,500,1000);
  hHT->Sumw2();
  hHT_wt = new TH1D("hHT_wt","HT;H_{T} [GeV];Events",20,500,1000);
  hHT_wt->Sumw2();
  hMET = new TH1D("hMET","MET;#slash{E}_{T} [GeV];Events",24,200,800);
  hMET->Sumw2();
  hMET_wt = new TH1D("hMET_wt","MET;#slash{E}_{T} [GeV];Events",24,200,800);
  hMET_wt->Sumw2();
  hNJets30 = new TH1D("hNJets30","NJets30;N_{jets}(p_{T} > 30);Events",10,1,10);
  hNJets30->Sumw2();
  hNJets30_wt = new TH1D("hNJets30_wt","NJets30;N_{jets}(p_{T} > 30);Events",10,1,10);
  hNJets30_wt->Sumw2();
  hNJets50 = new TH1D("hNJets50","NJets50;N_{jets}(p_{T} > 50);Events",10,1,10);
  hNJets50->Sumw2();
  hNJets50_wt = new TH1D("hNJets50_wt","NJets50;N_{jets}(p_{T} > 50);Events",10,1,10);
  hNJets50_wt->Sumw2();
  hNbJets = new TH1D("hNbJets","NbJets;N_{bjets};Events",5, 0, 5);
  hNbJets->Sumw2();
  hNbJets_wt = new TH1D("hNbJets_wt","NbJets;N_{bjets};Events",5, 0, 5);
  hNbJets_wt->Sumw2();
  hNTops = new TH1D("hNTops","NTops;N_{tops};Events",5, 0, 5);
  hNTops->Sumw2();
  hNTops_wt = new TH1D("hNTops_wt","NTops;N_{tops};Events",5, 0, 5);
  hNTops_wt->Sumw2();
  hMT2 = new TH1D("hMT2","MT2;M_{T2}[GeV];Events",24,200,800);
  hMT2->Sumw2();
  hMT2_wt = new TH1D("hMT2_wt","MT2;M_{T2}[GeV];Events",24,200,800);
  hMT2_wt->Sumw2();
  hSearchBins_wt = new TH1D("hSearchBins_wt", "SearchBins;Search bin;Events",nSB,0, nSB);
  hSearchBins_wt->Sumw2();
  hSearchBins = new TH1D("hSearchBins", "SearchBins;Search bin;Events",nSB,0,nSB);
  hSearchBins->Sumw2();


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
bool passIsoTrks1(const TLorentzVector isoTrksLVec, const double isoTrksIso, const double isoTrksMtw, const int isoTrkspdgId);
