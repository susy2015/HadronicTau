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
#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "TauResponse.h"

static const int nSB = 37; //We use 37 serach bins depending on Nbjet, Ntop, met and MT2 value.
static const double mtWErr = 0.3;
static const double effErr = 0.1;
static const double BmistagErr = 0.5;
static const double pdfErr = 0.04;
static const double scaleErr = 0.04;
static const double recoErr = 0.01;
static const double isoErr = 0.01;
static const double isotrkErr = 0.1;

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
  
  TH1D *hPredYields;
  TH1D *hPredYields_wt;
  
  TH1D *hPredAccStatUp;
  TH1D *hPredAccStatLow;
  TH1D *hPredAccPDFUp;
  TH1D *hPredAccPDFLow;
  TH1D *hPredAccScaleUp;
  TH1D *hPredAccScaleLow;
  TH1D *hPredmtWMetjecUp;
  TH1D *hPredmtWMetjecLow;
  TH1D *hPredmtWMetjerUp;
  TH1D *hPredmtWMetjerLow;
  TH1D *hPredmtWStatUp;
  TH1D *hPredmtWStatLow;
  TH1D *hPredtaumuStatUp;
  TH1D *hPredtaumuStatLow;
  TH1D *hPredMuRecoStatUp;
  TH1D *hPredMuRecoStatLow;
  TH1D *hPredMuRecoTPUp;
  TH1D *hPredMuRecoTPLow;
  TH1D *hPredMuIsoStatUp;
  TH1D *hPredMuIsoStatLow;
  TH1D *hPredMuIsoTPUp;
  TH1D *hPredMuIsoTPLow;
  TH1D *hPredIsoTrkEffStatUp;
  TH1D *hPredIsoTrkEffStatLow;
  TH1D *hPredIsoTrkEffTPUp;
  TH1D *hPredIsoTrkEffTPLow;
  TH1D *hPredTemplateUp;
  TH1D *hPredTemplateLow;
  TH1D *hPredBmistagUp;
  TH1D *hPredBmistagLow;
  TH1D *hweight;
  const TString title = "Hadronic-Tau Prediction";
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_Systematic"+index+".root";
  oFile = new TFile(filename, "recreate");

  hPredYields_wt = new TH1D("hPredYields_wt", title+";search bin;Events",nSB,0,37);
  hPredYields_wt->Sumw2();
  hPredYields = new TH1D("hPredYields", title+";search bin;Events",nSB,0,37);
  hPredYields->Sumw2();

  hPredAccStatUp = new TH1D("hPredAccStatUp", title+";search bin;Events",nSB,0,37);
  hPredAccStatLow = new TH1D("hPredAccStatLow", title+";search bin;Events",nSB,0,37);
  hPredAccPDFUp = new TH1D("hPredAccPDFUp", title+";search bin;Events",nSB,0,37);
  hPredAccPDFLow = new TH1D("hPredAccPDFLow", title+";search bin;Events",nSB,0,37);
  hPredAccScaleUp = new TH1D("hPredAccScaleUp", title+";search bin;Events",nSB,0,37);
  hPredAccScaleLow = new TH1D("hPredAccScaleLow", title+";search bin;Events",nSB,0,37);
  hPredmtWStatUp = new TH1D("hPredmtWStatUp", title+";search bin;Events",nSB,0,37);
  hPredmtWStatLow = new TH1D("hPredmtWStatLow", title+";search bin;Events",nSB,0,37);
  hPredmtWMetjecUp = new TH1D("hPredmtWMetjecUp", title+";search bin;Events",nSB,0,37);
  hPredmtWMetjecLow = new TH1D("hPredmtWMetjecLow", title+";search bin;Events",nSB,0,37);
  hPredmtWMetjerUp = new TH1D("hPredmtWMetjerUp", title+";search bin;Events",nSB,0,37);
  hPredmtWMetjerLow = new TH1D("hPredmtWMetjerLow", title+";search bin;Events",nSB,0,37);
  hPredMuRecoStatUp = new TH1D("hPredMuRecoStatUp", title+";search bin;Events",nSB,0,37);
  hPredMuRecoStatLow = new TH1D("hPredMuRecoStatLow", title+";search bin;Events",nSB,0,37);
  hPredMuRecoTPUp = new TH1D("hPredMuRecoTPUp", title+";search bin;Events",nSB,0,37);
  hPredMuRecoTPLow = new TH1D("hPredMuRecoTPLow", title+";search bin;Events",nSB,0,37);
  hPredMuIsoStatUp = new TH1D("hPredMuIsoStatUp", title+";search bin;Events",nSB,0,37);
  hPredMuIsoStatLow = new TH1D("hPredMuIsoStatLow", title+";search bin;Events",nSB,0,37);
  hPredMuIsoTPUp = new TH1D("hPredMuIsoTPUp", title+";search bin;Events",nSB,0,37);
  hPredMuIsoTPLow = new TH1D("hPredMuIsoTPLow", title+";search bin;Events",nSB,0,37);
  hPredIsoTrkEffStatUp = new TH1D("hPredIsoTrkEffStatUp", title+";search bin;Events",nSB,0,37);
  hPredIsoTrkEffStatLow = new TH1D("hPredIsoTrkEffStatLow", title+";search bin;Events",nSB,0,37);
  hPredIsoTrkEffTPUp = new TH1D("hPredIsoTrkEffTPUp", title+";search bin;Events",nSB,0,37);
  hPredIsoTrkEffTPLow = new TH1D("hPredIsoTrkEffTPLow", title+";search bin;Events",nSB,0,37);
  hPredTemplateUp = new TH1D("hPredTemplateUp", title+";search bin;Events",nSB,0,37);
  hPredTemplateLow = new TH1D("hPredTemplateLow", title+";search bin;Events",nSB,0,37);
  hPredtaumuStatUp = new TH1D("hPredtaumuStatUp", title+";search bin;Events",nSB,0,37);
  hPredtaumuStatLow = new TH1D("hPredtaumuStatLow",title+";search bin;Events",nSB,0,37);
  hPredBmistagUp = new TH1D("hPredBmistagUp", title+";search bin;Events",nSB,0,37);
  hPredBmistagLow = new TH1D("hPredBmistagLow", title+";search bin;Events",nSB,0,37);

  hweight = new TH1D("hweight", "Total wight", 20,0,1);
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

double calcMT(const TLorentzVector &objLVec, const TLorentzVector &metLVec){

  const double objMass = objLVec.M(), objPt = objLVec.Pt(), objPx = objLVec.Px(), objPy = objLVec.Py();

  const double met = metLVec.Pt(), metphi = metLVec.Phi();

  double mt = sqrt( objMass*objMass + 2*( met*sqrt(objMass*objMass + objPt*objPt) -( met*cos(metphi)*objPx + met*sin(metphi)*objPy ) ) );

  return mt;

}
bool passIsoTrks1(const TLorentzVector isoTrksLVec, const double isoTrksIso, const double isoTrksMtw, const int isoTrkspdgId);
