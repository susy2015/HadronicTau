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

static const int nSB = 59; //We use 59 serach bins depending on Nbjet, Ntop, met and MT2 value.

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

  vector<TH1*> hPredYields;
  vector<TH1*> hPredYields_wt;

  const TString title = "Hadronic-Tau Closure Test";
  TString name1 = "hPredYields_wt_";
  TString name = "hPredYields_";
  //int msMax = 1200;//stop mass highest range(T2tt)
  //int msMin = 150;//stop mass lowest range(T2tt)
  //int msMax = 1125;//stop mass highest range(T2tb)
  //int msMin = 350;//stop mass lowest range(T2tb)
  //int mlMax = 650;//lsp mass highest range
  int msMax = 2300;//gluino mass highest range(T1tttt)
  int msMin = 600;//gluino mass lowest range(T1tttt)
  int mlMax = 1600;//lsp mass highest range
  // int msMax = 2000;//gluino mass highest range(T1ttbb)
  //int msMin = 600;//gluino mass lowest range(T1ttbb)
  //int mlMax = 1450;//lsp mass highest range
  //int msMax = 1700;//gluino mass highest range(T5ttcc,T5ttttDM, T5ttttdgen)
  //int msMin = 600;//gluino mass lowest range(T5ttcc,T5ttttDM, T5ttttdgen)
  //int mlMax = 1375;//lsp mass highest range

  int find_idx(const double m1, const double m2);
  std::map<std::pair<double, double>, int> histidx;
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_SignalCont"+index+".root";
  oFile = new TFile(filename, "recreate");
  int idx = 0; 
  for(int ms = msMin; ms<=msMax; ms = ms+25){
    TString name1cp = name1;
    TString namecp = name;
    name1cp+=ms;
    name1cp+="_";
    namecp+=ms;
    namecp+="_";
    for(int ml = 0; ml<=mlMax; ml = ml+25){
      if(ml>ms)continue;
      TString name1cpy = name1cp;
      TString namecpy = namecp;
      name1cpy+=ml;
      namecpy+=ml;
      
      hPredYields_wt.push_back(new TH1D(name1cpy, title+";search bin;Events",nSB,0,nSB));
      hPredYields_wt.back()->Sumw2();
      hPredYields.push_back(new TH1D(namecpy, title+";search bin;Events",nSB,0,nSB));
      hPredYields.back()->Sumw2();
      
      std::pair<double, double> msl((double)ms, (double)ml);
      histidx[msl] = idx;
      idx++;
    }
  }
}

int BaseHistgram::find_idx(const double m1, const double m2){
  int idx = 0;
  std::pair<double, double> msl(m1, m2);
  idx = histidx[msl];
  return idx;
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

