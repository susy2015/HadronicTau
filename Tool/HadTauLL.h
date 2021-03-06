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

static const int nSB = 84; //We use nSB serach bins depending on Nbjet, Ntop, met and MT2 value.

using namespace std;
BaselineVessel *ExpBaseline = 0;
double Lumiscale = 1.0;
double EventWeight = 1.0;
class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;
  TH1D *hMET_tau;
  TH1D *hNbJets_tau;
  TH1D *hNTops_tau;
  TH1D *hMT2_tau;
  TH1D *hNJets_tau;
  TH1D *hHT_tau;
  TH1D *hYields_tau;
  TH1D *hdPhi0_tau;
  TH1D *hdPhi1_tau;
  TH1D *hdPhi2_tau;

  TH1D *hYields_tau_no_corr_SF;
  TH1D *hYields_tau_bSF;
  TH1D *hYields_tau_isrWght;

  TH1D *hYields_Veto_tau;
  TH1D *hYields_Veto_tau_SF;
  TH1D *hYields_Pass_tau;

  TH1D *hMET_LL;
  TH1D *hNbJets_LL;
  TH1D *hNTops_LL;
  TH1D *hMT2_LL;
  TH1D *hNJets_LL;
  TH1D *hHT_LL;
  TH1D *hYields_LL;
  TH1D *hdPhi0_LL;
  TH1D *hdPhi1_LL;
  TH1D *hdPhi2_LL;

  TH1D *hYields_LL_no_corr_SF;
  TH1D *hYields_LL_bSF;
  TH1D *hYields_LL_isrWght;

  TH1D *hYields_Veto_LL;
  TH1D *hYields_Veto_LL_SF;
  TH1D *hYields_Pass_LL;

  TH1D *hYields_tau_bSFup;
  TH1D *hYields_tau_bSFdown;
  TH1D *hYields_LL_bSFup;
  TH1D *hYields_LL_bSFdown;

  TH1D *hYields_tau_isrup;
  TH1D *hYields_tau_isrdown;
  TH1D *hYields_LL_isrup;
  TH1D *hYields_LL_isrdown;

  TH1D *hYields_Veto_tau_SFup;
  TH1D *hYields_Veto_tau_SFdn;
  TH1D *hYields_Veto_LL_SFup;
  TH1D *hYields_Veto_LL_SFdn;

  TH1D *hYields_tau_scaleUncup;
  TH1D *hYields_tau_scaleUncdn;
  TH1D *hYields_tau_pdfUncup;
  TH1D *hYields_tau_pdfUnccen;
  TH1D *hYields_tau_pdfUncdn;

  TH1D *hYields_tau_metMagUp;
  TH1D *hYields_tau_metMagDn;
  TH1D *hYields_tau_metPhiUp;
  TH1D *hYields_tau_metPhiDn;
  TH1D *hYields_tau_jecUp;
  TH1D *hYields_tau_jecDn;

  TH1D *hYields_LL_scaleUncup;
  TH1D *hYields_LL_scaleUncdn;
  TH1D *hYields_LL_pdfUncup;
  TH1D *hYields_LL_pdfUnccen;
  TH1D *hYields_LL_pdfUncdn;

  TH1D *hYields_LL_metMagUp;
  TH1D *hYields_LL_metMagDn;
  TH1D *hYields_LL_metPhiUp;
  TH1D *hYields_LL_metPhiDn;
  TH1D *hYields_LL_jecUp;
  TH1D *hYields_LL_jecDn;

  const TString title = "HadTau MC";
  const TString title_LL = "LL MC";
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_HadTauLL"+index+".root";
  oFile = new TFile(filename, "recreate");

  hMET_tau = new TH1D("hMET_tau",title+";met [GeV];Events",24,250.,850.);
  hMET_tau->Sumw2();
  hNbJets_tau = new TH1D("hNbJets_tau",title+";N_{bjets};Events",4, 1, 5);
  hNbJets_tau->Sumw2();
  hNTops_tau = new TH1D("hNTops_tau",title+";N_{tops};Events",4, 1, 5);
  hNTops_tau->Sumw2();
  hMT2_tau = new TH1D("hMT2_tau",title+";M_{T2}[GeV];Events",28,200,900);
  hMT2_tau->Sumw2();
  hYields_tau = new TH1D("hYields_tau", title+";search bin;Events",nSB,0,nSB);
  hYields_tau->Sumw2();
  hNJets_tau = new TH1D("hNJets_tau",title+";N_{jets};Events",6 ,4,10);
  hNJets_tau->Sumw2();
  hHT_tau = new TH1D("hHT_tau",title+";H_{T} [GeV];Events",68,300.,2000.);
  hHT_tau->Sumw2();  
  hdPhi0_tau = new TH1D("hdPhi0_tau", title+";dPhi0;Events", 16, 0, 3.2);
  hdPhi0_tau->Sumw2();
  hdPhi1_tau = new TH1D("hdPhi1_tau", title+";dPhi1;Events", 16, 0, 3.2);
  hdPhi1_tau->Sumw2();
  hdPhi2_tau = new TH1D("hdPhi2_tau", title+";dPhi2;Events", 16, 0, 3.2);
  hdPhi2_tau->Sumw2();

  hYields_tau_no_corr_SF = new TH1D("hYields_tau_no_corr_SF", title+";search bin;Events",nSB,0,nSB); hYields_tau_no_corr_SF->Sumw2();
  hYields_tau_bSF = new TH1D("hYields_tau_bSF", title+";search bin;Events",nSB,0,nSB); hYields_tau_bSF->Sumw2();
  hYields_tau_isrWght = new TH1D("hYields_tau_isrWght", title+";search bin;Events",nSB,0,nSB); hYields_tau_isrWght->Sumw2();

  hYields_Veto_tau = new TH1D("hYields_Veto_tau", title+";search bin;Events",nSB,0,nSB);
  hYields_Veto_tau->Sumw2();
  hYields_Veto_tau_SF = new TH1D("hYields_Veto_tau_SF", title+";search bin;Events",nSB,0,nSB);
  hYields_Veto_tau_SF->Sumw2();
  hYields_Pass_tau = new TH1D("hYields_Pass_tau", title+";search bin;Events",nSB,0,nSB);
  hYields_Pass_tau->Sumw2();
  
  hMET_LL = new TH1D("hMET_LL",title_LL+";met [GeV];Events",24,250.,850.);
  hMET_LL->Sumw2();
  hNbJets_LL = new TH1D("hNbJets_LL",title_LL+";N_{bjets};Events",4, 1, 5);
  hNbJets_LL->Sumw2();
  hNTops_LL = new TH1D("hNTops_LL",title_LL+";N_{tops};Events",4, 1, 5);
  hNTops_LL->Sumw2();
  hMT2_LL = new TH1D("hMT2_LL",title_LL+";M_{T2}[GeV];Events",28,200,900);
  hMT2_LL->Sumw2();
  hYields_LL = new TH1D("hYields_LL", title_LL+";search bin;Events",nSB,0,nSB);
  hYields_LL->Sumw2();
  hNJets_LL = new TH1D("hNJets_LL",title_LL+";N_{jets};Events",6 ,4,10);
  hNJets_LL->Sumw2();
  hHT_LL = new TH1D("hHT_LL",title_LL+";H_{T} [GeV];Events",68,300.,2000.);
  hHT_LL->Sumw2();  
  hdPhi0_LL = new TH1D("hdPhi0_LL", title_LL+";dPhi0;Events", 16, 0, 3.2);
  hdPhi0_LL->Sumw2();
  hdPhi1_LL = new TH1D("hdPhi1_LL", title_LL+";dPhi1;Events", 16, 0, 3.2);
  hdPhi1_LL->Sumw2();
  hdPhi2_LL = new TH1D("hdPhi2_LL", title_LL+";dPhi2;Events", 16, 0, 3.2);
  hdPhi2_LL->Sumw2();
   
  hYields_LL_no_corr_SF = new TH1D("hYields_LL_no_corr_SF", title+";search bin;Events",nSB,0,nSB); hYields_LL_no_corr_SF->Sumw2();
  hYields_LL_bSF = new TH1D("hYields_LL_bSF", title+";search bin;Events",nSB,0,nSB); hYields_LL_bSF->Sumw2();
  hYields_LL_isrWght = new TH1D("hYields_LL_isrWght", title+";search bin;Events",nSB,0,nSB); hYields_LL_isrWght->Sumw2();

  hYields_Veto_LL = new TH1D("hYields_Veto_LL", title+";search bin;Events",nSB,0,nSB);
  hYields_Veto_LL->Sumw2();
  hYields_Veto_LL_SF = new TH1D("hYields_Veto_LL_SF", title+";search bin;Events",nSB,0,nSB);
  hYields_Veto_LL_SF->Sumw2();
  hYields_Pass_LL = new TH1D("hYields_Pass_LL", title+";search bin;Events",nSB,0,nSB);
  hYields_Pass_LL->Sumw2();
  
  hYields_tau_bSFup = new TH1D("hYields_tau_bSFup", title+";search bin;Events",nSB,0,nSB);
  hYields_tau_bSFup->Sumw2();
  hYields_tau_bSFdown = new TH1D("hYields_tau_bSFdown", title+";search bin;Events",nSB,0,nSB);
  hYields_tau_bSFdown->Sumw2();
  hYields_LL_bSFup = new TH1D("hYields_LL_bSFup", title_LL+";search bin;Events",nSB,0,nSB);
  hYields_LL_bSFup->Sumw2();
  hYields_LL_bSFdown = new TH1D("hYields_LL_bSFdown", title_LL+";search bin;Events",nSB,0,nSB);
  hYields_LL_bSFdown->Sumw2();

  hYields_tau_isrup = new TH1D("hYields_tau_isrup", title+" isrup;search bin;Events",nSB,0,nSB);
  hYields_tau_isrup->Sumw2();
  hYields_tau_isrdown = new TH1D("hYields_tau_isrdown", title+" isrdown;search bin;Events",nSB,0,nSB);
  hYields_tau_isrdown->Sumw2();

  hYields_LL_isrup = new TH1D("hYields_LL_isrup", title_LL+" isrup;search bin;Events",nSB,0,nSB);
  hYields_LL_isrup->Sumw2();
  hYields_LL_isrdown = new TH1D("hYields_LL_isrdown", title_LL+" isrdown;search bin;Events",nSB,0,nSB);
  hYields_LL_isrdown->Sumw2();

  hYields_Veto_tau_SFup = new TH1D("hYields_Veto_tau_SFup", title+" SFup;search bin;Events",nSB,0,nSB); hYields_Veto_tau_SFup->Sumw2();
  hYields_Veto_tau_SFdn = new TH1D("hYields_Veto_tau_SFdn", title+" SFdn;search bin;Events",nSB,0,nSB); hYields_Veto_tau_SFdn->Sumw2();
  hYields_Veto_LL_SFup = new TH1D("hYields_Veto_LL_SFup", title_LL+" SFup;search bin;Events",nSB,0,nSB); hYields_Veto_LL_SFup->Sumw2();
  hYields_Veto_LL_SFdn = new TH1D("hYields_Veto_LL_SFdn", title_LL+" SFdn;search bin;Events",nSB,0,nSB); hYields_Veto_LL_SFdn->Sumw2();

  hYields_tau_scaleUncup = new TH1D("hYields_tau_scaleUncup", title+" scaleUncup;search bin;Events",nSB,0,nSB); hYields_tau_scaleUncup->Sumw2();
  hYields_tau_scaleUncdn = new TH1D("hYields_tau_scaleUncdn", title+" scaleUncdn;search bin;Events",nSB,0,nSB); hYields_tau_scaleUncdn->Sumw2();
  hYields_tau_pdfUncup = new TH1D("hYields_tau_pdfUncup", title+" pdfUncup;search bin;Events",nSB,0,nSB); hYields_tau_pdfUncup->Sumw2();
  hYields_tau_pdfUnccen = new TH1D("hYields_tau_pdfUnccen", title+" pdfUnccen;search bin;Events",nSB,0,nSB); hYields_tau_pdfUnccen->Sumw2();
  hYields_tau_pdfUncdn = new TH1D("hYields_tau_pdfUncdn", title+" pdfUncdn;search bin;Events",nSB,0,nSB); hYields_tau_pdfUncdn->Sumw2();

  hYields_tau_metMagUp = new TH1D("hYields_tau_metMagUp", title+" metMagUp;search bin;Events",nSB,0,nSB); hYields_tau_metMagUp->Sumw2();
  hYields_tau_metMagDn = new TH1D("hYields_tau_metMagDn", title+" metMagDn;search bin;Events",nSB,0,nSB); hYields_tau_metMagDn->Sumw2();
  hYields_tau_metPhiUp = new TH1D("hYields_tau_metPhiUp", title+" metPhiUp;search bin;Events",nSB,0,nSB); hYields_tau_metPhiUp->Sumw2();
  hYields_tau_metPhiDn = new TH1D("hYields_tau_metPhiDn", title+" metPhiDn;search bin;Events",nSB,0,nSB); hYields_tau_metPhiDn->Sumw2();
  hYields_tau_jecUp = new TH1D("hYields_tau_jecUp", title+" jecUp;search bin;Events",nSB,0,nSB); hYields_tau_jecUp->Sumw2();
  hYields_tau_jecDn = new TH1D("hYields_tau_jecDn", title+" jecDn;search bin;Events",nSB,0,nSB); hYields_tau_jecDn->Sumw2();

  hYields_LL_scaleUncup = new TH1D("hYields_LL_scaleUncup", title_LL+" scaleUncup;search bin;Events",nSB,0,nSB); hYields_LL_scaleUncup->Sumw2();
  hYields_LL_scaleUncdn = new TH1D("hYields_LL_scaleUncdn", title_LL+" scaleUncdn;search bin;Events",nSB,0,nSB); hYields_LL_scaleUncdn->Sumw2();
  hYields_LL_pdfUncup = new TH1D("hYields_LL_pdfUncup", title_LL+" pdfUncup;search bin;Events",nSB,0,nSB); hYields_LL_pdfUncup->Sumw2();
  hYields_LL_pdfUnccen = new TH1D("hYields_LL_pdfUnccen", title_LL+" pdfUnccen;search bin;Events",nSB,0,nSB); hYields_LL_pdfUnccen->Sumw2();
  hYields_LL_pdfUncdn = new TH1D("hYields_LL_pdfUncdn", title_LL+" pdfUncdn;search bin;Events",nSB,0,nSB); hYields_LL_pdfUncdn->Sumw2();

  hYields_LL_metMagUp = new TH1D("hYields_LL_metMagUp", title_LL+" metMagUp;search bin;Events",nSB,0,nSB); hYields_LL_metMagUp->Sumw2();
  hYields_LL_metMagDn = new TH1D("hYields_LL_metMagDn", title_LL+" metMagDn;search bin;Events",nSB,0,nSB); hYields_LL_metMagDn->Sumw2();
  hYields_LL_metPhiUp = new TH1D("hYields_LL_metPhiUp", title_LL+" metPhiUp;search bin;Events",nSB,0,nSB); hYields_LL_metPhiUp->Sumw2();
  hYields_LL_metPhiDn = new TH1D("hYields_LL_metPhiDn", title_LL+" metPhiDn;search bin;Events",nSB,0,nSB); hYields_LL_metPhiDn->Sumw2();
  hYields_LL_jecUp = new TH1D("hYields_LL_jecUp", title_LL+" jecUp;search bin;Events",nSB,0,nSB); hYields_LL_jecUp->Sumw2();
  hYields_LL_jecDn = new TH1D("hYields_LL_jecDn", title_LL+" jecDn;search bin;Events",nSB,0,nSB); hYields_LL_jecDn->Sumw2();
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

