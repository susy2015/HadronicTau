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

  TH1D *hYields_mu_bSFup;
  TH1D *hYields_mu_bSFdown;
  TH1D *hYields_el_bSFup;
  TH1D *hYields_el_bSFdown;

  TH1D *hYields_mu_isrup;
  TH1D *hYields_mu_isrdown;
  TH1D *hYields_el_isrup;
  TH1D *hYields_el_isrdown;

  TH1D *hYields_mu_SFup;
  TH1D *hYields_mu_SFdn;
  TH1D *hYields_mu_scaleUncup;
  TH1D *hYields_mu_scaleUncdn;
  TH1D *hYields_mu_pdfUncup;
  TH1D *hYields_mu_pdfUnccen;
  TH1D *hYields_mu_pdfUncdn;

  TH1D *hYields_mu_metMagUp;
  TH1D *hYields_mu_metMagDn;
  TH1D *hYields_mu_metPhiUp;
  TH1D *hYields_mu_metPhiDn;
  TH1D *hYields_mu_jecUp;
  TH1D *hYields_mu_jecDn;

  TH1D *hYields_el_SFup;
  TH1D *hYields_el_SFdn;
  TH1D *hYields_el_scaleUncup;
  TH1D *hYields_el_scaleUncdn;
  TH1D *hYields_el_pdfUncup;
  TH1D *hYields_el_pdfUnccen;
  TH1D *hYields_el_pdfUncdn;

  TH1D *hYields_el_metMagUp;
  TH1D *hYields_el_metMagDn;
  TH1D *hYields_el_metPhiUp;
  TH1D *hYields_el_metPhiDn;
  TH1D *hYields_el_jecUp;
  TH1D *hYields_el_jecDn;

  const TString title = "Muon CS";
  const TString title_el = "Electron CS";

};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_CS"+index+".root";
  oFile = new TFile(filename, "recreate");
 
  hMET_mu = new TH1D("hMET_mu",title+";met [GeV];Events",24,250.,850.);
  hMET_mu->Sumw2();
  hNbJets_mu = new TH1D("hNbJets_mu",title+";N_{bjets};Events",4, 1, 5);
  hNbJets_mu->Sumw2();
  hNTops_mu = new TH1D("hNTops_mu",title+";N_{tops};Events",4, 1, 5);
  hNTops_mu->Sumw2();
  hMT2_mu = new TH1D("hMT2_mu",title+";M_{T2}[GeV];Events",28,200,900);
  hMT2_mu->Sumw2();
  hYields_mu = new TH1D("hYields_mu", title+";search bin;Events",nSB,0,nSB);
  hYields_mu->Sumw2();
  hNJets_mu = new TH1D("hNJets_mu",title+";N_{jets};Events",6 ,4,10);
  hNJets_mu->Sumw2();
  hHT_mu = new TH1D("hHT_mu",title+";H_{T} [GeV];Events",68,300.,2000.);
  hHT_mu->Sumw2();  
  hdPhi0_mu = new TH1D("hdPhi0_mu", title+";dPhi0;Events", 16, 0, 3.2);
  hdPhi0_mu->Sumw2();
  hdPhi1_mu = new TH1D("hdPhi1_mu", title+";dPhi1;Events", 16, 0, 3.2);
  hdPhi1_mu->Sumw2();
  hdPhi2_mu = new TH1D("hdPhi2_mu", title+";dPhi2;Events", 16, 0, 3.2);
  hdPhi2_mu->Sumw2();
  
  hMET_el = new TH1D("hMET_el",title_el+";met [GeV];Events",24,250.,850.);
  hMET_el->Sumw2();
  hNbJets_el = new TH1D("hNbJets_el",title_el+";N_{bjets};Events",4, 1, 5);
  hNbJets_el->Sumw2();
  hNTops_el = new TH1D("hNTops_el",title_el+";N_{tops};Events",4, 1, 5);
  hNTops_el->Sumw2();
  hMT2_el = new TH1D("hMT2_el",title_el+";M_{T2}[GeV];Events",28,200,900);
  hMT2_el->Sumw2();
  hYields_el = new TH1D("hYields_el", title_el+";search bin;Events",nSB,0,nSB);
  hYields_el->Sumw2();
  hNJets_el = new TH1D("hNJets_el",title_el+";N_{jets};Events",6 ,4,10);
  hNJets_el->Sumw2();
  hHT_el = new TH1D("hHT_el",title_el+";H_{T} [GeV];Events",68,300.,2000.);
  hHT_el->Sumw2();  
  hdPhi0_el = new TH1D("hdPhi0_el", title_el+";dPhi0;Events", 16, 0, 3.2);
  hdPhi0_el->Sumw2();
  hdPhi1_el = new TH1D("hdPhi1_el", title_el+";dPhi1;Events", 16, 0, 3.2);
  hdPhi1_el->Sumw2();
  hdPhi2_el = new TH1D("hdPhi2_el", title_el+";dPhi2;Events", 16, 0, 3.2);
  hdPhi2_el->Sumw2();
  
  hYields_mu_bSFup = new TH1D("hYields_mu_bSFup", title+";search bin;Events",nSB,0,nSB);
  hYields_mu_bSFup->Sumw2();
  hYields_mu_bSFdown = new TH1D("hYields_mu_bSFdown", title+";search bin;Events",nSB,0,nSB);
  hYields_mu_bSFdown->Sumw2();
  hYields_el_bSFup = new TH1D("hYields_el_bSFup", title_el+";search bin;Events",nSB,0,nSB);
  hYields_el_bSFup->Sumw2();
  hYields_el_bSFdown = new TH1D("hYields_el_bSFdown", title_el+";search bin;Events",nSB,0,nSB);
  hYields_el_bSFdown->Sumw2();

  hYields_mu_isrup = new TH1D("hYields_mu_isrup", title+" isrup;search bin;Events",nSB,0,nSB);
  hYields_mu_isrup->Sumw2();
  hYields_mu_isrdown = new TH1D("hYields_mu_isrdown", title+" isrdown;search bin;Events",nSB,0,nSB);
  hYields_mu_isrdown->Sumw2();

  hYields_el_isrup = new TH1D("hYields_el_isrup", title_el+" isrup;search bin;Events",nSB,0,nSB);
  hYields_el_isrup->Sumw2();
  hYields_el_isrdown = new TH1D("hYields_el_isrdown", title_el+" isrdown;search bin;Events",nSB,0,nSB);
  hYields_el_isrdown->Sumw2();

  hYields_mu_SFup = new TH1D("hYields_mu_SFup", title+"  SFup;search bin; Events",nSB,0,nSB); hYields_mu_SFup->Sumw2();
  hYields_mu_SFdn = new TH1D("hYields_mu_SFdn", title+"  SFdn;search bin; Events",nSB,0,nSB); hYields_mu_SFdn->Sumw2();
  hYields_mu_scaleUncup = new TH1D("hYields_mu_scaleUncup", title+"  scaleUncup;search bin; Events",nSB,0,nSB); hYields_mu_scaleUncup->Sumw2();
  hYields_mu_scaleUncdn = new TH1D("hYields_mu_scaleUncdn", title+"  scaleUncdn;search bin; Events",nSB,0,nSB); hYields_mu_scaleUncdn->Sumw2();
  hYields_mu_pdfUncup = new TH1D("hYields_mu_pdfUncup", title+"  pdfUncup;search bin; Events",nSB,0,nSB); hYields_mu_pdfUncup->Sumw2();
  hYields_mu_pdfUnccen = new TH1D("hYields_mu_pdfUnccen", title+"  pdfUnccen;search bin; Events",nSB,0,nSB); hYields_mu_pdfUnccen->Sumw2();
  hYields_mu_pdfUncdn = new TH1D("hYields_mu_pdfUncdn", title+"  pdfUncdn;search bin; Events",nSB,0,nSB); hYields_mu_pdfUncdn->Sumw2();

  hYields_mu_metMagUp = new TH1D("hYields_mu_metMagUp", title+" metMagUp;search bin; Events",nSB,0,nSB); hYields_mu_metMagUp->Sumw2();
  hYields_mu_metMagDn = new TH1D("hYields_mu_metMagDn", title+" metMagDn;search bin; Events",nSB,0,nSB); hYields_mu_metMagDn->Sumw2();
  hYields_mu_metPhiUp = new TH1D("hYields_mu_metPhiUp", title+" metPhiUp;search bin; Events",nSB,0,nSB); hYields_mu_metPhiUp->Sumw2();
  hYields_mu_metPhiDn = new TH1D("hYields_mu_metPhiDn", title+" metPhiDn;search bin; Events",nSB,0,nSB); hYields_mu_metPhiDn->Sumw2();
  hYields_mu_jecUp = new TH1D("hYields_mu_jecUp", title+" jecUp;search bin; Events",nSB,0,nSB); hYields_mu_jecUp->Sumw2();
  hYields_mu_jecDn = new TH1D("hYields_mu_jecDn", title+" jecDn;search bin; Events",nSB,0,nSB); hYields_mu_jecDn->Sumw2();

  hYields_el_SFup = new TH1D("hYields_el_SFup", title_el+"  SFup;search bin; Events",nSB,0,nSB); hYields_el_SFup->Sumw2();
  hYields_el_SFdn = new TH1D("hYields_el_SFdn", title_el+"  SFdn;search bin; Events",nSB,0,nSB); hYields_el_SFdn->Sumw2();
  hYields_el_scaleUncup = new TH1D("hYields_el_scaleUncup", title_el+"  scaleUncup;search bin; Events",nSB,0,nSB); hYields_el_scaleUncup->Sumw2();
  hYields_el_scaleUncdn = new TH1D("hYields_el_scaleUncdn", title_el+"  scaleUncdn;search bin; Events",nSB,0,nSB); hYields_el_scaleUncdn->Sumw2();
  hYields_el_pdfUncup = new TH1D("hYields_el_pdfUncup", title_el+"  pdfUncup;search bin; Events",nSB,0,nSB); hYields_el_pdfUncup->Sumw2();
  hYields_el_pdfUnccen = new TH1D("hYields_el_pdfUnccen", title_el+"  pdfUnccen;search bin; Events",nSB,0,nSB); hYields_el_pdfUnccen->Sumw2();
  hYields_el_pdfUncdn = new TH1D("hYields_el_pdfUncdn", title_el+"  pdfUncdn;search bin; Events",nSB,0,nSB); hYields_el_pdfUncdn->Sumw2();

  hYields_el_metMagUp = new TH1D("hYields_el_metMagUp", title_el+" metMagUp;search bin; Events",nSB,0,nSB); hYields_el_metMagUp->Sumw2();
  hYields_el_metMagDn = new TH1D("hYields_el_metMagDn", title_el+" metMagDn;search bin; Events",nSB,0,nSB); hYields_el_metMagDn->Sumw2();
  hYields_el_metPhiUp = new TH1D("hYields_el_metPhiUp", title_el+" metPhiUp;search bin; Events",nSB,0,nSB); hYields_el_metPhiUp->Sumw2();
  hYields_el_metPhiDn = new TH1D("hYields_el_metPhiDn", title_el+" metPhiDn;search bin; Events",nSB,0,nSB); hYields_el_metPhiDn->Sumw2();
  hYields_el_jecUp = new TH1D("hYields_el_jecUp", title_el+" jecUp;search bin; Events",nSB,0,nSB); hYields_el_jecUp->Sumw2();
  hYields_el_jecDn = new TH1D("hYields_el_jecDn", title_el+" jecDn;search bin; Events",nSB,0,nSB); hYields_el_jecDn->Sumw2();

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
