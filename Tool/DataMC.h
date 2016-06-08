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

static const int nSB = 45; //We use 45 serach bins depending on Nbjet, Ntop, met and MT2 value.
static const int nTB = nSB + 2;// one extra bin for baseline and another bin for MT2 value less than 200 GeV 

using namespace std;

static BaselineVessel *ExpBaselineVessel;
void passBaselineFuncExp(NTupleReader& tr)
{
  (*ExpBaselineVessel)(tr);
}

AnaSamples::SampleSet        allSamples;
AnaSamples::SampleCollection allCollections(allSamples);

class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;

  TH1D *hMCMuonPt;
  TH1D *hMCMuonEta;
  TH1D *hMCJetPt;
  TH1D *hMCJetEta;
  TH1D *hMCJetPhi;
  TH1D *hMCJetMass;
  TH1D *hMC1stJetMass;
  TH1D *hMC2ndJetMass;
  TH1D *hMC3rdJetMass;
  TH1D *hMCmet;
  TH1D *hMCmetphi;
  TH1D *hMCmetphi_met;
  TH1D *hMCht;
  TH1D *hMCmt;
  TH1D *hMCnJet;
  TH1D *hDataMuonPt;
  TH1D *hDataMuonEta;
  TH1D *hDataJetPt;
  TH1D *hDataJetEta;
  TH1D *hDataJetPhi;
  TH1D *hDataJetMass;
  TH1D *hData1stJetMass;
  TH1D *hData2ndJetMass;
  TH1D *hData3rdJetMass;
  TH1D *hDatamet;
  TH1D *hDatametphi;
  TH1D *hDatametphi_met;
  TH1D *hDataht;
  TH1D *hDatamt;
  TH1D *hDatanJet;

  const TString title = "Hadronic-Tau Closure Test";
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_DataMC"+index+".root";
  oFile = new TFile(filename, "recreate");

  hMCMuonPt = new TH1D("hMCMuonPt", "Muon p_{T}", 50 , 0, 500);
  hMCMuonPt->Sumw2();
  hMCMuonEta = new TH1D("hMCMuonEta", "Muon #eta", 25 , -2.5, 2.5);
  hMCMuonEta->Sumw2();
  hMCJetPt = new TH1D("hMCJetPt", "Jet p_{T}", 50 , 0, 1000);
  hMCJetPt->Sumw2();
  hMCJetEta = new TH1D("hMCJetEta", "Jet #eta", 50 , -5, 5);
  hMCJetEta->Sumw2();
  hMCJetPhi = new TH1D("hMCJetPhi", "Jet #phi", 32 , -3.2, 3.2);
  hMCJetPhi->Sumw2();
  hMCJetMass = new TH1D("hMCJetMass", "Jet Mass", 50 , 0, 1000);
  hMCJetMass->Sumw2();
  hMC1stJetMass = new TH1D("hMC1stJetMass", "Leading Jet Mass", 50 , 0, 1000);
  hMC1stJetMass->Sumw2();
  hMC2ndJetMass = new TH1D("hMC2ndJetMass", " 2nd leading Jet Mass", 50 , 0, 1000);
  hMC2ndJetMass->Sumw2();
  hMC3rdJetMass = new TH1D("hMC3rdJetMass", "3rd leading Jet Mass", 50 , 0, 1000);
  hMC3rdJetMass->Sumw2();
  hMCmet = new TH1D("hMCmet", "met", 50 , 0, 800);
  hMCmet->Sumw2();
  hMCmetphi = new TH1D("hMCmetphi", "metphi", 32 , -3.2, 3.2);
  hMCmetphi->Sumw2();
  hMCmetphi_met = new TH1D("hMCmetphi_met", "metphi after met cut", 32 , -3.2, 3.2);
  hMCmetphi_met->Sumw2();
  hMCht = new TH1D("hMCht", "ht", 50 , 0, 1000);
  hMCht->Sumw2();
  hMCmt = new TH1D("hMCmt", "mt", 50 , 0, 500);
  hMCmt->Sumw2();
  hMCnJet = new TH1D("hMCnJet", "N_{jet}", 10 , 0, 10);
  hMCnJet->Sumw2();

  hDataMuonPt = new TH1D("hDataMuonPt", "Muon p_{T}", 50 , 0, 500);
  hDataMuonPt->Sumw2();
  hDataMuonEta = new TH1D("hDataMuonEta", "Muon #eta", 25 , -2.5, 2.5);
  hDataMuonEta->Sumw2();
  hDataJetPt = new TH1D("hDataJetPt", "Jet p_{T}", 50 , 0, 1000);
  hDataJetPt->Sumw2();
  hDataJetEta = new TH1D("hDataJetEta", "Jet #eta", 50 , -5, 5);
  hDataJetEta->Sumw2();
  hDataJetPhi = new TH1D("hDataJetPhi", "Jet #phi", 32 , -3.2, 3.2);
  hDataJetPhi->Sumw2();
  hDataJetMass = new TH1D("hDataJetMass", "Jet Mass", 50 , 0, 1000);
  hDataJetMass->Sumw2();
  hData1stJetMass = new TH1D("hData1stJetMass", "Leading Jet Mass", 50 , 0, 1000);
  hData1stJetMass->Sumw2();
  hData2ndJetMass = new TH1D("hData2ndJetMass", " 2nd leading Jet Mass", 50 , 0, 1000);
  hData2ndJetMass->Sumw2();
  hData3rdJetMass = new TH1D("hData3rdJetMass", "3rd leading Jet Mass", 50 , 0, 1000);
  hData3rdJetMass->Sumw2();
  hDatamet = new TH1D("hDatamet", "met", 50 , 0, 800);
  hDatamet->Sumw2();
  hDatametphi = new TH1D("hDatametphi", "metphi", 32, -3.2, 3.2);
  hDatametphi->Sumw2();
  hDatametphi_met = new TH1D("hDatametphi_met", "metphi after met cut", 32, -3.2, 3.2);
  hDatametphi_met->Sumw2();
  hDataht = new TH1D("hDataht", "ht", 50 , 0, 1000);
  hDataht->Sumw2();
  hDatamt = new TH1D("hDatamt", "mt", 50 , 0, 500);
  hDatamt->Sumw2();
  hDatanJet = new TH1D("hDatanJet", "N_{jet}", 10 , 0, 10);
  hDatanJet->Sumw2();
}


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
