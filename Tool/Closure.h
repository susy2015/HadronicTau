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

  TH1D *hPredHt;
  TH1D *hPredmet;
  TH1D *hPredmet_wt;
  TH1D *hPredNJets;
  TH1D *hPredNbJets;
  TH1D *hPredNbJets_wt;
  TH1D *hPredNTops;
  TH1D *hPredNTops_wt;
  TH1D *hPredMT2;
  TH1D *hPredMT2_wt;
  TH1D *hPredmTcomb;
  TH1D *hPreddPhi0;
  TH1D *hPreddPhi0_wt;
  TH1D *hPreddPhi1;
  TH1D *hPreddPhi1_wt;
  TH1D *hPreddPhi2;
  TH1D *hPreddPhi2_wt;
  
  TH1D *hPredYields;
  TH1D *hPredYields_wt;

  TH1D *hTrueHt;
  TH1D *hTruemet;
  TH1D *hTrueNJets;
  TH1D *hTrueNbJets;
  TH1D *hTrueNTops;
  TH1D *hTrueMT2;
  TH1D *hTruemTcomb;
  TH1D *hTrueYields;
  TH1D *hTruedPhi0;
  TH1D *hTruedPhi1;
  TH1D *hTruedPhi2;

  TH1D * hTruecutFlow;
  TH1D * hTrueFilter;

  TH1D *hcorrection;
  TH1D *hmtW;
  TH1D *hnomtW;
  TH1D *hmtW_wt;
  TH1D *hnomtW_wt;
  TH1D *htaumu;
  TH1D *htaumu_wt;
  TH1D *hnotaumu;
  TH1D *hnotaumu_wt;

  const TString title = "Hadronic-Tau Closure Test";
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_Closure"+index+".root";
  oFile = new TFile(filename, "recreate");
  hPredHt = new TH1D("hPredHt",title+";H_{T} [GeV];Events",25,200.,1000.);
  hPredHt->Sumw2();
  hPredmet = new TH1D("hPredmet",title+";met [GeV];Events",24,200.,800.);
  hPredmet->Sumw2();
  hPredmet_wt = new TH1D("hPredmet_wt",title+";met [GeV];Events",24,200.,800.);
  hPredmet_wt->Sumw2();
  hPredNJets = new TH1D("hPredNJets",title+";N_{jets};Events",6,4,10);
  hPredNJets->Sumw2();
  hPredNbJets = new TH1D("hPredNbJets",title+";N_{bjets};Events",4, 1, 5);
  hPredNbJets->Sumw2();
  hPredNbJets_wt = new TH1D("hPredNbJets_wt",title+";N_{bjets};Events",4, 1, 5);
  hPredNbJets_wt->Sumw2();
  hPredNTops = new TH1D("hPredNTops",title+";N_{tops};Events",4, 1, 5);
  hPredNTops->Sumw2();
  hPredNTops_wt = new TH1D("hPredNTops_wt",title+";N_{tops};Events",4, 1, 5);
  hPredNTops_wt->Sumw2();
  hPredMT2 = new TH1D("hPredMT2",title+";M_{T2}[GeV];Events",20,0,500);
  hPredMT2->Sumw2();
  hPredMT2_wt = new TH1D("hPredMT2_wt",title+";M_{T2}[GeV];Events",20,0,500);
  hPredMT2_wt->Sumw2();
  hPredmTcomb = new TH1D("hPredmTcomb",title+";M_{Tb}+0.5*M_{Tt}[GeV];Events",50,0,1000);
  hPredmTcomb->Sumw2();
  hPreddPhi0 = new TH1D("hPreddPhi0",title+";dPhi0;Events",16,0,3.2);
  hPreddPhi0->Sumw2();
  hPreddPhi0_wt = new TH1D("hPreddPhi0_wt",title+";dPhi0;Events",16,0,3.2);
  hPreddPhi0_wt->Sumw2();
  hPreddPhi1 = new TH1D("hPreddPhi1",title+";dPhi1;Events",16,0,3.2);
  hPreddPhi1->Sumw2();
  hPreddPhi1_wt = new TH1D("hPreddPhi1_wt",title+";dPhi1;Events",16,0,3.2);
  hPreddPhi1_wt->Sumw2();
  hPreddPhi2 = new TH1D("hPreddPhi2",title+";dPhi2;Events",16,0,3.2);
  hPreddPhi2->Sumw2();
  hPreddPhi2_wt = new TH1D("hPreddPhi2_wt",title+";dPhi2;Events",16,0,3.2);
  hPreddPhi2_wt->Sumw2();
  
  hPredYields_wt = new TH1D("hPredYields_wt", title+";search bin;Events",nSB,-0.5,44.5);
  hPredYields_wt->Sumw2();
  for(int bin = 1; bin <= hPredYields_wt->GetNbinsX(); ++bin) {
    TString label = "Bin ";
    label += bin;
    hPredYields_wt->GetXaxis()->SetBinLabel(bin,label);
  }
  hPredYields = new TH1D("hPredYields", title+";search bin;Events",nSB,-0.5,44.5);
  hPredYields->Sumw2();
  for(int bin = 1; bin <= hPredYields->GetNbinsX(); ++bin) {
    TString label = "Bin ";
    label += bin;
    hPredYields->GetXaxis()->SetBinLabel(bin,label);
  }

  hTrueHt = new TH1D("hTrueHt",title+";H_{T} [GeV];Events",25,200.,1000.);
  hTrueHt->Sumw2();
  hTruemet = new TH1D("hTruemet",title+";met [GeV];Events",24,200.,800.);
  hTruemet->Sumw2();
  hTrueNJets = new TH1D("hTrueNJets",title+";N_{jets};Events",6 ,4,10);
  hTrueNJets->Sumw2();
  hTrueNbJets = new TH1D("hTrueNbJets",title+";N_{bjets};Events",4, 1, 5);
  hTrueNbJets->Sumw2();
  hTrueNTops = new TH1D("hTrueNTops",title+";N_{tops};Events",4, 1, 5);
  hTrueNTops->Sumw2();
  hTrueMT2 = new TH1D("hTrueMT2",title+";M_{T2}[GeV];Events",20,0,500);
  hTrueMT2->Sumw2();
  hTruemTcomb = new TH1D("hTruemTcomb",title+";M_{Tb}+0.5*M_{Tt}[GeV];Events",50,0,1000);
  hTruemTcomb->Sumw2();
  hTrueYields = new TH1D("hTrueYields", title+";search bin;Events",nSB,-0.5,44.5);
  hTrueYields->Sumw2();
  for(int bin = 1; bin <= hTrueYields->GetNbinsX(); ++bin) {
    TString label = "Bin ";
    label += bin;
    hTrueYields->GetXaxis()->SetBinLabel(bin,label);
  }

  hTruedPhi0 = new TH1D("hTruedPhi0", title+";dPhi0;Events", 16, 0, 3.2);
  hTruedPhi0->Sumw2();
  hTruedPhi1 = new TH1D("hTruedPhi1", title+";dPhi1;Events", 16, 0, 3.2);
  hTruedPhi1->Sumw2();
  hTruedPhi2 = new TH1D("hTruedPhi2", title+";dPhi2;Events", 16, 0, 3.2);
  hTruedPhi2->Sumw2();

  hTruecutFlow = new TH1D("hTruecutFlow", "cut flow table", 20, 0, 20);
  hTruecutFlow->SetBit(TH1::kCanRebin);                                                                                            
  hTruecutFlow->Sumw2();

  hTrueFilter = new TH1D("hTrueFilter", "Event filter table", 3, 0, 3);
  hTrueFilter->SetBit(TH1::kCanRebin);
  hTrueFilter->Sumw2();

  hcorrection = new TH1D("hcorrection", "hcorrection", 20, 0.5, 1.5);
  hcorrection->Sumw2();
  hmtW = new TH1D("hmtW", "mtW correction;Search bin;Events", nTB, -0.5, 46.5);
  hmtW->Sumw2();
  hnomtW = new TH1D("hnomtW", "mtW correction;Search bin;Events", nTB, -0.5, 46.5);
  hnomtW->Sumw2();
  hmtW_wt = new TH1D("hmtW_wt", "mtW correction;Search bin;Events", nTB, -0.5, 46.5);
  hmtW_wt->Sumw2();
  hnomtW_wt = new TH1D("hnomtW_wt", "mtW correction;Search bin;Events", nTB, -0.5, 46.5);
  hnomtW_wt->Sumw2();
  htaumu = new TH1D("htaumu", "taumu contamination;Search bin;Events", nTB, -0.5, 46.5);
  htaumu->Sumw2();
  hnotaumu = new TH1D("hnotaumu", "taumu contamination;Search bin;Events", nTB, -0.5, 46.5);
  hnotaumu->Sumw2();
  htaumu_wt = new TH1D("htaumu_wt", "taumu contamination;Search bin;Events", nTB, -0.5, 46.5);
  htaumu_wt->Sumw2();
  hnotaumu_wt = new TH1D("hnotaumu_wt", "taumu contamination;Search bin;Events", nTB, -0.5, 46.5);
  hnotaumu_wt->Sumw2();
}

/*bool FillChain(TChain *chain, const TString &inputFileList)
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
}*/

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
//chain = new TChain("stopTreeMaker/AUX");
//chain->Add("root://cmsxrootd-site.fnal.gov//store/user/lpcsusyhad/Spring15_74X_Oct_2015_Ntp_v2X///WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Spring15DR74_Asympt25ns_Ntp_v2_WJetsToLNu_HT-100To200/150928_200340/0000/stopFlatNtuples_1.root"); 
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
