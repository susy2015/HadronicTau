#ifndef SYSTEMATICS_H
#define SYSTEMATICS_H

#include <iostream>
#include <vector>
#include "TH1.h"
#include "TFile.h"
#include "CommonShared.h"

class Systematics
{
  public:
  //Constructor
  // expStr = "tau", "LL"
  // csStr = "mu", "el"
    Systematics(const TString expStr="tau", const TString csStr="mu", const TString filename_CS = "Mix_CS.root", const TString filename_HadTauLL = "Mix_HadTauLL.root", const bool do_mergeBins = false);

  //Member Data
    std::vector<double> ISRsys_up;
    std::vector<double> ISRsys_dn;

    std::vector<double> bSFsys_up;
    std::vector<double> bSFsys_dn;

    std::vector<double> scaleUncsys_up;
    std::vector<double> scaleUncsys_dn;

    std::vector<double> pdfUncsys_up;
    std::vector<double> pdfUncsys_cen;
    std::vector<double> pdfUncsys_dn;

    std::vector<double> metMagsys_up;
    std::vector<double> metMagsys_dn;

    std::vector<double> metPhisys_up;
    std::vector<double> metPhisys_dn;

    std::vector<double> jecsys_up;
    std::vector<double> jecsys_dn;

    TFile *file_CS;
    TFile *file_HadTauLL;

    const bool do_mergeBins_;
};

Systematics::Systematics(const TString expStr, const TString csStr, const TString filename_CS, const TString filename_HadTauLL, const bool do_mergeBins) : 
   do_mergeBins_(do_mergeBins)
{
  std::cout<<"expStr : "<<expStr<<"  csStr : "<<csStr<<"  filename_CS : "<<filename_CS<<"  filename_HadTauLL : "<<filename_HadTauLL<<"  do_mergeBins_ : "<<do_mergeBins_<<std::endl;
// adjust merge bins before calculating ANY ratio

  file_CS = new TFile(filename_CS);
  file_HadTauLL = new TFile(filename_HadTauLL);

  TH1D *hCS = (TH1D*)file_CS->Get("hYields_"+csStr);
  TH1D *hExp = (expStr=="tau" || expStr=="LL") ? (TH1D*)file_HadTauLL->Get("hYields_"+expStr) : expStr=="comb"? (TH1D*)file_HadTauLL->Get("hYields_tau") : nullptr;
  if( hExp == nullptr ){ std::cout<<"illigal expStr : "<<expStr<<std::endl; return; }
  if( expStr == "comb" ) hExp->Add((TH1D*)file_HadTauLL->Get("hYields_LL"));
  if( do_mergeBins_ ) predSpec::adjustBins_merge(hCS);
  if( do_mergeBins_ ) predSpec::adjustBins_merge(hExp);
  TH1D *hRatio_Exp_CS = (TH1D*)hExp->Clone("Ratio");
  hRatio_Exp_CS->Divide(hCS);

  TH1D *hCS_pdfUncsys_cen = (TH1D*)file_CS->Get("hYields_"+csStr+"_pdfUnccen");
  TH1D *hExp_pdfUncsys_cen = (expStr=="tau" || expStr=="LL") ? (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_pdfUnccen") : expStr=="comb"? (TH1D*)file_HadTauLL->Get("hYields_tau_pdfUnccen") : nullptr;
  if( expStr == "comb" ) hExp_pdfUncsys_cen->Add((TH1D*)file_HadTauLL->Get("hYields_LL_pdfUnccen"));
  if( do_mergeBins_ ) predSpec::adjustBins_merge(hCS_pdfUncsys_cen);
  if( do_mergeBins_ ) predSpec::adjustBins_merge(hExp_pdfUncsys_cen);
  TH1D *hRatio_Exp_CS_pdfUncsys_cen = (TH1D*)hExp_pdfUncsys_cen->Clone("Ratio");
  hRatio_Exp_CS_pdfUncsys_cen->Divide(hCS_pdfUncsys_cen);

  const int kNDists = 2; // up and dn(down) 
  TH1D* hCS_ISR[kNDists];
  TH1D* hExp_ISR[kNDists];
  TH1D* hRatio_Exp_CS_ISR[kNDists];

  TH1D* hCS_bSF[kNDists];
  TH1D* hExp_bSF[kNDists];
  TH1D* hRatio_Exp_CS_bSF[kNDists];

  TH1D* hCS_scaleUncsys[kNDists];
  TH1D* hExp_scaleUncsys[kNDists];
  TH1D* hRatio_Exp_CS_scaleUncsys[kNDists];

  TH1D* hCS_pdfUncsys[kNDists];
  TH1D* hExp_pdfUncsys[kNDists];
  TH1D* hRatio_Exp_CS_pdfUncsys[kNDists];

  TH1D* hCS_metMagsys[kNDists];
  TH1D* hExp_metMagsys[kNDists];
  TH1D* hRatio_Exp_CS_metMagsys[kNDists];

  TH1D* hCS_metPhisys[kNDists];
  TH1D* hExp_metPhisys[kNDists];
  TH1D* hRatio_Exp_CS_metPhisys[kNDists];

  TH1D* hCS_jecsys[kNDists];
  TH1D* hExp_jecsys[kNDists];
  TH1D* hRatio_Exp_CS_jecsys[kNDists];

  for(unsigned int i = 0; i < kNDists; ++i)
  {
    TString name = "";
    if(      i == 0 ) name = "up";
    else if( i == 1 ) name = "down";
    // Get histograms from file                                                                                                                
    hCS_ISR[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_isr"+name);
    hExp_ISR[i] = (expStr=="tau" || expStr=="LL") ? (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_isr"+name) : expStr=="comb"? (TH1D*)file_HadTauLL->Get("hYields_tau_isr"+name) : nullptr;
    if( expStr == "comb" ) hExp_ISR[i]->Add((TH1D*)file_HadTauLL->Get("hYields_LL_isr"+name));
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hCS_ISR[i]);
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hExp_ISR[i]);
    //Ratio                                                                                                                                   
    hRatio_Exp_CS_ISR[i] = (TH1D*)hExp_ISR[i]->Clone("Ratio");
    hRatio_Exp_CS_ISR[i]->Divide(hCS_ISR[i]);

// bSF
    hCS_bSF[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_bSF"+name);
    hExp_bSF[i] = (expStr=="tau" || expStr=="LL") ? (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_bSF"+name) : expStr=="comb"? (TH1D*)file_HadTauLL->Get("hYields_tau_bSF"+name) : nullptr;
    if( expStr == "comb" ) hExp_bSF[i]->Add((TH1D*)file_HadTauLL->Get("hYields_LL_bSF"+name));
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hCS_bSF[i]);
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hExp_bSF[i]);
    //Ratio                                                                                                                                   
    hRatio_Exp_CS_bSF[i] = (TH1D*)hExp_bSF[i]->Clone("Ratio");
    hRatio_Exp_CS_bSF[i]->Divide(hCS_bSF[i]);

    if(      i == 0 ) name = "up";
    else if( i == 1 ) name = "dn";

    hCS_scaleUncsys[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_scaleUnc"+name);
    hExp_scaleUncsys[i] = (expStr=="tau" || expStr=="LL") ? (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_scaleUnc"+name) : expStr=="comb"? (TH1D*)file_HadTauLL->Get("hYields_tau_scaleUnc"+name) : nullptr;
    if( expStr == "comb" ) hExp_scaleUncsys[i]->Add((TH1D*)file_HadTauLL->Get("hYields_LL_scaleUnc"+name));
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hCS_scaleUncsys[i]);
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hExp_scaleUncsys[i]);
    hRatio_Exp_CS_scaleUncsys[i] = (TH1D*)hExp_scaleUncsys[i]->Clone("Ratio");
    hRatio_Exp_CS_scaleUncsys[i]->Divide(hCS_scaleUncsys[i]);

    hCS_pdfUncsys[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_pdfUnc"+name);
    hExp_pdfUncsys[i] = (expStr=="tau" || expStr=="LL") ? (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_pdfUnc"+name) : expStr=="comb"? (TH1D*)file_HadTauLL->Get("hYields_tau_pdfUnc"+name) : nullptr;
    if( expStr == "comb" ) hExp_pdfUncsys[i]->Add((TH1D*)file_HadTauLL->Get("hYields_LL_pdfUnc"+name));
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hCS_pdfUncsys[i]);
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hExp_pdfUncsys[i]);
//    hRatio_Exp_CS_pdfUncsys[i] = (TH1D*)hExp_pdfUncsys[i]->Clone("Ratio");
//    hRatio_Exp_CS_pdfUncsys[i]->Divide(hCS_pdfUncsys[i]);

    if(      i == 0 ) name = "Up";
    else if( i == 1 ) name = "Dn";

    hCS_metMagsys[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_metMag"+name);
    hExp_metMagsys[i] = (expStr=="tau" || expStr=="LL") ? (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_metMag"+name) : expStr=="comb"? (TH1D*)file_HadTauLL->Get("hYields_tau_metMag"+name) : nullptr;
    if( expStr == "comb" ) hExp_metMagsys[i]->Add((TH1D*)file_HadTauLL->Get("hYields_LL_metMag"+name));
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hCS_metMagsys[i]);
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hExp_metMagsys[i]);
    hRatio_Exp_CS_metMagsys[i] = (TH1D*)hExp_metMagsys[i]->Clone("Ratio");
    hRatio_Exp_CS_metMagsys[i]->Divide(hCS_metMagsys[i]);

    hCS_metPhisys[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_metPhi"+name);
    hExp_metPhisys[i] = (expStr=="tau" || expStr=="LL") ? (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_metPhi"+name) : expStr=="comb"? (TH1D*)file_HadTauLL->Get("hYields_tau_metPhi"+name) : nullptr;
    if( expStr == "comb" ) hExp_metPhisys[i]->Add((TH1D*)file_HadTauLL->Get("hYields_LL_metPhi"+name));
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hCS_metPhisys[i]);
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hExp_metPhisys[i]);
    hRatio_Exp_CS_metPhisys[i] = (TH1D*)hExp_metPhisys[i]->Clone("Ratio");
    hRatio_Exp_CS_metPhisys[i]->Divide(hCS_metPhisys[i]);

    hCS_jecsys[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_jec"+name);
    hExp_jecsys[i] = (expStr=="tau" || expStr=="LL") ? (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_jec"+name) : expStr=="comb"? (TH1D*)file_HadTauLL->Get("hYields_tau_jec"+name) : nullptr;
    if( expStr == "comb" ) hExp_jecsys[i]->Add((TH1D*)file_HadTauLL->Get("hYields_LL_jec"+name));
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hCS_jecsys[i]);
    if( do_mergeBins_ ) predSpec::adjustBins_merge(hExp_jecsys[i]);
    hRatio_Exp_CS_jecsys[i] = (TH1D*)hExp_jecsys[i]->Clone("Ratio");
    hRatio_Exp_CS_jecsys[i]->Divide(hCS_jecsys[i]);
  }
  double sum_CS_pdfUncsys_up = 0, sum_CS_pdfUncsys_cen = 0, sum_CS_pdfUncsys_dn = 0;
  for(unsigned int i=1; i<hCS_pdfUncsys_cen->GetNbinsX(); i++)
  {
     sum_CS_pdfUncsys_cen += hCS_pdfUncsys_cen->GetBinContent(i);
     sum_CS_pdfUncsys_up += hCS_pdfUncsys[0]->GetBinContent(i);
     sum_CS_pdfUncsys_dn += hCS_pdfUncsys[1]->GetBinContent(i);
  }
  double sum_Exp_pdfUncsys_up = 0, sum_Exp_pdfUncsys_cen = 0, sum_Exp_pdfUncsys_dn = 0;
  for(unsigned int i=1; i<hExp_pdfUncsys_cen->GetNbinsX(); i++)
  {
     sum_Exp_pdfUncsys_cen += hExp_pdfUncsys_cen->GetBinContent(i);
     sum_Exp_pdfUncsys_up += hExp_pdfUncsys[0]->GetBinContent(i);
     sum_Exp_pdfUncsys_dn += hExp_pdfUncsys[1]->GetBinContent(i);
  }
  hCS_pdfUncsys[0]->Scale(sum_CS_pdfUncsys_cen/sum_CS_pdfUncsys_up);
  hCS_pdfUncsys[1]->Scale(sum_CS_pdfUncsys_cen/sum_CS_pdfUncsys_dn);

  hExp_pdfUncsys[0]->Scale(sum_Exp_pdfUncsys_cen/sum_Exp_pdfUncsys_up);
  hExp_pdfUncsys[1]->Scale(sum_Exp_pdfUncsys_cen/sum_Exp_pdfUncsys_dn);

  for(unsigned int i = 0; i < kNDists; ++i)
  {
     hRatio_Exp_CS_pdfUncsys[i] = (TH1D*)hExp_pdfUncsys[i]->Clone("Ratio");
     hRatio_Exp_CS_pdfUncsys[i]->Divide(hCS_pdfUncsys[i]);
  }
/*
  hRatio_Exp_CS_pdfUncsys_cen->Scale(1./sum_pdfUncsys_cen);
  hRatio_Exp_CS_pdfUncsys[0]->Scale(1./sum_pdfUncsys_up);
  hRatio_Exp_CS_pdfUncsys[1]->Scale(1./sum_pdfUncsys_dn);
*/
  for(unsigned i=1; i<= hRatio_Exp_CS->GetNbinsX();i++)
  {
    ISRsys_up.push_back(std::fabs(hRatio_Exp_CS->GetBinContent(i)-hRatio_Exp_CS_ISR[0]->GetBinContent(i)));
    ISRsys_dn.push_back(std::fabs(hRatio_Exp_CS->GetBinContent(i)-hRatio_Exp_CS_ISR[1]->GetBinContent(i)));

    bSFsys_up.push_back(std::fabs(hRatio_Exp_CS->GetBinContent(i)-hRatio_Exp_CS_bSF[0]->GetBinContent(i)));
    bSFsys_dn.push_back(std::fabs(hRatio_Exp_CS->GetBinContent(i)-hRatio_Exp_CS_bSF[1]->GetBinContent(i)));

    scaleUncsys_up.push_back(std::fabs(hRatio_Exp_CS->GetBinContent(i)-hRatio_Exp_CS_scaleUncsys[0]->GetBinContent(i)));
    scaleUncsys_dn.push_back(std::fabs(hRatio_Exp_CS->GetBinContent(i)-hRatio_Exp_CS_scaleUncsys[1]->GetBinContent(i)));

    pdfUncsys_up.push_back(std::fabs(hRatio_Exp_CS_pdfUncsys_cen->GetBinContent(i)-hRatio_Exp_CS_pdfUncsys[0]->GetBinContent(i)));
    pdfUncsys_dn.push_back(std::fabs(hRatio_Exp_CS_pdfUncsys_cen->GetBinContent(i)-hRatio_Exp_CS_pdfUncsys[1]->GetBinContent(i)));
    pdfUncsys_cen.push_back(hRatio_Exp_CS_pdfUncsys_cen->GetBinContent(i));

    metMagsys_up.push_back(std::fabs(hRatio_Exp_CS->GetBinContent(i)-hRatio_Exp_CS_metMagsys[0]->GetBinContent(i)));
    metMagsys_dn.push_back(std::fabs(hRatio_Exp_CS->GetBinContent(i)-hRatio_Exp_CS_metMagsys[1]->GetBinContent(i)));

    metPhisys_up.push_back(std::fabs(hRatio_Exp_CS->GetBinContent(i)-hRatio_Exp_CS_metPhisys[0]->GetBinContent(i)));
    metPhisys_dn.push_back(std::fabs(hRatio_Exp_CS->GetBinContent(i)-hRatio_Exp_CS_metPhisys[1]->GetBinContent(i)));

    jecsys_up.push_back(std::fabs(hRatio_Exp_CS->GetBinContent(i)-hRatio_Exp_CS_jecsys[0]->GetBinContent(i)));
    jecsys_dn.push_back(std::fabs(hRatio_Exp_CS->GetBinContent(i)-hRatio_Exp_CS_jecsys[1]->GetBinContent(i)));
  }

}
#endif
