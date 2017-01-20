#ifndef SYSTEMATICS_H
#define SYSTEMATICS_H

#include <iostream>
#include <vector>
#include "TH1.h"
#include "TFile.h"

class Systematics
{
  public:
  //Constructor
  // expStr = "tau", "LL"
  // csStr = "mu", "el"
    Systematics(const TString expStr="tau", const TString csStr="mu", const TString filename_CS = "Mix_CS.root", const TString filename_HadTauLL = "Mix_HadTauLL.root");

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
};

Systematics::Systematics(const TString expStr, const TString csStr, const TString filename_CS, const TString filename_HadTauLL)
{
  std::cout<<"expStr : "<<expStr<<"  csStr : "<<csStr<<"  filename_CS : "<<filename_CS<<"  filename_HadTauLL : "<<filename_HadTauLL<<std::endl;

  file_CS = new TFile(filename_CS);
  file_HadTauLL = new TFile(filename_HadTauLL);

  TH1D *hCS = (TH1D*)file_CS->Get("hYields_"+csStr);
  TH1D *hExp = (TH1D*)file_HadTauLL->Get("hYields_"+expStr);

  TH1D *hRatio_Exp_CS = (TH1D*)hExp->Clone("Ratio");
  hRatio_Exp_CS->Divide(hCS);

  TH1D *hCS_pdfUncsys_cen = (TH1D*)file_CS->Get("hYields_"+csStr+"_pdfUnccen");
  TH1D *hExp_pdfUncsys_cen = (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_pdfUnccen");
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
    hExp_ISR[i] = (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_isr"+name);
    //Ratio                                                                                                                                   
    hRatio_Exp_CS_ISR[i] = (TH1D*)hExp_ISR[i]->Clone("Ratio");
    hRatio_Exp_CS_ISR[i]->Divide(hCS_ISR[i]);

// bSF
    hCS_bSF[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_bSF"+name);
    hExp_bSF[i] = (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_bSF"+name);
    //Ratio                                                                                                                                   
    hRatio_Exp_CS_bSF[i] = (TH1D*)hExp_bSF[i]->Clone("Ratio");
    hRatio_Exp_CS_bSF[i]->Divide(hCS_bSF[i]);

    if(      i == 0 ) name = "up";
    else if( i == 1 ) name = "dn";

    hCS_scaleUncsys[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_scaleUnc"+name);
    hExp_scaleUncsys[i] = (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_scaleUnc"+name);
    hRatio_Exp_CS_scaleUncsys[i] = (TH1D*)hExp_scaleUncsys[i]->Clone("Ratio");
    hRatio_Exp_CS_scaleUncsys[i]->Divide(hCS_scaleUncsys[i]);

    hCS_pdfUncsys[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_pdfUnc"+name);
    hExp_pdfUncsys[i] = (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_pdfUnc"+name);
    hRatio_Exp_CS_pdfUncsys[i] = (TH1D*)hExp_pdfUncsys[i]->Clone("Ratio");
    hRatio_Exp_CS_pdfUncsys[i]->Divide(hCS_pdfUncsys[i]);

    if(      i == 0 ) name = "Up";
    else if( i == 1 ) name = "Dn";

    hCS_metMagsys[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_metMag"+name);
    hExp_metMagsys[i] = (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_metMag"+name);
    hRatio_Exp_CS_metMagsys[i] = (TH1D*)hExp_metMagsys[i]->Clone("Ratio");
    hRatio_Exp_CS_metMagsys[i]->Divide(hCS_metMagsys[i]);

    hCS_metPhisys[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_metPhi"+name);
    hExp_metPhisys[i] = (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_metPhi"+name);
    hRatio_Exp_CS_metPhisys[i] = (TH1D*)hExp_metPhisys[i]->Clone("Ratio");
    hRatio_Exp_CS_metPhisys[i]->Divide(hCS_metPhisys[i]);

    hCS_jecsys[i] = (TH1D*)file_CS->Get("hYields_"+csStr+"_jec"+name);
    hExp_jecsys[i] = (TH1D*)file_HadTauLL->Get("hYields_"+expStr+"_jec"+name);
    hRatio_Exp_CS_jecsys[i] = (TH1D*)hExp_jecsys[i]->Clone("Ratio");
    hRatio_Exp_CS_jecsys[i]->Divide(hCS_jecsys[i]);
  }
  double sum_pdfUncsys_up = 0, sum_pdfUncsys_cen = 0, sum_pdfUncsys_dn = 0;
  for(unsigned int i=1; i<hRatio_Exp_CS_pdfUncsys_cen->GetNbinsX(); i++)
  {
     sum_pdfUncsys_cen += hRatio_Exp_CS_pdfUncsys_cen->GetBinContent(i);
     sum_pdfUncsys_up += hRatio_Exp_CS_pdfUncsys[0]->GetBinContent(i);
     sum_pdfUncsys_dn += hRatio_Exp_CS_pdfUncsys[1]->GetBinContent(i);
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
