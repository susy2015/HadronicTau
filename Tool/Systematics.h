#ifndef SYSTEMATICS_H
#define SYSTEMATICS_H

#include <iostream>
#include <vector>
#include "TH1.h"
#include "TFile.h"

class Systematics{
 public:
  //Constructor
  Systematics();

  //Member Data
  std::vector<double> ISRsys_up_mu;
  std::vector<double> ISRsys_down_mu;
  std::vector<double> ISRsys_up_el;
  std::vector<double> ISRsys_down_el;

  TFile *file_CS_ISR;
  TFile *file_LL_ISR;
};

Systematics::Systematics()
{

  file_CS_ISR = new TFile("Mixed_CS_ISR.root");
  file_LL_ISR = new TFile("Mixed_HadTauLL_ISR.root");
  
  
  TH1D *hCS_mu = (TH1D*)file_CS_ISR->Get("hYields_mu");
  TH1D *hCS_ele = (TH1D*)file_CS_ISR->Get("hYields_el");
  TH1D *hHadtau = (TH1D*)file_LL_ISR->Get("hYields_tau");
  TH1D *hLL = (TH1D*)file_LL_ISR->Get("hYields_LL");
  TH1D *hRatio_Hadtau_mu = (TH1D*)hHadtau->Clone("Ratio");
  hRatio_Hadtau_mu->Divide(hCS_mu);
  TH1D *hRatio_LL_mu = (TH1D*)hLL->Clone("Ratio");
  hRatio_LL_mu->Divide(hCS_mu);
  TH1D *hRatio_Hadtau_ele = (TH1D*)hHadtau->Clone("Ratio");
  hRatio_Hadtau_ele->Divide(hCS_ele);
  TH1D *hRatio_LL_ele = (TH1D*)hLL->Clone("Ratio");
  hRatio_LL_ele->Divide(hCS_ele);


  const unsigned int kNDists = 2;
  TH1* hCSISR_mu[kNDists];
  TH1* hHadtauISR[kNDists];
  TH1* hCSISR_ele[kNDists];
  TH1* hLLISR[kNDists];
  TH1* hRatio_Hadtau_muISR[kNDists];
  TH1* hRatio_LL_muISR[kNDists];
  TH1* hRatio_Hadtau_eleISR[kNDists];
  TH1* hRatio_LL_eleISR[kNDists];

  for(unsigned int i = 0; i < kNDists; ++i) {
    TString name = "";
    if(      i == 0 ) name = "up";
    else if( i == 1 ) name = "down";
    // Get histograms from file                                                                                                                
    hCSISR_mu[i] = (TH1D*)file_CS_ISR->Get("hYields_mu_isr"+name);
    hCSISR_ele[i] = (TH1D*)file_CS_ISR->Get("hYields_el_isr"+name);
    hHadtauISR[i] = (TH1D*)file_LL_ISR->Get("hYields_tau_isr"+name);
    hLLISR[i] = (TH1D*)file_LL_ISR->Get("hYields_LL_isr"+name);
    //Ratio                                                                                                                                   
    hRatio_Hadtau_muISR[i] = static_cast<TH1*>(hHadtauISR[i]->Clone("Ratio"));
    hRatio_Hadtau_muISR[i]->Divide(hCSISR_mu[i]);
    hRatio_LL_muISR[i] = static_cast<TH1*>(hLLISR[i]->Clone("Ratio"));
    hRatio_LL_muISR[i]->Divide(hCSISR_mu[i]);
    hRatio_Hadtau_eleISR[i] = static_cast<TH1*>(hHadtauISR[i]->Clone("Ratio"));
    hRatio_Hadtau_eleISR[i]->Divide(hCSISR_ele[i]);
    hRatio_LL_eleISR[i] = static_cast<TH1*>(hLLISR[i]->Clone("Ratio"));
    hRatio_LL_eleISR[i]->Divide(hCSISR_ele[i]);
  }
  /*
  std::cout << std::setprecision(3) << std::fixed;
  std::cout<<"Muon CS"<<std::endl;
  std::cout<<"Bin \t\t Ratio \t\tRatio Up \tRatio Down \t SysUp\t\t SysDown"<<std::endl;
  */
  for(unsigned i=1; i<= hRatio_Hadtau_mu->GetNbinsX();i++)
    {
      ISRsys_up_mu.push_back(std::fabs(hRatio_Hadtau_mu->GetBinContent(i)-hRatio_Hadtau_muISR[0]->GetBinContent(i)));
      ISRsys_down_mu.push_back(std::fabs(hRatio_Hadtau_mu->GetBinContent(i)-hRatio_Hadtau_muISR[1]->GetBinContent(i)));

    }
  /*
  std::cout<<std::endl;
  std::cout<<"Electron CS"<<std::endl;
  std::cout<<"Bin \t\t Ratio \t\tRatio Up \tRatio Down \t SysUp\t\t SysDown"<<std::endl;
  */
  for(unsigned i=1; i<= hRatio_Hadtau_ele->GetNbinsX();i++)
    {

      ISRsys_up_el.push_back(std::fabs(hRatio_Hadtau_ele->GetBinContent(i)-hRatio_Hadtau_eleISR[0]->GetBinContent(i)));
      ISRsys_down_el.push_back(std::fabs(hRatio_Hadtau_ele->GetBinContent(i)-hRatio_Hadtau_eleISR[1]->GetBinContent(i)));
    }


}
#endif
