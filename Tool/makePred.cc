#include "Systematics.h"
#include "CommonShared.h"
#include "Math/QuantFuncMathCore.h"
#include <fstream>

const double mc_lumi = 36813.714859265;
const double data_lumi = 36813.714859265;

const double scale_mc = data_lumi/mc_lumi;

TH1* hTF_tau_mu = 0, * hTF_tau_ele = 0, * hTF_LL_mu = 0, * hTF_LL_ele = 0, * hTF_sum_mu = 0, * hTF_sum_ele = 0;
TH1* hTF_tau_mu_SFup = 0, * hTF_tau_ele_SFup = 0, * hTF_LL_mu_SFup = 0, * hTF_LL_ele_SFup = 0, * hTF_sum_mu_SFup = 0, * hTF_sum_ele_SFup = 0;
TH1* hTF_tau_mu_SFdn = 0, * hTF_tau_ele_SFdn = 0, * hTF_LL_mu_SFdn = 0, * hTF_LL_ele_SFdn = 0, * hTF_sum_mu_SFdn = 0, * hTF_sum_ele_SFdn = 0;
TH1* hYields_lepSF_ratio_tau = 0, * hYields_lepSF_ratio_LL = 0, * hYields_lepSF_ratio_sum = 0;
TH1* hYields_lepSF_ratio_tau_SFup = 0, * hYields_lepSF_ratio_LL_SFup = 0, * hYields_lepSF_ratio_sum_SFup = 0;
TH1* hYields_lepSF_ratio_tau_SFdn = 0, * hYields_lepSF_ratio_LL_SFdn = 0, * hYields_lepSF_ratio_sum_SFdn = 0;
TH1* hYields_Data_mu = 0, * hYields_Data_ele = 0;

TH1* hYields_MC_mu =0, * hYields_MC_ele =0;
TH1* hYields_MC_mu_SFup =0, * hYields_MC_mu_SFdn =0, * hYields_MC_ele_SFup = 0, * hYields_MC_ele_SFdn = 0;
TH1* hYields_tau =0, * hYields_LL =0, * hYields_sum =0;
TH1* hYields_Veto_tau =0, * hYields_Veto_tau_SF =0, * hYields_Pass_tau =0;
TH1* hYields_Veto_LL =0, * hYields_Veto_LL_SF =0, * hYields_Pass_LL =0;
TH1* hYields_Veto_tau_SFup = 0, *hYields_Veto_tau_SFdn =0;
TH1* hYields_Veto_LL_SFup = 0, *hYields_Veto_LL_SFdn =0;

// key = "hadtau" or "lostle"
bool relSys_for_ZERO = true;
bool doVector_dataCard = false;
// Do adJustBins_merge before calculating ANY ratio or summation and others...
// Only adjustBins_merge for TF factor calcuation and systematics -- do NOT do this on data!
bool do_mergeBins = false;
void combCS_pred(const std::string key = "hadtau", const std::string sel_CS_str = "comb");

const std::string filename_CS = "Mix_CS.root", filename_HadTauLL = "Mix_HadTauLL.root";

const bool enable_prtTFfactors = true, prtExtraInfo = true;

void prtTFfactors();

void predFromCSonHadtauLL();

void predLLTry();

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr <<"Please give 1 arguments "<<" comb card or not"<<std::endl;
    std::cerr <<" ./makePred comb[mu, ele]" << std::endl;
    return -1;
  }

  const char *card_type_str= argv[1];

  TFile *file_Data_CS = new TFile("Data_MET_CS.root");
  TFile *file_MC_CS = new TFile("Mix_CS.root");
  TFile *file_HadTauLL = new TFile("Mix_HadTauLL.root");

  hYields_MC_mu = (TH1D*) file_MC_CS->Get("hYields_mu"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_MC_mu);
  hYields_MC_ele = (TH1D*) file_MC_CS->Get("hYields_el"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_MC_ele);

  hYields_MC_mu_SFup = (TH1D*) file_MC_CS->Get("hYields_mu_SFup"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_MC_mu_SFup);
  hYields_MC_ele_SFup = (TH1D*) file_MC_CS->Get("hYields_el_SFup"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_MC_ele_SFup);

  hYields_MC_mu_SFdn = (TH1D*) file_MC_CS->Get("hYields_mu_SFdn"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_MC_mu_SFdn);
  hYields_MC_ele_SFdn = (TH1D*) file_MC_CS->Get("hYields_el_SFdn"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_MC_ele_SFdn);

  hYields_tau = (TH1D*) file_HadTauLL->Get("hYields_tau");  if( do_mergeBins ) predSpec::adjustBins_merge(hYields_tau);
  hYields_LL = (TH1D*) file_HadTauLL->Get("hYields_LL"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_LL);
  hYields_sum = static_cast<TH1*>(hYields_tau->Clone("hYields_sum")); hYields_sum->Add(hYields_LL); //if( do_mergeBins ) predSpec::adjustBins_merge(hYields_sum);

  hYields_Veto_tau = (TH1D*) file_HadTauLL->Get("hYields_Veto_tau"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_Veto_tau);
  hYields_Veto_tau_SF = (TH1D*) file_HadTauLL->Get("hYields_Veto_tau_SF"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_Veto_tau_SF);
  hYields_Pass_tau = (TH1D*) file_HadTauLL->Get("hYields_Pass_tau"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_Pass_tau);
  hYields_Veto_tau_SFup = (TH1D*) file_HadTauLL->Get("hYields_Veto_tau_SFup"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_Veto_tau_SFup);
  hYields_Veto_tau_SFdn = (TH1D*) file_HadTauLL->Get("hYields_Veto_tau_SFdn"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_Veto_tau_SFdn);

  hYields_Veto_LL = (TH1D*) file_HadTauLL->Get("hYields_Veto_LL"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_Veto_LL);
  hYields_Veto_LL_SF = (TH1D*) file_HadTauLL->Get("hYields_Veto_LL_SF"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_Veto_LL_SF);
  hYields_Pass_LL = (TH1D*) file_HadTauLL->Get("hYields_Pass_LL"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_Pass_LL);
  hYields_Veto_LL_SFup = (TH1D*) file_HadTauLL->Get("hYields_Veto_LL_SFup"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_Veto_LL_SFup);
  hYields_Veto_LL_SFdn = (TH1D*) file_HadTauLL->Get("hYields_Veto_LL_SFdn"); if( do_mergeBins ) predSpec::adjustBins_merge(hYields_Veto_LL_SFdn);

// Do NOT adjustBins_merge for data!
  hYields_Data_mu = (TH1D*) file_Data_CS->Get("hYields_mu");
  hYields_Data_ele = (TH1D*) file_Data_CS->Get("hYields_el");

  hYields_MC_mu->Scale(scale_mc);
  hYields_MC_ele->Scale(scale_mc);

  hYields_MC_mu_SFup->Scale(scale_mc);
  hYields_MC_ele_SFup->Scale(scale_mc);

  hYields_MC_mu_SFdn->Scale(scale_mc);
  hYields_MC_ele_SFdn->Scale(scale_mc);

  hYields_tau->Scale(scale_mc);
  hYields_LL->Scale(scale_mc);
  hYields_sum->Scale(scale_mc);

  hYields_Veto_tau->Scale(scale_mc);
  hYields_Veto_tau_SF->Scale(scale_mc);
  hYields_Pass_tau->Scale(scale_mc);
  hYields_Veto_tau_SFup->Scale(scale_mc);
  hYields_Veto_tau_SFdn->Scale(scale_mc);

  hYields_Veto_LL->Scale(scale_mc);
  hYields_Veto_LL_SF->Scale(scale_mc);
  hYields_Pass_LL->Scale(scale_mc);
  hYields_Veto_LL_SFup->Scale(scale_mc);
  hYields_Veto_LL_SFdn->Scale(scale_mc);

// (pass+veto-veto_SF)/pass --> to be directly applied on the TF
  hYields_lepSF_ratio_tau = static_cast<TH1*>(hYields_Pass_tau->Clone("hYields_lepSF_ratio_tau"));
  hYields_lepSF_ratio_tau->Add(hYields_Veto_tau);
  hYields_lepSF_ratio_tau->Add(hYields_Veto_tau_SF, -1);
  hYields_lepSF_ratio_tau->Divide(hYields_Pass_tau);

  hYields_lepSF_ratio_LL = static_cast<TH1*>(hYields_Pass_LL->Clone("hYields_lepSF_ratio_LL"));
  hYields_lepSF_ratio_LL->Add(hYields_Veto_LL);
  hYields_lepSF_ratio_LL->Add(hYields_Veto_LL_SF, -1);
  hYields_lepSF_ratio_LL->Divide(hYields_Pass_LL);

  hYields_lepSF_ratio_sum = static_cast<TH1*>(hYields_Pass_tau->Clone("hYields_lepSF_ratio_sum"));
  hYields_lepSF_ratio_sum->Add(hYields_Pass_LL);
  TH1* hYields_lepSF_ratio_sum_tmp = (TH1*) hYields_lepSF_ratio_sum->Clone("hYields_lepSF_ratio_sum_tmp"); // sum_tmp is the sum of Pass_tau and Pass_LL
  hYields_lepSF_ratio_sum->Add(hYields_Veto_tau);
  hYields_lepSF_ratio_sum->Add(hYields_Veto_LL);
  hYields_lepSF_ratio_sum->Add(hYields_Veto_tau_SF, -1);
  hYields_lepSF_ratio_sum->Add(hYields_Veto_LL_SF, -1);
  hYields_lepSF_ratio_sum->Divide(hYields_lepSF_ratio_sum_tmp); 

  hYields_lepSF_ratio_tau_SFup = static_cast<TH1*>(hYields_Pass_tau->Clone("hYields_lepSF_ratio_tau_SFup"));
  hYields_lepSF_ratio_tau_SFup->Add(hYields_Veto_tau);
  hYields_lepSF_ratio_tau_SFup->Add(hYields_Veto_tau_SFup, -1);
  hYields_lepSF_ratio_tau_SFup->Divide(hYields_Pass_tau);

  hYields_lepSF_ratio_tau_SFdn = static_cast<TH1*>(hYields_Pass_tau->Clone("hYields_lepSF_ratio_tau_SFdn"));
  hYields_lepSF_ratio_tau_SFdn->Add(hYields_Veto_tau);
  hYields_lepSF_ratio_tau_SFdn->Add(hYields_Veto_tau_SFdn, -1);
  hYields_lepSF_ratio_tau_SFdn->Divide(hYields_Pass_tau);

  hYields_lepSF_ratio_LL_SFup = static_cast<TH1*>(hYields_Pass_LL->Clone("hYields_lepSF_ratio_LL_SFup"));
  hYields_lepSF_ratio_LL_SFup->Add(hYields_Veto_LL);
  hYields_lepSF_ratio_LL_SFup->Add(hYields_Veto_LL_SFup, -1);
  hYields_lepSF_ratio_LL_SFup->Divide(hYields_Pass_LL);

  hYields_lepSF_ratio_LL_SFdn = static_cast<TH1*>(hYields_Pass_LL->Clone("hYields_lepSF_ratio_LL_SFdn"));
  hYields_lepSF_ratio_LL_SFdn->Add(hYields_Veto_LL);
  hYields_lepSF_ratio_LL_SFdn->Add(hYields_Veto_LL_SFdn, -1);
  hYields_lepSF_ratio_LL_SFdn->Divide(hYields_Pass_LL);

  hYields_lepSF_ratio_sum_SFup = static_cast<TH1*>(hYields_Pass_tau->Clone("hYields_lepSF_ratio_sum_SFup"));
  hYields_lepSF_ratio_sum_SFup->Add(hYields_Pass_LL);
  TH1* hYields_lepSF_ratio_sum_SFup_tmp = (TH1*) hYields_lepSF_ratio_sum_SFup->Clone("hYields_lepSF_ratio_sum_SFup_tmp"); // sum_SFup_tmp is the sum_SFup of Pass_tau and Pass_LL
  hYields_lepSF_ratio_sum_SFup->Add(hYields_Veto_tau);
  hYields_lepSF_ratio_sum_SFup->Add(hYields_Veto_LL);
  hYields_lepSF_ratio_sum_SFup->Add(hYields_Veto_tau_SFup, -1);
  hYields_lepSF_ratio_sum_SFup->Add(hYields_Veto_LL_SFup, -1);
  hYields_lepSF_ratio_sum_SFup->Divide(hYields_lepSF_ratio_sum_SFup_tmp); 

  hYields_lepSF_ratio_sum_SFdn = static_cast<TH1*>(hYields_Pass_tau->Clone("hYields_lepSF_ratio_sum_SFdn"));
  hYields_lepSF_ratio_sum_SFdn->Add(hYields_Pass_LL);
  TH1* hYields_lepSF_ratio_sum_SFdn_tmp = (TH1*) hYields_lepSF_ratio_sum_SFdn->Clone("hYields_lepSF_ratio_sum_SFdn_tmp"); // sum_SFdn_tmp is the sum_SFdn of Pass_tau and Pass_LL
  hYields_lepSF_ratio_sum_SFdn->Add(hYields_Veto_tau);
  hYields_lepSF_ratio_sum_SFdn->Add(hYields_Veto_LL);
  hYields_lepSF_ratio_sum_SFdn->Add(hYields_Veto_tau_SFdn, -1);
  hYields_lepSF_ratio_sum_SFdn->Add(hYields_Veto_LL_SFdn, -1);
  hYields_lepSF_ratio_sum_SFdn->Divide(hYields_lepSF_ratio_sum_SFdn_tmp); 

// Derive the translation factors
  hTF_tau_mu = static_cast<TH1*>(hYields_tau->Clone("hTF_tau_mu")); hTF_tau_mu->Divide(hYields_MC_mu);
  hTF_tau_ele = static_cast<TH1*>(hYields_tau->Clone("hTF_tau_ele")); hTF_tau_ele->Divide(hYields_MC_ele);
  hTF_LL_mu = static_cast<TH1*>(hYields_LL->Clone("hTF_LL_mu")); hTF_LL_mu->Divide(hYields_MC_mu);
  hTF_LL_ele = static_cast<TH1*>(hYields_LL->Clone("hTF_LL_ele")); hTF_LL_ele->Divide(hYields_MC_ele);
  hTF_sum_mu = static_cast<TH1*>(hYields_sum->Clone("hTF_sum_mu")); hTF_sum_mu->Divide(hYields_MC_mu);
  hTF_sum_ele = static_cast<TH1*>(hYields_sum->Clone("hTF_sum_ele")); hTF_sum_ele->Divide(hYields_MC_ele);

// Derive the translation factors
  hTF_tau_mu_SFup = static_cast<TH1*>(hYields_tau->Clone("hTF_tau_mu_SFup")); hTF_tau_mu_SFup->Divide(hYields_MC_mu_SFup);
  hTF_tau_ele_SFup = static_cast<TH1*>(hYields_tau->Clone("hTF_tau_ele_SFup")); hTF_tau_ele_SFup->Divide(hYields_MC_ele_SFup);
  hTF_LL_mu_SFup = static_cast<TH1*>(hYields_LL->Clone("hTF_LL_mu_SFup")); hTF_LL_mu_SFup->Divide(hYields_MC_mu_SFup);
  hTF_LL_ele_SFup = static_cast<TH1*>(hYields_LL->Clone("hTF_LL_ele_SFup")); hTF_LL_ele_SFup->Divide(hYields_MC_ele_SFup);
  hTF_sum_mu_SFup = static_cast<TH1*>(hYields_sum->Clone("hTF_sum_mu_SFup")); hTF_sum_mu_SFup->Divide(hYields_MC_mu_SFup);
  hTF_sum_ele_SFup = static_cast<TH1*>(hYields_sum->Clone("hTF_sum_ele_SFup")); hTF_sum_ele_SFup->Divide(hYields_MC_ele_SFup);

// Derive the translation factors
  hTF_tau_mu_SFdn = static_cast<TH1*>(hYields_tau->Clone("hTF_tau_mu_SFdn")); hTF_tau_mu_SFdn->Divide(hYields_MC_mu_SFdn);
  hTF_tau_ele_SFdn = static_cast<TH1*>(hYields_tau->Clone("hTF_tau_ele_SFdn")); hTF_tau_ele_SFdn->Divide(hYields_MC_ele_SFdn);
  hTF_LL_mu_SFdn = static_cast<TH1*>(hYields_LL->Clone("hTF_LL_mu_SFdn")); hTF_LL_mu_SFdn->Divide(hYields_MC_mu_SFdn);
  hTF_LL_ele_SFdn = static_cast<TH1*>(hYields_LL->Clone("hTF_LL_ele_SFdn")); hTF_LL_ele_SFdn->Divide(hYields_MC_ele_SFdn);
  hTF_sum_mu_SFdn = static_cast<TH1*>(hYields_sum->Clone("hTF_sum_mu_SFdn")); hTF_sum_mu_SFdn->Divide(hYields_MC_mu_SFdn);
  hTF_sum_ele_SFdn = static_cast<TH1*>(hYields_sum->Clone("hTF_sum_ele_SFdn")); hTF_sum_ele_SFdn->Divide(hYields_MC_ele_SFdn);

  if( enable_prtTFfactors ) prtTFfactors();
  if( prtExtraInfo ){ predFromCSonHadtauLL(); predLLTry(); }

  combCS_pred("hadtau", card_type_str);
  combCS_pred("lostle", card_type_str);

  return 1;
}

// sel_CS_str : "comb" -- combined   "mu" for muon only   "ele" for electron only
void combCS_pred(const std::string key, const std::string sel_CS_str)
{
  TH1 * hTF_local_mu = key == "hadtau" ? (TH1*) hTF_tau_mu->Clone() : (TH1*) hTF_LL_mu->Clone();
  TH1 * hTF_local_ele = key == "hadtau" ? (TH1*) hTF_tau_ele->Clone() : (TH1*) hTF_LL_ele->Clone();
  TH1 * hYields_lepSF_ratio_local = key == "hadtau" ? (TH1*) hYields_lepSF_ratio_tau->Clone() : (TH1*) hYields_lepSF_ratio_LL->Clone();
  Systematics systCalc_mu = key == "hadtau"? Systematics("tau", "mu", filename_CS.c_str(), filename_HadTauLL.c_str(), do_mergeBins) : Systematics("LL", "mu", filename_CS.c_str(), filename_HadTauLL.c_str(), do_mergeBins);
  Systematics systCalc_el = key == "hadtau"? Systematics("tau", "el", filename_CS.c_str(), filename_HadTauLL.c_str(), do_mergeBins) : Systematics("LL", "el", filename_CS.c_str(), filename_HadTauLL.c_str(), do_mergeBins);

  TH1 * hTF_local_mu_SFup = key == "hadtau" ? (TH1*) hTF_tau_mu_SFup->Clone() : (TH1*) hTF_LL_mu_SFup->Clone();
  TH1 * hTF_local_ele_SFup = key == "hadtau" ? (TH1*) hTF_tau_ele_SFup->Clone() : (TH1*) hTF_LL_ele_SFup->Clone();
  TH1 * hYields_lepSF_ratio_local_SFup = key == "hadtau" ? (TH1*) hYields_lepSF_ratio_tau_SFup->Clone() : (TH1*) hYields_lepSF_ratio_LL_SFup->Clone();

  TH1 * hTF_local_mu_SFdn = key == "hadtau" ? (TH1*) hTF_tau_mu_SFdn->Clone() : (TH1*) hTF_LL_mu_SFdn->Clone();
  TH1 * hTF_local_ele_SFdn = key == "hadtau" ? (TH1*) hTF_tau_ele_SFdn->Clone() : (TH1*) hTF_LL_ele_SFdn->Clone();
  TH1 * hYields_lepSF_ratio_local_SFdn = key == "hadtau" ? (TH1*) hYields_lepSF_ratio_tau_SFdn->Clone() : (TH1*) hYields_lepSF_ratio_LL_SFdn->Clone();

  std::vector<double> systErr_rel_ISR_up_muVec, systErr_rel_ISR_dn_muVec, systErr_rel_ISR_up_eleVec, systErr_rel_ISR_dn_eleVec;
  for(unsigned int is=0; is<systCalc_mu.ISRsys_up.size(); is++)
  {
    const double rel_sys_up_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.ISRsys_up[is]/hTF_local_mu->GetBinContent(is+1);
    const double rel_sys_dn_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.ISRsys_dn[is]/hTF_local_mu->GetBinContent(is+1);
    const double rel_sys_up_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.ISRsys_up[is]/hTF_local_ele->GetBinContent(is+1);
    const double rel_sys_dn_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.ISRsys_dn[is]/hTF_local_ele->GetBinContent(is+1);

    systErr_rel_ISR_up_muVec.push_back(rel_sys_up_mu);
    systErr_rel_ISR_dn_muVec.push_back(rel_sys_dn_mu);
    systErr_rel_ISR_up_eleVec.push_back(rel_sys_up_el);
    systErr_rel_ISR_dn_eleVec.push_back(rel_sys_dn_el);
  }

  std::vector<double> systErr_rel_bTag_up_muVec, systErr_rel_bTag_dn_muVec, systErr_rel_bTag_up_eleVec, systErr_rel_bTag_dn_eleVec;
  for(unsigned int is=0; is<systCalc_mu.bSFsys_up.size(); is++)
  {
    const double rel_sys_up_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.bSFsys_up[is]/hTF_local_mu->GetBinContent(is+1);
    const double rel_sys_dn_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.bSFsys_dn[is]/hTF_local_mu->GetBinContent(is+1);
    const double rel_sys_up_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.bSFsys_up[is]/hTF_local_ele->GetBinContent(is+1);
    const double rel_sys_dn_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.bSFsys_dn[is]/hTF_local_ele->GetBinContent(is+1);

    systErr_rel_bTag_up_muVec.push_back(rel_sys_up_mu);
    systErr_rel_bTag_dn_muVec.push_back(rel_sys_dn_mu);
    systErr_rel_bTag_up_eleVec.push_back(rel_sys_up_el);
    systErr_rel_bTag_dn_eleVec.push_back(rel_sys_dn_el);
  }

  std::vector<double> systErr_rel_scaleUnc_up_muVec, systErr_rel_scaleUnc_dn_muVec, systErr_rel_scaleUnc_up_eleVec, systErr_rel_scaleUnc_dn_eleVec;
  for(unsigned int is=0; is<systCalc_mu.scaleUncsys_up.size(); is++)
  {
    const double rel_sys_up_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.scaleUncsys_up[is]/hTF_local_mu->GetBinContent(is+1);
    const double rel_sys_dn_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.scaleUncsys_dn[is]/hTF_local_mu->GetBinContent(is+1);
    const double rel_sys_up_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.scaleUncsys_up[is]/hTF_local_ele->GetBinContent(is+1);
    const double rel_sys_dn_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.scaleUncsys_dn[is]/hTF_local_ele->GetBinContent(is+1);

    systErr_rel_scaleUnc_up_muVec.push_back(rel_sys_up_mu);
    systErr_rel_scaleUnc_dn_muVec.push_back(rel_sys_dn_mu);
    systErr_rel_scaleUnc_up_eleVec.push_back(rel_sys_up_el);
    systErr_rel_scaleUnc_dn_eleVec.push_back(rel_sys_dn_el);
  }

  std::vector<double> systErr_rel_pdfUnc_up_muVec, systErr_rel_pdfUnc_dn_muVec, systErr_rel_pdfUnc_up_eleVec, systErr_rel_pdfUnc_dn_eleVec;
  for(unsigned int is=0; is<systCalc_mu.pdfUncsys_up.size(); is++)
  {
    double rel_sys_up_mu = systCalc_mu.pdfUncsys_cen[is] == 0 ? 0 : systCalc_mu.pdfUncsys_up[is]/systCalc_mu.pdfUncsys_cen[is];
    double rel_sys_dn_mu = systCalc_mu.pdfUncsys_cen[is] == 0 ? 0 : systCalc_mu.pdfUncsys_dn[is]/systCalc_mu.pdfUncsys_cen[is];
    double rel_sys_up_el = systCalc_el.pdfUncsys_cen[is] == 0 ? 0 : systCalc_el.pdfUncsys_up[is]/systCalc_el.pdfUncsys_cen[is];
    double rel_sys_dn_el = systCalc_el.pdfUncsys_cen[is] == 0 ? 0 : systCalc_el.pdfUncsys_dn[is]/systCalc_el.pdfUncsys_cen[is];
/*
    const double rel_sys_up_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.pdfUncsys_up[is]/hTF_local_mu->GetBinContent(is+1);
    const double rel_sys_dn_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.pdfUncsys_dn[is]/hTF_local_mu->GetBinContent(is+1);
    const double rel_sys_up_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.pdfUncsys_up[is]/hTF_local_ele->GetBinContent(is+1);
    const double rel_sys_dn_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.pdfUncsys_dn[is]/hTF_local_ele->GetBinContent(is+1);
*/
    if( rel_sys_up_mu >= 1.0 ) rel_sys_up_mu = 0.0;
    if( rel_sys_dn_mu >= 1.0 ) rel_sys_dn_mu = 0.0;
    if( rel_sys_up_el >= 1.0 ) rel_sys_up_el = 0.0;
    if( rel_sys_dn_el >= 1.0 ) rel_sys_dn_el = 0.0;
    systErr_rel_pdfUnc_up_muVec.push_back(rel_sys_up_mu);
    systErr_rel_pdfUnc_dn_muVec.push_back(rel_sys_dn_mu);
    systErr_rel_pdfUnc_up_eleVec.push_back(rel_sys_up_el);
    systErr_rel_pdfUnc_dn_eleVec.push_back(rel_sys_dn_el);
  }

  std::vector<double> systErr_rel_metMag_up_muVec, systErr_rel_metMag_dn_muVec, systErr_rel_metMag_up_eleVec, systErr_rel_metMag_dn_eleVec;
  for(unsigned int is=0; is<systCalc_mu.metMagsys_up.size(); is++)
  {
    double rel_sys_up_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.metMagsys_up[is]/hTF_local_mu->GetBinContent(is+1);
    double rel_sys_dn_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.metMagsys_dn[is]/hTF_local_mu->GetBinContent(is+1);
    double rel_sys_up_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.metMagsys_up[is]/hTF_local_ele->GetBinContent(is+1);
    double rel_sys_dn_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.metMagsys_dn[is]/hTF_local_ele->GetBinContent(is+1);

    if( rel_sys_up_mu >= 1.0 ) rel_sys_up_mu = 0.0;
    if( rel_sys_dn_mu >= 1.0 ) rel_sys_dn_mu = 0.0;
    if( rel_sys_up_el >= 1.0 ) rel_sys_up_el = 0.0;
    if( rel_sys_dn_el >= 1.0 ) rel_sys_dn_el = 0.0;
    systErr_rel_metMag_up_muVec.push_back(rel_sys_up_mu);
    systErr_rel_metMag_dn_muVec.push_back(rel_sys_dn_mu);
    systErr_rel_metMag_up_eleVec.push_back(rel_sys_up_el);
    systErr_rel_metMag_dn_eleVec.push_back(rel_sys_dn_el);
  }

  std::vector<double> systErr_rel_metPhi_up_muVec, systErr_rel_metPhi_dn_muVec, systErr_rel_metPhi_up_eleVec, systErr_rel_metPhi_dn_eleVec;
  for(unsigned int is=0; is<systCalc_mu.metPhisys_up.size(); is++)
  {
    const double rel_sys_up_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.metPhisys_up[is]/hTF_local_mu->GetBinContent(is+1);
    const double rel_sys_dn_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.metPhisys_dn[is]/hTF_local_mu->GetBinContent(is+1);
    const double rel_sys_up_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.metPhisys_up[is]/hTF_local_ele->GetBinContent(is+1);
    const double rel_sys_dn_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.metPhisys_dn[is]/hTF_local_ele->GetBinContent(is+1);

    systErr_rel_metPhi_up_muVec.push_back(rel_sys_up_mu);
    systErr_rel_metPhi_dn_muVec.push_back(rel_sys_dn_mu);
    systErr_rel_metPhi_up_eleVec.push_back(rel_sys_up_el);
    systErr_rel_metPhi_dn_eleVec.push_back(rel_sys_dn_el);
  }

  std::vector<double> systErr_rel_jec_up_muVec, systErr_rel_jec_dn_muVec, systErr_rel_jec_up_eleVec, systErr_rel_jec_dn_eleVec;
  for(unsigned int is=0; is<systCalc_mu.jecsys_up.size(); is++)
  {
    double rel_sys_up_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.jecsys_up[is]/hTF_local_mu->GetBinContent(is+1);
    double rel_sys_dn_mu = hTF_local_mu->GetBinContent(is+1) == 0 ? 0 : systCalc_mu.jecsys_dn[is]/hTF_local_mu->GetBinContent(is+1);
    double rel_sys_up_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.jecsys_up[is]/hTF_local_ele->GetBinContent(is+1);
    double rel_sys_dn_el = hTF_local_ele->GetBinContent(is+1) == 0 ? 0 : systCalc_el.jecsys_dn[is]/hTF_local_ele->GetBinContent(is+1);

    if( rel_sys_up_mu >= 1.0 ) rel_sys_up_mu = 0.0;
    if( rel_sys_dn_mu >= 1.0 ) rel_sys_dn_mu = 0.0;
    if( rel_sys_up_el >= 1.0 ) rel_sys_up_el = 0.0;
    if( rel_sys_dn_el >= 1.0 ) rel_sys_dn_el = 0.0;
    systErr_rel_jec_up_muVec.push_back(rel_sys_up_mu);
    systErr_rel_jec_dn_muVec.push_back(rel_sys_dn_mu);
    systErr_rel_jec_up_eleVec.push_back(rel_sys_up_el);
    systErr_rel_jec_dn_eleVec.push_back(rel_sys_dn_el);
  }

  std::cout<< std::setprecision(3)<<std::fixed;
  std::cout<<"\n"<<key.c_str()<<" predictions"<<std::endl;
  std::cout<<"Bin \t mu CS \t ele CS \t TF("<<key.c_str()<<"/CS)@mu \t TF("<<key.c_str()<<"/CS)@ele \t\t Pred("<<key.c_str()<<")@mu \t\t\t Pred("<<key.c_str()<<")@ele \t\t\t Pred("<<key.c_str()<<")@avg"<<std::endl;
// Combine mu and ele together
// TF_mu/2.0 * CS_mu + TF_ele/2.0 * CS_ele
//      A    * CS_mu +       B    * CS_ele
// Syst : sqrt( dA*dA * CS_mu*CS_mu + dB*dB * CS_ele*CS_ele )
// Stat : ?
// Using Poisson errors: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars
  const double alpha = 1 - 0.6827; // the "68.27%" intervals
// All errors stored here are relative
  std::vector<double> cached_rateVec, cached_stat_upErrVec, cached_stat_dnErrVec, cached_systErr_TF_statVec;
  std::vector<double> cached_stat_abs_upErrVec, cached_stat_abs_dnErrVec;
  std::vector<double> dummyVec;
  std::vector<double> cached_systErr_ISR_upErrVec, cached_systErr_ISR_dnErrVec;
  std::vector<double> cached_systErr_bTag_upErrVec, cached_systErr_bTag_dnErrVec;
  std::vector<double> cached_systErr_scaleUnc_upErrVec, cached_systErr_scaleUnc_dnErrVec;
  std::vector<double> cached_systErr_pdfUnc_upErrVec, cached_systErr_pdfUnc_dnErrVec;
  std::vector<double> cached_systErr_metMag_upErrVec, cached_systErr_metMag_dnErrVec;
  std::vector<double> cached_systErr_metPhi_upErrVec, cached_systErr_metPhi_dnErrVec;
  std::vector<double> cached_systErr_jec_upErrVec, cached_systErr_jec_dnErrVec;
  std::vector<double> cached_systErr_SF_upErrVec, cached_systErr_SF_dnErrVec;

  std::vector<double> cached_fin_TF_muVec, cached_fin_TFerr_muVec;
  std::vector<double> cached_fin_TF_eleVec, cached_fin_TFerr_eleVec;
  double sum_pred_local_avg = 0, sum_pred_local_from_mu = 0, sum_pred_local_from_ele = 0;
  for(unsigned int i=1; i<=hTF_local_mu->GetNbinsX(); i++)
  {
    const double lepSF_ratio_local = hYields_lepSF_ratio_local->GetBinContent(i);
    const double lepSF_ratio_local_SFup = hYields_lepSF_ratio_local_SFup->GetBinContent(i);
    const double lepSF_ratio_local_SFdn = hYields_lepSF_ratio_local_SFdn->GetBinContent(i);

// mu
    const double cont_TF_local_to_mu = hTF_local_mu->GetBinContent(i), err_TF_local_to_mu = hTF_local_mu->GetBinError(i);
    const double fin_TF_local_to_mu = cont_TF_local_to_mu * lepSF_ratio_local, finErr_TF_local_to_mu = err_TF_local_to_mu * lepSF_ratio_local;
    cached_fin_TF_muVec.push_back(fin_TF_local_to_mu); cached_fin_TFerr_muVec.push_back(finErr_TF_local_to_mu);

    const double cont_TF_local_to_mu_SFup = hTF_local_mu_SFup->GetBinContent(i);
    const double fin_TF_local_to_mu_SFup = cont_TF_local_to_mu_SFup * lepSF_ratio_local_SFup;

    const double cont_TF_local_to_mu_SFdn = hTF_local_mu_SFdn->GetBinContent(i);
    const double fin_TF_local_to_mu_SFdn = cont_TF_local_to_mu_SFdn * lepSF_ratio_local_SFdn;

    const double data_mu = hYields_Data_mu->GetBinContent(i);
    const double data_mu_dn_bound = (data_mu ==0 )? 0. : (ROOT::Math::gamma_quantile(alpha/2, data_mu, 1.0));
    const double data_mu_up_bound = ROOT::Math::gamma_quantile_c(alpha/2, data_mu+1.0, 1.0);
    const double data_mu_dn_err = data_mu - data_mu_dn_bound;
    const double data_mu_up_err = data_mu_up_bound - data_mu;

    const double pred_local_from_mu = data_mu * fin_TF_local_to_mu;
    const double pred_local_from_mu_stat_dn_err = data_mu_dn_err * fin_TF_local_to_mu;
    const double pred_local_from_mu_stat_up_err = data_mu_up_err * fin_TF_local_to_mu;
    const double pred_local_from_mu_syst = data_mu * finErr_TF_local_to_mu;
    const double pred_local_from_mu_syst_ISR_up_err = data_mu * fin_TF_local_to_mu * systErr_rel_ISR_up_muVec[i-1];
    const double pred_local_from_mu_syst_ISR_dn_err = data_mu * fin_TF_local_to_mu * systErr_rel_ISR_dn_muVec[i-1];
    const double pred_local_from_mu_syst_bTag_up_err = data_mu * fin_TF_local_to_mu * systErr_rel_bTag_up_muVec[i-1];
    const double pred_local_from_mu_syst_bTag_dn_err = data_mu * fin_TF_local_to_mu * systErr_rel_bTag_dn_muVec[i-1];

    const double pred_local_from_mu_syst_scaleUnc_up_err = data_mu * fin_TF_local_to_mu * systErr_rel_scaleUnc_up_muVec[i-1];
    const double pred_local_from_mu_syst_scaleUnc_dn_err = data_mu * fin_TF_local_to_mu * systErr_rel_scaleUnc_dn_muVec[i-1];

    const double pred_local_from_mu_syst_pdfUnc_up_err = data_mu * fin_TF_local_to_mu * systErr_rel_pdfUnc_up_muVec[i-1];
    const double pred_local_from_mu_syst_pdfUnc_dn_err = data_mu * fin_TF_local_to_mu * systErr_rel_pdfUnc_dn_muVec[i-1];

    const double pred_local_from_mu_syst_metMag_up_err = data_mu * fin_TF_local_to_mu * systErr_rel_metMag_up_muVec[i-1];
    const double pred_local_from_mu_syst_metMag_dn_err = data_mu * fin_TF_local_to_mu * systErr_rel_metMag_dn_muVec[i-1];

    const double pred_local_from_mu_syst_metPhi_up_err = data_mu * fin_TF_local_to_mu * systErr_rel_metPhi_up_muVec[i-1];
    const double pred_local_from_mu_syst_metPhi_dn_err = data_mu * fin_TF_local_to_mu * systErr_rel_metPhi_dn_muVec[i-1];

    const double pred_local_from_mu_syst_jec_up_err = data_mu * fin_TF_local_to_mu * systErr_rel_jec_up_muVec[i-1];
    const double pred_local_from_mu_syst_jec_dn_err = data_mu * fin_TF_local_to_mu * systErr_rel_jec_dn_muVec[i-1];

    const double pred_local_from_mu_syst_SF_up_err = data_mu * std::abs(fin_TF_local_to_mu_SFup - fin_TF_local_to_mu);
    const double pred_local_from_mu_syst_SF_dn_err = data_mu * std::abs(fin_TF_local_to_mu_SFdn - fin_TF_local_to_mu);

// el 
    const double cont_TF_local_to_ele = hTF_local_ele->GetBinContent(i), err_TF_local_to_ele = hTF_local_ele->GetBinError(i);
    const double fin_TF_local_to_ele = cont_TF_local_to_ele * lepSF_ratio_local, finErr_TF_local_to_ele = err_TF_local_to_ele * lepSF_ratio_local;
    cached_fin_TF_eleVec.push_back(fin_TF_local_to_ele); cached_fin_TFerr_eleVec.push_back(finErr_TF_local_to_ele);

    const double cont_TF_local_to_ele_SFup = hTF_local_ele_SFup->GetBinContent(i);
    const double fin_TF_local_to_ele_SFup = cont_TF_local_to_ele_SFup * lepSF_ratio_local_SFup;

    const double cont_TF_local_to_ele_SFdn = hTF_local_ele_SFdn->GetBinContent(i);
    const double fin_TF_local_to_ele_SFdn = cont_TF_local_to_ele_SFdn * lepSF_ratio_local_SFdn;

    const double data_ele = hYields_Data_ele->GetBinContent(i);
    const double data_ele_dn_bound = (data_ele ==0 )? 0. : (ROOT::Math::gamma_quantile(alpha/2, data_ele, 1.0));
    const double data_ele_up_bound = ROOT::Math::gamma_quantile_c(alpha/2, data_ele+1.0, 1.0);
    const double data_ele_dn_err = data_ele - data_ele_dn_bound;
    const double data_ele_up_err = data_ele_up_bound - data_ele;

    const double pred_local_from_ele = data_ele * fin_TF_local_to_ele; 
    const double pred_local_from_ele_stat_dn_err = data_ele_dn_err * fin_TF_local_to_ele;
    const double pred_local_from_ele_stat_up_err = data_ele_up_err * fin_TF_local_to_ele;
    const double pred_local_from_ele_syst = data_ele * finErr_TF_local_to_ele;
    const double pred_local_from_ele_syst_ISR_up_err = data_ele * fin_TF_local_to_ele * systErr_rel_ISR_up_eleVec[i-1];
    const double pred_local_from_ele_syst_ISR_dn_err = data_ele * fin_TF_local_to_ele * systErr_rel_ISR_dn_eleVec[i-1];
    const double pred_local_from_ele_syst_bTag_up_err = data_ele * fin_TF_local_to_ele * systErr_rel_bTag_up_eleVec[i-1];
    const double pred_local_from_ele_syst_bTag_dn_err = data_ele * fin_TF_local_to_ele * systErr_rel_bTag_dn_eleVec[i-1];

    const double pred_local_from_ele_syst_scaleUnc_up_err = data_ele * fin_TF_local_to_ele * systErr_rel_scaleUnc_up_eleVec[i-1];
    const double pred_local_from_ele_syst_scaleUnc_dn_err = data_ele * fin_TF_local_to_ele * systErr_rel_scaleUnc_dn_eleVec[i-1];

    const double pred_local_from_ele_syst_pdfUnc_up_err = data_ele * fin_TF_local_to_ele * systErr_rel_pdfUnc_up_eleVec[i-1];
    const double pred_local_from_ele_syst_pdfUnc_dn_err = data_ele * fin_TF_local_to_ele * systErr_rel_pdfUnc_dn_eleVec[i-1];

    const double pred_local_from_ele_syst_metMag_up_err = data_ele * fin_TF_local_to_ele * systErr_rel_metMag_up_eleVec[i-1];
    const double pred_local_from_ele_syst_metMag_dn_err = data_ele * fin_TF_local_to_ele * systErr_rel_metMag_dn_eleVec[i-1];

    const double pred_local_from_ele_syst_metPhi_up_err = data_ele * fin_TF_local_to_ele * systErr_rel_metPhi_up_eleVec[i-1];
    const double pred_local_from_ele_syst_metPhi_dn_err = data_ele * fin_TF_local_to_ele * systErr_rel_metPhi_dn_eleVec[i-1];

    const double pred_local_from_ele_syst_jec_up_err = data_ele * fin_TF_local_to_ele * systErr_rel_jec_up_eleVec[i-1];
    const double pred_local_from_ele_syst_jec_dn_err = data_ele * fin_TF_local_to_ele * systErr_rel_jec_dn_eleVec[i-1];

    const double pred_local_from_ele_syst_SF_up_err = data_ele * std::abs(fin_TF_local_to_ele_SFup - fin_TF_local_to_ele);
    const double pred_local_from_ele_syst_SF_dn_err = data_ele * std::abs(fin_TF_local_to_ele_SFdn - fin_TF_local_to_ele);

// avg
    const double pred_local_avg = sel_CS_str == "comb" ? 0.5*pred_local_from_mu + 0.5*pred_local_from_ele : sel_CS_str == "mu" ? pred_local_from_mu : sel_CS_str == "ele" ? pred_local_from_ele : 0;
    const double pred_local_avg_syst = sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst*pred_local_from_mu_syst + pred_local_from_ele_syst*pred_local_from_ele_syst)
                                     : sel_CS_str == "mu" ? pred_local_from_mu_syst : sel_CS_str == "ele" ? pred_local_from_ele_syst : 0;
    const double pred_local_avg_stat_up_err = sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_stat_up_err*pred_local_from_mu_stat_up_err + pred_local_from_ele_stat_up_err*pred_local_from_ele_stat_up_err)
                                            : sel_CS_str == "mu" ? pred_local_from_mu_stat_up_err : sel_CS_str == "ele" ? pred_local_from_ele_stat_up_err : 0;
    const double pred_local_avg_stat_dn_err = sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_stat_dn_err*pred_local_from_mu_stat_dn_err + pred_local_from_ele_stat_dn_err*pred_local_from_ele_stat_dn_err)
                                            : sel_CS_str == "mu" ? pred_local_from_mu_stat_dn_err : sel_CS_str == "ele" ? pred_local_from_ele_stat_dn_err : 0;
    const double pred_local_avg_syst_ISR_up_err = 
                              sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_ISR_up_err*pred_local_from_mu_syst_ISR_up_err + pred_local_from_ele_syst_ISR_up_err*pred_local_from_ele_syst_ISR_up_err)
                              : sel_CS_str == "mu" ? pred_local_from_mu_syst_ISR_up_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_ISR_up_err : 0;
    const double pred_local_avg_syst_ISR_dn_err = 
                              sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_ISR_dn_err*pred_local_from_mu_syst_ISR_dn_err + pred_local_from_ele_syst_ISR_dn_err*pred_local_from_ele_syst_ISR_dn_err)
                              : sel_CS_str == "mu" ? pred_local_from_mu_syst_ISR_dn_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_ISR_dn_err : 0;

    const double pred_local_avg_syst_bTag_up_err = 
                              sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_bTag_up_err*pred_local_from_mu_syst_bTag_up_err + pred_local_from_ele_syst_bTag_up_err*pred_local_from_ele_syst_bTag_up_err)
                              : sel_CS_str == "mu" ? pred_local_from_mu_syst_bTag_up_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_bTag_up_err : 0;
    const double pred_local_avg_syst_bTag_dn_err = 
                              sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_bTag_dn_err*pred_local_from_mu_syst_bTag_dn_err + pred_local_from_ele_syst_bTag_dn_err*pred_local_from_ele_syst_bTag_dn_err)
                              : sel_CS_str == "mu" ? pred_local_from_mu_stat_dn_err : sel_CS_str == "ele" ? pred_local_from_ele_stat_dn_err : 0;

    const double pred_local_avg_syst_scaleUnc_up_err = 
              sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_scaleUnc_up_err*pred_local_from_mu_syst_scaleUnc_up_err + pred_local_from_ele_syst_scaleUnc_up_err*pred_local_from_ele_syst_scaleUnc_up_err)
              : sel_CS_str == "mu" ? pred_local_from_mu_syst_scaleUnc_up_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_scaleUnc_up_err : 0;
    const double pred_local_avg_syst_scaleUnc_dn_err = 
              sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_scaleUnc_dn_err*pred_local_from_mu_syst_scaleUnc_dn_err + pred_local_from_ele_syst_scaleUnc_dn_err*pred_local_from_ele_syst_scaleUnc_dn_err)
              : sel_CS_str == "mu" ? pred_local_from_mu_syst_scaleUnc_dn_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_scaleUnc_dn_err : 0;

    const double pred_local_avg_syst_pdfUnc_up_err = 
                      sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_pdfUnc_up_err*pred_local_from_mu_syst_pdfUnc_up_err + pred_local_from_ele_syst_pdfUnc_up_err*pred_local_from_ele_syst_pdfUnc_up_err)
                      : sel_CS_str == "mu" ? pred_local_from_mu_syst_pdfUnc_up_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_pdfUnc_up_err : 0;
    const double pred_local_avg_syst_pdfUnc_dn_err = 
                      sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_pdfUnc_dn_err*pred_local_from_mu_syst_pdfUnc_dn_err + pred_local_from_ele_syst_pdfUnc_dn_err*pred_local_from_ele_syst_pdfUnc_dn_err)
                      : sel_CS_str == "mu" ? pred_local_from_mu_syst_pdfUnc_dn_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_pdfUnc_dn_err : 0;

    const double pred_local_avg_syst_metMag_up_err = 
                     sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_metMag_up_err*pred_local_from_mu_syst_metMag_up_err + pred_local_from_ele_syst_metMag_up_err*pred_local_from_ele_syst_metMag_up_err)
                     : sel_CS_str == "mu" ? pred_local_from_mu_syst_metMag_up_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_metMag_up_err : 0;
    const double pred_local_avg_syst_metMag_dn_err = 
                     sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_metMag_dn_err*pred_local_from_mu_syst_metMag_dn_err + pred_local_from_ele_syst_metMag_dn_err*pred_local_from_ele_syst_metMag_dn_err)
                     : sel_CS_str == "mu" ? pred_local_from_mu_syst_metMag_dn_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_metMag_dn_err : 0;

    const double pred_local_avg_syst_metPhi_up_err = 
                     sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_metPhi_up_err*pred_local_from_mu_syst_metPhi_up_err + pred_local_from_ele_syst_metPhi_up_err*pred_local_from_ele_syst_metPhi_up_err)
                     : sel_CS_str == "mu" ? pred_local_from_mu_syst_metPhi_up_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_metPhi_up_err : 0;
    const double pred_local_avg_syst_metPhi_dn_err = 
                     sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_metPhi_dn_err*pred_local_from_mu_syst_metPhi_dn_err + pred_local_from_ele_syst_metPhi_dn_err*pred_local_from_ele_syst_metPhi_dn_err)
                     : sel_CS_str == "mu" ? pred_local_from_mu_syst_metPhi_dn_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_metPhi_dn_err : 0;

    const double pred_local_avg_syst_jec_up_err = 
                     sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_jec_up_err*pred_local_from_mu_syst_jec_up_err + pred_local_from_ele_syst_jec_up_err*pred_local_from_ele_syst_jec_up_err)
                     : sel_CS_str == "mu" ? pred_local_from_mu_syst_jec_up_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_jec_up_err : 0;
    const double pred_local_avg_syst_jec_dn_err = 
                     sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_jec_dn_err*pred_local_from_mu_syst_jec_dn_err + pred_local_from_ele_syst_jec_dn_err*pred_local_from_ele_syst_jec_dn_err)
                     : sel_CS_str == "mu" ? pred_local_from_mu_syst_jec_dn_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_jec_dn_err : 0;

    const double pred_local_avg_syst_SF_up_err = 
                     sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_SF_up_err*pred_local_from_mu_syst_SF_up_err + pred_local_from_ele_syst_SF_up_err*pred_local_from_ele_syst_SF_up_err)
                     : sel_CS_str == "mu" ? pred_local_from_mu_syst_SF_up_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_SF_up_err : 0;
    const double pred_local_avg_syst_SF_dn_err = 
                     sel_CS_str == "comb" ? 0.5*sqrt(pred_local_from_mu_syst_SF_dn_err*pred_local_from_mu_syst_SF_dn_err + pred_local_from_ele_syst_SF_dn_err*pred_local_from_ele_syst_SF_dn_err)
                     : sel_CS_str == "mu" ? pred_local_from_mu_syst_SF_dn_err : sel_CS_str == "ele" ? pred_local_from_ele_syst_SF_dn_err : 0;

    const double pred_local_avg_rel_stat_up_err = pred_local_avg == 0 ? 0 : pred_local_avg_stat_up_err/pred_local_avg;
    const double pred_local_avg_rel_stat_dn_err = pred_local_avg == 0 ? 0 : pred_local_avg_stat_dn_err/pred_local_avg;
// if pred_local_avg = 0 then force data_mu =1 and data_ele =1 to calculate the relative syst unc.
    if( relSys_for_ZERO )
    {
      const double enable_mu = (sel_CS_str == "mu" || sel_CS_str == "comb") ? 1.0 : 0.0;
      const double enable_ele = (sel_CS_str == "ele" || sel_CS_str == "comb") ? 1.0 : 0.0;
      if( enable_mu == 0 && enable_ele ==0 ){ std::cout<<"both enable_mu and enable_ele is ZERO?"<<std::endl; return; }
      const double pred_local_avg_rel_syst = pred_local_avg == 0 ? sqrt(enable_mu*finErr_TF_local_to_mu*enable_mu*finErr_TF_local_to_mu + enable_ele*finErr_TF_local_to_ele*enable_ele*finErr_TF_local_to_ele)/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst/pred_local_avg;
      const double pred_local_avg_rel_syst_ISR_up_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_ISR_up_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_ISR_up_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_ISR_up_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_ISR_up_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_ISR_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_ISR_dn_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_ISR_dn_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_ISR_dn_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_ISR_dn_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_ISR_dn_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_ISR_dn_err/pred_local_avg;
      const double pred_local_avg_rel_syst_bTag_up_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_bTag_up_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_bTag_up_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_bTag_up_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_bTag_up_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_bTag_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_bTag_dn_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_bTag_dn_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_bTag_dn_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_bTag_dn_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_bTag_dn_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_bTag_dn_err/pred_local_avg;
  
      const double pred_local_avg_rel_syst_scaleUnc_up_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_scaleUnc_up_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_scaleUnc_up_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_scaleUnc_up_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_scaleUnc_up_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_scaleUnc_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_scaleUnc_dn_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_scaleUnc_dn_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_scaleUnc_dn_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_scaleUnc_dn_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_scaleUnc_dn_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_scaleUnc_dn_err/pred_local_avg;
  
      const double pred_local_avg_rel_syst_pdfUnc_up_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_pdfUnc_up_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_pdfUnc_up_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_pdfUnc_up_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_pdfUnc_up_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_pdfUnc_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_pdfUnc_dn_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_pdfUnc_dn_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_pdfUnc_dn_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_pdfUnc_dn_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_pdfUnc_dn_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_pdfUnc_dn_err/pred_local_avg;
  
      const double pred_local_avg_rel_syst_metMag_up_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_metMag_up_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_metMag_up_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_metMag_up_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_metMag_up_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_metMag_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_metMag_dn_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_metMag_dn_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_metMag_dn_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_metMag_dn_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_metMag_dn_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_metMag_dn_err/pred_local_avg;
  
      const double pred_local_avg_rel_syst_metPhi_up_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_metPhi_up_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_metPhi_up_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_metPhi_up_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_metPhi_up_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_metPhi_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_metPhi_dn_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_metPhi_dn_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_metPhi_dn_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_metPhi_dn_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_metPhi_dn_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_metPhi_dn_err/pred_local_avg;
  
      const double pred_local_avg_rel_syst_jec_up_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_jec_up_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_jec_up_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_jec_up_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_jec_up_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_jec_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_jec_dn_err = pred_local_avg == 0 ? sqrt(enable_mu*fin_TF_local_to_mu * systErr_rel_jec_dn_muVec[i-1]*enable_mu*fin_TF_local_to_mu * systErr_rel_jec_dn_muVec[i-1] + enable_ele*fin_TF_local_to_ele * systErr_rel_jec_dn_eleVec[i-1]*enable_ele*fin_TF_local_to_ele * systErr_rel_jec_dn_eleVec[i-1])/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_jec_dn_err/pred_local_avg;
  
      const double pred_local_avg_rel_syst_SF_up_err = pred_local_avg == 0 ? sqrt( std::abs(enable_mu*fin_TF_local_to_mu_SFup - enable_mu*fin_TF_local_to_mu)*std::abs(enable_mu*fin_TF_local_to_mu_SFup - enable_mu*fin_TF_local_to_mu) + std::abs(enable_ele*fin_TF_local_to_ele_SFup - enable_ele*fin_TF_local_to_ele) * std::abs(enable_ele*fin_TF_local_to_ele_SFup - enable_ele*fin_TF_local_to_ele) )/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_SF_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_SF_dn_err = pred_local_avg == 0 ? sqrt( std::abs(enable_mu*fin_TF_local_to_mu_SFdn - enable_mu*fin_TF_local_to_mu)*std::abs(enable_mu*fin_TF_local_to_mu_SFdn - enable_mu*fin_TF_local_to_mu) + std::abs(enable_ele*fin_TF_local_to_ele_SFdn - enable_ele*fin_TF_local_to_ele) * std::abs(enable_ele*fin_TF_local_to_ele_SFdn - enable_ele*fin_TF_local_to_ele) )/(enable_mu*fin_TF_local_to_mu + enable_ele*fin_TF_local_to_ele) : pred_local_avg_syst_SF_dn_err/pred_local_avg;

      cached_rateVec.push_back(pred_local_avg); cached_stat_upErrVec.push_back(pred_local_avg_rel_stat_up_err); cached_stat_dnErrVec.push_back(pred_local_avg_rel_stat_dn_err);
      cached_stat_abs_upErrVec.push_back(pred_local_avg_stat_up_err); cached_stat_abs_dnErrVec.push_back(pred_local_avg_stat_dn_err);
      cached_systErr_TF_statVec.push_back(pred_local_avg_rel_syst);
      cached_systErr_ISR_upErrVec.push_back(pred_local_avg_rel_syst_ISR_up_err); cached_systErr_ISR_dnErrVec.push_back(pred_local_avg_rel_syst_ISR_dn_err);
      cached_systErr_bTag_upErrVec.push_back(pred_local_avg_rel_syst_bTag_up_err); cached_systErr_bTag_dnErrVec.push_back(pred_local_avg_rel_syst_bTag_dn_err);
      cached_systErr_scaleUnc_upErrVec.push_back(pred_local_avg_rel_syst_scaleUnc_up_err); cached_systErr_scaleUnc_dnErrVec.push_back(pred_local_avg_rel_syst_scaleUnc_dn_err);
      cached_systErr_pdfUnc_upErrVec.push_back(pred_local_avg_rel_syst_pdfUnc_up_err); cached_systErr_pdfUnc_dnErrVec.push_back(pred_local_avg_rel_syst_pdfUnc_dn_err);
      cached_systErr_metMag_upErrVec.push_back(pred_local_avg_rel_syst_metMag_up_err); cached_systErr_metMag_dnErrVec.push_back(pred_local_avg_rel_syst_metMag_dn_err);
      cached_systErr_metPhi_upErrVec.push_back(pred_local_avg_rel_syst_metPhi_up_err); cached_systErr_metPhi_dnErrVec.push_back(pred_local_avg_rel_syst_metPhi_dn_err);
      cached_systErr_jec_upErrVec.push_back(pred_local_avg_rel_syst_jec_up_err); cached_systErr_jec_dnErrVec.push_back(pred_local_avg_rel_syst_jec_dn_err);
      cached_systErr_SF_upErrVec.push_back(pred_local_avg_rel_syst_SF_up_err); cached_systErr_SF_dnErrVec.push_back(pred_local_avg_rel_syst_SF_dn_err);
      dummyVec.push_back(0.0);
    }else
    {
      const double pred_local_avg_rel_syst = pred_local_avg == 0 ? 0 : pred_local_avg_syst/pred_local_avg;
      const double pred_local_avg_rel_syst_ISR_up_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_ISR_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_ISR_dn_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_ISR_dn_err/pred_local_avg;
      const double pred_local_avg_rel_syst_bTag_up_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_bTag_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_bTag_dn_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_bTag_dn_err/pred_local_avg;
  
      const double pred_local_avg_rel_syst_scaleUnc_up_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_scaleUnc_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_scaleUnc_dn_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_scaleUnc_dn_err/pred_local_avg;
  
      const double pred_local_avg_rel_syst_pdfUnc_up_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_pdfUnc_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_pdfUnc_dn_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_pdfUnc_dn_err/pred_local_avg;
  
      const double pred_local_avg_rel_syst_metMag_up_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_metMag_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_metMag_dn_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_metMag_dn_err/pred_local_avg;
  
      const double pred_local_avg_rel_syst_metPhi_up_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_metPhi_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_metPhi_dn_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_metPhi_dn_err/pred_local_avg;
  
      const double pred_local_avg_rel_syst_jec_up_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_jec_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_jec_dn_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_jec_dn_err/pred_local_avg;
  
      const double pred_local_avg_rel_syst_SF_up_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_SF_up_err/pred_local_avg;
      const double pred_local_avg_rel_syst_SF_dn_err = pred_local_avg == 0 ? 0 : pred_local_avg_syst_SF_dn_err/pred_local_avg;
  
      cached_rateVec.push_back(pred_local_avg); cached_stat_upErrVec.push_back(pred_local_avg_rel_stat_up_err); cached_stat_dnErrVec.push_back(pred_local_avg_rel_stat_dn_err);
      cached_stat_abs_upErrVec.push_back(pred_local_avg_stat_up_err); cached_stat_abs_dnErrVec.push_back(pred_local_avg_stat_dn_err);
      cached_systErr_TF_statVec.push_back(pred_local_avg_rel_syst);
      cached_systErr_ISR_upErrVec.push_back(pred_local_avg_rel_syst_ISR_up_err); cached_systErr_ISR_dnErrVec.push_back(pred_local_avg_rel_syst_ISR_dn_err);
      cached_systErr_bTag_upErrVec.push_back(pred_local_avg_rel_syst_bTag_up_err); cached_systErr_bTag_dnErrVec.push_back(pred_local_avg_rel_syst_bTag_dn_err);
      cached_systErr_scaleUnc_upErrVec.push_back(pred_local_avg_rel_syst_scaleUnc_up_err); cached_systErr_scaleUnc_dnErrVec.push_back(pred_local_avg_rel_syst_scaleUnc_dn_err);
      cached_systErr_pdfUnc_upErrVec.push_back(pred_local_avg_rel_syst_pdfUnc_up_err); cached_systErr_pdfUnc_dnErrVec.push_back(pred_local_avg_rel_syst_pdfUnc_dn_err);
      cached_systErr_metMag_upErrVec.push_back(pred_local_avg_rel_syst_metMag_up_err); cached_systErr_metMag_dnErrVec.push_back(pred_local_avg_rel_syst_metMag_dn_err);
      cached_systErr_metPhi_upErrVec.push_back(pred_local_avg_rel_syst_metPhi_up_err); cached_systErr_metPhi_dnErrVec.push_back(pred_local_avg_rel_syst_metPhi_dn_err);
      cached_systErr_jec_upErrVec.push_back(pred_local_avg_rel_syst_jec_up_err); cached_systErr_jec_dnErrVec.push_back(pred_local_avg_rel_syst_jec_dn_err);
      cached_systErr_SF_upErrVec.push_back(pred_local_avg_rel_syst_SF_up_err); cached_systErr_SF_dnErrVec.push_back(pred_local_avg_rel_syst_SF_dn_err);
      dummyVec.push_back(0.0);
    } 
    std::cout<<i-1<<"\t"<<std::setw(6)<<data_mu<<"\t"<<std::setw(6)<<data_ele
             <<"\t\t"<<std::setw(6)<<fin_TF_local_to_mu<<" +- "<<std::setw(6)<<finErr_TF_local_to_mu
             <<"\t"<<std::setw(6)<<fin_TF_local_to_ele<<" +- "<<std::setw(6)<<finErr_TF_local_to_ele
             <<"\t"<<std::setw(6)<<pred_local_from_mu<<" + "<<std::setw(6)<<pred_local_from_mu_stat_up_err<<" - "<<std::setw(6)<<pred_local_from_mu_stat_dn_err<<" +- "<<std::setw(6)<<pred_local_from_mu_syst
             <<"\t"<<std::setw(6)<<pred_local_from_ele<<" + "<<std::setw(6)<<pred_local_from_ele_stat_up_err<<" - "<<std::setw(6)<<pred_local_from_ele_stat_dn_err<<" +- "<<std::setw(6)<<pred_local_from_ele_syst
             <<"\t"<<std::setw(6)<<pred_local_avg<<" + "<<std::setw(6)<<pred_local_avg_stat_up_err<<" - "<<std::setw(6)<<pred_local_avg_stat_dn_err<<" +- "<<std::setw(6)<<pred_local_avg_syst
             <<std::endl;
    sum_pred_local_avg += pred_local_avg;
    sum_pred_local_from_mu += pred_local_from_mu;
    sum_pred_local_from_ele += pred_local_from_ele;
  }

  std::cout<<"\nsum_pred_local_from_mu : "<<sum_pred_local_from_mu<<"  sum_pred_local_from_ele : "<<sum_pred_local_from_ele<<"  sum_pred_local_avg : "<<sum_pred_local_avg<<std::endl<<std::endl;

  TString vecIdStr, divStr, cloPre, cloAft;
  int setWnum;
  if ( doVector_dataCard ){
     vecIdStr = "vector<double>";
     divStr = ",";
     cloPre = "{";
     cloAft = "};";
     setWnum = 0;
  }else{
     vecIdStr = "";
     divStr = " ";
     cloPre = "";
     cloAft = "";
     setWnum = 7;
  }

  std::ofstream of_local(key+"_"+sel_CS_str+".txt", std::ofstream::out);
  char tmpstr[200];
  of_local<< std::setprecision(3)<<std::fixed<<std::left;
  const int nTotChn = hTF_local_mu->GetNbinsX();
  of_local<<"luminosity = "<<data_lumi<<std::endl;
  of_local<<"channels = "<<nTotChn<<std::endl;
  of_local<<"sample = "<<key.c_str()<<"\n"<<std::endl;
  of_local<<std::setw(21)<<"channel = "; for(int i=0; i<nTotChn; i++){ sprintf(tmpstr, "bin%d", i+1); of_local<<std::setw(7)<<tmpstr<<" "; } of_local<<"\n"<<std::endl;
  of_local<<std::setw(21)<<vecIdStr+" rate    = "; of_local<<cloPre; for(auto it : cached_rateVec){ if(it != cached_rateVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it;} of_local<<cloAft; of_local<<"\n"<<std::endl;

  of_local<<"#"<<std::setw(21)<<vecIdStr+" stat_unc_abs_up = "; of_local<<cloPre; for(auto it : cached_stat_abs_upErrVec){ if(it != cached_stat_abs_upErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;
  of_local<<"#"<<std::setw(21)<<vecIdStr+" stat_unc_abs_dn = "; of_local<<cloPre; for(auto it : cached_stat_abs_dnErrVec){ if(it != cached_stat_abs_dnErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;

  of_local<<std::setw(21)<<vecIdStr+" stat_unc_up = "; of_local<<cloPre; for(auto it : cached_stat_upErrVec){ if(it != cached_stat_upErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;
  of_local<<std::setw(21)<<vecIdStr+" stat_unc_dn = "; of_local<<cloPre; for(auto it : cached_stat_dnErrVec){ if(it != cached_stat_dnErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<"\n"<<std::endl;

  of_local<<std::setw(21)<<vecIdStr+" syst_unc_TF_stat_up = "; of_local<<cloPre; for(auto it : cached_systErr_TF_statVec){ if(it != cached_systErr_TF_statVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;
  of_local<<std::setw(21)<<vecIdStr+" syst_unc_TF_stat_dn = "; of_local<<cloPre; for(auto it : cached_systErr_TF_statVec){ if(it != cached_systErr_TF_statVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<"\n"<<std::endl;

  of_local<<std::setw(21)<<vecIdStr+" syst_unc_ISR_up = "; of_local<<cloPre; for(auto it : cached_systErr_ISR_upErrVec){ if(it != cached_systErr_ISR_upErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;
  of_local<<std::setw(21)<<vecIdStr+" syst_unc_ISR_dn = "; of_local<<cloPre; for(auto it : cached_systErr_ISR_dnErrVec){ if(it != cached_systErr_ISR_dnErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<"\n"<<std::endl;

  of_local<<std::setw(21)<<vecIdStr+" syst_unc_bTag_up = "; of_local<<cloPre; for(auto it : cached_systErr_bTag_upErrVec){ if(it != cached_systErr_bTag_upErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;
  of_local<<std::setw(21)<<vecIdStr+" syst_unc_bTag_dn = "; of_local<<cloPre; for(auto it : cached_systErr_bTag_dnErrVec){ if(it != cached_systErr_bTag_dnErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<"\n"<<std::endl;

//  of_local<<std::setw(21)<<"syst_unc_scaleUnc_up = "; for(auto it: cached_systErr_scaleUnc_upErrVec){ of_local<<std::setw(7)<<it<<" "; } of_local<<std::endl;
//  of_local<<std::setw(21)<<"syst_unc_scaleUnc_dn = "; for(auto it: cached_systErr_scaleUnc_dnErrVec){ of_local<<std::setw(7)<<it<<" "; } of_local<<"\n"<<std::endl;

  of_local<<std::setw(21)<<vecIdStr+" syst_unc_pdfUnc_up = "; of_local<<cloPre; for(auto it : cached_systErr_pdfUnc_upErrVec){ if(it != cached_systErr_pdfUnc_upErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;
  of_local<<std::setw(21)<<vecIdStr+" syst_unc_pdfUnc_dn = "; of_local<<cloPre; for(auto it : cached_systErr_pdfUnc_dnErrVec){ if(it != cached_systErr_pdfUnc_dnErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<"\n"<<std::endl;

  of_local<<std::setw(21)<<vecIdStr+" syst_unc_metMag_up = "; of_local<<cloPre; for(auto it : cached_systErr_metMag_upErrVec){ if(it != cached_systErr_metMag_upErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;
  of_local<<std::setw(21)<<vecIdStr+" syst_unc_metMag_dn = "; of_local<<cloPre; for(auto it : cached_systErr_metMag_dnErrVec){ if(it != cached_systErr_metMag_dnErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<"\n"<<std::endl;

  of_local<<std::setw(21)<<vecIdStr+" syst_unc_metPhi_up = "; of_local<<cloPre; for(auto it : cached_systErr_metPhi_upErrVec){ if(it != cached_systErr_metPhi_upErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;
  of_local<<std::setw(21)<<vecIdStr+" syst_unc_metPhi_dn = "; of_local<<cloPre; for(auto it : cached_systErr_metPhi_dnErrVec){ if(it != cached_systErr_metPhi_dnErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<"\n"<<std::endl;

  of_local<<std::setw(21)<<vecIdStr+" syst_unc_jec_up = "; of_local<<cloPre; for(auto it : cached_systErr_jec_upErrVec){ if(it != cached_systErr_jec_upErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;
  of_local<<std::setw(21)<<vecIdStr+" syst_unc_jec_dn = "; of_local<<cloPre; for(auto it : cached_systErr_jec_dnErrVec){ if(it != cached_systErr_jec_dnErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<"\n"<<std::endl;

  of_local<<std::setw(21)<<vecIdStr+" syst_unc_SF_up = "; of_local<<cloPre; for(auto it : cached_systErr_SF_upErrVec){ if(it != cached_systErr_SF_upErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;
  of_local<<std::setw(21)<<vecIdStr+" syst_unc_SF_dn = "; of_local<<cloPre; for(auto it : cached_systErr_SF_dnErrVec){ if(it != cached_systErr_SF_dnErrVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<"\n"<<std::endl;

  of_local<<"## note : for special usage "<<std::endl;
  of_local<<"#"<<std::setw(21)<<vecIdStr+" fin_TF_to_mu = "; of_local<<cloPre; for(auto it : cached_fin_TF_muVec){ if(it != cached_fin_TF_muVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;
  of_local<<"#"<<std::setw(21)<<vecIdStr+" fin_TFerr_to_mu = "; of_local<<cloPre; for(auto it : cached_fin_TFerr_muVec){ if(it != cached_fin_TFerr_muVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<"\n"<<std::endl;

  of_local<<"#"<<std::setw(21)<<vecIdStr+" fin_TF_to_ele = "; of_local<<cloPre; for(auto it : cached_fin_TF_eleVec){ if(it != cached_fin_TF_eleVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<std::endl;
  of_local<<"#"<<std::setw(21)<<vecIdStr+" fin_TFerr_to_ele = "; of_local<<cloPre; for(auto it : cached_fin_TFerr_eleVec){ if(it != cached_fin_TFerr_eleVec.back()) of_local<<std::setw(setWnum)<<it<<divStr+" "; else of_local<<std::setw(setWnum)<<it; } of_local<<cloAft; of_local<<"\n"<<std::endl;

  of_local.close();

//  for(unsigned int ic=0; ic<cached_systErr_pdfUnc_upErrVec.size(); ic++){
//     std::cout<<"ic : "<<ic<<"  pdfUnc_up : "<<cached_systErr_pdfUnc_upErrVec[ic]<<"  pdfUnc_dn : "<<cached_systErr_pdfUnc_dnErrVec[ic]<<std::endl;
//  }
}

void predLLTry()
{
  std::cout<< std::setprecision(3)<<std::fixed;
  std::cout<<"\nLL Predictions"<<std::endl;
  std::cout<<"Bin \t mu CS \t ele CS \t TF(LL/CS)@mu \t TF(LL/CS)@ele \t\t Pred(LL)@mu \t\t\t Pred(LL)@ele \t\t\t Pred(LL)@avg"<<std::endl;
// Combine mu and ele together
// TF_mu/2.0 * CS_mu + TF_ele/2.0 * CS_ele
//      A    * CS_mu +       B    * CS_ele
// Syst : sqrt( dA*dA * CS_mu*CS_mu + dB*dB * CS_ele*CS_ele )
// Stat : ?
  for(unsigned int i=1; i<=hTF_tau_mu->GetNbinsX(); i++)
  {
    const double lepSF_ratio_tau = hYields_lepSF_ratio_tau->GetBinContent(i);
    const double lepSF_ratio_LL = hYields_lepSF_ratio_LL->GetBinContent(i);

// tau
    const double cont_TF_tau_to_mu = hTF_tau_mu->GetBinContent(i), err_TF_tau_to_mu = hTF_tau_mu->GetBinError(i);
    const double fin_TF_tau_to_mu = cont_TF_tau_to_mu * lepSF_ratio_tau, finErr_TF_tau_to_mu = err_TF_tau_to_mu * lepSF_ratio_tau;

    const double cont_TF_LL_to_mu = hTF_LL_mu->GetBinContent(i), err_TF_LL_to_mu = hTF_LL_mu->GetBinError(i);
    const double fin_TF_LL_to_mu = cont_TF_LL_to_mu * lepSF_ratio_LL, finErr_TF_LL_to_mu = err_TF_LL_to_mu * lepSF_ratio_LL;
    
    const double data_mu = hYields_Data_mu->GetBinContent(i);

    const double pred_tau_from_mu = data_mu * fin_TF_tau_to_mu; 
    const double pred_tau_from_mu_stat = data_mu ==0 ? 1.8 * fin_TF_tau_to_mu : sqrt(data_mu) * fin_TF_tau_to_mu;
    const double pred_tau_from_mu_syst = data_mu * finErr_TF_tau_to_mu;

    const double pred_LL_from_mu = data_mu * fin_TF_LL_to_mu; 
    const double pred_LL_from_mu_stat = data_mu ==0 ? 1.8 * fin_TF_LL_to_mu : sqrt(data_mu) * fin_TF_LL_to_mu;
    const double pred_LL_from_mu_syst = data_mu * finErr_TF_LL_to_mu;

// LL 
    const double cont_TF_tau_to_ele = hTF_tau_ele->GetBinContent(i), err_TF_tau_to_ele = hTF_tau_ele->GetBinError(i);
    const double fin_TF_tau_to_ele = cont_TF_tau_to_ele * lepSF_ratio_tau, finErr_TF_tau_to_ele = err_TF_tau_to_ele * lepSF_ratio_tau;

    const double cont_TF_LL_to_ele = hTF_LL_ele->GetBinContent(i), err_TF_LL_to_ele = hTF_LL_ele->GetBinError(i);
    const double fin_TF_LL_to_ele = cont_TF_LL_to_ele * lepSF_ratio_LL, finErr_TF_LL_to_ele = err_TF_LL_to_ele * lepSF_ratio_LL;
    
    const double data_ele = hYields_Data_ele->GetBinContent(i);

    const double pred_tau_from_ele = data_ele * fin_TF_tau_to_ele; 
    const double pred_tau_from_ele_stat = data_ele ==0 ? 1.8 * fin_TF_tau_to_ele : sqrt(data_ele) * fin_TF_tau_to_ele;
    const double pred_tau_from_ele_syst = data_ele * finErr_TF_tau_to_ele;

    const double pred_LL_from_ele = data_ele * fin_TF_LL_to_ele; 
    const double pred_LL_from_ele_stat = data_ele ==0 ? 1.8 * fin_TF_LL_to_ele : sqrt(data_ele) * fin_TF_LL_to_ele;
    const double pred_LL_from_ele_syst = data_ele * finErr_TF_LL_to_ele;

// avg
    const double pred_LL_avg = 0.5*pred_LL_from_mu + 0.5*pred_LL_from_ele;
    const double pred_LL_avg_syst = 0.5*sqrt(pred_LL_from_mu_syst*pred_LL_from_mu_syst + pred_LL_from_ele_syst*pred_LL_from_ele_syst);
    const double pred_LL_avg_stat = 0.5*sqrt(pred_LL_from_mu_stat*pred_LL_from_mu_stat + pred_LL_from_ele_stat*pred_LL_from_ele_stat);

    std::cout<<i-1<<"\t"<<std::setw(6)<<data_mu<<"\t"<<std::setw(6)<<data_ele
             <<"\t\t"<<std::setw(6)<<fin_TF_LL_to_mu<<" +- "<<std::setw(6)<<finErr_TF_LL_to_mu
             <<"\t"<<std::setw(6)<<fin_TF_LL_to_ele<<" +- "<<std::setw(6)<<finErr_TF_LL_to_ele
             <<"\t"<<std::setw(6)<<pred_LL_from_mu<<" +- "<<std::setw(6)<<pred_LL_from_mu_stat<<" +- "<<std::setw(6)<<pred_LL_from_mu_syst
             <<"\t"<<std::setw(6)<<pred_LL_from_ele<<" +- "<<std::setw(6)<<pred_LL_from_ele_stat<<" +- "<<std::setw(6)<<pred_LL_from_ele_syst
             <<"\t"<<std::setw(6)<<pred_LL_avg<<" +- "<<std::setw(6)<<pred_LL_avg_stat<<" +- "<<std::setw(6)<<pred_LL_avg_syst
             <<std::endl;
  }
}

void predFromCSonHadtauLL()
{
  std::cout<< std::setprecision(3)<<std::fixed;
  std::cout<<"\nPredictions using Muon CS"<<std::endl;
  std::cout<<"Bin \t Data CS \t TF(Hadtau/CS) \t\t TF(LL/CS) \t\t\t Pred(Hadtau) \t\t\t Pred(LL)"<<std::endl;
  for(unsigned int i=1; i<=hTF_tau_mu->GetNbinsX(); i++)
  {
    const double cont_TF_tau_to_mu = hTF_tau_mu->GetBinContent(i), err_TF_tau_to_mu = hTF_tau_mu->GetBinError(i);
    const double lepSF_ratio_tau = hYields_lepSF_ratio_tau->GetBinContent(i);
    const double fin_TF_tau_to_mu = cont_TF_tau_to_mu * lepSF_ratio_tau, finErr_TF_tau_to_mu = err_TF_tau_to_mu * lepSF_ratio_tau;

    const double cont_TF_LL_to_mu = hTF_LL_mu->GetBinContent(i), err_TF_LL_to_mu = hTF_LL_mu->GetBinError(i);
    const double lepSF_ratio_LL = hYields_lepSF_ratio_LL->GetBinContent(i);
    const double fin_TF_LL_to_mu = cont_TF_LL_to_mu * lepSF_ratio_LL, finErr_TF_LL_to_mu = err_TF_LL_to_mu * lepSF_ratio_LL;

    const double data_CS = hYields_Data_mu->GetBinContent(i);

    const double pred_tau = data_CS * fin_TF_tau_to_mu;
    const double pred_tau_stat = data_CS ==0 ? 1.8 * fin_TF_tau_to_mu : sqrt(data_CS) * fin_TF_tau_to_mu;
    const double pred_tau_syst = data_CS * finErr_TF_tau_to_mu;

    const double pred_LL = data_CS * fin_TF_LL_to_mu;
    const double pred_LL_stat = data_CS ==0 ? 1.8 * fin_TF_LL_to_mu : sqrt(data_CS) * fin_TF_LL_to_mu;
    const double pred_LL_syst = data_CS * finErr_TF_LL_to_mu;

    std::cout<<i-1<<"\t"<<std::setw(6)<<data_CS<<"\t\t"<<std::setw(6)<<fin_TF_tau_to_mu<<" +- "<<std::setw(6)<<finErr_TF_tau_to_mu<<"\t"<<std::setw(6)<<fin_TF_LL_to_mu<<" +- "<<std::setw(6)<<finErr_TF_LL_to_mu<<"\t"<<std::setw(6)<<pred_tau<<" +- "<<std::setw(6)<<pred_tau_stat<<" +- "<<std::setw(6)<<pred_tau_syst<<"\t"<<std::setw(6)<<pred_LL<<" +- "<<std::setw(6)<<pred_LL_stat<<" +- "<<std::setw(6)<<pred_LL_syst<<std::endl;
  }

  std::cout<< std::setprecision(3)<<std::fixed;
  std::cout<<"\nPredictions using Electron CS"<<std::endl;
  std::cout<<"Bin \t Data CS \t TF(Hadtau/CS) \t\t TF(LL/CS) \t\t\t Pred(Hadtau) \t\t\t Pred(LL)"<<std::endl;
  for(unsigned int i=1; i<=hTF_tau_ele->GetNbinsX(); i++)
  {
    const double cont_TF_tau_to_ele = hTF_tau_ele->GetBinContent(i), err_TF_tau_to_ele = hTF_tau_ele->GetBinError(i);
    const double lepSF_ratio_tau = hYields_lepSF_ratio_tau->GetBinContent(i);
    const double fin_TF_tau_to_ele = cont_TF_tau_to_ele * lepSF_ratio_tau, finErr_TF_tau_to_ele = err_TF_tau_to_ele * lepSF_ratio_tau;

    const double cont_TF_LL_to_ele = hTF_LL_ele->GetBinContent(i), err_TF_LL_to_ele = hTF_LL_ele->GetBinError(i);
    const double lepSF_ratio_LL = hYields_lepSF_ratio_LL->GetBinContent(i);
    const double fin_TF_LL_to_ele = cont_TF_LL_to_ele * lepSF_ratio_LL, finErr_TF_LL_to_ele = err_TF_LL_to_ele * lepSF_ratio_LL;

    const double data_CS = hYields_Data_ele->GetBinContent(i);

    const double pred_tau = data_CS * fin_TF_tau_to_ele;
    const double pred_tau_stat = data_CS ==0 ? 1.8 * fin_TF_tau_to_ele : sqrt(data_CS) * fin_TF_tau_to_ele;
    const double pred_tau_syst = data_CS * finErr_TF_tau_to_ele;

    const double pred_LL = data_CS * fin_TF_LL_to_ele;
    const double pred_LL_stat = data_CS ==0 ? 1.8 * fin_TF_LL_to_ele : sqrt(data_CS) * fin_TF_LL_to_ele;
    const double pred_LL_syst = data_CS * finErr_TF_LL_to_ele;

    std::cout<<i-1<<"\t"<<std::setw(6)<<data_CS<<"\t\t"<<std::setw(6)<<fin_TF_tau_to_ele<<" +- "<<std::setw(6)<<finErr_TF_tau_to_ele<<"\t"<<std::setw(6)<<fin_TF_LL_to_ele<<" +- "<<std::setw(6)<<finErr_TF_LL_to_ele<<"\t"<<std::setw(6)<<pred_tau<<" +- "<<std::setw(6)<<pred_tau_stat<<" +- "<<std::setw(6)<<pred_tau_syst<<"\t"<<std::setw(6)<<pred_LL<<" +- "<<std::setw(6)<<pred_LL_stat<<" +- "<<std::setw(6)<<pred_LL_syst<<std::endl;
  }
}

void prtTFfactors()
{
  std::cout << std::setprecision(3) << std::fixed;
  std::cout<<"\nMuon CS"<<std::endl;
  std::cout<<"Bin \t\t Hadtau \t\t LL \t\t\t sum \t\t\t Muon CS \t\t\t TF(Hadtau/CS) \t\t\t TF(LL/CS) \t\t\t TF(sum/CS) \t\t Data"<<std::endl;
  for(unsigned i=1; i<= hTF_tau_mu->GetNbinsX();i++)
  {
    const double pctErr_TF_tau_to_mu = hTF_tau_mu->GetBinContent(i) == 0? 0.0 : hTF_tau_mu->GetBinError(i)/hTF_tau_mu->GetBinContent(i);
    const double pctErr_TF_LL_to_mu = hTF_LL_mu->GetBinContent(i) == 0? 0.0 : hTF_LL_mu->GetBinError(i)/hTF_LL_mu->GetBinContent(i);
    const double pctErr_TF_sum_to_mu = hTF_sum_mu->GetBinContent(i) == 0? 0.0 : hTF_sum_mu->GetBinError(i)/hTF_sum_mu->GetBinContent(i);

    const double lepSF_ratio_tau = hYields_lepSF_ratio_tau->GetBinContent(i);
    const double lepSF_ratio_LL = hYields_lepSF_ratio_LL->GetBinContent(i);
    const double lepSF_ratio_sum = hYields_lepSF_ratio_sum->GetBinContent(i);

    std::cout<<i-1<<"\t"<<std::setw(6)<<hYields_tau->GetBinContent(i)<<" +- "<<std::setw(6)<<hYields_tau->GetBinError(i)<<"\t"<<std::setw(6)<<hYields_LL->GetBinContent(i)<<" +- "<<std::setw(6)<<hYields_LL->GetBinError(i)<<"\t"<<std::setw(6)<<hYields_sum->GetBinContent(i)<<" +- "<<std::setw(6)<<hYields_sum->GetBinError(i)<<"\t"<<std::setw(6)<<hYields_MC_mu->GetBinContent(i)<<" +- "<<std::setw(6)<<hYields_MC_mu->GetBinError(i)<<"\t"<<std::setw(6)<<hTF_tau_mu->GetBinContent(i)<<" +- "<<std::setw(6)<<hTF_tau_mu->GetBinError(i)<<" ("<<pctErr_TF_tau_to_mu*100.0<<"%)"<<" x "<<std::setw(5)<<lepSF_ratio_tau<<"\t"<<std::setw(6)<<hTF_LL_mu->GetBinContent(i)<<" +- "<<std::setw(6)<<hTF_LL_mu->GetBinError(i)<<" ("<<pctErr_TF_LL_to_mu*100.0<<"%)"<<" x "<<std::setw(5)<<lepSF_ratio_LL<<"\t"<<std::setw(6)<<hTF_sum_mu->GetBinContent(i)<<" +- "<<std::setw(6)<<hTF_sum_mu->GetBinError(i)<<" ("<<pctErr_TF_sum_to_mu*100.0<<"%)"<<" x "<<std::setw(5)<<lepSF_ratio_sum<<"\t"<<std::setw(6)<<hYields_Data_mu->GetBinContent(i)<<std::endl;
  }

  std::cout << std::setprecision(3) << std::fixed;
  std::cout<<"\nElectron CS"<<std::endl;
  std::cout<<"Bin \t\t Hadtau \t\t LL \t\t\t sum \t\t\t Electron CS \t\t\t TF(Hadtau/CS) \t\t\t TF(LL/CS) \t\t\t TF(sum/CS) \t\t Data"<<std::endl;
  for(unsigned i=1; i<= hTF_tau_ele->GetNbinsX();i++)
  {
    const double pctErr_TF_tau_to_ele = hTF_tau_ele->GetBinContent(i) == 0? 0.0 : hTF_tau_ele->GetBinError(i)/hTF_tau_ele->GetBinContent(i);
    const double pctErr_TF_LL_to_ele = hTF_LL_ele->GetBinContent(i) == 0? 0.0 : hTF_LL_ele->GetBinError(i)/hTF_LL_ele->GetBinContent(i);
    const double pctErr_TF_sum_to_ele = hTF_sum_ele->GetBinContent(i) == 0? 0.0 : hTF_sum_ele->GetBinError(i)/hTF_sum_ele->GetBinContent(i);

    const double lepSF_ratio_tau = hYields_lepSF_ratio_tau->GetBinContent(i);
    const double lepSF_ratio_LL = hYields_lepSF_ratio_LL->GetBinContent(i);
    const double lepSF_ratio_sum = hYields_lepSF_ratio_sum->GetBinContent(i);

    std::cout<<i-1<<"\t"<<std::setw(6)<<hYields_tau->GetBinContent(i)<<" +- "<<std::setw(6)<<hYields_tau->GetBinError(i)<<"\t"<<std::setw(6)<<hYields_LL->GetBinContent(i)<<" +- "<<std::setw(6)<<hYields_LL->GetBinError(i)<<"\t"<<std::setw(6)<<hYields_sum->GetBinContent(i)<<" +- "<<std::setw(6)<<hYields_sum->GetBinError(i)<<"\t"<<std::setw(6)<<hYields_MC_ele->GetBinContent(i)<<" +- "<<std::setw(6)<<hYields_MC_ele->GetBinError(i)<<"\t"<<std::setw(6)<<hTF_tau_ele->GetBinContent(i)<<" +- "<<std::setw(6)<<hTF_tau_ele->GetBinError(i)<<" ("<<pctErr_TF_tau_to_ele*100.0<<"%)"<<" x "<<std::setw(5)<<lepSF_ratio_tau<<"\t"<<std::setw(6)<<hTF_LL_ele->GetBinContent(i)<<" +- "<<std::setw(6)<<hTF_LL_ele->GetBinError(i)<<" ("<<pctErr_TF_LL_to_ele*100.0<<"%)"<<" x "<<std::setw(5)<<lepSF_ratio_LL<<"\t"<<std::setw(6)<<hTF_sum_ele->GetBinContent(i)<<" +- "<<std::setw(6)<<hTF_sum_ele->GetBinError(i)<<" ("<<pctErr_TF_sum_to_ele*100.0<<"%)"<<" x "<<std::setw(5)<<lepSF_ratio_sum<<"\t"<<std::setw(6)<<hYields_Data_ele->GetBinContent(i)<<std::endl;
  }
}

