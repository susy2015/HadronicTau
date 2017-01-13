const double mc_lumi = 8000.0;
const double data_lumi = 36352.970569733;

const double scale_mc = data_lumi/mc_lumi;

void makePred()
{
  TFile *file_Data_CS = new TFile("Data_MET_CS.root");
  TFile *file_MC_CS = new TFile("Mix_CS.root");
  TFile *file_HadTauLL = new TFile("Mix_HadTauLL.root");

  TH1* hYields_MC_mu = (TH1D*) file_MC_CS->Get("hYields_mu");
  TH1* hYields_MC_ele = (TH1D*) file_MC_CS->Get("hYields_el");

  TH1* hYields_tau = (TH1D*) file_HadTauLL->Get("hYields_tau");
  TH1* hYields_LL = (TH1D*) file_HadTauLL->Get("hYields_LL");
  TH1* hYields_sum = static_cast<TH1*>(hYields_tau->Clone("hYields_sum")); hYields_sum->Add(hYields_LL);

  TH1* hYields_Veto_tau = (TH1D*) file_HadTauLL->Get("hYields_Veto_tau");
  TH1* hYields_Veto_tau_SF = (TH1D*) file_HadTauLL->Get("hYields_Veto_tau_SF");
  TH1* hYields_Pass_tau = (TH1D*) file_HadTauLL->Get("hYields_Pass_tau");

  TH1* hYields_Veto_LL = (TH1D*) file_HadTauLL->Get("hYields_Veto_LL");
  TH1* hYields_Veto_LL_SF = (TH1D*) file_HadTauLL->Get("hYields_Veto_LL_SF");
  TH1* hYields_Pass_LL = (TH1D*) file_HadTauLL->Get("hYields_Pass_LL");

  TH1* hYields_Data_mu = (TH1D*) file_Data_CS->Get("hYields_mu");
  TH1* hYields_Data_ele = (TH1D*) file_Data_CS->Get("hYields_el");

  hYields_MC_mu->Scale(scale_mc);
  hYields_MC_ele->Scale(scale_mc);

  hYields_tau->Scale(scale_mc);
  hYields_LL->Scale(scale_mc);
  hYields_sum->Scale(scale_mc);

  hYields_Veto_tau->Scale(scale_mc);
  hYields_Veto_tau_SF->Scale(scale_mc);
  hYields_Pass_tau->Scale(scale_mc);

  hYields_Veto_LL->Scale(scale_mc);
  hYields_Veto_LL_SF->Scale(scale_mc);
  hYields_Pass_LL->Scale(scale_mc);

// (pass+veto-veto_SF)/pass --> to be directly applied on the TF
  TH1* hYields_lepSF_ratio_tau = static_cast<TH1*>(hYields_Pass_tau->Clone("hYields_lepSF_ratio_tau"));
  hYields_lepSF_ratio_tau->Add(hYields_Veto_tau);
  hYields_lepSF_ratio_tau->Add(hYields_Veto_tau_SF, -1);
  hYields_lepSF_ratio_tau->Divide(hYields_Pass_tau);

  TH1* hYields_lepSF_ratio_LL = static_cast<TH1*>(hYields_Pass_LL->Clone("hYields_lepSF_ratio_LL"));
  hYields_lepSF_ratio_LL->Add(hYields_Veto_LL);
  hYields_lepSF_ratio_LL->Add(hYields_Veto_LL_SF, -1);
  hYields_lepSF_ratio_LL->Divide(hYields_Pass_LL);

  TH1* hYields_lepSF_ratio_sum = static_cast<TH1*>(hYields_Pass_tau->Clone("hYields_lepSF_ratio_sum"));
  hYields_lepSF_ratio_sum->Add(hYields_Pass_LL);
  TH1* hYields_lepSF_ratio_sum_tmp = (TH1*) hYields_lepSF_ratio_sum->Clone("hYields_lepSF_ratio_sum_tmp"); // sum_tmp is the sum of Pass_tau and Pass_LL
  hYields_lepSF_ratio_sum->Add(hYields_Veto_tau);
  hYields_lepSF_ratio_sum->Add(hYields_Veto_LL);
  hYields_lepSF_ratio_sum->Add(hYields_Veto_tau_SF, -1);
  hYields_lepSF_ratio_sum->Add(hYields_Veto_LL_SF, -1);
  hYields_lepSF_ratio_sum->Divide(hYields_lepSF_ratio_sum_tmp); 

// Derive the translation factors
  TH1* hTF_tau_mu = static_cast<TH1*>(hYields_tau->Clone("hTF_tau_mu")); hTF_tau_mu->Divide(hYields_MC_mu);
  TH1* hTF_tau_ele = static_cast<TH1*>(hYields_tau->Clone("hTF_tau_ele")); hTF_tau_ele->Divide(hYields_MC_ele);
  TH1* hTF_LL_mu = static_cast<TH1*>(hYields_LL->Clone("hTF_LL_mu")); hTF_LL_mu->Divide(hYields_MC_mu);
  TH1* hTF_LL_ele = static_cast<TH1*>(hYields_LL->Clone("hTF_LL_ele")); hTF_LL_ele->Divide(hYields_MC_ele);
  TH1* hTF_sum_mu = static_cast<TH1*>(hYields_sum->Clone("hTF_sum_mu")); hTF_sum_mu->Divide(hYields_MC_mu);
  TH1* hTF_sum_ele = static_cast<TH1*>(hYields_sum->Clone("hTF_sum_ele")); hTF_sum_ele->Divide(hYields_MC_ele);

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

  std::vector<double> systErr_rel_ISR_up_muVec = {0.0223449, 0.0144005, 0.0153194, 0.0523595, 0.0118566, 0.0191526, 0.0273686, 0.0379992, 0.0157679, 0.0286615, 0.0156814, 0.0165899, 0.0233174, 0.00300184, 0.00584308, 0.0295252, 0.0717616, 0.0340199, 0.123465, 0.0388341, 0.0180426, 0.0177415, 0.0247206, 0.0308921, 0.00367983, 0.0079269, 0.0162445, 0.042758, 0.00660086, 0.0467009, 0.025531, 0.021758, 0.00524293, 0.00275218, 0.0298734, 0.0920234, 0.0678468, 0.0203841, 0.0221972, 0.0153352, 0.00937534, 0.0238859, 0.000558504, 0.0725679, 0.00191255, 0.0191745, 0.014791, 0.0835808, 0.0393838, 0.0187466, 0.0555218, 0.0358415, 0.0201543, 0.0389601, 0.040124, 0.0104884, 0.0365455, 0.00843715, 0.0502504, 0.0257902, 0.0202889, 0.0857592, 0.026684, 0.0560776, 0.0794008, 0.0377195, 0.0557219, 0.0160422, 0.0851351, 0.0161418, 0.05234, 0.0317891, 0.0212015, 0.0296471, 0.00768832, 0.00365471, 0.0134244, 0.127825, 0.0132602, 0.00752821, 0.0309767, 0.111493, 0.0735374, 0.0647696};
  std::vector<double> systErr_rel_ISR_dn_muVec = {0.0289214, 0.0241403, 0.0228128, 0.0720088, 0.0346087, 0.025797, 0.0438907, 0.056684, 0.0266344, 0.0408251, 0.0241304, 0.0295604, 0.0380576, 0.0473047, 0.00559874, 0.0412971, 0.093301, 0.0692335, 0.232483, 0.0597023, 0.0219244, 0.0230148, 0.0340433, 0.0455982, 0.000891312, 0.00515911, 0.0230902, 0.0655584, 0.0179195, 0.0617011, 0.0381435, 0.0296835, 0.00769949, 0.0024411, 0.0326997, 0.119734, 0.0991401, 0.0276846, 0.030658, 0.0220447, 0.0128859, 0.0386158, 0.00134284, 0.121772, 0.000662684, 0.0331461, 0.0257871, 0.148419, 0.0612226, 0.0311553, 0.0787899, 0.0591165, 0.0345857, 0.0653338, 0.0676296, 0.0156884, 0.0658163, 0.00458177, 0.0709808, 0.0381576, 0.0382583, 0.174665, 0.039951, 0.0804977, 0.134029, 0.057846, 0.0931922, 0.00175686, 0.129962, 0.0304459, 0.07831, 0.0505875, 0.0478878, 0.054331, 0.0141949, 0.00960203, 0.026401, 0.237861, 0.0184043, 0.0207878, 0.0597065, 0.212839, 0.13722, 0.130556};
  std::vector<double> systErr_rel_ISR_up_eleVec = {0.0225375, 0.0156323, 0.0373153, 0.0551565, 0.0216432, 0.0210539, 0.0180684, 0.0208353, 0.0213983, 0.0211316, 0.00672371, 0.0222339, 0.0327263, 0.0447201, 0.0051285, 0.0353762, 0.0321294, 0.00585544, 0.0316591, 0.0466683, 0.00666844, 0.0159263, 0.0217629, 0.0522447, 0.0324078, 0.0484707, 0.0198368, 0.0370599, 0.00393814, 0.0714938, 0.0248919, 0.0324219, 0.0775123, 0.0394969, 0.0187723, 0.064946, 0.0150259, 0.0233281, 0.0183106, 0.0122795, 0.00650126, 0.0519377, 0.0208866, 0.0474768, 0.0095227, 0.00670275, 0.0616357, 0.037756, 0.0415645, 0.0379646, 0.0104365, 0.0156736, 0.0237569, 0.0293648, 0.0898318, 0.0109936, 0.0317469, 0.11277, 0.0430595, 0.0298639, 0.0696875, 0.0761212, 0.0192987, 0.0293796, 0.102137, 0.0227234, 0.0826642, 0.0196719, 0.0270137, 0.0297348, 0.0397054, 0.0474895, 0.00394319, 0.0519013, 0.0635913, 0.0312723, 0.0501955, 0.169421, 0.0541806, 0.0337054, 0.0257037, 0.0811609, 0.036813, 0.012029};
  std::vector<double> systErr_rel_ISR_dn_eleVec = {0.0287617, 0.0273989, 0.0571614, 0.0750553, 0.0240476, 0.0288056, 0.0318027, 0.0341299, 0.0333132, 0.0308798, 0.0110599, 0.0359333, 0.0491028, 0.0944098, 0.011577, 0.0491549, 0.0485072, 0.015243, 0.0644796, 0.0648999, 0.0135061, 0.0210496, 0.0314689, 0.0705422, 0.0437992, 0.0773629, 0.0275093, 0.0557447, 0.00617299, 0.102722, 0.0369173, 0.0428255, 0.0952644, 0.0605482, 0.0346014, 0.0993139, 0.0319285, 0.0312642, 0.0251675, 0.0175068, 0.00822826, 0.077911, 0.0337727, 0.088492, 0.0162791, 0.0106498, 0.103137, 0.0573127, 0.0641788, 0.0574306, 0.0218164, 0.0171318, 0.039727, 0.0505629, 0.131371, 0.0179098, 0.0514278, 0.148043, 0.0618133, 0.0434104, 0.105552, 0.108192, 0.0284595, 0.0453019, 0.164564, 0.0326755, 0.145216, 0.0227579, 0.0485468, 0.0502088, 0.0613045, 0.070245, 0.0207897, 0.0895137, 0.101396, 0.0630546, 0.0910642, 0.277629, 0.0998254, 0.0763006, 0.0508281, 0.171823, 0.0768417, 0.0307887};

  std::vector<double> systErr_rel_bTag_up_muVec = {0.000,0.001,0.003,0.005,0.002,0.001,0.002,0.001,0.006,0.001,0.003,0.001,0.007,0.006,0.016,0.002,0.003,0.004,0.016,0.014,0.004,0.002,0.010,0.003,0.010,0.013,0.001,0.011,0.000,0.026,0.006,0.003,0.014,0.041,0.048,0.024,0.015,0.004,0.002,0.002,0.007,0.007,0.010,0.006,0.004,0.005,0.010,0.027,0.000,0.001,0.019,0.012,0.000,0.007,0.001,0.008,0.005,0.014,0.003,0.004,0.001,0.017,0.006,0.006,0.017,0.003,0.013,0.013,0.015,0.003,0.011,0.006,0.010,0.004,0.014,0.021,0.010,0.048,0.035,0.054,0.014,0.005,0.002,0.016};
  std::vector<double> systErr_rel_bTag_dn_muVec = {0.000,0.001,0.003,0.005,0.002,0.001,0.002,0.001,0.007,0.001,0.003,0.002,0.007,0.008,0.016,0.002,0.003,0.004,0.015,0.012,0.002,0.002,0.011,0.004,0.011,0.013,0.001,0.012,0.001,0.028,0.006,0.002,0.014,0.039,0.051,0.026,0.016,0.004,0.002,0.002,0.009,0.007,0.010,0.006,0.004,0.005,0.010,0.029,0.000,0.002,0.020,0.011,0.000,0.006,0.001,0.008,0.005,0.014,0.004,0.005,0.001,0.018,0.007,0.006,0.019,0.004,0.014,0.014,0.017,0.002,0.012,0.006,0.010,0.005,0.015,0.023,0.009,0.052,0.036,0.057,0.014,0.005,0.002,0.015};
  std::vector<double> systErr_rel_bTag_up_eleVec = {0.000, 0.001, 0.001, 0.001, 0.005, 0.001, 0.003, 0.004, 0.007, 0.000, 0.003, 0.007, 0.020, 0.020, 0.008, 0.003, 0.007, 0.006, 0.022, 0.006, 0.000, 0.003, 0.006, 0.010, 0.013, 0.020, 0.001, 0.005, 0.011, 0.027, 0.000, 0.006, 0.014, 0.002, 0.040, 0.016, 0.031, 0.005, 0.003, 0.008, 0.000, 0.007, 0.011, 0.021, 0.010, 0.009, 0.007, 0.033, 0.001, 0.007, 0.014, 0.000, 0.002, 0.003, 0.003, 0.010, 0.009, 0.007, 0.004, 0.010, 0.004, 0.023, 0.004, 0.004, 0.015, 0.000, 0.019, 0.003, 0.025, 0.006, 0.008, 0.018, 0.021, 0.006, 0.003, 0.025, 0.012, 0.046, 0.015, 0.059, 0.004, 0.004, 0.010, 0.008};
  std::vector<double> systErr_rel_bTag_dn_eleVec = {0.000, 0.001, 0.001, 0.001, 0.005, 0.001, 0.003, 0.004, 0.007, 0.000, 0.003, 0.007, 0.021, 0.019, 0.008, 0.003, 0.007, 0.006, 0.021, 0.008, 0.001, 0.003, 0.006, 0.010, 0.013, 0.021, 0.001, 0.006, 0.012, 0.029, 0.000, 0.006, 0.014, 0.002, 0.041, 0.017, 0.033, 0.006, 0.004, 0.009, 0.000, 0.007, 0.011, 0.023, 0.010, 0.009, 0.007, 0.034, 0.001, 0.007, 0.014, 0.001, 0.002, 0.003, 0.002, 0.010, 0.009, 0.007, 0.004, 0.011, 0.004, 0.024, 0.004, 0.004, 0.015, 0.000, 0.019, 0.004, 0.027, 0.007, 0.008, 0.019, 0.021, 0.007, 0.002, 0.027, 0.010, 0.051, 0.015, 0.062, 0.004, 0.006, 0.012, 0.008};

  std::cout<< std::setprecision(3)<<std::fixed;
  std::cout<<"\nHadtau Predictions"<<std::endl;
  std::cout<<"Bin \t mu CS \t ele CS \t TF(Hadtau/CS)@mu \t TF(Hadtau/CS)@ele \t\t Pred(tau)@mu \t\t\t\t Pred(tau)@ele \t\t\t\t Pred(tau)@avg"<<std::endl;
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
    const double data_mu_dn_bound = (data_mu ==0 )? 0. : (ROOT::Math::gamma_quantile(alpha/2, data_mu, 1.0));
    const double data_mu_up_bound = ROOT::Math::gamma_quantile_c(alpha/2, data_mu+1.0, 1.0);
    const double data_mu_dn_err = data_mu - data_mu_dn_bound;
    const double data_mu_up_err = data_mu_up_bound - data_mu;

    const double pred_tau_from_mu = data_mu * fin_TF_tau_to_mu;
    const double pred_tau_from_mu_stat_dn_err = data_mu_dn_err * fin_TF_tau_to_mu;
    const double pred_tau_from_mu_stat_up_err = data_mu_up_err * fin_TF_tau_to_mu;
    const double pred_tau_from_mu_syst = data_mu * finErr_TF_tau_to_mu;
    const double pred_tau_from_mu_syst_ISR_up_err = data_mu * fin_TF_tau_to_mu * systErr_rel_ISR_up_muVec[i-1];
    const double pred_tau_from_mu_syst_ISR_dn_err = data_mu * fin_TF_tau_to_mu * systErr_rel_ISR_dn_muVec[i-1];
    const double pred_tau_from_mu_syst_bTag_up_err = data_mu * fin_TF_tau_to_mu * systErr_rel_bTag_up_muVec[i-1];
    const double pred_tau_from_mu_syst_bTag_dn_err = data_mu * fin_TF_tau_to_mu * systErr_rel_bTag_dn_muVec[i-1];

    const double pred_LL_from_mu = data_mu * fin_TF_LL_to_mu; 
    const double pred_LL_from_mu_stat = data_mu ==0 ? 1.8 * fin_TF_LL_to_mu : sqrt(data_mu) * fin_TF_LL_to_mu;
    const double pred_LL_from_mu_syst = data_mu * finErr_TF_LL_to_mu;

// LL 
    const double cont_TF_tau_to_ele = hTF_tau_ele->GetBinContent(i), err_TF_tau_to_ele = hTF_tau_ele->GetBinError(i);
    const double fin_TF_tau_to_ele = cont_TF_tau_to_ele * lepSF_ratio_tau, finErr_TF_tau_to_ele = err_TF_tau_to_ele * lepSF_ratio_tau;

    const double cont_TF_LL_to_ele = hTF_LL_ele->GetBinContent(i), err_TF_LL_to_ele = hTF_LL_ele->GetBinError(i);
    const double fin_TF_LL_to_ele = cont_TF_LL_to_ele * lepSF_ratio_LL, finErr_TF_LL_to_ele = err_TF_LL_to_ele * lepSF_ratio_LL;
    
    const double data_ele = hYields_Data_ele->GetBinContent(i);
    const double data_ele_dn_bound = (data_ele ==0 )? 0. : (ROOT::Math::gamma_quantile(alpha/2, data_ele, 1.0));
    const double data_ele_up_bound = ROOT::Math::gamma_quantile_c(alpha/2, data_ele+1.0, 1.0);
    const double data_ele_dn_err = data_ele - data_ele_dn_bound;
    const double data_ele_up_err = data_ele_up_bound - data_ele;

    const double pred_tau_from_ele = data_ele * fin_TF_tau_to_ele; 
    const double pred_tau_from_ele_stat_dn_err = data_ele_dn_err * fin_TF_tau_to_ele;
    const double pred_tau_from_ele_stat_up_err = data_ele_up_err * fin_TF_tau_to_ele;
    const double pred_tau_from_ele_syst = data_ele * finErr_TF_tau_to_ele;
    const double pred_tau_from_ele_syst_ISR_up_err = data_ele * fin_TF_tau_to_ele * systErr_rel_ISR_up_eleVec[i-1];
    const double pred_tau_from_ele_syst_ISR_dn_err = data_ele * fin_TF_tau_to_ele * systErr_rel_ISR_dn_eleVec[i-1];
    const double pred_tau_from_ele_syst_bTag_up_err = data_ele * fin_TF_tau_to_ele * systErr_rel_bTag_up_eleVec[i-1];
    const double pred_tau_from_ele_syst_bTag_dn_err = data_ele * fin_TF_tau_to_ele * systErr_rel_bTag_dn_eleVec[i-1];

    const double pred_LL_from_ele = data_ele * fin_TF_LL_to_ele; 
    const double pred_LL_from_ele_stat = data_ele ==0 ? 1.8 * fin_TF_LL_to_ele : sqrt(data_ele) * fin_TF_LL_to_ele;
    const double pred_LL_from_ele_syst = data_ele * finErr_TF_LL_to_ele;

// avg
    const double pred_tau_avg = 0.5*pred_tau_from_mu + 0.5*pred_tau_from_ele;
    const double pred_tau_avg_syst = 0.5*sqrt(pred_tau_from_mu_syst*pred_tau_from_mu_syst + pred_tau_from_ele_syst*pred_tau_from_ele_syst);
    const double pred_tau_avg_stat_up_err = 0.5*sqrt(pred_tau_from_mu_stat_up_err*pred_tau_from_mu_stat_up_err + pred_tau_from_ele_stat_up_err*pred_tau_from_ele_stat_up_err);
    const double pred_tau_avg_stat_dn_err = 0.5*sqrt(pred_tau_from_mu_stat_dn_err*pred_tau_from_mu_stat_dn_err + pred_tau_from_ele_stat_dn_err*pred_tau_from_ele_stat_dn_err);
    const double pred_tau_avg_syst_ISR_up_err = 0.5*sqrt(pred_tau_from_mu_syst_ISR_up_err*pred_tau_from_mu_syst_ISR_up_err + pred_tau_from_ele_syst_ISR_up_err*pred_tau_from_ele_syst_ISR_up_err);
    const double pred_tau_avg_syst_ISR_dn_err = 0.5*sqrt(pred_tau_from_mu_syst_ISR_dn_err*pred_tau_from_mu_syst_ISR_dn_err + pred_tau_from_ele_syst_ISR_dn_err*pred_tau_from_ele_syst_ISR_dn_err);
    const double pred_tau_avg_syst_bTag_up_err = 0.5*sqrt(pred_tau_from_mu_syst_bTag_up_err*pred_tau_from_mu_syst_bTag_up_err + pred_tau_from_ele_syst_bTag_up_err*pred_tau_from_ele_syst_bTag_up_err);
    const double pred_tau_avg_syst_bTag_dn_err = 0.5*sqrt(pred_tau_from_mu_syst_bTag_dn_err*pred_tau_from_mu_syst_bTag_dn_err + pred_tau_from_ele_syst_bTag_dn_err*pred_tau_from_ele_syst_bTag_dn_err);

    const double pred_tau_avg_rel_stat_up_err = pred_tau_avg == 0 ? 0 : pred_tau_avg_stat_up_err/pred_tau_avg;
    const double pred_tau_avg_rel_stat_dn_err = pred_tau_avg == 0 ? 0 : pred_tau_avg_stat_dn_err/pred_tau_avg;
// if pred_tau_avg = 0 then force data_mu =1 and data_ele =1 to calculate the relative syst unc.
    const double pred_tau_avg_rel_syst = pred_tau_avg == 0 ? sqrt(finErr_TF_tau_to_mu*finErr_TF_tau_to_mu + finErr_TF_tau_to_ele*finErr_TF_tau_to_ele)/(fin_TF_tau_to_mu + fin_TF_tau_to_ele) : pred_tau_avg_syst/pred_tau_avg;
    const double pred_tau_avg_rel_syst_ISR_up_err = pred_tau_avg == 0 ? sqrt(fin_TF_tau_to_mu * systErr_rel_ISR_up_muVec[i-1]*fin_TF_tau_to_mu * systErr_rel_ISR_up_muVec[i-1] + fin_TF_tau_to_ele * systErr_rel_ISR_up_eleVec[i-1]*fin_TF_tau_to_ele * systErr_rel_ISR_up_eleVec[i-1])/(fin_TF_tau_to_mu + fin_TF_tau_to_ele) : pred_tau_avg_syst_ISR_up_err/pred_tau_avg;
    const double pred_tau_avg_rel_syst_ISR_dn_err = pred_tau_avg == 0 ? sqrt(fin_TF_tau_to_mu * systErr_rel_ISR_dn_muVec[i-1]*fin_TF_tau_to_mu * systErr_rel_ISR_dn_muVec[i-1] + fin_TF_tau_to_ele * systErr_rel_ISR_dn_eleVec[i-1]*fin_TF_tau_to_ele * systErr_rel_ISR_dn_eleVec[i-1])/(fin_TF_tau_to_mu + fin_TF_tau_to_ele) : pred_tau_avg_syst_ISR_dn_err/pred_tau_avg;
    const double pred_tau_avg_rel_syst_bTag_up_err = pred_tau_avg == 0 ? sqrt(fin_TF_tau_to_mu * systErr_rel_bTag_up_muVec[i-1]*fin_TF_tau_to_mu * systErr_rel_bTag_up_muVec[i-1] + fin_TF_tau_to_ele * systErr_rel_bTag_up_eleVec[i-1]*fin_TF_tau_to_ele * systErr_rel_bTag_up_eleVec[i-1])/(fin_TF_tau_to_mu + fin_TF_tau_to_ele) : pred_tau_avg_syst_bTag_up_err/pred_tau_avg;
    const double pred_tau_avg_rel_syst_bTag_dn_err = pred_tau_avg == 0 ? sqrt(fin_TF_tau_to_mu * systErr_rel_bTag_dn_muVec[i-1]*fin_TF_tau_to_mu * systErr_rel_bTag_dn_muVec[i-1] + fin_TF_tau_to_ele * systErr_rel_bTag_dn_eleVec[i-1]*fin_TF_tau_to_ele * systErr_rel_bTag_dn_eleVec[i-1])/(fin_TF_tau_to_mu + fin_TF_tau_to_ele) : pred_tau_avg_syst_bTag_dn_err/pred_tau_avg;

    cached_rateVec.push_back(pred_tau_avg); cached_stat_upErrVec.push_back(pred_tau_avg_rel_stat_up_err); cached_stat_dnErrVec.push_back(pred_tau_avg_rel_stat_dn_err);
    cached_stat_abs_upErrVec.push_back(pred_tau_avg_stat_up_err); cached_stat_abs_dnErrVec.push_back(pred_tau_avg_stat_dn_err);
    cached_systErr_TF_statVec.push_back(pred_tau_avg_rel_syst);
    cached_systErr_ISR_upErrVec.push_back(pred_tau_avg_rel_syst_ISR_up_err); cached_systErr_ISR_dnErrVec.push_back(pred_tau_avg_rel_syst_ISR_dn_err);
    cached_systErr_bTag_upErrVec.push_back(pred_tau_avg_rel_syst_bTag_up_err); cached_systErr_bTag_dnErrVec.push_back(pred_tau_avg_rel_syst_bTag_dn_err);
    dummyVec.push_back(0.0);

    std::cout<<i-1<<"\t"<<std::setw(6)<<data_mu<<"\t"<<std::setw(6)<<data_ele
             <<"\t\t"<<std::setw(6)<<fin_TF_tau_to_mu<<" +- "<<std::setw(6)<<finErr_TF_tau_to_mu
             <<"\t"<<std::setw(6)<<fin_TF_tau_to_ele<<" +- "<<std::setw(6)<<finErr_TF_tau_to_ele
             <<"\t"<<std::setw(6)<<pred_tau_from_mu<<" + "<<std::setw(6)<<pred_tau_from_mu_stat_up_err<<" - "<<std::setw(6)<<pred_tau_from_mu_stat_dn_err<<" +- "<<std::setw(6)<<pred_tau_from_mu_syst
             <<"\t"<<std::setw(6)<<pred_tau_from_ele<<" + "<<std::setw(6)<<pred_tau_from_ele_stat_up_err<<" - "<<std::setw(6)<<pred_tau_from_ele_stat_dn_err<<" +- "<<std::setw(6)<<pred_tau_from_ele_syst
             <<"\t"<<std::setw(6)<<pred_tau_avg<<" + "<<std::setw(6)<<pred_tau_avg_stat_up_err<<" - "<<std::setw(6)<<pred_tau_avg_stat_dn_err<<" +- "<<std::setw(6)<<pred_tau_avg_syst
             <<std::endl;
  }

  std::ofstream of_hadtau("hadtau.txt", std::ofstream::out);
  char tmpstr[200];
  of_hadtau<< std::setprecision(3)<<std::fixed<<std::left;
  const int nTotChn = hTF_tau_mu->GetNbinsX();
  of_hadtau<<"luminosity = "<<data_lumi<<std::endl;
  of_hadtau<<"channels = "<<nTotChn<<std::endl;
  of_hadtau<<"sample = hadtau\n"<<std::endl;
  of_hadtau<<std::setw(21)<<"channel = "; for(int i=0; i<nTotChn; i++){ sprintf(tmpstr, "bin%d", i+1); of_hadtau<<std::setw(7)<<tmpstr<<" "; } of_hadtau<<"\n"<<std::endl;
  of_hadtau<<std::setw(21)<<"rate    = "; for(auto it : cached_rateVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau<<std::setw(21)<<"#stat_unc_abs_up = "; for(auto it : cached_stat_abs_upErrVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"#stat_unc_abs_dn = "; for(auto it : cached_stat_abs_dnErrVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;

  of_hadtau<<std::setw(21)<<"stat_unc_up = "; for(auto it : cached_stat_upErrVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"stat_unc_dn = "; for(auto it : cached_stat_dnErrVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;
// use the closure unc for the syst from TF stat
  of_hadtau<<std::setw(21)<<"syst_unc_closure_up = "; for(auto it : cached_systErr_TF_statVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;  
  of_hadtau<<std::setw(21)<<"syst_unc_closure_dn = "; for(auto it : cached_systErr_TF_statVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau<<std::setw(21)<<"syst_unc_temp_up = "; for(auto it: cached_systErr_ISR_upErrVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"syst_unc_temp_dn = "; for(auto it: cached_systErr_ISR_dnErrVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau<<std::setw(21)<<"syst_unc_mistag_up = "; for(auto it: cached_systErr_bTag_upErrVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"syst_unc_mistag_dn = "; for(auto it: cached_systErr_bTag_dnErrVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau<<std::setw(21)<<"syst_unc_pdf_up = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"syst_unc_pdf_dn = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau<<std::setw(21)<<"syst_unc_MCScale_up = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"syst_unc_MCScale_dn = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau<<std::setw(21)<<"syst_unc_Mt_up = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"syst_unc_Mt_dn = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau<<std::setw(21)<<"syst_unc_taumu_up = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"syst_unc_taumu_dn = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau<<std::setw(21)<<"syst_unc_isotrk_up = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"syst_unc_isotrk_dn = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau<<std::setw(21)<<"syst_unc_mureco_up = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"syst_unc_mureco_dn = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau<<std::setw(21)<<"syst_unc_muiso_up = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"syst_unc_muiso_dn = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau<<std::setw(21)<<"syst_unc_llovr_up = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"syst_unc_llovr_dn = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau<<std::setw(21)<<"syst_unc_trg_up = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<std::endl;
  of_hadtau<<std::setw(21)<<"syst_unc_trg_dn = "; for(auto it: dummyVec){ of_hadtau<<std::setw(7)<<it<<" "; } of_hadtau<<"\n"<<std::endl;

  of_hadtau.close();
/*
  std::cout<< std::setprecision(3)<<std::fixed;
  std::cout<<"\nPredictions"<<std::endl;
  std::cout<<"Bin \t mu CS \t ele CS \t Pred(tau)@mu \t\t\t Pred(tau)@ele \t\t\t Pred(LL)@mu \t\t\t Pred(LL)@ele"<<std::endl;
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

    std::cout<<i-1<<"\t"<<std::setw(6)<<data_mu<<"\t"<<std::setw(6)<<data_ele
             <<"\t"<<std::setw(6)<<pred_tau_from_mu<<" +- "<<std::setw(6)<<pred_tau_from_mu_stat<<" +- "<<std::setw(6)<<pred_tau_from_mu_syst
             <<"\t"<<std::setw(6)<<pred_tau_from_ele<<" +- "<<std::setw(6)<<pred_tau_from_ele_stat<<" +- "<<std::setw(6)<<pred_tau_from_ele_syst
             <<"\t"<<std::setw(6)<<pred_LL_from_mu<<" +- "<<std::setw(6)<<pred_LL_from_mu_stat<<" +- "<<std::setw(6)<<pred_LL_from_mu_syst
             <<"\t"<<std::setw(6)<<pred_LL_from_ele<<" +- "<<std::setw(6)<<pred_LL_from_ele_stat<<" +- "<<std::setw(6)<<pred_LL_from_ele_syst
             <<std::endl;
  }
*/

}
