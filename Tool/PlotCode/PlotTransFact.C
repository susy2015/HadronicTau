const double mc_lumi = 8000.0;
const double data_lumi = 36352.970569733;

const double scale_mc = data_lumi/mc_lumi;

void PlotTransFact(){
  
  TFile *file_CS = new TFile("Mix_CS.root");
  TFile *file_LL = new TFile("Mix_HadTauLL.root");
  const unsigned int kNDists = 7;
  TH1* hCS_mu[kNDists];
  TH1* hHadtau[kNDists];
  TH1* hCS_ele[kNDists];
  TH1* hLL[kNDists];
  TH1* hRatio_Hadtau_mu[kNDists];
  TH1* hRatio_LL_mu[kNDists];
  TH1* hRatio_Hadtau_ele[kNDists];
  TH1* hRatio_LL_ele[kNDists];
  for(unsigned int i = 0; i < kNDists; ++i) {
    TString name = "";
    if(      i == 0 ) name = "MET";
    else if( i == 1 ) name = "NbJets";
    else if( i == 2 ) name = "NTops";
    else if( i == 3) name = "MT2";
    else if( i == 4) name = "NJets";
    else if( i == 5) name = "HT";
    else if (i == 6) name = "Yields";
    
    // Get histograms from file
    hCS_mu[i] = (TH1D*)file_CS->Get("h"+name+"_mu");
    hCS_ele[i] = (TH1D*)file_CS->Get("h"+name+"_el");
    hHadtau[i] = (TH1D*)file_LL->Get("h"+name+"_tau");
    hLL[i] = (TH1D*)file_LL->Get("h"+name+"_LL");
    //scaling
    hCS_mu[i]->Scale(scale_mc);
    hCS_ele[i]->Scale(scale_mc);
    //setting style
    hCS_mu[i]->SetTitle("Muon CS comparison for "+name);
    hCS_mu[i]->GetYaxis()->SetTitleOffset(0.7);
    hCS_mu[i]->GetYaxis()->SetTitleSize(0.06);
    hCS_mu[i]->GetYaxis()->SetTitleFont(42);
    hCS_mu[i]->GetYaxis()->SetLabelSize(0.045);
    hCS_mu[i]->GetYaxis()->SetLabelFont(42);
    hCS_mu[i]->SetMarkerStyle(20);
    hCS_mu[i]->SetMarkerColor(kRed);
    hCS_mu[i]->SetLineColor(hCS_mu[i]->GetMarkerColor());
    hCS_mu[i]->SetMarkerSize(0.9);
    hHadtau[i]->SetLineColor(kBlack);
    hCS_mu[i]->SetStats(0);
    hHadtau[i]->SetStats(0);
    hCS_ele[i]->SetTitle("Electron CS comparison for "+name);
    hCS_ele[i]->GetYaxis()->SetTitleOffset(0.7);
    hCS_ele[i]->GetYaxis()->SetTitleSize(0.06);
    hCS_ele[i]->GetYaxis()->SetTitleFont(42);
    hCS_ele[i]->GetYaxis()->SetLabelSize(0.045);
    hCS_ele[i]->GetYaxis()->SetLabelFont(42);
    hCS_ele[i]->SetMarkerStyle(20);
    hCS_ele[i]->SetMarkerColor(kRed);
    hCS_ele[i]->SetLineColor(hCS_ele[i]->GetMarkerColor());
    hCS_ele[i]->SetMarkerSize(0.9);
    hLL[i]->SetLineColor(kBlack);
    hCS_ele[i]->SetStats(0);
    hLL[i]->SetStats(0);

    //Ratio
    hRatio_Hadtau_mu[i] = static_cast<TH1*>(hHadtau[i]->Clone("Ratio"));
    hRatio_Hadtau_mu[i]->Divide(hCS_mu[i]);
    hRatio_Hadtau_mu[i]->GetYaxis()->SetTitle("#frac{Hadtau}{CS}");
    hRatio_Hadtau_mu[i]->GetYaxis()->SetRangeUser(0,2);
    hRatio_Hadtau_mu[i]->SetTitle("");
    hRatio_Hadtau_mu[i]->SetStats(0);
    hRatio_Hadtau_mu[i]->GetYaxis()->SetTitleSize(0.1);
    hRatio_Hadtau_mu[i]->GetYaxis()->SetTitleOffset(0.4);
    hRatio_Hadtau_mu[i]->GetYaxis()->SetLabelSize(0.10);
    hRatio_Hadtau_mu[i]->GetXaxis()->SetTitleSize(0.13);
    hRatio_Hadtau_mu[i]->GetXaxis()->SetTitleOffset(0.7);
    hRatio_Hadtau_mu[i]->GetXaxis()->SetLabelSize(0.1);

    hRatio_LL_mu[i] = static_cast<TH1*>(hLL[i]->Clone("Ratio"));
    hRatio_LL_mu[i]->Divide(hCS_mu[i]);
    hRatio_LL_mu[i]->GetYaxis()->SetTitle("#frac{LL}{CS}");
    hRatio_LL_mu[i]->GetYaxis()->SetRangeUser(0,2);
    hRatio_LL_mu[i]->SetTitle("");
    hRatio_LL_mu[i]->SetStats(0);
    hRatio_LL_mu[i]->GetYaxis()->SetTitleSize(0.1);
    hRatio_LL_mu[i]->GetYaxis()->SetTitleOffset(0.4);
    hRatio_LL_mu[i]->GetYaxis()->SetLabelSize(0.10);
    hRatio_LL_mu[i]->GetXaxis()->SetTitleSize(0.13);
    hRatio_LL_mu[i]->GetXaxis()->SetTitleOffset(0.7);
    hRatio_LL_mu[i]->GetXaxis()->SetLabelSize(0.1);

    hRatio_Hadtau_ele[i] = static_cast<TH1*>(hHadtau[i]->Clone("Ratio"));
    hRatio_Hadtau_ele[i]->Divide(hCS_ele[i]);
    hRatio_Hadtau_ele[i]->GetYaxis()->SetTitle("#frac{Hadtau}{CS}");
    hRatio_Hadtau_ele[i]->GetYaxis()->SetRangeUser(0,2);
    hRatio_Hadtau_ele[i]->SetTitle("");
    hRatio_Hadtau_ele[i]->SetStats(0);
    hRatio_Hadtau_ele[i]->GetYaxis()->SetTitleSize(0.1);
    hRatio_Hadtau_ele[i]->GetYaxis()->SetTitleOffset(0.4);
    hRatio_Hadtau_ele[i]->GetYaxis()->SetLabelSize(0.10);
    hRatio_Hadtau_ele[i]->GetXaxis()->SetTitleSize(0.13);
    hRatio_Hadtau_ele[i]->GetXaxis()->SetTitleOffset(0.7);
    hRatio_Hadtau_ele[i]->GetXaxis()->SetLabelSize(0.1);

    hRatio_LL_ele[i] = static_cast<TH1*>(hLL[i]->Clone("Ratio"));
    hRatio_LL_ele[i]->Divide(hCS_ele[i]);
    hRatio_LL_ele[i]->GetYaxis()->SetTitle("#frac{LL}{CS}");
    hRatio_LL_ele[i]->GetYaxis()->SetRangeUser(0,2);
    hRatio_LL_ele[i]->SetTitle("");
    hRatio_LL_ele[i]->SetStats(0);
    hRatio_LL_ele[i]->GetYaxis()->SetTitleSize(0.1);
    hRatio_LL_ele[i]->GetYaxis()->SetTitleOffset(0.4);
    hRatio_LL_ele[i]->GetYaxis()->SetLabelSize(0.10);
    hRatio_LL_ele[i]->GetXaxis()->SetTitleSize(0.13);
    hRatio_LL_ele[i]->GetXaxis()->SetTitleOffset(0.7);
    hRatio_LL_ele[i]->GetXaxis()->SetLabelSize(0.1);

    /*
    // Create legend                                                                                                             
    TLegend* leg_mu = new TLegend(0.6,0.75,0.9,0.89);
    leg_mu->SetBorderSize(0);
    leg_mu->SetFillColor(0);
    leg_mu->SetFillStyle(0);
    leg_mu->SetTextFont(42);
    leg_mu->SetTextSize(0.05);
    leg_mu->AddEntry(hData_mu[i],"Data","L");
    leg_mu->AddEntry(hMC_mu[i],"MC","P");

    TLegend* leg_ele = new TLegend(0.6,0.75,0.9,0.89);
    leg_ele->SetBorderSize(0);
    leg_ele->SetFillColor(0);
    leg_ele->SetFillStyle(0);
    leg_ele->SetTextFont(42);
    leg_ele->SetTextSize(0.05);
    leg_ele->AddEntry(hData_ele[i],"Data","L");
    leg_ele->AddEntry(hMC_ele[i],"MC","P");

    //Draw
    TCanvas* can_mu = new TCanvas(name+"_mu",name,800,600);
    TPad *pad1_mu = new TPad("pad1_mu", "pad1_mu", 0, 0.3, 1, 1.0);
    pad1_mu->SetBottomMargin(0); // Upper and lower plot are joined                                                                      
    pad1_mu->Draw();             // Draw the upper pad: pad1                                                                    
    pad1_mu->cd();
    hMC_mu[i]->Draw("PE1");
    hData_mu[i]->Draw("histsame");
    leg_mu->Draw("same");
    can_mu->cd();
    TPad *pad2_mu = new TPad("pad2_mu", "pad2_mu", 0, 0.01, 1, 0.3);
    pad2_mu->SetTopMargin(0);
    pad2_mu->SetBottomMargin(0.2);
    pad2_mu->Draw();
    pad2_mu->cd();
    hRatio_mu[i]->Draw("PE1");
    can_mu->SaveAs(name+"_muCS.pdf");

    TCanvas* can_ele = new TCanvas(name+"_ele",name,800,600);
    TPad *pad1_ele = new TPad("pad1_ele", "pad1_ele", 0, 0.3, 1, 1.0);
    pad1_ele->SetBottomMargin(0); // Upper and lower plot are joined                                                                      
    pad1_ele->Draw();             // Draw the upper pad: pad1                                                                    
    pad1_ele->cd();
    hMC_ele[i]->Draw("PE1");
    hData_ele[i]->Draw("histsame");
    leg_ele->Draw("same");
    can_ele->cd();
    TPad *pad2_ele = new TPad("pad2_ele", "pad2_ele", 0, 0.01, 1, 0.3);
    pad2_ele->SetTopMargin(0);
    pad2_ele->SetBottomMargin(0.2);
    pad2_ele->Draw();
    pad2_ele->cd();
    hRatio_ele[i]->Draw("PE1");
    can_ele->SaveAs(name+"_eleCS.pdf");
    */
  }
  std::cout << std::setprecision(2) << std::fixed;
  std::cout<<"Muon CS"<<std::endl;
  std::cout<<"Bin \t\t Hadtua \t\t LL \t\t Muon CS \t\t\t Ratio(Hadtau/CS) \t\t Ratio(LL/CS)"<<std::endl;
  for(unsigned i=1; i<= hRatio_Hadtau_mu[6]->GetNbinsX();i++)
    {
      std::cout<<i<<"\t\t"<<hHadtau[6]->GetBinContent(i)<<" ± "<<hHadtau[6]->GetBinError(i)<<"\t\t"<<hLL[6]->GetBinContent(i)<<" ± "<<hLL[6]->GetBinError(i)<<" \t "<<hCS_mu[6]->GetBinContent(i)<<" ± "<<hCS_mu[6]->GetBinError(i)<<"  \t\t  "<<hRatio_Hadtau_mu[6]->GetBinContent(i)<<" ± "<<hRatio_Hadtau_mu[6]->GetBinError(i)<<"  \t\t  "<<hRatio_LL_mu[6]->GetBinContent(i)<<" ± "<<hRatio_LL_mu[6]->GetBinError(i)<<std::endl;
    }
  std::cout<<std::endl;
 std::cout<<"Electron CS"<<std::endl; 
 std::cout<<"Bin \t\t Hadtua \t\t LL \t\t Electron CS \t\t\t Ratio(Hadtau/CS) \t\t Ratio(LL/CS)"<<std::endl;
 for(unsigned i=1; i<= hRatio_Hadtau_ele[6]->GetNbinsX();i++)
    {
      std::cout<<i<<"\t\t"<<hHadtau[6]->GetBinContent(i)<<" ± "<<hHadtau[6]->GetBinError(i)<<"\t\t"<<hLL[6]->GetBinContent(i)<<" ± "<<hLL[6]->GetBinError(i)<<" \t "<<hCS_ele[6]->GetBinContent(i)<<" ± "<<hCS_ele[6]->GetBinError(i)<<"  \t\t  "<<hRatio_Hadtau_ele[6]->GetBinContent(i)<<" ± "<<hRatio_Hadtau_ele[6]->GetBinError(i)<<"  \t\t  "<<hRatio_LL_ele[6]->GetBinContent(i)<<" ± "<<hRatio_LL_ele[6]->GetBinError(i)<<std::endl;
    }
}
