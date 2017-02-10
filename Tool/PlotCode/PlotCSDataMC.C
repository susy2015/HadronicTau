const double mc_lumi = 36813.714859265;
const double data_lumi = 36813.714859265;

const double scale_mc = data_lumi/mc_lumi;

const bool doshape = true;

void PlotCSDataMC()
{
  
  TFile *file_Data = new TFile("Data_MET_CS.root");
  TFile *file_MC = new TFile("Mix_CS.root");
//  const unsigned int kNDists = 28;
  const unsigned int kNDists = 7;
  TH1* hMC_mu[kNDists];
  TH1* hData_mu[kNDists];
  TH1* hMC_ele[kNDists];
  TH1* hData_ele[kNDists];
  TH1* hRatio_mu[kNDists];
  TH1* hRatio_ele[kNDists];
  for(unsigned int i = 0; i < kNDists; ++i) 
  {
    TString name = "";
    if(      i == 0 ) name = "MET";
    else if( i == 1 ) name = "NbJets";
    else if( i == 2 ) name = "NTops";
    else if( i == 3 ) name = "MT2";
    else if( i == 4 ) name = "NJets";
    else if( i == 5 ) name = "HT";
    else if (i == 6 ) name = "Yields";
    else if (i == 7 ) name = "CS_MT2_2D_nb1_nt1";
    else if (i == 8 ) name = "CS_met_2D_nb1_nt1";
    else if (i == 9 ) name = "CS_nJets_2D_nb1_nt1";
    else if (i == 10) name = "CS_recoTopPt_2D_nb1_nt1";
    else if (i == 11) name = "CS_MT2_2D_nb2_nt1";
    else if (i == 12) name = "CS_met_2D_nb2_nt1";
    else if (i == 13) name = "CS_nJets_2D_nb2_nt1";
    else if (i == 14) name = "CS_recoTopPt_2D_nb2_nt1";
    else if (i == 15) name = "CS_MT2_2D_nb1_nt2";
    else if (i == 16) name = "CS_met_2D_nb1_nt2";
    else if (i == 17) name = "CS_nJets_2D_nb1_nt2";
    else if (i == 18) name = "CS_recoTopPt_2D_nb1_nt2";
    else if (i == 19) name = "CS_MT2_2D_nb2_nt2";
    else if (i == 20) name = "CS_met_2D_nb2_nt2";
    else if (i == 21) name = "CS_nJets_2D_nb2_nt2";
    else if (i == 22) name = "CS_recoTopPt_2D_nb2_nt2";
    else if (i == 23) name = "CS_HT_2D_nb3_nt1";
    else if (i == 24) name = "CS_HT_2D_nb3_nt2";
    else if (i == 25) name = "CS_HT_2D_nb3_nt3";
    else if (i == 26) name = "CS_HT_2D_nb1_nt3";
    else if (i == 27) name = "CS_HT_2D_nb2_nt3";
    
    // Get histograms from file
    if( i>=7 && i<=kNDists-1 )
    {
       hMC_mu[i] = (TH1D*)file_MC->Get("mu"+name);
       hData_mu[i] = (TH1D*)file_Data->Get("mu"+name);
       hMC_ele[i] = (TH1D*)file_MC->Get("ele"+name);
       hData_ele[i] = (TH1D*)file_Data->Get("ele"+name);
    }else
    {
       hMC_mu[i] = (TH1D*)file_MC->Get("h"+name+"_mu");
       hData_mu[i] = (TH1D*)file_Data->Get("h"+name+"_mu");
       hMC_ele[i] = (TH1D*)file_MC->Get("h"+name+"_el");
       hData_ele[i] = (TH1D*)file_Data->Get("h"+name+"_el");
    }
//    hMC_mu[i]->Draw();
    //scaling
    hMC_mu[i]->Scale(scale_mc);
    hMC_ele[i]->Scale(scale_mc);
    if( doshape ){
       std::cout<<"\n"<<name<<" -->  hMC_mu[i]->Integral : "<<hMC_mu[i]->Integral()<<"  hData_mu[i]->Integral : "<<hData_mu[i]->Integral()<<std::endl;
       std::cout<<name<<" -->  hMC_ele[i]->Integral : "<<hMC_ele[i]->Integral()<<"  hData_ele[i]->Integral : "<<hData_ele[i]->Integral()<<std::endl;
       hMC_mu[i]->Scale(hData_mu[i]->Integral()/hMC_mu[i]->Integral());
       hMC_ele[i]->Scale(hData_ele[i]->Integral()/hMC_ele[i]->Integral());
    }
    //setting style
    hMC_mu[i]->SetTitle("Data_MC Shape comparison for "+name);
    hMC_mu[i]->GetYaxis()->SetTitleOffset(0.7);
    hMC_mu[i]->GetYaxis()->SetTitleSize(0.06);
    hMC_mu[i]->GetYaxis()->SetTitleFont(42);
    hMC_mu[i]->GetYaxis()->SetLabelSize(0.045);
    hMC_mu[i]->GetYaxis()->SetLabelFont(42);
    hMC_mu[i]->SetMarkerStyle(20);
    hMC_mu[i]->SetMarkerColor(kRed);
    hMC_mu[i]->SetLineColor(hMC_mu[i]->GetMarkerColor());
    hMC_mu[i]->SetMarkerSize(0.9);
    hData_mu[i]->SetLineColor(kBlack);
    hMC_mu[i]->SetStats(0);
    hData_mu[i]->SetStats(0);
    hMC_ele[i]->SetTitle("Data_MC Shape comparison for "+name);
    hMC_ele[i]->GetYaxis()->SetTitleOffset(0.7);
    hMC_ele[i]->GetYaxis()->SetTitleSize(0.06);
    hMC_ele[i]->GetYaxis()->SetTitleFont(42);
    hMC_ele[i]->GetYaxis()->SetLabelSize(0.045);
    hMC_ele[i]->GetYaxis()->SetLabelFont(42);
    hMC_ele[i]->SetMarkerStyle(20);
    hMC_ele[i]->SetMarkerColor(kRed);
    hMC_ele[i]->SetLineColor(hMC_ele[i]->GetMarkerColor());
    hMC_ele[i]->SetMarkerSize(0.9);
    hData_ele[i]->SetLineColor(kBlack);
    hMC_ele[i]->SetStats(0);
    hData_ele[i]->SetStats(0);

    //Ratio
    hRatio_mu[i] = static_cast<TH1*>(hData_mu[i]->Clone("Ratio"));
    hRatio_mu[i]->Divide(hMC_mu[i]);
    hRatio_mu[i]->GetYaxis()->SetTitle("#frac{Data}{MC}");
    hRatio_mu[i]->GetYaxis()->SetRangeUser(0,2);
    hRatio_mu[i]->SetTitle("");
    hRatio_mu[i]->SetStats(0);
    hRatio_mu[i]->GetYaxis()->SetTitleSize(0.1);
    hRatio_mu[i]->GetYaxis()->SetTitleOffset(0.4);
    hRatio_mu[i]->GetYaxis()->SetLabelSize(0.10);
    hRatio_mu[i]->GetXaxis()->SetTitleSize(0.13);
    hRatio_mu[i]->GetXaxis()->SetTitleOffset(0.7);
    hRatio_mu[i]->GetXaxis()->SetLabelSize(0.1);
    hRatio_mu[i]->SetLineWidth(2); hRatio_mu[i]->SetMarkerSize(2);

    hRatio_ele[i] = static_cast<TH1*>(hData_ele[i]->Clone("Ratio"));
    hRatio_ele[i]->Divide(hMC_ele[i]);
    hRatio_ele[i]->GetYaxis()->SetTitle("#frac{Data}{MC}");
    hRatio_ele[i]->GetYaxis()->SetRangeUser(0,2);
    hRatio_ele[i]->SetTitle("");
    hRatio_ele[i]->SetStats(0);
    hRatio_ele[i]->GetYaxis()->SetTitleSize(0.1);
    hRatio_ele[i]->GetYaxis()->SetTitleOffset(0.4);
    hRatio_ele[i]->GetYaxis()->SetLabelSize(0.10);
    hRatio_ele[i]->GetXaxis()->SetTitleSize(0.13);
    hRatio_ele[i]->GetXaxis()->SetTitleOffset(0.7);
    hRatio_ele[i]->GetXaxis()->SetLabelSize(0.1);
    hRatio_ele[i]->SetLineWidth(2); hRatio_ele[i]->SetMarkerSize(2);

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
    hMC_mu[i]->SetMaximum(hMC_mu[i]->GetMaximum()*1.25);
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
    if( doshape )
    {
      TF1 * fline = new TF1("line", "pol0", hRatio_mu[i]->GetBinLowEdge(1), hRatio_mu[i]->GetBinLowEdge(hRatio_mu[i]->GetNbinsX()) + hRatio_mu[i]->GetBinWidth(hRatio_mu[i]->GetNbinsX()));
      fline->SetParameter(0, 1);
      fline->SetLineStyle(2);
      fline->Draw("same");
    }
    can_mu->SaveAs(name+"_muCS.pdf");

    TCanvas* can_ele = new TCanvas(name+"_ele",name,800,600);
    TPad *pad1_ele = new TPad("pad1_ele", "pad1_ele", 0, 0.3, 1, 1.0);
    pad1_ele->SetBottomMargin(0); // Upper and lower plot are joined                                                                      
    pad1_ele->Draw();             // Draw the upper pad: pad1                                                                    
    pad1_ele->cd();
    hMC_ele[i]->SetMaximum(hMC_ele[i]->GetMaximum()*1.25);
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
    if( doshape )
    {
      TF1 * fline = new TF1("line", "pol0", hRatio_ele[i]->GetBinLowEdge(1), hRatio_ele[i]->GetBinLowEdge(hRatio_ele[i]->GetNbinsX()) + hRatio_ele[i]->GetBinWidth(hRatio_ele[i]->GetNbinsX()));
      fline->SetParameter(0, 1);
      fline->SetLineStyle(2);
      fline->Draw("same");
    }
    can_ele->SaveAs(name+"_eleCS.pdf");
    
  }
  if( !doshape )
  {
    std::cout << std::setprecision(2) << std::fixed;
    std::cout<<"Muon CS"<<std::endl;
    std::cout<<"Bin \t\t\t Data \t\t\t MC \t\t\t Ratio(Data/MC)"<<std::endl;
    for(unsigned i=1; i<= hRatio_mu[6]->GetNbinsX();i++)
    {
      std::cout<<i<<"\t\t\t"<<hData_mu[6]->GetBinContent(i)<<" ± "<<hData_mu[6]->GetBinError(i)<<" \t "<<hMC_mu[6]->GetBinContent(i)<<" ± "<<hMC_mu[6]->GetBinError(i)<<"  \t\t  "<<hRatio_mu[6]->GetBinContent(i)<<" ± "<<hRatio_mu[6]->GetBinError(i)<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<"Electron CS"<<std::endl;
    for(unsigned i=1; i<= hRatio_ele[6]->GetNbinsX();i++)
    {
      std::cout<<i<<"\t\t\t"<<hData_ele[6]->GetBinContent(i)<<" ± "<<hData_ele[6]->GetBinError(i)<<" \t "<<hMC_ele[6]->GetBinContent(i)<<" ± "<<hMC_ele[6]->GetBinError(i)<<"  \t\t  "<<hRatio_ele[6]->GetBinContent(i)<<" ± "<<hRatio_ele[6]->GetBinError(i)<<std::endl;
    }
  }
}
