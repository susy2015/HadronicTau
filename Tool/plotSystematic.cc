#include "SusyAnaTools/Tools/tdrstyle.h"
#include "SusyAnaTools/Tools/searchBins.h"

#include <stdio.h>      /* printf */
#include <math.h>       /* sqrt */

#include <sstream>
#include <iostream>
#include <fstream>

#include "TStopwatch.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"

#include "TPrincipal.h"

#include "THStack.h"
#include "TLine.h"

#include "TLatex.h"

#include "TF1.h"

#include "TFile.h"

#include "TGraphAsymmErrors.h"

bool doZoom = false;

SearchBins * sb =0;

int main(int argc, char* argv[])
{
   if (argc < 2)
   {
       std::cerr <<"Please give 1 arguments "<<"dataCardFile"<<std::endl;
       std::cerr <<"Valid configurations are " << std::endl;
       std::cerr <<"./PlotSystematic hadtau.txt" << std::endl;
       return -1;
   }

   const std::string dataCardFile = argv[1];
   std::cout<<"\nReading in "<<dataCardFile.c_str()<<" ..."<<std::endl;

   sb = new SearchBins();

   int nTotBins = sb->nSearchBins();

   setTDRStyle();

   std::map<std::string, std::vector<double>> parsingMap;

   std::ifstream fin(dataCardFile.c_str());

   std::string line;
   std::string key, dump, item;
   while(std::getline(fin, line) )
   {
      std::stringstream ss(line);
      ss >> key; ss >> dump; // key and "=" in dump
      if( key.find("#") != std::string::npos ){
         while(key.find_last_of("#") == key.size()-1)
         {
            key = dump; ss >> dump;
         }
         key = key.substr(key.find_last_of("#")+1);
      }
      if(   key == "channels" || key == "luminosity" || key == "channel" || key == "sample" || key == "note" ) continue;
      while(ss >> item)
      {
         double strToDouble = std::atof(item.c_str());
         if(strToDouble != strToDouble ) strToDouble = 0.0; // check if the number is NAN
         parsingMap[key].push_back(strToDouble);
      }
   }

   std::map<std::string, std::vector<TH1D*>> systHistMap;
   std::map<std::string, std::string> keyDisptMap;
   std::map<std::string, int> keyColorMap;
   TH1D * h1_rate =0, * h1_stat_up =0, * h1_stat_dn =0;
   std::vector<double> rateVec, stat_upVec, stat_dnVec, syst_upVec, syst_dnVec, sum_upVec, sum_dnVec;
   TH1D * h1_TF_to_mu =0, * h1_TF_to_ele =0;
   int nBins = -1;
   // first loop to complete the histogram of the systematics mapper -> one key (one source of systematics) could have
   // multiple (dn/up -> 2 ) uncerntainties 
   for(auto miter : parsingMap)
   {
      std::string key;
      if(miter.first.find("_up") == miter.first.size()-3 || miter.first.find("_dn") == miter.first.size()-3 )
      {
         key = miter.first.substr(0, miter.first.size()-3); // strip out the "_up" or "_dn"
      }else
      {
         key = miter.first;
      }
      std::string dispt;
      std::stringstream ss(key);
      while(std::getline(ss, dump, '_'))
      {
         if(dump == "syst" || dump == "unc" || dump == "abs") continue;
         dispt += dump + " ";
      }
      if( nBins == -1 )
      {
         nBins = miter.second.size();
         std::cout<<"nBins : "<<nBins<<std::endl;
         h1_rate = new TH1D("rate", "rate", nBins, 0, nBins); h1_rate->Sumw2();
         h1_stat_up = new TH1D("stat_up", "stat_up", nBins, 0, nBins); h1_stat_up->Sumw2();
         h1_stat_dn = new TH1D("stat_dn", "stat_dn", nBins, 0, nBins); h1_stat_dn->Sumw2();
         h1_TF_to_mu = new TH1D("TF_to_mu", "TF_to_mu", nBins, 0, nBins); h1_TF_to_mu->Sumw2();
         h1_TF_to_ele = new TH1D("TF_to_ele", "TF_to_ele", nBins, 0, nBins); h1_TF_to_ele->Sumw2();

         if( nTotBins != nBins){ std::cout<<"ERROR   nBins("<<nBins<<") from data card is NOT the same as nTotBins("<<nTotBins<<") from SearchBins class!"<<std::endl; return 0; }
      }

      if( key == "rate" ) { rateVec = miter.second; for(unsigned int iv=0; iv<miter.second.size(); iv++) h1_rate->Fill(iv, miter.second[iv]); }
      if( miter.first == "stat_unc_abs_up" ) { stat_upVec = miter.second; for(unsigned int iv=0; iv<miter.second.size(); iv++) h1_stat_up->Fill(iv, miter.second[iv]); }
      if( miter.first == "stat_unc_abs_dn" ) { stat_dnVec = miter.second; for(unsigned int iv=0; iv<miter.second.size(); iv++) h1_stat_dn->Fill(iv, miter.second[iv]); }

      if( key == "fin_TF_to_mu" ){ for(unsigned int iv=0; iv<miter.second.size(); iv++) h1_TF_to_mu->SetBinContent(iv+1, miter.second[iv]); }
      if( key == "fin_TFerr_to_mu" ){ for(unsigned int iv=0; iv<miter.second.size(); iv++) h1_TF_to_mu->SetBinError(iv+1, miter.second[iv]); }
      if( key == "fin_TF_to_ele" ){ for(unsigned int iv=0; iv<miter.second.size(); iv++) h1_TF_to_ele->SetBinContent(iv+1, miter.second[iv]); }
      if( key == "fin_TFerr_to_ele" ){ for(unsigned int iv=0; iv<miter.second.size(); iv++) h1_TF_to_ele->SetBinError(iv+1, miter.second[iv]); }

      // skip stat unc and pred as they will be shown in a seperate histogram
      if( key == "rate" || key == "stat_unc_abs" || key == "stat_unc" || key.find("fin_TF") != std::string::npos ) continue;
      keyDisptMap[key] = dispt;
      systHistMap[key].emplace_back(new TH1D(miter.first.c_str(), "Systematics", nBins, 0, nBins));
      for(unsigned int iv=0; iv<miter.second.size(); iv++)
      {
         const double tofill = miter.second[iv] >=1.00? 0 : miter.second[iv];
         systHistMap[key].back()->Fill(iv, tofill);
      }
   }

   std::map<std::string, std::vector<TH1D*>> systAccHistSqrtMap = systHistMap;
   std::vector<std::string> cachedKeyVec;
   int nType = -1;
   for(auto miter : systHistMap)
   {
      if( !cachedKeyVec.empty() )
      {
         if( miter.second.size() != systHistMap[cachedKeyVec.back()].size() ){ std::cout<<"mis-matching in vector<TH1D*> size!?"<<std::endl; return 0; }
         nType = (int)miter.second.size();
         for(int ih=0; ih<nType; ih++)
         {
            for(int ib=0; ib<nBins; ib++)
            {
               double pre_sqrt = systAccHistSqrtMap[cachedKeyVec.back()][ih]->GetBinContent(ib+1);
               double cur_cen = miter.second[ih]->GetBinContent(ib+1);
               double cur_sqrt = sqrt(cur_cen*cur_cen + pre_sqrt*pre_sqrt);
               systAccHistSqrtMap[miter.first][ih]->SetBinContent(ib+1, cur_sqrt); 
            }
         }
      }
      cachedKeyVec.push_back(miter.first);
   }

   for(auto miter : keyDisptMap )
   {
      if( miter.first == "syst_unc_ISR" ) keyColorMap[miter.first] = kTeal -7;
      if( miter.first == "syst_unc_SF" ) keyColorMap[miter.first] = kAzure-4;
      if( miter.first == "syst_unc_TF_stat" ) keyColorMap[miter.first] = kMagenta-2;
      if( miter.first == "syst_unc_bTag" ) keyColorMap[miter.first] = kOrange+2;
      if( miter.first == "syst_unc_jec" ) keyColorMap[miter.first] = kAzure+2;
      if( miter.first == "syst_unc_metMag" ) keyColorMap[miter.first] = kBlue;
      if( miter.first == "syst_unc_metPhi" ) keyColorMap[miter.first] = kYellow;
      if( miter.first == "syst_unc_pdfUnc" ) keyColorMap[miter.first] = kRed;
   }

   std::vector<double> xbin, xbinup, xbindn;
   for(int ib=0; ib<nBins; ib++){ xbin.push_back(0.5+ib); xbinup.push_back(0.50); xbindn.push_back(0.50); }
   for(int ib=0; ib<nBins; ib++)
   {
      h1_rate->SetBinError(ib+1, 0);
      const double stat_up = stat_upVec[ib];
      const double syst_up = systAccHistSqrtMap.rbegin()->second[1]->GetBinContent(ib+1) * rateVec[ib];
      double sum_up = sqrt(stat_up*stat_up + syst_up*syst_up);
      syst_upVec.push_back(syst_up);
      sum_upVec.push_back(sum_up); 
 
      const double stat_dn = stat_dnVec[ib];
      double syst_dn = systAccHistSqrtMap.rbegin()->second[1]->GetBinContent(ib+1) * rateVec[ib];
      double sum_dn = sqrt(stat_dn*stat_dn + syst_dn*syst_dn);
      if( sum_dn > rateVec[ib] ) sum_dn = rateVec[ib];
      if( syst_dn > rateVec[ib] ) syst_dn = rateVec[ib];
      syst_dnVec.push_back(syst_dn);
      sum_dnVec.push_back(sum_dn);
   }

   sb->print_searchBinsPred_latex(rateVec, stat_upVec, stat_dnVec, syst_upVec, syst_dnVec, "& $\\tauh$ Prediction \\\\");

// Start making plots
   tdrStyle->SetOptTitle(1);
   tdrStyle->SetTitleFont(42,"");
   tdrStyle->SetTitleColor(1);
   tdrStyle->SetTitleTextColor(1);
   tdrStyle->SetTitleFillColor(0);
   tdrStyle->SetTitleFontSize(0.1);
   tdrStyle->SetTitleAlign(13);
   tdrStyle->SetTitleX(0.6);
   tdrStyle->SetTitleH(0.05);
   tdrStyle->SetTitleBorderSize(0);
   tdrStyle->SetTitleAlign(13);
   tdrStyle->SetTitleX(0.19);
   tdrStyle->SetTitleH(0.038);
                                                                         
     //  For the frame
   tdrStyle->SetFrameBorderMode(0);
   tdrStyle->SetFrameBorderSize(1);
   tdrStyle->SetFrameFillColor(kBlack);
   tdrStyle->SetFrameFillStyle(0);
   tdrStyle->SetFrameLineColor(kBlack);
   tdrStyle->SetFrameLineStyle(0);
   tdrStyle->SetFrameLineWidth(1);
 
      //  For the Pad
   tdrStyle->SetPadBorderMode(0);
   tdrStyle->SetPadColor(kWhite);
   tdrStyle->SetPadGridX(false);
   tdrStyle->SetPadGridY(false);
   tdrStyle->SetGridColor(0);
   tdrStyle->SetGridStyle(3);
   tdrStyle->SetGridWidth(1);
 
      //  For the axis
   tdrStyle->SetAxisColor(1,"XYZ");
   tdrStyle->SetTickLength(0.03,"XYZ");
   tdrStyle->SetNdivisions(505,"XYZ");
   tdrStyle->SetPadTickX(1);
   tdrStyle->SetPadTickY(1);
   tdrStyle->SetStripDecimals(kFALSE);
 
   tdrStyle->SetLabelSize(0.050, "XYZ");
 
   tdrStyle->SetPadTopMargin(0.12); tdrStyle->SetPadBottomMargin(0.20);
   tdrStyle->SetPadLeftMargin(0.10); tdrStyle->SetPadRightMargin(0.05);
 
   tdrStyle->SetHistLineWidth(1);
 
   tdrStyle->SetPaintTextFormat("4.2f");

   tdrStyle->SetTitleXOffset(5.50); tdrStyle->SetTitleYOffset(6.50);

   const double width = 400 * nType, height = 400;
   TCanvas * cd = new TCanvas("cd", "cd", width, height);
   cd->Divide(nType, 1);
   cd->Print("systUnc.pdf[");
   for(auto miter : systAccHistSqrtMap)
   {
      for(int it=0; it<nType; it++)
      {
         cd->cd(it+1);
         miter.second[it]->SetStats(kFALSE);
         miter.second[it]->SetTitle(miter.first.c_str());
         miter.second[it]->Draw("hist");
      }
      cd->Print("systUnc.pdf");
   }
   cd->Print("systUnc.pdf]");

// 
   TCanvas * c1 = new TCanvas("c1", "c1", 800, 800);
   TLegend* leg = new TLegend(0.30,0.60,0.6,0.85);
   leg->SetFillColor(0);
   leg->SetTextSize(0.035);
   leg->SetTextFont(62);
   leg->SetBorderSize(0);
   std::map<std::string, std::vector<TH1D*>>::reverse_iterator rit;
   for(int it=0; it<nType; it++)
   {
     leg->Clear();
     for(rit = systAccHistSqrtMap.rbegin(); rit != systAccHistSqrtMap.rend(); ++rit)
     {
        rit->second[it]->SetFillColor(keyColorMap[rit->first]);
        rit->second[it]->SetLineColor(keyColorMap[rit->first]);
        rit->second[it]->SetMarkerColor(keyColorMap[rit->first]);
        if( it == 0 ) rit->second[it]->SetTitle("Systematics Down; Search Bins; ");
        if( it == 1 ) rit->second[it]->SetTitle("Systematics Up; Search Bins; ");
        if( cachedKeyVec.back() == rit->first )
        {
           rit->second[it]->GetYaxis()->SetRangeUser(0, 2.0);
           rit->second[it]->Draw("hist");
        }else
        {
           rit->second[it]->Draw("same hist");  
        }
        leg->AddEntry(rit->second[it], keyDisptMap[rit->first].c_str(), "F");
     }
     leg->Draw();
     c1->Print(TString("stack_syst" + std::to_string(it) + ".pdf"));
   }

//
   TLegend* leg_TF = new TLegend(0.30,0.75,0.6,0.85);
   leg_TF->SetFillColor(0);
   leg_TF->SetTextSize(0.035);
   leg_TF->SetTextFont(62);
   leg_TF->SetBorderSize(0);
   c1->cd();
   h1_TF_to_ele->GetYaxis()->SetRangeUser(0, h1_TF_to_ele->GetMaximum()*1.2);
   h1_TF_to_ele->SetLineColor(kRed);
   h1_TF_to_ele->SetMarkerColor(kRed);
   h1_TF_to_ele->Draw();
   c1->Print("ele_TF.pdf");

   c1->cd();
   h1_TF_to_mu->GetYaxis()->SetRangeUser(0, h1_TF_to_mu->GetMaximum()*1.2);
   h1_TF_to_mu->SetLineColor(kBlue);
   h1_TF_to_mu->SetMarkerColor(kBlue);
   h1_TF_to_mu->Draw();
   c1->Print("mu_TF.pdf");

   leg_TF->Clear();
   h1_TF_to_mu->GetYaxis()->SetRangeUser(0, h1_TF_to_mu->GetMaximum()*1.5);
   h1_TF_to_mu->SetTitle("Translation Factor; Search Bins; ");
   h1_TF_to_mu->Draw();
   h1_TF_to_ele->Draw("same");
   leg_TF->AddEntry(h1_TF_to_mu, "TF for #mu CS", "L"); 
   leg_TF->AddEntry(h1_TF_to_ele, "TF for e CS", "L");
   leg_TF->Draw();
   c1->Print("comp_TF.pdf");

/*
   leg->Clear();
   TCanvas *cr = new TCanvas("cr", "cr", 800, 600);
   cr->cd();
   TPad * pad1_cr = new TPad("pad1_cr", "pad1_cr", 0, 0.3, 1, 1.0);
   pad1_cr->SetBottomMargin(0); // Upper and lower plot are joined                                                                      
   pad1_cr->Draw();             // Draw the upper pad: pad1                                                                    
   pad1_cr->cd();
   h1_TF_to_mu->SetTitle("Translation Factor; Search Bins; ");
   h1_TF_to_mu->Draw();
   h1_TF_to_ele->Draw("same");
   leg->AddEntry(h1_TF_to_mu, "TF for #mu CS", "L"); 
   leg->AddEntry(h1_TF_to_ele, "TF for e CS", "L");
   cr->cd();
   TPad * pad2_cr = new TPad("pad2_cr", "pad2_cr", 0, 0.01, 1, 0.3);
   pad2_cr->SetTopMargin(0);
   pad2_cr->SetBottomMargin(0.2);
   pad2_cr->Draw();
   pad2_cr->cd();
   TH1D * ratio_TF_mu_OVER_e = (TH1D*) h1_TF_to_mu->Clone("ratio");
   ratio_TF_mu_OVER_e->Divide(h1_TF_to_ele);
   ratio_TF_mu_OVER_e->GetYaxis()->SetTitle("#frac{#mu}{e}");
   ratio_TF_mu_OVER_e->GetYaxis()->SetRangeUser(0, 2.);
   ratio_TF_mu_OVER_e->SetTitle("");
   ratio_TF_mu_OVER_e->SetStats(0);
   ratio_TF_mu_OVER_e->Draw();
   cr->Print("comp_TF.pdf");
*/

// 
   TGraphAsymmErrors * gr_stat_AsymErr = new TGraphAsymmErrors(nBins, &xbin[0], &rateVec[0], &xbindn[0], &xbinup[0], &stat_dnVec[0], &stat_upVec[0]);
   gr_stat_AsymErr->SetFillColor(kBlue+2); gr_stat_AsymErr->SetFillStyle(3244); gr_stat_AsymErr->SetLineWidth(1); gr_stat_AsymErr->SetLineColor(0);
   TGraphAsymmErrors * gr_sum_AsymErr = new TGraphAsymmErrors(nBins, &xbin[0], &rateVec[0], &xbindn[0], &xbinup[0], &sum_dnVec[0], &sum_upVec[0]);
   gr_sum_AsymErr->SetFillColor(kGray+2); gr_sum_AsymErr->SetFillStyle(3244); gr_sum_AsymErr->SetLineWidth(1); gr_sum_AsymErr->SetLineColor(0); gr_sum_AsymErr->SetMarkerSize(0); gr_sum_AsymErr->SetMarkerColor(0);

  // Create legend
   TLegend* leg2 = new TLegend(0.58,0.75,0.93,0.87);
   leg2->SetBorderSize(1);
   leg2->SetLineColor(1);
   leg2->SetLineWidth(2);
   leg2->SetFillColor(0);
   leg2->SetTextFont(42);
   leg2->SetTextSize(0.035);
   leg2->SetHeader("Hadronic #tau Background");
   leg2->AddEntry(h1_rate,"Data Prediction.","p");

   double legendX1 = .58;
   double legendX2 = .88;
   double legendY1 = .65;
   double legendY2 = .73;
   TLegend* catLeg_unc = new TLegend(legendX1,legendY1,legendX2,legendY2);
   catLeg_unc->SetTextSize(0.030);
   catLeg_unc->SetTextFont(42);

   catLeg_unc->AddEntry(gr_sum_AsymErr, "Bkg. Syst. Unc.", "F");
   catLeg_unc->AddEntry(gr_stat_AsymErr, "Bkg. Stat. Unc.", "F");

   if( doZoom ) h1_rate->GetYaxis()->SetRangeUser(0, 80);
   else h1_rate->GetYaxis()->SetRangeUser(0, h1_rate->GetMaximum()*1.15); 
   h1_rate->SetStats(kFALSE);
   h1_rate->SetTitle(";Search Bins;");
   h1_rate->SetMarkerSize(0.8);
   h1_rate->SetMarkerStyle(20);
   h1_rate->Draw("p");

   gr_sum_AsymErr->Draw("2");
   gr_stat_AsymErr->Draw("2");
   leg2->Draw("same");

   catLeg_unc->SetFillColor(kWhite);
   catLeg_unc->SetBorderSize(0);
   catLeg_unc->Draw();
    
   c1->Print("pred.pdf");

   return 1;
}
