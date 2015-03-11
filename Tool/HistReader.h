#ifndef HIST_READER_H
#define HIST_READER_H

#include <exception>
#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TString.h"


// Get histogram from ROOT file.
class HistReader {
public:
  static TH1* get(const char *fileName, const TString &histName, const TString &xTitle = "", const TString &yTitle = "");
};


TH1* HistReader::get(const char *fileName, const TString &histName, const TString &xTitle, const TString &yTitle) {
  TFile *file = new TFile(fileName,"READ");
  TH1* h = NULL;
  file->GetObject(histName,h);
  if( h == NULL ) {
    std::cerr << "\n\nERROR: Histogram '" << histName << "' not found\n\n" << std::endl;
    throw std::exception();
  } else {
    h->SetDirectory(0);
    h->UseCurrentStyle();
    if( xTitle != "" ) h->GetXaxis()->SetTitle(xTitle);
    if( yTitle != "" ) h->GetYaxis()->SetTitle(yTitle);
  }
  file->Close();

  return h;
}

#endif
