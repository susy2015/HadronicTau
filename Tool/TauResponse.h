#ifndef TAU_RESPONSE_H
#define TAU_RESPONSE_H

#include <iostream>
#include <exception>

#include "TH1.h"
#include "TString.h"

#include "HistReader.h"


// Tau-reponse templates and related information

class TauResponse {
public:
  // Definition of the pt bins for which separate
  // response templates are determined
  static unsigned int nBins() { return 4; }
  static unsigned int ptBin(double pt);

  // The acceptance inside which the response
  // distributions are valid. Corresponds to the
  // muon acceptance in the analysis.
  static double ptMin() { return ptMin(0); }
  static double etaMax() { return 2.4; }

  // Pt bin edges
  static double ptMin(unsigned int ptBin);
  static double ptMax(unsigned int ptBin);

  // Name of the response-template histogram
  // for the given pt bin
  static TString name(unsigned int ptBin);


  // Constructor. Takes name of file with response histograms
  // as input
  TauResponse(const char *fileName);

  // Get random number sampled from response distribution
  // The response distribution is chosen according to
  // the specified pt
  double getRandom(double pt) const { return resp_.at(ptBin(pt))->GetRandom(); }



private:
  static void checkPtBin(unsigned int ptBin);

  std::vector<TH1*> resp_;	//  the response distribtuons
};


TauResponse::TauResponse(const char *fileName) {
  for(unsigned int i = 0; i < nBins(); ++i) {
    resp_.push_back(HistReader::get(fileName,TauResponse::name(i)));
  }
}


unsigned int TauResponse::ptBin(double pt) {
  if( pt < ptMin() ) {
    std::cerr << "\n\nERROR in TauResponse::ptBin" << std::endl;
    std::cerr << "  No response available for pt = " << pt << " < " << ptMin() << std::endl;
    throw std::exception();
  }

  unsigned int bin = 0;
  if( pt > 30. )  bin = 1;
  if( pt > 50. )  bin = 2;
  if( pt > 100. ) bin = 3;

  return bin;
}


double TauResponse::ptMin(unsigned int ptBin) {
  checkPtBin(ptBin);
  double pt = 20.;
  if(      ptBin == 1 ) pt = 30.;
  else if( ptBin == 2 ) pt = 50.;
  else if( ptBin == 3 ) pt = 100.;

  return pt;
}


double TauResponse::ptMax(unsigned int ptBin) {
  checkPtBin(ptBin);
  double pt = 30.;
  if(      ptBin == 1 ) pt = 50.;
  else if( ptBin == 2 ) pt = 100.;
  else if( ptBin == 3 ) pt = 10000.;

  return pt;
}


TString TauResponse::name(unsigned int ptBin) {
  checkPtBin(ptBin);
  TString name = "hTauResp_";
  name += ptBin;

  return name;
}


void TauResponse::checkPtBin(unsigned int ptBin) {
  if( ptBin > 3 ) {
    std::cerr << "\n\nERROR in TauResponse: pt bin " << ptBin << " out of binning" << std::endl;
    throw std::exception();
  }
}

#endif
