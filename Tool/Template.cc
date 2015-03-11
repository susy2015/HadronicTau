#include <iostream>
#include <algorithm>
#include <cstring>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>


#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/customize.h"
#include "TStopwatch.h"
#include "TString.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Template.h"

#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TChain.h"

#include "TauResponse.h"
#include "utils.h"




// === Main Function ===================================================
int main(int argc, char* argv[]) {


  if (argc < 2)
    {
      std::cerr <<"Please give 2 arguments " << "inputList " << " " << "outputFileName" << std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Template List_ttbar.txt HadTau_TauResponseTemplates.root" << std::endl;
      return -1;
    }
  const char *inputFileList = argv[1];
  const char *outFileName = argv[2];

  TChain *fChain = new TChain("Hadtau");
  if(!FillChain(fChain, inputFileList))
    {
      std::cerr << "Cannot get the tree " << std::endl;
    }

  NTupleReader tr(fChain);

  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(outFileName);

  
  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;

  // Loop over the events (tree entries)
  int k = 0;
  while(tr.getNextEvent()){
    k++;
    std::vector<double>jetspt = tr.getVec<double>("jetspt");
    std::vector<double>jetseta = tr.getVec<double>("jetseta");
    std::vector<double>jetsphi = tr.getVec<double>("jetsphi");
    std::vector<double>gentaupt = tr.getVec<double>("gentaupt");
    std::vector<double>gentaueta = tr.getVec<double>("gentaueta");
    std::vector<double>gentauphi = tr.getVec<double>("gentauphi");
    int hadtauflag = tr.getVar<int>("hadtauflag");

    // Select only events where the W decayed into a hadronically decaying tau
    if( !(hadtauflag == 1) ) continue;
    // Kinematic variables of generator-level tau
    const double genTauPt = gentaupt.at(0);
    const double genTauEta = gentaueta.at(0);
    const double genTauPhi = gentauphi.at(0);

    // Use only events where the tau is inside the muon acceptance (pt>20, eta>2.1)because lateron we will apply the response to muon+jet events
    if( genTauPt < TauResponse::ptMin() ) continue;
    if( fabs(genTauEta) > TauResponse::etaMax() ) continue;




    // Do the matching
    int tauJetIdx = -1;		// Will store the index of the jet matched to the tau
    const float deltaRMax = genTauPt < 50. ? 0.2 : 0.1; // Increase deltaRMax at low pt to maintain high-enought matching efficiency
    const unsigned nObj = jetspt.size();
    if( !utils::findMatchedObject(tauJetIdx,genTauEta,genTauPhi,jetseta,jetsphi,nObj,deltaRMax) ) continue;
    

    // Fill histogram with relative visible energy of the tau
    // ("tau response template") for hadronically decaying taus
    for(unsigned jetIdx = 0; jetIdx < jetspt.size(); ++jetIdx) {	// Loop over reco jets
      // Select tau jet
      if( jetIdx == tauJetIdx ) {
	// Get the response pt bin for the tau
	const unsigned int ptBin = TauResponse::ptBin(genTauPt);
	// Fill the corresponding response template
	const double tauJetPt = jetspt.at(jetIdx);
	myBaseHistgram.hTauResp.at(ptBin)->Fill( tauJetPt / genTauPt );
	myBaseHistgram.hTauvisible.at(ptBin)->Fill(tauJetPt);
	myBaseHistgram.hTaugenerated.at(ptBin)->Fill(genTauPt);
	break;		// End the jet loop once the tau jet has been found
      }
    }	// End of loop over reco jets
  } // End of loop over tree entries

  std::cout << "Total Event: "<<k<<std::endl;
  // Normalize the response distributions to get the probability density
   for(unsigned int i = 0; i <myBaseHistgram. hTauResp.size(); ++i) {
    if(myBaseHistgram. hTauResp.at(i)->Integral("width") > 0. ) {
      myBaseHistgram.hTauResp.at(i)->Scale(1./myBaseHistgram.hTauResp.at(i)->Integral("width"));
    }
  }


  // --- Save the Histograms to File -----------------------------------
   (myBaseHistgram.oFile)->Write();

   return 0;
}
