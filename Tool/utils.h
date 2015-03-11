#include "TVector2.h"

// Collection of geometry-related helper functions.
namespace utils {
  // Return delta-phi between -Pi and Pi
  float deltaPhi(float phi1, float phi2) {
    return TVector2::Phi_mpi_pi(phi1-phi2);
  }

  // Return delta-R
  float deltaR(float eta1, float eta2, float phi1, float phi2) {
    float dphi = deltaPhi(phi1,phi2);
    float deta = eta1 - eta2;

    return sqrt( deta*deta + dphi*dphi );
  }

  // Find index 'matchedObjIdx' of the obj that is closest in deltaR
  // around the vector (eta,phi). Returns
  //  - true  : if an obj has been found.
  //  - false : otherwise. In that case, 'matchedObjIdx == -1'
  // If deltaRMax is specified, it is in addition required  that 
  // the objects are closer than deltaRMax.
  bool findMatchedObject(int &matchedObjIdx, double eta, double phi, std::vector<double> jetseta, std::vector<double> jetsphi, unsigned nObj, float deltaRMax) {
    matchedObjIdx = -1;
    double deltaRMin = 100000.;
    for(unsigned objIdx = 0; objIdx < nObj; ++objIdx) { // Loop over objects
      const double dr = deltaR(eta,jetseta.at(objIdx),phi,jetsphi.at(objIdx));
      if( dr < deltaRMin ) {
	deltaRMin = dr;
	matchedObjIdx = objIdx;
      }
    } // End of loop over objects

    bool match = false;
    if( deltaRMin < deltaRMax ) {
      match = true;
    } else {
      matchedObjIdx = -1;
    }

    return match;
  }


  double mt(double pt, double phi, double met, double metPhi) {
    const double dPhi = deltaPhi(phi,metPhi);                                
    return sqrt( 2.*pt*met*(1.-cos(dPhi)) );                                       
  }   
}
