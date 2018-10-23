#include <iostream>
#include <vector>

#include "PrepSystematics.h"

PrepSystematics::PrepSystematics(const bool isusegenmet) :
  usegenmet(isusegenmet)
{

}

void PrepSystematics::systematicPrep(NTupleReader & tr){

   const std::vector<TLorentzVector>& jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
   const std::vector<double> & recoJetsBtag = tr.getVec<double>("recoJetsBtag_0");
   const std::vector<double>& recoJetsJecUnc = tr.getVec<double>("recoJetsJecUnc");

   const std::vector<double>& metMagUpVec = tr.getVec<double>("metMagUp");
   const std::vector<double>& metMagDnVec = tr.getVec<double>("metMagDown");
   const std::vector<double>& metPhiUpVec = tr.getVec<double>("metPhiUp");
   const std::vector<double>& metPhiDnVec = tr.getVec<double>("metPhiDown");

   const double genmet = tr.hasVar("genmet") ? tr.getVar<double>("genmet") : 0;
   const double genmetphi = tr.hasVar("genmetphi") ? tr.getVar<double>("genmetphi") : -999;

   const double met = ( usegenmet && tr.hasVar("genmet") )? genmet : tr.getVar<double>("met");
   const double metphi = ( usegenmet && tr.hasVar("genmetphi") )? genmetphi : tr.getVar<double>("metphi");

   std::vector<TLorentzVector> *jetLVecUp = new std::vector<TLorentzVector>;
   std::vector<TLorentzVector> *jetLVecDn = new std::vector<TLorentzVector>;

   std::vector<double> *recoJetsBtagUp = new std::vector<double>;
   std::vector<double> *recoJetsBtagDn = new std::vector<double>;

//      std::vector<double> *dmetMag = new std::vector<double>;
//      std::vector<double> *dmetPhi = new std::vector<double>;

   double metUp = 0.0, metDn = 99999.0;
   double metPhiUp = metphi, metPhiDn = metphi;

   for(int iMet = 0; iMet < metMagUpVec.size(); ++iMet){
      metUp = std::max(metUp, metMagUpVec[iMet]);
      metDn = std::min(metDn, metMagDnVec[iMet]);

      if( TVector2::Phi_mpi_pi(metPhiUpVec[iMet] - metphi) > 0 && metPhiUp < metPhiUpVec[iMet] ) metPhiUp = metPhiUpVec[iMet];
      if( TVector2::Phi_mpi_pi(metPhiUpVec[iMet] - metphi) < 0 && metPhiDn > metPhiUpVec[iMet] ) metPhiDn = metPhiUpVec[iMet];

      if( TVector2::Phi_mpi_pi(metPhiDnVec[iMet] - metphi) > 0 && metPhiUp < metPhiDnVec[iMet] ) metPhiUp = metPhiDnVec[iMet];
      if( TVector2::Phi_mpi_pi(metPhiDnVec[iMet] - metphi) < 0 && metPhiDn > metPhiDnVec[iMet] ) metPhiDn = metPhiDnVec[iMet];

//         dmetMag->push_back((metMagUp[iMet] - met)/met);
//         dmetMag->push_back((metMagDown[iMet] - met)/met);

//         dmetPhi->push_back(TVector2::Phi_mpi_pi(metPhiUp[iMet] - metphi));
//         dmetPhi->push_back(TVector2::Phi_mpi_pi(metPhiDown[iMet] - metphi));
   }

   std::vector<double> tmpjetPtUp, tmpjetPtDn;

   for(int iJet = 0; iJet < jetsLVec.size(); ++iJet){
      tmpjetPtUp.push_back( jetsLVec[iJet].Pt() * (1 + recoJetsJecUnc[iJet]) );
      tmpjetPtDn.push_back( jetsLVec[iJet].Pt() * (1 - recoJetsJecUnc[iJet]) );
   }

   std::vector<size_t> ptIdxUp, ptIdxDn;
   stdindexSort::argsort(tmpjetPtUp.begin(), tmpjetPtUp.end(), std::greater<double>(), ptIdxUp);
   stdindexSort::argsort(tmpjetPtDn.begin(), tmpjetPtDn.end(), std::greater<double>(), ptIdxDn);
   for(unsigned int ip=0; ip<ptIdxUp.size(); ip++){
      unsigned int idxMapped = ptIdxUp[ip];
      jetLVecUp->push_back( jetsLVec[idxMapped] * (1 + recoJetsJecUnc[idxMapped]) );
      recoJetsBtagUp->push_back( recoJetsBtag[idxMapped] );
   }

   for(unsigned int ip=0; ip<ptIdxDn.size(); ip++){
      unsigned int idxMapped = ptIdxDn[ip];
      jetLVecDn->push_back( jetsLVec[idxMapped] * (1 - recoJetsJecUnc[idxMapped]) );
      recoJetsBtagDn->push_back( recoJetsBtag[idxMapped] );
   }

   tr.registerDerivedVar("met_metMagUp", metUp);
   tr.registerDerivedVar("met_metMagDn", metDn);

   tr.registerDerivedVar("metphi_metPhiUp", metPhiUp);
   tr.registerDerivedVar("metphi_metPhiDn", metPhiDn);

   tr.registerDerivedVec("jetLVec_jecUp", jetLVecUp);
   tr.registerDerivedVec("jetLVec_jecDn", jetLVecDn);

   tr.registerDerivedVec("recoJetsBtag_jecUp", recoJetsBtagUp);
   tr.registerDerivedVec("recoJetsBtag_jecDn", recoJetsBtagDn);
}
