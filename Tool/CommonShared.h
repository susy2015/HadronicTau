#ifndef COMMONSHARED_H
#define COMMONSHARED_H

#include "SusyAnaTools/Tools/baselineDef.h"
#include "SusyAnaTools/Tools/PDFUncertainty.h"

namespace systBaseline
{
   const std::string spec_jecUp = "jecUp";
   BaselineVessel * SRblv_jecUp =0;

   const std::string spec_jecDn = "jecDn";
   BaselineVessel * SRblv_jecDn =0;

   const std::string spec_metMagUp = "metMagUp";
   BaselineVessel * SRblv_metMagUp =0;

   const std::string spec_metMagDn = "metMagDn";
   BaselineVessel * SRblv_metMagDn =0;

   const std::string spec_metPhiUp = "metPhiUp";
   BaselineVessel * SRblv_metPhiUp =0;

   const std::string spec_metPhiDn = "metPhiDn";
   BaselineVessel * SRblv_metPhiDn =0;

   PDFUncertainty * pdfScale =0;
}

namespace predSpec
{
   const std::vector<std::vector<double>> mergeBins_dummy;
// Merge bins to fix the big TF issue with limited MC stat. unc --> merge bins for both mu and ele channels
// Current method is the following:
// Only combine neighbor bins. Do NOT support combine more than 2 bins (if true, why bins are designed so poor to have >=3 combined?)
// In the two bins, assume one bin has good stat. (the reason to combine) and the other has poor stat.
// Then add stat. from both bins to the bin with poor stat. and do NOT touch the bin with good stat.
   const std::vector<std::vector<double>> mergeBins = {
                                                        {18, 19},
                                                        {34, 35}
                                                      };
   void adjustBins_merge(TH1 * hist)
   {
      if( mergeBins.empty() ) return;
 
// no checking if nXbins is consistent with total search bins or not...
      const int nXbins = hist->GetXaxis()->GetNbins();
      for(auto pair : mergeBins)
      {
         double sum =0, sumErr =0;
         bool OKtoModif = true;
         int toModif_binIdx = -1;
         for(auto perbin : pair)
         {
            if( perbin > nXbins-1 ){ OKtoModif = false; break; }// no way to do since bin index to merge out-of-range
            if( toModif_binIdx == -1 )
            {
               toModif_binIdx = perbin;
            }else
            {
               if( hist->GetBinContent(perbin+1) < hist->GetBinContent(toModif_binIdx+1) ) toModif_binIdx = perbin;
            }
            sum += hist->GetBinContent(perbin+1); // bin counting starts from 0
            sumErr += hist->GetBinError(perbin+1) * hist->GetBinError(perbin+1);
         }
         sumErr = sqrt(sumErr);
         if( OKtoModif )
         {
            hist->SetBinContent(toModif_binIdx+1, sum);
            hist->SetBinError(toModif_binIdx+1, sumErr);
         }
      }
   }
}
#endif
