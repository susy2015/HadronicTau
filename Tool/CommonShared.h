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
#endif
