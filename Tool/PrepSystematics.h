#ifndef PREPSYSTEMATICS_H
#define PREPSYSTEMATICS_H

#include <iostream>
#include <vector>

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/baselineDef.h"

class PrepSystematics
{
   public:
   PrepSystematics(const bool isusegenmet = false);

   void operator()(NTupleReader& tr)
   {
      systematicPrep(tr);
   }

   private:

   const bool usegenmet;

   void systematicPrep(NTupleReader& tr);
};

#endif
