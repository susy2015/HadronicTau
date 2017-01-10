#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$2/src/SusyAnaTools/Tools/obj:$2/src/TopTagger/TopTagger/test:$2/src/opencv/lib

cd $2/src
eval `scramv1 runtime -sh`

cd ${_CONDOR_SCRATCH_DIR}

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$2/src/opencv/lib/

xrdcp root://cmseos.fnal.gov/$(echo $5 | sed 's|/eos/uscms||') .

./HadTauLL $1 -1 $3 $4 "condor"

rm $(echo $5 | sed 's|.*/||')
