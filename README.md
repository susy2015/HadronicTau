Check out analysis recipe as described in https://github.com/susy2015/SusyAnaTools

#HadronicTau code
```
cd $CMSSW_BASE/src
git clone -b Moriond2017 git@github.com:susy2015/HadronicTau.git
scram b -j9
cd HadronicTau/Tool/
make
```
#To check the control sample distribution, run CS

##Running interactively
```
./CS TTbarSingleLepT -1 0 1
```
TTbarSingleLepT is sample name, -1 is the number of events to be run, 0 and 1 are start and end root files respectively.

It will create a root file with naming convention TTbarSingleLepT_CS0.root(o for start file) which contains all the histograms.

##Running in condor
```
cd condorCS
voms-proxy-init
./condorSubmit.py -d TTbarSingleLep -n 10
```
TTbarSingleLep: sample set, 10 is number of files per job.

#To check the hadtau and lostlepton yield distribution, run HadTauLL
```
./HadTauLL TTbarSingleLepT -1 0 1
```
Condor jobs can be submitted in similar way as CS. One similar condor directory needs to be created.