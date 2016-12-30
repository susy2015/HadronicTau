#!/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_8/external/slc6_amd64_gcc491/bin/python
####!${SRT_CMSSW_RELEASE_BASE_SCRAMRTDEL}/external/${SCRAM_ARCH}/bin/python

from samples import SampleCollection
from os import system
import optparse 

submitFile = """universe = vanilla
Executable = $ENV(CMSSW_BASE)/src/HadronicTau/Tool/condorCS/goMakePlots.sh
Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $ENV(CMSSW_BASE)/src/HadronicTau/Tool/CS, $ENV(CMSSW_BASE)/src/HadronicTau/Tool/condorCS/goMakePlots.sh, $ENV(CMSSW_BASE)/lib/$ENV(SCRAM_ARCH)/librecipeAUXOxbridgeMT2.so, $ENV(CMSSW_BASE)/src/HadronicTau/Tool/TopTagger.cfg, $ENV(CMSSW_BASE)/src/HadronicTau/Tool/TrainingOutput_dR20_pt30_depth14_2016_Dec2.model
Output = logs/makePlots_$(Process).stdout
Error = logs/makePlots_$(Process).stderr
Log = logs/makePlots_$(Process).log
notify_user = ${LOGNAME}@FNAL.GOV
x509userproxy = $ENV(X509_USER_PROXY)
"""

parser = optparse.OptionParser("usage: %prog [options]\n")

parser.add_option ('-n', dest='numfile', type='int', default = 5, help="number of files per job")
parser.add_option ('-d', dest='datasets', type='string', default = '', help="List of datasets 'TTbarSingleLep'")
#parser.add_option ('-f', dest='subsample', type='string', default = '', help="Nmae of subsample 'WJetsToLNu_HT_2500toInf'")


options, args = parser.parse_args()

nFilesPerJob = options.numfile
#subsamplename = options.subsample

fileParts = [submitFile]

sc = SampleCollection()

datasets = []

if options.datasets:
    datasets = options.datasets.split(',')
else:
    print "No dataset pecified"
    exit(0)

for ds in datasets:
    ds = ds.strip()

    for s, n in sc.sampleList(ds):
        print n
        print s
       # if n!= subsamplename:
           # continue
        f = open(s)
        if not f == None:
            count = 0
            for l in f:
                if '.root' in l and not 'failed' in l:
                    count = count + 1
            for startFileNum in xrange(0, count, nFilesPerJob):
                fileParts.append("Arguments = %s $ENV(CMSSW_BASE) %i %i %s\nQueue\n\n"%(n, startFileNum, nFilesPerJob, s))
            f.close()

fout = open("condorCS_submit.txt", "w")
fout.write(''.join(fileParts))
fout.close()

system('mkdir -p logs')
system("echo 'condor_submit condorCS_submit.txt'")
system('condor_submit condorCS_submit.txt')
