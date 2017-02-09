#!/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_8/external/slc6_amd64_gcc491/bin/python
####!${SRT_CMSSW_RELEASE_BASE_SCRAMRTDEL}/external/${SCRAM_ARCH}/bin/python

from samples import SampleCollection
from os import system, environ
import re
import optparse 

submitFile = """universe = vanilla
Executable = $ENV(CMSSW_BASE)/src/HadronicTau/Tool/condorHadTauLL/goMakePlots.sh
Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $ENV(CMSSW_BASE)/src/HadronicTau/Tool/HadTauLL, $ENV(CMSSW_BASE)/src/HadronicTau/Tool/condorHadTauLL/goMakePlots.sh, $ENV(CMSSW_BASE)/lib/$ENV(SCRAM_ARCH)/librecipeAUXOxbridgeMT2.so, $ENV(CMSSW_BASE)/src/HadronicTau/Tool/TopTagger.cfg, $ENV(CMSSW_BASE)/src/HadronicTau/Tool/TrainingOutput_dR20_pt30_depth14_2016_Dec2.model, $ENV(CMSSW_BASE)/src/SusyAnaTools/Tools/CSVFiles/CSVv2_Moriond17_B_H.csv, $ENV(CMSSW_BASE)/src/SusyAnaTools/Tools/data/allINone_bTagEff.root, $ENV(CMSSW_BASE)/src/SusyAnaTools/Tools/ISR_Root_Files/allINone_ISRJets.root, $ENV(CMSSW_BASE)/src/SusyAnaTools/Tools/ISR_Root_Files/ISRWeights.root, $ENV(CMSSW_BASE)/src/SusyAnaTools/Tools/LeptonSF_Root_Files/allINone_leptonSF_Moriond17.root, $ENV(CMSSW_BASE)/src/SusyAnaTools/Tools/puppiSoftdropResol.root, TRANS_TAR_BALL_CMSSW
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

rel_base = environ['CMSSW_BASE']
rel_name = rel_base.split('/')[-1]
tar_file_name = rel_name+'.tar.gz'

cache_all_dir = rel_base+'/src/SusyAnaTools/Tools/cache_all.sh'
system('sh '+cache_all_dir)

tar_command = 'tar --exclude-caches-all -czf '+tar_file_name + ' -C '+rel_base+'/.. ' + rel_name
system(tar_command)
system('mv '+rel_name+'.tar.gz '+rel_base+'/src/SusyAnaTools/Tools/condor/')

submitFile = re.sub(r'TRANS_TAR_BALL_CMSSW', '$ENV(CMSSW_BASE)/src/SusyAnaTools/Tools/condor/'+tar_file_name, submitFile)

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
                fileParts.append("Arguments = %s %s %i %i %s\nQueue\n\n"%(n, rel_name, startFileNum, nFilesPerJob, s))
            f.close()

fout = open("condorHadTauLL_submit.txt", "w")
fout.write(''.join(fileParts))
fout.close()

system('mkdir -p logs')
system("echo 'condor_submit condorHadTauLL_submit.txt'")
system('condor_submit condorHadTauLL_submit.txt')
