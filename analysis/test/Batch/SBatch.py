#!/usr/bin/env python 
# catGetDatasetInfo v7-4-4 # to make dataset lists
# sed -i 's/^\/store/root:\/\/cms-xrdr.sdfarm.kr:1094\/\/xrd\/store/g' *

import os, json, array, sys
import numpy as np
from math import ceil       
username = os.environ['USER']

analysis = 'TT'

#pythonCfg = 'tth2mu.py'
pythonCfg = 'tth2muC.py'

RunFiles = [
    #         # 'WMinusH_HToMuMu',
    #         # 'WPlusH_HToMuMu',
    #         # 'ZH_HToMuMu',
    #         # 'VBF_HToMuMu',
    #         # 'GG_HToMuMu',
    #           'ttH',
    #         # "WWTo2L2Nu",
    #         # "WZTo3LNu_amcatnlo",
    #         # "WZTo2LQQ",
    #         # "ZZTo2L2Nu",
    #         # "ZZTo2L2Q",
    #          "ZZTo4L_powheg",
    #          "WWW",
    #          "WWZ",
    #          "WZZ",
    #          "ZZZ",
    #          "TTZToLLNuNu",
              "ZZ",
              "WZ",
              "WW",
#   #           "ttWToLNu",
    #         # "SingleTop_tW_noHadron",
    #         # "SingleTbar_tW_noHadron",
    #         # "SingleTop_tW",
    #         # "SingleTbar_tW",
#   #           "TTJets_DiLept",
#   #           "TTJets_DiLept_Tune4",
    #          'TTJets_aMC', 
    #          'DYJets',
#   #           'DYJets_MG_10to50',
#   #           'DYJets_MG2',
#   #           'DYJets_2J',
#   #           'DYJets_1J',
#   #           'DYJets_0J',
#   #           'DYJets_10to50', 
    #          'SingleMuon_Run2016B',
    #          'SingleMuon_Run2016Bv1',
    #          'SingleMuon_Run2016C',
    #          'SingleMuon_Run2016D',
    #          'SingleMuon_Run2016E',
    #          'SingleMuon_Run2016F',
    #          'SingleMuon_Run2016G',
    #          'SingleMuon_Run2016H',
    #         # 'SingleMuon_Run2016H_v3',
              ]
SetDir = "test"              
datadir = '{}/src/nano/analysis/data/dataset/'.format(os.environ['CMSSW_BASE'])
#version = os.environ["CMSSW_VERSION"]


for i in RunFiles:
    datasetName = i
    fileList = datadir + 'dataset_' + datasetName + '.txt'
    jobName = analysis+'_'+datasetName 

    Dirname = "{}/src/nano/analysis/test/Batch/{}/".format(os.environ['CMSSW_BASE'],jobName)
    DirnameJDS = "{}/src/nano/analysis/test/Batch/".format(os.environ['CMSSW_BASE'])
    if os.path.isdir(Dirname):
        print "ERROR: output directory already existing."
        sys.exit()
    else: os.makedirs(Dirname)

    Dirname_ = "{}/src/nano/analysis/test/Results/{}/{}/".format(os.environ['CMSSW_BASE'],SetDir,datasetName)
    if not os.path.isdir(Dirname_):
        os.makedirs(Dirname_)

    Dirname_ = "%s/src/nano/analysis/test/Batch/NanoAOD/"%(os.environ['CMSSW_BASE'])
    if not os.path.isdir(Dirname_):
        os.makedirs(Dirname_)

    files = np.array([])
    for f in open(fileList).readlines():
        f = f.strip()
        f = f.strip('\',"')
        if len(f) < 5: continue
        if '#' == f[0] or '.root' != f[-5:]: continue
        files = np.append(files,[f])
    nFiles = len(files)     
    maxFiles = 10
    nSection = int(ceil(1.0*nFiles/maxFiles))
    count = 0
    for section in range(nSection):
        begin = section*maxFiles
        end = min(begin + maxFiles, nFiles)
        FileNames = files[begin:end]
        FileNamesStr = " ".join(str(i) for i in FileNames)

        print "@@ Writing run script..."
        jds = "%ssubmit.jds" %Dirname 
        fout = open(jds, "w")
        print>>fout, "# Job description file for condor job"
        print>>fout, """executable = {0}/src/nano/analysis/test/Batch/tth2mu.sh
universe   = vanilla

log = condor.log

getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = {1}/job_{2}.log
error = {1}/job_{2}.err
transfer_input_files = NanoAOD
queue""" .format(os.environ['CMSSW_BASE'],Dirname, count)
        fout.close()
        count += 1 
        #jobName = analysis+'_'+datasetName
        subBatch = "condor_submit -batch-name {} -append 'arguments={} {} {}' {}".format(datasetName ,SetDir ,datasetName ,FileNamesStr ,jds)
        #print createbatch
        print subBatch 
            
        os.system(subBatch)
