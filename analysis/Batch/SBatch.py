#!/usr/bin/env python 
# catGetDatasetInfo v7-4-4 # to make dataset lists
# sed -i 's/^\/store/root:\/\/cms-xrdr.sdfarm.kr:1094\/\/xrd\/store/g' *

import os, json, array
import numpy as np
from math import ceil       
username = os.environ['USER']

analysis = 'tth2mu'
#analysis = 'TtbarDiLeptonAnalyzer'
pythonCfg = 'tth2mu.py'
#analysis=analysis+'Silver'
RunFiles = [
              'WMinusH_HToMuMu',
              'WPlusH_HToMuMu',
              'ZH_HToMuMu',
              'VBF_HToMuMu',
              'GG_HToMuMu',
              "WWTo2L2Nu",
              "WZTo3LNu_amcatnlo",
              "WZTo2LQQ",
              "ZZTo2L2Nu",
              "ZZTo2L2Q",
              "ZZTo4L",
              "WWW",
              "WWZ",
              "WZZ",
              "ZZZ",
#              "ttZToLLNuNu",
#              "ttWToLNu",
              "SingleTop_tW_noHadron",
              "SingleTbar_tW_noHadron",
              "SingleTop_tW",
              "SingleTbar_tW",
#              "TTJets_DiLept",
#              "TTJets_DiLept_Tune4",
              'TTJets_aMC', 
              'DYJets',
#              'DYJets_MG_10to50',
#              'DYJets_MG2',
#              'DYJets_2J',
#              'DYJets_1J',
#              'DYJets_0J',
#              'DYJets_10to50', 
              'SingleMuon_Run2016B',
              'SingleMuon_Run2016C',
              'SingleMuon_Run2016D',
              'SingleMuon_Run2016E',
              'SingleMuon_Run2016F',
              'SingleMuon_Run2016G',
              'SingleMuon_Run2016H',
            #  'SingleMuon_Run2016H_v3',
              ]
datadir = os.environ["CMSSW_BASE"]+'/nano/nanoAOD/data/dataset/'
#version = os.environ["CMSSW_VERSION"]


count = 0
for i in RunFiles:
    datasetName = i
    fileList = datadir + 'dataset_' + datasetName + '.txt'
    files = np.array([])
    for f in open(fileList).readlines():
        f = f.strip()
        f = f.strip('\',"')
        if len(f) < 5: continue
        if '#' == f[0] or '.root' != f[-5:]: continue
        files = np.append(files,[f])
    nFiles = len(files)     
    maxFiles = 20
    nSection = int(ceil(1.0*nFiles/maxFiles))
    for section in range(nSection):
        begin = section*maxFiles
        end = min(begin + maxFiles, nFiles)
        FileNames = files[begin:end]
        FileNamesStr = " ".join(str(i) for i in FileNames)
        #jobName = analysis+'_'+datasetName
        subBatch = "condor_submit -batch-name %s -append 'arguments=%s %s' tth2mu_cfg.jds" %(datasetName ,datasetName,FileNamesStr)
        #print createbatch
        print subBatch 
            
        os.system(subBatch)
