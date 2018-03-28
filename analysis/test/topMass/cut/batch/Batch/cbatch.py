#!/usr/bin/env python 
# catGetDatasetInfo v7-4-4 # to make dataset lists
# sed -i 's/^\/store/root:\/\/cms-xrdr.sdfarm.kr:1094\/\/xrd\/store/g' *

import os, json, array
import numpy as np
from math import ceil       
username = os.environ['USER']

analysis = 'copt'
#analysis = 'TtbarDiLeptonAnalyzer'
pythonCfg = 'copt.py'
#analysis=analysis+'Silver'
RunFiles = [
              "hadron",
              ]
datadir = os.environ["CMSSW_BASE"]+'/src/nano/analysis/data/dataset/'
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
        subBatch = "condor_submit -batch-name %s -append 'arguments=%s %s' copt_cfg.jds" %(datasetName ,datasetName,FileNamesStr)
        #print createbatch
        print subBatch 
            
        os.system(subBatch)
