#!/usr/bin/env python
import os,json,sys,shutil,time,getopt
#import CATTools.CatProducer.catDefinitions_cfi as cat

def submitjob(requestName, dataset, globalTag, lumiMask, submit):
    print 'v'*80
    print "creating job"
    print dataset

    isMC = False
    if 'MC' in requestName:
        isMC = True
    
    dataset = dataset.strip()
    if isMC :
        label = dataset.split("/")[1]
    else :
        label = dataset.split("/")[1]+"_"+dataset.split("/")[2]

    if requestName:
        dataRequestName = '%s_%s'%(requestName,label)
        outputDatasetTag = '%s_%s'%(requestName,dataset.split("/")[2])

    if submit:
        sendjob = "crab submit Data.outputDatasetTag='%s' JobType.psetName='%s.py' General.requestName='%s' Data.inputDataset='%s'"%(outputDatasetTag,requestName,dataRequestName,dataset)
    else :
        sendjob = "crab submit --dryrun Data.outputDatasetTag='%s' JobType.psetName='%s.py' General.requestName='%s' Data.inputDataset='%s'"%(outputDatasetTag,requestName,dataRequestName,dataset)

    print sendjob
    print "submiting job"
    os.system(sendjob)
    print '^'*80

    time.sleep(5)
    
submitBlock = None
requestName = ""
datasets = []
inputFile =None
submit = False
lumiMask =""
globalTag =""

try:
    opts, args = getopt.getopt(sys.argv[1:],"hsi:n:l:g:b:",["requestName","inputFile","lumiMask","globalTag","submitBlock"])
except getopt.GetoptError:          
    print 'Usage : ./submitCrab3.py -n <requestName> -i <inputFile> -l <lumiMask> -g <globalTag>'
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./submitCrab3.py -n <requestName> -i <inputFile> -l <lumiMask> -g <globalTag>'
        sys.exit()
    elif opt in ("-n", "--requestName"):
        requestName = arg
    elif opt in ("-s"):
        submit = True
    elif opt in ("-i", "--inputFile"):
        inputFile = arg
        if os.path.isfile(inputFile):
            lines = open(inputFile)
            datasets = lines.readlines()
        else:
            datasets.append(inputFile)
    elif opt in ("-l", "--lumiMask"):
        lumiMask = arg
    elif opt in ("-b", "--submitBlock"):
        submitBlock = arg
    elif opt in ("-g", "--globalTag"):
        globalTag = arg

if requestName == "" :
    print "requestName(-n) is mandantory"
    sys.exit(-1)

if inputFile is None:
    datasets = json.load(open("%s/src/cat/NanoAOD/data/dataset/dataset.json"%os.environ['CMSSW_BASE']))
    for d in datasets:
        dataset = d['DataSetName']
        if len( dataset ) == 0: continue

        isMC = False
        if 'MC' in requestName:
            isMC = True
            if d['type'] == 'Data':
                continue            
        if 'RD' in requestName:
            if d['type'] != 'Data':
                continue            

        if os.path.exists('crab_%s_%s' % (requestName, dataset.split('/')[1])): continue
        if os.path.exists('crab_%s_%s_%s' % (requestName, dataset.split('/')[1], dataset.split('/')[2])): continue
        if len( d['path']) == 0:
            #print d['path'], len( d['path'])
            submitjob(requestName, dataset, None,None, submit)
        
        #if submitBlock == '1' and 'QCD' in dataset:
        #    continue
        #if submitBlock == '2' and 'QCD' not in dataset:
        #    continue

else:
    for dataset in datasets:
        if len(dataset) < 10:
            continue
        if dataset.startswith("#"):
            continue
        submitjob(requestName, dataset, globalTag, lumiMask, submit)

if not submit:
    print "Dry run, not submitting job and only printing crab3 command"
    print "Add -s to submit job"
    print 'Usage : ./submitCrab3.py -n <requestName> -i <inputFile> -s'
