#!/bin/bash 

RESULT_DIR="NanoAOD"
DIR=`hostname`

/bin/mkdir -p $RESULT_DIR/$DIR

echo "hostname" >> $RESULT_DIR/$DIR/`hostname`.txt
echo "date" >> $RESULT_DIR/$DIR/`hostname`.txt
source /cvmfs/cms.cern.ch/cmsset_default.sh >> $RESULT_DIR/$DIR/`hostname`.txt
python /cms/scratch/daniel/CMSSW_8_0_26_patch1/src/CATTools/CatAnalyzer/test/Nano_AOD/tth2mu.py $@
echo "Done." >> $RESULT_DIR/$DIR/`hostname`.txt

