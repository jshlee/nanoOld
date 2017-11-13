import FWCore.ParameterSet.Config as cms

def customise(process):
    process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
    process.load('cat.NanoAOD.cmesons_cff')
    process.nanoAOD_step += process.cmesonTables
    
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
    process.MessageLogger.cerr.FwkSummary.reportEvery = cms.untracked.int32(1000)
    process.NANOAODSIMoutput.fileName = cms.untracked.string('nanoAOD.root')
    
    return(process)
