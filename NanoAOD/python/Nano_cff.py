import FWCore.ParameterSet.Config as cms

def customiseMuons(process):
    # additional variables needed for h2mu
    process.muonTable.variables.globalMu = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('isGlobalMuon'),
        expr = cms.string("isGlobalMuon()"),
        mcOnly = cms.bool(False),
        precision = cms.int32(-1),
        type = cms.string('bool')
    )
    
    return(process)
    
def customise(process):
    customiseMuons(process)
    
    process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
    process.load('cat.NanoAOD.cmesons_cff')
    process.nanoAOD_step += process.cmesonTables
    
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
    process.MessageLogger.cerr.FwkSummary.reportEvery = cms.untracked.int32(1000)
    fileName = cms.untracked.string('nanoAOD.root')
    if hasattr(process, 'NANOAODSIMoutput'):          
        process.NANOAODSIMoutput.fileName = fileName
    else:
        process.NANOAODoutput.fileName = fileName
    
    return(process)
