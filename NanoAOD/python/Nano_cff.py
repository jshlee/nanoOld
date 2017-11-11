import FWCore.ParameterSet.Config as cms

def customise(process):
    process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
    process.load('cat.NanoAOD.cmesons_cff')
    process.nanoAOD_step += process.cmesonTables
    return(process)
