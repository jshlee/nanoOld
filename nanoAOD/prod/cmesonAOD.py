# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: run2_2016MC -s NANO -n -1 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --filein file:ttbar_mc.root --conditions auto:run2_mc --era Run2_2016,run2_miniAOD_80XLegacy --customise nano/nanoAOD/nano_cff.customise
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('NANO',eras.Run2_2016,eras.run2_miniAOD_80XLegacy)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.maxEvents.input = cms.untracked.int32(1000)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:aod_ttbar.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('run2_2016MC nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('nanoAOD.root'),
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppies
makePuppies( process );
process.load('PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff')
process.load('PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi')
process.load('PhysicsTools.PatAlgos.slimming.packedPFCandidates_cfi')
process.load('PhysicsTools.PatAlgos.slimming.slimmedJets_cfi')
process.packedPFCandidates.inputVertices = cms.InputTag("offlinePrimaryVertices")
process.packedPFCandidates.PuppiNoLepSrc = cms.InputTag("puppi")
process.slimmedJets.src = cms.InputTag("patJets")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('nano.nanoAOD.cmesons_cff')
process.cmesonTable.vertexLabel = cms.InputTag("offlinePrimaryVertices")
process.cmesonTable.mcLabel  = cms.InputTag("genParticles")

process.p = cms.Path(process.makePatJets
                         +process.primaryVertexAssociation+process.puppi
                         +process.packedPFCandidates
                         +process.slimmedJets
                         +process.cmesonTables+process.cmesonCandidateTable)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.p,process.NANOAODSIMoutput_step)
