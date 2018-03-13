import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

##################### Tables for final output and docs ##########################
ldTable = cms.EDProducer("longDecayProducer",
  jetLabel = cms.InputTag("slimmedJets"),
  vertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
  mcLabel  = cms.InputTag("prunedGenParticles"),
  packed = cms.InputTag("packedGenParticles"),
  pruned = cms.InputTag("prunedGenParticles"),
  applySoftLeptonCut = cms.bool(True),
  doFullMCMatching = cms.bool(False),
  # -- cuts on initial track collection --
  # Track normalized Chi2 <
  tkChi2Cut = cms.double(100),
  # Number of valid hits on track >=
  tkNHitsCut = cms.int32(0),
  # Pt of track >
  tkPtCut = cms.double(0.),
  # Track impact parameter significance >
  tkIPSigXYCut = cms.double(100000),
  tkIPSigZCut = cms.double(100000),
  
  # -- cuts on the vertex --
  # Vertex chi2 <
  vtxChi2Cut = cms.double(10),
  # XY decay distance significance >
  vtxDecaySigXYCut = cms.double(-1),
  # XYZ decay distance significance >
  vtxDecaySigXYZCut = cms.double(-1.),
  
  # -- miscellaneous cuts --
  # POCA distance between tracks <
  tkDCACut = cms.double(100),
  # cos(angleXY) between x and p of V0 candidate >
  cosThetaXYCut = cms.double(100),
  # cos(angleXYZ) between x and p of V0 candidate >
  cosThetaXYZCut = cms.double(100),
)

longDecayCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("ldTable"),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("longDecay"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(True), 
    variables = cms.PSet(P4Vars,
        x   = Var("vx()", float, doc = "secondary vertex X position, in cm",precision=14),
        y   = Var("vy()", float, doc = "secondary vertex Y position, in cm",precision=14),
        z   = Var("vz()", float, doc = "secondary vertex Z position, in cm",precision=14),
        pdgId=Var("pdgId()", int, doc = "pdgId"),
    ),
)
longDecayCandidateTable.variables.pt.precision=14
longDecayCandidateTable.variables.phi.precision=14
longDecayCandidateTable.variables.eta.precision=14
longDecayCandidateTable.variables.mass.precision=14

ldTables = cms.Sequence(ldTable+longDecayCandidateTable)
