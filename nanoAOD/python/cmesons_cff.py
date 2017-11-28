import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

##################### Tables for final output and docs ##########################
cmesonTable = cms.EDProducer("CMesonProducer",
  jetLabel = cms.InputTag("slimmedJets"),
  vertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
  mcLabel  = cms.InputTag("prunedGenParticles"),
  applySoftLeptonCut = cms.bool(True),
  # -- cuts on initial track collection --
  # Track normalized Chi2 <
  tkChi2Cut = cms.double(100),
  # Number of valid hits on track >=
  tkNHitsCut = cms.int32(0),
  # Pt of track >
  tkPtCut = cms.double(0.),
  # Track impact parameter significance >
  tkIPSigXYCut = cms.double(100),
  tkIPSigZCut = cms.double(100),
  
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

cmesonCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("cmesonTable"),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("cmeson"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(True), 
    variables = cms.PSet(P4Vars,
        x   = Var("vx()", float, doc = "secondary vertex X position, in cm",precision=10),
        y   = Var("vy()", float, doc = "secondary vertex Y position, in cm",precision=10),
        z   = Var("vz()", float, doc = "secondary vertex Z position, in cm",precision=14),
        chi2= Var("vertexChi2()", float, doc = "chi2",precision=10),
        ndof= Var("vertexNdof()", int, doc = "number of degrees of freedom",precision=8),
        pdgId=Var("pdgId()", int, doc = "pdgId",precision=8),
    ),
)
cmesonCandidateTable.variables.pt.precision=10
cmesonCandidateTable.variables.phi.precision=12

#before cross linking
cmesonSequence = cms.Sequence()
#after cross linkining
cmesonTables = cms.Sequence(cmesonTable+cmesonCandidateTable)

