import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

##################### Tables for final output and docs ##########################

cmesonTable = cms.EDProducer("CMesonProducer",
  jetLabel = cms.InputTag("slimmedJets"),
  vertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
  #mcLabel  = cms.InputTag("prunedGenParticles"),
  maxNumPFCand = cms.int32(9999),
  maxDeltaR = cms.double(0.2),
  d0MassCut = cms.double(0.5),
  d0MassWindow = cms.double(0.05),
  applySoftLeptonCut = cms.bool(True)
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
        ndof   = Var("vertexNdof()", float, doc = "number of degrees of freedom",precision=8),
        chi2   = Var("vertexNormalizedChi2()", float, doc = "reduced chi2, i.e. chi/ndof",precision=8),
    ),
)
cmesonCandidateTable.variables.pt.precision=10
cmesonCandidateTable.variables.phi.precision=12

#before cross linking
cmesonSequence = cms.Sequence()
#after cross linkining
cmesonTables = cms.Sequence(cmesonTable+cmesonCandidateTable)

