import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

##################### Tables for final output and docs ##########################
cmesonTable = cms.EDProducer("CMesonProducer",
  jetLabel = cms.InputTag("slimmedJets"),
  vertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
  #mcLabel  = cms.InputTag("prunedGenParticles"),
  applySoftLeptonCut = cms.bool(True),
  # -- cuts on initial track collection --
  # Track normalized Chi2 <
  tkChi2Cut = cms.double(10.),
  # Number of valid hits on track >=
  tkNHitsCut = cms.int32(3),
  # Pt of track >
  tkPtCut = cms.double(0.35),
  # Track impact parameter significance >
  tkIPSigXYCut = cms.double(2.),
  tkIPSigZCut = cms.double(-1.),
  
  # -- cuts on the vertex --
  # Vertex chi2 <
  vtxChi2Cut = cms.double(6.63),
  # XY decay distance significance >
  vtxDecaySigXYCut = cms.double(15.),
  # XYZ decay distance significance >
  vtxDecaySigXYZCut = cms.double(-1.),
  
  # -- miscellaneous cuts --
  # POCA distance between tracks <
  tkDCACut = cms.double(1.),
  # check if either track has a hit radially inside the vertex position minus this number times the sigma of the vertex fit
  # note: Set this to -1 to disable this cut, which MUST be done if you want to run V0Producer on the AOD track collection!
  innerHitPosCut = cms.double(4.),
  # cos(angleXY) between x and p of V0 candidate >
  cosThetaXYCut = cms.double(0.998),
  # cos(angleXYZ) between x and p of V0 candidate >
  cosThetaXYZCut = cms.double(-2.),
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
        ndof= Var("vertexNdof()", float, doc = "number of degrees of freedom",precision=8),
        chi2= Var("vertexNormalizedChi2()", float, doc = "reduced chi2, i.e. chi/ndof",precision=8),
    ),
)
cmesonCandidateTable.variables.pt.precision=10
cmesonCandidateTable.variables.phi.precision=12

#before cross linking
cmesonSequence = cms.Sequence()
#after cross linkining
cmesonTables = cms.Sequence(cmesonTable+cmesonCandidateTable)

