/*
From cmssw RecoVertex/V0Producer/src/V0Fitter.cc
*/

#include<memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

using namespace edm;
using namespace std;
using namespace reco;

class CMesonProducer : public edm::stream::EDProducer<> {
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;
    
public:
  explicit CMesonProducer(const edm::ParameterSet & iConfig);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  void produce( edm::Event&, const edm::EventSetup& ) override;

  reco::VertexCompositeCandidate fit(vector<const pat::PackedCandidate*> cands, int pdgId, float &dca);
  
  float getDistance(int dim, reco::VertexCompositeCandidate vertex,reco::Vertex pv);
  float getDistanceSig(int dim, reco::VertexCompositeCandidate vertex,reco::Vertex pv);
  SVector3 getDistanceVector(int dim, reco::VertexCompositeCandidate vertex,reco::Vertex pv);

  edm::EDGetTokenT<edm::View<pat::Jet> > jetSrc_;
  edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
  edm::ESHandle<TransientTrackBuilder> trackBuilder_;

  const float gPionMass = 0.1396;
  const float gKaonMass = 0.4937;  
  const float jpsiMin_ = 2.5;
  const float jpsiMax_ = 3.4;
  const float D0Min_   = 1.7;
  const float D0Max_   = 2.0;
  const float DstarDiffMin_   = 0.14;
  const float DstarDiffMax_   = 0.16;
  //unsigned int maxNumPFCand_;
  bool applyCuts_;
  // cuts on initial track selection
  float tkChi2Cut_;
  int tkNHitsCut_;
  float tkPtCut_;
  float tkIPSigXYCut_;
  float tkIPSigZCut_;
  // cuts on the vertex
  float vtxChi2Cut_;
  float vtxDecaySigXYCut_;
  float vtxDecaySigXYZCut_;
  // miscellaneous cuts
  float tkDCACut_;
  float innerHitPosCut_;
  float cosThetaXYCut_;
  float cosThetaXYZCut_;  
};

CMesonProducer::CMesonProducer(const edm::ParameterSet & iConfig) :
  jetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetLabel"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  applyCuts_(iConfig.getParameter<bool>("applySoftLeptonCut"))  
{
  // cuts on initial track selection
  tkChi2Cut_ = iConfig.getParameter<double>("tkChi2Cut");
  tkNHitsCut_ = iConfig.getParameter<int>("tkNHitsCut");
  tkPtCut_ = iConfig.getParameter<double>("tkPtCut");
  tkIPSigXYCut_ = iConfig.getParameter<double>("tkIPSigXYCut");
  tkIPSigZCut_ = iConfig.getParameter<double>("tkIPSigZCut");   
  // cuts on vertex
  vtxChi2Cut_ = iConfig.getParameter<double>("vtxChi2Cut");
  vtxDecaySigXYZCut_ = iConfig.getParameter<double>("vtxDecaySigXYZCut");
  vtxDecaySigXYCut_ = iConfig.getParameter<double>("vtxDecaySigXYCut");
  // miscellaneous cuts
  tkDCACut_ = iConfig.getParameter<double>("tkDCACut");
  innerHitPosCut_ = iConfig.getParameter<double>("innerHitPosCut");
  cosThetaXYCut_ = iConfig.getParameter<double>("cosThetaXYCut");
  cosThetaXYZCut_ = iConfig.getParameter<double>("cosThetaXYZCut");
  
  produces<nanoaod::FlatTable>("cmeson");
  produces<reco::VertexCompositeCandidateCollection>();
}

void
CMesonProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);
  reco::Vertex pv = recVtxs->at(0);

  Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByToken(jetSrc_, jetHandle);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder_);

  auto cmVxCand = std::make_unique<reco::VertexCompositeCandidateCollection>();
  std::vector<float> dca,lxy,l3D,jetDR,legDR,diffMass;
  dca.clear();lxy.clear();l3D.clear();jetDR.clear();legDR.clear();diffMass.clear();
  std::vector<int> pid;pid.clear();
  
  for (const pat::Jet & aPatJet : *jetHandle){
    if (aPatJet.pt()< 30 or abs(aPatJet.eta())>3 ) continue;

    vector< const pat::PackedCandidate * > jetDaughters, softlepCands;
        
    for( unsigned int idx = 0 ; idx < aPatJet.numberOfDaughters() ; ++idx) {
      const pat::PackedCandidate * dauCand ( dynamic_cast<const pat::PackedCandidate*>(aPatJet.daughter(idx)));
      if ( dauCand->charge() == 0 ) continue;
      //if ( dauCand->pt() <1 ) continue;
      jetDaughters.emplace_back( dauCand );
      if ( abs(dauCand->pdgId()) == 11  || abs(dauCand->pdgId())==13) softlepCands.emplace_back(dauCand);
    }
    unsigned int dau_size = jetDaughters.size();
    if ( dau_size < 2 ) continue;

    // dont need to sort since all combinations are done.
    //sort(jetDaughters.begin(), jetDaughters.end(), [](const pat::PackedCandidate * a, const pat::PackedCandidate * b) {return a->pt() > b->pt(); }); 

    // if ( dau_size > maxNumPFCand_ ) dau_size = maxNumPFCand_;
    // jetDaughters.resize( dau_size );
    for ( unsigned int lep1_idx = 0 ; lep1_idx< dau_size-1 ; ++lep1_idx) {
      const pat::PackedCandidate* lep1Cand = jetDaughters[lep1_idx];
      if ( abs(lep1Cand->pdgId()) != 13 and abs(lep1Cand->pdgId()) != 11) continue;
      
      for ( unsigned int lep2_idx = lep1_idx+1 ; lep2_idx< dau_size ; ++lep2_idx) {
        const pat::PackedCandidate* lep2Cand = jetDaughters[lep2_idx];
	if ( abs(lep2Cand->pdgId()) != 13 and abs(lep2Cand->pdgId()) != 11) continue;

        int pdgMul = lep1Cand->pdgId() * lep2Cand->pdgId();
        if ( pdgMul != -121 and pdgMul != -169 ) continue; 

	float dca_Jpsi = -9;
	reco::VertexCompositeCandidate JpsiCand = this->fit({lep1Cand, lep2Cand}, 443, dca_Jpsi);
	
        // if ( abs(lep1Cand->pdgId() ) == 13 && abs(lep2Cand->pdgId()) == 13 ) {
        //   int lep1ID = (int)lep1Cand->isStandAloneMuon() + (int)lep1Cand->isGlobalMuon()*2;
        //   int lep2ID = (int)lep2Cand->isStandAloneMuon() + (int)lep2Cand->isGlobalMuon()*2;
        //   JpsiCand.setLeptonID( lep1ID, lep2ID );
        // }
        // else JpsiCand.setLeptonID( -1, -1 );
  
        if ( JpsiCand.mass() < jpsiMin_ || JpsiCand.mass() > jpsiMax_ ) continue;

	cmVxCand->emplace_back(JpsiCand);      
	dca.emplace_back(dca_Jpsi);
	lxy.emplace_back(getDistance(2,JpsiCand,pv));
	l3D.emplace_back(getDistance(3,JpsiCand,pv));
	jetDR.emplace_back(reco::deltaR( JpsiCand, aPatJet));
	legDR.emplace_back(reco::deltaR( *lep1Cand, *lep2Cand));
	diffMass.emplace_back(0);
	pid.emplace_back(443);
	
      }
    }

    if ( applyCuts_ && softlepCands.size()==0  ) continue;
    
    for ( unsigned int pion_idx = 0 ; pion_idx< dau_size ; ++pion_idx) {
      const pat::PackedCandidate * pionCand = jetDaughters[pion_idx];
      if ( abs(pionCand->pdgId()) == 13 or abs(pionCand->pdgId()) == 11) continue;
      
      for ( unsigned int kaon_idx = 0 ; kaon_idx< dau_size ; ++kaon_idx) {
        if ( pion_idx == kaon_idx ) continue;
	pat::PackedCandidate kaonCand(*jetDaughters[kaon_idx]);	
        if ( abs(kaonCand.pdgId()) == 13 or abs(kaonCand.pdgId()) == 11) continue;
        if ( pionCand->charge() * kaonCand.charge() != -1 ) continue;
	
        kaonCand.setMass(gKaonMass);
	
	float dca_D0 = -9;
	reco::VertexCompositeCandidate D0Cand = fit({pionCand, &kaonCand}, 421, dca_D0);
	
        if ( applyCuts_ ) D0Cand.addDaughter( *softlepCands[0] );

        if ( D0Cand.mass() < D0Min_ || D0Cand.mass() > D0Max_ ) continue;

	cmVxCand->emplace_back(D0Cand);
	dca.emplace_back(dca_D0);
	lxy.emplace_back(getDistance(2,D0Cand,pv));
	l3D.emplace_back(getDistance(3,D0Cand,pv));
	jetDR.emplace_back(reco::deltaR( D0Cand, aPatJet));
	legDR.emplace_back(reco::deltaR( *pionCand, kaonCand));
	diffMass.emplace_back(0);
	pid.emplace_back(421);
	
        if ( dau_size < 3 ) continue;

	for( unsigned int extra_pion_idx = 0 ;  extra_pion_idx < dau_size ; ++extra_pion_idx) {
	  if ( extra_pion_idx== pion_idx || extra_pion_idx == kaon_idx) continue;
	  const pat::PackedCandidate * pion2Cand = jetDaughters[extra_pion_idx];
	  if ( abs(pion2Cand->pdgId()) != 211) continue;
	    
	  float dca_Dstar = -9;
	  reco::VertexCompositeCandidate DstarCand = fit({pionCand,&kaonCand,pion2Cand}, pion2Cand->charge()*413, dca_Dstar);

	  if ( applyCuts_ ) DstarCand.addDaughter( *softlepCands[0] );
	    
	  float diffMass_Dstar = DstarCand.mass() - D0Cand.mass();
	  if ( diffMass_Dstar < DstarDiffMin_ || diffMass_Dstar > DstarDiffMax_ ) continue;
	    
	  cmVxCand->emplace_back(DstarCand);
	  dca.emplace_back(dca_Dstar);
	  lxy.emplace_back(getDistance(2,DstarCand,pv));
	  l3D.emplace_back(getDistance(3,DstarCand,pv));
	  jetDR.emplace_back(reco::deltaR( DstarCand, aPatJet));
	  legDR.emplace_back(reco::deltaR( D0Cand, *pion2Cand));      
	  diffMass.emplace_back(diffMass_Dstar);
	  pid.emplace_back(413);
	      
	}
      }
    }
  }
  
  auto cmesonTable = std::make_unique<nanoaod::FlatTable>(cmVxCand->size(),"cmeson",false);
  // For SV we fill from here only stuff that cannot be created with the SimpleFlatTableProducer 
  cmesonTable->addColumn<float>("lxy",lxy,"2D decay length in cm",nanoaod::FlatTable::FloatColumn,10);
  cmesonTable->addColumn<float>("l3D",l3D,"3D decay length in cm",nanoaod::FlatTable::FloatColumn,10);
  cmesonTable->addColumn<float>("dca",dca,"distance of closest approach cm",nanoaod::FlatTable::FloatColumn,10);
  cmesonTable->addColumn<float>("jetDR",jetDR,"DR between jet",nanoaod::FlatTable::FloatColumn,10);
  cmesonTable->addColumn<float>("legDR",legDR,"DR between leg",nanoaod::FlatTable::FloatColumn,10);
  cmesonTable->addColumn<float>("diffMass",diffMass,"diffMass",nanoaod::FlatTable::FloatColumn,10);
  cmesonTable->addColumn<int>("pid",pid,"pid of vertex cand",nanoaod::FlatTable::IntColumn,8);
 
  iEvent.put(std::move(cmesonTable),"cmeson");
  iEvent.put(std::move(cmVxCand));
  
}

reco::VertexCompositeCandidate CMesonProducer::fit(vector<const pat::PackedCandidate*> cands, int pdgId, float &dca)
{
  int charge = 0;
  vector<reco::TransientTrack> transientTracks;
  math::XYZTLorentzVector lv;
  for (auto trk : cands){
    if (trk->bestTrack() == nullptr) continue;
    const TransientTrack transientTrack = trackBuilder_->build(trk->bestTrack());
    transientTracks.emplace_back(transientTrack);
    charge += trk->charge();
    lv += trk->p4();    
  }

  if (transientTracks.size() < 2) return reco::VertexCompositeCandidate();
  
  KalmanVertexFitter m_kvf(true);
  
  TransientVertex tv = m_kvf.vertex(transientTracks);      
  if (!tv.isValid()) return reco::VertexCompositeCandidate();
  
  reco::Vertex theVtx = tv;
  // loose cut on chi2
  if (theVtx.normalizedChi2() > 10.) return reco::VertexCompositeCandidate();
  
  reco::Particle::Point vtx(theVtx.x(), theVtx.y(), theVtx.z());  
  const reco::Vertex::CovarianceMatrix vtxCov(theVtx.covariance());
    
  // impactPointTSCP DCA
  dca = -1;
  ClosestApproachInRPhi cApp;
  reco::TransientTrack tt1 = tv.refittedTracks()[0], tt2 = tv.refittedTracks()[1];
  cApp.calculate(tt1.impactPointTSCP().theState(), tt2.impactPointTSCP().theState());  
  if (!cApp.status()) return reco::VertexCompositeCandidate();
  
  dca = cApp.distance();

  // the POCA should at least be in the sensitive volume
  GlobalPoint cxPt = cApp.crossingPoint();
  if (sqrt(cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y()) > 120. || std::abs(cxPt.z()) > 300.)
    return reco::VertexCompositeCandidate();

  // the tracks should at least point in the same quadrant
  TrajectoryStateClosestToPoint const & posTSCP = tt1.trajectoryStateClosestToPoint(cxPt);
  TrajectoryStateClosestToPoint const & negTSCP = tt2.trajectoryStateClosestToPoint(cxPt);
  if (!posTSCP.isValid() || !negTSCP.isValid()) return reco::VertexCompositeCandidate();
  if (posTSCP.momentum().dot(negTSCP.momentum()) < 0) return reco::VertexCompositeCandidate();
  
  

  
  reco::VertexCompositeCandidate secVert(charge, lv, vtx, vtxCov, theVtx.chi2(), theVtx.ndof(), pdgId);
  for (auto trk : cands){
    if (trk->bestTrack() == nullptr) continue;
    secVert.addDaughter( *trk );    
  }
  
  return secVert;
}

CMesonProducer::SVector3
CMesonProducer::getDistanceVector(int dim, reco::VertexCompositeCandidate vertex,reco::Vertex pv)
{
  float z = 0.;
  if (dim == 3) z = vertex.vz() - pv.position().z();
  SVector3 distanceVector(vertex.vx() - pv.position().x(),
			  vertex.vy() - pv.position().y(),
			  z);
  return distanceVector;
}

float CMesonProducer::getDistance(int dim, reco::VertexCompositeCandidate vertex,reco::Vertex pv)
{
  SVector3 distanceVec = getDistanceVector(dim, vertex, pv);
  return ROOT::Math::Mag(distanceVec);
}

float CMesonProducer::getDistanceSig(int dim, reco::VertexCompositeCandidate vertex,reco::Vertex pv)
{
  SMatrixSym3D totalCov = vertex.vertexCovariance() + pv.covariance();
  SVector3 distVecXYZ = getDistanceVector(dim, vertex, pv);
  float distMagXYZ = ROOT::Math::Mag(distVecXYZ);
  float sigmaDistMagXYZ = sqrt(ROOT::Math::Similarity(totalCov, distVecXYZ)) / distMagXYZ;

  return sigmaDistMagXYZ;
}

void
CMesonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(CMesonProducer);
