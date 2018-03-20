#ifndef HadronProducer_H
#define HadronProducer_H
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

#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "TLorentzVector.h"

//#define debugMode

using namespace edm;
using namespace std;

class HadronProducer : public edm::stream::EDProducer<> {
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

public:
  struct hadronCandidate {
    reco::VertexCompositeCandidate vcc;
    pat::Jet jet;
    int nJet, nDau;
    float diffMass, lxy, lxySig, l3D, l3DSig, dca, angleXY, angleXYZ;
  };
  
  explicit HadronProducer(const edm::ParameterSet & iConfig);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  void produce( edm::Event&, const edm::EventSetup& ) override;

  reco::VertexCompositeCandidate fit(vector<reco::Candidate*>& cands,
				     reco::Vertex& pv, int pdgId,
				     float &dca, float &angleXY, float &angleXYZ);
    
  SVector3 getDistanceVector(int dim, reco::VertexCompositeCandidate& vertex,reco::Vertex& pv);
  pair<float, float> getDistance(int dim, reco::VertexCompositeCandidate& vertex,reco::Vertex& pv);
  
  vector<hadronCandidate> findJPsiCands(vector<reco::Candidate*> &leptons, reco::Vertex& pv, int nJet, const pat::Jet & aPatJet);
  vector<hadronCandidate> findD0Cands(vector<reco::Candidate*> &chargedHads, reco::Vertex& pv, int nJet, const pat::Jet & aPatJet);
  vector<hadronCandidate> findDStarCands(vector<HadronProducer::hadronCandidate>& d0cands, vector<reco::Candidate*> &chargedHads,
					 reco::Vertex& pv, int nJet, const pat::Jet & aPatJet);
  edm::EDGetTokenT<edm::View<pat::Jet> > jetLabel_;
  edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
  edm::ESHandle<TransientTrackBuilder> trackBuilder_;

  const int pion_pdgId_ = 211, kaon_pdgId_ = 321, proton_pdgId_ = 2122;
  const float pion_m_ = 0.1396, kaon_m_ = 0.4937, proton_m_ = 0.938272;

  const int jpsi_pdgId_ = 443, d0_pdgId_ = 421, dstar_pdgId_ = 413;
  const float jpsi_m_ = 3.096, d0_m_ = 1.865, dstar_m_ = 2.010;

  const int kshort_pdgId_ = 310, lambda_pdgId_ = 3122;
  const float kshort_m_ = 0.4976, lambda_m_ = 1.11568;
  
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
  float cosThetaXYCut_;
  float cosThetaXYZCut_;  
};


HadronProducer::HadronProducer(const edm::ParameterSet & iConfig) :
  jetLabel_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetLabel"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel")))
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
  cosThetaXYCut_ = iConfig.getParameter<double>("cosThetaXYCut");
  cosThetaXYZCut_ = iConfig.getParameter<double>("cosThetaXYZCut");
  
  produces<nanoaod::FlatTable>("had");
  produces<reco::VertexCompositeCandidateCollection>();
  produces<vector<pat::Jet>>("jet");
}

reco::VertexCompositeCandidate HadronProducer::fit(vector<reco::Candidate*>& cands,
						   reco::Vertex& pv, int pdgId,
						   float &dca, float &angleXY, float &angleXYZ)
{
  int charge = 0;
  vector<reco::TransientTrack> transientTracks;
  for (auto trk : cands){
    //if (trk->bestTrack() == nullptr) continue;
    const reco::TransientTrack transientTrack = trackBuilder_->build(trk->bestTrack());
    transientTracks.emplace_back(transientTrack);
    charge += trk->charge();
  }

  if (transientTracks.size() < 2) return reco::VertexCompositeCandidate();

  // impactPointTSCP DCA
  // measure distance between tracks at their closest approach
  dca = -1;
  ClosestApproachInRPhi cApp;
  reco::TransientTrack tt1 = transientTracks[0], tt2 = transientTracks[1];
  cApp.calculate(tt1.impactPointTSCP().theState(), tt2.impactPointTSCP().theState());  
  if (!cApp.status()) return reco::VertexCompositeCandidate();
  
  dca = cApp.distance();

  if (dca > tkDCACut_) return reco::VertexCompositeCandidate();
  
  // the POCA should at least be in the sensitive volume
  GlobalPoint cxPt = cApp.crossingPoint();
  if (sqrt(cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y()) > 120. || abs(cxPt.z()) > 300.)
    return reco::VertexCompositeCandidate();

  // the tracks should at least point in the same quadrant
  TrajectoryStateClosestToPoint const & posTSCP = tt1.trajectoryStateClosestToPoint(cxPt);
  TrajectoryStateClosestToPoint const & negTSCP = tt2.trajectoryStateClosestToPoint(cxPt);
  if (!posTSCP.isValid() || !negTSCP.isValid()) return reco::VertexCompositeCandidate();
  if (posTSCP.momentum().dot(negTSCP.momentum()) < 0) return reco::VertexCompositeCandidate();
  
  KalmanVertexFitter m_kvf(true);
  
  TransientVertex tv = m_kvf.vertex(transientTracks);
  if (!tv.isValid()) return reco::VertexCompositeCandidate();
  
  reco::Vertex theVtx = tv;
  // loose cut on chi2
  if (theVtx.normalizedChi2() > vtxChi2Cut_) return reco::VertexCompositeCandidate();

  // if (theVtx.normalizedChi2() < 0)
  //   std::cout << "-ve vtx: " << theVtx.normalizedChi2() << " chi2/ndof" << theVtx.chi2() << " / " << theVtx.ndof()
  // 	      << ". trknc2 " << tt1.normalizedChi2() << " " << tt2.normalizedChi2() << " dca:" << dca << " " << tv.trackWeight(tt1) << " " << tv.trackWeight(tt2)
  // 	      << ". (x,y,z) " << cxPt.x() << ", " << cxPt.y() << ", " << cxPt.z() << " " << theVtx.x() << ", " << theVtx.y() << ", " << theVtx.z()
  // 	      << std::endl;
  // else
  //   std::cout << "+ve vtx: " << theVtx.normalizedChi2() << " chi2/ndof" << theVtx.chi2() << " / " << theVtx.ndof()
  // 	      << ". trknc2 " << tt1.normalizedChi2() << " " << tt2.normalizedChi2() << " dca:" << dca << " " << tv.trackWeight(tt1) << " " << tv.trackWeight(tt2)
  // 	      << ". (x,y,z) " << cxPt.x() << ", " << cxPt.y() << ", " << cxPt.z() << " " << theVtx.x() << ", " << theVtx.y() << ", " << theVtx.z()
  // 	      << std::endl;    
  
  GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());

  math::XYZTLorentzVector tlv;
  GlobalVector totalP;
  int i = 0;
  for (auto trk : tv.refittedTracks()){
    TrajectoryStateClosestToPoint const & tscp = trk.trajectoryStateClosestToPoint(vtxPos);
    GlobalVector mom = tscp.momentum();
    double mass = cands[i]->mass();
    double energy = sqrt(mom.mag2() + mass*mass);    
    const math::XYZTLorentzVector lv(mom.x(), mom.y(), mom.z(), energy);
    totalP += mom;
    tlv += lv;
    ++i;
  }
  math::XYZPoint referencePos = pv.position();

  // 2D pointing angle
  double dx = theVtx.x()-referencePos.x();
  double dy = theVtx.y()-referencePos.y();
  double px = totalP.x();
  double py = totalP.y();
  angleXY = (dx*px+dy*py)/(sqrt(dx*dx+dy*dy)*sqrt(px*px+py*py));
  if (angleXY > cosThetaXYCut_) return reco::VertexCompositeCandidate();
  
  // 3D pointing angle
  double dz = theVtx.z()-referencePos.z();
  double pz = totalP.z();
  angleXYZ = (dx*px+dy*py+dz*pz)/(sqrt(dx*dx+dy*dy+dz*dz)*sqrt(px*px+py*py+pz*pz));
  if (angleXYZ > cosThetaXYZCut_) return reco::VertexCompositeCandidate();

  reco::Particle::Point vtx(theVtx.x(), theVtx.y(), theVtx.z());
  const reco::Vertex::CovarianceMatrix vtxCov(theVtx.covariance());

  reco::VertexCompositeCandidate secVert(charge, tlv, vtx, vtxCov, theVtx.chi2(), theVtx.ndof(), pdgId);
  for (auto trk : cands){    
    secVert.addDaughter(*trk);   
  }
  return secVert;
}


HadronProducer::SVector3
HadronProducer::getDistanceVector(int dim, reco::VertexCompositeCandidate& vertex,reco::Vertex& pv)
{
  float z = 0.;
  if (dim == 3) z = vertex.vz() - pv.position().z();
  SVector3 distanceVector(vertex.vx() - pv.position().x(),
			  vertex.vy() - pv.position().y(),
			  z);
  return distanceVector;
}

pair<float, float> HadronProducer::getDistance(int dim, reco::VertexCompositeCandidate& vertex,reco::Vertex& pv)
{
  SMatrixSym3D totalCov = vertex.vertexCovariance() + pv.covariance();
  SVector3 distVecXYZ = getDistanceVector(dim, vertex, pv);
  float distMagXYZ = ROOT::Math::Mag(distVecXYZ);
  float sigmaDistMagXYZ = sqrt(ROOT::Math::Similarity(totalCov, distVecXYZ)) / distMagXYZ;

  return make_pair(distMagXYZ,sigmaDistMagXYZ);
}

void
HadronProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HadronProducer);
#endif
