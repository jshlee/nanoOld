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
  void beginStream(edm::StreamID) override {};
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override {};

  reco::VertexCompositeCandidate fit(vector<const pat::PackedCandidate*> cands, int pdgId, float &dca);
  
  float get2Ddistance(reco::VertexCompositeCandidate vertex,reco::Vertex pv);
  float get3Ddistance(reco::VertexCompositeCandidate vertex,reco::Vertex pv);
  
  edm::EDGetTokenT<edm::View<pat::Jet> > jetSrc_;
  edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
  edm::ESHandle<TransientTrackBuilder> trackBuilder_;

  const float gJpsiMass = 3.096;
  const float gPionMass = 0.1396;
  const float gKaonMass = 0.4937;
  const float gD0Mass   = 1.86480;
  const float gDstarD0DiffMass = 0.145;
  float d0MassWindow_, maxDeltaR_, d0MassCut_;
  unsigned int maxNumPFCand_;
  bool applyCuts_;
  

};
//bool pTComp( const reco::Candidate* a, const reco::Candidate* b) { return a->pt() > b->pt();  }

CMesonProducer::CMesonProducer(const edm::ParameterSet & iConfig) :
  jetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetLabel"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel")))
{
  produces<nanoaod::FlatTable>("cmeson");
  produces<reco::VertexCompositeCandidateCollection>();
  
  maxNumPFCand_ = iConfig.getParameter<int>("maxNumPFCand");
  d0MassWindow_ = iConfig.getParameter<double>("d0MassWindow");
  d0MassCut_ = iConfig.getParameter<double>("d0MassCut");
  maxDeltaR_  = iConfig.getParameter<double>("maxDeltaR");
  applyCuts_ = iConfig.getParameter<bool>("applySoftLeptonCut");
}
void
CMesonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
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
    vector<TransientTrack> tracks;

    reco::VertexCompositeCandidate bestJpsi, bestD0, bestDstar;
    
    float dca_Jpsi=-9, legDR_Jpsi;
    float dca_D0=-9, legDR_D0;
    float dca_Dstar=-9, legDR_Dstar, diffMass_Dstar;
    
    for( unsigned int idx = 0 ; idx < aPatJet.numberOfDaughters() ; idx++) {
      const pat::PackedCandidate * dauCand ( dynamic_cast<const pat::PackedCandidate*>(aPatJet.daughter(idx)));
      if ( dauCand->charge() ==0 ) continue;
      //if ( dauCand->pt() <1 ) continue;
      jetDaughters.push_back( dauCand );
      if ( abs(dauCand->pdgId()) == 11  || abs(dauCand->pdgId())==13) softlepCands.push_back(dauCand);
    }
    unsigned int dau_size = jetDaughters.size();
    if ( dau_size < 2 ) continue;

    sort(jetDaughters.begin(), jetDaughters.end(), [](const pat::PackedCandidate * a, const pat::PackedCandidate * b) {return a->pt() > b->pt(); }); 

    // if ( dau_size > maxNumPFCand_ ) dau_size = maxNumPFCand_;
    // jetDaughters.resize( dau_size );
    for ( unsigned int lep1_idx = 0 ; lep1_idx< dau_size-1 ; lep1_idx++) {
      for ( unsigned int lep2_idx = lep1_idx+1 ; lep2_idx< dau_size ; lep2_idx++) {
	
        const pat::PackedCandidate* lep1Cand = jetDaughters[lep1_idx];
        const pat::PackedCandidate* lep2Cand = jetDaughters[lep2_idx];

        int pdgMul = lep1Cand->pdgId() * lep2Cand->pdgId();
        if ( pdgMul != -121 && pdgMul != -169 ) continue; 

	vector<const pat::PackedCandidate*> jpsiCands;
	jpsiCands.push_back(lep1Cand);
	jpsiCands.push_back(lep2Cand);
	
	float dca = 0;
	reco::VertexCompositeCandidate JpsiCand = fit(jpsiCands, 443, dca);
	
        // if ( abs(lep1Cand->pdgId() ) == 13 && abs(lep2Cand->pdgId()) == 13 ) {
        //   int lep1ID = (int)lep1Cand->isStandAloneMuon() + (int)lep1Cand->isGlobalMuon()*2;
        //   int lep2ID = (int)lep2Cand->isStandAloneMuon() + (int)lep2Cand->isGlobalMuon()*2;
        //   JpsiCand.setLeptonID( lep1ID, lep2ID );
        // }
        // else JpsiCand.setLeptonID( -1, -1 );
  
        if ( JpsiCand.mass() < 2 || JpsiCand.mass() > 4 ) continue;
        if ( abs(bestJpsi.mass() - gJpsiMass) > abs(JpsiCand.mass()-gJpsiMass) ){
	  bestJpsi = JpsiCand;
	  dca_Jpsi = dca;
	  legDR_Jpsi = reco::deltaR( *lep1Cand, *lep2Cand);
	}
      }
    }
    if ( dca_Jpsi > -1 ){
      cmVxCand->push_back(bestD0);      
      dca.push_back(dca_Jpsi);
      lxy.push_back(get2Ddistance(bestJpsi,pv));
      l3D.push_back(get3Ddistance(bestJpsi,pv));
      jetDR.push_back(reco::deltaR( bestJpsi, aPatJet));
      legDR.push_back(legDR_Jpsi);
      diffMass.push_back(0);
      pid.push_back(443);
    }
    if ( applyCuts_ && softlepCands.size()==0  ) continue;
    
    for ( unsigned int pion_idx = 0 ; pion_idx< dau_size ; pion_idx++) {
      for ( unsigned int kaon_idx = 0 ; kaon_idx< dau_size ; kaon_idx++) {
        if ( pion_idx == kaon_idx ) continue;
        const pat::PackedCandidate * pionCand = jetDaughters[pion_idx];
	pat::PackedCandidate kaonCand(*jetDaughters[kaon_idx]);	
        kaonCand.setMass(gKaonMass);

        if ( abs(pionCand->pdgId()) == 13 or abs(pionCand->pdgId()) == 11) continue;
        if ( abs(kaonCand.pdgId()) == 13 or abs(kaonCand.pdgId()) == 11) continue;
        if ( pionCand->charge() * kaonCand.charge() != -1 ) continue;

        auto D0 = pionCand->p4() + kaonCand.p4();
        if ( abs(D0.M() - gD0Mass) > d0MassCut_) continue;

	vector<const pat::PackedCandidate*> d0Cands;
	d0Cands.push_back(pionCand);
	d0Cands.push_back(&kaonCand);
	float dca = 0;
	reco::VertexCompositeCandidate D0Cand = fit(d0Cands, 421, dca);
	
        if ( applyCuts_ ) D0Cand.addDaughter( *softlepCands[0] );

        if ( abs(bestD0.mass() - gD0Mass)> abs(D0Cand.mass()-gD0Mass) ){
	  bestD0 = D0Cand;
	  dca_D0 = dca;
	  legDR_D0 = reco::deltaR( *pionCand, kaonCand);
	}

        if ( dau_size < 3 ) continue;

        if ( abs( D0.M() - gD0Mass) < d0MassWindow_ ) {
          for( unsigned int extra_pion_idx = 0 ;  extra_pion_idx < dau_size ; extra_pion_idx++) {
            if ( extra_pion_idx== pion_idx || extra_pion_idx == kaon_idx) continue;
            const pat::PackedCandidate * pion2Cand = jetDaughters[extra_pion_idx];
            if ( abs(pion2Cand->pdgId()) != 211) continue;

	    vector<const pat::PackedCandidate*> dstarCands;
	    dstarCands.push_back(pionCand);
	    dstarCands.push_back(&kaonCand);
	    dstarCands.push_back(pion2Cand);
	    
	    float dcaDstar = 0;
	    reco::VertexCompositeCandidate DstarCand = fit(dstarCands, pion2Cand->charge()*413, dcaDstar);

            if ( applyCuts_ ) D0Cand.addDaughter( *softlepCands[0] );
	    
            float diffMass = DstarCand.mass() - D0.M();
            if ( diffMass > 0.2) continue;
            if ( abs( diffMass_Dstar - gDstarD0DiffMass) > abs( diffMass - gDstarD0DiffMass) ){
	      bestDstar = DstarCand;
	      dca_Dstar = dcaDstar;
	      legDR_Dstar = reco::deltaR( D0Cand, *pion2Cand);
	      diffMass_Dstar = diffMass;
	    }
          }
        }
      }
    }
    if ( dca_D0 > -1 ){
      cmVxCand->push_back(bestD0);
      dca.push_back(dca_D0);
      lxy.push_back(get2Ddistance(bestD0,pv));
      l3D.push_back(get3Ddistance(bestD0,pv));
      jetDR.push_back(reco::deltaR( bestD0, aPatJet));
      legDR.push_back(legDR_D0);
      diffMass.push_back(0);
      pid.push_back(421);
    }
    if ( dca_Dstar > -1 ){
      cmVxCand->push_back(bestDstar);
      dca.push_back(dca_Dstar);
      lxy.push_back(get2Ddistance(bestDstar,pv));
      l3D.push_back(get3Ddistance(bestDstar,pv));
      jetDR.push_back(reco::deltaR( bestDstar, aPatJet));
      legDR.push_back(legDR_Dstar);      
      diffMass.push_back(diffMass_Dstar);
      pid.push_back(413);
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
    transientTracks.push_back(transientTrack);
    charge += trk->charge();
    lv += trk->p4();    
  }

  if (transientTracks.size() < 2) return reco::VertexCompositeCandidate();
  
  KalmanVertexFitter m_kvf(true);
  
  TransientVertex tv = m_kvf.vertex(transientTracks);      
  if (!tv.isValid()) return reco::VertexCompositeCandidate();
  
  reco::Vertex theVtx = tv;
  reco::Particle::Point vtx(theVtx.x(), theVtx.y(), theVtx.z());  
  const reco::Vertex::CovarianceMatrix vtxCov(theVtx.covariance());
  
  // for (auto refittedTrk : tv.refittedTracks()) {
  //   lv += refittedTrk.track().p4();
  // }
  reco::VertexCompositeCandidate secVert(charge, lv, vtx, vtxCov, theVtx.chi2(), theVtx.ndof(), pdgId);
  for (auto trk : cands){
    if (trk->bestTrack() == nullptr) continue;
    secVert.addDaughter( *trk );    
  }  
  
  // impactPointTSCP DCA
  dca = -1;
  ClosestApproachInRPhi cApp_impactPointTSCP;
  reco::TransientTrack tt1 = tv.refittedTracks()[0], tt2 = tv.refittedTracks()[1];
  cApp_impactPointTSCP.calculate(tt1.impactPointTSCP().theState(), tt2.impactPointTSCP().theState());  
  if ( cApp_impactPointTSCP.status() )
    dca = cApp_impactPointTSCP.distance();
  
  return secVert;
}

float CMesonProducer::get2Ddistance(reco::VertexCompositeCandidate vertex,reco::Vertex pv)
{
  SVector3 distanceVectorXY(vertex.vx() - pv.position().x(),
			    vertex.vy() - pv.position().y(), 0.);
  return ROOT::Math::Mag(distanceVectorXY);
}
float CMesonProducer::get3Ddistance(reco::VertexCompositeCandidate vertex,reco::Vertex pv)
{ 
  SVector3 distanceVector3D(vertex.vx() - pv.position().x(),
			    vertex.vy() - pv.position().y(),
			    vertex.vz() - pv.position().z());
  return ROOT::Math::Mag(distanceVector3D);
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
