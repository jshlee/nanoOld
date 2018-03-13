/*
  From cmssw RecoVertex/V0Producer/src/V0Fitter.cc
  matching based on Validation/RecoVertex/src/V0Validator.cc
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

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "TLorentzVector.h"

//#define debugMode

using namespace edm;
using namespace std;

class longDecayProducer : public edm::stream::EDProducer<>
{
public:
  explicit longDecayProducer(const edm::ParameterSet & iConfig);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  void produce( edm::Event&, const edm::EventSetup& ) override;

  reco::LeafCandidate getCandidate(const TrackingParticle* tp) {
    return reco::LeafCandidate( tp->charge(), tp->p4(), tp->vertex(), tp->pdgId(), tp->status() );
  };
  
  edm::EDGetTokenT<edm::View<pat::Jet> > jetLabel_;
  edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genLabel_;
  edm::EDGetTokenT<TrackingVertexCollection> trackingVertexLabel_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleLabel_;
  
  edm::ESHandle<TransientTrackBuilder> trackBuilder_;

  const float gPionMass = 0.1396;
  const float gKaonMass = 0.4937;
  const float gJpsiMass = 3.096;
  const float gD0Mass = 1.865;
  const float gDstarMass = 2.010;
  const float gProtonMass = 0.938272;
  const int pdgId_Kp = 321;
  const int pdgId_Jpsi = 443;
  const int pdgId_D0 = 421;
  const int pdgId_Dstar = 413;
  const int pdgId_KS = 310;
  const int pdgId_Lambda = 3122;
  const int pdgId_p = 2122;
  const float jpsiMin_ = 2.5;
  const float jpsiMax_ = 3.4;
  const float D0Min_   = 1.7;
  const float D0Max_   = 2.0;
  const float DstarDiffMin_   = 0.14;
  const float DstarDiffMax_   = 0.16;
  const float KSMin_   = 0.43;
  const float KSMax_   = 0.57;
  const float LambdaMin_   = 0.9;
  const float LambdaMax_   = 1.16;
  //unsigned int maxNumPFCand_;
  bool doFullMatch_;
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
  float cosThetaXYCut_;
  float cosThetaXYZCut_;  
};


longDecayProducer::longDecayProducer(const edm::ParameterSet & iConfig) :
  // jetLabel_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetLabel"))),
  // vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  genLabel_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genLabel"))),
  trackingVertexLabel_(consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("trackingVertexLabel"))),
  trackingParticleLabel_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticleLabel"))),
  // packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
  //mcSrc_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("mcLabel"))),
  doFullMatch_(iConfig.getParameter<bool>("doFullMCMatching")),
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
  cosThetaXYCut_ = iConfig.getParameter<double>("cosThetaXYCut");
  cosThetaXYZCut_ = iConfig.getParameter<double>("cosThetaXYZCut");
  
  produces<nanoaod::FlatTable>("longDecay");
  produces<std::vector<reco::LeafCandidate> >();
}

void
longDecayProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //bool runOnMC = !iEvent.isRealData();
  
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder_);
  
  // Handle<reco::VertexCollection> recVtxs;
  // iEvent.getByToken(vertexLabel_,recVtxs);
  // reco::Vertex pv = recVtxs->at(0);

  // Handle<edm::View<pat::Jet> > jetHandle;
  // iEvent.getByToken(jetLabel_, jetHandle);

  auto candidates = make_unique<std::vector<reco::LeafCandidate>>();
  vector<int> imother;
  
  edm::Handle<TrackingVertexCollection> trackingVertexs;
  iEvent.getByToken(trackingVertexLabel_, trackingVertexs);
  edm::Handle<TrackingParticleCollection> trackingParticles;
  iEvent.getByToken(trackingParticleLabel_, trackingParticles);

  for (auto const& trackVertex : *trackingVertexs.product()) {
    if (trackVertex.eventId().bunchCrossing() != 0)
      continue;  // Consider only in-time events
    if (trackVertex.nDaughterTracks() < 2) continue;  // Keep only V0 vertices
    
    for (TrackingVertex::tp_iterator source = trackVertex.sourceTracks_begin();
         source != trackVertex.sourceTracks_end(); ++source) {
      auto decayTrk = source->get();
      if (decayTrk->pdgId() != 310) continue;

      candidates->push_back(getCandidate(decayTrk));
      imother.push_back(1);
      
      cout << " # decayTrk = " << " pdg = " << decayTrk->pdgId() << ", pt = " << decayTrk->p4().Pt()<<", vert = "<< trackVertex.position() <<", inVolume "<<trackVertex.inVolume()<< endl;
      for (unsigned int i = 0; i < trackVertex.nDaughterTracks(); ++i){	
	auto dau = trackVertex.daughterTracks().at(i).get();

	candidates->push_back(getCandidate(dau));
	imother.push_back(1);
	
	cout << " dau = " << " pdg = " << dau->pdgId() << ", pt = " << dau->p4().Pt() << endl;
      
      }
    }
  }

  Handle<edm::View<reco::GenParticle> > genParticles;
  iEvent.getByToken(genLabel_,genParticles);
  int i = 0;
  for (auto gen : *genParticles){
    if (gen.pdgId() != 5122) continue;    
    cout << " # i = " << i<< " pdg = " << gen.pdgId() << ", pt = " << gen.pt()
  	 << ", status = " << gen.status()
  	 << ", daus = " << gen.numberOfDaughters()
  	 << endl;
    candidates->push_back(gen);
    imother.push_back(1);
    
    for (size_t j = 0; j < gen.numberOfDaughters(); ++j){
      auto dau = gen.daughterRefVector().at(j);
      candidates->push_back(*dau);
      imother.push_back(1);
      cout << " dau = " << i<< " pdg = " << dau->pdgId()
	   << endl;
    }
    ++i;
  }
  

  auto hadTable = make_unique<nanoaod::FlatTable>(candidates->size(),"longDecay",false);
  hadTable->addColumn<int>("mother",imother,"index of mother",nanoaod::FlatTable::IntColumn);
  
  iEvent.put(move(hadTable),"longDecay");
  iEvent.put(move(candidates));
}

void
longDecayProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(longDecayProducer);
