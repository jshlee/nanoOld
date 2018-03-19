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

#include "longDecayProducer.h"
//#define debugMode

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
  vector<int> isKsFromTsb;
  vector<uint8_t> inVol;
  vector<uint8_t> isKsFromTop;

  edm::Handle<TrackingVertexCollection> trackingVertexs;
  iEvent.getByToken(trackingVertexLabel_, trackingVertexs);
  edm::Handle<TrackingParticleCollection> trackingParticles;
  iEvent.getByToken(trackingParticleLabel_, trackingParticles);

  for (auto const& trackVertex : *trackingVertexs.product()) {

    if (trackVertex.eventId().bunchCrossing() != 0) continue;  // Consider only in-time events
    if (trackVertex.nDaughterTracks() < 2) continue;  // Keep only V0 vertices
    
    for (TrackingVertex::tp_iterator source = trackVertex.sourceTracks_begin(); source != trackVertex.sourceTracks_end(); ++source) {
      auto decayTrk = source->get();
      if (decayTrk->pdgId() != 310) continue; //&& decayTrk->pdgId() != 3122) continue;

      cout << "##################################start##################################" << endl;

      candidates->push_back(getCandidate(decayTrk));
      imother.push_back(1);
      inVol.push_back(trackVertex.inVolume());

      int count = 0;
      int KsFromQuark = 0;
      bool KsFromTquark = false;

      motherTracking(trackVertex, decayTrk, count, KsFromQuark, KsFromTquark, isKsFromTsb, isKsFromTop);

      // Check KS-Pion Matching through SimTrk
      for (unsigned int i = 0; i < trackVertex.nDaughterTracks(); ++i){	
	auto dau = trackVertex.daughterTracks().at(i).get();

	candidates->push_back(getCandidate(dau));
	imother.push_back(0);
        inVol.push_back(trackVertex.inVolume());
        isKsFromTop.push_back(KsFromTquark);
	isKsFromTsb.push_back(KsFromQuark);
	cout << "ev final ===> KsFromTquark : " << KsFromTquark << " KsFromQuark : " << KsFromQuark << endl;
        cout << "ev final size : " << "candidate : " << candidates->size() << " imother : " << imother.size() << " inVol : " << inVol.size() << " isKsFromTop : " << isKsFromTop.size() << " isKsFromTsb : " << isKsFromTsb.size() << endl;


//	cout << " dau = " << " pdg = " << dau->pdgId() << ", pt = " << dau->p4().Pt() << endl;
      
      }
      cout << "###################################end###################################" << endl;
    }
  }
/*
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
  
*/
  auto hadTable = make_unique<nanoaod::FlatTable>(candidates->size(),"longDecay",false);
  hadTable->addColumn<int>("mother",imother,"index of mother",nanoaod::FlatTable::IntColumn);
  hadTable->addColumn<int>("isKsFromTsb",isKsFromTsb,"track from t->s/b",nanoaod::FlatTable::IntColumn);
  hadTable->addColumn<uint8_t>("isKsFromTop",isKsFromTop,"track from top",nanoaod::FlatTable::UInt8Column);
  hadTable->addColumn<uint8_t>("inVol",inVol,"track in volume",nanoaod::FlatTable::UInt8Column); 
 
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
