#ifndef HadTruthProducer_H
#define HadTruthProducer_H
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

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "TLorentzVector.h"

using namespace edm;
using namespace std;

class HadTruthProducer : public edm::stream::EDProducer<>
{
public:
  explicit HadTruthProducer(const edm::ParameterSet & iConfig);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce( edm::Event&, const edm::EventSetup& ) override;

  reco::LeafCandidate getCandidate(const TrackingParticle* tp) {
    return reco::LeafCandidate( tp->charge(), tp->p4(), tp->vertex(), tp->pdgId(), tp->status() );
  };
  bool isKsFrom(const reco::GenParticle* particle, int pdgId, int count,int & KsFromQuark, bool & KsFromTop);
  bool isHadFrom(const reco::GenParticleRef &particle, int pdgId, int count, int & hadFromQuark, bool & hadFromTop);
  void motherTracking(const TrackingVertex trackVertex, const TrackingParticle *decayTrk, int count, int & KsFromQuark, bool & KsFromTop, std::vector<int> & isKsFromTsb, vector<bool> & isKsFromTop);

  int trackingVertex_pdgId(const TrackingVertex* tv);
  const reco::GenParticleRef getMother(const TrackingParticleRef& tp);

  edm::EDGetTokenT<reco::RecoToSimCollection> recoRecoToSim_;
  edm::EDGetTokenT<reco::SimToRecoCollection> recoSimToReco_;
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection > hadronCands_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genLabel_;
  edm::EDGetTokenT<TrackingVertexCollection> trackingVertexLabel_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleLabel_;
};

bool HadTruthProducer::isKsFrom(const reco::GenParticle* particle, int pdgId, int count,int & KsFromQuark, bool & KsFromTop){
  if(abs(particle->pdgId()) == pdgId && particle->status() == 62){
//    cout << " pdg : " << particle->pdgId() << " , status : " << particle->status() << endl;
//    auto ndau = particle->numberOfDaughters();
//    auto dau1 = particle->daughter(0);
    auto dau2 = particle->daughter(1);

//    cout << " dau1 pdg : " << dau1->pdgId() << " , status : " << dau1->status() << endl;
//    cout << " dau2 pdg : " << dau2->pdgId() << " , status : " << dau2->status() << endl;
//    cout << " ndau : " << ndau << endl;
    KsFromQuark = dau2->pdgId();
    return KsFromTop = true;
  }
  const reco::GenParticleRefVector& mothers = particle->motherRefVector(); 
/*
  if(mothers.empty()) //cout << "[" << count << "] empty!" << endl;
  if(mothers.isAvailable()) //cout << "[" << count << "] Available!" << endl;
  if(mothers.isNull()) //cout << "[" << count << "] Null!" << endl;
  if(mothers.isTransient()) //cout << "[" << count << "] Transient!" << endl;

  //cout << "[" << count << "] mothers capacity : " << mothers.capacity() << endl;
*/
  count = count + 1;
  for(reco::GenParticleRefVector::const_iterator im = mothers.begin(); im != mothers.end(); ++im){
    const reco::GenParticle& part = **im;
//    //cout << "[" << count << "] Momther pdgId : " << part.pdgId() << endl;
    if( isKsFrom( &part, pdgId, count, KsFromQuark, KsFromTop) ){
      return KsFromTop = true;
    }
  }
  return KsFromTop = false;
}

bool HadTruthProducer::isHadFrom(const reco::GenParticleRef& particle, int pdgId, int count,int & hadFromQuark, bool & hadFromTop){
  if(abs(particle->pdgId()) == pdgId && particle->status() == 62){
//    cout << " pdg : " << particle->pdgId() << " , status : " << particle->status() << endl;
//    auto ndau = particle->numberOfDaughters();
//    auto dau1 = particle->daughter(0);
    auto dau2 = particle->daughter(1);
//    cout << " dau1 pdg : " << dau1->pdgId() << " , status : " << dau1->status() << endl;
//    cout << " dau2 pdg : " << dau2->pdgId() << " , status : " << dau2->status() << endl;
//    cout << " ndau : " << ndau << endl;
    hadFromQuark = dau2->pdgId();
    return hadFromTop = true;
  }
//  cout << count << " pdgId : " << particle->pdgId() << " nMother : " << particle->numberOfMothers() << endl;
  count = count + 1;
  for(unsigned int im = 0; im < particle->numberOfMothers(); ++im){
    const reco::GenParticleRef& mothers = particle->motherRef(im);
//    cout << count << " ==> " << im << " th mother ===> mother " << mothers->numberOfMothers() << endl; 
    if( isHadFrom( mothers, pdgId, count, hadFromQuark, hadFromTop) ){
      return hadFromTop = true;
    }
  }
  return hadFromTop = false;
}

void HadTruthProducer::motherTracking(const TrackingVertex trackVertex, const TrackingParticle *decayTrk, int count, int & KsFromQuark, bool & KsFromTop, std::vector<int> & isKsFromTsb, vector<bool> & isKsFromTop){
  //cout << " ["<< count << "] decayTrk = " << " pdg = " << decayTrk->pdgId() << ", pt = " << decayTrk->p4().Pt() <<", vert = "<< trackVertex.position() <<", inVolume "<<trackVertex.inVolume()<< endl;
  //cout << " ["<< count << "] isGenEmpty = " << decayTrk->genParticles().empty() << endl;
  //cout << " ["<< count << "] isGenNull = " << decayTrk->genParticles().isNull() << endl;
  //cout << " ["<< count << "] isGenAvail = " << decayTrk->genParticles().isAvailable() << endl;
  //cout << " ["<< count << "] Capacity = " << decayTrk->genParticles().capacity() << endl;

  if(!decayTrk->genParticles().empty()){
    for(TrackingParticle::genp_iterator igen = decayTrk->genParticle_begin(); igen != decayTrk->genParticle_end(); ++igen){
      auto gen = igen->get();
      if(count != 0 || decayTrk->pdgId() == 310 || decayTrk->pdgId() == 3122) {
        isKsFromTop.push_back(isKsFrom(gen,6,count,KsFromQuark,KsFromTop));
        isKsFromTsb.push_back(KsFromQuark);
      }
    }
  }
  else{
    count = count + 1;
    //cout << "there is no gen info ===> go parent vertex...(decayTrk->pv)" << endl;
    auto pv = decayTrk->parentVertex().get();
    for (TrackingVertex::tp_iterator pr = pv->sourceTracks_begin(); pr != pv->sourceTracks_end(); ++pr) {
      auto decayTrk2 = pr->get();
      motherTracking(*pv, decayTrk2, count, KsFromQuark, KsFromTop, isKsFromTsb, isKsFromTop);
    }
  }
}

void
HadTruthProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(HadTruthProducer);
#endif
