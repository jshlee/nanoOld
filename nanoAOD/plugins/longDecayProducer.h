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

bool isKsFrom(const reco::GenParticle* particle, int pdgId, int count,int & KsFromQuark, bool & KsFromTquark){
  if(abs(particle->pdgId()) == pdgId && particle->status() == 62){
    cout << " pdg : " << particle->pdgId() << " , status : " << particle->status() << endl;
    auto dau1 = particle->daughter(0);
    auto dau2 = particle->daughter(1);
    auto ndau = particle->numberOfDaughters();

    cout << " dau1 pdg : " << dau1->pdgId() << " , status : " << dau1->status() << endl;
    cout << " dau2 pdg : " << dau2->pdgId() << " , status : " << dau2->status() << endl;
    cout << " ndau : " << ndau << endl;
    KsFromQuark = dau2->pdgId();
    return KsFromTquark = true;
  }
  const reco::GenParticleRefVector& mothers = particle->motherRefVector(); 
/*
  if(mothers.empty()) cout << "[" << count << "] empty!" << endl;
  if(mothers.isAvailable()) cout << "[" << count << "] Available!" << endl;
  if(mothers.isNull()) cout << "[" << count << "] Null!" << endl;
  if(mothers.isTransient()) cout << "[" << count << "] Transient!" << endl;

  cout << "[" << count << "] mothers capacity : " << mothers.capacity() << endl;
*/
  count = count + 1;
  for(reco::GenParticleRefVector::const_iterator im = mothers.begin(); im != mothers.end(); ++im){
    const reco::GenParticle& part = **im;
//    cout << "[" << count << "] Momther pdgId : " << part.pdgId() << endl;
    if( isKsFrom( &part, pdgId, count, KsFromQuark, KsFromTquark) ){
      return KsFromTquark = true;
    }
  }
  return KsFromTquark = false;
}

void motherTracking(const TrackingVertex trackVertex, const TrackingParticle *decayTrk, int count, int & KsFromQuark, bool & KsFromTquark, std::vector<int> & isKsFromTsb, vector<uint8_t> & isKsFromTop){
  cout << " ["<< count << "] decayTrk = " << " pdg = " << decayTrk->pdgId() << ", pt = " << decayTrk->p4().Pt() <<", vert = "<< trackVertex.position() <<", inVolume "<<trackVertex.inVolume()<< endl;
  cout << " ["<< count << "] isGenEmpty = " << decayTrk->genParticles().empty() << endl;
  cout << " ["<< count << "] isGenNull = " << decayTrk->genParticles().isNull() << endl;
  cout << " ["<< count << "] isGenAvail = " << decayTrk->genParticles().isAvailable() << endl;
  cout << " ["<< count << "] Capacity = " << decayTrk->genParticles().capacity() << endl;

  if(!decayTrk->genParticles().empty()){
    for(TrackingParticle::genp_iterator igen = decayTrk->genParticle_begin(); igen != decayTrk->genParticle_end(); ++igen){
      auto gen = igen->get();
      if(count != 0 || decayTrk->pdgId() == 310) {
        isKsFromTop.push_back(isKsFrom(gen,6,count,KsFromQuark,KsFromTquark));
        isKsFromTsb.push_back(KsFromQuark);
      }
    }
  }
  else{
    count = count + 1;
    cout << "there is no gen info ===> go parent vertex...(decayTrk->pv)" << endl;
    auto pv = decayTrk->parentVertex().get();
    for (TrackingVertex::tp_iterator pr = pv->sourceTracks_begin(); pr != pv->sourceTracks_end(); ++pr) {
      auto decayTrk2 = pr->get();
      motherTracking(*pv, decayTrk2, count, KsFromQuark, KsFromTquark, isKsFromTsb, isKsFromTop);
    }
  }
}
