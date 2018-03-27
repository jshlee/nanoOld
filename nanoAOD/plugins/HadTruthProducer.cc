#include "HadTruthProducer.h"
//#define debugMode

HadTruthProducer::HadTruthProducer(const edm::ParameterSet & iConfig) : recoRecoToSim_(consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("recoRecoToSim"))), recoSimToReco_(consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("recoSimToReco"))),
hadronCands_(consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("hadronCands"))),
genLabel_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genLabel"))),
trackingVertexLabel_(consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("trackingVertexLabel"))),
trackingParticleLabel_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticleLabel")))
{
  produces<nanoaod::FlatTable>("hadTruth");
  produces<nanoaod::FlatTable>("genHadron");
  produces<std::vector<reco::LeafCandidate> >();
}

void
HadTruthProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  Handle<reco::VertexCompositeCandidateCollection> hadronCands;
  iEvent.getByToken(hadronCands_,hadronCands);

  Handle<reco::RecoToSimCollection> recotosimH;
  iEvent.getByToken(recoRecoToSim_, recotosimH);
  auto recotosim = *recotosimH.product();
  Handle<reco::SimToRecoCollection> simtorecoH;
  iEvent.getByToken(recoSimToReco_, simtorecoH);
  auto simtoreco = *simtorecoH.product();
  
  edm::Handle<TrackingVertexCollection> trackingVertexs;
  iEvent.getByToken(trackingVertexLabel_, trackingVertexs);
  edm::Handle<TrackingParticleCollection> trackingParticles;
  iEvent.getByToken(trackingParticleLabel_, trackingParticles);

  vector<int> nmatchedv;

  vector<int> isHadFromTsb;
  vector<bool> isHadFromTop;
  
  for (reco::VertexCompositeCandidateCollection::const_iterator cand = hadronCands->begin();
       cand != hadronCands->end(); cand++) {

    int count = 0;
    int hadFromQuark = 0;
    bool hadFromTop = false;

    // for dstar and lambdaB, need to match with grand mother
//    cout <<"################################################# "<<endl;
//    cout <<"had " << cand->pdgId()<<endl;
//    cout <<"had pt = "<< cand->pt() << ", eta = "<< cand->eta() << ", pid = "<< cand->pdgId() << ", m = "<< cand->mass() <<endl;
    reco::GenParticleRef trueHad;
    
    int numberOfDaughters = cand->numberOfDaughters();
    int nmatched = 0;
    for (int ndau =0; ndau < numberOfDaughters; ++ndau){
      auto rcCand = dynamic_cast<const reco::RecoChargedCandidate*>(cand->daughter(ndau));
//      cout <<" dau pid " << rcCand->pdgId()<<endl;
      RefToBase<reco::Track> track(rcCand->track());
      if (recotosim.find(track) != recotosim.end()) {
	
	TrackingParticleRef tpref = recotosim[track].begin()->first;
//	cout <<" matched dau pid " << tpref->pdgId()<<endl;	
	if (rcCand->pdgId() == tpref->pdgId()){
	  auto mother = getMother(tpref);
	  if (mother.isNull()){
	    continue;
	  }
//	  cout <<"mum pt = "<< mother->pt() << ", eta = "<< mother->eta() << ", pid = "<< mother->pdgId()<<endl;
	  
	  if (trueHad.isNull()){
	    trueHad = mother;
	  }
	  if (mother != trueHad){
	    continue;
	  }
	  if (abs(mother->pdgId()) == cand->pdgId()){	    
	    nmatched++;
	  }
	}
      }
    }
    nmatchedv.push_back(nmatched);
    if (nmatched){
//      cout <<"nmatched "<< nmatched<<endl;
      cout <<"trueHad "<< trueHad->pdgId() <<endl;
    }
    if (nmatched == 2){
      cout << "nmatched is 2 ===> mother tracking" << endl;
      isHadFrom(trueHad,6,count,hadFromQuark,hadFromTop);
    }
    isHadFromTsb.push_back(hadFromQuark);
    isHadFromTop.push_back(hadFromTop);
  }
  
  auto hadTruthTable = make_unique<nanoaod::FlatTable>(hadronCands->size(),"hadTruth",false);
  hadTruthTable->addColumn<int>("nMatched",nmatchedv,"no. of dau match",nanoaod::FlatTable::IntColumn);
  hadTruthTable->addColumn<int>("isHadFromTsb",isHadFromTsb,"no. of dau match",nanoaod::FlatTable::IntColumn);
  hadTruthTable->addColumn<bool>("isHadFromTop",isHadFromTop,"Hadron from top",nanoaod::FlatTable::BoolColumn);  
  iEvent.put(move(hadTruthTable),"hadTruth");
  
  auto candidates = make_unique<std::vector<reco::LeafCandidate>>();
  vector<int> imother;
  vector<int> isKsFromTsb;
  vector<bool> inVol;
  vector<bool> isKsFromTop;

  for (auto const& trackVertex : *trackingVertexs.product()) {

    if (trackVertex.eventId().bunchCrossing() != 0) continue;  // Consider only in-time events
    if (trackVertex.nDaughterTracks() < 2) continue;  // Keep only V0 vertices
    
    for (TrackingVertex::tp_iterator source = trackVertex.sourceTracks_begin(); source != trackVertex.sourceTracks_end(); ++source) {
      auto decayTrk = source->get();
      if (decayTrk->pdgId() != 310) continue; //&& decayTrk->pdgId() != 3122) continue;

      //cout << "##################################start##################################" << endl;

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
	//cout << "ev final ===> KsFromTquark : " << KsFromTquark << " KsFromQuark : " << KsFromQuark << endl;
	//cout << "ev final size : " << "candidate : " << candidates->size() << " imother : " << imother.size() << " inVol : " << inVol.size() << " isKsFromTop : " << isKsFromTop.size() << " isKsFromTsb : " << isKsFromTsb.size() << endl;


	//	//cout << " dau = " << " pdg = " << dau->pdgId() << ", pt = " << dau->p4().Pt() << endl;
      
      }
      //cout << "###################################end###################################" << endl;
    }
  }
  /*
    Handle<edm::View<reco::GenParticle> > genParticles;
    iEvent.getByToken(genLabel_,genParticles);
    int i = 0;
    for (auto gen : *genParticles){
    if (gen.pdgId() != 5122) continue;    
    //cout << " # i = " << i<< " pdg = " << gen.pdgId() << ", pt = " << gen.pt()
    << ", status = " << gen.status()
    << ", daus = " << gen.numberOfDaughters()
    << endl;
    candidates->push_back(gen);
    imother.push_back(1);
    
    for (size_t j = 0; j < gen.numberOfDaughters(); ++j){
    auto dau = gen.daughterRefVector().at(j);
    candidates->push_back(*dau);
    imother.push_back(1);
    //cout << " dau = " << i<< " pdg = " << dau->pdgId()
    << endl;
    }
    ++i;
    }
  
  */
  auto genHadTable = make_unique<nanoaod::FlatTable>(candidates->size(),"genHadron",false);
  genHadTable->addColumn<int>("mother",imother,"index of mother",nanoaod::FlatTable::IntColumn);
  genHadTable->addColumn<int>("isKsFromTsb",isKsFromTsb,"track from t->s/b",nanoaod::FlatTable::IntColumn);
  genHadTable->addColumn<bool>("isKsFromTop",isKsFromTop,"track from top",nanoaod::FlatTable::BoolColumn);
  genHadTable->addColumn<bool>("inVol",inVol,"track in volume",nanoaod::FlatTable::BoolColumn); 
  
  iEvent.put(move(genHadTable),"genHadron");
  iEvent.put(move(candidates));
}

int HadTruthProducer::trackingVertex_pdgId(const TrackingVertex* tv)
{
  cout <<"tv->nSourceTracks() "<< tv->nSourceTracks() <<endl;
  cout <<"tv->nSourceTracks() "<< tv->nSourceTracks() <<endl;
  for (TrackingVertex::tp_iterator source = tv->sourceTracks_begin(); source != tv->sourceTracks_end(); ++source) {
    return source->get()->pdgId();
  }  
  return 0;
}

const reco::GenParticleRef HadTruthProducer::getMother(const TrackingParticleRef& tp)
{
  const TrackingVertexRef& tv = tp->parentVertex();
  if (tv->nSourceTracks()){
    for (TrackingVertex::tp_iterator source = tv->sourceTracks_begin(); source != tv->sourceTracks_end(); ++source) {
      auto mothers = source->get()->genParticles();
      if (!mothers.empty()){
	reco::GenParticleRefVector::const_iterator im = mothers.begin();
	return *im;
      }
    }
  }

  if (!tp->genParticles().empty()){
    auto genpart = tp->genParticles()[0];
    const reco::GenParticleRefVector& mothers = genpart->motherRefVector();
    if (!mothers.empty()){
      reco::GenParticleRefVector::const_iterator im = mothers.begin();
      return *im;
    }
  }
  cout<<"no match to mother "<<endl;
  return reco::GenParticleRef();
}
