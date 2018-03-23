#include "HadTruthProducer.h"
//#define debugMode

HadTruthProducer::HadTruthProducer(const edm::ParameterSet & iConfig) : recoRecoToSim_(consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("recoRecoToSim"))), recoSimToReco_(consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("recoSimToReco"))),
hadronCands_(consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("hadronCands"))),
genLabel_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genLabel"))),
trackingVertexLabel_(consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("trackingVertexLabel"))),
trackingParticleLabel_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticleLabel")))
{
  produces<nanoaod::FlatTable>("hadTruth");
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
  
  for (reco::VertexCompositeCandidateCollection::const_iterator cand = hadronCands->begin();
       cand != hadronCands->end(); cand++) {

    cout <<"pid of had " << cand->pdgId()<<endl;
    
    for (size_t ndau =0; ndau < cand->numberOfDaughters(); ++ndau){
      auto rcCand = dynamic_cast<const reco::RecoChargedCandidate*>(cand->daughter(ndau));
      RefToBase<reco::Track> track(rcCand->track());
      if (recotosim.find(track) != recotosim.end()) {
	
	TrackingParticleRef tpref = recotosim[track].begin()->first;
	// do matchin
      }
      // cout <<"pid of had " << track <<endl;
      // for (unsigned int daughter = 0; daughter < 2; ++daughter) {
      //   if (simtoreco.find(
      // 				   gen_vertex.daughterTracks()[daughter]) !=
      // 	  simtoreco.end()) {
      // 	if (!simtoreco[gen_vertex.daughterTracks()[daughter]].empty()) {
      // 	  candidateEff[daughter] = 1;  // Found a daughter track
      // 	  reco_daughter[daughter] =
      // 	    simtoreco[gen_vertex.daughterTracks()[daughter]]
      // 	    .begin()
      // 	    ->first.castTo<reco::TrackRef>();
      // 	}
      //   }
    }
  }
  
  auto candidates = make_unique<std::vector<reco::LeafCandidate>>();
  vector<int> imother;
  vector<int> isKsFromTsb;
  vector<uint8_t> inVol;
  vector<uint8_t> isKsFromTop;

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
  auto hadTable = make_unique<nanoaod::FlatTable>(candidates->size(),"hadTruth",false);
  hadTable->addColumn<int>("mother",imother,"index of mother",nanoaod::FlatTable::IntColumn);
  hadTable->addColumn<int>("isKsFromTsb",isKsFromTsb,"track from t->s/b",nanoaod::FlatTable::IntColumn);
  hadTable->addColumn<uint8_t>("isKsFromTop",isKsFromTop,"track from top",nanoaod::FlatTable::UInt8Column);
  hadTable->addColumn<uint8_t>("inVol",inVol,"track in volume",nanoaod::FlatTable::UInt8Column); 
  
  iEvent.put(move(hadTable),"hadTruth");
  iEvent.put(move(candidates));
}
