#include "HadronProducer.h"
//#define debugMode
using namespace edm;
using namespace std;

void
HadronProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder_);
  
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);
  reco::Vertex pv = recVtxs->at(0);

  Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByToken(jetLabel_, jetHandle);

  vector<hadronCandidate> hadronCandidates;
  
  math::XYZPoint primaryVertexPoint = pv.position();
  int njet = 0;
  for (const pat::Jet & aPatJet : *jetHandle){
    if (aPatJet.pt() < 30 or abs(aPatJet.eta()) > 3 ) continue;

    vector<const reco::Candidate*> chargedHadrons, leptons;
    
    for( unsigned int idx = 0 ; idx < aPatJet.numberOfDaughters() ; ++idx) {
      const reco::Candidate* dauCand ( dynamic_cast<const reco::Candidate*>(aPatJet.daughter(idx)));
#ifdef debugMode
      cout <<"dau pt = "<< dauCand->pt() << ", eta = "<< dauCand->eta() << ", pid = "<< dauCand->pdgId()<<endl;
#endif
      
      if ( dauCand->charge() == 0 ) continue;
      
      //if ( dauCand->pt() < tkPtCut_ ) continue;

      if (dauCand->bestTrack() == nullptr) continue;
      
      //const reco::Track * trk = dauCand->bestTrack();      
      // if (trk->normalizedChi2() > tkChi2Cut_) continue;
      // if (trk->numberOfValidHits() < tkNHitsCut_) continue;

      // float ipsigXY = std::abs(trk->dxy(primaryVertexPoint)/trk->dxyError());
      // //if (ipsigXY > tkIPSigXYCut_) continue;      
      // float ipsigZ = std::abs(trk->dz(primaryVertexPoint)/trk->dzError());
      //if (ipsigZ > tkIPSigZCut_) continue;
      
      if ( abs(dauCand->pdgId()) == 11  || abs(dauCand->pdgId())==13)
	leptons.emplace_back(dauCand);
      else
	chargedHadrons.emplace_back( dauCand );

    }
    unsigned int dau_size = chargedHadrons.size() + leptons.size();
    if ( dau_size < 2 ) continue;

    if (leptons.size() > 1){
      auto jpsiCands = findJPsiCands(leptons, pv, njet, aPatJet);      
      hadronCandidates.insert(hadronCandidates.end(), jpsiCands.begin(), jpsiCands.end());
    }
	
    ++njet;
  }


  // saving all variables
  auto had_cands = make_unique<reco::VertexCompositeCandidateCollection>();
  vector<int> had_nJet, had_nDau;
  vector<float> had_jetDR, had_legDR, had_diffMass;
  vector<float> had_lxy, had_lxySig, had_l3D, had_l3DSig, had_dca, had_angleXY, had_angleXYZ;
  vector<float> had_dau1_chi2, had_dau1_nHits, had_dau1_pt, had_dau1_ipsigZ, had_dau1_ipsigXY;
  vector<float> had_dau2_chi2, had_dau2_nHits, had_dau2_pt, had_dau2_ipsigZ, had_dau2_ipsigXY;
  
  for (auto cand: hadronCandidates){
    had_cands->push_back(cand.vcc);
    
    const reco::Track* dau1 = cand.vcc.daughter(0)->bestTrack();
    had_dau1_chi2.push_back(dau1->normalizedChi2());
    had_dau1_nHits.push_back(dau1->numberOfValidHits());
    had_dau1_pt.push_back(dau1->pt());
    had_dau1_ipsigZ.push_back(std::abs(dau1->dz(primaryVertexPoint)/dau1->dzError()));
    had_dau1_ipsigXY.push_back(std::abs(dau1->dxy(primaryVertexPoint)/dau1->dxyError()));

    const reco::Track* dau2 = cand.vcc.daughter(1)->bestTrack();
    had_dau2_chi2.push_back(dau2->normalizedChi2());
    had_dau2_nHits.push_back(dau2->numberOfValidHits());
    had_dau2_pt.push_back(dau2->pt());
    had_dau2_ipsigZ.push_back(std::abs(dau2->dz(primaryVertexPoint)/dau2->dzError()));
    had_dau2_ipsigXY.push_back(std::abs(dau2->dxy(primaryVertexPoint)/dau2->dxyError()));
    
    had_nJet.push_back(cand.nJet);
    had_nDau.push_back(cand.nDau);
    had_jetDR.push_back(cand.jetDR);
    had_legDR.push_back(cand.legDR);
    had_diffMass.push_back(cand.diffMass);
    had_lxy.push_back(cand.lxy);
    had_lxySig.push_back(cand.lxySig);
    had_l3D.push_back(cand.l3D);
    had_l3DSig.push_back(cand.l3DSig);
    had_dca.push_back(cand.dca);
    had_angleXY.push_back(cand.angleXY);
    had_angleXYZ.push_back(cand.angleXYZ);
  }

  
  auto had_table = make_unique<nanoaod::FlatTable>(had_cands->size(),"had",false);
  had_table->addColumn<int>("nJet",had_nJet,"nJet of vertex cand",nanoaod::FlatTable::IntColumn);
  had_table->addColumn<int>("nDau",had_nDau,"nDau of vertex cand",nanoaod::FlatTable::IntColumn);

  had_table->addColumn<float>("jetDR",had_jetDR,"DR between jet",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("legDR",had_legDR,"DR between daugthers",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("diffMass",had_diffMass,"mass difference",nanoaod::FlatTable::FloatColumn);

  had_table->addColumn<float>("lxy",had_lxy,"2D decay length in cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("lxySig",had_lxySig,"2D decay length sig in cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("l3D",had_l3D,"3D decay length in cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("l3DSig",had_l3DSig,"3D decay length sig in cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dca",had_dca,"distance of closest approach cm",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("angleXY",had_angleXY,"2D angle between vertex and tracks",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("angleXYZ",had_angleXYZ,"3D angle between vertex and tracks",nanoaod::FlatTable::FloatColumn);
  
  had_table->addColumn<float>("dau1_chi2",had_dau1_chi2,"dau1 chi2/ndof",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau1_nHits",had_dau1_nHits,"dau1 nHits",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau1_pt",had_dau1_pt,"dau1 Pt",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau1_ipsigXY",had_dau1_ipsigXY,"dau1 ipsigXY",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau1_ipsigZ",had_dau1_ipsigZ,"dau1 ipsigZ",nanoaod::FlatTable::FloatColumn);

  had_table->addColumn<float>("dau2_chi2",had_dau2_chi2,"dau2 chi2/ndof",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau2_nHits",had_dau2_nHits,"dau2 nHits",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau2_pt",had_dau2_pt,"dau2 Pt",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau2_ipsigXY",had_dau2_ipsigXY,"dau2 ipsigXY",nanoaod::FlatTable::FloatColumn);
  had_table->addColumn<float>("dau2_ipsigZ",had_dau2_ipsigZ,"dau2 ipsigZ",nanoaod::FlatTable::FloatColumn);

  iEvent.put(move(had_table),"had");
  iEvent.put(move(had_cands));
  
}

vector<HadronProducer::hadronCandidate> HadronProducer::findJPsiCands(vector<const reco::Candidate*> &leptons, reco::Vertex& pv, int nJet, const pat::Jet & aPatJet)
{
  vector<hadronCandidate> hadrons;
  // find jpsi to mumu or ee
  for (auto lep1Cand : leptons){
    //if (lep1Cand->pdgId() > 0) continue; 
    for (auto lep2Cand : leptons){
      // if (lep2Cand->pdgId() < 0) continue; 

      // int pdgMul = lep1Cand->pdgId() * lep2Cand->pdgId();
      // if ( pdgMul != -121 and pdgMul != -169 ) continue; 

      vector<const reco::Candidate*> cands{lep1Cand, lep2Cand};

      hadronCandidate hc;

      reco::VertexCompositeCandidate cand = fit(cands, pv, jpsi_pdgId_,
						hc.dca, hc.angleXY, hc.angleXYZ);

      if (cand.numberOfDaughters() < 2) continue;
      // if (abs(cand.mass() - jpsi_m_) > 0.5) continue;
      
      hc.vcc = cand;
      
      auto d2 = getDistance(2,cand,pv);
      hc.lxy = d2.first;      
      hc.lxySig = d2.second;
      auto d3 = getDistance(3,cand,pv);	
      hc.l3D = d3.first;
      hc.l3DSig = d3.second;

      hc.nJet = nJet;
      hc.nDau = 2;
      hc.diffMass = -9;

      hc.jetDR = reco::deltaR( cand, aPatJet);
      hc.legDR = reco::deltaR( *lep1Cand, *lep2Cand);
      
      hadrons.emplace_back(hc);
    }
  }
  return hadrons;
}
