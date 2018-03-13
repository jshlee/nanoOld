#define nanoAnalysis_cxx
#include "nano/analysis/src/nanoAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cstdlib>

using namespace std;
/*
To compile:
g++ `root-config --cflags --glibs` -lEG nanoAnalysis.cc -o nanoAnalysis
*/

void nanoAnalysis::SetOutput(string outputName)
{
  m_output = TFile::Open(outputName.c_str(), "recreate");

  m_tree = new TTree("events", "events");
  MakeBranch(m_tree);
  
  h_Event_Tot = new TH1D("Event_total", "Event_total" ,1,0,1);
  h_genweights = new TH1D("genweight", "genweight" , 1,0,1);
  h_weight = new TH1D("weight", "weight", 1,0,1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
}

void nanoAnalysis::MakeBranch(TTree* t)
{
  t->Branch("Event_No", &b_Event_No, "Event_No/I");
  t->Branch("Step", &b_Step, "Step/I");
  t->Branch("Dilep", "TLorentzVector", &b_Dilep);
  t->Branch("Mu1", "TLorentzVector", &b_Mu1);
  t->Branch("Mu2", "TLorentzVector", &b_Mu2);
  t->Branch("lep1", "TLorentzVector", &b_lep1);
  t->Branch("lep2", "TLorentzVector", &b_lep2);
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("weight", &b_weight, "weight/F");
  t->Branch("PV_npvs", &PV_npvs, "PV_npvs/I");
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("charge", &b_charge, "charge/I");
  t->Branch("nlep", &b_nlep, "nlep/I");
  t->Branch("nmuon", &b_nmuon, "nmuon/I");
  t->Branch("nelec", &b_nelec, "nelec/I");
  t->Branch("njet", &b_njet, "njet/I");
  t->Branch("nbjet", &b_nbjet, "nbjet/I");
  t->Branch("trig_m", &b_trig_m, "trig_m/O");
  t->Branch("trig_m2", &b_trig_m2, "trig_m2/O");
  t->Branch("trig_e", &b_trig_e, "trig_e/O");
  t->Branch("trig_mm", &b_trig_mm, "trig_mm/O");
  t->Branch("trig_em", &b_trig_em, "trig_em/O");
  t->Branch("trig_ee", &b_trig_ee, "trig_ee/O");
  t->Branch("Met", &b_Met, "Met/F");
  t->Branch("Met_phi", &b_Met_phi, "Met_phi/F");
  t->Branch("CSVv2", &b_CSVv2);
  t->Branch("FL", &b_FL, "FL/I");
  t->Branch("FH", &b_FH, "FH/I");
  t->Branch("SL", &b_SL, "SL/I");
  t->Branch("csvweight", "std::vector<float>", &b_csvweights);
  t->Branch("btagweight", &b_btagweight, "btagweight/F");
  t->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
  t->Branch("mueffweight_up", &b_mueffweight_up, "mueffweight_up/F");
  t->Branch("mueffweight_dn", &b_mueffweight_dn, "mueffweight_dn/F");
}

void nanoAnalysis::ResetBranch()
{
  b_Event_No = 0;
  b_Step = 0;
  b_Dilep.SetPtEtaPhiM(0,0,0,0);
  b_Mu1.SetPtEtaPhiM(0,0,0,0);
  b_Mu2.SetPtEtaPhiM(0,0,0,0);
  b_lep1.SetPtEtaPhiM(0,0,0,0);
  b_lep2.SetPtEtaPhiM(0,0,0,0);
  b_Mu_tlv.clear();
  b_El_tlv.clear();
  b_Jet_tlv.clear();
  b_bJet_tlv.clear();
  b_CSVv2.clear();
  b_Event_Total = 1;
  b_channel = -1;
  b_nlep = 0; b_nmuon = 0; b_nelec = 0;
  b_charge = 0;
  b_Met_phi = 0; b_Met = 0;
  b_FL = 0; b_SL = 0; b_FH = 0;
}

void nanoAnalysis::LoadModules(pileUpTool* pileUp, lumiTool* lumi, RoccoR* rocCor)
{
  //Get Modules
  m_rocCor = rocCor;
  m_lumi = lumi;
  m_pileUp = pileUp;
}

bool nanoAnalysis::Analysis()
{
  h_Event_Tot->Fill(0.5, b_Event_Total);
  h_cutFlow->Fill(0);
  //Run for MC
  if(m_isMC){
    Int_t nvtx = Pileup_nTrueInt;
    b_puweight = m_pileUp->getWeight(nvtx);
    b_genweight = genWeight;
    h_genweights->Fill(0.5, b_genweight);
    b_weight = b_genweight * b_puweight;
    h_weight->Fill(0.5, b_weight);
  }
  else {
    b_puweight = 1;
    b_genweight = 0;
    if(!(m_lumi->LumiCheck(run, luminosityBlock))) return false;
  }
  h_cutFlow->Fill(1);

  if (abs(PV_z) >= 24.) return false;
  if (PV_npvs == 0) return false;
  if (PV_ndof < 4) return false;
  h_cutFlow->Fill(2);
  
  auto looseMuons = LooseMuonSelection();
  auto looseElecs = LooseElectronSelection(looseMuons);
  auto looseJets = LooseJetSelection(looseMuons);
  auto looseBJets = LooseBJetSelection(looseMuons);
  
  auto muons = MuonSelection();
  auto elecs = ElectronSelection();
  auto jets = JetSelection(looseMuons);
  auto bjets = BJetSelection(looseMuons);
  
  if (looseMuons.size() + looseElecs.size() < 2) return false;

  b_trig_m = HLT_IsoTkMu24 || HLT_IsoMu24;
  b_trig_e = HLT_Ele27_WPTight_Gsf;  
  b_trig_mm = HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL
    || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  b_trig_em = HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL
    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
  b_trig_ee = HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  
  // make all the variables that you need to save here
  TParticle mu1;
  TParticle mu2; 
  b_Met = PuppiMET_pt;
  b_Met_phi = PuppiMET_phi;
  Bool_t IsoMu24 = false;
  Bool_t IsoTkMu24 = false;
  b_nlep = looseMuons.size() + looseElecs.size();
  b_nmuon = looseMuons.size();
  b_nelec = looseElecs.size();
  b_njet = looseJets.size();
  b_nbjet = looseBJets.size();
  

  // do ttH analysis steps here
  // ttH uses loose leptons
//////////////////////////////////////////// TTH -> MuMu Start //////////////////////////////////////////////////////
  for (UInt_t i = 0; i < nTrigObj; ++i){
    if (TrigObj_id[i] != 13) continue;
    if (TrigObj_pt[i] < 24) continue;
    Int_t bits = TrigObj_filterBits[i];
    if (bits & 0x2) IsoMu24 = true;
    if (bits & 0x8) IsoTkMu24 = true;  
  }
  if ((IsoMu24 || IsoTkMu24))
  {
    b_trig_m2 = true;
  }
  
  for(UInt_t i = 0; i < looseMuons.size(); i++)
  { 
    if( (b_Mu_tlv[0].Pt() > 26) || (b_Mu_tlv[0].Pt() > 26) )
    { 
      if( ( looseMuons[0].GetPdgCode() * looseMuons[i].GetPdgCode() ) < 0 )
      { 
        b_Mu1 = b_Mu_tlv[0];
        b_Mu2 = b_Mu_tlv[i];
        b_charge = 1;
        mu1 = looseMuons[0];
        mu2 = looseMuons[i];
        break;
      }
    }
  }
  b_Dilep = b_Mu1 + b_Mu2;
 
  /////////////////////// Fully Leptonic //////////////////
  if (looseMuons.size() + looseElecs.size() == 4)
  {

    Int_t mulpdg = -1;

    if (looseMuons.size() == 2)
    {
      looseElecs[0].Momentum(b_lep1);
      looseElecs[1].Momentum(b_lep2);
      mulpdg = looseElecs[0].GetPdgCode()*looseElecs[1].GetPdgCode();
      b_channel = CH_ELEL;
    }

    else if (looseMuons.size() == 3)
    {
      for (UInt_t i = 1; i < looseMuons.size(); i++)
      {
        if (looseMuons[i].Pt() != b_Mu2.Pt())
        {
          mulpdg = looseMuons[i].GetPdgCode()*looseElecs[0].GetPdgCode();
          looseMuons[i].Momentum(b_lep1);
          break;
        }
      }
      looseElecs[0].Momentum(b_lep2);
      b_channel = CH_MUEL;
    }

    else if(looseMuons.size() == 4)
    {
      for (UInt_t i = 1; i < looseMuons.size(); i++)
      {
        if (looseMuons[i].Pt() != b_Mu2.Pt())
        {
          mulpdg = looseMuons[i].GetPdgCode()*looseMuons[3].GetPdgCode();
          looseMuons[i].Momentum(b_lep1);
         break;
        }
      }
      looseMuons[3].Momentum(b_lep2);
      b_channel = CH_MUMU;
    }

    if (mulpdg < 0)
    {
      b_FL = 1;
    }
  }
  ////////////////////// Fully Hadronic //////////////////
  if (looseElecs.size() + looseMuons.size() == 0 && looseJets.size() >= 4)
  {
     if ((looseBJets.size() == 1)||(looseBJets.size() == 2))
     {
       b_FH = 4;
     }
  }
  if (looseElecs.size() + looseMuons.size() == 0 && looseJets.size() >= 3)
  {
     if ((looseBJets.size() == 1)||(looseBJets.size() == 2))
     {
       b_FH = 3;
     }
  }
  
  if (looseElecs.size() + looseMuons.size() == 0 && looseJets.size() >= 2)
  {
     if ((looseBJets.size() == 1)||(looseBJets.size() == 2))
     {
       b_FH = 2;
     }
  }

  /////////////////////// Semi Leptonic //////////////////
  if (looseMuons.size() + looseElecs.size() == 3)
  {
    if (looseBJets.size() == 1)
    {
      if (looseJets.size() >= 2)
      {
        b_SL = 1;
      }
    }
    else if (looseBJets.size() == 2)
    {
      if (looseJets.size() >=1)
      {
        b_SL = 1;
      }
    }
  }
  
  b_mueffweight = m_muonSF.getScaleFactor(mu1, 13, 0)*m_muonSF.getScaleFactor(mu2, 13, 0);
  b_mueffweight_up = m_muonSF.getScaleFactor(mu1, 13, +1)*m_muonSF.getScaleFactor(mu2, 13, +1);
  b_mueffweight_dn = m_muonSF.getScaleFactor(mu1, 13, -1)*m_muonSF.getScaleFactor(mu2, 13, -1);
  b_Event_No = 1;

//////////////////////////////////////////// TTH -> MuMu END //////////////////////////////////////////////////////

  // do top analysis steps here

  
  return true;
}


void nanoAnalysis::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //Prepare for new loop
    ResetBranch();
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    bool keep = Analysis();
    if (keep){
      m_tree->Fill();
    }
  }
  
}

int main(Int_t argc, Char_t** argv)
{
  string env = getenv("CMSSW_BASE");
  string username = getenv("USER");
  RoccoR* rocCor = new RoccoR(env+"/src/nano/analysis/data/rcdata.2016.v3/");
  lumiTool* lumi = new lumiTool(env+"/src/nano/analysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt");
  pileUpTool* pileUp = new pileUpTool();

  if(argc != 1)
  {
    //string dirName = env+("/src/nano/analysis/h2mu/Results/")+argv[1]+"/"+argv[2];
    string dirName = "root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/"+username+"/nanoAOD/"+std::string(argv[1])+"/"+std::string(argv[2]);
    string temp = argv[2];
    Bool_t isMC = false;
    Size_t found = temp.find("Run");
    if(found == string::npos) isMC = true;
    for(Int_t i = 3; i < argc; i++)
    {
      TFile *f = TFile::Open(argv[i], "read");

      TTree *tree;
      f->GetObject("Events", tree);

      temp = argv[i];
      found = temp.find_last_of('/');
      string outPutName = dirName+temp.substr(found);
      nanoAnalysis t(tree, isMC);
      t.LoadModules(pileUp, lumi, rocCor);
      t.SetOutput(outPutName);
      t.Loop();
    }
  }
  else
  {
    TFile *f = TFile::Open("root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/nanoAOD/run2_2016v3/ttHToMuMu_M125_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6/180125_131219/0000/nanoAOD_982.root", "read");
    
    TTree *tree;
    f->GetObject("Events", tree);
    
    nanoAnalysis t(tree, true);
    t.LoadModules(pileUp, lumi, rocCor);
    t.SetOutput("test.root");
    t.Loop();
  }

  return 0;
}

//Object Selections
vector<TParticle> nanoAnalysis::LooseMuonSelection()
{
  vector<TParticle> muons;
  for(UInt_t i = 0; i < nMuon; i++)
  {
    if (!Muon_trackerMu[i]) continue;
    if (!Muon_globalMu[i]) continue;
    if (!Muon_tightId[i]) continue;
    if (Muon_pfRelIso04_all[i] > 0.25) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    mom = mom * roccoR(mom, Muon_charge[i], Muon_genPartIdx[i], Muon_nTrackerLayers[i]);
    
    if (mom.Pt() < 10) continue;
    if (abs(mom.Eta() > 2.4)) continue;
   
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);

    b_Mu_tlv.push_back(mom);
    muons.push_back(muon);
  }
  return muons;
}

vector<TParticle> nanoAnalysis::LooseElectronSelection(vector<TParticle> leptons)
{
  vector<TParticle> electrons;
  for(UInt_t i = 0; i < nElectron; i++)
  {
    if ( Electron_pt[i] < 10) continue;
    if (abs(Electron_eta[i]) > 2.5 ) continue; //<~~~~~~~~~~~~~~ Higgs Electron pt == 10; Higgs Electron eta > 2.5  
    //if( Electron_pfRelIso03_all[i] > 0.15 || Electron_cutBased[i] < 3 ) continue; //<~~~~~~~~~~~~~~ Top doesn't use isolation? 
    if (Electron_cutBased[i] < 3) continue;
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( std::abs(el_scEta) > 1.4442 &&  std::abs(el_scEta) < 1.566 ) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);

    if (hasOverLap(mom, leptons)) continue;
    
    auto elec = TParticle();
    elec.SetPdgCode(11*Electron_charge[i]*-1);
    elec.SetMomentum(mom);

    b_El_tlv.push_back(mom);
    electrons.push_back(elec);
  }
  return electrons;
}
vector<TParticle> nanoAnalysis::LooseJetSelection(vector<TParticle> leptons)
{
  vector<TParticle> jets;
  float Jet_SF_CSV[19] = {1.0,};
  for(UInt_t i = 0; i < nJet; i++)
  {
    if (Jet_pt[i] < 30) continue;
    if (abs(Jet_eta[i]) > 4.7) continue; 

    if (Jet_jetId[i] < 1) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    
    if (hasOverLap(mom, leptons)) continue;
    
    auto jet = TParticle();
    jet.SetMomentum(mom);

    jets.push_back(jet);
    for (UInt_t iu = 0; iu < 19; iu++)
    {
      Jet_SF_CSV[iu] *= m_btagSF.getSF(jet, Jet_btagCSVV2[i], Jet_hadronFlavour[i], iu);
    }
  }
  for (UInt_t i =0; i<19; i++) b_csvweights.push_back(Jet_SF_CSV[i]);
  b_btagweight = Jet_SF_CSV[0];
  return jets;
}

vector<TParticle> nanoAnalysis::LooseBJetSelection(vector<TParticle> leptons)
{
  vector<TParticle> bJets;
  for (UInt_t i = 0; i < nJet; i++)
  {
    if (Jet_btagCSVV2[i] < 0.8484) continue;
    if (Jet_pt[i] < 20) continue;
    if (abs(Jet_eta[i]) > 2.4) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);

    if (hasOverLap(mom, leptons)) continue;
    
    auto bjet = TParticle();
    bjet.SetMomentum(mom);
    bJets.push_back(bjet);
    b_CSVv2.push_back(Jet_btagCSVV2[i]);
  }
  return bJets;
}

// for top selection
vector<TParticle> nanoAnalysis::MuonSelection()
{
  vector<TParticle> muons;
  for(UInt_t i = 0; i < nMuon; i++)
  {
    if (Muon_pt[i] < 10) continue;
    if (abs(Muon_eta[i]) > 2.4) continue;
    if (!Muon_trackerMu[i]) continue;
    if (!Muon_globalMu[i]) continue;
    if (!Muon_tightId[i]) continue;
    if (Muon_pfRelIso04_all[i] > 0.25) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    mom = mom * roccoR(mom, Muon_charge[i], Muon_genPartIdx[i], Muon_nTrackerLayers[i]);
   
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);
    muons.push_back(muon);
  }
  return muons;
}

vector<TParticle> nanoAnalysis::ElectronSelection()
{
  vector<TParticle> electrons;
  for(UInt_t i = 0; i < nElectron; i++)
  {
    if ( Electron_pt[i] < 10) continue;
    if (abs(Electron_eta[i]) > 2.5 ) continue; //<~~~~~~~~~~~~~~ Higgs Electron pt == 10; Higgs Electron eta > 2.5  
    //if( Electron_pfRelIso03_all[i] > 0.15 || Electron_cutBased[i] < 3 ) continue; //<~~~~~~~~~~~~~~ Top doesn't use isolation? 
    if (Electron_cutBased[i] < 3) continue;
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( std::abs(el_scEta) > 1.4442 &&  std::abs(el_scEta) < 1.566 ) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);

    auto elec = TParticle();
    elec.SetPdgCode(11*Electron_charge[i]*-1);
    elec.SetMomentum(mom);

    electrons.push_back(elec);
  }
  return electrons;
}

bool nanoAnalysis::hasOverLap(TLorentzVector cand, vector<TParticle> objects)
{
  for (auto obj: objects){
    TLorentzVector mom;
    obj.Momentum(mom);
    if (cand.DeltaR(mom) < 0.4){
      return true;
    }
  }
  return false;
}

vector<TParticle> nanoAnalysis::JetSelection(vector<TParticle> leptons)
{
  vector<TParticle> jets;
  float Jet_SF_CSV[19] = {1.0,};
  for(UInt_t i = 0; i < nJet; i++)
  {
    if (Jet_pt[i] < 30) continue;
    if (abs(Jet_eta[i]) > 4.7) continue; 

    if (Jet_jetId[i] < 1) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    
    if (hasOverLap(mom, leptons)) continue;
    
    auto jet = TParticle();
    jet.SetMomentum(mom);

    jets.push_back(jet);
    for (UInt_t iu = 0; iu < 19; iu++)
    {
      Jet_SF_CSV[iu] *= m_btagSF.getSF(jet, Jet_btagCSVV2[i], Jet_hadronFlavour[i], iu);
    }
  }
  for (UInt_t i =0; i<19; i++) b_csvweights.push_back(Jet_SF_CSV[i]);
  b_btagweight = Jet_SF_CSV[0];
  return jets;
}

vector<TParticle> nanoAnalysis::BJetSelection(vector<TParticle> leptons)
{
  vector<TParticle> bJets;
  for (UInt_t i = 0; i < nJet; i++)
  {
    if (Jet_btagCSVV2[i] < 0.8484) continue;
    if (Jet_pt[i] < 20) continue;
    if (abs(Jet_eta[i]) > 2.4) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);

    if (hasOverLap(mom, leptons)) continue;
    
    auto bjet = TParticle();
    bjet.SetMomentum(mom);
    bJets.push_back(bjet);
  }
  return bJets;
}

Double_t nanoAnalysis::roccoR(TLorentzVector m, int &q, int &nGen, int &nTrackerLayers)
{
  Float_t u1 = gRandom->Rndm();
  Float_t u2 = gRandom->Rndm();
  if (!m_isMC){
    return m_rocCor->kScaleDT(q, m.Pt(), m.Eta(), m.Phi(), 0, 0);
  }
  else {
    if (nGen > -1){
      return m_rocCor->kScaleFromGenMC(q, m.Pt(), m.Eta(), m.Phi(),
				       nTrackerLayers, GenPart_pt[nGen],
				       u1, 0, 0);
    }
    else
      return m_rocCor->kScaleAndSmearMC(q, m.Pt(), m.Eta(), m.Phi(),
					nTrackerLayers, u1, u2, 0, 0);
  }
  return 1.0;
}