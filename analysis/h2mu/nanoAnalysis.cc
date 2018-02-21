#define nanoAnalysis_cxx
#include "nanoAnalysis.h"
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

  ALL = new TTree("nEvent", "nEvent");
  MakeBranch(ALL);
  
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
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("PV_npvs", &PV_npvs, "PV_npvs/I");
  // t->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
  // t->Branch("mueffweight_up", &b_mueffweight_up, "mueffweight_up/F");
  // t->Branch("mueffweight_dn", &b_mueffweight_dn, "mueffweight_dn/F");
}

void nanoAnalysis::ResetBranch()
{
  b_Event_No = 0;
  b_Step = 0;
  b_Dilep.SetPtEtaPhiM(0,0,0,0);
  b_Mu1.SetPtEtaPhiM(0,0,0,0);
  b_Mu2.SetPtEtaPhiM(0,0,0,0);
  b_Mu_tlv.clear();
  b_El_tlv.clear();
  b_Jet_tlv.clear();
  b_bJet_tlv.clear();
  b_Event_Total = 1;
}

void nanoAnalysis::LoadModules(pileUpTool* pileUp, lumiTool* lumi, RoccoR* rocCor)
{
  //Get Modules
  m_rocCor = rocCor;
  m_lumi = lumi;
  m_pileUp = pileUp;
}

void nanoAnalysis::Analysis()
{
  h_Event_Tot->Fill(0.5, b_Event_Total);
  h_cutFlow->Fill(0);
  b_Step = 0;
  //Run for MC
  if(m_isMC)
  {
    Int_t nvtx = Pileup_nTrueInt;
    b_puweight = m_pileUp->getWeight(nvtx);
    b_genweight = genWeight;
    h_genweights->Fill(0.5, b_genweight);
    b_weight = b_genweight * b_puweight;
    h_weight->Fill(0.5, b_weight);
  } else
  {
    b_puweight = 1;
    b_genweight = 0;
    if(!(m_lumi->LumiCheck(run, luminosityBlock))) return;
  }
  b_Step = 1;
  h_cutFlow->Fill(1);

  if (abs(PV_z) >= 24.) return;
  if (PV_npvs == 0) return;
  if (PV_ndof < 4) return;
  b_Step = 2;
  h_cutFlow->Fill(2);

  auto muons = MuonSelection();
  auto electrons = ElectronSelection();
  auto jets = JetSelection();
  auto bjets = BtaggedSelection();
  
  if(muons.size() < 2) return;
  b_Step = 3;
  h_cutFlow->Fill(3);

  TParticle mu1;
  TParticle mu2;

  Bool_t charge = false;
  for(Int_t i = 0; i < muons.size(); i++)
  {
    if( (b_Mu_tlv[0].Pt() > 26) || (b_Mu_tlv[0].Pt() > 26) )
    {
      if( ( muons[0].GetPdgCode() * muons[i].GetPdgCode() ) < 0 )
      {
        b_Mu1 = b_Mu_tlv[0];
        b_Mu2 = b_Mu_tlv[i];
        charge = true;
        mu1 = muons[0];
        mu2 = muons[i];
        break;
      }
    }
  }
  if(!charge) return;
  b_Step = 4;
  h_cutFlow->Fill(4);
  b_Dilep = b_Mu1 + b_Mu2;
  if (b_Dilep.M() < 12.) return;

  // b_mueffweight = m_muonSF.getScaleFactor(mu1, 13, 0)*m_muonSF.getScaleFactor(mu2, 13, 0);
  // b_mueffweight_up = m_muonSF.getScaleFactor(mu1, 13, +1)*m_muonSF.getScaleFactor(mu2, 13, +1);
  // b_mueffweight_dn = m_muonSF.getScaleFactor(mu1, 13, -1)*m_muonSF.getScaleFactor(mu2, 13, -1);

  b_Step = 5;
  h_cutFlow->Fill(5);

  Bool_t IsoMu24 = false;
  Bool_t IsoTkMu24 = false;
  for (Int_t i = 0; i < nTrigObj; ++i){
    if (TrigObj_id[i] != 13) continue;
    if (TrigObj_pt[i] < 24) continue;
    Int_t bits = TrigObj_filterBits[i];
    if (bits & 0x2) IsoMu24 = true;
    if (bits & 0x8) IsoTkMu24 = true;	
  }
  if (!(IsoMu24 || IsoTkMu24)) return;
  b_Step = 6;
  h_cutFlow->Fill(6);

  b_Event_No = 1;
}


void nanoAnalysis::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  float nPassTrig = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //Prepare for new loop
    ResetBranch();
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    Analysis();
    ALL->Fill();
  }
}

int main(Int_t argc, Char_t** argv)
{
  string env = getenv("CMSSW_BASE");
  RoccoR* rocCor = new RoccoR(env+"/src/nano/analysis/data/rcdata.2016.v3/");
  lumiTool* lumi = new lumiTool(env+"/src/nano/analysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt");
  pileUpTool* pileUp = new pileUpTool();

  if(argc != 1)
  {
    string dirName = env+("/src/nano/analysis/h2mu/Results/")+argv[1]+"/"+argv[2];
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

//Object Selection
vector<TParticle> nanoAnalysis::MuonSelection()
{
  vector<TParticle> muons;
  for(UInt_t i = 0; i < nMuon; i++)
  {
    if (!Muon_trackerMu[i] || !Muon_globalMu[i] || !Muon_tightId[i] || Muon_pfRelIso04_all[i] > 0.25) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    mom = mom * roccoR(mom, Muon_charge[i], Muon_genPartIdx[i], Muon_nTrackerLayers[i]);

    if (mom.Pt() < 20) continue;
    else if (abs(mom.Eta()) > 2.4) continue;
    
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);

    b_Mu_tlv.push_back(mom);
    muons.push_back(muon);
  }
  return muons;
}

vector<TParticle> nanoAnalysis::ElectronSelection()
{
  vector<TParticle> electrons;
  Bool_t flag;
  for(UInt_t i = 0; i < nElectron; i++)
  {
    flag = false;
    if( Electron_pt[i] < 10 || abs(Electron_eta[i]) > 2.5 ) continue;
    if( Electron_pfRelIso03_all[i] > 0.15 || Electron_cutBased[i] < 3 ) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);

    for(UInt_t j = 0; j < b_Mu_tlv.size(); j++)
    {
      if(b_Mu_tlv[j].DeltaR(mom) < 0.3)
      {
        flag = true;
        break;
      }
    }
    if(flag) continue;

    auto elec = TParticle();
    elec.SetPdgCode(11*Electron_charge[i]*-1);
    elec.SetMomentum(mom);

    b_El_tlv.push_back(mom);
    electrons.push_back(elec);
  }
  return electrons;
}

vector<TParticle> nanoAnalysis::JetSelection()
{
  vector<TParticle> jets;
  Bool_t flag;
  for(UInt_t i = 0; i < nJet; i++)
  {
    flag = false;
    auto JET_LOOSE = (1<<0);
    if(Jet_pt[i] < 30 || abs(Jet_eta[i]) > 4.7) continue;
    if( Jet_jetId[i] & JET_LOOSE == 0 ) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);

    for(UInt_t j = 0; j < b_Mu_tlv.size(); j++)
    {
      if(b_Mu_tlv[j].DeltaR(mom) < 0.4)
      {
        flag = true;
        break;
      }
    }
    if(flag) continue;

    auto jet = TParticle();
    jet.SetMomentum(mom);

    b_Jet_tlv.push_back(mom);
    jets.push_back(jet);
  }
  return jets;
}

vector<TParticle> nanoAnalysis::BtaggedSelection()
{
  vector<TParticle> bJets;
  for (UInt_t i = 0; i < nJet; i++)
  {
    if (Jet_btagCSVV2[i] < 0.8484) continue;
    if (Jet_pt[i] < 20 || abs(Jet_eta[i]) > 2.4) continue;
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    auto bjet = TParticle();
    bjet.SetMomentum(mom);

    b_bJet_tlv.push_back(mom);
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
