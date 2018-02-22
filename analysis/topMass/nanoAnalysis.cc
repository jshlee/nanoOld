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

vector<TParticle> nanoAnalysis::muonSelection()
{
  vector<TParticle> muons; 
  for (UInt_t i = 0; i < nMuon; ++i){
    if (!Muon_tightId[i]) continue;
    if (Muon_pt[i] < 20) continue;
    if (std::abs(Muon_eta[i]) > 2.4) continue;
    if (Muon_pfRelIso04_all[i] > 0.15) continue;

    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);

    //mom *= roccoR(mom, Muon_charge[i], Muon_genPartIdx[i], Muon_nTrackerLayers[i]);
    
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);

    muons.push_back(muon);
  }
  return muons;
}


vector<TParticle> nanoAnalysis::elecSelection()
{
  vector<TParticle> elecs; 
  for (UInt_t i = 0; i < nElectron; ++i){
    if (Electron_pt[i] < 20) continue;
    if (std::abs(Electron_eta[i]) > 2.4) continue;
    if (Electron_cutBased[i] < 3) continue;
    if (Electron_pfRelIso03_all[i] > 0.0571) continue;
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( std::abs(el_scEta) > 1.4442 &&  std::abs(el_scEta) < 1.566 ) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);
   
    auto elec = TParticle();
    elec.SetPdgCode(11*Electron_charge[i]*-1);
    elec.SetMomentum(mom);
    elecs.push_back(elec);
  }
  return elecs;
}


vector<TParticle> nanoAnalysis::jetSelection()
{
  vector<TParticle> jets; 

  for (UInt_t i = 0; i < nJet; ++i){
    if (Jet_pt[i] < 30) continue;
    if (std::abs(Jet_eta[i]) > 2.4) continue;
    if (Jet_jetId[i] < 1) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps){
        if (mom.TLorentzVector::DeltaR(lep) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    auto jet = TParticle();
    jet.SetMomentum(mom);
    jets.push_back(jet);
    idxs.push_back(i);
  }
  return jets;
}


vector<TParticle> nanoAnalysis::bjetSelection()
{
  vector<TParticle> bjets; 
  for (int idx : idxs ){
    if (Jet_btagCSVV2[idx] < 0.8484) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[idx], Jet_eta[idx], Jet_phi[idx], Jet_mass[idx]);
    auto bjet = TParticle();
    bjet.SetMomentum(mom);
    bjets.push_back(bjet);
  }
  return bjets;
}


void nanoAnalysis::analysis()
{
  //h_nevents->Fill(0.5);  //ori
  h_cutFlow->Fill(0);

  //Run for MC
  if (m_isMC) {
    Int_t nvtx = Pileup_nTrueInt;
    b_puweight = m_pileUp->getWeight(nvtx);
      
    b_genweight = genWeight;
    h_genweights->Fill(0.5, b_genweight);
    b_weight = b_genweight * b_puweight;
    h_weights->Fill(0.5, b_weight); //ori
  }
  else
  {
    b_puweight = 1;
    b_genweight = 0;
    if (!(m_lumi->LumiCheck(run, luminosityBlock))) return;
  }
  h_nevents->Fill(0.5,b_genweight*b_puweight); //try
    
  h_cutFlow->Fill(1);

  if (std::abs(PV_z) >= 24.) return;
  //if (PV_npvs == 0) return;
  if (PV_ndof < 4) return;

  h_cutFlow->Fill(2);

  auto muons = muonSelection();
  auto elecs = elecSelection();

  if(muons.size()+ elecs.size() != 2) return;
  h_cutFlow->Fill(3);

  int mulpdg;
  if (muons.size() == 2 && muons[0].Pt() > 25 && muons[1].Pt() > 25 ) {
      muons[0].Momentum(b_lep1);
      muons[1].Momentum(b_lep2);
      mulpdg = muons[0].GetPdgCode()*muons[1].GetPdgCode();
      b_channel = CH_MUMU;
  }
  
  if (muons.size() == 1 && elecs.size() == 1 && muons[0].Pt() > 25 && elecs[0].Pt() > 25){
      muons[0].Momentum(b_lep1);
      elecs[0].Momentum(b_lep2);
      mulpdg = muons[0].GetPdgCode()*elecs[0].GetPdgCode();
      b_channel = CH_MUEL;
  }
/*
  if (elecs.size() != 2 && elecs.size() == 1 && elecs.pt() > 25){
      elecs[0].Momentum(b_lep1);
      elecs[1].Momentum(b_lep2);
      mulpdg = elecs[0].GetPdgCode()*elecs[i].GetPdgCode();
      b_channel = CH_EL;
  }
*/
  if (elecs.size() == 2 && elecs[0].Pt() > 25){
      elecs[0].Momentum(b_lep1);
      elecs[1].Momentum(b_lep2);
      mulpdg = elecs[0].GetPdgCode()*elecs[1].GetPdgCode();
      b_channel = CH_ELEL;
  }

  

  //vector<TLorentzVector> recoleps;
  recoleps.push_back(b_lep1);
  recoleps.push_back(b_lep2);

  b_dilep = b_lep1 + b_lep2;
  
  auto jets = jetSelection();
  auto bjets = bjetSelection();
  
  if (b_dilep.M() < 20.) return;
  if (mulpdg > 0 ) return;
  b_step1 = true;
  b_step = 1;
  h_cutFlow->Fill(4);


  if (b_channel != CH_MUEL && 76 < b_dilep.M() && b_dilep.M() < 106) return;
  b_step2 = true;
  b_step = 2;
  h_cutFlow->Fill(5);
  /*
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
  */
   
  if (b_channel == CH_MUMU){
    if ((!HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL) && (!HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ)
      && (!HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL) && (!HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ)) return;
  }

  if (b_channel == CH_MUEL){
    if ((!HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL) && (!HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL)
      && (!HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && (!HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ)) return;
  }
  if (b_channel == CH_ELEL){
    if (!HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) return;
  }
  b_met = MET_pt;
  b_njet = jets.size();
  b_nbjet = bjets.size();

  if (b_channel != CH_MUEL && b_met < 40) return;
  b_step3 = true;
  b_step = 3;
  h_cutFlow->Fill(6);

  if (b_njet < 2) return;
  b_step4 = true;
  b_step = 4;
  h_cutFlow->Fill(7);

  if (b_nbjet < 1) return;
  b_step5 = true;
  b_step = 5;
  h_cutFlow->Fill(8);

}

void nanoAnalysis::LoadModules(pileUpTool* pileUp, lumiTool* lumi)
{
  //m_rocCor = rocCor;
  m_lumi = lumi;
  m_pileUp = pileUp;
}


void nanoAnalysis::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  float nPassTrig = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    std::cout << jentry << "\n";
    //Prepare for new loop
    resetBranch();
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;

    
    analysis();

    m_tree->Fill();

  }
  //for sync
  //for (int x=1; x<10; x++){
  //  cout << h_cutFlow->GetBinContent(x) << endl;
  //}
}
/*
int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " data set txt file" << std::endl;
    return 1;
  }

  TChain * chain = new TChain("Events","");
  string filename;
  ifstream infile;
  infile.open(argv[1]);

  while(!infile.eof()){
    getline(infile,filename);
    if (filename.find('#') != string::npos) continue;
    if (filename.empty()) continue;

    cout<< "Opening File: "<<filename << endl;
    chain->Add(TString(filename+"/Events"));

  }
  infile.close();


  nanoAnalysis t(chain, true);
  t.setOutput("tree.root");
  t.Loop();

  return 0;

}
*/

int main(int argc, char* argv[])
{
  if(argc != 1)
  {
    std::string env = std::getenv("CMSSW_BASE");
    std::cout << "start\n";
    //RoccoR* rocCor = new RoccoR(env+"/src/nano/analysis/data/rcdata.2016.v3/");
    //std::cout << "roccor loaded\n";
    lumiTool* lumi = new lumiTool(env+"/src/nano/analysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt");
    std::cout << "lumiTool loaded\n";
    pileUpTool* pileUp = new pileUpTool();
    std::cout << "pileup loaded\n";
    std::string dirName = env+("/src/nano/analysis/topMass/Results/")+argv[1];
    std::string temp = argv[1];
    Bool_t isMC = false;
    Size_t found = temp.find("Run");
    if(found == std::string::npos) isMC = true;
    for(Int_t i = 2; i < argc; i++)
    {
      TFile *f = TFile::Open(argv[i], "read");
      std::cout << argv[i] << "\n";                          
      TTree *tree;                  
      f->GetObject("Events", tree);
      temp = argv[i];   
      found = temp.find_last_of('/');
      std::string outPutName = dirName+temp.substr(found);
      nanoAnalysis t(tree, isMC);
      t.LoadModules(pileUp, lumi);  
      t.setOutput(outPutName);
      t.Loop();
    }
  }
  else
  {
    TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v2/DoubleMuon/Run2016B-18Apr2017_ver1-v1/171219_063744/0000/nanoAOD_23.root", "read");

    TTree *tree;
    f->GetObject("Events", tree);

    nanoAnalysis t(tree);
    t.setOutput("test.root");
    t.Loop();
  }
  return 0;
}




void nanoAnalysis::setOutput(std::string outputName)
{
  m_output = TFile::Open(outputName.c_str(), "recreate");

  m_tree = new TTree("event", "event");
  MakeBranch(m_tree);

  
  h_nevents = new TH1D("nevents", "nevents", 1, 0, 1);
  h_genweights = new TH1D("genweight", "genweight", 1, 0, 1);
  h_weights = new TH1D("weight", "weight", 1, 0, 1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
}

void nanoAnalysis::MakeBranch(TTree* t)
{
  t->Branch("nvertex", &b_nvertex, "nvertex/I");
  t->Branch("step", &b_step, "step/I");
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("njet", &b_njet, "njet/I");
  t->Branch("nbjet", &b_nbjet, "nbjet/I");
  t->Branch("step1", &b_step1, "step1/O");
  t->Branch("step2", &b_step2, "step2/O");
  t->Branch("step3", &b_step3, "step3/O");
  t->Branch("step4", &b_step4, "step4/O");
  t->Branch("step5", &b_step5, "step5/O");
  t->Branch("step6", &b_step6, "step6/O");
  t->Branch("step7", &b_step7, "step7/O");

  m_tree->Branch("lep1", "TLorentzVector", &b_lep1);
  m_tree->Branch("lep1_pid", &b_lep1_pid, "lep1_pid/I");    
  m_tree->Branch("lep2", "TLorentzVector", &b_lep2);
  m_tree->Branch("lep2_pid", &b_lep2_pid, "lep2_pid/I");    
  t->Branch("dilep", "TLorentzVector", &b_dilep);
  //m_tree->Branch("jet1", "TLorentzVector", &b_jet1);
  //m_tree->Branch("jet1_CSVInclV2", &b_jet1_CSVInclV2, "jet1_CSVInclV2/F");
  //m_tree->Branch("jet2", "TLorentzVector", &b_jet2);
  //m_tree->Branch("jet2_CSVInclV2", &b_jet2_CSVInclV2, "jet2_CSVInclV2/F");

  t->Branch("met", &b_met, "met/F");
  t->Branch("weight", &b_weight, "weight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("PV_npvs", &PV_npvs, "PV_npvs/I");
  t->Branch("ncmeson", &ncmeson, "ncmeson/I");
  t->Branch("cmeson_dca", &cmeson_dca, "cmeson_dca/F");
  t->Branch("cmeson_angleXY", &cmeson_angleXY, "cmeson_angleXY/F");
  t->Branch("cmeson_lxy", &cmeson_lxy, "cmeson_lxy/F");
  t->Branch("cmeson_l3D", &cmeson_l3D, "cmeson_l3D/F");
  t->Branch("cmeson_jetDR", &cmeson_jetDR, "cmeson_jetDR/F");
  t->Branch("cmeson_legDR", &cmeson_legDR, "cmeson_legDR/F");
}


void nanoAnalysis::resetBranch()
{
  b_lep1.SetPtEtaPhiM(0,0,0,0);
  b_lep2.SetPtEtaPhiM(0,0,0,0);
  b_dilep.SetPtEtaPhiM(0,0,0,0);
  //b_jet1.SetPtEtaPhiM(0,0,0,0);
  //b_jet2.SetPtEtaPhiM(0,0,0,0);
  recoleps.clear();
  idxs.clear();
  b_lep1_pid = 0; b_lep2_pid = 0;
  b_jet1_CSVInclV2 = -1; b_jet2_CSVInclV2 = -1;
   

  b_nvertex = 0; b_step = -1; b_channel = 0; b_njet = 0; b_nbjet = 0;
  b_step1 = 0; b_step2 = 0; b_step3 = 0; b_step4 = 0; b_step5 = 0; b_step6 = 0; b_step7 = 0;
  b_met = 0; b_weight = 1; b_genweight = 1; b_puweight = 1;
}



Double_t nanoAnalysis::roccoR(TLorentzVector m, int q, int nGen, int nTrackerLayers)
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
