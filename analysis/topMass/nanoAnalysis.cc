#define nanoAnalysis_cxx
#include "nanoAnalysis.h"
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
    if (!Muon_trackerMu[i] || !Muon_globalMu[i] || !Muon_tightId[i]) continue;
    if (Muon_pt[i] < 20) continue;
    if (std::abs(Muon_eta[i]) > 2.4) continue;
    if (Muon_pfRelIso04_all[i] > 0.25) continue;

    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);

    mom *= roccoR(mom, Muon_charge[i], Muon_genPartIdx[i], Muon_nTrackerLayers[i]);
    
    auto muon = TParticle();
    muon.SetPdgCode(13*Muon_charge[i]*-1);
    muon.SetMomentum(mom);

    // int genN = Muon_genPartIdx[i];
    // cout << " genN "<< genN<< endl;
    // cout << " muon pid "<< muon.GetPdgCode()
    // 	 << " pt "<< muon.Pt()
    // 	 << " eta "<< muon.Eta()
    // 	 << " gen pid "<< GenPart_pdgId[genN]
    // 	 << " pt "<< GenPart_pt[genN]
    // 	 << " eta "<< GenPart_eta[genN]
    // 	 << endl;
    muons.push_back(muon);
  }
  return muons;
}

void nanoAnalysis::analysis()
{
  h_nevents->Fill(0.5);
  h_cutFlow->Fill(0);

  //Run for MC
  if (m_isMC) {
    Int_t nvtx = Pileup_nTrueInt;
    if (nvtx < hist_mc->GetNbinsX()) b_puweight = m_puWeightCalculator->getWeight(nvtx);
    else b_puweight = 1;
      
    b_genweight = genWeight;
    h_genweights->Fill(0.5, b_genweight);
    b_weight = b_genweight * b_puweight;
    h_weights->Fill(0.5, b_weight);
    
    mcAnalysis();
  }
  else {
    b_puweight = 1;
    b_genweight = 0;
    if (!lumiCheck()) return;
  }
    
  b_step = 1;
  h_cutFlow->Fill(1);

  if (std::abs(PV_z) >= 24.) return;
  if (PV_npvs == 0) return;
  if (PV_ndof < 4) return;

  b_step = 2;
  h_cutFlow->Fill(2);

  auto muons = muonSelection();
    
  if(muons.size() < 2) return;
  b_step = 3;
  h_cutFlow->Fill(3);

  if (muons[0].GetPdgCode()*muons[1].GetPdgCode() > 0 ) return;
  b_step = 4;
  h_cutFlow->Fill(4);

  muons[0].Momentum(b_lep1);
  muons[1].Momentum(b_lep2);
  
  b_dilep = b_lep1 + b_lep2;

  if (b_dilep.M() < 12.) return;
  b_step = 5;
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
  b_step = 6;
  h_cutFlow->Fill(6);
}

void nanoAnalysis::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  float nPassTrig = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //Prepare for new loop
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;

    resetBranch();
    
    analysis();

    m_tree->Fill();

  }
}

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

void nanoAnalysis::setOutput(std::string outputName)
{
  m_output = TFile::Open(outputName.c_str(), "recreate");

  m_tree = new TTree("event", "event");
  m_tree->Branch("nvertex", &b_nvertex, "nvertex/I");
  m_tree->Branch("step", &b_step, "step/I");
  m_tree->Branch("channel", &b_channel, "channel/I");
  m_tree->Branch("njet", &b_njet, "njet/I");
  m_tree->Branch("nbjet", &b_nbjet, "nbjet/I");
  m_tree->Branch("step1", &b_step1, "step1/O");
  m_tree->Branch("step2", &b_step2, "step2/O");
  m_tree->Branch("step3", &b_step3, "step3/O");
  m_tree->Branch("step4", &b_step4, "step4/O");
  m_tree->Branch("step5", &b_step5, "step5/O");
  m_tree->Branch("step6", &b_step6, "step6/O");
  m_tree->Branch("step7", &b_step7, "step7/O");

  m_tree->Branch("lep1", "TLorentzVector", &b_lep1);
  m_tree->Branch("lep1_pid", &b_lep1_pid, "lep1_pid/I");    
  m_tree->Branch("lep2", "TLorentzVector", &b_lep2);
  m_tree->Branch("lep2_pid", &b_lep2_pid, "lep2_pid/I");    
  m_tree->Branch("dilep", "TLorentzVector", &b_dilep);
  m_tree->Branch("jet1", "TLorentzVector", &b_jet1);
  m_tree->Branch("jet1_CSVInclV2", &b_jet1_CSVInclV2, "jet1_CSVInclV2/F");
  m_tree->Branch("jet2", "TLorentzVector", &b_jet2);
  m_tree->Branch("jet2_CSVInclV2", &b_jet2_CSVInclV2, "jet2_CSVInclV2/F");

  m_tree->Branch("met", &b_met, "met/F");
  m_tree->Branch("weight", &b_weight, "weight/F");
  m_tree->Branch("puweight", &b_puweight, "puweight/F");
  m_tree->Branch("genweight", &b_genweight, "genweight/F");
  
  h_nevents = new TH1D("nevents", "nevents", 1, 0, 1);
  h_genweights = new TH1D("genweight", "genweight", 1, 0, 1);
  h_weights = new TH1D("weight", "weight", 1, 0, 1);
  h_cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
}

void nanoAnalysis::resetBranch()
{
  b_lep1.SetPtEtaPhiM(0,0,0,0);
  b_lep2.SetPtEtaPhiM(0,0,0,0);
  b_dilep.SetPtEtaPhiM(0,0,0,0);
  b_jet1.SetPtEtaPhiM(0,0,0,0);
  b_jet2.SetPtEtaPhiM(0,0,0,0);
  
  b_lep1_pid = -1; b_lep2_pid = -1;
  b_jet1_CSVInclV2 = -1; b_jet2_CSVInclV2 = -1;

  b_nvertex = -1; b_step = -1; b_channel = -1; b_njet = -1; b_nbjet = -1;
  b_step1 = 0; b_step2 = 0; b_step3 = 0; b_step4 = 0; b_step5 = 0; b_step6 = 0; b_step7 = 0;
  b_met = -1; b_weight = -1; b_genweight = -1; b_puweight = -1;
}

Bool_t nanoAnalysis::lumiCheck()
{
  if ( lumiMap.find(run) == lumiMap.end() ) {
    return false;
  } else {
    for (UInt_t i = 0; i < lumiMap[run].size(); i++){
      if(lumiMap[run][i][0] <= luminosityBlock && lumiMap[run][i][1] >= luminosityBlock) return true;
    }
    return false;
  }
}

void nanoAnalysis::mcAnalysis()
{
 
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
