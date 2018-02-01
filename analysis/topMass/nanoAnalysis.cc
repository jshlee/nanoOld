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
  ResetBranch();
  Event_Tot->Fill(0.5, Event_Total);
  cutFlow->Fill(0);

  //Run for MC
  if (isMC) {
    Int_t nvtx = Pileup_nTrueInt;
    if (nvtx < hist_mc->GetNbinsX()) puweight = puWeightCalculator->getWeight(nvtx);
    else puweight = 1;
      
    genweight = genWeight;
    genweights->Fill(0.5, genweight);
    b_weight = genweight * puweight;
    weight->Fill(0.5, b_weight);
      
    fillMcBranch();
  }
  else {
    puweight = 1;
    genweight = 0;
    if (!LumiCheck()) return;
  }
    
  Step = 1;
  cutFlow->Fill(1);

  if (std::abs(PV_z) >= 24.) return;
  if (PV_npvs == 0) return;
  if (PV_ndof < 4) return;

  Step = 2;
  cutFlow->Fill(2);

  auto muons = muonSelection();
    
  if(muons.size() < 2) return;
  Step = 3;
  cutFlow->Fill(3);

  if (muons[0].GetPdgCode()*muons[1].GetPdgCode() > 0 ) return;
  Step = 4;
  cutFlow->Fill(4);

  muons[0].Momentum(Mu1);
  muons[1].Momentum(Mu2);
    
  Dilep = Mu1 + Mu2;

  if (Dilep.M() < 12.) return;
  Step = 5;
  cutFlow->Fill(5);

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
  Step = 6;
  cutFlow->Fill(6);

  Event_No = 1;
  ALL->Fill();
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
    analysis();
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
  t.SetOutput("tree.root");
  t.Loop();

  return 0;
}

void nanoAnalysis::SetOutput(std::string outputName)
{
  output = TFile::Open(outputName.c_str(), "recreate");

  ALL = new TTree("nEvent", "nEvent");
  MakeBranch(ALL);
  
  Event_Tot = new TH1D("Event_total", "Event_total" ,1,0,1);
  genweights = new TH1D("genweight", "genweight" , 1,0,1);
  weight = new TH1D("weight", "weight", 1,0,1);
  cutFlow = new TH1D("cutflow", "cutflow", 11, -0.5, 10.5);
}

void nanoAnalysis::MakeBranch(TTree* t)
{
  t->Branch("Event_No", &Event_No, "Event_No/I");
  t->Branch("Step", &Step, "Step/I");
  t->Branch("Dilep", "TLorentzVector", &Dilep);
  t->Branch("Mu1", "TLorentzVector", &Mu1);
  t->Branch("Mu2", "TLorentzVector", &Mu2);
  t->Branch("Nu_Mu", &Nu_Mu, "Nu_Mu/I");
  t->Branch("Mu_Pt", &Mu_Pt);
  t->Branch("Mu_Eta", &Mu_Eta);
  // t->Branch("Nu_El", &Nu_El, "Nu_El/I");
  // t->Branch("El_Pt", &El_Pt);
  // t->Branch("El_Eta", &El_Eta);
  // t->Branch("Nu_Jet", &Nu_Jet, "Nu_Jet/I");
  // t->Branch("Jet_Pt", &Jet_Pt);
  // t->Branch("Jet_Eta", &Jet_Eta);
  // t->Branch("Nu_BJet", &Nu_BJet, "Nu_BJet/I");
  t->Branch("genweight", &genweight, "genweight/F");
  t->Branch("puweight", &puweight, "puweight/F");
  t->Branch("PV_npvs", &PV_npvs, "PV_npvs/I");
}

void nanoAnalysis::ResetBranch()
{
  Event_No = 0;
  Step = 0;
  Dilep.SetPtEtaPhiM(0,0,0,0);
  GenLep1.SetPtEtaPhiM(0,0,0,0);
  GenLep2.SetPtEtaPhiM(0,0,0,0);
  Mu1.SetPtEtaPhiM(0,0,0,0);
  Mu2.SetPtEtaPhiM(0,0,0,0);
  Nu_Mu = 0;
  Mu_Pt.clear();
  Mu_Eta.clear();
  Mu_Phi.clear();
  Mu_M.clear();
  Mu_Charge.clear();
  Event_Total = 1;
}

inline Bool_t nanoAnalysis::LumiCheck()
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

void nanoAnalysis::fillMcBranch()
{
  for (UInt_t i = 0; i < nGenDressedLepton; i++){
    if(std::abs(GenDressedLepton_pdgId[i]) != 13) break;
    Bool_t bosonSample = false;
    Bool_t isFromBoson = false;
    for(UInt_t j = 0; j < nGenPart; j++)
      {
	if(GenPart_genPartIdxMother[j] == 23 || GenPart_genPartIdxMother[j] == 25)
	  {
	    bosonSample = true;
	    isFromBoson = true;
	  }
      }
    if (isFromBoson)
      {
	if(GenDressedLepton_pdgId[i] == 13)
	  GenLep1.SetPtEtaPhiM(GenDressedLepton_pt[i], GenDressedLepton_eta[i], GenDressedLepton_phi[i], GenDressedLepton_mass[i]);
	else 
	  GenLep2.SetPtEtaPhiM(GenDressedLepton_pt[i], GenDressedLepton_eta[i], GenDressedLepton_phi[i], GenDressedLepton_mass[i]);
      }
  }
}

Double_t nanoAnalysis::roccoR(TLorentzVector m, int q, int nGen, int nTrackerLayers)
{
  Float_t u1 = gRandom->Rndm();
  Float_t u2 = gRandom->Rndm();
  if (!isMC){
    return rocCor->kScaleDT(q, m.Pt(), m.Eta(), m.Phi(), 0, 0);
  }
  else {
    if (nGen > -1){
      return rocCor->kScaleFromGenMC(q, m.Pt(), m.Eta(), m.Phi(),
				     nTrackerLayers, GenPart_pt[nGen],
				     u1, 0, 0);
    }
    else
      return rocCor->kScaleAndSmearMC(q, m.Pt(), m.Eta(), m.Phi(),
				      nTrackerLayers, u1, u2, 0, 0);
  }
  return 1.0;
}
