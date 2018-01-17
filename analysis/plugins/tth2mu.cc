#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "../scripts/WeightCalculatorFromHistogram.cc"
#include "../src/RoccoR.cc"

#include <vector>
#include <cmath>

class TTH2MuAnalyzer {
  
  public:
    TTH2MuAnalyzer(std::string outputFile, std::string evnName);
    ~TTH2MuAnalyzer();

    void LoadLumiMap(std::map<unsigned int, std::vector<std::array<unsigned int, 2>>> rMap);
    void AnalyzeForMC(std::string inputFile);
    void AnalyzeForData(std::string inputFile);
    void EndOfAnalyze();

  private:
    TFile* f;

    std::map<unsigned int, std::vector<std::array<unsigned int, 2>>> lumiMap;
    std::string env;

    //Trees
    TTree* ALL;
    TTree* Cat[10];

    //histogram
    TH1D* Event_Tot;
    TH1D* genweights;
    TH1D* weight;

    //Variables
    TLorentzVector Dilep;
    TLorentzVector Mu1;
    TLorentzVector Mu2;
    TLorentzVector GenLep1;
    TLorentzVector GenLep2;

    std::vector<float> Mu_Pt;
    std::vector<float> Mu_Eta;
    std::vector<float> Mu_Charge;
    std::vector<float> Mu_Phi;
    std::vector<float> Mu_M;

    std::vector<float> El_Pt;
    std::vector<float> El_Eta;

    std::vector<float> Jet_Pt;
    std::vector<float> Jet_Eta;
    std::vector<float> Jet_CSVV2;
    std::vector<float> Jet_M;
    std::vector<float> Jet_Phi;

    int Event_No;
    int Event_Total;
    int Nu_Mu;
    int Nu_El;
    int Nu_Jet;
    int Nu_BJet;
    int Nu_NonBJet;
    float genweight;
    float puweight;
    float b_weight;
    // float GJson;

    //Calculator
    WeightCalculatorFromHistogram* puWeightCalculator;
    RoccoR* rocCor;

    TH1D* hist_data;
    TH1D* hist_mc;

    //flag for discrimination
    bool isMC;

    //private method
    void MakeBranch(TTree* t);
    void LoadModule();

    inline bool MuonSelection(
      float mu_pt,
      float mu_eta,
      float mu_phi,
      float mu_m,
      float mu_iso,
      int mu_charge,
      bool mu_id,
      int nTrack,
      bool Tracker,
      bool Global);
    inline bool ElecSelection(
      float elec_pt,
      float elec_eta,
      float elec_phi,
      float elec_m,
      float elec_iso,
      int elec_id,
      std::vector<TLorentzVector> mu_p);
    inline bool JetSelection(
      float jet_pt,
      float jet_eta,
      float jet_phi,
      float jet_m,
      int jet_id,
      float jet_b,
      std::vector<TLorentzVector> mu_p);
    inline bool BtaggedSelection(
      float jet_pt,
      float jet_eta,
      float jet_CSVV2);
    inline bool LumiCheck(unsigned int run, unsigned int lumiBlock);
    inline double MuScaleFactor(int mu_charge, float mu_pt, float mu_eta, float mu_phi, int nTrack);
};

TTH2MuAnalyzer::TTH2MuAnalyzer(std::string output, std::string envName){
  env = envName;
  LoadModule();
  f = TFile::Open(output.c_str(), "recreate");
  ALL = new TTree("nEvent", "nEvent");
  MakeBranch(ALL);
  for(int i = 0 ; i  < 10; i++){
    char catName[20];
    sprintf(catName, "Cat%d", i+1);
    Cat[i] = new TTree(catName, catName);
    MakeBranch(Cat[i]);
  }

  Event_Tot = new TH1D("Event_total", "Event_total" ,1,0,1);
  genweights = new TH1D("genweight", "genweight" , 1,0,1);
  weight = new TH1D("weight", "weight", 1,0,1);
}

TTH2MuAnalyzer::~TTH2MuAnalyzer(){
  EndOfAnalyze();
}

void TTH2MuAnalyzer::EndOfAnalyze(){
  f->Write();
  f->Close();
}

void TTH2MuAnalyzer::LoadLumiMap(std::map<unsigned int, std::vector<std::array<unsigned int, 2>>> rMap)
{
  lumiMap = rMap;
}

void TTH2MuAnalyzer::MakeBranch(TTree* t)
{
  t->Branch("Event_No", &Event_No, "Event_No/I");
  t->Branch("Dilep", "TLorentzVector", &Dilep);
  t->Branch("Mu1", "TLorentzVector", &Mu1);
  t->Branch("Mu2", "TLorentzVector", &Mu2);
  t->Branch("Nu_Mu", &Nu_Mu, "Nu_Mu/I");
  t->Branch("Mu_Pt", &Mu_Pt);
  t->Branch("Mu_Eta", &Mu_Eta);
  t->Branch("Nu_El", &Nu_El, "Nu_El/I");
  t->Branch("El_Pt", &El_Pt);
  t->Branch("El_Eta", &El_Eta);
  t->Branch("Nu_Jet", &Nu_Jet, "Nu_Jet/I");
  t->Branch("Jet_Pt", &Jet_Pt);
  t->Branch("Jet_Eta", &Jet_Eta);
  t->Branch("Nu_BJet", &Nu_BJet, "Nu_BJet/I");
  t->Branch("genweight", &genweight, "genweight/F");
  t->Branch("puweight", &puweight, "puweight/F");
}

void TTH2MuAnalyzer::LoadModule(){
  //Load Rochester and puWeightCalculator
  std::string temp = env+"/src/nano/analysis/data/pu_root/pileup_profile_Spring16.root";
  hist_mc = (TH1D*)TFile::Open(temp.c_str())->Get("pu_mc");
  hist_mc->SetDirectory(0);
  temp = env+"/src/nano/analysis/data/pu_root/PileupData_GoldenJSON_Full2016.root";
  hist_data = (TH1D*)TFile::Open(temp.c_str())->Get("pileup");
  hist_data->SetDirectory(0);
  puWeightCalculator = new WeightCalculatorFromHistogram(hist_mc, hist_data, true, true, false);
  temp = env+"/src/nano/analysis/data/rcdata.2016.v3/";
  rocCor = new RoccoR(temp);
}

inline double TTH2MuAnalyzer::MuScaleFactor(int mu_charge, float mu_pt, float mu_eta, float mu_phi, int nTrack)
{
  float scaleFactor = 1.0;
  float u1 = gRandom->Rndm();
  float u2 = gRandom->Rndm();
  if (!isMC)
  {
    scaleFactor = rocCor->kScaleDT(mu_charge, mu_pt, mu_eta, mu_phi, 0, 0);
  } else
  {
    if (mu_pt == GenLep1.Pt()) scaleFactor = rocCor->kScaleFromGenMC(mu_charge, mu_pt, mu_eta, mu_phi, nTrack, GenLep1.Pt(), u1, 0, 0);
    else if (mu_pt == GenLep2.Pt()) scaleFactor = rocCor->kScaleFromGenMC(mu_charge, mu_pt, mu_eta, mu_phi, nTrack, GenLep2.Pt(), u1, 0, 0);
    else scaleFactor = rocCor->kScaleAndSmearMC(mu_charge, mu_pt, mu_eta, mu_phi, nTrack, u1, u2, 0, 0);
  }
  return scaleFactor;
}

inline bool TTH2MuAnalyzer::MuonSelection(
  float mu_pt,
  float mu_eta,
  float mu_phi,
  float mu_m,
  float mu_iso,
  int mu_charge,
  bool mu_id,
  int nTrack,
  bool Tracker,
  bool Global )
{
  if (!Tracker || !Global || !mu_id) return false;
  TLorentzVector m;
  TLorentzVector mu;
  m.SetPtEtaPhiM(mu_pt, mu_eta, mu_phi, mu_m);
  mu = m * MuScaleFactor(mu_charge, mu_pt, mu_eta, mu_phi, nTrack);

  if (mu.Pt() < 20) return false;
  else if (std::abs(mu.Eta()) > 2.4) return false;
  else if (mu_iso > 0.25) return false;
  Mu_Pt.push_back(mu.Pt());
  Mu_Eta.push_back(mu.Eta());
  Mu_Charge.push_back(mu_charge);
  Mu_Phi.push_back(mu.Phi());
  Mu_M.push_back(mu.M());
  return true;
}

inline bool TTH2MuAnalyzer::ElecSelection(
  float elec_pt,
  float elec_eta,
  float elec_phi,
  float elec_m,
  float elec_iso,
  int elec_id,
  std::vector<TLorentzVector> mu_p )
{
  if (elec_pt < 10.) return false;
  else if (std::abs(elec_eta) > 2.5) return false;
  else if ((std::abs(elec_eta) > 1.4442) && (std::abs(elec_eta) < 1.566)) return false;
  else if (elec_iso > 0.15) return false;
  else if (elec_id < 3) return false;

  TLorentzVector ElP4;
  ElP4.SetPtEtaPhiM(elec_pt, elec_eta, elec_phi, elec_m);
  for(unsigned int i = 0; i < mu_p.size(); i++)
  {
    if(mu_p[i].DeltaR(ElP4) < 0.4) return false;
  }
  El_Pt.push_back(elec_pt);
  El_Eta.push_back(elec_eta);
  return true;
}

inline bool TTH2MuAnalyzer::JetSelection(
  float jet_pt,
  float jet_eta,
  float jet_phi,
  float jet_m,
  int jet_id,
  float jet_b,
  std::vector<TLorentzVector> mu_p )
{
  int JET_LOOSE = (1<<0);
  if (jet_pt < 30.) return false;
  else if (std::abs(jet_eta) > 4.7) return false;
  else if ((jet_id & JET_LOOSE) == 0) return false;
  TLorentzVector JetP4;
  JetP4.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_m);
  for(unsigned int i = 0; i < mu_p.size(); i++)
  {
    if(mu_p[i].DeltaR(JetP4) < 0.4) return false;
  }
  Jet_Pt.push_back(jet_pt);
  Jet_Eta.push_back(jet_eta);
  Jet_CSVV2.push_back(jet_b);
  Jet_M.push_back(jet_m);
  Jet_Phi.push_back(jet_phi);
  return true;
}

inline bool TTH2MuAnalyzer::BtaggedSelection(
  float jet_pt,
  float jet_eta,
  float jet_CSVV2 )
{
  if (jet_pt < 20.) return false;
  else if(std::abs(jet_eta) > 2.4) return false;
  else if(jet_CSVV2 < 0.848) return false;
  else return true;
}

inline bool TTH2MuAnalyzer::LumiCheck(unsigned int run, unsigned int lumiBlock)
{
  if ( lumiMap.find(run) == lumiMap.end() ) {
    return false;
  } else {
    for (unsigned int i = 0; i < lumiMap[run].size(); i++){
      if(lumiMap[run][i][0] <= lumiBlock && lumiMap[run][i][1] >= lumiBlock) return true;
    }
    return false;
  }
}

void TTH2MuAnalyzer::AnalyzeForMC(std::string inputFile)
{
  isMC = true;
  TFile* inFile = TFile::Open(inputFile.c_str(), "read");

  TTreeReader event("Events", inFile);
  TTreeReaderValue<unsigned int> rNMuon(event, "nMuon");
  TTreeReaderValue<unsigned int> rNElectron(event, "nElectron");
  TTreeReaderValue<unsigned int> rNJet(event, "nJet");
  //Muon Reader
  TTreeReaderArray<float> rMuon_Pt(event, "Muon_pt");
  TTreeReaderArray<float> rMuon_Eta(event, "Muon_eta");
  TTreeReaderArray<float> rMuon_Phi(event, "Muon_phi");
  TTreeReaderArray<float> rMuon_Mass(event, "Muon_mass");
  TTreeReaderArray<float> rMuon_pfRelIson04_all(event, "Muon_pfRelIso04_all");
  TTreeReaderArray<int> rMuon_Charge(event, "Muon_charge");
  TTreeReaderArray<bool> rMuon_MediumId(event, "Muon_mediumId");
  TTreeReaderArray<int> rMuon_nTrackerLayers(event, "Muon_nTrackerLayers");
  TTreeReaderArray<bool> rMuon_TrackerMu(event, "Muon_trackerMu");
  TTreeReaderArray<bool> rMuon_GlobalMu(event, "Muon_globalMu");
  //Electron Reader
  TTreeReaderArray<float> rElectron_Pt(event, "Electron_pt");
  TTreeReaderArray<float> rElectron_Eta(event, "Electron_eta");
  TTreeReaderArray<float> rElectron_Phi(event, "Electron_phi");
  TTreeReaderArray<float> rElectron_Mass(event, "Electron_mass");
  TTreeReaderArray<float> rElectron_pfRelIso03_all(event, "Electron_pfRelIso03_all");
  TTreeReaderArray<int> rElectron_CutBased(event, "Electron_cutBased");
  //Jet Reader
  TTreeReaderArray<float> rJet_Pt(event, "Jet_pt");
  TTreeReaderArray<float> rJet_Eta(event, "Jet_eta");
  TTreeReaderArray<float> rJet_Phi(event, "Jet_phi");
  TTreeReaderArray<float> rJet_Mass(event, "Jet_mass");
  TTreeReaderArray<int> rJet_JetId(event, "Jet_jetId");
  TTreeReaderArray<float> rJet_btagCSVV2(event, "Jet_btagCSVV2");
  //Primary Vertice Reader
  TTreeReaderValue<float> rPV_Z(event, "PV_z");
  TTreeReaderValue<int> rPV_npvs(event, "PV_npvs");
  TTreeReaderValue<float> rPV_ndof(event, "PV_ndof");
  //Triger Reader
  TTreeReaderValue<bool> rHLT_IsoMu24(event, "HLT_IsoMu24");
  TTreeReaderValue<bool> rHLT_IsoTkMu24(event, "HLT_IsoTkMu24");
  //Only For MC
  TTreeReaderValue<int> rPileup_nTrueInt(event, "Pileup_nTrueInt");
  TTreeReaderValue<float> rGenWeight(event, "genWeight");
  
  TTreeReaderValue<unsigned int> rNGenDressedLepton(event, "nGenDressedLepton");
  TTreeReaderArray<int> rGenDerssedLepton_pdgId(event, "GenDressedLepton_pdgId");
  TTreeReaderArray<int> rGenPart_genPartIdxMother(event, "GenPart_genPartIdxMother");
  TTreeReaderArray<float> rGenDressedLepton_pt(event, "GenDressedLepton_pt");
  TTreeReaderArray<float> rGenDerssedLepton_eta(event, "GenDressedLepton_eta");
  TTreeReaderArray<float> rGenDerssedLepton_phi(event, "GenDressedLepton_phi");
  TTreeReaderArray<float> rGenDerssedLepton_mass(event, "GenDressedLepton_mass");

  while(event.Next())
  {
    Mu_Pt.clear();
    Mu_Eta.clear();
    Mu_Charge.clear();
    Mu_Phi.clear();
    Mu_M.clear();

    El_Pt.clear();
    El_Eta.clear();

    Jet_Pt.clear();
    Jet_Eta.clear();
    Jet_CSVV2.clear();
    Jet_M.clear();
    Jet_Phi.clear();

    Dilep.SetPtEtaPhiM(0,0,0,0);
    GenLep1.SetPtEtaPhiM(0,0,0,0);
    GenLep2.SetPtEtaPhiM(0,0,0,0);
    Mu1.SetPtEtaPhiM(0,0,0,0);
    Mu2.SetPtEtaPhiM(0,0,0,0);

    Event_Total = 1;
    Event_Tot->Fill(0.5, Event_Total);
    
    Nu_Mu = 0;
    Nu_El = 0;
    Nu_Jet = 0;
    Nu_BJet = 0;

    int nvtx = int(*rPileup_nTrueInt);
    if (nvtx < hist_mc->GetNbinsX()) puweight = puWeightCalculator->getWeight(nvtx);
    else puweight = 1;
    
    genweight = *rGenWeight;
    genweights->Fill(0.5, genweight);
    b_weight = genweight * puweight;
    weight->Fill(0.5, b_weight);
    for(unsigned int i = 0; i < *rNGenDressedLepton; i++)
    {
      if(std::abs(rGenDerssedLepton_pdgId[i]) != 13) break;
      bool bosonSample = false;
      bool isFromBoson = false;
      for(int GenPart_genPartIdxMother : rGenPart_genPartIdxMother)
      {
        if(GenPart_genPartIdxMother == 23 || GenPart_genPartIdxMother == 25)
        {
          bosonSample = true;
          isFromBoson = true;
        }
      }
      if (isFromBoson)
      {
        if(rGenDerssedLepton_pdgId[i] == 13) GenLep1.SetPtEtaPhiM(rGenDressedLepton_pt[i], rGenDerssedLepton_eta[i], rGenDerssedLepton_phi[i], rGenDerssedLepton_mass[i]);
        else GenLep2.SetPtEtaPhiM(rGenDressedLepton_pt[i], rGenDerssedLepton_eta[i], rGenDerssedLepton_phi[i], rGenDerssedLepton_mass[i]);
      }
    }

    for(unsigned int i = 0; i < *rNMuon; i++)
    {
      if( MuonSelection(rMuon_Pt[i], rMuon_Eta[i], rMuon_Phi[i], rMuon_Mass[i], rMuon_pfRelIson04_all[i], rMuon_Charge[i], rMuon_MediumId[i], rMuon_nTrackerLayers[i], rMuon_TrackerMu[i], rMuon_GlobalMu[i]) )
      {
        Nu_Mu++;
      }
    }
    if (Nu_Mu < 2) continue;

    std::vector<TLorentzVector> Mu_P4;
    for(int i = 0; i < Nu_Mu; i++)
    {
      TLorentzVector MuP4;
      MuP4.SetPtEtaPhiM(Mu_Pt[i], Mu_Eta[i], Mu_Phi[i], Mu_M[i]);
    }

    for(unsigned int i = 0; i < *rNElectron; i++)
    {
      if( ElecSelection(rElectron_Pt[i], rElectron_Eta[i], rElectron_Phi[i], rElectron_Mass[i], rElectron_pfRelIso03_all[i], rElectron_CutBased[i], Mu_P4) )
      {
        Nu_El++;
      }
    }

    for(unsigned int i = 0; i < *rNJet; i++)
    {
      if( JetSelection(rJet_Pt[i], rJet_Eta[i], rJet_Phi[i], rJet_Mass[i], rJet_JetId[i], rJet_btagCSVV2[i], Mu_P4) )
      {
        Nu_Jet++;
      }
    }

    for(int i = 0; i < Nu_Jet; i++)
    {
      if ( BtaggedSelection(Jet_Pt[i], Jet_Eta[i], Jet_CSVV2[i]) )
      {
        Nu_BJet++;
      }
    }

    if ( !(*rHLT_IsoMu24 || *rHLT_IsoTkMu24) ) continue;

    if ( std::abs(*rPV_Z) >= 24. ) continue;
    if ( *rPV_npvs == 0 ) continue;
    if ( *rPV_ndof < 4 ) continue;

    bool Charge = false;
    for (int i = 0; i < Nu_Mu; i++)
    {
      if ( (Mu_Pt[0] > 26) || (Mu_Pt[i] > 26) )
      {
        if ( Mu_Charge[0] * Mu_Charge[i] < 0 )
        {
          Mu1.SetPtEtaPhiM(Mu_Pt[0], Mu_Eta[0], Mu_Phi[0], Mu_M[0]);
          Mu2.SetPtEtaPhiM(Mu_Pt[i], Mu_Eta[i], Mu_Phi[i], Mu_M[i]);
          Charge = true;
          break;
        }
      }
    }
    if (!Charge) continue;

    TLorentzVector Dilep_ = Mu1 + Mu2;
    Dilep.SetPtEtaPhiM(Dilep_.Pt(), Dilep_.Eta(), Dilep_.Phi(), Dilep_.M());
    if (Dilep.M() < 12.) continue;

    if( Nu_BJet == 1 )
    {
      if ( Nu_El == 1 )
      {
        if( Nu_Mu > 2 ) continue;
        else Cat[0]->Fill();
      }
      if( Nu_El == 2 )
      {
        if( Nu_Mu > 2 ) continue;
        else Cat[1]->Fill();
      }
      if( Nu_Mu == 3 )
      {
        if( Nu_El > 0 ) continue;
        else Cat[2]->Fill();
      }
      if( Nu_Mu == 4 )
      {
        if( Nu_El > 0 )continue;
        else Cat[3]->Fill();
      }
      if ( (Nu_Jet - Nu_BJet) == 4 ) Cat[4]->Fill();
    }
    else if( Nu_BJet == 2 )
    {
      if ( Nu_El == 1 )
      {
        if( Nu_Mu > 2 ) continue;
        else Cat[5]->Fill();
      }
      if( Nu_El == 2 )
      {
        if( Nu_Mu > 2 ) continue;
        else Cat[6]->Fill();
      }
      if( Nu_Mu == 3 )
      {
        if( Nu_El > 0 ) continue;
        else Cat[7]->Fill();
      }
      if( Nu_Mu == 4 )
      {
        if( Nu_El > 0 )continue;
        else Cat[8]->Fill();
      }
      if ( (Nu_Jet - Nu_BJet) == 4 ) Cat[9]->Fill();
    }

    Event_No = 1;
    ALL->Fill();
  }
}

void TTH2MuAnalyzer::AnalyzeForData(std::string inputFile)
{
  isMC = false;
  TFile* inFile = TFile::Open(inputFile.c_str(), "read");

  TTreeReader event("Events", inFile);
  TTreeReaderValue<unsigned int> rRun(event, "run");
  TTreeReaderValue<unsigned int> rLumiBlock(event, "luminosityBlock");

  TTreeReaderValue<unsigned int> rNMuon(event, "nMuon");
  TTreeReaderValue<unsigned int> rNElectron(event, "nElectron");
  TTreeReaderValue<unsigned int> rNJet(event, "nJet");
  //Muon Reader
  TTreeReaderArray<float> rMuon_Pt(event, "Muon_pt");
  TTreeReaderArray<float> rMuon_Eta(event, "Muon_eta");
  TTreeReaderArray<float> rMuon_Phi(event, "Muon_phi");
  TTreeReaderArray<float> rMuon_Mass(event, "Muon_mass");
  TTreeReaderArray<float> rMuon_pfRelIson04_all(event, "Muon_pfRelIso04_all");
  TTreeReaderArray<int> rMuon_Charge(event, "Muon_charge");
  TTreeReaderArray<bool> rMuon_MediumId(event, "Muon_mediumId");
  TTreeReaderArray<int> rMuon_nTrackerLayers(event, "Muon_nTrackerLayers");
  TTreeReaderArray<bool> rMuon_TrackerMu(event, "Muon_trackerMu");
  TTreeReaderArray<bool> rMuon_GlobalMu(event, "Muon_globalMu");
  //Electron Reader
  TTreeReaderArray<float> rElectron_Pt(event, "Electron_pt");
  TTreeReaderArray<float> rElectron_Eta(event, "Electron_eta");
  TTreeReaderArray<float> rElectron_Phi(event, "Electron_phi");
  TTreeReaderArray<float> rElectron_Mass(event, "Electron_mass");
  TTreeReaderArray<float> rElectron_pfRelIso03_all(event, "Electron_pfRelIso03_all");
  TTreeReaderArray<int> rElectron_CutBased(event, "Electron_cutBased");
  //Jet Reader
  TTreeReaderArray<float> rJet_Pt(event, "Jet_pt");
  TTreeReaderArray<float> rJet_Eta(event, "Jet_eta");
  TTreeReaderArray<float> rJet_Phi(event, "Jet_phi");
  TTreeReaderArray<float> rJet_Mass(event, "Jet_mass");
  TTreeReaderArray<int> rJet_JetId(event, "Jet_jetId");
  TTreeReaderArray<float> rJet_btagCSVV2(event, "Jet_btagCSVV2");
  //Primary Vertice Reader
  TTreeReaderValue<float> rPV_Z(event, "PV_z");
  TTreeReaderValue<int> rPV_npvs(event, "PV_npvs");
  TTreeReaderValue<float> rPV_ndof(event, "PV_ndof");
  //Triger Reader
  TTreeReaderValue<bool> rHLT_IsoMu24(event, "HLT_IsoMu24");
  TTreeReaderValue<bool> rHLT_IsoTkMu24(event, "HLT_IsoTkMu24");

  while(event.Next())
  {
    Mu_Pt.clear();
    Mu_Eta.clear();
    Mu_Charge.clear();
    Mu_Phi.clear();
    Mu_M.clear();

    El_Pt.clear();
    El_Eta.clear();

    Jet_Pt.clear();
    Jet_Eta.clear();
    Jet_CSVV2.clear();
    Jet_M.clear();
    Jet_Phi.clear();

    Dilep.SetPtEtaPhiM(0,0,0,0);
    GenLep1.SetPtEtaPhiM(0,0,0,0);
    GenLep2.SetPtEtaPhiM(0,0,0,0);
    Mu1.SetPtEtaPhiM(0,0,0,0);
    Mu2.SetPtEtaPhiM(0,0,0,0);

    Event_Total = 1;
    Event_Tot->Fill(0.5, Event_Total);
    
    Nu_Mu = 0;
    Nu_El = 0;
    Nu_Jet = 0;
    Nu_BJet = 0;

    puweight = 1;

    if(!LumiCheck(*rRun, *rLumiBlock)) continue;
    
    for(unsigned int i = 0; i < *rNMuon; i++)
    {
      if( MuonSelection(rMuon_Pt[i], rMuon_Eta[i], rMuon_Phi[i], rMuon_Mass[i], rMuon_pfRelIson04_all[i], rMuon_Charge[i], rMuon_MediumId[i], rMuon_nTrackerLayers[i], rMuon_TrackerMu[i], rMuon_GlobalMu[i]))
      {
        Nu_Mu++;
      }
    }
    if (Nu_Mu < 2) continue;

    std::vector<TLorentzVector> Mu_P4;
    for(int i = 0; i < Nu_Mu; i++)
    {
      TLorentzVector MuP4;
      MuP4.SetPtEtaPhiM(Mu_Pt[i], Mu_Eta[i], Mu_Phi[i], Mu_M[i]);
    }

    for(unsigned int i = 0; i < *rNElectron; i++)
    {
      if( ElecSelection(rElectron_Pt[i], rElectron_Eta[i], rElectron_Phi[i], rElectron_Mass[i], rElectron_pfRelIso03_all[i], rElectron_CutBased[i], Mu_P4) )
      {
        Nu_El++;
      }
    }

    for(unsigned int i = 0; i < *rNJet; i++)
    {
      if( JetSelection(rJet_Pt[i], rJet_Eta[i], rJet_Phi[i], rJet_Mass[i], rJet_JetId[i], rJet_btagCSVV2[i], Mu_P4) )
      {
        Nu_Jet++;
      }
    }

    for(int i = 0; i < Nu_Jet; i++)
    {
      if ( BtaggedSelection(Jet_Pt[i], Jet_Eta[i], Jet_CSVV2[i]) )
      {
        Nu_BJet++;
      }
    }

    if ( !(*rHLT_IsoMu24 || *rHLT_IsoTkMu24) ) continue;

    if (std::abs(*rPV_Z) >= 24.) continue;
    if (*rPV_npvs == 0) continue;
    if (*rPV_ndof < 4) continue;


    bool Charge = false;
    for (int i = 0; i < Nu_Mu; i++)
    {
      if ( (Mu_Pt[0] > 26) || (Mu_Pt[i] > 26) )
      {
        if ( Mu_Charge[0] * Mu_Charge[i] < 0 )
        {
          Mu1.SetPtEtaPhiM(Mu_Pt[0], Mu_Eta[0], Mu_Phi[0], Mu_M[0]);
          Mu2.SetPtEtaPhiM(Mu_Pt[i], Mu_Eta[i], Mu_Phi[i], Mu_M[i]);
          Charge = true;
          break;
        }
      }
    }
    if (!Charge) continue;

    TLorentzVector Dilep_ = Mu1 + Mu2;
    Dilep.SetPtEtaPhiM(Dilep_.Pt(), Dilep_.Eta(), Dilep_.Phi(), Dilep_.M());
    if (Dilep.M() < 12.) continue;

    if( Nu_BJet == 1 )
    {
      if ( Nu_El == 1 )
      {
        if( Nu_Mu > 2 ) continue;
        else Cat[0]->Fill();
      }
      if( Nu_El == 2 )
      {
        if( Nu_Mu > 2 ) continue;
        else Cat[1]->Fill();
      }
      if( Nu_Mu == 3 )
      {
        if( Nu_El > 0 ) continue;
        else Cat[2]->Fill();
      }
      if( Nu_Mu == 4 )
      {
        if( Nu_El > 0 )continue;
        else Cat[3]->Fill();
      }
      if ( (Nu_Jet - Nu_BJet) == 4 ) Cat[4]->Fill();
    }
    else if( Nu_BJet == 2 )
    {
      if ( Nu_El == 1 )
      {
        if( Nu_Mu > 2 ) continue;
        else Cat[5]->Fill();
      }
      if( Nu_El == 2 )
      {
        if( Nu_Mu > 2 ) continue;
        else Cat[6]->Fill();
      }
      if( Nu_Mu == 3 )
      {
        if( Nu_El > 0 ) continue;
        else Cat[7]->Fill();
      }
      if( Nu_Mu == 4 )
      {
        if( Nu_El > 0 )continue;
        else Cat[8]->Fill();
      }
      if ( (Nu_Jet - Nu_BJet) == 4 ) Cat[9]->Fill();
    }

    Event_No = 1;
    ALL->Fill();  
  }
}