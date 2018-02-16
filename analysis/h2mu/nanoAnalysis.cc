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
g++ `root-config --cflags --glibs` nanoAnalysis.cc -o nanoAnalysis
*/

void nanoAnalysis::SetOutput(std::string outputName)
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
  t->Branch("Nu_Mu", &b_Nu_Mu, "Nu_Mu/I");
  t->Branch("Mu_Pt", &b_Mu_Pt);
  t->Branch("Mu_Eta", &b_Mu_Eta);
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("PV_npvs", &PV_npvs, "PV_npvs/I");
}

void nanoAnalysis::ResetBranch()
{
  b_Event_No = 0;
  b_Step = 0;
  b_Dilep.SetPtEtaPhiM(0,0,0,0);
  b_Mu1.SetPtEtaPhiM(0,0,0,0);
  b_Mu2.SetPtEtaPhiM(0,0,0,0);
  b_Nu_Mu = 0;
  for(UInt_t i = 0; i < 5; i++)
  {
    b_Mu_Pt[i] = 0;
    b_Mu_Eta[i] = 0;
    b_Mu_Phi[i] = 0;
    b_Mu_M[i] = 0;
    b_Mu_Charge[i] = 0;  
  }
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

  if (std::abs(PV_z) >= 24.) return;
  if (PV_npvs == 0) return;
  if (PV_ndof < 4) return;

  b_Step = 2;
  h_cutFlow->Fill(2);

  MuonSelection();
  
  if(b_Nu_Mu < 2) return;
  b_Step = 3;
  h_cutFlow->Fill(3);

  Bool_t charge = false;
  for(Int_t i = 0; i < b_Nu_Mu; i++)
  {
    if( (b_Mu_Pt[0] > 26) || (b_Mu_Pt[i] > 26) )
    {
      if( ( b_Mu_Charge[0] * b_Mu_Charge[i] ) < 0 )
      {
        b_Mu1.SetPtEtaPhiM(b_Mu_Pt[0], b_Mu_Eta[0], b_Mu_Phi[0], b_Mu_M[0]);
        b_Mu2.SetPtEtaPhiM(b_Mu_Pt[i], b_Mu_Eta[i], b_Mu_Phi[i], b_Mu_M[i]);
        charge = true;
        break;
      }
    }
  }
  if(!charge) return;
  b_Step = 4;
  h_cutFlow->Fill(4);

  b_Dilep = b_Mu1 + b_Mu2;
  if (b_Dilep.M() < 12.) return;
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
  if(argc != 1)
  {
    std::string env = std::getenv("CMSSW_BASE");
    RoccoR* rocCor = new RoccoR(env+"/src/nano/analysis/data/rcdata.2016.v3/");
    lumiTool* lumi = new lumiTool(env+"/src/nano/analysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt");
    pileUpTool* pileUp = new pileUpTool();
    std::string dirName = env+("/src/nano/analysis/h2mu/Results/Nano_C_Test/")+argv[1];
    std::string temp = argv[1];
    Bool_t isMC = false;
    Size_t found = temp.find("Run");
    if(found == std::string::npos) isMC = true;
    for(Int_t i = 2; i < argc; i++)
    {
      TFile *f = TFile::Open(argv[i], "read");
    
      TTree *tree;
      f->GetObject("Events", tree);
      temp = argv[i];
      found = temp.find_last_of('/');
      std::string outPutName = dirName+temp.substr(found);
      nanoAnalysis t(tree, isMC);
      t.LoadModules(pileUp, lumi, rocCor);
      t.SetOutput(outPutName);
      t.Loop();
    }
  }
  else
  {
    TFile *f = TFile::Open("/xrootd/store/group/nanoAOD/run2_2016v3/SingleMuon/Run2016D-18Apr2017-v1/180117_173704/0000/nanoAOD_74.root", "read");
    
    TTree *tree;
    f->GetObject("Events", tree);
    
    nanoAnalysis t(tree);
    t.SetOutput("test.root");
    t.Loop();
  }

  return 0;
}

//Object Selection
void nanoAnalysis::MuonSelection()
{
  for(UInt_t i = 0; i < nMuon; i++)
  {
    if (!Muon_trackerMu[i] || !Muon_globalMu[i] || !Muon_tightId[i] || Muon_pfRelIso04_all[i] > 0.25) continue;
    TLorentzVector m;
    m.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    m = m * roccoR(m, Muon_charge[i], Muon_genPartIdx[i], Muon_nTrackerLayers[i]);
    if (m.Pt() < 20) continue;
    else if (std::abs(m.Eta()) > 2.4) continue;
    b_Mu_Pt[b_Nu_Mu] = m.Pt();
    b_Mu_Eta[b_Nu_Mu] = m.Eta();
    b_Mu_Phi[b_Nu_Mu] = m.Phi();
    b_Mu_M[b_Nu_Mu] = m.M();
    b_Mu_Charge[b_Nu_Mu] = Muon_charge[i];
    b_Nu_Mu++;
  }
  return;
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