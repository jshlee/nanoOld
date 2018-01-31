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

inline Bool_t nanoAnalysis::MuonSelection(const UInt_t i)
{
  if (!Muon_trackerMu[i] || !Muon_globalMu[i] || !Muon_tightId[i]) return false;
  TLorentzVector m;
  TLorentzVector mu;
  m.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
  mu = m * MuScaleFactor(i);

  if (mu.Pt() < 20) return false;
  else if (std::abs(mu.Eta()) > 2.4) return false;
  else if (Muon_pfRelIso04_all[i] > 0.25) return false;
  Mu_Pt.push_back(mu.Pt());
  Mu_Eta.push_back(mu.Eta());
  Mu_Charge.push_back(Muon_charge[i]);
  Mu_Phi.push_back(mu.Phi());
  Mu_M.push_back(mu.M());
  return true;
}

inline Double_t nanoAnalysis::MuScaleFactor(const UInt_t i)
{
  Float_t scaleFactor = 1.0;
  Float_t u1 = gRandom->Rndm();
  Float_t u2 = gRandom->Rndm();
  if (!isMC)
  {
    scaleFactor = rocCor->kScaleDT(Muon_charge[i], Muon_pt[i], Muon_eta[i], Muon_phi[i], 0, 0);
  } 
  else
  {
    if (Muon_pt[i] == GenLep1.Pt())      scaleFactor = rocCor->kScaleFromGenMC(Muon_charge[i], 
                                                                               Muon_pt[i], 
                                                                               Muon_eta[i],
                                                                               Muon_phi[i], 
                                                                               Muon_nTrackerLayers[i], 
                                                                               GenLep1.Pt(), 
                                                                               u1, 0, 0);
    else if (Muon_pt[i] == GenLep2.Pt()) scaleFactor = rocCor->kScaleFromGenMC(Muon_charge[i], 
                                                                               Muon_pt[i], 
                                                                               Muon_eta[i], 
                                                                               Muon_phi[i], 
                                                                               Muon_nTrackerLayers[i], 
                                                                               GenLep2.Pt(), 
                                                                               u1, 0, 0);
    else                            scaleFactor = rocCor->kScaleAndSmearMC(Muon_charge[i], 
                                                                           Muon_pt[i], 
                                                                           Muon_eta[i], 
                                                                           Muon_phi[i], 
                                                                           Muon_nTrackerLayers[i], 
                                                                           u1, u2, 0, 0);
  }
  return scaleFactor;
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
    Event_Tot->Fill(0.5, Event_Total);
    cutFlow->Fill(0);

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //Run for MC
    if(isMC)
    {
      Int_t nvtx = Pileup_nTrueInt;
      if (nvtx < hist_mc->GetNbinsX()) puweight = puWeightCalculator->getWeight(nvtx);
      else puweight = 1;
      
      genweight = genWeight;
      genweights->Fill(0.5, genweight);
      b_weight = genweight * puweight;
      weight->Fill(0.5, b_weight);
      for(UInt_t i = 0; i < nGenDressedLepton; i++)
      {
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
    } else
    {
      puweight = 1;
      genweight = 0;
      if(!LumiCheck()) {continue;}
    }
    Step = 1;
    cutFlow->Fill(1);

    if (std::abs(PV_z) >= 24.) {continue;}
    if (PV_npvs == 0) {continue;}
    if (PV_ndof < 4) {continue;}

    Step = 2;
    cutFlow->Fill(2);


    // if (Cut(ientry) < 0) continue;
    for (UInt_t i = 0; i < nMuon; i++)
    {
      if(MuonSelection(i)) Nu_Mu++;
    }
    if(Nu_Mu < 2) {continue;}
    Step = 3;
    cutFlow->Fill(3);

    Bool_t charge = false;
    for(Int_t i = 0; i < Nu_Mu; i++)
    {
      if( (Mu_Pt[0] > 26) || (Mu_Pt[i] > 26) )
      {
        if( ( Mu_Charge[0] * Mu_Charge[i] ) < 0 )
        {
          Mu1.SetPtEtaPhiM(Mu_Pt[0], Mu_Eta[0], Mu_Phi[0], Mu_M[0]);
          Mu2.SetPtEtaPhiM(Mu_Pt[i], Mu_Eta[i], Mu_Phi[i], Mu_M[i]);
          charge = true;
          break;
        }
      }
    }
    if(!charge) {continue;}
    Step = 4;
    cutFlow->Fill(4);

    TLorentzVector Dilep_ = Mu1 + Mu2;
    Dilep.SetPtEtaPhiM(Dilep_.Pt(), Dilep_.Eta(), Dilep_.Phi(), Dilep_.M());
    if (Dilep.M() < 12.) {continue;}
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
    if (!(IsoMu24 || IsoTkMu24)) {continue;}
    Step = 6;
    cutFlow->Fill(6);

    Event_No = 1;
    ALL->Fill();
  }
   
  // cout <<"pass hlt      "<< nPassTrig <<endl;
  // cout <<"total hlt     "<< nentries <<endl;
  // cout <<"pass hlt frac "<< nPassTrig/float(nentries) <<endl;
}

int main(Int_t argc, Char_t** argv)
{
  if(argc != 0)
  {
    std::string env = std::getenv("CMSSW_BASE");
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
