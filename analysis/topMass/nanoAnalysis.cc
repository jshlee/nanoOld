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
  t->Branch("nvertex", &b_nvertex, "nvertex/I");
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("Step", &b_Step, "Step/I");
  t->Branch("Dilep", "TLorentzVector", &b_Dilep);
  //t->Branch("Quadlep", "TLorentzVector", &b_Quadlep);
  t->Branch("Mu1", "TLorentzVector", &b_Mu1);
  t->Branch("Mu2", "TLorentzVector", &b_Mu2);
  t->Branch("lep1", "TLorentzVector", &b_lep1);
  t->Branch("lep2", "TLorentzVector", &b_lep2);
  t->Branch("genweight", &b_genweight, "genweight/F");
  t->Branch("puweight", &b_puweight, "puweight/F");
  t->Branch("weight", &b_weight, "weight/F");
 // t->Branch("PV_npvs", &PV_npvs, "PV_npvs/I");
  t->Branch("channel", &b_channel, "channel/I");
  t->Branch("charge", &b_charge, "charge/I");
  //t->Branch("FL", &b_FL, "FL/I");
  //t->Branch("FH", &b_FH, "FH/I");
  //t->Branch("SL", &b_SL, "SL/I");
  t->Branch("nlep", &b_nlep, "nlep/I");
  t->Branch("nmuon", &b_nmuon, "nmuon/I");
  t->Branch("nelec", &b_nelec, "nelec/I");
  t->Branch("njet", &b_njet, "njet/I");
  t->Branch("nbjet", &b_nbjet, "nbjet/I");
  t->Branch("Met_pt", &b_Met_pt, "Met_pt/F");
  t->Branch("Met_phi", &b_Met_phi, "Met_phi/F");
  t->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
  t->Branch("btagweight", &b_btagweight, "b_tagweight/F");
  t->Branch("b_b_Jet_CSVV2", &b_b_Jet_CSVV2);
  t->Branch("mueffweight_up", &b_mueffweight_up, "mueffweight_up/F");
  t->Branch("mueffweight_dn", &b_mueffweight_dn, "mueffweight_dn/F");
  t->Branch("eleffweight", &b_eleffweight, "eleffweight/F");
  t->Branch("eleffweight_up", &b_eleffweight_up, "eleffweight_up/F");
  t->Branch("eleffweight_dn", &b_eleffweight_dn, "eleffweight_dn/F");
}

void nanoAnalysis::ResetBranch()
{
  b_nvertex = 0;
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
  b_b_Jet_CSVV2.clear(); 
  b_Event_Total = 1;
  b_channel = -1;
  b_nlep = 0; b_nmuon = 0; b_nelec = 0; b_njet = 0; b_nbjet = 0;
  b_charge = 0;  
  b_Met_pt = 0; b_Met_phi = 0;
  b_mueffweight = 1; b_btagweight = 1;
  recoleps.clear();
  // b_mueffweight = 1; b_mueffweight_up = 1; b_mueffweight_down = 1;
 // b_btagweight = 1; b_btagweight_up = 1; b_btagweight_down = 1;
 // b_eleffweight = 1; 
  //b_eleffweight_up = 1; b_eleffweight_down = 1;
 // MET = 0; 

}

void nanoAnalysis::LoadModules(pileUpTool* pileUp, lumiTool* lumi, RoccoR* rocCor)
{
  //Get Modules
  m_rocCor = rocCor;
  m_lumi = lumi;
  m_pileUp = pileUp;
  m_btagSF.initCSVWeight();
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
  // b_Step = 1;
  // h_cutFlow->Fill(1);

  if (abs(PV_z) >= 24.) return;
  if (PV_npvs == 0) return;
  if (PV_ndof < 4) return;
  //b_Step = 2;
  //h_cutFlow->Fill(2);
  

  auto muons = MuonSelection();
  auto elecs = ElectronSelection();
  auto jets = JetSelection();
  auto bjets = BtaggedSelection();
  
  b_Met_pt = MET_pt;
  b_Met_phi = MET_phi;
  b_njet = jets.size();
  b_nbjet = bjets.size();
  b_nlep = muons.size() + elecs.size(); 
  b_nmuon = muons.size(); 
  b_nelec = elecs.size(); 

  if(b_nmuon +b_nelec != 2)
  {  
     return;
  }

  int mulpdg;
  if (muons.size() == 2){
      recolep1 = muons[0];
      recolep2 = muons[1];
      mulpdg = muons[0].GetPdgCode()*muons[1].GetPdgCode();
      b_channel = CH_MUMU;
  }

  if (muons.size() == 1 && elecs.size() == 1){
      recolep1 = muons[0];
      recolep2 = elecs[0];
      mulpdg = muons[0].GetPdgCode()*elecs[0].GetPdgCode();
      b_channel = CH_MUEL;
  }
  if (elecs.size() == 2){
      recolep1 = elecs[0];
      recolep2 = elecs[1];
      mulpdg = elecs[0].GetPdgCode()*elecs[1].GetPdgCode();
      b_channel = CH_ELEL;
  }
  
  recolep1.Momentum(b_lep1);
  recolep2.Momentum(b_lep2);

  recoleps.push_back(b_lep1);
  recoleps.push_back(b_lep2);

  b_Dilep = b_lep1 + b_lep2;
  
  if (b_channel == CH_MUMU){
    if (m_isMC){
      if (!(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ     //mumu
         || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ //
         || HLT_IsoMu24 || HLT_IsoTkMu24)) return; //mu
    }
    if (m_isDL){
      if (!(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ     //mumu
         || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ)) return;
    }
    if (m_isSL_m){
      if ((HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ        //mumu
         || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ)  //
         || !(HLT_IsoMu24 || HLT_IsoTkMu24)) return;  //mu
    }
  }

  if (b_channel == CH_MUEL){
    if (m_isMC){
      if (!( HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL        //emu//
          || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ  //
          || HLT_IsoMu24 || HLT_IsoTkMu24       //mu
          || HLT_Ele27_WPTight_Gsf)) return;    //e
    }
    if (m_isDL){
      if (!(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL         //emu
         || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ)) return;
    }
    if (m_isSL_e) {
      if ((HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL         //emu
         || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ) //
         || !(HLT_Ele27_WPTight_Gsf)      //e 
         || (HLT_IsoMu24 || HLT_IsoTkMu24)) return;   //mu
    }
    if (m_isSL_m) {
      if ((HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL         //emu
         || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ) //
         || (HLT_Ele27_WPTight_Gsf)      //e 
         || !(HLT_IsoMu24 || HLT_IsoTkMu24)) return;   //mu
    }
  }

  if (b_channel == CH_ELEL){
    if (m_isMC){
      if (!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ //ee
          || HLT_Ele27_WPTight_Gsf)) return;  //e
    }
    if (m_isDL){
      if (!HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) return;  //ee
    }
    if (m_isSL_e){
      if (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ  //ee
         || !(HLT_Ele27_WPTight_Gsf)) return;        //e
    }
  } 


 


  b_mueffweight = m_muonSF.getScaleFactor(recolep1, 13, 0)*m_muonSF.getScaleFactor(recolep2, 13, 0);
  b_mueffweight_up = m_muonSF.getScaleFactor(recolep1, 13, +1)*m_muonSF.getScaleFactor(recolep2, 13, +1);
  b_mueffweight_dn = m_muonSF.getScaleFactor(recolep1, 13, -1)*m_muonSF.getScaleFactor(recolep2, 13, -1);

  b_eleffweight    = m_elecSF.getScaleFactor(recolep1, 11, 0)*m_elecSF.getScaleFactor(recolep2, 11,  0);
  b_eleffweight_up = m_elecSF.getScaleFactor(recolep1, 11, +1)*m_elecSF.getScaleFactor(recolep2, 11, +1);
  b_eleffweight_dn = m_elecSF.getScaleFactor(recolep1, 11, -1)*m_elecSF.getScaleFactor(recolep2, 11, -1);
 
  //TParticle mu1;
  //TParticle mu2;
  /*
  for(UInt_t i = 0; i < muons.size(); i++)
  {
    if( ( muons[0].GetPdgCode() * muons[i].GetPdgCode() ) < 0 )
    {
      b_Mu1 = b_Mu_tlv[0];
      b_Mu2 = b_Mu_tlv[i];
      b_charge = 1;
      mu1 = muons[0];
      mu2 = muons[i];
      break;
    }    
  }*/
/*
  if (b_charge==0)
  {
     keep = true;
     return;
  }*/
  //b_Step = 5;
  //h_cutFlow->Fill(5);
  
  if (mulpdg > 0 || b_Dilep.M() < 20.)
  {
     keep = true;
     return;
  }
  b_Step = 1;
  h_cutFlow->Fill(1);
  

/*
  b_Met_pt = MET_pt;
  b_Met_phi = MET_phi;
  b_njet = jets.size();
  b_nbjet = bjets.size();
  b_nlep = muons.size() + elecs.size(); 
  b_nmuon = muons.size(); 
  b_nelec = elecs.size(); 
  */

  if (b_njet == 0 || b_nbjet == 0)
  {
     return;
  }
  /*
  /////////////////////// Full Leptonic //////////////////
  if (muons.size() + elecs.size() == 4)        
  {

    Int_t mulpdg = -1;

    if (muons.size() == 2)
    {
      elecs[0].Momentum(b_lep1);
      elecs[1].Momentum(b_lep2);
      mulpdg = elecs[0].GetPdgCode()*elecs[1].GetPdgCode();
      b_channel = CH_ELEL;
    }

    else if (muons.size() == 3)
    {
      for (UInt_t i = 1; i < muons.size(); i++)
      {
        if (muons[i].Pt() != b_Mu2.Pt())
        {
          mulpdg = muons[i].GetPdgCode()*elecs[0].GetPdgCode();
          muons[i].Momentum(b_lep1);
          break;
        }
      }
      elecs[0].Momentum(b_lep2);
      b_channel = CH_MUEL;
    }
    
    else if(muons.size() == 4)
    {
      for (UInt_t i = 1; i < muons.size(); i++)
      {
        if (muons[i].Pt() != b_Mu2.Pt())
        {
          mulpdg = muons[i].GetPdgCode()*muons[3].GetPdgCode();
          muons[i].Momentum(b_lep1);
         break;
        }
      }
      muons[3].Momentum(b_lep2);
      b_channel = CH_MUMU;
    } 
     
    b_Quadlep = b_lep1 + b_lep2 + b_Mu1 + b_Mu2;
       
    if (mulpdg < 0)
    {
      b_FL = 1;         
    }
  }
  
  /////////////////////// Full Hadronic //////////////////
  if (elecs.size() + muons.size() == 0 && jets.size() >= 3)
  {
     if ((bjets.size() == 1)||(bjets.size() == 2))
     {
       b_FH = 1; 
     }
  }
  if (elecs.size() + muons.size() == 0 && jets.size() >= 4)
  {
     if ((bjets.size() == 1)||(bjets.size() == 2))
     {
       b_FH = 2; 
     }
  }
  if (elecs.size() + muons.size() == 0 && jets.size() >= 2)
  {
     if ((bjets.size() == 1)||(bjets.size() == 2))
     {
       b_FH = 1; 
     }
  }

  /////////////////////// Semi Leptonic //////////////////
  if (muons.size() + elecs.size() == 3)
  {
    if (b_nbjet == 1)
    {
      if (b_njet >= 2)
      {
        b_SL = 1;   
      }
    }
    else if (b_nbjet == 2)
    {
      if (b_njet >=1)
      {
        b_SL = 1;   
      }
    }
  }
  */
  if (b_channel != CH_MUEL && 76. < b_Dilep.M() && b_Dilep.M() < 106.)
  {
     keep = true;
     return;
  }
  b_Step = 2;
  h_cutFlow->Fill(2);
  
  if (b_channel != CH_MUEL && b_Met_pt < 40.)
  {
     keep = true;
     return;
  }
  b_Step = 3;
  h_cutFlow->Fill(3);

  if (b_njet < 2)
  {
     keep = true;
     return;
  }
  b_Step = 4;
  h_cutFlow->Fill(4);
    
  if (b_nbjet < 1)
  {
     keep = true;
     return;
  }
  b_Step = 5;
  h_cutFlow->Fill(5);

  if (!(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ))
  {
     keep = true;
     return;    
  }
  b_Step = 6;
  h_cutFlow->Fill(6);

  b_Event_No = 1;
  keep = true;
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
    keep = false; 
    Analysis();
    if (keep)
    {
      ALL->Fill();
    }
  }
}

int main(Int_t argc, Char_t** argv)
{
  string env = getenv("CMSSW_BASE");
  RoccoR* rocCor = new RoccoR(env+"/src/nano/analysis/data/rcdata.2016.v3/");
  lumiTool* lumi = new lumiTool(env+"/src/nano/analysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt");
  pileUpTool* pileUp = new pileUpTool();
  
  //csvWeight.initCSVWeight(false, "csvv2");
//  bTagWeightL.init(3, "csvv2", BTagEntry::OP_LOOSE , 1);
//  bTagWeightM.init(3, "csvv2", BTagEntry::OP_MEDIUM , 1);
//  bTagWeightT.init(3, "csvv2", BTagEntry::OP_TIGHT , 1);
  
  if(argc != 1)
  {
    string dirName = env+("/src/nano/analysis/topMass/Results/")+argv[1];
    //string dirName = "root://cms-xrdr.sdfarm.kr:1094////xrd/store/user/daniellee/Nano/"+std::string(argv[1])+"/"+std::string(argv[2]);
    string temp = argv[1];
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
    if (!Muon_trackerMu[i] || !Muon_globalMu[i] || !Muon_tightId[i] || Muon_pfRelIso04_all[i] > 0.15) continue;   //<~~~~~~~~~~~~~ Higgs Muon pf isolation > 0.25;  top > 0.15
    //if (!Muon_tightId[i]) continue;
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
    //mom = mom * roccoR(mom, Muon_charge[i], Muon_genPartIdx[i], Muon_nTrackerLayers[i]);

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
    if( Electron_pt[i] < 20 || abs(Electron_eta[i]) > 2.4 ) continue;                                  //<~~~~~~~~~~~~~~ Higgs Electron pt == 10; Higgs Electron eta > 2.5  
    //if( Electron_pfRelIso03_all[i] > 0.15 || Electron_cutBased[i] < 3 ) continue;                      //<~~~~~~~~~~~~~~ Top doesn't use isolation? 
    if (Electron_cutBased[i] < 3) continue;
    float el_scEta = Electron_deltaEtaSC[i] + Electron_eta[i];
    if ( std::abs(el_scEta) > 1.4442 &&  std::abs(el_scEta) < 1.566 ) continue;
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
  float Jet_SF_CSV[19] = {1.0, };
  Bool_t flag;
  for(UInt_t i = 0; i < nJet; i++)
  {
    flag = false;
   // auto JET_LOOSE = (1<<0);
    if(Jet_pt[i] < 30 || abs(Jet_eta[i]) > 2.4) continue;                       //<~~~~~~~~~~~~ Higgs Jet eta > 4.7;
   // if(Jet_pt[i] < 30 || abs(Jet_eta[i]) > 4.7) continue;                       //<~~~~~~~~~~~~ Higgs Jet eta > 4.7;
   //if( Jet_jetId[i] & JET_LOOSE == 0 ) continue;
    if (Jet_jetId[i] < 1) continue;
    //if (Jet_btagCSVV2[i] < 0.8484) continue; 
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
   for (UInt_t iu = 0; iu < 19; iu++)
    {
      Jet_SF_CSV[iu] *= m_btagSF.getSF(jet, Jet_btagCSVV2[i], Jet_hadronFlavour[i], iu);
    } 
  }
  for (UInt_t i =0; i<19; i++) b_csvweights.push_back(Jet_SF_CSV[i]);
  b_btagweight = Jet_SF_CSV[0];
  return jets;
}

vector<TParticle> nanoAnalysis::BtaggedSelection()
{
  vector<TParticle> bJets;
  for (UInt_t i = 0; i < nJet; i++)
  {
    if (Jet_jetId[i] < 1) continue;
    if (Jet_btagCSVV2[i] < 0.8484) continue;
    if (Jet_pt[i] < 30 || abs(Jet_eta[i]) > 2.4) continue;           
    
    TLorentzVector mom;
    mom.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
    bool hasOverLap = false;
    for (auto lep : recoleps){
        if (mom.TLorentzVector::DeltaR(lep) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;  
    auto bjet = TParticle();
    bjet.SetMomentum(mom);    
    b_b_Jet_CSVV2.push_back(Jet_btagCSVV2[i]); 
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
