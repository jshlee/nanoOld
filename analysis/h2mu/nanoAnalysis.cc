#define nanoAnalysis_cxx
#include "nanoAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;
/*
To compile:
g++ `root-config --cflags --glibs` nanoAnalysis.cc -o nanoAnalysis
*/

void nanoAnalysis::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   float nPassTrig = 0;
   float nPassTrig1 = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (HLT_IsoMu24 || HLT_IsoTkMu24)
	nPassTrig++;

      bool IsoMu24 = false;
      bool IsoTkMu24 = false;
      for (int i = 0; i < nTrigObj; ++i){
	if (TrigObj_id[i] != 13) continue;
	if (TrigObj_pt[i] < 24) continue;
	int bits = TrigObj_filterBits[i];
	if (bits & 0x2) IsoMu24 = true;
	if (bits & 0x8) IsoTkMu24 = true;	
      }
      if (IsoMu24 || IsoTkMu24)
	nPassTrig1++;
      // if (Cut(ientry) < 0) continue;
   }
   
   cout <<"pass hlt full "<< nPassTrig1 <<endl;
   cout <<"pass hlt      "<< nPassTrig <<endl;
   cout <<"total hlt     "<< nentries <<endl;
   cout <<"pass hlt frac "<< nPassTrig/float(nentries) <<endl;
   cout <<"pass hlt frac "<< nPassTrig1/float(nentries) <<endl;
}

int main()
{

  TFile *f = new TFile("/home/nanoAOD/run2_2016v3/SingleMuon/Run2016D-18Apr2017-v1/180117_173704/0000/nanoAOD_74.root");
  
  TTree *tree;
  f->GetObject("Events", tree);
  
  nanoAnalysis t(tree);
  t.Loop();


  return 0;
}
