#define nanoAnalysis_cxx
#include "nanoAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

using namespace std;
/*
To compile:
g++ `root-config --cflags --glibs` nanoAnalysis.cc -o nanoAnalysis
To run:
nanoAnalysis  /cms/ldap_home/jlee/nanoAOD/src/nano/nanoAOD/data/dataset/dataset_SingleMuon_Run2016H.txt
*/

void nanoAnalysis::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   float nPassTrig = 0;
   float nPassIsoMu24 = 0;
   float nPassIsoTkMu24 = 0;
      
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (HLT_IsoMu24) nPassIsoMu24++;
      if (HLT_IsoTkMu24) nPassIsoTkMu24++;
      
      if (HLT_IsoMu24 || HLT_IsoTkMu24) nPassTrig++;
      
      // if (Cut(ientry) < 0) continue;
   }
   
   cout <<"pass hlt      "<< nPassTrig <<endl;
   cout <<"pass IsoMu24      "<< nPassIsoMu24 <<endl;
   cout <<"pass IsoTkMu24    "<< nPassIsoTkMu24 <<endl;
   cout <<"total hlt     "<< nentries <<endl;
   cout <<"pass hlt frac "<< nPassTrig/float(nentries) <<endl;
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


  nanoAnalysis t(chain);
  t.Loop();    
  
  return 0;
}
