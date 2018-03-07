import ROOT, math, os, copy, getopt, sys, array, math, glob
from ROOT import * 
from array import array

# TTree
f = ROOT.TFile("tmva_class_example.root", "recreate")
sigTree = ROOT.TTree("TreeS", "TreeS")
bkgTree = ROOT.TTree("TreeB", "TreeB")

# Variables
#event_no = array("i",[0])
cme_dca = array("f",[0])
cme_angleXY = array("f",[0])
cme_angleXYZ = array("f",[0])
cme_trk_normalizedChi2 = array("f",[0])
cme_trk_pt = array("f",[0])
cme_trk_ipsigXY = array("f",[0])
cme_trk_ipsigZ = array("f",[0])
cme_lxy = array("f",[0])
cme_lxySig = array("f",[0])
cme_l3D = array("f",[0])
cme_l3DSig = array("f",[0])
cme_jetDR = array("f",[0])
cme_legDR = array("f",[0])
cme_diffMass = array("f",[0])
cme_nJet = array("i",[0])
cme_mcMatch = array("i",[0])
cme_trk_nHits = array("i",[0])
cme_chi2 = array("f",[0])
cme_eta = array("f",[0])
cme_mass = array("f",[0])
cme_pt = array("f",[0])
cme_x = array("f",[0])
cme_y = array("f",[0])
cme_z = array("f",[0])
cme_ndof = array("i",[0])
cme_pdgId = array("i",[0])


# Branches
#sigTree.Branch("event_no", event_no, "event_no/I")
sigTree.Branch("cme_dca", cme_dca, "cme_dca/F")
sigTree.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
sigTree.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
sigTree.Branch("cme_trk_normalizedChi2", cme_trk_normalizedChi2, "cme_trk_normalizedChi2/F")
sigTree.Branch("cme_trk_pt", cme_trk_pt, "cme_trk_pt/F")
sigTree.Branch("cme_trk_ipsigXY", cme_trk_ipsigXY, "cme_trk_ipsigXY/F")
sigTree.Branch("cme_trk_ipsigZ", cme_trk_ipsigZ, "cme_trk_ipsigZ/F")
sigTree.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
sigTree.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
sigTree.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
sigTree.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
sigTree.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
sigTree.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
sigTree.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
sigTree.Branch("cme_nJet", cme_nJet, "cme_nJet/F")
sigTree.Branch("cme_mcMatch", cme_mcMatch, "cme_mcMatch/F")
sigTree.Branch("cme_trk_nHits", cme_trk_nHits, "cme_trk_nHits/F")
sigTree.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
sigTree.Branch("cme_eta", cme_eta, "cme_eta/F")
sigTree.Branch("cme_mass", cme_mass, "cme_mass/F")
sigTree.Branch("cme_pt", cme_pt, "cme_pt/F")
sigTree.Branch("cme_x", cme_x, "cme_x/F")
sigTree.Branch("cme_y", cme_y, "cme_y/F")
sigTree.Branch("cme_z", cme_z, "cme_z/F")
sigTree.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
sigTree.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")

#bkgTree.Branch("event_no", event_no, "event_no/I")
bkgTree.Branch("cme_dca", cme_dca, "cme_dca/F")
bkgTree.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
bkgTree.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
bkgTree.Branch("cme_trk_normalizedChi2", cme_trk_normalizedChi2, "cme_trk_normalizedChi2/F")
bkgTree.Branch("cme_trk_pt", cme_trk_pt, "cme_trk_pt/F")
bkgTree.Branch("cme_trk_ipsigXY", cme_trk_ipsigXY, "cme_trk_ipsigXY/F")
bkgTree.Branch("cme_trk_ipsigZ", cme_trk_ipsigZ, "cme_trk_ipsigZ/F")
bkgTree.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
bkgTree.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
bkgTree.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
bkgTree.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
bkgTree.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
bkgTree.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
bkgTree.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
bkgTree.Branch("cme_nJet", cme_nJet, "cme_nJet/F")
bkgTree.Branch("cme_mcMatch", cme_mcMatch, "cme_mcMatch/F")
bkgTree.Branch("cme_trk_nHits", cme_trk_nHits, "cme_trk_nHits/F")
bkgTree.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
bkgTree.Branch("cme_eta", cme_eta, "cme_eta/F")
bkgTree.Branch("cme_mass", cme_mass, "cme_mass/F")
bkgTree.Branch("cme_pt", cme_pt, "cme_pt/F")
bkgTree.Branch("cme_x", cme_x, "cme_x/F")
bkgTree.Branch("cme_y", cme_y, "cme_y/F")
bkgTree.Branch("cme_z", cme_z, "cme_z/F")
bkgTree.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
bkgTree.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")


events = ROOT.TChain("Events")
events.Add("/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0000/nanoAOD_*.root")
for iev, event in enumerate(events):

    if event.ncmeson == 0:
        continue
        
    for k in range(event.ncmeson):
                    
        cme_dca[0] = event.cmeson_dca[k]
        cme_angleXY[0] = event.cmeson_angleXY[k]               
        cme_angleXYZ[0] = event.cmeson_angleXYZ[k]               
        cme_trk_normalizedChi2[0] = event.cmeson_trk_normalizedChi2[k]               
        cme_trk_pt[0] = event.cmeson_trk_pt[k]               
        cme_trk_ipsigXY[0] = event.cmeson_trk_ipsigXY[k]               
        cme_trk_ipsigZ[0] = event.cmeson_trk_ipsigZ[k]               
        cme_lxy[0] = event.cmeson_lxy[k]               
        cme_lxySig[0] = event.cmeson_lxySig[k]               
        cme_l3D[0] = event.cmeson_l3D[k]               
        cme_l3DSig[0] = event.cmeson_l3DSig[k]               
        cme_jetDR[0] = event.cmeson_jetDR[k]               
        cme_legDR[0] = event.cmeson_legDR[k]               
        cme_diffMass[0] = event.cmeson_diffMass[k]               
        cme_nJet[0] = event.cmeson_nJet[k]               
        cme_mcMatch[0] = event.cmeson_mcMatch[k]               
        cme_trk_nHits[0] = event.cmeson_trk_nHits[k]               
        cme_chi2[0] = event.cmeson_chi2[k]               
        cme_eta[0] = event.cmeson_eta[k]               
        cme_pt[0] = event.cmeson_pt[k]               
        cme_x[0] = event.cmeson_x[k]               
        cme_y[0] = event.cmeson_y[k]               
        cme_z[0] = event.cmeson_z[k]               
        cme_ndof[0] = event.cmeson_ndof[k]               
        cme_pdgId[0] = event.cmeson_pdgId[k]               

        if cme_mcMatch[0] > 1:
            sigTree.Fill()
        else:
            bkgTree.Fill()
    
f.Write()
f.Close()



        


