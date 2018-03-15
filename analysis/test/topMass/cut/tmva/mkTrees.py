import ROOT, math, os, copy, getopt, sys, array, math, glob
from ROOT import * 
from array import array
import numpy as np

# TTree
f = ROOT.TFile("tmva_class_example.root", "recreate")
allTree = ROOT.TTree("All", "All")
sigTree0 = ROOT.TTree("TreeS0", "TreeS0")
sigTree1 = ROOT.TTree("TreeS1", "TreeS1")
bkgTree0 = ROOT.TTree("TreeB0", "TreeB0")
bkgTree1 = ROOT.TTree("TreeB1", "TreeB1")

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


# Making Branches
allTree.Branch("cme_dca", cme_dca, "cme_dca/F")
allTree.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
allTree.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
allTree.Branch("cme_trk_normalizedChi2", cme_trk_normalizedChi2, "cme_trk_normalizedChi2/F")
allTree.Branch("cme_trk_pt", cme_trk_pt, "cme_trk_pt/F")
allTree.Branch("cme_trk_ipsigXY", cme_trk_ipsigXY, "cme_trk_ipsigXY/F")
allTree.Branch("cme_trk_ipsigZ", cme_trk_ipsigZ, "cme_trk_ipsigZ/F")
allTree.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
allTree.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
allTree.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
allTree.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
allTree.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
allTree.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
allTree.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
allTree.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
allTree.Branch("cme_mcMatch", cme_mcMatch, "cme_mcMatch/I")
allTree.Branch("cme_trk_nHits", cme_trk_nHits, "cme_trk_nHits/I")
allTree.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
allTree.Branch("cme_eta", cme_eta, "cme_eta/F")
allTree.Branch("cme_mass", cme_mass, "cme_mass/F")
allTree.Branch("cme_pt", cme_pt, "cme_pt/F")
allTree.Branch("cme_x", cme_x, "cme_x/F")
allTree.Branch("cme_y", cme_y, "cme_y/F")
allTree.Branch("cme_z", cme_z, "cme_z/F")
allTree.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
allTree.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")

#sigTree.Branch("event_no", event_no, "event_no/I")
sigTree0.Branch("cme_dca", cme_dca, "cme_dca/F")
sigTree0.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
sigTree0.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
sigTree0.Branch("cme_trk_normalizedChi2", cme_trk_normalizedChi2, "cme_trk_normalizedChi2/F")
sigTree0.Branch("cme_trk_pt", cme_trk_pt, "cme_trk_pt/F")
sigTree0.Branch("cme_trk_ipsigXY", cme_trk_ipsigXY, "cme_trk_ipsigXY/F")
sigTree0.Branch("cme_trk_ipsigZ", cme_trk_ipsigZ, "cme_trk_ipsigZ/F")
sigTree0.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
sigTree0.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
sigTree0.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
sigTree0.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
sigTree0.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
sigTree0.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
sigTree0.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
sigTree0.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
sigTree0.Branch("cme_mcMatch", cme_mcMatch, "cme_mcMatch/I")
sigTree0.Branch("cme_trk_nHits", cme_trk_nHits, "cme_trk_nHits/I")
sigTree0.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
sigTree0.Branch("cme_eta", cme_eta, "cme_eta/F")
sigTree0.Branch("cme_mass", cme_mass, "cme_mass/F")
sigTree0.Branch("cme_pt", cme_pt, "cme_pt/F")
sigTree0.Branch("cme_x", cme_x, "cme_x/F")
sigTree0.Branch("cme_y", cme_y, "cme_y/F")
sigTree0.Branch("cme_z", cme_z, "cme_z/F")
sigTree0.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
sigTree0.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")

sigTree1.Branch("cme_dca", cme_dca, "cme_dca/F")
sigTree1.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
sigTree1.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
sigTree1.Branch("cme_trk_normalizedChi2", cme_trk_normalizedChi2, "cme_trk_normalizedChi2/F")
sigTree1.Branch("cme_trk_pt", cme_trk_pt, "cme_trk_pt/F")
sigTree1.Branch("cme_trk_ipsigXY", cme_trk_ipsigXY, "cme_trk_ipsigXY/F")
sigTree1.Branch("cme_trk_ipsigZ", cme_trk_ipsigZ, "cme_trk_ipsigZ/F")
sigTree1.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
sigTree1.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
sigTree1.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
sigTree1.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
sigTree1.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
sigTree1.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
sigTree1.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
sigTree1.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
sigTree1.Branch("cme_mcMatch", cme_mcMatch, "cme_mcMatch/I")
sigTree1.Branch("cme_trk_nHits", cme_trk_nHits, "cme_trk_nHits/I")
sigTree1.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
sigTree1.Branch("cme_eta", cme_eta, "cme_eta/F")
sigTree1.Branch("cme_mass", cme_mass, "cme_mass/F")
sigTree1.Branch("cme_pt", cme_pt, "cme_pt/F")
sigTree1.Branch("cme_x", cme_x, "cme_x/F")
sigTree1.Branch("cme_y", cme_y, "cme_y/F")
sigTree1.Branch("cme_z", cme_z, "cme_z/F")
sigTree1.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
sigTree1.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")

#bkgTree.Branch("event_no", event_no, "event_no/I")
bkgTree0.Branch("cme_dca", cme_dca, "cme_dca/F")
bkgTree0.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
bkgTree0.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
bkgTree0.Branch("cme_trk_normalizedChi2", cme_trk_normalizedChi2, "cme_trk_normalizedChi2/F")
bkgTree0.Branch("cme_trk_pt", cme_trk_pt, "cme_trk_pt/F")
bkgTree0.Branch("cme_trk_ipsigXY", cme_trk_ipsigXY, "cme_trk_ipsigXY/F")
bkgTree0.Branch("cme_trk_ipsigZ", cme_trk_ipsigZ, "cme_trk_ipsigZ/F")
bkgTree0.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
bkgTree0.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
bkgTree0.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
bkgTree0.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
bkgTree0.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
bkgTree0.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
bkgTree0.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
bkgTree0.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
bkgTree0.Branch("cme_mcMatch", cme_mcMatch, "cme_mcMatch/I")
bkgTree0.Branch("cme_trk_nHits", cme_trk_nHits, "cme_trk_nHits/I")
bkgTree0.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
bkgTree0.Branch("cme_eta", cme_eta, "cme_eta/F")
bkgTree0.Branch("cme_mass", cme_mass, "cme_mass/F")
bkgTree0.Branch("cme_pt", cme_pt, "cme_pt/F")
bkgTree0.Branch("cme_x", cme_x, "cme_x/F")
bkgTree0.Branch("cme_y", cme_y, "cme_y/F")
bkgTree0.Branch("cme_z", cme_z, "cme_z/F")
bkgTree0.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
bkgTree0.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")

bkgTree1.Branch("cme_dca", cme_dca, "cme_dca/F")
bkgTree1.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
bkgTree1.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
bkgTree1.Branch("cme_trk_normalizedChi2", cme_trk_normalizedChi2, "cme_trk_normalizedChi2/F")
bkgTree1.Branch("cme_trk_pt", cme_trk_pt, "cme_trk_pt/F")
bkgTree1.Branch("cme_trk_ipsigXY", cme_trk_ipsigXY, "cme_trk_ipsigXY/F")
bkgTree1.Branch("cme_trk_ipsigZ", cme_trk_ipsigZ, "cme_trk_ipsigZ/F")
bkgTree1.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
bkgTree1.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
bkgTree1.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
bkgTree1.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
bkgTree1.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
bkgTree1.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
bkgTree1.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
bkgTree1.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
bkgTree1.Branch("cme_mcMatch", cme_mcMatch, "cme_mcMatch/I")
bkgTree1.Branch("cme_trk_nHits", cme_trk_nHits, "cme_trk_nHits/I")
bkgTree1.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
bkgTree1.Branch("cme_eta", cme_eta, "cme_eta/F")
bkgTree1.Branch("cme_mass", cme_mass, "cme_mass/F")
bkgTree1.Branch("cme_pt", cme_pt, "cme_pt/F")
bkgTree1.Branch("cme_x", cme_x, "cme_x/F")
bkgTree1.Branch("cme_y", cme_y, "cme_y/F")
bkgTree1.Branch("cme_z", cme_z, "cme_z/F")
bkgTree1.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
bkgTree1.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")

events = ROOT.TChain("Events")

events.Add("/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0000/nanoAOD_8*.root")
events.Add("/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0001/nanoAOD_12*.root")
events.Add("/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0002/nanoAOD_22*.root")
events.Add("/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0003/nanoAOD_37*.root")
events.Add("/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0004/nanoAOD_40*.root")
events.Add("/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0005/nanoAOD_55*.root")
events.Add("/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0006/nanoAOD_66*.root")
events.Add("/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0007/nanoAOD_71*.root")
events.Add("/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0008/nanoAOD_87*.root")

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

        #print event.cmeson_lxySig[0]

        if np.isnan(event.cmeson_lxySig[0])and event.cmeson_pdgId[0]==421:
            print event.cmeson_lxySig[0]
            print event.cmeson_pdgId[0]




        #if cme_mcMatch[0] > 0:
        #    sigTree0.Fill()
        #    #print "sigTree0 : ", cme_mcMatch[0]
        #if cme_mcMatch[0] > 1:
        #    #print "sigTree1 : ", cme_mcMatch[0]
        #    sigTree1.Fill()
        #if cme_mcMatch[0] < 0:
        #    bkgTree0.Fill()
        #    #print "bkgTree0 : ", cme_mcMatch[0]
        #if cme_mcMatch[0] < 1:
        #    bkgTree1.Fill()
        #    #print "bkgTree1 : ", cme_mcMatch[0]
        ##print cme_mcMatch[0] 
        #allTree.Fill()
    
f.Write()
f.Close()



        


