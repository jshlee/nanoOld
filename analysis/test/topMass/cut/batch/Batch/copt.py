import ROOT, os, getopt, sys, array, math, glob, json
from ROOT import *
from array import array


## Make TTree ##
FileArg = sys.argv
tempdir = FileArg[1]
Dirname = "/cms/scratch/seulgi/nanoAOD/src/nano/analysis/test/topMass/cut/batch/Results/%s/" %tempdir
if not os.path.isdir(Dirname):
cattree = Dirname + temp

f = ROOT.TFile(cattree, "recreate")
D0sig = ROOT.TTree("DO_sig", "D0_sig")
D0bkg = ROOT.TTree("D0_bkg", "D0_bkg")
Dstarsig = ROOT.TTree("Dstar_sig", "Dstar_sig")
Dstarbkg = ROOT.TTree("Dstar_bkg", "Dstar_bkg")
JPsisig = ROOT.TTree("JPsi_sig", "JPsi_sig")
JPsibkg = ROOT.TTree("JPsi_bkg", "JPsi_bkg")

## Variables ##
cme_dca = array("f",[0])
cme_angleXY = array("f",[0])
cme_angleXYZ = array("f",[0])
cme_lxy = array("f",[0])
cme_lxySig = array("f",[0])
cme_l3D = array("f",[0])
cme_l3DSig = array("f",[0])
cme_jetDR = array("f",[0])
cme_legDR = array("f",[0])
cme_diffMass = array("f",[0])
cme_nJet = array("i",[0])
cmeTruth_nMatched = array("i",[0])
cme_chi2 = array("f",[0])
cme_eta = array("f",[0])
cme_phi = array("f",[0])
cme_mass = array("f",[0])
cme_pt = array("f",[0])
cme_x = array("f",[0])
cme_y = array("f",[0])
cme_z = array("f",[0])
cme_ndof = array("i",[0])
cme_pdgId = array("i",[0])

cme_jet_btagCMVA = array("f", [0])
cme_jet_btagCSVV2 = array("f", [0])
cme_jet_btagDeepB = array("f", [0])
cme_jet_btagDeepC = array("f", [0])
cme_dau1_chi2 = array("f", [0])
cme_dau1_ipsigXY = array("f", [0])
cme_dau1_ipsigZ = array("f", [0])
cme_dau1_nHits = array("f", [0])
cme_dau1_pt = array("f", [0])
cme_dau2_chi2 = array("f", [0])
cme_dau2_ipsigXY = array("f", [0])
cme_dau2_ipsigZ = array("f", [0])
cme_dau2_nHits = array("f", [0])
cme_dau2_pt = array("f", [0])

# Branche ##
D0sig.Branch("cme_dca", cme_dca, "cme_dca/F")
D0sig.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
D0sig.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
D0sig.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
D0sig.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
D0sig.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
D0sig.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
D0sig.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
D0sig.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
D0sig.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
D0sig.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
D0sig.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
D0sig.Branch("cme_eta", cme_eta, "cme_eta/F")
D0sig.Branch("cme_mass", cme_mass, "cme_mass/F")
D0sig.Branch("cme_phi", cme_phi, "cme_phi/F")
D0sig.Branch("cme_pt", cme_pt, "cme_pt/F")
D0sig.Branch("cme_x", cme_x, "cme_x/F")
D0sig.Branch("cme_y", cme_y, "cme_y/F")
D0sig.Branch("cme_z", cme_z, "cme_z/F")
D0sig.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
D0sig.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")
D0sig.Branch("cme_jet_btagCMVA", cme_jet_btagCMVA, "cme_jet_btagCMVA/F")
D0sig.Branch("cme_jet_btagCSVV2", cme_jet_btagCSVV2, "cme_jet_btagCSVV2/F")
D0sig.Branch("cme_jet_btagDeepB", cme_jet_btagDeepB, "cme_jet_btagDeepB/F")
D0sig.Branch("cme_jet_btagDeepC", cme_jet_btagDeepC, "cme_jet_btagDeepC/F")
D0sig.Branch("cme_nDau", cme_nDau, "cme_nDau/I")
D0sig.Branch("cme_dau1_chi2", cme_dau1_chi2, "cme_dau1_chi2/F")
D0sig.Branch("cme_dau1_pt", cme_dau1_pt, "cme_dau1_pt/F")
D0sig.Branch("cme_dau1_ipsigXY", cme_dau1_ipsigXY, "cme_dau1_ipsigXY/F")
D0sig.Branch("cme_dau1_ipsigZ", cme_dau1_ipsigZ, "cme_dau1_ipsigZ/F")
D0sig.Branch("cme_dau1_nHits", cme_dau1_nHits, "cme_dau1_nHits/I")
D0sig.Branch("cme_dau2_chi2", cme_dau2_chi2, "cme_dau2_chi2/F")
D0sig.Branch("cme_dau2_pt", cme_dau2_pt, "cme_dau2_pt/F")
D0sig.Branch("cme_dau2_ipsigXY", cme_dau2_ipsigXY, "cme_dau2_ipsigXY/F")
D0sig.Branch("cme_dau2_ipsigZ", cme_dau2_ipsigZ, "cme_dau2_ipsigZ/F")
D0sig.Branch("cme_dau2_nHits", cme_dau2_nHits, "cme_dau2_nHits/I")

D0bkg.Branch("cme_dca", cme_dca, "cme_dca/F")
D0bkg.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
D0bkg.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
D0bkg.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
D0bkg.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
D0bkg.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
D0bkg.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
D0bkg.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
D0bkg.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
D0bkg.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
D0bkg.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
D0bkg.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
D0bkg.Branch("cme_eta", cme_eta, "cme_eta/F")
D0bkg.Branch("cme_mass", cme_mass, "cme_mass/F")
D0bkg.Branch("cme_phi", cme_phi, "cme_phi/F")
D0bkg.Branch("cme_pt", cme_pt, "cme_pt/F")
D0bkg.Branch("cme_x", cme_x, "cme_x/F")
D0bkg.Branch("cme_y", cme_y, "cme_y/F")
D0bkg.Branch("cme_z", cme_z, "cme_z/F")
D0bkg.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
D0bkg.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")
D0bkg.Branch("cme_jet_btagCMVA", cme_jet_btagCMVA, "cme_jet_btagCMVA/F")
D0bkg.Branch("cme_jet_btagCSVV2", cme_jet_btagCSVV2, "cme_jet_btagCSVV2/F")
D0bkg.Branch("cme_jet_btagDeepB", cme_jet_btagDeepB, "cme_jet_btagDeepB/F")
D0bkg.Branch("cme_jet_btagDeepC", cme_jet_btagDeepC, "cme_jet_btagDeepC/F")
D0bkg.Branch("cme_nDau", cme_nDau, "cme_nDau/I")
D0bkg.Branch("cme_dau1_chi2", cme_dau1_chi2, "cme_dau1_chi2/F")
D0bkg.Branch("cme_dau1_pt", cme_dau1_pt, "cme_dau1_pt/F")
D0bkg.Branch("cme_dau1_ipsigXY", cme_dau1_ipsigXY, "cme_dau1_ipsigXY/F")
D0bkg.Branch("cme_dau1_ipsigZ", cme_dau1_ipsigZ, "cme_dau1_ipsigZ/F")
D0bkg.Branch("cme_dau1_nHits", cme_dau1_nHits, "cme_dau1_nHits/I")
D0bkg.Branch("cme_dau2_chi2", cme_dau2_chi2, "cme_dau2_chi2/F")
D0bkg.Branch("cme_dau2_pt", cme_dau2_pt, "cme_dau2_pt/F")
D0bkg.Branch("cme_dau2_ipsigXY", cme_dau2_ipsigXY, "cme_dau2_ipsigXY/F")
D0bkg.Branch("cme_dau2_ipsigZ", cme_dau2_ipsigZ, "cme_dau2_ipsigZ/F")

#Dstarsig.Branch("cme_dca", cme_dca, "cme_dca/F")
#Dstarsig.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
#Dstarsig.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
#Dstarsig.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
#Dstarsig.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
#Dstarsig.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
#Dstarsig.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
#Dstarsig.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
#Dstarsig.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
#Dstarsig.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
#Dstarsig.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
#Dstarsig.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
#Dstarsig.Branch("cme_eta", cme_eta, "cme_eta/F")
#Dstarsig.Branch("cme_mass", cme_mass, "cme_mass/F")
#Dstarsig.Branch("cme_phi", cme_phi, "cme_phi/F")
#Dstarsig.Branch("cme_pt", cme_pt, "cme_pt/F")
#Dstarsig.Branch("cme_x", cme_x, "cme_x/F")
#Dstarsig.Branch("cme_y", cme_y, "cme_y/F")
#Dstarsig.Branch("cme_z", cme_z, "cme_z/F")
#Dstarsig.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
#Dstarsig.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")
#Dstarsig.Branch("cme_jet_btagCMVA", cme_jet_btagCMVA, "cme_jet_btagCMVA/F")
#Dstarsig.Branch("cme_jet_btagCSVV2", cme_jet_btagCSVV2, "cme_jet_btagCSVV2/F")
#Dstarsig.Branch("cme_jet_btagDeepB", cme_jet_btagDeepB, "cme_jet_btagDeepB/F")
#Dstarsig.Branch("cme_jet_btagDeepC", cme_jet_btagDeepC, "cme_jet_btagDeepC/F")
#Dstarsig.Branch("cme_nDau", cme_nDau, "cme_nDau/I")
#Dstarsig.Branch("cme_dau1_chi2", cme_dau1_chi2, "cme_dau1_chi2/F")
#Dstarsig.Branch("cme_dau1_pt", cme_dau1_pt, "cme_dau1_pt/F")
#Dstarsig.Branch("cme_dau1_ipsigXY", cme_dau1_ipsigXY, "cme_dau1_ipsigXY/F")
#Dstarsig.Branch("cme_dau1_ipsigZ", cme_dau1_ipsigZ, "cme_dau1_ipsigZ/F")
#Dstarsig.Branch("cme_dau1_nHits", cme_dau1_nHits, "cme_dau1_nHits/I")
#Dstarsig.Branch("cme_dau2_chi2", cme_dau2_chi2, "cme_dau2_chi2/F")
#Dstarsig.Branch("cme_dau2_pt", cme_dau2_pt, "cme_dau2_pt/F")
#Dstarsig.Branch("cme_dau2_ipsigXY", cme_dau2_ipsigXY, "cme_dau2_ipsigXY/F")
#Dstarsig.Branch("cme_dau2_ipsigZ", cme_dau2_ipsigZ, "cme_dau2_ipsigZ/F")
#
#Dstarbkg.Branch("cme_dca", cme_dca, "cme_dca/F")
#Dstarbkg.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
#Dstarbkg.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
#Dstarbkg.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
#Dstarbkg.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
#Dstarbkg.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
#Dstarbkg.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
#Dstarbkg.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
#Dstarbkg.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
#Dstarbkg.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
#Dstarbkg.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
#Dstarbkg.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
#Dstarbkg.Branch("cme_eta", cme_eta, "cme_eta/F")
#Dstarbkg.Branch("cme_mass", cme_mass, "cme_mass/F")
#Dstarbkg.Branch("cme_phi", cme_phi, "cme_phi/F")
#Dstarbkg.Branch("cme_pt", cme_pt, "cme_pt/F")
#Dstarbkg.Branch("cme_x", cme_x, "cme_x/F")
#Dstarbkg.Branch("cme_y", cme_y, "cme_y/F")
#Dstarbkg.Branch("cme_z", cme_z, "cme_z/F")
#Dstarbkg.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
#Dstarbkg.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")
#Dstarbkg.Branch("cme_jet_btagCMVA", cme_jet_btagCMVA, "cme_jet_btagCMVA/F")
#Dstarbkg.Branch("cme_jet_btagCSVV2", cme_jet_btagCSVV2, "cme_jet_btagCSVV2/F")
#Dstarbkg.Branch("cme_jet_btagDeepB", cme_jet_btagDeepB, "cme_jet_btagDeepB/F")
#Dstarbkg.Branch("cme_jet_btagDeepC", cme_jet_btagDeepC, "cme_jet_btagDeepC/F")
#Dstarbkg.Branch("cme_nDau", cme_nDau, "cme_nDau/I")
#Dstarbkg.Branch("cme_dau1_chi2", cme_dau1_chi2, "cme_dau1_chi2/F")
#Dstarbkg.Branch("cme_dau1_pt", cme_dau1_pt, "cme_dau1_pt/F")
#Dstarbkg.Branch("cme_dau1_ipsigXY", cme_dau1_ipsigXY, "cme_dau1_ipsigXY/F")
#Dstarbkg.Branch("cme_dau1_ipsigZ", cme_dau1_ipsigZ, "cme_dau1_ipsigZ/F")
#Dstarbkg.Branch("cme_dau1_nHits", cme_dau1_nHits, "cme_dau1_nHits/I")
#Dstarbkg.Branch("cme_dau2_chi2", cme_dau2_chi2, "cme_dau2_chi2/F")
#Dstarbkg.Branch("cme_dau2_pt", cme_dau2_pt, "cme_dau2_pt/F")
#Dstarbkg.Branch("cme_dau2_ipsigXY", cme_dau2_ipsigXY, "cme_dau2_ipsigXY/F")
#Dstarbkg.Branch("cme_dau2_ipsigZ", cme_dau2_ipsigZ, "cme_dau2_ipsigZ/F")

JPsisig.Branch("cme_dca", cme_dca, "cme_dca/F")
JPsisig.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
JPsisig.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
JPsisig.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
JPsisig.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
JPsisig.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
JPsisig.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
JPsisig.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
JPsisig.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
JPsisig.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
JPsisig.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
JPsisig.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
JPsisig.Branch("cme_eta", cme_eta, "cme_eta/F")
JPsisig.Branch("cme_mass", cme_mass, "cme_mass/F")
JPsisig.Branch("cme_phi", cme_phi, "cme_phi/F")
JPsisig.Branch("cme_pt", cme_pt, "cme_pt/F")
JPsisig.Branch("cme_x", cme_x, "cme_x/F")
JPsisig.Branch("cme_y", cme_y, "cme_y/F")
JPsisig.Branch("cme_z", cme_z, "cme_z/F")
JPsisig.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
JPsisig.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")
JPsisig.Branch("cme_jet_btagCMVA", cme_jet_btagCMVA, "cme_jet_btagCMVA/F")
JPsisig.Branch("cme_jet_btagCSVV2", cme_jet_btagCSVV2, "cme_jet_btagCSVV2/F")
JPsisig.Branch("cme_jet_btagDeepB", cme_jet_btagDeepB, "cme_jet_btagDeepB/F")
JPsisig.Branch("cme_jet_btagDeepC", cme_jet_btagDeepC, "cme_jet_btagDeepC/F")
JPsisig.Branch("cme_nDau", cme_nDau, "cme_nDau/I")
JPsisig.Branch("cme_dau1_chi2", cme_dau1_chi2, "cme_dau1_chi2/F")
JPsisig.Branch("cme_dau1_pt", cme_dau1_pt, "cme_dau1_pt/F")
JPsisig.Branch("cme_dau1_ipsigXY", cme_dau1_ipsigXY, "cme_dau1_ipsigXY/F")
JPsisig.Branch("cme_dau1_ipsigZ", cme_dau1_ipsigZ, "cme_dau1_ipsigZ/F")
JPsisig.Branch("cme_dau1_nHits", cme_dau1_nHits, "cme_dau1_nHits/I")
JPsisig.Branch("cme_dau2_chi2", cme_dau2_chi2, "cme_dau2_chi2/F")
JPsisig.Branch("cme_dau2_pt", cme_dau2_pt, "cme_dau2_pt/F")
JPsisig.Branch("cme_dau2_ipsigXY", cme_dau2_ipsigXY, "cme_dau2_ipsigXY/F")
JPsisig.Branch("cme_dau2_ipsigZ", cme_dau2_ipsigZ, "cme_dau2_ipsigZ/F")

JPsibkg.Branch("cme_dca", cme_dca, "cme_dca/F")
JPsibkg.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
JPsibkg.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
JPsibkg.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
JPsibkg.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
JPsibkg.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
JPsibkg.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
JPsibkg.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
JPsibkg.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
JPsibkg.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
JPsibkg.Branch("cme_nJet", cme_nJet, "cme_nJet/I")
JPsibkg.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
JPsibkg.Branch("cme_eta", cme_eta, "cme_eta/F")
JPsibkg.Branch("cme_mass", cme_mass, "cme_mass/F")
JPsibkg.Branch("cme_phi", cme_phi, "cme_phi/F")
JPsibkg.Branch("cme_pt", cme_pt, "cme_pt/F")
JPsibkg.Branch("cme_x", cme_x, "cme_x/F")
JPsibkg.Branch("cme_y", cme_y, "cme_y/F")
JPsibkg.Branch("cme_z", cme_z, "cme_z/F")
JPsibkg.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
JPsibkg.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")
JPsibkg.Branch("cme_jet_btagCMVA", cme_jet_btagCMVA, "cme_jet_btagCMVA/F")
JPsibkg.Branch("cme_jet_btagCSVV2", cme_jet_btagCSVV2, "cme_jet_btagCSVV2/F")
JPsibkg.Branch("cme_jet_btagDeepB", cme_jet_btagDeepB, "cme_jet_btagDeepB/F")
JPsibkg.Branch("cme_jet_btagDeepC", cme_jet_btagDeepC, "cme_jet_btagDeepC/F")
JPsibkg.Branch("cme_nDau", cme_nDau, "cme_nDau/I")
JPsibkg.Branch("cme_dau1_chi2", cme_dau1_chi2, "cme_dau1_chi2/F")
JPsibkg.Branch("cme_dau1_pt", cme_dau1_pt, "cme_dau1_pt/F")
JPsibkg.Branch("cme_dau1_ipsigXY", cme_dau1_ipsigXY, "cme_dau1_ipsigXY/F")
JPsibkg.Branch("cme_dau1_ipsigZ", cme_dau1_ipsigZ, "cme_dau1_ipsigZ/F")
JPsibkg.Branch("cme_dau1_nHits", cme_dau1_nHits, "cme_dau1_nHits/I")
JPsibkg.Branch("cme_dau2_chi2", cme_dau2_chi2, "cme_dau2_chi2/F")
JPsibkg.Branch("cme_dau2_pt", cme_dau2_pt, "cme_dau2_pt/F")
JPsibkg.Branch("cme_dau2_ipsigXY", cme_dau2_ipsigXY, "cme_dau2_ipsigXY/F")
JPsibkg.Branch("cme_dau2_ipsigZ", cme_dau2_ipsigZ, "cme_dau2_ipsigZ/F")


for i, hadFile in enumerate(FileArg[2:]):
    InFile = TNetXNGFile(Nfile)
    events = InFile.Get("Events")

    
    for iev, event in enumerate(events):
        for k in range(event.nhad):
            
            if event.hadTruth_nMatched[k] != 2 or hadTruth_nMatched[k] != 0: continue
            cme_dca[0] = event.had_dca[k]
            cme_angleXY[0] = event.had_angleXY[k]
            cme_angleXYZ[0] = event.had_angleXYZ[k]
            cme_lxy[0] = event.had_lxy[k]
            cme_lxySig[0] = event.had_lxySig[k]
            cme_l3D[0] = event.had_l3D[k]
            cme_l3DSig[0] = event.had_l3DSig[k]
            cme_jetDR[0] = event.had_jetDR[k]
            cme_legDR[0] = event.had_legDR[k]
            cme_diffMass[0] = event.had_diffMass[k]
            cme_nJet[0] = event.had_nJet[k]
            cmeTruth_nMatched[0] = event.hadTruth_nMatched[k]
            cme_chi2[0] = event.had_chi2[k]
            cme_mass[0] = event.had_mass[k]
            cme_eta[0] = event.had_eta[k]
            cme_phi[0] = event.had_phi[k]
            cme_pt[0] = event.had_pt[k]
            cme_x[0] = event.had_x[k]
            cme_y[0] = event.had_y[k]
            cme_z[0] = event.had_z[k]
            cme_ndof[0] = event.had_ndof[k]
            cme_pdgId[0] = event.had_pdgId[k]
            cme_jet_btagCMVA[0] = event.had_jet_btagCMVA[k]
            cme_jet_btagCSVV2[0] = event.had_jet_btagCSVV2[k]
            cme_jet_btagDeepB[0] = event.had_jet_btagDeepB[k]
            cme_jet_btagDeepC[0] = event.had_jet_btagDeepC[k]
            cme_dau1_chi2[0] = event.had_dau1_chi2[k]
            cme_dau1_ipsigXY[0] = event.had_dau1_ipsigXY[k]
            cme_dau1_ipsigZ[0] = event.had_dau1_ipsigZ[k]
            cme_dau1_nHits[0] = event.had_dau1_nHits[k]
            cme_dau1_pt[0] = event.had_dau1_pt[k]
            cme_dau2_chi2[0] = event.had_dau2_chi2[k]
            cme_dau2_ipsigXY[0] = event.had_dau2_ipsigXY[k]
            cme_dau2_ipsigZ[0] = event.had_dau2_ipsigZ[k]
            cme_dau2_nHits[0] = event.had_dau2_nHits[k]
            cme_dau2_pt[0] = event.had_dau2_pt[k]

            if cmeTruth_nMatched[0] == 2:
                if cme_pdgId[0] == 421:
                    D0sig.Fill()
                elif cme_pdgId[0] == 443:
                    Jpsisig.Fill()
            
            else:
                if cme_pdgId[0] == 421:
                    D0bkg.Fill()
                elif cme_pdgId[0] == 443:
                    Jpsibkg.Fill()
                




