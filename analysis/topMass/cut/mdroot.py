import ROOT, math, os, copy, getopt, sys, array, math, glob
from ROOT import * 
from array import array


# TTree
f = ROOT.TFile("1.root", "recreate")
Cme = ROOT.TTree("Event_C", "Event_C")


# Variables
event_no = array("i",[0])
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
Cme.Branch("event_no", event_no, "event_no/I")
Cme.Branch("cme_dca", cme_dca, "cme_dca/F")
Cme.Branch("cme_angleXY", cme_angleXY, "cme_angleXY/F")
Cme.Branch("cme_angleXYZ", cme_angleXYZ, "cme_angleXYZ/F")
Cme.Branch("cme_trk_normalizedChi2", cme_trk_normalizedChi2, "cme_trk_normalizedChi2/F")
Cme.Branch("cme_trk_pt", cme_trk_pt, "cme_trk_pt/F")
Cme.Branch("cme_trk_ipsigXY", cme_trk_ipsigXY, "cme_trk_ipsigXY/F")
Cme.Branch("cme_trk_ipsigZ", cme_trk_ipsigZ, "cme_trk_ipsigZ/F")
Cme.Branch("cme_lxy", cme_lxy, "cme_lxy/F")
Cme.Branch("cme_lxySig", cme_lxySig, "cme_lxySig/F")
Cme.Branch("cme_l3D", cme_l3D, "cme_l3D/F")
Cme.Branch("cme_l3DSig", cme_l3DSig, "cme_l3DSig/F")
Cme.Branch("cme_jetDR", cme_jetDR, "cme_jetDR/F")
Cme.Branch("cme_legDR", cme_legDR, "cme_legDR/F")
Cme.Branch("cme_diffMass", cme_diffMass, "cme_diffMass/F")
Cme.Branch("cme_nJet", cme_nJet, "cme_nJet/F")
Cme.Branch("cme_mcMatch", cme_mcMatch, "cme_mcMatch/F")
Cme.Branch("cme_trk_nHits", cme_trk_nHits, "cme_trk_nHits/F")
Cme.Branch("cme_chi2", cme_chi2, "cme_chi2/F")
Cme.Branch("cme_eta", cme_eta, "cme_eta/F")
Cme.Branch("cme_mass", cme_mass, "cme_mass/F")
Cme.Branch("cme_pt", cme_pt, "cme_pt/F")
Cme.Branch("cme_x", cme_x, "cme_x/F")
Cme.Branch("cme_y", cme_y, "cme_y/F")
Cme.Branch("cme_z", cme_z, "cme_z/F")
Cme.Branch("cme_ndof", cme_ndof, "cme_ndof/I")
Cme.Branch("cme_pdgId", cme_pdgId, "cme_pdgId/I")



cmesondir0 = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180222_144032/0000/"
cmesondir1 = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180222_144032/0001/"
#cmesondir2 = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180222_144032/0002/"
cmesondir3 = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180222_144032/0003/"
#cmesondir4 = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180222_144032/0004/"
#cmesondir5 = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180222_144032/0005/"
#
#samples = [cmesondir0, cmesondir1, cmesondir2, cmesondir3, cmesondir4, cmesondir5]
samples = [cmesondir0, cmesondir1, cmesondir3]

#hlsit = []
ncme = 0
event_no[0] = 0
for j, sampledir in enumerate(samples):
    #nevents = 0
    filelist = [l for l in os.listdir(sampledir) if "root" in l]
    
    #cme_dca[0] = 0
    #cme_angleXY[0] = 0             
    #cme_angleXYZ[0] = 0            
    #cme_trk_normalizedChi2[0] = 0               
    #cme_trk_pt[0] = 0               
    #cme_trk_ipsigXY[0] = 0               
    #cme_trk_ipsigZ[0] = 0              
    #cme_lxy[0] = 0         
    #cme_lxySig[0] = 0              
    #cme_l3D[0] = 0         
    #cme_l3DSig[0] = 0              
    #cme_jetDR[0] = 0         
    #cme_legDR[0] = 0             
    #cme_diffMass[0] = 0              
    #cme_nJet[0] = 0         
    #cme_mcMatch[0] = 0              
    #cme_trk_nHits[0] = 0               
    #cme_chi2[0] = 0             
    #cme_eta[0] = 0       
    #cme_pt[0] = 0        
    #cme_x[0] = 0           
    #cme_y[0] = 0             
    #cme_z[0] = 0            
    #cme_ndof[0] = 0               
    #cme_pdgId[0] = 0
    #
    
    for i, fileName in enumerate(filelist):
        if i == 2: break
        print fileName
        inFile = ROOT.TFile(sampledir+fileName)
        events = inFile.Get("Events")
        #nevents += events.GetEntries()

        #if i == 2: break

        for iev, event in enumerate(events):

            #if event.ncmeson == 0:
            #    continue

            if event.ncmeson > 0:
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

                    Cme.Fill()


            ncme += event.ncmeson
        #Cme.Fill()
        inFile.Close()
        print ncme 
    
f.Write()
f.Close()






