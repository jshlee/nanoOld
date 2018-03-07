import ROOT, math, os, copy, getopt, sys, array, math, glob
from ROOT import * 
from array import array




## Print TTree
f = ROOT.TFile(cmeson, "recreate")
Cme = ROOT.TTree("nEvent", "nEvent")
#Cme.Branch("cmeson", "TLorentzVector", cmeson)
Cme.Branch("cmeson_dca", cmeson_dca)
Cme.Branch("cmeson_angleXY", cmeson_angleXY)
Cme.Branch("cmeson_angleXYZ", cmeson_angleXYZ)
Cme.Branch("cmeson_trk_normalizedChi2", cmeson_trk_normalizedChi2)
Cme.Branch("cmeson_trk_pt", cmeson_trk_pt)
Cme.Branch("cmeson_trk_ipsigXY", cmeson_trk_ipsigXY)
Cme.Branch("cmeson_trk_ipsigZ", cmeson_trk_ipsigZ)
Cme.Branch("cmeson_lxy", cmeson_lxy)
Cme.Branch("cmeson_lxySig", cmeson_lxySig)
Cme.Branch("cmeson_l3D", cmeson_l3D)
Cme.Branch("cmeson_l3DSig", cmeson_l3DSig)
Cme.Branch("cmeson_jetDR", cmeson_jetDR)
Cme.Branch("cmeson_legDR", cmeson_legDR)
Cme.Branch("cmeson_diffMass", cmeson_diffMass)
Cme.Branch("cmeson_nJet", cmeson_nJet)
Cme.Branch("cmeson_mcMatch", cmeson_mcMatch)
Cme.Branch("cmeson_trk_nHits", cmeson_trk_nHits)
Cme.Branch("cmeson_chi", cmeson_Chi)
Cme.Branch("cmeson_eta", cmeson_eta)
Cme.Branch("cmeson_mass", cmeson_mass)
Cme.Branch("cmeson_pt", cmeson_pt)
Cme.Branch("cmeson_x", cmeson_x)
Cme.Branch("cmeson_y", cmeson_y)
Cme.Branch("cmeson_z", cmeson_z)
Cme.Branch("cmeson_ndof", cmeson_ndof)
Cme.Branch("cmeson_pdgId", cmeson_pdgId)



def D0Selection (event):
    d0s = []
    for i in range(event.ncmeson):
        if event.cmeson_pdgId[i] != 421: continue
        d0 = ROOT.TLorentzVector()
        d0.SetPtEtaPhiM(event.cmeson_pt[i], event.cmeson_eta[i], event.cmeson_phi[i], event.cmeson_mass[i])
        d0s.append(d0)
    return d0s

inFile = ROOT.TFile("/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180124_054937/0000/nanoAOD_157.root")

events = inFile.Get("Events")
nevents += events.GetEntries()




'''
cmesondir0 = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180124_054937/0000/"

samples = [cmesondir0]

hlsit = []

for j, sampledir in enumerate(samples):
    nevents = 0
    filelist = [l for l in os.listdir(samples) if "root" in l]
    for i, fileName in enumerate():
        print fileName
        inFile = ROOT.TFile(samledir+fileName):
        events = inFile.Get("Events")
        nevents += events.GetEntries()p

        if i == 3: break

        

        if ncmeson == 1: 
            Cme.Fill()

f.Write()
f.Close()

'''

