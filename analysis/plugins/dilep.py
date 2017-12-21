import ROOT, math, os, copy, getopt, sys, array, math, glob
from ROOT import *
from array import array


## --- Variables --- ##
#dilep = ROOT.TLorentzVector()
#mu1 = ROOT.TLorentzVector()
#mu2 = ROOT.TLorentzVector()
#    
#jet_CSVV2 = ROOT.std.vector('float')()
#
#Event_N = array("i",[0])
#N_mu = array("i",[0])
#N_el = array("i",[0])
#N_jet = array("i",[0])
#N_bjet = array("i",[0])


def selMuon(event):
    charges =[]
    muons = []
    for i in range(event.nMuon):
        if abs(event.Muon_pt[i]) < 20 : continue
        if abs(event.Muon_eta[i]) > 2.4 : continue
        if not (event.Muon_tightId[i]) : continue
        if event.Muon_pfRelIso04_all[i] > 0.15 : continue

        mu = ROOT.TLorentzVector(event.Muon_pt[i], event.Muon_eta[i], event.Muon_phi[i], event.Muon_mass[i])
        muons.append(mu)
        charges.append(event.Muon_charge[i])
    return muons, charges

def selElec(event):
    charges =[]
    elecs = []
    for i in range(event.nElectron):
        if abs(event.Electron_pt[i]) < 20 : continue
        if abs(event.Electron_eta[i]) > 2.4 : continue
        if not (event.Electron_cutBased[i]) : continue
        if event.Electron_pfRelIso03_all[i] > 0.12 : continue

        elec = ROOT.TLorentzVector(event.Electron_pt[i], event.Electron_eta[i], event.Electron_phi[i], event.Electron_mass[i])
        elecs.append(elec)
        charges.append(event.Electron_charge[i])
    return elecs, charges

def selJet(event):
    jets = []
    for i in range(event.nJet):
        if event.Jet_pt[i] < 30: continue
        if abs(event.Jet_eta[i]) > 2.4: continue
        if not event.Jet_jetId[i] : continue
        """
        bool hasOverLap = false
        for lep in recolep:
            if (deltaR(jet.p4(),lep.p4()) < 0.4) hasOverLap = true;
        if (hasOverLap) continue;
        """

        jet = ROOT.TLorentzVector(event.Jet_pt[i], event.Jet_eta[i], event.Jet_phi[i], 0)
        jets.append(jet)
    return jets

def selBJet(event):
    bjets = []
    for i in range(event.nJet):
        if event.Jet_pt[i] < 20: continue
        if event.Jet_btagCSVV2[i] < 0.848: continue

        bjet = ROOT.TLorentzVector(event.Jet_pt[i], event.Jet_eta[i], event.Jet_phi[i], 0)
        bjets.append(bjet)
    return bjets

def pickD0(event):
    d0s = []
    for i in range(event.ncmeson):
        if event.cmeson_pid[i] != 421: continue
        if event.cmeson_lxy[i] < 0.1: continue
        if event.cmeson_l3D[i] < 0.2: continue

        d0 = ROOT.TLorentzVector()
        d0.SetPtEtaPhiM(event.cmeson_pt[i], event.cmeson_eta[i], event.cmeson_phi[i], event.cmeson_mass[i])
        d0s.append(d0)
    return d0s

    
def pickDstar(event):
    dstars = []
    for i in range(event.ncmeson):
        if event.cmeson_pid[i] != 413: continue
        if event.cmeson_lxy[i] < 0.1: continue
        if event.cmeson_l3D[i] < 0.2: continue

        dstar = ROOT.TLorentzVector(event.cmeson_pt[i], event.cmeson_eta[i], event.cmeson_phi[i], event.cmeson_mass[i])
        dstars.append(dstar)
    return dstars

def pickJpsi(event):
    jpsis = []
    for i in range(event.ncmeson):
        if event.cmeson_pid[i] != 443: continue
#        if event.cmeson_lxy[i] < 0.1: continue
#        if event.cmeson_l3D[i] < 0.2: continue

        jpsi = ROOT.TLorentzVector(event.cmeson_pt[i], event.cmeson_eta[i], event.cmeson_phi[i], event.cmeson_mass[i])
        jpsis.append(jpsi)
    return jpsis

#tstack = ROOT.THStack("lepton mass", "lepton mass")
tstack = ROOT.THStack("cmeson mass", "cmeson mass")

datalumi = 38*1000

ttbardir = "/xrootd/store/group/nanoAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/run2_2016MC_NANO_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/171111_222906/0000/"
dydir10to50 = "/xrootd/store/group/nanoAOD/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/run2_2016MC_NANO_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/171111_220746/0000/"
dydir1 = "/xrootd/store/group/nanoAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/run2_2016MC_NANO_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/171111_220522/0000/"
dydir2 = "/xrootd/store/group/nanoAOD/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/run2_2016MC_NANO_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/171111_220746/0001/"

#samples = [ttbardir, dydir10to50, dydir1]
samples = [ttbardir]
xsec = [831.76, 18610, 6025.2]
#samples = [ttbardir, dydir]

hlist = []
for j, sampledir in enumerate(samples):
    #h = copy.deepcopy(ROOT.TH1D("lepton mass", "lepton mass", 60, 20, 320))
    h = copy.deepcopy(ROOT.TH1D("cmeson mass", "cmeson mass", 30, 0, 60))
    nevents =0
    filelist = [l for l in os.listdir(sampledir) if "root" in l]
    for i, fileName in enumerate(filelist):
        print fileName
        inFile = ROOT.TFile(sampledir+fileName)
        events = inFile.Get("Events")
        nevents += events.GetEntries()

        if i == 3 : break

        for iev, event in enumerate(events):
            #if iev == 1000: break
            selectedMuons, mucharges = selMuon(event)
            #selectedElecs, eleccharges = selElec(event)
            
            #if len(selectedMuons) or len(selectedElecs) < 1: continue
            #recoleps = selectedMuons+selectedElecs
            #charges = mucharges+eleccharges
            recoleps = selectedMuons
            charges = mucharges

            #selectedJets = selJet(event, recoleps)
            selectedBJets = selBJet(event)
            selectedJets = selJet(event)
      
            selecD0 = pickD0(event)
            #D0 = selecD0[0]+selecD0[1]
            #print len(selecD0)

            if len(selectedMuons) != 2 : continue
            dilep = recoleps[0]+recoleps[1]

            #step1
            if abs(dilep.M()) < 20 : continue
            if charges[0]*charges[1] > 0 : continue
            step1 = True
            #h.Fill(abs(dilep.M()))

            #step2
            if 76 < abs(dilep.M()) < 106 : continue
            step2 = True
            #h.Fill(abs(dilep.M()))
            
            #step3
            if event.MET_pt < 40 : continue
            step3 = True
            #h.Fill(abs(dilep.M()))

            #step4
            if len(selectedJets) < 2 : continue
            step4 = True
            #h.Fill(abs(dilep.M()))

            #step5
            if len(selectedBJets) < 1 : continue
            step5 = True
            #h.Fill(abs(dilep.M()))
                   
            ## --- Begin D0 --- ##
            if len(selecD0) < 1 : continue
            k = len(selecD0)
            for m in range(k):
                D0 = selecD0[m]
            h.Fill(abs(D0.M()))

 


        inFile.Close()

    h.SetTitle(sampledir)
    scale = datalumi*xsec[j]/float(nevents)
    print scale
    h.Scale(scale)
    hlist.append(h)

"""
mumudir = "/xrootd/store/user/djeon/NanoTest/DoubleMuon/DoubleMuon2017C_2/171109_205048/0000/"
rdhist = copy.deepcopy(ROOT.TH1D("rdhist", "rdhist", 60, 20, 320))
for i, fileName in enumerate(os.listdir(mumudir)):
    inFile = ROOT.TFile(mumudir+fileName)
    events = inFile.Get("Events")
    for iev, event in enumerate(events):
        selectedMuons, mucharges = selMuon(event)
        selectedElecs, eleccharges = selElec(event)

        recoleps = selectedMuons
        charges = mucharges
        if len(selectedMuons) != 2 : continue
        dilep = recoleps[0]+recoleps[1]

        #step1
        if abs(dilep.M()) < 20 : continue
        if charges[0]*charges[1] > 0 : continue
        step1 = True
        rdhist.Fill(abs(dilep.M()))

    inFile.Close()
"""

outFile = ROOT.TFile("out2.root", "RECREATE")
#ALL = ROOT.TTree("nEvent", "nEvent")
#Cat1 = ROOT.TTree("Cat1", "Cat1")
#Cat2 = ROOT.TTree("Cat2", "Cat2")
#Cat3 = ROOT.TTree("Cat3", "Cat3")
#Cat4 = ROOT.TTree("Cat4", "Cat4")
#
#ALL.Branch("nEvent",nEvent, "nEvent/I")
#ALL.Branch("dilep", "TLorentzVector", dilep)
#ALL.Branch("mu1", "TLorentzVector", mu1)
#ALL.Branch("mu2", "TLorentzVector", mu2)
#ALL.Branch("N_mu", N_mu, "N_mu/I")
#ALL.Branch("N_el", N_el, "N_el/I")
#ALL.Branch("N_jet", N_jet, "N_jet/I")
#ALL.Branch("N_bjet", N_bjet, "N_bjet/I")
#
#
#Cat1.Branch("nEvent",nEvent, "nEvent/I")
#Cat1.Branch("dilep", "TLorentzVector", dilep)
#Cat1.Branch("mu1", "TLorentzVector", mu1)
#Cat1.Branch("mu2", "TLorentzVector", mu2)
#Cat1.Branch("N_mu", N_mu, "N_mu/I")
#Cat1.Branch("N_el", N_el, "N_el/I")
#Cat1.Branch("N_jet", N_jet, "N_jet/I")
#Cat1.Branch("N_bjet", N_bjet, "N_bjet/I")
#
#
#Cat2.Branch("nEvent",nEvent, "nEvent/I")
#Cat2.Branch("dilep", "TLorentzVector", dilep)
#Cat2.Branch("mu1", "TLorentzVector", mu1)
#Cat2.Branch("mu2", "TLorentzVector", mu2)
#Cat2.Branch("N_mu", N_mu, "N_mu/I")
#Cat2.Branch("N_el", N_el, "N_el/I")
#Cat2.Branch("N_jet", N_jet, "N_jet/I")
#Cat2.Branch("N_bjet", N_bjet, "N_bjet/I")
#
#
#Cat3.Branch("nEvent",nEvent, "nEvent/I")
#Cat3.Branch("dilep", "TLorentzVector", dilep)
#Cat3.Branch("mu1", "TLorentzVector", mu1)
#Cat3.Branch("mu2", "TLorentzVector", mu2)
#Cat3.Branch("N_mu", N_mu, "N_mu/I")
#Cat3.Branch("N_el", N_el, "N_el/I")
#Cat3.Branch("N_jet", N_jet, "N_jet/I")
#Cat3.Branch("N_bjet", N_bjet, "N_bjet/I")
#
#
#Cat4.Branch("nEvent",nEvent, "nEvent/I")
#Cat4.Branch("Event_N",Event_N, "Event_N/I")
#Cat4.Branch("dilep", "TLorentzVector", dilep)
#Cat4.Branch("mu1", "TLorentzVector", mu1)
#Cat4.Branch("mu2", "TLorentzVector", mu2)
#Cat4.Branch("N_mu", N_mu, "N_mu/I")
#Cat4.Branch("N_el", N_el, "N_el/I")
#Cat4.Branch("N_jet", N_jet, "N_jet/I")
#Cat4.Branch("N_bjet", N_bjet, "N_bjet/I")


outFile.cd()

for i, h in enumerate(hlist):
    h.SetLineColor(2+i*2)
    h.SetFillColor(2+i*2)
    if "DY" in h.GetTitle():
        h.SetLineColor(4)
        h.SetFillColor(4)
    tstack.Add(h)
    h.Write()

tstack.Write()

canv = ROOT.TCanvas()
canv.SetLogy()
tstack.Draw("hist")
#rdhist.Draw("samee1")
canv.Print("nana.png")

outFile.Write()
outFile.Close()
