import ROOT, os, getopt, sys, array, math, glob
from ROOT import * 
from array import array

### Rochester ###
ROOT.gROOT.LoadMacro("/cms/scratch/daniel/CMSSW_8_0_26_patch1/src/CATTools/CatAnalyzer/src/RoccoR.cc+")
roc = ROOT.std.string("/cms/scratch/daniel/CMSSW_8_0_26_patch1/src/CATTools/CatAnalyzer/data/rcdata.2016.v3/")
rocCor = ROOT.RoccoR(roc)

### Make TTREE ### 
FileArg = sys.argv
print FileArg
tempdir = FileArg[1]
Dirname = "/cms/scratch/daniel/CMSSW_8_0_26_patch1/src/CATTools/CatAnalyzer/test/Nano_AOD/Results/%s/"%tempdir 
if not os.path.isdir(Dirname):
    os.makedirs(Dirname)

temp = FileArg[2].split('/').pop()
cattree = Dirname+temp

#print cattree
f = ROOT.TFile(cattree, "recreate")
ALL = ROOT.TTree("nEvent", "nEvent")
Cat1 = ROOT.TTree("Cat1", "Cat1")
Cat2 = ROOT.TTree("Cat2", "Cat2")
Cat3 = ROOT.TTree("Cat3", "Cat3")
Cat4 = ROOT.TTree("Cat4", "Cat4")
Cat5 = ROOT.TTree("Cat5", "Cat5")
Cat6 = ROOT.TTree("Cat6", "Cat6")
Cat7 = ROOT.TTree("Cat7", "Cat7")
Cat8 = ROOT.TTree("Cat8", "Cat8")
Cat9 = ROOT.TTree("Cat9", "Cat9")
Cat10 = ROOT.TTree("Cat10", "Cat10")

### Variables ###
Dilep = ROOT.TLorentzVector()
Mu1 = ROOT.TLorentzVector()
Mu2 = ROOT.TLorentzVector()
GenLep1 = ROOT.TLorentzVector()
GenLep2 = ROOT.TLorentzVector()


Mu_Pt = ROOT.std.vector('float')()
Mu_Eta = ROOT.std.vector('float')()
Mu_Charge = ROOT.std.vector('float')()
Mu_Phi = ROOT.std.vector('float')()
Mu_M = ROOT.std.vector('float')()

El_Pt = ROOT.std.vector('float')()
El_Eta = ROOT.std.vector('float')()

Jet_Pt = ROOT.std.vector('float')()
Jet_Eta = ROOT.std.vector('float')()
Jet_CSVV2 = ROOT.std.vector('float')()
Jet_M = ROOT.std.vector('float')()
Jet_Phi = ROOT.std.vector('float')()

Event_No = array("i",[0])
Event_Total = array("i",[0])
Nu_Mu = array("i",[0])
Nu_El = array("i",[0])
Nu_Jet = array("i",[0])
Nu_BJet = array("i",[0])
Nu_NonBJet = array("i",[0])
genweight = array("f",[0])

Event_Tot = ROOT.TH1D("Event_total", "Event_total" ,1,0,1)

### Branches ###
ALL.Branch("Event_No", Event_No, "Event_No/I")
ALL.Branch("Dilep", "TLorentzVector", Dilep)
ALL.Branch("Mu1", "TLorentzVector", Mu1)
ALL.Branch("Mu2", "TLorentzVector", Mu2)
ALL.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
ALL.Branch("Mu_Pt", Mu_Pt)
ALL.Branch("Mu_Eta", Mu_Eta)
ALL.Branch("Nu_El", Nu_El, "Nu_El/I")
ALL.Branch("El_Pt", El_Pt)
ALL.Branch("El_Eta", El_Eta)
ALL.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
ALL.Branch("Jet_Pt", Jet_Pt)
ALL.Branch("Jet_Eta", Jet_Eta)
ALL.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
ALL.Branch("genweight", genweight, "genweight/F")

Cat1.Branch("Event_No", Event_No, "Event_No/I")
Cat1.Branch("Dilep", "TLorentzVector", Dilep)
Cat1.Branch("Mu1", "TLorentzVector", Mu1)
Cat1.Branch("Mu2", "TLorentzVector", Mu2)
Cat1.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
Cat1.Branch("Mu_Pt", Mu_Pt)
Cat1.Branch("Mu_Eta", Mu_Eta)
Cat1.Branch("Nu_El", Nu_El, "Nu_El/I")
Cat1.Branch("El_Pt", El_Pt)
Cat1.Branch("El_Eta", El_Eta)
Cat1.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
Cat1.Branch("Jet_Pt", Jet_Pt)
Cat1.Branch("Jet_Eta", Jet_Eta)
Cat1.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
Cat1.Branch("genweight", genweight, "genweight/F")

Cat2.Branch("Event_No", Event_No, "Event_No/I")
Cat2.Branch("Dilep", "TLorentzVector", Dilep)
Cat2.Branch("Mu1", "TLorentzVector", Mu1)
Cat2.Branch("Mu2", "TLorentzVector", Mu2)
Cat2.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
Cat2.Branch("Mu_Pt", Mu_Pt)
Cat2.Branch("Mu_Eta", Mu_Eta)
Cat2.Branch("Nu_El", Nu_El, "Nu_El/I")
Cat2.Branch("El_Pt", El_Pt)
Cat2.Branch("El_Eta", El_Eta)
Cat2.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
Cat2.Branch("Jet_Pt", Jet_Pt)
Cat2.Branch("Jet_Eta", Jet_Eta)
Cat2.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
Cat2.Branch("genweight", genweight, "genweight/F")

Cat3.Branch("Event_No", Event_No, "Event_No/I")
Cat3.Branch("Dilep", "TLorentzVector", Dilep)
Cat3.Branch("Mu1", "TLorentzVector", Mu1)
Cat3.Branch("Mu2", "TLorentzVector", Mu2)
Cat3.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
Cat3.Branch("Mu_Pt", Mu_Pt)
Cat3.Branch("Mu_Eta", Mu_Eta)
Cat3.Branch("Nu_El", Nu_El, "Nu_El/I")
Cat3.Branch("El_Pt", El_Pt)
Cat3.Branch("El_Eta", El_Eta)
Cat3.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
Cat3.Branch("Jet_Pt", Jet_Pt)
Cat3.Branch("Jet_Eta", Jet_Eta)
Cat3.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
Cat3.Branch("genweight", genweight, "genweight/F")

Cat4.Branch("Event_No", Event_No, "Event_No/I")
Cat4.Branch("Dilep", "TLorentzVector", Dilep)
Cat4.Branch("Mu1", "TLorentzVector", Mu1)
Cat4.Branch("Mu2", "TLorentzVector", Mu2)
Cat4.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
Cat4.Branch("Mu_Pt", Mu_Pt)
Cat4.Branch("Mu_Eta", Mu_Eta)
Cat4.Branch("Nu_El", Nu_El, "Nu_El/I")
Cat4.Branch("El_Pt", El_Pt)
Cat4.Branch("El_Eta", El_Eta)
Cat4.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
Cat4.Branch("Jet_Pt", Jet_Pt)
Cat4.Branch("Jet_Eta", Jet_Eta)
Cat4.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
Cat4.Branch("genweight", genweight, "genweight/F")

Cat5.Branch("Event_No", Event_No, "Event_No/I")
Cat5.Branch("Dilep", "TLorentzVector", Dilep)
Cat5.Branch("Mu1", "TLorentzVector", Mu1)
Cat5.Branch("Mu2", "TLorentzVector", Mu2)
Cat5.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
Cat5.Branch("Mu_Pt", Mu_Pt)
Cat5.Branch("Mu_Eta", Mu_Eta)
Cat5.Branch("Nu_El", Nu_El, "Nu_El/I")
Cat5.Branch("El_Pt", El_Pt)
Cat5.Branch("El_Eta", El_Eta)
Cat5.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
Cat5.Branch("Jet_Pt", Jet_Pt)
Cat5.Branch("Jet_Eta", Jet_Eta)
Cat5.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
Cat5.Branch("genweight", genweight, "genweight/F")

Cat6.Branch("Event_No", Event_No, "Event_No/I")
Cat6.Branch("Dilep", "TLorentzVector", Dilep)
Cat6.Branch("Mu1", "TLorentzVector", Mu1)
Cat6.Branch("Mu2", "TLorentzVector", Mu2)
Cat6.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
Cat6.Branch("Mu_Pt", Mu_Pt)
Cat6.Branch("Mu_Eta", Mu_Eta)
Cat6.Branch("Nu_El", Nu_El, "Nu_El/I")
Cat6.Branch("El_Pt", El_Pt)
Cat6.Branch("El_Eta", El_Eta)
Cat6.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
Cat6.Branch("Jet_Pt", Jet_Pt)
Cat6.Branch("Jet_Eta", Jet_Eta)
Cat6.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
Cat6.Branch("genweight", genweight, "genweight/F")

Cat7.Branch("Event_No", Event_No, "Event_No/I")
Cat7.Branch("Dilep", "TLorentzVector", Dilep)
Cat7.Branch("Mu1", "TLorentzVector", Mu1)
Cat7.Branch("Mu2", "TLorentzVector", Mu2)
Cat7.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
Cat7.Branch("Mu_Pt", Mu_Pt)
Cat7.Branch("Mu_Eta", Mu_Eta)
Cat7.Branch("Nu_El", Nu_El, "Nu_El/I")
Cat7.Branch("El_Pt", El_Pt)
Cat7.Branch("El_Eta", El_Eta)
Cat7.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
Cat7.Branch("Jet_Pt", Jet_Pt)
Cat7.Branch("Jet_Eta", Jet_Eta)
Cat7.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
Cat7.Branch("genweight", genweight, "genweight/F")

Cat8.Branch("Event_No", Event_No, "Event_No/I")
Cat8.Branch("Dilep", "TLorentzVector", Dilep)
Cat8.Branch("Mu1", "TLorentzVector", Mu1)
Cat8.Branch("Mu2", "TLorentzVector", Mu2)
Cat8.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
Cat8.Branch("Mu_Pt", Mu_Pt)
Cat8.Branch("Mu_Eta", Mu_Eta)
Cat8.Branch("Nu_El", Nu_El, "Nu_El/I")
Cat8.Branch("El_Pt", El_Pt)
Cat8.Branch("El_Eta", El_Eta)
Cat8.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
Cat8.Branch("Jet_Pt", Jet_Pt)
Cat8.Branch("Jet_Eta", Jet_Eta)
Cat8.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
Cat8.Branch("genweight", genweight, "genweight/F")

Cat9.Branch("Event_No", Event_No, "Event_No/I")
Cat9.Branch("Dilep", "TLorentzVector", Dilep)
Cat9.Branch("Mu1", "TLorentzVector", Mu1)
Cat9.Branch("Mu2", "TLorentzVector", Mu2)
Cat9.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
Cat9.Branch("Mu_Pt", Mu_Pt)
Cat9.Branch("Mu_Eta", Mu_Eta)
Cat9.Branch("Nu_El", Nu_El, "Nu_El/I")
Cat9.Branch("El_Pt", El_Pt)
Cat9.Branch("El_Eta", El_Eta)
Cat9.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
Cat9.Branch("Jet_Pt", Jet_Pt)
Cat9.Branch("Jet_Eta", Jet_Eta)
Cat9.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
Cat9.Branch("genweight", genweight, "genweight/F")

Cat10.Branch("Event_No", Event_No, "Event_No/I")
Cat10.Branch("Dilep", "TLorentzVector", Dilep)
Cat10.Branch("Mu1", "TLorentzVector", Mu1)
Cat10.Branch("Mu2", "TLorentzVector", Mu2)
Cat10.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
Cat10.Branch("Mu_Pt", Mu_Pt)
Cat10.Branch("Mu_Eta", Mu_Eta)
Cat10.Branch("Nu_El", Nu_El, "Nu_El/I")
Cat10.Branch("El_Pt", El_Pt)
Cat10.Branch("El_Eta", El_Eta)
Cat10.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
Cat10.Branch("Jet_Pt", Jet_Pt)
Cat10.Branch("Jet_Eta", Jet_Eta)
Cat10.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
Cat10.Branch("genweight", genweight, "genweight/F")

def MuScaleFactor (mu_charge, mu_pt, mu_eta, mu_phi, nTrack):
    scaleFactor = 1.0
    u1 = ROOT.gRandom.Rndm()
    u2 = ROOT.gRandom.Rndm()
    if "Run" in FileArg[1]:
        scaleFactor = rocCor.kScaleDT(mu_charge, mu_pt, mu_eta, mu_phi, 0, 0)
    else: 
        if mu_pt == GenLep1.Pt():
            scaleFactor = rocCor.kScaleFromGenMC(mu_charge, mu_pt, mu_eta, mu_phi, nTrack, GenLep1.Pt(), u1, 0, 0);
        if mu_pt == GenLep2.Pt():
            scaleFactor = rocCor.kScaleFromGenMC(mu_charge, mu_pt, mu_eta, mu_phi, nTrack, GenLep2.Pt(), u1, 0, 0);
        else:
            scaleFactor = rocCor.kScaleAndSmearMC(mu_charge, mu_pt, mu_eta, mu_phi, nTrack, u1, u2, 0, 0);   
    
    return scaleFactor 

def MuonSelection (mu_pt , mu_eta, mu_phi, mu_m, mu_iso, mu_charge, mu_id, nTrack):
    m = ROOT.TLorentzVector()
    mu = ROOT.TLorentzVector()
    m.SetPtEtaPhiM(mu_pt, mu_eta, mu_phi, mu_m)
    mu = m * MuScaleFactor(mu_charge, mu_pt, mu_eta, mu_phi, nTrack)

    if mu.Pt() < 20: return False 
    if abs(mu.Eta()) > 2.4: return False 
    if mu_iso > 0.25: return False
    if not mu_id: return False 
    Mu_Pt.push_back(mu.Pt())
    Mu_Eta.push_back(mu.Eta())
    Mu_Charge.push_back(mu_charge)
    Mu_Phi.push_back(mu.Phi())
    Mu_M.push_back(mu.M())
    return True

def ElecSelection (elec_pt, elec_eta, elec_phi, elec_m, elec_iso, elec_id, mu_p):
    if elec_pt < 10: return False   
    if abs(elec_eta) > 2.5: return False 
    if (abs(elec_eta > 1.4442) and (abs(elec_eta) < 1.566)): return False 
    if elec_iso > 0.15: return False 
    if elec_id < 3: return False 
    ElP4 = ROOT.TLorentzVector()
    ElP4.SetPtEtaPhiM(elec_pt, elec_eta, elec_phi, elec_m)
    for i, mu in enumerate (mu_p): 
        if mu.DeltaR(ElP4) < 0.4:
            return False 
    El_Pt.push_back(elec_pt)
    El_Eta.push_back(elec_eta)
    return True

def JetSelection (jet_pt, jet_eta, jet_phi, jet_m, jet_id, jet_b, mu_p):
    JET_LOOSE = (1<<0)
    if jet_pt < 30: return False  
    if abs(jet_eta) > 4.7: return False 
    if (jet_id & JET_LOOSE) == 0: return False 
    JetP4 = ROOT.TLorentzVector()
    JetP4.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_m)
    for i, mu in enumerate (mu_p):
        if mu.DeltaR(JetP4) < 0.4:
            return False         
    Jet_Pt.push_back(jet_pt)
    Jet_Eta.push_back(jet_eta)
    Jet_CSVV2.push_back(jet_b)
    Jet_M.push_back(jet_m)
    Jet_Phi.push_back(jet_phi)
    return True 

def BtaggedSelection (Jet_Pt, Jet_Eta, Jet_CSVV2):
    if Jet_Pt < 20: return False 
    if abs(Jet_Eta) > 2.4: return False 
    if Jet_CSVV2 < 0.848: return False 
    return True 

### Open Root File ###
#filelist = [#"/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016B-18Apr2017_ver2-v1/171112_160845/0000/run2_2016RD_NANO_*",
#            "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016B-18Apr2017_ver2-v1/171112_160845/0001/run2_2016RD_NANO_*",
#            "root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016C-18Apr2017-v1/171112_161105/0000/run2_2016RD_NANO*",
#            "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016D-18Apr2017-v1/171112_161322/0000/run2_2016RD_NANO*",
#            "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016E-18Apr2017-v1/171112_161612/0000/run2_2016RD_NANO_*",
#            "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016F-18Apr2017-v2/171112_161842/0000/run2_2016RD_NANO_*",
#            "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016G-18Apr2017-v1/171112_162118/0000/run2_2016RD_NANO_*",
#            "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016G-18Apr2017-v1/171112_162118/0001/run2_2016RD_NANO_*",
#            "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016H-18Apr2017-v1/171112_162332/0000/run2_2016RD_NANO_*",
#            "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016H-18Apr2017-v1/171112_162332/0001/run2_2016RD_NANO_*"
#            ]
#NanoFiles = []       
#for i, Nfile in enumerate(filelist):            
#    NanoFiles = NanoFiles + glob.glob(Nfile)
#    print NanoFiles

for i,Nfile in enumerate(FileArg[2:]):
#for i,Nfile in enumerate(NanoFiles):
    CurFile = TNetXNGFile(Nfile)
    Tree = CurFile.Get("Events")
    
    for ive, event in enumerate(Tree):
    ### Clear Vectors ### 
        Mu_Pt.clear()
        Mu_Eta.clear()
        Mu_Charge.clear()
        Mu_Phi.clear()
        Mu_M.clear()

        El_Pt.clear()
        El_Eta.clear()

        Jet_Pt.clear()
        Jet_Eta.clear()
        Jet_CSVV2.clear()
        Jet_M.clear()
        Jet_Phi.clear()
        
    ### Start Event Loop ###
                  
      ### Object Selection ########################################################################################################################
        ### Muon Selection ### 
        Event_Total[0] = 1
        Event_Tot.Fill(Event_Total[0])
        NuMu = 0
        NuEl = 0
        NuJet = 0
        NuBJet = 0
        
        ### Weights ###
        if "Run" not in FileArg[1]:
            genweight[0] = event.genWeight 

        ### Generated Lepton ###
            for i in range(event.nGenDressedLepton): 
                if abs(event.GenDressedLepton_pdgId[i]) != 13 : 
                    break    
            
                bosonSample = False
                isfromBoson = False 
                for k in range(event.nGenPart):
                    if (event.GenPart_genPartIdxMother[k] == 23 or event.GenPart_genPartIdxMother[k] == 25):
                        bosonSample = True 
                        isfromBoson = True 
                
                if isfromBoson == True:    
                    if event.GenDressedLepton_pdgId[i] == 13:
                        GenLep1.SetPtEtaPhiM(event.GenDressedLepton_pt[i], event.GenDressedLepton_eta[i], event.GenDressedLepton_phi[i], event.GenDressedLepton_mass[i])
                    else: 
                        GenLep2.SetPtEtaPhiM(event.GenDressedLepton_pt[i], event.GenDressedLepton_eta[i], event.GenDressedLepton_phi[i], event.GenDressedLepton_mass[i])
        
        ### Muon Selection ###        
        if event.nMuon > 0:
            for i in range(event.nMuon):
                if MuonSelection(event.Muon_pt[i], event.Muon_eta[i], event.Muon_phi[i], event.Muon_mass[i], event.Muon_pfRelIso04_all[i], event.Muon_charge[i], event.Muon_mediumId[i], event.Muon_nTrackerLayers[i]):
                    NuMu += 1
                Nu_Mu[0] = NuMu    
        if Nu_Mu < 2:
            continue

        ### Muon TLorentzVector ###    
        Mu_P4 = []            
        for i in range(NuMu):
            MuP4 = ROOT.TLorentzVector()
            MuP4.SetPtEtaPhiM(Mu_Pt[i], Mu_Eta[i], Mu_Phi[i], Mu_M[i])
            Mu_P4.append(MuP4)

        ### Election Selection ###  
        if event.nElectron > 0:
            for k in range(event.nElectron):
                if ElecSelection(event.Electron_pt[k], event.Electron_eta[k], event.Electron_phi[k], event.Electron_mass[k], event.Electron_pfRelIso03_all[k], event.Electron_cutBased[k], Mu_P4):
                    NuEl += 1 
            Nu_El[0] = NuEl

        ### Jet Selection ###
        if event.nJet > 0:
            for j in range(event.nJet):
                if JetSelection(event.Jet_pt[j], event.Jet_eta[j], event.Jet_phi[j], event.Jet_mass[j], event.Jet_jetId[j], event.Jet_btagCSVV2[j], Mu_P4):
                    NuJet += 1
            Nu_Jet[0] = NuJet
            
        ## B-Tagged Jet Selection ###
            for l in xrange(NuJet):
                if BtaggedSelection(Jet_Pt[l], Jet_Eta[l], Jet_CSVV2[l]):
                    NuBJet += 1
                Nu_BJet[0] = NuBJet    

      ### Event Selection ############################################################################################################################
        ### Muon Triggers ###
        if not event.HLT_IsoMu24 and not event.HLT_IsoTkMu24:
            continue   
        
        ### Primary Vertex ###
        if abs(event.PV_z) >= 24:
            continue 
        if event.PV_npvs == 0:
            continue 
        if event.PV_ndof < 4:
            continue 

        ### Muon With Opp Charge ###
        Charge = False
        for i in xrange(NuMu):   
            if ((Mu_Pt[0] > 26) or (Mu_Pt[i] > 26)): 
                if Mu_Charge[0] * Mu_Charge[i] < 0:
                    Mu1.SetPtEtaPhiM(Mu_Pt[0],Mu_Eta[0],Mu_Phi[0],Mu_M[0])
                    Mu2.SetPtEtaPhiM(Mu_Pt[i],Mu_Eta[i],Mu_Phi[i],Mu_M[i])
                    Charge = True
                    break
        if Charge == False: 
            continue
        
        Dilep_ = Mu1 + Mu2 
        Dilep.SetPtEtaPhiM(Dilep_.Pt(),Dilep_.Eta(),Dilep_.Phi(),Dilep_.M())
        if Dilep.M() < 12:
            continue 
      ### CATEGORY ####################################################################################################################################  
        
        ### CAT 1 1-BJet ###
        if NuBJet == 1: 
            ### Category 1: 1 Electron ### 
            if NuEl == 1:
                if NuMu > 2:
                    continue 
                else:
                    Cat1.Fill()
            ### Category 2: 2 Electrons ###    
            if NuEl == 2:
                if NuMu > 2:
                    continue 
                else:
                    Cat2.Fill()

            ### Category 3: 3 Muons ###
            if NuMu == 3:     
                if NuEl > 0:
                    continue
                else:
                    Cat3.Fill()
            ### Category 3: 4 Muons ###
            if NuMu == 4:    
                if NuEl > 0:
                    continue 
                else:
                    Cat4.Fill()
            if NuJet - NuBJet == 4:
                Cat5.Fill()

        ### CAT 2 2-BJet ###
        if NuBJet == 2: 
            ### Category 1: 1 Electron ### 
            if NuEl == 1:
                if NuMu > 2:
                    continue 
                else:
                    Cat6.Fill()
            ### Category 2: 2 Electrons ###    
            if NuEl == 2:
                if NuMu > 2:
                    continue 
                else:
                    Cat7.Fill()

            ### Category 3: 3 Muons ###
            if NuMu == 3:     
                if NuEl > 0:
                    continue
                else:
                    Cat8.Fill()
            ### Category 3: 4 Muons ###
            if NuMu == 4:    
                if NuEl > 0:
                    continue 
                else:
                    Cat9.Fill()
            ### Category 
            if NuJet - NuBJet == 4:
                Cat10.Fill()

        Event_No[0] = 1             
        ALL.Fill()    
      

f.Write()
f.Close()
