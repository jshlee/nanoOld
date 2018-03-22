import ROOT

#fdir = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0000/"
#fdir = "/xrootd/store/group/nanoAOD/run2_2016v3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180125_131546/0000/"

h = ROOT.TH1D("dr", "dr", 100, 0, 3)
h2 = ROOT.TH1D("mass sig", "mass sig", 100, 1.7, 2.0)
h3 = ROOT.TH1D("mass bkg", "mass bkg", 100, 1.7, 2.0)

h_genPt = ROOT.TH1D("gen pt", "gen pt", 100, 0, 100)
h_recoPt = ROOT.TH1D("reco pt", "reco pt", 100, 0, 100)
h_genEta = ROOT.TH1D("gen eta", "gen eta", 100, -3, 3)
h_recoEta = ROOT.TH1D("reco eta", "reco eta", 100, -3, 3)


events = ROOT.TChain("Events")

#events.Add("/xrootd/store/group/nanoAOD/run2_2016v3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180125_131546/0000/nanoAOD_9*.root")
events.Add("/xrootd/store/user/jlee/tsW_13TeV_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v4_hadAOD/nanoAOD_*.root")


#for ifile in range(1, 2):
#    f = ROOT.TFile(fdir+"nanoAOD_%d.root"%ifile)
#    print ifile
for iev, e in enumerate(events):
    if iev%100 == 0 : print "iev: ", iev
    #print "iev: ", iev
    for i in xrange(e.nhad):
        if e.had_pdgId[i] != 421: continue
        if e.had_pt[i] < 5: continue
        if abs(e.had_eta[i]) > 2.4: continue
        h3.Fill(e.had_mass[i])

    for i in xrange(e.nGenPart):
        if abs(e.GenPart_pdgId[i]) != 421: continue
        if e.GenPart_pt[i] < 5: continue
        if abs(e.GenPart_eta[i]) > 2.4: continue
        if abs(e.GenPart_status[i]) == 1: continue
        tlv_gen  = ROOT.TLorentzVector()
        tlv_gen.SetPtEtaPhiM(e.GenPart_pt[i], e.GenPart_eta[i], e.GenPart_phi[i], 1.86484)
        dr = 999
        mindrIdx = 999
        gim = e.GenPart_genPartIdxMother[i]
        #print e.ncmeson

        for m in xrange(e.nGenPart):
            #print "m, gim: ", m,gim
            if gim == -1: break
            if abs(e.GenPart_pdgId[gim]) == 6:
                for j in xrange(e.nhad):
                    #print "j: ", j
                    if e.had_pdgId[j] != 421: continue
                    if e.had_pt[j] < 5: continue
                    if abs(e.had_eta[j]) > 2.4: continue
                    
                    tlv_cme = ROOT.TLorentzVector()
                    tlv_cme.SetPtEtaPhiM(e.had_pt[j], e.had_eta[j], e.had_phi[j], e.had_mass[j])
                    if tlv_gen.DeltaR(tlv_cme) < dr:
                        mindrIdx = j
                        #print "mind: ",mindrIdx
                        #print "j: ", j
                        dr = tlv_gen.DeltaR(tlv_cme)
                    #print "dr: mind:" , dr, mindrIdx
                #print "dr2: , mind2: ", dr, mindrIdx
                if dr == 999: continue
                h.Fill(dr)
                if dr < 1.0 and mindrIdx!=999:
                    #print "min: " ,mindrIdx
                    h2.Fill(e.had_mass[mindrIdx])
                    
                    h_recoPt.Fill(e.had_pt[mindrIdx])
                    h_recoEta.Fill(e.had_eta[mindrIdx])
                h_genPt.Fill(e.GenPart_pt[i])
                h_genEta.Fill(e.GenPart_eta[i])
                break
            if abs(e.GenPart_pdgId[gim]) != 6:
                ppd = e.GenPart_genPartIdxMother[gim]
                gim = ppd
    
c = ROOT.TCanvas()
h.Draw()
c.Print("dr.png")

#hMP.Draw()
#c.Print("genmotherId.png")

h3.Draw()
h2.Draw("same")
h2.SetFillColor(2)
h3.SetFillColor(4)
h2.SetLineColor(2)
h3.SetLineColor(4)
h3.SetMinimum(0)
c.Print("m.png")

h_genPt.Draw()
h_recoPt.Draw("same")
h_genPt.SetLineColor(2)
c.Print("pt.png")

h_genEta.Draw()
h_recoEta.Draw("same")
h_genEta.SetLineColor(2)
c.Print("eta.png")
