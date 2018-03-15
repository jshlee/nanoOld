import ROOT

#fdir = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0000/"
#fdir = "/xrootd/store/group/nanoAOD/run2_2016v3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180125_131546/0000/"
fdir = "/cms/scratch/seulgi/nanoAOD/src/nano/analysis/topmass/Results/results_merged/"

h = ROOT.TH1D("dr", "dr", 100, 0, 3)
h2 = ROOT.TH1D("mass sig", "mass sig", 100, 1.7, 2)
h3 = ROOT.TH1D("mass bkg", "mass bkg", 100, 1.7, 2)

h_genPt = ROOT.TH1D("gen pt", "gen pt", 100, 0, 100)
h_recoPt = ROOT.TH1D("reco pt", "reco pt", 100, 0, 100)
h_genEta = ROOT.TH1D("gen eta", "gen eta", 100, -3, 3)
h_recoEta = ROOT.TH1D("reco eta", "reco eta", 100, -3, 3)

for ifile in range(1, 4):
    f = ROOT.TFile(fdir+"topmass_ch_%d.root"%ifile)
    for iev, e in enumerate(f.event):
        if iev%100 == 0 : print "iev: ", iev
        for i in xrange(e.ncmeson):
            if e.cmeson_pdgId[i] != 421: continue
            if e.cmeson_pt[i] < 5: continue
            if abs(e.cmeson_eta[i]) > 2.4: continue
            #tlv_cm = ROOT.TLorentzVector()
            #tlv_cm.SetPtEtaPhiM(e.cmeson_pt[i], e.cmeson_eta[i], e.cmeson_phi[i], e.cmeson_mass[i])
            h3.Fill(e.cmeson_mass[i])
            matched = False
            dr = 999
            for j in xrange(e.nGenPart):
                if abs(e.GenPart_pdgId[j]) != 421: continue
                if e.GenPart_pt[j] < 5: continue
                if abs(e.GenPart_eta[j]) > 2.4: continue
                tlv_gen  = ROOT.TLorentzVector()
                tlv_gen.SetPtEtaPhiM(e.GenPart_pt[j], e.GenPart_eta[j], e.GenPart_phi[j], 1.86484)
                tlv_reco = ROOT.TLorentzVector()
                tlv_reco.SetPtEtaPhiM(e.cmeson_pt[i], e.cmeson_eta[i], e.cmeson_phi[i], e.cmeson_mass[i])
                if tlv_gen.DeltaR(tlv_reco) < dr:
                    mindrIdx = j
                    dr = tlv_gen.DeltaR(tlv_reco)

            if dr == 999: continue
            h.Fill(dr)
            if dr < 0.1:
                matched = True
                h2.Fill(e.cmeson_mass[i])
                h_genPt.Fill(e.GenPart_pt[mindrIdx])
                h_recoPt.Fill(e.cmeson_pt[i])
                h_genEta.Fill(e.GenPart_eta[mindrIdx])
                h_recoEta.Fill(e.cmeson_eta[i])

c = ROOT.TCanvas()
h.Draw()
c.Print("dr.png")

h3.Draw()
h2.Draw("same")
h2.SetFillColor(2)
h3.SetFillColor(4)
h2.SetLineColor(2)
h3.SetLineColor(4)
h3.SetMinimum(0)
c.Print("m.png")

h_recoPt.Draw()
h_genPt.Draw("same")
h_genPt.SetLineColor(2)
c.Print("pt.png")

h_recoEta.Draw()
h_genEta.Draw("same")
h_genEta.SetLineColor(2)
c.Print("eta.png")
