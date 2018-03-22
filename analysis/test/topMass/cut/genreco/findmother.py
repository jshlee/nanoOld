import ROOT

#fdir = "/xrootd/store/user/jlee/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180307_072216/0000/"
fdir = "/xrootd/store/group/nanoAOD/run2_2016v3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180125_131546/0000/"

h = ROOT.TH1D("genmother", "genmother", 100, 0, 10)

for ifile in range(1, 2):
    f = ROOT.TFile(fdir+"nanoAOD_%d.root"%ifile)
    print ifile
    for iev, e in enumerate(f.Events):
        if iev%100 == 0 : print "iev: ", iev

        for i in xrange(e.nGenPart):
            print "real i: ", i
            if abs(e.GenPart_pdgId[i]) != 421: continue
            if e.GenPart_pt[i] < 5: continue
            if abs(e.GenPart_eta[i]) > 2.4: continue
            if abs(e.GenPart_status[i]) == 1: continue
            #if e.GenPart_genPartIdxMother[i] == -1: continue
            gim = e.GenPart_genPartIdxMother[i]
            print "gim_i: ",gim

            for m in xrange(e.nGenPart):
                print "gim: ", gim
                if gim == -1: break
                if abs(e.GenPart_pdgId[gim]) == 6:
                    h.Fill(e.GenPart_pdgId[gim])
                    print "last pdgId: ", e.GenPart_pdgId[gim]
                    break
                if abs(e.GenPart_pdgId[gim]) != 6:
                    print "not top -- ", e.GenPart_pdgId[gim]
                    ppd = e.GenPart_genPartIdxMother[gim]
                    print ppd
                    gim = ppd
            #print e.GenPart_pdgId[pdg]
            #h.Fill(e.GenPart_pdgId[pdg]) 
            #h.Fill(e.GenPart_genPartIdxMother[i])

c = ROOT.TCanvas()
h.Draw()
c.Print("genmother.png")
