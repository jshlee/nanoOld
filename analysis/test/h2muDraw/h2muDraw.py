#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys
from CATTools.CatAnalyzer.histoHelper import *
from CATTools.CatAnalyzer.histoHelper_h2mu import *
from ROOT import TLorentzVector
#import DYestimation
ROOT.gROOT.SetBatch(True)

datalumi = 36814
version = os.environ['CMSSW_VERSION']
CMS_lumi.lumi_sqrtS = "%.2f fb^{-1}, #sqrt{s} = 13 TeV 25ns "%(float(datalumi)/1000)

minirootdir = "/xrootd/store/user/daniellee/cattreeMiniRoc/CMSSW_8_0_26_patch1/results_merged/h2muAnalyzer_"
#minirootdir = "/cms/scratch/daniel/CMSSW_8_0_26_patch1/src/CATTools/CatAnalyzer/test/Nano_AOD/ROOT/"
nanorootdir = "/cms/scratch/daniel/CMSSW_8_0_26_patch1/src/CATTools/CatAnalyzer/test/Nano_AOD/ROOT/"
#nanorootdir = "/xrootd/store/user/daniellee/cattreeMiniC/CMSSW_8_0_26_patch1/results_merged/h2muAnalyzer_"

minifilelist = [
                'DYJets',
               # 'GG_HToMuMu',
               # 'NanoGluGlu_HToMuMu',
               # 'NanoDYJets_noRoC',
               # 'SingleMuon2_Run2016B',
               ]


nanofilelist = [
               # 'GG_HToMuMu',
               #'NanoGluGlu_HToMuMu',
               'NanoDYJets_ROC',
               # 'NanoSingleB',
               ]

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

minicut = 'dilep.M()>60&&step>4'
#minicut = 'Dilep.M()>60'
nanocut = 'Dilep.M()>60'
#nanocut = 'dilep.M()>60&&step>4'

miniplotvar = "dilep.M()"
#miniplotvar = "Dilep.M()"
nanoplotvar = "Dilep.M()"
#nanoplotvar = "lep1.Pt()"
binning = [150,50,200]
x_name = 'dilep_mass'
y_name = 'Events'
dolog = True
f_name = 'dilep_mass_C'
print "File Name : %s" %f_name
minitname = "cattree/nom"
#minitname = "nEvent"
nanotname = "nEvent"
#nanotname = "cattree/nom"

minirdfname = minirootdir + minifilelist[0] +".root"
nanordfname = nanorootdir + nanofilelist[0] +".root"
print "minirdname: %s\n tname: %s\n binning: %s\n plotvar: %s\n cut: %s\n"%(minirdfname, minitname, binning, miniplotvar, minicut)
print "nanordname: %s\n tname: %s\n binning: %s\n plotvar: %s\n cut: %s\n"%(nanordfname, nanotname, binning, nanoplotvar, nanocut)

minird = makeTH1(minirdfname, minitname, 'CatTuple', binning, miniplotvar, minicut)
nanord = makeTH1(nanordfname, nanotname, 'NanoAOD', binning, nanoplotvar, nanocut)


def drawRatio(name, cmsLumi, datamini, datanano, x_name, y_name, doLog=False, doRatio=True, ratioRange=0.45, legx=0.68, legfontsize=0.030):
    leg = ROOT.TLegend(legx,0.68,legx+0.2,0.91)
    leg.SetBorderSize(0)
    #leg.SetNColumns(2)
    leg.SetTextSize(legfontsize)
    leg.SetTextFont(42)
    leg.SetLineColor(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.AddEntry(datamini,"CatTuple","lp")
    leg.AddEntry(datanano,"NanoAOD","lp")
    hratio = datanano.Clone("hratio")
    hratio.Reset()

    hratio.Divide(datamini,datanano,1.,1.,"B")

    tdrstyle.setTDRStyle()

    setDefTH1Style(datamini, x_name, y_name)
    datamini.SetName('data')
    datamini.SetMaximum(datamini.GetMaximum()*1.8)
    if doLog:
        datamini.SetMaximum(datamini.GetMaximum()*100)
        #data.SetMinimum(10**-3)
    else:
        datamini.GetYaxis().SetTitleSize(0.04)
        datamini.GetYaxis().SetLabelSize(0.024)
        datamini.GetYaxis().SetTitleOffset(1.35)
    
    ratio_fraction = 0 
    if doRatio:
        ratio_fraction = 0.3
        datamini.GetXaxis().SetLabelSize(0)
        datamini.GetXaxis().SetTitleSize(0)
        setDefTH1Style(hratio, x_name, "Cat/Nano")
        hratio.GetYaxis().CenterTitle()
        hratio.GetYaxis().SetNdivisions(5)

    canv = makeCanvas(name, doRatio)
    pads=[canv]
    pads = rootplotcore.divide_canvas(canv, ratio_fraction)

    pads[0].cd()
    setMargins(pads[0], doRatio)
    if doLog:
        pads[0].SetLogy()
    
    datamini.SetMarkerColor(4)
    datanano.SetMarkerColor(2)
    datamini.Draw()
    datanano.Draw("same")
    datamini.Draw("esamex0")
    leg.Draw("same")

    pads[0].Update()

    if doRatio:
        pads[1].cd()
        pads[1].SetGridy()
        setMargins(pads[1], doRatio)
        hratio.SetLineColor(1)
        hratio.Draw("e")
        hratio.SetMaximum(1.+ratioRange)
        hratio.SetMinimum(1.-ratioRange)

    for p in pads:
        p.RedrawAxis()
        p.Modified()
        p.Update()

    canv.cd()


#iPos = 0 # in frame
    iPos = 11 # out frame
    if( iPos==0 ):
        cmsLumi.relPosX = 0.1
    cmsLumi.CMS_lumi(pads[0], 0, iPos)

    canv.Modified()
    canv.Update()
    return copy.deepcopy(canv)

canv = drawRatio(f_name, CMS_lumi, minird, nanord , x_name, y_name, dolog)
canv.SaveAs(f_name+".png")

