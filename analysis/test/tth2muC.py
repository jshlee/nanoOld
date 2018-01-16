import ROOT, os, getopt, sys, array, math, glob, json
from ROOT import * 
from array import array


### Make TTREE ### 
FileArg = sys.argv
print FileArg
tempdir = FileArg[1]
Dirname = "/cms/scratch/yckang/nanoAOD/src/nano/analysis/test/Results/Nano_C_Test/%s/"%tempdir
if not os.path.isdir(Dirname):
    os.makedirs(Dirname)

temp = FileArg[2].split('/').pop()
cattree = Dirname+temp

GJsonF = open("/cms/scratch/daniel/nanoAOD/src/nano/analysis/data/GoldenJson.txt")
Gfile = json.load(GJsonF)
GJsonF.seek(0)
json_hold = GJsonF.read()

lumiMap = ROOT.std.map('unsigned int', 'std::vector<std::array<unsigned int, 2>>')()

for i in Gfile :
    lumiVector = ROOT.std.vector('std::array<unsigned int, 2>')()
    for j in range(len(Gfile[i])) :
        lumiArray = ROOT.std.array('unsigned int', 2)()
        lumiArray[0] = int(Gfile[i][j][0])
        lumiArray[1] = int(Gfile[i][j][1])
        lumiVector.push_back(lumiArray)
    lumiMap[int(i)] = lumiVector

ROOT.gROOT.LoadMacro("/cms/scratch/yckang/nanoAOD/src/nano/analysis/plugins/tth2mu.cc+")
analyzer = ROOT.TTH2MuAnalyzer(cattree);

for i,Nfile in enumerate(FileArg[2:]):
    if "Run" not in FileArg[1]:
        analyzer.AnalyzeForMC(Nfile)
    else:
        analyzer.AnalyzeForData(Nfile)

analyzer.EndOfAnalyze()