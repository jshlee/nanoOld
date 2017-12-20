import os
import pprint
pp = pprint.pprint
filename = 'SingleMuon2_Run2016H.txt'
f = open(filename , "w")
a = []
pathlist=  [#"/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016B-18Apr2017_ver2-v1/171112_160845/0000/",
            #"/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016B-18Apr2017_ver2-v1/171112_160845/0001/",
      #      "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016C-18Apr2017-v1/171112_161105/0000/",
      #      "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016D-18Apr2017-v1/171112_161322/0000/",
      #      "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016E-18Apr2017-v1/171112_161612/0000/",
      #      "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016F-18Apr2017-v2/171112_161842/0000/",
      #      "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016G-18Apr2017-v1/171112_162118/0000/",
      #      "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016G-18Apr2017-v1/171112_162118/0001/",
            "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016H-18Apr2017-v1/171112_162332/0000/",
            "/xrootd/store/group/nanoAOD/SingleMuon/run2_2016RD_NANO_Run2016H-18Apr2017-v1/171112_162332/0001/"
            ]   

for path in pathlist:
    temp = os.listdir(path)
    temp = map(lambda p: os.path.join(path, p), temp) 
    a +=  temp

for i in a:
    newstr = str.replace(i, '/xrootd','root://cms-xrdr.sdfarm.kr:1094///xrd')
    f.write(newstr + "\n")

pp(newstr)
f.close()    
