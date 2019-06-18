import subprocess as sp
import shlex
import json
import shutil
import os 
from runUtils import checkProxy
import argparse
import ROOT as R

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', dest='location', help='Get samples from local or DAS', choices=["das","private"], default = "das")
    parser.add_argument('--tagMap', dest='update', help='update tagMapping.json', action = "store_true")
    args = parser.parse_args()

    if args.location == "private": 
        # VBFHToTauTau_M125_13TeV_powheg_pythia8_PUMoriond17_94X_mcRun2_asymptotic_v3_v2
        getfromPrivate("/dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190618/",update=args.update)
        # getfromPrivate("/dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190424/",update=args.update)
        #getfromPrivate("/dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190318/")
        #getfromPrivate("/dpm/oeaw.ac.at/home/cms/store/user/jaandrej/")
    else:
        getfromDAS()

def getNumberOfEvents (files) :
    nevents = 0
    for f in files:
        ch = R.TChain("Events")
        ch.Add(f)
        nevents += ch.GetEntries()
    return nevents

def getfromPrivate(location,update):

    if not os.path.exists( "samples" ):
        os.mkdir( "samples" )

    for tier in ["samples/data","samples/mc"]:
        if not os.path.exists( tier ):
            os.mkdir(tier)
        folders = getPrivateFiles(location + tier )
        k = folders.keys()
        k.sort()
        for folder in k:
            if folders[folder]:
                # sub = folder.replace(location+tier,"").split("/")[1]
                # if not os.path.exists( "/".join([tier, sub]) ):
                #     os.mkdir( "/".join([ tier, sub]) )
                with open("tagMapping.json","r+") as FSO:
                        tagmap = json.load(FSO)

                isMC = "_Run2018" not in folder 
                print isMC
                if update and isMC:
                    print folder
                    print "No of events: "
                    Nevents = getNumberOfEvents(files=folders[folder])
                    print Nevents
                    
                    
                    with open("tagMapping.json","r+") as FSO:
                        tagmap = json.load(FSO)
                    tagmap[folder.split('/')[-1]]["nevents_from_DAS"] = Nevents
                    tagmap[folder.split('/')[-1]]["nevents"] = Nevents
                    
                    with open("tagMapping.json", "w") as jsonFile:
                        json.dump(tagmap, jsonFile,indent=4,sort_keys=True)
                
                if isMC :
                    if (tagmap[folder.split('/')[-1]]["fromDAS"] == 'False') :
                        sub = folder.replace(location+tier,"").split("/")[1]
                        if not os.path.exists( "/".join([tier, sub]) ):
                            os.mkdir( "/".join([ tier, sub]) )

                        with open( "/".join([tier,sub, folder.replace(location+tier,"").split("/")[2]]) + ".txt","w") as FSO:
                            FSO.write( "\n".join( folders[folder] ) )
                else :
                    sub = folder.replace(location+tier,"").split("/")[1]
                    if not os.path.exists( "/".join([tier, sub]) ):
                        os.mkdir( "/".join([ tier, sub]) )

                    with open( "/".join([tier,sub, folder.replace(location+tier,"").split("/")[2]]) + ".txt","w") as FSO:
                        FSO.write( "\n".join( folders[folder] ) )

def getPrivateFiles(source):

    proc = sp.Popen(shlex.split("dpns-ls -R {0}".format(source) ), stdout=sp.PIPE, stderr=sp.PIPE, shell=False)
    (out, err) = proc.communicate()
    folder = {}
    current = ""
    for line in out.splitlines():
        if ":" in line:
            folder[line.replace(":","")] = []
            current = line.replace(":","")
        if current and ".root" in line:
            folder[current].append("/".join(["root://hephyse.oeaw.ac.at",current, line]) )

    return folder

def getfromDAS():

    if not checkProxy(): sys.exit()

    if os.path.exists("samples"):
        shutil.rmtree("samples")

    with open("sample_collection.json","r") as FSO:
        config = json.load(FSO)

    for frmt in config:
        f = config[frmt]

        for run in f["run"]:
            r = f["run"][run]

            for sample in r["samples"]:
                links = r["samples"][sample]

                for link in links:
                    folder = "/".join(["samples",frmt, sample])
                    filename = buildFileName( link, run, r["creation"] )
                    content = getDASquery(link)

                    content = content.replace("/store/","root://cms-xrd-global.cern.ch//store/")


                    writeFile(content, filename, folder)
                    # print getDASquery(link,"file")


def getDASquery(link, info="file"):

    proc = sp.Popen(shlex.split('dasgoclient --query="{0} dataset={1}" '.format(info,link) ), stdout=sp.PIPE, stderr=sp.PIPE, shell=False)
    (out, err) = proc.communicate()

    return out

def writeFile( content, filename, folder ):
    if not os.path.exists(folder):
        os.makedirs(folder)

    with open("/".join([folder,filename]), "w" ) as FSO:
        FSO.write(content)

def buildFileName(link, run, creation):

    parts = link.split("/")

    ext = ""
    if "_ext" in parts[2]:
        ext = "_ext" + parts[2].split("_ext")[1].split("-")[0]

    return "_".join([parts[1].replace("-","_"),run, creation]) + ext + ".txt"



if __name__ == '__main__':
    main()

# /dpm/oeaw.ac.at/home/cms/store/mc/RunIIFall17NanoAOD/VBFHToTauTau_M125_13TeV_powheg_pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/
# /dpm/oeaw.ac.at/home/cms/store/data/Run2017B/SingleMuon/NANOAOD/31Mar2018-v1
