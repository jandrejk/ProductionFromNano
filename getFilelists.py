import subprocess as sp
import shlex
import json
import shutil
import os 
from runUtils import checkProxy
import argparse
import ROOT as R
import re
from xsec import xsec, extensions
from terminal_colors import bcolors
extension_matcher = {}
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', dest='location', help='Get samples from local or DAS', choices=["das","private"], default = "das")
    parser.add_argument('-dbs', dest='dbs_instance', help='where are the  samples - storage site', choices=["global","phys01","phys02","phys03"],default="global")
    parser.add_argument('--tagMap', dest='update', help='update tagMapping.json', action = "store_true")
    parser.add_argument('-coll', dest='collection', help='which sample collection to use', choices=["sample_collection","sample_collection_nanoAODv5"], default = "sample_collection")
    args = parser.parse_args()
    
    if args.location == "private": 
        # VBFHToTauTau_M125_13TeV_powheg_pythia8_PUMoriond17_94X_mcRun2_asymptotic_v3_v2
        getfromPrivate("/dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190618/",update=args.update)
        # getfromPrivate("/dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190424/",update=args.update)
        #getfromPrivate("/dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190318/")
        #getfromPrivate("/dpm/oeaw.ac.at/home/cms/store/user/jaandrej/")
    else:
        getfromDAS(coll=args.collection,USERstorage=args.dbs_instance)

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

def getfromDAS(coll,USERstorage):

    if not checkProxy(): sys.exit()
    
    if os.path.exists("samples"):
        shutil.rmtree("samples")

    with open("{0}.json".format(coll),"r") as FSO:
        config = json.load(FSO)
        
    for frmt in config:
        f = config[frmt]

        for run in f["run"]:
            r = f["run"][run]

            with open("tagMapping.json","r+") as FSO:
                tagmap = json.load(FSO)
                    
            for sample in r["samples"]:
                links = r["samples"][sample]
                folder = "/".join(["samples",frmt, sample])
                
                for link in links:
                    filename = buildFileName( link, run, r["creation"] )
                    content = getDASquery(link,USERstorage=USERstorage)
                    content = content.replace("/store/","root://cms-xrd-global.cern.ch//store/")
                    numberOfEventsInEachFile = getDASquery(link=link, USERstorage=USERstorage,info='file',grep='file.nevents')
                    l1 = content.split("  \n")[:-1]
                    l2 = numberOfEventsInEachFile.split("  \n")[:-1]
                    
                    output_string = ""
                    for i,l in enumerate(l1) :
                        output_string += "{name}, {nevents} \n".format(name=l,nevents=l2[i])

                    print filename
                    writeFile(output_string.rstrip(), filename, folder) # strip away "\n" at the end of output_string


                    if "mc" in folder : # it is a MC sample and we are going to update the tagMapping.json
                        Nevents = re.search(r'"nevents":(\d+)', string=getDASquery(link=link,USERstorage=USERstorage,info="summary",grep="summary")).group(1)
                        Nevents = int(Nevents)
                        file_key = filename.split(".txt")[0] # strip away .txt enxtension
                        
                        if file_key not in tagmap.keys() :
                            print bcolors.OKBLUE + "creating entry in dict" + bcolors.ENDC
                            tagmap[file_key] = {}

                        tagmap[file_key]["fromDAS"] = True
                        tagmap[file_key]["nevents"] = Nevents
                        tagmap[file_key]["nevents_from_DAS"] = Nevents
                        tagmap[file_key]["putag"] = "pileup"
                        tagmap[file_key]["xsec"] = getXSec(name=filename)
                        CheckExtension(file_key=file_key, N=Nevents)
                        
                        
            for k in extension_matcher.keys() :
                tagmap[k]["nevents"] = extensions[extension_matcher[k]]
                                
            with open("tagMapping.json", "w") as jsonFile:
                json.dump(tagmap, jsonFile,indent=4,sort_keys=True)

def CheckExtension(file_key,N) :
    for key in extensions.keys() :
        if re.findall(pattern= '{0}'.format(key), string=file_key) :
            print "Found extension for {0}".format(file_key)
            extension_matcher[file_key] = key
            extensions[key] += N


def getXSec(name) :
    for key in xsec.keys() :
        expression = '{0}'
        if re.findall(pattern= expression.format(key), string=name) :
            return xsec[key]
    print bcolors.FAIL + "Failure: \n {name} did not match any x-section key in xsec.py and a value of -999 was assigned".format(name=name) + bcolors.ENDC
    return -999



def getDASquery(link, USERstorage, info="file", grep='file.name'):
    proc = sp.Popen(shlex.split('dasgoclient --query="{0} dataset={1} instance=prod/{2} | grep {3}" '.format(info,link,USERstorage,grep) ), stdout=sp.PIPE, stderr=sp.PIPE, shell=False)
    
    (out, err) = proc.communicate()

    return out

def writeFile( content, filename, folder ):
    if not os.path.exists(folder):
        os.makedirs(folder)

    with open("/".join([folder,filename]), "w" ) as FSO:
        FSO.write(content)

def buildFileName(link, run, creation):

    parts = link.split("/")
    expr = "(_v[0-9]-)|(_v[0-9]_)"
   
    ext = re.search(pattern=expr, string=link).group(1)[:-1]
    ext = ext.replace("-","_")
    if "_ext" in parts[2]:
        ext += "_ext" + parts[2].split("_ext")[1].split("-")[0]
    elif "_ext1" in parts[2]:
        ext += "_ext1" + parts[2].split("_ext1")[1].split("-")[0]
    elif "_ext2" in parts[2]:
        ext += "_ext2" + parts[2].split("_ext2")[1].split("-")[0]

    return "_".join([parts[1].replace("-","_"),run, creation]) + ext + ".txt"



if __name__ == '__main__':
    main()

# /dpm/oeaw.ac.at/home/cms/store/mc/RunIIFall17NanoAOD/VBFHToTauTau_M125_13TeV_powheg_pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/
# /dpm/oeaw.ac.at/home/cms/store/data/Run2017B/SingleMuon/NANOAOD/31Mar2018-v1
