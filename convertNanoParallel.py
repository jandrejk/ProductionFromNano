#!/usr/bin/env python

import os
import sys
import threading

from ROOT import gSystem, TChain, TSystem, TFile, TString, vector, TFileCollection, edm

# from PSet import process
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Config as cms
import json
from glob import glob

def getLumisToRun(JSON):
    if JSON == "": return vector( 'edm::LuminosityBlockRange' )()

    vlumis = vector( 'edm::LuminosityBlockRange' )()
    myList = LumiList.LumiList (filename = "/".join(["utils/cert_data",JSONfile]) ).getCMSSWString().split(',')
    lumisToProcess = cms.untracked.VLuminosityBlockRange( myList )


    for BlockRange in lumisToProcess:
        Block = BlockRange.split('-')

        startRun =       int(Block[0].split(':')[0])
        startLumiBlock = int(Block[0].split(':')[1])

        if len(Block) > 1:
            endRun =       int(Block[1].split(':')[0])
            endLumiBlock = int(Block[1].split(':')[1])
        else:
            endRun = startRun
            endLumiBlock = endLumiBlock

        vlumis.push_back( edm.LuminosityBlockRange( edm.LuminosityBlockID(startRun, startLumiBlock),
                                                    edm.LuminosityBlockID(endRun, endLumiBlock) ) )
    return vlumis
#######################################################################################################
def CopyFile(file_path,system,check_event) :
    aROOTFile = TFile.Open(file_path)
    aTree = aROOTFile.Get("Events")
    remoteEvts = aTree.GetEntries()   


    print "Remote file {file} has {n} number of events".format(file=file_path,n=remoteEvts)
    
    aROOTFile.Close()

    if system == "lxbatch" or system == "condor":
        os.system("xrdcp {0} {1}".format(file_path, file_path.split("/")[-1] ) )
        local_filename = file_path.split("/")[-1]

    aROOTFile = TFile.Open(local_filename)
    aTree = aROOTFile.Get("Events")

    if check_event > 0:
        aTree = aTree.CopyTree("event == {0}".format(check_event) )

    entries = aTree.GetEntries()

    if not entries:
        print "file is empty. Aborting"
        exit(1)
    if not remoteEvts == entries and not check_event > 0:
        print "File was not copied properly. Aborting"
        exit(2)
    else:
        print "Copying successful! TTree entries: ", entries

    aROOTFile.Close()
    return local_filename


#######################################################################################################

with open("configBall.json","r") as FSO:
    configBall = json.load(FSO)

aFile =            configBall["file"]
channel =          str(configBall["channel"])
systShift =        str(configBall["systShift"])
JSONfile =         str(configBall["certJson"])
nevents =          int(configBall["nevents"])
check_event =      int(configBall["check_event"])

print "-"*30
for i,k in configBall.items():
    if i == "file": continue
    print i+" "*(15 - len(i)),k
print "-"*30

if not "root://" in aFile: aFile = "file://" + aFile
print "Job file name is: {0}".format(aFile)


if "files" in configBall :
    nameOfFilesToBeMerged = []
    print "several files need to to be copied and merged before compiling"
    for name in configBall["files"] :
        print name
        nameOfFilesToBeMerged.append(CopyFile(file_path=name,system=str( configBall["system"] ),check_event=check_event))
    print "copying of all files sucessfull. Merging with the command:"
    hadd_cmd = "python haddnano.py {output} {input}".format(output=aFile,input=" ".join(nameOfFilesToBeMerged))
    print hadd_cmd
    os.system(hadd_cmd)
    aROOTFile = TFile.Open(aFile)
    aTree = aROOTFile.Get("Events")
    
else :
    fname = CopyFile(file_path=aFile, system=str( configBall["system"] ),check_event=check_event)
    aROOTFile = TFile.Open(fname)
    aTree = aROOTFile.Get("Events")


print "Compiling...."


# Some system have problem runnig compilation (missing glibc-static library?).
# First we try to compile, and only then we start time consuming cmssw

gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libZZMatrixElementMELA.so");

assert gSystem.CompileMacro('HTTEvent.cxx','k')
assert gSystem.CompileMacro('utils/TauTriggerSFs2017/src/TauTriggerSFs2017.cc','k')
assert gSystem.CompileMacro('HTXSClassification.cc','k')
assert gSystem.CompileMacro('EventWriter.C','k')
#status *= gSystem.CompileMacro('NanoEventsSkeleton.C') #RECOMPILE IF IT CHANGES!
assert gSystem.CompileMacro('NanoEventsSkeleton.C','k')

gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libTauAnalysisClassicSVfit.so')
gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libTauAnalysisSVfitTF.so')
gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libHTT-utilitiesRecoilCorrections.so')

assert gSystem.CompileMacro('HTauTauTreeFromNanoBase.C','k')
from ROOT import HTTParticle, HTTAnalysis

if channel=='mt':
    assert gSystem.CompileMacro('HMuTauhTreeFromNano.C','k')
    from ROOT import HMuTauhTreeFromNano as Ntuplizer

if channel=='et':
    assert gSystem.CompileMacro('HElTauhTreeFromNano.C','k')
    from ROOT import HElTauhTreeFromNano as Ntuplizer

if channel=='tt':
    assert gSystem.CompileMacro('HTauhTauhTreeFromNano.C','k')
    from ROOT import HTauhTauhTreeFromNano as Ntuplizer



vlumis = getLumisToRun(JSONfile)

HTTParticle.corrType = getattr(HTTAnalysis, systShift )

prefix = "-".join([channel, systShift])
Ntuplizer(  aTree, vlumis, prefix).Loop(nevents,check_event)


# for f in glob('*'):
#     if not f in glob(prefix + '*.root') and not f in glob("log*.txt"):
#         os.remove(f)


exit(0)
