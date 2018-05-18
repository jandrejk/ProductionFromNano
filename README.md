# ProductionFromNano

Tools for running HEPHY ntuple production, based on WAW ntuple production, from NanoAOD.

Port of [tools by A. Kalinowski, M. Bluj et al. based on KLUB/LLR trees](https://github.com/akalinow/Production.git) to CMS NanoAOD

---

* NanoEventsSkeleton.{h,C}: interface to NanoAOD ntuples produced from NanoAOD with MakeClass tool, it needs be updated after each modification of NanoAOD format.
* HTauTauTreeFromNanoBase.{h,C}: base class to translate to WAW format
* HElTauhTreeFromNano.{h,C}: specialization for the e+tau channel
* HMuTauhTreeFromNano.{h,C}: specialization for the mu+tau channel
* HTauhTauhTreeFromNano.{h,C}: specialization for the di-tau channel
* HTTEvent.{h,cxx}: definition of WAW analysis classes
* PropertyEnum.h, TriggerEnum.h, FilterEnum.h: definition of enums, (re)generated by the tool
* AnalysisEnums.h, SelectionBitsEnum.h: definition of enums
* convertNano.py: script to run conversion - currently not functional, use run.py and convertNanoParallel.py (but contains important lines to run on data)
* run.py: script to run NanoAOD-to-SyncNtuple conversion in parallel, calls convertNanoParallel.py
* convertNanoParallel.py: script to run conversion
* ParameterConfig.cc: defines cut values etc; not fully implemented in the code yet (only parts are used)
* Missing: production tools, need be taken modified from old WAW repo

---

Installation recipe for CMSSW_9_4_4
```
scram project -n CMSSW_9_4_4_fromNano CMSSW CMSSW_9_4_4
cd CMSSW_9_4_4_fromNano/src/
cmsenv
# NanoAOD and tools 
git cms-addpkg PhysicsTools/NanoAOD #not mandatory, but it initializes git for CMSSW which is already useful
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools # not used for now, but can be in future, e.g. JES?
# SVFit
git clone https://github.com/svfit/ClassicSVfit.git TauAnalysis/ClassicSVfit -b release_2018Mar20
git clone https://github.com/svfit/SVfitTF.git TauAnalysis/SVfitTF
# MET recoil corrections
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections
# production tools based on WawTools from NanoAOD
git clone https://github.com/mflechl/ProductionFromNano.git WawTools/NanoAODTools -b SM2017ML
# This is just needed to get rid of warnings
git cms-addpkg FWCore/MessageLogger
cat FWCore/MessageLogger/interface/MessageDrop.h | sed s#CMS_THREAD_SAFE##g > FWCore/MessageLogger/interface/MessageDrop.h2
mv FWCore/MessageLogger/interface/MessageDrop.h2 FWCore/MessageLogger/interface/MessageDrop.h 
# compile
scram b -j 4
```

---
How to run
```
cd $CMSSW_BASE/WawTools/NanoAODTools
./run.py
```
This neeeds to be extended for larger productions, and to work with non-local files.

---
Important ressources

* [https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html](https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html)
* [https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD)
* [https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py](https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py)



---
Release notes:
* 16.01.2018, M.Bluj, initial version for 80X (2016) inputs with CMSSW_9_4_2
* 26.02.2018, M.Bluj, update to NanoAOD of 05Feb2018 production of 2016 data wiht CMSSW_9_4_4
* 18.05.2018, M.Flechl, produce sync ntuple directly; produce e-tau ntuples; several fixes for the sync
