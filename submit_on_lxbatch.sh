export X509_USER_PROXY=${rundir}/proxy/x509_proxy

echo "---------------------"
echo "Grid certificate 1"
voms-proxy-info --all
echo "---------------------"

echo "Good night..."
date
sleep ${sleeping}
echo "Good morning..."

echo "---------------------"
echo "Current dir: `pwd`"
ls -l
echo "---------------------"

scram project -n CMSSW_9_4_4_fromNano CMSSW CMSSW_9_4_4
cd CMSSW_9_4_4_fromNano/src/
eval `scramv1 runtime -sh`
# NanoAOD and tools 
git cms-addpkg PhysicsTools/NanoAOD
git cms-addpkg FWCore/MessageLogger
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
git clone https://github.com/svfit/ClassicSVfit.git TauAnalysis/ClassicSVfit -b release_2018Mar20
git clone https://github.com/svfit/SVfitTF.git TauAnalysis/SVfitTF
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections
git clone https://github.com/MarkusSpanring/ProductionFromNano.git WawTools/NanoAODTools -b SM2017ML
cat FWCore/MessageLogger/interface/MessageDrop.h | sed s#CMS_THREAD_SAFE##g > FWCore/MessageLogger/interface/MessageDrop.h2
mv FWCore/MessageLogger/interface/MessageDrop.h2 FWCore/MessageLogger/interface/MessageDrop.h 
# compile
scram b -j 4
cd WawTools/NanoAODTools
cp ${rundir}/configBall.json .

echo "---------------------"
echo "Current dir: `pwd`"
ls -l
echo "---------------------"

./convertNanoParallel.py

echo "---------------------"
echo "Current dir: `pwd`"
ls -l
echo "---------------------"
chmod 777 ${channel}*root
mv -f ${channel}*root ${outdir}
