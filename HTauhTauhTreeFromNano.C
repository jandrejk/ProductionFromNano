#define HTauhTauhTreeFromNano_cxx
#include "HTauhTauhTreeFromNano.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>

/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauhTauhTreeFromNano::pairSelection(unsigned int iPair){
  std::cout<< "bool HTauhTauhTreeFromNano::pairSelection(unsigned int iPair)"<<std::endl;
  ///Baseline+post sync selection as on
  ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2016#Baseline_tau_h_tau_h
  ///Indices for multiplexed ID variables taken from   LLRHiggsTauTau/NtupleProducer/plugins/
  ///HTauTauNtuplizer.cc, MuFiller.cc, TauFiller.cc, EleFiller.cc

  if(httPairs_.empty()) return false;
  std::cout<<1<<std::endl;

  int pdgIdLeg1 = httPairs_[iPair].getLeg1().getPDGid();
  std::cout<<2<<std::endl;
  int pdgIdLeg2 = httPairs_[iPair].getLeg2().getPDGid();
  std::cout<<3<<std::endl;
  
  if( std::abs(pdgIdLeg1)!=15 || std::abs(pdgIdLeg2)!=15 ) return 0;
  std::cout<<4<<std::endl;
  
  unsigned int indexLeg1 = httPairs_[iPair].getIndexLeg1();
  std::cout<<5<<std::endl;
  unsigned int indexLeg2 = httPairs_[iPair].getIndexLeg2();
  std::cout<<6<<std::endl;
  
  //MB sort taus within the pair
  double pt_1 = httLeptonCollection[indexLeg1].getP4().Pt();
  std::cout<<7<<std::endl;
  double pt_2 = httLeptonCollection[indexLeg2].getP4().Pt();
  std::cout<<8<<std::endl;
  if(pt_2>pt_1){//tau with higher-Pt first
    unsigned int indexLegTmp = indexLeg1;
    std::cout<<9<<std::endl;
    indexLeg1 = indexLeg2;
    std::cout<<10<<std::endl;
    indexLeg2 = indexLegTmp;
    std::cout<<11<<std::endl;
  }

  TLorentzVector tau1P4 = httLeptonCollection[indexLeg1].getP4();
  std::cout<<12<<std::endl;
  TLorentzVector tau2P4 = httLeptonCollection[indexLeg2].getP4();
  std::cout<<13<<std::endl;

  bool tauBaselineSelection1 =  httLeptonCollection[indexLeg1].isFullHadLeadTau();
  std::cout<<14<<std::endl;
  bool tauBaselineSelection2 = httLeptonCollection[indexLeg2].isFullHadSubTau();
  std::cout<<15<<std::endl;

  bool baselinePair = tau1P4.DeltaR(tau2P4) > 0.5;
  std::cout<<16<<std::endl;
  bool boolAntiEleLeg1 = ( (int)httLeptonCollection[indexLeg1].getProperty(PropertyEnum::idAntiEle) & 0x1) == 0x1;   //Vloose AntiEle Id
  std::cout<<17<<std::endl;
  bool boolAntiEleLeg2 = ( (int)httLeptonCollection[indexLeg2].getProperty(PropertyEnum::idAntiEle) & 0x1) == 0x1;   //Vloose AntiEle Id
  std::cout<<18<<std::endl;
  bool boolAntiMuLeg1  = ( (int)httLeptonCollection[indexLeg1].getProperty(PropertyEnum::idAntiMu)  & 0x1) == 0x1;   //Loose AntiMu Id
  std::cout<<19<<std::endl;
  bool boolAntiMuLeg2  = ( (int)httLeptonCollection[indexLeg2].getProperty(PropertyEnum::idAntiMu)  & 0x1) == 0x1;   //Loose AntiMu Id
  std::cout<<20<<std::endl;
  bool boolAntiLeptonId = boolAntiEleLeg1 && boolAntiMuLeg1 && boolAntiEleLeg2 && boolAntiMuLeg2;
  std::cout<<21<<std::endl;

  bool boolExtraElectronVeto = thirdLeptonVeto(indexLeg1, indexLeg2, 11);
  std::cout<<22<<std::endl;
  bool boolExtraMuonVeto     = thirdLeptonVeto(indexLeg1, indexLeg2, 13);
  std::cout<<23<<std::endl;
  bool boolThirdLeptonVeto   = boolExtraMuonVeto || boolExtraElectronVeto;
  std::cout<<24<<std::endl;
  
  httEvent->clearSelectionWord();
  std::cout<<25<<std::endl;
  httEvent->setSelectionBit(SelectionBitsEnum::diMuonVeto,0); //only set explicitly for mutau
  std::cout<<26<<std::endl;
  httEvent->setSelectionBit(SelectionBitsEnum::diElectronVeto,0); //only set explicitly for etau  
  std::cout<<27<<std::endl;
  httEvent->setSelectionBit(SelectionBitsEnum::diLeptonVeto, 0);
  std::cout<<28<<std::endl;
  httEvent->setSelectionBit(SelectionBitsEnum::antiLeptonId, boolAntiLeptonId);
  std::cout<<29<<std::endl;
  httEvent->setSelectionBit(SelectionBitsEnum::extraMuonVeto,boolExtraMuonVeto);
  std::cout<<30<<std::endl;
  httEvent->setSelectionBit(SelectionBitsEnum::extraElectronVeto, boolExtraElectronVeto);
  std::cout<<31<<std::endl;
  httEvent->setSelectionBit(SelectionBitsEnum::thirdLeptonVeto, boolThirdLeptonVeto);
  std::cout<<32<<std::endl;

  return tauBaselineSelection1 && tauBaselineSelection2 && baselinePair
         && true;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
unsigned int HTauhTauhTreeFromNano::bestPair(std::vector<unsigned int> &pairIndexes){
  std::cout<<"unsigned int HTauhTauhTreeFromNano::bestPair(std::vector<unsigned int> &pairIndexes)"<<std::endl;
  unsigned int bestIndex = 9999;
  std::cout<<1<<std::endl;
  ///Pairs are already sorted during the ntuple creation?
  double iso_1=std::numeric_limits<double>::infinity(), iso_2=std::numeric_limits<double>::infinity(), pt_1=-1, pt_2=-1;
  std::cout<<2<<std::endl;
  if(pairIndexes.size()) {
    std::cout<<3<<std::endl;
    //return pairIndexes[0];//MB
    for(unsigned int ii=0;ii<2*pairIndexes.size();++ii){
      std::cout<<4<<std::endl;
      unsigned int i=(ii<pairIndexes.size()?ii:ii-pairIndexes.size());
      std::cout<<5<<std::endl;
      unsigned int iPair = pairIndexes[i];
      std::cout<<6<<std::endl;
      unsigned int indexLeg1 = httPairs_[iPair].getIndexLeg1();
      std::cout<<7<<std::endl;
      unsigned int indexLeg2 = httPairs_[iPair].getIndexLeg2();
      std::cout<<8<<std::endl;
      if(ii>=pairIndexes.size()){//invert legs
        std::cout<<9<<std::endl;
        indexLeg1 = httPairs_[iPair].getIndexLeg2();
        std::cout<<10<<std::endl;
        indexLeg2 = httPairs_[iPair].getIndexLeg1();
        std::cout<<11<<std::endl;
      }
      double pt_1_i = httLeptonCollection[indexLeg1].getP4().Pt();
      std::cout<<12<<std::endl;
      double pt_2_i = httLeptonCollection[indexLeg2].getP4().Pt();
      std::cout<<13<<std::endl;
      //MB: More isolated for MVAIso means higher value so inverted here to keep standard convention in comparison
      double iso_1_i = -httLeptonCollection[indexLeg1].getProperty(PropertyEnum::rawMVAoldDM2017v2);
      std::cout<<14<<std::endl;
      double iso_2_i = -httLeptonCollection[indexLeg2].getProperty(PropertyEnum::rawMVAoldDM2017v2);
      std::cout<<15<<std::endl;
      if(iso_1_i>iso_1) continue;
      std::cout<<16<<std::endl;
      if(iso_1_i==iso_1 && pt_1_i<pt_1) continue;
      std::cout<<17<<std::endl;
      if(iso_2_i>iso_2) continue;
      std::cout<<18<<std::endl;
      if(iso_2_i==iso_2 && pt_2_i<pt_2) continue;
      std::cout<<19<<std::endl;
      bestIndex = iPair;
      std::cout<<20<<std::endl;
      iso_1 = iso_1_i;
      iso_2 = iso_2_i;
      pt_1 = pt_1_i;
      pt_2 = pt_2_i;
      std::cout<<21<<std::endl;
    }
  }
  /*
  if(pairIndexes.size() && bestIndex!=pairIndexes[0]){
    unsigned int iPair = pairIndexes[0];
    unsigned int indexLeg1 = httPairs_[iPair].getIndexLeg1();
    unsigned int indexLeg2 = httPairs_[iPair].getIndexLeg2();
    double pt_1_i = httLeptonCollection[indexLeg1].getP4().Pt();
    double pt_2_i = httLeptonCollection[indexLeg2].getP4().Pt();
    std::cout<<"Pair sorting: "<<std::endl
             <<"best index = "<<bestIndex<<", index[0] = "<<pairIndexes[0]<<std::endl
             <<"\tiso1[best]="<<-iso_1<<", iso1[0]="<<httLeptonCollection[indexLeg1].getProperty(PropertyEnum::rawMVAoldDM)<<std::endl
             <<"\tpt1[best]="<<pt_1<<", pt1[0]="<<pt_1_i<<std::endl
             <<"\tiso2[best]="<<-iso_2<<", iso1[0]="<<httLeptonCollection[indexLeg2].getProperty(PropertyEnum::rawMVAoldDM)<<std::endl
             <<"\tpt2[best]="<<pt_2<<", pt1[0]="<<pt_2_i<<std::endl;
  }
  */

  return bestIndex;
};
/////////////////////////////////////////////////
/////////////////////////////////////////////////
