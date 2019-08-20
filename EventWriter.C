#include "EventWriter.h"

const float DEF = -10.;
const float DEFWEIGHT = 1.0;
const bool DEFFLAG = 0.;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fill(HTTEvent *ev, HTTJetCollection *jets, std::vector<HTTParticle> leptons, HTTPair *pair){
    // Make sure that current Collections are filled
    jets->setCurrentUncertShift( "", true );
    pair->setCurrentMETShift("");
  
    channel = pair->getFinalState();

    leg1 = pair->getLeg1();
    leg2 = pair->getLeg2();

    leg1P4=leg1.getP4();
    leg2P4=leg2.getP4();
  
    gen_match_1=leg1.getProperty(PropertyEnum::mc_match);
    gen_match_2=leg2.getProperty(PropertyEnum::mc_match);

    pdg1=std::abs(leg1.getProperty(PropertyEnum::pdgId));
    pdg2=std::abs(leg2.getProperty(PropertyEnum::pdgId));

    runID=ev->getRunId();

    lumiBlock=ev->getLSId();
    eventNr=ev->getEventId();
  
    npv=ev->getNPV();
    npu=ev->getNPU();
    rho=ev->getRho();
  
    fillLeg1Branches();
    
    fillLeg2Branches();
    
    fillJetBranches(jets);

    fillPairBranches(pair, jets);

    fillAdditionalLeptons( leptons, pair );

    fillMELA(jets);
  
    //////////////////////////////////////////////////////////////////  
 
    // Make sure that for 2016 badMuon filters are added
    Flag_goodVertices =                       ev->getFilter(FilterEnum::Flag_goodVertices);
    Flag_globalTightHalo2016Filter =          ev->getFilter(FilterEnum::Flag_globalTightHalo2016Filter);
    Flag_globalSuperTightHalo2016Filter =     ev->getFilter(FilterEnum::Flag_globalSuperTightHalo2016Filter);
    Flag_HBHENoiseFilter =                    ev->getFilter(FilterEnum::Flag_HBHENoiseFilter);
    Flag_HBHENoiseIsoFilter =                 ev->getFilter(FilterEnum::Flag_HBHENoiseIsoFilter);
    Flag_EcalDeadCellTriggerPrimitiveFilter = ev->getFilter(FilterEnum::Flag_EcalDeadCellTriggerPrimitiveFilter);
    Flag_BadPFMuonFilter =                    ev->getFilter(FilterEnum::Flag_BadPFMuonFilter);
    Flag_BadChargedCandidateFilter =          ev->getFilter(FilterEnum::Flag_BadChargedCandidateFilter);
    Flag_eeBadScFilter =                      ev->getFilter(FilterEnum::Flag_eeBadScFilter);
    Flag_ecalBadCalibFilter =                 ev->getFilter(FilterEnum::Flag_ecalBadCalibFilter);
    Flag_METFilters =                         ev->getFilter(FilterEnum::Flag_METFilters);    

    flagMETFilter = Flag_goodVertices
                    && Flag_globalSuperTightHalo2016Filter
                    && Flag_HBHENoiseFilter
                    && Flag_HBHENoiseIsoFilter
                    && Flag_EcalDeadCellTriggerPrimitiveFilter
                    && Flag_BadPFMuonFilter
                    && Flag_BadChargedCandidateFilter
                    && Flag_ecalBadCalibFilter;

    if(!isMC) flagMETFilter &= Flag_eeBadScFilter; // only for data
  
    //////////////////////////////////////////////////////////////////  

    // Fill vetos 
    extramuon_veto=ev->checkSelectionBit(SelectionBitsEnum::extraMuonVeto);
    extraelec_veto=ev->checkSelectionBit(SelectionBitsEnum::extraElectronVeto);
    diMuonVeto= ev->checkSelectionBit(SelectionBitsEnum::diMuonVeto);
    diElectronVeto= ev->checkSelectionBit(SelectionBitsEnum::diElectronVeto);
    dilepton_veto = ev->checkSelectionBit(SelectionBitsEnum::diLeptonVeto);

    passesDiMuonVeto=!diMuonVeto;           // silly to pass a veto...
    passesDiElectronVeto=!diElectronVeto;   // kept for backward compatibility

    passesThirdLepVeto=!ev->checkSelectionBit(SelectionBitsEnum::thirdLeptonVeto); 
    passesTauLepVetos = ev->checkSelectionBit(SelectionBitsEnum::antiLeptonId);

    //////////////////////////////////////////////////////////////////

    // Fill stxs related branches
    htxs_stage1cat = ev->getStage1Cat();
    std::vector<double> THU = ev->getTHU_uncertainties();
    THU_ggH_Mu =    THU[0];
    THU_ggH_Res =   THU[1];
    THU_ggH_Mig01 = THU[2];
    THU_ggH_Mig12 = THU[3];
    THU_ggH_VBF2j = THU[4];
    THU_ggH_VBF3j = THU[5];
    THU_ggH_PT60 =  THU[6];
    THU_ggH_PT120 = THU[7];
    THU_ggH_qmtop = THU[8];

    NNLO_ggH_weight = ev->getNNLO_ggH_weight();
  
    //////////////////////////////////////////////////////////////////
    //this is quite slow, calling the function for each trigger item...
    if ( channel == HTTAnalysis::MuTau ) //mu-tau
    {
        trg_singlemuon_22        = leg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu22) && pt_1 > 23;
        trg_singlemuon_22_eta2p1 = leg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu22_eta2p1) && pt_1 > 23;

        trg_singlemuonTk_22        = leg1.hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22) && pt_1 > 23;
        trg_singlemuonTk_22_eta2p1 = leg1.hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22_eta2p1) && pt_1 > 23;

        trg_singlemu_22 = trg_singlemuon_22 || trg_singlemuon_22_eta2p1 || trg_singlemuonTk_22 || trg_singlemuonTk_22_eta2p1;
        
        trg_crossmuon_mu19tau20 =  leg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu19_eta2p1_LooseIsoPFTau20) 
                                  && leg2.hasTriggerMatch(TriggerEnum::HLT_IsoMu19_eta2p1_LooseIsoPFTau20)
                                  && pt_1 > 20 && pt_2 > 32
                                  && abs(eta_1) < 2.1 && abs(eta_2) < 2.1;

        trg_crossmuon_mu19tau20_singleL1 =  leg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1) 
                                            && leg2.hasTriggerMatch(TriggerEnum::HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1)
                                            && pt_1 > 20 && pt_2 > 32
                                            && abs(eta_1) < 2.1 && abs(eta_2) < 2.1;

        trg_crossmu_mu19tau20 = trg_crossmuon_mu19tau20 || trg_crossmuon_mu19tau20_singleL1;

    }else if ( channel == HTTAnalysis::EleTau )
    {
        trg_singleelectron_25_eta2p1 =    leg1.hasTriggerMatch(TriggerEnum::HLT_Ele25_eta2p1_WPTight_Gsf)  && pt_1 > 26;
        
    } else if ( channel == HTTAnalysis::TauTau )
    {
        trg_doubletau_35_mediso_eta2p1 = ( leg1.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg) && leg2.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg) ) && pt_1 > 40 && pt_2 > 40;
        
        trg_doubletau_35_medCombiso_eta2p1 = ( leg1.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg) && leg2.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg) ) && pt_1 > 40 && pt_2 > 40;
        trg_doubletau_35 = trg_doubletau_35_mediso_eta2p1 || trg_doubletau_35_medCombiso_eta2p1;
    }
    
    TLorentzVector ll=ev->getGenBosonP4(false);
    TLorentzVector llvis=ev->getGenBosonP4(true);
    gen_Mll=ll.M();
    gen_ll_px=ll.Px();
    gen_ll_py=ll.Py();
    gen_ll_pz=ll.Pz();
    gen_vis_Mll=llvis.M();
    gen_ll_vis_px=llvis.Px();
    gen_ll_vis_py=llvis.Py();
    gen_ll_vis_pz=llvis.Pz();

    topPtReweightWeightRun2=ev->getTopPtReWeight(false);
    topPtReweightWeightRun1=ev->getTopPtReWeight(true);
   
    if(ev->getSampleType() == HTTEvent::DY || ev->getSampleType() == HTTEvent::DYLowM ) {
        w->var("z_gen_mass")->setVal(  ll.M()  );
        w->var("z_gen_pt")->setVal( ll.Pt() );
        zPtReweightWeight=w->function("zptmass_weight_nom")->getVal();
    }

    if(isMC)
    {
        xsec = ev->getXsec();
        genWeight = ev->getMCWeight();
        genNEventsWeight = ev->getGenNEventsWeight();
        puWeight = ev->getPUWeight();
        lumiWeight=1000*xsec*genNEventsWeight;
        NUP=ev->getLHEnOutPartons();
        fillStitchingWeight( ev->getSampleType() );
        fillScalefactors();
        fillLeptonFakeRateWeights();

    } else
    {
        xsec = DEFWEIGHT;
        genWeight = DEFWEIGHT;        
        genNEventsWeight = DEFWEIGHT;
        puWeight = DEFWEIGHT;
        weight = DEFWEIGHT;
        NUP=-1;
        lumiWeight=DEFWEIGHT;
    }
   
    // Make sure this weight is up to date. 
    // Trigger sf need to be applied according to what triggers are used
    weight = puWeight*stitchedWeight*genWeight*eleTauFakeRateWeight*muTauFakeRateWeight*idisoweight_1*idisoweight_2*sf_trk*sf_reco;
   
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillStitchingWeight(HTTEvent::sampleTypeEnum sampleType)
{
    if(sampleType == HTTEvent::DY)
    {
        if(NUP == 0) stitchedWeight = 0.03941307;
        if(NUP == 1) stitchedWeight = 0.01132200;
        if(NUP == 2) stitchedWeight = 0.01175009;
        if(NUP == 3) stitchedWeight = 0.01203609;
        if(NUP == 4) stitchedWeight = 0.00980626;


    } else if(sampleType == HTTEvent::WJets)
    {
        if(NUP == 0) stitchedWeight = 0.70788322;
        if(NUP == 1) stitchedWeight = 0.16372267;
        if(NUP == 2) stitchedWeight = 0.04846574;
        if(NUP == 3) stitchedWeight = 0.01574309;
        if(NUP == 4) stitchedWeight = 0.01585509;

    } else stitchedWeight=lumiWeight;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillLeptonFakeRateWeights()
{
    eleTauFakeRateWeight = 1.0;
    muTauFakeRateWeight = 1.0;
    antilep_tauscaling = 1.0;

    //values taken from here: https://indico.cern.ch/event/803335/contributions/3359969/attachments/1829820/2996253/TauPOG_HTT_workshop_20190415_v0.pdf
    // values taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingLegacyRun2#mu-%3Etau%20FR
    if(channel == HTTAnalysis::MuTau)
    {

        if((gen_match_2 == 1 || gen_match_2 == 3) && againstElectronVLooseMVA6_2 > 0.5 )
        {
            if( std::abs(eta_2) < 1.460 )       eleTauFakeRateWeight *= 1.21;
            else if ( std::abs(eta_2) > 1.559 ) eleTauFakeRateWeight *= 1.38;
        }
        if((gen_match_2 == 2 || gen_match_2 == 4) && againstMuonTight3_2 > 0.5 )
        {
            if( std::abs(eta_2) < 0.4 )       muTauFakeRateWeight *= 1.47;
            else if ( std::abs(eta_2) < 0.8 ) muTauFakeRateWeight *= 1.55;
            else if ( std::abs(eta_2) < 1.2 ) muTauFakeRateWeight *= 1.33;
            else if ( std::abs(eta_2) < 1.7 ) muTauFakeRateWeight *= 1.72;
            else if ( std::abs(eta_2) < 2.3 ) muTauFakeRateWeight *= 2.5;
        }        
    }
    if(channel == HTTAnalysis::EleTau)
    {
        if((gen_match_2 == 1 || gen_match_2 == 3) && againstElectronTightMVA6_2 > 0.5 )
        {
            if( std::abs(eta_2) < 1.460 )       eleTauFakeRateWeight *= 1.40;
            else if ( std::abs(eta_2) > 1.559 ) eleTauFakeRateWeight *= 1.90;
        }
        if((gen_match_2 == 2 || gen_match_2 == 4) && againstMuonLoose3_2> 0.5 )
        {
            if( std::abs(eta_2) < 0.4 )       muTauFakeRateWeight *= 1.22;
            else if ( std::abs(eta_2) < 0.8 ) muTauFakeRateWeight *= 1.12;
            else if ( std::abs(eta_2) < 1.2 ) muTauFakeRateWeight *= 1.26;
            else if ( std::abs(eta_2) < 1.7 ) muTauFakeRateWeight *= 1.22;
            else if ( std::abs(eta_2) < 2.3 ) muTauFakeRateWeight *= 2.39;
        }        
    }
    if(channel == HTTAnalysis::TauTau)
    {

        if((gen_match_1 == 1 || gen_match_1 == 3) && againstElectronVLooseMVA6_1 > 0.5 )
        {
            if( std::abs(eta_1) < 1.460 )       eleTauFakeRateWeight *= 1.21;
            else if ( std::abs(eta_1) > 1.559 ) eleTauFakeRateWeight *= 1.38;
        }
        if((gen_match_1 == 2 || gen_match_1 == 4) && againstMuonLoose3_1> 0.5 )
        {
            if( std::abs(eta_1) < 0.4 )       muTauFakeRateWeight *= 1.22;
            else if ( std::abs(eta_1) < 0.8 ) muTauFakeRateWeight *= 1.12;
            else if ( std::abs(eta_1) < 1.2 ) muTauFakeRateWeight *= 1.26;
            else if ( std::abs(eta_1) < 1.7 ) muTauFakeRateWeight *= 1.22;
            else if ( std::abs(eta_1) < 2.3 ) muTauFakeRateWeight *= 2.39;
        } 

        if((gen_match_2 == 1 || gen_match_2 == 3) && againstElectronVLooseMVA6_2 > 0.5 )
        {
            if( std::abs(eta_2) < 1.460 )       eleTauFakeRateWeight *= 1.21;
            else if ( std::abs(eta_2) > 1.559 ) eleTauFakeRateWeight *= 1.38;
        }
        if((gen_match_2 == 2 || gen_match_2 == 4) && againstMuonLoose3_2> 0.5 )
        {
            if( std::abs(eta_2) < 0.4 )       muTauFakeRateWeight *= 1.22;
            else if ( std::abs(eta_2) < 0.8 ) muTauFakeRateWeight *= 1.12;
            else if ( std::abs(eta_2) < 1.2 ) muTauFakeRateWeight *= 1.26;
            else if ( std::abs(eta_2) < 1.7 ) muTauFakeRateWeight *= 1.22;
            else if ( std::abs(eta_2) < 2.3 ) muTauFakeRateWeight *= 2.39;
        }        
    }
    antilep_tauscaling = eleTauFakeRateWeight * muTauFakeRateWeight;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillScalefactors()
{
    // from https://github.com/CMS-HTT/CorrectionsWorkspace/tree/2017_17NovReRecoData_Fall17MC
    singleTriggerSFLeg1 = DEFWEIGHT;
    singleTriggerSFLeg2 = DEFWEIGHT;

    double s1_data = DEFWEIGHT;
    double s1_mc   = DEFWEIGHT;
    double x1_data = DEFWEIGHT;
    double x1_mc   = DEFWEIGHT;
    double x2_data = DEFWEIGHT;
    double x2_mc   = DEFWEIGHT;

    isoWeight_1 = DEFWEIGHT;
    isoWeight_2 = DEFWEIGHT;
    idWeight_1 = DEFWEIGHT;
    idWeight_2 = DEFWEIGHT;

    idisoweight_1 = DEFWEIGHT;
    idisoweight_2 = DEFWEIGHT;

    // from https://github.com/truggles/TauTriggerSFs2017/tree/tauTriggers2017_MCv2_PreReMiniaod
    xTriggerSFLeg1 = DEFWEIGHT;
    xTriggerSFLeg2 = DEFWEIGHT;

    sf_trk =DEFWEIGHT;
    sf_reco=DEFWEIGHT;
    sf_SingleOrCrossTrigger = DEFWEIGHT;
    sf_SingleXorCrossTrigger = DEFWEIGHT;
    sf_SingleTrigger = DEFWEIGHT;
    sf_DoubleTauTight = DEFWEIGHT;
    sf_DoubleTauVTight = DEFWEIGHT;


    trk_sf = DEFWEIGHT;
    trigweight_1 = DEFWEIGHT;
    trigweight_2 = DEFWEIGHT;
    trigweightVTight_1 = DEFWEIGHT;
    trigweightVTight_2 = DEFWEIGHT;      

    idisoweight_1 = DEFWEIGHT;
    idisoweight_2 = DEFWEIGHT;

    // FIXME needs to be changed for 2016!
    if( channel == HTTAnalysis::MuTau )
    {
        w->var("m_pt")->setVal(  pt_1  );
        w->var("m_eta")->setVal( eta_1 );
        
        // w->var("t_pt")->setVal(  pt_2  );
        // w->var("t_eta")->setVal( eta_2 );
        // w->var("t_dm")->setVal( decayMode_2 );
        
        idisoweight_1 = w->function("m_idiso_desy_ratio")->getVal();

        if(pt_1 > 23){
            trigweight_1 = w->function("m_trgIsoMu22_desy_ratio")->getVal();
        }
        // else if(pt_1 < 23){
        //     trigweight_1 = w->function("m_trgMu19leg_eta2p1_desy_ratio")->getVal();
        //     trigweight_2 = w->function("t_genuine_TightIso_mt_ratio")->getVal();
        // }
        
        trk_sf =  w->function("m_trk_ratio")->getVal();
        sf_SingleTrigger = trigweight_1;
    }

    if( channel == HTTAnalysis::EleTau )
    {

        w->var("e_pt")->setVal(  pt_1  );
        w->var("e_eta")->setVal( eta_1 );


        idisoweight_1 = w->function("e_idiso_desy_ratio")->getVal();

        if( std::abs(eta_2) < 2.1 ){
            trigweight_1 = w->function("e_trgEle25leg_desy_ratio")->getVal();
        }
       

        trk_sf = w->function("e_trk_ratio")->getVal();
        

        sf_SingleTrigger = trigweight_1;
    }

    if( channel == HTTAnalysis::TauTau )
    {



        w->var("t_pt")->setVal( pt_1 );
        w->var("t_eta")->setVal( eta_1 );
        w->var("t_dm")->setVal( decayMode_1 );

        if( gen_match_1 == 5){
            trigweight_1 = w->function("t_genuine_TightIso_tt_ratio")->getVal();
            trigweightVTight_1 = w->function("t_genuine_VTightIso_tt_ratio")->getVal();
        }
        else{
          trigweight_1 = w->function("t_fake_TightIso_tt_ratio")->getVal();
          trigweightVTight_1 = w->function("t_fake_VTightIso_tt_ratio")->getVal();
        }

         w->var("t_pt")->setVal( pt_2 );
        w->var("t_eta")->setVal( eta_2 );
        w->var("t_dm")->setVal( decayMode_2 );
        
        if(gen_match_2 == 5){
            trigweight_2 = w->function("t_genuine_TightIso_tt_ratio")->getVal();
            trigweightVTight_2 = w->function("t_genuine_VTightIso_tt_ratio")->getVal();
        }
        else{
            trigweight_2 = w->function("t_fake_TightIso_tt_ratio")->getVal();
            trigweightVTight_2 = w->function("t_fake_VTightIso_tt_ratio")->getVal();
        }

      
        sf_DoubleTauVTight = trigweight_1 * trigweight_2;

    }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillLeg1Branches()
{
    UChar_t bitmask;

    pt_1=leg1P4.Pt();
    phi_1=leg1P4.Phi();
    eta_1=leg1P4.Eta();
    eta_SC_1=eta_1+leg1.getProperty(PropertyEnum::deltaEtaSC);
    m_1=leg1P4.M();
    q_1=leg1.getCharge();
    d0_1=leg1.getProperty(PropertyEnum::dxy);
    dZ_1=leg1.getProperty(PropertyEnum::dz);


    if      (pdg1==15) iso_1=leg1.getProperty( HTTEvent::usePropertyFor.at("tauIsolation") );
    else if (pdg1==13) iso_1=leg1.getProperty( HTTEvent::usePropertyFor.at("muonIsolation") );
    else if (pdg1==11) iso_1=leg1.getProperty( HTTEvent::usePropertyFor.at("electronIsolation") );;

    // bitmask=leg1.getProperty(PropertyEnum::idAntiEle);

    againstElectronMVA6_1 = leg1.getProperty(PropertyEnum::idAntiEle);
    againstElectronVLooseMVA6_1 =(againstElectronMVA6_1 & 0x1 )>0;
    againstElectronLooseMVA6_1  =(againstElectronMVA6_1 & 0x2 )>0;
    againstElectronMediumMVA6_1 =(againstElectronMVA6_1 & 0x4 )>0;
    againstElectronTightMVA6_1  =(againstElectronMVA6_1 & 0x8 )>0;
    againstElectronVTightMVA6_1 =(againstElectronMVA6_1 & 0x10)>0;

    againstMuon3_1=leg1.getProperty(PropertyEnum::idAntiMu);
    againstMuonLoose3_1=(againstMuon3_1 & 0x1)>0;
    againstMuonTight3_1=(againstMuon3_1 & 0x2)>0;

    byCombinedIsolationDeltaBetaCorrRaw3Hits_1=leg1.getProperty(PropertyEnum::rawIso);

    byIsolationMVA3oldDMwLTraw_1 =leg1.getProperty(HTTEvent::usePropertyFor.at("tauID"));

    // Tau ID based on BDT trained with 2017 data
    byIsolationMVArun2017v2DBoldDMwLTraw2017_1=leg1.getProperty(HTTEvent::usePropertyFor.at("tauIsolation")); // Raw value
    byIsolationMVArun2017v2DBoldDMwLT2017_1=leg1.getProperty(HTTEvent::usePropertyFor.at("tauID"));    // Bitmask

    byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x1)>0;
    byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1  = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x2)>0;
    byLooseIsolationMVArun2017v2DBoldDMwLT2017_1   = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x4)>0;
    byMediumIsolationMVArun2017v2DBoldDMwLT2017_1  = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x8)>0;
    byTightIsolationMVArun2017v2DBoldDMwLT2017_1   = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x10)>0;
    byVTightIsolationMVArun2017v2DBoldDMwLT2017_1  = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x20)>0;
    byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1  = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x40)>0;
    // std::cout<<"byIsolationMVArun2017v2DBoldDMwLTraw2017_1 : "<<byIsolationMVArun2017v2DBoldDMwLTraw2017_1 << std::endl;
    // std::cout<<"byIsolationMVArun2017v2DBoldDMwLT2017_1 : "<<byIsolationMVArun2017v2DBoldDMwLT2017_1 << std::endl;

    // Tau ID based on DeepTauIDv2 against jets
    byIsolationDeepTau2017v2VSjet_raw_1=leg1.getProperty(HTTEvent::usePropertyFor.at("tauIsolationDeepVSjet")); // Raw value
    byIsolationDeepTau2017v2VSjet_1=leg1.getProperty(HTTEvent::usePropertyFor.at("tauIDDeepVSjet"));    // Bitmask
    // std::cout<<"byIsolationDeepTau2017v2VSjet_raw_1 : "<<byIsolationDeepTau2017v2VSjet_raw_1 << std::endl;
    // std::cout<<"byIsolationDeepTau2017v2VSjet_1 : "<<byIsolationDeepTau2017v2VSjet_1 << std::endl;

    byVVVLooseIsolationDeepTau2017v2VSjet_1 = (byIsolationDeepTau2017v2VSjet_1 & 0x1)>0;
    byVVLooseIsolationDeepTau2017v2VSjet_1  = (byIsolationDeepTau2017v2VSjet_1 & 0x2)>0;
    byVLooseIsolationDeepTau2017v2VSjet_1   = (byIsolationDeepTau2017v2VSjet_1 & 0x4)>0;
    byLooseIsolationDeepTau2017v2VSjet_1  = (byIsolationDeepTau2017v2VSjet_1 & 0x8)>0;
    byMediumIsolationDeepTau2017v2VSjet_1   = (byIsolationDeepTau2017v2VSjet_1 & 0x10)>0;
    byTightIsolationDeepTau2017v2VSjet_1  = (byIsolationDeepTau2017v2VSjet_1 & 0x20)>0;
    byVTightIsolationDeepTau2017v2VSjet_1  = (byIsolationDeepTau2017v2VSjet_1 & 0x40)>0; // 4*16 = 64   
    byVVTightIsolationDeepTau2017v2VSjet_1  = (byIsolationDeepTau2017v2VSjet_1 & 0x80)>0; //8*16     = 128

    // Tau ID based on DeepTauIDv2 against electrons
    byIsolationDeepTau2017v2VSe_raw_1=leg1.getProperty(HTTEvent::usePropertyFor.at("tauIsolationDeepVSele")); // Raw value
    byIsolationDeepTau2017v2VSe_1=leg1.getProperty(HTTEvent::usePropertyFor.at("tauIDDeepVSele"));    // Bitmask

    byVVVLooseIsolationDeepTau2017v2VSe_1 = (byIsolationDeepTau2017v2VSe_1 & 0x1)>0;
    byVVLooseIsolationDeepTau2017v2VSe_1  = (byIsolationDeepTau2017v2VSe_1 & 0x2)>0;
    byVLooseIsolationDeepTau2017v2VSe_1   = (byIsolationDeepTau2017v2VSe_1 & 0x4)>0;
    byLooseIsolationDeepTau2017v2VSe_1  = (byIsolationDeepTau2017v2VSe_1 & 0x8)>0;
    byMediumIsolationDeepTau2017v2VSe_1   = (byIsolationDeepTau2017v2VSe_1 & 0x10)>0;
    byTightIsolationDeepTau2017v2VSe_1  = (byIsolationDeepTau2017v2VSe_1 & 0x20)>0;
    byVTightIsolationDeepTau2017v2VSe_1  = (byIsolationDeepTau2017v2VSe_1 & 0x40)>0; // 4*16 = 64   
    byVVTightIsolationDeepTau2017v2VSe_1  = (byIsolationDeepTau2017v2VSe_1 & 0x80)>0; //8*16     = 128

    // Tau ID based on DeepTauIDv2 against muons
    byIsolationDeepTau2017v2VSmu_raw_1=leg1.getProperty(HTTEvent::usePropertyFor.at("tauIsolationDeepVSmu")); // Raw value
    byIsolationDeepTau2017v2VSmu_1=leg1.getProperty(HTTEvent::usePropertyFor.at("tauIDDeepVSmu"));    // Bitmask

    byVLooseIsolationDeepTau2017v2VSmu_1   = (byIsolationDeepTau2017v2VSmu_1 & 0x1)>0;
    byLooseIsolationDeepTau2017v2VSmu_1  = (byIsolationDeepTau2017v2VSmu_1 & 0x2)>0;
    byMediumIsolationDeepTau2017v2VSmu_1   = (byIsolationDeepTau2017v2VSmu_1 & 0x4)>0;
    byTightIsolationDeepTau2017v2VSmu_1  = (byIsolationDeepTau2017v2VSmu_1 & 0x8)>0;
   
    chargedIsoPtSum_1=leg1.getProperty(PropertyEnum::chargedIso);
    neutralIsoPtSum_1=leg1.getProperty(PropertyEnum::neutralIso);
    puCorrPtSum_1=leg1.getProperty(PropertyEnum::puCorr);
    decayModeFindingOldDMs_1=leg1.getProperty(PropertyEnum::idDecayMode);
    decayModeFindingNewDMs_1=leg1.getProperty(PropertyEnum::idDecayModeNewDMs);
    decayMode_1=leg1.getProperty(PropertyEnum::decayMode);

    if (pdg1==13) id_m_loose_1=1; //already filtered at NanoAOD production
    id_m_medium_1=leg1.getProperty(PropertyEnum::mediumId);
    id_m_tight_1=leg1.getProperty(PropertyEnum::tightId);
    bitmask=leg1.getProperty(PropertyEnum::highPtId);
    id_m_highpt_1=(bitmask & 0x2)>0;

    int intmask=leg1.getProperty(PropertyEnum::cutBased);
    id_e_cut_veto_1=intmask>=1;
    id_e_cut_loose_1=intmask>=2;
    id_e_cut_medium_1=intmask>=3;
    id_e_cut_tight_1=intmask>=4;
    id_e_mva_nt_loose_1=DEF;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillLeg2Branches()
{
    UChar_t bitmask;

    pt_2=leg2P4.Pt();
    phi_2=leg2P4.Phi();
    eta_2=leg2P4.Eta();
    m_2=leg2P4.M();
    q_2=leg2.getCharge();
    d0_2=leg2.getProperty(PropertyEnum::dxy);
    dZ_2=leg2.getProperty(PropertyEnum::dz);


    if      (pdg2==15) iso_2=leg2.getProperty( HTTEvent::usePropertyFor.at("tauIsolation") );
    else if (pdg2==13) iso_2=leg2.getProperty( HTTEvent::usePropertyFor.at("muonIsolation") );
    else if (pdg2==11) iso_2=leg2.getProperty( HTTEvent::usePropertyFor.at("electronIsolation") );;


    againstElectronMVA6_2=leg2.getProperty(PropertyEnum::idAntiEle);
    againstElectronVLooseMVA6_2 = (againstElectronMVA6_2 & 0x1)>0;
    againstElectronLooseMVA6_2  = (againstElectronMVA6_2 & 0x2)>0;
    againstElectronMediumMVA6_2 = (againstElectronMVA6_2 & 0x4)>0;
    againstElectronTightMVA6_2  = (againstElectronMVA6_2 & 0x8)>0;
    againstElectronVTightMVA6_2 = (againstElectronMVA6_2 & 0x10)>0;

    againstMuon3_2=leg2.getProperty(PropertyEnum::idAntiMu);
    againstMuonLoose3_2=(againstMuon3_2 & 0x1)>0;
    againstMuonTight3_2=(againstMuon3_2 & 0x2)>0;

    byCombinedIsolationDeltaBetaCorrRaw3Hits_2=leg2.getProperty(PropertyEnum::rawIso);

    byIsolationMVA3oldDMwLTraw_2 =leg2.getProperty(PropertyEnum::rawMVAoldDM2017v2); //same as above!?

    byIsolationMVArun2017v2DBoldDMwLTraw2017_2=leg2.getProperty( HTTEvent::usePropertyFor.at("tauIsolation") ); // Raw value
    byIsolationMVArun2017v2DBoldDMwLT2017_2=leg2.getProperty( HTTEvent::usePropertyFor.at("tauID") );     // Bitmask

    byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x1)>0;
    byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2  = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x2)>0;
    byLooseIsolationMVArun2017v2DBoldDMwLT2017_2   = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x4)>0;;
    byMediumIsolationMVArun2017v2DBoldDMwLT2017_2  = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x8)>0;
    byTightIsolationMVArun2017v2DBoldDMwLT2017_2   = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x10)>0;
    byVTightIsolationMVArun2017v2DBoldDMwLT2017_2  = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x20)>0;
    byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2  = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x40)>0;
    // std::cout<<"byIsolationMVArun2017v2DBoldDMwLTraw2017_2 : "<<byIsolationMVArun2017v2DBoldDMwLTraw2017_2 << std::endl;
    // std::cout<<"byIsolationMVArun2017v2DBoldDMwLT2017_2 : "<<byIsolationMVArun2017v2DBoldDMwLT2017_2 << std::endl;

    // Tau ID based on DeepTauIDv2 against jets
    byIsolationDeepTau2017v2VSjet_raw_2=leg2.getProperty(HTTEvent::usePropertyFor.at("tauIsolationDeepVSjet")); // Raw value
    byIsolationDeepTau2017v2VSjet_2=leg2.getProperty(HTTEvent::usePropertyFor.at("tauIDDeepVSjet"));    // Bitmask
    // std::cout<<"byIsolationDeepTau2017v2VSjet_raw_2 : "<<byIsolationDeepTau2017v2VSjet_raw_2 << std::endl;
    // std::cout<<"byIsolationDeepTau2017v2VSjet_2 : "<<byIsolationDeepTau2017v2VSjet_2 << std::endl;
    
    byVVVLooseIsolationDeepTau2017v2VSjet_2 = (byIsolationDeepTau2017v2VSjet_2 & 0x1)>0;
    byVVLooseIsolationDeepTau2017v2VSjet_2  = (byIsolationDeepTau2017v2VSjet_2 & 0x2)>0;
    byVLooseIsolationDeepTau2017v2VSjet_2   = (byIsolationDeepTau2017v2VSjet_2 & 0x4)>0;
    byLooseIsolationDeepTau2017v2VSjet_2  = (byIsolationDeepTau2017v2VSjet_2 & 0x8)>0;
    byMediumIsolationDeepTau2017v2VSjet_2   = (byIsolationDeepTau2017v2VSjet_2 & 0x10)>0;
    byTightIsolationDeepTau2017v2VSjet_2  = (byIsolationDeepTau2017v2VSjet_2 & 0x20)>0;
    byVTightIsolationDeepTau2017v2VSjet_2  = (byIsolationDeepTau2017v2VSjet_2 & 0x40)>0; // 4*16 = 64   
    byVVTightIsolationDeepTau2017v2VSjet_2  = (byIsolationDeepTau2017v2VSjet_2 & 0x80)>0; //8*16     = 128

    // Tau ID based on DeepTauIDv2 against electrons
    byIsolationDeepTau2017v2VSe_raw_2=leg2.getProperty(HTTEvent::usePropertyFor.at("tauIsolationDeepVSele")); // Raw value
    byIsolationDeepTau2017v2VSe_2=leg2.getProperty(HTTEvent::usePropertyFor.at("tauIDDeepVSele"));    // Bitmask

    byVVVLooseIsolationDeepTau2017v2VSe_2 = (byIsolationDeepTau2017v2VSe_2 & 0x1)>0;
    byVVLooseIsolationDeepTau2017v2VSe_2  = (byIsolationDeepTau2017v2VSe_2 & 0x2)>0;
    byVLooseIsolationDeepTau2017v2VSe_2   = (byIsolationDeepTau2017v2VSe_2 & 0x4)>0;
    byLooseIsolationDeepTau2017v2VSe_2  = (byIsolationDeepTau2017v2VSe_2 & 0x8)>0;
    byMediumIsolationDeepTau2017v2VSe_2   = (byIsolationDeepTau2017v2VSe_2 & 0x10)>0;
    byTightIsolationDeepTau2017v2VSe_2  = (byIsolationDeepTau2017v2VSe_2 & 0x20)>0;
    byVTightIsolationDeepTau2017v2VSe_2  = (byIsolationDeepTau2017v2VSe_2 & 0x40)>0; // 4*16 = 64   
    byVVTightIsolationDeepTau2017v2VSe_2  = (byIsolationDeepTau2017v2VSe_2 & 0x80)>0; //8*16     = 128

    // Tau ID based on DeepTauIDv2 against muons
    byIsolationDeepTau2017v2VSmu_raw_2=leg2.getProperty(HTTEvent::usePropertyFor.at("tauIsolationDeepVSmu")); // Raw value
    byIsolationDeepTau2017v2VSmu_2=leg2.getProperty(HTTEvent::usePropertyFor.at("tauIDDeepVSmu"));    // Bitmask
    // std::cout<<"byIsolationDeepTau2017v2VSmu_raw_2 : "<<byIsolationDeepTau2017v2VSmu_raw_2 << std::endl;
    // std::cout<<"byIsolationDeepTau2017v2VSmu_2 : "<<byIsolationDeepTau2017v2VSmu_2 << std::endl;
    
    byVLooseIsolationDeepTau2017v2VSmu_2   = (byIsolationDeepTau2017v2VSmu_2 & 0x1)>0;
    byLooseIsolationDeepTau2017v2VSmu_2  = (byIsolationDeepTau2017v2VSmu_2 & 0x2)>0;
    byMediumIsolationDeepTau2017v2VSmu_2   = (byIsolationDeepTau2017v2VSmu_2 & 0x4)>0;
    byTightIsolationDeepTau2017v2VSmu_2  = (byIsolationDeepTau2017v2VSmu_2 & 0x8)>0;

    chargedIsoPtSum_2=leg2.getProperty(PropertyEnum::chargedIso);
    neutralIsoPtSum_2=leg2.getProperty(PropertyEnum::neutralIso);
    puCorrPtSum_2=leg2.getProperty(PropertyEnum::puCorr);
    decayModeFindingOldDMs_2=leg2.getProperty(PropertyEnum::idDecayMode);
    decayModeFindingNewDMs_2=leg2.getProperty(PropertyEnum::idDecayModeNewDMs);
    decayMode_2=leg2.getProperty(PropertyEnum::decayMode);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillJetBranches(HTTJetCollection *jets)
{
    for(unsigned int shift = 0; shift<jecShifts.size(); ++shift )
    {

        jets->setCurrentUncertShift( jecShifts[shift].second.first, jecShifts[shift].second.second);
        njets[shift]        = jets->getNJets(30);
        njetspt20[shift]    = jets->getNJets(20);
        njetingap[shift]    = jets->getNJetInGap(30);
        njetingap20[shift]  = jets->getNJetInGap(20);        


        if ( njetspt20[shift] >=1 )
        {
            jpt_1[shift]   = jets->getJet(0).Pt();
            jeta_1[shift]  = jets->getJet(0).Eta();
            jphi_1[shift]  = jets->getJet(0).Phi();
            if( strcmp(jecShifts[shift].first.c_str(), "") == 0 )
            {
                jm_1=   jets->getJet(0).M();
                jrawf_1=jets->getJet(0).getProperty(PropertyEnum::rawFactor);
                //jmva_1= jets->getJet(0).getProperty(PropertyEnum::puIdDisc);
                jcsv_1= jets->getJet(0).getProperty(PropertyEnum::btagDeepB);
            }
        }
        if ( njetspt20[shift] >=2 )
        {

            jpt_2[shift]   =  jets->getJet(1).Pt();
            jeta_2[shift]  = jets->getJet(1).Eta();
            jphi_2[shift]  = jets->getJet(1).Phi();
            
            mjj[shift]      = jets->getDiJetP4().M();
            dijetpt[shift]  = jets->getDiJetP4().Pt();
            dijetphi[shift] = jets->getDiJetP4().Phi();
            jdeta[shift]    = std::abs( jets->getJet(0).Eta()-jets->getJet(1).Eta() );
            jdphi[shift]    = jets->getJet(0).P4().DeltaPhi( jets->getJet(1).P4() );

            if( strcmp(jecShifts[shift].first.c_str(), "") == 0 )
            {
                jm_2=   jets->getJet(1).M();
                jrawf_2=jets->getJet(1).getProperty(PropertyEnum::rawFactor);
                //jmva_2= jets->getJet(1).getProperty(PropertyEnum::puIdDisc);
                jcsv_2= jets->getJet(1).getProperty(PropertyEnum::btagDeepB);
                jeta1eta2=jeta_1[shift]*jeta_2[shift];
                lep_etacentrality=TMath::Exp( -4/pow(jeta_1[shift]-jeta_2[shift],2) * pow( (eta_1-( jeta_1[shift]+jeta_2[shift] )*0.5), 2 ) );
            }
        }
    }
    jets->setCurrentUncertShift("",true);
    for(unsigned int shift = 0; shift<btagShifts.size(); ++shift )
    {
        jets->btagPromoteDemote( btagShifts[shift].second.first, btagShifts[shift].second.second );
        nbtag[shift]= jets->getNBtag();
        if (nbtag[shift] >= 1)
        {
            bpt_1[shift]=   jets->getBtagJet(0).Pt();
            beta_1[shift]=  jets->getBtagJet(0).Eta();
            bphi_1[shift]=  jets->getBtagJet(0).Phi();
            brawf_1[shift]= jets->getBtagJet(0).getProperty(PropertyEnum::rawFactor);
            //bmva_1[shift]=  jets->getBtagJet(0).getProperty(PropertyEnum::puIdDisc);
            bcsv_1[shift]=  jets->getBtagJet(0).getProperty(PropertyEnum::btagDeepB);
        }

        if (nbtag[shift] >= 2)
        {
            bpt_2[shift]=   jets->getBtagJet(1).Pt();
            beta_2[shift]=  jets->getBtagJet(1).Eta();
            bphi_2[shift]=  jets->getBtagJet(1).Phi();
            brawf_2[shift]= jets->getBtagJet(1).getProperty(PropertyEnum::rawFactor);
            //bmva_2[shift]=  jets->getBtagJet(1).getProperty(PropertyEnum::puIdDisc);
            bcsv_2[shift]=  jets->getBtagJet(1).getProperty(PropertyEnum::btagDeepB);
        }


        if (pdg1==15) gen_match_jetId_1=getGenMatch_jetId(leg1P4,jets);
        if (pdg2==15) gen_match_jetId_2=getGenMatch_jetId(leg2P4,jets);
    } 

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillPairBranches(HTTPair *pair, HTTJetCollection *jets)
{

    metcov00=pair->getMETMatrix().at(0);
    metcov01=pair->getMETMatrix().at(1);
    metcov10=pair->getMETMatrix().at(2);
    metcov11=pair->getMETMatrix().at(3);
    

    for(unsigned int shift = 0; shift<metShifts.size(); ++shift )
    {
        pair->setCurrentMETShift( metShifts[shift] );

        met[shift]   =pair->getMET().Mod();
        met_ex[shift]=pair->getMET().X();
        met_ey[shift]=pair->getMET().Y();
        metphi[shift]= TVector2::Phi_mpi_pi( pair->getMET().Phi() );

        mt_1[shift]=pair->getMTLeg1();
        mt_2[shift]=pair->getMTLeg2();

        m_sv[shift]=pair->getP4(HTTPair::SVFit).M();
        pt_sv[shift]=pair->getP4(HTTPair::SVFit).Pt();

        pt_tt[shift]=pair->getP4( HTTPair::Simple ).Pt();
        mt_tot[shift]=pair->getMTTOT();
        pt_sum[shift]=pt_1+pt_2+met[shift];
    }

    for(unsigned int shift = 0; shift<jecShifts.size(); ++shift )
    {
        jets->setCurrentUncertShift( jecShifts[shift].second.first, jecShifts[shift].second.second );
        htxs_reco_ggf[shift] = (int)getStage1Category(HiggsProdMode::GGF, pair->getP4( HTTPair::Simple ), jets);
        htxs_reco_vbf[shift] = (int)getStage1Category(HiggsProdMode::VBF, pair->getP4( HTTPair::Simple ), jets);
        //////////////////////////////////////////////////////////////////
        TLorentzVector ttjj; ttjj.SetPtEtaPhiM(-10,-10,-10,-10);
        if ( njetspt20[0] >=2 )
        {
            ttjj = pair->getP4( HTTPair::Simple ) + jets->getJet(0).P4() + jets->getJet(1).P4();
        }
        pt_ttjj[shift] = ttjj.Pt();
        m_ttjj[shift] = ttjj.M();
        //////////////////////////////////////////////////////////////////
    }
    jets->setCurrentUncertShift( "", true );
    pair->setCurrentMETShift( "" );
    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    pt_vis=pair->getPTVis();
    m_vis=pair->getMVis();
    dphi=leg1P4.DeltaPhi(leg2P4);
    dr_leptau=leg1P4.DeltaR(leg2P4);    
    //////////////////////////////////////////////////////////////////
    double zetaX = TMath::Cos(phi_1) + TMath::Cos(phi_2);
    double zetaY = TMath::Sin(phi_1) + TMath::Sin(phi_2);
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) {
      zetaX /= zetaR;
      zetaY /= zetaR;
    }
    pzetavis =  (leg1P4.Px() + leg2P4.Px())*zetaX + (leg1P4.Py() + leg2P4.Py())*zetaY;
    pzetamiss = met_ex[0]*zetaX + met_ey[0]*zetaY;
    dzeta = pzetamiss + (pzetavis - 1.85 * pzetavis);
    //////////////////////////////////////////////////////////////////
    double x1 = pt_1 / ( pt_1 + met[0] );
    double x2 = pt_2 / ( pt_2 + met[0] );
    if( TMath::Cos( dphi ) > -0.95
        && ( x1 > 0 && x1 < 1)
        && ( x2 > 0 && x2 < 1) )
    {
        m_coll =  m_vis / sqrt( x1 * x2 ) ;
    }
    else m_coll = DEF;
    //////////////////////////////////////////////////////////////////
    vector<TLorentzVector> objs;
    objs.push_back(leg1P4);
    objs.push_back(leg2P4);
    if ( njetspt20[0]>0 ) objs.push_back(jets->getJet(0).P4());
    if ( njetspt20[0]>1 ) objs.push_back(jets->getJet(1).P4());
    sphericity=calcSphericity(objs);  
    //////////////////////////////////////////////////////////////////
    TLorentzVector v0=leg1P4*( 1/sqrt(  pow(leg1P4.Px(),2)+pow(leg1P4.Py(),2)  ) ); //lep, normalized in transverse plane
    TLorentzVector v1=leg2P4*( 1/sqrt(  pow(leg2P4.Px(),2)+pow(leg2P4.Py(),2)  ) ); //tau, normalized in transverse plane
    float omega=v1.DeltaPhi(v0);
    float theta=-v0.Phi();
    float x=(     met_ex[0] * TMath::Sin(omega-theta)  - met_ey[0]*TMath::Cos(omega-theta)   ) / TMath::Sin(omega); //x coord in lep-tau system
    float y=(     met_ex[0] * TMath::Sin(theta)        + met_ey[0]*TMath::Cos(theta)         ) / TMath::Sin(omega); //y coord in lep-tau system
    met_centrality=( x+y ) / sqrt(x*x + y*y);
    //////////////////////////////////////////////////////////////////
 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillAdditionalLeptons( std::vector<HTTParticle> leptons, HTTPair *pair )
{
    unsigned int indexLeg1 = (unsigned int)pair->getIndexLeg1();
    unsigned int indexLeg2 = (unsigned int)pair->getIndexLeg2();

    UChar_t bitmask;

    for(unsigned int i = 0; i<leptons.size(); i++)
    {
        if( i == indexLeg1 || i == indexLeg2 || !leptons[i].isAdditionalLepton() ) continue;

        int lepPDGId = leptons[i].getPDGid();
        if( std::abs(lepPDGId) == 15 &&  channel == HTTAnalysis::TauTau ) continue;

        addlepton_p4.push_back( leptons[i].getP4() );
        addlepton_pt.push_back( leptons[i].getP4().Pt() );
        addlepton_eta.push_back( leptons[i].getP4().Eta() );
        addlepton_phi.push_back( leptons[i].getP4().Phi() );
        addlepton_m.push_back( leptons[i].getP4().M() );        
        addlepton_pdgId.push_back( leptons[i].getPDGid()*(-1) );
        addlepton_mc_match.push_back( leptons[i].getProperty(PropertyEnum::mc_match) );

        addlepton_tauAntiEle.push_back( (int)leptons[i].getProperty(PropertyEnum::idAntiEle) );
        addlepton_tauAntiMu.push_back( (int)leptons[i].getProperty(PropertyEnum::idAntiMu)  );

        addlepton_d0.push_back( leptons[i].getProperty(PropertyEnum::dxy) );
        addlepton_dZ.push_back( leptons[i].getProperty(PropertyEnum::dz) );
        addlepton_mt.push_back( leptons[i].getMT( pair->getMET() ) );

        if( std::abs(lepPDGId) == 13 )
        {
            nadditionalMu++;
            addmuon_pt.push_back( leptons[i].getP4().Pt() );
            addmuon_eta.push_back( leptons[i].getP4().Eta() );
            addmuon_phi.push_back( leptons[i].getP4().Phi() );
            addmuon_m.push_back( leptons[i].getP4().M() );
            addmuon_q.push_back(leptons[i].getProperty(PropertyEnum::charge));
            addmuon_iso.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("muonIsolation") ) );
            addmuon_gen_match.push_back(leptons[i].getProperty(PropertyEnum::mc_match) );

            addlepton_isoDeepTau.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("tauIsolationDeepVSmu") ) );
            addlepton_tauIDDeepTau.push_back(  (int)leptons[i].getProperty(HTTEvent::usePropertyFor.at("tauIDDeepVSmu")) );


            addlepton_iso.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("muonIsolation") ) );
            addlepton_mvis.push_back( ( leg2P4 + leptons[i].getP4() ).M() );

            addlepton_tauID.push_back( 0. );
            addlepton_tauDM.push_back( -1. );
            addlepton_tauCombIso.push_back( 0. );

        }else if( std::abs(lepPDGId) == 11 )
        {
            nadditionalEle++;
            addele_pt.push_back( leptons[i].getP4().Pt() );
            addele_eta.push_back( leptons[i].getP4().Eta() );
            addele_phi.push_back( leptons[i].getP4().Phi() );
            addele_m.push_back( leptons[i].getP4().M() );
            addele_q.push_back(leptons[i].getProperty(PropertyEnum::charge));
            addele_iso.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("electronIsolation") ) );
            addele_gen_match.push_back(leptons[i].getProperty(PropertyEnum::mc_match) );

            addlepton_isoDeepTau.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("tauIsolationDeepVSele") ) );
            addlepton_tauIDDeepTau.push_back(  (int)leptons[i].getProperty(HTTEvent::usePropertyFor.at("tauIDDeepVSele")) );
            

            addlepton_iso.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("electronIsolation") ) );
            addlepton_mvis.push_back( ( leg2P4 + leptons[i].getP4() ).M() );

            addlepton_tauID.push_back( 0. );
            addlepton_tauDM.push_back( -1. );
            addlepton_tauCombIso.push_back( 0. );


        }else if( std::abs(lepPDGId) == 15 )
        {
            nadditionalTau++;
            addtau_pt.push_back( leptons[i].getP4().Pt() );
            addtau_eta.push_back( leptons[i].getP4().Eta() );
            addtau_phi.push_back( leptons[i].getP4().Phi() );
            addtau_m.push_back( leptons[i].getP4().M() );
            addtau_q.push_back(leptons[i].getProperty(PropertyEnum::charge));
            addtau_byIsolationMVArun2v1DBoldDMwLTraw.push_back(leptons[i].getProperty( HTTEvent::usePropertyFor.at("tauIsolation") ));
            addtau_gen_match.push_back( leptons[i].getProperty(PropertyEnum::mc_match) );

            bitmask = leptons[i].getProperty(PropertyEnum::idAntiEle);
            bool antiEle = ( bitmask & 0x8) > 0;
            bitmask = leptons[i].getProperty(PropertyEnum::idAntiMu);
            bool antiMu  = (bitmask & 0x1 ) > 0;
            addtau_passesTauLepVetos.push_back( antiEle & antiMu );

            addtau_decayMode.push_back( leptons[i].getProperty(PropertyEnum::decayMode) );
            addtau_d0.push_back( leptons[i].getProperty(PropertyEnum::dxy) );
            addtau_dZ.push_back( leptons[i].getProperty(PropertyEnum::dz) );

            addtau_mt.push_back( leptons[i].getMT( pair->getMET() ) );
            addtau_mvis.push_back( ( leg1P4 + leptons[i].getP4() ).M() );

            addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits.push_back( leptons[i].getProperty(PropertyEnum::rawIso) );

            bitmask=leptons[i].getProperty(HTTEvent::usePropertyFor.at("tauID"));
            addtau_byVVLooseIsolationMVArun2v1DBoldDMwLT.push_back( (bitmask & 0x1)>0 );
            addtau_byVLooseIsolationMVArun2v1DBoldDMwLT.push_back( (bitmask & 0x2)>0 );
            addtau_byLooseIsolationMVArun2v1DBoldDMwLT.push_back( (bitmask & 0x4)>0 );
            addtau_byMediumIsolationMVArun2v1DBoldDMwLT.push_back( (bitmask & 0x8)>0 );
            addtau_byTightIsolationMVArun2v1DBoldDMwLT.push_back( (bitmask & 0x10)>0 );
            addtau_byVTightIsolationMVArun2v1DBoldDMwLT.push_back( (bitmask & 0x20)>0 );


            // Tau ID based on DeepTauIDv2 against jets
            bitmask=leptons[i].getProperty(HTTEvent::usePropertyFor.at("tauIDDeepVSjet"));    // Bitmask
    
            addtau_byVVVLooseIsolationDeepTau2017v2VSjet.push_back(  (bitmask & 0x1)>0 );
            addtau_byVVLooseIsolationDeepTau2017v2VSjet.push_back( (bitmask & 0x2)>0 );
            addtau_byVLooseIsolationDeepTau2017v2VSjet.push_back(  (bitmask & 0x4)>0 );
            addtau_byLooseIsolationDeepTau2017v2VSjet.push_back( (bitmask & 0x8)>0 );
            addtau_byMediumIsolationDeepTau2017v2VSjet.push_back(  (bitmask & 0x10)>0 );
            addtau_byTightIsolationDeepTau2017v2VSjet.push_back( (bitmask & 0x20)>0 );
            addtau_byVTightIsolationDeepTau2017v2VSjet.push_back( (bitmask & 0x40)>0 ); // 4*16 = 64   
            addtau_byVVTightIsolationDeepTau2017v2VSjet.push_back( (bitmask & 0x80)>0 ); //8*16     = 128

            addlepton_isoDeepTau.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("tauIsolationDeepVSjet") ) );
            addlepton_tauIDDeepTau.push_back(  (int)leptons[i].getProperty(HTTEvent::usePropertyFor.at("tauIDDeepVSjet")) );
                

            addlepton_iso.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("tauIsolation") ) );
            addlepton_mvis.push_back( ( leg1P4 + leptons[i].getP4() ).M() );

            addlepton_tauID.push_back(  (int)leptons[i].getProperty(HTTEvent::usePropertyFor.at("tauID")) );
            addlepton_tauDM.push_back(  (int)leptons[i].getProperty(PropertyEnum::decayMode) );
            addlepton_tauCombIso.push_back( leptons[i].getProperty(PropertyEnum::rawIso) );


        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillMELA(HTTJetCollection *jets)
{
  if(jets->getNJets(30) < 2) return;

  TLorentzVector tau1, tau2;
  tau1 = leg1P4;
  tau2 = leg2P4;

  float ch_1,ch_2;
  ch_1 = leg1.getCharge();
  ch_2 = leg2.getCharge();

  // Sanitize charge for application on same-sign events
  if (ch_1 * ch_2 > 0) {
    ch_2 = -ch_1;
  }

  TLorentzVector jet1, jet2;
  jet1 = jets->getJet(0).P4();
  jet2 = jets->getJet(1).P4();


  // Run MELA
  SimpleParticleCollection_t daughters;
  daughters.push_back(SimpleParticle_t(15 * ch_1, tau1));
  daughters.push_back(SimpleParticle_t(15 * ch_2, tau2));

  SimpleParticleCollection_t associated;
  associated.push_back(SimpleParticle_t(0, jet1));
  associated.push_back(SimpleParticle_t(0, jet2));

  SimpleParticleCollection_t associated2;
  associated2.push_back(SimpleParticle_t(0, jet2));
  associated2.push_back(SimpleParticle_t(0, jet1));

  mela->resetInputEvent();
  mela->setCandidateDecayMode(TVar::CandidateDecay_ff);
  mela->setInputEvent(&daughters, &associated, (SimpleParticleCollection_t *)0, false);

  // Hypothesis: SM Higgs
  mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
  mela->computeProdP(ME_vbf, false);
  mela->computeVBFAngles(ME_q2v1, ME_q2v2, ME_costheta1, ME_costheta2, ME_phi, ME_costhetastar, ME_phi1);

  // Hypothesis: Z + 2 jets
  // Compute the Hypothesis with flipped jets and sum them up for the discriminator.
  mela->setProcess(TVar::bkgZJets, TVar::MCFM, TVar::JJQCD);
  mela->computeProdP(ME_z2j_1, false);

  mela->resetInputEvent();
  mela->setInputEvent(&daughters, &associated2, (SimpleParticleCollection_t *)0, false);
  mela->computeProdP(ME_z2j_2, false);

  // Compute discriminator
  ME_D = ME_vbf / (ME_vbf + ME_z2j_1 + ME_z2j_2);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double EventWriter::calcSphericity(std::vector<TLorentzVector> p){

  TMatrixD S(3,3);

  std::vector<std::vector<double> > pvec;

  float denom=0;
  for (unsigned k=0; k<p.size(); k++){
    double dtmp[3]={p.at(k).Px(),p.at(k).Py(),p.at(k).Pz()};
    std::vector<double> vtmp(dtmp, dtmp+3);
    pvec.push_back(vtmp);
    denom+=dtmp[0]*dtmp[0]+dtmp[1]*dtmp[1]+dtmp[2]*dtmp[2];
  }

  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      float num=0;
      for (unsigned k=0; k<pvec.size(); k++){
        num+=pvec.at(k)[i]*pvec.at(k)[j];
      }
      S(i,j)=num/denom;
    }
  }
  return calcSphericityFromMatrix(S);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double EventWriter::calcSphericityFromMatrix(TMatrixD M) {

  //  TMatrixD M(3,3);
  //  M.SetMatrixArray(A);

  TMatrixDEigen V(M);
  TMatrixD Eval = V.GetEigenValues();

  double e1 = TMatrixDRow(Eval,0)(0);
  double e2 = TMatrixDRow(Eval,1)(1);
  double e3 = TMatrixDRow(Eval,2)(2);

  std::vector<double> evalvector;
  evalvector.push_back(e1);
  evalvector.push_back(e2);
  evalvector.push_back(e3);

  //sort eigenvalues to get the lowest two for use in sphericity (lowest to highest order)
  sort (evalvector.begin(), evalvector.end());

  //error-checking
  //this number should equal zero as the off-diagonal elements should all be zero in the eigenvalue matrix
  //returns error value of -1
  double check = TMatrixDRow(Eval,0)(1) + TMatrixDRow(Eval,0)(2) + TMatrixDRow(Eval,1)(0) + TMatrixDRow(Eval,1)(2) + TMatrixDRow(Eval,2)(0) + TMatrixDRow(Eval,2)(1);
  if (check != 0.0) {double err = -1; return err;}

  //for formula, see: http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node213.html
  
  double value = evalvector.at(0)+evalvector.at(1);
  double spher = 1.5*value;

  return spher;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int EventWriter::getGenMatch_jetId(TLorentzVector selObj, HTTJetCollection *jets){
  float minDR=1;
  int whichjet=0;

  for (unsigned i=0; i<(unsigned int)jets->getNJets(20); i++){
    TLorentzVector p4=jets->getJet(i).P4();
    if(p4.Pt() > 20 && fabs(p4.Eta() ) < 4.7 ){
      float tmpDR = calcDR( selObj.Eta(), selObj.Phi(), p4.Eta(), p4.Phi() );
      if( tmpDR < minDR ){
    minDR = tmpDR;
    whichjet=i;
      }
    }
  }

  if( minDR < 0.5 ) return jets->getJet(whichjet).getProperty(PropertyEnum::partonFlavour);
  return -99;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double EventWriter::calcDR(double eta1, double phi1, double eta2, double phi2){
  double deta = eta1-eta2;
  double dphi = TVector2::Phi_mpi_pi(phi1-phi2);
  return TMath::Sqrt( deta*deta+dphi*dphi );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::setDefault(){

    lumiWeight=DEFWEIGHT;
    runID=DEF;
    lumiBlock=DEF;
    eventNr=0; //unsigned
    entry=DEF;
    fileEntry=DEF;
    htxs_stage1cat=DEF;

    npv=DEF;
    npu=DEF;
    rho=DEF;
    gen_match_1=DEF;
    gen_match_2=DEF;
    gen_match_jetId_1=DEF;
    gen_match_jetId_2=DEF;
    NUP=DEF;
    //evtWeight=DEF;

    singleTriggerSFLeg1 = DEFWEIGHT;
    singleTriggerSFLeg2 = DEFWEIGHT;
    xTriggerSFLeg1 = DEFWEIGHT;
    xTriggerSFLeg2 = DEFWEIGHT;

    isoWeight_1 = DEFWEIGHT;
    isoWeight_2 = DEFWEIGHT;
    idWeight_1 = DEFWEIGHT;
    idWeight_2 = DEFWEIGHT;

    weight=DEFWEIGHT;
    puWeight=DEFWEIGHT;
    genWeight=DEFWEIGHT;
    idisoweight_1=DEFWEIGHT;
    idisoweight_2=DEFWEIGHT;
    sf_trk=DEFWEIGHT;
    sf_reco=DEFWEIGHT;
    sf_SingleOrCrossTrigger = DEFWEIGHT;
    sf_SingleXorCrossTrigger = DEFWEIGHT;
    sf_SingleTrigger = DEFWEIGHT;
    sf_DoubleTauTight = DEFWEIGHT;
    sf_DoubleTauVTight = DEFWEIGHT;
    trk_sf = DEFWEIGHT;
    trigweight_1 = DEFWEIGHT;
    trigweight_2 = DEFWEIGHT;
    trigweightVTight_1 = DEFWEIGHT;
    trigweightVTight_2 = DEFWEIGHT;      

    idisoweight_1 = DEFWEIGHT;
    idisoweight_2 = DEFWEIGHT;
    stitchedWeight=DEFWEIGHT;
    topPtReweightWeightRun2=DEFWEIGHT;
    topPtReweightWeightRun1=DEFWEIGHT;
    zPtReweightWeight=DEFWEIGHT;
    zPtReweightWeight1D=DEFWEIGHT;

    NNLO_ggH_weight=DEFWEIGHT;
    THU_ggH_Mu = DEFWEIGHT;
    THU_ggH_Res = DEFWEIGHT;
    THU_ggH_Mig01 = DEFWEIGHT;
    THU_ggH_Mig12 = DEFWEIGHT;
    THU_ggH_VBF2j = DEFWEIGHT;
    THU_ggH_VBF3j = DEFWEIGHT;
    THU_ggH_PT60 = DEFWEIGHT;
    THU_ggH_PT120 = DEFWEIGHT;
    THU_ggH_qmtop = DEFWEIGHT;
    
    eleTauFakeRateWeight=DEFWEIGHT;
    muTauFakeRateWeight=DEFWEIGHT;
    antilep_tauscaling=DEFWEIGHT;
//////////////////////////////////////////////////////////////////
    gen_Mll=DEF;
    gen_ll_px=DEF;
    gen_ll_py=DEF;
    gen_ll_pz=DEF;
    gen_vis_Mll=DEF;
    gen_ll_vis_px=DEF;
    gen_ll_vis_py=DEF;
    gen_ll_vis_pz=DEF;
    genJets=DEF;
    //////////////////////////////////////////////////////////////////  
    trg_singlemu_22=DEFFLAG;
    trg_crossmu_mu19tau20=DEFFLAG;
    trg_doubletau_35=DEFFLAG;
    trg_singlemuon_22=DEFFLAG;
    trg_singlemuon_22_eta2p1=DEFFLAG;
    trg_singlemuonTk_22=DEFFLAG;
    trg_singlemuonTk_22_eta2p1=DEFFLAG;
    trg_crossmuon_mu19tau20=DEFFLAG;
    trg_crossmuon_mu19tau20_singleL1=DEFFLAG;
    trg_singleelectron_25_eta2p1=DEFFLAG;
    trg_doubletau_35_mediso_eta2p1=DEFFLAG;
    trg_doubletau_35_medCombiso_eta2p1=DEFFLAG;
    Flag_goodVertices = DEFFLAG;
    Flag_globalTightHalo2016Filter = DEFFLAG;
    Flag_globalSuperTightHalo2016Filter = DEFFLAG;
    Flag_HBHENoiseFilter = DEFFLAG;
    Flag_HBHENoiseIsoFilter = DEFFLAG;
    Flag_EcalDeadCellTriggerPrimitiveFilter = DEFFLAG;
    Flag_BadPFMuonFilter = DEFFLAG;
    Flag_BadChargedCandidateFilter = DEFFLAG;
    Flag_eeBadScFilter = DEFFLAG;
    Flag_ecalBadCalibFilter = DEFFLAG;
    Flag_METFilters = DEFFLAG;

    //////////////////////////////////////////////////////////////////  
    pt_1=DEF;
    phi_1=DEF;
    eta_1=DEF;
    eta_SC_1=DEF;
    m_1=DEF;
    q_1=DEF;
    d0_1=DEF;
    dZ_1=DEF;
    iso_1=DEF;
    againstElectronLooseMVA6_1=DEF;
    againstElectronMediumMVA6_1=DEF;
    againstElectronTightMVA6_1=DEF;
    againstElectronVLooseMVA6_1=DEF;
    againstElectronVTightMVA6_1=DEF;
    againstMuonLoose3_1=DEF;
    againstMuonTight3_1=DEF;
    byCombinedIsolationDeltaBetaCorrRaw3Hits_1=DEF;
    
    byIsolationMVA3newDMwoLTraw_1=DEF;
    byIsolationMVA3oldDMwoLTraw_1=DEF;
    byIsolationMVA3newDMwLTraw_1=DEF;
    byIsolationMVA3oldDMwLTraw_1=DEF;
    
    byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1=DEF;
    byLooseIsolationMVArun2017v2DBoldDMwLT2017_1=DEF;
    byMediumIsolationMVArun2017v2DBoldDMwLT2017_1=DEF;
    byTightIsolationMVArun2017v2DBoldDMwLT2017_1=DEF;
    byVTightIsolationMVArun2017v2DBoldDMwLT2017_1=DEF;
    byVVVLooseIsolationDeepTau2017v2VSjet_1=DEF;
    byVVLooseIsolationDeepTau2017v2VSjet_1=DEF; 
    byVLooseIsolationDeepTau2017v2VSjet_1=DEF;  
    byLooseIsolationDeepTau2017v2VSjet_1=DEF; 
    byMediumIsolationDeepTau2017v2VSjet_1=DEF;  
    byTightIsolationDeepTau2017v2VSjet_1=DEF; 
    byVTightIsolationDeepTau2017v2VSjet_1=DEF; 
    byVVTightIsolationDeepTau2017v2VSjet_1=DEF; 
    byVVVLooseIsolationDeepTau2017v2VSe_1=DEF;
    byVVLooseIsolationDeepTau2017v2VSe_1=DEF; 
    byVLooseIsolationDeepTau2017v2VSe_1=DEF;  
    byLooseIsolationDeepTau2017v2VSe_1=DEF; 
    byMediumIsolationDeepTau2017v2VSe_1=DEF;  
    byTightIsolationDeepTau2017v2VSe_1=DEF; 
    byVTightIsolationDeepTau2017v2VSe_1=DEF; 
    byVVTightIsolationDeepTau2017v2VSe_1=DEF; 
    byVLooseIsolationDeepTau2017v2VSmu_1=DEF;  
    byLooseIsolationDeepTau2017v2VSmu_1=DEF; 
    byMediumIsolationDeepTau2017v2VSmu_1=DEF;  
    byTightIsolationDeepTau2017v2VSmu_1=DEF; 
    
    byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byLooseIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byMediumIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byTightIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byVTightIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byVVVLooseIsolationDeepTau2017v2VSjet_2=DEF;
    byVVLooseIsolationDeepTau2017v2VSjet_2=DEF; 
    byVLooseIsolationDeepTau2017v2VSjet_2=DEF;  
    byLooseIsolationDeepTau2017v2VSjet_2=DEF; 
    byMediumIsolationDeepTau2017v2VSjet_2=DEF;  
    byTightIsolationDeepTau2017v2VSjet_2=DEF; 
    byVTightIsolationDeepTau2017v2VSjet_2=DEF; 
    byVVTightIsolationDeepTau2017v2VSjet_2=DEF; 
    byVVVLooseIsolationDeepTau2017v2VSe_2=DEF;
    byVVLooseIsolationDeepTau2017v2VSe_2=DEF; 
    byVLooseIsolationDeepTau2017v2VSe_2=DEF;  
    byLooseIsolationDeepTau2017v2VSe_2=DEF; 
    byMediumIsolationDeepTau2017v2VSe_2=DEF;  
    byTightIsolationDeepTau2017v2VSe_2=DEF; 
    byVTightIsolationDeepTau2017v2VSe_2=DEF; 
    byVVTightIsolationDeepTau2017v2VSe_2=DEF; 
    byVLooseIsolationDeepTau2017v2VSmu_2=DEF;  
    byLooseIsolationDeepTau2017v2VSmu_2=DEF; 
    byMediumIsolationDeepTau2017v2VSmu_2=DEF;  
    byTightIsolationDeepTau2017v2VSmu_2=DEF; 
    
    chargedIsoPtSum_1=DEF;
    neutralIsoPtSum_1=DEF;
    puCorrPtSum_1=DEF;
    decayModeFindingOldDMs_1=DEF;
    decayModeFindingNewDMs_1=DEF;
    decayMode_1=DEF;
    id_e_mva_nt_loose_1=DEF;
    id_m_loose_1=DEF;
    id_m_medium_1=DEF;
    id_m_tight_1=DEF;
    id_m_highpt_1=DEF;
    id_e_cut_veto_1=DEF;
    id_e_cut_loose_1=DEF;
    id_e_cut_medium_1=DEF;
    id_e_cut_tight_1=DEF;
    
    //////////////////////////////////////////////////////////////////
    pt_2=DEF;
    phi_2=DEF;
    eta_2=DEF;
    m_2=DEF;
    q_2=DEF;
    d0_2=DEF;
    dZ_2=DEF;
    iso_2=DEF;
    againstElectronLooseMVA6_2=DEF;
    againstElectronMediumMVA6_2=DEF;
    againstElectronTightMVA6_2=DEF;
    againstElectronVLooseMVA6_2=DEF;
    againstElectronVTightMVA6_2=DEF;
    againstMuonLoose3_2=DEF;
    againstMuonTight3_2=DEF;
    byCombinedIsolationDeltaBetaCorrRaw3Hits_2=DEF;
    byIsolationMVA3newDMwoLTraw_2=DEF;
    byIsolationMVA3oldDMwoLTraw_2=DEF;
    byIsolationMVA3newDMwLTraw_2=DEF;
    byIsolationMVA3oldDMwLTraw_2=DEF;
    byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byLooseIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byMediumIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byTightIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byVTightIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    chargedIsoPtSum_2=DEF;
    neutralIsoPtSum_2=DEF;
    puCorrPtSum_2=DEF;
    decayModeFindingOldDMs_2=DEF;
    decayModeFindingNewDMs_2=DEF;
    decayMode_2=DEF;
    //////////////////////////////////////////////////////////////////
    for(unsigned int shift = 0; shift<jecShifts.size(); ++shift )
    {
        met[shift]=DEF;
        metphi[shift]=DEF;
        met_ex[shift]=DEF;
        met_ey[shift]=DEF;
        m_sv[shift]=DEF;
        pt_sv[shift]=DEF;
        pt_tt[shift]=DEF;
        pt_ttjj[shift]=DEF;
        m_ttjj[shift]=DEF;
        mt_1[shift]=DEF;
        mt_2[shift]=DEF;
        mt_tot[shift]=DEF;
        pt_sum[shift]=DEF;

        njets[shift]=DEF;
        njetspt20[shift]=DEF;
        njetingap[shift]=0;
        njetingap20[shift]=0;
        mjj[shift]=DEF;
        jdeta[shift]=DEF;
        dijetpt[shift]=DEF;
        dijetphi[shift]=DEF;
        jdphi[shift]=DEF;
        jpt_1[shift]=DEF;
        jpt_2[shift]=DEF;
        jeta_1[shift]=DEF;
        jeta_2[shift]=DEF;
        jphi_1[shift]=DEF;
        jphi_2[shift]=DEF;
        htxs_reco_ggf[shift]=0;
        htxs_reco_vbf[shift]=0;
    }

    jm_1=DEF;
    jm_2=DEF;
    jrawf_1=DEF;
    jrawf_2=DEF;
    jmva_1=DEF;
    jmva_2=DEF;
    jcsv_1=DEF;
    jcsv_2=DEF;

    for(int i =0 ; i < 5; ++i){
        nbtag[i]=DEF;
        bpt_1[i]=DEF;
        beta_1[i]=DEF;
        bphi_1[i]=DEF;
        brawf_1[i]=DEF;
        bmva_1[i]=DEF;
        bcsv_1[i]=DEF;
        bpt_2[i]=DEF;
        beta_2[i]=DEF;
        bphi_2[i]=DEF;
        brawf_2[i]=DEF;
        bmva_2[i]=DEF;
        bcsv_2[i]=DEF;
    }
    //////////////////////////////////////////////////////////////////
    ME_vbf = DEF;
    ME_q2v1 = DEF;
    ME_q2v2 = DEF;
    ME_costheta1 = DEF;
    ME_costheta2 = DEF;
    ME_phi = DEF;
    ME_costhetastar = DEF;
    ME_phi1 = DEF;
    ME_z2j_1 = DEF;
    ME_z2j_2 = DEF;
    ME_D = DEF;
    //////////////////////////////////////////////////////////////////  
    metcov00=DEF;
    metcov01=DEF;
    metcov10=DEF;
    metcov11=DEF;
    //////////////////////////////////////////////////////////////////
    passesTauLepVetos=DEF;
    passesThirdLepVeto=DEF;
    diMuonVeto=DEF;
    diElectronVeto=DEF;
    dilepton_veto=DEF;
    extramuon_veto=DEF;
    extraelec_veto=DEF;
    //////////////////////////////////////////////////////////////////
    pzetavis=DEF;
    pzetamiss=DEF;
    dzeta=DEF;
    //////////////////////////////////////////////////////////////////
    pt_vis=DEF;
    m_vis=DEF;
    m_coll=DEF;
    dphi=DEF;
    //////////////////////////////////////////////////////////////////
    dr_leptau=DEF;
    jeta1eta2=DEF;
    met_centrality=DEF;
    lep_etacentrality=DEF;
    sphericity=DEF;
    //////////////////////////////////////////////////////////////////

    addlepton_p4.clear();
    addlepton_pt.clear();
    addlepton_eta.clear();
    addlepton_phi.clear();
    addlepton_m.clear();    
    addlepton_iso.clear();
    addlepton_pdgId.clear();
    addlepton_mc_match.clear();
    addlepton_d0.clear();
    addlepton_dZ.clear();
    addlepton_mt.clear();
    addlepton_mvis.clear();
    addlepton_tauDM.clear();
    addlepton_tauCombIso.clear();
    addlepton_tauID.clear();
    addlepton_tauAntiEle.clear();
    addlepton_tauAntiMu.clear();
    addlepton_isoDeepTau.clear();
    addlepton_tauIDDeepTau.clear();

    //////////////////////////////////////////////////////////////////
    nadditionalMu = 0;
    addmuon_pt.clear();
    addmuon_eta.clear();
    addmuon_phi.clear();
    addmuon_m.clear();
    addmuon_q.clear();
    addmuon_iso.clear();
    addmuon_gen_match.clear();
    //////////////////////////////////////////////////////////////////
    nadditionalEle = 0;
    addele_pt.clear();
    addele_eta.clear();
    addele_phi.clear();
    addele_m.clear();
    addele_q.clear();
    addele_iso.clear();
    addele_gen_match.clear();
    //////////////////////////////////////////////////////////////////
    nadditionalTau = 0;
    addtau_pt.clear();
    addtau_eta.clear();
    addtau_phi.clear();
    addtau_m.clear();
    addtau_q.clear();
    addtau_byIsolationMVArun2v1DBoldDMwLTraw.clear();
    addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits.clear();
    addtau_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
    addtau_byTightCombinedIsolationDeltaBetaCorr3Hits.clear();
    addtau_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
    addtau_byVVLooseIsolationMVArun2v1DBoldDMwLT.clear();
    addtau_byVLooseIsolationMVArun2v1DBoldDMwLT.clear();
    addtau_byLooseIsolationMVArun2v1DBoldDMwLT.clear();
    addtau_byMediumIsolationMVArun2v1DBoldDMwLT.clear();
    addtau_byTightIsolationMVArun2v1DBoldDMwLT.clear();
    addtau_byVTightIsolationMVArun2v1DBoldDMwLT.clear();

    addtau_byVVVLooseIsolationDeepTau2017v2VSjet.clear();
    addtau_byVVLooseIsolationDeepTau2017v2VSjet.clear();
    addtau_byVLooseIsolationDeepTau2017v2VSjet.clear();
    addtau_byLooseIsolationDeepTau2017v2VSjet.clear();
    addtau_byMediumIsolationDeepTau2017v2VSjet.clear();
    addtau_byTightIsolationDeepTau2017v2VSjet.clear();
    addtau_byVTightIsolationDeepTau2017v2VSjet.clear();
    addtau_byVVTightIsolationDeepTau2017v2VSjet.clear();

    addtau_passesTauLepVetos.clear();
    addtau_decayMode.clear();
    addtau_d0.clear();
    addtau_dZ.clear();
    addtau_gen_match.clear();
    addtau_mt.clear();
    addtau_mvis.clear();
    //////////////////////////////////////////////////////////////////

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::initTree(TTree *t, vector< pair< string, pair<string,bool> > > jecShifts_,  bool isMC_, bool isSync_, vector<string> metShifts_ ){

    isMC=isMC_;
    isSync=isSync_;

    jecShifts = jecShifts_;
    metShifts = metShifts_;

    btagShifts.clear();
    btagShifts.push_back( make_pair( "",           make_pair("central","central") ) );
    if(isMC){
        btagShifts.push_back( make_pair( "MistagUp",   make_pair("up","central") ) );
        btagShifts.push_back( make_pair( "MistagDown", make_pair("down","central") ) );
        btagShifts.push_back( make_pair( "BtagUp",     make_pair("central","up") ) );
        btagShifts.push_back( make_pair( "BtagDown",   make_pair("central","down") ) );
    }

    const int erg_tev = 13;
    const float mPOLE = 125.6;
    TVar::VerbosityLevel verbosity = TVar::SILENT;
    mela = new Mela(erg_tev, mPOLE, verbosity);

    // used in 2017 but put inside the RooWorkspace below for 2018 and also for 2017 later
    // tauTrigSFTight = new TauTriggerSFs2017("utils/TauTriggerSFs2017/data/tauTriggerEfficiencies2017_New.root", "utils/TauTriggerSFs2017/data/tauTriggerEfficiencies2017.root","tight","MVA");
    // tauTrigSFVTight = new TauTriggerSFs2017("utils/TauTriggerSFs2017/data/tauTriggerEfficiencies2017_New.root", "utils/TauTriggerSFs2017/data/tauTriggerEfficiencies2017.root","vtight","MVA");

    TFile wsp("utils/CorrectionWorkspaces/htt_scalefactors_2016_v2.root");
    w = (RooWorkspace*)wsp.Get("w");



    if(!isSync){
        t->Branch("dr_leptau", &dr_leptau);

        t->Branch("jeta1eta2", &jeta1eta2);
        t->Branch("met_centrality", &met_centrality);
        t->Branch("lep_etacentrality", &lep_etacentrality);
        t->Branch("sphericity", &sphericity);

        t->Branch("addlepton_p4", &addlepton_p4);
        t->Branch("addlepton_pt", &addlepton_pt);
        t->Branch("addlepton_eta", &addlepton_eta);
        t->Branch("addlepton_phi", &addlepton_phi);
        t->Branch("addlepton_m", &addlepton_m);
        t->Branch("addlepton_iso", &addlepton_iso);
        t->Branch("addlepton_pdgId", &addlepton_pdgId);
        t->Branch("addlepton_mc_match", &addlepton_mc_match);
        t->Branch("addlepton_d0", &addlepton_d0);
        t->Branch("addlepton_dZ", &addlepton_dZ);
        t->Branch("addlepton_mt", &addlepton_mt);
        t->Branch("addlepton_mvis", &addlepton_mvis);
        t->Branch("addlepton_tauCombIso", &addlepton_tauCombIso);
        t->Branch("addlepton_tauID", &addlepton_tauID);
        t->Branch("addlepton_tauDM", &addlepton_tauDM);
        t->Branch("addlepton_tauAntiEle", &addlepton_tauAntiEle);
        t->Branch("addlepton_tauAntiMu", &addlepton_tauAntiMu);
        t->Branch("addlepton_isoDeepTau", &addlepton_isoDeepTau);
        t->Branch("addlepton_tauIDDeepTau", &addlepton_tauIDDeepTau);


        t->Branch("nadditionalMu", &nadditionalMu);
        t->Branch("addmuon_pt", &addmuon_pt);
        t->Branch("addmuon_eta", &addmuon_eta);
        t->Branch("addmuon_phi", &addmuon_phi);
        t->Branch("addmuon_m", &addmuon_m);
        t->Branch("addmuon_q", &addmuon_q);
        t->Branch("addmuon_iso", &addmuon_iso);
        t->Branch("addmuon_gen_match", &addmuon_gen_match);

        t->Branch("nadditionalEle", &nadditionalEle);
        t->Branch("addele_pt", &addele_pt);
        t->Branch("addele_eta", &addele_eta);
        t->Branch("addele_phi", &addele_phi);
        t->Branch("addele_m", &addele_m);
        t->Branch("addele_q", &addele_q);
        t->Branch("addele_iso", &addele_iso);
        t->Branch("addele_gen_match", &addele_gen_match);

        t->Branch("nadditionalTau", &nadditionalTau);
        t->Branch("addtau_pt", &addtau_pt);
        t->Branch("addtau_eta", &addtau_eta);
        t->Branch("addtau_phi", &addtau_phi);
        t->Branch("addtau_m", &addtau_m);
        t->Branch("addtau_q", &addtau_q);
        t->Branch("addtau_byIsolationMVArun2v1DBoldDMwLTraw", &addtau_byIsolationMVArun2v1DBoldDMwLTraw);
        t->Branch("addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits", &addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits);
        t->Branch("addtau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &addtau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
        t->Branch("addtau_byTightCombinedIsolationDeltaBetaCorr3Hits", &addtau_byTightCombinedIsolationDeltaBetaCorr3Hits);
        t->Branch("addtau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &addtau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
        t->Branch("addtau_byVVLooseIsolationMVArun2v1DBoldDMwLT", &addtau_byVVLooseIsolationMVArun2v1DBoldDMwLT);
        t->Branch("addtau_byVLooseIsolationMVArun2v1DBoldDMwLT", &addtau_byVLooseIsolationMVArun2v1DBoldDMwLT);
        t->Branch("addtau_byLooseIsolationMVArun2v1DBoldDMwLT", &addtau_byLooseIsolationMVArun2v1DBoldDMwLT);
        t->Branch("addtau_byMediumIsolationMVArun2v1DBoldDMwLT", &addtau_byMediumIsolationMVArun2v1DBoldDMwLT);
        t->Branch("addtau_byTightIsolationMVArun2v1DBoldDMwLT", &addtau_byTightIsolationMVArun2v1DBoldDMwLT);
        t->Branch("addtau_byVTightIsolationMVArun2v1DBoldDMwLT", &addtau_byVTightIsolationMVArun2v1DBoldDMwLT);
      
        t->Branch("addtau_byVVVLooseIsolationDeepTau2017v2VSjet", &addtau_byVVVLooseIsolationDeepTau2017v2VSjet);
        t->Branch("addtau_byVVLooseIsolationDeepTau2017v2VSjet", &addtau_byVVLooseIsolationDeepTau2017v2VSjet);
        t->Branch("addtau_byVLooseIsolationDeepTau2017v2VSjet", &addtau_byVLooseIsolationDeepTau2017v2VSjet);
        t->Branch("addtau_byLooseIsolationDeepTau2017v2VSjet", &addtau_byLooseIsolationDeepTau2017v2VSjet);
        t->Branch("addtau_byMediumIsolationDeepTau2017v2VSjet", &addtau_byMediumIsolationDeepTau2017v2VSjet);
        t->Branch("addtau_byTightIsolationDeepTau2017v2VSjet", &addtau_byTightIsolationDeepTau2017v2VSjet);
        t->Branch("addtau_byVTightIsolationDeepTau2017v2VSjet", &addtau_byVTightIsolationDeepTau2017v2VSjet);
        t->Branch("addtau_byVVTightIsolationDeepTau2017v2VSjet", &addtau_byVVTightIsolationDeepTau2017v2VSjet);

        t->Branch("addtau_passesTauLepVetos", &addtau_passesTauLepVetos);
        t->Branch("addtau_decayMode", &addtau_decayMode);
        t->Branch("addtau_d0", &addtau_d0);
        t->Branch("addtau_dZ", &addtau_dZ);
        t->Branch("addtau_gen_match", &addtau_gen_match);
        t->Branch("addtau_mt", &addtau_mt);
        t->Branch("addtau_mvis", &addtau_mvis);
    }

    for(unsigned int shift = 0; shift<metShifts.size(); ++shift )
    {
        t->Branch( ("met"+metShifts[shift]).c_str(),      &met[shift]);
        t->Branch( ("metphi"+metShifts[shift]).c_str(),   &metphi[shift]);
        t->Branch( ("met_ex"+metShifts[shift]).c_str(),   &met_ex[shift]);
        t->Branch( ("met_ey"+metShifts[shift]).c_str(),   &met_ey[shift]);
        t->Branch( ("m_sv"+metShifts[shift]).c_str(),     &m_sv[shift]);
        t->Branch( ("pt_sv"+metShifts[shift]).c_str(),    &pt_sv[shift]);
        t->Branch( ("pt_tt"+metShifts[shift]).c_str(),    &pt_tt[shift]);
        t->Branch( ("pt_ttjj"+metShifts[shift]).c_str(),  &pt_ttjj[shift]);
        t->Branch( ("m_ttjj"+metShifts[shift]).c_str(),   &m_ttjj[shift]);
        t->Branch( ("pt_sum"+metShifts[shift]).c_str(),   &pt_sum[shift]);
        t->Branch( ("mt_1"+metShifts[shift]).c_str(),     &mt_1[shift]);
        t->Branch( ("mt_2"+metShifts[shift]).c_str(),     &mt_2[shift]);
        t->Branch( ("mt_tot"+metShifts[shift]).c_str(),   &mt_tot[shift]);
    }

    for(unsigned int shift = 0; shift<jecShifts.size(); ++shift )
    {
        t->Branch( ("njets"+jecShifts[shift].first).c_str() ,      &njets[shift]);
        t->Branch( ("njetspt20"+jecShifts[shift].first).c_str(),   &njetspt20[shift]);
        t->Branch( ("njetingap"+jecShifts[shift].first).c_str(),   &njetingap[shift]);
        t->Branch( ("njetingap20"+jecShifts[shift].first).c_str(), &njetingap20[shift]);
        t->Branch( ("dijetpt"+jecShifts[shift].first).c_str(),     &dijetpt[shift]);
        t->Branch( ("dijetphi"+jecShifts[shift].first).c_str(),    &dijetphi[shift]);
        t->Branch( ("jdphi"+jecShifts[shift].first).c_str(),       &jdphi[shift]);
        t->Branch( ("jdeta"+jecShifts[shift].first).c_str(),       &jdeta[shift]);
        t->Branch( ("mjj"+jecShifts[shift].first).c_str(),         &mjj[shift]);

        t->Branch( ("jpt_1"+jecShifts[shift].first).c_str(),       &jpt_1[shift]);
        t->Branch( ("jpt_2"+jecShifts[shift].first).c_str(),       &jpt_2[shift]);
        t->Branch( ("jeta_1"+jecShifts[shift].first).c_str(),      &jeta_1[shift]);
        t->Branch( ("jeta_2"+jecShifts[shift].first).c_str(),      &jeta_2[shift]);
        t->Branch( ("jphi_1"+jecShifts[shift].first).c_str(),      &jphi_1[shift]);
        t->Branch( ("jphi_2"+jecShifts[shift].first).c_str(),      &jphi_2[shift]);

        t->Branch( ("htxs_reco_ggf"+jecShifts[shift].first).c_str(), &htxs_reco_ggf[shift]);
        t->Branch( ("htxs_reco_vbf"+jecShifts[shift].first).c_str(), &htxs_reco_vbf[shift]);
    }

    for(unsigned int shift=0; shift < btagShifts.size(); ++shift )
    {
        t->Branch( ("nbtag"+btagShifts[shift].first).c_str(),   &nbtag[shift]);
        t->Branch( ("bpt_1"+btagShifts[shift].first).c_str(),   &bpt_1[shift]);
        t->Branch( ("bpt_2"+btagShifts[shift].first).c_str(),   &bpt_2[shift]);
        t->Branch( ("beta_1"+btagShifts[shift].first).c_str(),  &beta_1[shift]);
        t->Branch( ("beta_2"+btagShifts[shift].first).c_str(),  &beta_2[shift]);
        t->Branch( ("bphi_1"+btagShifts[shift].first).c_str(),  &bphi_1[shift]);
        t->Branch( ("bphi_2"+btagShifts[shift].first).c_str(),  &bphi_2[shift]);
        t->Branch( ("brawf_1"+btagShifts[shift].first).c_str(), &brawf_1[shift]);
        t->Branch( ("brawf_2"+btagShifts[shift].first).c_str(), &brawf_2[shift]);
        t->Branch( ("bmva_1"+btagShifts[shift].first).c_str(),  &bmva_1[shift]);
        t->Branch( ("bmva_2"+btagShifts[shift].first).c_str(),  &bmva_2[shift]);
        t->Branch( ("bcsv_1"+btagShifts[shift].first).c_str(),  &bcsv_1[shift]);
        t->Branch( ("bcsv_2"+btagShifts[shift].first).c_str(),  &bcsv_2[shift]);
    }


    t->Branch("ME_vbf", &ME_vbf);
    t->Branch("ME_q2v1", &ME_q2v1);
    t->Branch("ME_q2v2", &ME_q2v2);
    t->Branch("ME_costheta1", &ME_costheta1);
    t->Branch("ME_costheta2", &ME_costheta2);
    t->Branch("ME_phi", &ME_phi);
    t->Branch("ME_costhetastar", &ME_costhetastar);
    t->Branch("ME_phi1", &ME_phi1);
    t->Branch("ME_z2j_1", &ME_z2j_1);
    t->Branch("ME_z2j_2", &ME_z2j_2);
    t->Branch("ME_D", &ME_D);

    t->Branch("gen_Mll", &gen_Mll);
    t->Branch("genpX", &gen_ll_px);
    t->Branch("genpY", &gen_ll_py);
    t->Branch("genpZ", &gen_ll_pz);
    t->Branch("gen_vis_Mll", &gen_vis_Mll);
    t->Branch("vispX", &gen_ll_vis_px);
    t->Branch("vispY", &gen_ll_vis_py);
    t->Branch("vispZ", &gen_ll_vis_pz);
    t->Branch("npv", &npv);
    t->Branch("npu", &npu);
    t->Branch("rho", &rho);
    t->Branch("NUP", &NUP);

    t->Branch("flagMETFilter", &flagMETFilter);
    t->Branch("Flag_goodVertices", &Flag_goodVertices);
    t->Branch("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter);
    t->Branch("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter);
    t->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter);
    t->Branch("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter);
    t->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter);
    t->Branch("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter);
    t->Branch("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter);
    t->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter);
    t->Branch("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter);
    t->Branch("Flag_METFilters", &Flag_METFilters);

    t->Branch("gen_match_1", &gen_match_1);
    t->Branch("gen_match_2", &gen_match_2);
    t->Branch("gen_match_jetId_1", &gen_match_jetId_1);
    t->Branch("gen_match_jetId_2", &gen_match_jetId_2);
    t->Branch("genJets", &genJets);
    t->Branch("pdg_1", &pdg1);
    t->Branch("pdg_2", &pdg2);

    t->Branch("pt_1", &pt_1);
    t->Branch("phi_1", &phi_1);
    t->Branch("eta_1", &eta_1);
    t->Branch("eta_SC_1", &eta_SC_1);
    t->Branch("m_1", &m_1);
    t->Branch("q_1", &q_1);
    t->Branch("d0_1", &d0_1);
    t->Branch("dZ_1", &dZ_1);
    t->Branch("iso_1", &iso_1);
    t->Branch("againstElectronMVA6_1", &againstElectronMVA6_1);
    t->Branch("againstElectronLooseMVA6_1", &againstElectronLooseMVA6_1);
    t->Branch("againstElectronMediumMVA6_1", &againstElectronMediumMVA6_1);
    t->Branch("againstElectronTightMVA6_1", &againstElectronTightMVA6_1);
    t->Branch("againstElectronVLooseMVA6_1", &againstElectronVLooseMVA6_1);
    t->Branch("againstElectronVTightMVA6_1", &againstElectronVTightMVA6_1);
    t->Branch("againstMuon3_1", &againstMuon3_1);    
    t->Branch("againstMuonLoose3_1", &againstMuonLoose3_1);
    t->Branch("againstMuonTight3_1", &againstMuonTight3_1);
    t->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_1", &byCombinedIsolationDeltaBetaCorrRaw3Hits_1);
    t->Branch("byIsolationMVA3newDMwoLTraw_1", &byIsolationMVA3newDMwoLTraw_1);
    t->Branch("byIsolationMVA3oldDMwoLTraw_1", &byIsolationMVA3oldDMwoLTraw_1);
    t->Branch("byIsolationMVA3newDMwLTraw_1", &byIsolationMVA3newDMwLTraw_1);
    t->Branch("byIsolationMVA3oldDMwLTraw_1", &byIsolationMVA3oldDMwLTraw_1);
    t->Branch("byIsolationMVArun2017v2DBoldDMwLTraw2017_1", &byIsolationMVArun2017v2DBoldDMwLTraw2017_1);
    t->Branch("byIsolationMVArun2017v2DBoldDMwLT2017_1", &byIsolationMVArun2017v2DBoldDMwLT2017_1);   
    t->Branch("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1);
    t->Branch("byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1);
    t->Branch("byLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &byLooseIsolationMVArun2017v2DBoldDMwLT2017_1);
    t->Branch("byMediumIsolationMVArun2017v2DBoldDMwLT2017_1", &byMediumIsolationMVArun2017v2DBoldDMwLT2017_1);
    t->Branch("byTightIsolationMVArun2017v2DBoldDMwLT2017_1", &byTightIsolationMVArun2017v2DBoldDMwLT2017_1);
    t->Branch("byVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &byVTightIsolationMVArun2017v2DBoldDMwLT2017_1);
    t->Branch("byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1);
    
    t->Branch("byIsolationDeepTau2017v2VSjet_raw_1", &byIsolationDeepTau2017v2VSjet_raw_1);
    t->Branch("byIsolationDeepTau2017v2VSjet_1", &byIsolationDeepTau2017v2VSjet_1);
    t->Branch("byVVVLooseIsolationDeepTau2017v2VSjet_1", &byVVVLooseIsolationDeepTau2017v2VSjet_1);
    t->Branch("byVVLooseIsolationDeepTau2017v2VSjet_1", &byVVLooseIsolationDeepTau2017v2VSjet_1);
    t->Branch("byVLooseIsolationDeepTau2017v2VSjet_1", &byVLooseIsolationDeepTau2017v2VSjet_1);
    t->Branch("byLooseIsolationDeepTau2017v2VSjet_1", &byLooseIsolationDeepTau2017v2VSjet_1);
    t->Branch("byMediumIsolationDeepTau2017v2VSjet_1", &byMediumIsolationDeepTau2017v2VSjet_1);
    t->Branch("byTightIsolationDeepTau2017v2VSjet_1", &byTightIsolationDeepTau2017v2VSjet_1);
    t->Branch("byVTightIsolationDeepTau2017v2VSjet_1", &byVTightIsolationDeepTau2017v2VSjet_1);
    t->Branch("byVVTightIsolationDeepTau2017v2VSjet_1", &byVVTightIsolationDeepTau2017v2VSjet_1);

    t->Branch("byIsolationDeepTau2017v2VSe_raw_1", &byIsolationDeepTau2017v2VSe_raw_1);
    t->Branch("byIsolationDeepTau2017v2VSe_1", &byIsolationDeepTau2017v2VSe_1);
    t->Branch("byVVVLooseIsolationDeepTau2017v2VSe_1", &byVVVLooseIsolationDeepTau2017v2VSe_1);
    t->Branch("byVVLooseIsolationDeepTau2017v2VSe_1", &byVVLooseIsolationDeepTau2017v2VSe_1);
    t->Branch("byVLooseIsolationDeepTau2017v2VSe_1", &byVLooseIsolationDeepTau2017v2VSe_1);
    t->Branch("byLooseIsolationDeepTau2017v2VSe_1", &byLooseIsolationDeepTau2017v2VSe_1);
    t->Branch("byMediumIsolationDeepTau2017v2VSe_1", &byMediumIsolationDeepTau2017v2VSe_1);
    t->Branch("byTightIsolationDeepTau2017v2VSe_1", &byTightIsolationDeepTau2017v2VSe_1);
    t->Branch("byVTightIsolationDeepTau2017v2VSe_1", &byVTightIsolationDeepTau2017v2VSe_1);
    t->Branch("byVVTightIsolationDeepTau2017v2VSe_1", &byVVTightIsolationDeepTau2017v2VSe_1);

    t->Branch("byIsolationDeepTau2017v2VSmu_raw_1", &byIsolationDeepTau2017v2VSmu_raw_1);
    t->Branch("byIsolationDeepTau2017v2VSmu_1", &byIsolationDeepTau2017v2VSmu_1);
    t->Branch("byVLooseIsolationDeepTau2017v2VSmu_1", &byVLooseIsolationDeepTau2017v2VSmu_1);
    t->Branch("byLooseIsolationDeepTau2017v2VSmu_1", &byLooseIsolationDeepTau2017v2VSmu_1);
    t->Branch("byMediumIsolationDeepTau2017v2VSmu_1", &byMediumIsolationDeepTau2017v2VSmu_1);
    t->Branch("byTightIsolationDeepTau2017v2VSmu_1", &byTightIsolationDeepTau2017v2VSmu_1);
    
    t->Branch("chargedIsoPtSum_1", &chargedIsoPtSum_1);
    t->Branch("neutralIsoPtSum_1", &neutralIsoPtSum_1);
    t->Branch("puCorrPtSum_1", &puCorrPtSum_1);
    t->Branch("decayModeFindingOldDMs_1", &decayModeFindingOldDMs_1);
    t->Branch("decayModeFindingNewDMs_1", &decayModeFindingNewDMs_1);
    t->Branch("decayMode_1", &decayMode_1);
    t->Branch("id_e_mva_nt_loose_1", &id_e_mva_nt_loose_1);

    t->Branch("id_m_loose_1", &id_m_loose_1);
    t->Branch("id_m_medium_1", &id_m_medium_1);
    t->Branch("id_m_tight_1", &id_m_tight_1);
    t->Branch("id_m_highpt_1", &id_m_highpt_1);
    t->Branch("id_e_cut_veto_1", &id_e_cut_veto_1);
    t->Branch("id_e_cut_loose_1", &id_e_cut_loose_1);
    t->Branch("id_e_cut_medium_1", &id_e_cut_medium_1);
    t->Branch("id_e_cut_tight_1", &id_e_cut_tight_1);

    
    t->Branch("pt_2", &pt_2);
    t->Branch("phi_2", &phi_2);
    t->Branch("eta_2", &eta_2); 
    t->Branch("m_2", &m_2);
    t->Branch("q_2", &q_2);
    t->Branch("d0_2", &d0_2);
    t->Branch("dZ_2", &dZ_2);

    t->Branch("iso_2", &iso_2);
    t->Branch("againstElectronMVA6_2", &againstElectronMVA6_2);
    t->Branch("againstElectronLooseMVA6_2", &againstElectronLooseMVA6_2);
    t->Branch("againstElectronMediumMVA6_2", &againstElectronMediumMVA6_2);
    t->Branch("againstElectronTightMVA6_2", &againstElectronTightMVA6_2);
    t->Branch("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2);
    t->Branch("againstElectronVTightMVA6_2", &againstElectronVTightMVA6_2);
    t->Branch("againstMuon3_2", &againstMuon3_2);    
    t->Branch("againstMuonLoose3_2", &againstMuonLoose3_2);
    t->Branch("againstMuonTight3_2", &againstMuonTight3_2);
    t->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
    t->Branch("byIsolationMVA3newDMwoLTraw_2", &byIsolationMVA3newDMwoLTraw_2);
    t->Branch("byIsolationMVA3oldDMwoLTraw_2", &byIsolationMVA3oldDMwoLTraw_2);
    t->Branch("byIsolationMVA3newDMwLTraw_2", &byIsolationMVA3newDMwLTraw_2);
    t->Branch("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2);
    t->Branch("byIsolationMVArun2017v2DBoldDMwLTraw2017_2", &byIsolationMVArun2017v2DBoldDMwLTraw2017_2);
    t->Branch("byIsolationMVArun2017v2DBoldDMwLT2017_2", &byIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &byLooseIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byMediumIsolationMVArun2017v2DBoldDMwLT2017_2", &byMediumIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byTightIsolationMVArun2017v2DBoldDMwLT2017_2", &byTightIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &byVTightIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2);

    t->Branch("byIsolationDeepTau2017v2VSjet_raw_2", &byIsolationDeepTau2017v2VSjet_raw_2);
    t->Branch("byIsolationDeepTau2017v2VSjet_2", &byIsolationDeepTau2017v2VSjet_2);
    t->Branch("byVVVLooseIsolationDeepTau2017v2VSjet_2", &byVVVLooseIsolationDeepTau2017v2VSjet_2);
    t->Branch("byVVLooseIsolationDeepTau2017v2VSjet_2", &byVVLooseIsolationDeepTau2017v2VSjet_2);
    t->Branch("byVLooseIsolationDeepTau2017v2VSjet_2", &byVLooseIsolationDeepTau2017v2VSjet_2);
    t->Branch("byLooseIsolationDeepTau2017v2VSjet_2", &byLooseIsolationDeepTau2017v2VSjet_2);
    t->Branch("byMediumIsolationDeepTau2017v2VSjet_2", &byMediumIsolationDeepTau2017v2VSjet_2);
    t->Branch("byTightIsolationDeepTau2017v2VSjet_2", &byTightIsolationDeepTau2017v2VSjet_2);
    t->Branch("byVTightIsolationDeepTau2017v2VSjet_2", &byVTightIsolationDeepTau2017v2VSjet_2);
    t->Branch("byVVTightIsolationDeepTau2017v2VSjet_2", &byVVTightIsolationDeepTau2017v2VSjet_2);

    t->Branch("byIsolationDeepTau2017v2VSe_raw_2", &byIsolationDeepTau2017v2VSe_raw_2);
    t->Branch("byIsolationDeepTau2017v2VSe_2", &byIsolationDeepTau2017v2VSe_2);
    t->Branch("byVVVLooseIsolationDeepTau2017v2VSe_2", &byVVVLooseIsolationDeepTau2017v2VSe_2);
    t->Branch("byVVLooseIsolationDeepTau2017v2VSe_2", &byVVLooseIsolationDeepTau2017v2VSe_2);
    t->Branch("byVLooseIsolationDeepTau2017v2VSe_2", &byVLooseIsolationDeepTau2017v2VSe_2);
    t->Branch("byLooseIsolationDeepTau2017v2VSe_2", &byLooseIsolationDeepTau2017v2VSe_2);
    t->Branch("byMediumIsolationDeepTau2017v2VSe_2", &byMediumIsolationDeepTau2017v2VSe_2);
    t->Branch("byTightIsolationDeepTau2017v2VSe_2", &byTightIsolationDeepTau2017v2VSe_2);
    t->Branch("byVTightIsolationDeepTau2017v2VSe_2", &byVTightIsolationDeepTau2017v2VSe_2);
    t->Branch("byVVTightIsolationDeepTau2017v2VSe_2", &byVVTightIsolationDeepTau2017v2VSe_2);

    t->Branch("byIsolationDeepTau2017v2VSmu_raw_2", &byIsolationDeepTau2017v2VSmu_raw_2);
    t->Branch("byIsolationDeepTau2017v2VSmu_2", &byIsolationDeepTau2017v2VSmu_2);
    t->Branch("byVLooseIsolationDeepTau2017v2VSmu_2", &byVLooseIsolationDeepTau2017v2VSmu_2);
    t->Branch("byLooseIsolationDeepTau2017v2VSmu_2", &byLooseIsolationDeepTau2017v2VSmu_2);
    t->Branch("byMediumIsolationDeepTau2017v2VSmu_2", &byMediumIsolationDeepTau2017v2VSmu_2);
    t->Branch("byTightIsolationDeepTau2017v2VSmu_2", &byTightIsolationDeepTau2017v2VSmu_2);

    t->Branch("chargedIsoPtSum_2", &chargedIsoPtSum_2);
    t->Branch("neutralIsoPtSum_2", &neutralIsoPtSum_2);
    t->Branch("puCorrPtSum_2", &puCorrPtSum_2);
    t->Branch("decayModeFindingOldDMs_2", &decayModeFindingOldDMs_2);
    t->Branch("decayModeFindingNewDMs_2", &decayModeFindingNewDMs_2);
    t->Branch("decayMode_2", &decayMode_2);

    t->Branch("pzetavis", &pzetavis);
    t->Branch("pzetamiss", &pzetamiss);
    t->Branch("dzeta", &dzeta);
    t->Branch("m_vis", &m_vis);
    t->Branch("m_coll", &m_coll);
    t->Branch("pt_vis", &pt_vis);
    t->Branch("dphi", &dphi);

    t->Branch("passesTauLepVetos", &passesTauLepVetos);
    t->Branch("passesThirdLepVeto", &passesThirdLepVeto);
    t->Branch("passesDiMuonVeto", &passesDiMuonVeto);
    t->Branch("passesDiElectronVeto", &passesDiElectronVeto);
    t->Branch("diMuonVeto", &diMuonVeto);
    t->Branch("diElectronVeto", &diElectronVeto);    
    t->Branch("dilepton_veto", &dilepton_veto);
    t->Branch("extraelec_veto", &extraelec_veto);
    t->Branch("extramuon_veto", &extramuon_veto);

    t->Branch("metcov00", &metcov00);
    t->Branch("metcov01", &metcov01);
    t->Branch("metcov10", &metcov10);
    t->Branch("metcov11", &metcov11);
    t->Branch("htxs_stage1cat", &htxs_stage1cat);

    t->Branch("jm_1", &jm_1);
    t->Branch("jm_2", &jm_2);
    t->Branch("jrawf_1", &jrawf_1);
    t->Branch("jrawf_2", &jrawf_2);    
    t->Branch("jmva_1", &jmva_1);
    t->Branch("jmva_2",&jmva_2);    
    t->Branch("jcsv_1", &jcsv_1);
    t->Branch("jcsv_2",&jcsv_2);    

    t->Branch("weight", &weight);
    t->Branch("eventWeight", &weight);
    t->Branch("lumiWeight", &lumiWeight);
    t->Branch("puweight", &puWeight);
    t->Branch("genweight", &genWeight);
    t->Branch("xsec", &xsec);
    t->Branch("genNEventsWeight", &genNEventsWeight);

    t->Branch("singleTriggerSFLeg1",&singleTriggerSFLeg1);
    t->Branch("singleTriggerSFLeg2",&singleTriggerSFLeg2);
    t->Branch("xTriggerSFLeg1",&xTriggerSFLeg1);
    t->Branch("xTriggerSFLeg2",&xTriggerSFLeg2);

    t->Branch("isoWeight_1",&isoWeight_1);
    t->Branch("isoWeight_2",&isoWeight_2);
    t->Branch("idWeight_1",&idWeight_1);
    t->Branch("idWeight_2",&idWeight_2);

    t->Branch("idisoweight_1", &idisoweight_1);
    t->Branch("idisoweight_2", &idisoweight_2);
    t->Branch("sf_trk", &sf_trk);
    t->Branch("sf_reco", &sf_reco);
    t->Branch("sf_SingleOrCrossTrigger", &sf_SingleOrCrossTrigger);
    t->Branch("sf_SingleXorCrossTrigger", &sf_SingleXorCrossTrigger);
    t->Branch("sf_SingleTrigger", &sf_SingleTrigger);

    t->Branch("trk_sf", &trk_sf);
    t->Branch("trigweight_1", &trigweight_1);
    t->Branch("trigweight_2", &trigweight_2);
    t->Branch("trigweightVTight_1", &trigweightVTight_1);
    t->Branch("trigweightVTight_2", &trigweightVTight_2);
    t->Branch("idisoweight_1", &idisoweight_1);
    t->Branch("idisoweight_2", &idisoweight_2);
   

    t->Branch("sf_DoubleTauTight", &sf_DoubleTauTight);
    t->Branch("sf_DoubleTauVTight", &sf_DoubleTauVTight);

    t->Branch("stitchedWeight", &stitchedWeight);
    t->Branch("topPtReweightWeightRun2", &topPtReweightWeightRun2);
    t->Branch("topPtReweightWeightRun1", &topPtReweightWeightRun1);
    t->Branch("zPtReweightWeight", &zPtReweightWeight);
    t->Branch("zPtReweightWeight1D", &zPtReweightWeight1D);
    t->Branch("eleTauFakeRateWeight", &eleTauFakeRateWeight);
    t->Branch("muTauFakeRateWeight", &muTauFakeRateWeight);
    t->Branch("antilep_tauscaling", &antilep_tauscaling);

    t->Branch("NNLO_ggH_weight", &NNLO_ggH_weight);

    t->Branch("THU_ggH_Mu", &THU_ggH_Mu);
    t->Branch("THU_ggH_Res", &THU_ggH_Res);
    t->Branch("THU_ggH_Mig01", &THU_ggH_Mig01);
    t->Branch("THU_ggH_Mig12", &THU_ggH_Mig12);
    t->Branch("THU_ggH_VBF2j", &THU_ggH_VBF2j);
    t->Branch("THU_ggH_VBF3j", &THU_ggH_VBF3j);
    t->Branch("THU_ggH_PT60", &THU_ggH_PT60);
    t->Branch("THU_ggH_PT120", &THU_ggH_PT120);
    t->Branch("THU_ggH_qmtop", &THU_ggH_qmtop);
    
    t->Branch("trg_singlemu_22", &trg_singlemu_22);
    t->Branch("trg_crossmu_mu19tau20", &trg_crossmu_mu19tau20);
    t->Branch("trg_doubletau_35", &trg_doubletau_35);
    t->Branch("trg_singlemuon_22", &trg_singlemuon_22);
    t->Branch("trg_singlemuon_22_eta2p1", &trg_singlemuon_22_eta2p1);
    t->Branch("trg_singlemuonTk_22", &trg_singlemuonTk_22);
    t->Branch("trg_singlemuonTk_22_eta2p1", &trg_singlemuonTk_22_eta2p1);
    t->Branch("trg_crossmuon_mu19tau20", &trg_crossmuon_mu19tau20);
    t->Branch("trg_crossmuon_mu19tau20_singleL1", &trg_crossmuon_mu19tau20_singleL1);
    t->Branch("trg_singleelectron_25_eta2p1", &trg_singleelectron_25_eta2p1);
    t->Branch("trg_doubletau_35_mediso_eta2p1", &trg_doubletau_35_mediso_eta2p1);
    t->Branch("trg_doubletau_35_medCombiso_eta2p1", &trg_doubletau_35_medCombiso_eta2p1);

    t->Branch("fileEntry", &fileEntry);
    t->Branch("entry", &entry);
    t->Branch("run", &runID);
    t->Branch("lumi", &lumiBlock);
    t->Branch("evt", &eventNr);


}
