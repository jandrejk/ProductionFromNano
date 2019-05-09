struct TriggerData {
    std::string path_name;
    unsigned int leg1Id=0,         leg2Id=0;//0-undefined, 11-electron, 13-muon, 15-tau
    int          leg1BitMask=0,    leg2BitMask=0;//definition depends on Id, cf. PhysicsTools/NanoAOD/python/triggerObjects_cff.py
    float        leg1Pt=0.,        leg2Pt=0.;
    float        leg1L1Pt=0.,      leg2L1Pt=0.;
    float        leg1Eta=0.,       leg2Eta=0.;
    float        leg1OfflinePt=0., leg2OfflinePt=0.;
};

enum class TriggerEnum {
HLT_IsoMu24 = 0,
HLT_IsoMu27,
HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1,
HLT_Ele27_WPTight_Gsf,
HLT_Ele32_WPTight_Gsf,
HLT_Ele35_WPTight_Gsf,
HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1,
HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1,
HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg,
HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg,
HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg,
HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg,
NONE
};

vector<TriggerData>  getTriggerSettings(){
	vector<TriggerData> triggerBits_;
    TriggerData aTrgData;

    // 2017 94X  Filter

    // Electron
    // 0 = *CaloIdLTrackIdLIsoVL*TrackIso*Filter
    // 1 = hltEle*WPTight*TrackIsoFilter*
    // 2 = hltEle*WPLoose*TrackIsoFilter
    // 3 = *OverlapFilterIsoEle*PFTau*
    // 4 = hltEle*Ele*CaloIdLTrackIdLIsoVL*Filter
    // 5 = hltMu*TrkIsoVVL*Ele*CaloIdLTrackIdLIsoVL*Filter*
    // 6 = *OverlapFilterIsoEle*PFTau*
    // 7 = hltEle*Ele*Ele*CaloIdLTrackIdLDphiLeg*Filter
    // 8 = max(filter('hltL3fL1Mu*DoubleEG*Filtered*'),filter('hltMu*DiEle*CaloIdLTrackIdLElectronleg*Filter'))
    // 9 = max(filter('hltL3fL1DoubleMu*EG*Filter*'),filter('hltDiMu*Ele*CaloIdLTrackIdLElectronleg*Filter'))
            

    // Muon
    // 0 = *RelTrkIsoVVLFiltered0p4
    // 1 = hltL3crIso*Filtered0p07
    // 2 = *OverlapFilterIsoMu*PFTau*
    // 3 = max(filter('hltL3crIsoL1*SingleMu*Filtered0p07'),filter('hltL3crIsoL1sMu*Filtered0p07'))
    // 4 = hltDiMuon*Filtered*
    // 5 = hltMu*TrkIsoVVL*Ele*CaloIdLTrackIdLIsoVL*Filter*
    // 6 = hltOverlapFilterIsoMu*PFTau*
    // 7 = hltL3fL1TripleMu*
    // 8 = max(filter('hltL3fL1DoubleMu*EG*Filtered*'),filter('hltDiMu*Ele*CaloIdLTrackIdLElectronleg*Filter'))
    // 9 = max(filter('hltL3fL1Mu*DoubleEG*Filtered*'),filter('hltMu*DiEle*CaloIdLTrackIdLElectronleg*Filter'))

    // Tau
    // 0 = *LooseChargedIso*
    // 1 = *MediumChargedIso*
    // 2 = *TightChargedIso*
    // 3 = *TightOOSCPhotons*
    // 4 = *Hps*
    // 5 = hltSelectedPFTau*MediumChargedIsolationL1HLTMatched*
    // 6 = hltDoublePFTau*TrackPt1*ChargedIsolation*Dz02Reg
    // 7 = hltOverlapFilterIsoEle*PFTau*
    // 8 = hltOverlapFilterIsoMu*PFTau*
    // 9 = hltDoublePFTau*TrackPt1*ChargedIsolation*
           

    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_IsoMu24";
    triggerBits_.back().leg1Id=13;
    triggerBits_.back().leg1BitMask=(1<<3);
    // triggerBits_.back().leg1Pt=24;
    // triggerBits_.back().leg1L1Pt=22;
    // triggerBits_.back().leg1OfflinePt=25;

    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_IsoMu27";
    triggerBits_.back().leg1Id=13;
    triggerBits_.back().leg1BitMask=(1<<3);
    // triggerBits_.back().leg1Pt=27;
    // triggerBits_.back().leg1L1Pt=22; //22 or 25...
    // triggerBits_.back().leg1OfflinePt=28;

    // Mu tauh triggers
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1";
    triggerBits_.back().leg1Id=13;
    triggerBits_.back().leg1BitMask=(1<<1) + (1<<2); //iso+OL
    // triggerBits_.back().leg1Pt=20;
    // triggerBits_.back().leg1Eta=2.1;
    triggerBits_.back().leg1L1Pt=-1;
    //  triggerBits_.back().leg1L1Pt=18;
    // triggerBits_.back().leg1OfflinePt=20;
    triggerBits_.back().leg2Id=15;
    triggerBits_.back().leg2BitMask=(1<<0) + (1<<8); //looseChargedIso+OL
    // triggerBits_.back().leg2Pt=27;
    // triggerBits_.back().leg2Eta=2.1;
    //  triggerBits_.back().leg2L1Pt=20;
    triggerBits_.back().leg2L1Pt=-1;

    //Single e triggers
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_Ele27_WPTight_Gsf";
    triggerBits_.back().leg1Id=11;
    triggerBits_.back().leg1BitMask=(1<<1);
    // triggerBits_.back().leg1Pt=27;
    triggerBits_.back().leg1L1Pt=-1;

    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_Ele32_WPTight_Gsf";
    triggerBits_.back().leg1Id=11;
    triggerBits_.back().leg1BitMask=(1<<1);
    // triggerBits_.back().leg1Pt=32;
    triggerBits_.back().leg1L1Pt=-1;

    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_Ele35_WPTight_Gsf";
    triggerBits_.back().leg1Id=11;
    triggerBits_.back().leg1BitMask=(1<<1);
    // triggerBits_.back().leg1Pt=35;
    triggerBits_.back().leg1L1Pt=-1;


    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1";
    triggerBits_.back().leg1Id=11;
    triggerBits_.back().leg1BitMask=(1<<1)+(1<<3);
    // triggerBits_.back().leg1Pt=24;
    // triggerBits_.back().leg1Eta=2.1;
    triggerBits_.back().leg2Id=15;
    triggerBits_.back().leg2BitMask=(1<<0) + (1<<7); 
    // triggerBits_.back().leg2Pt=30;
    // triggerBits_.back().leg2Eta=2.1;


    // Single tauh triggers
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1";
    triggerBits_.back().leg1Id=15;
    triggerBits_.back().leg1BitMask=2;

    // tauh tauh triggers
    ///9th tau bit(1<<8) for di-tau dz filter  (should be OK 80X triggers)
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg";
    triggerBits_.back().leg1Id=15;
    triggerBits_.back().leg1BitMask= (1<<2) + (1<<3) + (1<<6); //TightChargedIso+photons+dz
    // triggerBits_.back().leg1Pt=35;
    // triggerBits_.back().leg1Eta=2.1;
    triggerBits_.back().leg2Id=15;
    triggerBits_.back().leg2BitMask= (1<<2) + (1<<3) + (1<<6);
    // triggerBits_.back().leg2Pt=35;
    // triggerBits_.back().leg2Eta=2.1;

    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg";
    triggerBits_.back().leg1Id=15;
    triggerBits_.back().leg1BitMask=(1<<1) + (1<<3) + (1<<6); //MediumChargedIso+photons+dz
    // triggerBits_.back().leg1Pt=40;
    // triggerBits_.back().leg1Eta=2.1;
    triggerBits_.back().leg2Id=15;
    triggerBits_.back().leg2BitMask=(1<<1) + (1<<3) + (1<<6);
    // triggerBits_.back().leg2Pt=40;
    // triggerBits_.back().leg2Eta=2.1;

    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg";
    triggerBits_.back().leg1Id=15;
    triggerBits_.back().leg1BitMask=(1<<2) + (1<<6);; //TightChargedIso+dz
    // triggerBits_.back().leg1Pt=40;
    // triggerBits_.back().leg1Eta=2.1;
    triggerBits_.back().leg2Id=15;
    triggerBits_.back().leg2BitMask=(1<<2) + (1<<6);;
    // triggerBits_.back().leg2Pt=40;
    // triggerBits_.back().leg2Eta=2.1;

    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg";
    triggerBits_.back().leg1Id=15;
    triggerBits_.back().leg1BitMask=(1<<1) + (1<<4);
    
    triggerBits_.back().leg2Id=15;
    triggerBits_.back().leg2BitMask=(1<<1) + (1<<4);

    return triggerBits_;
}