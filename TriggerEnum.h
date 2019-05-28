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
HLT_IsoMu22 = 0,
HLT_IsoMu22_eta2p1,
HLT_IsoTkMu22,
HLT_IsoTkMu22_eta2p1,
HLT_IsoMu19_eta2p1_LooseIsoPFTau20,
HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1,
HLT_Ele25_eta2p1_WPTight_Gsf,
HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg,
HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg,
NONE
};

vector<TriggerData>  getTriggerSettings(){
    vector<TriggerData> triggerBits_;
    TriggerData aTrgData;

    // These Filter bits are taken from cmssw/PhysicsTools/NanoAOD/python/triggerObjects_cff.py
    // See line 138-160 for era-dependent modification used for the 2016 data set

    // Electron
    // 0 = *CaloIdLTrackIdLIsoVL*TrackIso*Filter
    // 1 = hltEle*WPTight*TrackIsoFilter*
    // 2 = hltEle*WPLoose*TrackIsoFilter
    // 3 = *OverlapFilter*IsoEle*PFTau*

    // Muon
    // 0 = *RelTrkIso*Filtered0p4
    // 1 = hltL3cr*IsoFiltered0p09
    // 2 = *OverlapFilter*IsoMu*PFTau*
    // 3 = hltL3f*IsoFiltered0p09
    // sel.qualityBitsDoc = cms.string("1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau, 8 = IsoTkMu")
    
    // Tau
    // 0 = *LooseIso*-*VLooseIso*
    // 1 = *Medium*Iso*
    // 2 = *VLooseIso*
    // 3 = 0
    // 4 = hltL2TauIsoFilter
    // 5 = *OverlapFilter*IsoMu*
    // 6 = *OverlapFilter*IsoEle*
    // 7 = *L1HLTMatched*
    // 8 = *Dz02*
    // sel.qualityBitsDoc = cms.string("1 = LooseIso, 2 = Medium(Comb)Iso, 4 = VLooseIso, 8 = None, 16 = L2p5 pixel iso, 32 = OverlapFilter IsoMu, 64 = OverlapFilter IsoEle, 128 = L1-HLT matched, 256 = Dz")
   
    // single-mu trigger with mu-filter to match: hltL3cr IsoL1sMu20L1f0L2f10QL3f22QL3trk IsoFiltered0p09 ->2
    // hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_IsoMu22";
    triggerBits_.back().leg1Id=13;
    // triggerBits_.back().leg1BitMask=(1<<1);    
    triggerBits_.back().leg1BitMask=0;    
    // single-mu trigger  with mu-filter to match: hltL3cr IsoL1sSingleMu20erL1f0L2f10QL3f22QL3trkIso Filtered0p09
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_IsoMu22_eta2p1";
    triggerBits_.back().leg1Id=13;
    // triggerBits_.back().leg1BitMask=(1<<1);
    triggerBits_.back().leg1BitMask=0;
    // single-mu trigger  with mu-filter to match: hltL3f L1sMu20L1f0Tkf22QL3trk IsoFiltered0p09
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_IsoTkMu22";
    triggerBits_.back().leg1Id=13;
    // triggerBits_.back().leg1BitMask=(1<<3);
    triggerBits_.back().leg1BitMask=0;
    // single-mu trigger  with mu-filter to match: hltL3 fL1sMu20erL1f0Tkf22QL3trkIso Filtered0p09
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_IsoTkMu22_eta2p1";
    triggerBits_.back().leg1Id=13;
    // triggerBits_.back().leg1BitMask=(1<<3);
    triggerBits_.back().leg1BitMask=0;


    // mu-tauh triggers with 
    // mu-filter   : hltL3cr IsoL1sMu18erTauJet20erL1f0L2f10QL3f19QL3trk IsoFiltered0p09 and hlt OverlapFilter IsoMu 19LooseIso PFTau 20
    // tauh-filter : hltPFTau20Track LooseIso AgainstMuon and hlt OverlapFilter IsoMu 19LooseIsoPFTau20
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_IsoMu19_eta2p1_LooseIsoPFTau20";
    triggerBits_.back().leg1Id=13;
    // triggerBits_.back().leg1BitMask=(1<<1) + (1<<2); //iso+OL (overlap)
    triggerBits_.back().leg1BitMask=0;
    triggerBits_.back().leg2Id=15;
    // triggerBits_.back().leg2BitMask=(1<<0) + (1<<5); //looseIso+OL
    triggerBits_.back().leg2BitMask=0;
    // mu-tauh triggers with 
    // mu-filter   : hltL3cr IsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trk IsoFiltered0p09 and hlt OverlapFilter Single IsoMu 19LooseIso PFTau 20
    // tauh-filter : hltPFTau20Track LooseIso AgainstMuon and hlt OverlapFilter Single IsoMu 19LooseIsoPFTau20
    triggerBits_.push_back(aTrgData); 
    triggerBits_.back().path_name="HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1";
    triggerBits_.back().leg1Id=13;
    // triggerBits_.back().leg1BitMask=(1<<1) + (1<<2); //iso+OL (overlap)
    triggerBits_.back().leg1BitMask=0;
    triggerBits_.back().leg2Id=15;
    // triggerBits_.back().leg2BitMask=(1<<0) + (1<<5); //looseIso+OL
    triggerBits_.back().leg2BitMask=0;
    
    
    // single-ele trigger with ele filter to match : hltEle 25er WPTight Gsf TrackIsoFilter
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_Ele25_eta2p1_WPTight_Gsf";
    triggerBits_.back().leg1Id=11;
    // triggerBits_.back().leg1BitMask=(1<<1);
    triggerBits_.back().leg1BitMask=0;
  

    // double-tauh trigger with tau filter to match: 
    // hltDoublePFTau35TrackPt1 Medium Isolation Dz02Reg
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg";
    triggerBits_.back().leg1Id=15;
    // triggerBits_.back().leg1BitMask=(1<<1) + (1<<8) //MediumIso+dz
    triggerBits_.back().leg1BitMask=0;
    triggerBits_.back().leg2Id=15;
    // triggerBits_.back().leg2BitMask=(1<<1) + (1<<8)
    triggerBits_.back().leg2BitMask=0;
    // double-tauh trigger with tau filter to match: 
    // hltDoublePFTau35TrackPt1 Medium CombinedIsolation Dz02 Reg
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg";
    triggerBits_.back().leg1Id=15;
    // triggerBits_.back().leg1BitMask=(1<<1) + (1<<8) //MediumIso+dz
    triggerBits_.back().leg1BitMask=0;
    triggerBits_.back().leg2Id=15;
    // triggerBits_.back().leg2BitMask=(1<<1) + (1<<8)
    triggerBits_.back().leg2BitMask=0;
    


    return triggerBits_;
}