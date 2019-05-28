#include "HTTEvent.h"
#include "HTXSClassification.h"
#include "TTree.h"
#include "TLorentzVector.h"
//#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFunctor.h"
#include "TFile.h"

#include "utils/TauTriggerSFs2017/interface/TauTriggerSFs2017.h"


#ifndef __EventWriter__
#define __EventWriter__

class EventWriter
{

 public:
  //ClassDef(EventWriter,0);

  HTTParticle leg1, leg2;
  HTTAnalysis::finalState channel;
  TLorentzVector leg1P4, leg2P4;

  bool isSync;
  bool isMC;

  Float_t lumiWeight;
  Int_t runID;
  Float_t lumiBlock;
  ULong64_t eventNr;
  int entry;
  Int_t fileEntry;

  int npv;
  float npu;
  float rho;
  int gen_match_1;
  int gen_match_2;
  float genPt_1;
  float genPt_2;
  int gen_match_jetId_1;
  int gen_match_jetId_2;
  int NUP;
  //float evtWeight;

  float singleTriggerSFLeg1;
  float singleTriggerSFLeg2;
  float xTriggerSFLeg1;
  float xTriggerSFLeg2;

  float isoWeight_1;
  float isoWeight_2;
  float idWeight_1;
  float idWeight_2;

  float weight;
  float puWeight;
  float genWeight;
  float xsec;
  float genNEventsWeight;
  float idisoweight_1;
  float idisoweight_2;
  float NNLO_ggH_weight;
  float THU_ggH_Mu;
  float THU_ggH_Res;
  float THU_ggH_Mig01;
  float THU_ggH_Mig12;
  float THU_ggH_VBF2j;
  float THU_ggH_VBF3j;
  float THU_ggH_PT60;
  float THU_ggH_PT120;
  float THU_ggH_qmtop;

  float sf_trk;
  float sf_reco;
  float sf_SingleOrCrossTrigger;
  float sf_SingleXorCrossTrigger;
  float sf_SingleTrigger;
  float sf_DoubleTauTight;
  float sf_DoubleTauVTight;

  float trk_sf;
  float trigweight_1;
  float trigweight_2;
  float trigweightVTight_1;
  float trigweightVTight_2;      



  float stitchedWeight;
  float topPtReweightWeightRun1;
  float topPtReweightWeightRun2;
  float zPtReweightWeight;
  float zPtReweightWeight1D;

  float gen_Mll;
  float gen_ll_px;
  float gen_ll_py;
  float gen_ll_pz;
  float gen_vis_Mll;
  float gen_ll_vis_px;
  float gen_ll_vis_py;
  float gen_ll_vis_pz;
  float gen_top_pt_1;
  float gen_top_pt_2;
  int genJets;
  //////////////////////////////////////////////////////////////////  
  int trg_singlemuon_22;
  int trg_singlemuon_22_eta2p1;
  int trg_singlemuonTk_22;
  int trg_singlemuonTk_22_eta2p1;
  int trg_crossmuon_mu19tau20;
  int trg_crossmuon_mu19tau20_singleL1;
  int trg_singleelectron_25_eta2p1;
  int trg_doubletau_35_mediso_eta2p1;
  int trg_doubletau_35_medCombiso_eta2p1;

  // int trg_singletau_leading;
  // int trg_singletau_trailing;
  // int trg_singlemuon_27;
  // int trg_singlemuon_24;
  // int trg_crossmuon_mu20tau27;
  // int trg_singleelectron_35;
  // int trg_singleelectron_32;
  // int trg_singleelectron_27;
  // int trg_crossele_ele24tau30;
  // int trg_doubletau_40_tightiso;
  // int trg_doubletau_40_mediso_tightid;
  // int trg_doubletau_35_tightiso_tightid;
  // int trg_doubletau_35_mediso_HPS;
  // int trg_doubletau;


  int flagMETFilter;
  int Flag_METFilters;
  int Flag_goodVertices;
  int Flag_globalTightHalo2016Filter;
  int Flag_globalSuperTightHalo2016Filter;
  int Flag_HBHENoiseFilter;
  int Flag_HBHENoiseIsoFilter;
  int Flag_EcalDeadCellTriggerPrimitiveFilter;
  int Flag_BadPFMuonFilter;
  int Flag_BadChargedCandidateFilter;
  int Flag_eeBadScFilter;
  int Flag_ecalBadCalibFilter;

  //////////////////////////////////////////////////////////////////  
  float pt_1;
  float phi_1;
  float eta_1;
  float eta_SC_1;
  float m_1;
  int q_1;
  float d0_1;
  float dZ_1;
  float iso_1;
  int againstElectronMVA6_1;  
  int againstElectronLooseMVA6_1;
  int againstElectronMediumMVA6_1;
  int againstElectronTightMVA6_1;
  int againstElectronVLooseMVA6_1;
  int againstElectronVTightMVA6_1;
  int againstMuon3_1;
  int againstMuonLoose3_1;
  int againstMuonTight3_1;

  float byCombinedIsolationDeltaBetaCorrRaw3Hits_1;

  float byIsolationMVA3newDMwoLTraw_1;
  float byIsolationMVA3oldDMwoLTraw_1;
  float byIsolationMVA3newDMwLTraw_1;
  float byIsolationMVA3oldDMwLTraw_1;

  float byIsolationMVArun2017v2DBoldDMwLTraw2017_1;
  int byIsolationMVArun2017v2DBoldDMwLT2017_1;
  int byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1;
  int byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1;
  int byLooseIsolationMVArun2017v2DBoldDMwLT2017_1;
  int byMediumIsolationMVArun2017v2DBoldDMwLT2017_1;
  int byTightIsolationMVArun2017v2DBoldDMwLT2017_1;
  int byVTightIsolationMVArun2017v2DBoldDMwLT2017_1;
  int byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1;

  float chargedIsoPtSum_1;
  float neutralIsoPtSum_1;
  float puCorrPtSum_1;
  int decayModeFindingOldDMs_1;
  int decayMode_1;
  float id_e_mva_nt_loose_1;
  float id_m_loose_1;
  float id_m_medium_1;
  float id_m_tight_1;
  float id_m_highpt_1;
  float id_e_cut_veto_1;
  float id_e_cut_loose_1;
  float id_e_cut_medium_1;
  float id_e_cut_tight_1;
  
  //////////////////////////////////////////////////////////////////
  float pt_2;
  float phi_2;
  float eta_2;
  float m_2;
  int q_2;
  float d0_2;
  float dZ_2;
  float iso_2;
  int againstElectronMVA6_2;  
  int againstElectronLooseMVA6_2;
  int againstElectronMediumMVA6_2;
  int againstElectronTightMVA6_2;
  int againstElectronVLooseMVA6_2;
  int againstElectronVTightMVA6_2;
  int againstMuon3_2;
  int againstMuonLoose3_2;
  int againstMuonTight3_2;

  float byCombinedIsolationDeltaBetaCorrRaw3Hits_2;

  float byIsolationMVA3newDMwoLTraw_2;
  float byIsolationMVA3oldDMwoLTraw_2;
  float byIsolationMVA3newDMwLTraw_2;
  float byIsolationMVA3oldDMwLTraw_2;

  float byIsolationMVArun2017v2DBoldDMwLTraw2017_2;
  int byIsolationMVArun2017v2DBoldDMwLT2017_2;  
  int byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2;  
  int byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2;
  int byLooseIsolationMVArun2017v2DBoldDMwLT2017_2;
  int byMediumIsolationMVArun2017v2DBoldDMwLT2017_2;
  int byTightIsolationMVArun2017v2DBoldDMwLT2017_2;
  int byVTightIsolationMVArun2017v2DBoldDMwLT2017_2;
  int byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2;

  float chargedIsoPtSum_2;
  float neutralIsoPtSum_2;
  float puCorrPtSum_2;
  int decayModeFindingOldDMs_2;
  int decayMode_2;
  //////////////////////////////////////////////////////////////////

  int njets[56];
  int njetspt20[56];
  int njetingap[56];
  int njetingap20[56];
  float mjj[56];
  float jdeta[56];
  float dijetpt[56];
  float dijetphi[56];
  float jdphi[56];
  float jpt_1[56];
  float jpt_2[56];
  float jeta_1[56];
  float jeta_2[56];
  float jphi_1[56];
  float jphi_2[56];
  float jm_1;
  float jm_2;  
  float jrawf_1;
  float jrawf_2;  
  float jmva_1;
  float jmva_2;  
  float jcsv_1;
  float jcsv_2;

  int nbtag[5];
  float bpt_1[5];
  float bpt_2[5];
  float beta_1[5];
  float beta_2[5];
  float bphi_1[5];
  float bphi_2[5];
  float brawf_1[5];
  float brawf_2[5];
  float bmva_1[5];
  float bmva_2[5];
  float bcsv_1[5];
  float bcsv_2[5];
  //////////////////////////////////////////////////////////////////
  float corrmet;
  float corrmet_ex;
  float corrmet_ey;
  float corrmetphi;
  float metcov00;
  float metcov01;
  float metcov10;
  float metcov11;
  //////////////////////////////////////////////////////////////////
  float met[56];
  float met_ex[56];
  float met_ey[56];
  float metphi[56];
  float mt_1[56];
  float mt_2[56];
  float mt_tot[56];
  float m_sv[56];
  float pt_sv[56];
  float pt_tt[56];
  float pt_ttjj[56];
  float m_ttjj[56];
  float pt_sum[56];
  int htxs_reco_ggf[56];
  int htxs_reco_vbf[56];
  int htxs_stage1cat;
  //////////////////////////////////////////////////////////////////
  float eleTauFakeRateWeight;
  float muTauFakeRateWeight;
  float antilep_tauscaling;
  //////////////////////////////////////////////////////////////////
  bool passesTauLepVetos;
  bool passesThirdLepVeto;
  bool passesDiMuonVeto;
  bool passesDiElectronVeto;
  bool diMuonVeto;
  bool diElectronVeto;
  bool dilepton_veto;
  bool extramuon_veto;
  bool extraelec_veto;
  //////////////////////////////////////////////////////////////////
  float pzetavis;
  float pzetamiss;
  float dzeta;
  //////////////////////////////////////////////////////////////////
  float pt_vis;
  float mt_3;
  float pfpt_tt;
  float m_vis;
  float m_coll;
  float dphi;
  //////////////////////////////////////////////////////////////////
  float dr_leptau;
  float jeta1eta2;
  float met_centrality;
  float lep_etacentrality;
  float sphericity;
  //////////////////////////////////////////////////////////////////
  int pdg1;
  int pdg2;

  EventWriter(){}  
  ~EventWriter(){}  

  void setDefault();
  void fill(HTTEvent *ev, HTTJetCollection *jets, std::vector<HTTParticle> leptons, HTTPair *pair);
  void initTree(TTree *t, vector< pair< string, pair<string,bool> > >jecShifts_, bool isMC_, bool isSync_, vector<string> metShifts);

  double calcSphericity(std::vector<TLorentzVector> p);
  double calcSphericityFromMatrix(TMatrixD M);

  int getGenMatch_jetId(TLorentzVector selObj, HTTJetCollection *jets);

  double calcDR(double eta1, double phi1, double eta2, double phi2);


  void fillScalefactors();
  void fillLeptonFakeRateWeights();
  void fillStitchingWeight(HTTEvent::sampleTypeEnum sampleType);
  void fillLeg1Branches();
  void fillLeg2Branches();
  void fillJetBranches(HTTJetCollection *jets);
  void fillPairBranches(HTTPair *pair, HTTJetCollection *jets);
  void fillAdditionalLeptons( std::vector<HTTParticle> leptons, HTTPair *pair);

  RooWorkspace *w;
  TauTriggerSFs2017 *tauTrigSFTight;
  TauTriggerSFs2017 *tauTrigSFVTight;
  vector< pair< string, pair<string,bool> > > jecShifts;
  vector< string > metShifts;
  vector< pair< string, pair<string, string> > >  btagShifts;

  vector<TLorentzVector> addlepton_p4;
  vector<double> addlepton_pt;
  vector<double> addlepton_eta;
  vector<double> addlepton_phi;
  vector<double> addlepton_m;
  vector<double> addlepton_iso;
  vector<int> addlepton_pdgId;
  vector<int> addlepton_mc_match;
  vector<double> addlepton_d0;
  vector<double> addlepton_dZ;
  vector<double> addlepton_mt;
  vector<double> addlepton_mvis;
  vector<double> addlepton_tauCombIso;
  vector<int> addlepton_tauID;
  vector<int> addlepton_tauDM;
  vector<int> addlepton_tauAntiEle;
  vector<int> addlepton_tauAntiMu;
  //////////////////////////////////////////////////////////////////
  int nadditionalMu;
  vector<double> addmuon_pt;
  vector<double> addmuon_eta;
  vector<double> addmuon_phi;
  vector<double> addmuon_m;
  vector<int> addmuon_q;
  vector<double> addmuon_iso;
  vector<int> addmuon_gen_match;
  //////////////////////////////////////////////////////////////////
  int nadditionalEle;
  vector<double> addele_pt;
  vector<double> addele_eta;
  vector<double> addele_phi;
  vector<double> addele_m;
  vector<int> addele_q;
  vector<double> addele_iso;
  vector<int> addele_gen_match;
  //////////////////////////////////////////////////////////////////
  int nadditionalTau;
  vector<double> addtau_pt;
  vector<double> addtau_eta;
  vector<double> addtau_phi;
  vector<double> addtau_m;
  vector<double> addtau_q;
  vector<double> addtau_byIsolationMVArun2v1DBoldDMwLTraw;
  vector<double> addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits;
  vector<int> addtau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
  vector<int> addtau_byTightCombinedIsolationDeltaBetaCorr3Hits;
  vector<int> addtau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
  vector<int> addtau_byVVLooseIsolationMVArun2v1DBoldDMwLT;
  vector<int> addtau_byVLooseIsolationMVArun2v1DBoldDMwLT;
  vector<int> addtau_byLooseIsolationMVArun2v1DBoldDMwLT;
  vector<int> addtau_byMediumIsolationMVArun2v1DBoldDMwLT;
  vector<int> addtau_byTightIsolationMVArun2v1DBoldDMwLT;
  vector<int> addtau_byVTightIsolationMVArun2v1DBoldDMwLT;

  vector<int> addtau_passesTauLepVetos;
  vector<int> addtau_decayMode;
  vector<double> addtau_d0;
  vector<double> addtau_dZ;
  vector<int> addtau_gen_match;
  vector<double> addtau_mt;
  vector<double> addtau_mvis;
  //////////////////////////////////////////////////////////////////

};
#endif
