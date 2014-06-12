//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 18 14:43:11 2014 by ROOT version 5.28/00b
// from TTree SimpleNtupleEoverP/SimpleNtupleEoverP
// found on file: input.root
//////////////////////////////////////////////////////////

#ifndef analysis_h
#define analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;

const int NELE_max = 30;
const int NGEN_max = 50;
const int NPHO_max = 200;


class analysis {
    public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain
    
    // Declaration of leaf types
    Long64_t        bxId;
    Long64_t        eventId;
    Int_t           lumiId;
    Int_t           runId;
    Int_t           timeStampHigh;
    Int_t           isW;
    Int_t           isZ;
    Int_t           PV_n;
    Float_t         PV_z;
    Float_t         PV_d0;
    Float_t         PV_SumPt;
    Float_t         rho;
    Int_t           NELE;
    Float_t         ELE_pt[NELE_max];   //[NELE]
    Float_t         ELE_eta[NELE_max];   //[NELE]
    Float_t         ELE_phi[NELE_max];   //[NELE]
    Float_t         ELE_px[NELE_max];   //[NELE]
    Float_t         ELE_py[NELE_max];   //[NELE]
    Float_t         ELE_pz[NELE_max];   //[NELE]
    Float_t         ELE_e[NELE_max];   //[NELE]
    Int_t           ELE_q[NELE_max];   //[NELE]
    Int_t           ELE_id[NELE_max];   //[NELE]
    Float_t         ELE_sigmaIetaIeta[NELE_max];   //[NELE]
    Float_t         ELE_DphiIn[NELE_max];   //[NELE]
    Float_t         ELE_DetaIn[NELE_max];   //[NELE]
    Float_t         ELE_HOverE[NELE_max];   //[NELE]
    Float_t         ELE_ooemoop[NELE_max];   //[NELE]
    Float_t         ELE_tkIso[NELE_max];   //[NELE]
    Float_t         ELE_emIso[NELE_max];   //[NELE]
    Float_t         ELE_hadIso[NELE_max];   //[NELE]
    Float_t         ELE_effAreaForIso[NELE_max];   //[NELE]
    Float_t         ELE_combIso[NELE_max];   //[NELE]
    Float_t         ELE_dxy[NELE_max];   //[NELE]
    Float_t         ELE_dz[NELE_max];   //[NELE]
    Float_t         ELE_scE_regression[NELE_max];   //[NELE]
    Float_t         ELE_scEtaWidth[NELE_max];   //[NELE]
    Float_t         ELE_scPhiWidth[NELE_max];   //[NELE]
    Float_t         ELE_scERaw[NELE_max];   //[NELE]
    Float_t         ELE_scE[NELE_max];   //[NELE]
    Float_t         ELE_e5x5[NELE_max];   //[NELE]
    Float_t         ELE_e1x5[NELE_max];   //[NELE]
    Float_t         ELE_e2x5Max[NELE_max];   //[NELE]
    Float_t         ELE_es[NELE_max];   //[NELE]
    Float_t         ELE_fbrem[NELE_max];   //[NELE]
    Int_t           psize;
    Float_t         pSC_energy[NPHO_max];   //[psize]
    Double_t        pEta[NPHO_max];   //[psize]
    Double_t        pPhi[NPHO_max];   //[psize]
    Double_t        pEsenergy[NPHO_max];   //[psize]
    Double_t        pRawenergy[NPHO_max];   //[psize]
    Float_t         pPt[NPHO_max];   //[psize]
    Float_t         p1x5_energy[NPHO_max];   //[psize]
    Float_t         p2x5_energy[NPHO_max];   //[psize]
    Float_t         p3x3_energy[NPHO_max];   //[psize]
    Float_t         p5x5_energy[NPHO_max];   //[psize]
    Float_t         pR9[NPHO_max];   //[psize]
    Bool_t          pGap[NPHO_max];   //[psize]
    Float_t         pSigmaIetaIeta[NPHO_max];   //[psize]
    Float_t         pHoverE[NPHO_max];   //[psize]
    Float_t         p_photonenergy[NPHO_max];   //[psize]
    Int_t           p_nCrystals[NPHO_max];   //[psize]
    Int_t           p_adc[NPHO_max][25][10];   //[psize]
    Int_t           p_gainId[NPHO_max][25][10];   //[psize]
    Float_t         p_intercalib[NPHO_max][25];   //[psize]
    Int_t           p_recieta[NPHO_max][25];   //[psize]
    Int_t           p_reciphi[NPHO_max][25];   //[psize]
    Float_t         p_receta[NPHO_max][25];   //[psize]
    Float_t         p_recphi[NPHO_max][25];   //[psize]
    Float_t         p_recenergy[NPHO_max][25];   //[psize]
    Float_t         p_rectime[NPHO_max][25];   //[psize]
    UInt_t          p_recflag[NPHO_max][25];   //[psize]
    Bool_t          p_recflags[NPHO_max][25][19];   //[psize]
    Int_t           NGEN;
    Float_t         GEN_pt[NGEN_max];   //[NGEN]
    Float_t         GEN_eta[NGEN_max];   //[NGEN]
    Float_t         GEN_phi[NGEN_max];   //[NGEN]
    Float_t         GEN_px[NGEN_max];   //[NGEN]
    Float_t         GEN_py[NGEN_max];   //[NGEN]
    Float_t         GEN_pz[NGEN_max];   //[NGEN]
    Float_t         GEN_e[NGEN_max];   //[NGEN]
    Float_t         GEN_mass[NGEN_max];   //[NGEN]
    Int_t           GEN_q[NGEN_max];   //[NGEN]
    Int_t           GEN_id[NGEN_max];   //[NGEN]
    Int_t           GEN_status[NGEN_max];   //[NGEN]
    Int_t           GEN_parent[NGEN_max];   //[NGEN]
    Float_t         ele1_charge;
    Float_t         ele1_p;
    Float_t         ele1_pt;
    Float_t         ele1_eta;
    Float_t         ele1_phi;
    Int_t           ele1_isEB;
    Int_t           ele1_isEBEEGap;
    Int_t           ele1_isEBEtaGap;
    Int_t           ele1_isEBPhiGap;
    Int_t           ele1_isEEDeeGap;
    Int_t           ele1_isEERingGap;
    Int_t           ele1_isTrackerDriven;
    Float_t         ele1_sigmaIetaIeta;
    Float_t         ele1_DphiIn;
    Float_t         ele1_DetaIn;
    Float_t         ele1_HOverE;
    Float_t         ele1_tkIso;
    Float_t         ele1_emIso;
    Float_t         ele1_hadIso;
    Float_t         ele1_scERaw;
    Float_t         ele1_scEtRaw;
    Float_t         ele1_scE;
    Float_t         ele1_scEt;
    Float_t         ele1_scE_regression;
    Float_t         ele1_scEerr_regression;
    Float_t         ele1_scE_regression_PhotonTuned;
    Float_t         ele1_scEerr_regression_PhotonTuned;
    Float_t         ele1_scERaw_PUcleaned;
    Float_t         ele1_es;
    Float_t         ele1_scLaserCorr;
    Float_t         ele1_scCrackCorr;
    Float_t         ele1_scLocalContCorr;
    Float_t         ele1_scEta;
    Float_t         ele1_scPhi;
    Float_t         ele1_scLocalEta;
    Float_t         ele1_scLocalPhi;
    Float_t         ele1_scEtaWidth;
    Float_t         ele1_scPhiWidth;
    Float_t         ele1_scEtaWidth_PUcleaned;
    Float_t         ele1_scPhiWidth_PUcleaned;
    Float_t         ele1_fCorrection_PUcleaned;
    Float_t         ele1_fEta;
    Float_t         ele1_fEtaCorr;
    Float_t         ele1_tkP;
    Float_t         ele1_tkPt;
    Float_t         ele1_fbrem;
    Float_t         ele1_dxy_PV;
    Float_t         ele1_dz_PV;
    Float_t         ele1_sigmaP;
    Float_t         ele1_eSeedBC;
    Float_t         ele1_e5x5;
    Float_t         ele1_e3x3;
    Float_t         ele1_scNxtal;
    Int_t           ele1_bcN;
    Float_t         ele1_5x5LaserCorr;
    Float_t         ele1_3x3LaserCorr;
    Float_t         ele1_seedE;
    Float_t         ele1_seedLaserAlpha;
    Float_t         ele1_seedLaserCorr;
    Float_t         ele1_seedICConstant;
    Int_t           ele1_seedIeta;
    Int_t           ele1_seedIphi;
    Int_t           ele1_seedIx;
    Int_t           ele1_seedIy;
    Int_t           ele1_seedZside;
    Float_t         ele1_EOverP;
    Int_t           ele1_nRecHits;
    vector<float>   *ele1_recHit_E;
    vector<int>     *ele1_recHit_flag;
    vector<int>     *ele1_recHit_hashedIndex;
    vector<int>     *ele1_recHit_ietaORix;
    vector<int>     *ele1_recHit_iphiORiy;
    vector<int>     *ele1_recHit_zside;
    vector<float>   *ele1_recHit_laserCorrection;
    vector<float>   *ele1_recHit_Alpha;
    vector<float>   *ele1_recHit_ICConstant;
    Float_t         ele2_charge;
    Float_t         ele2_p;
    Float_t         ele2_pt;
    Float_t         ele2_eta;
    Float_t         ele2_phi;
    Int_t           ele2_isEB;
    Int_t           ele2_isEBEEGap;
    Int_t           ele2_isEBEtaGap;
    Int_t           ele2_isEBPhiGap;
    Int_t           ele2_isEEDeeGap;
    Int_t           ele2_isEERingGap;
    Int_t           ele2_isTrackerDriven;
    Float_t         ele2_sigmaIetaIeta;
    Float_t         ele2_DphiIn;
    Float_t         ele2_DetaIn;
    Float_t         ele2_HOverE;
    Float_t         ele2_tkIso;
    Float_t         ele2_emIso;
    Float_t         ele2_hadIso;
    Float_t         ele2_dxy_PV;
    Float_t         ele2_dz_PV;
    Float_t         ele2_sigmaP;
    Float_t         ele2_scERaw;
    Float_t         ele2_scEtRaw;
    Float_t         ele2_scE;
    Float_t         ele2_scEt;
    Float_t         ele2_scE_regression;
    Float_t         ele2_scEerr_regression;
    Float_t         ele2_scE_regression_PhotonTuned;
    Float_t         ele2_scEerr_regression_PhotonTuned;
    Float_t         ele2_scERaw_PUcleaned;
    Float_t         ele2_es;
    Float_t         ele2_scLaserCorr;
    Float_t         ele2_scCrackCorr;
    Float_t         ele2_scLocalContCorr;
    Float_t         ele2_scEta;
    Float_t         ele2_scPhi;
    Float_t         ele2_scLocalEta;
    Float_t         ele2_scLocalPhi;
    Float_t         ele2_scEtaWidth;
    Float_t         ele2_scPhiWidth;
    Float_t         ele2_scEtaWidth_PUcleaned;
    Float_t         ele2_scPhiWidth_PUcleaned;
    Float_t         ele2_fCorrection_PUcleaned;
    Float_t         ele2_fEta;
    Float_t         ele2_fEtaCorr;
    Float_t         ele2_tkP;
    Float_t         ele2_tkPt;
    Float_t         ele2_fbrem;
    Float_t         ele2_eSeedBC;
    Float_t         ele2_e5x5;
    Float_t         ele2_e3x3;
    Float_t         ele2_scNxtal;
    Int_t           ele2_bcN;
    Float_t         ele2_5x5LaserCorr;
    Float_t         ele2_3x3LaserCorr;
    Float_t         ele2_seedE;
    Float_t         ele2_seedLaserAlpha;
    Float_t         ele2_seedLaserCorr;
    Float_t         ele2_seedICConstant;
    Int_t           ele2_seedIeta;
    Int_t           ele2_seedIphi;
    Int_t           ele2_seedIx;
    Int_t           ele2_seedIy;
    Int_t           ele2_seedZside;
    Float_t         ele2_EOverP;
    Int_t           ele2_nRecHits;
    vector<float>   *ele2_recHit_E;
    vector<int>     *ele2_recHit_flag;
    vector<int>     *ele2_recHit_hashedIndex;
    vector<int>     *ele2_recHit_ietaORix;
    vector<int>     *ele2_recHit_iphiORiy;
    vector<int>     *ele2_recHit_zside;
    vector<float>   *ele2_recHit_laserCorrection;
    vector<float>   *ele2_recHit_Alpha;
    vector<float>   *ele2_recHit_ICConstant;
    Float_t         met_et;
    Float_t         met_phi;
    Float_t         ele1Met_mt;
    Float_t         ele1Met_Dphi;
    Float_t         ele1ele2_m;
    Float_t         ele1ele2_scM;
    Float_t         ele1ele2_scM_regression;
    Float_t         PUit_TrueNumInteractions;
    Int_t           PUit_NumInteractions;
    Float_t         PUit_zpositions;
    Float_t         PUit_sumpT_lowpT;
    Float_t         PUit_sumpT_highpT;
    Float_t         PUit_ntrks_lowpT;
    Float_t         PUit_ntrks_highpT;
    Float_t         PUoot_early_TrueNumInteractions;
    Int_t           PUoot_early;
    Float_t         PUoot_early_zpositions;
    Float_t         PUoot_early_sumpT_lowpT;
    Float_t         PUoot_early_sumpT_highpT;
    Float_t         PUoot_early_ntrks_lowpT;
    Float_t         PUoot_early_ntrks_highpT;
    Float_t         PUoot_late_TrueNumInteractions;
    Int_t           PUoot_late;
    Float_t         PUoot_late_zpositions;
    Float_t         PUoot_late_sumpT_lowpT;
    Float_t         PUoot_late_sumpT_highpT;
    Float_t         PUoot_late_ntrks_lowpT;
    Float_t         PUoot_late_ntrks_highpT;
    Float_t         mcV_E;
    Float_t         mcV_Px;
    Float_t         mcV_Py;
    Float_t         mcV_Pz;
    Int_t           mcV_Charge;
    Int_t           mcV_PdgId;
    Float_t         mcF1_fromV_E;
    Float_t         mcF1_fromV_Px;
    Float_t         mcF1_fromV_Py;
    Float_t         mcF1_fromV_Pz;
    Int_t           mcF1_fromV_Charge;
    Int_t           mcF1_fromV_PdgId;
    Float_t         mcF2_fromV_E;
    Float_t         mcF2_fromV_Px;
    Float_t         mcF2_fromV_Py;
    Float_t         mcF2_fromV_Pz;
    Int_t           mcF2_fromV_Charge;
    Int_t           mcF2_fromV_PdgId;
    
    // List of branches
    TBranch        *b_bxId;   //!
    TBranch        *b_eventId;   //!
    TBranch        *b_lumiId;   //!
    TBranch        *b_runId;   //!
    TBranch        *b_timeStampHigh;   //!
    TBranch        *b_isW;   //!
    TBranch        *b_isZ;   //!
    TBranch        *b_PV_n;   //!
    TBranch        *b_PV_z;   //!
    TBranch        *b_PV_d0;   //!
    TBranch        *b_PV_SumPt;   //!
    TBranch        *b_rho;   //!
    TBranch        *b_NELE;   //!
    TBranch        *b_ELE_pt;   //!
    TBranch        *b_ELE_eta;   //!
    TBranch        *b_ELE_phi;   //!
    TBranch        *b_ELE_px;   //!
    TBranch        *b_ELE_py;   //!
    TBranch        *b_ELE_pz;   //!
    TBranch        *b_ELE_e;   //!
    TBranch        *b_ELE_q;   //!
    TBranch        *b_ELE_id;   //!
    TBranch        *b_ELE_sigmaIetaIeta;   //!
    TBranch        *b_ELE_DphiIn;   //!
    TBranch        *b_ELE_DetaIn;   //!
    TBranch        *b_ELE_HOverE;   //!
    TBranch        *b_ELE_ooemoop;   //!
    TBranch        *b_ELE_tkIso;   //!
    TBranch        *b_ELE_emIso;   //!
    TBranch        *b_ELE_hadIso;   //!
    TBranch        *b_ELE_effAreaForIso;   //!
    TBranch        *b_ELE_combIso;   //!
    TBranch        *b_ELE_dxy;   //!
    TBranch        *b_ELE_dz;   //!
    TBranch        *b_ELE_scE_regression;   //!
    TBranch        *b_ELE_scEtaWidth;   //!
    TBranch        *b_ELE_scPhiWidth;   //!
    TBranch        *b_ELE_scERaw;   //!
    TBranch        *b_ELE_scE;   //!
    TBranch        *b_ELE_e5x5;   //!
    TBranch        *b_ELE_e1x5;   //!
    TBranch        *b_ELE_e2x5Max;   //!
    TBranch        *b_ELE_es;   //!
    TBranch        *b_ELE_fbrem;   //!
    TBranch        *b_psize;   //!
    TBranch        *b_pSC_energy;   //!
    TBranch        *b_pEta;   //!
    TBranch        *b_pPhi;   //!
    TBranch        *b_pEsenergy;   //!
    TBranch        *b_pRawenergy;   //!
    TBranch        *b_pPt;   //!
    TBranch        *b_p1x5_energy;   //!
    TBranch        *b_p2x5_energy;   //!
    TBranch        *b_p3x3_energy;   //!
    TBranch        *b_p5x5_energy;   //!
    TBranch        *b_pR9;   //!
    TBranch        *b_pGap;   //!
    TBranch        *b_pSigmaIetaIeta;   //!
    TBranch        *b_pHoverE;   //!
    TBranch        *b_p_photonenergy;   //!
    TBranch        *b_p_nCrystals;   //!
    TBranch        *b_p_adc;   //!
    TBranch        *b_p_gainId;   //!
    TBranch        *b_p_intercalib;   //!
    TBranch        *b_p_recieta;   //!
    TBranch        *b_p_reciphi;   //!
    TBranch        *b_p_receta;   //!
    TBranch        *b_p_recphi;   //!
    TBranch        *b_p_recenergy;   //!
    TBranch        *b_p_rectime;   //!
    TBranch        *b_p_recflag;   //!
    TBranch        *b_p_recflags;   //!
    TBranch        *b_NGEN;   //!
    TBranch        *b_GEN_pt;   //!
    TBranch        *b_GEN_eta;   //!
    TBranch        *b_GEN_phi;   //!
    TBranch        *b_GEN_px;   //!
    TBranch        *b_GEN_py;   //!
    TBranch        *b_GEN_pz;   //!
    TBranch        *b_GEN_e;   //!
    TBranch        *b_GEN_mass;   //!
    TBranch        *b_GEN_q;   //!
    TBranch        *b_GEN_id;   //!
    TBranch        *b_GEN_status;   //!
    TBranch        *b_GEN_parent;   //!
    TBranch        *b_ele1_charge;   //!
    TBranch        *b_ele1_p;   //!
    TBranch        *b_ele1_pt;   //!
    TBranch        *b_ele1_eta;   //!
    TBranch        *b_ele1_phi;   //!
    TBranch        *b_ele1_isEB;   //!
    TBranch        *b_ele1_isEBEEGap;   //!
    TBranch        *b_ele1_isEBEtaGap;   //!
    TBranch        *b_ele1_isEBPhiGap;   //!
    TBranch        *b_ele1_isEEDeeGap;   //!
    TBranch        *b_ele1_isEERingGap;   //!
    TBranch        *b_ele1_isTrackerDriven;   //!
    TBranch        *b_ele1_sigmaIetaIeta;   //!
    TBranch        *b_ele1_DphiIn;   //!
    TBranch        *b_ele1_DetaIn;   //!
    TBranch        *b_ele1_HOverE;   //!
    TBranch        *b_ele1_tkIso;   //!
    TBranch        *b_ele1_emIso;   //!
    TBranch        *b_ele1_hadIso;   //!
    TBranch        *b_ele1_scERaw;   //!
    TBranch        *b_ele1_scEtRaw;   //!
    TBranch        *b_ele1_scE;   //!
    TBranch        *b_ele1_scEt;   //!
    TBranch        *b_ele1_scE_regression;   //!
    TBranch        *b_ele1_scEerr_regression;   //!
    TBranch        *b_ele1_scE_regression_PhotonTuned;   //!
    TBranch        *b_ele1_scEerr_regression_PhotonTuned;   //!
    TBranch        *b_ele1_scERaw_PUcleaned;   //!
    TBranch        *b_ele1_es;   //!
    TBranch        *b_ele1_scLaserCorr;   //!
    TBranch        *b_ele1_scCrackCorr;   //!
    TBranch        *b_ele1_scLocalContCorr;   //!
    TBranch        *b_ele1_scEta;   //!
    TBranch        *b_ele1_scPhi;   //!
    TBranch        *b_ele1_scLocalEta;   //!
    TBranch        *b_ele1_scLocalPhi;   //!
    TBranch        *b_ele1_scEtaWidth;   //!
    TBranch        *b_ele1_scPhiWidth;   //!
    TBranch        *b_ele1_scEtaWidth_PUcleaned;   //!
    TBranch        *b_ele1_scPhiWidth_PUcleaned;   //!
    TBranch        *b_ele1_fCorrection_PUcleaned;   //!
    TBranch        *b_ele1_fEta;   //!
    TBranch        *b_ele1_fEtaCorr;   //!
    TBranch        *b_ele1_tkP;   //!
    TBranch        *b_ele1_tkPt;   //!
    TBranch        *b_ele1_fbrem;   //!
    TBranch        *b_ele1_dxy_PV;   //!
    TBranch        *b_ele1_dz_PV;   //!
    TBranch        *b_ele1_sigmaP;   //!
    TBranch        *b_ele1_eSeedBC;   //!
    TBranch        *b_ele1_e5x5;   //!
    TBranch        *b_ele1_e3x3;   //!
    TBranch        *b_ele1_scNxtal;   //!
    TBranch        *b_ele1_bcN;   //!
    TBranch        *b_ele1_5x5LaserCorr;   //!
    TBranch        *b_ele1_3x3LaserCorr;   //!
    TBranch        *b_ele1_seedE;   //!
    TBranch        *b_ele1_seedLaserAlpha;   //!
    TBranch        *b_ele1_seedLaserCorr;   //!
    TBranch        *b_ele1_seedICConstant;   //!
    TBranch        *b_ele1_seedIeta;   //!
    TBranch        *b_ele1_seedIphi;   //!
    TBranch        *b_ele1_seedIx;   //!
    TBranch        *b_ele1_seedIy;   //!
    TBranch        *b_ele1_seedZside;   //!
    TBranch        *b_ele1_EOverP;   //!
    TBranch        *b_ele1_nRecHits;   //!
    TBranch        *b_ele1_recHit_E;   //!
    TBranch        *b_ele1_recHit_flag;   //!
    TBranch        *b_ele1_recHit_hashedIndex;   //!
    TBranch        *b_ele1_recHit_ietaORix;   //!
    TBranch        *b_ele1_recHit_iphiORiy;   //!
    TBranch        *b_ele1_recHit_zside;   //!
    TBranch        *b_ele1_recHit_laserCorrection;   //!
    TBranch        *b_ele1_recHit_Alpha;   //!
    TBranch        *b_ele1_recHit_ICConstant;   //!
    TBranch        *b_ele2_charge;   //!
    TBranch        *b_ele2_p;   //!
    TBranch        *b_ele2_pt;   //!
    TBranch        *b_ele2_eta;   //!
    TBranch        *b_ele2_phi;   //!
    TBranch        *b_ele2_isEB;   //!
    TBranch        *b_ele2_isEBEEGap;   //!
    TBranch        *b_ele2_isEBEtaGap;   //!
    TBranch        *b_ele2_isEBPhiGap;   //!
    TBranch        *b_ele2_isEEDeeGap;   //!
    TBranch        *b_ele2_isEERingGap;   //!
    TBranch        *b_ele2_isTrackerDriven;   //!
    TBranch        *b_ele2_sigmaIetaIeta;   //!
    TBranch        *b_ele2_DphiIn;   //!
    TBranch        *b_ele2_DetaIn;   //!
    TBranch        *b_ele2_HOverE;   //!
    TBranch        *b_ele2_tkIso;   //!
    TBranch        *b_ele2_emIso;   //!
    TBranch        *b_ele2_hadIso;   //!
    TBranch        *b_ele2_dxy_PV;   //!
    TBranch        *b_ele2_dz_PV;   //!
    TBranch        *b_ele2_sigmaP;   //!
    TBranch        *b_ele2_scERaw;   //!
    TBranch        *b_ele2_scEtRaw;   //!
    TBranch        *b_ele2_scE;   //!
    TBranch        *b_ele2_scEt;   //!
    TBranch        *b_ele2_scE_regression;   //!
    TBranch        *b_ele2_scEerr_regression;   //!
    TBranch        *b_ele2_scE_regression_PhotonTuned;   //!
    TBranch        *b_ele2_scEerr_regression_PhotonTuned;   //!
    TBranch        *b_ele2_scERaw_PUcleaned;   //!
    TBranch        *b_ele2_es;   //!
    TBranch        *b_ele2_scLaserCorr;   //!
    TBranch        *b_ele2_scCrackCorr;   //!
    TBranch        *b_ele2_scLocalContCorr;   //!
    TBranch        *b_ele2_scEta;   //!
    TBranch        *b_ele2_scPhi;   //!
    TBranch        *b_ele2_scLocalEta;   //!
    TBranch        *b_ele2_scLocalPhi;   //!
    TBranch        *b_ele2_scEtaWidth;   //!
    TBranch        *b_ele2_scPhiWidth;   //!
    TBranch        *b_ele2_scEtaWidth_PUcleaned;   //!
    TBranch        *b_ele2_scPhiWidth_PUcleaned;   //!
    TBranch        *b_ele2_fCorrection_PUcleaned;   //!
    TBranch        *b_ele2_fEta;   //!
    TBranch        *b_ele2_fEtaCorr;   //!
    TBranch        *b_ele2_tkP;   //!
    TBranch        *b_ele2_tkPt;   //!
    TBranch        *b_ele2_fbrem;   //!
    TBranch        *b_ele2_eSeedBC;   //!
    TBranch        *b_ele2_e5x5;   //!
    TBranch        *b_ele2_e3x3;   //!
    TBranch        *b_ele2_scNxtal;   //!
    TBranch        *b_ele2_bcN;   //!
    TBranch        *b_ele2_5x5LaserCorr;   //!
    TBranch        *b_ele2_3x3LaserCorr;   //!
    TBranch        *b_ele2_seedE;   //!
    TBranch        *b_ele2_seedLaserAlpha;   //!
    TBranch        *b_ele2_seedLaserCorr;   //!
    TBranch        *b_ele2_seedICConstant;   //!
    TBranch        *b_ele2_seedIeta;   //!
    TBranch        *b_ele2_seedIphi;   //!
    TBranch        *b_ele2_seedIx;   //!
    TBranch        *b_ele2_seedIy;   //!
    TBranch        *b_ele2_seedZside;   //!
    TBranch        *b_ele2_EOverP;   //!
    TBranch        *b_ele2_nRecHits;   //!
    TBranch        *b_ele2_recHit_E;   //!
    TBranch        *b_ele2_recHit_flag;   //!
    TBranch        *b_ele2_recHit_hashedIndex;   //!
    TBranch        *b_ele2_recHit_ietaORix;   //!
    TBranch        *b_ele2_recHit_iphiORiy;   //!
    TBranch        *b_ele2_recHit_zside;   //!
    TBranch        *b_ele2_recHit_laserCorrection;   //!
    TBranch        *b_ele2_recHit_Alpha;   //!
    TBranch        *b_ele2_recHit_ICConstant;   //!
    TBranch        *b_met_et;   //!
    TBranch        *b_met_phi;   //!
    TBranch        *b_ele1Met_mt;   //!
    TBranch        *b_ele1Met_Dphi;   //!
    TBranch        *b_ele1ele2_m;   //!
    TBranch        *b_ele1ele2_scM;   //!
    TBranch        *b_ele1ele2_scM_regression;   //!
    TBranch        *b_PUit_TrueNumInteractions;   //!
    TBranch        *b_PUit_NumInteractions;   //!
    TBranch        *b_PUit_zpositions;   //!
    TBranch        *b_PUit_sumpT_lowpT;   //!
    TBranch        *b_PUit_sumpT_highpT;   //!
    TBranch        *b_PUit_ntrks_lowpT;   //!
    TBranch        *b_PUit_ntrks_highpT;   //!
    TBranch        *b_PUoot_early_TrueNumInteractions;   //!
    TBranch        *b_PUoot_early;   //!
    TBranch        *b_PUoot_early_zpositions;   //!
    TBranch        *b_PUoot_early_sumpT_lowpT;   //!
    TBranch        *b_PUoot_early_sumpT_highpT;   //!
    TBranch        *b_PUoot_early_ntrks_lowpT;   //!
    TBranch        *b_PUoot_early_ntrks_highpT;   //!
    TBranch        *b_PUoot_late_TrueNumInteractions;   //!
    TBranch        *b_PUoot_late;   //!
    TBranch        *b_PUoot_late_zpositions;   //!
    TBranch        *b_PUoot_late_sumpT_lowpT;   //!
    TBranch        *b_PUoot_late_sumpT_highpT;   //!
    TBranch        *b_PUoot_late_ntrks_lowpT;   //!
    TBranch        *b_PUoot_late_ntrks_highpT;   //!
    TBranch        *b_mcV_E;   //!
    TBranch        *b_mcV_Px;   //!
    TBranch        *b_mcV_Py;   //!
    TBranch        *b_mcV_Pz;   //!
    TBranch        *b_mcV_Charge;   //!
    TBranch        *b_mcV_PdgId;   //!
    TBranch        *b_mcF1_fromV_E;   //!
    TBranch        *b_mcF1_fromV_Px;   //!
    TBranch        *b_mcF1_fromV_Py;   //!
    TBranch        *b_mcF1_fromV_Pz;   //!
    TBranch        *b_mcF1_fromV_Charge;   //!
    TBranch        *b_mcF1_fromV_PdgId;   //!
    TBranch        *b_mcF2_fromV_E;   //!
    TBranch        *b_mcF2_fromV_Px;   //!
    TBranch        *b_mcF2_fromV_Py;   //!
    TBranch        *b_mcF2_fromV_Pz;   //!
    TBranch        *b_mcF2_fromV_Charge;   //!
    TBranch        *b_mcF2_fromV_PdgId;   //!
    
    analysis(TTree *tree=0);
    virtual ~analysis();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual void     LoopPho();
    virtual void     LoopEle();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
};

#endif

