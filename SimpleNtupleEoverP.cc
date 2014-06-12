// -*- C++ -*-
//
// Package:   SimpleNtupleEoverP
// Class:     SimpleNtupleEoverP
//
#include "Calibration/EcalCalibNtuple/plugins/SimpleNtupleEoverP.h"
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
#include "PhysicsTools/NtupleUtils/interface/readJSONFile.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsMCRcd.h"


#define PI 3.141592654
#define TWOPI 6.283185308


// MY STUFF
const int NELE_max = 100;
int NELE=-99;
float ELE_pt[NELE_max];
float ELE_eta[NELE_max];
float ELE_phi[NELE_max];
float ELE_px[NELE_max];
float ELE_py[NELE_max];
float ELE_pz[NELE_max];
float ELE_e[NELE_max];
int ELE_q[NELE_max];
int ELE_id[NELE_max];
float ELE_sigmaIetaIeta[NELE_max];
float ELE_DphiIn[NELE_max];
float ELE_DetaIn[NELE_max];
float ELE_HOverE[NELE_max];
float ELE_ooemoop[NELE_max];
float ELE_tkIso[NELE_max];
float ELE_emIso[NELE_max];
float ELE_hadIso[NELE_max];
float ELE_effAreaForIso[NELE_max];
float ELE_combIso[NELE_max];
float ELE_dxy[NELE_max];
float ELE_dz[NELE_max];
float ELE_scE_regression[NELE_max];
float ELE_scEtaWidth[NELE_max];
float ELE_scPhiWidth[NELE_max];
float ELE_scERaw[NELE_max];
float ELE_scE[NELE_max];
float ELE_e5x5[NELE_max];
float ELE_e1x5[NELE_max];
float ELE_e2x5Max[NELE_max];
float ELE_es[NELE_max];
float ELE_fbrem[NELE_max];
int ELE_adc[NELE_max][25][10]; 
int ELE_gainId[NELE_max][25][10];
float ELE_intercalib[NELE_max][25];
int ELE_recieta[NELE_max][25];
int ELE_reciphi[NELE_max][25];
float ELE_receta[NELE_max][25];
float ELE_recphi[NELE_max][25];
float  ELE_recenergy[NELE_max][25];
float  ELE_rectime[NELE_max][25];
uint32_t ELE_recflag[NELE_max][25];
bool ELE_recflags[NELE_max][25][19];

const int doPhotons=0;
const int NPHO_max=1;//100;
int pho_size;
double p_esenergy[NPHO_max];
double p_rawenergy[NPHO_max];
float p_ecalenergy[NPHO_max];
float p_energy_1x5[NPHO_max];
float p_energy_2x5[NPHO_max];
float p_energy_3x3[NPHO_max];
float p_energy_5x5[NPHO_max];
float p_pt[NPHO_max];
float p_r9[NPHO_max];
bool p_isGap[NPHO_max];
float p_sigmaIetaIeta[NPHO_max];
float p_hOverE[NPHO_max];
double p_eta[NPHO_max];
double p_phi[NPHO_max];
//-----------------------------------
int p_adc[NPHO_max][25][10];
int p_gainId[NPHO_max][25][10];
float p_intercalib[NPHO_max][25];
int p_recieta[NPHO_max][25];
int p_reciphi[NPHO_max][25];
float p_receta[NPHO_max][25];
float p_recphi[NPHO_max][25];
float  p_recenergy[NPHO_max][25];
float  p_rectime[NPHO_max][25];
uint32_t p_recflag[NPHO_max][25];
bool p_recflags[NPHO_max][25][19];
float p_photonenergy[NPHO_max];
int p_nCrystals[NPHO_max];

const int NGEN_max = 1000;
int NGEN=-99;
float GEN_pt[NGEN_max];
float GEN_eta[NGEN_max];
float GEN_phi[NGEN_max];
float GEN_px[NGEN_max];
float GEN_py[NGEN_max];
float GEN_pz[NGEN_max];
float GEN_e[NGEN_max];
float GEN_mass[NGEN_max];
int GEN_q[NGEN_max];
int GEN_id[NGEN_max];
int GEN_status[NGEN_max];
int GEN_parent[NGEN_max];

// ---------------------------------------

SimpleNtupleEoverP::SimpleNtupleEoverP(const edm::ParameterSet& iConfig)
{
	std::cout<< ">>> SimpleNtupleEoverP::SimpleNtupleEoverP begin <<<" << std::endl;
	
	
	
	// Initialize parameters for cluster corrections
	std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::read parameters <<<" << std::endl;
	InitializeParams(params);
	
	edm::Service<TFileService> fs ;
	outTree_  = fs -> make <TTree>("SimpleNtupleEoverP","SimpleNtupleEoverP"); 
	
	eventType_    = iConfig.getUntrackedParameter<int>("eventType",1);
	dataRun_      = iConfig.getParameter<std::string>("dataRun");
	jsonFileName_ = iConfig.getParameter<std::string>("jsonFileName");
	digiCollection_EB_   = iConfig.getParameter<edm::InputTag>("digiCollection_EB");
	digiCollection_EE_   = iConfig.getParameter<edm::InputTag>("digiCollection_EE");
	recHitCollection_EB_ = iConfig.getParameter<edm::InputTag>("recHitCollection_EB");
	recHitCollection_EE_ = iConfig.getParameter<edm::InputTag>("recHitCollection_EE");
	BSTag_               = iConfig.getParameter<edm::InputTag>("theBeamSpotTag");
	SRFlagCollection_EB_ = iConfig.getParameter<edm::InputTag>("SRFlagCollection_EB");
	SRFlagCollection_EE_ = iConfig.getParameter<edm::InputTag>("SRFlagCollection_EE");
	MCPileupTag_         = iConfig.getParameter<edm::InputTag>("MCPileupTag");
	MCtruthTag_          = iConfig.getParameter<edm::InputTag>("MCtruthTag");
	PVTag_               = iConfig.getParameter<edm::InputTag>("PVTag");
	rhoTag_              = iConfig.getParameter<edm::InputTag>("rhoTag");
	//EleTag_              = iConfig.getParameter<edm::InputTag>("EleTag");
    //EleTag_ = edm::InputTag("gsfElectrons","","reRECO");
    EleTag_ = edm::InputTag("gsfElectrons","","RECO");
	PFMetTag_            = iConfig.getParameter<edm::InputTag>("PFMetTag");
	conversionsInputTag_ = iConfig.getParameter<edm::InputTag>("conversionsInputTag");
	
	
	//---- flags ----
	jsonFlag_         = iConfig.getUntrackedParameter<bool>("jsonFlag",        false);
	saveMCPU_         = iConfig.getUntrackedParameter<bool>("saveMCPU",        false);
	verbosity_        = iConfig.getUntrackedParameter<bool>("verbosity",       false);
	doWZSelection_    = iConfig.getUntrackedParameter<bool>("doWZSelection",   false);
	applyCorrections_ = iConfig.getUntrackedParameter<bool>("applyCorrections",false);
	dataFlag_         = iConfig.getUntrackedParameter<bool>("dataFlag",         true);
	saveRecHitMatrix_ = iConfig.getUntrackedParameter<bool>("saveRecHitMatrix",false);
	saveFbrem_        = iConfig.getUntrackedParameter<bool>("saveFbrem",       false);  
	saveMCInfo_       = iConfig.getUntrackedParameter<bool>("saveMCInfo",      false);  
	eventNaiveId_     = 0;
	mcAnalysisZW_     = NULL;
	
	
	//---- Initialize tree branches ----
	
	// event variables
	std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::set event branches <<<" << std::endl;
	outTree_ -> Branch("bxId",         &bxId,                  "bxId/L");
	outTree_ -> Branch("eventId",      &eventId,            "eventId/L");
	outTree_ -> Branch("lumiId",       &lumiId,              "lumiId/I");
	outTree_ -> Branch("runId",        &runId,                "runId/I");
	outTree_ -> Branch("timeStampHigh",&timeStampHigh,"timeStampHigh/I");
	outTree_ -> Branch("isW", &isW, "isW/I");
	outTree_ -> Branch("isZ", &isZ, "isZ/I");  
	
	outTree_ -> Branch("PV_n",    &PV_n,        "PV_n/I");
	outTree_ -> Branch("PV_z",    &PV_z,        "PV_z/F");
	outTree_ -> Branch("PV_d0",   &PV_d0,      "PV_d0/F");
	outTree_ -> Branch("PV_SumPt",&PV_SumPt,"PV_SumPt/F");
	
	outTree_ -> Branch("rho",   &rho,  "rho/F");
	
	
	// MY STUFF:
	for (int i=0; i<NELE_max; i++) {
		ELE_pt[i] = -99;
		ELE_eta[i] = -99;
		ELE_phi[i] = -99;
		ELE_px[i] = -9999;
		ELE_py[i] = -9999;
		ELE_pz[i] = -9999;
		ELE_e[i] = -9999;
		ELE_q[i] = -99;
		ELE_id[i] = -99;
		ELE_sigmaIetaIeta[i]= -9999;
		ELE_DphiIn[i]= -9999;
		ELE_DetaIn[i]= -9999;
		ELE_HOverE[i]= -9999;
		ELE_ooemoop[i]= -9999;
		ELE_tkIso[i]= -9999;
		ELE_emIso[i]= -9999;
		ELE_hadIso[i]= -9999;
		ELE_effAreaForIso[i]= -9999;
		ELE_combIso[i]= -9999;
		ELE_dxy[i]= -9999;
		ELE_dz[i]= -9999;
		ELE_scE_regression[i]= -9999;
		ELE_scEtaWidth[i]= -9999;
		ELE_scPhiWidth[i]= -9999;
		ELE_scERaw[i]= -9999;
		ELE_scE[i]= -9999;
		ELE_e5x5[i]= -9999;
		ELE_e1x5[i]= -9999;
		ELE_e2x5Max[i]= -9999;
		ELE_es[i]= -9999;
		ELE_fbrem[i]= -9999;
		for (int j=0; j<25; j++) {
            ELE_intercalib[i][j] = -9999;
            ELE_recieta[i][j] = -9999;
            ELE_reciphi[i][j] = -9999;
            ELE_receta[i][j] = -9999;
            ELE_recphi[i][j] = -9999;
            ELE_recenergy[i][j] = -9999;
            ELE_rectime[i][j] = -9999;
            ELE_recflag[i][j] = -9999;
			for (int k=0; k<10; k++) {
				ELE_adc[i][j][k] = -9999;
				ELE_gainId[i][j][k] = -9999;
			}
            for (int k=0; k<19; k++) {
                ELE_recflags[i][j][k] = -9999;
			}
		}
	}
	std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::set MY STUFF branches (ELE) <<<" << std::endl;
	outTree_ -> Branch("NELE",&NELE,"NELE/I");
	outTree_ -> Branch("ELE_pt",&ELE_pt,"ELE_pt[NELE]/F");
	outTree_ -> Branch("ELE_eta",&ELE_eta,"ELE_eta[NELE]/F");
	outTree_ -> Branch("ELE_phi",&ELE_phi,"ELE_phi[NELE]/F");
	outTree_ -> Branch("ELE_px",&ELE_px,"ELE_px[NELE]/F");
	outTree_ -> Branch("ELE_py",&ELE_py,"ELE_py[NELE]/F");
	outTree_ -> Branch("ELE_pz",&ELE_pz,"ELE_pz[NELE]/F");
	outTree_ -> Branch("ELE_e",&ELE_e,"ELE_e[NELE]/F");
	outTree_ -> Branch("ELE_q",&ELE_q,"ELE_q[NELE]/I");
	outTree_ -> Branch("ELE_id",&ELE_id,"ELE_id[NELE]/I");

	outTree_ -> Branch("ELE_sigmaIetaIeta",&ELE_sigmaIetaIeta,"ELE_sigmaIetaIeta[NELE]/F");
	outTree_ -> Branch("ELE_DphiIn",&ELE_DphiIn,"ELE_DphiIn[NELE]/F");
	outTree_ -> Branch("ELE_DetaIn",&ELE_DetaIn,"ELE_DetaIn[NELE]/F");
	outTree_ -> Branch("ELE_HOverE",&ELE_HOverE,"ELE_HOverE[NELE]/F");
	outTree_ -> Branch("ELE_ooemoop",&ELE_ooemoop,"ELE_ooemoop[NELE]/F");
	outTree_ -> Branch("ELE_tkIso",&ELE_tkIso,"ELE_tkIso[NELE]/F");
	outTree_ -> Branch("ELE_emIso",&ELE_emIso,"ELE_emIso[NELE]/F");
	outTree_ -> Branch("ELE_hadIso",&ELE_hadIso,"ELE_hadIso[NELE]/F");
	outTree_ -> Branch("ELE_effAreaForIso",&ELE_effAreaForIso,"ELE_effAreaForIso[NELE]/F");
	outTree_ -> Branch("ELE_combIso",&ELE_combIso,"ELE_combIso[NELE]/F");
	outTree_ -> Branch("ELE_dxy",&ELE_dxy,"ELE_dxy[NELE]/F");
	outTree_ -> Branch("ELE_dz",&ELE_dz,"ELE_dz[NELE]/F");
	outTree_ -> Branch("ELE_scE_regression",&ELE_scE_regression,"ELE_scE_regression[NELE]/F");
	outTree_ -> Branch("ELE_scEtaWidth",&ELE_scEtaWidth,"ELE_scEtaWidth[NELE]/F");
	outTree_ -> Branch("ELE_scPhiWidth",&ELE_scPhiWidth,"ELE_scPhiWidth[NELE]/F");
	outTree_ -> Branch("ELE_scERaw",&ELE_scERaw,"ELE_scERaw[NELE]/F");
	outTree_ -> Branch("ELE_scE",&ELE_scE,"ELE_scE[NELE]/F");
	outTree_ -> Branch("ELE_e5x5",&ELE_e5x5,"ELE_e5x5[NELE]/F");
	outTree_ -> Branch("ELE_e1x5",&ELE_e1x5,"ELE_e1x5[NELE]/F");
	outTree_ -> Branch("ELE_e2x5Max",&ELE_e2x5Max,"ELE_e2x5Max[NELE]/F");
	outTree_ -> Branch("ELE_es",&ELE_es,"ELE_es[NELE]/F");
	outTree_ -> Branch("ELE_fbrem",&ELE_fbrem,"ELE_fbrem[NELE]/F");

    /*
     // decreasing the size, commenting out for H->gg samples
	outTree_ -> Branch("ELE_adc",ELE_adc,"ELE_adc[NELE][25][10]/I");
	outTree_ -> Branch("ELE_gainId",ELE_gainId,"ELE_gainId[NELE][25][10]/I");
    outTree_ -> Branch("ELE_intercalib",ELE_intercalib,"ELE_intercalib[NELE][25]/F");
    outTree_ -> Branch("ELE_recieta",ELE_recieta,"ELE_recieta[NELE][25]/I");
    outTree_ -> Branch("ELE_reciphi",ELE_reciphi,"ELE_reciphi[NELE][25]/I");
    outTree_ -> Branch("ELE_receta",ELE_receta,"ELE_receta[NELE][25]/F");
    outTree_ -> Branch("ELE_recphi",ELE_recphi,"ELE_recphi[NELE][25]/F");
    outTree_ -> Branch("ELE_recenergy",ELE_recenergy,"ELE_recenergy[NELE][25]/F");
    outTree_ -> Branch("ELE_rectime",ELE_rectime,"ELE_rectime[NELE][25]/F");
    outTree_ -> Branch("ELE_recflag",ELE_recflag,"ELE_recflag[NELE][25]/i");
    outTree_ -> Branch("ELE_recflags",ELE_recflags,"ELE_recflags[NELE][25][19]/O"); */


    outTree_->Branch("psize",&pho_size,"psize/I");
    outTree_->Branch("pSC_energy",p_ecalenergy,"pSC_energy[psize]/F");
    outTree_->Branch("pEta",p_eta,"pEta[psize]/D");
    outTree_->Branch("pPhi",p_phi,"pPhi[psize]/D");
    outTree_->Branch("pEsenergy",p_esenergy,"pEsenergy[psize]/D");
    outTree_->Branch("pRawenergy",p_rawenergy,"pRawenergy[psize]/D");
    outTree_->Branch("pPt",p_pt,"pPt[psize]/F");
    outTree_->Branch("p1x5_energy",p_energy_1x5,"p1x5_energy[psize]/F");
    outTree_->Branch("p2x5_energy",p_energy_2x5,"p2x5_energy[psize]/F");
    outTree_->Branch("p3x3_energy",p_energy_3x3,"p3x3_energy[psize]/F");
    outTree_->Branch("p5x5_energy",p_energy_5x5,"p5x5_energy[psize]/F");
    outTree_->Branch("pR9",p_r9,"pR9[psize]/F");
    outTree_->Branch("pGap",p_isGap,"pGap[psize]/O");
    outTree_->Branch("pSigmaIetaIeta",p_sigmaIetaIeta,"pSigmaIetaIeta[psize]/F");
    outTree_->Branch("pHoverE",p_hOverE,"pHoverE[psize]/F");

    outTree_->Branch("p_photonenergy",p_photonenergy,"p_photonenergy[psize]/F");
    outTree_->Branch("p_nCrystals",p_nCrystals,"p_nCrystals[psize]/I");
    
	outTree_ -> Branch("p_adc",p_adc,"p_adc[psize][25][10]/I");
	outTree_ -> Branch("p_gainId",p_gainId,"p_gainId[psize][25][10]/I");
    outTree_ -> Branch("p_intercalib",p_intercalib,"p_intercalib[psize][25]/F");
    outTree_ -> Branch("p_recieta",p_recieta,"p_recieta[psize][25]/I");
    outTree_ -> Branch("p_reciphi",p_reciphi,"p_reciphi[psize][25]/I");
    outTree_ -> Branch("p_receta",p_receta,"p_receta[psize][25]/F");
    outTree_ -> Branch("p_recphi",p_recphi,"p_recphi[psize][25]/F");
    outTree_ -> Branch("p_recenergy",p_recenergy,"p_recenergy[psize][25]/F");
    outTree_ -> Branch("p_rectime",p_rectime,"p_rectime[psize][25]/F");
    outTree_ -> Branch("p_recflag",p_recflag,"p_recflag[psize][25]/i");
    outTree_ -> Branch("p_recflags",p_recflags,"p_recflags[psize][25][19]/O");

    for (int i=0; i<NPHO_max; i++) {
		p_ecalenergy[i] = -99;
		p_eta[i] = -99;
		p_phi[i] = -99;
		p_esenergy[i] = -9999;
		p_rawenergy[i] = -9999;
		p_pt[i] = -9999;
		p_energy_1x5[i] = -9999;
		p_energy_2x5[i] = -99;
		p_energy_3x3[i] = -99;
		p_energy_5x5[i]= -9999;
		p_r9[i]= -9999;
		p_isGap[i]= -9999;
		p_sigmaIetaIeta[i]= -9999;
		p_hOverE[i]= -9999;
		p_photonenergy[i]= -9999;
		p_nCrystals[i]= -9999;
		
		for (int j=0; j<25; j++) {
            p_intercalib[i][j] = -9999;
            p_recieta[i][j] = -9999;
            p_reciphi[i][j] = -9999;
            p_receta[i][j] = -9999;
            p_recphi[i][j] = -9999;
            p_recenergy[i][j] = -9999;
            p_rectime[i][j] = -9999;
            p_recflag[i][j] = -9999;
			for (int k=0; k<10; k++) {
				p_adc[i][j][k] = -9999;
				p_gainId[i][j][k] = -9999;
			}
            for (int k=0; k<19; k++) {
                p_recflags[i][j][k] = -9999;
			}
		}
	}

    

	// GEN PARTICLES
	for (int i=0; i<NGEN_max; i++) {
		GEN_pt[i] = -99;
		GEN_eta[i] = -99;
		GEN_phi[i] = -99;
		GEN_px[i] = -9999;
		GEN_py[i] = -9999;
		GEN_pz[i] = -9999;
		GEN_e[i] = -9999;
		GEN_mass[i] = -9999;
		GEN_q[i] = -99;
		GEN_id[i] = -99;
		GEN_status[i] = -99;
		GEN_parent[i] = -99;
	}
	std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::set MY STUFF branches (GEN) <<<" << std::endl;
	outTree_ -> Branch("NGEN",&NGEN,"NGEN/I");
	outTree_ -> Branch("GEN_pt",&GEN_pt,"GEN_pt[NGEN]/F");
	outTree_ -> Branch("GEN_eta",&GEN_eta,"GEN_eta[NGEN]/F");
	outTree_ -> Branch("GEN_phi",&GEN_phi,"GEN_phi[NGEN]/F");
	outTree_ -> Branch("GEN_px",&GEN_px,"GEN_px[NGEN]/F");
	outTree_ -> Branch("GEN_py",&GEN_py,"GEN_py[NGEN]/F");
	outTree_ -> Branch("GEN_pz",&GEN_pz,"GEN_pz[NGEN]/F");
	outTree_ -> Branch("GEN_e",&GEN_e,"GEN_e[NGEN]/F");
	outTree_ -> Branch("GEN_mass",&GEN_mass,"GEN_mass[NGEN]/F");
	outTree_ -> Branch("GEN_q",&GEN_q,"GEN_q[NGEN]/I");
	outTree_ -> Branch("GEN_id",&GEN_id,"GEN_id[NGEN]/I");
	outTree_ -> Branch("GEN_status",&GEN_status,"GEN_status[NGEN]/I");
	outTree_ -> Branch("GEN_parent",&GEN_parent,"GEN_parent[NGEN]/I");
	//--------------------------------------------------------------
	
	
	
	// ele1 variables
	std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::set ele1 branches <<<" << std::endl;
	outTree_ -> Branch("ele1_charge",&ele1_charge,"ele1_charge/F");
	outTree_ -> Branch("ele1_p",     &ele1_p,          "ele1_p/F");
	outTree_ -> Branch("ele1_pt",    &ele1_pt,        "ele1_pt/F");
	outTree_ -> Branch("ele1_eta",   &ele1_eta,      "ele1_eta/F");
	outTree_ -> Branch("ele1_phi",   &ele1_phi,      "ele1_phi/F");
	
	outTree_ -> Branch("ele1_isEB",       &ele1_isEB,              "ele1_isEB/I");
	outTree_ -> Branch("ele1_isEBEEGap",  &ele1_isEBEEGap,    "ele1_isEBEEGap/I");
	outTree_ -> Branch("ele1_isEBEtaGap", &ele1_isEBEtaGap,  "ele1_isEBEtaGap/I");
	outTree_ -> Branch("ele1_isEBPhiGap", &ele1_isEBPhiGap,  "ele1_isEBPhiGap/I");
	outTree_ -> Branch("ele1_isEEDeeGap", &ele1_isEEDeeGap,  "ele1_isEEDeeGap/I");
	outTree_ -> Branch("ele1_isEERingGap",&ele1_isEERingGap,"ele1_isEERingGap/I");
	
	outTree_ -> Branch("ele1_isTrackerDriven",&ele1_isTrackerDriven,"ele1_isTrackerDriven/I");
	outTree_ -> Branch("ele1_sigmaIetaIeta",  &ele1_sigmaIetaIeta,    "ele1_sigmaIetaIeta/F");
	outTree_ -> Branch("ele1_DphiIn",         &ele1_DphiIn,                  "ele1_DphiIn/F");
	outTree_ -> Branch("ele1_DetaIn",         &ele1_DetaIn,                  "ele1_DetaIn/F");
	outTree_ -> Branch("ele1_HOverE",         &ele1_HOverE,                  "ele1_HOverE/F");
	outTree_ -> Branch("ele1_tkIso",          &ele1_tkIso,                    "ele1_tkIso/F");
	outTree_ -> Branch("ele1_emIso",          &ele1_emIso,                    "ele1_emIso/F");
	outTree_ -> Branch("ele1_hadIso",         &ele1_hadIso,                  "ele1_hadIso/F");
	
	outTree_ -> Branch("ele1_scERaw",     &ele1_scERaw,         "ele1_scERaw/F");
	outTree_ -> Branch("ele1_scEtRaw",    &ele1_scEtRaw,       "ele1_scEtRaw/F");
	outTree_ -> Branch("ele1_scE",    &ele1_scE,         "ele1_scE/F");
	outTree_ -> Branch("ele1_scEt",       &ele1_scEt,             "ele1_scEt/F");
	outTree_ -> Branch("ele1_scE_regression",    &ele1_scE_regression,   "ele1_scE_regression/F");
	outTree_ -> Branch("ele1_scEerr_regression",       &ele1_scEerr_regression,         "ele1_scEerr_regression/F");
	outTree_ -> Branch("ele1_scE_regression_PhotonTuned",    &ele1_scE_regression_PhotonTuned ,   "ele1_scE_regression_PhotonTuned/F");
	outTree_ -> Branch("ele1_scEerr_regression_PhotonTuned",       &ele1_scEerr_regression_PhotonTuned ,         "ele1_scEerr_regression_PhotonTuned/F");
	outTree_ -> Branch("ele1_scERaw_PUcleaned",    &ele1_scERaw_PUcleaned,   "ele1_scERaw_PUcleaned/F");
	outTree_ -> Branch("ele1_es",       &ele1_es,  "ele1_es/F");
	outTree_ -> Branch("ele1_scLaserCorr",       &ele1_scLaserCorr,         "ele1_scLaserCorr/F");
	outTree_ -> Branch("ele1_scCrackCorr",   &ele1_scCrackCorr,        "ele1_scCrackCorr/F");
	outTree_ -> Branch("ele1_scLocalContCorr",    &ele1_scLocalContCorr,         "ele1_scLocalContCorr/F");
	
	outTree_ -> Branch("ele1_scEta",  &ele1_scEta,     "ele1_scEta/F");
	outTree_ -> Branch("ele1_scPhi",  &ele1_scPhi,     "ele1_scPhi/F");
	outTree_ -> Branch("ele1_scLocalEta",    &ele1_scLocalEta,         "ele1_scLocalEta/F");
	outTree_ -> Branch("ele1_scLocalPhi",    &ele1_scLocalPhi,         "ele1_scLocalPhi/F");
	outTree_ -> Branch("ele1_scEtaWidth",    &ele1_scEtaWidth,         "ele1_scEtaWidth/F");
	outTree_ -> Branch("ele1_scPhiWidth",    &ele1_scPhiWidth,         "ele1_scPhiWidth/F");
	outTree_ -> Branch("ele1_scEtaWidth_PUcleaned",       &ele1_scEtaWidth_PUcleaned,         "ele1_scEtaWidth_PUcleaned/F");
	outTree_ -> Branch("ele1_scPhiWidth_PUcleaned",    &ele1_scPhiWidth_PUcleaned,   "ele1_scPhiWidth_PUcleaned/F");
	outTree_ -> Branch("ele1_fCorrection_PUcleaned",       &ele1_fCorrection_PUcleaned,         "ele1_fCorrection_PUcleaned/F");
	
	outTree_ -> Branch("ele1_fEta",       &ele1_fEta,   "ele1_fEta/F");
	outTree_ -> Branch("ele1_fEtaCorr",   &ele1_fEtaCorr,         "ele1_fEtaCorr/F");
	outTree_ -> Branch("ele1_tkP",        &ele1_tkP,   "ele1_tkP/F");
	outTree_ -> Branch("ele1_tkPt",       &ele1_tkPt,  "ele1_tkPt/F");
	outTree_ -> Branch("ele1_fbrem",       &ele1_fbrem,  "ele1_fbrem/F");
	
	outTree_ -> Branch("ele1_dxy_PV",        &ele1_dxy_PV,   "ele1_dxy_PV/F");
	outTree_ -> Branch("ele1_dz_PV",       &ele1_dz_PV,  "ele1_dz_PV/F");
	outTree_ -> Branch("ele1_sigmaP",       &ele1_sigmaP,  "ele1_sigmaP/F");
	
	outTree_ -> Branch("ele1_eSeedBC",       &ele1_eSeedBC,  "ele1_eSeedBC/F");
	outTree_ -> Branch("ele1_e5x5",       &ele1_e5x5,  "ele1_e5x5/F");
	outTree_ -> Branch("ele1_e3x3",       &ele1_e3x3,  "ele1_e3x3/F");
	outTree_ -> Branch("ele1_scNxtal",       &ele1_scNxtal,  "ele1_scNxtal/F");
	outTree_ -> Branch("ele1_bcN",           &ele1_bcN,  "ele1_bcN/I");
	outTree_ -> Branch("ele1_5x5LaserCorr",       &ele1_5x5LaserCorr,  "ele1_5x5LaserCorr/F");
	outTree_ -> Branch("ele1_3x3LaserCorr",       &ele1_3x3LaserCorr,  "ele1_3x3LaserCorr/F");
	
	outTree_ -> Branch("ele1_seedE",       &ele1_seedE,  "ele1_seedE/F");
	outTree_ -> Branch("ele1_seedLaserAlpha",       &ele1_seedLaserAlpha,  "ele1_seedLaserAlpha/F");
	outTree_ -> Branch("ele1_seedLaserCorr",       &ele1_seedLaserCorr,  "ele1_seedLaserCorr/F");
	outTree_ -> Branch("ele1_seedICConstant",       &ele1_seedICConstant,  "ele1_seedICConstant/F");
	outTree_ -> Branch("ele1_seedIeta",       &ele1_seedIeta,  "ele1_seedIeta/I");
	outTree_ -> Branch("ele1_seedIphi",       &ele1_seedIphi,  "ele1_seedIphi/I");
	outTree_ -> Branch("ele1_seedIx",       &ele1_seedIx,  "ele1_seedIx/I");
	outTree_ -> Branch("ele1_seedIy",       &ele1_seedIy,  "ele1_seedIy/I");
	outTree_ -> Branch("ele1_seedZside",       &ele1_seedZside,  "ele1_seedZside/I");
	outTree_ -> Branch("ele1_EOverP",       &ele1_EOverP,  "ele1_EOverP/F");
	
	outTree_ -> Branch("ele1_nRecHits",&ele1_nRecHits,"ele1_nRecHits/I");
	outTree_ -> Branch("ele1_recHit_E",              "std::vector<float>",&ele1_recHit_E);
	outTree_ -> Branch("ele1_recHit_flag",           "std::vector<int>",  &ele1_recHit_flag);
	outTree_ -> Branch("ele1_recHit_hashedIndex",    "std::vector<int>",  &ele1_recHit_hashedIndex);
	outTree_ -> Branch("ele1_recHit_ietaORix",       "std::vector<int>",  &ele1_recHit_ietaORix);
	outTree_ -> Branch("ele1_recHit_iphiORiy",       "std::vector<int>",  &ele1_recHit_iphiORiy);
	outTree_ -> Branch("ele1_recHit_zside",          "std::vector<int>",  &ele1_recHit_zside);
	outTree_ -> Branch("ele1_recHit_laserCorrection","std::vector<float>",&ele1_recHit_laserCorrection);
	outTree_ -> Branch("ele1_recHit_Alpha",          "std::vector<float>",&ele1_recHit_Alpha);
	outTree_ -> Branch("ele1_recHit_ICConstant",     "std::vector<float>",&ele1_recHit_ICConstant);
	
	if(saveRecHitMatrix_)
	{
		outTree_ -> Branch("ele1_recHitMatrix_E",              "std::vector<float>",&ele1_recHitMatrix_E);
		outTree_ -> Branch("ele1_recHitMatrix_flag",           "std::vector<int>",  &ele1_recHitMatrix_flag);
		outTree_ -> Branch("ele1_recHitMatrix_hashedIndex",    "std::vector<int>",  &ele1_recHitMatrix_hashedIndex);
		outTree_ -> Branch("ele1_recHitMatrix_ietaORix",       "std::vector<int>",  &ele1_recHitMatrix_ietaORix);
		outTree_ -> Branch("ele1_recHitMatrix_iphiORiy",       "std::vector<int>",  &ele1_recHitMatrix_iphiORiy);
		outTree_ -> Branch("ele1_recHitMatrix_zside",          "std::vector<int>",  &ele1_recHitMatrix_zside);
		outTree_ -> Branch("ele1_recHitMatrix_laserCorrection","std::vector<float>",&ele1_recHitMatrix_laserCorrection);
		outTree_ -> Branch("ele1_recHitMatrix_ICConstant",     "std::vector<float>",&ele1_recHitMatrix_ICConstant);
		outTree_ -> Branch("ele1_recHitMatrix_samples",        "std::vector<float>",&ele1_recHitMatrix_samples);
	}
	
	
	// ele2 variables
	std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::set ele2 branches <<<" << std::endl;
	outTree_ -> Branch("ele2_charge",&ele2_charge,"ele2_charge/F");
	outTree_ -> Branch("ele2_p",     &ele2_p,          "ele2_p/F");
	outTree_ -> Branch("ele2_pt",    &ele2_pt,        "ele2_pt/F");
	outTree_ -> Branch("ele2_eta",   &ele2_eta,      "ele2_eta/F");
	outTree_ -> Branch("ele2_phi",   &ele2_phi,      "ele2_phi/F");
	
	outTree_ -> Branch("ele2_isEB",       &ele2_isEB,  "ele2_isEB/I");
	outTree_ -> Branch("ele2_isEBEEGap",       &ele2_isEBEEGap,  "ele2_isEBEEGap/I");
	outTree_ -> Branch("ele2_isEBEtaGap",       &ele2_isEBEtaGap,  "ele2_isEBEtaGap/I");
	outTree_ -> Branch("ele2_isEBPhiGap",       &ele2_isEBPhiGap,  "ele2_isEBPhiGap/I");
	outTree_ -> Branch("ele2_isEEDeeGap",       &ele2_isEEDeeGap,  "ele2_isEEDeeGap/I");
	outTree_ -> Branch("ele2_isEERingGap",       &ele2_isEERingGap,  "ele2_isEERingGap/I");
	
	outTree_ -> Branch("ele2_isTrackerDriven",&ele2_isTrackerDriven,"ele2_isTrackerDriven/I");
	outTree_ -> Branch("ele2_sigmaIetaIeta",  &ele2_sigmaIetaIeta,    "ele2_sigmaIetaIeta/F");
	outTree_ -> Branch("ele2_DphiIn",         &ele2_DphiIn,                  "ele2_DphiIn/F");
	outTree_ -> Branch("ele2_DetaIn",         &ele2_DetaIn,                  "ele2_DetaIn/F");
	outTree_ -> Branch("ele2_HOverE",         &ele2_HOverE,                  "ele2_HOverE/F");
	outTree_ -> Branch("ele2_tkIso",          &ele2_tkIso,                    "ele2_tkIso/F");
	outTree_ -> Branch("ele2_emIso",          &ele2_emIso,                    "ele2_emIso/F");
	outTree_ -> Branch("ele2_hadIso",         &ele2_hadIso,                  "ele2_hadIso/F");
	
	outTree_ -> Branch("ele2_dxy_PV",        &ele2_dxy_PV,   "ele2_dxy_PV/F");
	outTree_ -> Branch("ele2_dz_PV",       &ele2_dz_PV,  "ele2_dz_PV/F");
	outTree_ -> Branch("ele2_sigmaP",       &ele2_sigmaP,  "ele2_sigmaP/F");
	
	outTree_ -> Branch("ele2_scERaw",     &ele2_scERaw,         "ele2_scERaw/F");
	outTree_ -> Branch("ele2_scEtRaw",    &ele2_scEtRaw,       "ele2_scEtRaw/F");
	outTree_ -> Branch("ele2_scE",    &ele2_scE,         "ele2_scE/F");
	outTree_ -> Branch("ele2_scEt",       &ele2_scEt,             "ele2_scEt/F");
	outTree_ -> Branch("ele2_scE_regression",    &ele2_scE_regression,   "ele2_scE_regression/F");
	outTree_ -> Branch("ele2_scEerr_regression",       &ele2_scEerr_regression,         "ele2_scEerr_regression/F");
	outTree_ -> Branch("ele2_scE_regression_PhotonTuned",    &ele2_scE_regression_PhotonTuned ,   "ele2_scE_regression_PhotonTuned/F");
	outTree_ -> Branch("ele2_scEerr_regression_PhotonTuned",&ele2_scEerr_regression_PhotonTuned ,"ele2_scEerr_regression_PhotonTuned/F");
	outTree_ -> Branch("ele2_scERaw_PUcleaned",    &ele2_scERaw_PUcleaned,   "ele2_scERaw_PUcleaned/F");
	outTree_ -> Branch("ele2_es",       &ele2_es,  "ele2_es/F");
	outTree_ -> Branch("ele2_scLaserCorr",       &ele2_scLaserCorr,         "ele2_scLaserCorr/F");
	outTree_ -> Branch("ele2_scCrackCorr",   &ele2_scCrackCorr,        "ele2_scCrackCorr/F");
	outTree_ -> Branch("ele2_scLocalContCorr",    &ele2_scLocalContCorr,         "ele2_scLocalContCorr/F");
	
	outTree_ -> Branch("ele2_scEta",  &ele2_scEta,     "ele2_scEta/F");
	outTree_ -> Branch("ele2_scPhi",  &ele2_scPhi,     "ele2_scPhi/F");
	outTree_ -> Branch("ele2_scLocalEta",    &ele2_scLocalEta,         "ele2_scLocalEta/F");
	outTree_ -> Branch("ele2_scLocalPhi",    &ele2_scLocalPhi,         "ele2_scLocalPhi/F");
	outTree_ -> Branch("ele2_scEtaWidth",    &ele2_scEtaWidth,         "ele2_scEtaWidth/F");
	outTree_ -> Branch("ele2_scPhiWidth",    &ele2_scPhiWidth,         "ele2_scPhiWidth/F");
	outTree_ -> Branch("ele2_scEtaWidth_PUcleaned",       &ele2_scEtaWidth_PUcleaned,         "ele2_scEtaWidth_PUcleaned/F");
	outTree_ -> Branch("ele2_scPhiWidth_PUcleaned",    &ele2_scPhiWidth_PUcleaned,   "ele2_scPhiWidth_PUcleaned/F");
	outTree_ -> Branch("ele2_fCorrection_PUcleaned",       &ele2_fCorrection_PUcleaned,         "ele2_fCorrection_PUcleaned/F");
	
	outTree_ -> Branch("ele2_fEta",       &ele2_fEta,   "ele2_fEta/F");
	outTree_ -> Branch("ele2_fEtaCorr",   &ele2_fEtaCorr,         "ele2_fEtaCorr/F");
	outTree_ -> Branch("ele2_tkP",        &ele2_tkP,   "ele2_tkP/F");
	outTree_ -> Branch("ele2_tkPt",       &ele2_tkPt,  "ele2_tkPt/F");
	outTree_ -> Branch("ele2_fbrem",       &ele2_fbrem,  "ele2_fbrem/F");
	
	outTree_ -> Branch("ele2_eSeedBC",    &ele2_eSeedBC,  "ele2_eSeedBC/F");
	outTree_ -> Branch("ele2_e5x5",       &ele2_e5x5,  "ele2_e5x5/F");
	outTree_ -> Branch("ele2_e3x3",       &ele2_e3x3,  "ele2_e3x3/F");
	outTree_ -> Branch("ele2_scNxtal",       &ele2_scNxtal,  "ele2_scNxtal/F");
	outTree_ -> Branch("ele2_bcN",           &ele2_bcN,  "ele2_bcN/I");
	outTree_ -> Branch("ele2_5x5LaserCorr",       &ele2_5x5LaserCorr,  "ele2_5x5LaserCorr/F");
	outTree_ -> Branch("ele2_3x3LaserCorr",       &ele2_3x3LaserCorr,  "ele2_3x3LaserCorr/F");
	
	outTree_ -> Branch("ele2_seedE",         &ele2_seedE,                  "ele2_seedE/F");
	outTree_ -> Branch("ele2_seedLaserAlpha",&ele2_seedLaserAlpha,"ele2_seedLaserAlpha/F");
	outTree_ -> Branch("ele2_seedLaserCorr", &ele2_seedLaserCorr,  "ele2_seedLaserCorr/F");
	outTree_ -> Branch("ele2_seedICConstant",&ele2_seedICConstant,"ele2_seedICConstant/F");
	outTree_ -> Branch("ele2_seedIeta",      &ele2_seedIeta,            "ele2_seedIeta/I");
	outTree_ -> Branch("ele2_seedIphi",      &ele2_seedIphi,            "ele2_seedIphi/I");
	outTree_ -> Branch("ele2_seedIx",        &ele2_seedIx,                "ele2_seedIx/I");
	outTree_ -> Branch("ele2_seedIy",        &ele2_seedIy,                "ele2_seedIy/I");
	outTree_ -> Branch("ele2_seedZside",     &ele2_seedZside,          "ele2_seedZside/I");
	outTree_ -> Branch("ele2_EOverP",        &ele2_EOverP,                "ele2_EOverP/F");
	
	outTree_ -> Branch("ele2_nRecHits",&ele2_nRecHits,"ele2_nRecHits/I");
	outTree_ -> Branch("ele2_recHit_E",              "std::vector<float>",&ele2_recHit_E);
	outTree_ -> Branch("ele2_recHit_flag",           "std::vector<int>",  &ele2_recHit_flag);
	outTree_ -> Branch("ele2_recHit_hashedIndex",    "std::vector<int>",  &ele2_recHit_hashedIndex);
	outTree_ -> Branch("ele2_recHit_ietaORix",       "std::vector<int>",  &ele2_recHit_ietaORix);
	outTree_ -> Branch("ele2_recHit_iphiORiy",       "std::vector<int>",  &ele2_recHit_iphiORiy);
	outTree_ -> Branch("ele2_recHit_zside",          "std::vector<int>",  &ele2_recHit_zside);
	outTree_ -> Branch("ele2_recHit_laserCorrection","std::vector<float>",&ele2_recHit_laserCorrection);
	outTree_ -> Branch("ele2_recHit_Alpha",          "std::vector<float>",&ele2_recHit_Alpha);
	outTree_ -> Branch("ele2_recHit_ICConstant",     "std::vector<float>",&ele2_recHit_ICConstant);
	
	if(saveRecHitMatrix_)
	{
		outTree_ -> Branch("ele2_recHitMatrix_E",              "std::vector<float>",&ele2_recHitMatrix_E);
		outTree_ -> Branch("ele2_recHitMatrix_flag",           "std::vector<int>",  &ele2_recHitMatrix_flag);
		outTree_ -> Branch("ele2_recHitMatrix_hashedIndex",    "std::vector<int>",  &ele2_recHitMatrix_hashedIndex);
		outTree_ -> Branch("ele2_recHitMatrix_ietaORix",       "std::vector<int>",  &ele2_recHitMatrix_ietaORix);
		outTree_ -> Branch("ele2_recHitMatrix_iphiORiy",       "std::vector<int>",  &ele2_recHitMatrix_iphiORiy);
		outTree_ -> Branch("ele2_recHitMatrix_zside",          "std::vector<int>",  &ele2_recHitMatrix_zside);
		outTree_ -> Branch("ele2_recHitMatrix_laserCorrection","std::vector<float>",&ele2_recHitMatrix_laserCorrection);
		outTree_ -> Branch("ele2_recHitMatrix_ICConstant",     "std::vector<float>",&ele2_recHitMatrix_ICConstant);
		outTree_ -> Branch("ele2_recHitMatrix_samples",        "std::vector<float>",&ele2_recHitMatrix_samples);
	}
	
	
	// met variables
	std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::set met branches <<<" << std::endl;
	outTree_ -> Branch("met_et",      &met_et,            "met_et/F");
	outTree_ -> Branch("met_phi",     &met_phi,          "met_phi/F");
	outTree_ -> Branch("ele1Met_mt",  &ele1Met_mt,    "ele1Met_mt/F");
	outTree_ -> Branch("ele1Met_Dphi",&ele1Met_Dphi,"ele1Met_Dphi/F");
	
	
	// di-electron variables
	std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::set dielectron branches <<<" << std::endl;
	outTree_ -> Branch("ele1ele2_m",             &ele1ele2_m,                          "ele1ele2_m/F");
	outTree_ -> Branch("ele1ele2_scM",           &ele1ele2_scM,                      "ele1ele2_scM/F");
	outTree_ -> Branch("ele1ele2_scM_regression",&ele1ele2_scM_regression,"ele1ele2_scM_regression/F");
	
	
	// fbrem variables  
	if(saveFbrem_)
	{
		std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::set fBrem branches <<<" << std::endl;
		outTree_ -> Branch("ele1_inner_p",  &ele1_inner_p,     "ele1_inner_p/F");
		outTree_ -> Branch("ele1_inner_x",  &ele1_inner_x,     "ele1_inner_x/F");
		outTree_ -> Branch("ele1_inner_y",  &ele1_inner_y,     "ele1_inner_y/F");
		outTree_ -> Branch("ele1_inner_z",  &ele1_inner_z,     "ele1_inner_z/F");
		outTree_ -> Branch("ele1_outer_p",  &ele1_outer_p,     "ele1_outer_p/F");
		outTree_ -> Branch("ele1_outer_x",  &ele1_outer_x,     "ele1_outer_x/F");
		outTree_ -> Branch("ele1_outer_y",  &ele1_outer_y,     "ele1_outer_y/F");
		outTree_ -> Branch("ele1_outer_z",  &ele1_outer_z,     "ele1_outer_z/F");
		outTree_ -> Branch("ele1_tangent_n",&ele1_tangent_n, "ele1_tangent_n/I");
		outTree_ -> Branch("ele1_tangent_p",    "std::vector<float>",&ele1_tangent_p);
		outTree_ -> Branch("ele1_tangent_x",    "std::vector<float>",&ele1_tangent_x);
		outTree_ -> Branch("ele1_tangent_y",    "std::vector<float>",&ele1_tangent_y);
		outTree_ -> Branch("ele1_tangent_z",    "std::vector<float>",&ele1_tangent_z);
		outTree_ -> Branch("ele1_tangent_dP",   "std::vector<float>",&ele1_tangent_dP);
		outTree_ -> Branch("ele1_tangent_dPerr","std::vector<float>",&ele1_tangent_dPerr);
		
		outTree_ -> Branch("ele2_inner_p",  &ele2_inner_p,    "ele2_inner_p/F");
		outTree_ -> Branch("ele2_inner_x",  &ele2_inner_x,    "ele2_inner_x/F");
		outTree_ -> Branch("ele2_inner_y",  &ele2_inner_y,    "ele2_inner_y/F");
		outTree_ -> Branch("ele2_inner_z",  &ele2_inner_z,    "ele2_inner_z/F");
		outTree_ -> Branch("ele2_outer_p",  &ele2_outer_p,    "ele2_outer_p/F");
		outTree_ -> Branch("ele2_outer_x",  &ele2_outer_x,    "ele2_outer_x/F");
		outTree_ -> Branch("ele2_outer_y",  &ele2_outer_y,    "ele2_outer_y/F");
		outTree_ -> Branch("ele2_outer_z",  &ele2_outer_z,    "ele2_outer_z/F");
		outTree_ -> Branch("ele2_tangent_n",&ele2_tangent_n,"ele2_tangent_n/I");
		outTree_ -> Branch("ele2_tangent_p",    "std::vector<float>",&ele2_tangent_p);
		outTree_ -> Branch("ele2_tangent_x",    "std::vector<float>",&ele2_tangent_x);
		outTree_ -> Branch("ele2_tangent_y",    "std::vector<float>",&ele2_tangent_y);
		outTree_ -> Branch("ele2_tangent_z",    "std::vector<float>",&ele2_tangent_z);
		outTree_ -> Branch("ele2_tangent_dP",   "std::vector<float>",&ele2_tangent_dP);
		outTree_ -> Branch("ele2_tangent_dPerr","std::vector<float>",&ele2_tangent_dPerr);
	}
	
	
	// PU variables
	if(saveMCPU_)
	{
		std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::set MC PU branches <<<" << std::endl;
		outTree_ -> Branch("PUit_TrueNumInteractions",&PUit_TrueNumInteractions,"PUit_TrueNumInteractions/F");
		outTree_ -> Branch("PUit_NumInteractions",    &PUit_NumInteractions,        "PUit_NumInteractions/I");
		outTree_ -> Branch("PUit_zpositions",         &PUit_zpositions,                  "PUit_zpositions/F");
		outTree_ -> Branch("PUit_sumpT_lowpT",        &PUit_sumpT_lowpT,                "PUit_sumpT_lowpT/F");
		outTree_ -> Branch("PUit_sumpT_highpT",       &PUit_sumpT_highpT,              "PUit_sumpT_highpT/F");
		outTree_ -> Branch("PUit_ntrks_lowpT",        &PUit_ntrks_lowpT,                "PUit_ntrks_lowpT/F");
		outTree_ -> Branch("PUit_ntrks_highpT",       &PUit_ntrks_highpT,              "PUit_ntrks_highpT/F");
		
		outTree_ -> Branch("PUoot_early_TrueNumInteractions",&PUoot_early_TrueNumInteractions,"PUoot_early_TrueNumInteractions/F");
		outTree_ -> Branch("PUoot_early",                    &PUoot_early,                                        "PUoot_early/I");
		outTree_ -> Branch("PUoot_early_zpositions",         &PUoot_early_zpositions,                  "PUoot_early_zpositions/F");
		outTree_ -> Branch("PUoot_early_sumpT_lowpT",        &PUoot_early_sumpT_lowpT,                "PUoot_early_sumpT_lowpT/F");
		outTree_ -> Branch("PUoot_early_sumpT_highpT",       &PUoot_early_sumpT_highpT,              "PUoot_early_sumpT_highpT/F");
		outTree_ -> Branch("PUoot_early_ntrks_lowpT",        &PUoot_early_ntrks_lowpT,                "PUoot_early_ntrks_lowpT/F");
		outTree_ -> Branch("PUoot_early_ntrks_highpT",       &PUoot_early_ntrks_highpT,              "PUoot_early_ntrks_highpT/F");
		
		outTree_ -> Branch("PUoot_late_TrueNumInteractions",&PUoot_late_TrueNumInteractions,"PUoot_late_TrueNumInteractions/F");
		outTree_ -> Branch("PUoot_late",                    &PUoot_late,                                        "PUoot_late/I");
		outTree_ -> Branch("PUoot_late_zpositions",         &PUoot_late_zpositions,                  "PUoot_late_zpositions/F");
		outTree_ -> Branch("PUoot_late_sumpT_lowpT",        &PUoot_late_sumpT_lowpT,                "PUoot_late_sumpT_lowpT/F");
		outTree_ -> Branch("PUoot_late_sumpT_highpT",       &PUoot_late_sumpT_highpT,              "PUoot_late_sumpT_highpT/F");
		outTree_ -> Branch("PUoot_late_ntrks_lowpT",        &PUoot_late_ntrks_lowpT,                "PUoot_late_ntrks_lowpT/F");
		outTree_ -> Branch("PUoot_late_ntrks_highpT",       &PUoot_late_ntrks_highpT,              "PUoot_late_ntrks_highpT/F");
	}
	
	
	// MC info
	if(saveMCInfo_)
	{
		std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::set MC branches <<<" << std::endl;
		outTree_ -> Branch("mcV_E",     &mcV_E,           "mcV_E/F");
		outTree_ -> Branch("mcV_Px",    &mcV_Px,         "mcV_Px/F");
		outTree_ -> Branch("mcV_Py",    &mcV_Py,         "mcV_Py/F");
		outTree_ -> Branch("mcV_Pz",    &mcV_Pz,         "mcV_Pz/F");
		outTree_ -> Branch("mcV_Charge",&mcV_Charge, "mcV_Charge/I");
		outTree_ -> Branch("mcV_PdgId", &mcV_PdgId,   "mcV_PdgId/I");
		
		outTree_ -> Branch("mcF1_fromV_E",     &mcF1_fromV_E,          "mcF1_fromV_E/F");
		outTree_ -> Branch("mcF1_fromV_Px",    &mcF1_fromV_Px,        "mcF1_fromV_Px/F");
		outTree_ -> Branch("mcF1_fromV_Py",    &mcF1_fromV_Py,        "mcF1_fromV_Py/F");
		outTree_ -> Branch("mcF1_fromV_Pz",    &mcF1_fromV_Pz,        "mcF1_fromV_Pz/F");
		outTree_ -> Branch("mcF1_fromV_Charge",&mcF1_fromV_Charge,"mcF1_fromV_Charge/I");
		outTree_ -> Branch("mcF1_fromV_PdgId", &mcF1_fromV_PdgId,  "mcF1_fromV_PdgId/I");
		
		outTree_ -> Branch("mcF2_fromV_E",     &mcF2_fromV_E,          "mcF2_fromV_E/F");
		outTree_ -> Branch("mcF2_fromV_Px",    &mcF2_fromV_Px,        "mcF2_fromV_Px/F");
		outTree_ -> Branch("mcF2_fromV_Py",    &mcF2_fromV_Py,        "mcF2_fromV_Py/F");
		outTree_ -> Branch("mcF2_fromV_Pz",    &mcF2_fromV_Pz,        "mcF2_fromV_Pz/F");
		outTree_ -> Branch("mcF2_fromV_Charge",&mcF2_fromV_Charge,"mcF2_fromV_Charge/I");
		outTree_ -> Branch("mcF2_fromV_PdgId", &mcF2_fromV_PdgId,  "mcF2_fromV_PdgId/I");
	}
	
	
	std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::get cluster corrections <<<" << std::endl;
	EcalClusterCrackCorrection = EcalClusterFunctionFactory::get()->create("EcalClusterCrackCorrection", iConfig);
	EcalClusterLocalContCorrection = EcalClusterFunctionFactory::get()->create("EcalClusterLocalContCorrection", iConfig);
	
	
	// JSON file map 
	if( jsonFlag_ == true )
	{
		std::cout<< ">>>>>> SimpleNtupleEoverP::SimpleNtupleEoverP::read json file <<<" << std::endl;
		jsonMap_ = readJSONFile(jsonFileName_);
	}
	
	std::cout<< ">>> SimpleNtupleEoverP::SimpleNtupleEoverP end <<<" << std::endl;
}

// --------------------------------------------------------------------



SimpleNtupleEoverP::~SimpleNtupleEoverP ()
{
	std::cout<< ">>> SimpleNtupleEoverP::~SimpleNtupleEoverP <<< analyzed " <<  eventNaiveId_ << " events" << std::endl;
	// save tree
	// outTree_ -> Write();
}

// -----------------------------------------------------------------------------------------



void SimpleNtupleEoverP::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if( verbosity_ )
		std::cout<< ">>> SimpleNtupleEoverP::analyze begin <<<" << std::endl;
	
	
	++eventNaiveId_;
	
	bool isGoodEvent = false;
	
	// event variables 
	bxId = iEvent.bunchCrossing();
	runId = iEvent.id().run();
	lumiId = iEvent.luminosityBlock();
	eventId = iEvent.id().event();
	timeStampHigh = (int)(iEvent.time().value() >> 32);
	
	isW = -1;
	isZ = -1;
	
	PUit_TrueNumInteractions= -99.;
	PUit_NumInteractions= -99;
	PUit_zpositions = -99.;
	PUit_sumpT_lowpT = -99.;
	PUit_sumpT_highpT = -99.;
	PUit_ntrks_lowpT = -99.;
	PUit_ntrks_highpT = -99.;
	
	PUoot_early_TrueNumInteractions = -99.;
	PUoot_early = -99;
	PUoot_early_zpositions = -99.;
	PUoot_early_sumpT_lowpT = -99.;
	PUoot_early_sumpT_highpT = -99.;
	PUoot_early_ntrks_lowpT = -99.;
	PUoot_early_ntrks_highpT = -99.;
	
	PUoot_late_TrueNumInteractions = -99.;
	PUoot_late = -99;
	PUoot_late_zpositions = -99.;
	PUoot_late_sumpT_lowpT = -99.;
	PUoot_late_sumpT_highpT = -99.; 
	PUoot_late_ntrks_lowpT = -99.;
	PUoot_late_ntrks_highpT = -99.;
	
	if(saveMCInfo_){
		mcV_E = -9999.;
		mcV_Px = -9999.;
		mcV_Py = -9999.;
		mcV_Pz = -9999.;
		mcV_Charge = -99;
		mcV_PdgId = -99;
		
		mcF1_fromV_E = -9999.;
		mcF1_fromV_Px = -9999.;
		mcF1_fromV_Py = -9999.;
		mcF1_fromV_Pz = -9999.;
		mcF1_fromV_Charge = -99;
		mcF1_fromV_PdgId = -99;
		mcF2_fromV_E = -9999.;
		mcF2_fromV_Px = -9999.;
		mcF2_fromV_Py = -9999.;
		mcF2_fromV_Pz = -9999.;
		mcF2_fromV_Charge = -99;
		mcF2_fromV_PdgId = -99;
	}
	
	// electron variables  
	ele1_charge =-99.;
	ele1_p =-99.;
	ele1_pt =-99.;
	ele1_eta =-99.;
	ele1_phi=-99.;
	ele1_isTrackerDriven=-99;
	
	ele1_sigmaIetaIeta =-99.;
	ele1_DphiIn =-99.;
	ele1_DetaIn =-99.;
	ele1_HOverE =-99.;
	ele1_tkIso =-99.;
	ele1_emIso =-99.;
	ele1_hadIso =-99.;
	
	ele1_dxy_PV=-99.;
	ele1_dz_PV=-99.;
	ele1_sigmaP = -99.;
	ele1_effAreaForIso=-99.;
	
	ele1_scERaw =-99.;
	ele1_scEtRaw = -99.;
	ele1_scEt = -99.;
	ele1_scLocalEta =-99.;
	ele1_scLocalPhi = -99.;
	ele1_scEtaWidth = -99.;
	ele1_scPhiWidth = -99.;
	ele1_scCrackCorr =-99.;
	ele1_scLocalContCorr =-99.;
	ele1_scE =-99.;
	ele1_scEta =-99.;
	ele1_scPhi =-99.;
	ele1_scLaserCorr =-99.;
	ele1_scE_regression =-99.;
	ele1_scEerr_regression=-99.;
	ele1_scE_regression_PhotonTuned=-99.;
	ele1_scEerr_regression_PhotonTuned=-99.;
	ele1_scERaw_PUcleaned=-99.;
	ele1_scEtaWidth_PUcleaned=-99.;
	ele1_scPhiWidth_PUcleaned=-99.;
	ele1_fCorrection_PUcleaned=-99.;
	
	ele1_fEta=-99.;
	ele1_fEtaCorr=-99.;
	ele1_tkP=-99.;
	ele1_tkPt=-99.;
	ele1_fbrem=-99.;
	
	
	ele1_eSeedBC=-99.;
	ele1_e5x5=-99.;
	ele1_e3x3=-99.;
	ele1_scNxtal=-99.;
	ele1_bcN    =-99;
	ele1_5x5LaserCorr=-99.;
	ele1_3x3LaserCorr=-99.;
	ele1_es=-99.;
	
	ele1_seedE=-99.;
	ele1_seedLaserAlpha=-99.;
	ele1_seedLaserCorr=-99.;
	ele1_seedICConstant=-99.;
	ele1_seedIeta=-9999;
	ele1_seedIphi=-9999;
	ele1_seedIx=-9999;
	ele1_seedIy=-9999;
	ele1_seedZside=-9999;
	ele1_EOverP=-99.;
	
	ele1_recHit_E.clear();
	ele1_recHit_flag.clear();
	ele1_recHit_hashedIndex.clear();
	ele1_recHit_ietaORix.clear();
	ele1_recHit_iphiORiy.clear();
	ele1_recHit_zside.clear();
	ele1_recHit_laserCorrection.clear();
	ele1_recHit_Alpha.clear();
	ele1_recHit_ICConstant.clear();
	ele1_nRecHits = -9999;
	
	ele1_recHitMatrix_E.clear();
	ele1_recHitMatrix_flag.clear();
	ele1_recHitMatrix_hashedIndex.clear();
	ele1_recHitMatrix_ietaORix.clear();
	ele1_recHitMatrix_iphiORiy.clear();
	ele1_recHitMatrix_zside.clear();
	ele1_recHitMatrix_laserCorrection.clear();
	ele1_recHitMatrix_ICConstant.clear();
	ele1_recHitMatrix_samples.clear();
	
	ele1_isEB= -9999;
	ele1_isEBEEGap= -9999;
	ele1_isEBEtaGap= -9999;
	ele1_isEBPhiGap= -9999;
	ele1_isEEDeeGap= -9999;
	ele1_isEERingGap= -9999;
	
	ele1_eRegrInput_rawE = -99.;
	ele1_eRegrInput_r9 = -99.;
	ele1_eRegrInput_eta = -99.;
	ele1_eRegrInput_phi= -99.;
	ele1_eRegrInput_r25= -99.;
	ele1_eRegrInput_etaW= -99.;
	ele1_eRegrInput_phiW= -99.;
	ele1_eRegrInput_rho= -99.;
	ele1_eRegrInput_Deta_bC_sC= -99.;
	ele1_eRegrInput_Dphi_bC_sC= -99.;
	ele1_eRegrInput_bCE_Over_sCE= -99.;
	ele1_eRegrInput_e3x3_Over_bCE= -99.;
	ele1_eRegrInput_e5x5_Over_bCE= -99.;
	ele1_eRegrInput_sigietaieta_bC1= -99.;
	ele1_eRegrInput_sigiphiiphi_bC1= -99.;
	ele1_eRegrInput_sigietaiphi_bC1= -99.;
	ele1_eRegrInput_bEMax_Over_bCE= -99.;
	ele1_eRegrInput_bE2nd_Over_bCE= -99.;
	ele1_eRegrInput_bEtop_Over_bCE= -99.;
	ele1_eRegrInput_bEbot_Over_bCE= -99.;
	
	ele1_eRegrInput_bEleft_Over_bCE= -99.;
	ele1_eRegrInput_bEright_Over_bCE= -99.;
	ele1_eRegrInput_be2x5max_Over_bCE= -99.;
	ele1_eRegrInput_be2x5top_Over_bCE= -99.;
	ele1_eRegrInput_be2x5bottom_Over_bCE= -99.;
	ele1_eRegrInput_be2x5left_Over_bCE= -99.;
	ele1_eRegrInput_be2x5right_Over_bCE= -99.;
	
	ele1_eRegrInput_seedbC_eta= -99.;
	ele1_eRegrInput_seedbC_phi= -99.;
	ele1_eRegrInput_seedbC_eta_p5= -99.;
	ele1_eRegrInput_seedbC_phi_p2= -99.;
	ele1_eRegrInput_seedbC_bieta= -99.;
	ele1_eRegrInput_seedbC_phi_p20= -99.;
	ele1_eRegrInput_seedbC_etacry= -99.;
	ele1_eRegrInput_seedbC_phicry= -99.;
	
	ele1_eRegrInput_ESoSC= -99.;
	ele1_eRegrInput_nPV= -99.;
	ele1_eRegrInput_SCsize= -99.;
	
	
	ele2_charge =-99.;
	ele2_p =-99.;
	ele2_pt =-99.;
	ele2_eta =-99.;
	ele2_phi=-99.;
	ele2_isTrackerDriven=-99;
	
	ele2_sigmaIetaIeta =-99.;
	ele2_DphiIn =-99.;
	ele2_DetaIn =-99.;
	ele2_HOverE =-99.;
	ele2_tkIso =-99.;
	ele2_emIso =-99.;
	ele2_hadIso =-99.;
	
	ele2_dxy_PV=-99.;
	ele2_dz_PV=-99.;
	ele2_sigmaP = -99.;
	ele2_effAreaForIso=-99.;
	
	ele2_scERaw =-99.;
	ele2_scEtRaw = -99.;
	ele2_scEt = -99.;
	ele2_scLocalEta =-99.;
	ele2_scLocalPhi = -99.;
	ele2_scEtaWidth = -99.;
	ele2_scPhiWidth = -99.;
	ele2_scCrackCorr =-99.;
	ele2_scLocalContCorr =-99.;
	ele2_scE =-99.;
	ele2_scEta =-99.;
	ele2_scPhi =-99.;
	ele2_scLaserCorr =-99.;
	ele2_scE_regression =-99.;
	ele2_scEerr_regression=-99.;
	ele2_scE_regression_PhotonTuned=-99.;
	ele2_scEerr_regression_PhotonTuned=-99.;
	ele2_scERaw_PUcleaned=-99.;
	ele2_scEtaWidth_PUcleaned=-99.;
	ele2_scPhiWidth_PUcleaned=-99.;
	ele2_fCorrection_PUcleaned=-99.;
	
	ele2_fEta=-99.;
	ele2_fEtaCorr=-99.;
	ele2_tkP=-99.;
	ele2_tkPt=-99.;
	ele2_fbrem=-99.;
	
	ele2_eSeedBC=-99.;
	ele2_e5x5=-99.;
	ele2_e3x3=-99.;
	ele2_scNxtal=-99.;
	ele2_bcN    =-99;
	ele2_5x5LaserCorr=-99.;
	ele2_3x3LaserCorr=-99.;
	ele2_es=-99.;
	
	ele2_seedE=-99.;
	ele2_seedLaserAlpha=-99.;
	ele2_seedLaserCorr=-99.;
	ele2_seedICConstant=-99.;
	ele2_seedIeta=-9999;
	ele2_seedIphi=-9999;
	ele2_seedIx=-9999;
	ele2_seedIy=-9999;
	ele2_seedZside=-9999;
	ele2_EOverP=-99.;
	
	ele2_recHit_E.clear();
	ele2_recHit_flag.clear();
	ele2_recHit_hashedIndex.clear();
	ele2_recHit_ietaORix.clear();
	ele2_recHit_iphiORiy.clear();
	ele2_recHit_zside.clear();
	ele2_recHit_laserCorrection.clear();
	ele2_recHit_Alpha.clear();
	ele2_recHit_ICConstant.clear();
	ele2_nRecHits = -9999;
	
	ele2_recHitMatrix_E.clear();
	ele2_recHitMatrix_flag.clear();
	ele2_recHitMatrix_hashedIndex.clear();
	ele2_recHitMatrix_ietaORix.clear();
	ele2_recHitMatrix_iphiORiy.clear();
	ele2_recHitMatrix_zside.clear();
	ele2_recHitMatrix_laserCorrection.clear();
	ele2_recHitMatrix_ICConstant.clear();
	ele2_recHitMatrix_samples.clear();
	
	ele2_isEB= -9999;
	ele2_isEBEEGap= -9999;
	ele2_isEBEtaGap= -9999;
	ele2_isEBPhiGap= -9999;
	ele2_isEEDeeGap= -9999;
	ele2_isEERingGap= -9999;
	
	ele2_eRegrInput_rawE = -99.;
	ele2_eRegrInput_r9 = -99.;
	ele2_eRegrInput_eta = -99.;
	ele2_eRegrInput_phi= -99.;
	ele2_eRegrInput_r25= -99.;
	ele2_eRegrInput_etaW= -99.;
	ele2_eRegrInput_phiW= -99.;
	ele2_eRegrInput_rho= -99.;
	ele2_eRegrInput_Deta_bC_sC= -99.;
	ele2_eRegrInput_Dphi_bC_sC= -99.;
	ele2_eRegrInput_bCE_Over_sCE= -99.;
	ele2_eRegrInput_e3x3_Over_bCE= -99.;
	ele2_eRegrInput_e5x5_Over_bCE= -99.;
	ele2_eRegrInput_sigietaieta_bC1= -99.;
	ele2_eRegrInput_sigiphiiphi_bC1= -99.;
	ele2_eRegrInput_sigietaiphi_bC1= -99.;
	ele2_eRegrInput_bEMax_Over_bCE= -99.;
	ele2_eRegrInput_bE2nd_Over_bCE= -99.;
	ele2_eRegrInput_bEtop_Over_bCE= -99.;
	ele2_eRegrInput_bEbot_Over_bCE= -99.;
	
	ele2_eRegrInput_bEleft_Over_bCE= -99.;
	ele2_eRegrInput_bEright_Over_bCE= -99.;
	ele2_eRegrInput_be2x5max_Over_bCE= -99.;
	ele2_eRegrInput_be2x5top_Over_bCE= -99.;
	ele2_eRegrInput_be2x5bottom_Over_bCE= -99.;
	ele2_eRegrInput_be2x5left_Over_bCE= -99.;
	ele2_eRegrInput_be2x5right_Over_bCE= -99.;
	
	ele2_eRegrInput_seedbC_eta= -99.;
	ele2_eRegrInput_seedbC_phi= -99.;
	ele2_eRegrInput_seedbC_eta_p5= -99.;
	ele2_eRegrInput_seedbC_phi_p2= -99.;
	ele2_eRegrInput_seedbC_bieta= -99.;
	ele2_eRegrInput_seedbC_phi_p20= -99.;
	ele2_eRegrInput_seedbC_etacry= -99.;
	ele2_eRegrInput_seedbC_phicry= -99.;
	
	ele2_eRegrInput_ESoSC= -99.;
	ele2_eRegrInput_nPV= -99.;
	ele2_eRegrInput_SCsize= -99.;
	
	if(saveFbrem_){
		ele1_inner_p = -9999.;
		ele1_inner_x = -9999.;
		ele1_inner_y = -9999.;
		ele1_inner_z = -9999.;
		ele1_outer_p = -9999.;
		ele1_outer_x = -9999.;
		ele1_outer_y = -9999.;
		ele1_outer_z = -9999.;
		ele1_tangent_n = -1;
		ele1_tangent_p.clear();
		ele1_tangent_x.clear();
		ele1_tangent_y.clear();
		ele1_tangent_z.clear();
		ele1_tangent_dP.clear();
		ele1_tangent_dPerr.clear();
		
		ele2_inner_p = -9999.;
		ele2_inner_x = -9999.;
		ele2_inner_y = -9999.;
		ele2_inner_z = -9999.;
		ele2_outer_p = -9999.;
		ele2_outer_x = -9999.;
		ele2_outer_y = -9999.;
		ele2_outer_z = -9999.;
		ele2_tangent_n = -1;
		ele2_tangent_p.clear();
		ele2_tangent_x.clear();
		ele2_tangent_y.clear();
		ele2_tangent_z.clear();
		ele2_tangent_dP.clear();
		ele2_tangent_dPerr.clear();
	}
	met_et=-99.;
	met_phi=-99.;
	
	ele1Met_mt=-99.;
	ele1Met_Dphi=-99.;
	
	
	// di-electron variables
	ele1ele2_m=-99.;
	ele1ele2_scM=-99.;
	ele1ele2_scM_regression=-99.;
	
	// Accept event from json file
	if( jsonFlag_ == true )
	{
		bool skipEvent = false;
		if(AcceptEventByRunAndLumiSection(runId,lumiId,jsonMap_) == false) skipEvent = true;
		if( (jsonFlag_ == true) && (skipEvent == true) ) return;
	}
	
	//************* VERTEXES
	fillPVInfo (iEvent,iSetup) ;
	
	//************* RHO for ISO
	fillRhoInfo(iEvent,iSetup) ; 
	
	//************* MC PU
	if (saveMCPU_) fillMCPUInfo (iEvent, iSetup);
	
	//************ MC INFO
	if(saveMCInfo_) {
		edm::Handle<reco::GenParticleCollection> genParticles;
		iEvent.getByLabel(MCtruthTag_, genParticles);
		eventType_ = 0;
		mcAnalysisZW_ = new MCDumperZW(genParticles, eventType_, verbosity_);
		fillMCInfo(iEvent, iSetup);
	}
	
	// MY STUFF (GEN): 
	{
		Handle<reco::GenParticleCollection> genCollection;
		if(!(iEvent.getByLabel("genParticles" ,genCollection)))
		{
			std::cout << "gen particles not found" << std::endl;
			return;
		}
		iEvent.getByLabel("genParticles",genCollection);
		//int n_daughters=0;
		NGEN = 0;
		for(reco::GenParticleCollection::const_iterator it=genCollection->begin(); it!=genCollection->end() ; ++it)
		{
            if (NGEN<NGEN_max) {
                if ((abs(it->pdgId())==11 && it->status()==1 && it->p4().E()>10.0 && fabs(it->p4().Eta())<4.5)//ele
                    //|| abs(it->pdgId())==13//mu
                    //|| abs(it->pdgId())==15//tau
                    || (doPhotons && (abs(it->pdgId())==22 && it->p4().E()>10.0 && fabs(it->p4().Eta())<4.5 && it->p4().Pt()>10.0 && it->status()==1))//photon
                    //|| abs(it->pdgId())==23//Z0
                    //|| abs(it->pdgId())==24//24 = W+
                    //|| abs(it->pdgId())==25//H
                    ) {
                    
                    GEN_pt[NGEN] = it->p4().Pt();
                    GEN_eta[NGEN] = it->p4().Eta();
                    GEN_phi[NGEN] = it->p4().Phi();
                    GEN_px[NGEN] = it->p4().Px();
                    GEN_py[NGEN] = it->p4().Py();
                    GEN_pz[NGEN] = it->p4().Pz();
                    GEN_e[NGEN] = it->p4().E();
                    GEN_mass[NGEN] = it->mass();
                    GEN_q[NGEN] = it->charge();
                    GEN_id[NGEN] = it->pdgId();
                    GEN_status[NGEN] = it->status();
                    int parID = -99;
                    const Candidate * mom = it->mother();
                    parID = mom->pdgId();
                    if (parID == GEN_id[NGEN]) {
                        const Candidate * mom2 = mom->mother();
                        parID = mom2->pdgId();
                        if (parID == GEN_id[NGEN]) {
                            const Candidate * mom3 = mom2->mother();
                            parID = mom3->pdgId();
                            if (parID == GEN_id[NGEN]) {
                                const Candidate * mom4 = mom3->mother();
                                parID = mom4->pdgId();
                                if (parID == GEN_id[NGEN]) {
                                    const Candidate * mom5 = mom4->mother();
                                    parID = mom5->pdgId();
                                }
                            }
                        }
                    }
                    GEN_parent[NGEN] = parID;
                    NGEN++; 
                }
            }
        }
	}
	//--------------------------------------------------
	
    /* // disabled for smaller size, May 2014
	// from Michael
    edm::ESHandle<CaloGeometry> geomHndl;
    iSetup.get<CaloGeometryRecord>().get(geomHndl);
    CaloGeometry const* geometry = geomHndl.product();
    
    edm::ESHandle<EcalIntercalibConstantsMC> pIcal;
    iSetup.get<EcalIntercalibConstantsMCRcd>().get(pIcal);
    const EcalIntercalibConstantsMC* ical = pIcal.product();
    const EcalIntercalibConstantMap& icalMap = ical->getMap();
    EcalIntercalibConstantMap::const_iterator icalit;
    
	Handle<EBDigiCollection> ebDigis;
	//if(!(iEvent.getByLabel("selectDigi","selectedEcalEBDigiCollection",ebDigis)))
	if(!(iEvent.getByLabel("simEcalDigis","ebDigis",ebDigis)))
	{
		std::cout << "EB digi collection not found" << std::endl;
		return;
	}
	Handle<EEDigiCollection> eeDigis;
	//if(!(iEvent.getByLabel("selectDigi","selectedEcalEEDigiCollection",eeDigis)))
	if(!(iEvent.getByLabel("simEcalDigis","eeDigis",eeDigis)))	
	{
		std::cout << "EE digi collection not found" << std::endl;
        return;
	} 
	Handle<EERecHitCollection> eeRecHits;
    if(!(iEvent.getByLabel("ecalRecHit","EcalRecHitsEE",eeRecHits)))
    {
        std::cout << "EE RecHits collection not found" << std::endl;
        return;
    }
    Handle<EBRecHitCollection> ebRecHits;
    if(!(iEvent.getByLabel("ecalRecHit","EcalRecHitsEB",ebRecHits)))
        // if(!(iEvent.getByLabel("selectDigi","selectedEcalEEDigiCollection",eeDigis)))
    {
        std::cout << "EB RecHits collection not found" << std::endl;
        return;
    }
 */
	/*
	// from this code (different DIGI collection), from Michael - is limited (?), "selectedDigi"
	//=========== EB DIGI HITS
	edm::Handle<EBDigiCollection> ebDigis;
	if(1){
		iEvent.getByLabel (digiCollection_EB_, ebDigis) ;
		if (! (ebDigis.isValid ()) ) {
			std::cerr << "EcalValidation::analyze -->  ebDigis not found" << std::endl;
		}
	}
	//=========== EE DIGI HITS
	edm::Handle<EEDigiCollection> eeDigis;
	if(1){
		iEvent.getByLabel (digiCollection_EE_, eeDigis) ;
		if (! (eeDigis.isValid ()) ) {
			std::cerr << "EcalValidation::analyze -->  eeDigis not found" << std::endl;
		}
	}
	*/
	
    //************* PHOTONS
    Handle<reco::PhotonCollection> photonCollection;
    if(!(iEvent.getByLabel("photons" ,photonCollection)))
    {
        std::cout << "photons not found" << std::endl;
        return;
    }
    iEvent.getByLabel("photons",photonCollection);
    pho_size=0;
    /*
    for(int i=0;i<200;i++)
    {
        p_ecalenergy[i]=-100;
        p_eta[i]=-100;
        p_phi[i]=-100;
        p_r9[i]=-100;
        p_energy_1x5[i]=-100;
        p_energy_2x5[i]=-100;
        p_energy_3x3[i]=-100;
        p_energy_5x5[i]=-100;
        p_esenergy[i]=-100;
        p_rawenergy[i]=-100;
        p_pt[i]=-100;
        p_isGap[i]=0;
        p_hOverE[i]=-100;
        p_sigmaIetaIeta[i]=-100;
    }
     */
    int p_size = 0;
    for(reco::PhotonCollection::const_iterator it=photonCollection->begin(); it!=photonCollection->end() ; ++it)
    {
        if (doPhotons && it->pt()>10) {
            p_eta[pho_size] = it->eta();
            p_phi[pho_size] = it->phi();
            p_ecalenergy[pho_size] = it->superCluster()->energy();
            p_r9[pho_size] = it->r9();
            p_energy_1x5[pho_size] = it->e1x5();
            p_energy_2x5[pho_size] = it->e2x5();
            p_energy_3x3[pho_size] = it->e3x3();
            p_energy_5x5[pho_size] = it->e5x5();
            p_esenergy[pho_size] = it->superCluster()->preshowerEnergy();
            p_rawenergy[pho_size]= it->superCluster()->rawEnergy();
            p_pt[pho_size]= it->pt();
            p_isGap[pho_size] = it->isEBEEGap();
            p_hOverE[pho_size]= it->hadronicOverEm();
            p_sigmaIetaIeta[pho_size]= it->sigmaIetaIeta();
            //const reco::CaloClusterPtr phoSeed = it->superCluster()->seed();
            p_photonenergy[pho_size]= it->getCorrectedEnergy(it->getCandidateP4type());
            p_nCrystals[pho_size]= it->superCluster()->size();
            
            /* // disabled for smaller size, May 2014
            // 2013-10-16 update from Michael
            DetId id = it->superCluster()->seed()->seed();
            if(id.subdetId()==EcalBarrel)
            {
                id = EBDetId::offsetBy(id,-2,-2); //set id to lowest ieta, iphi in 5x5
            }
            else
            {
                id = EEDetId::offsetBy(id,-2,-2); //set id to lowest ix, iy in 5x5
            }
            for(int i=0 ; i<25 ; i++)
            {
                //shift id through the 5x5 around the central xtal id
                if(id.subdetId()==EcalBarrel)
                {
                    EBDetId nId;
                    nId =(EBDetId) EBDetId::offsetBy(id,i%5,i/5); //increments in eta every 1 step, in phi every 5 steps from lower left corner of 5x5
                    if(nId)
                    {
                        for(EBDigiCollection::const_iterator ebdigiItr = ebDigis->begin();ebdigiItr!=ebDigis->end();++ebdigiItr)
                        {
                            if(ebdigiItr->id()!=nId)
                            continue;
                            if(ebdigiItr!=ebDigis->end()) //if find() fails, it returns last digi
                            {
                                EBDataFrame dataframe = (*ebdigiItr);
                                for(int j=0 ; j<10 ; j++)
                                {
                                    p_adc[p_size][i][j] = dataframe.sample(j).adc();
                                    p_gainId[p_size][i][j] = dataframe.sample(j).gainId();
                                }
                            }
                        }
                        for(EBRecHitCollection::const_iterator ebhitItr = ebRecHits->begin();ebhitItr!=ebRecHits->end();ebhitItr++)
                        {
                            if(ebhitItr->id()!=nId)
                            continue;
                            if(ebhitItr!=ebRecHits->end()) //if find() fails, it returns last digi
                            {
                                icalit = icalMap.find(nId);
                                EcalIntercalibConstant icalconst = 1;
                                if( icalit!=icalMap.end() )
                                {
                                    icalconst = (*icalit);
                                }
                                else
                                {
                                    edm::LogError("EcalRecHitError") << "No intercalib const found for xtal " << nId.rawId() << "! something wrong with EcalIntercalibConstants in your DB? ";
                                }
                                p_recieta[p_size][i] = nId.ieta();
                                p_reciphi[p_size][i] = nId.iphi();
                                
                                GlobalPoint const& position = geometry->getPosition(nId);
                                p_receta[p_size][i] = position.eta();
                                p_recphi[p_size][i] = position.phi();
                                
                                
                                p_intercalib[p_size][i] = icalconst;
                                p_recenergy[p_size][i] = ebhitItr->energy();
                                p_rectime[p_size][i] = ebhitItr->time();
                                p_recflag[p_size][i] = ebhitItr->flags();
                                for(int j =0;j<19;j++) //see ecalrechit.cc for enumeration of flags
                                p_recflags[p_size][i][j] = ebhitItr->checkFlag(j);
                            }
                        }
                    }
                }
                else
                {
                    EEDetId nId;
                    nId = (EEDetId) EEDetId::offsetBy(id,i%5,i/5); //increments in eta every 1 step, in phi every 5 steps from lower left corner of 5x5
                    if(nId)
                    {
                        //loops through entire digicollection
                        for(EEDigiCollection::const_iterator eedigiItr = eeDigis->begin();eedigiItr!=eeDigis->end();++eedigiItr)
                        {
                            if(eedigiItr->id()!=nId) //if not the crystal we are looking for, continue
                            continue;
                            if(eedigiItr!=eeDigis->end())
                            {
                                EEDataFrame dataframe = (*eedigiItr);
                                for(int j=0 ; j<10 ; j++)
                                {
                                    p_adc[p_size][i][j] = dataframe.sample(j).adc();
                                    p_gainId[p_size][i][j] = dataframe.sample(j).gainId();
                                }
                            }
                        }
                        for(EERecHitCollection::const_iterator eehitItr = eeRecHits->begin();eehitItr!=eeRecHits->end();++eehitItr)
                        {
                            if(eehitItr->id()!=nId)
                            continue;
                            if(eehitItr!=eeRecHits->end()) //if find() fails, it returns last digi
                            {
                                icalit = icalMap.find(nId);
                                EcalIntercalibConstant icalconst = 1;
                                if( icalit!=icalMap.end() )
                                {
                                    icalconst = (*icalit);
                                }
                                else
                                {
                                    edm::LogError("EcalRecHitError") << "No intercalib const found for xtal " << nId.rawId() << "! something wrong with EcalIntercalibConstants in your DB? ";
                                }
                                p_recieta[p_size][i] = nId.ix();
                                p_reciphi[p_size][i] = nId.iy();
                                
                                GlobalPoint const& position = geometry->getPosition(nId);
                                p_receta[p_size][i] = position.eta();
                                p_recphi[p_size][i] = position.phi();
                                
                                p_intercalib[p_size][i] = icalconst;
                                p_recenergy[p_size][i] = eehitItr->energy();
                                p_rectime[p_size][i] = eehitItr->time();
                                p_recflag[p_size][i] = eehitItr->flags();
                                for(int j =0;j<19;j++) //see ecalrechit.cc for enumeration of flags
                                p_recflags[p_size][i][j] = eehitItr->checkFlag(j);
                            }
                        }
                    }
                }
            }
             */
            p_size++;
            
            pho_size++;
        }
    }
    
	
	//************* ELECTRONS
	Handle<View<reco::GsfElectron> > electronHandle;
	iEvent.getByLabel(EleTag_,electronHandle);
	View<reco::GsfElectron> electrons = *electronHandle;
	
	///---- get the number of the electron in the event to know if it's a W or a Z ----
	//if(doWZSelection_){
	if (1) {
		
		int nEleTight=0,nEleMedium=0,nEleLoose=0;
		if (nEleLoose==nEleTight);
		eleIts_.clear();
		
		
		// MY STUFF
		NELE = 0;
		NELE = electrons.size();
		int ELE_id_ = -99;
		//---------------------				
		
		int p_size = 0;
		for ( unsigned int iEle=0; iEle<electrons.size(); ++iEle ){
			// MY STUFF
			reco::GsfElectron electron = electrons.at(iEle);
            if ((electron.p4()).P()>5) {
                ELE_id_ = -99;
                if (LooseEle(iEvent,iSetup,iEle)) ELE_id_=2;
                if (MediumEle(iEvent,iSetup,iEle)) ELE_id_=1;
                if (TightEle(iEvent,iSetup,iEle)) ELE_id_=0;
                
                ELE_pt[p_size] = electron.pt();
                ELE_eta[p_size] = electron.eta();
                ELE_phi[p_size] = electron.phi();
                ELE_px[p_size] = electron.px();
                ELE_py[p_size] = electron.py();
                ELE_pz[p_size] = electron.pz();
                ELE_e[p_size] = (electron.p4()).P();//electron.P();
                ELE_q[p_size] = electron.charge();
                ELE_id[p_size] = ELE_id_;//electron.id();
                
                ELE_sigmaIetaIeta[p_size] = electron.sigmaIetaIeta();
                ELE_DphiIn[p_size] = electron.deltaPhiSuperClusterTrackAtVtx();
                ELE_DetaIn[p_size] = electron.deltaEtaSuperClusterTrackAtVtx();
                ELE_HOverE[p_size] = electron.hadronicOverEm();
                ELE_ooemoop[p_size] = (1.0/electron.ecalEnergy() - electron.eSuperClusterOverP()/electron.ecalEnergy());
                ELE_tkIso[p_size] = electron.dr03TkSumPt();
                ELE_emIso[p_size] = electron.dr03EcalRecHitSumEt();
                ELE_hadIso[p_size] = electron.dr03HcalDepth1TowerSumEt()+electron.dr03HcalDepth2TowerSumEt();
                ELE_effAreaForIso[p_size] = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, electron.superCluster()->eta(), ElectronEffectiveArea::kEleEAData2011);
                ELE_combIso[p_size] = ELE_tkIso[p_size] +std::max(ELE_emIso[p_size] + ELE_hadIso[p_size] - rho*ELE_effAreaForIso[p_size],float(0.));
                reco::GsfTrackRef eleTrack  = electron.gsfTrack();
                ELE_dxy[p_size] = eleTrack->dxy (PVPoint_);
                ELE_dz[p_size] = eleTrack->dz (PVPoint_);
                Handle< double > rhoHandle;
                iEvent.getByLabel(rhoTag_,rhoHandle);
                //std::pair<double,double> cor = ecorr_.CorrectedEnergyWithErrorV3(electron,*hVertexProduct,*rhoHandle,lazyTools,iSetup);
                //ELE_scE_regression[p_size] = cor.first;
                reco::SuperClusterRef scRef = electron.superCluster();
                ELE_scEtaWidth[p_size] = scRef->etaWidth();
                ELE_scPhiWidth[p_size] = scRef->phiWidth();
                ELE_scERaw[p_size] = scRef->rawEnergy();
                ELE_scE[p_size] = scRef->energy();
                ELE_e5x5[p_size] = electron.e5x5();
                ELE_e1x5[p_size] = electron.e1x5();
                ELE_e2x5Max[p_size] = electron.e2x5Max();
                ELE_es[p_size] = scRef->preshowerEnergy();
                ELE_fbrem[p_size] = electron.fbrem();
                //ELE_[p_size] =
                // ---------------
                
                /* // disabled for smaller size, May 2014 
                // 2013-10-16 update from Michael
                DetId id = electron.superCluster()->seed()->seed();
                if(id.subdetId()==EcalBarrel)
                {
                    //EBDetId p_id = it->superCluster()->seed()->seed();
                    id = EBDetId::offsetBy(id,-2,-2); //set id to lowest ieta, iphi in 5x5
                }
                else
				//if(id.subdetId()==EcalEndcap)
                {
                    id = EEDetId::offsetBy(id,-2,-2); //set id to lowest ix, iy in 5x5
                }
                for(int i=0 ; i<25 ; i++)
                {
                    //shift id through the 5x5 around the central xtal id
                    if(id.subdetId()==EcalBarrel)
                    {
                        EBDetId nId;
                        nId =(EBDetId) EBDetId::offsetBy(id,i%5,i/5); //increments in eta every 1 step, in phi every 5 steps from lower left corner of 5x5
                        //std::cout << "eta " << i%5 << " phi " << i/5 << std::endl;
                        if(nId)
                        {
                            for(EBDigiCollection::const_iterator ebdigiItr = ebDigis->begin();ebdigiItr!=ebDigis->end();++ebdigiItr)
                            {
                                if(ebdigiItr->id()!=nId)
								continue;
                                //EBDigiCollection::const_iterator ebdigiItr = ebDigis->find(nId.rawId()); //failing here for EB-
                                if(ebdigiItr!=ebDigis->end()) //if find() fails, it returns last digi
                                {
                                    EBDataFrame dataframe = (*ebdigiItr);
                                    for(int j=0 ; j<10 ; j++)
                                    {
                                        ELE_adc[p_size][i][j] = dataframe.sample(j).adc();
                                        ELE_gainId[p_size][i][j] = dataframe.sample(j).gainId();
                                    }
                                }
                                //   std::cout <<"EB energy " << it->e5x5() << " gain " <<  dataframe.sample(5).gainId()<< " adc " << dataframe.sample(5).adc() << std::endl;
                            }
                            for(EBRecHitCollection::const_iterator ebhitItr = ebRecHits->begin();ebhitItr!=ebRecHits->end();ebhitItr++)
                            {
                                if(ebhitItr->id()!=nId)
                                continue;
                                //EBDigiCollection::const_iterator ebdigiItr = ebDigis->find(nId.rawId()); //failing here for EB-
                                if(ebhitItr!=ebRecHits->end()) //if find() fails, it returns last digi
                                {
                                    
                                    //std::cout << "about to get intercalib map" << std::endl;
                                    icalit = icalMap.find(nId);
                                    EcalIntercalibConstant icalconst = 1;
                                    if( icalit!=icalMap.end() )
                                    {
                                        icalconst = (*icalit);
                                    }
                                    else
                                    {
                                        edm::LogError("EcalRecHitError") << "No intercalib const found for xtal " << nId.rawId() << "! something wrong with EcalIntercalibConstants in your DB? ";
                                    }
                 
                                    //std::cout << "got intercalib map" << std::endl;
                                    ELE_recieta[p_size][i] = nId.ieta();
                                    ELE_reciphi[p_size][i] = nId.iphi();
                                    
                                    GlobalPoint const& position = geometry->getPosition(nId);
                                    ELE_receta[p_size][i] = position.eta();
                                    ELE_recphi[p_size][i] = position.phi();
                                    
                                    
                                    ELE_intercalib[p_size][i] = icalconst;
                                    ELE_recenergy[p_size][i] = ebhitItr->energy();
                                    ELE_rectime[p_size][i] = ebhitItr->time();
                                    ELE_recflag[p_size][i] = ebhitItr->flags();
                                    for(int j =0;j<19;j++) //see ecalrechit.cc for enumeration of flags
                                    ELE_recflags[p_size][i][j] = ebhitItr->checkFlag(j);
                                }
                            }
                        }
                    }
                    else
                    {
                        EEDetId nId;
                        nId = (EEDetId) EEDetId::offsetBy(id,i%5,i/5); //increments in eta every 1 step, in phi every 5 steps from lower left corner of 5x5
                        if(nId)
                        {
                            //loops through entire digicollection
                            for(EEDigiCollection::const_iterator eedigiItr = eeDigis->begin();eedigiItr!=eeDigis->end();++eedigiItr)
                            {
                                if(eedigiItr->id()!=nId) //if not the crystal we are looking for, continue
								continue;
                                //EEDigiCollection::const_iterator eedigiItr = eeDigis->find(nId.rawId());
                                if(eedigiItr!=eeDigis->end())
                                {
                                    EEDataFrame dataframe = (*eedigiItr);
                                    for(int j=0 ; j<10 ; j++)
                                    {
                                        ELE_adc[p_size][i][j] = dataframe.sample(j).adc();
                                        ELE_gainId[p_size][i][j] = dataframe.sample(j).gainId();
                                    }
                                }
                            }
                            for(EERecHitCollection::const_iterator eehitItr = eeRecHits->begin();eehitItr!=eeRecHits->end();++eehitItr)
                            {
                                if(eehitItr->id()!=nId)
                                continue;
                                //EBDigiCollection::const_iterator ebdigiItr = ebDigis->find(nId.rawId()); //failing here for EB-
                                if(eehitItr!=eeRecHits->end()) //if find() fails, it returns last digi
                                {
                                    icalit = icalMap.find(nId);
                                    EcalIntercalibConstant icalconst = 1;
                                    if( icalit!=icalMap.end() )
                                    {
                                        icalconst = (*icalit);
                                    }
                                    else
                                    {
                                        edm::LogError("EcalRecHitError") << "No intercalib const found for xtal " << nId.rawId() << "! something wrong with EcalIntercalibConstants in your DB? ";
                                    }
                 
                                    //std::cout << "got intercalib map" << std::endl;
                                    ELE_recieta[p_size][i] = nId.ix();
                                    ELE_reciphi[p_size][i] = nId.iy();
                                    
                                    GlobalPoint const& position = geometry->getPosition(nId);
                                    ELE_receta[p_size][i] = position.eta();
                                    ELE_recphi[p_size][i] = position.phi();
                                    
                                    ELE_intercalib[p_size][i] = icalconst;
                                    ELE_recenergy[p_size][i] = eehitItr->energy();
                                    ELE_rectime[p_size][i] = eehitItr->time();
                                    ELE_recflag[p_size][i] = eehitItr->flags();
                                    for(int j =0;j<19;j++) //see ecalrechit.cc for enumeration of flags
                                    ELE_recflags[p_size][i][j] = eehitItr->checkFlag(j);
                                }
                            }
                        }
                    }
                }
                 */
                p_size++;
            }

			
		}
		/*
		if( nEleTight < 1 ) return ;
		if( nEleTight > 2 ) return ;
		if( nEleMedium > 1) return ;
		if( nEleLoose > 0 ) return ;
		*/
		///---- check if the event is good----
		
		if( (nEleTight == 1) && (nEleMedium == 0) ){
			isW=1;isZ=0;
			std::map<float,int>::const_iterator mapIt = eleIts_.begin();
			fillEleInfo ( iEvent, iSetup, mapIt->second, "ele1" ); 
			fillMetInfo (iEvent, iSetup);
			
		}
		
		if( (nEleTight == 2) || (nEleTight == 1 && nEleMedium == 1) ){
			isW=0;isZ=1;
			std::map<float,int>::const_iterator mapIt = eleIts_.begin();
			fillEleInfo ( iEvent, iSetup, mapIt->second, "ele1" );
			mapIt++;
			fillEleInfo ( iEvent, iSetup, mapIt->second, "ele2" ); 
			fillDoubleEleInfo (iEvent, iSetup);
			fillMetInfo (iEvent, iSetup);
			
		}
		
		
		if( (nEleTight == 1) && (nEleMedium == 0) ) isGoodEvent = myWselection ( iEvent, iSetup); 
		if( (nEleTight == 2) || (nEleTight == 1 && nEleMedium == 1) ) isGoodEvent = myZselection ( iEvent, iSetup); 
		///---- save the entry of the tree only if W/Z event ----
		//if ( isGoodEvent )   outTree_ -> Fill();
		if ( isGoodEvent );
		//if (NELE>0 || ) outTree_ -> Fill();
		
	}
    if (NELE>0 || pho_size>0) outTree_ -> Fill();
    
    
	/*
	else{
		int nEle = electrons.size();
		if ( nEle == 1 ) { isW = 1; isZ = 0; }
		if ( nEle == 2 ) { isW = 0; isZ = 1; }
		
		if ( isW == 1 ) fillEleInfo ( iEvent, iSetup, 0, "ele1" ); 
		
		
		if ( isZ == 1 ) { 
			fillEleInfo ( iEvent, iSetup, 0, "ele1" ); 
			fillEleInfo ( iEvent, iSetup, 1, "ele2" );
			fillDoubleEleInfo (iEvent, iSetup);
		}
		
		fillMetInfo (iEvent, iSetup);
		
		if (isW==1 || isZ==1) outTree_ -> Fill();
	}*/
	
	
	if( verbosity_ )
		std::cout<< ">>> SimpleNtupleEoverP::analyze end <<<" << std::endl;
}



//-----------------------------------------------------------------------------------------
bool  SimpleNtupleEoverP::TightEle (const edm::Event & iEvent, const edm::EventSetup & iESetup, const int &iEle){
	
	//************* ELECTRONS
	Handle<View<reco::GsfElectron> > electronHandle;
	iEvent.getByLabel(EleTag_,electronHandle);
	View<reco::GsfElectron> electrons = *electronHandle;
	
	//********* CONVERSION TOOLS
	edm::Handle<reco::ConversionCollection> conversions_h;
	iEvent.getByLabel(conversionsInputTag_, conversions_h);
	
	bool isTightEle = false; 
	
	reco::GsfElectron electron = electrons.at(iEle);
	
	edm::InputTag EleBad = edm::InputTag("gsfElectrons");
	
	if(EleTag_== EleBad && !electron.ecalDriven()) return false; 
	
	float pt = electron.pt();
	float eta = electron.eta();
	
	float tkIso  = electron.dr03TkSumPt();
	float emIso  = electron.dr03EcalRecHitSumEt(); 
	float hadIso = electron.dr03HcalDepth1TowerSumEt()+electron.dr03HcalDepth2TowerSumEt();
	float AEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, electron.superCluster()->eta(), ElectronEffectiveArea::kEleEAData2011);
	
	// default
	float combIso = tkIso +std::max(emIso + hadIso - rho*AEff,float(0.));
	
	// 2011-like test
	//float combIso = tkIso + emIso + hadIso - rho*3.14159*0.3*0.3;
	
	int isEB = electron.isEB();
	float sigmaIetaIeta = electron.sigmaIetaIeta();
	float DetaIn        = electron.deltaEtaSuperClusterTrackAtVtx();
	float DphiIn        = electron.deltaPhiSuperClusterTrackAtVtx();
	float HOverE        = electron.hadronicOverEm();
	float ooemoop       = (1.0/electron.ecalEnergy() - electron.eSuperClusterOverP()/electron.ecalEnergy());
	
	int mishits             = electron.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
	int nAmbiguousGsfTracks = electron.ambiguousGsfTracksSize();
	
	
	reco::GsfTrackRef eleTrack  = electron.gsfTrack() ;
	float dxy           = eleTrack->dxy(PVPoint_);  
	float dz            = eleTrack->dz (PVPoint_);
	
	edm::Handle<reco::BeamSpot> BSHandle;
	iEvent.getByLabel(BSTag_, BSHandle);
	const reco::BeamSpot BS = *BSHandle;
	
	
	bool isConverted = ConversionTools::hasMatchedConversion(electron, conversions_h, BS.position());
	
	// default
	if(  (pt > 20.) && (fabs(eta) < 2.5) &&
	   ( ( (isEB == 1) && (fabs(DetaIn)  < 0.004) ) || ( (isEB == 0) && (fabs(DetaIn)  < 0.007) ) ) &&
	   ( ( (isEB == 1) && (fabs(DphiIn)  < 0.060) ) || ( (isEB == 0) && (fabs(DphiIn)  < 0.030) ) ) &&
	   ( ( (isEB == 1) && (sigmaIetaIeta < 0.010) ) || ( (isEB == 0) && (sigmaIetaIeta < 0.030) ) ) &&
	   ( ( (isEB == 1) && (HOverE        < 0.120) ) || ( (isEB == 0) && (HOverE        < 0.100) ) ) &&
	   ( ( (isEB == 1) && (fabs(ooemoop) < 0.050) ) || ( (isEB == 0) && (fabs(ooemoop) < 0.050) ) ) &&
	   ( ( (isEB == 1) && (fabs(dxy)     < 0.020) ) || ( (isEB == 0) && (fabs(dxy)     < 0.020) ) ) &&
	   ( ( (isEB == 1) && (fabs(dz)      < 0.100) ) || ( (isEB == 0) && (fabs(dz)      < 0.100) ) ) &&
	   ( ( (isEB == 1) && (!isConverted) ) || ( (isEB == 0) && (!isConverted) ) ) &&
	   ( ( (isEB == 1) && (combIso/pt    < 0.070) ) || ( (isEB == 0) && (combIso/pt    < 0.060) ) ) &&
	   ( mishits == 0 ) &&
	   ( nAmbiguousGsfTracks == 0 ) )
		isTightEle=true;
	
	// 2011-like test
	//if(  (pt > 20.) && (fabs(eta) < 2.5) &&
	//     ( ( (isEB == 1) && (fabs(DetaIn)  < 0.004) ) || ( (isEB == 0) && (fabs(DetaIn)  < 0.007) ) ) &&
	//     ( ( (isEB == 1) && (fabs(DphiIn)  < 0.060) ) || ( (isEB == 0) && (fabs(DphiIn)  < 0.030) ) ) &&
	//     ( ( (isEB == 1) && (sigmaIetaIeta < 0.010) ) || ( (isEB == 0) && (sigmaIetaIeta < 0.030) ) ) &&
	//     //( ( (isEB == 1) && (HOverE        < 0.040) ) || ( (isEB == 0) && (HOverE        < 0.025) ) ) &&
	//     ( ( (isEB == 1) && (combIso/pt    < 0.070) ) || ( (isEB == 0) && (combIso/pt    < 0.060) ) ) &&
	//     ( mishits == 0 ) &&
	//     ( nAmbiguousGsfTracks == 0 ) )
	//  isTightEle=true;
	
	return isTightEle;
	
}

bool  SimpleNtupleEoverP::MediumEle (const edm::Event & iEvent, const edm::EventSetup & iESetup, const int &iEle){
	
	//************* ELECTRONS
	Handle<View<reco::GsfElectron> > electronHandle;
	iEvent.getByLabel(EleTag_,electronHandle);
	View<reco::GsfElectron> electrons = *electronHandle;
	
	//********* CONVERSION TOOLS
	edm::Handle<reco::ConversionCollection> conversions_h;
	iEvent.getByLabel(conversionsInputTag_, conversions_h);
	
	bool isMediumEle = false; 
	
	reco::GsfElectron electron = electrons.at(iEle);
	
	edm::InputTag EleBad = edm::InputTag("gsfElectrons");
	
	if(EleTag_== EleBad && !electron.ecalDriven()) return false; 
	
	float pt = electron.pt();
	float eta = electron.eta();
	
	float tkIso  = electron.dr03TkSumPt();
	float emIso  = electron.dr03EcalRecHitSumEt(); 
	float hadIso = electron.dr03HcalDepth1TowerSumEt()+electron.dr03HcalDepth2TowerSumEt();
	
	float AEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, electron.superCluster()->eta(), ElectronEffectiveArea::kEleEAData2011);
	
	// default
	float combIso = tkIso +std::max(emIso + hadIso - rho*AEff,float(0.));
	
	// 2011-like test
	//float combIso = tkIso + emIso + hadIso - rho*3.14159*0.3*0.3;
	
	int isEB = electron.isEB();
	float sigmaIetaIeta = electron.sigmaIetaIeta();
	float DetaIn        = electron.deltaEtaSuperClusterTrackAtVtx();
	float DphiIn        = electron.deltaPhiSuperClusterTrackAtVtx();
	float HOverE        = electron.hadronicOverEm();
	float ooemoop       = (1.0/electron.ecalEnergy() - electron.eSuperClusterOverP()/electron.ecalEnergy());
	
	int mishits             = electron.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
	int nAmbiguousGsfTracks = electron.ambiguousGsfTracksSize();
	
	
	reco::GsfTrackRef eleTrack  = electron.gsfTrack() ;
	float dxy           = eleTrack->dxy(PVPoint_);  
	float dz            = eleTrack->dz (PVPoint_);
	
	edm::Handle<reco::BeamSpot> BSHandle;
	iEvent.getByLabel(BSTag_, BSHandle);
	const reco::BeamSpot BS = *BSHandle;
	
	
	bool isConverted = ConversionTools::hasMatchedConversion(electron, conversions_h, BS.position());
	
	// default
	if(  (pt > 12.) && (fabs(eta) < 2.5) &&
	   ( ( (isEB == 1) && (fabs(DetaIn)  < 0.004) ) || ( (isEB == 0) && (fabs(DetaIn)  < 0.007) ) ) &&
	   ( ( (isEB == 1) && (fabs(DphiIn)  < 0.060) ) || ( (isEB == 0) && (fabs(DphiIn)  < 0.030) ) ) &&
	   ( ( (isEB == 1) && (sigmaIetaIeta < 0.010) ) || ( (isEB == 0) && (sigmaIetaIeta < 0.030) ) ) &&
	   ( ( (isEB == 1) && (HOverE        < 0.120) ) || ( (isEB == 0) && (HOverE        < 0.100) ) ) &&
	   ( ( (isEB == 1) && (fabs(ooemoop) < 0.050) ) || ( (isEB == 0) && (fabs(ooemoop) < 0.050) ) ) &&
	   ( ( (isEB == 1) && (fabs(dxy)     < 0.020) ) || ( (isEB == 0) && (fabs(dxy)     < 0.020) ) ) &&
	   ( ( (isEB == 1) && (fabs(dz)      < 0.100) ) || ( (isEB == 0) && (fabs(dz)      < 0.100) ) ) &&
	   ( ( (isEB == 1) && (!isConverted) ) || ( (isEB == 0) && (!isConverted) ) ) &&
	   ( ( (isEB == 1) && (combIso/pt    < 0.070) ) || ( (isEB == 0) && (combIso/pt    < 0.060) ) ) &&
	   ( mishits == 0 ) &&
	   ( nAmbiguousGsfTracks == 0 )      
	   )
		isMediumEle=true;
	
	// 2011-like test
	//if(  (pt > 12.) && (fabs(eta) < 2.5) &&
	//     ( ( (isEB == 1) && (fabs(DetaIn)  < 0.004) ) || ( (isEB == 0) && (fabs(DetaIn)  < 0.007) ) ) &&
	//     ( ( (isEB == 1) && (fabs(DphiIn)  < 0.060) ) || ( (isEB == 0) && (fabs(DphiIn)  < 0.030) ) ) &&
	//     ( ( (isEB == 1) && (sigmaIetaIeta < 0.010) ) || ( (isEB == 0) && (sigmaIetaIeta < 0.030) ) ) &&
	//     //( ( (isEB == 1) && (HOverE        < 0.040) ) || ( (isEB == 0) && (HOverE        < 0.025) ) ) &&
	//     ( ( (isEB == 1) && (combIso/pt    < 0.070) ) || ( (isEB == 0) && (combIso/pt    < 0.060) ) ) &&
	//     ( mishits == 0 ) &&
	//     ( nAmbiguousGsfTracks == 0 )      
	//       )
	//   isMediumEle=true;
	
	
	return isMediumEle;
}
// ----------------------------------------------------------------------------------------
bool SimpleNtupleEoverP::LooseEle (const edm::Event & iEvent, const edm::EventSetup & iESetup,const int &iEle){
	
	//************* ELECTRONS
	Handle<View<reco::GsfElectron> > electronHandle;
	iEvent.getByLabel(EleTag_,electronHandle);
	View<reco::GsfElectron> electrons = *electronHandle;
	
	//********* CONVERSION TOOLS
	edm::Handle<reco::ConversionCollection> conversions_h;
	iEvent.getByLabel(conversionsInputTag_, conversions_h);
	
	bool isLooseEle = false; 
	
	reco::GsfElectron electron = electrons.at(iEle);
	
	edm::InputTag EleBad = edm::InputTag("gsfElectrons");
	
	if(EleTag_== EleBad && !electron.ecalDriven()) return false; 
	
	float pt = electron.pt();
	float eta = electron.eta();
	
	float tkIso  = electron.dr03TkSumPt();
	float emIso  = electron.dr03EcalRecHitSumEt(); 
	float hadIso = electron.dr03HcalDepth1TowerSumEt()+electron.dr03HcalDepth2TowerSumEt();
	
	float AEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, electron.superCluster()->eta(), ElectronEffectiveArea::kEleEAData2011);
	
	// default
	float combIso = tkIso +std::max(emIso + hadIso - rho*AEff,float(0.));
	
	// 2011-like test
	//float combIso = tkIso + emIso + hadIso - rho*3.14159*0.3*0.3;      
	
	int isEB = electron.isEB();
	float sigmaIetaIeta = electron.sigmaIetaIeta();
	float DetaIn        = electron.deltaEtaSuperClusterTrackAtVtx();
	float DphiIn        = electron.deltaPhiSuperClusterTrackAtVtx();
	float HOverE        = electron.hadronicOverEm();
	float ooemoop       = (1.0/electron.ecalEnergy() - electron.eSuperClusterOverP()/electron.ecalEnergy());
	
	int mishits             = electron.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
	int nAmbiguousGsfTracks = electron.ambiguousGsfTracksSize();
	
	
	reco::GsfTrackRef eleTrack  = electron.gsfTrack() ;
	float dxy           = eleTrack->dxy(PVPoint_);  
	float dz            = eleTrack->dz (PVPoint_);
	
	edm::Handle<reco::BeamSpot> BSHandle;
	iEvent.getByLabel(BSTag_, BSHandle);
	const reco::BeamSpot BS = *BSHandle;
	
	
	bool isConverted = ConversionTools::hasMatchedConversion(electron, conversions_h, BS.position());
	
	// default
	if(   (pt > 10.) && (fabs(eta) < 2.5) &&
       ( ( (isEB == 1) && (fabs(DetaIn)  < 0.007) ) || ( (isEB == 0) && (fabs(DetaIn)  < 0.009) ) ) &&
       ( ( (isEB == 1) && (fabs(DphiIn)  < 0.150) ) || ( (isEB == 0) && (fabs(DphiIn)  < 0.100) ) ) &&
       ( ( (isEB == 1) && (sigmaIetaIeta < 0.010) ) || ( (isEB == 0) && (sigmaIetaIeta < 0.030) ) ) &&
       ( ( (isEB == 1) && (HOverE        < 0.120) ) || ( (isEB == 0) && (HOverE        < 0.100) ) ) &&
       ( ( (isEB == 1) && (fabs(ooemoop) < 0.050) ) || ( (isEB == 0) && (fabs(ooemoop) < 0.050) ) ) &&
       ( ( (isEB == 1) && (fabs(dxy)     < 0.020) ) || ( (isEB == 0) && (fabs(dxy)     < 0.020) ) ) &&
       ( ( (isEB == 1) && (fabs(dz)      < 0.200) ) || ( (isEB == 0) && (fabs(dz)      < 0.200) ) ) &&
       ( ( (isEB == 1) && (!isConverted) ) || ( (isEB == 0) && (!isConverted) ) ) &&
       ( mishits == 0 ) &&
       ( nAmbiguousGsfTracks == 0 ) &&
       ( ( (isEB == 1 ) && (combIso/pt    < 0.150) ) || ( (isEB == 0 ) && (combIso/pt    < 0.100) ) ) )
		isLooseEle=true;;
	
	// 2011-like test
	//if(   (pt > 10.) && (fabs(eta) < 2.5) &&
	//      ( ( (isEB == 1) && (fabs(DetaIn)  < 0.007) ) || ( (isEB == 0) && (fabs(DetaIn)  < 0.010) ) ) &&
	//      ( ( (isEB == 1) && (fabs(DphiIn)  < 0.800) ) || ( (isEB == 0) && (fabs(DphiIn)  < 0.700) ) ) &&
	//      ( ( (isEB == 1) && (sigmaIetaIeta < 0.010) ) || ( (isEB == 0) && (sigmaIetaIeta < 0.030) ) ) &&
	//      ( ( (isEB == 1) && (HOverE        < 0.150) ) || ( (isEB == 0) && (HOverE        < 0.070) ) ) &&
	//      ( mishits == 0 ) &&
	//      ( nAmbiguousGsfTracks == 0 ) &&
	//      ( ( (isEB == 1 ) && (combIso/pt    < 0.150) ) || ( (isEB == 0 ) && (combIso/pt    < 0.100) ) ) )
	//  isLooseEle=true;;
	
	
	return isLooseEle;
}

// -----------------------------------------------------------------------------------------

bool SimpleNtupleEoverP::myWselection (const edm::Event & iEvent, const edm::EventSetup & iSetup)
{ 
	// default
	float combIso = ele1_tkIso +std::max(ele1_emIso + ele1_hadIso - rho*ele1_effAreaForIso,float(0.));
	
	// 2011-like test
	//float combIso = ele1_tkIso + ele1_emIso + ele1_hadIso - rho*3.14159*0.3*0.3;
	
	if( ele1_pt < 30. ) return false;
	
	// default
	if( ( ele1_isEB == 1 ) && ( combIso/ele1_pt > 0.05 ) ) return false;
	if( ( ele1_isEB == 1 ) && ( fabs(ele1_DetaIn) > 0.004 ) ) return false;
	if( ( ele1_isEB == 1 ) && ( fabs(ele1_DphiIn) > 0.030 ) ) return false;
	if( ( ele1_isEB == 1 ) && ( ele1_sigmaIetaIeta > 0.010 ) ) return false;
	if( ( ele1_isEB == 1 ) && ( ele1_HOverE > 0.120 ) ) return false;
	if( ( ele1_isEB == 1 ) && ( fabs(ele1_ooemoop) > 0.050 ) ) return false;
	
	if( ( ele1_isEB == 0 ) && ( combIso/ele1_pt > 0.035 ) ) return false;
	if( ( ele1_isEB == 0 ) && ( fabs(ele1_DetaIn) > 0.005 ) ) return false;
	if( ( ele1_isEB == 0 ) && ( fabs(ele1_DphiIn) > 0.020 ) ) return false;
	if( ( ele1_isEB == 0 ) && ( ele1_sigmaIetaIeta > 0.030 ) ) return false;
	if( ( ele1_isEB == 0 ) && ( ele1_HOverE > 0.100 ) ) return false;
	if( ( ele1_isEB == 0 ) && ( fabs(ele1_ooemoop) > 0.050 ) ) return false;
	
	// 2011-like test
	//if( ( ele1_isEB == 1 ) && ( combIso/ele1_pt > 0.04 ) ) return false;
	//if( ( ele1_isEB == 1 ) && ( fabs(ele1_DetaIn) > 0.004 ) ) return false;
	//if( ( ele1_isEB == 1 ) && ( fabs(ele1_DphiIn) > 0.030 ) ) return false;
	//if( ( ele1_isEB == 1 ) && ( ele1_HOverE > 0.25 ) ) return false;
	
	//if( ( ele1_isEB == 0 ) && ( combIso/ele1_pt > 0.03 ) ) return false;
	//if( ( ele1_isEB == 0 ) && ( fabs(ele1_DetaIn) > 0.005 ) ) return false;
	//if( ( ele1_isEB == 0 ) && ( fabs(ele1_DphiIn) > 0.020 ) ) return false;
	//if( ( ele1_isEB == 0 ) && ( ele1_HOverE > 0.025 ) ) return false;
	
	if( met_et       < 25.00 ) return false;
	if( ele1Met_mt   < 50.00 ) return false;
	if( ele1Met_Dphi <  1.57 ) return false;
	
	return true;
	
}

// -----------------------------------------------------------------------------------------

bool SimpleNtupleEoverP::myZselection (const edm::Event & iEvent, const edm::EventSetup & iSetup)
{ 
	if( met_et     >  40. ) return false;
	if( ele1ele2_m <  60. ) return false;
	if( ele1ele2_m > 120. ) return false;
	if( (ele1_charge * ele2_charge) != -1. ) return false;
	
	return true;
	
}

// -----------------------------------------------------------------------------------------

void SimpleNtupleEoverP::fillPVInfo (const edm::Event & iEvent, const edm::EventSetup & iESetup){
	
	edm::Handle<reco::VertexCollection> vertexes;
	iEvent.getByLabel(PVTag_, vertexes);
	PV_n = vertexes -> size();
	
	reco::Vertex PV;
	bool PVfound = (vertexes -> size() != 0);
	if(PVfound){
		
		PV = vertexes->at(0);
		PV_z = PV.z();
		PV_d0 = PV.position().Rho();
	}
	else {
		//creating a dummy PV
		PV_z=-999.;
		PV_d0=-999.;
	}  
	
	math::XYZPoint PVPoint(PV.position().x(), PV.position().y(), PV.position().z());
	PVPoint_=PVPoint;
	
}
//------------------------------------------------------------------------------------------------------------
void SimpleNtupleEoverP::fillRhoInfo(const edm::Event & iEvent, const edm::EventSetup & iSetup) 
{
	
	Handle< double > rhoHandle;
	iEvent.getByLabel(rhoTag_,rhoHandle);
	rho = *rhoHandle;  
}

// ------------------------------------------------------------------------------------------------------------

void SimpleNtupleEoverP::fillEleInfo(const edm::Event & iEvent, const edm::EventSetup & iSetup, const int iEle, const std::string eleName)
{
	if( verbosity_ )
		std::cout<< ">>> SimpleNtupleEoverP::fillEleInfo start <<<" << std::endl;
	
	//*********** TRACKER GEOMETRY                             
	edm::ESHandle<TrackerGeometry> pDD ;
	iSetup.get<TrackerDigiGeometryRecord> ().get(pDD);
	
	//*********** MAGNETIC FIELD                               
	edm::ESHandle<MagneticField> theMagField ;
	iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
	
	//*********** CALO TOPOLOGY
	edm::ESHandle<CaloTopology> pTopology;
	iSetup.get<CaloTopologyRecord>().get(pTopology);
	const CaloTopology *topology = pTopology.product();
	
	//*********** IC CONSTANTS
	edm::ESHandle<EcalIntercalibConstants> theICConstants;
	iSetup.get<EcalIntercalibConstantsRcd>().get(theICConstants);
	const EcalIntercalibConstantMap& ICMap = theICConstants->getMap();
	
	//*********** ADCToGeV
	edm::ESHandle<EcalADCToGeVConstant> theADCToGeV;
	iSetup.get<EcalADCToGeVConstantRcd>().get(theADCToGeV);
	
	//*********** LASER ALPHAS
	edm::ESHandle<EcalLaserAlphas> theEcalLaserAlphas;
	iSetup.get<EcalLaserAlphasRcd>().get(theEcalLaserAlphas);
	const EcalLaserAlphaMap* theEcalLaserAlphaMap = theEcalLaserAlphas.product();
	
	//*********** LASER CORRECTION
	edm::ESHandle<EcalLaserDbService> theLaser;
	iSetup.get<EcalLaserDbRecord>().get(theLaser);
	
	//*********** EB DIGI HITS
	edm::Handle<EBDigiCollection> ebDigis;
	if(saveRecHitMatrix_){
		iEvent.getByLabel (digiCollection_EB_, ebDigis) ;
		if (! (ebDigis.isValid ()) ) {
			std::cerr << "EcalValidation::analyze -->  ebDigis not found" << std::endl;
		}
	}
	//*********** EE DIGI HITS
	edm::Handle<EEDigiCollection> eeDigis;
	if(saveRecHitMatrix_){
		iEvent.getByLabel (digiCollection_EE_, eeDigis) ;
		if (! (eeDigis.isValid ()) ) {
			std::cerr << "EcalValidation::analyze -->  eeDigis not found" << std::endl;
		}
	}
	//*********** EB REC HITS
	edm::Handle<EcalRecHitCollection> recHitsEB;
	iEvent.getByLabel( recHitCollection_EB_, recHitsEB );
	const EcalRecHitCollection* theBarrelEcalRecHits = recHitsEB.product () ;
	if ( ! recHitsEB.isValid() ) {
		std::cerr << "SimpleNtupleEoverP::analyze --> recHitsEB not found" << std::endl; 
	}
	
	//*********** EE REC HITS
	edm::Handle<EcalRecHitCollection> recHitsEE;
	iEvent.getByLabel( recHitCollection_EE_, recHitsEE );
	const EcalRecHitCollection* theEndcapEcalRecHits = recHitsEE.product () ;
	if ( ! recHitsEE.isValid() ) {
		std::cerr << "SimpleNtupleEoverP::analyze --> recHitsEE not found" << std::endl; 
	}
	
	//*********** EB SR FLAGS
	edm::Handle<EBSrFlagCollection> SRFlagsEB;
	std::vector<EcalTrigTowerDetId> TTIdList;
	
	if( saveRecHitMatrix_ )
	{
		//    std::cout << ""
		iEvent.getByLabel(SRFlagCollection_EB_, SRFlagsEB );
		for(EBSrFlagCollection::const_iterator it = SRFlagsEB->begin(); it != SRFlagsEB->end(); ++it)
		{
			const int flag = it->value();
			if( flag != EcalSrFlag::SRF_FULL ) continue;
			const EcalTrigTowerDetId TTId = it->id();
			TTIdList.push_back(TTId);
			//std::cout << "flag: " << flag << "   TTId: " << TTId << "   iEta: " << TTId.ieta() << "   iPhi: " << TTId.iphi() << std::endl;
		}
	} 
	
	//*********** EE SR FLAGS
	edm::Handle<EESrFlagCollection> SRFlagsEE;
	std::vector<EcalScDetId> SCIdList;
	if( saveRecHitMatrix_ )
	{
		iEvent.getByLabel(SRFlagCollection_EE_, SRFlagsEE );
		for(EESrFlagCollection::const_iterator it = SRFlagsEE->begin(); it != SRFlagsEE->end(); ++it)
		{
			const int flag = it->value();
			if( flag != EcalSrFlag::SRF_FULL ) continue;
			const EcalScDetId SCId = it->id();
			SCIdList.push_back(SCId);
			//std::cout << "flag: " << flag << "   SCId: " << SCId << "   iEta: " << SCId.ieta() << "   iPhi: " << SCId.iphi() << std::endl;
		}
	}      
	
	//************* ELECTRONS
	Handle<View<reco::GsfElectron> > electronHandle;
	iEvent.getByLabel(EleTag_,electronHandle);
	View<reco::GsfElectron> electrons = *electronHandle;
	
	//************* VERTEX COLLECTION
	edm::Handle<reco::VertexCollection> hVertexProduct;
	iEvent.getByLabel("offlinePrimaryVerticesWithBS",hVertexProduct);
	
	//************* CLUSTER LAZY TOOLS
	EcalClusterLazyTools lazyTools(iEvent,iSetup,recHitCollection_EB_,recHitCollection_EE_); 
	
	//************* REGRESSION
	if( !ecorr_.IsInitialized() ){
		//ecorr_.Initialize(iSetup,"gbrv3ele_52x.root");
		//ecorr_.Initialize(iSetup,"/afs/cern.ch/user/b/bendavid/cmspublic/regweights52xV3/gbrv3ele_52x.root");
		ecorr_.Initialize(iSetup,"/afs/crc.nd.edu/user/a/adrozdet/Private/ANALYSIS/CMSSW_6_1_2/src/UserCode/Bicocca/Calibration/EcalCalibNtuple/test/gbrv3ele_52x.root");
	}
	if( !ecorrPho_.IsInitialized()){
		//ecorrPho_.Initialize(iSetup,"gbrv3ph_52x.root");
		//ecorrPho_.Initialize(iSetup,"/afs/cern.ch/user/b/bendavid/cmspublic/regweights52xV3/gbrv3ph_52x.root");
		ecorrPho_.Initialize(iSetup,"/afs/crc.nd.edu/user/a/adrozdet/Private/ANALYSIS/CMSSW_6_1_2/src/UserCode/Bicocca/Calibration/EcalCalibNtuple/test/gbrv3ph_52x.root");
	} 
	
	//************* CLUSTER PU CLEANING TOOLS
	EcalClusterPUCleaningTools cleaningTools(iEvent, iSetup, recHitCollection_EB_, recHitCollection_EE_); 
	
	
	
	// Take the correct ele
	bool printOut = false;
	reco::GsfElectron electron = electrons.at(iEle);
	
	if ( eleName == "ele1")
	{
		edm::InputTag EleBad = edm::InputTag("gsfElectrons");
		
		if(EleTag_== EleBad && !electron.ecalDriven()) return ; 
		
		ele1=electron.p4();
		ele1_charge=electron.charge();
		ele1_p=ele1.P();
		ele1_pt=ele1.Pt();
		ele1_eta=ele1.eta();
		ele1_phi=ele1.phi();
		
		ele1_isEB=electron.isEB();
		ele1_isEBEEGap=electron.isEBEEGap();
		ele1_isEBEtaGap=electron.isEBEtaGap();
		ele1_isEBPhiGap=electron.isEBPhiGap();
		ele1_isEEDeeGap=electron.isEEDeeGap();
		ele1_isEERingGap=electron.isEERingGap();
		
		ele1_sigmaIetaIeta=electron.sigmaIetaIeta();
		ele1_DphiIn=electron.deltaPhiSuperClusterTrackAtVtx();
		ele1_DetaIn=electron.deltaEtaSuperClusterTrackAtVtx();
		ele1_HOverE=electron.hadronicOverEm();
		ele1_ooemoop  = (1.0/electron.ecalEnergy() - electron.eSuperClusterOverP()/electron.ecalEnergy());
		ele1_tkIso=electron.dr03TkSumPt();
		ele1_emIso=electron.dr03EcalRecHitSumEt();
		ele1_hadIso=electron.dr03HcalDepth1TowerSumEt()+electron.dr03HcalDepth2TowerSumEt();
		
		ele1_effAreaForIso = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, electron.superCluster()->eta(), ElectronEffectiveArea::kEleEAData2011);
		
		reco::GsfTrackRef eleTrack  = electron.gsfTrack();
		ele1_dxy_PV = eleTrack->dxy (PVPoint_);
		ele1_dz_PV = eleTrack->dz (PVPoint_);
		ele1_sigmaP =electron.corrections().trackMomentumError;
		
		reco::SuperClusterRef scRef = electron.superCluster();
		const edm::Ptr<reco::CaloCluster>& seedCluster = scRef->seed();
		
		double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
		double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());
		
		Handle< double > rhoHandle;
		iEvent.getByLabel(rhoTag_,rhoHandle);
		std::pair<double,double> cor = ecorr_.CorrectedEnergyWithErrorV3(electron,*hVertexProduct,*rhoHandle,lazyTools,iSetup);
		
		ele1_scERaw=scRef->rawEnergy();
		ele1_scEtRaw=scRef->rawEnergy()*(Rt/R);
		ele1_scEt=scRef->energy()*(Rt/R);
		ele1_scEtaWidth=scRef->etaWidth();
		ele1_scPhiWidth=scRef->phiWidth();
		ele1_scE=scRef->energy();
		ele1_scEta=scRef->eta();
		ele1_scPhi=scRef->phi();
		ele1_scE_regression=cor.first;
		ele1_scEerr_regression = cor.second;
		
		std::pair<double,double> corPho = ecorrPho_.CorrectedEnergyWithErrorV3(electron,*hVertexProduct,*rhoHandle,lazyTools,iSetup);
		ele1_scE_regression_PhotonTuned=corPho.first;
		ele1_scEerr_regression_PhotonTuned = corPho.second;
		
		
		
		EcalClusterLocal ecalLocalCoord;
		float bcLocalEta, bcLocalPhi, bcThetatilt, bcPhitilt;  
		int bcIeta, bcIphi;
		bcLocalEta = 0;
		
		if ( electron.isEB() )
			ecalLocalCoord.localCoordsEB(*seedCluster,iSetup,bcLocalEta,bcLocalPhi,bcIeta,bcIphi,bcThetatilt,bcPhitilt);  
		if ( electron.isEE() )
			ecalLocalCoord.localCoordsEE(*seedCluster,iSetup,bcLocalEta,bcLocalPhi,bcIeta,bcIphi,bcThetatilt,bcPhitilt);
		
		
		ele1_scLocalEta=bcLocalEta;
		ele1_scLocalPhi=bcLocalPhi;
		
		
		// crack correction variables and local containment corrections
		EcalClusterCrackCorrection -> init(iSetup);
		EcalClusterLocalContCorrection -> init(iSetup);
		double crackcor = 1.;
		double localContCorr = 1.;
		
		for(reco::CaloCluster_iterator cIt = electron.superCluster()->clustersBegin();
			cIt != electron.superCluster()->clustersEnd(); ++cIt)
		{
			const reco::CaloClusterPtr cc = *cIt; 
			crackcor *= ( (electron.superCluster()->rawEnergy() + (*cIt)->energy()*(EcalClusterCrackCorrection->getValue(*cc)-1.)) / electron.superCluster()->rawEnergy() );
		}
		localContCorr = EcalClusterLocalContCorrection->getValue(*electron.superCluster(), 1) ;
		
		ele1_scCrackCorr=crackcor;
		ele1_scLocalContCorr=localContCorr;
		
		
		reco::SuperCluster cleanedSC   = cleaningTools.CleanedSuperCluster(0.02, *scRef, iEvent );
		reco::CaloClusterPtr myseed = (*scRef).seed();
		if (  !((myseed->seed()).rawId()) || (myseed->seed()).rawId()==0 )
		{
			ele1_scERaw_PUcleaned=-9999.;
			ele1_scEtaWidth_PUcleaned=-9999.;
			ele1_scPhiWidth_PUcleaned=-9999.;
			ele1_fCorrection_PUcleaned=-9999.;
		}
		else
		{
			ele1_scERaw_PUcleaned=cleanedSC.energy();
			ele1_scEtaWidth_PUcleaned=cleanedSC.etaWidth();
			ele1_scPhiWidth_PUcleaned=cleanedSC.phiWidth();   
			float fCorrCleaned = fClusterCorrections(cleanedSC.energy() + scRef->preshowerEnergy(), cleanedSC.eta(),cleanedSC.phiWidth()/cleanedSC.etaWidth(),params)/(cleanedSC.energy()+ scRef->preshowerEnergy());
			ele1_fCorrection_PUcleaned=fCorrCleaned;
		}
		
		
		ele1_fEta = scRef->energy()/scRef->rawEnergy();
		ele1_tkP = electron.trackMomentumAtVtx().R();
		ele1_tkPt=electron.trackMomentumAtVtx().Rho();
		ele1_fbrem=electron.fbrem();
		ele1_e5x5=electron.e5x5();
		ele1_eSeedBC=(scRef->seed())->energy();
		ele1_es=scRef->preshowerEnergy();
		
		float E3x3 = 0;
		
		if( electron.isEB() ){
			E3x3 = EcalClusterTools::e3x3( *scRef, theBarrelEcalRecHits, topology);
		}
		if( electron.isEE() ){
			E3x3 = EcalClusterTools::e3x3( *scRef, theEndcapEcalRecHits, topology);
		}
		
		ele1_e3x3=E3x3;
		ele1_scNxtal = scRef->hitsAndFractions().size();
		ele1_bcN = electron.basicClustersSize();
		
		float energy=0.;
		int ieta=0;
		int iphi=0;
		int ix=0;
		int iy=0;
		int zside=0;
		float seedICConstant = -1.;
		float seedLaserAlpha = -1.;
		float seedLaserCorrection = -1.;
		
		if(electron.isEB())
		{
			if( printOut )
				std::cout << "*** EB ***" << std::endl;
			
			std::pair<DetId, float> id = EcalClusterTools::getMaximum(seedCluster->hitsAndFractions(), theBarrelEcalRecHits);
			
			// flag
			EcalRecHitCollection::const_iterator it = theBarrelEcalRecHits->find(id.first);
			
			if( it != theBarrelEcalRecHits->end() )
			{
				const EcalRecHit& rh = (*it);
				energy = rh.energy();
				ieta = (EBDetId(id.first)).ieta();
				iphi = (EBDetId(id.first)).iphi();
				ix = -999;
				iy = -999;
				zside = 0;
			}
			
			// intercalib constant
			EcalIntercalibConstantMap::const_iterator ICMapIt = ICMap.find(EBDetId(id.first));
			if( ICMapIt != ICMap.end() )
				seedICConstant = *ICMapIt;
			
			// laser alphas
			EcalLaserAlphaMap::const_iterator italpha = theEcalLaserAlphaMap->find(id.first);
			if( italpha != theEcalLaserAlphaMap->end() )
				seedLaserAlpha = (*italpha);
			
			// laser correction
			seedLaserCorrection = theLaser->getLaserCorrection(EBDetId(id.first), iEvent.time());
		}
		
		else
		{
			if( printOut )
				std::cout << "*** EE ***" << std::endl;
			
			std::pair<DetId, float> id = EcalClusterTools::getMaximum(seedCluster->hitsAndFractions(), theEndcapEcalRecHits);
			
			// flag - OutOfTime
			EcalRecHitCollection::const_iterator it = theEndcapEcalRecHits->find(id.first);
			
			if( it != theEndcapEcalRecHits->end() )
			{
				const EcalRecHit& rh = (*it);
				energy = rh.energy();
				ix = (EEDetId(id.first)).ix();
				iy = (EEDetId(id.first)).iy();
				ieta = -999;
				iphi = -999;
				zside = (EEDetId(id.first)).zside();
			}
			
			// intercalib constant
			EcalIntercalibConstantMap::const_iterator ICMapIt = ICMap.find(EEDetId(id.first));
			if( ICMapIt != ICMap.end() )
				seedICConstant = *ICMapIt;
			
			// laser alphas
			EcalLaserAlphaMap::const_iterator italpha = theEcalLaserAlphaMap->find(id.first);
			if( italpha != theEcalLaserAlphaMap->end() )
				seedLaserAlpha = (*italpha);
			
			// laser correction
			seedLaserCorrection = theLaser->getLaserCorrection(EEDetId(id.first), iEvent.time());
		}
		
		ele1_seedE=energy;
		ele1_seedLaserAlpha=seedLaserAlpha;
		ele1_seedLaserCorr=seedLaserCorrection;
		ele1_seedICConstant=seedICConstant;
		ele1_seedIeta =ieta;
		ele1_seedIphi = iphi;
		ele1_seedIx=ix;
		ele1_seedIy=iy;
		ele1_seedZside=zside;
		ele1_EOverP=electron.eSuperClusterOverP();
		
		// rechit variables
		int numRecHit = 0;
		float sumRecHitE = 0.;
		float sumLaserCorrectionRecHitE = 0.;
		float sumRecHitE5x5 = 0.;
		float sumLaserCorrectionRecHitE5x5 = 0.;
		float sumRecHitE3x3 = 0.;
		float sumLaserCorrectionRecHitE3x3 = 0.;
		
		const std::vector<std::pair<DetId,float> >& hits = scRef->hitsAndFractions();
		
		if( printOut )
		{
			std::cout << "runId: " << iEvent.id().run() 
			<< std::fixed
			<< "   electron eta: " << std::setprecision(2) << std::setw(5) << electron.eta()
			<< "   electron phi: " << std::setprecision(2) << std::setw(5) << electron.phi()
			<< "   SC energy: "    << std::setprecision(2) << std::setw(6) << scRef -> energy()
			<< std::endl;
		} 
		
		
		if(saveRecHitMatrix_)
		{
			float theLaserCorrection = -1.;
			float theICCorrection = -1.;
			
			if(electron.isEB())
			{
				DetId seedId = EcalClusterTools::getMaximum( scRef->hitsAndFractions(), theBarrelEcalRecHits ).first;    
				//save the matrix in case of eleSeed
				std::vector<DetId> rectangle =  EcalClusterTools::matrixDetId(topology, seedId, -9, 9, -9, 9);
				
				int it = 0;
				for(std::vector<DetId>::const_iterator itr = rectangle.begin(); itr != rectangle.end(); ++itr)
				{
					++it;
					EcalRecHitCollection::const_iterator itrRecHit = theBarrelEcalRecHits->find(*itr) ;
					if(itrRecHit == theBarrelEcalRecHits->end()) continue;
					
					// fill recHit variables
					EBDetId barrelId(*itr);
					EcalTrigTowerDetId towerId = barrelId.tower();
					
					// laser correction
					theLaserCorrection = theLaser->getLaserCorrection(barrelId, iEvent.time());
					// IC correction
					EcalIntercalibConstantMap::const_iterator ICMapIt = ICMap.find(barrelId);
					theICCorrection = *ICMapIt;
					
					int SRFlag = 3;
					std::vector<EcalTrigTowerDetId>::iterator TTIdListIt = std::find(TTIdList.begin(),TTIdList.end(),towerId);
					if( TTIdListIt == TTIdList.end() ) SRFlag = 1;
					
					bool digiFound = false;
					for(EBDigiCollection::const_iterator digiItr = ebDigis->begin(); digiItr != ebDigis->end(); ++digiItr)
					{
						if(digiItr->id() != barrelId )continue;
						digiFound = true;
						EcalDataFrame df = *digiItr;
						for(int iSample = 0; iSample < 10; ++iSample)
							ele1_recHitMatrix_samples.push_back(df.sample(iSample).adc());
					}
					if( digiFound == false ) continue;
                    
					ele1_recHitMatrix_E.push_back(itrRecHit->energy());
					ele1_recHitMatrix_flag.push_back(SRFlag*1000+itrRecHit->recoFlag());
					ele1_recHitMatrix_hashedIndex.push_back(barrelId.hashedIndex());
					ele1_recHitMatrix_ietaORix.push_back(barrelId.ieta());
					ele1_recHitMatrix_iphiORiy.push_back(barrelId.iphi());
					ele1_recHitMatrix_zside.push_back(0);
					ele1_recHitMatrix_laserCorrection.push_back(theLaserCorrection);
					ele1_recHitMatrix_ICConstant.push_back(theICCorrection);
				}
			}
			
			else
			{
				DetId seedId = EcalClusterTools::getMaximum( scRef->hitsAndFractions(), theEndcapEcalRecHits ).first;    
				//save the matrix in case of eleSeed
				std::vector<DetId> rectangle =  EcalClusterTools::matrixDetId(topology, seedId, -9, 9, -9, 9);
				
				for(std::vector<DetId>::const_iterator itr = rectangle.begin(); itr != rectangle.end(); ++itr)
				{
					EcalRecHitCollection::const_iterator itrRecHit = theEndcapEcalRecHits->find(*itr) ;
					if(itrRecHit == theEndcapEcalRecHits->end()) continue;
					
					// fill recHit variables
					EEDetId endcapId(*itr);
					EcalScDetId scId = endcapId.sc();
					
					// laser correction
					theLaserCorrection = theLaser->getLaserCorrection(endcapId, iEvent.time());
					// IC correction
					EcalIntercalibConstantMap::const_iterator ICMapIt = ICMap.find(endcapId);
					theICCorrection = *ICMapIt;
					
					int SRFlag = 3;
					std::vector<EcalScDetId>::iterator SCIdListIt = std::find(SCIdList.begin(),SCIdList.end(),scId);
					if( SCIdListIt == SCIdList.end() ) SRFlag = 1;
					
					bool digiFound = false;
					for(EEDigiCollection::const_iterator digiItr = eeDigis->begin(); digiItr != eeDigis->end(); ++digiItr)
					{
						if(digiItr->id() != endcapId) continue;
						digiFound = true;
						EcalDataFrame df = *digiItr;
						for(int iSample = 0; iSample < 10; ++iSample)            
							ele1_recHitMatrix_samples.push_back(df.sample(iSample).adc());
					}
					if( digiFound == false ) continue;
					
					ele1_recHitMatrix_E.push_back(itrRecHit->energy());
					ele1_recHitMatrix_flag.push_back(SRFlag*1000+itrRecHit->recoFlag());
					ele1_recHitMatrix_hashedIndex.push_back(endcapId.hashedIndex());
					ele1_recHitMatrix_ietaORix.push_back(endcapId.ix());
					ele1_recHitMatrix_iphiORiy.push_back(endcapId.iy());
					ele1_recHitMatrix_zside.push_back(0);
					ele1_recHitMatrix_laserCorrection.push_back(theLaserCorrection);
					ele1_recHitMatrix_ICConstant.push_back(theICCorrection);
				}
			}
		}
		
		for(std::vector<std::pair<DetId,float> >::const_iterator rh = hits.begin(); rh!=hits.end(); ++rh)
		{
			float rhLaserCorrection = 1.;
			float rhICCorrection = 1.;
			float theLaserCorrection = -1.;
			float theICCorrection = -1.;
			float theAlpha = -1.;
            
			if ((*rh).first.subdetId()== EcalBarrel)
			{
				EBRecHitCollection::const_iterator itrechit = theBarrelEcalRecHits->find((*rh).first);
				if (itrechit==theBarrelEcalRecHits->end()) continue;
				EBDetId barrelId (itrechit->id ()); 
				++numRecHit;
				
				// laser correction
				theLaserCorrection = theLaser->getLaserCorrection(barrelId, iEvent.time());
				if ( applyCorrections_ ) rhLaserCorrection = theLaserCorrection;
				// IC correction
				EcalIntercalibConstantMap::const_iterator ICMapIt = ICMap.find(barrelId);
				theICCorrection = *ICMapIt;
				if ( applyCorrections_ ) rhICCorrection = theICCorrection;
				// Alpha
				EcalLaserAlphaMap::const_iterator italpha = theEcalLaserAlphaMap->find(itrechit->id());
				if( italpha != theEcalLaserAlphaMap->end() )
					theAlpha = (*italpha);
				sumRecHitE += itrechit->energy() * rhLaserCorrection * rhICCorrection ;
				sumLaserCorrectionRecHitE += theLaserCorrection * itrechit->energy() * rhLaserCorrection * rhICCorrection;
				// check if rh is inside the 5x5 matrix
				if ( fabs(barrelId.ieta() - ele1_seedIeta) < 3 && fabs(barrelId.iphi() - ele1_seedIphi) < 3 ) {
					sumRecHitE5x5 += itrechit->energy() * rhLaserCorrection * rhICCorrection ;
					sumLaserCorrectionRecHitE5x5 += theLaserCorrection * itrechit->energy() * rhLaserCorrection * rhICCorrection;
				}
				// check if rh is inside the 3x3 matrix
				if ( fabs(barrelId.ieta() - ele1_seedIeta) < 1 && fabs(barrelId.iphi() - ele1_seedIphi) < 1 ) {
					sumRecHitE3x3 += itrechit->energy() * rhLaserCorrection * rhICCorrection ;
					sumLaserCorrectionRecHitE3x3 += theLaserCorrection * itrechit->energy() * rhLaserCorrection * rhICCorrection ;
				}
				// fill recHit variables
				ele1_recHit_E.push_back(itrechit->energy() * rhLaserCorrection);
				ele1_recHit_flag.push_back(itrechit->recoFlag());
				ele1_recHit_hashedIndex.push_back(barrelId.hashedIndex());
				ele1_recHit_ietaORix.push_back(barrelId.ieta());
				ele1_recHit_iphiORiy.push_back(barrelId.iphi());
				ele1_recHit_zside.push_back(0);
				ele1_recHit_laserCorrection.push_back(theLaserCorrection);
				ele1_recHit_ICConstant.push_back(theICCorrection);
				ele1_recHit_Alpha.push_back(theAlpha);
				
				if( printOut && itrechit->energy() > 1. )
				{
					std::cout << std::fixed
					<< "    recHitLC: "    << std::setprecision(6) << std::setw(8) << theLaserCorrection
					<< "    recHitIC: "    << std::setprecision(6) << std::setw(8) << theICCorrection
					<< std::endl;
				}
			}
			
			if( (*rh).first.subdetId()== EcalEndcap )
			{
				EERecHitCollection::const_iterator itrechit = theEndcapEcalRecHits->find((*rh).first);
				if (itrechit==theEndcapEcalRecHits->end()) continue;
				EEDetId endcapId (itrechit->id ()); 
				++numRecHit;
				
				// laser correction
				theLaserCorrection = theLaser->getLaserCorrection(endcapId, iEvent.time());
				if ( applyCorrections_ ) rhLaserCorrection = theLaserCorrection;
				// IC correction
				EcalIntercalibConstantMap::const_iterator ICMapIt = ICMap.find(endcapId);
				theICCorrection = *ICMapIt;
				if ( applyCorrections_ ) rhICCorrection = theICCorrection;
				
				sumRecHitE += itrechit->energy() * rhLaserCorrection * rhICCorrection ;
				sumLaserCorrectionRecHitE += theLaserCorrection * itrechit->energy() * rhLaserCorrection * rhICCorrection;
				// check if rh is inside the 5x5 matrix
				if ( fabs(endcapId.ix() - ele1_seedIx) < 3 && fabs(endcapId.iy() - ele1_seedIy) < 3 ) {
					sumRecHitE5x5 += itrechit->energy() * rhLaserCorrection * rhICCorrection ;
					sumLaserCorrectionRecHitE5x5 += theLaserCorrection * itrechit->energy() * rhLaserCorrection * rhICCorrection;
				}
				// check if rh is inside the 3x3 matrix
				if ( fabs(endcapId.ix() - ele1_seedIx) < 1 && fabs(endcapId.iy() - ele1_seedIy) < 1 ) {
					sumRecHitE3x3 += itrechit->energy() * rhLaserCorrection * rhICCorrection ;
					sumLaserCorrectionRecHitE3x3 += theLaserCorrection * itrechit->energy() * rhLaserCorrection * rhICCorrection ;
				}
				// fill recHit variables
				ele1_recHit_E.push_back(itrechit->energy() * rhLaserCorrection);
				ele1_recHit_flag.push_back(itrechit->recoFlag());
				ele1_recHit_hashedIndex.push_back(endcapId.hashedIndex());
				ele1_recHit_ietaORix.push_back(endcapId.ix());
				ele1_recHit_iphiORiy.push_back(endcapId.iy());
				ele1_recHit_zside.push_back(endcapId.zside());
				ele1_recHit_laserCorrection.push_back(theLaserCorrection);
				ele1_recHit_ICConstant.push_back(theICCorrection);
				
				if( printOut && itrechit->energy() > 1. )
				{
					std::cout << std::fixed
					<< "    recHitLC: "    << std::setprecision(6) << std::setw(8) << theLaserCorrection
					<< "    recHitIC: "    << std::setprecision(6) << std::setw(8) << theICCorrection
					<< std::endl;
				}
			}
			
		}
		
		ele1_nRecHits = numRecHit; 
		
		ele1_scLaserCorr = sumLaserCorrectionRecHitE/sumRecHitE;
		if ( applyCorrections_ ) ele1_scE = scRef->energy()*ele1_fEta*sumRecHitE;
		else ele1_scE = scRef->energy();
		
		ele1_fEtaCorr = fClusterCorrections(sumRecHitE+ele1_es,ele1_scEta,scRef->phiWidth()/scRef->etaWidth(),params)/fClusterCorrections(scRef->rawEnergy()+ele1_es,ele1_scEta,scRef->phiWidth()/scRef->etaWidth(),params);
		
		ele1_5x5LaserCorr = sumLaserCorrectionRecHitE5x5/sumRecHitE5x5;
		
		ele1_3x3LaserCorr = sumLaserCorrectionRecHitE3x3/sumRecHitE3x3;
		
		/// add regression input variables
		reco::SuperClusterRef s = electron.superCluster();
		reco::CaloClusterPtr b = s->seed(); //seed  basic cluster
		reco::CaloClusterPtr b2;
		reco::CaloClusterPtr bclast;
		reco::CaloClusterPtr bclast2;
		bool isbarrel =  b->hitsAndFractions().at(0).first.subdetId()==EcalBarrel;
		
		
		ele1_eRegrInput_rawE = s->rawEnergy();
		ele1_eRegrInput_r9 = lazyTools.e3x3(*b)/s->rawEnergy();
		ele1_eRegrInput_eta = s->eta();
		ele1_eRegrInput_phi = s->phi();
		ele1_eRegrInput_r25 = lazyTools.e5x5(*b)/s->rawEnergy();
		ele1_eRegrInput_hoe = electron.hcalOverEcal();
		ele1_eRegrInput_etaW = s->etaWidth() ;
		ele1_eRegrInput_phiW = s->phiWidth() ;
		ele1_eRegrInput_SCsize = s->clustersSize() ;
		ele1_eRegrInput_rho = *rhoHandle;
		ele1_eRegrInput_nPV = hVertexProduct->size();
		
		//seed basic cluster variables
		double bemax = lazyTools.eMax(*b);
		double be2nd = lazyTools.e2nd(*b);
		double betop = lazyTools.eTop(*b);
		double bebottom = lazyTools.eBottom(*b);
		double beleft = lazyTools.eLeft(*b);
		double beright = lazyTools.eRight(*b);
		
		double be2x5max = lazyTools.e2x5Max(*b);
		double be2x5top = lazyTools.e2x5Top(*b);
		double be2x5bottom = lazyTools.e2x5Bottom(*b);
		double be2x5left = lazyTools.e2x5Left(*b);
		double be2x5right = lazyTools.e2x5Right(*b);
		
		ele1_eRegrInput_Deta_bC_sC = b->eta()-s->eta();
		ele1_eRegrInput_Dphi_bC_sC = reco::deltaPhi(b->phi(),s->phi());
		ele1_eRegrInput_bCE_Over_sCE = b->energy()/b->energy();
		ele1_eRegrInput_e3x3_Over_bCE = lazyTools.e3x3(*b)/b->energy();
		ele1_eRegrInput_e5x5_Over_bCE = lazyTools.e5x5(*b)/b->energy();
		ele1_eRegrInput_sigietaieta_bC1 = sqrt(lazyTools.localCovariances(*b)[0]);
		ele1_eRegrInput_sigiphiiphi_bC1 = sqrt(lazyTools.localCovariances(*b)[2]);
		ele1_eRegrInput_sigietaiphi_bC1 = lazyTools.localCovariances(*b)[1];
		ele1_eRegrInput_bEMax_Over_bCE = bemax/b->energy();
		ele1_eRegrInput_bE2nd_Over_bCE = be2nd/b->energy();
		ele1_eRegrInput_bEtop_Over_bCE = betop/b->energy();
		ele1_eRegrInput_bEbot_Over_bCE = bebottom/b->energy();
		ele1_eRegrInput_bEleft_Over_bCE = beleft/b->energy();
		ele1_eRegrInput_bEright_Over_bCE = beright/b->energy();
		ele1_eRegrInput_be2x5max_Over_bCE = be2x5max/b->energy();
		ele1_eRegrInput_be2x5top_Over_bCE = be2x5top/b->energy();
		ele1_eRegrInput_be2x5bottom_Over_bCE = be2x5bottom/b->energy();
		ele1_eRegrInput_be2x5left_Over_bCE = be2x5left/b->energy();
		ele1_eRegrInput_be2x5right_Over_bCE = be2x5right/b->energy();
		
		
		
		if( isbarrel )
		{
			// seed cluster
			float betacry, bphicry, bthetatilt, bphitilt;
			int bieta, biphi;
			EcalClusterLocal _ecalLocal;
			_ecalLocal.localCoordsEB(*b,iSetup,betacry,bphicry,bieta,biphi,bthetatilt,bphitilt);
			
			ele1_eRegrInput_seedbC_eta = bieta;
			ele1_eRegrInput_seedbC_phi = biphi;
			ele1_eRegrInput_seedbC_eta_p5 = bieta%5;
			ele1_eRegrInput_seedbC_phi_p2 = biphi%2;
			ele1_eRegrInput_seedbC_bieta = (TMath::Abs(bieta)<=25)*(bieta%25) + (TMath::Abs(bieta)>25)*((bieta-25*TMath::Abs(bieta)/bieta)%20);
			ele1_eRegrInput_seedbC_phi_p20 = biphi%20;
			ele1_eRegrInput_seedbC_etacry = betacry;
			ele1_eRegrInput_seedbC_phicry = bphicry;
			ele1_eRegrInput_ESoSC = -99. ;
		}
		else
		{
			ele1_eRegrInput_ESoSC = s->preshowerEnergy()/s->rawEnergy();
			ele1_eRegrInput_seedbC_eta = -99.;
			ele1_eRegrInput_seedbC_phi = -99.;
			ele1_eRegrInput_seedbC_eta_p5 = -99.;
			ele1_eRegrInput_seedbC_phi_p2 = -99.;
			ele1_eRegrInput_seedbC_bieta = -99.;
			ele1_eRegrInput_seedbC_phi_p20 = -99.;
			ele1_eRegrInput_seedbC_etacry = -99.;
			ele1_eRegrInput_seedbC_phicry =-99.;
		}
		
		if(saveFbrem_){
			reco::GsfTrackRef eleTrack  = electron.gsfTrack();
			reco::TrackRef eleNTrack = electron.closestTrack();
			GlobalPoint outPos(eleTrack->extra()->outerPosition().x(), eleTrack->extra()->outerPosition().y(), eleTrack->extra()->outerPosition().z());
			GlobalPoint innPos(eleTrack->extra()->innerPosition().x(), eleTrack->extra()->innerPosition().y(), eleTrack->extra()->innerPosition().z());
			std::vector<reco::GsfTangent> eleTangent = eleTrack->gsfExtra()->tangents();
			int numberOfValidHits_Trk = 0;
			if(!eleNTrack.isNull())
				for(unsigned int trkNH = 0; trkNH < eleNTrack->extra()->recHitsSize(); ++trkNH){
					if(!eleNTrack->extra()->recHit(trkNH)->isValid()) continue;
					DetId id = eleNTrack->extra()->recHit(trkNH)->geographicalId();
					const GeomDet* det = pDD->idToDet(id);
					
					if(eleTrack->extra()->seedRef().isNull()) continue;
					edm::RefToBase<TrajectorySeed> seed = eleTrack->extra()->seedRef();
					ElectronSeedRef elseed=seed.castTo<ElectronSeedRef>();
					TrajectoryStateOnSurface t = trajectoryStateTransform::transientState(elseed->startingState(), &(det->surface()), &(*theMagField));
					if(!t.isValid()) continue;
					
					GlobalPoint hitPos = t.globalPosition();
					if(hitPos.x() == ele1_inner_x && hitPos.y() == ele1_inner_y && hitPos.z() == ele1_inner_z){
						ele1_inner_p = (sqrt(eleTrack->extra()->innerMomentum().Mag2()) );
						ele1_inner_x = innPos.x();
						ele1_inner_y = innPos.y();
						ele1_inner_z = innPos.z();
					}
					if(hitPos.x() == ele1_outer_x && hitPos.y() == ele1_outer_y && hitPos.z() == ele1_outer_z){
						ele1_outer_p = (sqrt(eleTrack->extra()->outerMomentum().Mag2()) );
						ele1_outer_x = outPos.x();
						ele1_outer_y = outPos.y();
						ele1_outer_z = outPos.z();
					}
					for(unsigned int pp=0; pp<eleTangent.size(); ++pp ){
						GlobalPoint tangPos( eleTangent.at(pp).position().x(),
											eleTangent.at(pp).position().y(),
											eleTangent.at(pp).position().z() );
						if(hitPos.x() != tangPos.x() && hitPos.y() != tangPos.y() && hitPos.z() != tangPos.z()) continue;
						float tangMom = sqrt(eleTangent.at(pp).momentum().Mag2());
						ele1_tangent_p.push_back(tangMom);
						ele1_tangent_x.push_back(tangPos.x());
						ele1_tangent_y.push_back(tangPos.y());
						ele1_tangent_z.push_back(tangPos.z());
						ele1_tangent_dP.push_back(eleTangent.at(pp).deltaP().value());
						ele1_tangent_dPerr.push_back(eleTangent.at(pp).deltaP().error());
						++numberOfValidHits_Trk;
					}
					//      ele1_tangent_n = eleTangent.size();                                                          
				}
			ele1_tangent_n = numberOfValidHits_Trk;
		}
	}
	
	if ( eleName == "ele2" )
	{
		edm::InputTag EleBad = edm::InputTag("gsfElectrons");
		
		if(EleTag_== EleBad && !electron.ecalDriven()) return ; 
		
		ele2=electron.p4();
		ele2_charge=electron.charge();
		ele2_p=ele2.P();
		ele2_pt=ele2.Pt();
		ele2_eta=ele2.eta();
		ele2_phi=ele2.phi();
		
		ele2_isEB=electron.isEB();
		ele2_isEBEEGap=electron.isEBEEGap();
		ele2_isEBEtaGap=electron.isEBEtaGap();
		ele2_isEBPhiGap=electron.isEBPhiGap();
		ele2_isEEDeeGap=electron.isEEDeeGap();
		ele2_isEERingGap=electron.isEERingGap();
		
		ele2_sigmaIetaIeta=electron.sigmaIetaIeta();
		ele2_DphiIn=electron.deltaPhiSuperClusterTrackAtVtx();
		ele2_DetaIn=electron.deltaEtaSuperClusterTrackAtVtx();
		ele2_HOverE=electron.hadronicOverEm();
		ele2_ooemoop=  (1.0/electron.ecalEnergy() - electron.eSuperClusterOverP()/electron.ecalEnergy());
		
		ele2_tkIso=electron.dr03TkSumPt();
		ele2_emIso=electron.dr03EcalRecHitSumEt();
		ele2_hadIso=electron.dr03HcalDepth1TowerSumEt()+electron.dr03HcalDepth2TowerSumEt();
		
		ele2_effAreaForIso = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, electron.superCluster()->eta(), ElectronEffectiveArea::kEleEAData2011);
		
		reco::GsfTrackRef eleTrack  = electron.gsfTrack() ;
		ele2_dxy_PV = eleTrack->dxy (PVPoint_);
		ele2_dz_PV = eleTrack->dz (PVPoint_);
		ele2_sigmaP =electron.corrections().trackMomentumError;
		
		
		reco::SuperClusterRef scRef = electron.superCluster();
		const edm::Ptr<reco::CaloCluster>& seedCluster = scRef->seed();
		
		double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
		double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());
		
		Handle< double > rhoHandle;
		iEvent.getByLabel(rhoTag_,rhoHandle);
		std::pair<double,double> cor = ecorr_.CorrectedEnergyWithErrorV3(electron,*hVertexProduct,*rhoHandle,lazyTools,iSetup);
		
		
		ele2_scERaw=scRef->rawEnergy();
		ele2_scEtRaw=scRef->rawEnergy()*(Rt/R);
		ele2_scEt=scRef->energy()*(Rt/R);
		ele2_scEtaWidth=scRef->etaWidth();
		ele2_scPhiWidth=scRef->phiWidth();
		ele2_scE=scRef->energy();
		ele2_scEta=scRef->eta();
		ele2_scPhi=scRef->phi();
		ele2_scE_regression=cor.first;
		ele2_scEerr_regression = cor.second;
		
		std::pair<double,double> corPho = ecorrPho_.CorrectedEnergyWithErrorV3(electron,*hVertexProduct,*rhoHandle,lazyTools,iSetup);
		ele2_scE_regression_PhotonTuned=corPho.first;
		ele2_scEerr_regression_PhotonTuned = corPho.second;
		
		
		EcalClusterLocal ecalLocalCoord;
		float bcLocalEta, bcLocalPhi, bcThetatilt, bcPhitilt;  
		int bcIeta, bcIphi;
		bcLocalEta = 0;
		
		if ( electron.isEB() )
			ecalLocalCoord.localCoordsEB(*seedCluster,iSetup,bcLocalEta,bcLocalPhi,bcIeta,bcIphi,bcThetatilt,bcPhitilt);  
		if ( electron.isEE() )
			ecalLocalCoord.localCoordsEE(*seedCluster,iSetup,bcLocalEta,bcLocalPhi,bcIeta,bcIphi,bcThetatilt,bcPhitilt);
		
		
		ele2_scLocalEta=bcLocalEta;
		ele2_scLocalPhi=bcLocalPhi;
		
		
		// crack correction variables and local containment corrections
		EcalClusterCrackCorrection -> init(iSetup);
		EcalClusterLocalContCorrection -> init(iSetup);
		double crackcor = 1.;
		double localContCorr = 1.;
		
		for(reco::CaloCluster_iterator cIt = electron.superCluster()->clustersBegin();
			cIt != electron.superCluster()->clustersEnd(); ++cIt)
		{
			const reco::CaloClusterPtr cc = *cIt; 
			crackcor *= ( (electron.superCluster()->rawEnergy() + (*cIt)->energy()*(EcalClusterCrackCorrection->getValue(*cc)-1.)) / electron.superCluster()->rawEnergy() );
		}
		localContCorr = EcalClusterLocalContCorrection->getValue(*electron.superCluster(), 1) ;
		
		ele2_scCrackCorr=crackcor;
		ele2_scLocalContCorr=localContCorr;
		
		
		reco::SuperCluster cleanedSC   = cleaningTools.CleanedSuperCluster(0.02, *scRef, iEvent );
		reco::CaloClusterPtr myseed = (*scRef).seed();
		if (  !((myseed->seed()).rawId()) || (myseed->seed()).rawId()==0 )
		{
			ele2_scERaw_PUcleaned=-9999.;
			ele2_scEtaWidth_PUcleaned=-9999.;
			ele2_scPhiWidth_PUcleaned=-9999.;
			ele2_fCorrection_PUcleaned=-9999.;
		}
		else
		{
			ele2_scERaw_PUcleaned=cleanedSC.energy();
			ele2_scEtaWidth_PUcleaned=cleanedSC.etaWidth();
			ele2_scPhiWidth_PUcleaned=cleanedSC.phiWidth();   
			float fCorrCleaned = fClusterCorrections(cleanedSC.energy() + scRef->preshowerEnergy(), cleanedSC.eta(),cleanedSC.phiWidth()/cleanedSC.etaWidth(),params)/(cleanedSC.energy()+ scRef->preshowerEnergy());
			ele2_fCorrection_PUcleaned=fCorrCleaned;
		}
		
		
		ele2_fEta = scRef->energy()/scRef->rawEnergy();
		ele2_tkP = electron.trackMomentumAtVtx().R();
		ele2_tkPt=electron.trackMomentumAtVtx().Rho();
		ele2_fbrem=electron.fbrem();
		ele2_e5x5=electron.e5x5();
		ele2_eSeedBC=(scRef->seed())->energy();
		ele2_es=scRef->preshowerEnergy();
		
		float E3x3 = 0;
		
		if( electron.isEB() ){
			E3x3 = EcalClusterTools::e3x3( *scRef, theBarrelEcalRecHits, topology);
		}
		
		if ( electron.isEE() ){
			E3x3 = EcalClusterTools::e3x3( *scRef, theEndcapEcalRecHits, topology);
		}
		
		ele2_e3x3=E3x3;
		ele2_scNxtal = scRef->hitsAndFractions().size();
		ele2_bcN = electron.basicClustersSize();
		
		float energy=0.;
		int ieta=0;
		int iphi=0;
		int ix=0;
		int iy=0;
		int zside=0;
		float seedICConstant = -1.;
		float seedLaserAlpha = -1.;
		float seedLaserCorrection = -1.;
		
		if(electron.isEB())
		{
			std::pair<DetId, float> id = EcalClusterTools::getMaximum(seedCluster->hitsAndFractions(), theBarrelEcalRecHits);
			
			// flag
			EcalRecHitCollection::const_iterator it = theBarrelEcalRecHits->find(id.first);
			
			if( it != theBarrelEcalRecHits->end() )
			{
				const EcalRecHit& rh = (*it);
				energy = rh.energy();
				ieta = (EBDetId(id.first)).ieta();
				iphi = (EBDetId(id.first)).iphi();
				ix = -999;
				iy = -999;
				zside = 0;
			}
			
			// intercalib constant
			EcalIntercalibConstantMap::const_iterator ICMapIt = ICMap.find(EBDetId(id.first));
			if( ICMapIt != ICMap.end() )
				seedICConstant = *ICMapIt;
			
			// laser alphas
			EcalLaserAlphaMap::const_iterator italpha = theEcalLaserAlphaMap->find(id.first);
			if( italpha != theEcalLaserAlphaMap->end() )
				seedLaserAlpha = (*italpha);
			
			// laser correction
			seedLaserCorrection = theLaser->getLaserCorrection(EBDetId(id.first), iEvent.time());
		}
		
		else
		{
			std::pair<DetId, float> id = EcalClusterTools::getMaximum(seedCluster->hitsAndFractions(), theEndcapEcalRecHits);
			
			// flag - OutOfTime
			EcalRecHitCollection::const_iterator it = theEndcapEcalRecHits->find(id.first);
			
			if( it != theEndcapEcalRecHits->end() )
			{
				const EcalRecHit& rh = (*it);
				energy = rh.energy();
				ix = (EEDetId(id.first)).ix();
				iy = (EEDetId(id.first)).iy();
				ieta = -999;
				iphi = -999;
				zside = (EEDetId(id.first)).zside();
			}
			
			// intercalib constant
			EcalIntercalibConstantMap::const_iterator ICMapIt = ICMap.find(EEDetId(id.first));
			if( ICMapIt != ICMap.end() )
				seedICConstant = *ICMapIt;
			
			// laser alphas
			EcalLaserAlphaMap::const_iterator italpha = theEcalLaserAlphaMap->find(id.first);
			if( italpha != theEcalLaserAlphaMap->end() )
				seedLaserAlpha = (*italpha);
			
			// laser correction
			seedLaserCorrection = theLaser->getLaserCorrection(EEDetId(id.first), iEvent.time());
		}
		
		ele2_seedE=energy;
		ele2_seedLaserAlpha=seedLaserAlpha;
		ele2_seedLaserCorr=seedLaserCorrection;
		ele2_seedICConstant=seedICConstant;
		ele2_seedIeta =ieta;
		ele2_seedIphi = iphi;
		ele2_seedIx=ix;
		ele2_seedIy=iy;
		ele2_seedZside=zside;
		ele2_EOverP=electron.eSuperClusterOverP();
		
		// rechit variables
		int numRecHit = 0;
		float sumRecHitE = 0.;
		float sumLaserCorrectionRecHitE = 0.;
		float sumRecHitE5x5 = 0.;
		float sumLaserCorrectionRecHitE5x5 = 0.;
		float sumRecHitE3x3 = 0.;
		float sumLaserCorrectionRecHitE3x3 = 0.;
		
		bool printOut = false;
		const std::vector<std::pair<DetId,float> >& hits = scRef->hitsAndFractions();
		
		if( printOut )
		{
			std::cout << "runId: " << iEvent.id().run() 
			<< std::fixed
			<< "   electron eta: " << std::setprecision(2) << std::setw(5) << electron.eta()
			<< "   electron phi: " << std::setprecision(2) << std::setw(5) << electron.phi()
			<< "   SC energy: "    << std::setprecision(2) << std::setw(6) << scRef -> energy()
			<< std::endl;
		} 
		
		if(saveRecHitMatrix_)
		{
			float theLaserCorrection = -1.;
			float theICCorrection = -1.;
			
			if(electron.isEB())
			{
				DetId seedId = EcalClusterTools::getMaximum( scRef->hitsAndFractions(), theBarrelEcalRecHits ).first;    
				//save the matrix in case of eleSeed
				std::vector<DetId> rectangle =  EcalClusterTools::matrixDetId(topology, seedId, -9, 9, -9, 9);
				
				for(std::vector<DetId>::const_iterator itr = rectangle.begin(); itr != rectangle.end(); ++itr)
				{
					EcalRecHitCollection::const_iterator itrRecHit = theBarrelEcalRecHits->find(*itr) ;
					if(itrRecHit == theBarrelEcalRecHits->end()) continue;
					
					// fill recHit variables
					EBDetId barrelId(*itr);
					EcalTrigTowerDetId towerId = barrelId.tower();
					
					// laser correction
					theLaserCorrection = theLaser->getLaserCorrection(barrelId, iEvent.time());
					// IC correction
					EcalIntercalibConstantMap::const_iterator ICMapIt = ICMap.find(barrelId);
					theICCorrection = *ICMapIt;
					
					int SRFlag = 3;
					std::vector<EcalTrigTowerDetId>::iterator TTIdListIt = std::find(TTIdList.begin(),TTIdList.end(),towerId);
					if( TTIdListIt == TTIdList.end() ) SRFlag = 1;
					
					bool digiFound = false;
					for(EBDigiCollection::const_iterator digiItr = ebDigis->begin(); digiItr != ebDigis->end(); ++digiItr)
					{
						if(digiItr->id() != barrelId )continue;
						digiFound = true;
						EcalDataFrame df = *digiItr;
						for(int iSample = 0; iSample < 10; ++iSample)
							ele2_recHitMatrix_samples.push_back(df.sample(iSample).adc());
					}
					if( digiFound == false ) continue;
					
					ele2_recHitMatrix_E.push_back(itrRecHit->energy());
					ele2_recHitMatrix_flag.push_back(SRFlag*1000+itrRecHit->recoFlag());
					ele2_recHitMatrix_hashedIndex.push_back(barrelId.hashedIndex());
					ele2_recHitMatrix_ietaORix.push_back(barrelId.ieta());
					ele2_recHitMatrix_iphiORiy.push_back(barrelId.iphi());
					ele2_recHitMatrix_zside.push_back(0);
					ele2_recHitMatrix_laserCorrection.push_back(theLaserCorrection);
					ele2_recHitMatrix_ICConstant.push_back(theICCorrection);
					
				}
			}
			
			else
			{
				DetId seedId = EcalClusterTools::getMaximum( scRef->hitsAndFractions(), theEndcapEcalRecHits ).first;    
				//save the matrix in case of eleSeed
				std::vector<DetId> rectangle =  EcalClusterTools::matrixDetId(topology, seedId, -9, 9, -9, 9);
				
				for(std::vector<DetId>::const_iterator itr = rectangle.begin(); itr != rectangle.end(); ++itr)
				{
					EcalRecHitCollection::const_iterator itrRecHit = theEndcapEcalRecHits->find(*itr) ;
					if(itrRecHit == theEndcapEcalRecHits->end()) continue;
					
					// fill recHit variables
					EEDetId endcapId(*itr);
					EcalScDetId scId = endcapId.sc();
					
					// laser correction
					theLaserCorrection = theLaser->getLaserCorrection(endcapId, iEvent.time());
					// IC correction
					EcalIntercalibConstantMap::const_iterator ICMapIt = ICMap.find(endcapId);
					theICCorrection = *ICMapIt;
					
					int SRFlag = 3;
					std::vector<EcalScDetId>::iterator SCIdListIt = std::find(SCIdList.begin(),SCIdList.end(),scId);
					if( SCIdListIt == SCIdList.end() ) SRFlag = 1;
					
					bool digiFound = false;
					for(EEDigiCollection::const_iterator digiItr = eeDigis->begin(); digiItr != eeDigis->end(); ++digiItr)
					{
						if(digiItr->id() != endcapId) continue;
						digiFound = true;
						EcalDataFrame df = *digiItr;
						for(int iSample = 0; iSample < 10; ++iSample)
							ele2_recHitMatrix_samples.push_back(df.sample(iSample).adc());
					}
					if( digiFound == false ) continue;
					
					ele2_recHitMatrix_E.push_back(itrRecHit->energy());
					ele2_recHitMatrix_flag.push_back(SRFlag*1000+itrRecHit->recoFlag());
					ele2_recHitMatrix_hashedIndex.push_back(endcapId.hashedIndex());
					ele2_recHitMatrix_ietaORix.push_back(endcapId.ix());
					ele2_recHitMatrix_iphiORiy.push_back(endcapId.iy());
					ele2_recHitMatrix_zside.push_back(0);
					ele2_recHitMatrix_laserCorrection.push_back(theLaserCorrection);
					ele2_recHitMatrix_ICConstant.push_back(theICCorrection);
					
				}
			}
		}	
		
		for(std::vector<std::pair<DetId,float> >::const_iterator rh = hits.begin(); rh!=hits.end(); ++rh)
		{
			float rhLaserCorrection = 1.;
			float rhICCorrection = 1.;
			float theLaserCorrection = -1.;
			float theICCorrection = -1.;
			float theAlpha = -1.;
            
			if ((*rh).first.subdetId()== EcalBarrel)
			{
				EBRecHitCollection::const_iterator itrechit = theBarrelEcalRecHits->find((*rh).first);
				if (itrechit==theBarrelEcalRecHits->end()) continue;
				EBDetId barrelId (itrechit->id ()); 
				++numRecHit;
				
				// laser correction
				theLaserCorrection = theLaser->getLaserCorrection(barrelId, iEvent.time());
				if ( applyCorrections_ ) rhLaserCorrection = theLaserCorrection;
				// IC correction
				EcalIntercalibConstantMap::const_iterator ICMapIt = ICMap.find(barrelId);
				theICCorrection = *ICMapIt;
				if ( applyCorrections_ ) rhICCorrection = theICCorrection;
				// Alpha
				EcalLaserAlphaMap::const_iterator italpha = theEcalLaserAlphaMap->find(itrechit->id());
				if( italpha != theEcalLaserAlphaMap->end() )
					theAlpha = (*italpha);
				sumRecHitE += itrechit->energy() * rhLaserCorrection * rhICCorrection ;
				sumLaserCorrectionRecHitE += theLaserCorrection * itrechit->energy() * rhLaserCorrection * rhICCorrection;
				// check if rh is inside the 5x5 matrix
				if ( fabs(barrelId.ieta() - ele2_seedIeta) < 3 && fabs(barrelId.iphi() - ele2_seedIphi) < 3 ) {
					sumRecHitE5x5 += itrechit->energy() * rhLaserCorrection * rhICCorrection ;
					sumLaserCorrectionRecHitE5x5 += theLaserCorrection * itrechit->energy() * rhLaserCorrection * rhICCorrection;
				}
				// check if rh is inside the 3x3 matrix
				if ( fabs(barrelId.ieta() - ele2_seedIeta) < 1 && fabs(barrelId.iphi() - ele2_seedIphi) < 1 ) {
					sumRecHitE3x3 += itrechit->energy() * rhLaserCorrection * rhICCorrection ;
					sumLaserCorrectionRecHitE3x3 += theLaserCorrection * itrechit->energy() * rhLaserCorrection * rhICCorrection ;
				}
				// fill recHit variables
				ele2_recHit_E.push_back(itrechit->energy() * rhLaserCorrection);
				ele2_recHit_flag.push_back(itrechit->recoFlag());
				ele2_recHit_hashedIndex.push_back(barrelId.hashedIndex());
				ele2_recHit_ietaORix.push_back(barrelId.ieta());
				ele2_recHit_iphiORiy.push_back(barrelId.iphi());
				ele2_recHit_zside.push_back(0);
				ele2_recHit_laserCorrection.push_back(theLaserCorrection);
				ele2_recHit_ICConstant.push_back(theICCorrection);
				ele2_recHit_Alpha.push_back(theAlpha);
				
				if( printOut && itrechit->energy() > 1. )
				{
					std::cout << std::fixed
					<< "    recHitLC: "    << std::setprecision(6) << std::setw(8) << theLaserCorrection
					<< "    recHitIC: "    << std::setprecision(6) << std::setw(8) << theICCorrection
					<< std::endl;
				}
			}
			
			if ((*rh).first.subdetId()== EcalEndcap)
			{
				EERecHitCollection::const_iterator itrechit = theEndcapEcalRecHits->find((*rh).first);
				if (itrechit==theEndcapEcalRecHits->end()) continue;
				EEDetId endcapId (itrechit->id ()); 
				++numRecHit;
				
				// laser correction
				theLaserCorrection = theLaser->getLaserCorrection(endcapId, iEvent.time());
				if ( applyCorrections_ ) rhLaserCorrection = theLaserCorrection;
				// IC correction
				EcalIntercalibConstantMap::const_iterator ICMapIt = ICMap.find(endcapId);
				theICCorrection = *ICMapIt;
				if ( applyCorrections_ ) rhICCorrection = theICCorrection;
				
				sumRecHitE += itrechit->energy() * rhLaserCorrection * rhICCorrection ;
				sumLaserCorrectionRecHitE += theLaserCorrection * itrechit->energy() * rhLaserCorrection * rhICCorrection;
				// check if rh is inside the 5x5 matrix
				if ( fabs(endcapId.ix() - ele2_seedIx) < 3 && fabs(endcapId.iy() - ele2_seedIy) < 3 ) {
					sumRecHitE5x5 += itrechit->energy() * rhLaserCorrection * rhICCorrection ;
					sumLaserCorrectionRecHitE5x5 += theLaserCorrection * itrechit->energy() * rhLaserCorrection * rhICCorrection;
				}
				// check if rh is inside the 3x3 matrix
				if ( fabs(endcapId.ix() - ele2_seedIx) < 1 && fabs(endcapId.iy() - ele2_seedIy) < 1 ) {
					sumRecHitE3x3 += itrechit->energy() * rhLaserCorrection * rhICCorrection ;
					sumLaserCorrectionRecHitE3x3 += theLaserCorrection * itrechit->energy() * rhLaserCorrection * rhICCorrection ;
				}
				// fill recHit variables
				ele2_recHit_E.push_back(itrechit->energy() * rhLaserCorrection);
				ele2_recHit_flag.push_back(itrechit->recoFlag());
				ele2_recHit_hashedIndex.push_back(endcapId.hashedIndex());
				ele2_recHit_ietaORix.push_back(endcapId.ix());
				ele2_recHit_iphiORiy.push_back(endcapId.iy());
				ele2_recHit_zside.push_back(endcapId.zside());
				ele2_recHit_laserCorrection.push_back(theLaserCorrection);
				ele2_recHit_ICConstant.push_back(theICCorrection);
				
				if( printOut && itrechit->energy() > 1. )
				{
					std::cout << std::fixed
					<< "    recHitLC: "    << std::setprecision(6) << std::setw(8) << theLaserCorrection
					<< "    recHitIC: "    << std::setprecision(6) << std::setw(8) << theICCorrection
					<< std::endl;
				}
			}
		}
		
		ele2_nRecHits = numRecHit;
		ele2_scLaserCorr = sumLaserCorrectionRecHitE/sumRecHitE;
		if ( applyCorrections_ ) ele2_scE = scRef->energy()*ele2_fEta*sumRecHitE;
		else ele2_scE = scRef->energy();
		
		ele2_fEtaCorr = fClusterCorrections(sumRecHitE+ele2_es,ele2_scEta,scRef->phiWidth()/scRef->etaWidth(),params)/fClusterCorrections(scRef->rawEnergy()+ele2_es,ele2_scEta,scRef->phiWidth()/scRef->etaWidth(),params);
        
		ele2_5x5LaserCorr = sumLaserCorrectionRecHitE5x5/sumRecHitE5x5;
		
		ele2_3x3LaserCorr = sumLaserCorrectionRecHitE3x3/sumRecHitE3x3;
		
		/// add regression input variables
		reco::SuperClusterRef s = electron.superCluster();
		reco::CaloClusterPtr b = s->seed(); //seed  basic cluster
		reco::CaloClusterPtr b2;
		reco::CaloClusterPtr bclast;
		reco::CaloClusterPtr bclast2;
		bool isbarrel =  b->hitsAndFractions().at(0).first.subdetId()==EcalBarrel;
		
		
		ele2_eRegrInput_rawE = s->rawEnergy();
		ele2_eRegrInput_r9 = lazyTools.e3x3(*b)/s->rawEnergy();
		ele2_eRegrInput_eta = s->eta();
		ele2_eRegrInput_phi = s->phi();
		ele2_eRegrInput_r25 = lazyTools.e5x5(*b)/s->rawEnergy();
		ele2_eRegrInput_hoe = electron.hcalOverEcal();
		ele2_eRegrInput_etaW = s->etaWidth() ;
		ele2_eRegrInput_phiW = s->phiWidth() ;
		ele2_eRegrInput_SCsize = s->clustersSize() ;
		ele2_eRegrInput_rho = *rhoHandle;
		ele2_eRegrInput_nPV = hVertexProduct->size();
		
		//seed basic cluster variables
		double bemax = lazyTools.eMax(*b);
		double be2nd = lazyTools.e2nd(*b);
		double betop = lazyTools.eTop(*b);
		double bebottom = lazyTools.eBottom(*b);
		double beleft = lazyTools.eLeft(*b);
		double beright = lazyTools.eRight(*b);
		
		double be2x5max = lazyTools.e2x5Max(*b);
		double be2x5top = lazyTools.e2x5Top(*b);
		double be2x5bottom = lazyTools.e2x5Bottom(*b);
		double be2x5left = lazyTools.e2x5Left(*b);
		double be2x5right = lazyTools.e2x5Right(*b);
		
		ele2_eRegrInput_Deta_bC_sC = b->eta()-s->eta();
		ele2_eRegrInput_Dphi_bC_sC = reco::deltaPhi(b->phi(),s->phi());
		ele2_eRegrInput_bCE_Over_sCE = b->energy()/b->energy();
		ele2_eRegrInput_e3x3_Over_bCE = lazyTools.e3x3(*b)/b->energy();
		ele2_eRegrInput_e5x5_Over_bCE = lazyTools.e5x5(*b)/b->energy();
		ele2_eRegrInput_sigietaieta_bC1 = sqrt(lazyTools.localCovariances(*b)[0]);
		ele2_eRegrInput_sigiphiiphi_bC1 = sqrt(lazyTools.localCovariances(*b)[2]);
		ele2_eRegrInput_sigietaiphi_bC1 = lazyTools.localCovariances(*b)[1];
		ele2_eRegrInput_bEMax_Over_bCE = bemax/b->energy();
		ele2_eRegrInput_bE2nd_Over_bCE = be2nd/b->energy();
		ele2_eRegrInput_bEtop_Over_bCE = betop/b->energy();
		ele2_eRegrInput_bEbot_Over_bCE = bebottom/b->energy();
		ele2_eRegrInput_bEleft_Over_bCE = beleft/b->energy();
		ele2_eRegrInput_bEright_Over_bCE = beright/b->energy();
		ele2_eRegrInput_be2x5max_Over_bCE = be2x5max/b->energy();
		ele2_eRegrInput_be2x5top_Over_bCE = be2x5top/b->energy();
		ele2_eRegrInput_be2x5bottom_Over_bCE = be2x5bottom/b->energy();
		ele2_eRegrInput_be2x5left_Over_bCE = be2x5left/b->energy();
		ele2_eRegrInput_be2x5right_Over_bCE = be2x5right/b->energy();
		
		
		
		if( isbarrel)
		{
			// seed cluster
			float betacry, bphicry, bthetatilt, bphitilt;
			int bieta, biphi;
			EcalClusterLocal _ecalLocal;
			_ecalLocal.localCoordsEB(*b,iSetup,betacry,bphicry,bieta,biphi,bthetatilt,bphitilt);
			
			ele2_eRegrInput_seedbC_eta = bieta;
			ele2_eRegrInput_seedbC_phi = biphi;
			ele2_eRegrInput_seedbC_eta_p5 = bieta%5;
			ele2_eRegrInput_seedbC_phi_p2 = biphi%2;
			ele2_eRegrInput_seedbC_bieta = (TMath::Abs(bieta)<=25)*(bieta%25) + (TMath::Abs(bieta)>25)*((bieta-25*TMath::Abs(bieta)/bieta)%20);
			ele2_eRegrInput_seedbC_phi_p20 = biphi%20;
			ele2_eRegrInput_seedbC_etacry = betacry;
			ele2_eRegrInput_seedbC_phicry = bphicry;
			ele2_eRegrInput_ESoSC = -99. ;
		}
		else
		{
			ele2_eRegrInput_ESoSC = s->preshowerEnergy()/s->rawEnergy();
			ele2_eRegrInput_seedbC_eta = -99.;
			ele2_eRegrInput_seedbC_phi = -99.;
			ele2_eRegrInput_seedbC_eta_p5 = -99.;
			ele2_eRegrInput_seedbC_phi_p2 = -99.;
			ele2_eRegrInput_seedbC_bieta = -99.;
			ele2_eRegrInput_seedbC_phi_p20 = -99.;
			ele2_eRegrInput_seedbC_etacry = -99.;
			ele2_eRegrInput_seedbC_phicry =-99.;
		}
		
		if(saveFbrem_){
			reco::GsfTrackRef eleTrack  = electron.gsfTrack();
			reco::TrackRef eleNTrack = electron.closestTrack();
			GlobalPoint outPos(eleTrack->extra()->outerPosition().x(), eleTrack->extra()->outerPosition().y(), eleTrack->extra()->outerPosition().z());
			GlobalPoint innPos(eleTrack->extra()->innerPosition().x(), eleTrack->extra()->innerPosition().y(), eleTrack->extra()->innerPosition().z());
			std::vector<reco::GsfTangent> eleTangent = eleTrack->gsfExtra()->tangents();
			int numberOfValidHits_Trk = 0;
			if(!eleNTrack.isNull())
				for(unsigned int trkNH = 0; trkNH < eleNTrack->extra()->recHitsSize(); ++trkNH){
					if(!eleNTrack->extra()->recHit(trkNH)->isValid()) continue;
					DetId id = eleNTrack->extra()->recHit(trkNH)->geographicalId();
					const GeomDet* det = pDD->idToDet(id);
					if(eleTrack->extra()->seedRef().isNull()) continue;
					edm::RefToBase<TrajectorySeed> seed = eleTrack->extra()->seedRef();
					ElectronSeedRef elseed=seed.castTo<ElectronSeedRef>();
					TrajectoryStateOnSurface t = trajectoryStateTransform::transientState(elseed->startingState(), &(det->surface()), &(*theMagField));
					if(!t.isValid()) continue;
					
					GlobalPoint hitPos = t.globalPosition();
					
					if(hitPos.x() == ele2_inner_x && hitPos.y() == ele2_inner_y && hitPos.z() == ele2_inner_z){
						ele2_inner_p = (sqrt(eleTrack->extra()->innerMomentum().Mag2()) );
						ele2_inner_x = innPos.x();
						ele2_inner_y = innPos.y();
						ele2_inner_z = innPos.z();
					}
					if(hitPos.x() == ele2_outer_x && hitPos.y() == ele2_outer_y && hitPos.z() == ele2_outer_z){
						ele2_outer_p = (sqrt(eleTrack->extra()->outerMomentum().Mag2()) );
						ele2_outer_x = outPos.x();
						ele2_outer_y = outPos.y();
						ele2_outer_z = outPos.z();
					}
					
					
					for(unsigned int pp=0; pp<eleTangent.size(); ++pp ){
						GlobalPoint tangPos( eleTangent.at(pp).position().x(),
											eleTangent.at(pp).position().y(),
											eleTangent.at(pp).position().z() );
						
						if(hitPos.x() != tangPos.x() && hitPos.y() != tangPos.y() && hitPos.z() != tangPos.z()) continue;
						float tangMom = sqrt(eleTangent.at(pp).momentum().Mag2());
						ele2_tangent_p.push_back(tangMom);
						ele2_tangent_x.push_back(tangPos.x());
						ele2_tangent_y.push_back(tangPos.y());
						ele2_tangent_z.push_back(tangPos.z());
						ele2_tangent_dP.push_back(eleTangent.at(pp).deltaP().value());
						ele2_tangent_dPerr.push_back(eleTangent.at(pp).deltaP().error());
						++numberOfValidHits_Trk;
					}
				}
			ele2_tangent_n = numberOfValidHits_Trk;
		}
	}
	
	
	if( verbosity_ )
		std::cout<< ">>> SimpleNtupleEoverP::fillEleInfo end <<<" << std::endl;
}

// -----------------------------------------------------------------------------------------

void SimpleNtupleEoverP::fillMetInfo(const edm::Event & iEvent, const edm::EventSetup & iSetup){
	//std::cout << "SimpleNtupleCalib::fillPFMetInfo" << std::endl;
	
	//*********** MET
	edm::Handle<edm::View<reco::MET> > PFmetHandle;
	iEvent.getByLabel(PFMetTag_,PFmetHandle);
	View<reco::MET>  mets = *PFmetHandle;
	reco::MET MET = mets.at(0);
	
	met = MET.p4();
	
	float cx0,cx1;
	float cy0,cy1; 
	
	if( dataFlag_ == true && (dataRun_=="2012A" || dataRun_=="2012B")){
		
		if( dataRun_ == "2012A" ){
			cx1 = +2.65299e-01; cx0 = +3.54233e-01;
			cy1 = -1.66425e-01; cy0 = +1.88923e-01;
		}
		if( dataRun_ == "2012B" ){
			cx1 = +2.65299e-01; cx0 = +3.54233e-01;
			cy1 = -1.66425e-01; cy0 = +1.88923e-01;
		}
		
		
		float metx = met.px();
		float mety = met.py();
		
		metx -= (cx0 + cx1*PV_n);
		mety -= (cy0 + cy1*PV_n);
		
		met.SetPxPyPzE(metx,mety,0,sqrt(metx*metx+mety*mety));
	}
	
	p_met = &met;
	met_et = p_met->Et();
	met_phi = p_met->phi();
	
	ele1Met_mt = sqrt( 2. * ele1_pt * met_et * ( 1. - cos( deltaPhi(ele1_phi,met_phi) ) ) );
	ele1Met_Dphi = deltaPhi(ele1_phi,met_phi);
	
	
}

//----------------------------------------------------------------------------------------------

void SimpleNtupleEoverP::fillDoubleEleInfo(const edm::Event & iEvent, const edm::EventSetup & iSetup){
	
	ele1ele2_m = (ele1 + ele2).mass();
	
	ROOT::Math::PtEtaPhiEVector ele1_sc(ele1_scE*sin(2*atan(exp(-1.*ele1_eta))),ele1_eta,ele1_phi,ele1_scE);
	ROOT::Math::PtEtaPhiEVector ele2_sc(ele2_scE*sin(2*atan(exp(-1.*ele2_eta))),ele2_eta,ele2_phi,ele2_scE);
	ele1ele2_scM = (ele1_sc + ele2_sc).mass();
	
	ROOT::Math::PtEtaPhiEVector ele1_sc_regression(ele1_scE_regression*sin(2*atan(exp(-1.*ele1_eta))),ele1_eta,ele1_phi,ele1_scE_regression);
	ROOT::Math::PtEtaPhiEVector ele2_sc_regression(ele2_scE_regression*sin(2*atan(exp(-1.*ele2_eta))),ele2_eta,ele2_phi,ele2_scE_regression);
	ele1ele2_scM_regression = (ele1_sc_regression + ele2_sc_regression).mass();
	
}


double SimpleNtupleEoverP::deltaPhi(const double& phi1, const double& phi2){ 
	double deltaphi = fabs(phi1 - phi2);
	if (deltaphi > 3.141592654) deltaphi = 6.283185308 - deltaphi;
	return deltaphi;
}

// -----------------------------------------------------------------------------------------
void SimpleNtupleEoverP::fillMCPUInfo (const edm::Event & iEvent, const edm::EventSetup & iSetup) 
{
	//std::cout << "SimpleNtupleCalib::fillMCPUInfo" << std::endl;
	
	edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
	iEvent.getByLabel(MCPileupTag_, PupInfo);
	
	
	// loop on BX
	// loop on BX
	std::vector<PileupSummaryInfo>::const_iterator PVI;
	for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
	{
		std::vector<float> temp_mc_PU_zpositions   = PVI->getPU_zpositions();
		std::vector<float> temp_mc_PU_sumpT_lowpT  = PVI->getPU_sumpT_lowpT();
		std::vector<float> temp_mc_PU_sumpT_highpT = PVI->getPU_sumpT_highpT();
		std::vector<int> temp_mc_PU_ntrks_lowpT    = PVI->getPU_ntrks_lowpT();
		std::vector<int> temp_mc_PU_ntrks_highpT   = PVI->getPU_ntrks_highpT();
		
		// in-time pileup
		if( PVI->getBunchCrossing() == 0 )
		{
			PUit_TrueNumInteractions = PVI->getTrueNumInteractions();
			PUit_NumInteractions = PVI->getPU_NumInteractions();    
			
			for(std::vector<float>::const_iterator it = temp_mc_PU_zpositions.begin(); it < temp_mc_PU_zpositions.end(); ++it)
				PUit_zpositions = *it;
			
			for(std::vector<float>::const_iterator it = temp_mc_PU_sumpT_lowpT.begin(); it < temp_mc_PU_sumpT_lowpT.end(); ++it)
				PUit_sumpT_lowpT=*it;
			
			for(std::vector<float>::const_iterator it = temp_mc_PU_sumpT_highpT.begin(); it < temp_mc_PU_sumpT_highpT.end(); ++it)
				PUit_sumpT_highpT = *it;
			
			for(std::vector<int>::const_iterator it = temp_mc_PU_ntrks_lowpT.begin(); it < temp_mc_PU_ntrks_lowpT.end(); ++it)
				PUit_ntrks_lowpT = *it;
			
			for(std::vector<int>::const_iterator it = temp_mc_PU_ntrks_highpT.begin(); it < temp_mc_PU_ntrks_highpT.end(); ++it)
				PUit_ntrks_highpT = *it;
		}
		
		// out-of-time pileup
		else
		{
			if (PVI->getBunchCrossing() < 0)
			{
				PUoot_early_TrueNumInteractions = PVI->getTrueNumInteractions();
				PUoot_early = PVI->getPU_NumInteractions();
				
				for(std::vector<float>::const_iterator it = temp_mc_PU_zpositions.begin(); it < temp_mc_PU_zpositions.end(); ++it)
					PUoot_early_zpositions = *it;
				
				for(std::vector<float>::const_iterator it = temp_mc_PU_sumpT_lowpT.begin(); it < temp_mc_PU_sumpT_lowpT.end(); ++it)
					PUoot_early_sumpT_lowpT = *it;
				
				for(std::vector<float>::const_iterator it = temp_mc_PU_sumpT_highpT.begin(); it < temp_mc_PU_sumpT_highpT.end(); ++it)
					PUoot_early_sumpT_highpT = *it;
				
				for(std::vector<int>::const_iterator it = temp_mc_PU_ntrks_lowpT.begin(); it < temp_mc_PU_ntrks_lowpT.end(); ++it)
					PUoot_early_ntrks_lowpT = *it;
				
				for(std::vector<int>::const_iterator it = temp_mc_PU_ntrks_highpT.begin(); it < temp_mc_PU_ntrks_highpT.end(); ++it)
					PUoot_early_ntrks_highpT = *it;
			}
			else
			{
				PUoot_late_TrueNumInteractions = PVI->getTrueNumInteractions();
				PUoot_late = PVI->getPU_NumInteractions();
				
				for(std::vector<float>::const_iterator it = temp_mc_PU_zpositions.begin(); it < temp_mc_PU_zpositions.end(); ++it)
					PUoot_late_zpositions = *it;
				
				for(std::vector<float>::const_iterator it = temp_mc_PU_sumpT_lowpT.begin(); it < temp_mc_PU_sumpT_lowpT.end(); ++it)
					PUoot_late_sumpT_lowpT = *it;
				
				for(std::vector<float>::const_iterator it = temp_mc_PU_sumpT_highpT.begin(); it < temp_mc_PU_sumpT_highpT.end(); ++it)
					PUoot_late_sumpT_highpT = *it;
				
				for(std::vector<int>::const_iterator it = temp_mc_PU_ntrks_lowpT.begin(); it < temp_mc_PU_ntrks_lowpT.end(); ++it)
					PUoot_late_ntrks_lowpT = *it;
				
				for(std::vector<int>::const_iterator it = temp_mc_PU_ntrks_highpT.begin(); it < temp_mc_PU_ntrks_highpT.end(); ++it)
					PUoot_late_ntrks_highpT = *it;          
			}
		}
	} // loop on BX
	
	
}


// -----------------------------------------------------------------------------------------
void SimpleNtupleEoverP::fillMCInfo (const edm::Event & iEvent, const edm::EventSetup & iESetup)
{
	//  std::cout << "SimpleNtupleEoverP::fillMCDecayInfo" << std::endl; 
	
	if(mcAnalysisZW_->isValid())
    {
		// Vboson Z 23  W +/-24
		mcV_E = mcAnalysisZW_->mcV()->p4().E();
		mcV_Px = mcAnalysisZW_->mcV()->p4().Px();
		mcV_Px = mcAnalysisZW_->mcV()->p4().Py();
		mcV_Px = mcAnalysisZW_->mcV()->p4().Pz();
		mcV_Charge = mcAnalysisZW_->mcV()->charge();
		mcV_PdgId = mcAnalysisZW_->mcV()->pdgId();
		
		// lepton1 pT > lepton2 pT
		//      mcAnalysisZW_->PrintEventInfo();
		mcF1_fromV_E = mcAnalysisZW_->mcF1_fromV()->p4().E();
		mcF1_fromV_Px = mcAnalysisZW_->mcF1_fromV()->p4().Px();
		mcF1_fromV_Py = mcAnalysisZW_->mcF1_fromV()->p4().Py();
		mcF1_fromV_Pz = mcAnalysisZW_->mcF1_fromV()->p4().Pz();
		mcF1_fromV_Charge = mcAnalysisZW_->mcF1_fromV()->charge();
		mcF1_fromV_PdgId = mcAnalysisZW_->mcF1_fromV()->pdgId();
		
		mcF2_fromV_E = mcAnalysisZW_->mcF2_fromV()->p4().E();
		mcF2_fromV_Px = mcAnalysisZW_->mcF2_fromV()->p4().Px();
		mcF2_fromV_Py = mcAnalysisZW_->mcF2_fromV()->p4().Py();
		mcF2_fromV_Pz = mcAnalysisZW_->mcF2_fromV()->p4().Pz();
		mcF2_fromV_Charge = mcAnalysisZW_->mcF2_fromV()->charge();
		mcF2_fromV_PdgId = mcAnalysisZW_->mcF2_fromV()->pdgId();
		
    }
	
}

//----------------------------------------------------------------------------------------------
DEFINE_FWK_MODULE(SimpleNtupleEoverP);
