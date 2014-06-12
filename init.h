#ifndef init_h
#define init_h

#include <iostream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include </usr/include/math.h>
#include <ctime>

#include "TRandom3.h"
#include "TLorentzVector.h"

#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>

#include "analysis.h"
#include "utils.h"
#include "root.h"
#include "postInitialResultsAnalysis.h"

// ------------------------------------ NAMESPACES
using namespace std;

template <typename T> string tostr(const T& t) { ostringstream os; os<<t; return os.str(); }

const float pi=3.14159265359;
//
// ==================================================================================== main.cpp
//
const int doMainAnalysis             = 0;
const int doMainAnalysisPhoton       = 1;
const int doPostInitialResultsGraphs = 0;
//
// ==================================================================================== analysis.cc
//
const float correctionEnergy = 1.0;

const int DEBUG                             =0;

const int doSignalProfileFit                =0;

const int useOnlyDebugEventsMatchedToPerfect=1;
const int useDebugPerfectBecomeBadEvents    =1;

const int nDebugEvents=0;
const float pedestalDebugThreshold=190;
const float resolutionDebugThresholdMin=-0.1;
const float resolutionDebugThresholdMax=2.8;

const int nObj = 6;
const string objName[nObj] = {"GEN_Prompt","REC_Prompt","GEN_NotPrompt","REC_NotPrompt","GEN_mother_e","REC_not_matched"};
const string objNameLong[nObj] = {"GEN promt electron","REC prompt electron", "GEN non-promt electron","REC non-prompt electron","GEN electron as parent","REC not matched electron"};
const string objNameLong1[nObj] = {"GEN","REC","GEN","REC","GEN","REC"};
const string objNameLong2[nObj] = {"promt electron","prompt electron", "non-promt electron","non-prompt electron","parent e electron","not matched electron"};
const int objLineColor[nObj] = {1,2,1,2,1,2};
const int objLineStyle[nObj] = {1,1,1,1,1,1};
const int nLev = 4;
const string levName[nLev] = {"all","GenPreSel","andAllPromptGenMatched","andSoftIDReco"};
const float GenPtCut=10;
const float EtaMin=0;
const float EtaMax=2.5;
const float EtaGapMin=1.444, EtaGapMax=1.566;
const int takeOutEtaGap = 1;
const float cutDR = 0.03;//0.03;//0.2;//0.003;
const int nVar = 8;
const string varName[nVar] = {"Pt","Eta","Phi","E","PedestalBins","PedestalSlope","DefaultAveragePedestal","dR"};
const string varNameLong[nVar] = {"p_{T} [GeV]","#eta","#phi [rad]","Energy [GeV]",
	"0,1,2 time bins (pedestal) [ADC counts]",
	"pedestal slope [ADC counts]",
	"default average pedestal [ADC counts]",
	"min(#DeltaR(#eta,#phi)_{GEN,RECO})"
};
const int usualNBins = 100;
const int nBins[nVar] = {usualNBins,usualNBins,usualNBins,usualNBins,350,200,250,100};
const float minVar[nVar] = {0,-5.5,-1.0*pi,0,150,-10,150,0};
const float maxVar[nVar] = {50,5.5,1.0*pi,100,500,10,400,0.1};
const int varYFlagLog[nVar] = {0,0,0,0,1,1,1,1};
const float varMinValue[nVar] = {0,0,0,0.0,0.9,0.9,0.9,0.9};

const float adc_f_i_EB[10] = {0,0,0,0.0123648,0.761847,0.999833,0.886227,0.672032,0.473137,0.318696};//https://cmssdt.cern.ch/SDT/doxygen/CMSSW_3_9_7/doc/html/d8/dfe/EBShape_8cc-source.html
const float adc_f_i_EE[10] = {0,0,0,0.114335, 0.745677,0.986438,0.885063,0.677502,0.484776,0.339289};
//
// ==================================================================================== postInitialResultsAnalysis.cc
//
const int nFilesRes = 4;
const string inputResultsFiles[nFilesRes] = {
	/*"backup_15_44_56_27_09_2013_dE_over_E_resolution_promptMatchEffic_GenPtCut15GeV_0to0.5Eta_ALL/results.out",
	 "backup_15_55_54_27_09_2013_dE_over_E_resolution_promptMatchEffic_GenPtCut15GeV_0.5to1.479Eta_ALL/results.out",
	 "backup_16_02_31_27_09_2013_dE_over_E_resolution_promptMatchEffic_GenPtCut15GeV_1.479to3Eta_ALL/results.out"*/
	/*"backup_2013_10_03_14_09_40_dE_over_E__e5x5_as_E__resolution_promptMatchEffic_GenPtCut15GeV_0to0.5Eta_ALL/results.out",
	 "backup_2013_10_03_14_13_56_dE_over_E__e5x5_as_E__resolution_promptMatchEffic_GenPtCut15GeV_0.5to1.444Eta_ALL/results.out",
	 "backup_2013_10_03_14_22_16_dE_over_E__e5x5_as_E__resolution_promptMatchEffic_GenPtCut15GeV_1.566to2.5Eta_ALL/results.out"*/
	/*"backup_2013_10_04_16_18_03_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0to0.5Eta_ALL_v2/results.out",
	"backup_2013_10_04_16_22_08_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0.5to1.444Eta_ALL_v2/results.out",
	"backup_2013_10_04_16_30_34_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_1.566to2.5Eta_ALL_v2/results.out"*/
	/*"backup_2013_10_07_14_39_11_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0to0.5Eta_ALL_v3/results.out",
	"backup_2013_10_07_14_50_02_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0.5to1.444Eta_ALL_v3/results.out",
	"backup_2013_10_07_14_59_15_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_1.566to2.5Eta_ALL_v3/results.out"*/
	/*"backup_2013_10_07_20_53_47_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0to0.5Eta_ALL_v3_ReDoneInputs/results.out",
	"backup_2013_10_07_21_13_46_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0.5to1.444Eta_ALL_v3_ReDoneInputs/results.out",
	"backup_2013_10_07_21_23_33_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_1.566to2.5Eta_ALL_v3_ReDoneInputs/results.out"*/
	/*"backup_2013_10_08_00_22_03_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0to0.5Eta_ALL_v4_GausSigmaAdded/results.out",
	"backup_2013_10_08_00_35_06_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0.5to1.444Eta_ALL_v4_GausSigmaAdded/results.out",
	"backup_2013_10_08_00_51_56_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_1.566to2.5Eta_ALL_v4_GausSigmaAdded/results.out",
	"backup_2013_10_08_09_04_17_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0to2.5Eta_ALL_v4_GausSigmaAdded/results.out"*/
	/*"backup_2013_10_08_11_22_04_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0to0.5Eta_ALL_v5_HZZ2e_HZZ4e_eff_added/results.out",
	"backup_2013_10_08_13_01_29_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0.5to1.444Eta_ALL_v5_HZZ2e_HZZ4e_eff_added/results.out",
	"backup_2013_10_08_13_11_55_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_1.566to2.5Eta_ALL_v5_HZZ2e_HZZ4e_eff_added/results.out",
	"backup_2013_10_08_11_10_07_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0to2.5Eta_ALL_v5_HZZ2e_HZZ4e_eff_added/results.out"*/
	/*"backup_2013_10_08_14_25_37_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0to0.5Eta_ALL_v6_correctedEnergy/results.out",
	"backup_2013_10_08_14_42_24_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0.5to1.444Eta_ALL_v6_correctedEnergy/results.out",
	"backup_2013_10_08_14_56_05_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_1.566to2.5Eta_ALL_v6_correctedEnergy/results.out",
	"backup_2013_10_08_15_04_54_ratioE__e5x5_plus_es_as_RecoE__moreParamsOutputGenPt10_0to2.5Eta_ALL_v6_correctedEnergy/results.out"*/
	/*"backup_2013_10_23_17_37_54_allPU_allAge_fakeRatesInputs_0to0.5eta_OK/results.out",
	"backup_2013_10_23_17_51_40_allPU_allAge_fakeRatesInputs_0.5to1.444eta_OK/results.out",
	"backup_2013_10_23_18_08_39_allPU_allAge_fakeRatesInputs_1.566to2.5eta_OK/results.out",
	"backup_2013_10_23_17_22_52_allPU_allAge_fakeRatesInputs_0to2.5eta_OK/results.out"*/
    "backup_2013_11_11_14_14_31______MODIFIED_PATH_TO_ROOT____PF_electrons_test____allPU_allAge_fakeRatesInputs_0to0.5eta_OK/results.out",
    "backup_2013_11_11_14_19_37______MODIFIED_PATH_TO_ROOT____PF_electrons_test____allPU_allAge_fakeRatesInputs_0.5to1.444eta_OK/results.out",
    "backup_2013_11_11_14_24_35______MODIFIED_PATH_TO_ROOT____PF_electrons_test____allPU_allAge_fakeRatesInputs_1.566to2.5eta_OK/results.out",
    "backup_2013_11_11_14_29_21______MODIFIED_PATH_TO_ROOT____PF_electrons_test____allPU_allAge_fakeRatesInputs_0to2.5eta_OK/results.out"
};
const string resultsRootFiles = "resultsRootFiles.txt";
const int nVarRes = 7+6;
const int nEta = 4;
const int nPU = 3;
const int nLumi = 6;
//const string nameVar[nVarRes] = {"ERatioMean","ERatioSigma","Efficiency"};
//const string nameVarLong[nVarRes] = {"mean of (E_{REC}-E_{GEN})/E_{GEN} [GeV]","#sigma of (E_{REC}-E_{GEN})/E_{GEN} [GeV]","efficiency"};
const string nameVar[nVarRes] = {"ERatioMean","ERatioSigma","ERatioFWHM","Efficiency","Efficiency10","Efficiency20","ERatioSigmaGaus",
	"EffHZZ2e2mu","EffHZZ2e2mu10","EffHZZ2e2mu20",
	"EffHZZ4e","EffHZZ4e10","EffHZZ4e20"
};
const string nameVarLong[nVarRes] = {
	"mean of E_{REC}/E_{GEN} [GeV]",
	"#sigma of E_{REC}/E_{GEN} [GeV]",
	"FWHM of E_{REC}/E_{GEN} [GeV]",
	"efficiency",
	"efficiency 10%",
	"efficiency 20%",
	"Gaus fit #sigma of E_{REC}/E_{GEN} peak [GeV]",
	"HZZ2e2mu efficiency",
	"HZZ2e2mu efficiency 10%",
	"HZZ2e2mu efficiency 20%",
	"HZZ4e efficiency",
	"HZZ4e efficiency 10%",
	"HZZ4e efficiency 20%"	
};
const string nameEta[nEta] = {"eta0to0.5","eta0.5to1.444","eta1.566to2.5","eta0to2.5"};
const string namePU[nPU] = {"PU0","PU70","PU140"};
const string nameLumi[nLumi] = {"0fb-1","150fb-1","300fb-1","500fb-1","1Kfb-1","3Kfb-1"};
const int colorsLumi[nLumi] = {1,2,4,6,8,49};
const int colors[nPU] = {1,2,4};
const float rangeMinVar[nVarRes] = {0.8,0.1,0.0,0.35,0.,0.,0.0, 0.,0.,0.,0.,0.,0.};
const float rangeMaxVar[nVarRes] = {1.0,0.3,0.3,1.00,1.0,1.0,0.15, 1.0,1.0,1.0,1.0,1.0,1.0};



#endif // #ifdef init_cxx
