#ifndef root_h
#define root_h

#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TString.h>

#include "TList.h"
#include "TKey.h"
#include "TObject.h"
#include "TCollection.h"

#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "THStack.h"


#include </usr/include/math.h>
#include <string>

#include "init.h"


using namespace std;

void setTDRStyle1(int gridFlagX, int gridFlagY, int logFlagX, int logFlagY);

const int nFigTypes = 3;
const string figType[nFigTypes] = {".eps",".png",".C"};
const int saveThis[nFigTypes] = {0,1,1};
const string plotsPrefix = "PLOTS/";
const string plotsPrefixC = "PLOTS/C/";
const string histoPrefix = "histo";


#endif // #ifdef root_cxx
