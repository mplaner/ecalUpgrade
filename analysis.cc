#define analysis_cxx
#include "init.h"
#include <TString.h>

using namespace std;

const int debug=0;

void analysis::LoopEle()
{
    int t = system("ls -l input.root | awk '{print $11}' | sed 's%\/% %g' | sed 's%.root%%' | awk '{print $4}' | sed 's%\_% %g' | awk '{print $11}' > PUname.txt");
    ifstream myfilePUname;
	myfilePUname.open("PUname.txt",std::ifstream::in);
    string PUmode;
    myfilePUname >> PUmode;
    cout << "PU mode is: " << PUmode << endl;
    myfilePUname.close();
    const int nPUoptions = 3;
    const string PUoptions[nPUoptions] = {"NoPileUp","PU70bx25","PU140Bx25"};
    const float cutR9forPUoptions[nPUoptions] = {0.94,0.94,0.94};//{0.8,0.8,0.8};//0.92,0.85,0.8};
    const float dRcut=0.15;//dR cut: <0.15, see comment above and: backup_2013_12_05_14_47_54_PHOTON__dR_cut_definition
    const float cutFBrem=0.2;
    const float barrelEtaMax = 1.444;
    const float endcapEtaMin = 1.566;
    const float endcapEtaMax = 3.0;
    const float genPtMin = 30;
    float tmpR9cut=-0.8;
    for (int i=0; i<nPUoptions; i++) {
        if (PUmode == PUoptions[i]) tmpR9cut=cutR9forPUoptions[i];
    }
    const float cutR9=tmpR9cut;
    cout << "PU mode is: " << PUmode << ". R9 cut to be used: " << cutR9 << endl;
    const int nEta = 6;
    const float etaMin[nEta] = {0,0.5,1.0,1.566,2.0,2.5};
    const float etaMax[nEta] = {0.5,1.0,1.444,2.0,2.5,3.0};
    //const int etaColor[nEta] = {1,2,4,6,8,9};//{20,21,22,23,24};
    
    vector<float> nGenEle[nEta], nRecEle[nEta], resolEle[nEta];
    
    TH1F *h_recFBrem;
	h_recFBrem = new TH1F("h_recFBrem","",100,0.,1.);
    
    if (fChain == 0) return; Long64_t nentries = fChain->GetEntriesFast(); Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        //for (Long64_t jentry=0; jentry<10;jentry++) {
        Long64_t ientry = LoadTree(jentry); if (ientry < 0) break; nb = fChain->GetEntry(jentry);   nbytes += nb; // if (Cut(ientry) < 0) continue;
        
        //const int cNELE=NELE, cNGEN=NGEN;
        
        TLorentzVector r[NELE_max], g[NGEN_max];
        for (int i=0; i<NGEN; i++) {
            float x,y,z;
            x = GEN_px[i];
            y = GEN_py[i];
            z = GEN_pz[i];
            g[i].SetXYZT(x,y,z,sqrt(x*x+y*y+z*z));
        }
        for (int j=0; j<NELE; j++) {
            float eta,phi,pt,e;
            pt  = ELE_pt[j];
            eta = ELE_eta[j];
            phi = ELE_phi[j];
            e   = ELE_e[j];
            r[j].SetPtEtaPhiE(pt,eta,phi,e);
        }
        
        int gen_taken[NGEN]; for (int i=0; i<NGEN; i++) gen_taken[i]=0;
        int rec_taken[NELE]; for (int j=0; j<NELE; j++) rec_taken[j]=0;
        for (int i=0; i<NGEN; i++) {
            if (abs(GEN_id[i])==11 && GEN_status[i]==1 && (GEN_parent[i]==25 || GEN_parent[i]==23)
                && fabs(g[i].Eta())<endcapEtaMax
                && g[i].Eta()>barrelEtaMax*(-1.0) // TEMP special for 2013 samples bug fix: no EE-
                && (fabs(g[i].Eta())<barrelEtaMax || fabs(g[i].Eta())>endcapEtaMin)
                && g[i].Pt()>genPtMin && !gen_taken[i]) {
                float etaBin=g[i].Eta();
                for (int k=0; k<nEta; k++) {
                    if (fabs(etaBin) >= etaMin[k] && fabs(etaBin) < etaMax[k] ) {
                        nGenEle[k].push_back(etaBin);
                    }
                }
                float mindr=100000;
                int selected[NELE]; for (int j=0; j<NELE; j++) selected[j]=0;
                int nselected=0;
                int selected_j=-1;
                for (int j=0; j<NELE; j++) {
                    float tmpdr = r[j].DeltaR(g[i]);
                    if (!rec_taken[j] && tmpdr<dRcut && tmpdr<mindr) {
                        selected[j]=1;
                        nselected++;
                    }
                }
                if (nselected>1) {
                    float minde=100000;
                    for (int j=0; j<NELE; j++) {
                        if (selected[j] && fabs(g[i].E()-r[j].E())<minde) {
                            minde=fabs(g[i].E()-r[j].E());
                            selected_j=j;
                        }
                    }
                } else {
                    for (int j=0; j<NELE; j++) {
                        if (selected[j]) selected_j=j;
                    }
                }
                if (selected_j>=0) {
                    int j=selected_j;
                    gen_taken[i]=1;
                    rec_taken[j]=1;
                    float etaBin=r[j].Eta();
                    for (int k=0; k<nEta; k++) {
                        if (fabs(etaBin) >= etaMin[k] && fabs(etaBin) < etaMax[k] ) {
                            h_recFBrem->Fill(ELE_fbrem[j]);
                            nRecEle[k].push_back(etaBin);
                            if (ELE_fbrem[j]<cutFBrem) {
                                resolEle[k].push_back(ELE_e5x5[j]/g[i].E());
                            }
                        }
                    }
                }
            }
        }
    }
    
    if (1) {
        int logOY=0;
        setTDRStyle1(1,1,0,logOY);//gridX/Y, logX/Y
        gStyle->TStyle::SetOptStat(0);
        gStyle->TStyle::SetOptFit(0);
        TCanvas *c1 = new TCanvas("c1","myPlots",0,0,600,600);
        c1->cd(1);
        
        {
            c1->SetLogy(1);
			TString name1 = "FBrem for REC electrons";
			const TString xTitle = name1;
			h_recFBrem->GetXaxis()->SetTitle(xTitle);
			const TString yTitle = "Number of entries";
			h_recFBrem->GetYaxis()->SetTitle(yTitle);
			h_recFBrem->GetYaxis()->SetTitleOffset(2.0);
			h_recFBrem->GetXaxis()->SetNdivisions(505);
			h_recFBrem->GetYaxis()->SetNdivisions(505);
			h_recFBrem->SetLineColor(1);
			h_recFBrem->SetLineStyle(1);
			h_recFBrem->SetLineWidth(2);
			h_recFBrem->SetMarkerStyle(21);
			h_recFBrem->SetMarkerColor(1);
			h_recFBrem->SetMarkerSize(0.1);
			h_recFBrem->SetFillColor(1);
			h_recFBrem->SetFillStyle(0);
			//h_recFBrem->SetMinimum(0);
			h_recFBrem->Draw();
			for (int ii=0; ii<nFigTypes; ii++) {
				if (saveThis[ii]) {
					TString PP = plotsPrefix;
					if (ii==2) PP = plotsPrefixC;
					TString name1;
					name1 = PP + histoPrefix
					+ "_h_recFBrem"
					+ figType[ii];
					TString ss = name1;
					c1->SaveAs(ss);
				}
			}
			c1->Clear();
            c1->SetLogy(0);
		}

        
    }
    
    ofstream myfile;
	myfile.open("results.out2",std::ofstream::app);
    
    for (int j=0; j<nEta; j++) {
        myfile << "resultsForElectrons " << fabs(etaMin[j]-etaMax[j])/2.0 + etaMin[j];
        float ratio=0;
        if (nGenEle[j].size()) ratio=float(nRecEle[j].size())/float(nGenEle[j].size());
        myfile << " " << nGenEle[j].size()
        << " " << nRecEle[j].size()
        << " " << ratio
        << " " << defineSmallestIntervalToContain(0.68,resolEle[j])/2.0//resolEle[j]
        << endl;
    }
    
    myfile.close();
    
}//end: analysis::LoopEle()

void analysis::LoopPho()
{
    /*
    const int testN=7;
    float a1[testN] = {23,45,12,78,78,9,64};
    float a2[testN] = {0.23,0.45,0.12,0.78,0.39,0.9,0.64};
    vector<float> v1,v2;
    for (int i=0; i<testN; i++) {
        v1.push_back(a1[i]);
        v2.push_back(a2[i]);
    }
    for (int i=0; i<testN; i++) cout << v1[i] << " ";
    cout << endl;
    for (int i=0; i<testN; i++) cout << v2[i] << " ";
    cout << endl;
    sortKeepAligned(v1,v2);
    for (int i=0; i<testN; i++) cout << v1[i] << " ";
    cout << endl;
    for (int i=0; i<testN; i++) cout << v2[i] << " ";
    cout << endl;
    */
    //exit(1);

    
    int t = system("ls -l input.root | awk '{print $11}' | sed 's%\/% %g' | sed 's%.root%%' | awk '{print $4}' | sed 's%\_% %g' | awk '{print $11}' > PUname.txt");
    ifstream myfilePUname;
	myfilePUname.open("PUname.txt",std::ifstream::in);
    string PUmode;
    myfilePUname >> PUmode;
    cout << "PU mode is: " << PUmode << endl;
    myfilePUname.close();
    const int nPUoptions = 3;
    const string PUoptions[nPUoptions] = {"NoPileUp","PU70bx25","PU140Bx25"};
    const float cutR9forPUoptions[nPUoptions] = {0.94,0.94,0.94};//{0.8,0.8,0.8};//0.92,0.85,0.8};
    /*
     NOTES:
     dR cut: relaxed to 0.15 to have ~100% efficiency to find some cluster close to the original GEN photon, but then will apply inv. mass constrain
     
    */
    const float dRcut=0.15;//dR cut: <0.15, see comment above and: backup_2013_12_05_14_47_54_PHOTON__dR_cut_definition
    const float barrelEtaMax = 1.444;
    const float endcapEtaMin = 1.566;
    const float endcapEtaMax = 3.0;
    const float genPtMin = 10;
    float tmpR9cut=-0.8;
    for (int i=0; i<nPUoptions; i++) {
        if (PUmode == PUoptions[i]) tmpR9cut=cutR9forPUoptions[i];
    }
    const float cutR9=tmpR9cut;
    cout << "PU mode is: " << PUmode << ". R9 cut to be used: " << cutR9 << endl;
    
    /*
    const int nEta = 13;//+2;
    const float etaMin[nEta] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.444,1.566,1.7,1.9,2.1,2.3};//,0,1.566};
    const float etaMax[nEta] = {0.2,0.4,0.6,0.8,1.0,1.2,1.444,1.566,1.7,1.9,2.1,2.3,2.5};//,1.444,2.5};
    const int etaColor[nEta] = {20,21,22,23,24,25,26,27,28,41,42,43,44};
     */
    const int nEta = 6;
    const float etaMin[nEta] = {0,0.5,1.0,1.566,2.0,2.5};
    const float etaMax[nEta] = {0.5,1.0,1.444,2.0,2.5,3.0};
    const int etaColor[nEta] = {1,2,4,6,8,9};//{20,21,22,23,24};
    const int nVar=4+6+3+1+3+1;//=18
    const string varName[nVar] = {"ET","Esc","e5x5","e3x3","genEta","genPhi","genEt","recEta","recPhi","r9","resolEsc","resolE5x5","resolE3x3","sigmaIEtaIEta","recFakesEsc","recFakesEt","recFakesEta","nCrystals"};
    const string varNameXaxis[nVar] = {"REC #gamma E_{T} [GeV]","REC #gamma E_{SC} [GeV]","REC #gamma energy in 5x5 [GeV]","REC #gamma energy in 3x3 [GeV]",
        "GEN #gamma #eta","GEN #gamma #phi [rad]","GEN #gamma E_{T} [GeV]","REC #gamma #eta","REC #gamma #phi [rad]",
        "REC r9","E_{SC,REC}/E_{GEN}","E5x5_{REC}/E_{GEN}","E3x3_{REC}/E_{GEN}","#sigma(i#eta,i#eta)",
        "fake: REC #gamma E_{SC} [GeV]","fake: REC #gamma E_{T} [GeV]","fake: REC #gamma #eta","nCrystals"};
    const int stdNBins=200;
    const int nBins[nVar] = {stdNBins,stdNBins,stdNBins,stdNBins,62,20,stdNBins,62,20,stdNBins,stdNBins,stdNBins,stdNBins,stdNBins,stdNBins,stdNBins,stdNBins,stdNBins};
    const float minVar[nVar] = {0,0,0,0,        -3.1,-pi,0,   -3.1,-pi,0  ,0,0,0,0,  0,0,-3.1,    0};
    const float maxVar[nVar] = {200,200,250,200, 3.1, pi,200,  3.1, pi,1.1,2,2,2,0.1,200,200,3.1, 1000};
    const float minVal[nVar] = {-0,-0,-0,-0,-0,0,-0,-0,0,-0,-0,-0,-0,-0,-0,-0,-0,-0};
    const int setMinVal[nVar] = {0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0};
    TH1F *histo1d[nVar][nEta];
    vector<float> v_histo1d[nVar][nEta];
    for (int i=0; i<nVar; i++) {
		for (int j=0; j<nEta; j++) {
            string etaName = "eta_" + tostr(etaMin[j]) + "_to_" + tostr(etaMax[j]);
            TString name1;
            name1 = histoPrefix
            + "_PHOTON_" + varName[i]
            + "_" + etaName;
            const TString name = name1;
            const TString title = name;
            const TString xTitle = varNameXaxis[i];//name1;
            //const int NB = nBins[k];
            float maxVarTmp = maxVar[i];
            if (i==2) maxVarTmp = (etaMax[j]-0.5)*maxVar[i] + maxVar[i];
            histo1d[i][j] = new TH1F(name, title, nBins[i], minVar[i], maxVarTmp);
            histo1d[i][j]->GetXaxis()->SetTitle(xTitle);
            const TString yTitle = "Number of entries";
            histo1d[i][j]->GetYaxis()->SetTitle(yTitle);
            histo1d[i][j]->GetYaxis()->SetTitleOffset(2.0);
            histo1d[i][j]->GetXaxis()->SetNdivisions(505);
            histo1d[i][j]->GetYaxis()->SetNdivisions(505);
            histo1d[i][j]->SetLineColor(etaColor[j]);
            histo1d[i][j]->SetLineStyle(1);
            histo1d[i][j]->SetLineWidth(2);
            histo1d[i][j]->SetMarkerStyle(21);
            histo1d[i][j]->SetMarkerColor(1);
            histo1d[i][j]->SetMarkerSize(0.1);
            histo1d[i][j]->SetFillColor(1);
            histo1d[i][j]->SetFillStyle(0);
            if (setMinVal[i]) histo1d[i][j]->SetMinimum(minVal[i]);
        }
    }
    TH1F *GEN_REC_match_dr;
	GEN_REC_match_dr = new TH1F("GEN_REC_match_dr","",100,0.,10.);
    TH1F *GEN_REC_match_dr_zoom;
	GEN_REC_match_dr_zoom = new TH1F("GEN_REC_match_dr_zoom","",100,0.,0.3);
    TH1F *nREC_in_dR_hist;
	nREC_in_dR_hist = new TH1F("nREC_in_dR_hist","",10,0.,10.);
    
    TH1F *histo1d_tmp_pR9 = new TH1F("tmp_pR9", "tmp_pR9", 150, 0., 1.5);
    
    const int nEtaPairs = 3;
    const float etaMinPairs[nEtaPairs] = {0,0,1.566};
    const float etaMaxPairs[nEtaPairs] = {3.0,1.444,3.0};
    const int etaColorPairs[nEtaPairs] = {1,2,4};//{20,21,22};
    const int nVarPairs=4;
    const string varNamePairs[nVarPairs] = {"mGen","mRecSc","mRecE5x5","mRecE3x3"};
    //const int stdNBins=200;
    const int nBinsPairs[nVarPairs] = {stdNBins,stdNBins,stdNBins,stdNBins};
    const float minVarPairs[nVarPairs] = {0,0,0,0};
    const float maxVarPairs[nVarPairs] = {200,200,200,200};
    TH1F *histo1dPairs[nVarPairs][nEtaPairs];
    for (int i=0; i<nVarPairs; i++) {
		for (int j=0; j<nEtaPairs; j++) {
            string etaName = "eta_" + tostr(etaMinPairs[j]) + "_to_" + tostr(etaMaxPairs[j]);
            TString name1;
            name1 = histoPrefix
            + "_PHOTON_" + varNamePairs[i]
            + "_" + etaName;
            const TString name = name1;
            const TString title = name;
            const TString xTitle = name1;
            //const int NB = nBins[k];
            histo1dPairs[i][j] = new TH1F(name, title, nBinsPairs[i], minVarPairs[i], maxVarPairs[i]);
            histo1dPairs[i][j]->GetXaxis()->SetTitle(xTitle);
            const TString yTitle = "Number of entries";
            histo1dPairs[i][j]->GetYaxis()->SetTitle(yTitle);
            histo1dPairs[i][j]->GetYaxis()->SetTitleOffset(2.0);
            histo1dPairs[i][j]->GetXaxis()->SetNdivisions(505);
            histo1dPairs[i][j]->GetYaxis()->SetNdivisions(505);
            histo1dPairs[i][j]->SetLineColor(etaColorPairs[j]);
            histo1dPairs[i][j]->SetLineStyle(1);
            histo1dPairs[i][j]->SetLineWidth(2);
            histo1dPairs[i][j]->SetMarkerStyle(21);
            histo1dPairs[i][j]->SetMarkerColor(1);
            histo1dPairs[i][j]->SetMarkerSize(0.1);
            histo1dPairs[i][j]->SetFillColor(1);
            histo1dPairs[i][j]->SetFillStyle(0);
            //histo1dPairs[i][j]->SetMinimum(...);
        }
    }
    
    
    vector<float> v_m1,v_m2,v_m3;
    
    if (fChain == 0) return; Long64_t nentries = fChain->GetEntriesFast(); Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //for (Long64_t jentry=0; jentry<10;jentry++) {
        Long64_t ientry = LoadTree(jentry); if (ientry < 0) break; nb = fChain->GetEntry(jentry);   nbytes += nb; // if (Cut(ientry) < 0) continue;
        
        // re-define r9:
        for (int i=0;i<psize;i++) {
            pR9[i]=p3x3_energy[i]/p5x5_energy[i];
        }
            
        
        int rec_pair_found=0;
        int REC_in_dR[psize];
        int REC_prompt[psize];
        for (int i=0;i<psize;i++) {
            REC_in_dR[i]=0;
            REC_prompt[i]=0;
        }
        
        int GEN_promptPhoton[NGEN];
        int GEN_matched[NGEN];
        for (int i=0; i<NGEN; i++) {
            GEN_promptPhoton[i]=0;
            TLorentzVector g;
            float x,y,z;
            x = GEN_px[i];
            y = GEN_py[i];
            z = GEN_pz[i];
            g.SetXYZT(x,y,z,sqrt(x*x+y*y+z*z));
			if (abs(GEN_id[i])==22 && GEN_status[i]==1 && GEN_parent[i]==25
                && fabs(g.Eta())<endcapEtaMax
                && g.Eta() > (-1.0)*barrelEtaMax // TEMP to take out EE- because of the bug in 2013 samples production for aging model
                && (fabs(g.Eta())<barrelEtaMax || fabs(g.Eta())>endcapEtaMin)
                && g.Pt()>genPtMin) {
                GEN_promptPhoton[i] = 1;
            }
        }
        int nGenPrompt=0;
        int nGenPromptBar=0;
        int nGenPromptEnd=0;
        for (int i=0; i<NGEN; i++) {
            if (GEN_promptPhoton[i]) nGenPrompt++;
        }
        if (nGenPrompt>2) cout << "WARNING: nGenPrompt>2, it's: " << nGenPrompt << endl; // checked on 05/12/13: no such WARNINGS, even without pt and eta cuts

        //first photon
        for (int i1=0; i1<NGEN; i1++) {
            GEN_matched[i1]=0;
            if (GEN_promptPhoton[i1]) {
                float dRmin=10000;
                int ind=-1;
                TLorentzVector g;
                float x,y,z;
                x = GEN_px[i1];
                y = GEN_py[i1];
                z = GEN_pz[i1];
                g.SetXYZT(x,y,z,sqrt(x*x+y*y+z*z));
                for (int i=0;i<psize;i++) {
                    TLorentzVector r;
                    float eta,phi,pt,e;
                    pt  = pPt[i];
                    eta = pEta[i];
                    phi = pPhi[i];
                    e   = pSC_energy[i];
                    r.SetPtEtaPhiE(pt,eta,phi,e);
                    float tmpdr = r.DeltaR(g);
                    if (tmpdr<dRcut) {
                        REC_in_dR[i] = 1;
                    }
                    if (tmpdr<dRmin) {
                        dRmin = tmpdr;
                        ind=i;
                    }
                }
                GEN_REC_match_dr->Fill(dRmin);
                //GEN_REC_match_dr_zoom->Fill(dRmin);//FIX: will fill in below, only for selected pair...
                if (ind>=0) {
                    if (dRmin<dRcut) {
                        GEN_matched[i1]=1;
                        REC_prompt[ind]=1;
                    }
                }
            }
        }
        //second photon:
        for (int i1=0; i1<NGEN; i1++) {
            if (GEN_promptPhoton[i1] && !GEN_matched[i1]) {
                float dRmin=10000;
                int ind=-1;
                TLorentzVector g;
                float x,y,z;
                x = GEN_px[i1];
                y = GEN_py[i1];
                z = GEN_pz[i1];
                g.SetXYZT(x,y,z,sqrt(x*x+y*y+z*z));
                for (int i=0;i<psize;i++) {
                    if (!REC_prompt[i]) {
                        TLorentzVector r;
                        float eta,phi,pt,e;
                        pt  = pPt[i];
                        eta = pEta[i];
                        phi = pPhi[i];
                        e   = pSC_energy[i];
                        r.SetPtEtaPhiE(pt,eta,phi,e);
                        float tmpdr = r.DeltaR(g);
                        if (tmpdr<dRcut) {
                            REC_in_dR[i] = 1;
                        }
                        if (tmpdr<dRmin) {
                            dRmin = tmpdr;
                            ind=i;
                        }
                    }
                    GEN_REC_match_dr->Fill(dRmin);
                    //GEN_REC_match_dr_zoom->Fill(dRmin);//FIX: will fill in below, only for selected pair...
                }
                if (ind>=0) {
                    if (dRmin<dRcut) {
                        GEN_matched[i1]=1;
                        REC_prompt[ind]=1;
                    }
                }
            }
        }
        
        if (nGenPrompt==2) {
            int count=0;
            for (int i=0;i<psize;i++) {
                if (REC_in_dR[i]) count++;
            }
            nREC_in_dR_hist->Fill(count);
        }
        
        int gind1=-1,gind2=-1;
        int rind1=-1,rind2=-1;
        TLorentzVector g1,g2,r1,r2;
        
        for (int i1=0; i1<NGEN; i1++) {
            if (GEN_promptPhoton[i1]) {
                TLorentzVector g;
                float x,y,z;
                x = GEN_px[i1];
                y = GEN_py[i1];
                z = GEN_pz[i1];
                g.SetXYZT(x,y,z,sqrt(x*x+y*y+z*z));
                if (gind1==-1) {
                    g1=g;
                    gind1=i1;
                }
                if (gind1!=-1 && gind2==-1 && i1!=gind1) {
                    g2=g;
                    gind2=i1;
                }
                float etaBin=g.Eta();
                for (int j=0; j<nEta; j++) {
                    histo1d[4][j]->Fill(g.Eta(),1.0/float(nEta));
                    if (fabs(etaBin) >= etaMin[j] && fabs(etaBin) < etaMax[j] ) {
                        histo1d[5][j]->Fill(g.Phi());
                        histo1d[6][j]->Fill(g.Pt());
                        v_histo1d[6][j].push_back(g.Pt());
                    }
                }
            }
        }
        
        /*
        float tmpdm=10000000;
        float truem=invMass(g1.E(),g1.Eta(),g1.Phi(),g2.E(),g2.Eta(),g2.Phi());
        for (int i=0;i<psize;i++) {
            if (REC_prompt[i]) {
                for (int ii=0;ii<psize;ii++) {
                    if (REC_prompt[ii]) {
                        if (i!=ii) {
                            TLorentzVector r11;
                            float eta,phi,pt,e;
                            pt  = pPt[i];
                            eta = pEta[i];
                            phi = pPhi[i];
                            e   = pSC_energy[i];
                            r11.SetPtEtaPhiE(pt,eta,phi,e);
                            TLorentzVector r22;
                            //float eta,phi,pt,e;
                            pt  = pPt[ii];
                            eta = pEta[ii];
                            phi = pPhi[ii];
                            e   = pSC_energy[ii];
                            r22.SetPtEtaPhiE(pt,eta,phi,e);
                            float m1=invMass(r11.E(),r11.Eta(),r11.Phi(),r22.E(),r22.Eta(),r22.Phi());
                            if (fabs(m1-truem)<tmpdm) {
                                rind1=i;
                                rind2=ii;
                                r1=r11;
                                r2=r22;
                            }
                        }
                    }
                }
            }
        }
        */
        
        // actually, there is more "per photon" flexibility if full match and unambiguous selection is done on per-photon basis (otherwise efficiency of individual photon reco is artif. lowered)
        // the algo below also constructed in such a way that gind1 corresponds to rind1 and gind2 to rind2 (for resolution, etc.). Selection by energy is based on E_SC.
        for (int i=0;i<psize;i++) {
            if (REC_prompt[i] && pR9[i]>cutR9) {
                float tmpde=100000000;
                TLorentzVector r;
                float eta,phi,pt,e;
                pt  = pPt[i];
                eta = pEta[i];
                phi = pPhi[i];
                e   = pSC_energy[i];
                r.SetPtEtaPhiE(pt,eta,phi,e);
                float tmpdr = r.DeltaR(g1);
                //GEN_REC_match_dr_zoom->Fill(tmpdr);
                if (tmpdr<dRcut && fabs(pSC_energy[i]-g1.E())<tmpde) {
                    tmpde = fabs(pSC_energy[i]-g1.E());
                    r1=r;
                    rind1=i;
                }
            }
        }
        for (int i=0;i<psize;i++) {
            if (REC_prompt[i] && i!=rind1 && pR9[i]>cutR9) {
                float tmpde=100000000;
                TLorentzVector r;
                float eta,phi,pt,e;
                pt  = pPt[i];
                eta = pEta[i];
                phi = pPhi[i];
                e   = pSC_energy[i];
                r.SetPtEtaPhiE(pt,eta,phi,e);
                float tmpdr = r.DeltaR(g2);
                if (tmpdr<dRcut && fabs(pSC_energy[i]-g2.E())<tmpde) {
                    tmpde = fabs(pSC_energy[i]-g2.E());
                    r2=r;
                    rind2=i;
                }
            }
        }
        /*
        int flagFake=0;
        for (int i=0;i<psize;i++) {
            if (i!=rind1 && i!=rind2 && pR9[i]>cutR9 && pHoverE[i]<0.05) {
                flagFake=1;
                TLorentzVector r;
                float eta,phi,pt,e;
                pt  = pPt[i];
                eta = pEta[i];
                phi = pPhi[i];
                e   = pSC_energy[i];
                r.SetPtEtaPhiE(pt,eta,phi,e);
                histo1d[14][0]->Fill(pPt[i]);
                histo1d[15][0]->Fill(pSC_energy[i]);
                histo1d[16][0]->Fill(pEta[i]);
                cout << "Fake " <<pPt[i] << " " << pEta[i] << " " << pPhi[i] << endl;
            }
        }
        if (flagFake) {
            if (gind1>=0) cout << "GEN PHO 1 " << g1.Pt() << " " << g1.Eta() << " " << g1.Phi() << endl;
            if (gind2>=0) cout << "GEN PHO 2 " << g2.Pt() << " " << g2.Eta() << " " << g2.Phi() << endl;
            if (rind1>=0) cout << "prompt PHO 1 " << pPt[rind1] << " " << pEta[rind1] << " " << pPhi[rind1] << endl;
            if (rind2>=0) cout << "prompt PHO 2 " << pPt[rind2] << " " << pEta[rind2] << " " << pPhi[rind2] << endl;
            cout << "----------------------------------------------------------------------------------------\n";
        }
        */
        
        /* // OLD: was picking up all prompt REC photons (see above), so, could have been more than 2
        for (int i=0;i<psize;i++) {
            if (REC_prompt[i]) {
                //float etaBin=pEta[i];
                TLorentzVector r;
                float eta,phi,pt,e;
                pt  = pPt[i];
                eta = pEta[i];
                phi = pPhi[i];
                e   = pSC_energy[i];
                r.SetPtEtaPhiE(pt,eta,phi,e);
                if (rind1==-1) {
                    r1=r;
                    rind1=i;
                }
                if (rind1!=-1 && rind2==-1 && i!=rind1) {
                    r2=r;
                    rind2=i;
                }
                float etaBin=r.Eta();
                for (int j=0; j<nEta; j++) {
                    histo1d[7][j]->Fill(r.Eta(),1.0/float(nEta));
                    if (fabs(etaBin) >= etaMin[j] && fabs(etaBin) < etaMax[j] ) {
                        histo1d[0][j]->Fill(pPt[i]);
                        histo1d[1][j]->Fill(pSC_energy[i]);
                        histo1d[2][j]->Fill(p5x5_energy[i]);
                        histo1d[3][j]->Fill(p3x3_energy[i]);
                        histo1d[8][j]->Fill(r.Phi());
                        histo1d[9][j]->Fill(pR9[i]);
                    }
                }
            }
        }
        */
        
        
        for (int i=0;i<psize;i++) {
            if (i==rind1 || i==rind2) {
                histo1d_tmp_pR9->Fill(pR9[i]);
                TLorentzVector r=r1, g=g1;
                if (i==rind2) {r=r2; g=g2;}
                float tmpdr = r.DeltaR(g);
                GEN_REC_match_dr_zoom->Fill(tmpdr);
                float etaBin=r.Eta();
                for (int j=0; j<nEta; j++) {
                    histo1d[7][j]->Fill(r.Eta(),1.0/float(nEta));
                    if (fabs(etaBin) >= etaMin[j] && fabs(etaBin) < etaMax[j] ) {
                        histo1d[0][j]->Fill(pPt[i]);
                        histo1d[1][j]->Fill(pRawenergy[i]);//pSC_energy[i]);
                        v_histo1d[1][j].push_back(pRawenergy[i]);//pSC_energy[i]);
                        histo1d[2][j]->Fill(p5x5_energy[i]);
                        v_histo1d[2][j].push_back(p5x5_energy[i]);
                        histo1d[3][j]->Fill(p3x3_energy[i]);
                        v_histo1d[3][j].push_back(p3x3_energy[i]);
                        histo1d[8][j]->Fill(r.Phi());
                        histo1d[9][j]->Fill(pR9[i]);
                        float genE=g1.E();
                        if (i==rind2) genE=g2.E();
                        histo1d[10][j]->Fill(pRawenergy[i]/genE);//pSC_energy[i]/genE);
                        v_histo1d[10][j].push_back(pRawenergy[i]/genE);//pSC_energy[i]/genE);
                        histo1d[11][j]->Fill(p5x5_energy[i]/genE);
                        v_histo1d[11][j].push_back(p5x5_energy[i]/genE);
                        histo1d[12][j]->Fill(p3x3_energy[i]/genE);
                        v_histo1d[12][j].push_back(p3x3_energy[i]/genE);
                        histo1d[13][j]->Fill(pSigmaIetaIeta[i]);

                        histo1d[17][j]->Fill(p_nCrystals[i]);
                        v_histo1d[17][j].push_back(p_nCrystals[i]);
                    }
                }
            }
        }
        
        
        // WORK with identified PAIRS:
        if (gind1!=-1 && gind2!=-1) {
            /* histo1dPairs[0][0]->Fill((g1+g2).M());
            if ((g1+g2).M()<5) {
                cout << g1.Eta() << " " << g1.Phi() << " " << g1.Pt() << " "
                << g2.Eta() << " " << g2.Phi() << " " << g2.Pt()
                << endl;
            } */
            for (int j=0; j<nEtaPairs; j++) {
                if (fabs(g1.Eta()) >= etaMinPairs[j] && fabs(g1.Eta()) < etaMaxPairs[j]
                    && fabs(g2.Eta()) >= etaMinPairs[j] && fabs(g2.Eta()) < etaMaxPairs[j]
                    ) {
                    //histo1dPairs[0][j]->Fill((g1+g2).M());
                    histo1dPairs[0][j]->Fill(invMass(g1.E(),g1.Eta(),g1.Phi(),g2.E(),g2.Eta(),g2.Phi()));
                }
            }
        }
        if (rind1!=-1 && rind2!=-1) {
            for (int j=0; j<nEtaPairs; j++) {
                if (fabs(r1.Eta()) >= etaMinPairs[j] && fabs(r1.Eta()) < etaMaxPairs[j]
                    && fabs(r2.Eta()) >= etaMinPairs[j] && fabs(r2.Eta()) < etaMaxPairs[j]
                    ) {
                    //histo1dPairs[1][j]->Fill((r1+r2).M());
                    float m1=invMass(r1.E(),r1.Eta(),r1.Phi(),r2.E(),r2.Eta(),r2.Phi());
                    histo1dPairs[1][j]->Fill(m1);
                    v_m1.push_back(m1);
                    /*
                    TLorentzVector r1_e5x5;
                    {
                        float eta,phi,pt,e;
                        pt  = pPt[rind1];
                        eta = pEta[rind1];
                        phi = pPhi[rind1];
                        e   = p5x5_energy[rind1];
                        r1_e5x5.SetPtEtaPhiE(pt,eta,phi,e);
                    }
                    TLorentzVector r2_e5x5;
                    {
                        float eta,phi,pt,e;
                        pt  = pPt[rind2];
                        eta = pEta[rind2];
                        phi = pPhi[rind2];
                        e   = p5x5_energy[rind2];
                        r2_e5x5.SetPtEtaPhiE(pt,eta,phi,e);
                    }
                    //histo1dPairs[2][j]->Fill((r1_e5x5+r2_e5x5).M());*/
                    float m2 = invMass(p5x5_energy[rind1],r1.Eta(),r1.Phi(),p5x5_energy[rind2],r2.Eta(),r2.Phi());
                    float m3 = invMass(p3x3_energy[rind1],r1.Eta(),r1.Phi(),p3x3_energy[rind2],r2.Eta(),r2.Phi());
                    histo1dPairs[2][j]->Fill(m2);
                    histo1dPairs[3][j]->Fill(m3);
                    v_m2.push_back(m2);
                    v_m3.push_back(m3);
                }
            }
        }
        
        
    }//MAIN LOOP: end
    
    histo1d[14][0]->Scale(1.0/float(nentries));
    histo1d[15][0]->Scale(1.0/float(nentries));
    histo1d[16][0]->Scale(1.0/float(nentries));
    
    /* sort(v_m1.begin(),v_m1.end());
    for (int i=0; i<v_m1.size(); i++) {
        cout << v_m1[i] << " ";
    }
    cout << endl; */
    
    ofstream myfile;
	myfile.open("results.out2",std::ofstream::app);
    
    { // energy resolutions vs. E inside particular eta slice
        myfile << "---------------- START: energy resolutions vs. E inside particular eta slice ------------------\n";
        const int nSteps=5;
        for (int j=0; j<nEta; j++) {
            myfile << "EtaSlice: " << etaMin[j] << ".." << etaMax[j] << endl;
            vector<float> v1 = v_histo1d[2][j];//1: E_SC (maybe raw energy); 2: 5x5; 3: 3x3
            vector<float> v2 = v_histo1d[11][j];//11 = e5x5/genE; 12 = e3x3/genE
            //cout << "1" << endl << flush;
            sortKeepAligned(v1,v2);
            //cout << "2" << endl << flush;
            const int n=v_histo1d[2][j].size();//v1.size();
            myfile << "size is: " << n << endl;
            const int n1=int(float(n)/float(nSteps));
            myfile << "size n1 is: " << n1 << endl;
            //cout << "3" << endl << flush;
            for (int i=0; i<nSteps; i++) {
                vector<float> w1,w2;
                int minI=i*n1;
                int maxI=(1+i)*n1;
                if (i==nSteps-1 || maxI>n) maxI = n;
                myfile << "index interval: " << minI << " " << maxI << endl;
                float meanX=0;
                float meanXcount=0;
                for (int ii=minI; ii<maxI; ii++) {
                    w1.push_back(v1[ii]);
                    meanX+=v1[ii];
                    meanXcount++;
                    w2.push_back(v2[ii]);
                }
                float sigma = 0.5*defineSmallestIntervalToContain(0.68,w2);
                meanX=meanX/meanXcount;
                myfile << "E interval: " << setprecision(3) << w1[0] << ".." << w1[w1.size()-1] << " mean energy: " << meanX << " resolution sigma: " << sigma << endl;
                myfile << "Excel2 " << etaMin[j] << " " << etaMax[j]
                << " " << fabs(etaMin[j]-etaMax[j])/2.0 + etaMin[j]
                << " " << w1[0] << " " << w1[w1.size()-1]
                << " " << meanX << " " << sigma
                << " " << n1
                << endl;
            }
            //cout << endl;
        }
        cout << "---------------- END: energy resolutions vs. E inside particular eta slice ------------------\n";
    }
    
    myfile << "68%_interval_for_Esc_(inv._mass): " << defineSmallestIntervalToContain(0.68,v_m1)/2.0 << endl;
    myfile << "68%_interval_for_E5x5_(inv._mass): " << defineSmallestIntervalToContain(0.68,v_m2)/2.0 << endl;
    myfile << "68%_interval_for_E3x3_(inv._mass): " << defineSmallestIntervalToContain(0.68,v_m3)/2.0 << endl;
    
    for (int i=0; i<nVar; i++) {
        if (i>=10 && i<13) {
            if (i==10) {
                myfile << "ETA_INTERVALS:";
                for (int j=0; j<nEta; j++) {
                    myfile << " " << etaMin[j] << ".." << etaMax[j];
                }
                myfile << endl;
                myfile << "ETA_MID_BINS:";
                for (int j=0; j<nEta; j++) {
                    myfile << " " << fabs(etaMin[j]-etaMax[j])/2.0 + etaMin[j];
                }
                myfile << endl;
            }
            myfile << varName[i] << "_68%_interval_for_individual_photons_vs._eta:";
            for (int j=0; j<nEta; j++) {
                myfile << " " << setprecision(3) << defineSmallestIntervalToContain(0.68,v_histo1d[i][j])/2.0;
            }
            myfile << endl;
        }
    }
    for (int i=0; i<nVar; i++) {
        if (i==10) {
            myfile << "ETA_INTERVALS:";
            for (int j=0; j<nEta; j++) {
                myfile << " " << etaMin[j] << ".." << etaMax[j];
            }
            myfile << endl;
            myfile << "ETA_MID_BINS:";
            for (int j=0; j<nEta; j++) {
                myfile << " " << fabs(etaMin[j]-etaMax[j])/2.0 + etaMin[j];
            }
            myfile << endl;
            myfile << "Reconstruction_efficiency_for_individual_photons_vs._eta:";
            for (int j=0; j<nEta; j++) {
                if (v_histo1d[6][j].size()!=0) myfile << " " << setprecision(3) << float(v_histo1d[10][j].size())/float(v_histo1d[6][j].size());
                else myfile << " nan";
            }
            myfile << endl;
            myfile << "GEN_counts_for_reconstruction_efficiency_for_individual_photons_vs._eta:";
            for (int j=0; j<nEta; j++) {
                myfile << " " << setprecision(3) << v_histo1d[6][j].size();
            }
            myfile << endl;
            myfile << "REC_counts_for_reconstruction_efficiency_for_individual_photons_vs._eta:";
            for (int j=0; j<nEta; j++) {
                myfile << " " << setprecision(3) << v_histo1d[10][j].size();
            }
            myfile << endl;
        }
    }
    for (int j=0; j<nEta; j++) {
        myfile << "forExcel " << fabs(etaMin[j]-etaMax[j])/2.0 + etaMin[j];
        for (int i=0; i<nVar; i++) {
            if (i>=10 && i<13) {
                myfile << " " << setprecision(3) << defineSmallestIntervalToContain(0.68,v_histo1d[i][j])/2.0;
            }
        }
        for (int i=0; i<nVar; i++) {
            if (i==10) {
                //myfile << " " << setprecision(3) << defineSmallestIntervalToContain(0.68,v_histo1d[i][j]);
                if (v_histo1d[6][j].size()!=0) myfile << " " << setprecision(3) << float(v_histo1d[10][j].size())/float(v_histo1d[6][j].size());
                else myfile << " nan";
                myfile << " " << setprecision(3) << v_histo1d[6][j].size();
                myfile << " " << setprecision(3) << v_histo1d[10][j].size();
            }
        }
        myfile << endl;
    }
    for (int j=0; j<nEta; j++) {
        myfile << "NXTALS " << fabs(etaMin[j]-etaMax[j])/2.0 + etaMin[j];
        for (int i=0; i<nVar; i++) {
            if (i==17) {
                myfile << " " << setprecision(3) << defineAverage(v_histo1d[i][j]);
            }
        }
        myfile << endl;
    }
    histo1d_tmp_pR9->Scale(1.0/float(nentries));
    for (int i=0; i<150; i++) {
        myfile << "r9_var_printout " << histo1d_tmp_pR9->GetBinCenter(i+1) << " " << histo1d_tmp_pR9->GetBinContent(i+1) << endl;
    }

    /* for (int j=0; j<nEta; j++) {
        myfile << fabs(etaMin[j]-etaMax[j])/2.0 + etaMin[j];
        for (int i=0; i<nVar; i++) {
            if (i==10) {
                //myfile << " " << setprecision(3) << defineSmallestIntervalToContain(0.68,v_histo1d[i][j]);
                if (v_histo1d[6][j].size()!=0) myfile << " " << setprecision(3) << float(v_histo1d[10][j].size())/float(v_histo1d[6][j].size());
                else myfile << " nan";
                myfile << " " << setprecision(3) << v_histo1d[6][j].size();
                myfile << " " << setprecision(3) << v_histo1d[10][j].size();
            }
        }
        myfile << endl;
    }*/
    
    
    myfile.close();
    
    
    if (1) {
        
        int logOY=0;
        setTDRStyle1(1,1,0,logOY);//gridX/Y, logX/Y
        gStyle->TStyle::SetOptStat(0);
        gStyle->TStyle::SetOptFit(0);
        TCanvas *c1 = new TCanvas("c1","myPlots",0,0,600,600);
        c1->cd(1);

        {
            gStyle->TStyle::SetOptStat(111110);
            for (int i=0; i<nVar; i++) {
                if (i<4 || i==13) {//ET and all energy variables
                    for (int j=0; j<nEta; j++) {
                        string etaName = "#eta#in[" + tostr(etaMin[j]) + ";" + tostr(etaMax[j]) + "]: ";
                        string etaNameF = "eta_" + tostr(etaMin[j]) + "_to_" + tostr(etaMax[j]);
                        const TString xTitle = etaName + varNameXaxis[i];
                        histo1d[i][j]->GetXaxis()->SetTitle(xTitle);
                        histo1d[i][j]->Draw();
                        for (int ii=0; ii<nFigTypes; ii++) {
                            if (saveThis[ii]) {
                                TString PP = plotsPrefix;
                                if (ii==2) PP = plotsPrefixC;
                                TString name1;
                                name1 = PP + histoPrefix
                                + "_PHOTON_" + varName[i] + etaNameF
                                + figType[ii];
                                TString ss = name1;
                                c1->SaveAs(ss);
                            }
                        }
                        c1->Clear();
                    }
                }
            }
            gStyle->TStyle::SetOptStat(0);
        }
        
        
        {
            for (int i=0; i<nVarPairs; i++) {
                for (int j=0; j<nEtaPairs; j++) {
                    if (j==0) histo1dPairs[i][j]->Draw();
                    else histo1dPairs[i][j]->Draw("same");
                }
                for (int ii=0; ii<nFigTypes; ii++) {
                    if (saveThis[ii]) {
                        TString PP = plotsPrefix;
                        if (ii==2) PP = plotsPrefixC;
                        TString name1;
                        name1 = PP + histoPrefix
                        + "_PHOTON_PAIR_" + varNamePairs[i]
                        + figType[ii];
                        TString ss = name1;
                        c1->SaveAs(ss);
                    }
                }
                c1->Clear();
            }
        }
        
        {
            for (int i=0; i<nVar; i++) {
                for (int j=0; j<nEta; j++) {
                    if (j==0) histo1d[i][j]->Draw();
                    else histo1d[i][j]->Draw("same");
                }
                for (int ii=0; ii<nFigTypes; ii++) {
                    if (saveThis[ii]) {
                        TString PP = plotsPrefix;
                        if (ii==2) PP = plotsPrefixC;
                        TString name1;
                        name1 = PP + histoPrefix
                        + "_PHOTON_" + varName[i]
                        + figType[ii];
                        TString ss = name1;
                        c1->SaveAs(ss);
                    }
                }
                c1->Clear();
            }
        }
        
        {
            c1->SetLogy(1);
            string namedr = tostr(dRcut);
			TString name1 = "\# of REC photons within " + namedr;
			const TString xTitle = name1;
			nREC_in_dR_hist->GetXaxis()->SetTitle(xTitle);
			const TString yTitle = "Number of entries";
			nREC_in_dR_hist->GetYaxis()->SetTitle(yTitle);
			nREC_in_dR_hist->GetYaxis()->SetTitleOffset(2.0);
			nREC_in_dR_hist->GetXaxis()->SetNdivisions(505);
			nREC_in_dR_hist->GetYaxis()->SetNdivisions(505);
			nREC_in_dR_hist->SetLineColor(1);
			nREC_in_dR_hist->SetLineStyle(1);
			nREC_in_dR_hist->SetLineWidth(2);
			nREC_in_dR_hist->SetMarkerStyle(21);
			nREC_in_dR_hist->SetMarkerColor(1);
			nREC_in_dR_hist->SetMarkerSize(0.1);
			nREC_in_dR_hist->SetFillColor(1);
			nREC_in_dR_hist->SetFillStyle(0);
			//nREC_in_dR_hist->SetMinimum(0);
			nREC_in_dR_hist->Draw();
			for (int ii=0; ii<nFigTypes; ii++) {
				if (saveThis[ii]) {
					TString PP = plotsPrefix;
					if (ii==2) PP = plotsPrefixC;
					TString name1;
					name1 = PP + histoPrefix
					+ "_nREC_in_dR"
					+ figType[ii];
					TString ss = name1;
					c1->SaveAs(ss);
				}
			}
			c1->Clear();
            c1->SetLogy(0);
		}
        
        {
            c1->SetLogy(1);
			TString name1 = "min(dR(#eta,#phi)) for all GEN and REC photons";
			const TString xTitle = name1;
			GEN_REC_match_dr->GetXaxis()->SetTitle(xTitle);
			const TString yTitle = "Number of entries";
			GEN_REC_match_dr->GetYaxis()->SetTitle(yTitle);
			GEN_REC_match_dr->GetYaxis()->SetTitleOffset(2.0);
			GEN_REC_match_dr->GetXaxis()->SetNdivisions(505);
			GEN_REC_match_dr->GetYaxis()->SetNdivisions(505);
			GEN_REC_match_dr->SetLineColor(1);
			GEN_REC_match_dr->SetLineStyle(1);
			GEN_REC_match_dr->SetLineWidth(2);
			GEN_REC_match_dr->SetMarkerStyle(21);
			GEN_REC_match_dr->SetMarkerColor(1);
			GEN_REC_match_dr->SetMarkerSize(0.1);
			GEN_REC_match_dr->SetFillColor(1);
			GEN_REC_match_dr->SetFillStyle(0);
			//GEN_REC_match_dr->SetMinimum(0);
			GEN_REC_match_dr->Draw();
			for (int ii=0; ii<nFigTypes; ii++) {
				if (saveThis[ii]) {
					TString PP = plotsPrefix;
					if (ii==2) PP = plotsPrefixC;
					TString name1;
					name1 = PP + histoPrefix
					+ "_GEN_REC_mindR"
					+ figType[ii];
					TString ss = name1;
					c1->SaveAs(ss);
				}
			}
			c1->Clear();
            c1->SetLogy(0);
		}
        
		{
            c1->SetLogy(1);
			TString name1 = "min(dR(#eta,#phi)) for all GEN and REC photons";
			const TString xTitle = name1;
			GEN_REC_match_dr_zoom->GetXaxis()->SetTitle(xTitle);
			const TString yTitle = "Number of entries";
			GEN_REC_match_dr_zoom->GetYaxis()->SetTitle(yTitle);
			GEN_REC_match_dr_zoom->GetYaxis()->SetTitleOffset(2.0);
			GEN_REC_match_dr_zoom->GetXaxis()->SetNdivisions(505);
			GEN_REC_match_dr_zoom->GetYaxis()->SetNdivisions(505);
			GEN_REC_match_dr_zoom->SetLineColor(1);
			GEN_REC_match_dr_zoom->SetLineStyle(1);
			GEN_REC_match_dr_zoom->SetLineWidth(2);
			GEN_REC_match_dr_zoom->SetMarkerStyle(21);
			GEN_REC_match_dr_zoom->SetMarkerColor(1);
			GEN_REC_match_dr_zoom->SetMarkerSize(0.1);
			GEN_REC_match_dr_zoom->SetFillColor(1);
			GEN_REC_match_dr_zoom->SetFillStyle(0);
			//GEN_REC_match_dr_zoom->SetMinimum(0);
			GEN_REC_match_dr_zoom->Draw();
			for (int ii=0; ii<nFigTypes; ii++) {
				if (saveThis[ii]) {
					TString PP = plotsPrefix;
					if (ii==2) PP = plotsPrefixC;
					TString name1;
					name1 = PP + histoPrefix
					+ "_GEN_REC_mindR_zoom"
					+ figType[ii];
					TString ss = name1;
					c1->SaveAs(ss);
				}
			}
			c1->Clear();
            c1->SetLogy(0);
		}
        
    }
    
} //void analysis::LoopPho()








void analysis::Loop()
{
    //... last time in: backup_2014_03_18_17_11_07_before_new_photon_analysis_with_FIXED_INPUT
    // NOW: working on smaller version of the input ROOT trees, no ADC details for electrons, etc.
}// END: void analysis::Loop()




















// -------------------------------------------------------------------------------------                                                                                           
// ------------- don't go there...                                                                                                                                                 
// -------------------------------------------------------------------------------------                                                                                           

//#ifdef analysis_cxx
analysis::analysis(TTree *tree)
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0) {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("input.root");
        if (!f) {
            f = new TFile("input.root");
            f->cd("input.root:/simpleNtupleEoverP");
        }
        tree = (TTree*)gDirectory->Get("SimpleNtupleEoverP");
        
    }
    Init(tree);
}

analysis::~analysis()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t analysis::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t analysis::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (!fChain->InheritsFrom(TChain::Class()))  return centry;
    TChain *chain = (TChain*)fChain;
    if (chain->GetTreeNumber() != fCurrent) {
        fCurrent = chain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void analysis::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).
    
    // Set object pointer
    ele1_recHit_E = 0;
    ele1_recHit_flag = 0;
    ele1_recHit_hashedIndex = 0;
    ele1_recHit_ietaORix = 0;
    ele1_recHit_iphiORiy = 0;
    ele1_recHit_zside = 0;
    ele1_recHit_laserCorrection = 0;
    ele1_recHit_Alpha = 0;
    ele1_recHit_ICConstant = 0;
    ele2_recHit_E = 0;
    ele2_recHit_flag = 0;
    ele2_recHit_hashedIndex = 0;
    ele2_recHit_ietaORix = 0;
    ele2_recHit_iphiORiy = 0;
    ele2_recHit_zside = 0;
    ele2_recHit_laserCorrection = 0;
    ele2_recHit_Alpha = 0;
    ele2_recHit_ICConstant = 0;
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);
    
    fChain->SetBranchAddress("bxId", &bxId, &b_bxId);
    fChain->SetBranchAddress("eventId", &eventId, &b_eventId);
    fChain->SetBranchAddress("lumiId", &lumiId, &b_lumiId);
    fChain->SetBranchAddress("runId", &runId, &b_runId);
    fChain->SetBranchAddress("timeStampHigh", &timeStampHigh, &b_timeStampHigh);
    fChain->SetBranchAddress("isW", &isW, &b_isW);
    fChain->SetBranchAddress("isZ", &isZ, &b_isZ);
    fChain->SetBranchAddress("PV_n", &PV_n, &b_PV_n);
    fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
    fChain->SetBranchAddress("PV_d0", &PV_d0, &b_PV_d0);
    fChain->SetBranchAddress("PV_SumPt", &PV_SumPt, &b_PV_SumPt);
    fChain->SetBranchAddress("rho", &rho, &b_rho);
    fChain->SetBranchAddress("NELE", &NELE, &b_NELE);
    fChain->SetBranchAddress("ELE_pt", ELE_pt, &b_ELE_pt);
    fChain->SetBranchAddress("ELE_eta", ELE_eta, &b_ELE_eta);
    fChain->SetBranchAddress("ELE_phi", ELE_phi, &b_ELE_phi);
    fChain->SetBranchAddress("ELE_px", ELE_px, &b_ELE_px);
    fChain->SetBranchAddress("ELE_py", ELE_py, &b_ELE_py);
    fChain->SetBranchAddress("ELE_pz", ELE_pz, &b_ELE_pz);
    fChain->SetBranchAddress("ELE_e", ELE_e, &b_ELE_e);
    fChain->SetBranchAddress("ELE_q", ELE_q, &b_ELE_q);
    fChain->SetBranchAddress("ELE_id", ELE_id, &b_ELE_id);
    fChain->SetBranchAddress("ELE_sigmaIetaIeta", ELE_sigmaIetaIeta, &b_ELE_sigmaIetaIeta);
    fChain->SetBranchAddress("ELE_DphiIn", ELE_DphiIn, &b_ELE_DphiIn);
    fChain->SetBranchAddress("ELE_DetaIn", ELE_DetaIn, &b_ELE_DetaIn);
    fChain->SetBranchAddress("ELE_HOverE", ELE_HOverE, &b_ELE_HOverE);
    fChain->SetBranchAddress("ELE_ooemoop", ELE_ooemoop, &b_ELE_ooemoop);
    fChain->SetBranchAddress("ELE_tkIso", ELE_tkIso, &b_ELE_tkIso);
    fChain->SetBranchAddress("ELE_emIso", ELE_emIso, &b_ELE_emIso);
    fChain->SetBranchAddress("ELE_hadIso", ELE_hadIso, &b_ELE_hadIso);
    fChain->SetBranchAddress("ELE_effAreaForIso", ELE_effAreaForIso, &b_ELE_effAreaForIso);
    fChain->SetBranchAddress("ELE_combIso", ELE_combIso, &b_ELE_combIso);
    fChain->SetBranchAddress("ELE_dxy", ELE_dxy, &b_ELE_dxy);
    fChain->SetBranchAddress("ELE_dz", ELE_dz, &b_ELE_dz);
    fChain->SetBranchAddress("ELE_scE_regression", ELE_scE_regression, &b_ELE_scE_regression);
    fChain->SetBranchAddress("ELE_scEtaWidth", ELE_scEtaWidth, &b_ELE_scEtaWidth);
    fChain->SetBranchAddress("ELE_scPhiWidth", ELE_scPhiWidth, &b_ELE_scPhiWidth);
    fChain->SetBranchAddress("ELE_scERaw", ELE_scERaw, &b_ELE_scERaw);
    fChain->SetBranchAddress("ELE_scE", ELE_scE, &b_ELE_scE);
    fChain->SetBranchAddress("ELE_e5x5", ELE_e5x5, &b_ELE_e5x5);
    fChain->SetBranchAddress("ELE_e1x5", ELE_e1x5, &b_ELE_e1x5);
    fChain->SetBranchAddress("ELE_e2x5Max", ELE_e2x5Max, &b_ELE_e2x5Max);
    fChain->SetBranchAddress("ELE_es", ELE_es, &b_ELE_es);
    fChain->SetBranchAddress("ELE_fbrem", ELE_fbrem, &b_ELE_fbrem);
    fChain->SetBranchAddress("psize", &psize, &b_psize);
    fChain->SetBranchAddress("pSC_energy", pSC_energy, &b_pSC_energy);
    fChain->SetBranchAddress("pEta", pEta, &b_pEta);
    fChain->SetBranchAddress("pPhi", pPhi, &b_pPhi);
    fChain->SetBranchAddress("pEsenergy", pEsenergy, &b_pEsenergy);
    fChain->SetBranchAddress("pRawenergy", pRawenergy, &b_pRawenergy);
    fChain->SetBranchAddress("pPt", pPt, &b_pPt);
    fChain->SetBranchAddress("p1x5_energy", p1x5_energy, &b_p1x5_energy);
    fChain->SetBranchAddress("p2x5_energy", p2x5_energy, &b_p2x5_energy);
    fChain->SetBranchAddress("p3x3_energy", p3x3_energy, &b_p3x3_energy);
    fChain->SetBranchAddress("p5x5_energy", p5x5_energy, &b_p5x5_energy);
    fChain->SetBranchAddress("pR9", pR9, &b_pR9);
    fChain->SetBranchAddress("pGap", pGap, &b_pGap);
    fChain->SetBranchAddress("pSigmaIetaIeta", pSigmaIetaIeta, &b_pSigmaIetaIeta);
    fChain->SetBranchAddress("pHoverE", pHoverE, &b_pHoverE);
    fChain->SetBranchAddress("p_photonenergy", p_photonenergy, &b_p_photonenergy);
    fChain->SetBranchAddress("p_nCrystals", p_nCrystals, &b_p_nCrystals);
    fChain->SetBranchAddress("p_adc", p_adc, &b_p_adc);
    fChain->SetBranchAddress("p_gainId", p_gainId, &b_p_gainId);
    fChain->SetBranchAddress("p_intercalib", p_intercalib, &b_p_intercalib);
    fChain->SetBranchAddress("p_recieta", p_recieta, &b_p_recieta);
    fChain->SetBranchAddress("p_reciphi", p_reciphi, &b_p_reciphi);
    fChain->SetBranchAddress("p_receta", p_receta, &b_p_receta);
    fChain->SetBranchAddress("p_recphi", p_recphi, &b_p_recphi);
    fChain->SetBranchAddress("p_recenergy", p_recenergy, &b_p_recenergy);
    fChain->SetBranchAddress("p_rectime", p_rectime, &b_p_rectime);
    fChain->SetBranchAddress("p_recflag", p_recflag, &b_p_recflag);
    fChain->SetBranchAddress("p_recflags", p_recflags, &b_p_recflags);
    fChain->SetBranchAddress("NGEN", &NGEN, &b_NGEN);
    fChain->SetBranchAddress("GEN_pt", GEN_pt, &b_GEN_pt);
    fChain->SetBranchAddress("GEN_eta", GEN_eta, &b_GEN_eta);
    fChain->SetBranchAddress("GEN_phi", GEN_phi, &b_GEN_phi);
    fChain->SetBranchAddress("GEN_px", GEN_px, &b_GEN_px);
    fChain->SetBranchAddress("GEN_py", GEN_py, &b_GEN_py);
    fChain->SetBranchAddress("GEN_pz", GEN_pz, &b_GEN_pz);
    fChain->SetBranchAddress("GEN_e", GEN_e, &b_GEN_e);
    fChain->SetBranchAddress("GEN_mass", GEN_mass, &b_GEN_mass);
    fChain->SetBranchAddress("GEN_q", GEN_q, &b_GEN_q);
    fChain->SetBranchAddress("GEN_id", GEN_id, &b_GEN_id);
    fChain->SetBranchAddress("GEN_status", GEN_status, &b_GEN_status);
    fChain->SetBranchAddress("GEN_parent", GEN_parent, &b_GEN_parent);
    fChain->SetBranchAddress("ele1_charge", &ele1_charge, &b_ele1_charge);
    fChain->SetBranchAddress("ele1_p", &ele1_p, &b_ele1_p);
    fChain->SetBranchAddress("ele1_pt", &ele1_pt, &b_ele1_pt);
    fChain->SetBranchAddress("ele1_eta", &ele1_eta, &b_ele1_eta);
    fChain->SetBranchAddress("ele1_phi", &ele1_phi, &b_ele1_phi);
    fChain->SetBranchAddress("ele1_isEB", &ele1_isEB, &b_ele1_isEB);
    fChain->SetBranchAddress("ele1_isEBEEGap", &ele1_isEBEEGap, &b_ele1_isEBEEGap);
    fChain->SetBranchAddress("ele1_isEBEtaGap", &ele1_isEBEtaGap, &b_ele1_isEBEtaGap);
    fChain->SetBranchAddress("ele1_isEBPhiGap", &ele1_isEBPhiGap, &b_ele1_isEBPhiGap);
    fChain->SetBranchAddress("ele1_isEEDeeGap", &ele1_isEEDeeGap, &b_ele1_isEEDeeGap);
    fChain->SetBranchAddress("ele1_isEERingGap", &ele1_isEERingGap, &b_ele1_isEERingGap);
    fChain->SetBranchAddress("ele1_isTrackerDriven", &ele1_isTrackerDriven, &b_ele1_isTrackerDriven);
    fChain->SetBranchAddress("ele1_sigmaIetaIeta", &ele1_sigmaIetaIeta, &b_ele1_sigmaIetaIeta);
    fChain->SetBranchAddress("ele1_DphiIn", &ele1_DphiIn, &b_ele1_DphiIn);
    fChain->SetBranchAddress("ele1_DetaIn", &ele1_DetaIn, &b_ele1_DetaIn);
    fChain->SetBranchAddress("ele1_HOverE", &ele1_HOverE, &b_ele1_HOverE);
    fChain->SetBranchAddress("ele1_tkIso", &ele1_tkIso, &b_ele1_tkIso);
    fChain->SetBranchAddress("ele1_emIso", &ele1_emIso, &b_ele1_emIso);
    fChain->SetBranchAddress("ele1_hadIso", &ele1_hadIso, &b_ele1_hadIso);
    fChain->SetBranchAddress("ele1_scERaw", &ele1_scERaw, &b_ele1_scERaw);
    fChain->SetBranchAddress("ele1_scEtRaw", &ele1_scEtRaw, &b_ele1_scEtRaw);
    fChain->SetBranchAddress("ele1_scE", &ele1_scE, &b_ele1_scE);
    fChain->SetBranchAddress("ele1_scEt", &ele1_scEt, &b_ele1_scEt);
    fChain->SetBranchAddress("ele1_scE_regression", &ele1_scE_regression, &b_ele1_scE_regression);
    fChain->SetBranchAddress("ele1_scEerr_regression", &ele1_scEerr_regression, &b_ele1_scEerr_regression);
    fChain->SetBranchAddress("ele1_scE_regression_PhotonTuned", &ele1_scE_regression_PhotonTuned, &b_ele1_scE_regression_PhotonTuned);
    fChain->SetBranchAddress("ele1_scEerr_regression_PhotonTuned", &ele1_scEerr_regression_PhotonTuned, &b_ele1_scEerr_regression_PhotonTuned);
    fChain->SetBranchAddress("ele1_scERaw_PUcleaned", &ele1_scERaw_PUcleaned, &b_ele1_scERaw_PUcleaned);
    fChain->SetBranchAddress("ele1_es", &ele1_es, &b_ele1_es);
    fChain->SetBranchAddress("ele1_scLaserCorr", &ele1_scLaserCorr, &b_ele1_scLaserCorr);
    fChain->SetBranchAddress("ele1_scCrackCorr", &ele1_scCrackCorr, &b_ele1_scCrackCorr);
    fChain->SetBranchAddress("ele1_scLocalContCorr", &ele1_scLocalContCorr, &b_ele1_scLocalContCorr);
    fChain->SetBranchAddress("ele1_scEta", &ele1_scEta, &b_ele1_scEta);
    fChain->SetBranchAddress("ele1_scPhi", &ele1_scPhi, &b_ele1_scPhi);
    fChain->SetBranchAddress("ele1_scLocalEta", &ele1_scLocalEta, &b_ele1_scLocalEta);
    fChain->SetBranchAddress("ele1_scLocalPhi", &ele1_scLocalPhi, &b_ele1_scLocalPhi);
    fChain->SetBranchAddress("ele1_scEtaWidth", &ele1_scEtaWidth, &b_ele1_scEtaWidth);
    fChain->SetBranchAddress("ele1_scPhiWidth", &ele1_scPhiWidth, &b_ele1_scPhiWidth);
    fChain->SetBranchAddress("ele1_scEtaWidth_PUcleaned", &ele1_scEtaWidth_PUcleaned, &b_ele1_scEtaWidth_PUcleaned);
    fChain->SetBranchAddress("ele1_scPhiWidth_PUcleaned", &ele1_scPhiWidth_PUcleaned, &b_ele1_scPhiWidth_PUcleaned);
    fChain->SetBranchAddress("ele1_fCorrection_PUcleaned", &ele1_fCorrection_PUcleaned, &b_ele1_fCorrection_PUcleaned);
    fChain->SetBranchAddress("ele1_fEta", &ele1_fEta, &b_ele1_fEta);
    fChain->SetBranchAddress("ele1_fEtaCorr", &ele1_fEtaCorr, &b_ele1_fEtaCorr);
    fChain->SetBranchAddress("ele1_tkP", &ele1_tkP, &b_ele1_tkP);
    fChain->SetBranchAddress("ele1_tkPt", &ele1_tkPt, &b_ele1_tkPt);
    fChain->SetBranchAddress("ele1_fbrem", &ele1_fbrem, &b_ele1_fbrem);
    fChain->SetBranchAddress("ele1_dxy_PV", &ele1_dxy_PV, &b_ele1_dxy_PV);
    fChain->SetBranchAddress("ele1_dz_PV", &ele1_dz_PV, &b_ele1_dz_PV);
    fChain->SetBranchAddress("ele1_sigmaP", &ele1_sigmaP, &b_ele1_sigmaP);
    fChain->SetBranchAddress("ele1_eSeedBC", &ele1_eSeedBC, &b_ele1_eSeedBC);
    fChain->SetBranchAddress("ele1_e5x5", &ele1_e5x5, &b_ele1_e5x5);
    fChain->SetBranchAddress("ele1_e3x3", &ele1_e3x3, &b_ele1_e3x3);
    fChain->SetBranchAddress("ele1_scNxtal", &ele1_scNxtal, &b_ele1_scNxtal);
    fChain->SetBranchAddress("ele1_bcN", &ele1_bcN, &b_ele1_bcN);
    fChain->SetBranchAddress("ele1_5x5LaserCorr", &ele1_5x5LaserCorr, &b_ele1_5x5LaserCorr);
    fChain->SetBranchAddress("ele1_3x3LaserCorr", &ele1_3x3LaserCorr, &b_ele1_3x3LaserCorr);
    fChain->SetBranchAddress("ele1_seedE", &ele1_seedE, &b_ele1_seedE);
    fChain->SetBranchAddress("ele1_seedLaserAlpha", &ele1_seedLaserAlpha, &b_ele1_seedLaserAlpha);
    fChain->SetBranchAddress("ele1_seedLaserCorr", &ele1_seedLaserCorr, &b_ele1_seedLaserCorr);
    fChain->SetBranchAddress("ele1_seedICConstant", &ele1_seedICConstant, &b_ele1_seedICConstant);
    fChain->SetBranchAddress("ele1_seedIeta", &ele1_seedIeta, &b_ele1_seedIeta);
    fChain->SetBranchAddress("ele1_seedIphi", &ele1_seedIphi, &b_ele1_seedIphi);
    fChain->SetBranchAddress("ele1_seedIx", &ele1_seedIx, &b_ele1_seedIx);
    fChain->SetBranchAddress("ele1_seedIy", &ele1_seedIy, &b_ele1_seedIy);
    fChain->SetBranchAddress("ele1_seedZside", &ele1_seedZside, &b_ele1_seedZside);
    fChain->SetBranchAddress("ele1_EOverP", &ele1_EOverP, &b_ele1_EOverP);
    fChain->SetBranchAddress("ele1_nRecHits", &ele1_nRecHits, &b_ele1_nRecHits);
    fChain->SetBranchAddress("ele1_recHit_E", &ele1_recHit_E, &b_ele1_recHit_E);
    fChain->SetBranchAddress("ele1_recHit_flag", &ele1_recHit_flag, &b_ele1_recHit_flag);
    fChain->SetBranchAddress("ele1_recHit_hashedIndex", &ele1_recHit_hashedIndex, &b_ele1_recHit_hashedIndex);
    fChain->SetBranchAddress("ele1_recHit_ietaORix", &ele1_recHit_ietaORix, &b_ele1_recHit_ietaORix);
    fChain->SetBranchAddress("ele1_recHit_iphiORiy", &ele1_recHit_iphiORiy, &b_ele1_recHit_iphiORiy);
    fChain->SetBranchAddress("ele1_recHit_zside", &ele1_recHit_zside, &b_ele1_recHit_zside);
    fChain->SetBranchAddress("ele1_recHit_laserCorrection", &ele1_recHit_laserCorrection, &b_ele1_recHit_laserCorrection);
    fChain->SetBranchAddress("ele1_recHit_Alpha", &ele1_recHit_Alpha, &b_ele1_recHit_Alpha);
    fChain->SetBranchAddress("ele1_recHit_ICConstant", &ele1_recHit_ICConstant, &b_ele1_recHit_ICConstant);
    fChain->SetBranchAddress("ele2_charge", &ele2_charge, &b_ele2_charge);
    fChain->SetBranchAddress("ele2_p", &ele2_p, &b_ele2_p);
    fChain->SetBranchAddress("ele2_pt", &ele2_pt, &b_ele2_pt);
    fChain->SetBranchAddress("ele2_eta", &ele2_eta, &b_ele2_eta);
    fChain->SetBranchAddress("ele2_phi", &ele2_phi, &b_ele2_phi);
    fChain->SetBranchAddress("ele2_isEB", &ele2_isEB, &b_ele2_isEB);
    fChain->SetBranchAddress("ele2_isEBEEGap", &ele2_isEBEEGap, &b_ele2_isEBEEGap);
    fChain->SetBranchAddress("ele2_isEBEtaGap", &ele2_isEBEtaGap, &b_ele2_isEBEtaGap);
    fChain->SetBranchAddress("ele2_isEBPhiGap", &ele2_isEBPhiGap, &b_ele2_isEBPhiGap);
    fChain->SetBranchAddress("ele2_isEEDeeGap", &ele2_isEEDeeGap, &b_ele2_isEEDeeGap);
    fChain->SetBranchAddress("ele2_isEERingGap", &ele2_isEERingGap, &b_ele2_isEERingGap);
    fChain->SetBranchAddress("ele2_isTrackerDriven", &ele2_isTrackerDriven, &b_ele2_isTrackerDriven);
    fChain->SetBranchAddress("ele2_sigmaIetaIeta", &ele2_sigmaIetaIeta, &b_ele2_sigmaIetaIeta);
    fChain->SetBranchAddress("ele2_DphiIn", &ele2_DphiIn, &b_ele2_DphiIn);
    fChain->SetBranchAddress("ele2_DetaIn", &ele2_DetaIn, &b_ele2_DetaIn);
    fChain->SetBranchAddress("ele2_HOverE", &ele2_HOverE, &b_ele2_HOverE);
    fChain->SetBranchAddress("ele2_tkIso", &ele2_tkIso, &b_ele2_tkIso);
    fChain->SetBranchAddress("ele2_emIso", &ele2_emIso, &b_ele2_emIso);
    fChain->SetBranchAddress("ele2_hadIso", &ele2_hadIso, &b_ele2_hadIso);
    fChain->SetBranchAddress("ele2_dxy_PV", &ele2_dxy_PV, &b_ele2_dxy_PV);
    fChain->SetBranchAddress("ele2_dz_PV", &ele2_dz_PV, &b_ele2_dz_PV);
    fChain->SetBranchAddress("ele2_sigmaP", &ele2_sigmaP, &b_ele2_sigmaP);
    fChain->SetBranchAddress("ele2_scERaw", &ele2_scERaw, &b_ele2_scERaw);
    fChain->SetBranchAddress("ele2_scEtRaw", &ele2_scEtRaw, &b_ele2_scEtRaw);
    fChain->SetBranchAddress("ele2_scE", &ele2_scE, &b_ele2_scE);
    fChain->SetBranchAddress("ele2_scEt", &ele2_scEt, &b_ele2_scEt);
    fChain->SetBranchAddress("ele2_scE_regression", &ele2_scE_regression, &b_ele2_scE_regression);
    fChain->SetBranchAddress("ele2_scEerr_regression", &ele2_scEerr_regression, &b_ele2_scEerr_regression);
    fChain->SetBranchAddress("ele2_scE_regression_PhotonTuned", &ele2_scE_regression_PhotonTuned, &b_ele2_scE_regression_PhotonTuned);
    fChain->SetBranchAddress("ele2_scEerr_regression_PhotonTuned", &ele2_scEerr_regression_PhotonTuned, &b_ele2_scEerr_regression_PhotonTuned);
    fChain->SetBranchAddress("ele2_scERaw_PUcleaned", &ele2_scERaw_PUcleaned, &b_ele2_scERaw_PUcleaned);
    fChain->SetBranchAddress("ele2_es", &ele2_es, &b_ele2_es);
    fChain->SetBranchAddress("ele2_scLaserCorr", &ele2_scLaserCorr, &b_ele2_scLaserCorr);
    fChain->SetBranchAddress("ele2_scCrackCorr", &ele2_scCrackCorr, &b_ele2_scCrackCorr);
    fChain->SetBranchAddress("ele2_scLocalContCorr", &ele2_scLocalContCorr, &b_ele2_scLocalContCorr);
    fChain->SetBranchAddress("ele2_scEta", &ele2_scEta, &b_ele2_scEta);
    fChain->SetBranchAddress("ele2_scPhi", &ele2_scPhi, &b_ele2_scPhi);
    fChain->SetBranchAddress("ele2_scLocalEta", &ele2_scLocalEta, &b_ele2_scLocalEta);
    fChain->SetBranchAddress("ele2_scLocalPhi", &ele2_scLocalPhi, &b_ele2_scLocalPhi);
    fChain->SetBranchAddress("ele2_scEtaWidth", &ele2_scEtaWidth, &b_ele2_scEtaWidth);
    fChain->SetBranchAddress("ele2_scPhiWidth", &ele2_scPhiWidth, &b_ele2_scPhiWidth);
    fChain->SetBranchAddress("ele2_scEtaWidth_PUcleaned", &ele2_scEtaWidth_PUcleaned, &b_ele2_scEtaWidth_PUcleaned);
    fChain->SetBranchAddress("ele2_scPhiWidth_PUcleaned", &ele2_scPhiWidth_PUcleaned, &b_ele2_scPhiWidth_PUcleaned);
    fChain->SetBranchAddress("ele2_fCorrection_PUcleaned", &ele2_fCorrection_PUcleaned, &b_ele2_fCorrection_PUcleaned);
    fChain->SetBranchAddress("ele2_fEta", &ele2_fEta, &b_ele2_fEta);
    fChain->SetBranchAddress("ele2_fEtaCorr", &ele2_fEtaCorr, &b_ele2_fEtaCorr);
    fChain->SetBranchAddress("ele2_tkP", &ele2_tkP, &b_ele2_tkP);
    fChain->SetBranchAddress("ele2_tkPt", &ele2_tkPt, &b_ele2_tkPt);
    fChain->SetBranchAddress("ele2_fbrem", &ele2_fbrem, &b_ele2_fbrem);
    fChain->SetBranchAddress("ele2_eSeedBC", &ele2_eSeedBC, &b_ele2_eSeedBC);
    fChain->SetBranchAddress("ele2_e5x5", &ele2_e5x5, &b_ele2_e5x5);
    fChain->SetBranchAddress("ele2_e3x3", &ele2_e3x3, &b_ele2_e3x3);
    fChain->SetBranchAddress("ele2_scNxtal", &ele2_scNxtal, &b_ele2_scNxtal);
    fChain->SetBranchAddress("ele2_bcN", &ele2_bcN, &b_ele2_bcN);
    fChain->SetBranchAddress("ele2_5x5LaserCorr", &ele2_5x5LaserCorr, &b_ele2_5x5LaserCorr);
    fChain->SetBranchAddress("ele2_3x3LaserCorr", &ele2_3x3LaserCorr, &b_ele2_3x3LaserCorr);
    fChain->SetBranchAddress("ele2_seedE", &ele2_seedE, &b_ele2_seedE);
    fChain->SetBranchAddress("ele2_seedLaserAlpha", &ele2_seedLaserAlpha, &b_ele2_seedLaserAlpha);
    fChain->SetBranchAddress("ele2_seedLaserCorr", &ele2_seedLaserCorr, &b_ele2_seedLaserCorr);
    fChain->SetBranchAddress("ele2_seedICConstant", &ele2_seedICConstant, &b_ele2_seedICConstant);
    fChain->SetBranchAddress("ele2_seedIeta", &ele2_seedIeta, &b_ele2_seedIeta);
    fChain->SetBranchAddress("ele2_seedIphi", &ele2_seedIphi, &b_ele2_seedIphi);
    fChain->SetBranchAddress("ele2_seedIx", &ele2_seedIx, &b_ele2_seedIx);
    fChain->SetBranchAddress("ele2_seedIy", &ele2_seedIy, &b_ele2_seedIy);
    fChain->SetBranchAddress("ele2_seedZside", &ele2_seedZside, &b_ele2_seedZside);
    fChain->SetBranchAddress("ele2_EOverP", &ele2_EOverP, &b_ele2_EOverP);
    fChain->SetBranchAddress("ele2_nRecHits", &ele2_nRecHits, &b_ele2_nRecHits);
    fChain->SetBranchAddress("ele2_recHit_E", &ele2_recHit_E, &b_ele2_recHit_E);
    fChain->SetBranchAddress("ele2_recHit_flag", &ele2_recHit_flag, &b_ele2_recHit_flag);
    fChain->SetBranchAddress("ele2_recHit_hashedIndex", &ele2_recHit_hashedIndex, &b_ele2_recHit_hashedIndex);
    fChain->SetBranchAddress("ele2_recHit_ietaORix", &ele2_recHit_ietaORix, &b_ele2_recHit_ietaORix);
    fChain->SetBranchAddress("ele2_recHit_iphiORiy", &ele2_recHit_iphiORiy, &b_ele2_recHit_iphiORiy);
    fChain->SetBranchAddress("ele2_recHit_zside", &ele2_recHit_zside, &b_ele2_recHit_zside);
    fChain->SetBranchAddress("ele2_recHit_laserCorrection", &ele2_recHit_laserCorrection, &b_ele2_recHit_laserCorrection);
    fChain->SetBranchAddress("ele2_recHit_Alpha", &ele2_recHit_Alpha, &b_ele2_recHit_Alpha);
    fChain->SetBranchAddress("ele2_recHit_ICConstant", &ele2_recHit_ICConstant, &b_ele2_recHit_ICConstant);
    fChain->SetBranchAddress("met_et", &met_et, &b_met_et);
    fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
    fChain->SetBranchAddress("ele1Met_mt", &ele1Met_mt, &b_ele1Met_mt);
    fChain->SetBranchAddress("ele1Met_Dphi", &ele1Met_Dphi, &b_ele1Met_Dphi);
    fChain->SetBranchAddress("ele1ele2_m", &ele1ele2_m, &b_ele1ele2_m);
    fChain->SetBranchAddress("ele1ele2_scM", &ele1ele2_scM, &b_ele1ele2_scM);
    fChain->SetBranchAddress("ele1ele2_scM_regression", &ele1ele2_scM_regression, &b_ele1ele2_scM_regression);
    fChain->SetBranchAddress("PUit_TrueNumInteractions", &PUit_TrueNumInteractions, &b_PUit_TrueNumInteractions);
    fChain->SetBranchAddress("PUit_NumInteractions", &PUit_NumInteractions, &b_PUit_NumInteractions);
    fChain->SetBranchAddress("PUit_zpositions", &PUit_zpositions, &b_PUit_zpositions);
    fChain->SetBranchAddress("PUit_sumpT_lowpT", &PUit_sumpT_lowpT, &b_PUit_sumpT_lowpT);
    fChain->SetBranchAddress("PUit_sumpT_highpT", &PUit_sumpT_highpT, &b_PUit_sumpT_highpT);
    fChain->SetBranchAddress("PUit_ntrks_lowpT", &PUit_ntrks_lowpT, &b_PUit_ntrks_lowpT);
    fChain->SetBranchAddress("PUit_ntrks_highpT", &PUit_ntrks_highpT, &b_PUit_ntrks_highpT);
    fChain->SetBranchAddress("PUoot_early_TrueNumInteractions", &PUoot_early_TrueNumInteractions, &b_PUoot_early_TrueNumInteractions);
    fChain->SetBranchAddress("PUoot_early", &PUoot_early, &b_PUoot_early);
    fChain->SetBranchAddress("PUoot_early_zpositions", &PUoot_early_zpositions, &b_PUoot_early_zpositions);
    fChain->SetBranchAddress("PUoot_early_sumpT_lowpT", &PUoot_early_sumpT_lowpT, &b_PUoot_early_sumpT_lowpT);
    fChain->SetBranchAddress("PUoot_early_sumpT_highpT", &PUoot_early_sumpT_highpT, &b_PUoot_early_sumpT_highpT);
    fChain->SetBranchAddress("PUoot_early_ntrks_lowpT", &PUoot_early_ntrks_lowpT, &b_PUoot_early_ntrks_lowpT);
    fChain->SetBranchAddress("PUoot_early_ntrks_highpT", &PUoot_early_ntrks_highpT, &b_PUoot_early_ntrks_highpT);
    fChain->SetBranchAddress("PUoot_late_TrueNumInteractions", &PUoot_late_TrueNumInteractions, &b_PUoot_late_TrueNumInteractions);
    fChain->SetBranchAddress("PUoot_late", &PUoot_late, &b_PUoot_late);
    fChain->SetBranchAddress("PUoot_late_zpositions", &PUoot_late_zpositions, &b_PUoot_late_zpositions);
    fChain->SetBranchAddress("PUoot_late_sumpT_lowpT", &PUoot_late_sumpT_lowpT, &b_PUoot_late_sumpT_lowpT);
    fChain->SetBranchAddress("PUoot_late_sumpT_highpT", &PUoot_late_sumpT_highpT, &b_PUoot_late_sumpT_highpT);
    fChain->SetBranchAddress("PUoot_late_ntrks_lowpT", &PUoot_late_ntrks_lowpT, &b_PUoot_late_ntrks_lowpT);
    fChain->SetBranchAddress("PUoot_late_ntrks_highpT", &PUoot_late_ntrks_highpT, &b_PUoot_late_ntrks_highpT);
    fChain->SetBranchAddress("mcV_E", &mcV_E, &b_mcV_E);
    fChain->SetBranchAddress("mcV_Px", &mcV_Px, &b_mcV_Px);
    fChain->SetBranchAddress("mcV_Py", &mcV_Py, &b_mcV_Py);
    fChain->SetBranchAddress("mcV_Pz", &mcV_Pz, &b_mcV_Pz);
    fChain->SetBranchAddress("mcV_Charge", &mcV_Charge, &b_mcV_Charge);
    fChain->SetBranchAddress("mcV_PdgId", &mcV_PdgId, &b_mcV_PdgId);
    fChain->SetBranchAddress("mcF1_fromV_E", &mcF1_fromV_E, &b_mcF1_fromV_E);
    fChain->SetBranchAddress("mcF1_fromV_Px", &mcF1_fromV_Px, &b_mcF1_fromV_Px);
    fChain->SetBranchAddress("mcF1_fromV_Py", &mcF1_fromV_Py, &b_mcF1_fromV_Py);
    fChain->SetBranchAddress("mcF1_fromV_Pz", &mcF1_fromV_Pz, &b_mcF1_fromV_Pz);
    fChain->SetBranchAddress("mcF1_fromV_Charge", &mcF1_fromV_Charge, &b_mcF1_fromV_Charge);
    fChain->SetBranchAddress("mcF1_fromV_PdgId", &mcF1_fromV_PdgId, &b_mcF1_fromV_PdgId);
    fChain->SetBranchAddress("mcF2_fromV_E", &mcF2_fromV_E, &b_mcF2_fromV_E);
    fChain->SetBranchAddress("mcF2_fromV_Px", &mcF2_fromV_Px, &b_mcF2_fromV_Px);
    fChain->SetBranchAddress("mcF2_fromV_Py", &mcF2_fromV_Py, &b_mcF2_fromV_Py);
    fChain->SetBranchAddress("mcF2_fromV_Pz", &mcF2_fromV_Pz, &b_mcF2_fromV_Pz);
    fChain->SetBranchAddress("mcF2_fromV_Charge", &mcF2_fromV_Charge, &b_mcF2_fromV_Charge);
    fChain->SetBranchAddress("mcF2_fromV_PdgId", &mcF2_fromV_PdgId, &b_mcF2_fromV_PdgId);
    Notify();
}

Bool_t analysis::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.
    
    return kTRUE;
}

void analysis::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t analysis::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
//#endif // #ifdef analysis_cxx



