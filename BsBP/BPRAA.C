#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
//#include "tnp_weight_lowptPbPb.h"



//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void BPRAA(){



	
	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);

	c->cd();

	c->SetLeftMargin(0.16);




//	TString InfileBs = "FinalFiles/BsPPCorrYieldPT.root";
	TString InfileBP = "../BP/EffAna/FinalFiles/BPPPCorrYieldPT.root";


	TFile * FileBP = new TFile(InfileBP.Data());
	TH1D * BPCross = (TH1D *) FileBP->Get("CorrDiffHisBin");
	BPCross->SetMarkerStyle(20);
	BPCross->SetMarkerSize(1);
	BPCross->SetMarkerColor(1);
	BPCross->SetLineColor(1);



	//B+ PbPb//
	
	const int NBins = 7;

	float BPXsecPPY[NBins];
	float BPXsecPPX[NBins] = {6,8.5,12.5,17.5,25,40,55};

	float BPXSecPPYErrUp[NBins];
	float BPXSecPPYErrDown[NBins];


	float BPXSecPPYErrUpPercent[NBins];
	float BPXSecPPYErrDownPercent[NBins];


	for(int i = 0; i < NBins; i++){

		BPXsecPPY[i] = BPCross->GetBinContent(i+1);
		BPXSecPPYErrUp[i] = BPCross->GetBinError(i+1);
		BPXSecPPYErrDown[i] = BPCross->GetBinError(i+1);
		BPXSecPPYErrUpPercent[i] = BPXSecPPYErrUp[i]/BPXsecPPY[i];
		BPXSecPPYErrDownPercent[i] = BPXSecPPYErrDown[i]/BPXsecPPY[i];
		
	}



	float BPXSecPPXErrUp[NBins] = {1,1.5,2.5,2.5,5,10,5};
	float BPXSecPPXErrDown[NBins] = {1,1.5,2.5,2.5,5,10,5};
	


	float BPXsecPbPbY[NBins] = {4.82132e+06/11.1,311668,270167,64384.4,208537/11.1,28700.6/11.1,7000.73/11.1};
	float BPXsecPbPbX[NBins] = {6,8.73,12.4,17.2,25,40,55};


	float BPXSecPbPbXErrUp[NBins] = {1,1.27,2.6,2.8,5,10,5};
	float BPXSecPbPbXErrDown[NBins] = {1,1.23,2.4,2.2,5,10,5};

	float BPXSecPbPbYErrUpPercent[NBins] = {0.278198,0.159,0.041,0.0654,0.0690334,0.104543,0.24575};
	float BPXSecPbPbYErrDownPercent[NBins] = {0.278198,0.145,0.0795,0.065,0.0690334,0.104543,0.24575};


	float BPXSecPbPbYErrUp[NBins];
	float BPXSecPbPbYErrDown[NBins];

	for(int i = 0; i < NBins; i++){

		BPXSecPbPbYErrUp[i] = BPXSecPbPbYErrUpPercent[i] * BPXsecPbPbY[i];
		BPXSecPbPbYErrDown[i] = BPXSecPbPbYErrDownPercent[i] * BPXsecPbPbY[i];

	}











	TH2D * HisEmpty = new TH2D("HisEmpty","",100,5,60,100,100.0,2000000);
	HisEmpty->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmpty->GetYaxis()->SetTitle("Cross Section (Or Production Yield)");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
	HisEmpty->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty->Draw();

	TGraphAsymmErrors *BPPPCrossGraph = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY,BPXSecPPXErrDown, BPXSecPPXErrUp,BPXSecPPYErrDown,BPXSecPPYErrUp);
	
	TGraphAsymmErrors *BPPbPbCrossGraph = new TGraphAsymmErrors(NBins, BPXsecPbPbX, BPXsecPbPbY,BPXSecPbPbXErrDown, BPXSecPbPbXErrUp,BPXSecPbPbYErrDown,BPXSecPbPbYErrUp);

	BPPbPbCrossGraph->SetLineColor(kGreen+2);
//	BPPbPbCrossGraph->SetFillColorAlpha(kGreen-9,0.5);
	BPPbPbCrossGraph->SetMarkerStyle(20);
	BPPbPbCrossGraph->SetMarkerSize(1);
	BPPbPbCrossGraph->SetMarkerColor(kGreen+2);

	BPPPCrossGraph->SetLineColor(kBlue+2);
//	BPPPCrossGraph->SetFillColorAlpha(kBlue-9,0.5);
	BPPPCrossGraph->SetMarkerStyle(21);
	BPPPCrossGraph->SetMarkerSize(1);
	BPPPCrossGraph->SetMarkerColor(kBlue+2);


	BPPbPbCrossGraph->Draw("epsame");	
	BPPPCrossGraph->Draw("epsame");	

	TLegend* leg = new TLegend(0.30,0.60,0.60,0.80,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(BPPbPbCrossGraph,"2018 PbPb 5.02 TeV","PL");
	leg->AddEntry(BPPPCrossGraph,"2017 pp 5.02 TeV","PL");
	leg->Draw("same");


	c->SaveAs("RAA/BPPbPbPPCross.png");

	




	//2015 Reference//




	TH2D * HisEmpty2 = new TH2D("HisEmpty2","",100,5,60,100,100.0,3000000);
	HisEmpty2->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmpty2->GetYaxis()->SetTitle("Cross Section (Or Production Yield)");
	HisEmpty2->GetXaxis()->CenterTitle();
	HisEmpty2->GetYaxis()->CenterTitle();
	HisEmpty2->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty2->Draw();

	BPPPCrossGraph->Draw("ep");


	const int NBins2015 = 5;
	float BPXsecPPX2015[NBins2015] = {8.5,12.5,17.5,25,40};
	float BPXSecPPXErrDown2015[NBins2015] = {1.5,2.5,2.5,5,10};
	float BPXSecPPXErrUp2015[NBins2015] = {1.5,2.5,2.5,5,10};
	
	float BPXsecPPY2015[NBins2015] = {2610000,744000,197000,46500,5300};
	float BPXSecPPYErrDown2015[NBins2015] = {170000,29000,9000,2400,500};
	float BPXSecPPYErrUp2015[NBins2015] = {170000,29000,9000,2400,500};


	TGraphAsymmErrors *BPPPCrossGraph2015 = new TGraphAsymmErrors(NBins2015, BPXsecPPX2015, BPXsecPPY2015,BPXSecPPXErrDown2015, BPXSecPPXErrUp2015,BPXSecPPYErrDown2015,BPXSecPPYErrUp2015);



	BPPPCrossGraph2015->SetLineColor(kOrange+2);
	BPPPCrossGraph2015->SetMarkerStyle(25);
	BPPPCrossGraph2015->SetMarkerSize(1);
	BPPPCrossGraph2015->SetMarkerColor(kOrange+2);
	BPPPCrossGraph2015->Draw("epSAME");


	TFile * finFONLL = new TFile("FONLLInput/BPFONLL.root");
	finFONLL->cd();
	TGraphAsymmErrors *BPFONLL = (TGraphAsymmErrors*) finFONLL->Get("gaeSigmaBplus");
	BPFONLL->SetLineColor(kRed+2);
	BPFONLL->SetMarkerStyle(20);
	BPFONLL->SetMarkerSize(1);
	BPFONLL->SetMarkerColor(kRed+2);
	BPFONLL->Draw("epSAME");



	TLegend* leg3 = new TLegend(0.50,0.55,0.80,0.80,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetTextSize(0.040);
	leg3->SetTextFont(42);
	leg3->SetFillStyle(0);
	leg3->SetLineWidth(3);
	leg3->AddEntry(BPPPCrossGraph,"2017 pp 5.02 TeV","PL");
	leg3->AddEntry(BPPPCrossGraph2015,"2015 pp 5.02 TeV","PL");
	leg3->AddEntry(BPFONLL,"FONLL Calculations","PL");
	leg3->Draw("same");



	c->SaveAs("RAA/BPCrossComp.png");
	c->SetLogy();

	c->SaveAs("RAA/BPCrossCompLog.png");



	//FONLL












	TCanvas * c2 = new TCanvas("c","c",600,600);

	c2->cd();

	c2->SetLeftMargin(0.16);



	TH2D * HisEmptyRAA = new TH2D("HisEmpty","",100,5,60,100,0,1.3);
	HisEmptyRAA->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmptyRAA->GetYaxis()->SetTitle("RAA = #frac{1}{TAA} #frac{dN_{PbPb}/dp_{T}}{d #sigma_{pp}/d p_{T}}");
	HisEmptyRAA->GetXaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->SetTitleOffset(1.8);
	
	HisEmptyRAA->Draw();

	


	float BPRAAY[NBins];
//	float BPRAAX[NBins] = {8.73,12.4,17.2,27.3};
	float BPRAAX[NBins] = {6,8.73,12.4,17.2,25,40,55};


//	float BPRAAXErrUp[NBins] = {1.27,2.6,2.8,22.7};
//	float BPRAAXErrDown[NBins] = {1.73,2.4,2.2,7.3};

	float BPRAAXErrUp[NBins] = {1,1.27,2.6,2.8,5,10,5};
	float BPRAAXErrDown[NBins] = {1,1.23,2.4,2.2,5,10,5};



	float BPRAAYErrUp[NBins];
	float BPRAAYErrDown[NBins];


	//2015 Data Sets
	
	//const int NBins2015 = 5;

	float BPRAAY2015[NBins2015] = {0.35,0.448,0.440,0.615,0.35};
	float BPRAAX2015[NBins2015] = {8.5,12.5,17.5,25,40};


	float BPRAAXErrUp2015[NBins2015] = {1.5,2.5,2.5,5,10};
	float BPRAAXErrDown2015[NBins2015] = {1.5,2.5,2.5,5,10};

	float BPRAAYErrUp2015[NBins2015] = {0.11,0.074,0.075,0.092,0.11};
	float BPRAAYErrDown2015[NBins2015] = {0.11,0.074,0.075,0.092,0.11};



	for(int i = 0; i < NBins; i++){

		BPRAAY[i] = BPXsecPbPbY[i] / BPXsecPPY[i];
	
		cout << "i = " << i <<  "    BPRAAY[i] = " << BPRAAY[i] << endl;

		BPRAAYErrUp[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYErrUpPercent[i] * BPXSecPbPbYErrUpPercent[i] + BPXSecPPYErrUpPercent[i] * BPXSecPPYErrUpPercent[i]);
		BPRAAYErrDown[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYErrDownPercent[i] * BPXSecPbPbYErrDownPercent[i] + BPXSecPPYErrDownPercent[i] * BPXSecPPYErrDownPercent[i]);
		

			
	}



	





	TGraphAsymmErrors *BPRAAGraph = new TGraphAsymmErrors(NBins, BPRAAX, BPRAAY,BPRAAXErrDown, BPRAAXErrUp,BPRAAYErrDown,BPRAAYErrUp);
	BPRAAGraph->SetName("BPRAAGraph");


	TGraphAsymmErrors *BPRAAGraph2015 = new TGraphAsymmErrors(NBins2015, BPRAAX2015, BPRAAY2015,BPRAAXErrDown2015, BPRAAXErrUp2015,BPRAAYErrDown2015,BPRAAYErrUp2015);



	BPRAAGraph->SetLineColor(kRed+2);
//	BPRAAGraph->SetFillColorAlpha(kRed+2,0.5);
	BPRAAGraph->SetMarkerStyle(20);
	BPRAAGraph->SetMarkerSize(1);

	BPRAAGraph->SetMarkerColor(kRed+2);

	BPRAAGraph2015->SetLineColor(kBlue+2);
//	BPRAAGraph->SetFillColorAlpha(kRed+2,0.5);
	BPRAAGraph2015->SetMarkerStyle(21);
	BPRAAGraph2015->SetMarkerSize(1);
	BPRAAGraph2015->SetMarkerColor(kBlue+2);
	

	

	BPRAAGraph->Draw("epSAME");
	BPRAAGraph2015->Draw("epSAME");


	TLegend* leg2 = new TLegend(0.30,0.75,0.60,0.90,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.040);
	leg2->SetTextFont(42);
	leg2->SetFillStyle(0);
	leg2->SetLineWidth(3);
	leg2->AddEntry(BPRAAGraph,"2018 PbPb + 2017 pp","PL");
	leg2->AddEntry(BPRAAGraph2015,"2015 PbPb + 2015 pp","PL");
	leg2->Draw("same");

	TLine * Unity = new TLine(5,1,60,1);
	Unity->SetLineWidth(2);
	Unity->SetLineStyle(2);
	Unity->SetLineColor(1);
	Unity->Draw("SAME");





	c2->SaveAs("RAA/BPRAA.png");
	

	TFile * fout = new TFile("OutFile/BPRAA.root","RECREATE");

	fout->cd();

	BPRAAGraph->Write();

	fout->Close();







}
