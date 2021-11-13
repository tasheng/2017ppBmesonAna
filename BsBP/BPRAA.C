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




	TString InfileBs = "../Bs/EffAna/FinalFiles/BsPPCorrYieldPT.root";
	TString InfileBP = "../BP/EffAna/FinalFiles/BPPPCorrYieldPT.root";

	TFile * FileBs = new TFile(InfileBs.Data());

	TH1D * BsCross = (TH1D *) FileBs->Get("CorrDiffHisBin");
	BsCross->SetMarkerStyle(20);
	BsCross->SetMarkerSize(1);
	BsCross->SetMarkerColor(1);
	BsCross->SetLineColor(1);



	BsCross->Draw("ep");

	c->SaveAs("Mult/BsCross.png");

	TFile * FileBP = new TFile(InfileBP.Data());
	TH1D * BPCross = (TH1D *) FileBP->Get("CorrDiffHisBin");
	BPCross->SetMarkerStyle(20);
	BPCross->SetMarkerSize(1);
	BPCross->SetMarkerColor(1);
	BPCross->SetLineColor(1);



	//B+ PbPb//
	
	const int NBins = 4;

	float BPXsecPPY[NBins];
	float BPXsecPPX[NBins] = {8.5,12.5,17.5,35};

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



	float BPXSecPPXErrUp[NBins] = {1.5,2.5,2.5,15};
	float BPXSecPPXErrDown[NBins] = {1.5,2.5,2.5,15};
	


	float BPXsecPbPbY[NBins] = {311668,270167,64384.4,7704.11};
	float BPXsecPbPbX[NBins] = {8.73,12.4,17.2,27.3};


	float BPXSecPbPbXErrUp[NBins] = {1.27,2.6,2.8,22.7};
	float BPXSecPbPbXErrDown[NBins] = {1.23,2.4,2.2,7.3};

	float BPXSecPbPbYErrUpPercent[NBins] = {0.159,0.041,0.0654,0.069};
	float BPXSecPbPbYErrDownPercent[NBins] = {0.145,0.0795,0.065,0.0526};


	float BPXSecPbPbYErrUp[NBins];
	float BPXSecPbPbYErrDown[NBins];

	for(int i = 0; i < NBins; i++){

		BPXSecPbPbYErrUp[i] = BPXSecPbPbYErrUpPercent[i] * BPXsecPbPbY[i];
		BPXSecPbPbYErrDown[i] = BPXSecPbPbYErrDownPercent[i] * BPXsecPbPbY[i];

	}











	TH2D * HisEmpty = new TH2D("HisEmpty","",100,7,50,100,100.0,2000000);
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


	c->SaveAs("BPRAA/BPPbPbPPCross.png");


	TCanvas * c2 = new TCanvas("c","c",600,600);

	c2->cd();

	c2->SetLeftMargin(0.16);


	TH2D * HisEmptyRAA = new TH2D("HisEmpty","",100,7,50,100,0,1);
	HisEmptyRAA->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmptyRAA->GetYaxis()->SetTitle("RAA = #frac{1}{TAA} #frac{dN_{PbPb}/dp_{T}}{d #sigma_{pp}/d p_{T}}");
	HisEmptyRAA->GetXaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->SetTitleOffset(1.8);
	
	HisEmptyRAA->Draw();








	float BPRAAY[NBins];
	float BPRAAX[NBins] = {8.73,12.4,17.2,27.3};


	float BPRAAXErrUp[NBins] = {1.27,2.6,2.8,22.7};
	float BPRAAXErrDown[NBins] = {1.73,2.4,2.2,7.3};

	float BPRAAYErrUp[NBins];
	float BPRAAYErrDown[NBins];


	//2015 Data Sets
	
	const int NBins2015 = 5;

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


	TLegend* leg2 = new TLegend(0.30,0.70,0.60,0.90,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.040);
	leg2->SetTextFont(42);
	leg2->SetFillStyle(0);
	leg2->SetLineWidth(3);
	leg2->AddEntry(BPRAAGraph,"2018 PbPb + 2017 pp","PL");
	leg2->AddEntry(BPRAAGraph2015,"2015 PbPb + 2015 pp","PL");
	leg2->Draw("same");




	c2->SaveAs("BPRAA/BPRAA.png");
	






}
