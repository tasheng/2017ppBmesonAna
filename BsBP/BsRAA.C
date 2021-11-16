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

void BsRAA(){




	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);

	c->cd();

	c->SetLeftMargin(0.16);




	TString InfileBs = "../Bs/EffAna/FinalFiles/BsPPCorrYieldPT.root";


	TFile * FileBs = new TFile(InfileBs.Data());

	TH1D * BsCross = (TH1D *) FileBs->Get("CorrDiffHisBin");
	BsCross->SetMarkerStyle(20);
	BsCross->SetMarkerSize(1);
	BsCross->SetMarkerColor(1);
	BsCross->SetLineColor(1);



	BsCross->Draw("ep");



	//B+ PbPb//
	
	const int NBins = 4;

	float BsXsecPPY[NBins];
	float BsXsecPPX[NBins] = {8.5,12.5,17.5,35};

	float BsXSecPPYErrUp[NBins];
	float BsXSecPPYErrDown[NBins];


	float BsXSecPPYErrUpPercent[NBins];
	float BsXSecPPYErrDownPercent[NBins];


	for(int i = 0; i < NBins; i++){

		BsXsecPPY[i] = BsCross->GetBinContent(i+1);
		BsXSecPPYErrUp[i] = BsCross->GetBinError(i+1);
		BsXSecPPYErrDown[i] = BsCross->GetBinError(i+1);
		BsXSecPPYErrUpPercent[i] = BsXSecPPYErrUp[i]/BsXsecPPY[i];
		BsXSecPPYErrDownPercent[i] = BsXSecPPYErrDown[i]/BsXsecPPY[i];
		
	}



	float BsXSecPPXErrUp[NBins] = {1.5,2.5,2.5,15};
	float BsXSecPPXErrDown[NBins] = {1.5,2.5,2.5,15};
	


	float BsXsecPbPbY[NBins] = {160432,75523.7,25354.5,2272.18};
	float BsXsecPbPbX[NBins] = {8.75,12.6,17.4,27.3};


	float BsXSecPbPbXErrUp[NBins] = {1.25,2.4,2.4,22.7};
	float BsXSecPbPbXErrDown[NBins] = {1.25,2.6,2.6,7.3};

	float BsXSecPbPbYErrUpPercent[NBins] = {0.513,0.224,0.216,0.216};
	float BsXSecPbPbYErrDownPercent[NBins] = {0.483,0.256,0.207,0.163};


	float BsXSecPbPbYErrUp[NBins];
	float BsXSecPbPbYErrDown[NBins];

	for(int i = 0; i < NBins; i++){

		BsXSecPbPbYErrUp[i] = BsXSecPbPbYErrUpPercent[i] * BsXsecPbPbY[i];
		BsXSecPbPbYErrDown[i] = BsXSecPbPbYErrDownPercent[i] * BsXsecPbPbY[i];

	}











	TH2D * HisEmpty = new TH2D("HisEmpty","",100,7,50,100,100.0,250000);
	HisEmpty->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
	HisEmpty->GetYaxis()->SetTitle("Cross Section (Or Production Yield)");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
	HisEmpty->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty->Draw();

	TGraphAsymmErrors *BPPPCrossGraph = new TGraphAsymmErrors(NBins, BsXsecPPX, BsXsecPPY,BsXSecPPXErrDown, BsXSecPPXErrUp,BsXSecPPYErrDown,BsXSecPPYErrUp);
	
	TGraphAsymmErrors *BsPbPbCrossGraph = new TGraphAsymmErrors(NBins, BsXsecPbPbX, BsXsecPbPbY,BsXSecPbPbXErrDown, BsXSecPbPbXErrUp,BsXSecPbPbYErrDown,BsXSecPbPbYErrUp);

	BsPbPbCrossGraph->SetLineColor(kGreen+2);
//	BsPbPbCrossGraph->SetFillColorAlpha(kGreen-9,0.5);
	BsPbPbCrossGraph->SetMarkerStyle(20);
	BsPbPbCrossGraph->SetMarkerSize(1);
	BsPbPbCrossGraph->SetMarkerColor(kGreen+2);

	BPPPCrossGraph->SetLineColor(kBlue+2);
//	BPPPCrossGraph->SetFillColorAlpha(kBlue-9,0.5);
	BPPPCrossGraph->SetMarkerStyle(21);
	BPPPCrossGraph->SetMarkerSize(1);
	BPPPCrossGraph->SetMarkerColor(kBlue+2);


	BsPbPbCrossGraph->Draw("epsame");	
	BPPPCrossGraph->Draw("epsame");	

	TLegend* leg = new TLegend(0.30,0.60,0.60,0.80,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(BsPbPbCrossGraph,"2018 PbPb 5.02 TeV","PL");
	leg->AddEntry(BPPPCrossGraph,"2017 pp 5.02 TeV","PL");
	leg->Draw("same");


	c->SaveAs("RAA/BsPbPbPPCross.png");


	TCanvas * c2 = new TCanvas("c","c",600,600);

	c2->cd();

	c2->SetLeftMargin(0.16);


	TH2D * HisEmptyRAA = new TH2D("HisEmpty","",100,7,50,100,0,5.5);
	HisEmptyRAA->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
	HisEmptyRAA->GetYaxis()->SetTitle("RAA = #frac{1}{TAA} #frac{dN_{PbPb}/dp_{T}}{d #sigma_{pp}/d p_{T}}");
	HisEmptyRAA->GetXaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->SetTitleOffset(1.8);
	HisEmptyRAA->GetXaxis()->SetTitleOffset(1.3);
	
	HisEmptyRAA->Draw();








	float BPRAAY[NBins];
	float BPRAAX[NBins] = {8.73,12.4,17.2,27.3};


	float BPRAAXErrUp[NBins] = {1.27,2.6,2.8,22.7};
	float BPRAAXErrDown[NBins] = {1.73,2.4,2.2,7.3};

	float BPRAAYErrUp[NBins];
	float BPRAAYErrDown[NBins];


	//2015 Data Sets
	
	const int NBins2015 = 2;

	float BPRAAY2015[NBins2015] = {1.51,0.87};
	float BPRAAX2015[NBins2015] = {11.0,32.5};


	float BPRAAXErrUp2015[NBins2015] = {4.0,17.5};
	float BPRAAXErrDown2015[NBins2015] = {4.0,17.5};

	float BPRAAYErrUp2015[NBins2015] = {0.61,0.30};
	float BPRAAYErrDown2015[NBins2015] = {0.61,0.30};



	for(int i = 0; i < NBins; i++){

		BPRAAY[i] = BsXsecPbPbY[i] / BsXsecPPY[i];
	
		cout << "i = " << i <<  "    BPRAAY[i] = " << BPRAAY[i] << endl;

		BPRAAYErrUp[i] = BPRAAY[i] * TMath::Sqrt(BsXSecPbPbYErrUpPercent[i] * BsXSecPbPbYErrUpPercent[i] + BsXSecPPYErrUpPercent[i] * BsXSecPPYErrUpPercent[i]);
		BPRAAYErrDown[i] = BPRAAY[i] * TMath::Sqrt(BsXSecPbPbYErrDownPercent[i] * BsXSecPbPbYErrDownPercent[i] + BsXSecPPYErrDownPercent[i] * BsXSecPPYErrDownPercent[i]);
		

			
	}



	





	TGraphAsymmErrors *BPRAAGraph = new TGraphAsymmErrors(NBins, BPRAAX, BPRAAY,BPRAAXErrDown, BPRAAXErrUp,BPRAAYErrDown,BPRAAYErrUp);

	BPRAAGraph->SetName("BsRAAGraph");


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

	TLine * Unity = new TLine(7,1,50,1);
	Unity->SetLineWidth(2);
	Unity->SetLineStyle(2);
	Unity->SetLineColor(1);
	Unity->Draw("SAME");





	c2->SaveAs("RAA/BsRAA.png");


	TFile * fout = new TFile("OutFile/BsRAA.root","RECREATE");

	fout->cd();

	BPRAAGraph->Write();

	fout->Close();


}
