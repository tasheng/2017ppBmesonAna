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
	TString InfileBP = "../../BP/EffAna/FinalFiles/BPPPCorrYieldPT.root";


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



	//Syst//



	float BPXSecPPYSystUp[NBins];
	float BPXSecPPYSystDown[NBins];



	float BPTrackingSyst[NBins] = {0.05,0.05,0.05,0.05,0.05,0.05,0.05};

	float BPMDDataSyst[NBins] = {0.208,0.00947,0.00458,0.0181,0.0490,0.0459,0.0123};
	float BPPDFSyst[NBins] = {0.0212,0.0122,0.0091,0.0411,0.0120,0.0150,0.0089};
	float BPPtShapeSyst[NBins] = {0.0096,0.00814,0.000408,0.000617,0.000833,0.00138,0.000640};

	float BPTnPSystDown[NBins] = {0.00454,0.00476,0.00382,0.00354,0.00429,0.00581,0.00613};
	float BPTnPSystUp[NBins] = {0.00454,0.00476,0.00382,0.00354,0.00429,0.00581,0.00613};

	float BPTotalSystDown[NBins];
	float BPTotalSystUp[NBins];

	for(int i = 0; i < NBins; i++){

		
	//	BPTotalSystDown[i] = TMath::Sqrt(BPMDDataSyst[i] * BPMDDataSyst[i] + BPPDFSyst[i] * BPPDFSyst[i] + BPPtShapeSyst[i] * BPPtShapeSyst[i] + BPTnPSystDown[i] * BPTnPSystDown[i]);
	//	BPTotalSystUp[i] = TMath::Sqrt(BPMDDataSyst[i] * BPMDDataSyst[i] + BPPDFSyst[i] * BPPDFSyst[i] + BPPtShapeSyst[i] * BPPtShapeSyst[i] + BPTnPSystUp[i] * BPTnPSystUp[i]);
		BPTotalSystDown[i] = TMath::Sqrt(BPTrackingSyst[i] * BPTrackingSyst[i] + BPMDDataSyst[i] * BPMDDataSyst[i] + BPPDFSyst[i] * BPPDFSyst[i] + BPPtShapeSyst[i] * BPPtShapeSyst[i] + BPTnPSystDown[i] * BPTnPSystDown[i]);
		BPTotalSystUp[i] = TMath::Sqrt(BPTrackingSyst[i] * BPTrackingSyst[i] + BPMDDataSyst[i] * BPMDDataSyst[i] + BPPDFSyst[i] * BPPDFSyst[i] + BPPtShapeSyst[i] * BPPtShapeSyst[i] + BPTnPSystUp[i] * BPTnPSystUp[i]);



	}
	
	for(int i = 0; i < NBins; i++){

		BPXSecPPYSystUp[i] = BPXsecPPY[i] * ( BPTotalSystUp[i]);
		BPXSecPPYSystDown[i] = BPXsecPPY[i] * (BPTotalSystDown[i] );
		//cout << "i = " << i << "     BPXSecPPYSystDown[i] = " << BPXSecPPYSystDown[i] << "  BPXSecPPYErrDown[i] =  " << BPXSecPPYErrDown[i] << endl;
		
	}


	//PbPb

	float BPXSecPbPbYSystUpPercent[NBins] = {0.3577,0.1404,0.1714,0.0775,0.0858,0.0715,0.1253};
	float BPXSecPbPbYSystDownPercent[NBins] = {0.3210,0.1359,0.1705,0.0761,0.0843,0.0699,0.1220};


	float BPXSecPbPbYSystUp[NBins];
	float BPXSecPbPbYSystDown[NBins];


	for(int i = 0; i < NBins; i++){

		BPXSecPbPbYSystDown[i] = (BPXSecPbPbYSystDownPercent[i]) * BPXsecPbPbY[i];
		BPXSecPbPbYSystUp[i] = (BPXSecPbPbYSystUpPercent[i]) * BPXsecPbPbY[i];

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


	TGraphAsymmErrors *BPPPCrossGraphSyst  = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY, BPXSecPPXErrDown, BPXSecPPXErrUp, BPXSecPPYSystDown,BPXSecPPYSystUp);
  	TGraphAsymmErrors *BPPbPbCrossGraphSyst    = new TGraphAsymmErrors(NBins, BPXsecPbPbX, BPXsecPbPbY, BPXSecPbPbXErrDown, BPXSecPbPbXErrUp, BPXSecPbPbYSystDown,BPXSecPbPbYSystUp);


	BPPPCrossGraphSyst->SetFillColorAlpha(kBlue-9,0.5);
	BPPPCrossGraphSyst->SetLineColor(kBlue-9);
	BPPbPbCrossGraphSyst->SetFillColorAlpha(kGreen-9,0.5);
	BPPbPbCrossGraphSyst->SetLineColor(kGreen-9);



	BPPbPbCrossGraph->Draw("epsame");	
	BPPPCrossGraph->Draw("epsame");	

	BPPPCrossGraphSyst->Draw("5same");	
	BPPbPbCrossGraphSyst->Draw("5same");	
	

	TLegend* leg = new TLegend(0.45,0.65,0.75,0.85,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(BPPbPbCrossGraph,"2018 PbPb 5.02 TeV","PL");
	leg->AddEntry(BPPPCrossGraph,"2017 pp 5.02 TeV","PL");
	leg->Draw("same");


	c->SaveAs("RAAPlots/BP/BPPbPbPPCross.png");

	c->SetLogy();

	c->SaveAs("RAAPlots/BP/BPPbPbPPCrossLog.png");






	//Fid and Non-Fid

	
	//Fid Non Fid

	TH2D * HisEmpty6 = new TH2D("HisEmpty6","",100,5,60,100,100.0,8000000);
	HisEmpty6->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmpty6->GetYaxis()->SetTitle("Cross Section (Or Production Yield)");
	HisEmpty6->GetXaxis()->CenterTitle();
	HisEmpty6->GetYaxis()->CenterTitle();
	HisEmpty6->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty6->Draw();

	

	TString InfileBP2 = "../Comparisons/NoFiducial/FinalFiles/BPPPCorrYieldPT.root";


	TFile * FileBP2 = new TFile(InfileBP2.Data());
	TH1D * BPCross2 = (TH1D *) FileBP2->Get("CorrDiffHisBin");
	BPCross2->SetMarkerStyle(20);
	BPCross2->SetMarkerSize(1);
	BPCross2->SetMarkerColor(1);
	BPCross2->SetLineColor(1);


	


	float BPXsecPPY2[NBins];


	float BPXSecPPYErrUpPercent2[NBins];
	float BPXSecPPYErrDownPercent2[NBins];
	float BPXSecPPYErrUp2[NBins];
	float BPXSecPPYErrDown2[NBins];

	for(int i = 0; i < NBins; i++){

		BPXsecPPY2[i] = BPCross2->GetBinContent(i+1);
		BPXSecPPYErrUp2[i] = BPCross2->GetBinError(i+1);
		BPXSecPPYErrDown2[i] = BPCross2->GetBinError(i+1);
		BPXSecPPYErrUpPercent2[i] = BPXSecPPYErrUp2[i]/BPXsecPPY2[i];
		BPXSecPPYErrDownPercent2[i] = BPXSecPPYErrDown2[i]/BPXsecPPY2[i];
		
	}

	TGraphAsymmErrors *BPPPCrossGraph2 = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY2,BPXSecPPXErrDown, BPXSecPPXErrUp,BPXSecPPYErrDown2,BPXSecPPYErrUp2);


	BPPPCrossGraph2->SetLineColor(kGreen+2);
//	BPPPCrossGraph->SetFillColorAlpha(kBlue-9,0.5);
	BPPPCrossGraph2->SetMarkerStyle(20);
	BPPPCrossGraph2->SetMarkerSize(1);
	BPPPCrossGraph2->SetMarkerColor(kGreen+2);

	BPPPCrossGraph2->Draw("ep");
	BPPPCrossGraph->Draw("epSAME");	

	TLegend* leg5 = new TLegend(0.45,0.55,0.75,0.80,NULL,"brNDC");
	leg5->SetBorderSize(0);
	leg5->SetTextSize(0.040);
	leg5->SetTextFont(42);
	leg5->SetFillStyle(0);
	leg5->SetLineWidth(3);
	leg5->AddEntry(BPPPCrossGraph,"Fiducial Region Used","PL");
	leg5->AddEntry(BPPPCrossGraph2,"No Fiducial Region","PL");
	leg5->Draw("same");





    c->SaveAs("RAAPlots/BP/BPFidOrNotComp.png");





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


	float BPRAAYSystUp[NBins] ;
	float BPRAAYSystDown[NBins];



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

	float BPRAAYSystUp2015[NBins2015] = {0.064,0.077,0.074,0.102,0.0539};
	float BPRAAYSystDown2015[NBins2015] = {0.064,0.077,0.074,0.102,0.0539};



	for(int i = 0; i < NBins; i++){

		BPRAAY[i] = BPXsecPbPbY[i] / BPXsecPPY[i];


		BPRAAYErrUp[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYErrUpPercent[i] * BPXSecPbPbYErrUpPercent[i] + BPXSecPPYErrUpPercent[i] * BPXSecPPYErrUpPercent[i]);
		BPRAAYErrDown[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYErrDownPercent[i] * BPXSecPbPbYErrDownPercent[i] + BPXSecPPYErrDownPercent[i] * BPXSecPPYErrDownPercent[i]);
		

		BPRAAYSystDown[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYSystDownPercent[i] * BPXSecPbPbYSystDownPercent[i] + BPTotalSystDown[i] * BPTotalSystDown[i]);
		BPRAAYSystUp[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYSystUpPercent[i] * BPXSecPbPbYSystUpPercent[i] + BPTotalSystUp[i] * BPTotalSystUp[i]);

		cout << "i = " << i <<  "    BPRAAY[i] = " << BPRAAY[i] << "   BPRAAYErrUp[i]  = "  << BPRAAYErrUp[i]/BPRAAY[i]   << "  BPRAAYErrDown[i] =   " << BPRAAYErrDown[i]/BPRAAY[i]  <<   "   BPRAAYSystUp[i]  = "  << BPRAAYSystUp[i]/BPRAAY[i]   << "  BPRAAYSystDown[i] =   " << BPRAAYSystDown[i]/BPRAAY[i]  << endl;

			
	}






	TGraphAsymmErrors *BPRAAGraph = new TGraphAsymmErrors(NBins, BPRAAX, BPRAAY,BPRAAXErrDown, BPRAAXErrUp,BPRAAYErrDown,BPRAAYErrUp);
	BPRAAGraph->SetName("BPRAAGraph");
	TGraphAsymmErrors *BPRAAGraphSyst = new TGraphAsymmErrors(NBins, BPRAAX, BPRAAY,BPRAAXErrDown, BPRAAXErrUp,BPRAAYSystDown,BPRAAYSystUp);



	TGraphAsymmErrors *BPRAAGraph2015 = new TGraphAsymmErrors(NBins2015, BPRAAX2015, BPRAAY2015,BPRAAXErrDown2015, BPRAAXErrUp2015,BPRAAYErrDown2015,BPRAAYErrUp2015);
	TGraphAsymmErrors *BPRAAGraphSyst2015 = new TGraphAsymmErrors(NBins2015, BPRAAX2015, BPRAAY2015,BPRAAXErrDown2015, BPRAAXErrUp2015,BPRAAYSystDown2015,BPRAAYSystUp2015);


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
	

	BPRAAGraphSyst->SetFillColorAlpha(kRed-9,0.5);
	BPRAAGraphSyst->SetLineColor(kRed-9);
	BPRAAGraphSyst2015->SetFillColorAlpha(kBlue-9,0.5);
	BPRAAGraphSyst2015->SetLineColor(kBlue-9);





	TLine * Unity = new TLine(5,1,60,1);
	Unity->SetLineWidth(2);
	Unity->SetLineStyle(2);
	Unity->SetLineColor(1);

	HisEmptyRAA->Draw();
	BPRAAGraph->Draw("ep");
	BPRAAGraphSyst->Draw("5same");
	Unity->Draw("SAME");
	
	c2->SaveAs("RAAPlots/BP/BPRAA.png");
	c2->SaveAs("RAAPlots/BP/BPRAA.pdf");


	HisEmptyRAA->Draw();
	BPRAAGraph->Draw("ep");
	BPRAAGraph->Draw("epSAME");
	BPRAAGraph2015->Draw("epSAME");
	BPRAAGraphSyst->Draw("5same");
	BPRAAGraphSyst2015->Draw("5same");


	TLegend* leg2 = new TLegend(0.30,0.75,0.60,0.90,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.040);
	leg2->SetTextFont(42);
	leg2->SetFillStyle(0);
	leg2->SetLineWidth(3);
	leg2->AddEntry(BPRAAGraph,"2018 PbPb + 2017 pp","PL");
	leg2->AddEntry(BPRAAGraph2015,"2015 PbPb + 2015 pp","PL");
	leg2->Draw("same");



	Unity->Draw("SAME");




	c2->SaveAs("RAAPlots/BP/BPRAACompairson.png");
	c2->SaveAs("RAAPlots/BP/BPRAACompairson.pdf");
	

	TFile * fout = new TFile("OutFile/BPRAA.root","RECREATE");

	fout->cd();

	BPRAAGraph->Write();

	fout->Close();







}
