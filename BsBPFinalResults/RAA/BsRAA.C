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




	TString InfileBs = "../../Bs/EffAna/FinalFiles/BsPPCorrYieldPT.root";
	//TString InfileBs = "BsPPCorrYieldPT.root";


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




	//Syst//



	float BsXSecPPYSystUp[NBins];
	float BsXSecPPYSystDown[NBins];




	float BsTrackingSyst[NBins] = {0.10,0.10,0.10,0.10};
	float BsMDDataSyst[NBins] = {0.0466,0.0870,0.0725,0.0495};
	float BsPtShapeSyst[NBins] = {0.0194,0.00607,0.00100,0.000312};	
	float BsPDFSyst[NBins] = {0.0545,0.012,0.0213,0.0395};

	float BsTnPSystDown[NBins] = {0.00412,0.00385,0.00376,0.00457};
	float BsTnPSystUp[NBins] = {0.00412,0.00385,0.00376,0.00457};

	float BsTotalSystDown[NBins];
	float BsTotalSystUp[NBins];

	for(int i = 0; i < NBins; i++){

		
		BsTotalSystDown[i] = TMath::Sqrt(BsTrackingSyst[i] * BsTrackingSyst[i] + BsMDDataSyst[i] * BsMDDataSyst[i] + BsPDFSyst[i] * BsPDFSyst[i] + BsPtShapeSyst[i] * BsPtShapeSyst[i] + BsTnPSystDown[i] * BsTnPSystDown[i]);
		BsTotalSystUp[i] = TMath::Sqrt(BsTrackingSyst[i] * BsTrackingSyst[i] +BsMDDataSyst[i] * BsMDDataSyst[i] + BsPDFSyst[i] * BsPDFSyst[i] + BsPtShapeSyst[i] * BsPtShapeSyst[i] + BsTnPSystUp[i] * BsTnPSystUp[i]);


	}

	
	for(int i = 0; i < NBins; i++){

		BsXSecPPYSystUp[i] = BsXsecPPY[i] * ( BsTotalSystUp[i]);
		BsXSecPPYSystDown[i] = BsXsecPPY[i] * (BsTotalSystDown[i] );
		//cout << "i = " << i << "     BsXSecPPYSystDown[i] = " << BsXSecPPYSystDown[i] << "  BsXSecPPYErrDown[i] =  " << BsXSecPPYErrDown[i] << endl;
		
	}



	//PbPb

	float BsXSecPbPbYSystUpPercent[NBins] = {0.4564,0.1482,0.1218,0.1647};
	float BsXSecPbPbYSystDownPercent[NBins] = {0.4564,0.1454,0.1210,0.1640};


	float BsXSecPbPbYSystUp[NBins];
	float BsXSecPbPbYSystDown[NBins];


	for(int i = 0; i < NBins; i++){

		BsXSecPbPbYSystDown[i] = (BsXSecPbPbYSystDownPercent[i]) * BsXsecPbPbY[i];
		BsXSecPbPbYSystUp[i] = (BsXSecPbPbYSystUpPercent[i]) * BsXsecPbPbY[i];

	}








	TH2D * HisEmpty = new TH2D("HisEmpty","",100,7,50,100,200.0,350000);
	HisEmpty->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
	HisEmpty->GetYaxis()->SetTitle("Cross Section (Or Production Yield)");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
	HisEmpty->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty->Draw();

	TGraphAsymmErrors *BsPPCrossGraph = new TGraphAsymmErrors(NBins, BsXsecPPX, BsXsecPPY,BsXSecPPXErrDown, BsXSecPPXErrUp,BsXSecPPYErrDown,BsXSecPPYErrUp);
	
	TGraphAsymmErrors *BPPbPbCrossGraph = new TGraphAsymmErrors(NBins, BsXsecPbPbX, BsXsecPbPbY,BsXSecPbPbXErrDown, BsXSecPbPbXErrUp,BsXSecPbPbYErrDown,BsXSecPbPbYErrUp);

	BPPbPbCrossGraph->SetLineColor(kGreen+2);
//	BPPbPbCrossGraph->SetFillColorAlpha(kGreen-9,0.5);
	BPPbPbCrossGraph->SetMarkerStyle(20);
	BPPbPbCrossGraph->SetMarkerSize(1);
	BPPbPbCrossGraph->SetMarkerColor(kGreen+2);

	BsPPCrossGraph->SetLineColor(kBlue+2);
//	BsPPCrossGraph->SetFillColorAlpha(kBlue-9,0.5);
	BsPPCrossGraph->SetMarkerStyle(21);
	BsPPCrossGraph->SetMarkerSize(1);
	BsPPCrossGraph->SetMarkerColor(kBlue+2);


	TGraphAsymmErrors *BPPPCrossGraphSyst  = new TGraphAsymmErrors(NBins, BsXsecPPX, BsXsecPPY, BsXSecPPXErrDown, BsXSecPPXErrUp, BsXSecPPYSystDown,BsXSecPPYSystUp);
  	TGraphAsymmErrors *BPPbPbCrossGraphSyst    = new TGraphAsymmErrors(NBins, BsXsecPbPbX, BsXsecPbPbY, BsXSecPbPbXErrDown, BsXSecPbPbXErrUp, BsXSecPbPbYSystDown,BsXSecPbPbYSystUp);


	BPPPCrossGraphSyst->SetFillColorAlpha(kBlue-9,0.5);
	BPPPCrossGraphSyst->SetLineColor(kBlue-9);
	BPPbPbCrossGraphSyst->SetFillColorAlpha(kGreen-9,0.5);
	BPPbPbCrossGraphSyst->SetLineColor(kGreen-9);




	BPPbPbCrossGraph->Draw("epsame");	
	BsPPCrossGraph->Draw("epsame");	


	BPPPCrossGraphSyst->Draw("5same");	
	BPPbPbCrossGraphSyst->Draw("5same");	

	TLegend* leg = new TLegend(0.45,0.65,0.75,0.85,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(BPPbPbCrossGraph,"2018 PbPb 5.02 TeV","PL");
	leg->AddEntry(BsPPCrossGraph,"2017 pp 5.02 TeV","PL");
	leg->Draw("same");


	c->SaveAs("RAAPlots/Bs/BsPbPbPPCross.png");


	c->SetLogy();

	c->SaveAs("RAAPlots/Bs/BsPbPbPPCrossLog.png");


	//2015 References



	TH2D * HisEmpty2 = new TH2D("HisEmpty2","",100,5,60,100,100.0,3000000);
	HisEmpty2->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
	HisEmpty2->GetYaxis()->SetTitle("Cross Section (Or Production Yield)");
	HisEmpty2->GetXaxis()->CenterTitle();
	HisEmpty2->GetYaxis()->CenterTitle();
	HisEmpty2->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty2->Draw();

	BsPPCrossGraph->Draw("ep");


	const int NBins2015PP = 3;
	float BsXsecPPX2015[NBins2015PP] = {11,17.5,35.0};
	float BsXSecPPXErrDown2015[NBins2015PP] = {4,2.5,15};
	float BsXSecPPXErrUp2015[NBins2015PP] = {4,2.5,15};
	
	float BsXsecPPY2015[NBins2015PP] = {316000,34100,3830};
	float BsXSecPPYErrDown2015[NBins2015PP] = {37000,6300,670};
	float BsXSecPPYErrUp2015[NBins2015PP] = {37000,6300,670};


	TGraphAsymmErrors *BsPPCrossGraph2015 = new TGraphAsymmErrors(NBins2015PP, BsXsecPPX2015, BsXsecPPY2015,BsXSecPPXErrDown2015, BsXSecPPXErrUp2015,BsXSecPPYErrDown2015,BsXSecPPYErrUp2015);



	BsPPCrossGraph2015->SetLineColor(kOrange+2);
	BsPPCrossGraph2015->SetMarkerStyle(25);
	BsPPCrossGraph2015->SetMarkerSize(1);
	BsPPCrossGraph2015->SetMarkerColor(kOrange+2);
	BsPPCrossGraph2015->Draw("epSAME");

/*
	TFile * finFONLL = new TFile("../BsFONLL.root");
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
	leg3->AddEntry(BsPPCrossGraph,"2017 pp 5.02 TeV","PL");
	leg3->AddEntry(BsPPCrossGraph2015,"2015 pp 5.02 TeV","PL");
	leg3->AddEntry(BPFONLL,"FONLL Calculations","PL");
	leg3->Draw("same");


*/
	c->SaveAs("RAAPlots/Bs/BsCrossComp.png");
	c->SetLogy();

	c->SaveAs("RAAPlots/Bs/BsCrossCompLog.png");

/*


	//Introduce Non-Fiducial
	//
	//
	TH2D * HisEmpty6 = new TH2D("HisEmpty6","",100,5,60,100,100.0,8000000);
	HisEmpty6->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmpty6->GetYaxis()->SetTitle("Cross Section (Or Production Yield)");
	HisEmpty6->GetXaxis()->CenterTitle();
	HisEmpty6->GetYaxis()->CenterTitle();
	HisEmpty6->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty6->Draw();

	

	TString InfileBP2 = "FinalFiles/BsPPCorrYieldPT2.root";


	TFile * FileBP2 = new TFile(InfileBP2.Data());
	TH1D * BsCross2 = (TH1D *) FileBP2->Get("CorrDiffHisBin");
	BsCross2->SetMarkerStyle(20);
	BsCross2->SetMarkerSize(1);
	BsCross2->SetMarkerColor(1);
	BsCross2->SetLineColor(1);


	


	float BsXsecPPY2[NBins];


	float BsXSecPPYErrUpPercent2[NBins];
	float BsXSecPPYErrDownPercent2[NBins];
	float BsXSecPPYErrUp2[NBins];
	float BsXSecPPYErrDown2[NBins];

	for(int i = 0; i < NBins; i++){

		BsXsecPPY2[i] = BsCross2->GetBinContent(i+1);
		BsXSecPPYErrUp2[i] = BsCross2->GetBinError(i+1);
		BsXSecPPYErrDown2[i] = BsCross2->GetBinError(i+1);
		BsXSecPPYErrUpPercent2[i] = BsXSecPPYErrUp2[i]/BsXsecPPY2[i];
		BsXSecPPYErrDownPercent2[i] = BsXSecPPYErrDown2[i]/BsXsecPPY2[i];
		
	}

	TGraphAsymmErrors *BsPPCrossGraph2 = new TGraphAsymmErrors(NBins, BsXsecPPX, BsXsecPPY2,BsXSecPPXErrDown, BsXSecPPXErrUp,BsXSecPPYErrDown2,BsXSecPPYErrUp2);


	BsPPCrossGraph2->SetLineColor(kGreen+2);
//	BsPPCrossGraph->SetFillColorAlpha(kBlue-9,0.5);
	BsPPCrossGraph2->SetMarkerStyle(20);
	BsPPCrossGraph2->SetMarkerSize(1);
	BsPPCrossGraph2->SetMarkerColor(kGreen+2);

	BsPPCrossGraph2->Draw("ep");
	BsPPCrossGraph->Draw("epSAME");	

	TLegend* leg5 = new TLegend(0.45,0.55,0.75,0.80,NULL,"brNDC");
	leg5->SetBorderSize(0);
	leg5->SetTextSize(0.040);
	leg5->SetTextFont(42);
	leg5->SetFillStyle(0);
	leg5->SetLineWidth(3);
	leg5->AddEntry(BsPPCrossGraph,"Fiducial Region Used","PL");
	leg5->AddEntry(BsPPCrossGraph2,"No Fiducial Region","PL");
	leg5->Draw("same");





	c->SaveAs("RAAPlots/Bs/BsCrossLog.png");

*/
//DONE


	TCanvas * c2 = new TCanvas("c","c",600,600);

	c2->cd();

	c2->SetLeftMargin(0.16);


	TLine * Unity = new TLine(7,1,50,1);
	Unity->SetLineWidth(2);
	Unity->SetLineStyle(2);
	Unity->SetLineColor(1);



	TH2D * HisEmptyRAA = new TH2D("HisEmpty","",100,7,50,100,0,5.5);
	HisEmptyRAA->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
	HisEmptyRAA->GetYaxis()->SetTitle("RAA = #frac{1}{TAA} #frac{dN_{PbPb}/dp_{T}}{d #sigma_{pp}/d p_{T}}");
	HisEmptyRAA->GetXaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->SetTitleOffset(1.8);
	HisEmptyRAA->GetXaxis()->SetTitleOffset(1.3);
	







	float BsRAAY[NBins];
	float BsRAAX[NBins] = {8.73,12.4,17.2,27.3};


	float BsRAAXErrUp[NBins] = {1.27,2.6,2.8,22.7};
	float BsRAAXErrDown[NBins] = {1.73,2.4,2.2,7.3};

	float BsRAAYErrUp[NBins];
	float BsRAAYErrDown[NBins];


	//2015 Data Sets
	
	const int NBins2015 = 2;

	float BsRAAY2015[NBins2015] = {1.51,0.87};
	float BsRAAX2015[NBins2015] = {11.0,32.5};


	float BsRAAXErrUp2015[NBins2015] = {4.0,17.5};
	float BsRAAXErrDown2015[NBins2015] = {4.0,17.5};

	float BsRAAYErrUp2015[NBins2015] = {0.61,0.30};
	float BsRAAYErrDown2015[NBins2015] = {0.61,0.30};


	float BsRAAYSystUp2015[NBins2015] = {0.50,0.17};
	float BsRAAYSystDown2015[NBins2015] = {0.50,0.17};


	float BsRAAYSystUp[NBins] ;
	float BsRAAYSystDown[NBins];


	for(int i = 0; i < NBins; i++){

		BsRAAY[i] = BsXsecPbPbY[i] / BsXsecPPY[i];
	
		cout << "i = " << i <<  "    BsRAAY[i] = " << BsRAAY[i] << endl;

		BsRAAYErrUp[i] = BsRAAY[i] * TMath::Sqrt(BsXSecPbPbYErrUpPercent[i] * BsXSecPbPbYErrUpPercent[i] + BsXSecPPYErrUpPercent[i] * BsXSecPPYErrUpPercent[i]);
		BsRAAYErrDown[i] = BsRAAY[i] * TMath::Sqrt(BsXSecPbPbYErrDownPercent[i] * BsXSecPbPbYErrDownPercent[i] + BsXSecPPYErrDownPercent[i] * BsXSecPPYErrDownPercent[i]);
		

		BsRAAYSystDown[i] = BsRAAY[i] * TMath::Sqrt(BsXSecPbPbYSystDownPercent[i] * BsXSecPbPbYSystDownPercent[i] + BsTotalSystDown[i] * BsTotalSystDown[i]);
		BsRAAYSystUp[i] = BsRAAY[i] * TMath::Sqrt(BsXSecPbPbYSystUpPercent[i] * BsXSecPbPbYSystUpPercent[i] + BsTotalSystUp[i] * BsTotalSystUp[i]);


		cout << "i = " << i <<  "    BPRAAY[i] = " << BsRAAY[i] << "   BPRAAYErrUp[i]  = "  << BsRAAYErrUp[i]/BsRAAY[i]   << "  BPRAAYErrDown[i] =   " << BsRAAYErrDown[i]/BsRAAY[i] <<   "   BPRAAYSystUp[i]  = "  << BsRAAYSystUp[i]/BsRAAY[i]   << "  BPRAAYSystDown[i] =   " << BsRAAYSystDown[i]/BsRAAY[i]  << endl;


			
	}



	





	TGraphAsymmErrors *BsRAAGraph = new TGraphAsymmErrors(NBins, BsRAAX, BsRAAY,BsRAAXErrDown, BsRAAXErrUp,BsRAAYErrDown,BsRAAYErrUp);

	BsRAAGraph->SetName("BsRAAGraph");
	TGraphAsymmErrors *BsRAAGraphSyst = new TGraphAsymmErrors(NBins, BsRAAX, BsRAAY,BsRAAXErrDown, BsRAAXErrUp,BsRAAYSystDown,BsRAAYSystUp);



	TGraphAsymmErrors *BsRAAGraph2015 = new TGraphAsymmErrors(NBins2015, BsRAAX2015, BsRAAY2015,BsRAAXErrDown2015, BsRAAXErrUp2015,BsRAAYErrDown2015,BsRAAYErrUp2015);
	TGraphAsymmErrors *BsRAAGraphSyst2015 = new TGraphAsymmErrors(NBins2015, BsRAAX2015, BsRAAY2015,BsRAAXErrDown2015, BsRAAXErrUp2015,BsRAAYSystDown2015,BsRAAYSystUp2015);





	BsRAAGraph->SetLineColor(kRed+2);
//	BsRAAGraph->SetFillColorAlpha(kRed+2,0.5);
	BsRAAGraph->SetMarkerStyle(20);
	BsRAAGraph->SetMarkerSize(1);

	BsRAAGraph->SetMarkerColor(kRed+2);

	BsRAAGraph2015->SetLineColor(kBlue+2);
//	BsRAAGraph->SetFillColorAlpha(kRed+2,0.5);
	BsRAAGraph2015->SetMarkerStyle(21);
	BsRAAGraph2015->SetMarkerSize(1);
	BsRAAGraph2015->SetMarkerColor(kBlue+2);




	BsRAAGraphSyst->SetFillColorAlpha(kRed-9,0.5);
	BsRAAGraphSyst->SetLineColor(kRed-9);
	BsRAAGraphSyst2015->SetFillColorAlpha(kBlue-9,0.5);
	BsRAAGraphSyst2015->SetLineColor(kBlue-9);


	HisEmptyRAA->Draw();
	BsRAAGraph->Draw("ep");
	BsRAAGraphSyst->Draw("5same");
	Unity->Draw("SAME");
	c2->SaveAs("RAAPlots/Bs/BsRAA.png");
	c2->SaveAs("RAAPlots/Bs/BsRAA.pdf");


	HisEmptyRAA->Draw();

	BsRAAGraph->Draw("ep");
	BsRAAGraph2015->Draw("epSAME");
	BsRAAGraph2015->Draw("epSAME");
	BsRAAGraphSyst->Draw("5same");
	BsRAAGraphSyst2015->Draw("5same");
	Unity->Draw("SAME");



	TLegend* leg2 = new TLegend(0.30,0.70,0.60,0.90,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.040);
	leg2->SetTextFont(42);
	leg2->SetFillStyle(0);
	leg2->SetLineWidth(3);
	leg2->AddEntry(BsRAAGraph,"2018 PbPb + 2017 pp","PL");
	leg2->AddEntry(BsRAAGraph2015,"2015 PbPb + 2015 pp","PL");
	leg2->Draw("same");






	c2->SaveAs("RAAPlots/Bs/BsRAAComparison.png");
	c2->SaveAs("RAAPlots/Bs/BsRAAComparison.pdf");


	TFile * fout = new TFile("OutFile/BsRAA.root","RECREATE");

	fout->cd();

	BsRAAGraph->Write();

	fout->Close();


}
