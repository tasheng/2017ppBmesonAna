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
//#include "uti.h"
//#include "parameters.h"

//#include "his.h"
using namespace std;

using std::cout;
using std::endl;


void ReweightBpt(int Opt){


	TString inputMC;
	TString inputFONLL;
	TString MethodLabel = "FONLL";
	TCut GenCut;

	TString ParName;


	if(Opt == 0){

		ParName = "BP";
		inputMC =  "../Unskimmed/NewOfficialMC/BPMC.root";
		inputFONLL = "FONLLFine/BPFONLLFine.root";
		// inputFONLL = "FONLLFine/fonllOutput_pp_Bplus_5p03TeV_yFid.root";
		GenCut = "TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0";

	}

	if(Opt == 1){

		ParName = "Bs";
		inputMC =  "../Unskimmed/NewOfficialMC/BsMC.root";
		inputFONLL = "FONLLFine/BPFONLLFine.root";
		// inputFONLL = "FONLLFine/fonllOutput_pp_Bplus_5p03TeV_yFid.root";
		GenCut = "(TMath::Abs(Gy)<2.4&&TMath::Abs(GpdgId)==531&&GisSignal>0)";

	}



	gStyle->SetOptStat(0); 
	TCanvas * c = new TCanvas("c","c",1800,600);
	c->Divide(3,1);
	TCut weighpthat = "(weight)";
	

	TFile *finMC = new TFile(inputMC.Data());
	TFile *finFONLL = new TFile(inputFONLL.Data());

	TTree * tMC = (TTree*) finMC->Get("Bfinder/ntGen");
	tMC->AddFriend("hltanalysis/HltTree");
	tMC->AddFriend("hiEvtAnalyzer/HiTree");  
	tMC->AddFriend("skimanalysis/HltTree");

	// const int nBinsReweight = 180;
	const int nBinsReweight = 220;
	float ptBinsReweight[nBinsReweight + 1];
	float BptMin = 5.0;

	for(int i = 0; i < nBinsReweight + 1; i++){

		ptBinsReweight[i] = BptMin + 0.25 * i;

	}


	

	TH1D * GptMC = new TH1D("GptMC","GptMC",nBinsReweight,ptBinsReweight);
		


	GptMC->GetYaxis()->SetTitleOffset(1.1);

	tMC->Project("GptMC","Gpt",TCut(weighpthat)*TCut(GenCut));
	GptMC->Sumw2();

	GptMC->GetXaxis()->SetTitle("Gen p_{T} (GeV)");
	GptMC->GetYaxis()->SetTitle("PYTHIA, #entries ");
	GptMC->SetTitle("");
	GptMC->Sumw2();
	GptMC->Scale(1.0/GptMC->Integral());
	GptMC->SetMarkerStyle(20);
	GptMC->SetMarkerColor(kBlack);
	GptMC->SetMarkerSize(1);
	c->cd(1);
	gPad->SetLogy();
	GptMC->Draw("ep");


	TGraphAsymmErrors* gaeBplusReference = (TGraphAsymmErrors*)finFONLL->Get("gaeSigmaBplus");

	TH1D * GptFONLL = new TH1D("GptFONLL","",nBinsReweight,ptBinsReweight);

	double x;
	double y;
	double yErr;
	for(int i=0;i<nBinsReweight;i++){

		gaeBplusReference->GetPoint(i,x,y);

		cout << "i = " << i << "  x = " << x << "  y = "  << y << endl;

		
		yErr = gaeBplusReference->GetErrorY(i);
		GptFONLL->SetBinContent(i+1,y);
		//GptFONLL->SetBinError(i+1,yReweightBptErr);
	}




	GptFONLL->GetYaxis()->SetTitleOffset(1.4);
	GptFONLL->GetXaxis()->SetTitle(Form("%s p_{T} (GeV)",MethodLabel.Data()));
	GptFONLL->GetYaxis()->SetTitle(Form("%s_pp, #entries ",MethodLabel.Data()));
	GptFONLL->Sumw2();
	GptFONLL->Scale(1.0/GptFONLL->Integral());
	GptFONLL->SetMarkerStyle(20);
	GptFONLL->SetMarkerColor(kBlack);
	GptFONLL->SetMarkerSize(1);
	c->cd(2);
	gPad->SetLogy();
	GptFONLL->Draw("ep");


	TH1D * BptRatio = (TH1D *) GptFONLL->Clone("BptRatio");
	BptRatio->GetYaxis()->SetTitleOffset(1.3); 
	BptRatio->Divide(GptMC);
	c->cd(3);



	TF1 *f1;

	//if(MethodLabel == "FONLL") f1 = new TF1("f1","[0]/(x*x*x)-[1]/(x*x) + [2]",5,120);

	//if(MethodLabel == "NLO") f1 = new TF1("f1","[0]/(x*x*x)+ [1] * x + [2]",5,120);


	if(MethodLabel =="FONLL") f1 =  new TF1("f1"," [0]/x**[1] + [2] + [3] * x",5,100);
	
	if(MethodLabel =="NLO") f1 =  new TF1("f1"," ([0]-[1]*x)*TMath::Exp(-[2]*x) + [3] + [4] *x ",5,300);
	
	if(MethodLabel == "DATAEXP") f1 =  new TF1("f1"," ([0]-[1]*x)*TMath::Exp(-[2]*x) + [3] + [4] *x ",5,300);


	f1->SetParLimits(0,0,10000);
	f1->SetParLimits(1,0,4);
	f1->SetParLimits(2,0,2);
	f1->SetParLimits(3,-1,1.0);
	if(MethodLabel =="NLO") 	f1->SetParLimits(4,0,0.1);

	//TF1 * f1 = new TF1("f1","([0]-[1]*x)*TMath::Exp(-[2]*x)+[3]",5,105);

	f1->SetParLimits(2,0,10);
	f1->SetParLimits(1,0,100);


	BptRatio->Fit(f1,"R");
	//BptRatio->Fit(f1,"L q m","",5,105);

	BptRatio->GetXaxis()->SetTitle("B_{s}^{0} p_{T}");
	BptRatio->GetYaxis()->SetTitle(Form("%s_pp/PYTHIA",MethodLabel.Data()));
	BptRatio->SetMarkerStyle(20);
	BptRatio->SetMarkerColor(kBlack);
	BptRatio->SetMarkerSize(1);
	BptRatio->SetMinimum(0);
	BptRatio->SetMaximum(2.5);

	BptRatio->Draw("ep");

	float Chi2 =  f1->GetChisquare();
	int nDOF = f1->GetNDF();
    float nChi2 =  Chi2/nDOF;
   // float nChi2 = f->GetChisquare();
	int nDigit_chi2BakerCousins = 2;
    int nDigit_nChi2 = 2;
	//nChi2 = roundToNdigit(nChi2);
    nDigit_chi2BakerCousins = 2;

	//cout << "Chi2 = " <<  Chi2 << "  nChi2 = " << nChi2 << endl;

	// nDigit_nChi2 = sigDigitAfterDecimal(nChi2);
    /*
	TLatex* texChi = new TLatex(0.35,0.70, Form("#chi^{2}/nDOF: %.*f/%d = %.*f", nDigit_chi2BakerCousins, Chi2, nDOF, nDigit_nChi2, nChi2));
	texChi->SetNDC();
	texChi->SetTextAlign(12);
	texChi->SetTextSize(0.04);
	texChi->SetTextFont(42);
	texChi->Draw("SAME");
	*/

	//GptFONLL->Draw("ep");
	//
	//
//	TString BptReweightFunc =  Form("%f+%f/(Bpt*Bpt)",f1->GetParameter(0),f1->GetParameter(1));
	//TString BptReweightFunc = Form("%f*TMath::Exp(-%f*Bpt)+%f/(Bpt*Bpt+ %f * Bpt + %f*%f)",f1->GetParameter(3),f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(4),f1->GetParameter(2),f1->GetParameter(2));
	TString BptReweightFunc;

	//if(MethodLabel == "FONLL") BptReweightFunc = Form("%f/(x*x*x)-%f/(x*x)+%f",f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2));
	//if(MethodLabel == "NLO") BptReweightFunc = Form("%f/(x*x*x)-%f*x+%f",f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2));

	//if(MethodLabel == "FONLL") BptReweightFunc = Form("(%f - %f*x)*TMath::Exp(-%f * x) + %f",f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),f1->GetParameter(3));
	if(MethodLabel == "FONLL") BptReweightFunc = Form("%f/x**(%f) + %f + %f * x",f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),f1->GetParameter(3));
	if(MethodLabel == "NLO") BptReweightFunc = Form("(%f - %f*x)*TMath::Exp(-%f * x) + %f + %f * x",f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),f1->GetParameter(3),f1->GetParameter(4));
	if(MethodLabel == "DATAEXP") BptReweightFunc = Form("(%f - %f*x)*TMath::Exp(-%f * x) + %f + %f * x",f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),f1->GetParameter(3),f1->GetParameter(4));


//	TString BptReweightFunc = Form("(%f-%f*x)*TMath::Exp(-%f*x)+%f",f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),f1->GetParameter(3));

	cout << "Bpt Func = " << BptReweightFunc.Data() << endl;


	TLine * Unity = new TLine(5,1,50,1);
	Unity->SetLineWidth(2);
	Unity->SetLineStyle(2);
	Unity->SetLineColor(kBlue + 2);
	Unity->Draw();

	c->SaveAs(Form("plotReweight/%sBptReweigt%s.png",ParName.Data(),MethodLabel.Data()));
	c->SaveAs(Form("plotReweight/%sBptReweigt%s.pdf",ParName.Data(),MethodLabel.Data()));


	ofstream foutResults(Form("ResultFile/ReweightBpt_%s.txt", ParName.Data()));
	foutResults	<< "([0]-x)*TMath::Exp(-[1]*x)+[2]" << endl;
	foutResults	<< "Fitting Results: " << BptReweightFunc.Data() << endl;
	foutResults << "Paramater 0 = " << f1->GetParameter(0) << endl;
	foutResults << "Paramater 0 Error = " << f1->GetParError(0) << endl;
	foutResults << "Paramater 1 = " << f1->GetParameter(1) << endl;
	foutResults << "Paramater 1 Error = " << f1->GetParError(1) << endl;
	foutResults << "Paramater 2 = " << f1->GetParameter(2) << endl;
	foutResults << "Paramater 2 Error = " << f1->GetParError(2) << endl;

	foutResults << "Paramater 3 = " << f1->GetParameter(3) << endl;
	foutResults << "Paramater 3 Error = " << f1->GetParError(3) << endl;
/*
	foutResults << "Paramater 4 = " << f1->GetParameter(4) << endl;
	foutResults << "Paramater 4 Error = " << f1->GetParError(4) << endl;
*/

}

