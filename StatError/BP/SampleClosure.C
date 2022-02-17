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
#include <vector>
#include <random>

//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void SampleClosure(int PtOpt){


	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	double ptWidth;
	double RawYield;
	double CorrectedYield;

	double NCorrMin;
	double NCorrMax;


	double AverageRawYield = 0;
	double RawYieldPullValue;
	double RawYieldError;



//	double NEvents = 1000;

	int NEvents;

	int ptBins[2];
	double BsRaw = 0;

	TString outfilefolder;



	if(PtOpt == 0){

		NEvents = 31;
		ptBins[0] = 5;
		ptBins[1] = 7;
		NCorrMin = 10000;
		NCorrMax = 20000;
		BsRaw = 2000;
	}

	if(PtOpt == 1){

		NEvents = 31;
		ptBins[0] = 7;
		ptBins[1] = 10;
		NCorrMin = 13000;
		NCorrMax = 20000;
		BsRaw = 2000;
	}


	if(PtOpt == 2){

		NEvents = 31;
		ptBins[0] = 10;
		ptBins[1] = 15;
		NCorrMin = 30000;
		NCorrMax = 40000;
		BsRaw = 2000;
	}


	
	if(PtOpt == 3){

		NEvents = 31;
		ptBins[0] = 15;
		ptBins[1] = 20;
		NCorrMin = 9000;
		NCorrMax = 16000;
		BsRaw = 2000;
	}



	
	if(PtOpt == 4){

		NEvents = 31;
		ptBins[0] = 20;
		ptBins[1] = 30;
		NCorrMin = 5000;
		NCorrMax = 12000;
		BsRaw = 2000;
	}



	
	if(PtOpt == 5){

		NEvents = 31;
		ptBins[0] = 30;
		ptBins[1] = 50;
		NCorrMin = 0;
		NCorrMax = 4000;
		BsRaw = 2000;
	}


	
	if(PtOpt == 6){

		NEvents = 31;
		ptBins[0] = 50;
		ptBins[1] = 60;
		NCorrMin = 0;
		NCorrMax = 3000;
		BsRaw = 2000;
	}




	NEvents = BsRaw;
	const int NCand = 10;




	TH1D * SampleSize = new TH1D("SampleSize","",100, 0 ,100);
	if(PtOpt == 0 ) SampleSize = new TH1D("SampleSize","",100, 8000 ,35000);
	if(PtOpt == -1 ) SampleSize = new TH1D("SampleSize","",100, 100 ,200);



	SampleSize->GetXaxis()->SetTitle("Pseudo Data Resamples Size");
	SampleSize->GetYaxis()->SetTitle("Counts");
	SampleSize->SetTitle("");
	SampleSize->GetXaxis()->CenterTitle();
	SampleSize->GetYaxis()->CenterTitle();




	const int NBinsEffDis = 150;



	TH1D * EffDis = new TH1D("EffDis","",NBinsEffDis,0, 100);

	if(PtOpt == 0)  EffDis = new TH1D("EffDis","",NBinsEffDis,0, 20);
	EffDis->GetXaxis()->SetTitle("<1/acc x eff> Distribution");
	EffDis->GetYaxis()->SetTitle("Counts");
	EffDis->SetTitle("");
	EffDis->GetXaxis()->CenterTitle();
	EffDis->GetYaxis()->CenterTitle();



	TH1D * CorrYieldDis = new TH1D("CorrYieldDis","",50,NCorrMin, NCorrMax);
	CorrYieldDis->GetXaxis()->SetTitle("Corrected Yield Distribution");
	CorrYieldDis->GetYaxis()->SetTitle("Counts");
	CorrYieldDis->SetTitle("");
	CorrYieldDis->GetXaxis()->CenterTitle();
	CorrYieldDis->GetYaxis()->CenterTitle();



	TH1D * RawYieldPlot = new TH1D("RawYieldPlot","",50,NEvents * 0.7, NEvents * 1.3);
	if(PtOpt == 0) RawYieldPlot = new TH1D("RawYieldPlot","",50,500, 4000);
	if(PtOpt == 1) RawYieldPlot = new TH1D("RawYieldPlot","",50,280,450);
	if(PtOpt == 2) RawYieldPlot = new TH1D("RawYieldPlot","",50,200,400);
	if(PtOpt == 3) RawYieldPlot = new TH1D("RawYieldPlot","",50,260,420);


	RawYieldPlot->GetXaxis()->SetTitle("Raw Yield Distribution");
	RawYieldPlot->GetYaxis()->SetTitle("Counts");
	RawYieldPlot->SetTitle("");
	RawYieldPlot->GetXaxis()->CenterTitle();
	RawYieldPlot->GetYaxis()->CenterTitle();
	


	TH1D * RawYieldErrorPlot = new TH1D("RawYieldErrorPlot","",50,0,sqrt(NEvents) * 1.2);
	if(PtOpt == 1)  RawYieldErrorPlot = new TH1D("RawYieldErrorPlot","",50,17, 23);
	RawYieldErrorPlot->GetXaxis()->SetTitle("Raw Yield Error Distribution");
	RawYieldErrorPlot->GetYaxis()->SetTitle("Counts");
	RawYieldErrorPlot->SetTitle("");
	RawYieldErrorPlot->GetXaxis()->CenterTitle();
	RawYieldErrorPlot->GetYaxis()->CenterTitle();


	TH1D * RawYieldPull = new TH1D("RawYieldPull","",50,-5, 5);
	RawYieldPull->GetXaxis()->SetTitle("Raw Yield Pull Distribution");
	RawYieldPull->GetYaxis()->SetTitle("Counts");
	RawYieldPull->SetTitle("");
	RawYieldPull->GetXaxis()->CenterTitle();
	RawYieldPull->GetYaxis()->CenterTitle();



	int NFiles = 1000;

	TString FitFile; 
	TString EffFile; 

	TFile * finFit;


	int NCounts;
	double TotalEfficiency; 
	double Efficiency;

	double trgtnp1;
	double trktnp1;
	double muidtnp1;

	double trgtnp2;
	double trktnp2;
	double muidtnp2;

	double tnptotal1;
	double tnptotal2;

	int DoLater = 5;


	TFile * finEff;

	
    finEff = new TFile("2DEffMaps/EffFineBDT.root");
	finEff->cd();
	TH2D * invEff2D = (TH2D *) finEff->Get("invEff2D");


	int XBin;
	int YBin;


	for(int q = 1; q < NFiles; q++){

		cout << "q = " << q << endl;

		EffFile = Form("DataResample/%d-%d/Data_%d.root",ptBins[0],ptBins[1],q);
	
		TFile  * inf = new TFile(EffFile.Data());
		inf->cd();

		TTree * EffInfoTree = (TTree * ) inf->Get("ntKp");


		
		Efficiency = 0;
		TotalEfficiency = 0;
		NCounts = 0;

		Int_t BsizeNew;
		Float_t BmassNew[NCand];
		Float_t BptNew[NCand];
		Float_t ByNew[NCand];
		Float_t BEffInv[NCand];

		Int_t Bmu1Type[NCand];
		Int_t Bmu2Type[NCand];

		Float_t Bmu1ptNew[NCand];
		Float_t Bmu2ptNew[NCand];
		Float_t Bmu1etaNew[NCand];
		Float_t Bmu2etaNew[NCand];

		EffInfoTree->SetBranchAddress("Bsize",&BsizeNew);
//		EffInfoTree->SetBranchAddress("BEffInv",BEffInv);
		EffInfoTree->SetBranchAddress("Bmass",BmassNew);
		EffInfoTree->SetBranchAddress("By",ByNew);
		EffInfoTree->SetBranchAddress("Bpt",BptNew);
		EffInfoTree->SetBranchAddress("Bmu1pt",Bmu1ptNew);
		EffInfoTree->SetBranchAddress("Bmu2pt",Bmu2ptNew);
		EffInfoTree->SetBranchAddress("Bmu1eta",Bmu1etaNew);
		EffInfoTree->SetBranchAddress("Bmu2eta",Bmu2etaNew);


	//	cout << "Pass 1 " << endl;


		NEvents = EffInfoTree->GetEntries();
		SampleSize->Fill(NEvents);




		for( int i = 0; i < NEvents; i++){

			EffInfoTree->GetEntry(i);
			for(int j = 0; j < 1; j++){

				XBin = invEff2D->GetXaxis()->FindBin( BptNew[j]);
				YBin = invEff2D->GetYaxis()->FindBin( TMath::Abs(ByNew[j]));
				BEffInv[j] = invEff2D->GetBinContent(XBin,YBin);

				if((BptNew[j] > ptBins[0] && BptNew[j] < ptBins[1] && TMath::Abs(BmassNew[j] - 5.27932) < 0.08  && ((BptNew[j] > 5 && BptNew[j] < 10 && ByNew[j] > 1.5 )||(BptNew[j] > 10))))
				{


	
					TotalEfficiency = BEffInv[j] + TotalEfficiency;
					//cout << "TotalEfficiency = " << TotalEfficiency << "  BEffInv[j] = " << BEffInv[j] << endl;

					NCounts = NCounts + 1;
				}

			}

		}

		//		cout << "TotalEfficiency = " << TotalEfficiency << "  NCounts = " << NCounts << endl;

		if(NCounts > 0) Efficiency = TotalEfficiency/NCounts;

		FitFile = Form("RawYields/%d-%d/yields_Bp_full_%d.root",ptBins[0],ptBins[1],q);
		finFit	= new TFile(FitFile.Data());

		TH1D * hisFit = (TH1D *) finFit->Get("hPt");
		RawYield = hisFit->GetBinContent(1) * (ptBins[1] - ptBins[0]);
		RawYieldError =  hisFit->GetBinError(1) * (ptBins[1] - ptBins[0]);


		CorrectedYield = RawYield * Efficiency;

		cout << "RawYield = " << RawYield << "   Efficiency = " << Efficiency  <<   "   CorrectedYield = " << CorrectedYield << endl;

		EffDis->Fill(Efficiency);

		CorrYieldDis->Fill(CorrectedYield);


		AverageRawYield = AverageRawYield + RawYield;

		RawYieldPullValue = (RawYield - BsRaw)/RawYieldError;

//		finEff->Close();
		finFit->Close();
		inf->Close();

		RawYieldPlot->Fill(RawYield);
		RawYieldErrorPlot->Fill(RawYieldError);
		RawYieldPull->Fill(RawYieldPullValue);
	}

	AverageRawYield = AverageRawYield/NFiles;

	cout << "Average Raw Yield is " << AverageRawYield << endl;

	double StatError = CorrYieldDis->GetRMS()/CorrYieldDis->GetMean();

	cout << "Statistical Error = " << 	StatError << endl;

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	CorrYieldDis->Draw();

	TLatex * texMeanOnly2 = new TLatex(0.88,0.77,Form("RMS/Mean = %.3f",StatError));
	texMeanOnly2->SetNDC();
	texMeanOnly2->SetTextAlign(32);
	texMeanOnly2->SetTextFont(42);
	texMeanOnly2->SetTextSize(0.03);
	texMeanOnly2->SetLineWidth(2);
	texMeanOnly2->Draw("SAME");

	c->SaveAs(Form("Plots/CorrYield_%d_%d_.png",ptBins[0],ptBins[1]));

	double SampleSizeMean = SampleSize->GetMean();


	SampleSize->Draw();

	TLatex * SampleMeanTex = new TLatex(0.88,0.77,Form("NEvent Mean = %.0f",SampleSizeMean));
	SampleMeanTex->SetNDC();
	SampleMeanTex->SetTextAlign(32);
	SampleMeanTex->SetTextFont(42);
	SampleMeanTex->SetTextSize(0.03);
	SampleMeanTex->SetLineWidth(2);
	SampleMeanTex->Draw("SAME");

	c->SaveAs(Form("Plots/SampleSize_%d_%d_.png",ptBins[0],ptBins[1]));


	RawYieldPull->Draw();
	TF1 *f4 = new TF1("f4","gaus",-3,3);
	RawYieldPull->Fit(f4,"R");


	double	MeanYieldPullFinal = f4->GetParameter(1);
	double	WidthYieldPullFinal = f4->GetParameter(2);

	double	MeanYieldPullFinalErr = f4->GetParError(1);
	double	WidthYieldPullFinalErr = f4->GetParError(2);

	TLatex * texMeanYield = new TLatex(0.88,0.77,Form("Mean = %.3f #pm %.3f",MeanYieldPullFinal,MeanYieldPullFinalErr));
	texMeanYield->SetNDC();
	texMeanYield->SetTextAlign(32);
	texMeanYield->SetTextFont(42);
	texMeanYield->SetTextSize(0.03);
	texMeanYield->SetLineWidth(2);
	texMeanYield->Draw("SAME");


	TLatex * texWidthYield = new TLatex(0.88,0.67,Form("Width = %.3f #pm %.3f",WidthYieldPullFinal,WidthYieldPullFinalErr));
	texWidthYield->SetNDC();
	texWidthYield->SetTextAlign(32);
	texWidthYield->SetTextFont(42);
	texWidthYield->SetTextSize(0.03);
	texWidthYield->SetLineWidth(2);
	texWidthYield->Draw("SAME");

	c->SaveAs(Form("PlotsRaw/RawYieldPull_%d_%d_.png",ptBins[0],ptBins[1]));
	//c->SaveAs(Form("PlotsRaw/pdfForGJ/RawYieldPull_%d_%d_.pdf",ptBins[0],ptBins[1]));

	TF1 *f3 = new TF1("f3","gaus",5000,3500);


	RawYieldPlot->Draw();
	RawYieldPlot->Fit(f3,"R");


	double	MeanYieldFinal = f3->GetParameter(1);
	double	WidthYieldFinal = f3->GetParameter(2);

	double	MeanYieldFinalErr = f3->GetParError(1);
	double	WidthYieldFinalErr = f3->GetParError(2);

	TLatex * texMeanRawYield = new TLatex(0.88,0.77,Form("Mean = %.3f #pm %.3f",MeanYieldFinal,MeanYieldFinalErr));
	texMeanRawYield->SetNDC();
	texMeanRawYield->SetTextAlign(32);
	texMeanRawYield->SetTextFont(42);
	texMeanRawYield->SetTextSize(0.03);
	texMeanRawYield->SetLineWidth(2);
	texMeanRawYield->Draw("SAME");


	TLatex * texWidthRawYield = new TLatex(0.88,0.67,Form("Width = %.3f #pm %.3f",WidthYieldFinal,WidthYieldFinalErr));
	texWidthRawYield->SetNDC();
	texWidthRawYield->SetTextAlign(32);
	texWidthRawYield->SetTextFont(42);
	texWidthRawYield->SetTextSize(0.03);
	texWidthRawYield->SetLineWidth(2);
	texWidthRawYield->Draw("SAME");

	c->SaveAs(Form("PlotsRaw/RawYieldPlot_%d_%d_.png",ptBins[0],ptBins[1]));
	//c->SaveAs(Form("PlotsRaw/pdfForGJ/RawYieldPlot_%d_%d_.pdf",ptBins[0],ptBins[1]));

	TF1 *f2 = new TF1("f2","gaus",100,1000);

	RawYieldErrorPlot->Draw();
	RawYieldErrorPlot->Fit(f2,"R");

	double	MeanYieldErrorFinal = f2->GetParameter(1);
	double	WidthYieldErrorFinal = f2->GetParameter(2);

	double	MeanYieldErrorFinalErr = f2->GetParError(1);
	double	WidthYieldErrorFinalErr = f2->GetParError(2);



	TLatex * texMeanRawYieldError = new TLatex(0.88,0.77,Form("Mean = %.3f #pm %.3f",MeanYieldErrorFinal,MeanYieldErrorFinalErr));
	texMeanRawYieldError->SetNDC();
	texMeanRawYieldError->SetTextAlign(32);
	texMeanRawYieldError->SetTextFont(42);
	texMeanRawYieldError->SetTextSize(0.03);
	texMeanRawYieldError->SetLineWidth(2);
	texMeanRawYieldError->Draw("SAME");


	TLatex * texWidthRawYieldError = new TLatex(0.88,0.67,Form("Width = %.3f #pm %.3f",WidthYieldErrorFinal,WidthYieldErrorFinalErr));
	texWidthRawYieldError->SetNDC();
	texWidthRawYieldError->SetTextAlign(32);
	texWidthRawYieldError->SetTextFont(42);
	texWidthRawYieldError->SetTextSize(0.03);
	texWidthRawYieldError->SetLineWidth(2);
	texWidthRawYieldError->Draw("SAME");

	c->SaveAs(Form("PlotsRaw/RawYieldPlotError_%d_%d_.png",ptBins[0],ptBins[1]));

	//Eff Error//


	EffDis->Draw();


	double EffMeanCal = EffDis->GetMean();
	double EffWidthCal = EffDis->GetRMS();
	double EffError = EffWidthCal/EffMeanCal;



	TLatex * texEff = new TLatex(0.88,0.77,Form("RMS/Mean = %.3f",EffError));
	texEff->SetNDC();
	texEff->SetTextAlign(32);
	texEff->SetTextFont(42);
	texEff->SetTextSize(0.03);
	texEff->SetLineWidth(2);
	texEff->Draw("SAME");


	c->SaveAs(Form("PlotsEff/Eff_%d_%d_.png",ptBins[0],ptBins[1]));


	int MeanBinEff =  EffDis->GetXaxis()->FindBin(EffMeanCal);
	double PercentileEff = 0.3414 * 2;
	double UpPercentileEff = 0;
	double DownPercentileEff = 0;
	int UpBinEff = MeanBinEff;
	int DownBinEff = MeanBinEff;


	cout << "Mean Bin = " << MeanBinEff << endl;
	cout << "Eff :: Performing Up Down Integration" << endl;


	while(UpPercentileEff < PercentileEff){
		UpBinEff = UpBinEff + 1;
		UpPercentileEff = EffDis->Integral(MeanBinEff,UpBinEff)/EffDis->Integral(MeanBinEff,200);
		cout << "Current Up Percentile Eff = " << UpPercentileEff <<  endl;
	}

	while(DownPercentileEff < PercentileEff){
		DownBinEff = DownBinEff - 1;
		DownPercentileEff = EffDis->Integral(DownBinEff,MeanBinEff)/EffDis->Integral(0,MeanBinEff);
		cout << "Current Down Percentile Eff = " << DownPercentileEff <<  endl;

	}

	double SystErrorUpEff = EffDis->GetBinCenter(UpBinEff);
	double SystErrorDownEff = EffDis->GetBinCenter(DownBinEff);


	cout << "Method 2:  SystErrorUp = " << SystErrorUpEff << "  SystErrorDown = " << SystErrorDownEff << endl;


	cout << "INFO : MeanEffForCal = " << EffMeanCal << "  RMS = " << EffWidthCal << endl;

	double PercentageUpEff = (SystErrorUpEff - EffMeanCal)/EffMeanCal;

	double PercentageDownEff = (EffMeanCal - SystErrorDownEff)/EffMeanCal;

	cout << "Method 1 Eff:  Percentile = " << EffWidthCal/EffMeanCal << endl;
	cout << "Method 2 Eff:  PercentileUp = " << PercentageUpEff << "  PercentileDown = " << PercentageDownEff << endl;

	EffDis->Draw();
	texEff->Draw("SAME");
	TLatex * texEff2 = new TLatex(0.88,0.67,Form("RMS/Mean = +%.3f and -%.3f",PercentageUpEff, PercentageDownEff));
	texEff2->SetNDC();
	texEff2->SetTextAlign(32);
	texEff2->SetTextFont(42);
	texEff2->SetTextSize(0.03);
	texEff2->SetLineWidth(2);
	texEff2->Draw("SAME");

	c->SaveAs(Form("PlotsEff/EffDis_%d_%d_Assym.png",ptBins[0],ptBins[1]));

	//Asymmetric Error//


	double MeanEffForCal = CorrYieldDis->GetMean();
	double Percentile = 0.3414 * 2;
	int MeanBin =  CorrYieldDis->GetXaxis()->FindBin(MeanEffForCal);
	int ZeroBin =  CorrYieldDis->GetXaxis()->FindBin(0.0);

	double UpPercentile = 0;
	double DownPercentile = 0;
	int UpBin = MeanBin;
	int DownBin = MeanBin;


	cout << "Mean Bin = " << MeanBin << "  Zero Bin = " << ZeroBin << endl;
	cout << "Performing Up Down Integration" << endl;


	while(UpPercentile < Percentile){
		UpBin = UpBin + 1;
		UpPercentile = CorrYieldDis->Integral(MeanBin,UpBin)/CorrYieldDis->Integral(MeanBin,50);
		cout << "Current Up Percentile = " << UpPercentile <<  endl;
	}

	while(DownPercentile < Percentile){
		DownBin = DownBin - 1;
		DownPercentile = CorrYieldDis->Integral(DownBin,MeanBin)/CorrYieldDis->Integral(0,MeanBin);
		cout << "Current Down Percentile = " << DownPercentile <<  endl;

	}

	double SystErrorUp = CorrYieldDis->GetBinCenter(UpBin);
	double SystErrorDown = CorrYieldDis->GetBinCenter(DownBin);


	cout << "Method 2:  SystErrorUp = " << SystErrorUp << "  SystErrorDown = " << SystErrorDown << endl;


	cout << "INFO : MeanEffForCal = " << MeanEffForCal << "  RMS = " << CorrYieldDis->GetRMS() << endl;

	double PercentageUp = (SystErrorUp - MeanEffForCal)/MeanEffForCal;

	double PercentageDown = (MeanEffForCal - SystErrorDown)/MeanEffForCal;

	cout << "Method 1:  Percentile = " << CorrYieldDis->GetRMS()/CorrYieldDis->GetMean() << endl;
	cout << "Method 2:  PercentileUp = " << PercentageUp << "  PercentileDown = " << PercentageDown << endl;

	CorrYieldDis->Draw();
	texMeanOnly2->Draw("SAME");
	TLatex * texMeanOnly3 = new TLatex(0.88,0.67,Form("RMS/Mean = +%.3f and -%.3f",PercentageUp, PercentageDown));
	texMeanOnly3->SetNDC();
	texMeanOnly3->SetTextAlign(32);
	texMeanOnly3->SetTextFont(42);
	texMeanOnly3->SetTextSize(0.03);
	texMeanOnly3->SetLineWidth(2);
	texMeanOnly3->Draw("SAME");

	c->SaveAs(Form("Plots/CorrYield_%d_%d_Assym.png",ptBins[0],ptBins[1]));

//	c->SaveAs(Form("Plots/pdfForGJ/CorrYield_%d_%d_Assym.pdf",ptBins[0],ptBins[1]));





//	c->SaveAs(Form("PlotsEff/pdfForGJ/EffDis_%d_%d_Assym.pdf",ptBins[0],ptBins[1]));





}
