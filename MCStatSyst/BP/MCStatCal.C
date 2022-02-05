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
//#include "tnp_weight_lowptPbPb.h"



//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void MCStatCal(){

	const int NFiles = 10000;
	const int NBins = 7;
	//const int NBins = 6;

	int TnP = 1;


	double BRchain = 6.02061e-5;

	//	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);


	TString FileName;

	FileName = "../../SkimmedSamples/BPData.root";
	TFile * fin = new TFile(FileName.Data());
	fin->cd();

	TTree * EffInfoTree = (TTree * ) fin->Get("ntKp");

	int NEvents = EffInfoTree->GetEntries();

	const int NCand = 10;

	Int_t BsizeNew;
	Int_t runNew;
	Int_t lumiNew;
	Int_t evtNew;
	Float_t BmassNew[NCand];
	Float_t BptNew[NCand];
	Float_t ByNew[NCand];
	Float_t BEff[NCand];
	Float_t BEffErr[NCand];


	Float_t BEffInv[NCand];
	Float_t BSelInv[NCand];
	Float_t BAccInv[NCand];

	Float_t BEffInvErr[NCand];
	Float_t BEffInv1D[NCand];
	Float_t BEffInvErr1D[NCand];

	Float_t BEffInvFit[NCand];
	Float_t BEffInvErrFit[NCand];

	Float_t BEffInvBDTWeighted[NCand];
	Float_t BEffInvErrBDTWeighted[NCand];



	Float_t BEffInvUp[NCand];
	Float_t BEffInvErrUp[NCand];
	Float_t BEffInvDown[NCand];
	Float_t BEffInvErrDown[NCand];

	EffInfoTree->SetBranchAddress("Bsize",&BsizeNew);
	EffInfoTree->SetBranchAddress("Bmass",BmassNew);
	EffInfoTree->SetBranchAddress("By",ByNew);
	EffInfoTree->SetBranchAddress("Bpt",BptNew);





	double ptBins[NBins + 1];

	int Counts[NBins];
	double SumCounts[NBins];
	double SumCountsErr[NBins];
	double NewEff[NBins];
	double NewEffErr[NBins];

	double NewEffReal[NBins];
	double NewEffRealErr[NBins];



	double SumCountsUp[NBins];
	double SumCountsErrUp[NBins];

	double SumCountsDown[NBins];
	double SumCountsErrDown[NBins];


	double SumCountsEff[NBins];
	double SumCountsEffErr[NBins];



	double SumCountsSyst[NBins];
	double SumCountsSystErr[NBins];
	double NewEffSyst[NBins];
	double NewEffSystErr[NBins];


	//	double CorrectionFactor[NBins];

	double NewEffUp[NBins];
	double NewEffErrUp[NBins];
	double NewEffDown[NBins];
	double NewEffErrDown[NBins];


	double lumi = 302.3;
	std::vector<double> ptbinsvec;
	std::vector<double> corrfactvec;

	std::vector<double> minvec;
	std::vector<double> maxvec;
	std::vector<double> meanvec;


	std::vector<double> minselvec;
	std::vector<double> maxselvec;



	std::vector<double> minaccvec;
	std::vector<double> maxaccvec;




	if(NBins == 1){


		ptbinsvec.push_back(10.0);
		ptbinsvec.push_back(50);
	}


	if(NBins == 3){

		ptbinsvec.push_back(5);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(50);

	}


	if(NBins == 4){

		ptbinsvec.push_back(5);
		ptbinsvec.push_back(10);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(50);



	}

	if(NBins == 6){


		ptbinsvec.push_back(5);
		ptbinsvec.push_back(7);		
		ptbinsvec.push_back(10);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(50);
		ptbinsvec.push_back(100);



	}



	if(NBins == 7){


		//		ptbinsvec.push_back(3);
		ptbinsvec.push_back(5);
		ptbinsvec.push_back(7);		
		ptbinsvec.push_back(10);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(30);

		ptbinsvec.push_back(50);
		ptbinsvec.push_back(60);




		minvec.push_back(45);
		minvec.push_back(22);
		minvec.push_back(11.3);
		minvec.push_back(5);
		minvec.push_back(3.5);
		minvec.push_back(2);
		minvec.push_back(1.5);


		minaccvec.push_back(4);
		minaccvec.push_back(3.5);
		minaccvec.push_back(2.1);
		minaccvec.push_back(1.57);
		minaccvec.push_back(1.57);
		minaccvec.push_back(1.57);
		minaccvec.push_back(1.57);


		minselvec.push_back(30);
		minselvec.push_back(10);
		minselvec.push_back(5);
		minselvec.push_back(0);
		minselvec.push_back(0);
		minselvec.push_back(0);
		minselvec.push_back(80);




		maxvec.push_back(60);
		maxvec.push_back(26);
		maxvec.push_back(13.5);
		maxvec.push_back(8);
		maxvec.push_back(5.5);
		maxvec.push_back(5);
		maxvec.push_back(8);


		maxaccvec.push_back(20);
		maxaccvec.push_back(20);
		maxaccvec.push_back(20);
		maxaccvec.push_back(20);
		maxaccvec.push_back(20);
		maxaccvec.push_back(20);
		maxaccvec.push_back(20);


		maxselvec.push_back(10);
		maxselvec.push_back(10);
		maxselvec.push_back(10);
		maxselvec.push_back(10);
		maxselvec.push_back(10);
		maxselvec.push_back(10);
		maxselvec.push_back(10);






		meanvec.push_back(53.538932);
		meanvec.push_back(23.943692);
		meanvec.push_back(12.038493);
		meanvec.push_back(5.7727937);
		meanvec.push_back(3.7848080);
		meanvec.push_back(3.0837052);
		meanvec.push_back(2.9782204);



	}



	if(NBins == 9){

		ptbinsvec.push_back(2);
		ptbinsvec.push_back(3);
		ptbinsvec.push_back(5);
		ptbinsvec.push_back(7);		
		ptbinsvec.push_back(10);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(30);
		ptbinsvec.push_back(50);
		ptbinsvec.push_back(100);



	}
	if(NBins == 10){

		ptbinsvec.push_back(1);
		ptbinsvec.push_back(2);
		ptbinsvec.push_back(3);
		ptbinsvec.push_back(5);
		ptbinsvec.push_back(7);		
		ptbinsvec.push_back(10);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(30);
		ptbinsvec.push_back(50);
		ptbinsvec.push_back(100);



	}



	double 	Mean[NBins];
	double 	Min[NBins];
	double 	Max[NBins];


	double 	MinAcc[NBins];
	double 	MaxAcc[NBins];

	double 	MinSel[NBins];
	double 	MaxSel[NBins];


	double SumCountsSel[NBins];
	double SumCountsAcc[NBins];

	double NewSel[NBins];
	double NewAcc[NBins];



	for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ptbinsvec[i];
	}


	for(int i = 0; i < NBins; i++){

		Mean[i] = meanvec[i];
		Min[i] = minvec[i];
		Max[i] = maxvec[i];

		MinAcc[i] = minaccvec[i];
		MaxAcc[i] = maxaccvec[i];

		MinSel[i] = minselvec[i];
		MaxSel[i] = maxselvec[i];


	}



	for(int i = 0; i < NBins; i++){
		Counts[i] = 0;
		SumCounts[i] = 0;
		SumCountsErr[i] = 0;
		SumCountsEff[i] = 0;
		SumCountsEffErr[i] = 0;
		SumCountsSyst[i] = 0;
		SumCountsSystErr[i] = 0;
		//	CorrectionFactor[i] = corrfactvec[i];
		SumCountsUp[i] = 0;
		SumCountsErrUp[i] = 0;
		SumCountsDown[i] = 0;
		SumCountsErrDown[i] = 0;

		SumCountsAcc[i] = 0;
		SumCountsSel[i] = 0;


	}


	TH1D * EffInvDistribution[NBins];
	TH1D * SelInvDistribution[NBins];
	TH1D * AccInvDistribution[NBins];

	for(int i = 0; i < NBins; i++){

		EffInvDistribution[i] = new TH1D("EffInvDistribution","",200,Min[i],Max[i]);
		EffInvDistribution[i]->GetXaxis()->SetTitle("<1/(acc x eff)>");
		EffInvDistribution[i]->GetYaxis()->SetTitle("Counts");
		EffInvDistribution[i]->SetTitle(Form("MC Smeared Distribution for %.0f < Bpt < %.0f",ptBins[i],ptBins[i+1]));

		EffInvDistribution[i]->GetXaxis()->CenterTitle();
		EffInvDistribution[i]->GetYaxis()->CenterTitle();
		EffInvDistribution[i]->GetYaxis()->SetTitleOffset(1.4);


		SelInvDistribution[i] = new TH1D("SelInvDistribution","",200,MinSel[i],MaxSel[i]);
		SelInvDistribution[i]->GetXaxis()->SetTitle("<1/(eff)>");
		SelInvDistribution[i]->GetYaxis()->SetTitle("Counts");
		SelInvDistribution[i]->SetTitle(Form("MC Smeared Distribution for %.0f < Bpt < %.0f",ptBins[i],ptBins[i+1]));

		SelInvDistribution[i]->GetXaxis()->CenterTitle();
		SelInvDistribution[i]->GetYaxis()->CenterTitle();
		SelInvDistribution[i]->GetYaxis()->SetTitleOffset(1.4);



		AccInvDistribution[i] = new TH1D("AccInvDistribution","",200,MinAcc[i],MaxAcc[i]);
		AccInvDistribution[i]->GetXaxis()->SetTitle("<1/(acc)>");
		AccInvDistribution[i]->GetYaxis()->SetTitle("Counts");
		AccInvDistribution[i]->SetTitle(Form("MC Smeared Distribution for %.0f < Bpt < %.0f",ptBins[i],ptBins[i+1]));

		AccInvDistribution[i]->GetXaxis()->CenterTitle();
		AccInvDistribution[i]->GetYaxis()->CenterTitle();
		AccInvDistribution[i]->GetYaxis()->SetTitleOffset(1.4);


	}



	int EtaBin;
	int PtBin;






	for(int q = 0; q < NFiles; q++){

		if(q%1 == 0)	cout << "Now Working On File " << q << endl;

		TFile * finEff;


		finEff = new TFile("OutFiles/GenStatSyst.root");

		finEff->cd();

		TH2D * invEff2D = (TH2D *) finEff->Get(Form("EffBptByInvTrial%d",q));
		TH2D * invSel2D = (TH2D *) finEff->Get(Form("SelBptByInvTrial%d",q));
		TH2D * invAcc2D = (TH2D *) finEff->Get(Form("AccBptByInvTrial%d",q));


		int XBin;
		int YBin;





		for( int i = 0; i < NEvents; i++){

			EffInfoTree->GetEntry(i);
			//MuonInfoTree->GetEntry(i);


			for(int j = 0; j < BsizeNew; j++){


				for(int k = 0; k < NBins; k++){

					if(BptNew[j] > ptBins[k] && BptNew[j] < ptBins[k+1] && TMath::Abs(BmassNew[j] - 5.27932) < 0.08 &&  TMath::Abs(ByNew[j]) < 2.4  && ((BptNew[j] > 5 && BptNew[j] < 10 && abs(ByNew[j]) > 1.5 )||(BptNew[j] > 10)))
					{


						XBin = invEff2D->GetXaxis()->FindBin( BptNew[j]);
						YBin = invEff2D->GetYaxis()->FindBin( TMath::Abs(ByNew[j]));
						BEffInv[j] = invEff2D->GetBinContent(XBin,YBin);
						BSelInv[j] = invSel2D->GetBinContent(XBin,YBin);						
						BAccInv[j] = invAcc2D->GetBinContent(XBin,YBin);





						if(BEffInv[j] > 0){
							SumCounts[k] = SumCounts[k] + BEffInv[j];
							SumCountsSel[k] = SumCountsSel[k] + BSelInv[j];
							SumCountsAcc[k] = SumCountsAcc[k] + BAccInv[j];

							Counts[k] = Counts[k] + 1;

						}
					}

				}


			}

		}


		for(int k = 0; k < NBins; k++){

			NewSel[k] = SumCountsSel[k]/Counts[k];
			NewAcc[k] = SumCountsAcc[k]/Counts[k];
			NewEff[k] = SumCounts[k]/Counts[k];
	
			//cout << "SumCounts[k] = " << SumCounts[k] << "    Counts[k]" << Counts[k]  << "   NewEff[k] = " << NewEff[k] << endl;

			EffInvDistribution[k]->Fill(NewEff[k]);
			SelInvDistribution[k]->Fill(NewSel[k]);			
			AccInvDistribution[k]->Fill(NewAcc[k]);



			Counts[k] = 0;
			SumCounts[k] = 0;
			SumCountsSel[k] = 0;
			SumCountsAcc[k] = 0;



		}

		finEff->Close();
	

	}


	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();




	TLine *l4[NBins];

	TLatex * texChi[NBins];


	//Efficiency//


	float RMS;
	float MeanValue;
	float Error;

	TH1D * MCSystStatHis = new TH1D("MCSystStatHis","",NBins,ptBins);
	MCSystStatHis->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	MCSystStatHis->GetYaxis()->SetTitle("MC Stat Syst (%)");
	MCSystStatHis->GetYaxis()->SetTitleOffset(1.4);
	MCSystStatHis->GetXaxis()->CenterTitle();
	MCSystStatHis->GetYaxis()->CenterTitle();

	MCSystStatHis->SetMarkerSize(1);
	MCSystStatHis->SetMarkerColor(1);
	MCSystStatHis->SetMarkerStyle(20);
	MCSystStatHis->SetLineColor(1);


	for(int i = 0; i < NBins; i++){

		EffInvDistribution[i]->SetMinimum(0);
		EffInvDistribution[i]->Draw();

		texChi[i] = new TLatex(0.10,0.85, Form("B^{+} GeV/c %.0f < p_{T} < %0.f GeV/c", ptBins[i],ptBins[i+1] ));
		texChi[i]->SetNDC();
		texChi[i]->SetTextAlign(12);
		texChi[i]->SetTextSize(0.06);
		texChi[i]->SetTextFont(42);
		texChi[i]->SetTextColor(1);

		texChi[i]->Draw("SAME");

		l4[i] = new TLine(Mean[i],0,Mean[i],EffInvDistribution[i]->GetMaximum());
		l4[i]->SetLineStyle(2);
		l4[i]->SetLineWidth(2);
		l4[i]->SetLineColor(2);
		l4[i]->Draw("SAME");

		c->SaveAs(Form("Plots/Eff/MCStatEff_%d.png",i));
		

		RMS = EffInvDistribution[i]->GetRMS();
		MeanValue = EffInvDistribution[i]->GetMean();
		Error = RMS/MeanValue;

		MCSystStatHis->SetBinContent(i+1,Error);
		MCSystStatHis->SetBinError(i+1,Error/100);
	
		cout << "<1/acc eff> Stat Syst: " <<  Error << endl;
	}


	MCSystStatHis->Draw("ep");
	c->SaveAs("Plots/MCStatSyst.png");



	//Selection//



	for(int i = 0; i < NBins; i++){



		SelInvDistribution[i]->SetMinimum(0);
		SelInvDistribution[i]->Draw();

		texChi[i] = new TLatex(0.10,0.85, Form("B^{+} GeV/c %.0f < p_{T} < %0.f GeV/c", ptBins[i],ptBins[i+1] ));
		texChi[i]->SetNDC();
		texChi[i]->SetTextAlign(12);
		texChi[i]->SetTextSize(0.06);
		texChi[i]->SetTextFont(42);
		texChi[i]->SetTextColor(1);

		texChi[i]->Draw("SAME");
		c->SaveAs(Form("Plots/Acc/MCStatAcc_%d.png",i));

		RMS = SelInvDistribution[i]->GetRMS();
		MeanValue = SelInvDistribution[i]->GetMean();
		Error = RMS/MeanValue;


		cout << "<1/eff> Stat Syst: " <<  Error << endl;
	}



	//Acceptance//

	for(int i = 0; i < NBins; i++){

		AccInvDistribution[i]->SetMinimum(0);
		AccInvDistribution[i]->Draw();

		texChi[i] = new TLatex(0.10,0.85, Form("B^{+} GeV/c %.0f < p_{T} < %0.f GeV/c", ptBins[i],ptBins[i+1] ));
		texChi[i]->SetNDC();
		texChi[i]->SetTextAlign(12);
		texChi[i]->SetTextSize(0.06);
		texChi[i]->SetTextFont(42);
		texChi[i]->SetTextColor(1);

		texChi[i]->Draw("SAME");
		c->SaveAs(Form("Plots/Sel/MCStatSel_%d.png",i));
	
	
		RMS = AccInvDistribution[i]->GetRMS();
		MeanValue = AccInvDistribution[i]->GetMean();
		Error = RMS/MeanValue;

		cout << "<1/acc> Stat Syst: " <<  Error << endl;


	}



}
