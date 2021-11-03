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

void ReAnaEffMulti(){

	const int NBins = 10;
	//const int NBins = 6;

	int TnP = 1;


	double BRchain = 6.02061e-5;

//	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	


	TString FileName;

	FileName = "/data/szhaozho/2017ppNewAna/BP/SkimmedOptCutBP.root";
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
	EffInfoTree->SetBranchAddress("BEffInv",BEffInv);
	EffInfoTree->SetBranchAddress("BEffInvErr",BEffInvErr);
	EffInfoTree->SetBranchAddress("BEff",BEff);
	EffInfoTree->SetBranchAddress("BEffInv1D",BEffInv1D);
	EffInfoTree->SetBranchAddress("BEffInvErr1D",BEffInvErr1D);

	EffInfoTree->SetBranchAddress("BEffInvFit",BEffInvFit);
	EffInfoTree->SetBranchAddress("BEffInvErrFit",BEffInvErrFit);

	EffInfoTree->SetBranchAddress("BEffInv",BEffInv);
	EffInfoTree->SetBranchAddress("BEffInvErr",BEffInvErr);

	EffInfoTree->SetBranchAddress("BEffInvBDTWeighted",BEffInvBDTWeighted);
	EffInfoTree->SetBranchAddress("BEffInvErrBDTWeighted",BEffInvErrBDTWeighted);



	EffInfoTree->SetBranchAddress("BEffInvUp",BEffInvUp);
	EffInfoTree->SetBranchAddress("BEffInvErrUp",BEffInvErrUp);

	EffInfoTree->SetBranchAddress("BEffInvDown",BEffInvDown);
	EffInfoTree->SetBranchAddress("BEffInvErrDown",BEffInvErrDown);


	/*
	   Int_t Bmu1Type[NCand];
	   Int_t Bmu2Type[NCand];





	   TTree * MuonInfoTree = (TTree * ) fin->Get("MuonInfoTree");

	   MuonInfoTree->SetBranchAddress("Bmu1Type",Bmu1Type);
	   MuonInfoTree->SetBranchAddress("Bmu2Type",Bmu2Type);


	   Float_t Bmu1etaNew[NCand];
	   Float_t Bmu2etaNew[NCand];

	   Float_t Bmu1ptNew[NCand];
	   Float_t Bmu2ptNew[NCand];


	   EffInfoTree->SetBranchAddress("Bmu1eta",Bmu1etaNew);
	   EffInfoTree->SetBranchAddress("Bmu2eta",Bmu2etaNew);

	   EffInfoTree->SetBranchAddress("Bmu1pt",Bmu1ptNew);
	   EffInfoTree->SetBranchAddress("Bmu2pt",Bmu2ptNew);

*/


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

		/*
		   corrfactvec.push_back(1.24759);
		   corrfactvec.push_back(1.05256);
		   corrfactvec.push_back(1.02614);
		   corrfactvec.push_back(1.01174);
		   */


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


		ptbinsvec.push_back(3);
		ptbinsvec.push_back(5);
		ptbinsvec.push_back(7);		
		ptbinsvec.push_back(10);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(50);
		ptbinsvec.push_back(100);



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
		ptbinsvec.push_back(0);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(25);
		ptbinsvec.push_back(30);
		ptbinsvec.push_back(35);
		ptbinsvec.push_back(40);		
		ptbinsvec.push_back(50);
		ptbinsvec.push_back(65);
		ptbinsvec.push_back(80);
		ptbinsvec.push_back(100);
		ptbinsvec.push_back(130);



	}


	for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ptbinsvec[i];
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




	}



	/*
	   const int NBins = 3;
	   double ptBins[NBins + 1] ={5,15,20,50};

	   int Counts[NBins]={0,0,0};
	   double SumCounts[NBins]={0,0,0};
	   double SumCountsErr[NBins]={0,0,0};
	   */


	/*
	   const int NBins = 4;
	   double ptBins[NBins + 1] ={5,10,15,20,50};

	   int Counts[NBins]={0,0,0,0};
	   double SumCounts[NBins]={0,0,0,0};
	   double SumCountsErr[NBins]={0,0,0,0};
	   */



	int EtaBin;
	int PtBin;




	double trgtnp1;
	double trktnp1;
	double muidtnp1;

	double trgtnp1systup;
	double trgtnp1systdown;
	double trgtnp1statup;
	double trgtnp1statdown;


	double trktnp1systup;
	double trktnp1systdown;
	double trktnp1statup;
	double trktnp1statdown;

	double muidtnp1systup;
	double muidtnp1systdown;
	double muidtnp1statup;
	double muidtnp1statdown;


	double tnptotal1;
	double tnptotal1up;
	double tnptotal1down;


	double tnptotal1systup;
	double tnptotal1systdown;
	double tnptotal1statup;
	double tnptotal1statdown;



	double trgtnp2;
	double trktnp2;
	double muidtnp2;

	double trgtnp2systup;
	double trgtnp2systdown;
	double trgtnp2statup;
	double trgtnp2statdown;


	double trktnp2systup;
	double trktnp2systdown;
	double trktnp2statup;
	double trktnp2statdown;

	double muidtnp2systup;
	double muidtnp2systdown;
	double muidtnp2statup;
	double muidtnp2statdown;


	double tnptotal2;
	double tnptotal2up;
	double tnptotal2down;

	double tnptotal2systup;
	double tnptotal2systdown;
	double tnptotal2statup;
	double tnptotal2statdown;



	double tnpabssystup;
	double tnpabssystdown;


	//TnP Syst DONE//

	TFile * RawYield = new TFile("RawYields/yields_Bp_binned_ptMulti.root");
	RawYield->cd();
	TH1D * hPt = (TH1D *) RawYield->Get("hPt");


	double RawCount;
	double RawCountErr;

	double CorrYield= 0;
	double CorrYieldErr= 0;

	double CorrYieldDiff[NBins];
	double CorrYieldDiffErr[NBins];




	//pt Binned Correction//
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	c->SetLogy();

	TH1D * CorrDiffHis = new TH1D("hPtSigma","",NBins,ptBins);




	CorrDiffHis->GetYaxis()->SetTitleOffset(1.3);
	CorrDiffHis->GetXaxis()->CenterTitle();
	CorrDiffHis->GetYaxis()->CenterTitle();

	TFile * fin1DEff = new TFile("NewEff2DMaps/EffFineBDT.root");
	fin1DEff->cd();
	TH1D * Eff1DHis = (TH1D * ) fin1DEff->Get("Eff1DHisMult");
	
	TH1D * CorrDiffHisBin = (TH1D * ) CorrDiffHis ->Clone("CorrDiffHisBin");

	CorrDiffHisBin->GetXaxis()->SetTitle("Event Multiplicity (N_{Trk})");
	CorrDiffHisBin->GetYaxis()->SetTitle("d #sigma/dN_{Trk} (pb)");

	CorrDiffHisBin->SetMarkerColor(kRed);
	CorrDiffHisBin->SetLineColor(kRed);	
	CorrDiffHisBin->SetMarkerSize(1);
	CorrDiffHisBin->SetMarkerStyle(20);

	float Eff1D[NBins];
	float Eff1DErr[NBins];

	for(int i = 0; i < NBins;i++){
		RawCount = hPt->GetBinContent(i+1);
		RawCountErr = hPt->GetBinError(i+1);
		Eff1D[i] = Eff1DHis->GetBinContent(i+1);
		Eff1DErr[i] = Eff1DHis->GetBinError(i+1);
	
//		CorrYieldDiff[i] = (RawCount *  Eff1D[i])/(BRchain*2* lumi);
//		CorrYieldDiffErr[i] = TMath::Sqrt((RawCountErr *  Eff1D[i]) *(RawCountErr  *  Eff1D[i]) + (RawCount *  Eff1DErr[i]) * (RawCount  *  Eff1DErr[i]))/(BRchain*2* lumi);
		CorrYieldDiff[i] = (RawCount /  Eff1D[i])/(BRchain*2* lumi);
		CorrYieldDiffErr[i] = TMath::Sqrt((RawCountErr /  Eff1D[i]) *(RawCountErr  /  Eff1D[i]) + (RawCount /Eff1D[i] *  Eff1DErr[i]) * (RawCount /Eff1D[i] *  Eff1DErr[i]))/(BRchain*2* lumi);

		CorrDiffHisBin->SetBinContent(i+1,CorrYieldDiff[i]);
		CorrDiffHisBin->SetBinError(i+1,CorrYieldDiffErr[i]);

	}
	
	CorrDiffHisBin->SetMaximum(5e6);
	
	CorrDiffHisBin->Draw("epSAME");


	//2015 B+ pp//


	const int NPt = 5;
	double Pt[NPt+1] = {7,10,15,20,30,50};
	double BPXSec[NPt] = {2610000,744000,197000,46500,5300};
	double StatErr[NPt] = {170000,29000,9000,2400,500};
	double SystErr[NPt] = {230000,59000,15000.3500,400};
	double TotalErr[NPt];

/*
	TH1D * BPXSecHis = new TH1D("BPXSecHis","",NPt,Pt);

	//BPXSecHis->SetTitle("2015 B^{+} pp Cross Section");
	BPXSecHis->SetMarkerColor(kGreen);
	BPXSecHis->SetLineColor(kGreen);
	
	BPXSecHis->SetMarkerSize(1);
	BPXSecHis->SetMarkerStyle(20);

	CorrDiffHisBin->GetXaxis()->SetTitle("Event Multiplicity (N_{Trk})");
	CorrDiffHisBin->GetYaxis()->SetTitle("d #sigma/dN_{Trk} (pb)");
	CorrDiffHisBin->GetYaxis()->SetTitleOffset(1.5);
	CorrDiffHisBin->GetXaxis()->CenterTitle();
	CorrDiffHisBin->GetYaxis()->CenterTitle();


	for(int i = 0; i < NPt; i++){

		TotalErr[i] = sqrt(StatErr[i] * StatErr[i] + SystErr[i] * SystErr[i]);

		BPXSecHis->SetBinContent(i+1,BPXSec[i]);
		BPXSecHis->SetBinError(i+1,TotalErr[i]);

	}

//	BPXSecHis->Draw("epSAME");
*/
	TLegend *leg = new TLegend(0.37,0.80,0.75,0.85,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);

//	leg->AddEntry(CorrDiffHis,"2017 pp","pl");
	leg->AddEntry(CorrDiffHisBin,"2017 pp Eff Binned Corr","pl");	
//	leg->AddEntry(BPXSecHis,"2015 pp","pl");
	
	//DONE 2015 B+ pp//

	leg->Draw("SAME");


	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.048);

	lat->DrawLatex(0.50,0.68,"B^{+} |y| < 2.4");

	c->SaveAs("FinalPlots/BPlusMulti.png");
	

	float InclBeautyXSec = 0;
	float InclBeautyXSecError = 0;
	
	float BXSecBin;
	float BXSecBinError;

	for(int i = 0; i < NBins;i++){
	
		BXSecBin = CorrDiffHisBin->GetBinContent(i+1) * (ptBins[i+1] - ptBins[i]);
		BXSecBinError = CorrDiffHisBin->GetBinError(i+1) * (ptBins[i+1] - ptBins[i]);

		InclBeautyXSec = InclBeautyXSec + BXSecBin;
		InclBeautyXSecError = InclBeautyXSecError + BXSecBinError * BXSecBinError;
	}

	InclBeautyXSec = InclBeautyXSec/100000;
	InclBeautyXSecError = sqrt(InclBeautyXSecError)/100000;

	cout << "----------------------------- BlockBuster !! --------------------------------------" << endl;

	cout << "INCLUSIVE pp to bb cross section at 5.02 TeV |y| < 2.4: " << InclBeautyXSec << " +- " << InclBeautyXSecError << "ub" << endl;
	cout << "----------------------------- BlockBuster !! --------------------------------------" << endl;




}
