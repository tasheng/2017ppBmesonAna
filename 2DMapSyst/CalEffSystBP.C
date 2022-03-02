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

void CalEffSystBP(){

	const int NBins = 7;
	//const int NBins = 6;

	int TnP = 1;


	double BRchain = 6.02061e-5;

//	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);


	TString FileName;

	FileName = "../SkimmedSamples/BPData.root";
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
	

	Float_t BEffInvTnPUp[NCand];
	Float_t BEffInvTnPDown[NCand];
	Float_t BEffInvBDT[NCand];
	Float_t BEffInvBpt[NCand];


	/*
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
	*/





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


	//Syst Collection
	double SumCountsTnPUpSyst[NBins];
	double SumCountsTnPDownSyst[NBins];
	double SumCountsBDTSyst[NBins];
	double SumCountsBptSyst[NBins];



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


//		ptbinsvec.push_back(3);
		ptbinsvec.push_back(5);
		ptbinsvec.push_back(7);		
		ptbinsvec.push_back(10);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(30);
		
		ptbinsvec.push_back(50);
		ptbinsvec.push_back(60);



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

		
		SumCountsTnPUpSyst[i] = 0;
		SumCountsTnPDownSyst[i] = 0;
		SumCountsBDTSyst[i] = 0;
		SumCountsBptSyst[i] = 0;



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


	TFile * finSyst2D = new TFile("../BP/EffAna/NewEff2DMaps/BPSyst2D.root");
	


	TH2D * invEff2D = (TH2D *) finSyst2D->Get("invEff2D");
	TH2D * invEff2DTnPSystUp = (TH2D *) finSyst2D->Get("invEff2DTnPSystUp");
	TH2D * invEff2DTnPSystDown = (TH2D *) finSyst2D->Get("invEff2DTnPSystDown");
	TH2D * invEff2DBDTSyst = (TH2D *) finSyst2D->Get("invEff2DBDTSyst");
	TH2D * invEff2DBptSyst = (TH2D *) finSyst2D->Get("invEff2DBptSyst");


	int XBin;
	int YBin;


	for( int i = 0; i < NEvents; i++){

		EffInfoTree->GetEntry(i);
		//MuonInfoTree->GetEntry(i);


		for(int j = 0; j < BsizeNew; j++){


			for(int k = 0; k < NBins; k++){

				//	if((BptNew[j] > ptBins[k] && BptNew[j] < ptBins[k+1] && TMath::Abs(BmassNew[j] - 5.27932) < 0.08  && ((BptNew[j] > 7 && BptNew[j] < 10 && ByNew[j] > 1.5 )||(BptNew[j] > 10)) && (Bmu1Type > -0.1 && Bmu2Type > -0.1)))
				if(BptNew[j] > ptBins[k] && BptNew[j] < ptBins[k+1] && TMath::Abs(BmassNew[j] - 5.27932) < 0.08 &&  TMath::Abs(ByNew[j]) < 2.4  && ((BptNew[j] > 5 && BptNew[j] < 10 && abs(ByNew[j]) > 1.5 )||(BptNew[j] > 10)))
				{


					XBin = invEff2D->GetXaxis()->FindBin( BptNew[j]);
					YBin = invEff2D->GetYaxis()->FindBin( TMath::Abs(ByNew[j]));
	
					BEffInv[j] = invEff2D->GetBinContent(XBin,YBin);


					BEffInvErr[j] = invEff2D->GetBinError(XBin,YBin);
					BEff[j] = 1.0/invEff2D->GetBinContent(XBin,YBin);

					BEffErr[j] = BEffInvErr[j]/(BEffInv[j] * BEffInv[j]);
	


					BEffInvTnPUp[j] = invEff2DTnPSystUp->GetBinContent(XBin,YBin);
					BEffInvTnPDown[j] = invEff2DTnPSystDown->GetBinContent(XBin,YBin);
					BEffInvBDT[j] = invEff2DBDTSyst->GetBinContent(XBin,YBin);
					BEffInvBpt[j] = invEff2DBptSyst->GetBinContent(XBin,YBin);


					if(BEffInv[j] > 0){
						SumCounts[k] = SumCounts[k] + BEffInv[j];
						SumCountsErr[k] = SumCountsErr[k] + BEffInvErr[j] * BEffInvErr[j];
						SumCountsEff[k] = SumCountsEff[k] + BEff[j];
						SumCountsEffErr[k] = SumCountsEffErr[k] + BEffErr[j] * BEffErr[j];
						SumCountsSyst[k] = 	SumCountsSyst[k]  + BEffInvBDTWeighted[j];
						SumCountsSystErr[k] = 	SumCountsSystErr[k]  + BEffInvErrBDTWeighted[j] * BEffInvErrBDTWeighted[j];

						SumCountsUp[k] = SumCountsUp[k] + BEffInvUp[j];
						SumCountsErrUp[k] = SumCountsErrUp[k] + BEffInvErrUp[j] * BEffInvErrUp[j];

						SumCountsDown[k] = SumCountsDown[k] + BEffInvDown[j];
						SumCountsErrDown[k] = SumCountsErrUp[k] + BEffInvErrDown[j] * BEffInvErrDown[j];


						SumCountsTnPUpSyst[k] = BEffInvTnPUp[j] + SumCountsTnPUpSyst[k];
						SumCountsTnPDownSyst[k] = BEffInvTnPDown[j] + SumCountsTnPDownSyst[k];

						SumCountsBDTSyst[k] = BEffInvBDT[j] + SumCountsBDTSyst[k];
						SumCountsBptSyst[k] = BEffInvBpt[j] + SumCountsBptSyst[k];

			
						Counts[k] = Counts[k] + 1;

						//cout << "SumCounts = " << SumCounts[k] << endl;

						//cout << "SumCountsUp = " << SumCountsUp[k] << endl;
						//cout << "Candidate: " << "   Bpt = " << BptNew[j] << "   Efficiency =  " << BEffInv[j] << "   Efficiency Error  = " << BEffInvErr[j] << endl;
					}
				}

			}


		}

	}



	double EffTnPUp[NBins];
	double EffTnPDown[NBins];
	double EffBDT[NBins];
	double EffBpt[NBins];


	for(int i = 0; i < NBins; i++){


		NewEff[i] = SumCounts[i]/Counts[i];
		NewEffErr[i] = TMath::Sqrt(SumCountsErr[i])/Counts[i];


		NewEffUp[i] = SumCountsUp[i]/Counts[i];
		NewEffErrUp[i] = TMath::Sqrt(SumCountsErrUp[i])/Counts[i];



		NewEffDown[i] = SumCountsDown[i]/Counts[i];
		NewEffErrDown[i] = TMath::Sqrt(SumCountsErrDown[i])/Counts[i];


		NewEffReal[i] = SumCountsEff[i]/Counts[i];
		NewEffRealErr[i] = TMath::Sqrt(SumCountsEffErr[i])/Counts[i];





		NewEffSyst[i] = SumCountsSyst[i]/Counts[i];
		NewEffSystErr[i] = TMath::Sqrt(SumCountsSystErr[i])/Counts[i];






		EffTnPUp[i] = SumCountsTnPUpSyst[i]/Counts[i];
		EffTnPDown[i] = SumCountsTnPDownSyst[i]/Counts[i];
		EffBDT[i] = SumCountsBDTSyst[i]/Counts[i];
		EffBpt[i] = SumCountsBptSyst[i]/Counts[i];








	}



	TH1D * Eff2DHis = new TH1D("Eff2DHis","",NBins,ptBins);

	Eff2DHis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff2DHis->GetYaxis()->SetTitle("<1/alpha #times #epsilon>");
	Eff2DHis->GetXaxis()->CenterTitle();	
	Eff2DHis->GetYaxis()->CenterTitle();
	Eff2DHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff2DHis->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Eff2DTnPUpSystHis = new TH1D("Eff2DTnPUpSystHis","",NBins,ptBins);

	Eff2DTnPUpSystHis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff2DTnPUpSystHis->GetYaxis()->SetTitle("<1/alpha #times #epsilon>");
	Eff2DTnPUpSystHis->GetXaxis()->CenterTitle();	
	Eff2DTnPUpSystHis->GetYaxis()->CenterTitle();
	Eff2DTnPUpSystHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff2DTnPUpSystHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff2DTnPDownSystHis = new TH1D("Eff2DTnPDownSystHis","",NBins,ptBins);

	Eff2DTnPDownSystHis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff2DTnPDownSystHis->GetYaxis()->SetTitle("<1/alpha #times #epsilon>");
	Eff2DTnPDownSystHis->GetXaxis()->CenterTitle();	
	Eff2DTnPDownSystHis->GetYaxis()->CenterTitle();
	Eff2DTnPDownSystHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff2DTnPDownSystHis->GetYaxis()->SetTitleOffset(1.5);



	TH1D * Eff2DBDTHis = new TH1D("Eff2DBDTHis","",NBins,ptBins);

	Eff2DBDTHis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff2DBDTHis->GetYaxis()->SetTitle("<1/alpha #times #epsilon>");
	Eff2DBDTHis->GetXaxis()->CenterTitle();	
	Eff2DBDTHis->GetYaxis()->CenterTitle();
	Eff2DBDTHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff2DBDTHis->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Eff2DBptHis = new TH1D("Eff2DBptHis","",NBins,ptBins);

	Eff2DBptHis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff2DBptHis->GetYaxis()->SetTitle("<1/alpha #times #epsilon>");
	Eff2DBptHis->GetXaxis()->CenterTitle();	
	Eff2DBptHis->GetYaxis()->CenterTitle();
	Eff2DBptHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff2DBptHis->GetYaxis()->SetTitleOffset(1.5);


	for(int i = 0; i < NBins; i++){

		Eff2DHis->SetBinContent(i+1, NewEff[i]);
		Eff2DHis->SetBinError(i+1, NewEffErr[i]);


		Eff2DTnPUpSystHis->SetBinContent(i+1, EffTnPUp[i]);		
		Eff2DTnPUpSystHis->SetBinError(i+1, NewEffErr[i]);

		Eff2DTnPDownSystHis->SetBinContent(i+1, EffTnPDown[i]);		
		Eff2DTnPDownSystHis->SetBinError(i+1, EffTnPUp[i]);
		

		Eff2DBDTHis->SetBinContent(i+1, EffBDT[i]);		
		Eff2DBDTHis->SetBinError(i+1, NewEffErr[i]);
	

		Eff2DBptHis->SetBinContent(i+1, EffBpt[i]);		
		Eff2DBptHis->SetBinError(i+1, NewEffErr[i]);


	}
	


	TFile * fout = new TFile("OutFiles/BPSyst2D.root","RECREATE");
	fout->cd();

	Eff2DHis->Write();
	Eff2DTnPUpSystHis->Write();
	Eff2DTnPDownSystHis->Write();
	Eff2DBDTHis->Write();
	Eff2DBptHis->Write();

	fout->Close();


}
