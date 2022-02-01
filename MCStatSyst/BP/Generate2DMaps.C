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
//#include "Rebin.h"



//#include "his.h"
using namespace std;

using std::cout;
using std::endl;



void Generate2DMaps(){
	
	float Factor = 4.0;

	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);



	const int yBinN = 5;
	double yBinning[yBinN+1] = {0.0,0.5, 1.0, 1.5,2.0, 2.4};




	double LowBinWidth = 0.5/Factor;
	int NLowBin = 10/LowBinWidth;
	//	int NLowBin = 5;

	double MidBinWidth = 1/Factor;
	int NMidBin = 10/MidBinWidth;
	double HighBinWidth = 1/Factor;
	int NHighBin = 40/HighBinWidth;
	const int BptBin = NHighBin + NMidBin + NLowBin;
	double BptBinning[BptBin + 1];


	//cout << "NLowBin = " << NLowBin << " NMidBin = " << NMidBin << " NHighBin = " << NHighBin << endl;

	for(int i = 0; i < NLowBin; i++){
		BptBinning[i] = 0 + i * LowBinWidth;
	}
	for(int i = 0; i < NMidBin; i++){
		BptBinning[i+NLowBin] = 10 + i * MidBinWidth;
	}
	for(int i = 0; i <  NHighBin+1; i++){
		BptBinning[i+NLowBin+NMidBin] = 20 + i * HighBinWidth;
	}



	TString FileName = "2DMaps/EffFineBDT.root";

	TFile * fin = new TFile(FileName.Data());
	fin->cd();

	TH2D * EffBptByInv = (TH2D * ) fin->Get("invEff2D");
	TH2D * hEffSelectionFine = (TH2D * ) fin->Get("invEffonly2D");
	TH2D * hEffAcc2DFine = (TH2D * ) fin->Get("invAcc2D");


	int NTrials = 10000;

	double EffInv2DMapCenter[BptBin][yBinN];
	double EffInv2DMapError[BptBin][yBinN];


	double EffSel2DMapCenter[BptBin][yBinN];
	double EffSel2DMapError[BptBin][yBinN];


	double Acc2DMapCenter[BptBin][yBinN];
	double Acc2DMapError[BptBin][yBinN];


	//loop through the 2D map histogram//



		for(int i = 0; i < BptBin; i++){

			for(int j = 0; j < yBinN; j++){

				EffInv2DMapCenter[i][j] = EffBptByInv->GetBinContent(i+1,j+1);
				EffInv2DMapError[i][j] = EffBptByInv->GetBinError(i+1,j+1);

				Acc2DMapCenter[i][j] = hEffAcc2DFine->GetBinContent(i+1,j+1);
				Acc2DMapError[i][j] = hEffAcc2DFine->GetBinError(i+1,j+1);

				EffSel2DMapCenter[i][j] = hEffSelectionFine->GetBinContent(i+1,j+1);
				EffSel2DMapError[i][j] = hEffSelectionFine->GetBinError(i+1,j+1);

	//			cout << "i = " << i << "   j =  "<< j << "  Acc2DMapCenter = " << Acc2DMapCenter[i][j] << "  Acc2DMapError = " << Acc2DMapError[i][j] << endl;
	//			cout << "i = " << i << "   j =  "<< j << "  EffSel2DMapCenter = " << EffSel2DMapCenter[i][j] << "  EffSel2DMapError = " << EffSel2DMapError[i][j] << endl;

			}

		}

	






	TFile * fout = new TFile("OutFiles/GenStatSyst.root","RECREATE");
	fout->cd();


	TH2D * EffBptByInvTrial;
	TH2D * SelBptByInvTrial;
	TH2D * AccBptByInvTrial;


	double GeneratedEffInv;
	double GeneratedSelInv;
	double GeneratedAccInv;


	for(int iTrial=0;iTrial < NTrials; iTrial++){

		cout << "Now working on Trial  " << iTrial <<  "  out of  " << NTrials  << endl;

		EffBptByInvTrial = new TH2D(Form("EffBptByInvTrial%d",iTrial),"",BptBin,BptBinning,yBinN,yBinning);
		SelBptByInvTrial = new TH2D(Form("SelBptByInvTrial%d",iTrial),"",BptBin,BptBinning,yBinN,yBinning);
		AccBptByInvTrial = new TH2D(Form("AccBptByInvTrial%d",iTrial),"",BptBin,BptBinning,yBinN,yBinning);

		for(int i = 0; i < BptBin; i++){

			for(int j = 0; j < yBinN; j++){
				GeneratedEffInv = gRandom->Gaus(EffInv2DMapCenter[i][j],EffInv2DMapError[i][j]);
				GeneratedSelInv = gRandom->Gaus(EffSel2DMapCenter[i][j],EffSel2DMapCenter[i][j]);
				GeneratedAccInv = gRandom->Gaus(Acc2DMapCenter[i][j],Acc2DMapError[i][j]);

				//cout << "GeneratedEffInv = " << GeneratedEffInv << endl;

				EffBptByInvTrial->SetBinContent(i+1,j+1,GeneratedEffInv);
				SelBptByInvTrial->SetBinContent(i+1,j+1,GeneratedSelInv);
				AccBptByInvTrial->SetBinContent(i+1,j+1,GeneratedAccInv);
			}

		}

		EffBptByInvTrial->Write();
		SelBptByInvTrial->Write();
		AccBptByInvTrial->Write();

	}

	fout->Close();

}
