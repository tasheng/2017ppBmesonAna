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
#include "tnp_weight_lowptPbPb.h"



//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void ReAnaEffSyst(int CentMin, int CentMax, int Weight ){

	int TnP = 1;
	
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);


	TString WeightName;
	/*

	   if(Weight == 16) WeightName = "NoWeight";
	   if(Weight == 1) WeightName = "FONLL";
	   if(Weight == 11) WeightName = "Linear";
	   if(Weight == 12) WeightName = "Quadratic";
	   if(Weight == 13) WeightName = "LInverse";
	   if(Weight == 14) WeightName = "LSqrt";
	   if(Weight == 15) WeightName = "LLog";
	   */

	if(Weight == 16) WeightName = "NoWeight";
	if(Weight == 1) WeightName = "FONLL";
	if(Weight == 11) WeightName = "Linear";
	if(Weight == 12) WeightName = "Quadratic";
	if(Weight == 13) WeightName = "LInverse";
	if(Weight == 14) WeightName = "LSqrt";
	if(Weight == 15) WeightName = "LLog";
	if(Weight == 0) WeightName = "NoTnP";


	TString FileName;
	/*
	   if(TnP == 0) FileName = Form("%dBptBins/EffInfo_%d_%d.root",NBptBins,CentMin,CentMax);
	   if(TnP == 1) FileName = Form("%dBptBinsTnP/EffInfo_%d_%d.root",NBptBins,CentMin,CentMax);
	   if(TnP == 2) FileName = Form("NonFiducial/EffInfo_%d_%d.root",CentMin,CentMax);
	   */

	FileName = Form("CheckSystNuno/%s/EffInfo_%d_%d.root",WeightName.Data(),CentMin,CentMax);

	TFile * fin = new TFile(FileName.Data());
	fin->cd();

	TTree * EffInfoTree = (TTree * ) fin->Get("EffInfoTree");
	//TTree * MuonInfoTree = (TTree * ) fin->Get("MuonInfoTree");


	int NEvents = EffInfoTree->GetEntries();

	const int NCand = 10;

	Int_t BsizeNew;
	Int_t runNew;
	Int_t lumiNew;
	Int_t evtNew;
	Float_t BmassNew[NCand];
	Float_t BptNew[NCand];
	Float_t ByNew[NCand];


	Int_t Bmu1Type[NCand]; 
	Int_t Bmu2Type[NCand]; 


	EffInfoTree->SetBranchAddress("BsizeNew",&BsizeNew);
	EffInfoTree->SetBranchAddress("BmassNew",BmassNew);
	EffInfoTree->SetBranchAddress("ByNew",ByNew);
	EffInfoTree->SetBranchAddress("BptNew",BptNew);


	TFile * finSystMap = new TFile(Form("10kTH2D/GenStatSyst_%d_%d.root",CentMin,CentMax));
	//TFile * finSystMap = new TFile(	Form("CheckSystNuno/%s/EffFine_%d_%d.root",WeightName.Data(),CentMin,CentMax));

	finSystMap->cd();

	TH2D * EffInv2DTrial;
	TH2D * SelInv2DTrial;
	TH2D * AccInv2DTrial;

	const int NTrials = 10000;

	//	const int NTrials = 1;

	TString outfileSyst;

	outfileSyst = Form("NunoSyst/%d-%d/AllTrials.root",CentMin,CentMax);



	TFile * finSyst = new TFile(outfileSyst.Data(),"RECREATE");
	finSyst->cd();
	TTree * NunoSyst = new TTree("NunoSyst","NunoSyst");


	Int_t BsizeNewSyst;
	Float_t BmassNewSyst[NCand];
	Float_t BptNewSyst[NCand];
	Float_t ByNewSyst[NCand];
	Float_t BEffInvSyst[NCand];
	Float_t BSelInvSyst[NCand];
	Float_t BAccInvSyst[NCand];

	Int_t Bmu1TypeSyst[NCand]; 
	Int_t Bmu2TypeSyst[NCand]; 


	NunoSyst->Branch("BsizeNewSyst",&BsizeNewSyst,"BsizeNewSyst/I");
	NunoSyst->Branch("BmassNewSyst",BmassNewSyst,"BmassNewSyst/F");
	NunoSyst->Branch("BptNewSyst",BptNewSyst,"BptNewSyst/F");
	NunoSyst->Branch("ByNewSyst",ByNewSyst,"ByNewSyst/F");
	NunoSyst->Branch("BEffInvSyst",BEffInvSyst,"BEffInvSyst/F");
	NunoSyst->Branch("BSelInvSyst",BSelInvSyst,"BSelInvSyst/F");
	NunoSyst->Branch("BAccInvSyst",BAccInvSyst,"BAccInvSyst/F");

	NunoSyst->Branch("Bmu1TypeSyst",Bmu1TypeSyst,"Bmu1TypeSyst/I");
	NunoSyst->Branch("Bmu2TypeSyst",Bmu2TypeSyst,"Bmu2TypeSyst/I");



	int XLOC;
	int YLOC;





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

	int EtaBin;
	int PtBin;

	double tnptotal1L2;
	double tnptotal1L3;
	double tnptotal2L2;
	double tnptotal2L3;





	Float_t Bmu1etaNew[NCand];
	Float_t Bmu2etaNew[NCand];

	Float_t Bmu1ptNew[NCand];
	Float_t Bmu2ptNew[NCand];


	EffInfoTree->SetBranchAddress("Bmu1etaNew",Bmu1etaNew);
	EffInfoTree->SetBranchAddress("Bmu2etaNew",Bmu2etaNew);

	EffInfoTree->SetBranchAddress("Bmu1ptNew",Bmu1ptNew);
	EffInfoTree->SetBranchAddress("Bmu2ptNew",Bmu2ptNew);


	for(int iTrial = 0; iTrial < NTrials;  iTrial++){

		cout << "Now Working on Trial = " << iTrial << " OUT OF " << NTrials << endl;

		EffInv2DTrial = (TH2D * ) finSystMap->Get(Form("EffBptByInvTrial%d",iTrial));
		//EffInv2DTrial = (TH2D * ) finSystMap->Get("EffBptByInv");
		SelInv2DTrial = (TH2D * ) finSystMap->Get(Form("SelBptByInvTrial%d",iTrial));
		AccInv2DTrial = (TH2D * ) finSystMap->Get(Form("AccBptByInvTrial%d",iTrial));



		for( int i = 0; i < NEvents; i++){

	//		MuonInfoTree->GetEntry(i);
			EffInfoTree->GetEntry(i);
			BsizeNewSyst = BsizeNew;

			for(int j = 0; j < BsizeNewSyst; j++){

				BmassNewSyst[j] = BmassNew[j];
				BptNewSyst[j] = BptNew[j];
				ByNewSyst[j] = ByNew[j];
				Bmu1TypeSyst[j] = Bmu1Type[j];
				Bmu2TypeSyst[j] = Bmu2Type[j];


				XLOC = EffInv2DTrial->GetXaxis()->FindBin(BptNewSyst[j]);
				YLOC = EffInv2DTrial->GetYaxis()->FindBin(ByNewSyst[j]);




				BEffInvSyst[j] = EffInv2DTrial->GetBinContent(XLOC,YLOC);
				BSelInvSyst[j] = SelInv2DTrial->GetBinContent(XLOC,YLOC);
				BAccInvSyst[j] = AccInv2DTrial->GetBinContent(XLOC,YLOC);
			



			}


			NunoSyst->Fill();
		}

	}

	finSyst->Write();
	finSyst->Close();

}
