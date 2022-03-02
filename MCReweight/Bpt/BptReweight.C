#ifndef __CINT__
#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TCut.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>
using namespace std;

using std::cout;
using std::endl;
#endif


void BptReweight(int Opt){

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	gStyle->SetOptStat(0);

	TString infileMC;
	TString infileData;

	TString PartName;

	if(Opt == 0) infileMC = "../../SkimmedSamples/BPMC.root";
	if(Opt == 1)  infileMC = "../../SkimmedSamples/BsMC.root";

	if(Opt == 0) infileData = "../../BP/RawYieldFits/ROOTfiles/yields_Bp_binned_pt.root";
	if(Opt == 1) infileData = "../../Bs/RawYieldFits/ROOTfiles/yields_Bs_binned_pt.root";


	if(Opt == 0) PartName = "BP";
	if(Opt == 1) PartName = "Bs";

	
	TString TreeName;

	if(Opt == 0) TreeName = "ntKp";
	if(Opt == 1) TreeName = "ntphi";

	const int NPtBins = (Opt==1)? 4:7;

	float PtBins[NPtBins +1];

	std::vector<float> PtVec;
	

	if(NPtBins == 4){

		PtVec.push_back(7);
		PtVec.push_back(10);
		PtVec.push_back(15);
		PtVec.push_back(20);
		PtVec.push_back(50);

	}

	if(NPtBins == 7){

		PtVec.push_back(5);		
		PtVec.push_back(7);
		PtVec.push_back(10);
		PtVec.push_back(15);
		PtVec.push_back(20);
		PtVec.push_back(30);		
		PtVec.push_back(50);
		PtVec.push_back(60);

	}


	for(int i = 0; i< NPtBins+1; i++){

		PtBins[i] = PtVec[i];

	}

	TFile * finMC = new TFile(infileMC.Data());
	finMC->cd();

	TTree * ntKp = (TTree *) finMC->Get(TreeName.Data());
	TH1D * BptMC = new TH1D("BptMC","",NPtBins,PtBins);
	ntKp->Project("BptMC","Bpt","(Bgen == 23333) && ((Bpt > 5 && Bpt < 10 && abs(By) > 1.5  && abs(By) < 2.4)||(Bpt > 10 && abs(By) > 1.5  && abs(By) < 2.4))");


	BptMC->SetMarkerSize(1);
	BptMC->SetMarkerStyle(20);
	BptMC->SetMarkerColor(kRed);
	BptMC->SetLineColor(kRed);

	BptMC->Scale(1.0/BptMC->Integral());


	TFile * finData = new TFile(infileData); 
	finData->cd();

	TH1D * BptData = (TH1D *) finData->Get("hPt");
	
	BptData->SetMarkerSize(1);
	BptData->SetMarkerStyle(20);
	BptData->SetMarkerColor(kBlue);
	BptData->SetLineColor(kBlue);



	BptData->Scale(1.0/BptData->Integral());



	BptData->Draw("ep");
	BptMC->Draw("epSAME");

	c->SaveAs(Form("DataMCCompSide_%s.png",PartName.Data()));

	TH1D * DataMCRatio = (TH1D *) BptData->Clone("DataMCRatio");
	
	DataMCRatio->Sumw2();
	DataMCRatio->Divide(BptMC);

	DataMCRatio->SetMarkerSize(1);
	DataMCRatio->SetMarkerStyle(20);
	DataMCRatio->SetMarkerColor(kBlack);
	DataMCRatio->SetLineColor(kBlack);

	DataMCRatio->Draw("ep");


	c->SaveAs(Form("DataMCRatio_%s.png",PartName.Data()));



}
