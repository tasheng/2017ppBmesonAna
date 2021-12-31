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


void FitBptShape(int Opt){

	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	gStyle->SetOptStat(0);

	TString BptWeightFile;
	TString OutName;

//	if(Opt == 0) BptWeightFile = "/Users/zhaozhongshi/Downloads/BPw.root"; 
//	if(Opt == 1)  BptWeightFile = "/Users/zhaozhongshi/Downloads/weights_Bs.root"; 


	if(Opt == 0) BptWeightFile = "../../BP/EffAna/BDTWeights/BPw.root"; 
	if(Opt == 1)  BptWeightFile = "../../Bs/EffAna/BDTWeights/Bsw.root"; 


	if(Opt == 0) OutName = "BPPtWeight.png"; 
	if(Opt == 1) OutName = "BsPtWeight.png"; 


	TFile * infile = new TFile(BptWeightFile.Data());
	infile->cd();


	TH1D * weights_Bpt = (TH1D *) infile->Get("weights_Bpt");

	weights_Bpt->SetMarkerSize(1);
	weights_Bpt->SetMarkerStyle(20);
	weights_Bpt->SetMarkerColor(kBlack);
	weights_Bpt->SetLineColor(kBlack);

	TF1 * f1 = new TF1("f1","[0]/(x*x) + [1] * TMath::Log(x) + [2]",0,50);

	f1->SetParLimits(0,1,100);
	f1->SetParLimits(1,0,10);
	f1->SetParLimits(2,-1.0,2);


	weights_Bpt->Fit(f1,"R");


	weights_Bpt->GetXaxis()->SetTitle("B p_{T} (GeV/c)");
	weights_Bpt->GetYaxis()->SetTitle("Data/MC");
	weights_Bpt->GetXaxis()->CenterTitle();	
	weights_Bpt->GetYaxis()->CenterTitle();


	weights_Bpt->Draw("ep");
	f1->Draw("SAME");
	

	c->SaveAs(OutName.Data());

	float p0 = f1->GetParameter(0);
	float p1 = f1->GetParameter(1);
	float p2 = f1->GetParameter(2);
	
	float Chi2ndf = f1->GetChisquare()/3;
	


	TString FuncForm = Form("%f/(x*x) +%f*TMath::Log(x) +%f",p0,p1,p2);


	cout <<  FuncForm.Data()  << "  Chi2ndf = " << Chi2ndf << endl;



}
