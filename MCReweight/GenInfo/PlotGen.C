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




//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void PlotGen(int Opt){


	TCanvas *c = new TCanvas();
	c->cd();

	gStyle->SetOptStat(0);



	TString infile;


	if(Opt == 0) infile = "../../UnskimmedSamples/OfficialMC/BPMC.root";
	if(Opt == 1)  infile = "/data/szhaozho/ppNewTMVA/CMSSW_10_3_2/src/Bs/ComputBDTOfficial/BsMC.root";

	TString Part;

	if(Opt == 0) Part = "BP";
	if(Opt == 1)  Part = "Bs";



	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TTree * ntGen = (TTree *) fin->Get("Bfinder/ntGen");
	
	ntGen->AddFriend("hiEvtAnalyzer/HiTree");

	TH1D * Gpt = new TH1D("Gpt","",100,0,100);
	Gpt->GetXaxis()->SetTitle("Generated B p_{T}");
	Gpt->GetYaxis()->SetTitle("Normalized Counts");
	Gpt->GetYaxis()->SetTitleOffset(1.4);
	Gpt->GetXaxis()->CenterTitle();
	Gpt->GetYaxis()->CenterTitle();

	Gpt->SetMarkerStyle(20);
	Gpt->SetMarkerSize(1);
	Gpt->SetMarkerColor(1);
	Gpt->SetLineColor(1);


	TH1D * JPsiPt = new TH1D("JPsiPt","",100,0,25);
	JPsiPt->GetXaxis()->SetTitle("Generated J/#psi p_{T}");
	JPsiPt->GetYaxis()->SetTitle("Normalized Counts");
	JPsiPt->GetYaxis()->SetTitleOffset(1.4);
	JPsiPt->GetXaxis()->CenterTitle();
	JPsiPt->GetYaxis()->CenterTitle();


	JPsiPt->SetMarkerStyle(20);
	JPsiPt->SetMarkerSize(1);
	JPsiPt->SetMarkerColor(1);
	JPsiPt->SetLineColor(1);


	TH1D * pthat = new TH1D("pthat","",100,0,100);
	pthat->GetXaxis()->SetTitle("pthat");
	pthat->GetYaxis()->SetTitle("Normalized Counts");
	pthat->GetYaxis()->SetTitleOffset(1.4);
	pthat->GetXaxis()->CenterTitle();
	pthat->GetYaxis()->CenterTitle();


	pthat->SetMarkerStyle(20);
	pthat->SetMarkerSize(1);
	pthat->SetMarkerColor(1);
	pthat->SetLineColor(1);


	TCut Weight = "weight"; 
	TCut GenCut;
	TCut JPsiGenCut;


	if(Opt == 0) GenCut = "(TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0)";
	if(Opt == 1) GenCut = "(TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==531 && GisSignal>0)";



	if(Opt == 0) JPsiGenCut = "(TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==443)";
	if(Opt == 1) JPsiGenCut = "(TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==443)";


	ntGen->Project("Gpt","Gpt",TCut(Weight) * TCut(GenCut));
	ntGen->Project("JPsiPt","Gpt",TCut(Weight) * TCut(JPsiGenCut));
	ntGen->Project("pthat","pthat",TCut(Weight));

	c->SetLogy();


	pthat->Scale(1.0/pthat->Integral());
	pthat->Draw("ep");
	c->SaveAs(Form("MCPlots/%spthat.png",Part.Data()));	


	Gpt->Scale(1.0/Gpt->Integral());
	Gpt->Draw("ep");
	c->SaveAs(Form("MCPlots/%sGpt.png",Part.Data()));	


	JPsiPt->Scale(1.0/JPsiPt->Integral());
	JPsiPt->Draw("ep");
	c->SaveAs(Form("MCPlots/%sJPsiPt.png",Part.Data()));	
	


}
