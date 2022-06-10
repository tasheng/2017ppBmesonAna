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


	if(Opt == 0) infile = "../../Unskimmed/NewOfficialMC/BPMC.root";
	if(Opt == 1) infile = "../../Unskimmed/NewOfficialMC/BsMC.root";

	TString Part;

	if(Opt == 0) Part = "BP";
	if(Opt == 1)  Part = "Bs";



	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TTree * ntGen = (TTree *) fin->Get("Bfinder/ntGen");
	
	ntGen->AddFriend("hiEvtAnalyzer/HiTree");

	TH1D * Gpt = new TH1D("Gpt","After reweighting",100,0,100);
	Gpt->GetXaxis()->SetTitle("Generated B p_{T}");
	Gpt->GetYaxis()->SetTitle("Normalized Counts");
	Gpt->GetYaxis()->SetTitleOffset(1.4);
	Gpt->GetXaxis()->CenterTitle();
	Gpt->GetYaxis()->CenterTitle();

	Gpt->SetMarkerStyle(20);
	Gpt->SetMarkerSize(1);
	Gpt->SetMarkerColor(1);
	Gpt->SetLineColor(1);


	TH1D * JPsiPt = new TH1D("JPsiPt","After reweighting",100,0,25);
	JPsiPt->GetXaxis()->SetTitle("Generated J/#psi p_{T}");
	JPsiPt->GetYaxis()->SetTitle("Normalized Counts");
	JPsiPt->GetYaxis()->SetTitleOffset(1.4);
	JPsiPt->GetXaxis()->CenterTitle();
	JPsiPt->GetYaxis()->CenterTitle();


	JPsiPt->SetMarkerStyle(20);
	JPsiPt->SetMarkerSize(1);
	JPsiPt->SetMarkerColor(1);
	JPsiPt->SetLineColor(1);


	TH1D * pthat = new TH1D("pthat","After reweighting",100,0,100);
	pthat->GetXaxis()->SetTitle("pthat");
	pthat->GetYaxis()->SetTitle("Normalized Counts");
	pthat->GetYaxis()->SetTitleOffset(1.4);
	pthat->GetXaxis()->CenterTitle();
	pthat->GetYaxis()->CenterTitle();


	pthat->SetMarkerStyle(20);
	pthat->SetMarkerSize(1);
	pthat->SetMarkerColor(1);
	pthat->SetLineColor(1);

	// TH1D * Gpt_noweight = new TH1D("Gpt_nw","Before reweighting",100,0,100);
  // Gpt->Copy(*Gpt_noweight);
	TH1D * Gpt_noweight = (TH1D*) Gpt->Clone("Gpt_nw");
  Gpt_noweight->SetTitle("Before reweighting");
  Gpt_noweight->SetLineColor(kRed);
  Gpt_noweight->SetMarkerColor(kRed);
	// TH1D * JPsiPt_noweight = new TH1D("JPsiPt_nw","Before reweighting",100,0,25);
  // JPsiPt->Copy(*JPsiPt_noweight);
	TH1D * JPsiPt_noweight = (TH1D*) JPsiPt->Clone("JPsiPt_nw");
  JPsiPt_noweight->SetTitle("Before reweighting");
  JPsiPt_noweight->SetLineColor(kRed);
  JPsiPt_noweight->SetMarkerColor(kRed);
	// TH1D * pthat_noweight = new TH1D("pthat_nw","Before reweighting",100,0,100);
  // pthat->Copy(*pthat_noweight);
	TH1D * pthat_noweight = (TH1D*) pthat->Clone("pthat_nw");
  pthat_noweight->SetTitle("Before reweighting");
  pthat_noweight->SetLineColor(kRed);
  pthat_noweight->SetMarkerColor(kRed);

	TCut Weight = "weight"; 
	TCut GenCut;
	TCut JPsiGenCut;


	if(Opt == 0) GenCut = "(TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0)";
	if(Opt == 1) GenCut = "(TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==531 && GisSignal>0)";



	if(Opt == 0) JPsiGenCut = "(TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==443)";
	if(Opt == 1) JPsiGenCut = "(TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==443)";

  // sanity check

	TH1D * weight = new TH1D("hweight","weight",40,0.8,1.2);
  ntGen->Draw("weight >> hweight", "weight < 1.4 && weight > 0.6");
  c->SaveAs("weight.png");
  return;

  // re-weighted
	// ntGen->Project("Gpt","Gpt",TCut(Weight) * TCut(GenCut));
	// ntGen->Project("JPsiPt","Gpt",TCut(Weight) * TCut(JPsiGenCut));
	// ntGen->Project("pthat","pthat",TCut(Weight));
  ntGen->Draw("pthat >> pthat", Weight);
  // unweighted
	// ntGen->Project("Gpt_nw","Gpt", TCut(GenCut));
	// ntGen->Project("JPsiPt_nw","Gpt", TCut(JPsiGenCut));
	// ntGen->Project("pthat_nw","pthat");
  ntGen->Draw("pthat >> pthat_nw");

	c->SetLogy();


	pthat_noweight->Scale(1.0/pthat_noweight->Integral());
	pthat_noweight->Draw("ep");
	c->SaveAs(Form("MCPlots/%spthat_bf.png",Part.Data()));	
	pthat->Scale(1.0/pthat->Integral());
	pthat->Draw("ep");
	c->SaveAs(Form("MCPlots/%spthat_af.png",Part.Data()));	
	// pthat_noweight->Scale(1.0/pthat_noweight->Integral());
	pthat_noweight->Draw("epsame");
  c->BuildLegend();
	c->SaveAs(Form("MCPlots/%spthat.png",Part.Data()));	


	Gpt->Scale(1.0/Gpt->Integral());
	Gpt->Draw("ep");
	Gpt_noweight->Scale(1.0/Gpt_noweight->Integral());
	Gpt_noweight->Draw("epsame");
  c->BuildLegend();
	c->SaveAs(Form("MCPlots/%sGpt.png",Part.Data()));	


	JPsiPt->Scale(1.0/JPsiPt->Integral());
	JPsiPt->Draw("ep");
	JPsiPt_noweight->Scale(1.0/JPsiPt_noweight->Integral());
	JPsiPt_noweight->Draw("epsame");
  c->BuildLegend();
	c->SaveAs(Form("MCPlots/%sJPsiPt.png",Part.Data()));	
	


}
