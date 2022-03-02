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


using namespace std;

using std::cout;
using std::endl;

void PlotEffSyst(int Opt){



	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	gStyle->SetOptStat(0);

	TString infile;


	if(Opt == 0) infile =  "../BP/EffAna/NewEff2DMaps/BPSyst.root";
	if(Opt == 1) infile =  "../Bs/EffAna/NewEff2DMaps/BsSyst.root";


	TString BmesonName;



	if(Opt == 0) BmesonName =  "BP";
	if(Opt == 1) BmesonName =  "Bs";


	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TH1D * Eff1DHis = (TH1D * ) fin->Get("Eff1DHis");
	TH1D * Eff1DHisTnPUp = (TH1D * ) fin->Get("Eff1DHisTnPUp");
	TH1D * Eff1DHisBpt = (TH1D * ) fin->Get("Eff1DHisBpt");
	TH1D * Eff1DHisBDT = (TH1D * ) fin->Get("Eff1DHisBDT");
	
	
	TH1D * TnPSyst = (TH1D *) Eff1DHisTnPUp->Clone("TnPSyst");
	TnPSyst->GetYaxis()->SetTitle("TnP Systematic Uncertainties (%)");
	TnPSyst->GetYaxis()->SetTitleOffset(1.3);

	TnPSyst->SetLineColor(kBlack);
	TnPSyst->SetMarkerColor(kBlack);	
	TnPSyst->Reset();


	TH1D * BptSyst = (TH1D *) Eff1DHisBpt->Clone("BptSyst");
	BptSyst->GetYaxis()->SetTitle("B-meson p_{T} Shape Systematic Uncertainties (%)");
	BptSyst->GetYaxis()->SetTitleOffset(1.3);
	
	BptSyst->SetLineColor(kBlack);
	BptSyst->SetMarkerColor(kBlack);	
	BptSyst->Reset();

	TH1D * BDTSyst = (TH1D *) Eff1DHisBDT->Clone("BDTSyst");
	BDTSyst->GetYaxis()->SetTitle("MC/Data Discrepancy Systematic Uncertainties (%)");
	BDTSyst->GetYaxis()->SetTitleOffset(1.3);
	
	BDTSyst->SetLineColor(kBlack);
	BDTSyst->SetMarkerColor(kBlack);	
	BDTSyst->Reset();








	float SystValue;
	float SystValueError;

	for(int i = 0; i < Eff1DHisTnPUp->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisTnPUp->GetBinContent(i+1) -  Eff1DHis->GetBinContent(i+1) )/Eff1DHis->GetBinContent(i+1) * 100;




		cout << "TnP Syst: " << SystValue << endl;
		
		TnPSyst->SetBinContent(i+1,SystValue);
		TnPSyst->SetBinError(i+1,0.001);

		

	}




	for(int i = 0; i < Eff1DHisTnPUp->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisBpt->GetBinContent(i+1) -  Eff1DHis->GetBinContent(i+1) )/Eff1DHis->GetBinContent(i+1) * 100;

		BptSyst->SetBinContent(i+1,SystValue);
		BptSyst->SetBinError(i+1,0.001);

		cout << "Bpt Syst: " << SystValue << endl;

			

	}




	for(int i = 0; i < Eff1DHisTnPUp->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisBDT->GetBinContent(i+1) -  Eff1DHis->GetBinContent(i+1) )/Eff1DHis->GetBinContent(i+1) * 100;

		BDTSyst->SetBinContent(i+1,SystValue);
		BDTSyst->SetBinError(i+1,0.001);

		cout << "BDT Syst: " << SystValue << endl;

			

	}











	TnPSyst->Draw("ep");

	c->SaveAs(Form("SystPlots/Pt/TnPSystRatio%s.png",BmesonName.Data()));
	
	BptSyst->Draw("ep");

	c->SaveAs(Form("SystPlots/Pt/BptSysRatiot%s.png",BmesonName.Data()));

	BDTSyst->Draw("ep");

	c->SaveAs(Form("SystPlots/Pt/BDTSystRatio%s.png",BmesonName.Data()));

	


	TH1D * Eff1DHisMult = (TH1D * ) fin->Get("Eff1DHisMult");
	TH1D * Eff1DHisTnPUpMult = (TH1D * ) fin->Get("Eff1DHisTnPUpMult");
	TH1D * Eff1DHisBptMult = (TH1D * ) fin->Get("Eff1DHisBptMult");
	TH1D * Eff1DHisBDTMult = (TH1D * ) fin->Get("Eff1DHisBDTMult");



	
	TH1D * TnPSystMult = (TH1D *) Eff1DHisTnPUp->Clone("TnPSyst");
	TnPSystMult->GetYaxis()->SetTitle("TnP Systematic Uncertainties (%)");
	TnPSystMult->GetYaxis()->SetTitleOffset(1.3);
	
	TnPSystMult->SetLineColor(kBlack);
	TnPSystMult->SetMarkerColor(kBlack);	
	TnPSystMult->Reset();


	TH1D * BptSystMult = (TH1D *) Eff1DHisBpt->Clone("BptSystMult");
	BptSystMult->GetYaxis()->SetTitle("B-meson p_{T} Shape Systematic Uncertainties (%)");
	BptSystMult->GetYaxis()->SetTitleOffset(1.3);
	
	BptSystMult->SetLineColor(kBlack);
	BptSystMult->SetMarkerColor(kBlack);	
	BptSystMult->Reset();

	TH1D * BDTSystMult = (TH1D *) Eff1DHisBDT->Clone("BDTSystMult");
	BDTSystMult->GetYaxis()->SetTitle("BDT Systematic Uncertainties (%)");
	BDTSystMult->GetYaxis()->SetTitleOffset(1.3);

	BDTSystMult->SetLineColor(kBlack);
	BDTSystMult->SetMarkerColor(kBlack);	
	BDTSystMult->Reset();



	for(int i = 0; i < Eff1DHisTnPUpMult->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisTnPUpMult->GetBinContent(i+1) -  Eff1DHisMult->GetBinContent(i+1) )/Eff1DHisMult->GetBinContent(i+1) * 100;




		cout << "TnP Syst: " << SystValue << endl;
		
		TnPSystMult->SetBinContent(i+1,SystValue);
		TnPSystMult->SetBinError(i+1,0.001);

		

	}




	for(int i = 0; i < Eff1DHisTnPUpMult->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisBptMult->GetBinContent(i+1) -  Eff1DHisMult->GetBinContent(i+1) )/Eff1DHisMult->GetBinContent(i+1) * 100;

		BptSystMult->SetBinContent(i+1,SystValue);
		BptSystMult->SetBinError(i+1,0.001);

		cout << "Bpt Syst: " << SystValue << endl;

			

	}




	for(int i = 0; i < Eff1DHisTnPUpMult->GetNbinsX(); i++){

		SystValue = abs(Eff1DHisBDTMult->GetBinContent(i+1) -  Eff1DHisMult->GetBinContent(i+1) )/Eff1DHisMult->GetBinContent(i+1) * 100;

		BDTSystMult->SetBinContent(i+1,SystValue);
		BDTSystMult->SetBinError(i+1,0.001);

		cout << "BDT Syst: " << SystValue << endl;

			

	}








	TnPSystMult->Draw("ep");

	c->SaveAs(Form("SystPlots/Mult/TnPSystRatio%sMult.png",BmesonName.Data()));



	BptSystMult->Draw("ep");
	c->SaveAs(Form("SystPlots/Mult/BptSystRatio%sMult.png",BmesonName.Data()));




	BDTSystMult->Draw("ep");
	c->SaveAs(Form("SystPlots/Mult/BDTSystRatio%sMult.png",BmesonName.Data()));



}


