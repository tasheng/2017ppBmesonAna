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


void PVZMC(int Opt){
	
	TString OutfileName;
	TString DataInfile;
	TString MCInfile;
	TString TTreeName;
	
	if(Opt == 0){
			
		 OutfileName = "BPPVZMCData.png";
		 DataInfile =  "../../Unskimmed/BPData.root";
		 MCInfile =  "../../Unskimmed/NewOfficialMC/BPMC.root";
		 TTreeName = "Bfinder/ntKp";

	}


	if(Opt == 1){

		 OutfileName = "BsPVZMCData.png";
		 DataInfile =  "../../Unskimmed/BsData.root";
		 MCInfile =  "../../Unskimmed/NewOfficialMC/BsMC.root";
		 TTreeName = "Bfinder/ntphi";

	}




	TCanvas * c = new TCanvas("c","c",1800,600);
	c->Divide(3,1);
	c->cd(1);
	gPad->SetLeftMargin(0.17);


	gStyle->SetOptStat(0);

	TFile * fin = new TFile(DataInfile.Data());

	fin->cd();

	TTree * ntKp = (TTree *) fin->Get(TTreeName.Data());

	ntKp->AddFriend("skimanalysis/HltTree");


	TH1D * DataPVzHis = new TH1D("DataPVzHis","",300,-30,30);

	ntKp->Project("DataPVzHis","PVz","pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1");

	DataPVzHis->GetXaxis()->SetTitle("PV_{z} (cm)");
	DataPVzHis->GetYaxis()->SetTitle("Number of Events");
	DataPVzHis->SetTitle("Data");
	
	DataPVzHis->GetXaxis()->CenterTitle();
	DataPVzHis->GetYaxis()->CenterTitle();
	DataPVzHis->GetYaxis()->SetTitleOffset(1.8);
	DataPVzHis->Scale(1.0/DataPVzHis->Integral());

	DataPVzHis->SetMarkerStyle(20);
	DataPVzHis->SetMarkerSize(1);
	DataPVzHis->SetMarkerColor(kBlue);
	DataPVzHis->SetLineColor(kBlue);
	

	

	DataPVzHis->Sumw2();
	
	TH1D * DataPVzHis2 = (TH1D *) DataPVzHis->Clone("DataPVzHis2");
	
	DataPVzHis->Draw("ep");

	TF1 * FittoData = new TF1("f1","gaus",-30,30);
	DataPVzHis->Fit(FittoData,"R");
	
	DataPVzHis->SetMaximum(0.015);

	DataPVzHis->Draw("ep");
	FittoData->Draw("SAME");

	float p0 =FittoData->GetParameter(0);
	float p1 =FittoData->GetParameter(1);
	float p2 =FittoData->GetParameter(2);


	c->cd(2);
	gPad->SetLeftMargin(0.18);

	
	TFile * finMC = new TFile(MCInfile.Data());

	finMC->cd();

	TTree * ntKpMC = (TTree *) finMC->Get(TTreeName.Data());

	ntKpMC->AddFriend("skimanalysis/HltTree");
	

	TH1D * MCPVzHis = new TH1D("MCPVzHis","",300,-30,30);

	ntKpMC->Project("MCPVzHis","PVz","pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1");
	
	MCPVzHis->GetXaxis()->SetTitle("PV_{z} (cm)");
	MCPVzHis->GetYaxis()->SetTitle("Number of Events");
	MCPVzHis->SetTitle("MC");
	MCPVzHis->GetXaxis()->CenterTitle();
	MCPVzHis->GetYaxis()->CenterTitle();
	MCPVzHis->GetYaxis()->SetTitleOffset(1.8);
	MCPVzHis->Scale(1.0/MCPVzHis->Integral());


	MCPVzHis->SetMarkerStyle(20);
	MCPVzHis->SetMarkerSize(1);
	MCPVzHis->SetMarkerColor(kGreen);
	MCPVzHis->SetLineColor(kGreen);

	MCPVzHis->Sumw2();
	MCPVzHis->Draw("ep");
	TF1 * FittoMC = new TF1("f1","gaus",-30,30);
	MCPVzHis->Fit(FittoMC,"R");


		
	MCPVzHis->SetMaximum(0.015);
	MCPVzHis->Draw("ep");
	FittoMC->Draw("SAME");
	

	float p3 =FittoMC->GetParameter(0);
	float p4 =FittoMC->GetParameter(1);
	float p5 =FittoMC->GetParameter(2);



	c->cd(3);
	gPad->SetLeftMargin(0.18);

	TString FinalFunc =  Form("(%f * TMath::Exp(-(PVz-%f)*(PVz-%f)/(2 * %f * %f)))/(%f * TMath::Exp(-(PVz-%f)*(PVz-%f)/(2 * %f * %f)))",p0,p1,p1,p2,p2,p3,p4,p4,p5,p5);

	TH1D * MCPVzHisWerighted = new TH1D("MCPVzHisWerighted","",300,-30,30);
	ntKpMC->Project("MCPVzHisWerighted","PVz",Form("(pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1)*(%s)",FinalFunc.Data()));

	MCPVzHisWerighted->GetXaxis()->SetTitle("PV_{z} (cm)");
	MCPVzHisWerighted->GetYaxis()->SetTitle("Number of Events");
	MCPVzHisWerighted->SetTitle("Reweighted MC + Data");
	
	MCPVzHisWerighted->GetXaxis()->CenterTitle();
	MCPVzHisWerighted->GetYaxis()->CenterTitle();
	MCPVzHisWerighted->GetYaxis()->SetTitleOffset(1.8);
	MCPVzHisWerighted->Scale(1.0/MCPVzHisWerighted->Integral());

	MCPVzHisWerighted->SetMarkerStyle(20);
	MCPVzHisWerighted->SetMarkerSize(1);
	MCPVzHisWerighted->SetMarkerColor(kGreen);
	MCPVzHisWerighted->SetLineColor(kGreen);

	MCPVzHisWerighted->Sumw2();
	
	MCPVzHisWerighted->SetMaximum(0.015);

	MCPVzHisWerighted->Draw("ep");

	DataPVzHis2->Draw("epSAME");
	
	c->SaveAs(OutfileName.Data());
	
  std::cout << FinalFunc << "\n";


}
