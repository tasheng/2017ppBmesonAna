#ifndef __CINT__
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
#endif


void PVzQuickPlot(){


	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();

	gStyle->SetOptStat(0);

	TFile * fin = new TFile("/data/szhaozho/2017ppSamplesNew/BDTOutput/AllMerge/BPDataAllBDT.root");

	fin->cd();

	TTree * ntKp = (TTree *) fin->Get("Bfinder/ntKp");

	ntKp->AddFriend("skimanalysis/HltTree");


	TH1D * DataPVzHis = new TH1D("DataPVzHis","",300,-30,30);

	ntKp->Project("DataPVzHis","PVz","HBHENoiseFilterResult == 1 && pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1");
	

	DataPVzHis->GetXaxis()->SetTitle("PV_{z} (cm)");
	DataPVzHis->GetYaxis()->SetTitle("Number of Events");
	DataPVzHis->GetXaxis()->CenterTitle();
	DataPVzHis->GetYaxis()->CenterTitle();
	DataPVzHis->GetYaxis()->SetTitleOffset(1.2);
	
	DataPVzHis->Draw("hist");
	c->SaveAs("DataPVz.png");

	TFile * finMC = new TFile("/data/szhaozho/2017ppSamplesNew/BDTOutput/AllMerge/BPMCAllBDT.root");

	finMC->cd();

	TTree * ntKpMC = (TTree *) finMC->Get("Bfinder/ntKp");

	ntKpMC->AddFriend("skimanalysis/HltTree");
	

	TH1D * MCPVzHis = new TH1D("MCPVzHis","",300,-30,30);

	ntKpMC->Project("MCPVzHis","PVz","HBHENoiseFilterResult == 1 && pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1");
	

	MCPVzHis->GetXaxis()->SetTitle("PV_{z} (cm)");
	MCPVzHis->GetYaxis()->SetTitle("Number of Events");
	MCPVzHis->GetXaxis()->CenterTitle();
	MCPVzHis->GetYaxis()->CenterTitle();
	MCPVzHis->GetYaxis()->SetTitleOffset(1.2);

	
	MCPVzHis->Draw("hist");
	c->SaveAs("MCPVz.png");


}
