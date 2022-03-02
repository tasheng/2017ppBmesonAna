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
//#include "parameters.h"


//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void Quick2DMap(){



	gStyle->SetOptStat(0);
	gStyle->SetPadRightMargin(10);

	TCanvas * cAll = new TCanvas("c","c",1800,600);
	cAll->SetLogz();
	cAll->Divide(3,1);
	TPad *p1 = cAll->cd(1); 
	p1->SetLogz();
	TPad *p2 = cAll->cd(2); 
	p2->SetLogz();
	TPad *p3 = cAll->cd(3); 
	p3->SetLogz();

/*
	
	TFile * finTnP1 =  new TFile(Form("TnP/TNP2D_Bs_Cent30-90.root"));
	TFile * finTnP2 =  new TFile(Form("TnP/TNP2D_Bs_Cent30-90.root"));
	TFile * finTnP3 =  new TFile(Form("TnP/TNP2D_Bs_Cent30-90.root"));

	TH2D * tnp_scale1 = (TH2D *) finTnP1->Get("tnp_scale");
	TH2D * tnp_scale2 = (TH2D *) finTnP1->Get("tnp_scale");
	TH2D * tnp_scale3 = (TH2D *) finTnP1->Get("tnp_scale");

*/


	TFile * fin1 = new TFile("2DMaps/EffFine_0_90.root");
	fin1->cd();
	TH2D * EffBptByInv1 = (TH2D *) fin1->Get("hEff2DInv2Shots");
	cAll->cd(1);
//	EffBptByInv1->Divide(tnp_scale1);

	EffBptByInv1->GetYaxis()->SetTitle("B |y|");
	EffBptByInv1->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
	EffBptByInv1->GetXaxis()->CenterTitle();
	EffBptByInv1->GetYaxis()->CenterTitle();
	EffBptByInv1->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv1->SetTitle("1/(Efficiency x Acceptance) 2D Map at Centrality 0 - 90");


	EffBptByInv1->Draw("COLZ");
		





	TFile * fin2 = new TFile("2DMaps/EffFine_0_30.root");
	fin2->cd();
	TH2D * EffBptByInv2 = (TH2D *) fin2->Get("hEff2DInv2Shots");
	cAll->cd(2);
//EffBptByInv2->Divide(tnp_scale2);

	EffBptByInv2->GetYaxis()->SetTitle("B |y|");
	EffBptByInv2->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
	EffBptByInv2->GetXaxis()->CenterTitle();
	EffBptByInv2->GetYaxis()->CenterTitle();
	EffBptByInv2->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv2->SetTitle("1/(Efficiency x Acceptance) 2D Map at Centrality 0 - 30");


	EffBptByInv2->Draw("COLZ");


	TFile * fin3 = new TFile("2DMaps/EffFine_30_90.root");
	cAll->cd(3);

	TH2D * EffBptByInv3 = (TH2D *) fin3->Get("hEff2DInv2Shots");
//	EffBptByInv3->Divide(tnp_scale3);


	EffBptByInv3->GetYaxis()->SetTitle("B |y|");
	EffBptByInv3->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
	EffBptByInv3->GetXaxis()->CenterTitle();
	EffBptByInv3->GetYaxis()->CenterTitle();
	EffBptByInv3->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv3->SetTitle("1/(Efficiency x Acceptance) 2D Map at Centrality 30 - 90");

	EffBptByInv3->Draw("COLZ");
	cAll->SaveAs("Eff2DMap/Eff_NoCorr.png");




	//Selection Eff//
	TH2D * hEffSelectionFine1 = (TH2D *) fin1->Get("hEffSelectionFine");
	cAll->cd(1);

//	hEffSelectionFine1->Divide(tnp_scale1);

	hEffSelectionFine1->GetYaxis()->SetTitle("B |y|");
	hEffSelectionFine1->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
	hEffSelectionFine1->GetXaxis()->CenterTitle();
	hEffSelectionFine1->GetYaxis()->CenterTitle();
	hEffSelectionFine1->GetYaxis()->SetTitleOffset(1.2);
	hEffSelectionFine1->SetTitle("1/Efficiency 2D Map at Centrality 0 - 90");


	hEffSelectionFine1->Draw("COLZ");

	TH2D * hEffSelectionFine2 = (TH2D *) fin2->Get("hEffSelectionFine");
//	hEffSelectionFine2->Divide(tnp_scale2);
	cAll->cd(2);

	hEffSelectionFine2->GetYaxis()->SetTitle("B |y|");
	hEffSelectionFine2->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
	hEffSelectionFine2->GetXaxis()->CenterTitle();
	hEffSelectionFine2->GetYaxis()->CenterTitle();
	hEffSelectionFine2->GetYaxis()->SetTitleOffset(1.2);
	hEffSelectionFine2->SetTitle("1/Efficiency 2D Map at Centrality 0 - 30");


	hEffSelectionFine2->Draw("COLZ");


	TH2D * hEffSelectionFine3 = (TH2D *) fin3->Get("hEffSelectionFine");
//	hEffSelectionFine3->Divide(tnp_scale3);
	cAll->cd(3);

	hEffSelectionFine3->GetYaxis()->SetTitle("B |y|");
	hEffSelectionFine3->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
	hEffSelectionFine3->GetXaxis()->CenterTitle();
	hEffSelectionFine3->GetYaxis()->CenterTitle();
	hEffSelectionFine3->GetYaxis()->SetTitleOffset(1.2);
	hEffSelectionFine3->SetTitle("1/Efficiency 2D Map at Centrality 30 - 90");

	hEffSelectionFine3->Draw("COLZ");
	cAll->SaveAs("Eff2DMap/Sel_NoCorr.png");

	//Selection Acc//
	
	TFile * finAcc1 = new TFile("ZZGJComparison/AccNew.root");
	finAcc1->cd();
	cAll->cd(1);

	TH2D * hEffAcc2DFine1 = (TH2D *) finAcc1->Get("hEffAcc2DFine");

	hEffAcc2DFine1->GetYaxis()->SetTitle("B |y|");
	hEffAcc2DFine1->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
	hEffAcc2DFine1->GetXaxis()->CenterTitle();
	hEffAcc2DFine1->GetYaxis()->CenterTitle();
	hEffAcc2DFine1->GetYaxis()->SetTitleOffset(1.2);
	hEffAcc2DFine1->SetTitle("1/Acceptance 2D Map at Centrality 0 - 90");


	hEffAcc2DFine1->Draw("COLZ");

	
	TFile * finAcc2 = new TFile("ZZGJComparison/AccNew.root");
	finAcc2->cd();
	cAll->cd(2);


	TH2D * hEffAcc2DFine2 = (TH2D *) finAcc2->Get("hEffAcc2DFine");

	hEffAcc2DFine2->GetYaxis()->SetTitle("B |y|");
	hEffAcc2DFine2->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
	hEffAcc2DFine2->GetXaxis()->CenterTitle();
	hEffAcc2DFine2->GetYaxis()->CenterTitle();
	hEffAcc2DFine2->GetYaxis()->SetTitleOffset(1.2);
	hEffAcc2DFine2->SetTitle("1/Acceptance 2D Map at Centrality 0 - 30");


	hEffAcc2DFine2->Draw("COLZ");


	TFile * finAcc3 = new TFile("ZZGJComparison/AccNew.root");
	finAcc3->cd();

	cAll->cd(3);

	TH2D * hEffAcc2DFine3 = (TH2D *) finAcc3->Get("hEffAcc2DFine");

	hEffAcc2DFine3->GetYaxis()->SetTitle("B |y|");
	hEffAcc2DFine3->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
	hEffAcc2DFine3->GetXaxis()->CenterTitle();
	hEffAcc2DFine3->GetYaxis()->CenterTitle();
	hEffAcc2DFine3->GetYaxis()->SetTitleOffset(1.2);
	hEffAcc2DFine3->SetTitle("1/Acceptance 2D Map at Centrality 30 - 90");

	hEffAcc2DFine3->Draw("COLZ");
	cAll->SaveAs("Eff2DMap/Acc_NoCorr.png");



	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	c->SetLogz();
	
	EffBptByInv1->Draw("COLZ");
	c->SaveAs("Eff2DMap/SingleBsApp/Eff2D_0_90.png");

	hEffSelectionFine1->Draw("COLZ");
	c->SaveAs("Eff2DMap/SingleBsApp/Sel2D_0_90.png");
	
	
	hEffAcc2DFine1->Draw("COLZ");
	c->SaveAs("Eff2DMap/SingleBsApp/Acc2D_0_90.png");


	//Corrected//

	TFile * fin1C = new TFile("2DMaps/EffFine_0_90.root");
	fin1C->cd();
	TH2D * EffBptByInv1C = (TH2D *) fin1C->Get("EffBptByInv");
	
	EffBptByInv1C->GetYaxis()->SetTitle("B |y|");
	EffBptByInv1C->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
	EffBptByInv1C->GetXaxis()->CenterTitle();
	EffBptByInv1C->GetYaxis()->CenterTitle();
	EffBptByInv1C->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv1C->SetTitle("Efficiency 2D Map at Centrality 0 - 90");
	//EffBptByInv1C->GetZaxis()->SetTitleOffset(1.1);
	EffBptByInv1C->GetZaxis()->SetLabelOffset(1);
	EffBptByInv1C->Draw("COLZ");

	c->SaveAs("Eff2DMap/Eff_0_90_Corr.png");

	TFile * fin2C = new TFile("2DMaps/EffFine_0_30.root");
	fin2C->cd();
	TH2D * EffBptByInv2C = (TH2D *) fin2C->Get("EffBptByInv");
	
	EffBptByInv2C->GetYaxis()->SetTitle("B |y|");
	EffBptByInv2C->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
	EffBptByInv2C->GetXaxis()->CenterTitle();
	EffBptByInv2C->GetYaxis()->CenterTitle();
	EffBptByInv2C->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv2C->SetTitle("Efficiency 2D Map at Centrality 0 - 30");


	EffBptByInv2C->Draw("COLZ");
	c->SaveAs("Eff2DMap/Eff_0_30_Corr.png");


	TFile * fin3C = new TFile("2DMaps/EffFine_30_90.root");
	fin3->cd();
	TH2D * EffBptByInv3C = (TH2D *) fin3C->Get("EffBptByInv");
	
	EffBptByInv3C->GetYaxis()->SetTitle("B |y|");
	EffBptByInv3C->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
	EffBptByInv3C->GetXaxis()->CenterTitle();
	EffBptByInv3C->GetYaxis()->CenterTitle();
	EffBptByInv3C->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv3C->SetTitle("Efficiency 2D Map at Centrality 30 - 90");

	EffBptByInv3C->Draw("COLZ");
	c->SaveAs("Eff2DMap/Eff_30_90_Corr.png");


	TCanvas * c15 = new TCanvas("c15","c15",600,600);
	c15->cd();
	//Draw New B Cand Distribution//

	TFile * finNew = new TFile("CheckSystNunoNew/NoTnP/EffInfo_0_90.root");
	finNew->cd();
	TH2D * ByBptHis = new TH2D("ByBptHis","",45,5,50,24,0,2.4);
	TTree * t = (TTree *) finNew->Get("EffInfoTree");
	t->Project("ByBptHis","ByNew:BptNew");
	
	ByBptHis->GetXaxis()->SetTitle("B_{s}^{0} p_{T} (GeV/c)");
	ByBptHis->GetYaxis()->SetTitle("B_{s}^{0} |y|");
	ByBptHis->GetXaxis()->CenterTitle();
	ByBptHis->GetYaxis()->CenterTitle();
	ByBptHis->GetYaxis()->SetTitleOffset(1.4);
	ByBptHis->Draw("COLZ");
	
	c15->SaveAs("Eff2DMap/hBptByData.png");
	


	//DONE


	TFile * finNewMC = new TFile("CheckSystNunoNew/NoTnP/EffInfo_MC_Closure-1.root");
	finNewMC->cd();
	TH2D * ByBptHisMC = new TH2D("ByBptHisMC","",45,5,50,24,0,2.4);
	TTree * t = (TTree *) finNewMC->Get("EffInfoTree");
	t->Project("ByBptHisMC","ByNew:BptNew");
	
	ByBptHisMC->GetXaxis()->SetTitle("B_{s}^{0} p_{T} (GeV/c)");
	ByBptHisMC->GetYaxis()->SetTitle("B_{s}^{0} |y|");
	ByBptHisMC->GetXaxis()->CenterTitle();
	ByBptHisMC->GetYaxis()->CenterTitle();
	ByBptHisMC->GetYaxis()->SetTitleOffset(1.4);
	ByBptHisMC->Draw("COLZ");
	
	c15->SaveAs("Eff2DMap/hBptByMC.png");




}
