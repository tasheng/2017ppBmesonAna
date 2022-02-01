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

void All2DMap(){



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


	TFile * fin1 = new TFile("NewEff2DMaps/EffFine_0_90.root");
	fin1->cd();
	TH2D * EffBptByInv1 = (TH2D *) fin1->Get("invEff2D");
	cAll->cd(1);
//	EffBptByInv1->Divide(tnp_scale1);

	EffBptByInv1->GetYaxis()->SetTitle("B^{+} |y|");
	EffBptByInv1->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	EffBptByInv1->GetXaxis()->CenterTitle();
	EffBptByInv1->GetYaxis()->CenterTitle();
	EffBptByInv1->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv1->SetTitle("1/(Efficiency x Acceptance) 2D Map at Centrality 0 - 90");


	EffBptByInv1->Draw("COLZ");
		





	TFile * fin2 = new TFile("NewEff2DMaps/EffFine_0_30.root");
	fin2->cd();
	TH2D * EffBptByInv2 = (TH2D *) fin2->Get("invEff2D");
	cAll->cd(2);
//EffBptByInv2->Divide(tnp_scale2);

	EffBptByInv2->GetYaxis()->SetTitle("B^{+} |y|");
	EffBptByInv2->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	EffBptByInv2->GetXaxis()->CenterTitle();
	EffBptByInv2->GetYaxis()->CenterTitle();
	EffBptByInv2->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv2->SetTitle("1/(Efficiency x Acceptance) 2D Map at Centrality 0 - 30");


	EffBptByInv2->Draw("COLZ");


	TFile * fin3 = new TFile("NewEff2DMaps/EffFine_30_90.root");
	cAll->cd(3);

	TH2D * EffBptByInv3 = (TH2D *) fin3->Get("invEff2D");
//	EffBptByInv3->Divide(tnp_scale3);


	EffBptByInv3->GetYaxis()->SetTitle("B^{+} |y|");
	EffBptByInv3->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	EffBptByInv3->GetXaxis()->CenterTitle();
	EffBptByInv3->GetYaxis()->CenterTitle();
	EffBptByInv3->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv3->SetTitle("1/(Efficiency x Acceptance) 2D Map at Centrality 30 - 90");

	EffBptByInv3->Draw("COLZ");
	cAll->SaveAs("BP2DMaps/Eff_Corr.png");




	//Selection Eff//
	TH2D * hEffSelectionFine1 = (TH2D *) fin1->Get("invEffonly2D");
	cAll->cd(1);

//	hEffSelectionFine1->Divide(tnp_scale1);

	hEffSelectionFine1->GetYaxis()->SetTitle("B^{+} |y|");
	hEffSelectionFine1->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	hEffSelectionFine1->GetXaxis()->CenterTitle();
	hEffSelectionFine1->GetYaxis()->CenterTitle();
	hEffSelectionFine1->GetYaxis()->SetTitleOffset(1.2);
	hEffSelectionFine1->SetTitle("1/Efficiency 2D Map at Centrality 0 - 90");


	hEffSelectionFine1->Draw("COLZ");

	TH2D * hEffSelectionFine2 = (TH2D *) fin2->Get("invEffonly2D");
//	hEffSelectionFine2->Divide(tnp_scale2);
	cAll->cd(2);

	hEffSelectionFine2->GetYaxis()->SetTitle("B^{+} |y|");
	hEffSelectionFine2->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	hEffSelectionFine2->GetXaxis()->CenterTitle();
	hEffSelectionFine2->GetYaxis()->CenterTitle();
	hEffSelectionFine2->GetYaxis()->SetTitleOffset(1.2);
	hEffSelectionFine2->SetTitle("1/Efficiency 2D Map at Centrality 0 - 30");


	hEffSelectionFine2->Draw("COLZ");


	TH2D * hEffSelectionFine3 = (TH2D *) fin3->Get("invEffonly2D");
//	hEffSelectionFine3->Divide(tnp_scale3);
	cAll->cd(3);

	hEffSelectionFine3->GetYaxis()->SetTitle("B^{+} |y|");
	hEffSelectionFine3->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	hEffSelectionFine3->GetXaxis()->CenterTitle();
	hEffSelectionFine3->GetYaxis()->CenterTitle();
	hEffSelectionFine3->GetYaxis()->SetTitleOffset(1.2);
	hEffSelectionFine3->SetTitle("1/Efficiency 2D Map at Centrality 30 - 90");

	hEffSelectionFine3->Draw("COLZ");
	cAll->SaveAs("BP2DMaps/Sel_Corr.png");

	//Selection Acc//

	TH2D * hEffAcc2DFine1 = (TH2D *) fin1->Get("invAcc2D");

	hEffAcc2DFine1->GetYaxis()->SetTitle("B^{+} |y|");
	hEffAcc2DFine1->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	hEffAcc2DFine1->GetXaxis()->CenterTitle();
	hEffAcc2DFine1->GetYaxis()->CenterTitle();
	hEffAcc2DFine1->GetYaxis()->SetTitleOffset(1.2);
	hEffAcc2DFine1->SetTitle("1/Acceptance 2D Map at Centrality 0 - 90");


	hEffAcc2DFine1->Draw("COLZ");

	


	TH2D * hEffAcc2DFine2 = (TH2D *) fin2->Get("invAcc2D");

	hEffAcc2DFine2->GetYaxis()->SetTitle("B^{+} |y|");
	hEffAcc2DFine2->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	hEffAcc2DFine2->GetXaxis()->CenterTitle();
	hEffAcc2DFine2->GetYaxis()->CenterTitle();
	hEffAcc2DFine2->GetYaxis()->SetTitleOffset(1.2);
	hEffAcc2DFine2->SetTitle("1/Acceptance 2D Map at Centrality 0 - 30");


	hEffAcc2DFine2->Draw("COLZ");



	TH2D * hEffAcc2DFine3 = (TH2D *) fin3->Get("invAcc2D");

	hEffAcc2DFine3->GetYaxis()->SetTitle("B^{+} |y|");
	hEffAcc2DFine3->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	hEffAcc2DFine3->GetXaxis()->CenterTitle();
	hEffAcc2DFine3->GetYaxis()->CenterTitle();
	hEffAcc2DFine3->GetYaxis()->SetTitleOffset(1.2);
	hEffAcc2DFine3->SetTitle("1/Acceptance 2D Map at Centrality 30 - 90");

	hEffAcc2DFine3->Draw("COLZ");
	cAll->SaveAs("BP2DMaps/Acc_Corr.png");



	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();
	c->SetLogz();
	
	EffBptByInv1->Draw("COLZ");
	c->SaveAs("BP2DMaps/SingleBPApp/Eff2D_0_90.png");

	hEffSelectionFine1->Draw("COLZ");
	c->SaveAs("BP2DMaps/SingleBPApp/Sel2D_0_90.png");
	

	hEffAcc2DFine1->Draw("COLZ");
	c->SaveAs("BP2DMaps/SingleBPApp/Acc2D_0_90.png");

	
	EffBptByInv2->Draw("COLZ");
	c->SaveAs("BP2DMaps/SingleBPApp/Eff2D_0_30.png");
	
	EffBptByInv3->Draw("COLZ");
	c->SaveAs("BP2DMaps/SingleBPApp/Eff2D_30_90.png");



	TH2D * EffBptByInv1Up = (TH2D *) fin1->Get("invEff2DSystUp");
	EffBptByInv1Up->GetYaxis()->SetTitle("B^{+} |y|");
	EffBptByInv1Up->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	EffBptByInv1Up->GetXaxis()->CenterTitle();
	EffBptByInv1Up->GetYaxis()->CenterTitle();
	EffBptByInv1Up->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv1Up->SetTitle("1/(Efficiency x Acceptance) 2D Map TnP Syst Up at Centrality 0 - 90");
	EffBptByInv1Up->Draw("COLZ");
	c->SaveAs("BP2DMaps/SingleBPApp/Eff2D_0_90_Up.png");


	TH2D * EffBptByInv1Down = (TH2D *) fin1->Get("invEff2DSystDown");
	EffBptByInv1Down->GetYaxis()->SetTitle("B^{+} |y|");
	EffBptByInv1Down->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	EffBptByInv1Down->GetXaxis()->CenterTitle();
	EffBptByInv1Down->GetYaxis()->CenterTitle();
	EffBptByInv1Down->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv1Down->SetTitle("1/(Efficiency x Acceptance) 2D Map TnP Syst Up at Centrality 0 - 90");
	EffBptByInv1Down->Draw("COLZ");
	c->SaveAs("BP2DMaps/SingleBPApp/Eff2D_0_90_Down.png");

	//0 - 30%//

	TH2D * EffBptByInv2Up = (TH2D *) fin2->Get("invEff2DSystUp");
	EffBptByInv2Up->GetYaxis()->SetTitle("B^{+} |y|");
	EffBptByInv2Up->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	EffBptByInv2Up->GetXaxis()->CenterTitle();
	EffBptByInv2Up->GetYaxis()->CenterTitle();
	EffBptByInv2Up->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv2Up->SetTitle("1/(Efficiency x Acceptance) 2D Map TnP Syst Up at Centrality 0 - 30");
	EffBptByInv2Up->Draw("COLZ");
	c->SaveAs("BP2DMaps/SingleBPApp/Eff2D_0_30_Up.png");


	TH2D * EffBptByInv2Down = (TH2D *) fin2->Get("invEff2DSystDown");
	EffBptByInv2Down->GetYaxis()->SetTitle("B^{+} |y|");
	EffBptByInv2Down->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	EffBptByInv2Down->GetXaxis()->CenterTitle();
	EffBptByInv2Down->GetYaxis()->CenterTitle();
	EffBptByInv2Down->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv2Down->SetTitle("1/(Efficiency x Acceptance) 2D Map TnP Syst Up at Centrality 0 - 30");
	EffBptByInv2Down->Draw("COLZ");
	c->SaveAs("BP2DMaps/SingleBPApp/Eff2D_0_30_Down.png");


	//30 - 90%//

	TH2D * EffBptByInv3Up = (TH2D *) fin3->Get("invEff2DSystUp");
	EffBptByInv3Up->GetYaxis()->SetTitle("B^{+} |y|");
	EffBptByInv3Up->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	EffBptByInv3Up->GetXaxis()->CenterTitle();
	EffBptByInv3Up->GetYaxis()->CenterTitle();
	EffBptByInv3Up->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv3Up->SetTitle("1/(Efficiency x Acceptance) 2D Map TnP Syst Up at Centrality 30 - 90");
	EffBptByInv3Up->Draw("COLZ");
	c->SaveAs("BP2DMaps/SingleBPApp/Eff2D_30_90_Up.png");


	TH2D * EffBptByInv3Down = (TH2D *) fin3->Get("invEff2DSystDown");
	EffBptByInv3Down->GetYaxis()->SetTitle("B^{+} |y|");
	EffBptByInv3Down->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	EffBptByInv3Down->GetXaxis()->CenterTitle();
	EffBptByInv3Down->GetYaxis()->CenterTitle();
	EffBptByInv3Down->GetYaxis()->SetTitleOffset(1.2);
	EffBptByInv3Down->SetTitle("1/(Efficiency x Acceptance) 2D Map TnP Syst Up at Centrality 30 - 90");
	EffBptByInv3Down->Draw("COLZ");
	c->SaveAs("BP2DMaps/SingleBPApp/Eff2D_30_90_Down.png");



}
