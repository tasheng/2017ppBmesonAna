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
//#include "tnp_weight_lowptPbPb.h"



//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void BsBP(){

	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);

	c->cd();

	TString InfileBs = "../Bs/EffAna/FinalFiles/BsPPCorrYieldPT.root";
	TString InfileBP = "../BP/EffAna/FinalFiles/BPPPCorrYieldPT.root";

	TFile * FileBs = new TFile(InfileBs.Data());

	TH1D * BsCross = (TH1D *) FileBs->Get("CorrDiffHisBin");
	BsCross->SetMarkerStyle(20);
	BsCross->SetMarkerSize(1);
	BsCross->SetMarkerColor(1);
	BsCross->SetLineColor(1);
	
	BsCross->Draw("ep");

	c->SaveAs("BsCross.png");

	TFile * FileBP = new TFile(InfileBP.Data());
	TH1D * BPCross = (TH1D *) FileBP->Get("CorrDiffHisBin");
	BPCross->SetMarkerStyle(20);
	BPCross->SetMarkerSize(1);
	BPCross->SetMarkerColor(1);
	BPCross->SetLineColor(1);
	
	BPCross->Draw("ep");

	c->SaveAs("BPCross.png");

	
	BsCross->SetLineColor(kBlue+2);
	BsCross->SetMarkerColor(kBlue+2);
	
	BPCross->SetLineColor(kGreen+2);
	BPCross->SetMarkerColor(kGreen+2);
	
	BPCross->Draw("ep");
	BsCross->Draw("epSAME");

	TLegend* leg2 = new TLegend(0.40,0.45,0.75,0.70,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.050);
	leg2->SetTextFont(42);
	leg2->SetFillStyle(0);
	leg2->SetLineWidth(3);
	leg2->AddEntry(BPCross,"B^{+}","PL");
	leg2->AddEntry(BsCross,"B^{0}_{s}","PL");
	leg2->Draw("same");


	c->SaveAs("BsBPCross.png");


	TH1D * BsBPRatio =  (TH1D * ) BsCross->Clone("BsBPRatio");

	BsBPRatio->GetXaxis()->SetTitle("B p_{T} (GeV/c)");
	BsBPRatio->GetYaxis()->SetTitle("B_{s}^{0}/B^{+}");
	BsBPRatio->SetTitle("B_{s}^{0}/B^{+} vs B p_{T}");

	BsBPRatio->SetMarkerStyle(20);
	BsBPRatio->SetMarkerSize(1);
	BsBPRatio->SetMarkerColor(kBlack);
	BsBPRatio->SetLineColor(kBlack);

	BsBPRatio->Sumw2();
	BPCross->Sumw2();




	BsBPRatio->Divide(BPCross);


	BsBPRatio->Draw("ep");

	TLine * BsBPFF = new TLine(2,0.255,100,0.255);
	BsBPFF->SetLineStyle(2);
	BsBPFF->SetLineWidth(2);
	BsBPFF->SetLineColor(2);
	BsBPFF->Draw("SAME");


	TLegend* leg = new TLegend(0.30,0.50,0.60,0.70,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(BsBPRatio,"Data Points","PL");
	leg->AddEntry(BsBPFF,"B^{0}_{s}/B^{+} Frag Frac","L");
	leg->Draw("same");


	c->SaveAs("BsBPRatio.png");


	TH1D * BsCross2D = (TH1D *) FileBs->Get("hPtSigma");

	TH1D * BPCross2D = (TH1D *) FileBP->Get("hPtSigma");


	TH1D * BsBPRatio2D =  (TH1D * ) BsCross2D->Clone("BsBPRatio2D");

	BsBPRatio2D->GetXaxis()->SetTitle("B p_{T} (GeV/c)");
	BsBPRatio2D->GetYaxis()->SetTitle("B_{s}^{0}/B^{+}");
	BsBPRatio2D->SetTitle("B_{s}^{0}/B^{+} vs B p_{T} with 2D Map");

	BsBPRatio2D->SetMarkerStyle(20);
	BsBPRatio2D->SetMarkerSize(1);
	BsBPRatio2D->SetMarkerColor(kBlack);
	BsBPRatio2D->SetLineColor(kBlack);

	BsBPRatio2D->Sumw2();
	BPCross2D->Sumw2();

	BsBPRatio2D->Divide(BPCross);


	BsBPRatio2D->Draw("ep");





	c->SaveAs("BsBPRatio2D.png");


	c->SetLogy();

	BPCross->Draw("ep");
	BsCross->Draw("epSAME");
	leg2->Draw("SAME");

	
	
	c->SaveAs("BsBPCrossLog.png");

}
