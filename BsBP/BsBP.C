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


	
	TFile * FileBP = new TFile(InfileBP.Data());
	TH1D * BPCross = (TH1D *) FileBP->Get("CorrDiffHisBin");
	

	TH1D * BsBPRatio =  (TH1D * ) BsCross->Clone("BsBPRatio");

	BsBPRatio->GetXaxis()->SetTitle("B p_{T} (GeV/c)");
	BsBPRatio->GetYaxis()->SetTitle("B_s^{0}/B^{+}");
	BsBPRatio->SetMarkerStyle(20);
	BsBPRatio->SetMarkerSize(1);
	BsBPRatio->SetMarkerColor(kBlack);
	BsBPRatio->SetLineColor(kBlack);

	BsBPRatio->Sumw2();
	BPCross->Sumw2();

	BsBPRatio->Divide(BPCross);
	

	BsBPRatio->Draw("ep");

	c->SaveAs("BsBPRatio.png");


}
