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



void  DrawEff(){

	gStyle->SetOptStat(0);
	
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();


	TFile *fin = new TFile("NewEff2DMaps/EffFineBDT.root");
	fin->cd();

	TH1D * Eff1DHis = (TH1D *) fin->Get("Eff1DHis");
	Eff1DHis->SetMarkerStyle(20);
	Eff1DHis->SetMarkerSize(1);
	Eff1DHis->SetMarkerColor(1);
	Eff1DHis->SetLineColor(1);
	
	Eff1DHis->Draw("ep");

	c->SaveAs("DrawEff/BsEff1DHis.png");

	TH1D * Eff1DHisMult = (TH1D *) fin->Get("Eff1DHisMult");
	Eff1DHisMult->SetMarkerStyle(20);
	Eff1DHisMult->SetMarkerSize(1);
	Eff1DHisMult->SetMarkerColor(1);
	Eff1DHisMult->SetLineColor(1);
	
	Eff1DHisMult->Draw("ep");

	c->SaveAs("DrawEff/BsEff1DHisMult.png");

}
