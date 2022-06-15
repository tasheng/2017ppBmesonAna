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
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>


using namespace std;

using std::cout;
using std::endl;



void  DrawEff(bool useTnp=true){

	gStyle->SetOptStat(0);
	
	TCanvas * c = new TCanvas("c","c",600,600);
	c->cd();


  TString inMap = "NewEff2DMaps/EffFineBDT.root";
  TString suffix = "";
  if (!useTnp) {
    inMap = "NewEff2DMaps/EffFineNoTnP.root";
    suffix = "_noTnP";
  }
	TFile *fin = new TFile(inMap);
	fin->cd();

	TH1D * Eff1DHis = (TH1D *) fin->Get("Eff1DHis");
	Eff1DHis->SetMarkerStyle(20);
	Eff1DHis->SetMarkerSize(1);
	Eff1DHis->SetMarkerColor(1);
	Eff1DHis->SetLineColor(1);
	
	Eff1DHis->Draw("ep");

  TString eff1dhis = TString::Format("DrawEff/BPEff1DHis%s.png", suffix.Data());
	c->SaveAs(eff1dhis);

	TH1D * Eff1DHisMult = (TH1D *) fin->Get("Eff1DHisMult");
	Eff1DHisMult->SetMarkerStyle(20);
	Eff1DHisMult->SetMarkerSize(1);
	Eff1DHisMult->SetMarkerColor(1);
	Eff1DHisMult->SetLineColor(1);
	
	Eff1DHisMult->Draw("ep");

  TString eff1dhismult = TString::Format("DrawEff/BPEff1DHisMult%s.png", suffix.Data());
	c->SaveAs(eff1dhismult);

  TH2D * invEff2DTnPSyst = (TH2D * ) fin->Get("invEff2D");
  invEff2DTnPSyst->GetYaxis()->SetTitle("B |y|");
  invEff2DTnPSyst->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
  invEff2DTnPSyst->GetXaxis()->CenterTitle();
  invEff2DTnPSyst->GetYaxis()->CenterTitle();
  invEff2DTnPSyst->GetYaxis()->SetTitleOffset(1.2);
  invEff2DTnPSyst->SetTitle("");

  invEff2DTnPSyst->Draw("COLZ");
  TString eff2dtnp = TString::Format("Eff2DMapTnP/Eff2D%s.png", suffix.Data());
  c->SaveAs(eff2dtnp);

  if (useTnp) {
    TH2D * invEff2DTnPSystUp = (TH2D * ) fin->Get("invEff2DTnPSystUp");
    invEff2DTnPSystUp->GetYaxis()->SetTitle("B |y|");
    invEff2DTnPSystUp->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
    invEff2DTnPSystUp->GetXaxis()->CenterTitle();
    invEff2DTnPSystUp->GetYaxis()->CenterTitle();
    invEff2DTnPSystUp->GetYaxis()->SetTitleOffset(1.2);
    invEff2DTnPSystUp->SetTitle("");

    invEff2DTnPSystUp->Draw("COLZ");
    TString eff2dtnpUp = TString::Format("Eff2DMapTnP/Eff2D_Up%s.png", suffix.Data());
    c->SaveAs(eff2dtnpUp);

    TH2D * invEff2DTnPSystDown = (TH2D * ) fin->Get("invEff2DTnPSystDown");
    invEff2DTnPSystDown->GetYaxis()->SetTitle("B |y|");
    invEff2DTnPSystDown->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
    invEff2DTnPSystDown->GetXaxis()->CenterTitle();
    invEff2DTnPSystDown->GetYaxis()->CenterTitle();
    invEff2DTnPSystDown->GetYaxis()->SetTitleOffset(1.2);
    invEff2DTnPSystDown->SetTitle("");

    invEff2DTnPSystDown->Draw("COLZ");
    TString eff2dtnpDown = TString::Format("Eff2DMapTnP/Eff2D_Down%s.png", suffix.Data());
    c->SaveAs(eff2dtnpDown);
  }

}
