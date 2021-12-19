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

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
//#include "tnp_weight_lowptPBsb.h"



//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void QuickComp(){


	gStyle->SetOptStat(0);
	TCanvas * c = new TCanvas("c","c",600,600);	
	c->cd();

	TFile * fin = new TFile("NewEff2DMaps/EffFineBDT.root");

	TFile * fin2 = new TFile("NewEff2DMaps/EffFineBDTNew.root");

	TH1D * Eff1DHis = (TH1D * ) fin->Get("Eff1DHis");
	TH1D * Eff1DHisNew = (TH1D * ) fin2->Get("Eff1DHis");

	Eff1DHis->SetMarkerSize(1);
	Eff1DHis->SetMarkerStyle(20);
	Eff1DHis->SetMarkerColor(kBlue);
	Eff1DHis->SetLineColor(kBlue);

	Eff1DHisNew->SetMarkerSize(1);
	Eff1DHisNew->SetMarkerStyle(20);
	Eff1DHisNew->SetMarkerColor(kGreen);
	Eff1DHisNew->SetLineColor(kGreen);	


	Eff1DHis->Draw("ep");
	Eff1DHisNew->Draw("epSAME");

	
	TLegend* leg = new TLegend(0.42,0.30,0.75,0.50,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(Eff1DHis,"No Rescaling of BsvpvErr","PL");
	leg->AddEntry(Eff1DHisNew,"Rescale of BsvpvErr","PL");	

	leg->Draw("same");

	c->SaveAs("EffComp.png");


}



