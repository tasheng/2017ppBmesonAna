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


void ModifyEff(){

	float MCGenMatched = 1;
	float MCGenMatchedErr = 1;
	
	float GenCounts;
	
	float EffModified;

	float EffModifiedErr;


	TString infile = "NewEff2DMaps/EffFineBDT.root";
	
	TFile * fin = new TFile(infile.Data());
	fin->cd();

	TH1D * Eff1DHis =  (TH1D * ) fin->Get("Eff1DHis");

	TH1D * Eff1DGENHis =  (TH1D * ) fin->Get("Eff1DGENHis");


	GenCounts = Eff1DGENHis->GetBinContent(1);


	EffModified = MCGenMatched/GenCounts;
	EffModifiedErr = MCGenMatchedErr/GenCounts;

	EffModified  = Eff1DHis->GetBinContent(2)/3;
	EffModifiedErr  = Eff1DHis->GetBinError(2)/3;

	Eff1DHis->SetBinContent(1,EffModified);
	Eff1DHis->SetBinError(1,EffModifiedErr);

	TFile * fout = new TFile("NewEff2DMaps/EffFineBDT2.root","RECREATE");
	fout->cd();
	Eff1DHis->Write();	
	fout->Close();

}
