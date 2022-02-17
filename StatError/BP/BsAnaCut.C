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




//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void BsAnaCut(int Opt){


	const int NBins = 4;
	float PtBins[NBins + 1] = {7,10,15,20,50};


	TString FileName;
	FileName="/data3/tasheng/CutSkim/BsData.root"; 


	TFile * fin = new TFile(FileName.Data());
	fin->cd();

	TTree * ntphi = (TTree * ) fin->Get("ntphi");

	float Bpt;
	ntphi->SetBranchAddress("Bpt",&Bpt);
	



	TString OutputFileName = Form("InputFiles/BsData_%d.root",Opt);

	TFile *fout = new TFile(OutputFileName.Data(),"RECREATE");


	TTree * ntphi_new = ntphi->CloneTree(0);
	

	int NEvents = ntphi->GetEntries();

	for(int i = 0; i < NEvents; i++){

		ntphi->GetEntry(i);
	
	
		if(Bpt > PtBins[Opt] && Bpt <  PtBins[Opt+1])	ntphi_new->Fill();

	}


	ntphi_new->Write();
	fout->Close();

}




