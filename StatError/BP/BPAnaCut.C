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

void BPAnaCut(int Opt){


	const int NBins = 7;
	float PtBins[NBins + 1] = {5,7,10,15,20,30,50,60};


	TString FileName;
	FileName="/data3/tasheng/CutSkim/BPData.root"; 


	TFile * fin = new TFile(FileName.Data());
	fin->cd();

	TTree * ntKp = (TTree * ) fin->Get("ntKp");

	float Bpt;
	ntKp->SetBranchAddress("Bpt",&Bpt);
	



	TString OutputFileName = Form("InputFiles/BPData_%d.root",Opt);

	TFile *fout = new TFile(OutputFileName.Data(),"RECREATE");


	TTree * ntKp_new = ntKp->CloneTree(0);
	

	int NEvents = ntKp->GetEntries();

	for(int i = 0; i < NEvents; i++){

		ntKp->GetEntry(i);
	
	
		if(Bpt > PtBins[Opt] && Bpt <  PtBins[Opt+1])	ntKp_new->Fill();

	}


	ntKp_new->Write();
	fout->Close();

}




