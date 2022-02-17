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
#include <vector>
#include <random>

//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void SpiltDataBoot(int PtOpt){


	int Factor = 1;

	double NBackground = 40;
	double Lambda = -1.0;

	TString Func = Form("%f * TMath::Exp(%f * x)",NBackground,Lambda);

	TF1 * bkgd = new TF1("bkgd",Func.Data(),5,6);






	TString infile;



	TString outfilefolder;

	cout << "Pass 0 " << endl;

	const int NBins = 7;
	int PtBin[NBins + 1] = {5,7,10,15,20,30,50,60};
	int Start[NBins] = {0,0,0,0,0,0,0};




	infile = Form("InputFiles/BPData_%d.root",PtOpt);

	TFile * fin = new TFile(infile.Data());
	fin->cd();

	cout << "Pass 0.5 " << endl;

	TTree * ntKp = (TTree * ) fin->Get("ntKp");
	//TTree * MuonInfoTree = (TTree * ) fin->Get("MuonInfoTree");

	int NEvents = ntKp->GetEntries();
	
//	cout << "Pass 0.6 " << endl;


	int BsizeNew;




	ntKp->SetBranchAddress("Bsize",&BsizeNew);


	



	cout << "Pass 1 " << endl;




	TString OutFileName;

	OutFileName = Form("DataResample/%d-%d/Data_0.root",PtBin[PtOpt],PtBin[PtOpt+1]);

	
	TFile * outf = new TFile(OutFileName.Data(),"RECREATE");;


	outf->cd();

	TTree* ntKp_New;
	TTree* ntKp_NewFit;
//	TTree* MuonInfoTree_New;

	int NSets = 1000;

	int Index = 0;
	int IndexPre = 0;

	cout << "Pass 2 " << endl;
	NEvents = ntKp->GetEntries();

	cout << "NEvent = " << NEvents << endl;
	ntKp_New = ntKp->CloneTree(0);
	ntKp_NewFit = ntKp->CloneTree(0);
	ntKp_NewFit->SetObject("ntKpFit","ntKpFit");
//	MuonInfoTree_New = MuonInfoTree->CloneTree(0);

//	cout << "NEvents = " << NEvents << endl;

	int SampleSize;
	std::vector<int> EventVec;


	for(int i = 0; i < NEvents; i++){

		EventVec.push_back(i);

	}
	

	TString f1Name = Form("TMath::Poisson(x,%d)",NEvents);
	TF1 *f1 = new TF1("f1",f1Name.Data(),0,1000000);

	TF1 *f2 = new TF1("f2","x",-0.5,NEvents-0.5);



	//		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

	

	for(int i = Start[PtOpt]; i < NSets; i++){

		cout << "Now Working on " << i << "  File" << endl;

		//std::shuffle(EventVec.begin(), EventVec.end(), std::default_random_engine(seed));

	

		SampleSize =  f1->GetRandom();

//		cout << "SampleSize = " << SampleSize << endl;

		for(int k = 0; k< SampleSize;k++){
			EventVec.push_back(f2->GetRandom());
		}


		outf = new TFile(Form("DataResample/%d-%d/Data_%d.root",PtBin[PtOpt],PtBin[PtOpt+1],i),"RECREATE");
		outf->cd();
		ntKp_New = ntKp->CloneTree(0);
		ntKp_NewFit = ntKp->CloneTree(0);
		ntKp_NewFit->SetObject("ntKpFit","ntKpFit");
//		MuonInfoTree_New = MuonInfoTree->CloneTree(0);


		for(int j = 0; j < SampleSize; j++){

			ntKp->GetEntry(EventVec[j]);
//			MuonInfoTree->GetEntry(EventVec[j]);


			ntKp_New->Fill();
//			MuonInfoTree_New->Fill();
			if(BsizeNew > 0) ntKp_NewFit->Fill();

		}


	

		outf->Write();
		outf->Close();
		EventVec.clear();
	}






}
