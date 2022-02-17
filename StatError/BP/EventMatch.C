#ifndef __CINT__
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
#endif


void EventMatch(int CentMin, int CentMax, int PtOpt){



	TString BType;

	if(PtOpt == -1) BType = "10-50";
	if(PtOpt == 0) BType = "7-10";
	if(PtOpt == 1) BType = "10-15";
	if(PtOpt == 2) BType = "15-20";
	if(PtOpt == 3) BType = "20-50";


	int EventMuon;
	int RunMuon;
	int LumiMuon;

	//	unsigned long int MuonKey;

	int EventData;
	int RunData;
	int LumiData;


	//int HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1;

	//	unsigned long int DataKey;

	TFile * MuonInfoFile = new TFile("/export/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/BsBP2018PbPb/FilePrep/InputFiles/BMuonInfo.root");



	TString infile = Form("CheckSystNuno/%d-%d/%s/EffInfo.root",CentMin,CentMax,BType.Data());

	TFile * fin = new TFile(infile.Data());


	TTree * root = (TTree *) MuonInfoFile->Get("Bfinder/root");
	root->SetBranchAddress("EvtInfo.RunNo",&RunMuon);
	root->SetBranchAddress("EvtInfo.EvtNo",&EventMuon);
	root->SetBranchAddress("EvtInfo.LumiNo",&LumiMuon);

	TTree * HltTree = (TTree *) MuonInfoFile->Get("hltanalysis/HltTree");
//	root->SetBranchAddress("HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1",&HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1);





	TTree * EffInfoTree = (TTree *) fin->Get("EffInfoTree");
	EffInfoTree->SetBranchAddress("runNew",&RunData);
	EffInfoTree->SetBranchAddress("evtNew",&EventData);
	EffInfoTree->SetBranchAddress("lumiNew",&LumiData);


	int NMuonEvent = root->GetEntries();
	const int NDataEvent = EffInfoTree->GetEntries();
//	const int NDataEvent = 1;

	int EventtoEntry[NDataEvent]; 

	for(int i = 0; i < NDataEvent; i++){

		EventtoEntry[i] = 0;


	}



	for(int i = 0; i < NDataEvent; i++){

		EffInfoTree->GetEntry(i);

		

		for(int j  = 0; j < NMuonEvent; j++){
			root->GetEntry(j);
	
			/*
			if(j%1000000==0){
				cout << "Now checking on entry " << j << endl;  
				
				cout << "RunMuon = " << RunMuon << "  RunData =  " << RunData <<  endl;
				cout << "EventMuon = " << EventMuon << "  EventData =  " << EventData << endl;
				cout << "LumiMuon = " << LumiMuon << "  LumiData =  " << LumiData << endl;

			}

			*/

			if(RunMuon == RunData && EventMuon == EventData && LumiMuon == LumiData ){

				cout << "FUCK BRO: Event " << i << "  Matched to Entry " << j << endl;
				EventtoEntry[i] = j;

			}

		}

	}

	TString outfile = Form("MatchedInfo/%d-%d/%s/EffInfo.root",CentMin,CentMax,BType.Data());

	TFile * fout = new TFile(outfile.Data(),"RECREATE");
	fout->cd();
	TTree * EffInfoTree_new = EffInfoTree->CloneTree(0);
	TTree * root_new = root->CloneTree(0);
	TTree * HltTree_new = HltTree->CloneTree(0);

	for(int i = 0; i < NDataEvent; i++){

		cout << "Event " << i << "   Corresponds to Entry  " << EventtoEntry[i] << endl;

		EffInfoTree->GetEntry(i);
		root->GetEntry( EventtoEntry[i]);
		HltTree->GetEntry( EventtoEntry[i]);


		EffInfoTree_new->Fill();
		root_new->Fill();
		HltTree_new->Fill();

	}
	fout->Write();

	fin->Close();

	MuonInfoFile->Close();

	fout->Close();



}
