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
#include "TTree.h"

using namespace std;

using std::cout;
using std::endl;
#endif



void BMuonMatch(int CentMin, int CentMax,int PtOpt){

	TString BType;

	if(PtOpt == -1) BType = "10-50";
	if(PtOpt == 0) BType = "7-10";
	if(PtOpt == 1) BType = "10-15";
	if(PtOpt == 2) BType = "15-20";
	if(PtOpt == 3) BType = "20-50";

	const int NCand = 1000;

	TString infile = Form("MatchedInfo/%d-%d/%s/EffInfo.root",CentMin,CentMax,BType.Data());
	
	TFile * fin = new TFile(infile.Data());
	fin->cd();
	
	TTree * root = (TTree *) fin->Get("root");
	TTree * EffInfoTree = (TTree *) fin->Get("EffInfoTree");

	int evt;
	int lumi;
	int run;



	float Bmu1pt;
	float Bmu1eta;
	float Bmu1phi;

	float Bmu2pt;
	float Bmu2eta;
	float Bmu2phi;

	int MuonSize; 

	float MuonInfoPt[NCand];
	float MuonInfoEta[NCand];
	float MuonInfoPhi[NCand];



	float MuTrgMatchFilterTrgObjPt[NCand];
	float MuTrgMatchFilterTrgObjEta[NCand];
	float MuTrgMatchFilterTrgObjPhi[NCand];

	float MuTrgMatchFilterTrgObjE[NCand];
	Bool_t MuonTrigger[NCand];

	root->SetBranchAddress("MuonInfo.size",&MuonSize);

	root->SetBranchAddress("MuonInfo.pt",MuonInfoPt);
	root->SetBranchAddress("MuonInfo.eta",MuonInfoEta);
	root->SetBranchAddress("MuonInfo.phi",MuonInfoPhi);


	root->SetBranchAddress("MuonInfo.MuTrgMatchFilterTrgObjE",MuTrgMatchFilterTrgObjE);
	root->SetBranchAddress("MuonInfo.MuTrgMatchFilterTrgObjPt",MuTrgMatchFilterTrgObjPt);
	root->SetBranchAddress("MuonInfo.MuTrgMatchFilterTrgObjEta",MuTrgMatchFilterTrgObjEta);
	root->SetBranchAddress("MuonInfo.MuTrgMatchFilterTrgObjPhi",MuTrgMatchFilterTrgObjPhi);
	root->SetBranchAddress("MuonInfo.isTriggered",MuonTrigger);




	


	EffInfoTree->SetBranchAddress("runNew",&run);
	EffInfoTree->SetBranchAddress("evtNew",&evt);
	EffInfoTree->SetBranchAddress("lumiNew",&lumi);



	EffInfoTree->SetBranchAddress("Bmu1ptNew",&Bmu1pt);
	EffInfoTree->SetBranchAddress("Bmu1etaNew",&Bmu1eta);
	EffInfoTree->SetBranchAddress("Bmu1phiNew",&Bmu1phi);


	EffInfoTree->SetBranchAddress("Bmu2ptNew",&Bmu2pt);
	EffInfoTree->SetBranchAddress("Bmu2etaNew",&Bmu2eta);
	EffInfoTree->SetBranchAddress("Bmu2phiNew",&Bmu2phi);


	int NEvents = EffInfoTree->GetEntries();


	TString OutputFileName = Form("FinalInput/%d-%d/%s/MuonInfo.root",CentMin,CentMax,BType.Data());

	TFile * fout = new TFile(OutputFileName.Data(),"RECREATE");	
	TTree * MuonInfoTree = new TTree("MuonInfoTree","MuonInfoTree");


	int evtMatched;
	int lumiMatched;
	int runMatched;


	double PtDiff = 0.3;
	double EtaDiff = 0.3;
	double PhiDiff = 0.2;



	float Bmu1ptMatched;
	float Bmu1etaMatched;
	float Bmu1phiMatched;

	float Bmu2ptMatched;
	float Bmu2etaMatched;
	float Bmu2phiMatched;


	int Bmu1Type;
	int Bmu2Type;

	int Muon1Matched = 0;
	int Muon2Matched = 0;


	MuonInfoTree->Branch("evt",&evtMatched,"evt/I");
	MuonInfoTree->Branch("run",&runMatched,"run/I");
	MuonInfoTree->Branch("lumi",&lumiMatched,"lumi/I");


	MuonInfoTree->Branch("Bmu1ptMatched",&Bmu1ptMatched,"Bmu1ptMatched/F");
	MuonInfoTree->Branch("Bmu1etaMatched",&Bmu1etaMatched,"Bmu1etaMatched/F");
	MuonInfoTree->Branch("Bmu1phiMatched",&Bmu1phiMatched,"Bmu1phiMatched/F");

	MuonInfoTree->Branch("Bmu2ptMatched",&Bmu2ptMatched,"Bmu2ptMatched/F");
	MuonInfoTree->Branch("Bmu2etaMatched",&Bmu2etaMatched,"Bmu2etaMatched/F");
	MuonInfoTree->Branch("Bmu2phiMatched",&Bmu2phiMatched,"Bmu2phiMatched/F");

	MuonInfoTree->Branch("Bmu1Type",&Bmu1Type,"Bmu1Type/I");
	MuonInfoTree->Branch("Bmu2Type",&Bmu2Type,"Bmu2Type/I");


	



	for(int i = 0; i < NEvents; i++){
		
		root->GetEntry(i);
		EffInfoTree->GetEntry(i);
		
		//Muon1Matched = 0;
		//Muon2Matched = 0;
	
		evtMatched = evt;
		lumiMatched = lumi;
		runMatched = run;

	

		for(int j = 0; j < MuonSize; j++){
		
		//	if( abs(Bmu1pt - MuTrgMatchFilterTrgObjPt[j * 3]) < PtDiff && abs(Bmu1eta - MuTrgMatchFilterTrgObjEta[j * 3]) < EtaDiff && (Bmu1phi - MuTrgMatchFilterTrgObjPhi[j* 3]) < PhiDiff ){
		
		

			if( Bmu1pt == MuonInfoPt[j] && Bmu1eta == MuonInfoEta[j] && Bmu1phi == MuonInfoPhi[j] ){

			
				Bmu1ptMatched = Bmu1pt;
				Bmu1etaMatched = Bmu1eta;
				Bmu1phiMatched = Bmu1phi;


				if(MuTrgMatchFilterTrgObjE[j * 3] < 0){
					Bmu1Type = -1;
					cout << "Muon 1 - Event Fucked Up : " << i << endl;
				}
				if(MuTrgMatchFilterTrgObjE[j * 3] > 0 && MuTrgMatchFilterTrgObjE[j * 3 + 1] > 0  && MuTrgMatchFilterTrgObjE[j * 3 + 2] < 0) Bmu1Type = 0;
				if(MuTrgMatchFilterTrgObjE[j * 3] > 0 && MuTrgMatchFilterTrgObjE[j * 3 + 1] > 0  && MuTrgMatchFilterTrgObjE[j * 3 + 2] > 0) Bmu1Type = 1;
				
							
				
			}
			
	

	//		if(abs(Bmu2pt - MuTrgMatchFilterTrgObjPt[j * 3]) < PtDiff && abs(Bmu2eta - MuTrgMatchFilterTrgObjEta[j * 3]) < EtaDiff && (Bmu2phi - MuTrgMatchFilterTrgObjPhi[j* 3]) < PhiDiff){
			


			if( Bmu2pt == MuonInfoPt[j] && Bmu2eta == MuonInfoEta[j] && Bmu2phi == MuonInfoPhi[j] ){


				Bmu2ptMatched = Bmu2pt;
				Bmu2etaMatched = Bmu2eta;
				Bmu2phiMatched = Bmu2phi;

				if(MuTrgMatchFilterTrgObjE[j * 3] < 0){
					Bmu2Type = -1;
					cout << "Muon 2 - Event Fucked Up : " << i << endl;
				}

				if(MuTrgMatchFilterTrgObjE[j * 3] > 0 && MuTrgMatchFilterTrgObjE[j * 3 + 1] > 0  && MuTrgMatchFilterTrgObjE[j * 3 + 2] < 0) Bmu2Type = 0;
				if(MuTrgMatchFilterTrgObjE[j * 3] > 0 && MuTrgMatchFilterTrgObjE[j * 3 + 1] > 0  && MuTrgMatchFilterTrgObjE[j * 3 + 2] > 0) Bmu2Type = 1;
			
			}


		}
		
		MuonInfoTree->Fill();

	}


	cout << "DONE BRO" << endl;
	fout->Write();
	fout->Close();




}
