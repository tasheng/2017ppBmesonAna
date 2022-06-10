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
#include "tnp_weight.h"

using namespace std;

using std::cout;
using std::endl;


void TnPWeight(int Opt){


	const int NCand = 8000;

	TString infile;
	TString TreeName;
	TString outfile;




	if(Opt == 0){


	//	infile = "../UnskimmedSamples/OfficialMC/BPMC.root";
		// infile = "../../bmva/TMVA/BP/sample/BPMC_5_60.root";
		// infile = "../../dat/BP_MC_all.root";
		infile = "../../dat/Unskimmed_gen/BP_MC_merged.root";

		
		TreeName = "Bfinder/ntKp";
		outfile = "BPTnPInfo.root";

	}
	if(Opt == 1){

		// infile = "../UnskimmedSamples/OfficialMC/BsMC.root";
		// infile = "../../bmva/TMVA/Bs/sample/BsMC_7_50.root";
		// infile = "../../dat/Bs_MC_all.root";
		infile = "../../dat/Unskimmed_gen/Bs_MC_merged.root";
		TreeName = "Bfinder/ntphi";
		//	outfile = "BsTnPInfo.root";
		outfile = "BsTnPInfo.root";

	}



	if(Opt == 2){

		infile = "/data/szhaozho/ppNewTMVA/CMSSW_10_3_2/src/ForMariaNew/BZ/BZMC.root";
		TreeName = "ntKstar";
		outfile = "BZTnPInfo.root";
	}






	float Bmu1pt[NCand];
	float Bmu2pt[NCand];

	float Bmu1eta[NCand];
	float Bmu2eta[NCand];

	int Bsize;

	int BsizeNew;


	TFile * fin = new TFile(infile.Data());
	fin->cd();

	cout << "Pass 1" << endl;


	cout << "TreeName = " << TreeName.Data() << endl;

	TTree * t = (TTree * ) fin->Get(TreeName.Data());

	cout << "Pass 1.5" << endl;



	t->SetBranchAddress("Bmu1pt",Bmu1pt);
	t->SetBranchAddress("Bmu2pt",Bmu2pt);
	t->SetBranchAddress("Bmu1eta",Bmu1eta);
	t->SetBranchAddress("Bmu2eta",Bmu2eta);
	t->SetBranchAddress("Bsize",&Bsize);


	cout << "Pass 2" << endl;


	int BsizeTnP;


	float TnPMu1Nominal[NCand];
	float TnPMu1StatError[NCand];
	float TnPMu1SystError[NCand];
	float TnPMu1Error[NCand];

	float TnPMu2Nominal[NCand];
	float TnPMu2StatError[NCand];
	float TnPMu2SystError[NCand];
	float TnPMu2Error[NCand];



	float TnPNominal[NCand];
	float TnPStatError[NCand];
	float TnPSystError[NCand];
	float TnPError[NCand];


	TFile * fout = new TFile(outfile.Data(),"RECREATE");
	fout->cd();

	TTree * TnPInfo = new TTree("TnPInfo","TnPInfo");
	//TnPInfo->Branch("BsizeTnP",BsizeTnP,"BsizeTnP/F");

	TnPInfo->Branch("BsizeNew",&BsizeNew,"BsizeNew/I");

	TnPInfo->Branch("TnPMu1Nominal",TnPMu1Nominal,"TnPMu1Nominal[BsizeNew]/F");
	TnPInfo->Branch("TnPMu1StatError",TnPMu1StatError,"TnPMu1StatError[BsizeNew]/F");
	TnPInfo->Branch("TnPMu1SystError",TnPMu1SystError,"TnPMu1SystError[BsizeNew]/F");
	TnPInfo->Branch("TnPMu1Error",TnPMu1Error,"TnPMu1Error[BsizeNew]/F");

	TnPInfo->Branch("TnPMu2Nominal",TnPMu2Nominal,"TnPMu2Nominal[BsizeNew]/F");
	TnPInfo->Branch("TnPMu2StatError",TnPMu2StatError,"TnPMu2StatError[BsizeNew]/F");
	TnPInfo->Branch("TnPMu2SystError",TnPMu2SystError,"TnPMu2SystError[BsizeNew]/F");
	TnPInfo->Branch("TnPMu2Error",TnPMu2Error,"TnPMu2Error[BsizeNew]/F");


	TnPInfo->Branch("TnPNominal",TnPNominal,"TnPNominal[BsizeNew]/F");
	TnPInfo->Branch("TnPStatError",TnPStatError,"TnPStatError[BsizeNew]/F");
	TnPInfo->Branch("TnPSystError",TnPSystError,"TnPSystError[BsizeNew]/F");
	TnPInfo->Branch("TnPError",TnPError,"TnPError[BsizeNew]/F");



	cout << "Pass 4" << endl;


	int NEvent = t->GetEntries();

	//NEvent = 10;

	for(int i = 0; i < NEvent; i++){

		if(i%100000 == 0) cout << "Event = " << i << endl;

		t->GetEntry(i);
		//MuonInfoTree->GetEntry(i);

		//		BsizeTnP = Bsize;

//		cout  << "i = " << i << "  Bmu1pt =" << Bmu1pt << "  Bmu1eta = " << Bmu1eta << endl;
		BsizeNew = Bsize;




		for(int j = 0; j < Bsize; j++){

		auto TnPSF1 = tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(Bmu1pt[j],Bmu1eta[j]);
		auto TnPSF2 = tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(Bmu2pt[j],Bmu2eta[j]);
		
		//cout << "Bmu1pt[j] = " << Bmu1pt[j] << " Bmu1eta[j] =  " << Bmu1eta[j] << endl;
		//cout << "Bmu2pt[j] = " << Bmu2pt[j] << " Bmu2eta[j] =  " << Bmu2eta[j] << endl;


		TnPMu1Nominal[j] =  std::get<0>(TnPSF1);
		TnPMu1StatError[j] = std::get<1>(TnPSF1);
		TnPMu1SystError[j] = std::get<2>(TnPSF1);
		TnPMu1Error[j] = std::get<3>(TnPSF1);


		TnPMu2Nominal[j] =  std::get<0>(TnPSF2);
		TnPMu2StatError[j] = std::get<1>(TnPSF2);
		TnPMu2SystError[j] = std::get<2>(TnPSF2);
		TnPMu2Error[j] = std::get<3>(TnPSF2);

//		cout << "TnPNominal = " << TnPNominal << "   TnPMu1Nominal = " << TnPMu1Nominal <<  "   TnPMu2Nominal = " << TnPMu2Nominal  << endl;

		TnPNominal[j] =  TnPMu1Nominal[j] * TnPMu2Nominal[j];
		TnPStatError[j] = TnPNominal[j] * TMath::Sqrt(TnPMu1StatError[j]/TnPMu1Nominal[j] * TnPMu1StatError[j]/TnPMu1Nominal[j] + TnPMu2StatError[j]/TnPMu2Nominal[j] * TnPMu2StatError[j]/TnPMu2Nominal[j]);
		TnPSystError[j] = TnPNominal[j] * TMath::Sqrt(TnPMu1SystError[j]/TnPMu1Nominal[j] * TnPMu1SystError[j]/TnPMu1Nominal[j] + TnPMu2SystError[j]/TnPMu2Nominal[j] * TnPMu2SystError[j]/TnPMu2Nominal[j]);
		TnPError[j] = TnPNominal[j] * TMath::Sqrt(TnPMu1Error[j]/TnPMu1Nominal[j] * TnPMu1Error[j]/TnPMu1Nominal[j] + TnPMu2Error[j]/TnPMu2Nominal[j] * TnPMu2Error[j]/TnPMu2Nominal[j]);



		}


		TnPInfo->Fill();	


	}

	fout->Write();
	fout->Close();

	fin->Close();
}


