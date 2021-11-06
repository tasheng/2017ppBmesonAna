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



void  MCEff(){

	gStyle->SetOptStat(0);

	int ptmin = 10;
	int ptmax = 50;

	//TString infile = "/data/szhaozho/2017ppSamplesNew/BDTOutput/AllMerge/BPMCAllBDT.root";
//	TString infile = "../../UnskimmedSamples/BsMC.root";
	TString infile = "../../UnskimmedSamples/BsMC.root";
	TFile * fin = new TFile(infile.Data());

	fin->cd();

	TTree * ntphi = (TTree * ) fin->Get("Bfinder/ntphi");
	//	TTree * BDT = (TTree * ) fin->Get("BDT");
	TTree * ntHi = (TTree * ) fin->Get("hiEvtAnalyzer/HiTree");
	TTree * ntSkim = (TTree * ) fin->Get("skimanalysis/HltTree");
	TTree * ntHlt = (TTree *) fin->Get("hltanalysis/HltTree");
	//	TTree * TnPInfo = (TTree * ) fin->Get("TnPInfo");
	//	TTree * CentWeightTree =	(TTree * ) fin->Get("CentWeightTree");
	TTree * ntGen = (TTree * ) fin->Get("Bfinder/ntGen");
	TTree * BDT1 = (TTree *) fin->Get("BDT_pt_3_5");
	TTree * BDT2 = (TTree *) fin->Get("BDT_pt_5_7");
	TTree * BDT3 = (TTree *) fin->Get("BDT_pt_7_10");
	TTree * BDT4 = (TTree *) fin->Get("BDT_pt_10_15");
	TTree * BDT5 = (TTree *) fin->Get("BDT_pt_15_20");
	TTree * BDT6 = (TTree *) fin->Get("BDT_pt_20_30");
	TTree * BDT9 = (TTree *) fin->Get("BDT_pt_30_50");
	
//	TTree * BDT7 = (TTree *) fin->Get("BDT_pt_2_3");
	TTree * BDT7 = (TTree *) fin->Get("BDT_pt_2_3");
	TTree * BDT8 = (TTree *) fin->Get("BDT_pt_1_2");
	TTree * root = (TTree * ) fin->Get("Bfinder/root");

	Int_t nMult;

	root->SetBranchAddress("EvtInfo.nMult",&nMult);

	const int NCand = 8000;


	int run;
	int lumi;
	int evt;
	int hiBin;
	Float_t PVz;
	Int_t pclusterCompatibilityFilter;
	Int_t pprimaryVertexFilter;
	Int_t phfCoincFilter2Th4;





	Int_t   Bsize;
	Float_t Btrk1Pt[NCand];
	Float_t Btrk2Pt[NCand];

	Float_t Btrk1PtErr[NCand];
	Float_t Btrk2PtErr[NCand];


	Float_t Bchi2cl[NCand];
	Float_t BsvpvDistance[NCand];
	Float_t BsvpvDisErr[NCand];
	Float_t Bpt[NCand];
	Float_t Btrk1Eta[NCand];
	Float_t Btrk2Eta[NCand];
	Float_t By[NCand];

	Bool_t Bmu1isTriggered[NCand];
	Bool_t Bmu2isTriggered[NCand];

	Float_t Bmass[NCand];


	Float_t Bmumumass[NCand];
	Float_t Bmu1eta[NCand];
	Float_t Bmu1pt[NCand];
	Float_t Bmu2eta[NCand];
	Float_t Bmu2pt[NCand];

	//	Float_t Bmu1phi[NCand];
	//	Float_t Bmu2phi[NCand];

	Bool_t Bmu1TMOneStationTight[NCand];
	Int_t Bmu1InPixelLayer[NCand];
	Int_t Bmu1InStripLayer[NCand];

	Bool_t Bmu2TMOneStationTight[NCand];	
	Int_t Bmu2InPixelLayer[NCand];
	Int_t Bmu2InStripLayer[NCand];


	Bool_t Bmu1isGlobalMuon[NCand];
	Bool_t Bmu2isGlobalMuon[NCand];


	Bool_t Bmu1isTrackerMuon[NCand];
	Bool_t Bmu2isTrackerMuon[NCand];

	Float_t Bmu1dxyPV[NCand];
	Float_t Bmu2dxyPV[NCand];

	Float_t Bmu1dzPV[NCand];
	Float_t Bmu2dzPV[NCand];

	Bool_t Btrk1highPurity[NCand];
	Bool_t Btrk2highPurity[NCand];

	Float_t Btktkmass[NCand];

	Float_t Btrk1PixelHit[NCand];
	Float_t Btrk2PixelHit[NCand];

	Float_t Btrk1StripHit[NCand];
	Float_t Btrk2StripHit[NCand];

	Float_t Btrk1Chi2ndf[NCand];
	Float_t Btrk2Chi2ndf[NCand];


	Float_t Btrk1nStripLayer[NCand];
	Float_t Btrk2nStripLayer[NCand];

	Float_t Btrk1nPixelLayer[NCand];
	Float_t Btrk2nPixelLayer[NCand];


	Float_t Bgen[NCand];



	//	Float_t pthatweight;

	Float_t pthat;
	Float_t weight;

	Float_t Bdtheta[NCand];


	Double_t BDT_pt_3_5[NCand];
	Double_t BDT_pt_5_7[NCand];
	Double_t BDT_pt_7_10[NCand];
	Double_t BDT_pt_10_15[NCand];
	Double_t BDT_pt_15_20[NCand];
	Double_t BDT_pt_20_30[NCand];
//	Double_t BDT_pt_2_3[NCand];
	Double_t BDT_pt_2_3[NCand];
	Double_t BDT_pt_1_2[NCand];
	Double_t BDT_pt_30_50[NCand];
	
	/*

	   BDT->SetBranchAddress("BDT_5_7",BDT_pt_5_7);
	   BDT->SetBranchAddress("BDT_7_10",BDT_pt_7_10);
	   BDT->SetBranchAddress("BDT_10_15",BDT_pt_10_15);
	   BDT->SetBranchAddress("BDT_15_20",BDT_pt_15_20);
	   BDT->SetBranchAddress("BDT_20_30",BDT_pt_20_30);
	   BDT->SetBranchAddress("BDT_30_40",BDT_pt_30_40);
	   BDT->SetBranchAddress("BDT_40_50",BDT_pt_40_50);
	   BDT->SetBranchAddress("BDT_50_60",BDT_pt_50_60);

	   BDT->SetBranchAddress("run",&run);
	   BDT->SetBranchAddress("evt",&evt);
	   BDT->SetBranchAddress("lumi",&lumi);

*/

	ntHi->SetBranchAddress("hiBin",&hiBin);
	ntHi->SetBranchAddress("pthat",&pthat);
	ntHi->SetBranchAddress("weight",&weight);


	int HBHENoiseFilterResult;
	int pPAprimaryVertexFilter;
	int pBeamScrapingFilter;

	ntSkim->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
	ntSkim->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
	ntSkim->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);


	ntphi->SetBranchAddress("Bsize",&Bsize);
	ntphi->SetBranchAddress("PVz",&PVz);
	ntphi->SetBranchAddress("Btrk1Pt",Btrk1Pt);
	ntphi->SetBranchAddress("Btrk2Pt",Btrk2Pt);
	ntphi->SetBranchAddress("Btrk1PtErr",Btrk1PtErr);
	ntphi->SetBranchAddress("Btrk2PtErr",Btrk2PtErr);
	ntphi->SetBranchAddress("Bchi2cl",Bchi2cl);
	ntphi->SetBranchAddress("BsvpvDistance",BsvpvDistance);
	ntphi->SetBranchAddress("BsvpvDisErr",BsvpvDisErr);
	ntphi->SetBranchAddress("Bpt",Bpt);
	ntphi->SetBranchAddress("By",By);
	ntphi->SetBranchAddress("Btrk1Eta",Btrk1Eta);
	ntphi->SetBranchAddress("Btrk2Eta",Btrk2Eta);
	ntphi->SetBranchAddress("Bmass",Bmass);
	ntphi->SetBranchAddress("Bdtheta",Bdtheta);



	ntphi->SetBranchAddress("Bmu1isTriggered",Bmu1isTriggered);
	ntphi->SetBranchAddress("Bmu2isTriggered",Bmu2isTriggered);


	ntphi->SetBranchAddress("Bmumumass",Bmumumass);
	ntphi->SetBranchAddress("Bmu1eta",Bmu1eta);
	ntphi->SetBranchAddress("Bmu2eta",Bmu2eta);
	ntphi->SetBranchAddress("Bmu1pt",Bmu1pt);
	ntphi->SetBranchAddress("Bmu2pt",Bmu2pt);

	//	ntphi->SetBranchAddress("Bmu1phi",Bmu1phi);
	//	ntphi->SetBranchAddress("Bmu2phi",Bmu2phi);

	ntphi->SetBranchAddress("Bmu1TMOneStationTight",Bmu1TMOneStationTight);
	ntphi->SetBranchAddress("Bmu1InPixelLayer",Bmu1InPixelLayer);
	ntphi->SetBranchAddress("Bmu1InStripLayer",Bmu1InStripLayer);

	ntphi->SetBranchAddress("Bmu2TMOneStationTight",Bmu2TMOneStationTight);
	ntphi->SetBranchAddress("Bmu2InPixelLayer",Bmu2InPixelLayer);
	ntphi->SetBranchAddress("Bmu2InStripLayer",Bmu2InStripLayer);


	ntphi->SetBranchAddress("Bmu1isGlobalMuon",Bmu1isGlobalMuon);
	ntphi->SetBranchAddress("Bmu2isGlobalMuon",Bmu2isGlobalMuon);

	ntphi->SetBranchAddress("Bmu1isTrackerMuon",Bmu1isTrackerMuon);
	ntphi->SetBranchAddress("Bmu2isTrackerMuon",Bmu2isTrackerMuon);


	//	Int_t HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1;


	//	ntHlt->SetBranchAddress("HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1",&HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1);



	ntphi->SetBranchAddress("Bmu1dxyPV",Bmu1dxyPV);
	ntphi->SetBranchAddress("Bmu2dxyPV",Bmu2dxyPV);
	ntphi->SetBranchAddress("Bmu1dzPV",Bmu1dzPV);
	ntphi->SetBranchAddress("Bmu2dzPV",Bmu2dzPV);


	ntphi->SetBranchAddress("Btrk1highPurity",Btrk1highPurity);
	ntphi->SetBranchAddress("Btrk2highPurity",Btrk2highPurity);

	ntphi->SetBranchAddress("Btktkmass",Btktkmass);


	ntphi->SetBranchAddress("Btrk1PixelHit",Btrk1PixelHit);
	ntphi->SetBranchAddress("Btrk2PixelHit",Btrk2PixelHit);
	ntphi->SetBranchAddress("Btrk1StripHit",Btrk1StripHit);
	ntphi->SetBranchAddress("Btrk2StripHit",Btrk2StripHit);



	ntphi->SetBranchAddress("Btrk1Chi2ndf",Btrk1Chi2ndf);
	ntphi->SetBranchAddress("Btrk2Chi2ndf",Btrk2Chi2ndf);

	ntphi->SetBranchAddress("Bgen",Bgen);


	ntphi->SetBranchAddress("Btrk1nStripLayer",Btrk1nStripLayer);	
	ntphi->SetBranchAddress("Btrk2nStripLayer",Btrk2nStripLayer);
	ntphi->SetBranchAddress("Btrk1nPixelLayer",Btrk1nPixelLayer);
	ntphi->SetBranchAddress("Btrk2nPixelLayer",Btrk2nPixelLayer);


	BDT1->SetBranchAddress("BDT_pt_3_5",BDT_pt_3_5);
	BDT2->SetBranchAddress("BDT_pt_5_7",BDT_pt_5_7);
	BDT3->SetBranchAddress("BDT_pt_7_10",BDT_pt_7_10);
	BDT4->SetBranchAddress("BDT_pt_10_15",BDT_pt_10_15);
	BDT5->SetBranchAddress("BDT_pt_15_20",BDT_pt_15_20);
	BDT6->SetBranchAddress("BDT_pt_20_30",BDT_pt_20_30);
	BDT9->SetBranchAddress("BDT_pt_30_50",BDT_pt_30_50);

	BDT7->SetBranchAddress("BDT_pt_2_3",BDT_pt_2_3);
	BDT8->SetBranchAddress("BDT_pt_1_2",BDT_pt_1_2);

	Bool_t Bmu1SoftMuID[NCand];
	Bool_t Bmu2SoftMuID[NCand];
	Bool_t Bmu1isAcc[NCand];
	Bool_t Bmu2isAcc[NCand];

	ntphi->SetBranchAddress("Bmu1SoftMuID",Bmu1SoftMuID);
	ntphi->SetBranchAddress("Bmu2SoftMuID",Bmu2SoftMuID);


	ntphi->SetBranchAddress("Bmu1isAcc",Bmu1isAcc);
	ntphi->SetBranchAddress("Bmu2isAcc",Bmu2isAcc);



	Int_t Gsize;
	Float_t Gy[NCand];
	Float_t Gpt[NCand];
	Int_t GisSignal[NCand];
	Int_t GcollisionId[NCand];
	Int_t GpdgId[NCand];

	Float_t Gmu1pt[NCand];
	Float_t Gmu1eta[NCand];
	Float_t Gmu1phi[NCand];

	Float_t Gmu2pt[NCand];
	Float_t Gmu2eta[NCand];
	Float_t Gmu2phi[NCand];

	Float_t Gtk1pt[NCand];
	Float_t Gtk1eta[NCand];
	Float_t Gtk1phi[NCand];

	Float_t Gtk2pt[NCand];
	Float_t Gtk2eta[NCand];
	Float_t Gtk2phi[NCand];

	ntGen->SetBranchAddress("Gsize",&Gsize);
	ntGen->SetBranchAddress("Gy",Gy);
	ntGen->SetBranchAddress("Gpt",Gpt);
	ntGen->SetBranchAddress("GisSignal",GisSignal);
	ntGen->SetBranchAddress("GcollisionId",GcollisionId);
	ntGen->SetBranchAddress("GpdgId",GpdgId);
	ntGen->SetBranchAddress("Gmu1pt",Gmu1pt);
	ntGen->SetBranchAddress("Gmu1eta",Gmu1eta);
	ntGen->SetBranchAddress("Gmu1phi",Gmu1phi);
	ntGen->SetBranchAddress("Gmu2pt",Gmu2pt);
	ntGen->SetBranchAddress("Gmu2eta",Gmu2eta);
	ntGen->SetBranchAddress("Gmu2phi",Gmu2phi);


	ntGen->SetBranchAddress("Gtk1pt",Gtk1pt);
	ntGen->SetBranchAddress("Gtk1eta",Gtk1eta);
	ntGen->SetBranchAddress("Gtk1phi",Gtk1phi);

	ntGen->SetBranchAddress("Gtk2pt",Gtk2pt);
	ntGen->SetBranchAddress("Gtk2eta",Gtk2eta);
	ntGen->SetBranchAddress("Gtk2phi",Gtk2phi);



	Int_t HLT_HIL1DoubleMu0_v1;
	ntHlt->SetBranchAddress("HLT_HIL1DoubleMu0_v1",&HLT_HIL1DoubleMu0_v1);


	double CentWeight;

	//CentWeightTree->SetBranchAddress("CentWeight",&CentWeight);


	double muid1[NCand];
	double muid2[NCand];
	double trk1[NCand];
	double trk2[NCand];
	double trg1[NCand];
	double trg2[NCand];

	/*
	   TnPInfo->SetBranchAddress("muid1",muid1);
	   TnPInfo->SetBranchAddress("muid2",muid2);
	   TnPInfo->SetBranchAddress("trk1",trk1);
	   TnPInfo->SetBranchAddress("trk2",trk2);
	   TnPInfo->SetBranchAddress("trg1",trg1);
	   TnPInfo->SetBranchAddress("trg2",trg2);
	   */




	//Syst Purpose//
	double muid1statup[NCand];
	double trk1statup[NCand];
	double trg1statup[NCand];

	double muid1statdown[NCand];
	double trk1statdown[NCand];
	double trg1statdown[NCand];


	double muid1systup[NCand];
	double trk1systup[NCand];
	double trg1systup[NCand];

	double muid1systdown[NCand];
	double trk1systdown[NCand];
	double trg1systdown[NCand];


	double muid2statup[NCand];
	double trk2statup[NCand];
	double trg2statup[NCand];

	double muid2statdown[NCand];
	double trk2statdown[NCand];
	double trg2statdown[NCand];


	double muid2systup[NCand];
	double trk2systup[NCand];
	double trg2systup[NCand];

	double muid2systdown[NCand];
	double trk2systdown[NCand];
	double trg2systdown[NCand];



	double muid1syst;
	double muid1stat;
	double muid2syst;
	double muid2stat;

	double trk1syst;
	double trk1stat;
	double trk2syst;
	double trk2stat;

	double trg1syst;
	double trg1stat;
	double trg2syst;
	double trg2stat;


	double tnptotal1syst;
	double tnptotal1stat;


	double tnptotal2syst;
	double tnptotal2stat;


	double tnptotal1err;
	double tnptotal2err;

	double tnptotalerr;


	/*



	   TnPInfo->SetBranchAddress("muid1statup",muid1statup);
	   TnPInfo->SetBranchAddress("trk1statup",trk1statup);
	   TnPInfo->SetBranchAddress("trg1statup",trg1statup);
	   TnPInfo->SetBranchAddress("muid1statdown",muid1statdown);
	   TnPInfo->SetBranchAddress("trk1statdown",trk1statdown);
	   TnPInfo->SetBranchAddress("trg1statdown",trg1statdown);
	   TnPInfo->SetBranchAddress("muid1systup",muid1systup);
	   TnPInfo->SetBranchAddress("trk1systup",trk1systup);
	   TnPInfo->SetBranchAddress("trg1systup",trg1systup);
	   TnPInfo->SetBranchAddress("muid1systdown",muid1systdown);
	   TnPInfo->SetBranchAddress("trk1systdown",trk1systdown);
	   TnPInfo->SetBranchAddress("trg1systdown",trg1systdown);



	   TnPInfo->SetBranchAddress("muid2statup",muid2statup);
	   TnPInfo->SetBranchAddress("trk2statup",trk2statup);
	   TnPInfo->SetBranchAddress("trg2statup",trg2statup);
	   TnPInfo->SetBranchAddress("muid2statdown",muid2statdown);
	   TnPInfo->SetBranchAddress("trk2statdown",trk2statdown);
	   TnPInfo->SetBranchAddress("trg2statdown",trg2statdown);
	   TnPInfo->SetBranchAddress("muid2systup",muid2systup);
	   TnPInfo->SetBranchAddress("trk2systup",trk2systup);
	   TnPInfo->SetBranchAddress("trg2systup",trg2systup);
	   TnPInfo->SetBranchAddress("muid2systdown",muid2systdown);
	   TnPInfo->SetBranchAddress("trk2systdown",trk2systdown);
	   TnPInfo->SetBranchAddress("trg2systdown",trg2systdown);
	   */




	int NEvents = ntphi->GetEntries();
	
	const int yBinN = 5;
	double yBinning[yBinN+1] = {0.0,0.5, 1.0, 1.5,2.0, 2.4};

	double LowBinWidth = 0.5;

	int NLowBin = 10/LowBinWidth;
	//	int NLowBin = 5;

	double MidBinWidth = 1;
	int NMidBin = 10/MidBinWidth;
	double HighBinWidth = 1;
	int NHighBin = 30/HighBinWidth;
	const int BptBin = NHighBin + NMidBin + NLowBin;
	double BptBinning[BptBin + 1];




	for(int i = 0; i < NLowBin; i++){
		BptBinning[i] = 0 + i * LowBinWidth;
	}
	for(int i = 0; i < NMidBin; i++){
		BptBinning[i+NLowBin] = 10 + i * MidBinWidth;
	}
	for(int i = 0; i <  NHighBin+1; i++){
		BptBinning[i+NLowBin+NMidBin] = 20 + i * HighBinWidth;
	}



	double PVzWeight;


	double EventWeight;
	double TnPWeight;
	double muidWeight;
	double trkWeight;
	double TotalWeight;
	double muidtrkWeight;


	double TotalWeightSystUp;
	double TotalWeightSystDown;

	double TotalWeightMuidUp;
	double TotalWeightMuidDown;
	double TotalWeightTrkUp;
	double TotalWeightTrkDown;
	double TotalWeightTrgUp;
	double TotalWeightTrgDown;


	double muid1total;
	double muid2total;

	double trk1total;
	double trk2total;

	double trg1total;
	double trg2total;


	double muidtotalerr;
	double trktotalerr;
	double trgtotalerr;


	TH2D * NoWeightHis = new TH2D("NoWeightHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * EvtWeightHis = new TH2D("EvtWeightHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * muidWeightHis = new TH2D("muidWeightHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * trkWeightHis = new TH2D("trkWeightHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * muidtrkWeightHis = new TH2D("muidtrkWeightHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHis = new TH2D("TnPWeightHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHisSystUp = new TH2D("TnPWeightHisSystUp","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHisSystDown = new TH2D("TnPWeightHisSystDown","",BptBin,BptBinning,yBinN,yBinning);



	TH2D * TnPWeightHisMuidUp = new TH2D("TnPWeightHisMuidUp","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHisMuidDown = new TH2D("TnPWeightHisMuidDown","",BptBin,BptBinning,yBinN,yBinning);

	TH2D * TnPWeightHisTrkUp = new TH2D("TnPWeightHisTrkUp","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHisTrkDown = new TH2D("TnPWeightHisTrkDown","",BptBin,BptBinning,yBinN,yBinning);


	TH2D * TnPWeightHisTrgUp = new TH2D("TnPWeightHisTrgUp","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * TnPWeightHisTrgDown = new TH2D("TnPWeightHisTrgDown","",BptBin,BptBinning,yBinN,yBinning);

	//Gen//
	TH2D * NoWeightGenHis = new TH2D("NoWeightGenHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * EvtWeightGenHis = new TH2D("EvtWeightGenHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * EvtWeightGenAccHis = new TH2D("EvtWeightGenAccHis","",BptBin,BptBinning,yBinN,yBinning);
	TH2D * NoWeightGenAccHis = new TH2D("NoWeightGenAccHis","",BptBin,BptBinning,yBinN,yBinning);



	TH1D * Bmu1ptHis = new TH1D("Bmu1ptHis","",200,0,50);
	Bmu1ptHis->GetXaxis()->SetTitle("Bmu1pt (GeV/c)");	
	Bmu1ptHis->GetYaxis()->SetTitle("Counts");
	Bmu1ptHis->GetXaxis()->CenterTitle();	
	Bmu1ptHis->GetYaxis()->CenterTitle();
	Bmu1ptHis->GetXaxis()->SetTitleOffset(1.2);	
	Bmu1ptHis->GetYaxis()->SetTitleOffset(1.5);



	const int NPtBins = 10;
	double PtBin[NPtBins + 1] = {0,1,2,3,5,7,10,15,20,50,100};
	

	const int NPtBins1D = 10;
	double  PtBin1D[NPtBins1D + 1] = {0,2,3,5,7,10,15,20,30,50,100};



	//	TH1D * Eff1DRECOHis = new TH1D("Eff1DRECOHis","",NPtBins,PtBin);
	TH1D * Eff1DRECOHis = new TH1D("Eff1DRECOHis","",NPtBins1D,PtBin1D);

	Eff1DRECOHis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff1DRECOHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOHis->GetXaxis()->CenterTitle();	
	Eff1DRECOHis->GetYaxis()->CenterTitle();
	Eff1DRECOHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOHis->GetYaxis()->SetTitleOffset(1.5);

	const int NMultiBin = 10;
	double  MultiBin1D[NMultiBin + 1] = {0,15,25,30,35,40,50,65,80,100,130};



	TH1D * Eff1DRECOMultHis = new TH1D("Eff1DRECOMultHis","",NMultiBin,MultiBin1D);

	Eff1DRECOMultHis->GetXaxis()->SetTitle("Multiplicity");
	Eff1DRECOMultHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DRECOMultHis->GetXaxis()->CenterTitle();	
	Eff1DRECOMultHis->GetYaxis()->CenterTitle();
	Eff1DRECOMultHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DRECOMultHis->GetYaxis()->SetTitleOffset(1.5);



//	TH1D * Eff1DGENHis = new TH1D("Eff1DGENHis","",NPtBins,PtBin);
	TH1D * Eff1DGENHis = new TH1D("Eff1DGENHis","",NPtBins1D,PtBin1D);
	
	Eff1DGENHis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	Eff1DGENHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DGENHis->GetXaxis()->CenterTitle();	
	Eff1DGENHis->GetYaxis()->CenterTitle();
	Eff1DGENHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DGENHis->GetYaxis()->SetTitleOffset(1.5);

	TH1D * Eff1DGENMultHis = new TH1D("Eff1DGENMultHis","",NMultiBin,MultiBin1D);

	Eff1DGENMultHis->GetXaxis()->SetTitle("Multiplicity");
	Eff1DGENMultHis->GetYaxis()->SetTitle("#alpha #times #epsilon");
	Eff1DGENMultHis->GetXaxis()->CenterTitle();	
	Eff1DGENMultHis->GetYaxis()->CenterTitle();
	Eff1DGENMultHis->GetXaxis()->SetTitleOffset(1.2);	
	Eff1DGENMultHis->GetYaxis()->SetTitleOffset(1.5);



	TH1D * Bmu2ptHis = new TH1D("Bmu2ptHis","",200,0,50);
	Bmu2ptHis->GetXaxis()->SetTitle("Bmu2pt (GeV/c)");	
	Bmu2ptHis->GetYaxis()->SetTitle("Counts");
	Bmu2ptHis->GetXaxis()->CenterTitle();	
	Bmu2ptHis->GetYaxis()->CenterTitle();
	Bmu2ptHis->GetXaxis()->SetTitleOffset(1.2);	
	Bmu2ptHis->GetYaxis()->SetTitleOffset(1.5);

	const int NPtbin = 200;
	double PtMin = 0;
	double PtMax = 100;
	double PtStep = (PtMax - PtMin)/NPtbin;

	const int NEtabin = 4;
	double Etabinning[NEtabin + 1] = {0,1.2,1.8,2.1,2.4};



	TH2D * Bmu1pteta = new TH2D("Bmu1pteta","",NPtbin,PtMin,PtMax,NEtabin,Etabinning);
	Bmu1pteta->GetXaxis()->SetTitle("Bmu1pt");	
	Bmu1pteta->GetYaxis()->SetTitle("Bmu1eta");
	Bmu1pteta->GetXaxis()->CenterTitle();	
	Bmu1pteta->GetYaxis()->CenterTitle();
	Bmu1pteta->GetXaxis()->SetTitleOffset(1.2);	
	Bmu1pteta->GetYaxis()->SetTitleOffset(1.5);



	TH2D * Bmu2pteta = new TH2D("Bmu2pteta","",NPtbin,PtMin,PtMax,NEtabin,Etabinning);
	Bmu2pteta->GetXaxis()->SetTitle("Bmu2pt");	
	Bmu2pteta->GetYaxis()->SetTitle("Bmu2eta");
	Bmu2pteta->GetXaxis()->CenterTitle();	
	Bmu2pteta->GetYaxis()->CenterTitle();
	Bmu2pteta->GetXaxis()->SetTitleOffset(1.2);	
	Bmu2pteta->GetYaxis()->SetTitleOffset(1.5);



	TH1D * Bmu1etaHis = new TH1D("Bmu1etaHis","",100,0,3);
	Bmu1etaHis->GetXaxis()->SetTitle("Bmu1eta");	
	Bmu1etaHis->GetYaxis()->SetTitle("Counts");
	Bmu1etaHis->GetXaxis()->CenterTitle();	
	Bmu1etaHis->GetYaxis()->CenterTitle();
	Bmu1etaHis->GetXaxis()->SetTitleOffset(1.2);	
	Bmu1etaHis->GetYaxis()->SetTitleOffset(1.5);




	TH1D * Bmu2etaHis = new TH1D("Bmu2etaHis","",100,0,3);
	Bmu2etaHis->GetXaxis()->SetTitle("Bmu2eta");	
	Bmu2etaHis->GetYaxis()->SetTitle("Counts");
	Bmu2etaHis->GetXaxis()->CenterTitle();	
	Bmu2etaHis->GetYaxis()->CenterTitle();
	Bmu2etaHis->GetXaxis()->SetTitleOffset(1.2);	
	Bmu2etaHis->GetYaxis()->SetTitleOffset(1.5);




	TH1D * Bmu1TrgSF = new TH1D("Bmu1TrgSF","",100,0.7,1.3);
	Bmu1TrgSF->GetXaxis()->SetTitle("Bmu1TrgSF");	
	Bmu1TrgSF->GetYaxis()->SetTitle("Counts");
	Bmu1TrgSF->GetXaxis()->CenterTitle();	
	Bmu1TrgSF->GetYaxis()->CenterTitle();
	Bmu1TrgSF->GetXaxis()->SetTitleOffset(1.2);	
	Bmu1TrgSF->GetYaxis()->SetTitleOffset(1.5);


	TH1D * Bmu2TrgSF = new TH1D("Bmu2TrgSF","",100,0.7,1.3);
	Bmu2TrgSF->GetXaxis()->SetTitle("Bmu2TrgSF");	
	Bmu2TrgSF->GetYaxis()->SetTitle("Counts");
	Bmu2TrgSF->GetXaxis()->CenterTitle();	
	Bmu2TrgSF->GetYaxis()->CenterTitle();
	Bmu2TrgSF->GetXaxis()->SetTitleOffset(1.2);	
	Bmu2TrgSF->GetYaxis()->SetTitleOffset(1.5);






	//	TFile * fout = new TFile("WeightCheck.root","RECREATE");

	TString outfileName = "NewEff2DMaps/EffFineBDT.root";
	TFile * fout = new TFile(outfileName.Data(),"RECREATE");

	fout->cd();


	//NEvents = 20;

	for(int i = 0; i < NEvents; i++){


		if(i%10000==0) cout << "Now Working on = " << i  <<  endl;

	

	Float_t Gtk1pt[NCand];
	Float_t Gtk1eta[NCand];
	Float_t Gtk1phi[NCand];
		ntphi->GetEntry(i);
		ntSkim->GetEntry(i);
		ntHi->GetEntry(i);
		ntHlt->GetEntry(i);
		//	BDT->GetEntry(i);
		//	CentWeightTree->GetEntry(i);
		//	TnPInfo->GetEntry(i);
		ntGen->GetEntry(i);

		//BDT//
		BDT1->GetEntry(i);
		BDT2->GetEntry(i);
		BDT3->GetEntry(i);
		BDT4->GetEntry(i);
		BDT5->GetEntry(i);
		BDT6->GetEntry(i);
		BDT7->GetEntry(i);
		BDT8->GetEntry(i);
	
		root->GetEntry(i);

		//cout << "Bsize = " << Bsize << endl;


		for(int j = 0; j < Bsize; j++){
			//cout << "j = " << j << endl;
			//	if( (Bgen[j] == 23333)  && (((BDT_pt_1_2[j]>0.04 && BsvpvDistance[j]/BsvpvDisErr[j] > 5 && Bchi2cl[j] > 0.05 && Bpt[j] > 1.0 && Bpt[j] < 2.0) || (BDT_pt_1_2[j]>0.04 && BsvpvDistance[j]/BsvpvDisErr[j] > 5 && Bchi2cl[j] > 0.05 && Bpt[j] > 0.0 && Bpt[j] < 1.0) )   || (((HBHENoiseFilterResult == 1 && pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1)  &&  (Bmu1isTriggered[j] == 1 && Bmu2isTriggered[j] == 1 ) &&  (Btrk1Pt[j] > 0.2 && Bchi2cl[j] > 0.05 && BsvpvDistance[j]/BsvpvDisErr[j] > 2.0 && Bpt[j] > 2 && abs(Btrk1Eta[j]-0.0) < 2.4  && (TMath::Abs(By[j])<2.4&&TMath::Abs(Bmumumass[j]-3.096916)<0.15&&((abs(Bmu1eta[j])<1.2&&Bmu1pt[j]>3.5)||(abs(Bmu1eta[j])>1.2&&abs(Bmu1eta[j])<2.1&&Bmu1pt[j]>(5.47-1.89*abs(Bmu1eta[j])))||(abs(Bmu1eta[j])>2.1&&abs(Bmu1eta[j])<2.4&&Bmu1pt[j]>1.5))&&((abs(Bmu2eta[j])<1.2&&Bmu2pt[j]>3.5)||(abs(Bmu2eta[j])>1.2&&abs(Bmu2eta[j])<2.1&&Bmu2pt[j]>(5.47-1.89*abs(Bmu2eta[j])))||(abs(Bmu2eta[j])>2.1&&abs(Bmu2eta[j])<2.4&&Bmu2pt[j]>1.5))&&Bmu1TMOneStationTight[j]&&Bmu2TMOneStationTight[j]&&Bmu1InPixelLayer[j]>0&&(Bmu1InPixelLayer[j]+Bmu1InStripLayer[j])>5&&Bmu2InPixelLayer[j]>0&&(Bmu2InPixelLayer[j]+Bmu2InStripLayer[j])>5&&Bmu1dxyPV[j]<0.3&&Bmu2dxyPV[j]<0.3&&Bmu1dzPV[j]<20&&Bmu2dzPV[j]<20&&Bmu1isTrackerMuon[j]&&Bmu2isTrackerMuon[j]&&Bmu1isGlobalMuon[j]&&Bmu2isGlobalMuon[j]&&Btrk1highPurity[j]&&abs(Btrk1Eta[j])<2.4&&Btrk1Pt[j]>0.2)  && (Btrk1PixelHit[j] + Btrk1StripHit[j] > 10) &&  (Btrk1PtErr[j]/Btrk1Pt[j] < 0.1)&& Btrk1Chi2ndf[j]/(Btrk1nStripLayer[j]+Btrk1nPixelLayer[j]) < 0.18   && (abs(PVz)<15))  &&( (Bpt[j] > 3 && Bpt[j] < 5 && BDT_pt_3_5[j] > 0.08) || (Bpt[j] > 5 && Bpt[j] < 7 && BDT_pt_5_7[j] > 0.06) || (Bpt[j] > 7 && Bpt[j] < 10 && BDT_pt_7_10[j] > 0.07) || (Bpt[j] > 10 && Bpt[j] < 15 && BDT_pt_10_15[j] > 0.08) || (Bpt[j] > 15 && Bpt[j] < 20 && BDT_pt_15_20[j] > 0.12)  || (Bpt[j] > 20 && Bpt[j] < 50 && BDT_pt_20_50[j] > 0.12) || (Bpt[j] > 50 && Bpt[j] < 100))) || ((HBHENoiseFilterResult == 1 && pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1) &&  (Bpt[j] > 2 && Bpt[j] < 3 && BDT_pt_2_3[j] > -0.12) ) ))){
			if( (Bgen[j] == 23333)  &&((   (((abs(Btktkmass[j]-1.019455)<0.015)&& BDT_pt_2_3[j] > -0.10 && TMath::Abs(Bmumumass[j]-3.096916)<0.15 && Bpt[j] > 0 && Bpt[j] < 5 && (abs(Btrk1Eta[j])<2.4 && abs(Btrk2Eta[j])<2.4 && Btrk1Pt[j]>0.0 && Btrk2Pt[j]>0.0) && Btrk1Pt[j] > 0.2 && Btrk2Pt[j] > 0.2  && Bchi2cl[j] > 0.05 && BsvpvDistance[j]/BsvpvDisErr[j] > 2.0)  && ( (Bpt[j] < 2 && Bpt[j] > 0 && BDT_pt_1_2[j] > -0.38 ) || (Bpt[j] < 3 && Bpt[j] > 2 && BDT_pt_2_3[j] > -0.05 ) || (Bpt[j] < 5 && Bpt[j] > 3 && BDT_pt_3_5[j] > -0.40)  ))) || ((  (Bpt[j] > 50) || (((Bpt[j] > 5 && Bpt[j] < 7 &&  BDT_pt_5_7[j] > -0.28) || (Bpt[j] > 7 && Bpt[j] < 10 &&  BDT_pt_7_10[j] > -0.25) || (Bpt[j] > 10 && Bpt[j] < 15 &&  BDT_pt_10_15[j] > -0.23) || (Bpt[j] > 15 && Bpt[j] < 20 &&  BDT_pt_15_20[j] > -0.20 ) || (Bpt[j] > 20 && Bpt[j] < 50 &&  BDT_pt_30_50[j] > -0.25 ) ) && (Btrk1Pt[j] > 0.0 && Btrk2Pt[j] > 0.0 && Bchi2cl[j] > 0.05 && BsvpvDistance[j]/BsvpvDisErr[j] > 2.0 && abs(Btrk1Eta[j]-0.0) < 2.4 && abs(Btrk2Eta[j]-0.0) < 2.4 && (TMath::Abs(By[j])<2.4&&TMath::Abs(Bmumumass[j]-3.096916)<0.15&&((abs(Bmu1eta[j])<1.2&&Bmu1pt[j]>3.5)||(abs(Bmu1eta[j])>1.2&&abs(Bmu1eta[j])<2.1&&Bmu1pt[j]>(5.47-1.89*abs(Bmu1eta[j])))||(abs(Bmu1eta[j])>2.1&&abs(Bmu1eta[j])<2.4&&Bmu1pt[j]>1.5))&&((abs(Bmu2eta[j])<1.2&&Bmu2pt[j]>3.5)||(abs(Bmu2eta[j])>1.2&&abs(Bmu2eta[j])<2.1&&Bmu2pt[j]>(5.47-1.89*abs(Bmu2eta[j])))||(abs(Bmu2eta[j])>2.1&&abs(Bmu2eta[j])<2.4&&Bmu2pt[j]>1.5))&&Bmu1TMOneStationTight[j]&&Bmu2TMOneStationTight[j]&&Bmu1InPixelLayer[j]>0&&(Bmu1InPixelLayer[j]+Bmu1InStripLayer[j])>5&&Bmu2InPixelLayer[j]>0&&(Bmu2InPixelLayer[j]+Bmu2InStripLayer[j])>5&&Bmu1dxyPV[j]<0.3&&Bmu2dxyPV[j]<0.3&&Bmu1dzPV[j]<20&&Bmu2dzPV[j]<20&&Bmu1isTrackerMuon[j]&&Bmu2isTrackerMuon[j]&&Bmu1isGlobalMuon[j]&&Bmu2isGlobalMuon[j]&&Btrk1highPurity[j]&&Btrk2highPurity[j]&&abs(Btrk1Eta[j])<2.4&&abs(Btrk2Eta[j])<2.4&&Btrk1Pt[j]>0.&&Btrk2Pt[j]>0.&&abs(Btktkmass[j]-1.019455)<0.015) && (abs(PVz)<15) && (Btrk1PixelHit[j] + Btrk1StripHit[j] > 10) && (Btrk2PixelHit[j] + Btrk2StripHit[j] > 10) && (Btrk1PtErr[j]/Btrk1Pt[j] < 0.1)&& (Btrk2PtErr[j]/Btrk2Pt[j] < 0.1) && Btrk1Chi2ndf[j]/(Btrk1nStripLayer[j]+Btrk1nPixelLayer[j]) < 0.18 && Btrk2Chi2ndf[j]/(Btrk2nStripLayer[j]+Btrk2nPixelLayer[j]) < 0.18   && (abs(PVz)<15))))))){
				PVzWeight = 1;
				//EventWeight = pthat * weight;
				EventWeight = weight;
				//Turn off TnP weight for now
				muid1[j] = 1;
				muid2[j] = 1;	Float_t Gtk1pt[NCand];
	Float_t Gtk1eta[NCand];
	Float_t Gtk1phi[NCand];
				trk1[j] = 1;
				trk2[j] = 1;
				trg1[j] = 1;
				trg2[j] = 1;
		
				TnPWeight = muid1[j] * trk1[j] * trg1[j] * muid2[j] * trk2[j] * trg2[j];

				muidWeight = EventWeight * muid1[j] * muid2[j];
				trkWeight = EventWeight * trk1[j] * trk2[j];

				muidtrkWeight = EventWeight * muid1[j] * muid2[j] * trk1[j] * trk2[j];
				TotalWeight = EventWeight * TnPWeight;




				NoWeightHis->Fill(Bpt[j],abs(By[j]),1);
				EvtWeightHis->Fill(Bpt[j],abs(By[j]),EventWeight);
				muidWeightHis->Fill(Bpt[j],abs(By[j]),muidWeight);
				trkWeightHis->Fill(Bpt[j],abs(By[j]),trkWeight);
				muidtrkWeightHis->Fill(Bpt[j],abs(By[j]),muidtrkWeight);
				TnPWeightHis->Fill(Bpt[j],abs(By[j]),TotalWeight);
	
				Eff1DRECOHis->Fill(Bpt[j],TotalWeight);
	
				Eff1DRECOMultHis->Fill(nMult,TotalWeight);


				//Now Everything is About the Error//
				if(muid1systup[j] >= muid1systdown[j]) muid1syst = muid1systup[j];
				if(muid1systdown[j] > muid1systup[j]) muid1syst = muid1systdown[j];
				if(muid2systup[j] >= muid2systdown[j]) muid2syst = muid2systup[j];
				if(muid2systdown[j] > muid2systup[j]) muid2syst = muid2systdown[j];

				if(muid1statup[j] >= muid1statdown[j]) muid1stat = muid1statup[j];
				if(muid1statdown[j] > muid1statup[j]) muid1stat = muid1statdown[j];
				if(muid2statup[j] >= muid2statdown[j]) muid2stat = muid2statup[j];
				if(muid2statdown[j] > muid2statup[j]) muid2stat = muid2statdown[j];

				//Trk
				if(trk1systup[j] >= trk1systdown[j]) trk1syst = trk1systup[j];
				if(trk1systdown[j] > trk1systup[j]) trk1syst = trk1systdown[j];
				if(trk2systup[j] >= trk2systdown[j]) trk2syst = trk2systup[j];
				if(trk2systdown[j] > trk2systup[j]) trk2syst = trk2systdown[j];

				if(trk1statup[j] >= trk1statdown[j]) trk1stat = trk1statup[j];
				if(trk1statdown[j] > trk1statup[j]) trk1stat = trk1statdown[j];
				if(trk2statup[j] >= trk2statdown[j]) trk2stat = trk2statup[j];
				if(trk2statdown[j] > trk2statup[j]) trk2stat = trk2statdown[j];



				//Trg
				if(trg1systup[j] >= trg1systdown[j]) trg1syst = trg1systup[j];
				if(trg1systdown[j] > trg1systup[j]) trg1syst = trg1systdown[j];
				if(trg2systup[j] >= trg2systdown[j]) trg2syst = trg2systup[j];
				if(trg2systdown[j] > trg2systup[j]) trg2syst = trg2systdown[j];

				if(trg1statup[j] >= trg1statdown[j]) trg1stat = trg1statup[j];
				if(trg1statdown[j] > trg1statup[j]) trg1stat = trg1statdown[j];
				if(trg2statup[j] >= trg2statdown[j]) trg2stat = trg2statup[j];
				if(trg2statdown[j] > trg2statup[j]) trg2stat = trg2statdown[j];


				tnptotal1syst = sqrt(muid1syst/muid1[j] * muid1syst/muid1[j] + trk1syst/trk1[j] * trk1syst/trk1[j] + trg1syst/trg1[j] * trg1syst/trg1[j]);
				tnptotal1stat = sqrt(muid1stat/muid1[j] * muid1stat/muid1[j] + trk1stat/trk1[j] * trk1stat/trk1[j] + trg1stat/trg1[j] * trg1stat/trg1[j]);


				tnptotal2syst = sqrt(muid2syst/muid2[j] * muid2syst/muid2[j] + trk2syst/trk2[j] * trk2syst/trk2[j] + trg2syst/trg2[j] * trg2syst/trg2[j]);
				tnptotal2stat = sqrt(muid2stat/muid2[j] * muid2stat/muid2[j] + trk2stat/trk2[j] * trk2stat/trk2[j] + trg2stat/trg2[j] * trg2stat/trg2[j]);


				tnptotal1err = sqrt(tnptotal1stat * tnptotal1stat + tnptotal1syst * tnptotal1syst);
				tnptotal2err = sqrt(tnptotal2stat * tnptotal2stat + tnptotal2syst * tnptotal2syst);	

				tnptotalerr = sqrt(tnptotal1err * tnptotal1err + tnptotal2err * tnptotal2err);

				//	cout << "tnptotalerr = " << tnptotalerr << endl;

				TotalWeightSystUp = TotalWeight * ( 1 + tnptotalerr);
				TotalWeightSystDown = TotalWeight * ( 1 - tnptotalerr);


				TnPWeightHisSystUp->Fill(Bpt[j],abs(By[j]),TotalWeightSystUp);
				TnPWeightHisSystDown->Fill(Bpt[j],abs(By[j]),TotalWeightSystDown);


				muid1total =  sqrt(muid1syst/muid1[j] * muid1syst/muid1[j] + muid1stat/muid1[j] * muid1stat/muid1[j]);
				muid2total =  sqrt(muid2syst/muid2[j] * muid2syst/muid2[j] + muid2stat/muid2[j] * muid2stat/muid2[j]);

				trk1total =  sqrt(trk1syst/trk1[j] * trk1syst/trk1[j] + trk1stat/trk1[j] * trk1stat/trk1[j]);
				trk2total =	 sqrt(trk2syst/trk2[j] * trk2syst/trk2[j] + trk2stat/trk2[j] * trk2stat/trk2[j]);

				trg1total =  sqrt(trg1syst/trg1[j] * trg1syst/trg1[j] + trg1stat/trg1[j] * trg1stat/trg1[j]);
				trg2total =	sqrt(trg2syst/trg2[j] * trg2syst/trg2[j] + trg2stat/trg2[j] * trg2stat/trg2[j]);


				muidtotalerr = sqrt(muid1total * muid1total + muid2total * muid2total);
				trktotalerr = sqrt(trk1total * trk1total + trk2total * trk2total);
				trgtotalerr = sqrt(trg1total * trg1total + trg2total * trg2total);

				//		cout << "muid1total = " << muid1total << "    muid2total = " <<muid2total << "  muidtotalerr = " << muidtotalerr  << endl;


				TotalWeightMuidUp = TotalWeight * ( 1 + muidtotalerr);
				TotalWeightMuidDown = TotalWeight * ( 1 - muidtotalerr);

				TotalWeightTrkUp = TotalWeight * ( 1 + trktotalerr);
				TotalWeightTrkDown = TotalWeight * ( 1 - trktotalerr);

				TotalWeightTrgUp = TotalWeight * ( 1 + trgtotalerr);
				TotalWeightTrgDown = TotalWeight * ( 1 - trgtotalerr);

				//For the Moment//

				TotalWeightMuidUp = 1;
				TotalWeightMuidDown = 1;


				TotalWeightTrkUp = 1;
				TotalWeightTrkDown = 1;

				TotalWeightTrgUp = 1;
				TotalWeightTrgDown = 1;

				TnPWeightHisMuidUp->Fill(Bpt[j],abs(By[j]),TotalWeightMuidUp);
				TnPWeightHisMuidDown->Fill(Bpt[j],abs(By[j]),TotalWeightMuidDown);


				//cout << "TotalWeightMuidUp = " << TotalWeightMuidUp  << "   TotalWeightMuidDown = " << TotalWeightMuidDown << endl;


				TnPWeightHisTrkUp->Fill(Bpt[j],abs(By[j]),TotalWeightTrkUp);
				TnPWeightHisTrkDown->Fill(Bpt[j],abs(By[j]),TotalWeightTrkDown);



				TnPWeightHisTrgUp->Fill(Bpt[j],abs(By[j]),TotalWeightTrgUp);
				TnPWeightHisTrgDown->Fill(Bpt[j],abs(By[j]),TotalWeightTrgDown);


				if(Bpt[j] < ptmax && Bpt[j] > ptmin){
					Bmu1ptHis->Fill(Bmu1pt[j],TotalWeight);
					Bmu2ptHis->Fill(Bmu2pt[j],TotalWeight);
					Bmu1etaHis->Fill(abs(Bmu1eta[j]),TotalWeight);
					Bmu2etaHis->Fill(abs(Bmu2eta[j]),TotalWeight);
					Bmu1TrgSF->Fill(trg1[j],TotalWeight);
					Bmu2TrgSF->Fill(trg2[j],TotalWeight);

					Bmu1pteta->Fill(Bmu1pt[j],abs(Bmu1eta[j]),TotalWeight);
					Bmu2pteta->Fill(Bmu2pt[j],abs(Bmu2eta[j]),TotalWeight);
				}

			}





		}



		}




		cout << "Now Loop Gen" << endl;



		for(int i = 0; i < NEvents; i++){


			ntGen->GetEntry(i);
			ntHi->GetEntry(i);
			//CentWeightTree->GetEntry(i);
			ntphi->GetEntry(i);
			root->GetEntry(i);

			PVzWeight = 1;

			CentWeight = 1;

			//	PVzWeight = (0.163562 * TMath::Exp(- 0.021039 * (PVz - 0.426587)*(PVz - 0.426587)))/(0.159629 * TMath::Exp(- 0.020014 * (PVz - 0.589381)*(PVz - 0.589381)));

			//	PVzWeight = (TMath::Gaus(PVz,0.432315,4.874300)/(sqrt(2*3.14159)*4.874300))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989));
			//	EventWeight = PVzWeight * CentWeight * pthat * weight;
			
			EventWeight = PVzWeight * CentWeight * weight;
			//cout << "GSize = " << Gsize << endl;


			for(int j = 0; j < Gsize; j++){



				//cout << "Gpt = " << Gpt[j] << "    GpdgId[j] = " << GpdgId[j]  << endl;





				if((TMath::Abs(Gy[j])<2.4 && TMath::Abs(GpdgId[j])==531 && GisSignal[j]>0 ) ){
					

					NoWeightGenHis->Fill(Gpt[j],abs(Gy[j]),1);
					EvtWeightGenHis->Fill(Gpt[j],abs(Gy[j]),EventWeight);
					Eff1DGENHis->Fill(Gpt[j],EventWeight);
					Eff1DGENMultHis->Fill(nMult,EventWeight);
				}

				if((TMath::Abs(Gy[j])<2.4 && TMath::Abs(GpdgId[j])==531 && GisSignal[j]==1 && GcollisionId[j]==0 && TMath::Abs(Gtk1eta[j])<2.4 && TMath::Abs(Gtk2eta[j])<2.4  && ((TMath::Abs(Gmu1eta[j])<1.2 && Gmu1pt[j]>3.5) || (TMath::Abs(Gmu1eta[j])>1.2 && TMath::Abs(Gmu1eta[j])<2.1 && Gmu1pt[j]>5.47-1.89*TMath::Abs(Gmu1eta[j])) || (TMath::Abs(Gmu1eta[j])>2.1 && TMath::Abs(Gmu2eta[j])<2.4 && Gmu1pt[j]>1.5)) && ((TMath::Abs(Gmu2eta[j])<1.2 && Gmu2pt[j]>3.5) || (TMath::Abs(Gmu2eta[j])>1.2 && TMath::Abs(Gmu2eta[j])<2.1 && Gmu2pt[j]>5.47-1.89*TMath::Abs(Gmu2eta[j])) || (TMath::Abs(Gmu2eta[j])>2.1 && TMath::Abs(Gmu2eta[j])<2.4 && Gmu2pt[j]>1.5)) && ((Gpt[j]>5 && Gpt[j]<10 && TMath::Abs(Gy[j])>1.5) || (Gpt[j]>10))) ){

					NoWeightGenAccHis->Fill(Gpt[j],abs(Gy[j]),EventWeight);
					EvtWeightGenAccHis->Fill(Gpt[j],abs(Gy[j]),EventWeight);

				}



			}

		}


		cout << "START MAKING HIS BRO" << endl;

		TH1D * Eff1DHis = (TH1D * ) Eff1DRECOHis->Clone("Eff1DHis");
		Eff1DHis->Sumw2();
		Eff1DGENHis->Sumw2();

		Eff1DHis->Divide(Eff1DGENHis);

		TH1D * Eff1DHisMult = (TH1D * ) Eff1DRECOMultHis->Clone("Eff1DHisMult");
		Eff1DHisMult->Sumw2();
		Eff1DGENMultHis->Sumw2();

		Eff1DHisMult->Divide(Eff1DGENMultHis);




		TH2D * invAcc2D = (TH2D * ) EvtWeightGenHis->Clone("invAcc2D");
		invAcc2D->Sumw2();
		invAcc2D->Divide(EvtWeightGenAccHis);

		TH2D * invEffonly2D = (TH2D * ) EvtWeightGenAccHis->Clone("invEffonly2D");
		invEffonly2D->Sumw2();
		invEffonly2D->Divide(TnPWeightHis);

		TH2D * invEff2D = (TH2D * ) EvtWeightGenHis->Clone("invEff2D");
		invEff2D->Sumw2();
		invEff2D->Divide(TnPWeightHis);

		//Systematics


		TH2D * invEff2DSystUp = (TH2D * ) EvtWeightGenHis->Clone("invEff2DSystUp");
		invEff2DSystUp->Sumw2();
		invEff2DSystUp->Divide(TnPWeightHisSystDown);

		TH2D * invEff2DSystDown = (TH2D * ) EvtWeightGenHis->Clone("invEff2DSystDown");
		invEff2DSystDown->Sumw2();
		invEff2DSystDown->Divide(TnPWeightHisSystUp);



		TH2D * invEff2DMuidUp = (TH2D * ) EvtWeightGenHis->Clone("invEff2DMuidUp");
		invEff2DMuidUp->Sumw2();
		invEff2DMuidUp->Divide(TnPWeightHisMuidDown);

		TH2D * invEff2DMuidDown = (TH2D * ) EvtWeightGenHis->Clone("invEff2DMuidDown");
		invEff2DMuidDown->Sumw2();
		invEff2DMuidDown->Divide(TnPWeightHisMuidUp);


		TH2D * invEff2DTrkUp = (TH2D * ) EvtWeightGenHis->Clone("invEff2DTrkUp");
		invEff2DTrkUp->Sumw2();
		invEff2DTrkUp->Divide(TnPWeightHisTrkDown);

		TH2D * invEff2DTrkDown = (TH2D * ) EvtWeightGenHis->Clone("invEff2DTrkDown");
		invEff2DTrkDown->Sumw2();
		invEff2DTrkDown->Divide(TnPWeightHisTrkUp);



		TH2D * invEff2DTrgUp = (TH2D * ) EvtWeightGenHis->Clone("invEff2DTrgUp");
		invEff2DTrgUp->Sumw2();
		invEff2DTrgUp->Divide(TnPWeightHisTrgDown);

		TH2D * invEff2DTrgDown = (TH2D * ) EvtWeightGenHis->Clone("invEff2DTrgDown");
		invEff2DTrgDown->Sumw2();
		invEff2DTrgDown->Divide(TnPWeightHisTrgUp);





		NoWeightHis->Write();
		EvtWeightHis->Write();
		muidWeightHis->Write();
		trkWeightHis->Write();
		muidtrkWeightHis->Write();
		TnPWeightHis->Write();

		NoWeightGenHis->Write();
		EvtWeightGenHis->Write();

		NoWeightGenAccHis->Write();
		EvtWeightGenAccHis->Write();

		invAcc2D->Write();
		invEffonly2D->Write();
		invEff2D->Write();


		invEff2DSystUp->Write();
		invEff2DSystDown->Write();

		invEff2DMuidUp->Write();
		invEff2DMuidDown->Write();

		invEff2DTrkUp->Write();
		invEff2DTrkDown->Write();

		invEff2DTrgUp->Write();
		invEff2DTrgDown->Write();
	

		Eff1DRECOHis->Write();
		Eff1DGENHis->Write();
		Eff1DHis->Write();


		Eff1DRECOMultHis->Write();
		Eff1DGENMultHis->Write();
		Eff1DHisMult->Write();


		TFile * fout2 = new TFile(Form("BPMuonInfoPlots_%d_%d.root",ptmin,ptmax),"RECREATE");
		fout2->cd();

		TCanvas *c = new TCanvas("c","c",600,600);
		c->cd();

		Eff1DHis->Draw("ep");
		c->SaveAs("1DEffPlots/Eff1DHis.png");

		for(int i = 0; i < NPtBins; i ++){

			cout << "pT =  " << PtBin[i] << " - " << PtBin[i+1] << "  Eff:  " << Eff1DHis->GetBinContent(i+1)   << "   Eff Error: " << Eff1DHis->GetBinError(i+1)<< endl;


		}

		Bmu1ptHis->Draw();
		c->SaveAs("MuonInfoPlots/Bmu1ptHis.png");

		Bmu2ptHis->Draw();
		c->SaveAs("MuonInfoPlots/Bmu2ptHis.png");


		Bmu1etaHis->Draw();
		c->SaveAs("MuonInfoPlots/Bmu1etaHis.png");

		Bmu2etaHis->Draw();
		c->SaveAs("MuonInfoPlots/Bmu2etaHis.png");

		Bmu1TrgSF->Draw();
		c->SaveAs("MuonInfoPlots/Bmu1TrgSF.png");

		Bmu2TrgSF->Draw();
		c->SaveAs("MuonInfoPlots/Bmu2TrgSF.png");

		Bmu1ptHis->Write();
		Bmu2ptHis->Write();
		Bmu1etaHis->Write();
		Bmu2etaHis->Write();
		Bmu1TrgSF->Write();
		Bmu2TrgSF->Write();

		Bmu1pteta->Write();
		Bmu2pteta->Write();
		fout2->Close();


		c->SetLogz();

		invEff2DSystUp->GetYaxis()->SetTitle("B |y|");
		invEff2DSystUp->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
		invEff2DSystUp->GetXaxis()->CenterTitle();
		invEff2DSystUp->GetYaxis()->CenterTitle();
		invEff2DSystUp->GetYaxis()->SetTitleOffset(1.2);
		invEff2DSystUp->SetTitle("");

		invEff2DSystUp->Draw("COLZ");
		c->SaveAs("Eff2DMapTnP/Eff2D_Up.png");




		invEff2DSystDown->GetYaxis()->SetTitle("B |y|");
		invEff2DSystDown->GetXaxis()->SetTitle("B_{s} p_{T} (GeV/c)");
		invEff2DSystDown->GetXaxis()->CenterTitle();
		invEff2DSystDown->GetYaxis()->CenterTitle();
		invEff2DSystDown->GetYaxis()->SetTitleOffset(1.2);
		invEff2DSystDown->SetTitle("");

		invEff2DSystDown->Draw("COLZ");
		c->SaveAs("Eff2DMapTnP/Eff2D_Down.png");


		fout->Close();
		fin->Close();







	}
