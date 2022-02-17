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


void ManualCutBootStrap(int CentMin, int CentMax, int Shape, int TwoShot, int PtOpt){


	int globalbin;

	TString FolderName;

	if(Shape == 16) FolderName = "NoWeight";
	if(Shape == 1) FolderName = "FONLL";
	if(Shape == 11) FolderName = "Linear";
	if(Shape == 12) FolderName = "Quadratic";
	if(Shape == 13) FolderName = "LInverse";
	if(Shape == 14) FolderName = "LSqrt";
	if(Shape == 15) FolderName = "LLog";
	if(Shape == 0) FolderName = "NoTnP";


	double ptahtCutValueUp;
	double ptahtCutValueDown;




	int BptLow;
	int BptHigh;

	int HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1;

	TString outfilept;


	if(PtOpt == -2){
		BptLow = 5;
		BptHigh = 50;
		outfilept = "5-50";
	}

	if(PtOpt == -1){
		BptLow = 10;
		BptHigh = 50;
		outfilept = "10-50";
	}

	if(PtOpt == 0){
		BptLow = 7;
		BptHigh = 10;
		outfilept = "7-10";

	}
	if(PtOpt == 1){
		BptLow = 10;
		BptHigh = 15;
		outfilept = "10-15";

	}
	if(PtOpt == 2){
		BptLow = 15;
		BptHigh = 20;
		outfilept = "15-20";

	}
	if(PtOpt == 3){
		BptLow = 20;
		BptHigh = 50;
		outfilept = "20-50";

	}



	TString FileName;
	//	FileName = "/export/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/2018Ana/Samples/FinalAnaSamples/PrivateMC-Data-Official-SemiFinal/Data_Bs_PbPb_TMVA_BDT_PbPb.root";

	FileName="/export/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/BPlusAna/InputFiles/ntKp_20190808_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_HF_and_MuonJSON_BDT.root";





	TFile * fin = new TFile(FileName.Data());
	TTree * ntphi = (TTree * ) fin->Get("Bfinder/ntKp");
	/*
	   TTree * BDT1 = (TTree * ) fin->Get("BDT_pt_5_10");
	   TTree * BDT2 = (TTree * ) fin->Get("BDT_pt_10_15");
	   TTree * BDT3 = (TTree * ) fin->Get("BDT_pt_15_20");
	   TTree * BDT4 = (TTree * ) fin->Get("BDT_pt_20_50");
	   */
	TTree * BDT = (TTree * ) fin->Get("BDT");

	TTree * ntHi = (TTree * ) fin->Get("hiEvtAnalyzer/HiTree");
	TTree * ntSkim = (TTree * ) fin->Get("skimanalysis/HltTree");
	TTree * ntHlt = (TTree *) fin->Get("hltanalysis/HltTree");




	TF1 * func;

	if(CentMin == 0 && CentMax == 90) func = new TF1("func","1/0.144708+TMath::Exp(-1.035696*(x-15.321432))+TMath::Exp(-0.204131*(x-30.289313))",5,50);

	if(CentMin == 0 && CentMax == 30) func = new TF1("func","1/0.131733+TMath::Exp(-1.051514*(x-15.468517))+TMath::Exp(-0.204272*(x-31.002324))",5,50);

	//	if(CentMin == 30 && CentMax == 40) func = new TF1("func","1/0.131733+TMath::Exp(-1.051514*(x-15.468517))+TMath::Exp(-0.204272*(x-31.002324))",5,50);
	//	if(CentMin == 40 && CentMax == 90) func = new TF1("func","1/0.131733+TMath::Exp(-1.051514*(x-15.468517))+TMath::Exp(-0.204272*(x-31.002324))",5,50);


	if(CentMin == 30 && CentMax == 90) func = new TF1("func","1/0.204307+TMath::Exp(-1.041731*(x-14.608514))+TMath::Exp(-0.206513*(x-27.599694))",5,50);





	const int NCand = 10000;



	int run;
	int lumi;
	int evt;
	int hiBin;
	Float_t PVz;
	Int_t pclusterCompatibilityFilter;
	Int_t pprimaryVertexFilter;
	Int_t phfCoincFilter2Th4;


	//TFile * finEff = new TFile(Form("CheckSystNuno/%s/EffFine_%d_%d.root",FolderName.Data(),CentMin,CentMax));

	TFile * finEff = new TFile(Form("NewEff2DMaps/EffFine_%d_%d.root",CentMin,CentMax));


	finEff->cd();


	TH2D * EffBptBy = (TH2D *) finEff->Get("invEff2D");
	TH2D * EffBptByInv = (TH2D *) finEff->Get("invEff2D");
	TH2D * EffBptByInvErr = (TH2D *) finEff->Get("invEff2D");


	TH2D * hEff2DInv2Shots = (TH2D *) finEff->Get("invEff2D");



	TH2D * EffBptByInvBDTWeighted = (TH2D *) finEff->Get("invEff2D");





	TFile * finTnP =  new TFile(Form("/export/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/2018Ana/BsRAA2015RunII/FromGJ/TnP/TNP2D_Bs_Cent%d-%d.root",CentMin,CentMax));

	finTnP->cd();

	TH2D * tnp_scale = (TH2D *) finTnP->Get("tnp_scale");
	TH2D * tnp_total_d = (TH2D *) finTnP->Get("tnp_total_d");
	TH2D * tnp_total_u = (TH2D *) finTnP->Get("tnp_total_u");

	//finTnP->Close();


	//	TFile * finSystWeight =  new TFile("/export/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/2018Ana/BsRAA2015RunII/FromGJ/weights.root");

	TFile * finSystWeight =  new TFile("/export/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/2018Ana/BsRAA2015RunII/FromGJ/weights_5bins.root");


	finSystWeight->cd();

	TH1D * weights_BDT_pt_5_10 = (TH1D *) finSystWeight->Get("weights_BDT_pt_5_10");
	TH1D * weights_BDT_pt_10_15 = (TH1D *) finSystWeight->Get("weights_BDT_pt_10_15");
	TH1D * weights_BDT_pt_15_20 = (TH1D *) finSystWeight->Get("weights_BDT_pt_15_20");
	TH1D * weights_BDT_pt_20_50 = (TH1D *) finSystWeight->Get("weights_BDT_pt_20_50");



	Int_t   Gsize;

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


	Float_t Bmass[NCand];


	Float_t Bmumumass[NCand];
	Float_t Bmu1eta[NCand];
	Float_t Bmu1pt[NCand];
	Float_t Bmu2eta[NCand];
	Float_t Bmu2pt[NCand];

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


	Double_t BDT_pt_5_7[NCand];
	Double_t BDT_pt_7_10[NCand];
	Double_t BDT_pt_10_15[NCand];
	Double_t BDT_pt_15_20[NCand];
	Double_t BDT_pt_20_30[NCand];
	Double_t BDT_pt_30_40[NCand];
	Double_t BDT_pt_40_50[NCand];
	Double_t BDT_pt_50_60[NCand];


	Bool_t Bmu1SoftMuID[NCand];
	Bool_t Bmu2SoftMuID[NCand];
	Bool_t Bmu1isAcc[NCand];
	Bool_t Bmu2isAcc[NCand];
	Bool_t Bmu1isTriggered[NCand];
	Bool_t Bmu2isTriggered[NCand];
	Float_t Bdtheta[NCand];




	ntphi->SetBranchAddress("Bmu1SoftMuID",Bmu1SoftMuID);
	ntphi->SetBranchAddress("Bmu2SoftMuID",Bmu2SoftMuID);


	ntphi->SetBranchAddress("Bmu1isAcc",Bmu1isAcc);
	ntphi->SetBranchAddress("Bmu2isAcc",Bmu2isAcc);


	ntphi->SetBranchAddress("Bmu1isTriggered",Bmu1isTriggered);
	ntphi->SetBranchAddress("Bmu2isTriggered",Bmu2isTriggered);

	ntphi->SetBranchAddress("Bdtheta",Bdtheta);





	Int_t GpdgId[NCand];
	Int_t GisSignal[NCand];
	Int_t GcollisionId[NCand];
	Float_t Gy[NCand];
	Float_t Gpt[NCand];

	Float_t Bmu1phi[NCand];
	Float_t Bmu2phi[NCand];

	ntHlt->SetBranchAddress("HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1",&HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1);


	ntHi->SetBranchAddress("hiBin",&hiBin);



	ntphi->SetBranchAddress("RunNo",&run);
	ntphi->SetBranchAddress("EvtNo",&evt);
	ntphi->SetBranchAddress("LumiNo",&lumi);



	ntSkim->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);
	ntSkim->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
	ntSkim->SetBranchAddress("phfCoincFilter2Th4",&phfCoincFilter2Th4);


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


	ntphi->SetBranchAddress("Bmumumass",Bmumumass);
	ntphi->SetBranchAddress("Bmu1eta",Bmu1eta);
	ntphi->SetBranchAddress("Bmu2eta",Bmu2eta);
	ntphi->SetBranchAddress("Bmu1pt",Bmu1pt);
	ntphi->SetBranchAddress("Bmu2pt",Bmu2pt);

	ntphi->SetBranchAddress("Bmu1phi",Bmu1phi);
	ntphi->SetBranchAddress("Bmu2phi",Bmu2phi);

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



	ntphi->SetBranchAddress("Btrk1nStripLayer",Btrk1nStripLayer);	
	ntphi->SetBranchAddress("Btrk2nStripLayer",Btrk2nStripLayer);
	ntphi->SetBranchAddress("Btrk1nPixelLayer",Btrk1nPixelLayer);
	ntphi->SetBranchAddress("Btrk2nPixelLayer",Btrk2nPixelLayer);

	BDT->SetBranchAddress("BDT_5_7",BDT_pt_5_7);
	BDT->SetBranchAddress("BDT_7_10",BDT_pt_7_10);
	BDT->SetBranchAddress("BDT_10_15",BDT_pt_10_15);
	BDT->SetBranchAddress("BDT_15_20",BDT_pt_15_20);
	BDT->SetBranchAddress("BDT_20_30",BDT_pt_20_30);
	BDT->SetBranchAddress("BDT_30_40",BDT_pt_30_40);
	BDT->SetBranchAddress("BDT_40_50",BDT_pt_40_50);
	BDT->SetBranchAddress("BDT_50_60",BDT_pt_50_60);

	TFile * fout = new TFile(Form("CheckSystNuno/%d-%d/%s/EffInfo.root",CentMin,CentMax,outfilept.Data()),"RECREATE");

	fout->cd();

	TTree* EffInfoTree = new TTree("EffInfoTree","EffInfoTree");	


	Int_t BsizeNew;
	Int_t runNew;
	Int_t lumiNew;
	Int_t evtNew;
	Int_t hiBinNew;

	Float_t BmassNew[NCand];
	Float_t BptNew[NCand];
	Float_t ByNew[NCand];
	Float_t BEff[NCand];
	Float_t BEffInv[NCand];
	Float_t BEffInvErr[NCand];

	Float_t BEffInvUp[NCand];
	Float_t BEffInvErrUp[NCand];
	Float_t BEffInvDown[NCand];
	Float_t BEffInvErrDown[NCand];



	Float_t BEffInvBDTWeighted[NCand];
	Float_t BEffInvErrBDTWeighted[NCand];


	Float_t BEff1D[NCand];
	Float_t BEffInv1D[NCand];
	Float_t BEffInvErr1D[NCand];

	Float_t BEffInvFit[NCand];
	Float_t BEffInvErrFit[NCand];

	Float_t BEffInvCent[NCand];
	Float_t BEffInvErrCent[NCand];



	Float_t TnPScale[NCand];
	Float_t TnPErrUp[NCand];
	Float_t TnPErrDown[NCand];


	


	int NumCand = 0;
	int BSizeCount;
	int iPass;

	int XLoc;
	int YLoc;

	int XLocTnP;
	int YLocTnP;
	int XLoc1D;
	int XLocCent;
	EffInfoTree->Branch("BsizeNew",&BsizeNew,"BsizeNew/I");
	EffInfoTree->Branch("runNew",&runNew,"runNew/I");
	EffInfoTree->Branch("evtNew",&evtNew,"evtNew/I");
	EffInfoTree->Branch("lumiNew",&lumiNew,"lumiNew/I");
	EffInfoTree->Branch("hiBinNew",&hiBinNew,"hiBinNew/I");


	EffInfoTree->Branch("BmassNew",BmassNew,"BmassNew/F");
	EffInfoTree->Branch("ByNew",ByNew,"ByNew/F");
	EffInfoTree->Branch("BptNew",BptNew,"BptNew/F");
	EffInfoTree->Branch("BEff",BEff,"BEff/F");
	EffInfoTree->Branch("BEffInv",BEffInv,"BEffInv/F");
	EffInfoTree->Branch("BEffInvErr",BEffInvErr,"BEffInvErr/F");


	EffInfoTree->Branch("BEffInvBDTWeighted",BEffInvBDTWeighted,"BEffInvBDTWeighted/F");
	EffInfoTree->Branch("BEffInvErrBDTWeighted",BEffInvErrBDTWeighted,"BEffInvErrBDTWeighted/F");



	EffInfoTree->Branch("BEff1D",BEff1D,"BEff1D/F");
	EffInfoTree->Branch("BEffInv1D",BEffInv1D,"BEffInv1D/F");
	EffInfoTree->Branch("BEffInvErr1D",BEffInvErr1D,"BEffInvErr1D/F");


	EffInfoTree->Branch("BEffInvFit",BEffInvFit,"BEffInvFit/F");
	EffInfoTree->Branch("BEffInvErrFit",BEffInvErrFit,"BEffInvErrFit/F");


	EffInfoTree->Branch("BEffInvCent",BEffInvCent,"BEffInvCent/F");
	EffInfoTree->Branch("BEffInvErrCent",BEffInvErrCent,"BEffInvErrCent/F");


	EffInfoTree->Branch("TnPScale",TnPScale,"TnPScale/F");
	EffInfoTree->Branch("TnPErrUp",TnPErrUp,"TnPErrUp/F");
	EffInfoTree->Branch("TnPErrDown",TnPErrDown,"TnPErrDown/F");

	EffInfoTree->Branch("BEffInvUp",BEffInvUp,"BEffInvUp/F");
	EffInfoTree->Branch("BEffInvErrUp",BEffInvErrUp,"BEffInvErrUp/F");

	EffInfoTree->Branch("BEffInvDown",BEffInvDown,"BEffInvDown/F");
	EffInfoTree->Branch("BEffInvErrDown",BEffInvErrDown,"BEffInvErrDown/F");


	Bool_t Bmu1isTriggeredNew[NCand];
	Bool_t Bmu2isTriggeredNew[NCand];

	EffInfoTree->Branch("Bmu1isTriggeredNew",Bmu1isTriggeredNew,"Bmu1isTriggeredNew/O");
	EffInfoTree->Branch("Bmu2isTriggeredNew",Bmu2isTriggeredNew,"Bmu2isTriggeredNew/O");


	Float_t Bmu1ptNew[NCand];
	Float_t Bmu2ptNew[NCand];
	Float_t Bmu1etaNew[NCand];
	Float_t Bmu2etaNew[NCand];
	Float_t Bmu1phiNew[NCand];
	Float_t Bmu2phiNew[NCand];

	EffInfoTree->Branch("Bmu1ptNew",Bmu1ptNew,"Bmu1ptNew/F");
	EffInfoTree->Branch("Bmu2ptNew",Bmu2ptNew,"Bmu2ptNew/F");
	EffInfoTree->Branch("Bmu1etaNew",Bmu1etaNew,"Bmu1etaNew/F");
	EffInfoTree->Branch("Bmu2etaNew",Bmu2etaNew,"Bmu2etaNew/F");
	EffInfoTree->Branch("Bmu1phiNew",Bmu1phiNew,"Bmu1phiNew/F");
	EffInfoTree->Branch("Bmu2phiNew",Bmu2phiNew,"Bmu2phiNew/F");


	int NEvents = ntphi->GetEntries();

	//int NEvents = 10000;

	int LocWeightX;


	int LocWeightY;

	double BDTSystemWeight;

	//NEvents = 100;


	//Copy a new tree//

	cout << "Copying Fit Tree for Framework" << endl;  

	TTree* EffInfoTree_NewFit;
	EffInfoTree_NewFit = EffInfoTree->CloneTree(0);
	EffInfoTree_NewFit->SetObject("EffInfoTreeFit","EffInfoTreeFit");


	//Done Copying//


	for(int i = 0; i < NEvents; i++){


		if(i%1000==0) std::cout<<std::setiosflags(std::ios::left)<<"  [ \033[1;36m"<<std::setw(10)<<i<<"\033[0m"<<" / "<<std::setw(10)<<ntphi->GetEntries()<<" ] "<<"\033[1;36m"<<Form("%.0f",100.*i/ntphi->GetEntries())<<"%\033[0m"<<"\r"<<std::flush;

		ntHlt->GetEntry(i);
		ntphi->GetEntry(i);
		ntSkim->GetEntry(i);
		ntHi->GetEntry(i);
		BDT->GetEntry(i);


		//	cout << "Bsize = " << Bsize << endl;






		BSizeCount = 0;
		iPass = 0;

		runNew = run;
		lumiNew = lumi;
		evtNew = evt;


		double lowstBDTCut = 0.32;







		if(CentMin == 30) lowstBDTCut = 0.32;
		if(CentMin == 0) lowstBDTCut = 0.32;


		for(int j =0; j < Bsize; j++){


			if((HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 == 1) && (Bpt[j] > BptLow && Bpt[j] < BptHigh) && (pprimaryVertexFilter == 1 && phfCoincFilter2Th4 ==1 && pclusterCompatibilityFilter ==1 && Btrk1Pt[j]>0.9 && Bpt[j]>5.0 && (BsvpvDistance[j]/BsvpvDisErr[j])>2.0 && Bchi2cl[j]>0.05 && TMath::Abs(Btrk1Eta[j])<2.4 && TMath::Abs(By[j])<2.4 && TMath::Abs(PVz)<15 && Bmass[j]>5 && Bmass[j]<6 && TMath::Abs(Bmumumass[j]-3.096900)<0.15 && Bmu1SoftMuID[j]==1 && Bmu2SoftMuID[j]==1 && Bmu1isAcc[j]==1 && Bmu2isAcc[j]==1 && Bmu1isTriggered[j]==1 && Bmu2isTriggered[j]==1 && (Btrk1PixelHit[j]+Btrk1StripHit[j])>=11 && (Btrk1Chi2ndf[j]/(Btrk1nStripLayer[j]+Btrk1nPixelLayer[j]))<0.18 && TMath::Abs(Btrk1PtErr[j]/Btrk1Pt[j])<0.1 && ((Bpt[j]>5 && Bpt[j]<7 && (BsvpvDistance[j]/BsvpvDisErr[j])>12 && cos(Bdtheta[j])>0.95) || (Bpt[j]>7 && Bpt[j]<10 && (BsvpvDistance[j]/BsvpvDisErr[j])>9.0 && cos(Bdtheta[j])>0.92) || (Bpt[j]>10 && Bpt[j]<60)) && ((Bpt[j]>5 && Bpt[j]<7 && BDT_pt_5_7[j]>0.02) || (Bpt[j]>7 && Bpt[j]<10 && BDT_pt_7_10[j]>0.03) || (Bpt[j]>10 && Bpt[j]<15 && BDT_pt_10_15[j]>0.09) || (Bpt[j]>15 && Bpt[j]<20 && BDT_pt_15_20[j]>0.07) || (Bpt[j]>20 && Bpt[j]<30 && BDT_pt_20_30[j]>0.10) || (Bpt[j]>30 && Bpt[j]<40 && BDT_pt_30_40[j]>0.16) || (Bpt[j]>40 && Bpt[j]<50 && BDT_pt_40_50[j]>0.20) || (Bpt[j]>50 && Bpt[j]<60 && BDT_pt_50_60[j]>0.27))) && hiBin >= CentMin * 2&& hiBin <= CentMax * 2 ){

				BmassNew[iPass] =  Bmass[j];
				BptNew[iPass] =  Bpt[j];
				ByNew[iPass] =  TMath::Abs(By[j]);

				XLoc = EffBptBy->GetXaxis()->FindBin(BptNew[iPass]);
				YLoc = EffBptBy->GetYaxis()->FindBin(ByNew[iPass]);


				XLocTnP = tnp_scale->GetXaxis()->FindBin(BptNew[iPass]);
				YLocTnP = tnp_scale->GetYaxis()->FindBin(ByNew[iPass]);
				TnPScale[iPass] = tnp_scale->GetBinContent(XLocTnP,YLocTnP);
				TnPErrUp[iPass] = tnp_total_u->GetBinContent(XLocTnP,YLocTnP);
				TnPErrDown[iPass] = tnp_total_d->GetBinContent(XLocTnP,YLocTnP);

				Bmu1isTriggeredNew[iPass] = Bmu1isTriggered[j];
				Bmu2isTriggeredNew[iPass] = Bmu2isTriggered[j];
				Bmu1ptNew[iPass] = Bmu1pt[j];
				Bmu2ptNew[iPass] = Bmu2pt[j];
				Bmu1etaNew[iPass] = Bmu1eta[j];
				Bmu2etaNew[iPass] = Bmu2eta[j];
				Bmu1phiNew[iPass] = Bmu1phi[j];
				Bmu2phiNew[iPass] = Bmu2phi[j];




				if(Shape > 0 && TwoShot == 0){

					BEff[iPass] =  EffBptBy->GetBinContent(XLoc,YLoc);
					BEffInv[iPass] =  EffBptByInv->GetBinContent(XLoc,YLoc) / TnPScale[iPass];
					BEffInvErr[iPass] =  EffBptByInvErr->GetBinContent(XLoc,YLoc) / TnPScale[iPass];

					BEffInvUp[iPass] =  EffBptByInv->GetBinContent(XLoc,YLoc) / (TnPScale[iPass] *(1 - TnPErrDown [iPass]));
					BEffInvErrUp[iPass] =  EffBptByInvErr->GetBinContent(XLoc,YLoc) / (TnPScale[iPass] *(1 - TnPErrDown [iPass]));
					BEffInvDown[iPass] =  EffBptByInv->GetBinContent(XLoc,YLoc) / (TnPScale[iPass] *(1 + TnPErrUp [iPass]));
					BEffInvErrDown[iPass] =  EffBptByInvErr->GetBinContent(XLoc,YLoc) / (TnPScale[iPass] *(1 + TnPErrUp [iPass]));


				}


				if(Shape == 0 && TwoShot == 0){

					BEff[iPass] =  EffBptBy->GetBinContent(XLoc,YLoc);
					BEffInv[iPass] =  EffBptByInv->GetBinContent(XLoc,YLoc);
					BEffInvErr[iPass] =  EffBptByInvErr->GetBinContent(XLoc,YLoc);

					//cout << "BEffInv[iPass]  = " << BEffInv[iPass]  << endl;

				}




				if(Shape > 0 && TwoShot == 1){

					BEff[iPass] =  hEff2DInv2Shots->GetBinContent(XLoc,YLoc);
					BEffInv[iPass] =  hEff2DInv2Shots->GetBinContent(XLoc,YLoc) / TnPScale[iPass];
					BEffInvErr[iPass] =  hEff2DInv2Shots->GetBinError(XLoc,YLoc) / TnPScale[iPass];

					BEffInvUp[iPass] =  hEff2DInv2Shots->GetBinContent(XLoc,YLoc) / (TnPScale[iPass] *(1 - TnPErrDown [iPass]));
					BEffInvErrUp[iPass] =  hEff2DInv2Shots->GetBinError(XLoc,YLoc) / (TnPScale[iPass] *(1 - TnPErrDown [iPass]));
					BEffInvDown[iPass] =  hEff2DInv2Shots->GetBinContent(XLoc,YLoc) / (TnPScale[iPass] *(1 + TnPErrUp [iPass]));
					BEffInvErrDown[iPass] =  hEff2DInv2Shots->GetBinError(XLoc,YLoc) / (TnPScale[iPass] *(1 + TnPErrUp [iPass]));

				}


				if(Shape == 0 && TwoShot == 1){

					BEff[iPass] =  hEff2DInv2Shots->GetBinContent(XLoc,YLoc);
					BEffInv[iPass] =  hEff2DInv2Shots->GetBinContent(XLoc,YLoc);
					BEffInvErr[iPass] =  hEff2DInv2Shots->GetBinError(XLoc,YLoc);
				}



				LocWeightX = EffBptByInvBDTWeighted->GetXaxis()->FindBin(BptNew[iPass]);
				LocWeightY = EffBptByInvBDTWeighted->GetYaxis()->FindBin(ByNew[iPass]);


				BEffInvBDTWeighted[iPass] =  EffBptByInvBDTWeighted->GetBinContent(LocWeightX,LocWeightY)/ TnPScale[iPass]; 
				BEffInvErrBDTWeighted[iPass] =  EffBptByInvBDTWeighted->GetBinError(LocWeightX,LocWeightY)/ TnPScale[iPass];









				BEffInvFit[iPass] = func->Eval(BptNew[iPass]);
				BEffInvErrFit[iPass] = 0;




				NumCand = NumCand + 1;

				BSizeCount = BSizeCount + 1;
				iPass = iPass + 1;

			}



		}




		BsizeNew = BSizeCount;
		hiBinNew = hiBin;



		if(BsizeNew > 0) EffInfoTree->Fill();

		if(BsizeNew > 0) EffInfoTree_NewFit->Fill();


	}


	fout->Write();
	fout->Close();




}
