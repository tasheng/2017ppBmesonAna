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

void PlotNPShapesNew(){
	


	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	TString infile;

	const int NBins = 9;

	int ptU[NBins] = {2,3,5,7,10,15,20,30,50};
	int ptL[NBins] = {1,2,3,5,7,10,15,20,30};
/*	
	int ptU[NBins] = {100};
	int ptL[NBins] = {0};	
*/
	for(int i = 0; i < NBins; i++){

	infile = Form("OutFile/BsNPStudies_%d_%d.root",ptL[i],ptU[i]);

	TFile *fin = new TFile(infile.Data());

	TH1D * BmassBs = (TH1D * ) fin->Get("BmassBs");
	TH1D * BmassBs_nosig = (TH1D * ) fin->Get("BmassBs_nosig");
	TH1D * BmassBsPiK = (TH1D * ) fin->Get("BmassBsPiK");
	TH1D * BmassBsKK = (TH1D * ) fin->Get("BmassBsKK");
	TH1D * BmassBsKKfake_binfo = (TH1D * ) fin->Get("BmassBsKKfake_binfo");
	TH1D * BmassBsXpipi = (TH1D * ) fin->Get("BmassBsXpipi");
	TH1D * BmassBsKKX = (TH1D * ) fin->Get("BmassBsKKX");
	TH1D * BmassB0KStar = (TH1D * ) fin->Get("BmassB0KStar");

	TCanvas *c = new TCanvas("c","c",600,600);
	c->cd();


	BmassBs->Sumw2();
	BmassBs_nosig->Sumw2();
	BmassBsPiK->Sumw2();
	BmassBsKK->Sumw2();
	BmassBsKKfake_binfo->Sumw2();
	BmassBsXpipi->Sumw2();
	BmassBsKKX->Sumw2();
	BmassB0KStar->Sumw2();
	/*
	BmassBs->Scale(1.0/BmassBs->Integral());
	BmassBs_nosig->Scale(1.0/BmassBs_nosig->Integral());
	BmassBsPiK->Scale(1.0/BmassBsPiK->Integral());
	BmassBsKK->Scale(1.0/BmassBsKK->Integral());
	BmassBsKKfake_binfo->Scale(1.0/BmassBsKKfake_binfo->Integral());
	BmassBsXpipi->Scale(1.0/BmassBsXpipi->Integral());
	BmassBsKKX->Scale(1.0/BmassBsKKX->Integral());
	*/

	BmassBs->SetMarkerStyle(20);
	BmassBs->SetMarkerColor(2);
	BmassBs->SetMarkerSize(1);
	BmassBs->SetLineColor(1);


	BmassBs_nosig->SetMarkerStyle(20);
	BmassBs_nosig->SetMarkerColor(3);
	BmassBs_nosig->SetMarkerSize(1);
	BmassBs_nosig->SetLineColor(1);


	BmassBsPiK->SetMarkerStyle(20);
	BmassBsPiK->SetMarkerColor(4);
	BmassBsPiK->SetMarkerSize(1);
	BmassBsPiK->SetLineColor(1);


	BmassBsKK->SetMarkerStyle(20);
	BmassBsKK->SetMarkerColor(1);
	BmassBsKK->SetMarkerSize(1);
	BmassBsKK->SetLineColor(1);


	BmassBsKKfake_binfo->SetMarkerStyle(20);
	BmassBsKKfake_binfo->SetMarkerColor(5);
	BmassBsKKfake_binfo->SetMarkerSize(1);
	BmassBsKKfake_binfo->SetLineColor(1);


	BmassBsXpipi->SetMarkerStyle(20);
	BmassBsXpipi->SetMarkerColor(6);
	BmassBsXpipi->SetMarkerSize(1);
	BmassBsXpipi->SetLineColor(1);


	BmassBsKKX->SetMarkerStyle(20);
	BmassBsKKX->SetMarkerColor(7);
	BmassBsKKX->SetMarkerSize(1);
	BmassBsKKX->SetLineColor(1);

	BmassB0KStar->SetMarkerStyle(20);
	BmassB0KStar->SetMarkerColor(29);
	BmassB0KStar->SetMarkerSize(1);
	BmassB0KStar->SetLineColor(1);




	BmassBs->GetYaxis()->SetTitle("Counts");
	BmassBs->GetYaxis()->SetTitleOffset(1.3);
	BmassBs->GetXaxis()->SetTitle("Bmass (#mu #mu K^{+} K^{-}) (GeV/c^{2})");
	BmassBs->GetXaxis()->CenterTitle();
	BmassBs->GetYaxis()->CenterTitle();

	BmassBs->SetTitle(Form("B_{s} Invariant Mass at %d < B p_{T} < %d",ptL[i],ptU[i]));




	BmassBs->Draw("ep");
	BmassBs_nosig->Draw("epSAME");
	BmassBsPiK->Draw("epSAME");
	BmassBsKK->Draw("epSAME");
	BmassBsKKfake_binfo->Draw("epSAME");
	BmassBsXpipi->Draw("epSAME");
	BmassBsKKX->Draw("epSAME");
	BmassB0KStar->Draw("epSAME");



	TLegend* leg = new TLegend(0.10,0.43,0.30,0.85,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(BmassBs,"Bs Signal","pl");
	leg->AddEntry(BmassBs_nosig,"Bs Inclusive Background","pl");
	leg->AddEntry(BmassBsPiK,"Bs -> J/#Psi  #pi K ","pl");
	leg->AddEntry(BmassBsKK,"Bs -> J/#Psi  K K (Not from #Phi meson resonance)","pl");
	leg->AddEntry(BmassBsKKfake_binfo,"B+ ->  J/#Psi  K+ K- (Random K-)","pl");
	leg->AddEntry(BmassBsXpipi,"X-> J/#Psi #pi #pi","pl");
	leg->AddEntry(BmassBsKKX,"Bs -> J/#Psi X K K","pl");
	leg->AddEntry(BmassB0KStar,"B^{0} -> J/#Psi K*","pl");

	leg->Draw("same");


	c->SaveAs(Form("NPNewPlots/NPBackground_%d.png",i));

	BmassBs_nosig->GetYaxis()->SetTitle("Counts");
	BmassBs_nosig->GetYaxis()->SetTitleOffset(1.3);
	BmassBs_nosig->GetXaxis()->SetTitle("Bmass (#mu #mu K^{+} K^{-}) (GeV/c^{2})");
	BmassBs_nosig->GetXaxis()->CenterTitle();
	BmassBs->GetYaxis()->CenterTitle();

	BmassBs_nosig->SetTitle(Form("B_{s} Invariant Mass at %d < B p_{T} < %d",ptL[i],ptU[i]));



	BmassBs_nosig->Draw("ep");
	BmassBsPiK->Draw("epSAME");
	BmassBsKK->Draw("epSAME");
	BmassBsKKfake_binfo->Draw("epSAME");
	BmassBsXpipi->Draw("epSAME");
	BmassBsKKX->Draw("epSAME");
	BmassB0KStar->Draw("epSAME");

	TH1D * BmassBsSub = (TH1D *) BmassBs_nosig->Clone("BmassBsSub");
	BmassBsSub->Add(BmassBsKK,-1);
	BmassBsSub->Add(BmassB0KStar,-1);
	BmassBsSub->Draw("epSAME");

	BmassBsSub->SetMarkerStyle(20);
	BmassBsSub->SetMarkerColor(2);
	BmassBsSub->SetMarkerSize(1);
	BmassBsSub->SetLineColor(1);



	TLegend* leg2 = new TLegend(0.28,0.53,0.40,0.85,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.04);
	leg2->SetTextFont(42);
	leg2->SetFillStyle(0);
	leg2->AddEntry(BmassBs_nosig,"Inclusive Background","pl");
	leg2->AddEntry(BmassBsPiK,"B_{s} -> J/#Psi  #pi K ","pl");
	leg2->AddEntry(BmassBsKK,"B_{s} -> J/#Psi  K K (No #Phi meson)","pl");
	leg2->AddEntry(BmassBsKKfake_binfo,"B^{+} ->  J/#Psi  K+ K- (Add Another K-)","pl");
	leg2->AddEntry(BmassBsXpipi,"X(3872) -> J/#Psi #pi #pi","pl");
	leg2->AddEntry(BmassBsKKX,"B_{s} -> J/#Psi h K K","pl");
	leg2->AddEntry(BmassB0KStar,"B^{0} -> J/#Psi K*","pl");
	leg2->AddEntry(BmassBsSub,"Inclusive Background - (KK) - (K*)","pl");


	leg2->Draw("same");

	//Evaluate NP to Sig Ratio//
	float BsMassPDG = 5.36688;
	float Width = 0.08;
	float SigLow = BsMassPDG - Width;
	float SigHigh = BsMassPDG + Width;
	
	int BinLow = BmassBs_nosig->GetXaxis()->FindBin(SigLow);
	int BinHigh = BmassBs_nosig->GetXaxis()->FindBin(SigHigh);

	float NPBkgd = BmassBs_nosig->Integral(BinLow,BinHigh);
	float Sig = BmassBs->Integral(BinLow,BinHigh);
	float Ratio = NPBkgd/Sig;

	cout << "NPBkgd = " <<  NPBkgd << "  Sig = " << Sig << "    Ratio = " << Ratio << endl;


	c->SaveAs(Form("NPNewPlots/NPBackgroundONLY_%d.png",i));

/*
	//BsBPB0//


	TH1D * BmassB0Inclsive = (TH1D * ) fin->Get("BmassB0Inclsive");
	TH1D * BmassBsInclsive = (TH1D * ) fin->Get("BmassBsInclsive");
	TH1D * BmassBPInclsive = (TH1D * ) fin->Get("BmassBPInclsive");
	TH1D * BmassOtherInclsive = (TH1D * ) fin->Get("BmassOtherInclsive");

	BmassB0Inclsive->Sumw2();
	BmassBsInclsive->Sumw2();
	BmassBPInclsive->Sumw2();
	BmassOtherInclsive->Sumw2();


	BmassB0Inclsive->SetMarkerStyle(20);
	BmassB0Inclsive->SetMarkerColor(1);
	BmassB0Inclsive->SetMarkerSize(1);
	BmassB0Inclsive->SetLineColor(1);

	BmassBsInclsive->SetMarkerStyle(20);
	BmassBsInclsive->SetMarkerColor(2);
	BmassBsInclsive->SetMarkerSize(1);
	BmassBsInclsive->SetLineColor(2);

	BmassBPInclsive->SetMarkerStyle(20);
	BmassBPInclsive->SetMarkerColor(3);
	BmassBPInclsive->SetMarkerSize(1);
	BmassBPInclsive->SetLineColor(3);

	BmassOtherInclsive->SetMarkerStyle(20);
	BmassOtherInclsive->SetMarkerColor(4);
	BmassOtherInclsive->SetMarkerSize(1);
	BmassOtherInclsive->SetLineColor(4);


	BmassBsInclsive->GetYaxis()->SetTitle("Counts");
	BmassBsInclsive->GetYaxis()->SetTitleOffset(1.3);
	BmassBsInclsive->GetXaxis()->SetTitle("Bmass (#mu #mu K^{+} K^{-}) (GeV/c^{2})");
	BmassBsInclsive->GetXaxis()->CenterTitle();
	BmassBsInclsive->GetYaxis()->CenterTitle();

	BmassBsInclsive->SetTitle(Form("B_{s} Invariant Mass at %d < B p_{T} < %d",ptL[i],ptU[i]));

	BmassBsInclsive->Draw("ep");
	BmassB0Inclsive->Draw("epSAME");
	BmassBPInclsive->Draw("epSAME");
	BmassOtherInclsive->Draw("epSAME");


	TLegend* leg3 = new TLegend(0.28,0.50,0.40,0.80,NULL,"brNDC");
	leg3->SetBorderSize(0);
	leg3->SetTextSize(0.04);
	leg3->SetTextFont(42);
	leg3->SetFillStyle(0);
	leg3->AddEntry(BmassBsInclsive,"B_{s} Inclusive Decay","pl");
	leg3->AddEntry(BmassB0Inclsive,"B^{0} Inclusive Decay","pl");
	leg3->AddEntry(BmassBPInclsive,"B^{+} Inclusive Decay","pl");
	leg3->AddEntry(BmassOtherInclsive,"Everything Else","pl");

	leg3->Draw("same");

	c->SaveAs(Form("BmesonBackground_%d.png",i));

	//BsBPB0//

	TH1D * B0KStarKPiMass = (TH1D * ) fin->Get("B0KStarKPiMass");
	TH1D * B0KStarPiKMass = (TH1D * ) fin->Get("B0KStarPiKMass");
	
	TH1D * B0KStarKPiMassInc = (TH1D * ) B0KStarKPiMass->Clone("B0KStarKPiMassInc");

	B0KStarKPiMassInc->Add(B0KStarPiKMass);
	B0KStarKPiMassInc->GetYaxis()->SetTitle("Counts");
	B0KStarKPiMassInc->GetXaxis()->SetTitle("#pi K Mass (GeV/c^{2})");

	B0KStarKPiMassInc->SetMarkerStyle(20);
	B0KStarKPiMassInc->SetMarkerColor(1);
	B0KStarKPiMassInc->SetMarkerSize(1);
	B0KStarKPiMassInc->SetLineColor(1);

	B0KStarKPiMassInc->GetYaxis()->SetTitleOffset(1.3);
	B0KStarKPiMassInc->GetXaxis()->CenterTitle();
	B0KStarKPiMassInc->GetYaxis()->CenterTitle();

	B0KStarKPiMassInc->Draw("ep");
	
	c->SaveAs(Form("B0KPiMass_%d.png",i));
*/
	}




}
