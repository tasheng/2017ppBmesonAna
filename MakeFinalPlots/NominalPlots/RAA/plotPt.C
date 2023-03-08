/*
   Macro to plot xsec and ratio vs pt for Bs and Bp)

Input: txt files in inputDir, with 10 columns:
ptmin, ptmax, central val, statUp, statDown, systUp, systDown, glbUp, glbDown, abscissae

Output: xsec vs pt, ratio vs pt.

*/
#include <Rtypes.h>
#include <TVirtualPad.h>
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <Riostream.h>

#include <map>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"

#include "TPaveStats.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "CMS_lumi.C"
#include "tdrstyle.C"
#include "DrawLHCb.C"

#include "auxiliaryPt.h"
#include "auxiliaryRef.h"

#include "theoryPrediction/drawTheory.h"
#include "theoryPrediction/uti.h"

//#include "outsideSource/lhcb.C"
#endif
using namespace std;


constexpr double xHigh = 60;
double yLow = 0.0;
double yHigh = 1.55;

const map<int, string> pname = {
    {0, "BP"},
    {1, "Bs"}
    };


bool drawLight = 1;
void adjustLegend(TLegend* l);


void plotPt(bool bSavePlots       = 1,
		bool bDoDebug         = 0, //  figure out if things are read properly
		bool whichPlot        = 1, //0 is x-sec, 1 is for ratio
		bool drawRef          = 1, //draw Ref (for ratio only
    bool drawBs = 1,
		const char* inputDir  = "dataSource", //inptu txt files
		const char* outputDir = "figs")// where the output figures will be
		{
      bool drawThm = (!drawBs);
      bool drawlhcb         = 1;
		gSystem->mkdir(Form("./%s/png",outputDir), kTRUE);
		gSystem->mkdir(Form("./%s/pdf",outputDir), kTRUE);

		//   gROOT->ProcessLine(".x lhcb.C");
		//  lhcb();
		//set the style
		setTDRStyle();

    if (drawBs) {
      yLow = 0.3;
      yHigh = 1.7;
      drawLight = 0;
    }


		//samples:
		const unsigned int nMes      = 2;
		const char* inputFileType[2] = {"RAA_pt", "RAA_pt"};
		const char* mesonName[nMes]  = {"Bs", "Bp"};

		Int_t endMes = nMes;

    std::vector<int> nlines = {4, 4};
		for (Int_t ib=0; ib<endMes; ib++){
			ifstream in;
			string inputFileName = Form("%s/%s_%s_New.txt",inputDir,inputFileType[whichPlot],mesonName[ib]);
			//if(whichPlot==1) inputFileName = Form("%s/%s.txt",inputDir,inputFileType[whichPlot]);
	//		if(whichPlot==1)  inputFileName = Form("%s/%s_%s.txt",inputDir,inputFileType[whichPlot],mesonName[ib]);		
			cout << "########## Input file name: " << inputFileName << endl;

			in.open(inputFileName.c_str());
			if (!in.is_open()) {
				cout << "input file " << inputFileName << " cannot be open" << endl;
				continue;
			}
			double x[20];

			//get first line, the header, and discard it
			string tmpstrg;
			getline(in,tmpstrg);//ignore first line/ the header

      cout << "start reading" << endl;
			unsigned nEntry=0;
			// while(in >> x[0] >> x[1] >> x[2] >> x[3] >> x[4] >> x[5] >> x[6] >> x[7] >> x[8] >> x[9])
			// {
			for (auto iline = 0; iline < nlines[ib]; ++iline)
			{
        in >> x[0] >> x[1] >> x[2] >> x[3] >> x[4] >> x[5] >> x[6] >> x[7] >> x[8] >> x[9];
				glbSystDown = x[8]*100;
				glbSystUp   = x[7]*100;

        cout << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ", " << x[4] << ", " << x[5] << ", " << x[6] << ", " << x[7] << ", " << x[8] << ", " << x[9] << endl;

					if(ib==0){//bs
						if(nEntry==0){
						binLow[nEntry]    = x[9];							
						glbSystUpBs  = x[8]*100;
						glbSystDownBs   = x[7]*100;
						bs_low[nEntry] = x[2]; //central value
						//stat uncert
						bs_low_yStatL[nEntry] = x[4]*x[2];
						bs_low_yStatH[nEntry] =  x[3]*x[2];
						//bin width
				

						bs_low_xErrL[nEntry] = x[9]-x[0];
						bs_low_xErrH[nEntry] = x[1]-x[9];
						//systm. uncert
						bs_low_ySystL[nEntry] = x[6]*x[2];
						bs_low_ySystH[nEntry] = x[5]*x[2];
						cout << "nEntry = " << nEntry << "   x[1] =  " << x[1] << "   x[9] = " << x[9] << "  bs_low_xErrH[0] =  " << bs_low_xErrH[0] << endl;

						if(bDoDebug){
							cout << "Central  \t statUncertL \t systUncerL " << endl;
							cout<< bs_low[nEntry] << "\t" << bs_low_yStatL[nEntry] << "\t" <<   bs_low_ySystL[nEntry] << endl;
							cout<<"==========================================="<<endl;
						}
						}
					}
					else{//bp
						//central value	
						if(nEntry < nBinsLowNew){
						binLowNew[nEntry]    = x[9];	
						glbSystUpBp = x[8]*100;
						glbSystDownBp   = x[7]*100;

						bpl_low[nEntry] = x[2];
						//stat uncert
						bpl_low_yStatL[nEntry] = x[4]*x[2];
						bpl_low_yStatH[nEntry] =  x[3]*x[2];
						//bin width
						bpl_low_xErrL[nEntry] = x[9]-x[0];
						bpl_low_xErrH[nEntry] = x[1]-x[9];
						//systm. uncert
						bpl_low_ySystL[nEntry] = x[6]*x[2];
						bpl_low_ySystH[nEntry] = x[5]*x[2];

						if(bDoDebug){
							cout << "Central  \t statUncertL  \t systUncertL " << endl;
							cout<< bpl_low[nEntry] << "\t" << bpl_low_yStatL[nEntry] << "\t" <<   bpl_low_ySystL[nEntry] << endl;
							cout<<"==========================================="<<endl;
						}
					}
					
				}
				// high-pt bins, |y|<2.4
					if(ib==0){//bs
						if(nEntry>0){
						//central value
						binHigh[nEntry-1] = x[9];										
						bs_high[nEntry-1] = x[2];
						//stat uncert
						bs_high_yStatL[nEntry-1] = x[4]*x[2];
						bs_high_yStatH[nEntry-1] = x[3]*x[2];
						//bin width
						bs_high_xErrL[nEntry-1] = x[9]-x[0];
						bs_high_xErrH[nEntry-1] = x[1]-x[9];
						//systm. uncert
						bs_high_ySystL[nEntry-1] = x[6]*x[2];
						bs_high_ySystH[nEntry-1] = x[5]*x[2];
						//	cout << "Bs Bro: " << " bs_high_ySystL = " << bs_high_ySystL[nEntry-1] << "  bs_high_ySystH = " << bs_high_ySystH[nEntry-1] << endl;

						// if(bDoDebug){
						// 	if(nEntry==1) cout << "Central  \t statUncertL  \t systUncertL "<< endl;
						// 	cout<< bs_high[nEntry-2] << "\t" << bs_high_yStatL[nEntry-2] << "\t" << bs_high_ySystL[nEntry-2] << endl;
						// }
						}
					}else{//bp
						if(nEntry >= nBinsLowNew){
              auto iEntryHigh = nEntry - nBinsLowNew;
						//central value
						binHighNew[iEntryHigh] = x[9];				
						bpl_high[iEntryHigh] = x[2];
						//stat uncert
						bpl_high_yStatL[iEntryHigh] = x[4]*x[2];
						bpl_high_yStatH[iEntryHigh] = x[3]*x[2];
						//bin width
						bpl_high_xErrL[iEntryHigh] = x[9]-x[0];
						bpl_high_xErrH[iEntryHigh] = x[1]-x[9];
						//systm. uncert
						bpl_high_ySystL[iEntryHigh] = x[6]*x[2];
						bpl_high_ySystH[iEntryHigh] = x[5]*x[2];

            cout << "ysyst:" << bpl_high_ySystL[iEntryHigh] << ", ratio: " <<
              bpl_high_ySystL[iEntryHigh] / bpl_high[iEntryHigh] << "\n";


						}
					
				}//high-pt bins
				nEntry++;
			}//reading input file line by line

			cout<<"@@@@@@@@@@@ Finished meson "<<mesonName[ib]<<endl;
			if(bDoDebug){
				if(ib==0) cout<<"Element_low = "<< binLow[0] <<"\t bs_low[0]= " << bs_low[0] << "\t statUncertL = "<<bs_low_yStatL[0]<<endl;

				if(ib==1)
				{
					for(int i=0; i<3; i++){
						cout<<"Element_high " << i << "\t binHigh[i] = "<< binHigh[i] <<"\t bs_high[i]= " << bs_high[i] << "\t statUncertL = "<<bs_high_yStatL[i]<<"\t systUncertL = "<<bs_high_ySystL[i]<< endl;
					}
				}
			}//bDebug
			in.close();//close input file

		}//for each meson,ib
//		bs_low_xErrH[0] = 1.0;
//		bs_low_xErrL[1] = 1.73;
		
//		bs_high_xErrH[4] = 5.0;


		//----------------------------------------------------------------
		// gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
		// Bs
		TGraphAsymmErrors *pgBs_low = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_yStatL,bs_low_yStatH);
		TGraphAsymmErrors *pgBs_high= new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high,
				bs_high_xErrL, bs_high_xErrH,
				bs_high_yStatL,bs_high_yStatH);
		TGraphAsymmErrors *pgBs_lowWhite = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_yStatL,bs_low_yStatH);   //AddWhite//
		

		// Bplus
		TGraphAsymmErrors *pgBpl_low = new TGraphAsymmErrors(nBinsLowNew, binLowNew, bpl_low,
				bpl_low_xErrL, bpl_low_xErrH,
				bpl_low_yStatL,bpl_low_yStatH);
		
		TGraphAsymmErrors *pgBpl_lowWhite = new TGraphAsymmErrors(nBinsLowNew, binLowNew, bpl_low,
				bpl_low_xErrL, bpl_low_xErrH,
				bpl_low_yStatL,bpl_low_yStatH);
		
		TGraphAsymmErrors *pgBpl_high= new TGraphAsymmErrors(nBinsHighNew,binHighNew,bpl_high,
				bpl_high_xErrL, bpl_high_xErrH,
				bpl_high_yStatL,bpl_high_yStatH);

		for(int i = 0; i < 5; i++){

			cout << "binLowHigh = " << binHighNew[i] << "   bs_high_xErrL = " << bs_high_xErrL[i] << "   bs_high_xErrH = " << bs_high_xErrH[i] << endl;
		}



		//Systmeatic uncertainty
		// Bs
		TGraphAsymmErrors *pgBs_syst_low    = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_ySystL,bs_low_ySystH);
		TGraphAsymmErrors *pgBs_syst_high   = new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high,
				bs_high_xErrL, bs_high_xErrH,
				bs_high_ySystL,bs_high_ySystH);


		// Bp
		TGraphAsymmErrors *pgBpl_syst_low    = new TGraphAsymmErrors(nBinsLowNew, binLowNew,  bpl_low,
				bpl_low_xErrL, bpl_low_xErrH,
				bpl_low_ySystL, bpl_low_ySystH);
		TGraphAsymmErrors *pgBpl_syst_high   = new TGraphAsymmErrors(nBinsHighNew,binHighNew, bpl_high,
				bpl_high_xErrL,bpl_high_xErrH,
				bpl_high_ySystL,bpl_high_ySystH);
		// pgBs_syst_high->Draw("APL");
		//==========================================
		//------------------------------------------
		TGraphAsymmErrors *pgRatio_low    = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_yStatL,bs_low_yStatH);

		TGraphAsymmErrors *pgRatio_lowWhite    = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_yStatL,bs_low_yStatH);

		TGraphAsymmErrors *pgRatio_high   = new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high,
				bs_high_xErrL, bs_high_xErrH,
				bs_high_yStatL,bs_high_yStatH);

		TGraphAsymmErrors *pgRatio_syst_low = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_ySystL,bs_low_ySystH);

		TGraphAsymmErrors *pgRatio_syst_lowWhite = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
				bs_low_xErrL, bs_low_xErrH,
				bs_low_ySystL,bs_low_ySystH);
		
		TGraphAsymmErrors *pgRatio_syst_high= new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high,
				bs_high_xErrL, bs_high_xErrH,
				bs_high_ySystL,bs_high_ySystH);

		//================== reference
		TGraphAsymmErrors * FragBand = new TGraphAsymmErrors(BandBin,BandX,BandY,BandXErr,BandXErr,BandYErr,BandYErr);
		FragBand->SetName("BandErr");
		FragBand->SetMarkerStyle(20);
		FragBand->SetMarkerSize(0.8);
		FragBand->SetFillColorAlpha(kGreen,0.5);
		FragBand->SetFillStyle(3004);
		FragBand->SetLineWidth(2);
		FragBand->SetLineColor(kGreen);

		TGraphAsymmErrors *LHCb7TeVRef = new TGraphAsymmErrors(LHCb7TeVNPoints, LHCb7TeVX, LHCb7TeVY,
				LHCb7TeVXErrDown, LHCb7TeVXErrUp,
				LHCb7TeVYErr,LHCb7TeVYErr);

		LHCb7TeVRef->SetName("LHCb7TeVRef");




		//==================
		// // **************** marker setup
		// Bs
		// marker style
		pgBs_low->SetMarkerStyle(markerLow[0]);
		pgBs_lowWhite->SetMarkerStyle(markerHigh[0]);
		pgBs_high->SetMarkerStyle(markerHigh[0]);

		pgBpl_lowWhite->SetMarkerStyle(markerHigh[1]);
		pgBpl_low->SetMarkerStyle(markerLow[1]);
		pgBpl_high->SetMarkerStyle(markerHigh[1]);
			

		pgRatio_low->SetMarkerStyle(markerRatio[0]);
		pgRatio_lowWhite->SetMarkerStyle(markerRatio[1]);

		pgRatio_high->SetMarkerStyle(markerRatio[1]);

		// marker size
		pgBs_low->SetMarkerSize(markerSizeLow[0]);
		pgBs_lowWhite->SetMarkerSize(markerSizeLow[0]);
		
		pgBs_high->SetMarkerSize(markerSizeHigh[0]);

		pgBpl_lowWhite->SetMarkerSize(markerSizeLow[1]*0.86);
		pgBpl_low->SetMarkerSize(markerSizeLow[1]);
		pgBpl_high->SetMarkerSize(markerSizeHigh[1]);

		pgRatio_low->SetMarkerSize(markerSizeRatio[0]);
		pgRatio_lowWhite->SetMarkerSize(markerSizeRatio[0]);
		
		pgRatio_high->SetMarkerSize(markerSizeRatio[1]);

		// marker color
		pgBs_low->SetMarkerColor(colorLow[0]);
		pgBs_high->SetMarkerColor(colorHigh[0]);
	
		pgBs_lowWhite->SetMarkerColor(kWhite);
		pgBpl_lowWhite->SetMarkerColor(kWhite);
	
		pgBs_low->SetMarkerColor(kBlue+2);
		pgBs_high->SetMarkerColor(kBlue+2);

		//cout << "Color = " << TColor::GetColor("#C32148") << endl; 

//		pgBs_low->SetMarkerColor(924);
		pgBpl_low->SetMarkerColor(colorLow[1]);
		pgBpl_high->SetMarkerColor(colorHigh[1]);
    if (!drawBs) {
      pgBpl_low->SetMarkerColor(kAzure-1);
      pgBpl_high->SetMarkerColor(kAzure-1);
    } else {
      pgBpl_low->SetMarkerColor(kGreen+3);
      pgBpl_high->SetMarkerColor(kGreen+3);
    }

		pgRatio_lowWhite->SetMarkerColor(kWhite);
		pgRatio_low->SetMarkerColor(kRed + 2);
		pgRatio_high->SetMarkerColor(kRed + 2);

		// line color
		pgBs_low->SetLineColor(colorLow[0]);
		pgBs_high->SetLineColor(colorHigh[0]);
		pgBs_low->SetLineColor(kBlue+2);	
		pgBs_high->SetLineColor(kBlue+2);

		pgBpl_low->SetLineColor(colorLow[1]);
		pgBpl_high->SetLineColor(colorHigh[1]);
		pgBpl_low->SetLineColor(kGreen+3);
		pgBpl_high->SetLineColor(kGreen+3);


		pgRatio_lowWhite->SetLineColor(kRed + 2);
		pgRatio_low->SetLineColor(kRed + 2);
		pgRatio_high->SetLineColor(kRed + 2);


		// systematic boxes
		pgBs_syst_low->SetFillColorAlpha(kBlue-9,0.5);
		pgBs_syst_high->SetFillColorAlpha(kBlue-9,0.5);
    if (!drawBs) {
      pgBpl_syst_low->SetFillColorAlpha(kAzure+7,0.5);
      pgBpl_syst_high->SetFillColorAlpha(kAzure+7,0.5);
    } else {
      pgBpl_syst_low->SetFillColorAlpha(kGreen-9,0.5);
      pgBpl_syst_high->SetFillColorAlpha(kGreen-9,0.5);
    }

		pgRatio_syst_low->SetFillColorAlpha(colorRatio[0],0.2);
		pgRatio_lowWhite->SetFillColor(kWhite);

		pgRatio_syst_high->SetFillColorAlpha(colorRatio[1],0.2);

		//Reference 
		LHCb7TeVRef->SetMarkerColor(kBlue);
		LHCb7TeVRef->SetLineColor(kBlue);
		LHCb7TeVRef->SetMarkerStyle(20);
		LHCb7TeVRef->SetMarkerSize(markerSizeLow[1]);


		//================
		//=============================================
		// read TAMU
		std::ifstream TAMUBsBP("outsideSource/TAMUPT.dat");
		const int NBinsTAMU = 150;
		double TAMUBsBPPt[NBinsTAMU];
		double TAMUBsBPPtErr[NBinsTAMU];
		double TAMUBsBPRatio[NBinsTAMU];
		double TAMUBsBPRatioErr[NBinsTAMU];
		for(int i = 0; i < NBinsTAMU; i++){
			TAMUBsBP >> TAMUBsBPPt[i] >> TAMUBsBPRatio[i];
		}
		for(int i = 0; i < NBinsTAMU; i++){
			TAMUBsBPPtErr[i] = 0.1;
			TAMUBsBPRatioErr[i] = 0.005;
			//     cout  << "TAMUBsBPPt[i] = "  <<  TAMUBsBPPt[i] << "     TAMUBsBPRatio[i] = " <<  TAMUBsBPRatio[i] << endl;
		}
		TGraphErrors* TAMUTheory = new TGraphErrors(NBinsTAMU,TAMUBsBPPt,TAMUBsBPRatio,TAMUBsBPPtErr,TAMUBsBPRatioErr);
		TAMUTheory->SetName("TAMUTheory");
		/*
		   TAMUTheory->SetMarkerStyle(20);
		   TAMUTheory->SetMarkerSize(0.0);
		   TAMUTheory->SetFillColor(kOrange+3);
		   TAMUTheory->SetFillStyle(3002);
		   TAMUTheory->SetLineColor(kOrange);
		   TAMUTheory->SetLineWidth(3);
		   */

		TAMUTheory->SetLineWidth(3);
//		TAMUTheory->SetMarkerColor(kOrange);
//		TAMUTheory->SetLineColor(kOrange);

		TAMUTheory->SetMarkerColor(kOrange + 1);   //Stronger TAMU Color
		TAMUTheory->SetLineColor(kOrange + 1);
		//=============================================

		// read CAO
		std::ifstream CAOBsBP("outsideSource/CAOPT.dat");
		const int NBinsCAO = 40;
		double CAOBsBPPt[NBinsCAO];
		double CAOBsBPPtErr[NBinsCAO];
		double CAOBsBPRatio[NBinsCAO];
		double CAOBsBPRatioErr[NBinsCAO];
		for(int i = 0; i < NBinsCAO; i++){
			CAOBsBP >> CAOBsBPPt[i] >> CAOBsBPRatio[i];
		}
		for(int i = 0; i < NBinsCAO; i++){
			CAOBsBPPtErr[i] = 0.1;
			CAOBsBPRatioErr[i] = 0.005;
			//     cout  << "CAOBsBPPt[i] = "  <<  CAOBsBPPt[i] << "     CAOBsBPRatio[i] = " <<  CAOBsBPRatio[i] << endl;
		}
		TGraphErrors* CAOTheory = new TGraphErrors(NBinsCAO,CAOBsBPPt,CAOBsBPRatio,CAOBsBPPtErr,CAOBsBPRatioErr);
		CAOTheory->SetName("CAOTheory");
		CAOTheory->SetMarkerStyle(20);
		CAOTheory->SetLineWidth(3);
		CAOTheory->SetMarkerColor(kGreen+1);
		CAOTheory->SetLineColor(kGreen+1);



		//-------------------------------------------
		TF1 *f4 = new TF1("f4","1",5, xHigh);
		f4->SetLineWidth(1);
		f4->SetLineColor(1);
		f4->SetLineStyle(1);
		f4->GetYaxis()->SetTitle(yAxName[whichPlot]);
		f4->GetXaxis()->SetTitle(xAxName[0]);
		f4->GetXaxis()->CenterTitle(kTRUE);
		f4->GetYaxis()->CenterTitle();
		if(whichPlot==1){
			f4->GetYaxis()->SetTitleSize(0.06*0.83);
			f4->GetXaxis()->SetTitleSize(0.06*0.83);
		
			f4->GetYaxis()->SetTitleOffset(1.40);
			f4->GetXaxis()->SetTitleOffset(1.20);


		}

		if(whichPlot==0){
			f4->GetYaxis()->SetTitleSize(0.06*0.80);
			f4->GetYaxis()->SetTitleOffset(1.5);
			f4->GetXaxis()->SetTitleOffset(1.05);
			f4->GetXaxis()->SetTitleSize(f4->GetXaxis()->GetTitleSize() * 0.77);
			f4->GetXaxis()->SetTitleOffset(1.18);
		
			f4->GetXaxis()->SetTitleSize(0.06*0.83);  //Unify Textsize
			f4->GetYaxis()->SetTitleSize(0.06*0.83);
			f4->GetYaxis()->SetTitleOffset(1.40);
			f4->GetXaxis()->SetTitleOffset(1.20);
		

			cout << "Offset = " << f4->GetYaxis()->GetTitleOffset() << endl;


		}

		f4->GetYaxis()->SetRangeUser(1e2,2e7);
		//if(whichPlot==1) f4->GetYaxis()->SetRangeUser(0.0,1.8);
		if(whichPlot==1) f4->GetYaxis()->SetRangeUser(0.0,4.0);
		if(whichPlot==0) f4->GetYaxis()->SetRangeUser(100.0,2000000); //1.8 -> 0.9
		
		f4->GetYaxis()->SetRangeUser(yLow, yHigh);
		f4->SetLineColor(0);
		//f4->GetXaxis()->SetNdivisions(-6);

		//---------------- general stuff
		TLatex *lat = new TLatex();
		lat->SetNDC();

		// // ##################################################### x-sec canvas
		TCanvas *pc1 = new TCanvas("pc1","pc1");
		f4->SetLineStyle(2);
		f4->Draw();// axis
		//pc1->SetBottomMargin(0.10);
	
	//	cout << "pc1->GetBottomMargin() = " << pc1->GetBottomMargin() << endl;
	
		pc1->SetBottomMargin(0.135); //Enlarge Margin

		CMS_lumi(pc1,19011,0);

		//Added to fix//
		pgBs_syst_high->SetLineColor(kBlue-9);
		pgBs_syst_low->SetLineColor(kBlue-9);
		pgBpl_syst_high->SetLineColor(kGreen-9);
		pgBpl_syst_low->SetLineColor(kGreen-9);

    // D0 comparison
    TFile fd0("dataSource/HEPData-ins1616207-v1-Table_1.root");
    auto d0Raa = fd0.Get<TH1F>("Table 1/Hist1D_y1");
    auto d0RaaStat = fd0.Get<TH1F>("Table 1/Hist1D_y1_e1");
    auto d0RaaSyst = fd0.Get<TH1F>("Table 1/Hist1D_y1_e2");
    auto d0RaaGraph = fd0.Get<TGraphAsymmErrors>("Table 1/Graph1D_y1");

    TFile fch("dataSource/HEPData-ins1496050-v2-Table_15.root");
    auto chRaa = fch.Get<TH1F>("Table 15/Hist1D_y1");
    auto chRaaStat = fch.Get<TH1F>("Table 15/Hist1D_y1_e1");
    auto chRaaSyst = fch.Get<TH1F>("Table 15/Hist1D_y1_e2");
    auto chRaaGraph = fch.Get<TGraphAsymmErrors>("Table 15/Graph1D_y1");

    auto gd0 = (TGraphAsymmErrors*) pgBpl_high->Clone();
    auto gd0Syst = (TGraphAsymmErrors*) pgBpl_syst_high->Clone();
    unsigned ngd0 = d0Raa->GetNbinsX();
    gd0->Set(ngd0);
    gd0Syst->Set(ngd0);
    for (auto i = 0; i < ngd0; ++i) {
      gd0->SetPoint(i, d0RaaGraph->GetPointX(i), d0RaaGraph->GetPointY(i));
      gd0Syst->SetPoint(i, d0RaaGraph->GetPointX(i), d0RaaGraph->GetPointY(i));
      gd0->SetPointError(i, d0RaaGraph->GetErrorX(i), d0RaaGraph->GetErrorX(i),
                         d0RaaStat->GetBinContent(i + 1), d0RaaStat->GetBinContent(i + 1));
      gd0Syst->SetPointError(i, d0RaaGraph->GetErrorX(i), d0RaaGraph->GetErrorX(i),
                             d0RaaSyst->GetBinContent(i + 1), d0RaaSyst->GetBinContent(i + 1));
    }

		gd0Syst->SetFillColorAlpha(kRed-9,0.5);
		gd0Syst->SetLineColor(kRed-9);
		gd0->SetLineColor(kRed+2);
		gd0->SetMarkerColor(kRed+2);

    auto gch = (TGraphAsymmErrors*) pgBpl_high->Clone();
    auto gchSyst = (TGraphAsymmErrors*) pgBpl_syst_high->Clone();
    unsigned ngch = chRaa->GetNbinsX();
    gch->Set(ngch);
    gchSyst->Set(ngch);
    for (auto i = 0; i < ngch; ++i) {
      gch->SetPoint(i, chRaaGraph->GetPointX(i), chRaaGraph->GetPointY(i));
      gchSyst->SetPoint(i, chRaaGraph->GetPointX(i), chRaaGraph->GetPointY(i));
      gch->SetPointError(i, chRaaGraph->GetErrorX(i), chRaaGraph->GetErrorX(i),
                         chRaaStat->GetBinContent(i + 1), chRaaStat->GetBinContent(i + 1));
      gchSyst->SetPointError(i, chRaaGraph->GetErrorX(i), chRaaGraph->GetErrorX(i),
                             chRaaSyst->GetBinContent(i + 1), chRaaSyst->GetBinContent(i + 1));
    }

		gchSyst->SetFillColorAlpha(kYellow-9,0.5);
		gchSyst->SetLineColor(kYellow-9);
		gch->SetLineColor(kYellow+2);
		gch->SetMarkerColor(kYellow+2);

		if(whichPlot==0)// x-section
		{
      TLegend *legendThm=new TLegend(0.18,0.64,0.60,0.915,"");
      TLegend *legendAds=new TLegend(0.71,0.64,0.95,0.78,"");
      if(drawThm){
        adjustLegend(legendThm);
        adjustLegend(legendAds);
        plotTheory();
        TGraphAsymmErrors* gThmDummy1 = new TGraphAsymmErrors();
        TGraphAsymmErrors* gThmDummy2 = new TGraphAsymmErrors();
        TGraphAsymmErrors* gThmDummy3 = new TGraphAsymmErrors();
        TGraphAsymmErrors* gThmDummy4 = new TGraphAsymmErrors();
        TGraphAsymmErrors* gThmDummy5 = new TGraphAsymmErrors();
        gThmDummy1->SetLineColor(kOrange+8);
        gThmDummy2->SetLineColor(kRed-4);
        gThmDummy3->SetLineColor(0);
        gThmDummy4->SetLineColor(0);
        gThmDummy5->SetLineColor(0);
        gThmDummy3->SetFillColorAlpha(kYellow+2,0.5);
        gThmDummy4->SetFillColorAlpha(kViolet-8,0.5);
        gThmDummy5->SetFillColorAlpha(kGreen-2,0.5);
        gThmDummy1->SetLineWidth(8.);
        gThmDummy2->SetLineWidth(8.);
        gThmDummy2->SetLineStyle(6);
        gThmDummy3->SetFillStyle(3344);
        gThmDummy4->SetFillStyle(3352);
        gThmDummy5->SetFillStyle(3325);
        TLegendEntry *ent_thm1 = legendThm->AddEntry(gThmDummy1,"TAMU","l");
        TLegendEntry *ent_thm2 = legendThm->AddEntry(gThmDummy2,"Djordjevic","l");
        TLegendEntry *ent_thm3 = legendThm->AddEntry(gThmDummy3,"CUJET3.0","f");
        TLegendEntry *ent_thm4 = legendThm->AddEntry(gThmDummy4,"AdS/CFT HH D(p)","f");
        //TLegendEntry *ent_thm5 = legendThm->AddEntry(gThmDummy5,"AdS/CFT HH D = const","f");
        TLegendEntry *ent_thm5 = legendThm->AddEntry(gThmDummy5,"AdS/CFT HH D=const","f");
        //legendAds->SetHeader("AdS/CFT HH");
        //TLegendEntry *ent_thm4 = legendAds->AddEntry(gThmDummy4,"D(p)","f");
        //TLegendEntry *ent_thm5 = legendAds->AddEntry(gThmDummy5,"D = const","f");
        ent_thm1->SetTextSize(0.038);
        ent_thm2->SetTextSize(0.038);
        ent_thm3->SetTextSize(0.038);
        ent_thm4->SetTextSize(0.038);
        ent_thm5->SetTextSize(0.038);
      }
      if(drawThm) {
        legendThm->Draw();
        //legendAds->Draw();
      }
			//gPad->SetLogy();

      // D0 comparison
      if (drawLight) {
        gd0Syst->Draw("5");
        gd0->Draw("P");
        gchSyst->Draw("5");
        gch->Draw("P");
      }


      if (drawBs) {
        pgBs_syst_high->SetLineWidth(1);

        pgBs_syst_high->Draw("5same");
        pgBs_syst_low->Draw("5same");
        //	pgBs_low->SetMarkerStyle(24);
        pgBs_lowWhite->Draw("P");
        pgBs_low->Draw("P");
        pgBs_high->Draw("P");
        pgBs_lowWhite->Draw("P");
      }

      pgBpl_syst_low->Draw("5same");
      pgBpl_syst_high->Draw("5same");

			
			pgBpl_lowWhite->Draw("P");
			pgBpl_low->Draw("P");

			pgBpl_high->Draw("P");
			pgBpl_lowWhite->Draw("P");

			TLine * Unity = new TLine(5,1,xHigh,1);
			Unity->SetLineWidth(2);
			Unity->SetLineStyle(2);
			Unity->SetLineColor(1);
		
			Unity->Draw("SAME");


		}
		else{//ratio plot
			if(drawRef){
				//		FragBand->Draw("5same");
				TAMUTheory->Draw("l");
				CAOTheory->Draw("l");
				LHCb7TeVRef->Draw("P");
			}

			pgRatio_syst_low->SetLineColor(colorRatio[0]);
			pgRatio_syst_high->SetLineColor(colorRatio[1]);

			pgRatio_syst_low->Draw("5same");
			pgRatio_low->Draw("P");
			pgRatio_lowWhite->Draw("P");
			
			pgRatio_syst_high->Draw("5same");
			pgRatio_high->Draw("P");


		}

		//supplemental info on plot:
		if(whichPlot==0){

			lat->SetTextFont(42);
			lat->SetTextSize(ltxSetTextSize2 * 1.3);
			lat->SetTextSize(ltxSetTextSize4); //Enlarge Labels

			cout << "ltxSetTextSize4 = " << ltxSetTextSize4 << endl;
			lat->DrawLatex(xsec_ltxText1_xStart + 0.30,xsec_ltxText1_yStart - 0.30+ 0.037,Form("global uncertainty: #pm %.1f%%",glbSystUpBp));

      if (drawThm) {
        double ShiftX = 0.05;
        double ShiftY = 0.13;
        double ysep = 0.035;
        double yoffset = 0.005;

        TLatex bold;
        bold.SetNDC();
        bold.SetTextAngle(0);
        bold.SetTextColor(kBlack);
        bold.SetTextFont(42);
        bold.SetTextAlign(31);
        bold.SetTextSize(ltxSetTextSize2);
        bold.SetTextFont(42);
        bold.SetTextSize(ltxSetTextSize2);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextAngle(0);
        latex.SetTextColor(kBlack);

        latex.SetTextFont(42);
        latex.SetTextAlign(11);
        latex.SetTextSize(ltxSetTextSize2);

        TLegend *legb = new TLegend(legXsec_xLowStart + ShiftX + 0.02,
                                    legXsec_y + 0.13 + yoffset * 4 - ysep * 0,
                                    legXsec_xLowEnd+ShiftX-0.10+0.06,
                                    legXsec_y + 0.13 - yoffset - ysep * 2,
                                    "                    ","brNDC");
        legb->SetBorderSize(0);
        legb->SetTextSize(ltxSetTextSize2);
        legb->SetTextFont(42);
        legb->AddEntry(pgBpl_low,"1.5 < |y| < 2.4","p");
        legb->AddEntry(pgBpl_high,"|y| < 2.4 ","p");
        legb->Draw();

        TLegend *legd = new TLegend(legXsec_xLowStart + ShiftX + 0.02,
                                    legXsec_y + 0.13 + yoffset * 4 - ysep * 3,
                                    legXsec_xLowEnd+ShiftX-0.10+0.06,
                                    legXsec_y + 0.13 - yoffset - ysep * 4,
                                    "                    ","brNDC");
        legd->SetBorderSize(0);
        legd->SetTextSize(ltxSetTextSize2);
        legd->SetTextFont(42);
        legd->AddEntry(gd0, "|y| < 1 ","p");
        legd->Draw();

        TLegend *legh = new TLegend(legXsec_xLowStart + ShiftX + 0.02,
                                    legXsec_y + 0.13 + yoffset * 4 - ysep * 5,
                                    legXsec_xLowEnd+ShiftX-0.10+0.06,
                                    legXsec_y + 0.13 - yoffset - ysep * 6,
                                    "                    ","brNDC");
        legh->SetBorderSize(0);
        legh->SetTextSize(ltxSetTextSize2);
        legh->SetTextFont(42);
        legh->AddEntry(gch, "|#eta| < 1 ","p");
        legh->Draw();

        latex.DrawLatex(legXsec_xLowStart+ShiftX+0.02 + 0.009,
                        legXsec_y + 0.13, "#bf{B^{+}}, Cent. 0-90%");
        latex.DrawLatex(legXsec_xLowStart+ShiftX+0.02 + 0.009,
                        legXsec_y + 0.13 - ysep * 3, "#bf{D^{0}}, Cent. 0-100%");
        latex.DrawLatex(legXsec_xLowStart+ShiftX+0.02 + 0.009,
                        legXsec_y + 0.13 - ysep * 5, "#bf{h^{+}}, Cent. 0-100%");

      } else {
        double shiftX = 0.05;
        double shiftY = -0.03;
        double ysep = 0.04;
        double yoffset = 0.005;

        TLatex latex;
        latex.SetNDC();
        latex.SetTextAngle(0);
        latex.SetTextColor(kBlack);

        latex.SetTextFont(42);
        latex.SetTextAlign(11);
        latex.SetTextSize(ltxSetTextSize2);

        TLegend *legb = new TLegend(legXsec_xLowStart + shiftX + 0.02,
                                    legXsec_y + shiftY + 0.13 + yoffset * 4 - ysep * 0,
                                    legXsec_xLowEnd+shiftX-0.10+0.06,
                                    legXsec_y + shiftY + 0.13 - yoffset - ysep * 2,
                                    "                    ","brNDC");
        legb->SetBorderSize(0);
        legb->SetTextSize(ltxSetTextSize2);
        legb->SetTextFont(42);
        legb->AddEntry(pgBpl_low,"1.5 < |y| < 2.4","p");
        legb->AddEntry(pgBpl_high,"|y| < 2.4 ","p");
        legb->Draw();

        TLegend *legbs = new TLegend(legXsec_xLowStart + shiftX + 0.02,
                                    legXsec_y + shiftY + 0.13 + yoffset * 4 - ysep * 3,
                                    legXsec_xLowEnd+shiftX-0.10+0.06,
                                    legXsec_y + shiftY + 0.13 - yoffset - ysep * 5,
                                    "                    ","brNDC");
        legbs->SetBorderSize(0);
        legbs->SetTextSize(ltxSetTextSize2);
        legbs->SetTextFont(42);
        legbs->AddEntry(pgBs_low,"1.5 < |y| < 2.4","p");
        legbs->AddEntry(pgBs_high,"|y| < 2.4 ","p");
        legbs->Draw();

        latex.DrawLatex(legXsec_xLowStart+shiftX+0.02 + 0.009,
                        legXsec_y + shiftY + 0.13, "#bf{B^{+}}, Cent. 0-90%");
        latex.DrawLatex(legXsec_xLowStart+shiftX+0.02 + 0.009,
                        legXsec_y + shiftY + 0.13 - ysep * 3, "#bf{B^{0}_{s}}, Cent. 0-90%");
      }


		}else{

			lat->SetTextSize(ltxSetTextSize2*1.7);
			lat->SetTextSize(ltxSetTextSize2*1.7/1.3 *1.17 * ltxSetTextSize4/ltxSetTextSize2); //Enlarge Labels
			lat->SetTextSize(ltxSetTextSize2*1.7/1.3 *1.17 * ltxSetTextSize4/ltxSetTextSize2 /1.05); //Enlarge Labels - Reduce by 8% to match
			

			lat->SetTextFont(42);
			//lat->DrawLatex(ratio_ltxText1_xStart,ratio_ltxText1_yStart,"#bf{#frac{B_{s}^{0}}{B^{+}}}");

			//cout << "ltxSetTextSize1 = " << ltxSetTextSize1 << endl;

			lat->SetTextFont(42);
			lat->SetTextSize(ltxSetTextSize2 * 1.3);
			lat->SetTextSize(ltxSetTextSize4*1.08/1.08); //Enlarge Labels
	
			//lat->DrawLatex(ratio_ltxText2_xStart-0.04,ratio_ltxText2_yStart,"Cent. 0-90%");
	//		lat->DrawLatex(ratio_ltxText2_xStart-0.04,ratio_ltxText2_yStart,"Cent. 0-90%");
	//		lat->DrawLatex(ratio_ltxText2_xStart-0.06,ratio_ltxText2_yStart,"Centrality 0-90%"); //Expand Cent.

			lat->SetTextFont(42);
			lat->SetTextSize(ltxSetTextSize2 * 1.3);
			lat->SetTextSize(ltxSetTextSize4*1.08); //Enlarge Labels
			lat->SetTextSize(ltxSetTextSize4*1.08/1.08); //Enlarge Labels - Reduce 8%
	
		//	lat->DrawLatex(legRatio_xLowStart-0.14,legRatio_y-0.15 - 0.30,Form("global uncert.: #pm %.1f %%",glbSystDown));
		//	lat->DrawLatex(legRatio_xLowStart-0.20,legRatio_y-0.15 - 0.30,Form("global uncertainty: #pm %.1f%%",glbSystDown));

			cout << "legRatio_xLowStart = " << legRatio_xLowStart << "   legRatio_y = " << legRatio_y << endl; 

			// legend
			TLegend *legRatio = new TLegend(legRatio_xLowStart-0.04,legRatio_y+0.08,legRatio_xLowEnd-0.04,legRatio_y+0.18,NULL,"brNDC");
			legRatio->SetBorderSize(0);
			legRatio->SetTextFont(42);
			legRatio->SetTextSize(ltxSetTextSize2 * 1.3);
			legRatio->SetTextSize(ltxSetTextSize4*1.08); //Enlarge Labels
			legRatio->SetTextSize(ltxSetTextSize4*1.08/1.08); //Enlarge Labels - Reduce by 6% to match

			
			legRatio->SetLineColor(1);
			legRatio->SetLineStyle(1);
			legRatio->SetLineWidth(1);
			legRatio->SetFillColor(19);
			legRatio->SetFillStyle(0);
			TLegendEntry *entry1 = legRatio->AddEntry("pgRatio_low","1.5 < |y| < 2.4","p");
			entry1->SetTextFont(42);
			entry1->SetMarkerStyle(markerRatio[0]);
			entry1->SetMarkerSize(1.7);
			entry1->SetMarkerColor(kRed+3);
			
			entry1->SetFillStyle(1001);

			TLegendEntry *entry2 = legRatio->AddEntry("pgRatio_high","|y| < 2.4","p");
			entry2->SetTextFont(42);
			entry2->SetMarkerStyle(markerRatio[1]);
			entry2->SetMarkerColor(kRed+3);
		
			entry2->SetMarkerSize(1.7);
			entry2->SetFillStyle(1001);

			legRatio->Draw();
			//---------------

		}

		// gPad->RedrawAxis();
		pc1->Update();

		if(bSavePlots)
		{
			if (whichPlot==0)
			{
				pc1->SaveAs(Form("%s/pdf/BRAA_vsPt_%s.pdf",outputDir, pname.at(drawBs).c_str()));
				pc1->SaveAs(Form("%s/png/BRAA_vsPt_%s.png",outputDir, pname.at(drawBs).c_str()));
			}else{
				pc1->SaveAs(Form("%s/pdf/ratio_vsPt_ref%d_%d.pdf",outputDir,drawRef,drawlhcb));
				pc1->SaveAs(Form("%s/png/ratio_vsPt_ref%d_%d.png",outputDir,drawRef,drawlhcb));
			}
		}

		}

void adjustLegend(TLegend* l){
	l->SetBorderSize(0);
	l->SetLineColor(0);
	l->SetFillColor(0);
	l->SetFillStyle(1001);
	l->SetTextFont(42);
	l->SetTextSize(0.04);
}
