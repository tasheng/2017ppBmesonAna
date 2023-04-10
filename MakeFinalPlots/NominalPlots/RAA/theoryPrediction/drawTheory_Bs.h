#include <iostream>
#include <fstream>
#include "TGraph.h"
Color_t colorTAMU_Bs = kOrange+8;
Style_t styleTAMU_Bs = 1001;
Color_t colorCUJET_Bs = kYellow+2;
Style_t styleCUJET_Bs = 3344;
using namespace std;
void plotTheory_Bs(double xThreshold = 10)
{
	const int n = 1000;
	Int_t nbin=0;
	Float_t temp;
	Float_t aCx[n],aCy[n],aCyl[n],aCyh[n];
	Float_t bCx[n],bCy[n],bCyl[n],bCyh[n];
	// TAMU
	ifstream getdata_tamu("theoryPrediction/theorypre_Bs/TAMU_20180509.txt");
	if(!getdata_tamu.is_open()) {
		cout<<"Opening the file fails: TAMU"<<endl;
		return;
	}
	nbin=0;
	while(!getdata_tamu.eof())
	{
		getdata_tamu>>aCx[nbin]>>aCy[nbin];
		nbin++;
	}
	Float_t* aTAMUB5TeVx = new Float_t[nbin];
	Float_t* aTAMUB5TeVxe = new Float_t[nbin];
	Float_t* aTAMUB5TeVy = new Float_t[nbin];
	Float_t* aTAMUB5TeVye = new Float_t[nbin];
	for(int i=0;i<nbin;i++)
	{
		aTAMUB5TeVx[i] = aCx[i];
		aTAMUB5TeVxe[i] = 0;
		aTAMUB5TeVy[i] = aCy[i];
		aTAMUB5TeVye[i] = 0.;
	}
	TGraphErrors* gTAMUB5TeV = new TGraphErrors(nbin, aTAMUB5TeVx, aTAMUB5TeVy, aTAMUB5TeVxe, aTAMUB5TeVye);
	gTAMUB5TeV->SetName("gTAMUB5TeV");
	gTAMUB5TeV->SetLineColor(colorTAMU_Bs);
	gTAMUB5TeV->SetFillColor(colorTAMU_Bs);
	gTAMUB5TeV->SetFillStyle(styleTAMU_Bs);
	gTAMUB5TeV->SetLineWidth(4);
	//gTAMUB5TeV->Draw("3 same");
	gTAMUB5TeV->Draw("L same");

	// CUJET
	ifstream getdata_cujet("theoryPrediction/theorypre_Bs/CUJET_20180510.dat");
	if(!getdata_cujet.is_open()) {
		cout<<"Opening the file fails: CUJET"<<endl;
		return;
	}
	nbin=0;
	while(!getdata_cujet.eof())
	{
		getdata_cujet>>aCx[nbin]>>aCyl[nbin]>>aCyh[nbin];
		nbin++;
	}
	Float_t* aCUJETB5TeVx = new Float_t[nbin];
	Float_t* aCUJETB5TeVxe = new Float_t[nbin];
	Float_t* aCUJETB5TeVy = new Float_t[nbin];
	Float_t* aCUJETB5TeVye = new Float_t[nbin];
	for(int i=0;i<nbin;i++)
	{
		aCUJETB5TeVx[i] = aCx[i];
		aCUJETB5TeVxe[i] = 0;
		aCUJETB5TeVy[i] = (aCyh[i]+aCyl[i])/2;
		aCUJETB5TeVye[i] = (aCyh[i]-aCyl[i])/2;
	}
	TGraphErrors* gCUJETB5TeV = new TGraphErrors(nbin, aCUJETB5TeVx, aCUJETB5TeVy, aCUJETB5TeVxe, aCUJETB5TeVye);
	gCUJETB5TeV->SetName("gCUJETB5TeV");
	gCUJETB5TeV->SetLineWidth(1);
	gStyle->SetHatchesLineWidth(3);
//	gCUJETB5TeV->SetLineColor(kRed-4);
	gCUJETB5TeV->SetFillColor(colorCUJET_Bs);
	gCUJETB5TeV->SetFillColorAlpha(colorCUJET_Bs, 0.5);
	gCUJETB5TeV->SetFillStyle(styleCUJET_Bs);
	gCUJETB5TeV->Draw("3 same");

  // remove points for x < threshold
  auto removePointsBelowX = []<typename T>(T* graph, double xThreshold) {
    for (int i=0; i<graph->GetN(); ++i) {
      double x = graph->GetX()[i];
      if (x < xThreshold) {
        graph->RemovePoint(i);
        i--; // need to decrement the counter to account for the removed point
      }
    }
  };
  removePointsBelowX(gCUJETB5TeV, xThreshold);
  removePointsBelowX(gTAMUB5TeV, xThreshold);
}
