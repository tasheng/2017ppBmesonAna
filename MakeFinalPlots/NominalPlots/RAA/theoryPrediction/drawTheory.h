#include <iostream>
#include <fstream>
#include "TGraph.h"
const int n = 1000;
Float_t aCx[n],aCy[n],aCyl[n],aCyh[n];
Float_t bCx[n],bCy[n],bCye[n],bDummy[n];
Int_t nbin=0;
Float_t temp;
Float_t cCx[n],cCy[n],cCye[n],cDummy[n];
Float_t dCx[n],dCy[n],dCyl[n],dCyh[n];
using namespace std;
void plotTheory(bool addpostPRL = false, double xThreshold = 10)
{
  // CUJET
  ifstream countlines("theoryPrediction/theorypre/CUJET_5.02TeV_B0_RAA_0-100.dat");
  int counter=0;
  while(!countlines.eof())
    {
      countlines>>temp>>temp>>temp;
      counter++;
    }
  countlines.close();
  countlines.clear();
  ifstream getcujet("theoryPrediction/theorypre/CUJET_5.02TeV_B0_RAA_0-100.dat");

  const int size = counter;
  double pt[size],RAAhigh[size],RAAlow[size],RAA[size],RAAerr[size];
  double zero[size];
  double width[size];

  for(int i=0;i<size;i++)
    {
      getcujet>>pt[i]>>RAAhigh[i]>>RAAlow[i];
      RAA[i] = (RAAhigh[i]+RAAlow[i])/2.;
      RAAerr[i] = (RAAhigh[i]-RAAlow[i])/2.;
      width[i]=0.5;
    }
  getcujet.close();

  TGraphAsymmErrors* gCUJETB5TeV = new TGraphAsymmErrors(size,pt,RAA,width,width,RAAerr,RAAerr);
  gCUJETB5TeV->SetName("gCUJETB5TeV");
  gCUJETB5TeV->SetLineWidth(1);
  gStyle->SetHatchesLineWidth(3);
  gCUJETB5TeV->SetLineColor(kRed-4);
  //gCUJETB5TeV->SetFillColor(kRed-4); //CWRv8
  //gCUJETB5TeV->SetFillColorAlpha(kRed-4, 0.5); //CWRv8
  // gCUJETB5TeV->SetFillColor(kYellow+2);
  // gCUJETB5TeV->SetFillColorAlpha(kYellow+2, 0.5);
  gCUJETB5TeV->SetFillColorAlpha(kYellow+2, 0.4);
  // gCUJETB5TeV->SetFillStyle(3344);
  gCUJETB5TeV->SetFillStyle(1001);
  gCUJETB5TeV->Draw("3 same");
  /*
  TGraph* gCUJETB5TeV = new TGraph(nbin,aCx,aCy);
  //gCUJETB5TeV->SetLineColor(4); // PAS color
  //gCUJETB5TeV->SetMarkerColor(4); // PAS color
  gCUJETB5TeV->SetLineColor(kRed-4); //paper color 20170224
  gCUJETB5TeV->SetMarkerColor(kRed-4); //paper color 20170224
  gCUJETB5TeV->SetLineWidth(3);
  gCUJETB5TeV->SetLineStyle(2);
  gCUJETB5TeV->SetMarkerSize(0.15);
  gCUJETB5TeV->Draw("c same");
  */
  //

  // AdS/CFT Constant
  //ifstream getadscft1("theoryPrediction/theorypre/ADSCFT_20170310/RAA_5TeV_0080_B0_DiffusionConstant_band.dat");
  //ifstream getadscft1("theoryPrediction/theorypre/ADSCFT_20170318/RAA_5TeV_0080_B0_DiffusionConstant_band.dat");
  ifstream getadscft1("theoryPrediction/theorypre/ADSCFT_20170320/RAA_5TeV_0080_B0_DiffusionConstant_band.dat");
  if(!getadscft1.is_open()) cout<<"Opening the file fails: ADS/CFT DiffusionConstant"<<endl;
  nbin=0;
  while(!getadscft1.eof())
    {
      getadscft1>>bCx[nbin]>>bCy[nbin]>>bDummy[nbin]>>bCye[nbin];
      nbin++;
    }
  bCx[nbin] = bCx[nbin-1];
  bCy[nbin] = bCy[nbin-1];
  bCye[nbin] = bCye[nbin-1];
  nbin++;
  Float_t* aADSCFT1B5TeVx = new Float_t[nbin];
  Float_t* aADSCFT1B5TeVxe = new Float_t[nbin];
  Float_t* aADSCFT1B5TeVy = new Float_t[nbin];
  Float_t* aADSCFT1B5TeVye = new Float_t[nbin];
  for(int i=0;i<nbin;i++)
    {
      aADSCFT1B5TeVx[i] = bCx[i];
      aADSCFT1B5TeVxe[i] = 0;
      aADSCFT1B5TeVy[i] = bCy[i];
      aADSCFT1B5TeVye[i] = bCye[i];
    }
  TGraphErrors* gADSCFT1B5TeV = new TGraphErrors(nbin, aADSCFT1B5TeVx, aADSCFT1B5TeVy, aADSCFT1B5TeVxe, aADSCFT1B5TeVye);
  gADSCFT1B5TeV->SetName("gADSCFT1B5TeV");
  gADSCFT1B5TeV->SetLineColor(kViolet-2);
  // gADSCFT1B5TeV->SetFillColor(kViolet-2);
  // gADSCFT1B5TeV->SetFillColorAlpha(kViolet-2,0.5);
  gADSCFT1B5TeV->SetFillColorAlpha(kViolet-2,0.4);
  // gADSCFT1B5TeV->SetFillStyle(3325);
  gADSCFT1B5TeV->SetFillStyle(1001);
  gADSCFT1B5TeV->Draw("4 same");

  // AdS/CFT Dynamic
  //ifstream getadscft2("theoryPrediction/theorypre/ADSCFT_20170310/RAA_5TeV_0080_B0_DiffusionDynamic_band.dat");
  //ifstream getadscft2("theoryPrediction/theorypre/ADSCFT_20170318/RAA_5TeV_0080_B0_DiffusionDynamic_band.dat");
  ifstream getadscft2("theoryPrediction/theorypre/ADSCFT_20170320/RAA_5TeV_0080_B0_DiffusionDynamic_band.dat");
  if(!getadscft2.is_open()) cout<<"Opening the file fails: ADS/CFT DiffusionDynamic"<<endl;
  nbin=0;
  while(!getadscft2.eof())
    {
      getadscft2>>cCx[nbin]>>cCy[nbin]>>bDummy[nbin]>>cCye[nbin];
      nbin++;
    }
  cCx[nbin] = cCx[nbin-1];
  cCy[nbin] = cCy[nbin-1];
  cCye[nbin] = cCye[nbin-1];
  nbin++;
  Float_t* aADSCFT2B5TeVx = new Float_t[nbin];
  Float_t* aADSCFT2B5TeVxe = new Float_t[nbin];
  Float_t* aADSCFT2B5TeVy = new Float_t[nbin];
  Float_t* aADSCFT2B5TeVye = new Float_t[nbin];
  for(int i=0;i<nbin;i++)
    {
      aADSCFT2B5TeVx[i] = cCx[i];
      aADSCFT2B5TeVxe[i] = 0;
      aADSCFT2B5TeVy[i] = cCy[i];
      aADSCFT2B5TeVye[i] = cCye[i];
    }
  TGraphErrors* gADSCFT2B5TeV = new TGraphErrors(nbin, aADSCFT2B5TeVx, aADSCFT2B5TeVy, aADSCFT2B5TeVxe, aADSCFT2B5TeVye);
  gADSCFT2B5TeV->SetName("gADSCFT2B5TeV");
  gADSCFT2B5TeV->SetLineColor(kGreen-8);
  // gADSCFT2B5TeV->SetFillColor(kGreen-8);
  gADSCFT2B5TeV->SetFillColorAlpha(kGreen-8,0.4);
  // gADSCFT2B5TeV->SetFillStyle(3352);
  gADSCFT2B5TeV->SetFillStyle(1001);
  gADSCFT2B5TeV->Draw("3 same");

  // Magdalena
  TGraph* gMagdalenaB5TeV = new TGraph("theoryPrediction/theorypre/Magdalena.dat");
  gMagdalenaB5TeV->SetName("gMagdalenaB5TeV");
  //gMagdalenaB5TeV->SetLineColor(kGreen+4); //CWRv8
  gMagdalenaB5TeV->SetLineColor(kRed-4);
  gMagdalenaB5TeV->SetLineWidth(8);
  gMagdalenaB5TeV->SetLineStyle(6);
  //gMagdalenaB5TeV->SetFillColor(kGreen+4);
  //gMagdalenaB5TeV->SetFillStyle(3004);
  //gMagdalenaB5TeV->Draw("f same");
  gMagdalenaB5TeV->Draw("c same");

  // TAMU
  //ifstream getdata("theoryPrediction/theorypre/TAMUminbiasB.txt");
  ifstream getdata("theoryPrediction/theorypre/TAMUminbiasB_20170419.txt");
  if(!getdata.is_open()) cout<<"Opening the file fails: TAMU"<<endl;
  nbin=0;
  while(!getdata.eof())
    {
      getdata>>aCx[nbin]>>aCyh[nbin]>>aCyl[nbin];
      if(aCx[nbin] != 0 && aCyh[nbin] != 0)
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
      aTAMUB5TeVy[i] = (aCyh[i]+aCyl[i])/2.;
      aTAMUB5TeVye[i] = abs(aCyh[i]-aCyl[i])/2.;
    }
  TGraphErrors* gTAMUB5TeV = new TGraphErrors(nbin, aTAMUB5TeVx, aTAMUB5TeVy, aTAMUB5TeVxe, aTAMUB5TeVye);
  gTAMUB5TeV->SetName("gTAMUB5TeV");
  //gTAMUB5TeV->SetLineColor(kViolet-6);
  //gTAMUB5TeV->SetFillColor(kViolet-6);
  gTAMUB5TeV->SetLineColor(kViolet+1);
  gTAMUB5TeV->SetFillColor(kViolet+1);
  gTAMUB5TeV->SetFillStyle(1001);
  gTAMUB5TeV->Draw("3 same");
  //gTAMUB5TeV->Draw("4 same");
  //gTAMUB5TeV->Draw("pe same");

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
  removePointsBelowX(gADSCFT1B5TeV, xThreshold);
  removePointsBelowX(gADSCFT2B5TeV, xThreshold);
  removePointsBelowX(gTAMUB5TeV, xThreshold);
  removePointsBelowX(gMagdalenaB5TeV, xThreshold);

  if(addpostPRL){
    ifstream polmodel("theoryPrediction/theorypre/Pol_20170921.txt");
    if(!polmodel.is_open()) cout<<"Opening the file fails: Pol's Model"<<endl;
    nbin=0;
    while(!polmodel.eof())
      {
        polmodel>>dCx[nbin]>>dCyh[nbin]>>dCyl[nbin];
        nbin++;
      }
    dCx[nbin-1] = dCx[nbin-2];
    dCyl[nbin-1] = dCyl[nbin-2];
    dCyh[nbin-1] = dCyh[nbin-2];
    Float_t* aPOLMODEL1B5TeVx = new Float_t[nbin];
    Float_t* aPOLMODEL1B5TeVxe = new Float_t[nbin];
    Float_t* aPOLMODEL1B5TeVy = new Float_t[nbin];
    Float_t* aPOLMODEL1B5TeVye = new Float_t[nbin];
    for(int i=0;i<nbin;i++)
      {
        aPOLMODEL1B5TeVx[i] = dCx[i];
        aPOLMODEL1B5TeVxe[i] = 0;
        aPOLMODEL1B5TeVy[i] = fabs(dCyh[i]+dCyl[i])/2;
        aPOLMODEL1B5TeVye[i] = fabs(dCyh[i]-dCyl[i])/2;
      }
    TGraphErrors* gPOLMODEL1B5TeV = new TGraphErrors(nbin, aPOLMODEL1B5TeVx, aPOLMODEL1B5TeVy, aPOLMODEL1B5TeVxe, aPOLMODEL1B5TeVye);
    gPOLMODEL1B5TeV->SetName("gPOLMODEL1B5TeV");
    gPOLMODEL1B5TeV->SetLineColor(kPink+9);
    gPOLMODEL1B5TeV->SetFillColor(kPink+9);
    gPOLMODEL1B5TeV->SetFillColorAlpha(kPink+9,0.5);
    gPOLMODEL1B5TeV->SetFillStyle(3306);
    gPOLMODEL1B5TeV->Draw("4 same");
      
  }

}
