// -*- C++ -*-
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   cut.h
 * Author: candan
 *
 * Created on August 31, 2017, 12:57 PM
 */

#ifndef CUTS_H
#define CUTS_H

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


//for pp //
TCut mycut="(hiBin < 181) && Btrk1Pt > 1.0 && Btrk2Pt > 1.0 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4 && abs(Btrk2Eta-0.0) < 2.4 && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&Btrk2highPurity&&abs(Btrk1Eta)<2.4&&abs(Btrk2Eta)<2.4&&Btrk1Pt>1.&&Btrk2Pt>1.&&abs(Btktkmass-1.019455)<0.015) && (abs(PVz)<15&&pclusterCompatibilityFilter&&pprimaryVertexFilter) && (Btrk1PixelHit + Btrk1StripHit > 10) && (Btrk2PixelHit + Btrk2StripHit > 10) && (Btrk1PtErr/Btrk1Pt < 0.1)&& (Btrk2PtErr/Btrk2Pt < 0.1) && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18 && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18&&((Bpt>5&&Bpt<10&&BDT_pt_5_10>0.32)||(Bpt>10&&Bpt<15&&BDT_pt_10_15> 0.29)||(Bpt>15&&Bpt<20&&BDT_pt_15_20>0.35)||(Bpt>20&&Bpt<50&&BDT_pt_20_50>0.33))&&abs(PVz)<15&&pclusterCompatibilityFilter&&pprimaryVertexFilter && phfCoincFilter2Th4)";

//Event InfoFinder
          TLorentzVector* b4P = new TLorentzVector;;
        
          b4P->SetPtEtaPhiM(BInfo->pt[j],BInfo->eta[j],BInfo->phi[j],BInfo->mass[j]);
            float By =b4P->Rapidity();
//1st muon
          b4P->SetPtEtaPhiM(MuonInfo->pt[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]],MuonInfo->eta[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]],MuonInfo->phi[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]],MUON_MASS);  
            Float_t mu1px,mu1py,mu1pz,mu1E;
            mu1px = b4P->Px();
            mu1py = b4P->Py();
            mu1pz = b4P->Pz();
            mu1E = b4P->E();


 //2nd muon 
           b4P->SetPtEtaPhiM(MuonInfo->pt[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]],MuonInfo->eta[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]],MuonInfo->phi[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]],MUON_MASS);
           Float_t mu2px,mu2py,mu2pz,mu2E;
           mu2px = b4P->Px();
           mu2py = b4P->Py();
           mu2pz = b4P->Pz();
           mu2E = b4P->E();    
     
//Bmumu mass      
           b4P->SetPxPyPzE(mu1px+mu2px,mu1py+mu2py,mu1pz+mu2pz,mu1E+mu2E); 
           float Bmumumass= b4P->Mag(); 
              
//Bmu1pt ve Bmu1eta cut 
       float Bmu1pt = MuonInfo->pt[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
       float Bmu1eta = MuonInfo->eta[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]; 
       float Bmu1phi = MuonInfo->phi[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
	   float Bmu2pt = MuonInfo->pt[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];	
       float Bmu2eta = MuonInfo->eta[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
       float Bmu2phi = MuonInfo->phi[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];


///Bmu1TMOneStationTight & Bmu2TMOneStationTight  boolen

      bool Bmu1TMOneStationTight = MuonInfo->TMOneStationTight[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];  
      bool Bmu2TMOneStationTight = MuonInfo->TMOneStationTight[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
        
 //Bmu1InPixelLayer && Bmu2InPixelLayer       
      int Bmu1InPixelLayer = MuonInfo->i_nPixelLayer[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]]; 
      int Bmu2InPixelLayer= MuonInfo->i_nPixelLayer[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];   

 //Bmu1InStripLayer && Bmu2InStripLayer 
      int Bmu1InStripLayer = MuonInfo->i_nStripLayer[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];                                                
	  int Bmu2InStripLayer = MuonInfo->i_nStripLayer[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]; 


 //Bmu1dxyPV && Bmu2dxyPV
      float Bmu1dxyPV = MuonInfo->dxyPV[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];                                                                        
	  float Bmu2dxyPV = MuonInfo->dxyPV[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];
 
 //Bmu1dzPV && Bmu2dzPV
      float Bmu1dzPV = MuonInfo->dzPV[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];                                                                       
	  float Bmu2dzPV = MuonInfo->dzPV[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];

//Bmu1isGlobalMuon && Bmu2isGlobalMuon
     bool Bmu1isGlobalMuon= MuonInfo->isGlobalMuon[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];                                                            
	 bool Bmu2isGlobalMuon= MuonInfo->isGlobalMuon[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]]; 
     
     bool Btrk1highPurity = TrackInfo->highPurity[BInfo->rftk1_index[j]]; 
     bool Btrk2highPurity = TrackInfo->highPurity[BInfo->rftk2_index[j]]; 
     //Bmass
     float Bmass = BInfo->mass[j];
     
     float Btrk1Pt= TrackInfo->pt[BInfo->rftk1_index[j]];  
     float Btrk2Pt = TrackInfo->pt[BInfo->rftk2_index[j]];

    
//Bchi2cl
     float Bchi2cl= TMath::Prob(BInfo->vtxchi2[j],BInfo->vtxdof[j]);
//(Bd0/Bd0Err)
     float Bd0= TMath::Sqrt((BInfo->vtxX[j]-EvtInfo->PVx)*(BInfo->vtxX[j]-EvtInfo->PVx)+(BInfo->vtxY[j]-EvtInfo->PVy)*(BInfo->vtxY[j]-EvtInfo->PVy));                                                                                                                                    
     float Bd0Err = TMath::Sqrt(BInfo->vtxXErr[j]*BInfo->vtxXErr[j]+BInfo->vtxYErr[j]*BInfo->vtxYErr[j]);


//Blxy
     float Blxy = ((BInfo->vtxX[j]-EvtInfo->PVx)*b4P->Px() + (BInfo->vtxY[j]-EvtInfo->PVy)*b4P->Py())/BInfo->pt[j]; 

//Btktkmass
//first track:

     b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk2_index[j]],TrackInfo->eta[BInfo->rftk2_index[j]],TrackInfo->phi[BInfo->rftk2_index[j]],KAON_MASS);       
     float  tk1px = b4P->Px();
     float  tk1py = b4P->Py();
     float  tk1pz = b4P->Pz();
     float  tk1E = b4P->E();

//2nd track
     b4P->SetPtEtaPhiM(TrackInfo->pt[BInfo->rftk1_index[j]],TrackInfo->eta[BInfo->rftk1_index[j]],TrackInfo->phi[BInfo->rftk1_index[j]],KAON_MASS);
     float   tk2px = b4P->Px();
     float   tk2py = b4P->Py();
     float   tk2pz = b4P->Pz();	
     float   tk2E = b4P->E();

//Btktkmass
    b4P->SetPxPyPzE(tk1px+tk2px,tk1py+tk2py,tk1pz+tk2pz,tk1E+tk2E);                     
    float Btktkmass = b4P->Mag();
      
//cos(Bdtheta)
    TVector3* bP = new TVector3;
    TVector3* bVtx = new TVector3;
    bP->SetPtEtaPhi(BInfo->pt[j],BInfo->eta[j]*0,BInfo->phi[j]);
    bVtx->SetXYZ(BInfo->vtxX[j]-EvtInfo->PVx,
    BInfo->vtxY[j]-EvtInfo->PVy,
    BInfo->vtxZ[j]*0-EvtInfo->PVz*0);
    float Bdtheta = bP->Angle(*bVtx);



//New Variables//
	float PVz = EvtInfo->PVz;


    float  Btrk1Dz1 = TrackInfo->dz1[BInfo->rftk1_index[j]];
    float  Btrk1Dxy1 = TrackInfo->dxy1[BInfo->rftk1_index[j]];
    float  Btrk1DxyError1= TrackInfo->dxyerror1[BInfo->rftk1_index[j]];
    float  Btrk1DzError1= TrackInfo->dzerror1[BInfo->rftk1_index[j]];

    float  Btrk2Dz1 = TrackInfo->dz1[BInfo->rftk2_index[j]];
    float  Btrk2Dxy1 = TrackInfo->dxy1[BInfo->rftk2_index[j]];
    float  Btrk2DxyError1= TrackInfo->d0error[BInfo->rftk2_index[j]];
    float  Btrk2DzError1= TrackInfo->dzerror[BInfo->rftk2_index[j]];


	float BsvpvDistance =  BInfo->svpvDistance[j];
	float BsvpvDisErr = BInfo->svpvDisErr[j];
	float Balpha  = BInfo->alpha[j];
	float Bpt = BInfo->pt[j];
	
	float Btrk1Eta = TrackInfo->eta[BInfo->rftk1_index[j]];
	float Btrk2Eta = TrackInfo->eta[BInfo->rftk2_index[j]];
	float Btrk1PtErr = TrackInfo->ptErr[BInfo->rftk1_index[j]];
	float Btrk2PtErr = TrackInfo->ptErr[BInfo->rftk2_index[j]];

    int  Btrk1PixelHit = TrackInfo->pixelhit[BInfo->rftk1_index[j]];
    int  Btrk1StripHit = TrackInfo->striphit[BInfo->rftk1_index[j]];
    int  Btrk2PixelHit = TrackInfo->pixelhit[BInfo->rftk2_index[j]];
    int  Btrk2StripHit = TrackInfo->striphit[BInfo->rftk2_index[j]];


	float Btrk1Chi2ndf = TrackInfo->chi2[BInfo->rftk1_index[j]]/TrackInfo->ndf[BInfo->rftk1_index[j]];
	float Btrk2Chi2ndf = TrackInfo->chi2[BInfo->rftk2_index[j]]/TrackInfo->ndf[BInfo->rftk2_index[j]];


	int Btrk1nPixelLayer = TrackInfo->nPixelLayer[BInfo->rftk1_index[j]];
    int Btrk1nStripLayer = TrackInfo->nStripLayer[BInfo->rftk1_index[j]];

	int Btrk2nPixelLayer = TrackInfo->nPixelLayer[BInfo->rftk2_index[j]];
    int Btrk2nStripLayer = TrackInfo->nStripLayer[BInfo->rftk2_index[j]];

    bool Bmu1isTrackerMuon = MuonInfo->isTrackerMuon[BInfo->uj_rfmu1_index[BInfo->rfuj_index[j]]];
    bool Bmu2isTrackerMuon = MuonInfo->isTrackerMuon[BInfo->uj_rfmu2_index[BInfo->rfuj_index[j]]];


	//Mass Recalculation//

	TLorentzVector* b4P1Pi = new TLorentzVector;
	TLorentzVector* b4P2Pi= new TLorentzVector; 
	TLorentzVector* b4P1K= new TLorentzVector;
	TLorentzVector* b4P2K= new TLorentzVector; 

	TLorentzVector* b4PiK= new TLorentzVector;
	TLorentzVector* b4KPi= new TLorentzVector; 
	TLorentzVector* b4PiPi= new TLorentzVector; 
	TLorentzVector* b4KK= new TLorentzVector; 

	TLorentzVector* b4BPiK= new TLorentzVector;
	TLorentzVector* b4BKPi= new TLorentzVector; 
	TLorentzVector* b4BPiPi= new TLorentzVector; 
	TLorentzVector* b4BKK= new TLorentzVector; 

	TLorentzVector* b4mumu= new TLorentzVector;

	TLorentzVector* mumu1 = new TLorentzVector;
	TLorentzVector* mumu2 = new TLorentzVector;



	double PiMass = 0.134976;
	double KMass = 0.493677;
	double MuMass = 0.10566;


	float  Btrk1Phi = TrackInfo->phi[BInfo->rftk1_index[j]];
	float  Btrk2Phi = TrackInfo->phi[BInfo->rftk2_index[j]];



	b4P1Pi->SetPtEtaPhiM(Btrk1Pt,Btrk1Eta,Btrk1Phi,PiMass);
	b4P2Pi->SetPtEtaPhiM(Btrk2Pt,Btrk2Eta,Btrk2Phi,PiMass);
	b4P1K->SetPtEtaPhiM(Btrk1Pt,Btrk1Eta,Btrk1Phi,KMass);
	b4P2K->SetPtEtaPhiM(Btrk2Pt,Btrk2Eta,Btrk2Phi,KMass);

	*b4PiPi = *b4P1Pi + *b4P2Pi;
	*b4PiK = *b4P1Pi + *b4P2K;
	*b4KPi = *b4P1K + *b4P2Pi;
	*b4KK = *b4P1K + *b4P2K;


	float PiPiMass = b4PiPi->M();
	float PiKMass = b4PiK->M();
	float KPiMass = b4KPi->M();
	float KKMass = b4KK->M();


	mumu1->SetPtEtaPhiM(Bmu1pt,Bmu1eta,Bmu1phi,MuMass);
	mumu2->SetPtEtaPhiM(Bmu2pt,Bmu2eta,Bmu2phi,MuMass);
	
	*b4mumu = *mumu1 + *mumu2;



	*b4BPiPi = *b4PiPi + *b4mumu;
	*b4BPiK = *b4PiK + *b4mumu;
	*b4BKPi = *b4KPi + *b4mumu;
	*b4BKK = *b4KK + *b4mumu;



	float PiPiBMass = b4BPiPi->M();
	float PiKBMass = b4BPiK->M();
	float KPiBMass = b4BKPi->M();
	float KKBMass = b4BKK->M();


	



#endif /* CUTS_H */

