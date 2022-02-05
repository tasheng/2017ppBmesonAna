# CMS 2017 pp Bs and B+ Analysis Codes

## Basic Organization of the Analysis Codes

These codes are made to analyze the CMS 2017 pp datasets and obtain the cross section of the B+ and Bs. Therefore, the organization of the codes is BP, Bs, and BsBP. They consist of raw yield extraction, efficiency correction, and cross section calculations. There are also folders named SkimmedSamples and UnskimmedSamples. These are the folders where you put your input files to perform the analysis

To obtain the codes, do:

git clone https://github.com/MYOMAO/2017ppBmesonAna.git


## Input Data and MC Samples (Ntuplized)

Here we assume we already have the data and official MC samples as well as the non-prompt J/psi samples. The detailed information to obtain the ntuplized data from CMS DAS (AOD data format) is documented on the data processing github:


### On EOS

You could find the samples on lxplus EOS (lxplus7.cern.ch):

/eos/cms/store/group/phys_heavyions/zshi/ForHenrique/FinalSamples/

You should copy the entire folder to your working directory

If you do 

ls /eos/cms/store/group/phys_heavyions/zshi/ForHenrique/FinalSamples/

you will see two folders

Skimmed  UnSkimmed

So basically, copy the files in Skimmed to SkimmedSamples and UnSkimmed to UnskimmedSamples. Then you will have all the files in the right place and ready to perform the analysis.

### On Grendel

Same as before on EOS, on grendel01.mit.edu

/data/szhaozho/2017ppSamplesFinal/

Do the same to copy to the folders to run the codes.


## Cuts Information

These cuts are presented as a string. There 2 ways to apply cuts. One is simply apply it as a string to project to a histogram. 

In this analysis codes, I apply them as along with the candidate ID: [j].  

The details of the cuts can be found at Section 5 of the analysis notes: AN-21-091.

### Quality Cuts - (Save as Pre-Filter)

B+: (((BDT_pt_1_2>0.04 && BsvpvDistance/BsvpvDisErr > 5 && Bchi2cl > 0.05 && Bpt > 1.0 && Bpt < 2.0) || (BDT_pt_1_2>0.04 && BsvpvDistance/BsvpvDisErr > 5 && Bchi2cl > 0.05 && Bpt > 0.0 && Bpt < 1.0) || (((HBHENoiseFilterResult == 1 && pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1) &&  (Bpt > 2 && Bpt < 3 && BDT_pt_2_3 > -0.12) )) ) ||  ( Bpt > 3 && ((HBHENoiseFilterResult == 1 && pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1)  &&  (Bmu1isTriggered == 1 && Bmu2isTriggered == 1 ) &&  (Btrk1Pt > 0.2 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0  && abs(Btrk1Eta-0.0) < 2.4  && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&abs(Btrk1Eta)<2.4&&Btrk1Pt>0.2)  && (Btrk1PixelHit + Btrk1StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1)&& Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18   && (abs(PVz)<15)))))


Bs: ((((abs(Btktkmass-1.019455)<0.015)&& BDT_pt_2_3 > -0.10 && TMath::Abs(Bmumumass-3.096916)<0.15 && Bpt > 0 && Bpt < 5 && (abs(Btrk1Eta)<2.4 && abs(Btrk2Eta)<2.4 && Btrk1Pt>0.0 && Btrk2Pt>0.0) && Btrk1Pt > 0.2 && Btrk2Pt > 0.2  && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0)  && ( (Bpt < 2 && Bpt > 0 && BDT_pt_1_2 > -0.38 ) || (Bpt < 3 && Bpt > 2 && BDT_pt_2_3 > -0.05 ) || (Bpt < 5 && Bpt > 3 && BDT_pt_3_5 > -0.40)  )))  ||  ( Bpt > 3 && ((HBHENoiseFilterResult == 1 && pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1 && (abs(PVz)<15))  &&  (Bmu1isTriggered == 1 && Bmu2isTriggered == 1 ) &&  (Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0)    && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon)  && ( Btrk1Pt > 0.2 && Btrk2Pt > 0.2 && abs(Btrk1Eta-0.0) < 2.4 && abs(Btrk2Eta-0.0) < 2.4  && Btrk1highPurity  && Btrk2highPurity  && Btrk1PixelHit + Btrk1StripHit > 10  && Btrk2PixelHit + Btrk2StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1)  &&  (Btrk2PtErr/Btrk2Pt < 0.1)    && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18   && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18 ))

Notices: we also include BDT for the low pT (pT < 5 GeV/c), which is irrelavent in this analysis. 

### BDT Cuts

The BDT cuts are optimized for each pT bins. The Bs is optmized: 1-2-3-5-7-10-15-20-30-50 and the B+ is optimized: 1-2-3-5-7-10-15-20-50-100. The structure is BDT_pt_X_Y:

B+: ( (Bpt > 3 && Bpt < 5 && BDT_pt_3_5 > 0.08) || (Bpt > 5 && Bpt < 7 && BDT_pt_5_7 > 0.06) || (Bpt > 7 && Bpt < 10 && BDT_pt_7_10 > 0.07) || (Bpt > 10 && Bpt < 15 && BDT_pt_10_15 > 0.08) || (Bpt > 15 && Bpt < 20 && BDT_pt_15_20 > 0.12)  || (Bpt > 20 && Bpt < 50 && BDT_pt_20_50 > 0.12) || (Bpt > 50 && Bpt < 100))

Bs: ( (Bpt > 5 && Bpt < 7 &&  BDT_pt_5_7 > -0.28) || (Bpt > 7 && Bpt < 10 &&  BDT_pt_7_10 > -0.30) || (Bpt > 10 && Bpt < 15 &&  BDT_pt_10_15 > -0.30) || (Bpt > 15 && Bpt < 20 &&  BDT_pt_15_20 > -0.25 ) || (Bpt > 20 && Bpt < 50 &&  BDT_pt_30_50 > -0.25 )  || (Bpt > 50) )

### Total Cuts

Total cuts is basically the Quality Cuts && BDT cuts. Here is the the total cuts:

B+:  (((BDT_pt_1_2>0.04 && BsvpvDistance/BsvpvDisErr > 5 && Bchi2cl > 0.05 && Bpt > 1.0 && Bpt < 2.0) || (BDT_pt_1_2>0.04 && BsvpvDistance/BsvpvDisErr > 5 && Bchi2cl > 0.05 && Bpt > 0.0 && Bpt < 1.0) )   || (((HBHENoiseFilterResult == 1 && pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1)  &&  (Bmu1isTriggered == 1 && Bmu2isTriggered == 1 ) &&  (Btrk1Pt > 0.2 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0 && Bpt > 2 && abs(Btrk1Eta-0.0) < 2.4  && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&abs(Btrk1Eta)<2.4&&Btrk1Pt>0.2)  && (Btrk1PixelHit + Btrk1StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1)&& Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18   && (abs(PVz)<15))  &&( (Bpt > 3 && Bpt < 5 && BDT_pt_3_5 > 0.08) || (Bpt > 5 && Bpt < 7 && BDT_pt_5_7 > 0.06) || (Bpt > 7 && Bpt < 10 && BDT_pt_7_10 > 0.07) || (Bpt > 10 && Bpt < 15 && BDT_pt_10_15 > 0.08) || (Bpt > 15 && Bpt < 20 && BDT_pt_15_20 > 0.12)  || (Bpt > 20 && Bpt < 50 && BDT_pt_20_50 > 0.12) || (Bpt > 50 && Bpt < 100))) || ((HBHENoiseFilterResult == 1 && pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1) &&  (Bpt > 2 && Bpt < 3 && BDT_pt_2_3 > -0.12) ) ))


Bs: ((((abs(Btktkmass-1.019455)<0.015)&& BDT_pt_2_3 > -0.10 && TMath::Abs(Bmumumass-3.096916)<0.15 && Bpt > 0 && Bpt < 5 && (abs(Btrk1Eta)<2.4 && abs(Btrk2Eta)<2.4 && Btrk1Pt>0.0 && Btrk2Pt>0.0) && Btrk1Pt > 0.2 && Btrk2Pt > 0.2  && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0)  && ( (Bpt < 2 && Bpt > 0 && BDT_pt_1_2 > -0.38 ) || (Bpt < 3 && Bpt > 2 && BDT_pt_2_3 > -0.05 ) || (Bpt < 5 && Bpt > 3 && BDT_pt_3_5 > -0.40)  )))  ||  ( ( (Bpt > 5 && Bpt < 7 &&  BDT_pt_5_7 > -0.28) || (Bpt > 7 && Bpt < 10 &&  BDT_pt_7_10 > -0.30) || (Bpt > 10 && Bpt < 15 &&  BDT_pt_10_15 > -0.30) || (Bpt > 15 && Bpt < 20 &&  BDT_pt_15_20 > -0.25 ) || (Bpt > 20 && Bpt < 50 &&  BDT_pt_30_50 > -0.25 )  || (Bpt > 50) ) && ((HBHENoiseFilterResult == 1 && pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1 && (abs(PVz)<15))  &&  (Bmu1isTriggered == 1 && Bmu2isTriggered == 1 ) &&  (Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0)    && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon)  && ( Btrk1Pt > 0.2 && Btrk2Pt > 0.2 && abs(Btrk1Eta-0.0) < 2.4 && abs(Btrk2Eta-0.0) < 2.4  && Btrk1highPurity  && Btrk2highPurity  && Btrk1PixelHit + Btrk1StripHit > 10  && Btrk2PixelHit + Btrk2StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1)  &&  (Btrk2PtErr/Btrk2Pt < 0.1)    && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18   && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18 ))

## Prerequisite - MC Reweighting

Before we start doing the analysis, we have some prerequisites. We need to perform global event level to reweight the MC to make it better match to the data. There are two reweight processes: PVz and Bpt. These are all included in the folder: MCReweight

Before reweihghting, we could produce some plots to show the MC performances in the GEN and RECO level. To do so, first get in the folder:

cd MCReweight/GenInfo

To obtain the Gen info for B+, simply do 

root -b -l -q PlotGen.C'(0)'

Here the argument 0 is for B+ and 1 for Bs. The plots are saved as:

Pthat Distribution: MCPlots/BPpthat.png 

Generated B+ pT distribution: MCPlots/BPGpt.png

Generated B+ J/psi distribution: MCPlots/BPJPsiPt.png

Same plots are generated for Bs when running 

root -b -l -q PlotGen.C'(1)'


### PVz Reweighting

To reweight the PVz, we simply produce the data and MC PVz distribution, fit them with Gaussian functions, and then define the weight as the ratio of Data Gaussian function to the MC Gaussian function. Finally, recheck the data-MC agreement by applying the weight to MC and compare it to data. To produce the stduies, simply do 

cd PVZ

root -b -l -q PVzMC.root'(0)'

Again, here the argument 0 is for B+ and 1 for Bs. The plots are saved as:

BPPVZMCData.png

There are 3 panels in the plots. The left is data, the middle is MC, and the right the PVz reweighted MC and comparison with data. Also, the Gaussian fit function parameters (constant, mean, and width) are printed out as well when running the codes.

Same for Bs as mentioned above. We simply need to change the argument from 0 to 1. 

root -b -l -q PVzMC.root'(1)'

We will get 

BsPVZMCData.png


### B pT Reweighting


Likewise, to run the B pT reweighting, we go to the folder Bpt

cd Bpt

The B pT shape in the data are already computed from raw yield extraction for each pT bin and saved at:

../../BP/RawYieldFits/ROOTfiles/yields_Bp_binned_pt.root

../../Bs/RawYieldFits/ROOTfiles/yields_Bs_binned_pt.root

The signal raw yield in MC is obtained by doing GenMatch selection. To compute the Bpt weight, we simply need to fit the raw yield ratio of data to MC as a function of Bpt with a function. To do so, run

root -b -l -q BptReweight.C'(0)'

Again, here the argument 0 is for B+ and 1 for Bs. The plots are saved as:

DataMCCompSide_BP.png and DataMCRatio_BP.png

The data/MC ratio as well as the fits function can be found at DataMCRatio_BP.png.

Same for Bs. We simply need to replace the argument 0 to 1 and run the codes

root -b -l -q BptReweight.C'(1)'


To obtain the fit function on the Bpt shap over a fine binning for B+, simply do

root -b -l -q FitBptShape.C'(0)'


Again, here the argument 0 stands for B+ and 1 stands for Bs.


After runing the codes, fit on data/MC ratio as a function of B+ pT plot is produced:

BPPtWeight.png


The fit function as well as the chi^2/ndf will also be printed out:

8.646057/(x*x) +0.425834*TMath::Log(x) +-0.148249 

Chi2ndf = 4.90555


Here the B+ has taken a function form: p0/x^2 + p1 * log(x) - p2

For Bs, just follow the same procedure to run 


root -b -l -q FitBptShape.C'(1)'

After runing the codes, fit on data/MC ratio as a function of Bs pT plot is produced:

BsPtWeight.png

Again, the fit function as well as the chi^2/ndf are printed out:

1.000000/(x*x) +0.463653*TMath::Log(x) +-0.238962 

Chi2ndf = 8.94302

Then, we could enter these fit functions to the code for Bs and B+ efficiency correction:

https://github.com/MYOMAO/2017ppBmesonAna/blob/master/BP/EffAna/MCEff.C#L790

https://github.com/MYOMAO/2017ppBmesonAna/blob/master/Bs/EffAna/MCEff.C#L982

At this point, the Bpt weight is also included. 

# Nominal Analysis


## Raw Yield Fit

The analysis consist of two parts: raw yield extraction and efficiency correction. To perform the signal raw yield extraction, we will use the skimmed files. To run the codes, for example, for B+, simply do:

cd BP/RawYieldFits

source doRoofit.sh

Here you can change the flags in doRoofit.sh

DOANALYSISPbPb_ROOFIT_FULL_BP                means the include pT bin

DOANALYSISPbPb_ROOFIT_BINNED_PT_BP           means the pT differential bin

DOANALYSISPbPb_ROOFIT_BINNED_Y_BP            means the rapidity bin

DOANALYSISPbPb_ROOFIT_BINNED_MULTI_BP        means the multiplicity bin

The binning of the fits can be modified in parametersNew.h

For pT differential binning, we can modify the 

double ptBins_bp[nBins_bp+1] = {7,10,15,20,50};


You can change the flags and perform the fits. The fit plots are saved plotFits/final_roofit/

The yield information are stored in the ROOTfiles/ folder

The pt differential file is saved as 

The procedure to obtain raw yield in Bs is basically the same as BP.

## Efficiency Correction

The next step is to obtain the efficiency correction. The steps are also quite straightforward. To Run the efficiency simply do:

cd BP/EffAna

root -b -l -q MCEff.C'(1,0)'

The first option means the application of Tag and Probe scale factor: 0 mean no TnP applied and 1 means with TnP applied

The second option means the rescaling of BsvpvDisErr: 0 mean no rescaling and 1 means with scale. The detailed of rescaling BsvpvDisErr can be found in the AN appenidx section C.

Currently, we find that the recaling does not significantly change the efficiency. Therefore, using 0 or 1 in the second option will not change the results. 

Inside the efficiency codes

You can find the binnings:

const int NPtBins1D = 4;
double  PtBin1D[NPtBins1D + 1] = {7,10,15,20,50};

This is the binning for 1D pT-binned efficiency correction

The corresponding histogram is: 

TH1D * Eff1DRECOHis = new TH1D("Eff1DRECOHis","",NPtBins1D,PtBin1D);

For multiplicity, we have 

const int NMultiBin = 10;
double  MultiBin1D[NMultiBin + 1] = {0,15,25,30,35,40,50,65,80,100,130};

The corresponding histogram is: 

TH1D * Eff1DRECOMultHis = new TH1D("Eff1DRECOMultHis","",NMultiBin,MultiBin1D);

We also have the method with the 2D map, which is the a TH2D of efficiency vs pT and |y|:

TH2D * invEff2D = (TH2D * ) EvtWeightGenHis->Clone("invEff2D");

These codes run on the unskimmed MC files (not flattened)

Once you finishing running the codes, the plots are saved in the root file 

EffFineBDT.root

The 1D and 2D efficiency histograms are saved inside the file

Again, Bs efficiency follows the same procedures

## Cross Section Calculations

The cross section calculations are also straightforward. Since we have obtained the raw yield and efficiency, to obtain the cross section, we simply run 

cd BP/EffAna/

For pT, simply run

root -b -l -q CrossSectionAna.C'(1)'

Here again, the first option 1 enables TnP and 0 disables TnP. 

For Multiplicity, simply run

root -b -l -q CrossSectionAnaMult.C'(1)'

Here again, the first option 1 enables TnP and 0 disables TnP. 


The output file of the cross section is saved at FinalFile/

FinalFile/BPPPCorrYieldPT.root and FinalFile/BPPPCorrYieldMult.root

Here, when we run the cross section, we need to make sure that the bin is indeed the same. For example, in here, NPtBins = 4, so in the CrossSectionAna.C, we make sure that 

const int NBins = 4;

Here, we again use two efficiency correction methods. The binned pT one is at 

TH1D * CorrDiffHisBin = (TH1D * ) CorrDiffHis ->Clone("CorrDiffHisBin");

The 2D map one is at:

TH1D * CorrDiffHis = new TH1D("hPtSigma","",NBins,ptBins);

Repeat this for Bs, we could get the cross section for Bs

At this point we get the Bs and BP cross sections as function of B meson pT and event multiplicity. We will proceed to the next step to obtain RAA and Bs/B+ ratios from these pp cross sections.


## Results Involving with Bs and B+ Cross Sections


To obtain more relevant physics results for the B-meson measurement to compare with theoretical calculations, we need to calculate the physical obserables such as Bs/B+ production yield ratio in pp and PbPb collisions and Nuclear modification factor RAA. The folder BsBPFinalResults contains the codes to produce the plots.


### Bs and B+ RAA

With the PbPb production yield and pp cross section, we can calculate the RAA. To obtain the RAA plots, first go to the folder BsBPFinalResults/RAA

cd BsBPFinalResults/RAA

To obatin B+ RAA, simply do:

Note that we use a fiducial region of |By| > 1.5 for B pT < 10 GeV/c.



### Bs/B+ 

To Run the Bs/B+ as a function of multiplicity. It is also very simple. We could go to 

cd BsBPFinalResults/BsBPRatio

And then run 

root -b -l -q PlotsBsBPRatio.C'(0)'

The second argument: 0 stands for Binned pT efficiency correction while 1 stands for 2D map efficiency correction

The plots are stored at Plots/Binned and Plots/2DMap

If we run 

root -b -l -q BPRAA.C

The B+ RAA plots are stored at RAAPlots/BP/

There are basically 3 types of plots:

1. The cross section of 2017 pp and production yield of 2018 PbPb, both plotted together, in linear scale: RAAPlots/BP/BPPbPbPPCross.png and log scale RAAPlots/BP/BPPbPbPPCrossLog.png. 

2. We also produce a comparison plot between the cross section with and without fiducial region: RAAPlots/BP/BPFidOrNotComp.png

3. The RAA, that is the ratio between production yield of 2018 PbPb to 2017 pp located at RAAPlots/BP/BPRAA.png. A comparison with 2015 B+ RAA is also made at RAAPlots/BP/BPRAACompairson.png.

A root file containing all these plots are also saved at OutFile/BPRAA.root

Like wise, to obatin the RAA of Bs, simply do 

root -b -l -q BsRAA.C

Same plots are saved at RAAPlots/BP/ and analysis histogram stored at OutFile/BsRAA.root


## Generate Comparison Plots 

To run the comparison of the 2017 pp Bs and B+ results with 2015 pp and FONLL calculations, go to the folder:

Caveat - fiducial region: for the full 2015 pp results, since the measurement |By| > 2.4 for B pT < 10 GeV/c unlike the 2018 PbPb where a fiducial region |By| > 2.4 for B pT < 10 GeV/c. Therefore, we have produced two sets of analysis. One with the fiducial region cut of |By| > 1.5 for B pT < 10 GeV/c. The other one is without the fiducial region (still keeping |By| < 2.4 for B pT < 10 GeV/c).

We can change the configuration to remove the fiducial region. Here we do not go through the details about that. I have produced two files one for Bs and one for B+ where the fiducial region is removed in order to compare with 2015 pp results. You can find the cross section files at:

BsBPFinalResults/Comparisons/NoFiducial/FinalFiles/


### Fiducial Region 

Since the 2015 pp reference does not have a fiducial region, in our comparison, we will only include FONLL calculations with a fiducial 1.5 < |By| < 2.4 selection for B pt < 10 GeV/c. For B pt > 10 GeV/c, we will simply use |By| < 2.4. We put all the codes under the folder: Comparisons/Fiducial

cd Comparisons/Fiducial


To get the comparison of B+ 2017 data with FONLL calculations, simply do 

root -b -l -q BPComparison.C 

The plots will be stored at Plots/BP/

We could look at the plot 

Plots/BP/BPCrossCompLog.png 

We can see there are two panels. The top panel is the comparison of the B+ pp cross section with FONLL cross section calculation with 2 methods: the binned pt (blue) and 2D map (orange). The bottom panel is their ratio and comparison with unity (red line). 

Similarly, we could run 

root -b -l -q BsComparison.C 

to obtain the comparison of Bs cross section with FONLL calculations.

Again, the files can be found at Plots/Bs just like B+ stated above

### Without Fiducial Region  (Optional)

This is optional, but if you also want to run the results without fiducial region, simply go to the folder


### Generating FONLL Files (Optional)

This is also optional, if you want to generate FONLL files, first obtain the calculations from the FONLL website:

http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html

Turn on the "Include PDFs uncertainties" and "uncertainties range from scale and masses" options.

Then, to generate the FONLL results. Then dump the results to the dat file

Then go to  

Run 

root -b -l -q 

to generate the FONLL with the desired binning defined in the code:


# Systematic Studies (Binned pT Method)

Now we have finished the nominal analysis and obtain the central values of the B-meson cross section. Since the analysis involves with two parts: signal raw yield extraction and efficiency correction. For the binned pT correctio, we have also obtained the statistical uncertainties. For the 2D Map method, an addition data boostraping approaching is needed to extract the statistical uncertainties of the results. The signal raw yield extraction ,which is called PDF variation, is the same for both methods. Henrique Legoinha is the contact person for signal extraction for everything involving with efficiency correction.

The following sections focus on the systematic uncertainties on the efficiency using the binned pT correction method. 

## TnP Systematics

One of the systematic uncertainties is related to muon. A data driven method, name Tag and Probe (TnP) is applied as the nominal of our analysis. Here we evaluate the TnP systematic uncertainties. To do so, fist obtain the TnP information for the MC samples under the folder MakeMCTnP. 

cd MakeMCTnP

To produce the TnPInfo tree with the scale factor and their uncertainties the based on the TnP header file, simply do

root -b -l -q TnPWeight.C'(0)'

Again, the argument 0 stands for B+ while the argument 1 stands for Bs.

Then, the output file is named "BPTnPInfo.root" (BsTnPInfo.root for Bs). 

Then merging the TnP input file with the MC file, we get the MC with TnPInfo tree that can be used in the analysis.

Now, to obtain the TnP systematic uncertainties, it turns out that it is also in the code MCEff.C. Thus, we simply run (not need to run again if you already run it):

cd BP/EffAna/

root -b -l -q MCEff.C'(1,0)'

We will save the TnP systematic in the file NewEff2DMap/BPSyst.root

When we open the file, we can see:

root -l NewEff2DMaps/BPSyst.root 
root [0] 
Attaching file NewEff2DMaps/BPSyst.root as _file0...
(TFile *) 0x113db70
root [1] .ls
TFile**		NewEff2DMaps/BPSyst.root	
 TFile*		NewEff2DMaps/BPSyst.root	
  KEY: TH1D	Eff1DHis;1	
  KEY: TH1D	Eff1DHisTnPUp;1	
  KEY: TH1D	Eff1DHisTnPDown;1	
  KEY: TH1D	Eff1DHisBpt;1	
  KEY: TH1D	Eff1DHisBDT;1	
  KEY: TH1D	Eff1DHisMult;1	
  KEY: TH1D	Eff1DHisTnPUpMult;1	
  KEY: TH1D	Eff1DHisTnPDownMult;1	
  KEY: TH1D	Eff1DHisBptMult;1	
  KEY: TH1D	Eff1DHisBDTMult;1	

In fact, the file also contains that the Splot weight variation for MC/Data Discrepancy systamtics as well as Bpt Systematics

Here, the histograms of Eff1DHis is the nominal, Eff1DHisTnPUp is the TnP variation up, and the Eff1DHisTnPDown is TnP variation down. Basically, we take the percent variation of TnP Up/Down to the nominal as the TnP systematic uncertainties

Finally, to obtain the TnP systematic uncertainties, we go to the folder SystStudies/

cd SystStudies

Then run the code:

root -b -l -q PlotEffSyst.C'(0)'

Again, the argument 0 stands for B+ while the argument 1 stands for Bs.

The plots of TnP Systematic for B+ is saved at SystPlots

Inside the folder SystPlots, it constain both B mesons:

BP/ and Bs/

Inside one of the Bmesons, it has both Pt and Multiplicity:

Pt/ and Mult/

Inside each of them, there two types of plots saved.

1. The comparison between the variated efficiency with the nominal efficiency. All of them overlay in a single plot with legend identifying them.

In this case, for TnP, the plot for Pt is located at: SystPlots/BP/Pt/TnPSystComp.png

2. The percent deviation of the variated efficiciency (again, quoted as systematic uncertainties)

In this case, for TnP, the plot for Pt is located at: SystPlots/BP/Pt/TnPSystRatio.png

We will collect the numbers here and enter them in the array for TnP systematic uncertainties in the systematic summary stage to plot the systematic uncertainties for each pT and multiplicity bin for B+ and Bs.

Likewise, to do this for Bs, simply run 

root -b -l -q PlotEffSyst.C'(0)'

Same plots will be produced under the folder 

SystPlots/Bs/

## MC-Data Discrepancy 

As mentioned above, the codes MCEff.C and PlotEffSyst.C includes also MC-Data discrepancy using Splot Techniques. Here, the prerequite is the input weight files: BPw.root and Bsw.root under the folder BP/EffAna/BDTWeights. They can be obtained from Splot studies done Henrique. The plots are saved as, for instance, for B+ Pt bins

SystPlots/BP/Pt/MCDataSystComp.png and SystPlots/BP/Pt/MCDataSystRatio.png


## Bpt Systematics

Again, as mentioned above, the codes MCEff.C and PlotEffSyst.C includes also Bpt shape. The plots are saved as, for instance, for B+ Pt bins

SystPlots/BP/Pt/BptSystComp.png and SystPlots/BP/Pt/BptSystRatio.png


## MC Stat Systematics 

It looks like in this case for binned pT correction, the way to calculate MC stat is undocumented in the previous analysis. So here for binned pT correction, we believe this uncertainty is small and decide to neglect it.  

## Generating Summary Plots

Now, with all systematic uncertainties in hand, we are ready to plot them for B+ and Bs pT and event multiplicity bins. The codes are located at: PlotSystSummary

cd PlotSystSummary

Inside the folder, you will find the header files:

BpSystValues.h and BsSystValues.h

If you open the header file, you will see

NPtBins = 7

NCentBins = 2

This stands for the pt and multiplicity binning 

So basically, we have 4 types of systematic uncertainties:

1. Global: they are basically the systematic uncertainties on lumi and Bs and B+ decay branching ratios

2. Signal: the systematic uncertainties related to (Henrique has this number)

3. Eff: the efficiency correction systematic uncertainties: basically the quadruture of Bpt shape and MC-Data discrepancy as well as the tracking efficiency systematic uncertainties above added up 

4. TnP: the efficiency correction based on muon TnP, which we have calculated above.

The total here is the quadrature sum of Eff, TnP, and signal extractions. 

**Currently, the numerical values of this table is not yet finalized. Once the systematic uncertainties studies are finalized, you can refresh the header files and rerun the codes to obtain the plot on the summary of systematic uncertainties. 


Anyhow, to plot the systematic uncertainties for pt, event multiplicity, as well as inclusive bins, for B+, we can do:

root -b -l -q PlotSystBp.C

The plots are saved at PlotSyst/Bp/

Pt: PlotSyst/Bp/BpPtSyst.png or PlotSyst/Bp/BpPtSyst.pdf

Mult: PlotSyst/Bp/BpCentSyst.png or PlotSyst/Bp/BpCentSyst.pdf

Inc: PlotSyst/Bp/BpIncSyst.png or PlotSyst/Bp/BpIncSyst.pdf

All plots in 1 canvas: PlotSyst/Bp/BpSysSumPlot.png or PlotSyst/Bp/BpSysSumPlot.pdf

Same for Bs, simply do 

root -b -l -q PlotSystBs.C




# Systematic Studies (2D Map Method)

Similar to binned pT method, 2D map efficiency correction also has the systematic uncertainties including:

PDF variation - same as the binned pT method since the efficiency correction does not affect signal raw yield extraction

MC-Data Disagreement

pT Shape Systematics

TnP Systematics

MC Statistics Systematics (New)

## MC-Data Disagreement, pT Shape, and TnP Systematics  

Basically, when we run the codes MCEff.C, we have computed the nominal 2D map as a function of B-meson pT and rapidity |y|: 1/(acc X eff) as well as the (1/acc X eff) varied by sPlot weight for MC-Data Disagreement, Bpt weight to acceount for pT Shape Systematics, and muon TnP scale factors to account for TnP Systematics. All these files are stored at 
## MC-Data Disagreement, pT Shape, and TnP Systematics  

NewEff2DMaps/BPSyst2D.root

NewEff2DMaps/BsSyst2D.root

for Bs and B+ respectfully under the folder EffAna/

To evaluate the <1/(acc X eff)> for nominal and varied cases, go to the folder 2DMapSyst

cd 2DMapSyst

Then run

root -b -l -q CalEffSystBP.C

for B+, a file BPSyst2D.root will be saved at

OutFiles/BPSyst2D.root

The file content is as follows:


 TFile*		OutFiles/BPSyst2D.root	
  KEY: TH1D	Eff2DHis;1	
  KEY: TH1D	Eff2DTnPUpSystHis;1	
  KEY: TH1D	Eff2DTnPDownSystHis;1	
  KEY: TH1D	Eff2DBDTHis;1	
  KEY: TH1D	Eff2DBptHis;1	

The nominal is Eff2DHis, TnP systematic varied up is Eff2DTnPUpSystHis, down is Eff2DTnPDownSystHis, the MC-Data discrepancy is Eff2DBDTHis, and the Bpt shape variation is Eff2DBptHis.


Finally, to evaluate the percent of systematic uncertainties, simply run:

root -b -l -q PlotEffSyst2D.C'(0)'

Here the argument 0 is for B+ and the argument 1 is for Bs,. 


The results will be printed out as percent:

TnP Syst: 0.463625

TnP Syst: 0.445255

TnP Syst: 0.371436

TnP Syst: 0.361883

TnP Syst: 0.425895

TnP Syst: 0.625206

TnP Syst: 0.611058

Bpt Syst: 0.0116134

Bpt Syst: 0.00430247

Bpt Syst: 0.000850791

Bpt Syst: 0.000144326

Bpt Syst: 0.000195505

Bpt Syst: 0.000757004

Bpt Syst: 0.00308549

BDT Syst: 2.08743

BDT Syst: 1.78165

BDT Syst: 0.769814

BDT Syst: 24.561

BDT Syst: 0.821754

BDT Syst: 2.62644

BDT Syst: 0

The plots comparing the varied efficiency to the nominal one are saved as 

SystPlots/BP/Pt/TnPSystComp.png

SystPlots/BP/Pt/MCDataSystComp.png

SystPlots/BP/Pt/BptSystComp.png


The computed systematic uncertainties plots are saved as 

SystPlots/BP/Pt/TnPSystRatio.png

SystPlots/BP/Pt/BptSysRatio.png

SystPlots/BP/Pt/MCDataSystRatio.png


Likewise for Bs, we will do the same as B+ to produce the systematic uncertainties 


## MC Stat Systematics

To evaluatiobn he MC Stat Systematics, just go to the folder


cd MCStatSyst/BP

Then run the code Generate2DMaps.C to generate 10k (default number than can be changed) of 2D maps based on the mean and error of each bin in the 2D map:

root -b -l -q Generate2DMaps.C

The 10k 2D maps will then be saved at:

OutFiles/GenStatSyst.root

with 10k 2D maps of 1/(eff X acc) name as EffBptByInvTrial for Total Efficiency
 
1/(eff) name as SelBptByInvTrial for Selection Efficiency

1/(acc) name as SelBptByInvTrial for Acceptance Efficiency

Next we run:

root -b -l -q MCStatCal.C

Here, we evaluate the <1/(eff X acc)> for each 2D 1/(eff X acc) maps based on the same dataset from the skimmed B+ data and fill each of <1/(eff X acc)> them into a histogram. 

The code will save the plots and also compare the nominal <1/(eff X acc)> (in a red vertical line) to the distribution of the histogram. Ideally, the distribution should be a perfectly symmetric Gaussian. The red vertical line will fall right at the mean of the distribution. If that is not the case, we will need to check the nominal value is correctly input in the code. For B+, since we have 7 bins, it will be at 

meanvec.push_back(53.538932);
	
meanvec.push_back(23.943692);
	
meanvec.push_back(12.038493);
	
meanvec.push_back(5.7727937);
	
meanvec.push_back(3.7848080);
	
meanvec.push_back(3.0837052);
	
meanvec.push_back(2.9782204);

Also make sure that the are input 2D map for the 10k generation are synchronous to the nominal 2D map.

If all above look good, the code MCStatCal.C will use mean/RMS as the MC stat systematics. They will be printed out as

<1/acc eff> Stat Syst: 0.0121105

<1/acc eff> Stat Syst: 0.00857881

<1/acc eff> Stat Syst: 0.00443399

<1/acc eff> Stat Syst: 0.00627536

<1/acc eff> Stat Syst: 0.00750213

<1/acc eff> Stat Syst: 0.0223612

<1/acc eff> Stat Syst: 0.0957221

for each pT bin.

The histigram for MC stat systematics in terms of percentage just like the MC-Data discrepancy, TnP, and Bpt shape systematics, is also saved as 

Plots/MCStatSyst.png




## Report the Systematic Uncertanties to the Final Results


###Cross Section

Finally, with all uncertainties values in hand, we can also add that to the final results in addition to the nominal cnetral values and the statistical uncertainties. To do so simply go back to the BsBPFinalResults/Comparisons/Fiducial/ and check out the codes: BPComparison.C and BsComparison.C


float BPMDDataSyst[NBins];

float BPPDFSyst[NBins];

float BPPtShapeSyst[NBins];

float BPTnPSystDown[NBins];
	
float BPTnPSystUp[NBins];


These are from line 153 to 161. (The number BPPtShapeSyst, the last pT bin 50 - 60 GeV/c of BPMDDataSyst, and the first 5 - 7 GeV/c, and last bins 50 - 60 GeV/c BPPDFSyst are no yet finalized). We will need to update them once we have obtained our final results on systematic uncertainties. 


Now run again the codes: 

cd BsBPFinalResults/Comparisons/Fiducial/

root -b -l -q BPComparison.C

We will get thes plots under the folder for B+ cross section including the total systematic uncertainties drawn in boxes with png and pdf formats: 


Plots/BP/BPCrossONLY.png

Plots/BP/BPCrossONLY.pdf

We also have the comparison plot, which also includes the total systematic uncertainties of B+ 2017 pp and 2018 PbPb drawn in boxes with png and pdf formats: 

Plots/BP/BPCrossCompLog.png

Plots/BP/BPCrossCompLog.pdf

At these point, all information about the measurements are reported in the figures. 

Similarly, to do the same for Bs, just run:

root -b -l -q BsComparison.C


Same plots are saved under the folder: Plots/Bs/


### RAA

Likewise, systematic uncertainties are also included to the RAA. To do so, go back to the RAA folder: BsBPFinalResults/Comparisons/RAA


cd BsBPFinalResults/Comparisons/RAA/


run root -b -l -q BPRAA.C

The output results of B+ RAA with systematic uncertainties as a function of B+ pT is saved at 

RAAPlots/BP/BPRAA.png

and comparision with 2015 B+ RAA results is saved at  

RAAPlots/BP/BPRAACompairson.png


To do this for Bs, simply run 


run root -b -l -q BsRAA.C


The same results for Bs save at RAAPlots/Bs



## Final Results Plots With Correct Style for Publication and Presentation 

In the end, with all the results available, we should make format plots that compile with the CMS rules of publication and presentation. To generate the plot with only 2017 pp and 2018 PbPb, go to the folder MakeFinalPlot/NominalPlots

cd MakeFinalPlots/NominalPlots

### Cross Section

To produce the pp cross section of Bs and B+ as a function of pT in the same plot, simply do:

cd CrossSection/

root -b -l -q plotPt'(1,1,0,1,1)'

The cross section plot files will be saved at figs/pdf/xsec_vsPt.pdf as pdf format and figs/png/xsec_vsPt.png as png

Here, the dataSource folder includes all in the information of the measurement

The files: corryield_pt_Bp_New.txt and corryield_pt_Bs_New.txt contain

pT bin lower and upper, central value of the cross section, statstistical uncertainties, systematic uncertainties, global systematic uncertanties, and abscissae for the bin center for B+ and Bs repectfully

The data information is current manually entered based on the analysis results. You can confirm this from Table 8 and Table 9 in the Results section summarizing the Bs and B+ cross section results the AN-21-091.


### Nuclear Modification Factor

To produce the RAA of Bs and B+ as a function of pT in the same plot, simply do:

cd RAA/

root -b -l -q plotPt'(1,1,0,1,1)'

The RAA plot files will be saved at figs/pdf/BRAA_vsPt.pdf as pdf format and figs/png/BRAA_vsPt.png as png

The files RAA_pt_Bp.txt and RAA_pt_Bs.txt contain 

pT bin lower and upper, central value of RAA, statstistical uncertainties, systematic uncertainties, global systematic uncertanties, and abscissae for the bin center

Again, the data information is current manually entered based on the analysis results. You can confirm this from Table 11 and Table 12 in the Results section summarizing the Bs and B+ cross section results the AN-21-091

Note that, here for the RAA calculations, the systematic uncertainties here assume the cancellation of Branching Ratio systematics. All other systematic uncertainties are added into quadruture. 

When the results have been updated, we will need to update the data in the dataSource folder. It would be benefit to produce a code to automatically generate the txt files in the dataSource to synchronize any changes in the analysis. 


### Cross Section Comparison Plots

To compare our pp results with 2015 pp and FONLL calculations, we can go to the folder MakeFinalPlots/ComparisonPlots

cd MakeFinalPlots/ComparisonPlots

The FONLL pp cross section calculations for Bs and B+ are saved at the folder FONLL/

B+: FONLLs/output-Bp.root 

Bs: FONLLs/output-Bs.root 


To compare the B+ cross section results, simply do:

root -b -l -q BPComparison.C

The files are saved at 

Plots/BP/

The comparison with 2015 results is saved as Plots/BP/BP2015CompLog.png as png format and Plots/BP/BP2015CompLog.pdf pdf format in log scale

The comparison with FONLL results is saved as Plots/BP/BPFONLLCompLog.png as png format and Plots/BP/BPFONLLCompLog.pdf pdf format in log scale

The files contains the comparison and the ratio between 2017 pp and the references.

Likewise for Bs, just do 

root -b -l -q BsComparison.C

The files are saved at 

Plots/Bs/

The comparison with 2015 results is saved as Plots/Bs/Bs2015CompLog.png as png format and Plots/Bs/Bs2015CompLog.pdf pdf format in log scale

The comparison with FONLL results is saved as Plots/Bs/BsFONLLCompLog.png as png format and Plots/Bs/BsFONLLCompLog.pdf pdf format in log scale


At this point, we have obtained the comparison with 2015 pp and FONLL results. Note that here the 2D map correction method is commented out. We could uncomment them to also include the 2D map efficiency correction medthod, which will be presented in orange, if needed later. 


### RAA Comparison Plots


The comparison of RAA results with 2015 are documented in the RAA Section under  BsBPFinalResults/Comparisons/RAA/ folder. 


# Tests

In addition to systematic uncertainties, there are several tests we need to do in order to validate our results. The tests are documented at the appendix after the systematic uncertainties section. 

## Bs NonPrompt Studies


To check the leakage b-hadron decay modes in the Bs signal region, we need to conduct a studies of different b hadrons decay using the inclusive non-prompt J/psi samples. Here, we have prepared everything: the BDT values are all computed, to conduct the studies. The BDT is computed on the root tree, a Bfinder produce tree where the tracks are taken into fully reconstructing the Bs via the decay channel Bs -> J/psi phi -> mu+ mu- K+ K- are all saved. Here, we have computed the BDT values for all pT bins and apply the analysis cuts to the samples.

The codes are located at: Bs/Tests/ppNonPrompStudies

cd Bs/Tests/ppNonPrompStudies


The files with all possible decay models are saved at:

[szhaozho@GRENDEL01 ppNonPrompStudies]$ ls OutFile
BsNPStudies_0_100.root  BsNPStudies_1_2.root    BsNPStudies_20_30.root  BsNPStudies_30_50.root  BsNPStudies_5_7.root
BsNPStudies_10_15.root  BsNPStudies_15_20.root  BsNPStudies_2_3.root    BsNPStudies_3_5.root    BsNPStudies_7_10.root

Here we have the studies of each pT bin as well as the inclusive pT bin as seen above 

Inside each file, the signal Bs decay as well as the non-signal and a decomposition of the decay models are saved. For instance:

[szhaozho@GRENDEL01 OutFile]$ root -l BsNPStudies_15_20.root 
root [0] 
Attaching file BsNPStudies_15_20.root as _file0...
(TFile *) 0x2b3d6a0
root [1] .ls
TFile**		BsNPStudies_15_20.root	
 TFile*		BsNPStudies_15_20.root	
  KEY: TH1D	Bmass;1	Bmass
  KEY: TH1D	BmassBs;1	BmassBs
  KEY: TH1D	BmassBs_nosig;1	BmassBs_nosig
  KEY: TH1D	BmassBsPiK;1	BmassBsPiK
  KEY: TH1D	BmassBsKK;1	BmassBsKK
  KEY: TH1D	BmassBsKKfake_binfo;1	BmassBsKKfake_binfo
  KEY: TH1D	BmassBsXpipi;1	BmassBsXpipi
  KEY: TH1D	BmassBsKKX;1	BmassBsKKX
  KEY: TH1D	BmassB0KStar;1	BmassB0KStar
  KEY: TH1D	BmassB0Hard;1	BmassB0Hard
  KEY: TH1D	BmassB0Leak;1	BmassB0Leak

The details of each decay mode is documented in the appendix C of the analysis note. 

Here, to obtain the plots, simply run the code:

root -b -l -q PlotNPShapesNew.C


The plots will be produced and saved at NPNewPlots/ as:

NPBackground_0.png  NPBackground_3.png  NPBackground_6.png  NPBackgroundONLY_0.png  NPBackgroundONLY_3.png  NPBackgroundONLY_6.png
NPBackground_1.png  NPBackground_4.png  NPBackground_7.png  NPBackgroundONLY_1.png  NPBackgroundONLY_4.png  NPBackgroundONLY_7.png
NPBackground_2.png  NPBackground_5.png  NPBackground_8.png  NPBackgroundONLY_2.png  NPBackgroundONLY_5.png  NPBackgroundONLY_8.png

NPBackground_5.png: is for pT bin of 10 - 15 GeV/c with both signal decay channel: Bs -> J/psi phi -> mu+ mu- K+ K- a in red plot together with all other non-prompt background channels 
NPBackgroundONLY_5.png: is for  pT bin of 10 - 15 GeV/c with all other non-prompt background channels as well as the red one of the inclusive non-prompt background channels subtracting the K* and KK channel, where peaking ocruing within our signal region: 5.36 +- 0.08 GeV/c^2

I have also saved one before cut name as: NPNewPlotsNoCut

Here, we can see that after applying the phi meson veto cut, we have significantly suppressed the non-prompt background to the degree of < 10% compared to the results without cuts. 

## Technical Support

The is the very first version of the instruction. It is still far from being complete. If you have any question, please feel free to email me: zzshi@mit.edu or ping me on Skype and Slack.
