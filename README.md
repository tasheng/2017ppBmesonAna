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


# Systematic Studies 


## TnP Systematics


## MC/Data Discrepancy 


## Bpt Systematics


## MC Stat Systematics 


## Generating Summary Plots


## Report the Systematic Uncertanties to the Final Results


# Tests




## Bs NonPrompt Studies




## Technical Support

The is the very first version of the instruction. It is still far from being complete. If you have any question, please feel free to email me: zzshi@mit.edu or ping me on Skype and Slack.
