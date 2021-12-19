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

### RAA and Bs/B+

To Run the RAA and Bs/B+. It is also very simple. We could go to 

cd BsBPFinalResults/BsBPRatio

And then run 

root -b -l -q PlotsBsBPRatio.C'(0,0)'

The first argument: 0 stands for the ratio as a function of pT while 1 stands for multiplicity

The second argument: 0 stands for Binned pT efficiency correction while 1 stands for 2D map efficiency correction

The plots are stored at Pt/

If we run 

root -b -l -q BPRAA.C

The B+ RAA plots are stored at BPRAA/



## Generate Comparison Plots 

To run the comparison of the 2017 pp Bs and B+ results with 2015 pp and FONLL calculations, go to the folder:

Caveat - fiducial region: for the full 2015 pp results, since the measurement |By| > 2.4 for B pT < 10 GeV/c unlike the 2018 PbPb where a fiducial region |By| > 2.4 for B pT < 10 GeV/c. Therefore, we have produced two sets of analysis 

We can change the configuration to remove the fiducial region. Here we do not go through the details about that. I have produced two files one for Bs and one for B+ where the fiducial region is removed in order to compare with 2015 pp results. You can find the cross section files at:




### Without Fiducial Region 





# Systematic Studies 




## Technical Support

The is the very first version of the instruction. It is still far from being complete. If you have any question, please feel free to email me: zzshi@mit.edu or ping me on Skype and Slack.
