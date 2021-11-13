# CMS 2017 Bs and B+ Analysis Codes

## Basic Organization of the Analysis Codes

These codes are made to analyze the CMS 2017 pp datasets and obtain the cross section of the B+ and Bs. Therefore, the organization of the codes is BP, Bs, and BsBP. They consist of raw yield extraction, efficiency correction, and cross section calculations. There are also folders named SkimmedSamples and UnskimmedSamples. These are the folders where you put your input files to perform the analysis

To obtain the codes, do:

git clone https://github.com/MYOMAO/2017ppBmesonAna.git


## Input Data and MC Samples

Here we assume you already have the data and official MC samples. You could find the samples on lxplus EOS:

/eos/cms/store/group/phys_heavyions/zshi/ForHenrique/

You should copy the entire folder to your working directory

if you do ls /eos/cms/store/group/phys_heavyions/zshi/ForHenrique/, you will see two folders

Skimmed  UnSkimmed

So basically, copy the files in Skimmed to SkimmedSamples and UnSkimmed to UnskimmedSamples. Then you will have all the files in the right place and ready to perform the analysis.


## Raw Yield Fit

The analysis consist of two parts: Raw yield extraction and efficiency correction. To perform the signal raw yield extraction, we will use the skimmed files. To run the codes, for example, for B+, simply do:

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


The procedure to obtain raw yield in Bs is basically the same as BP.

## Efficiency Correction
The next step is to obtain the efficiency correction. The steps are also quite straightforward. To Run the efficiency simply do:

cd BP/EffAna

root -b -l -q MCEff.C

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

These codes run on the unskimmed MC files (not flattened): 

Once you finishing running the codes, the plots are saved in the root file 

EffFineBDT.root

The 1D and 2D efficiency histograms are saved inside the file

Again, Bs efficiency follows the same procedures

## Cross Section Calculations

The cross section calculations are also straightforward. Since we have obtained the raw yield and efficiency, to obtain the cross section, we simply run 

cd BP/EffAna/

For pT, simply run

root -b -l -q CrossSectionAna.C


For Multiplicity, simply run

root -b -l -q CrossSectionAnaMult.C

The output file of the cross section is saved at FinalFile/

FinalFile/BPPPCorrYieldPT.root and FinalFile/BPPPCorrYieldMult.root

Here, when we run the cross section, we need to make sure that the bin is indeed the same. For example, in here, NPtBins = 4, so in the CrossSectionAna.C, we make sure that 

const int NBins = 4;

Here, we again use two efficiency correction methods. The binned pT one is at 

TH1D * CorrDiffHisBin = (TH1D * ) CorrDiffHis ->Clone("CorrDiffHisBin");

The 2D map one is at:

TH1D * CorrDiffHis = new TH1D("hPtSigma","",NBins,ptBins);

Repeat this for Bs, we could get the cross section for Bs

At this point we get the Bs and BP cross sections


## RAA and Bs/B+

To Run the RAA and Bs/B+. It is also very simple. We could go to 

cd BsBP/

And run 

root -b -l -q BsBP.C to get the BsBP.C ratio in Pt

The plots are stored at Pt/

If we run 

root -b -l -q BPRAA.C

The B+ RAA plots are stored at BPRAA/

root -b -l -q BsBPMult.C to get the BsBP.C ratio in Multiplicity 

The plots are stored at Mult/

## Technical Support

The is the very first version of the instruction. It is still far from being complete. If you have any question, please feel free to email me: zzshi@mit.edu or ping me on Skype and Slack.
