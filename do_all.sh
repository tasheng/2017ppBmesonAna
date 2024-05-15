#! /bin/bash

sync_with_main () {
    mkdir -p BP/EffAna/BDTWeights/
    mkdir -p Bs/EffAna/BDTWeights/
    rsync -a ../2017ppBmesonAna/BP/EffAna/BDTWeights/BPw.root BP/EffAna/BDTWeights/BPw.root
    rsync -a ../2017ppBmesonAna/Bs/EffAna/BDTWeights/Bsw.root Bs/EffAna/BDTWeights/Bsw.root
    rsync -a ../2017ppBmesonAna/NewBptStudies/ResultFile/BptWeight_BP.root NewBptStudies/ResultFile/BptWeight_BP.root
    rsync -a ../2017ppBmesonAna/NewBptStudies/ResultFile/BptWeight_Bs.root NewBptStudies/ResultFile/BptWeight_Bs.root
}

maketnp () {
    ## run TnP and attach them to the unskimmed MC files
    pushd MakeMCTnP/
    root -q -b -l TnPWeight.C'(0)' &
    root -q -b -l TnPWeight.C'(1)' &
    wait
    popd
}

# new bpt reweighting uses FONLL and doesn't depend on acceptance
# function >> MCEff.C
bptshape () {
    pushd NewBptStudies
    python ReweightY.py &
    wait
    root -q -b -l ReweightBpt.C'(0)' &
    root -q -b -l ReweightBpt.C'(1)' &
    wait
}

yield () {
    ## yield extraction
    pushd henri2022
    # roofitB.C contains 2 By cuts
    bash bpdoRoofit.sh &
    bash bsdoRoofit.sh &
    wait
    popd
}

yield () {
    ## yield extraction
    pushd henri2022
    # roofitB.C contains 2 By cuts
    bash bpdoRoofit.sh &
    bash bsdoRoofit.sh &
    wait
    popd
}

bpEff () {
    ## eff correction (2d map)
    # this is everything from yield to syst, other than MC syst
    pushd BP/EffAna
    # about 1hr
    # takes unskimmed MC sample as input
    echo "Takes BPw.root as input. Make sure it is up to date"
    ls -l BDTWeights/BPw.root
    # root -b -l -q MCEff.C+'(1,0)' 
    # CrossSectionAna.C contains 1 By cuts
    root -b -l -q CrossSectionAna.C'(1)'
    # root -b -l -q CrossSectionAnaMult.C'(1)'
    # >> BP/EffAna/FinalFiles/BPPPCorrYieldPT.root
    popd
}


bsEff () {
    pushd Bs/EffAna
    # about 1hr
    echo "Takes Bsw.root as input"
    ls -l BDTWeights/Bsw.root
    # root -b -l -q MCEff.C+'(1,0)' 

    root -b -l -q CrossSectionAna.C'(1)'
    # root -b -l -q CrossSectionAnaMult.C'(1)'
    # >> Bs/EffAna/FinalFiles/BsPPCorrYieldPT.root
    popd
}


nominal () {
    # central value
    cd BsBPFinalResults/BsBPRatio
    root -b -l -q PlotBsBPRatio.C'(1)'
    cd ../RAA
    root -b -l -q BPRAA.C
    root -b -l -q BsRAA.C
    cd ..

    cd Comparisons/Fiducial
    root -b -l -q BPComparison.C
    root -b -l -q BsComparison.C
    # root -b -l -q BPNewFidNoScale.C
    # root -b -l -q BsNewFidNoScale.C
    cd ../../..
}

syst () {
    pushd 2DMapSyst
    root -b -l -q CalEffSystBP.C # >> outfiles/bpsyst2d.root
    root -b -l -q CalEffSystBs.C
    root -b -l -q PlotEffSyst2D.C'(0)'
    root -b -l -q PlotEffSyst2D.C'(1)'
    # need to pass the above to cross calculation
    popd
}

## MC Stat Systematics
# takes 2D map eff as input
bpStat () {
    pushd MCStatSyst/BP
    root -b -l -q Generate2DMaps.C
    # more than 2hr
    root -b -l -q MCStatCal.C > mcstat.log
    popd
}

bsStat() {
    cd MCStatSyst/Bs
    root -b -l -q Generate2DMaps.C
    # ~1hr
    root -b -l -q MCStatCal.C > mcstat.log
    cd ../..
}

comp () {
    # get pdf variation errors
    python master.py

	    # Get pre-selection error
	#     python comppre.py
    python preeff.py |& tee presel_table.txt
 
   # comparison plot again
    pushd BsBPFinalResults/Comparisons/Fiducial/
    # << BP/EffAna/FinalFiles/BPPPCorrYieldPT.root
  
    root -b -l -q BPComparison.C
    root -b -l -q BsComparison.C

   python syst_table.py
    cd ../../RAA/

    root -b -l -q BPRAA.C
    root -b -l -q BsRAA.C
    popd
}

paperPlots () {
    # input
    pushd MakeFinalPlots/NominalPlots/CrossSection
    root -b -l -q plotPt.C'(1,1,0,1,1)'
    cd ../RAA
    root -b -l -q plotPt.C'(1,1,0,1,1,0)'
    root -b -l -q plotPt.C'(1,1,0,1,1,1)'
    root -b -l -q plotPt.C'(1,1,0,1,0,0)'
    root -b -l -q plotPt.C'(1,1,0,1,0,1)'
    root -b -l -q plotPt.C'(1,1,0,1,1,1,1)'
    popd
}
# maketnp

#bptshape



#yield
# wait

#sync_with_main

bpEff &
bsEff &
wait


syst
nominal


# bpStat&
# bsStat&
# wait

comp
paperPlots
