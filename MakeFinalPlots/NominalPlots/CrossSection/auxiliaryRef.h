#ifndef AUXILIARYREF
#define AUXILIARYREF


// fs/fu
// PDG
double BsFrag = 10.3;
double BsErr  = 0.5;
double BPFrag = 40.5;
double BPErr  = 0.6;
const int BandBin = 1;
//double BandY[BandBin] = {BsFrag/BPFrag};
//double BandY[BandBin] = {0.256};

double BandY[BandBin] = {0.244};


//double bin_min = 5;
//double bin_max = 50;
double bin_min = 4;
double bin_max = 25;




double BandX[BandBin] = {(bin_max+bin_min)/2};
double BandXErr[BandBin] = {(bin_max-bin_min)/2};


//forCent
double binCent_min = 0;
double binCent_max = 450;
double BandXCent[BandBin] = {(binCent_max+binCent_min)/2};
double BandXErrCent[BandBin] = {(binCent_max-binCent_min)/2};

//double BandYErr[BandBin] = { BsFrag/BPFrag*TMath::Sqrt((BsErr/BsFrag)*(BsErr/BsFrag)+(BPErr/BPFrag)*(BPErr/BPFrag)) };
double BandYErr[BandBin] = { 0.012 };



double BPBsFFR = 1.93;
double BPBsFFRErr = 0.16;

const int LHCb7TeVNPoints = 12;
double LHCb7TeVX[LHCb7TeVNPoints] = {1.25,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.73,12.555681053847229,18.378237750168154};
double LHCb7TeVXErrUp[LHCb7TeVNPoints] = {0.75,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.27,1.44,21.62};
double LHCb7TeVXErrDown[LHCb7TeVNPoints] = {0.75,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.73,0.56,4.38};

double LHCb7TeVY[LHCb7TeVNPoints] = {0.2262,0.2414,0.2281,0.2390,0.2390,0.2433,0.2452,0.2452,0.2186,0.2205,0.2243,0.2239};
double LHCb7TeVYErr[LHCb7TeVNPoints+1] = {0.0057,0.0057,0.0038,0.0038,0.0038,0.0057,0.0057,0.0057,0.0057,0.0057,0.0057};


#endif
