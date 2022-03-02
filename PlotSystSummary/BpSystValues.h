const int GlobalSystNBin = 1;
double GlobalX[GlobalSystNBin] = {3};
double GlobalY[GlobalSystNBin] = {0};

double GlobalXLow[GlobalSystNBin] = {1.5};
double GlobalXHigh[GlobalSystNBin] = {1.5};
double GlobalSystLow[GlobalSystNBin] = {0.0385};
double GlobalSystHigh[GlobalSystNBin] = {0.0385};




const int NPtBins = 7;
double PtBin[NPtBins] ={6,8.5,12.5,17.5,25,40,55};
double PtBinY[NPtBins] ={0,0,0,0,0,0,0};
double PtBinWidthHigh[NPtBins] ={1,1.5,2.5,2.5,5,10,5};
double PtBinWidthLow[NPtBins] ={1,1.5,2.5,2.5,5,10,5};


double TotalSystHigh[NPtBins] = {0.1404,0.1714,0.0775,0.0715};
double TotalSystLow[NPtBins] = {0.1359,0.1705,0.0761,0.0699};

double SignalSystHigh[NPtBins] = {0.0446153,0.0272526,0.0279989,0.0257498};
double SignalSystLow[NPtBins] = {0.0446153,0.0272526,0.0279989,0.0257498};

double EffSystHigh[NPtBins] ={0.111905,0.163718,0.0613144,0.0543553};
double EffSystLow[NPtBins] = {0.111905,0.163718,0.0613144,0.0543553};


double TnPSystHigh[NPtBins] = {0.0721,0.0392,0.0383,0.0387};
double TnPSystLow[NPtBins] = {0.0628,0.0323,0.0353,0.0356};



const int NCentBins = 10;
double CentBin[NCentBins] ={7.5,20,27.5,32.5,37.5,45,57.5,72.5,90,115};
double CentBinY[NCentBins] ={0,0,0,0,0,0,0,0,0,0};
double CentBinWidthHigh[NCentBins] ={7.5,5,2.5,2.5,2.5,5,7.5,7.5,10,15};
double CentBinWidthLow[NCentBins] ={7.5,5,2.5,2.5,2.5,5,7.5,7.5,10,15};


double GlobalXCent[GlobalSystNBin] = {-5};
double GlobalXCentLow[GlobalSystNBin] = {2.5};
double GlobalXCentHigh[GlobalSystNBin] = {2.5};

double GlobalSystLowCent[GlobalSystNBin] = {0.0292};
double GlobalSystHighCent[GlobalSystNBin] = {0.0292};


double TotalSystHighCent[NCentBins] = {0.1553,0.1189};
double TotalSystLowCent[NCentBins] = {0.1543,0.1177};

double SignalSystHighCent[NCentBins] = {0.0253372,0.0280745};
double SignalSystLowCent[NCentBins] = {0.0253372,0.0280745};

double EffSystHighCent[NCentBins] ={0.145451,0.100873};
double EffSystLowCent[NCentBins] = {0.145451,0.100873};


double TnPSystHighCent[NCentBins] = {0.0418,0.0414};
double TnPSystLowCent[NCentBins] = {0.0383,0.0380};














const int NInclusive = 1;
double InclusiveBin[NInclusive] ={45};
double InclusiveBinY[NInclusive] ={0};
double InclusiveBinWidthLow[NInclusive] ={45};
double InclusiveBinWidthHigh[NInclusive] ={45};


double GlobalSystLowInc[GlobalSystNBin] = {0.0292};
double GlobalSystHighInc[GlobalSystNBin] = {0.0292};

double TotalSystHighInc[NInclusive] = {0.1392};
double TotalSystLowInc[NInclusive] = {0.1382};

double SignalSystHighInc[NInclusive] = {0.0263483};
double SignalSystLowInc[NInclusive] = {0.0263483};

double EffSystHighInc[NInclusive] ={0.12768};
double EffSystLowInc[NInclusive] = {0.12768};


double TnPSystHighInc[NInclusive] = {0.0416};
double TnPSystLowInc[NInclusive] = {0.0381};





TH2D * HisEmptyPt = new TH2D("HisEmptyPt","",100,0,100,100,-1,1);
TH2D * HisEmptyCent = new TH2D("HisEmptyCent","",105,-10,120,100,-1,1);
TH2D * HisEmptyInc = new TH2D("HisEmptyInc","",105,-10,120,100,-1,1);

