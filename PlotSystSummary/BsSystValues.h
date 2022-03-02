const int GlobalSystNBin = 1;
double GlobalX[GlobalSystNBin] = {3};
double GlobalY[GlobalSystNBin] = {0};

double GlobalXLow[GlobalSystNBin] = {1.5};
double GlobalXHigh[GlobalSystNBin] = {1.5};
double GlobalSystLow[GlobalSystNBin] = {0.0792};
double GlobalSystHigh[GlobalSystNBin] = {0.0792};




const int NPtBins = 4;
double PtBin[NPtBins] ={8.5,12.5,17.5,35};
double PtBinY[NPtBins] ={0,0,0,0};
double PtBinWidthHigh[NPtBins] ={1.5,2.5,2.5,15};
double PtBinWidthLow[NPtBins] ={1.5,2.5,2.5,15};


double TotalSystHigh[NPtBins] = {0.4564,0.1482,0.1218,0.1647};
double TotalSystLow[NPtBins] = {0.4564,0.1454,0.1210,0.1640};

double SignalSystHigh[NPtBins] = {0.0122566,0.0356535,0.0177282,0.0639269};
double SignalSystLow[NPtBins] = {0.0122566,0.0356535,0.0177282,0.0639269};

double EffSystHigh[NPtBins] ={0.447535,0.130863,0.114543,0.146684};
double EffSystLow[NPtBins] = {0.447535,0.130863,0.114543,0.146684};


double TnPSystHigh[NPtBins] = {0.0888,0.0597,0.0373,0.0391};
double TnPSystLow[NPtBins] = {0.0749,0.0524,0.0346,0.0360};



const int NCentBins = 2;
double CentBin[NCentBins] ={15,60};
double CentBinY[NCentBins] ={0,0};
double CentBinWidthHigh[NCentBins] ={15,30};
double CentBinWidthLow[NCentBins] ={15,30};


double GlobalXCent[GlobalSystNBin] = {-5};
double GlobalXCentLow[GlobalSystNBin] = {2.5};
double GlobalXCentHigh[GlobalSystNBin] = {2.5};

double GlobalSystLowCent[GlobalSystNBin] = {0.0750};
double GlobalSystHighCent[GlobalSystNBin] = {0.0750};


double TotalSystHighCent[NCentBins] = {0.1368,0.1271};
double TotalSystLowCent[NCentBins] = {0.1352,0.1260};

double SignalSystHighCent[NCentBins] = {0.0252517,0.0324507};
double SignalSystLowCent[NCentBins] = {0.0252517,0.0324507};

double EffSystHighCent[NCentBins] ={0.123685,0.10888};
double EffSystLowCent[NCentBins] = {0.123685,0.10888};


double TnPSystHighCent[NCentBins] = {0.0552,0.0463};
double TnPSystLowCent[NCentBins] = {0.0485,0.0419};














const int NInclusive = 1;
double InclusiveBin[NInclusive] ={45};
double InclusiveBinY[NInclusive] ={0};
double InclusiveBinWidthLow[NInclusive] ={45};
double InclusiveBinWidthHigh[NInclusive] ={45};


double GlobalSystLowInc[GlobalSystNBin] = {0.0750};
double GlobalSystHighInc[GlobalSystNBin] = {0.0750};

double TotalSystHighInc[NInclusive] = {0.1300};
double TotalSystLowInc[NInclusive] = {0.1277};

double SignalSystHighInc[NInclusive] = {0.0228449};
double SignalSystLowInc[NInclusive] = {0.0228449};

double EffSystHighInc[NInclusive] ={0.113699};
double EffSystLowInc[NInclusive] = {0.113699};


double TnPSystHighInc[NInclusive] = {0.0529};
double TnPSystLowInc[NInclusive] = {0.0470};





TH2D * HisEmptyPt = new TH2D("HisEmptyPt","",100,0,100,100,-1,1);
TH2D * HisEmptyCent = new TH2D("HisEmptyCent","",105,-10,120,100,-1,1);
TH2D * HisEmptyInc = new TH2D("HisEmptyInc","",105,-10,120,100,-1,1);

