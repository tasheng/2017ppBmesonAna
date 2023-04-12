#include <vector>
#include <array>
#include <map>

// PbPb published binning
/* const std::vector<double> ptbinsvec = { */
/* 	/\* 5, 7, 10, 15, 20, 50, 60 *\/ */
/* 	7, 10, 15, 20, 50 */
/* }; */

const unsigned nptBins = 4;
const std::array<double, nptBins + 1> ptbinsvec = {
	7, 10, 15, 20, 50
};

// // const unsigned nptBins = 1;
// const unsigned nptBins = 2;
// const std::array<double, nptBins + 1> ptbinsvec = {
// 	5, 7, 10
// 	// 7, 10
// };
const std::vector<float> abscissae = {
8.5, 12.5, 17.5, 35
};

enum Tracking{
  loose = 0,
  standard = 1,
  tight = 2,
};

std::map<Tracking, double> ptErr{
  {Tracking::standard, 0.1},
  {Tracking::tight, 0.05},
  {Tracking::loose, 0.15}
};

std::map<Tracking, double> chi2Nlayer{
  {Tracking::standard, 0.18},
  {Tracking::tight, 0.15},
  {Tracking::loose, 0.18}
};

const double epsilon = 1e-7;

// 7 bins
/* float BPXsecPbPbY[NBins] = {4.82132e+06/11.1,311668,270167,64384.4,208537/11.1,28700.6/11.1,7000.73/11.1}; */
/* float BPXsecPbPbX[NBins] = {6,8.73,12.4,17.2,25,40,55}; */
/* float BPXSecPbPbXErrUp[NBins] = {1,1.27,2.6,2.8,5,10,5}; */
/* float BPXSecPbPbXErrDown[NBins] = {1,1.23,2.4,2.2,5,10,5}; */
/* float BPXSecPbPbYErrUpRatio[NBins] = {0.278198,0.159,0.041,0.0654,0.0690334,0.104543,0.24575}; */
/* float BPXSecPbPbYErrDownRatio[NBins] = {0.278198,0.145,0.0795,0.065,0.0690334,0.104543,0.24575}; */

const unsigned nptBinsPbPb = 4;

//BP_PBPB
float BPXsecPbPbY[nptBinsPbPb] = {311668, 270167, 64384.4, 7704};
float BPXsecPbPbX[nptBinsPbPb] = {8.73,12.4,17.2,27.3};
float BPXSecPbPbXErrUp[nptBinsPbPb] = {1.27,2.6,2.8,22.7};
float BPXSecPbPbXErrDown[nptBinsPbPb] = {1.23,2.4,2.2,7.3};

float BPXSecPbPbYErrUpRatio[nptBinsPbPb] =   {0.159, 0.041, 0.0654, 0.0691};
float BPXSecPbPbYErrDownRatio[nptBinsPbPb] = {0.145, 0.0795, 0.065, 0.0526};

float BPXSecPbPbYSystUpRatio[nptBinsPbPb]   = {0.1404, 0.1714, 0.0775, 0.0715};
float BPXSecPbPbYSystDownRatio[nptBinsPbPb] = {0.1359, 0.1705, 0.0761, 0.0698};


//Bs_PBPB
float BsXsecPbPbY[nptBinsPbPb] = {160432,75523.7,25354.5,2272.18};	
float BsXsecPbPbX[nptBinsPbPb] = {8.75,12.6,17.4,27.3};
float BsXSecPbPbXErrUp[nptBinsPbPb] = {1.25,2.4,2.4,22.7};
float BsXSecPbPbXErrDown[nptBinsPbPb] = {1.25,2.6,2.6,7.3};

float BsXSecPbPbYErrUpPercent[nptBinsPbPb] = {0.513,0.224,0.216,0.216};
float BsXSecPbPbYErrDownPercent[nptBinsPbPb] = {0.483,0.256,0.207,0.163};

float BsXSecPbPbYSystUpPercent[nptBinsPbPb] = {0.4564,0.1482,0.1218,0.1647};
float BsXSecPbPbYSystDownPercent[nptBinsPbPb] = {0.4564,0.1454,0.1210,0.1640};






// 2015 data


	float BsRAAY2015[2] = {1.51,0.87};
	float BsRAAX2015[2] = {11.0,32.5};
	float BsRAAXErrUp2015[2] = {4.0,17.5};
	float BsRAAXErrDown2015[2] = {4.0,17.5};
	float BsRAAYErrUp2015[2] = {0.61,0.30};
	float BsRAAYErrDown2015[2] = {0.61,0.30};
	float BsRAAYSystUp2015[2] = {0.50,0.17};
	float BsRAAYSystDown2015[2] = {0.50,0.17};

  const int NBins2015 = 5;
	float BPRAAY2015[NBins2015] = {0.35,0.448,0.440,0.615,0.35};
	float BPRAAX2015[NBins2015] = {8.5,12.5,17.5,25,40};
	float BPRAAXErrUp2015[NBins2015] = {1.5,2.5,2.5,5,10};
	float BPRAAXErrDown2015[NBins2015] = {1.5,2.5,2.5,5,10};
	float BPRAAYErrUp2015[NBins2015] = {0.11,0.074,0.075,0.092,0.11};
	float BPRAAYErrDown2015[NBins2015] = {0.11,0.074,0.075,0.092,0.11};
	float BPRAAYSystUp2015[NBins2015] = {0.064,0.077,0.074,0.102,0.0539};
	float BPRAAYSystDown2015[NBins2015] = {0.064,0.077,0.074,0.102,0.0539};