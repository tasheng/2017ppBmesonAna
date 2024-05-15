#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"
#include <iostream>

//
// Global variables
//
TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

// bool writeExtraText = true;//false;
bool writeExtraText = false;//false;
TString extraText   = "Preliminary";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.40 * (1 + 0.20);
float lumiTextOffset   = 0.10;
float cmsTextSize      = 0.55;
float cmsTextOffset    = 0.01;  // in case of ticks

float relPosX    = 0.045;
float relPosY    = 0.045;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "20.1 fb^{-1}";
TString lumi_8TeV  = "19.7 fb^{-1}";
TString lumi_7TeV  = " (1.7 nb^{-1})";
TString lumi_5TeV  = " (302.3 pb^{-1})";




bool drawLogo      = false;

void CMS_lumi( TPad* pad, int iPeriod=3, int iPosX=10 );
