#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>
#include "BpSystValues.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"


void PlotSystBp(){

	//Global Bin//


	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);




	HisEmptyPt->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmptyPt->GetYaxis()->SetTitle("Systematic Uncertainties");
	HisEmptyPt->GetXaxis()->CenterTitle();
	HisEmptyPt->GetYaxis()->CenterTitle();
	HisEmptyPt->GetYaxis()->SetTitleOffset(1.4);


	TGraphAsymmErrors * GlobalSystGraph = new TGraphAsymmErrors(GlobalSystNBin,GlobalX,GlobalY,GlobalXHigh,GlobalXLow,GlobalSystHigh,GlobalSystLow);
	GlobalSystGraph->SetLineWidth(2);
	GlobalSystGraph->SetLineColor(2);
	GlobalSystGraph->SetFillColor(0);


	TGraphAsymmErrors * TotalSyst = new TGraphAsymmErrors(NPtBins,PtBin,PtBinY,PtBinWidthHigh,PtBinWidthLow,TotalSystHigh,TotalSystLow);
	TotalSyst->SetLineWidth(2);
	TotalSyst->SetLineColor(1);
	TotalSyst->SetFillColor(0);


	TGraphAsymmErrors * SignalSyst = new TGraphAsymmErrors(NPtBins,PtBin,PtBinY,PtBinWidthHigh,PtBinWidthLow,SignalSystHigh,SignalSystLow);
	SignalSyst->SetLineWidth(2);
	SignalSyst->SetLineColor(6);
	SignalSyst->SetFillColor(0);



	TGraphAsymmErrors * EffSyst = new TGraphAsymmErrors(NPtBins,PtBin,PtBinY,PtBinWidthHigh,PtBinWidthLow,EffSystHigh,EffSystLow);
	EffSyst->SetLineWidth(2);
	EffSyst->SetLineColor(4);
	EffSyst->SetFillColor(0);

	TGraphAsymmErrors * TnPSyst = new TGraphAsymmErrors(NPtBins,PtBin,PtBinY,PtBinWidthHigh,PtBinWidthLow,TnPSystHigh,TnPSystLow);
	TnPSyst->SetLineWidth(2);
	TnPSyst->SetLineColor(3);
	TnPSyst->SetFillColor(0);





	TCanvas *c = new TCanvas("c","c",600,600);
	c->cd();
	HisEmptyPt->Draw();
	GlobalSystGraph->Draw("5same");
	TotalSyst->Draw("5same");
	EffSyst->Draw("5same");
	TnPSyst->Draw("5same");
	SignalSyst->Draw("5same");


	TLegend * leg = new TLegend(0.20,0.15,0.70,0.40,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(GlobalSystGraph,"Global Systematics","l");
	leg->AddEntry(TotalSyst,"Total Systematics","l");
	leg->AddEntry(SignalSyst,"Signal Extraction Systematics","l");
	leg->AddEntry(EffSyst,"Efficiency Correction Systematics","l");
	leg->AddEntry(TnPSyst,"Tag and Probe Systematics","l");
	leg->Draw("SAME");

	TLine *linePt = new TLine(0,0,100,0);
	linePt->SetLineStyle(2);
	linePt->SetLineWidth(2);
	linePt->Draw("SAME");


	c->SaveAs("PlotSyst/Bp/BpPtSyst.png");
	c->SaveAs("PlotSyst/Bp/BpPtSyst.pdf");

	HisEmptyCent->GetXaxis()->SetTitle("Centrality (%)");
	HisEmptyCent->GetYaxis()->SetTitle("Systematic Uncertainties");
	HisEmptyCent->GetXaxis()->CenterTitle();
	HisEmptyCent->GetYaxis()->CenterTitle();
	HisEmptyCent->GetYaxis()->SetTitleOffset(1.4);



	TGraphAsymmErrors * GlobalSystGraphCent = new TGraphAsymmErrors(GlobalSystNBin,GlobalXCent,GlobalY,GlobalXCentHigh,GlobalXCentLow,GlobalSystHighCent,GlobalSystLowCent);
	GlobalSystGraphCent->SetLineWidth(2);
	GlobalSystGraphCent->SetLineColor(2);
	GlobalSystGraphCent->SetFillColor(0);


	TGraphAsymmErrors * TotalSystCent = new TGraphAsymmErrors(NCentBins,CentBin,CentBinY,CentBinWidthHigh,CentBinWidthLow,TotalSystHighCent,TotalSystLowCent);
	TotalSystCent->SetLineWidth(2);
	TotalSystCent->SetLineColor(1);
	TotalSystCent->SetFillColor(0);


	TGraphAsymmErrors * SignalSystCent = new TGraphAsymmErrors(NCentBins,CentBin,CentBinY,CentBinWidthHigh,CentBinWidthLow,SignalSystHighCent,SignalSystLowCent);
	SignalSystCent->SetLineWidth(2);
	SignalSystCent->SetLineColor(6);
	SignalSystCent->SetFillColor(0);



	TGraphAsymmErrors * EffSystCent = new TGraphAsymmErrors(NPtBins,CentBin,CentBinY,CentBinWidthHigh,CentBinWidthLow,EffSystHighCent,EffSystLowCent);
	EffSystCent->SetLineWidth(2);
	EffSystCent->SetLineColor(4);
	EffSystCent->SetFillColor(0);

	TGraphAsymmErrors * TnPSystCent = new TGraphAsymmErrors(NPtBins,CentBin,CentBinY,CentBinWidthHigh,CentBinWidthLow,TnPSystHighCent,TnPSystLowCent);
	TnPSystCent->SetLineWidth(2);
	TnPSystCent->SetLineColor(3);
	TnPSystCent->SetFillColor(0);



	TCanvas *cCent = new TCanvas("cCent","cCent",600,600);
	cCent->cd();
	HisEmptyCent->Draw();
	GlobalSystGraphCent->Draw("5same");
	TotalSystCent->Draw("5same");
	EffSystCent->Draw("5same");
	TnPSystCent->Draw("5same");
	SignalSystCent->Draw("5same");
	

	TLegend * legCent = new TLegend(0.20,0.15,0.70,0.40,NULL,"brNDC");
	legCent->SetBorderSize(0);
	legCent->SetTextSize(0.04);
	legCent->SetTextFont(42);
	legCent->SetFillStyle(0);
	legCent->AddEntry(GlobalSystGraph,"Global Systematics","l");
	legCent->AddEntry(TotalSyst,"Total Systematics","l");
	legCent->AddEntry(SignalSyst,"Signal Extraction Systematics","l");
	legCent->AddEntry(EffSyst,"Efficiency Correction Systematics","l");
	legCent->AddEntry(TnPSyst,"Tag and Probe Systematics","l");
	legCent->Draw("SAME");


	TLine *lineCent = new TLine(-10,0,120,0);
	lineCent->SetLineStyle(2);
	lineCent->SetLineWidth(2);
	lineCent->Draw("SAME");



	cCent->SaveAs("PlotSyst/Bp/BpCentSyst.png");
	cCent->SaveAs("PlotSyst/Bp/BpCentSyst.pdf");



	HisEmptyInc->GetXaxis()->SetTitle("Centrality (%)");
	HisEmptyInc->GetYaxis()->SetTitle("Systematic Uncertainties");
	HisEmptyInc->GetXaxis()->CenterTitle();
	HisEmptyInc->GetYaxis()->CenterTitle();
	HisEmptyInc->GetYaxis()->SetTitleOffset(1.4);



	TGraphAsymmErrors * GlobalSystGraphInc = new TGraphAsymmErrors(GlobalSystNBin,GlobalXCent,GlobalY,GlobalXCentHigh,GlobalXCentLow,GlobalSystHighInc,GlobalSystLowInc);
	GlobalSystGraphInc->SetLineWidth(2);
	GlobalSystGraphInc->SetLineColor(2);
	GlobalSystGraphInc->SetFillColor(0);


	TGraphAsymmErrors * TotalSystInc = new TGraphAsymmErrors(NInclusive,InclusiveBin,InclusiveBinY,InclusiveBinWidthHigh,InclusiveBinWidthLow,TotalSystHighInc,TotalSystLowInc);
	TotalSystInc->SetLineWidth(2);
	TotalSystInc->SetLineColor(1);
	TotalSystInc->SetFillColor(0);


	TGraphAsymmErrors * SignalSystInc = new TGraphAsymmErrors(NInclusive,InclusiveBin,InclusiveBinY,InclusiveBinWidthHigh,InclusiveBinWidthLow,SignalSystHighInc,SignalSystLowInc);
	SignalSystInc->SetLineWidth(2);
	SignalSystInc->SetLineColor(6);
	SignalSystInc->SetFillColor(0);



	TGraphAsymmErrors * EffSystInc = new TGraphAsymmErrors(NInclusive,InclusiveBin,InclusiveBinY,InclusiveBinWidthHigh,InclusiveBinWidthLow,EffSystHighInc,EffSystLowInc);
	EffSystInc->SetLineWidth(2);
	EffSystInc->SetLineColor(4);
	EffSystInc->SetFillColor(0);

	TGraphAsymmErrors * TnPSystInc = new TGraphAsymmErrors(NInclusive,InclusiveBin,InclusiveBinY,InclusiveBinWidthHigh,InclusiveBinWidthLow,TnPSystHighInc,TnPSystLowInc);
	TnPSystInc->SetLineWidth(2);
	TnPSystInc->SetLineColor(3);
	TnPSystInc->SetFillColor(0);



	TCanvas *cInc = new TCanvas("cInc","cInc",600,600);
	cInc->cd();
	HisEmptyInc->Draw();
	GlobalSystGraphInc->Draw("5same");
	TotalSystInc->Draw("5same");
	EffSystInc->Draw("5same");
	TnPSystInc->Draw("5same");
	SignalSystInc->Draw("5same");
	

	TLegend * legInc = new TLegend(0.20,0.15,0.70,0.40,NULL,"brNDC");
	legInc->SetBorderSize(0);
	legInc->SetTextSize(0.04);
	legInc->SetTextFont(42);
	legInc->SetFillStyle(0);
	legInc->AddEntry(GlobalSystGraph,"Global Systematics","l");
	legInc->AddEntry(TotalSyst,"Total Systematics","l");
	legInc->AddEntry(SignalSyst,"Signal Extraction Systematics","l");
	legInc->AddEntry(EffSyst,"Efficiency Correction Systematics","l");
	legInc->AddEntry(TnPSyst,"Tag and Probe Systematics","l");
	legInc->Draw("SAME");


	TLine *lineInc = new TLine(-10,0,120,0);
	lineInc->SetLineStyle(2);
	lineInc->SetLineWidth(2);
	lineInc->Draw("SAME");



	cInc->SaveAs("PlotSyst/Bp/BpIncSyst.png");
	cInc->SaveAs("PlotSyst/Bp/BpIncSyst.pdf");







	TCanvas *cComb = new TCanvas("cComb","cComb",1800,600);
	cComb->Divide(3,1);
	cComb->cd(1);
	HisEmptyPt->Draw();
	GlobalSystGraph->Draw("5same");
	TotalSyst->Draw("5same");
	EffSyst->Draw("5same");
	TnPSyst->Draw("5same");
	SignalSyst->Draw("5same");
	linePt->Draw("SAME");
	leg->Draw("SAME");


	cComb->cd(2);
	HisEmptyCent->Draw();
	GlobalSystGraphCent->Draw("5same");
	TotalSystCent->Draw("5same");
	EffSystCent->Draw("5same");
	TnPSystCent->Draw("5same");
	SignalSystCent->Draw("5same");
	lineCent->Draw("SAME");
	legCent->Draw("SAME");

	cComb->cd(3);
	HisEmptyInc->Draw();
	GlobalSystGraphInc->Draw("5same");
	TotalSystInc->Draw("5same");
	EffSystInc->Draw("5same");
	TnPSystInc->Draw("5same");
	SignalSystInc->Draw("5same");
	lineInc->Draw("SAME");
	legInc->Draw("SAME");
	cComb->SaveAs("PlotSyst/Bp/BpSysSumPlot.png");
	cComb->SaveAs("PlotSyst/Bp/BpSysSumPlot.pdf");


	
}
