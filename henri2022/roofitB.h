#include "TAxis.h"
#include "uti.h"
#include "TSystem.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooHist.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include <RooBifurGauss.h>
#include <RooCmdArg.h>
#include <fstream>
#include <string>
#include <iomanip>
#include "RooMCStudy.h"
#include <RooMinuit.h>

void plot_jpsifit(RooWorkspace& w, RooAbsPdf* model, RooDataSet* ds, TString plotName,  bool with_sig);
void fix_parameters(RooWorkspace& w, TString pdfName, bool release=false);
void fit_jpsinp (RooWorkspace& w, int pti, int ptf, bool includeSignal=true);

// draw legend and suppress parameters
const bool drawLegend = false;
using namespace RooFit;
using namespace std;

#define BS_MASS 5.36682
#define BP_MASS 5.27915

float bkgd;
double Significance;
double real_significance;

double minhisto=5.;
double maxhisto=6.;
int nbinsmasshisto=50;

Int_t _count=0;
RooWorkspace* outputw = new RooWorkspace("w");
RooWorkspace* w_val= new RooWorkspace("w_vl");


RooFitResult *fit(TString variation, TString pdf,TString tree, TCanvas* c, TCanvas* cMC, RooDataSet* ds, RooDataSet* dsMC, RooDataHist* dh, RooRealVar* mass, RooPlot* &outframe, int ptmin, int ptmax, int isMC, TString npfit, RooWorkspace& w)
{
	cout<<"total data: "<<ds->numEntries()<<endl;
	TH1* h = dh->createHistogram("Bmass");
	h->Sumw2(kFALSE);
	h->SetBinErrorOption(TH1::kPoisson);
	h->SetMarkerSize(1);
	h->SetMarkerStyle(20);
	h->SetLineColor(1);
	h->SetLineWidth(2);
	RooPlot* frameMC = mass->frame();
	frameMC->SetTitle("");
	if(tree=="ntKp")frameMC->SetXTitle("m_{J/#psiK^{#pm}} (GeV/c^{2})");
	if(tree=="ntphi")frameMC->SetXTitle("m_{J/#psiK^{+}K^{-}} (GeV/c^{2})");

	frameMC->SetYTitle("Events / (20 MeV/c^{2})");
	frameMC->GetXaxis()->CenterTitle();
	frameMC->GetXaxis()->SetTitleOffset(1.0);
	frameMC->GetYaxis()->SetTitleOffset(2.);
	frameMC->GetXaxis()->SetTitleSize(0.055);
	frameMC->GetYaxis()->SetTitleSize(0.035);
	frameMC->GetXaxis()->SetTitleFont(42);
	frameMC->GetYaxis()->SetTitleFont(42);
	frameMC->GetXaxis()->SetLabelFont(42);
	frameMC->GetYaxis()->SetLabelFont(42);
	frameMC->GetXaxis()->SetLabelSize(0.035);
	frameMC->GetYaxis()->SetLabelSize(0.035);
	frameMC->SetStats(0);
	frameMC->GetXaxis()->SetNdivisions(-50205);

	// give better initial values to the Bmesons mass
	double init_mean = BS_MASS;
	double m1 = 5.35;
	double m2 = 5.38;
	if (tree == "ntKp"){
		init_mean = BP_MASS;
		m1 = 5.27;
		m2 = 5.29;}
	// give better initial values to the Bmesons mass

	RooRealVar meanMC(Form("meanMC%d_%s",_count,pdf.Data()),"",init_mean,5.25,5.4) ;
	RooRealVar sigma1MC(Form("sigma1MC%d_%s",_count,pdf.Data()),"",0.05,0.01,0.11) ;
	RooRealVar sigma2MC(Form("sigma2MC%d_%s",_count, pdf.Data()),"",0.03,0.01,0.06) ;
	RooRealVar sigma3MC(Form("sigma3MC%d_%s",_count, pdf.Data()),"",0.01,0.005,0.02) ;
	RooRealVar sigma4cbMC(Form("sigma4cbMC%d_%s",_count, pdf.Data()),"",0.0266,0.01,0.1) ;
	RooRealVar alphaMC(Form("alphaMC%d_%s",_count,pdf.Data()),"",4.,0,40);
	RooRealVar nMC(Form("nMC_%d_%s", _count, pdf.Data()),"",100,0,500);

	RooRealVar* scale;
	scale = new RooRealVar("scale","scale",1,0,2);

	RooProduct scaled_sigma1MC("scaled_sigma1MC","scaled_sigma1MC", RooArgList(*scale,sigma1MC));
	RooProduct scaled_sigma2MC("scaled_sigma2MC","scaled_sigma2MC", RooArgList(*scale,sigma2MC));
	RooProduct scaled_sigma3MC("scaled_sigma3MC","scaled_sigma3MC", RooArgList(*scale,sigma3MC));
	RooProduct scaled_sigma4cbMC("scaled_sigma4cbMC","scaled_sigma4cbMC", RooArgList(*scale,sigma4cbMC));

	RooGaussian sig1MC(Form("sig1MC%d_%s",_count,pdf.Data()),"",*mass,meanMC,scaled_sigma1MC);  
	RooGaussian sig2MC(Form("sig2MC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma2MC);  
	RooGaussian sig3MC(Form("sig3MC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma3MC);  
	RooCBShape  CBMC(Form("CBMC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma4cbMC, alphaMC, nMC);

	RooRealVar sig1fracMC(Form("sig1fracMC%d_%s",_count, pdf.Data()),"", 0.2, 0.001, 0.99);
	RooRealVar sig2fracMC(Form("sig2fracMC%d_%s",_count, pdf.Data()),"", 0.7, 0.001, 0.99);

	RooRealVar nsigMC(Form("nsigMC%d_%s",_count, pdf.Data()),"",0.5 * dsMC->sumEntries(), 0, 1.2 * dsMC->sumEntries());

	RooAddPdf* sigMC;
	if((variation=="" && pdf=="") || variation== "background" || (variation=="signal" && pdf=="fixed" )) sigMC = new RooAddPdf(Form("sigMC%d_%s",_count,pdf.Data()),"",RooArgList(sig1MC,sig2MC),sig1fracMC);
	if(variation=="signal" && pdf=="3gauss") sigMC = new RooAddPdf(Form("sigMC%d_%s",_count, pdf.Data()), "", RooArgList(sig1MC, sig2MC, sig3MC), RooArgList(sig1fracMC, sig2fracMC), true);
	if(variation=="signal" && pdf=="gauss_cb") sigMC = new RooAddPdf(Form("sigMC%d_%s",_count, pdf.Data()), "", RooArgList(sig1MC, CBMC), sig1fracMC);

	RooAddPdf* modelMC;
	if((variation=="signal" && (pdf=="gauss_cb"|| pdf=="3gauss"|| pdf=="fixed"))||variation=="background") modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(*sigMC),RooArgList(nsigMC));
	if(variation =="" && pdf=="") modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(*sigMC),RooArgList(nsigMC));

//////////ROOFIT ROOFIT ROOFIT  MC MC MC MC 
	cout << "Starting the signal MC fit" << endl;

	mass->setRange("signal",init_mean-0.2, init_mean+0.2);    //focus the MC fit to the signal region to prevent statistical flutuations
	scale->setConstant();
	RooFitResult * fitResultMC;
	fitResultMC = modelMC->fitTo(*dsMC,Save(), Range("signal"));

///////////ROOFIT ROOFIT ROOFIT MC MC MC MC

/*	std::cout<<"mean_MC= "<<meanMC.getVal()<<std::endl;
	std::cout<<"sigma1_MC= "<<sigma1MC.getVal()<<std::endl;
	std::cout<<"sigma2_MC= "<<sigma2MC.getVal()<<std::endl;
	std::cout<<"fraction1_MC= "<<sig1fracMC.getVal()<<std::endl;
	std::cout<<"fraction2_MC= "<<sig2fracMC.getVal()<<std::endl;
*/
	cMC->cd();
	dsMC->plotOn(frameMC,Name(Form("dsMC_cut%d",_count)),Binning(nbinsmasshisto),MarkerSize(1.55),MarkerStyle(20),LineColor(1),LineWidth(4));
	modelMC->plotOn(frameMC,Name(Form("sigMC%d_%s",_count, pdf.Data())),Components(*sigMC),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("L"),LineColor(kOrange-3),LineWidth(4));
	modelMC->plotOn(frameMC,Name(Form("sigFMC%d_%s",_count, pdf.Data())),Components(*sigMC),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("F"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(4));
	modelMC->plotOn(frameMC,Name(Form("modelMC%d_%s",_count, pdf.Data())),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("L"),LineColor(2),LineWidth(4));
	modelMC->paramOn(frameMC,Layout(0.62, 0.955, 0.9), Format("NEU",AutoPrecision(1)));
	frameMC->getAttText()->SetTextSize(0.02);
	double plot_minMC = 5.25;
	if(tree=="ntKp") plot_minMC = 5.1;
	frameMC->GetXaxis()->SetRangeUser(plot_minMC, 5.5);
	frameMC->Draw();
	cMC->RedrawAxis();
	TLatex* texB_pt;
	texB_pt = new TLatex(0.2,0.80, Form("%d < p_{T} <%d GeV",ptmin, ptmax));
	texB_pt->SetNDC();
	texB_pt->SetTextFont(42);
	texB_pt->SetTextSize(0.03);
	texB_pt->SetLineWidth(1);
	texB_pt->Draw();
	//cMC->SetLogy();

// FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp
if(npfit != "1" && variation=="" && pdf==""){ 

	// DEFINE MODEL to fit the non prompt background
	// MC Part. Reconstructed Background Model
	RooRealVar* m_nonprompt_scale=0;
	RooRealVar* m_nonprompt_shift=0;
	m_nonprompt_scale = new RooRealVar(Form("m_nonprompt_scale%d",_count), "m_nonprompt_scale",0.04, 0.02, 0.065);
  	m_nonprompt_shift = new RooRealVar(Form("m_nonprompt_shift%d",_count), "m_nonprompt_shift", 5.13425, 5.1, 5.25);
	RooGenericPdf* erfc = new RooGenericPdf(Form("erfc%d",_count), "erfc", Form("TMath::Erfc((Bmass-m_nonprompt_shift%d)/m_nonprompt_scale%d)",_count,_count), RooArgList(*mass, *m_nonprompt_scale, *m_nonprompt_shift));
	// MC Part. Reconstructed Background Model

	// MC Combinatorial Background Model
	RooRealVar* alpha_np; 
	if (ptmin==50){ alpha_np = new RooRealVar(Form("alpha_np%d", _count), "alpha_np",-0.1, -2, -0.05);}
	else{alpha_np = new RooRealVar(Form("alpha_np%d", _count), "alpha_np",-0.6, -3., 1.);}
	RooExponential* COMB_jpsi = new RooExponential(Form("COMB_jpsi%d",_count), "COMB_jpsi", *mass, *alpha_np);
	// MC Combinatorial Background Model

	RooRealVar jpsinp_fraction(Form("jpsinp_fraction%d",_count), "fraction", 0.35, 0.2, 1);
	RooAddPdf* m_jpsinp_cont = new RooAddPdf(Form("m_jpsinp_cont%d",_count), "model for jpsi nonprompt bg", RooArgList(*COMB_jpsi, *erfc), RooArgList(jpsinp_fraction));
	w.import(*m_jpsinp_cont);
	// DEFINE MODEL to fit the non prompt background
	
	fit_jpsinp(w, ptmin, ptmax);}
// FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp

	//FIT THE DATA FIT THE DATA FIT THE DATA
	scale->setConstant(false);	
	RooAddPdf* model;

	c->cd();
	RooPlot* frame = mass->frame("");
	TPad *p1 = new TPad("p1","p1",0.,0.28,1.,0.99);
	p1->SetBorderMode(1); 
	p1->SetFrameBorderMode(0); 
	p1->SetBorderSize(2);
	p1->SetBottomMargin(0.10);
	p1->Draw(); 

	TPad *p2 = new TPad("p2","p2",0.,0.10,1.,0.3); 
	p2->SetTopMargin(0.);    
	p2->SetBorderMode(0);
	p2->SetBorderSize(2); 
	p2->SetFrameBorderMode(0); 
	p2->SetTicks(1,1); 
	p2->Draw();

	double n_signal_initial = ds->sumEntries(TString::Format("abs(Bmass-%g)<0.05",init_mean)) - ds->sumEntries(TString::Format("abs(Bmass-%g)<0.10&&abs(Bmass-%g)>0.05",init_mean,init_mean));
	if(n_signal_initial<0){n_signal_initial=1;}
	p1->cd();

///////////////// SIGNAL FUNCTIONS

	RooRealVar nsig(Form("nsig%d",_count),"",n_signal_initial,-1,ds->sumEntries()*3);
	RooRealVar mean(Form("mean%d",_count),"",meanMC.getVal(),m1,m2) ;
	RooRealVar sigma1(Form("sigma1%d",_count),"",sigma1MC.getVal(),0.01,0.1) ;
	RooRealVar sigma2(Form("sigma2%d",_count),"",sigma2MC.getVal(),0.01,0.1) ;
	RooRealVar sigma3(Form("sigma3%d",_count),"",sigma3MC.getVal(),0.01,0.1) ;
	RooRealVar sigma4cb(Form("sigma4cb%d",_count),"",sigma4cbMC.getVal(),0.01,0.1) ;
	RooRealVar alpha(Form("alpha%d_%s",_count,pdf.Data()),"",alphaMC.getVal(),0,50);
	RooRealVar n(Form("n_%d_%s", _count, pdf.Data()),"",nMC.getVal(),0,500);
	RooProduct scaled_sigma1(Form("scaled_sigma1%d",_count),"scaled_sigma1", RooArgList(*scale,sigma1));
	RooProduct scaled_sigma2(Form("scaled_sigma2%d",_count),"scaled_sigma2", RooArgList(*scale,sigma2));
	RooProduct scaled_sigma3(Form("scaled_sigma3%d",_count),"scaled_sigma3", RooArgList(*scale,sigma3));
	RooProduct scaled_sigmacb(Form("scaled_sigmacb%d",_count),"scaled_sigmacb", RooArgList(*scale,sigma4cb));
	RooGaussian sig1(Form("sig1%d",_count),"",*mass,mean,scaled_sigma1);  
	RooGaussian sig2(Form("sig2%d",_count),"",*mass,mean,scaled_sigma2);  
	RooGaussian sig3(Form("sig3%d",_count),"",*mass,mean,scaled_sigma3);  
	RooCBShape  CB(Form("CB%d_%s",_count, pdf.Data()),"",*mass,mean,scaled_sigmacb, alpha, n);
	RooRealVar c1(Form("c1%d",_count),"",1.,0.,5.);
	RooRealVar sig1frac(Form("sig1frac%d",_count),"",sig1fracMC.getVal(),0.,1.);
	RooRealVar sig2frac(Form("sig2frac%d",_count),"",sig2fracMC.getVal(),0.,1.);

	RooAddPdf* sig;
	if(variation=="signal" && pdf=="3gauss") sig = new RooAddPdf(Form("sig%d",_count), "", RooArgList(sig1, sig2, sig3), RooArgList(sig1frac, sig2frac), true);
	if(variation=="signal" && pdf=="gauss_cb") sig = new RooAddPdf(Form("sig%d",_count),"",RooArgList(sig1, CB), sig1frac);
	if((variation=="" && pdf=="") || variation== "background" || (variation=="signal" && pdf=="fixed")) sig = new RooAddPdf(Form("sig%d",_count),"",RooArgList(sig1,sig2),sig1frac);
///////////////// SIGNAL FUNCTIONS

/////////////////  BACKGROUND FUNCTIONS
	RooRealVar nbkg(Form("nbkg%d",_count),"",ds->sumEntries() - n_signal_initial,0.,ds->sumEntries());
	RooRealVar a0(Form("a0%d",_count),"",-0.1,-5,5);
	RooRealVar a1(Form("a1%d",_count),"",0.01,-5,5);
	RooRealVar a2(Form("a2%d",_count),"",1,-5e3,5e3);
	RooPolynomial bkg_1st(Form("bkg%d",_count),"",*mass,RooArgSet(a0));
	RooPolynomial bkg_2nd(Form("bkg%d",_count),"",*mass,RooArgSet(a0,a1));
	RooPolynomial bkg_3rd(Form("bkg%d",_count),"",*mass,RooArgSet(a0,a1,a2));
	RooRealVar lambda(Form("lambda%d", _count), "lambda",-0.5, -3., 1.);
	RooExponential bkg(Form("bkg%d",_count),"",*mass,lambda);
	RooRealVar nbkg_part_r(Form("nbkg_part_r%d",_count),"",1,0,1e5);

	// B+ PEAKING AND PART. RECONSTRUCTED BACKGROUNDS
	RooProduct *nbkg_peaking;
	RooAbsPdf* jpsipi ;   
	RooAbsPdf* erfc ;   
	if(npfit != "1"){ 
	jpsipi = w.pdf("jpsipi");   
	erfc = w.pdf(Form("erfc%d",_count));   
	RooRealVar* jpsipi_to_signal_ratio = w.var("jpsipi_to_signal_ratio"); 
	nbkg_peaking = new RooProduct(Form("nbkg_peaking%d",_count), "number of jpsi pi with fixed ratio to n_signal", RooArgList(nsig, *jpsipi_to_signal_ratio));
				}
	// B+ PEAKING AND PART. RECONSTRUCTED BACKGROUNDS

/////////////////  BACKGROUND FUNCTIONS

//////////////// MODEL MODEL MODEL MODEL
/////////////////Bs Bs Bs Bs Bs Bs Bs Bs
	if((variation=="" && pdf=="") || (pdf=="mass_range")) model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg),RooArgList(nsig,nbkg));
	if(variation=="background" && pdf=="1st") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg_1st),RooArgList(nsig,nbkg));
	if(variation=="background" && pdf=="2nd") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg_2nd),RooArgList(nsig,nbkg));
	if(variation=="background" && pdf=="3rd") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg_3rd),RooArgList(nsig,nbkg));
	if(variation=="signal" && pdf=="1gauss") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(sig1,bkg),RooArgList(nsig,nbkg));
	if(variation=="signal" && (pdf=="3gauss"|| pdf=="fixed"|| pdf=="gauss_cb" )) model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig, bkg),RooArgList(nsig, nbkg));

/////////////////BP BP BP BP BP BP BP BP
    if(npfit != "1" && variation=="" && pdf=="") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(bkg,*sig, *erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
	if(npfit != "1" && pdf=="mass_range"){ model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg,*jpsipi),RooArgList(nsig,nbkg,*nbkg_peaking));}
    if(npfit != "1" && variation=="background" && pdf=="1st") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(bkg_1st,*sig,*erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
	if(npfit != "1" && variation=="background" && pdf=="2nd") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(bkg_2nd,*sig,*erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
	if(npfit != "1" && variation=="background" && pdf=="3rd") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(bkg_3rd,*sig,*erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
	if(npfit != "1" && variation=="signal" && pdf=="1gauss") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(bkg,sig1,*erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
	if(npfit != "1" && (variation=="signal" && (pdf=="3gauss"|| pdf=="fixed"|| pdf=="gauss_cb" ))) model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig, bkg, *erfc, *jpsipi),RooArgList(nsig, nbkg, nbkg_part_r,*nbkg_peaking));
//////////////// MODEL MODEL MODEL MODEL

cout << "AQUI5" << endl;

//////////////// SET PARAMETERS FROM MC FITS
	sigma1.setConstant();
	if(pdf!="1gauss"){
		sigma2.setConstant();
		sig1frac.setConstant();
	}
	if(variation=="signal" && pdf=="3gauss"){
		sigma3.setConstant();
		sig2frac.setConstant();  
	}
	if(variation=="signal" && pdf=="gauss_cb"){
		sigma4cb.setConstant();
		n.setConstant();
		alpha.setConstant();
	}
	if(variation=="signal" && pdf=="fixed") mean.setConstant();

//////////////// SET PARAMETERS FROM MC FITS

////// ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT

	mass->setRange("m_range", 5.19 , 6.);    //set a range to be used if pdf = mass_range
	mass->setRange("all", minhisto, maxhisto);    

  	TString fitRange = (pdf == "mass_range") ? "m_range" : "all";
	RooFitResult* fitResult = model->fitTo(*ds,Save(), Minos(),Extended(kTRUE), Range(fitRange));

////// ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT

	w_val->import(*model);
	w_val->import(nsig);

  ds->plotOn(frame, Name(Form("ds_cut%d", _count)), Binning(nbinsmasshisto), MarkerSize(1), MarkerStyle(20), MarkerColor(1), LineColor(1), LineWidth(2)); 
  
  if(npfit != "1")	{
	//TString option = (pdf == "mass_range")? "L" : "LF";
    //RooCmdArg drawRange = (pdf == "mass_range")? Range(fitRange) : RooCmdArg();
    model->plotOn(frame,  Name(Form("erfc%d",_count)) , Components(Form("erfc%d",_count)),Range(fitRange), Precision(1e-6), DrawOption("L"), FillStyle(3005),FillColor(kGreen+4), LineStyle(1), LineColor(kGreen+4), LineWidth(3));
  	model->plotOn(frame, RooFit::Name("B->J/#psi #pi"),Components("jpsipi"), NormRange("bmass"),LineColor(kPink+10), LineStyle(kDashed));
					}

	model->plotOn(frame, Name(Form("bkg%d", _count)) ,  Components(bkg), Range(fitRange), Precision(1e-6), DrawOption("L"),LineStyle(7), LineColor(4), LineWidth(3));
	model->plotOn(frame, Name(Form("model%d", _count)), Range(fitRange), Precision(1e-6), DrawOption("L"), LineColor(2), LineWidth(3));

	if(pdf!="1gauss") {model->plotOn(frame, Name(Form("sig%d", _count)),  Components(*sig), Range(fitRange), Precision(1e-6), DrawOption("LF"), FillStyle(3002), FillColor(kOrange-3), LineStyle(7), LineColor(kOrange-3), LineWidth(3));} 
	else {model->plotOn(frame, Name(Form("sig%d", _count)),  Components(sig1), NormRange(fitRange), Precision(1e-6), DrawOption("LF"), FillStyle(3002), FillColor(kOrange-3), LineStyle(7), LineColor(kOrange-3), LineWidth(3));}

  if(!drawLegend){model->paramOn(frame,Layout(0.6, 0.98, 0.68), Format("NEU",AutoPrecision(3)));} 

	frame->getAttText()->SetTextSize(0.00);
	frame->getAttFill()->SetFillStyle(0);
	frame->getAttLine()->SetLineWidth(0);
	frame->SetTitle("");
	frame->SetXTitle("");
	frame->SetYTitle("Events / (20 MeV/c^{2})");
	frame->GetYaxis()->CenterTitle();
	frame->GetYaxis()->SetTitleOffset(1.4);
	frame->GetYaxis()->SetTitleSize(0.05);
	frame->GetYaxis()->SetTitleFont(42);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetYaxis()->SetLabelSize(0.05);
	frame->SetStats(0);
	double plot_min = 5.25;
	if(tree=="ntKp") plot_min = 5.05;
	frame->GetXaxis()->SetRangeUser(plot_min,5.5);
	frame->GetXaxis()->SetNdivisions(-50205);	
	frame->Draw();

	RooHist* pull_hist = frame->pullHist(Form("ds_cut%d",_count),Form("model%d",_count));
	//  RooHist* pull_hist = frame->pullHist("Data","Fit");

	RooPlot* pull_plot = mass->frame();
	(pull_plot->GetXaxis())->SetRangeUser(minhisto,maxhisto); //maxhisto
	RooRealVar x("x","",1e3,0,1e6);
	x.setVal(0.);
	RooGenericPdf* line_ref = new RooGenericPdf("ref_0", "ref_0", "x", RooArgList(x));
	line_ref->plotOn(pull_plot, LineStyle(7), LineColor(13), LineWidth(2));  
	pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"XP");
	pull_plot->SetTitle("");

	if(tree=="ntKp")pull_plot->SetXTitle("m_{J/#psiK^{#pm}} (GeV/c^{2})");
	if(tree=="ntphi")pull_plot->SetXTitle("m_{J/#psiK^{+}K^{-}} (GeV/c^{2})");

	pull_plot->SetYTitle("Pull");
	pull_plot->GetYaxis()->SetTitleFont(42);  
	pull_plot->GetYaxis()->SetTitleSize(0.17);
	pull_plot->GetYaxis()->CenterTitle();
	pull_plot->GetYaxis()->SetLabelOffset(0.01);
	pull_plot->GetYaxis()->SetLabelSize(0.15);
	pull_plot->GetYaxis()->SetNdivisions(305);
	pull_plot->GetXaxis()->SetTitleSize(0.20);
	pull_plot->GetXaxis()->SetTitleOffset(1.0);
	pull_plot->GetXaxis()->CenterTitle();
	pull_plot->GetXaxis()->SetLabelFont(42);
	pull_plot->GetXaxis()->SetLabelOffset(0.01);
	pull_plot->GetXaxis()->SetLabelSize(0.15);
	pull_plot->GetXaxis()->SetTickLength(0.15);
	pull_plot->GetXaxis()->SetNdivisions(-50205);

	/*
	   mass->setRange("signal",5.2,5.6);
	   RooAbsReal* sigIntegral = sig->createIntegral(*mass,NormSet(*mass),Range("signal"));
	   double sigIntegralErr = sigIntegral->getPropagatedError(*fitResult);
	   RooAbsReal* bkgIntegral = bkg.createIntegral(*mass,NormSet(*mass),Range("signal"));
	   double bkgIntegralErr = bkgIntegral->getPropagatedError(*fitResult);
	   cout<<"nsig: "<<nsig.getVal()<<endl;
	//cout<<"nsig: "<<model->getParameters(*mass)->getRealValue("nsig")<<endl;
	cout<<"nsig error: "<<nsig.getError()<<endl;
	cout<<"sig integral: "<<sigIntegral->getVal()<<endl;
	cout<<"sig integral error: "<<sigIntegralErr<<endl;
	cout<<"nbkg: "<<nbkg.getVal()<<endl;
	//cout<<"nbkg: "<<model->getParameters(*mass)->getRealValue("nbkg")<<endl;
	cout<<"nbkg error: "<<nbkg.getError()<<endl;
	cout<<"bkg integral: "<<bkgIntegral->getVal()<<endl;
	cout<<"bkg integral error: "<<bkgIntegralErr<<endl;*/
	//cout<<"#expected events: "<<model->expectedEvents(RooArgSet(*mass))<<endl;
	//cout<<"model integral: "<<model->getNormIntegral(RooArgSet(*mass))->getVal()<<endl;
	//fitResult->Print("v");

	//RooHist* datahist = new RooHist();
	//datahist = frame->getHist("ds");
	//TGraphAsymmErrors* datagraph = static_cast<TGraphAsymmErrors*>(datahist);
	//RooCurve* modelcurve = new RooCurve();
	//modelcurve = frame->getCurve("model");

	Double_t yield = nsig.getVal();
	Double_t yieldErr = nsig.getError();
	TH1D* fh = (TH1D*)h->Clone("fh");
	double nexpected = model->expectedEvents(RooArgSet(*mass));
	double dataArr[nbinsmasshisto]; double dataErrArr[nbinsmasshisto]; double fitArr[nbinsmasshisto];
	double val = 0;
	for(int i = 0; i < nbinsmasshisto; i++){
		dataArr[i] = h->GetBinContent(i+1);
		dataErrArr[i] = h->GetBinError(i+1);
		mass->setVal(h->GetBinCenter(i+1));
		val = model->getVal(RooArgSet(*mass))*1./nbinsmasshisto*nexpected;
		fitArr[i] = val;
		fh->SetBinContent(i+1, fitArr[i]);
		fh->SetBinError(i+1, sqrt(fitArr[i]));
	}
	cout<<"frame->chiSquare: "<<frame->chiSquare(Form("model%d",_count),Form("ds_cut%d",_count), fitResult->floatParsFinal().getSize())<<endl;
	cout<<"pars= "<<fitResult->floatParsFinal().getSize()<<endl;
	cout<<"chi2: "<<frame->chiSquare(Form("model%d",_count),Form("ds%d",_count))<<endl;
	//RooChi2Var chi2_lowstat("chi2_lowstat","chi2",*model,*dh);
	//cout<<chi2_lowstat.getVal()<<endl;

	double chiRoo = frame->chiSquare(Form("model%d",_count),Form("ds_cut%d",_count), fitResult->floatParsFinal().getSize());
	double chi2Std = 0;
	double chi2Neyman = 0;
	double chi2Peason = 0;
	double chi2BakerCousins = 0;
	chi2Cal(dataArr, dataErrArr, fitArr, nbinsmasshisto, chi2Std, chi2Neyman, chi2Peason, chi2BakerCousins);
	printf("chi2 Standard: %f\n",chi2Std);
	printf("chi2 Neyman: %f\n",chi2Neyman);
	printf("chi2 Peason: %f\n",chi2Peason);
	printf("chi2 Baker & Cousins: %f\n",chi2BakerCousins);

	TLegend *leg = new TLegend(0.62,0.55,0.89,0.75,NULL,"brNDC"); 
	leg = new TLegend(0.65/1.20+0.10,0.42,0.89+0.10,0.90,NULL,"brNDC");
	
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(h,"Data","pe");
	leg->AddEntry(frame->findObject(Form("model%d",_count)),"Fit","l");
	leg->AddEntry(frame->findObject(Form("sig%d",_count)),"Signal","f");
	leg->AddEntry(frame->findObject(Form("bkg%d",_count)),"Background","l");
	if(npfit != "1") leg -> AddEntry(frame->findObject(Form("erfc%d",_count)),"B #rightarrow J/#psi X","f");
  	if (drawLegend) {leg -> Draw();}

	float posit = 0.58;
	int nDOF = nbinsmasshisto-(fitResult->floatParsFinal().getSize());
	float nChi2 = chi2BakerCousins/(nbinsmasshisto-(fitResult->floatParsFinal().getSize()));
	int nDigit_chi2BakerCousins = 2;
	int nDigit_nChi2 = 2;
	chi2BakerCousins = roundToNdigit(chi2BakerCousins);
	nChi2 = roundToNdigit(nChi2);
	nDigit_chi2BakerCousins = sigDigitAfterDecimal(chi2BakerCousins);
	nDigit_nChi2 = sigDigitAfterDecimal(nChi2);
	TLatex* texChi = new TLatex(0.65,0.64, Form("#chi^{2}/nDOF = %2f", chiRoo));
	texChi->SetNDC();
	texChi->SetTextAlign(12);
	texChi->SetTextSize(0.03);
	texChi->SetTextFont(42);

	nsig.setVal(0.);
	nsig.setConstant();
	RooFitResult* fitResult_nosig = model->fitTo(*ds,Save());
	RooAbsReal* nll_nosig = model->createNLL(*ds);
	double log_likelihood_nosig= nll_nosig->getVal();

	RooAbsReal* nll = model->createNLL(*ds);
	double log_likelihood= nll->getVal();
    real_significance = sqrt(2*(-log_likelihood+log_likelihood_nosig));
	std::cout<<"REAL SIGNIFICANCE= "<<real_significance<<std::endl;
	std::cout<<"***********************************************************************"<<std::endl;

	std::cout<<"sigma1_MC error hi= "<<sigma1MC.getAsymErrorHi()<<std::endl;
	std::cout<<"sigma1_MC error = "<<sigma1MC.getError()<<std::endl;
	std::cout<<"sigma2_MC error hi= "<<sigma2MC.getAsymErrorHi()<<std::endl;
	std::cout<<"sigma2_MC error = "<<sigma2MC.getError()<<std::endl;
	//RooFitResult* fitResult = model->fitTo(*ds,Save(),Minos());
	//ds->plotOn(frame,Name(Form("ds%d",_count)),Binning(nbinsmasshisto),MarkerSize(1.55),MarkerStyle(20),LineColor(1),LineWidth(4));

	// double width = 0.08;
	// double BmassH = BS_MASS + width;
	// double BmassL = BS_MASS - width;
	// mass->setRange("signal",BmassL,BmassH);
	RooAbsReal *bkgIntegral = bkg.createIntegral(*mass,NormSet(*mass),Range("signal"));
	// bkgIntegralErr = bkgIntegral->getPropagatedError(*fitResult);	
	cout<<"bkg integral: "<<bkgIntegral->getVal()<<endl;
	//  cout<<"bkg integral error: "<<bkgIntegralErr<<endl;

	bkgd = nbkg.getVal()*bkgIntegral->getVal();
	Significance = yield/sqrt(bkgd+yield);

	Double_t Significance =  real_significance;
	//	Double_t Significance =  yield/TMath::Sqrt(bkgd+yield);
	int nDigit_Significance = 3;
	Significance = roundToNdigit(Significance);
	nDigit_Significance = sigDigitAfterDecimal(Significance);
	TLatex* texSig = new TLatex(0.65,0.67,Form("Significance = %.*f", nDigit_Significance, Significance));
	cout<<"Significance = "<<Significance<<endl;
	texSig->SetNDC();
	texSig->SetTextFont(42);
	texSig->SetTextSize(0.03);
	texSig->SetLineWidth(2);

	cMC->cd();
//	leg->Draw("same");

	p2->cd();
	pull_plot->GetYaxis()->SetTitleOffset(0.4);
	pull_plot->SetYTitle("Pull");
	pull_plot->Draw();

	p1->cd();
	outframe = frame;
	outputw->import(*model);

	cout << "------------------------------------------------------------------------------------------------" << endl;

	Double_t yieldPrintErr = nsig.getError();
	Double_t yieldPrintErrUp = nsig.getAsymErrorHi();
	Double_t yieldPrintErrDown = -1 * nsig.getAsymErrorLo();

	cout << "yield Error = " << yieldPrintErr << "     yield Error Up = " << yieldPrintErrUp << "    yieldPrintErrDown = " << yieldPrintErrDown << endl;

	cout << "------------------------------------------------------------------------------------------------" << endl;

	RooAbsReal* RangeBakground = bkg.createIntegral(*mass,NormSet(*mass),Range("signal")); 

	//cout << "RangeBakground = " << RangeBakground << endl;

	double Calback = RangeBakground->getVal() * nbkg.getVal();
	double StatSig = yield/sqrt(yield + Calback);
	
	/*TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.045);
	lat->DrawLatex(0.63,0.55,Form("S = %.1f",yield));	
	lat->DrawLatex(0.63,0.50,Form("B = %.1f",Calback));	
	lat->DrawLatex(0.63,0.45,Form("S/#sqrt{S+B} = %.1f",StatSig));	
	lat->Draw("SAME");*/

	return fitResult;

} // END OF MAIN FITTING FUNCTION

void fit_jpsinp(RooWorkspace& w, int pti, int ptf, bool includeSignal) {

	RooDataSet* d_s = (RooDataSet*) w.data("jpsinp");
	if(pti>= 5 && ptf<=7){	     d_s = (RooDataSet*) d_s->reduce("BDT_pt_5_7 > 0.08");}
	else if(pti>= 7 && ptf<=10){ d_s = (RooDataSet*) d_s->reduce("BDT_pt_7_10 > 0.07");}
	else if(pti>= 10 && ptf<=15){d_s = (RooDataSet*) d_s->reduce("BDT_pt_10_15 > 0.0");}
	else if(pti>= 15 && ptf<=20){d_s = (RooDataSet*) d_s->reduce("BDT_pt_15_20 > 0.02");}
	else if(pti>= 20 && ptf<=50){d_s = (RooDataSet*) d_s->reduce("BDT_pt_20_50 > 0.04");}
	else {d_s->reduce("(BDT_pt_5_7 > 0.08 && Bpt >= 5 && Bpt < 7) || (BDT_pt_7_10 > 0.07 && Bpt >= 7 && Bpt < 10) || (BDT_pt_10_15 > 0.0 && Bpt >= 10 && Bpt < 15) || (BDT_pt_15_20 > 0.02 && Bpt >= 15 && Bpt < 20) || (BDT_pt_20_50 > 0.04 && Bpt >= 20 && Bpt < 50) || (Bpt >= 20 && Bpt < 50) ");}

  	// Apply y selections and pT bins
  	d_s = (RooDataSet*) d_s->reduce(Form("(Bpt>%d && Bpt < %d)&&((Bpt < 10 &&  abs(By) > 1.5 ) || (Bpt > 10))",pti , ptf) );
	// Get rid of B+ at gen level
	RooDataSet* ds_cont = (RooDataSet*) d_s->reduce("Bgen != 23333 && Bgen != 23335");

	// Crteate the necessary folders and define paths
  	gSystem->mkdir(Form("./results/BP/%i_%i", pti, ptf),true);
	TString jpsi_fit_plot = "./results/BP/" + TString::Format("%i_%i/np_fit_pt%i-%i%s.pdf", pti, ptf, pti, ptf, "");
  	TString jpsi_plot_with_sig = "./results/BP/" + TString::Format("%i_%i/np_fit_signal_pt%i-%i%s.pdf", pti, ptf, pti, ptf, "");

	//[START] FIX SHAPE (NP background)
	RooAbsPdf* m_jpsinp_cont = w.pdf(Form("m_jpsinp_cont%d",_count));
	// FIT
	auto cont_result = m_jpsinp_cont->fitTo(*ds_cont, Save());
	plot_jpsifit(w, m_jpsinp_cont, ds_cont, jpsi_fit_plot, false);
	fix_parameters(w, Form("m_jpsinp_cont%d",_count));
	//[END] FIX SHAPE (NP background) 

	// TEST THE FIT TEST THE FIT
	// Import from WorkSpace and define variables
	RooAbsPdf* erfc = w.pdf(Form("erfc%d",_count));
	RooAbsPdf* COMB_jpsi = w.pdf(Form("COMB_jpsi%d",_count));
	RooAbsPdf* signalnp = w.pdf("signalnp");    
	RooAbsPdf* jpsipi = w.pdf("jpsipi");    
  	RooRealVar* sigma1_np = w.var("sigma1np");          
	RooRealVar* jpsipi_to_signal_ratio = w.var("jpsipi_to_signal_ratio");  
	RooRealVar n_cont(Form("n_cont_np%d",_count), "n_cont_np", 1000, 0., (d_s->sumEntries())*2);
	RooRealVar n_erfc(Form("n_nonprompt%d",_count), "n_nonprompt", 1000, 0., (d_s->sumEntries())*2);
	RooRealVar nbin_signal(Form("n_signal_np%d",_count), "n_signal_np", 1000, 0., 150000);
	RooProduct* n_jpsipi = new RooProduct(Form("n_jpsipi_by_signal%d",_count), "number of jpsis pi with fixed ratio to n_signal", RooArgList(nbin_signal, *jpsipi_to_signal_ratio));
	
	// Unfix the signal with to better describe each pT bin peak
  	sigma1_np->setConstant(false);
	// BUILD TOTAL PDF
	RooAddPdf* model_inclusive = new RooAddPdf(Form("model_inclusive%d",_count), "NP with B+", RooArgList(*signalnp, *jpsipi, *erfc, *COMB_jpsi), RooArgList(nbin_signal, *n_jpsipi, n_erfc, n_cont));
    // FITFITFIT JUST TO CHECK 
    model_inclusive->fitTo(*d_s, Save(), Extended(), NumCPU(4));
	// Plot
    plot_jpsifit(w, model_inclusive, d_s, jpsi_plot_with_sig, true);
	// TEST THE FIT TEST THE FIT
}

void plot_jpsifit(RooWorkspace& w, RooAbsPdf* model, RooDataSet* ds, TString plotName, bool with_sig) {
  RooRealVar Bmass = *(w.var("Bmass"));
  Bmass.setRange("bmass", 5.0, 6.0);
  RooPlot* massframe = Bmass.frame(Title(" "));
  ds->plotOn(massframe, RooFit::Name("NP"), MarkerSize(0.9));
  model->plotOn(massframe, RooFit::Name("NP Fit"), NormRange("bmass"),LineColor(kRed), LineStyle(1), LineWidth(2));
  model->plotOn(massframe, RooFit::Name("peaking"),Components(Form("erfc%d",_count)), NormRange("bmass"), LineColor(kGreen+3), LineStyle(1), LineWidth(3), DrawOption("L"));
  model->plotOn(massframe, RooFit::Name("COMB_jpsi"),Components(Form("COMB_jpsi%d",_count)), NormRange("bmass"),LineColor(kBlue), LineStyle(kDashed));
  model->paramOn(massframe,  Layout(0.1, 0.65, 0.4),Format("NEU", AutoPrecision(1)));
  if (with_sig) {
    model->plotOn(massframe, RooFit::Name("signal"),Components("signalnp"), NormRange("bmass"),
                  LineColor(kOrange-3), LineStyle(1), LineWidth(3), FillStyle(3002),FillColor(kOrange-3), VLines(), DrawOption("LF"));
    model->plotOn(massframe, RooFit::Name("B->J/#psi #pi"),Components("jpsipi"), NormRange("bmass"),LineColor(kPink+10), LineStyle(kDashed));
  }

  TCanvas can_np;
  can_np.SetTitle("");
  massframe->GetXaxis()->SetRangeUser(5.00,5.5);
  massframe->Draw();

  TLatex txt;
  TLegend *leg = new TLegend (0.75, 0.55, 0.85, 0.9);
  leg->SetTextSize(0.04);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(massframe->findObject("NP"), "NP MC", "p)");
  leg->AddEntry(massframe->findObject("peaking"), "erfc", "l");
  leg->AddEntry(massframe->findObject("COMB_jpsi"), "combinatorial", "l");
  leg->AddEntry(massframe->findObject("NP Fit"),"Fit","l");
  // compare yields with gen particles
  if (with_sig) {
    leg->AddEntry(massframe->findObject("signal"), "signal", "f");
    leg->AddEntry(massframe->findObject("B->J/#psi #pi"), "B->J/#psi #pi", "l");
  				}
  leg->Draw();
  can_np.SaveAs(plotName);
				}

/* Fix or release the parameters for a given PDF */
void fix_parameters(RooWorkspace& w, TString pdfName, bool release) {
  RooAbsPdf* pdf = w.pdf(pdfName);
  RooAbsData* ds = w.data("data");
  RooArgSet* par_set = pdf->getParameters(*ds);
  auto itr = par_set->createIterator();
  bool toFix = ! release;
  std::string fix_or_float = (toFix)? "fix " : "float ";
  std::cout << fix_or_float << "parameters:";
  for (auto i = 0; i < par_set->getSize(); ++i) {
    RooRealVar* var = (RooRealVar*) itr->Next();
    var->setConstant(true);
    TString name = var->GetName();
    std::cout << name << ", ";
    var->setConstant(toFix);
  }
  std::cout << "\n";
}

void latex_table(std::string filename, int n_col, int n_lin, std::vector<std::string> col_name, std::vector<std::string> labels, 
		std::vector<std::vector<double> > numbers, std::string caption){
	
	std::ofstream file_check;
	std::ofstream file;

	//Begin Document
                                                                                    
	file.open(filename + ".tex");
	file_check.open(filename + "_check.tex");

	file_check << "\\documentclass{article}" << std::endl;
	//file << "\\usepackage[utf8]{inputenc}" << std::endl;     
	file_check << "\\usepackage{rotating}" << std::endl;                                                                                   
	// file_check << "\\usepackage{cancel}" << std::endl;
	file_check << "\\usepackage{geometry}" << std::endl;
	file_check << "\\usepackage{booktabs}" << std::endl;
	file_check << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm,}" << std::endl;

	file_check << "\\title{Bs/B+}" << std::endl;
	file_check << "\\author{Henrique Legoinha}" << std::endl;
	file_check << "\\date{June 2022}" << std::endl;
	file_check << "\\begin{document}" << std::endl;
	file_check << "\\maketitle" << std::endl;

	// Create table                                                                                                                                
	//file_check << "\\begin{table}[!h]" << std::endl;
	//file_check << "\\centering" << std::endl;                                                                                                         
	//setup table size                                                                                                                             
	std::string col="c";
	//std::setprecision(3);

	for(int i=1; i<n_col; i++)
		col+="|c";

	file_check << "\\begin{sidewaystable}"<< std::endl;
	file_check << "\\begin{tabular}{"+col+"}" << std::endl;
	file_check << "\\toprule" << std::endl;
	file << "\\begin{tabular}{"+col+"}" << std::endl;
	file << "\\toprule" << std::endl;

	for(int c=0; c<n_col-1; c++){
		file << col_name[c] << " & ";
		file_check << col_name[c] << " & ";
	}

	file << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;
	file_check << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;


	for(int i=1; i<n_lin; i++)
	{
		file << labels[i-1] << " & ";
		file_check << labels[i-1] << " & ";

		for(int c=1; c<n_col-1; c++){
			file << std::setprecision(3)<< numbers[c-1][i-1]<< " \\% & ";
			file_check << std::setprecision(3) << numbers[c-1][i-1]<< " \\% & ";
									}

		file << std::setprecision(3)<< numbers[n_col-2][i-1]<< " \\% \\\\" << std::endl;
		file_check << std::setprecision(3)<< numbers[n_col-2][i-1]<< " \\% \\\\" << std::endl; 
	}

	file << "\\bottomrule" << std::endl;
	file_check << "\\bottomrule" << std::endl;

	//End Table                                                                                                                                    
	file << "\\end{tabular}" << std::endl;
	file_check << "\\end{tabular}" << std::endl;
	file_check << "\\caption{"+caption+"}" << std::endl;
	file_check << "\\end{sidewaystable}"<< std::endl;

	//file_check << "\\end{table}" << std::endl;
	//End document                                                                                                                                 

	file_check << "\\end{document}" << std::endl;

	file.close();
	file_check.close();

	system(("pdflatex " + filename+ "_check.tex").c_str());
	system(("open " + filename + "_check.pdf").c_str());
}

void validate_fit(RooWorkspace* w, TString tree, TString variable, int full, int CentMin, int CentMax, int ptMin, int ptMax)
{
	std::cout << "Now Perform Check on Fit" << std::endl;
	RooRealVar Bmass = *(w->var("Bmass"));
	RooAbsPdf* model  = w->pdf(Form("model%d",_count));
	//RooDataSet* data = (RooDataSet*) w->data("data");

	//model->fitTo(*data);
	vector<RooRealVar> params;
	params.push_back(*(w->var(Form("nsig%d",_count))));
	params.push_back(*(w->var(Form("mean%d",_count))));

	/*
	//Retrieving Parameters//

	RooRealVar* lambda =  (w->var(Form("lambda%d", _count)));
	RooRealVar* nbkg   =  (w->var(Form("nbkg%d", _count)));
	RooRealVar * mean  = (w->var(Form("mean%d", _count)));
	RooRealVar * sigma1  = (w->var(Form("sigma1%d", _count)));
	RooRealVar * sigma2  = (w->var(Form("sigma2%d", _count)));
	//	RooRealVar * c1  = (w->var(Form("c1%d", _count)));
	RooRealVar * sig1frac = (w->var(Form("sig1frac%d", _count)));

	lambda->setVal(lambda->getVal());
	nbkg->setVal(nbkg->getVal());
	mean->setVal(mean->getVal());
	sigma1->setVal(sigma1->getVal());
	sigma2->setVal(sigma2->getVal());
	//	c1->setVal(c1->getVal());
	sig1frac->setVal(sig1frac->getVal());
	*/

	/*
	// test removing bakcground
	RooRealVar* lambda =  (w->var(Form("lambda%d", _count)));
	RooRealVar* nbkg   =  (w->var(Form("nbkg%d", _count)));

	//	double step = 0.002;

	//	double nkgdValue = 51.15;
	//	double lambdaValue = -2.0236;

	cout << "BRO lambda = " << lambda->getVal() << endl;
	cout << "BRIS back = " << nbkg->getVal() << endl;

	double nkgdValue = nbkg->getVal();
	double lambdaValue = lambda->getVal();

	nbkg->setVal(50);
	lambda->setVal(lambdaValue);
	*/
	//	RooRealVar* nbkg   =  (w->var(Form("nbkg%d", _count)));
	//	nbkg->setVal(50);

	double n_signal_init = params[0].getVal();
	double n_signal_error_init = params[0].getError();

	cout << "N_signal initial value: " << n_signal_init << endl;
	cout << "N_signal initial Error value: " << n_signal_error_init << endl;

	int params_size = params.size();

	cout << "params_size " << params_size << endl;

	RooMCStudy* mcstudy = new RooMCStudy(*model, Bmass,  Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));

	mcstudy->generateAndFit(5000);

	cout << "DONE Generate and Fit " << endl;


	TString XName[2] = {"nsig","mean"};

	vector<RooPlot*> framesPull, framesParam, framesError;

	for(int i = 0; i < params_size; ++i)
	{
		framesPull.push_back(mcstudy->plotPull(params.at(i),FrameBins(50),FrameRange(-5,5)));
		framesPull[i]->SetTitle("");
		framesParam.push_back(mcstudy->plotParam(params.at(i),FrameBins(50)));
		framesParam[i]->SetTitle("");
		framesError.push_back(mcstudy->plotError(params.at(i),FrameBins(50)));
		framesError[i]->SetTitle("");
	}

	vector<TGraph*> h1;
	vector<TGraph*> h2;
	vector<TGraph*> h3;

	for(int i = 0; i < params_size; ++i){
		h1.push_back(static_cast<TGraph*>(framesPull.at(i)->getObject(0)));
		h2.push_back(static_cast<TGraph*>(framesParam.at(i)->getObject(0)));
		h3.push_back(static_cast<TGraph*>(framesError.at(i)->getObject(0)));
	}

	gStyle->SetOptFit(0111);

	TCanvas* c_pull = new TCanvas("pulls", "pulls", 900, 800);

	gPad->SetLeftMargin(0.15);
	//gPad->SetBottomMargin(0.0);

	for(int i = 0; i < params_size; ++i){
		c_pull->cd();
		h1[i]->SetTitle("");
		h1[i]->Draw();
		c_pull->Update();
		h1[i]->Fit("gaus","","",-5,5);
		h1[i]->GetFunction("gaus")->SetLineColor(4);
		h1[i]->GetFunction("gaus")->SetLineWidth(5);
		h1[i]->GetXaxis()->SetTitle("Pull");
		h1[i]->GetYaxis()->SetTitle("");
		h1[i]->Draw("APsame");
	}


	TCanvas* c_params = new TCanvas("params", "params", 900, 800);
	for(int i = 0; i < params_size; ++i){
		c_params->cd();
		framesParam.at(i)->GetYaxis()->SetTitle("Toy MCs");
		framesParam.at(i)->GetYaxis()->SetTitleOffset(1.4);
		framesParam.at(i)->Draw();
	}


	/*TCanvas* c_params = new TCanvas("params", "params", 900, 800);

	  for(int i = 0; i < params_size; ++i){
	  c_params->cd();
	  h2[i]->SetTitle("");
	  h2[i]->Draw();
	  c_params->Update();
	  if(particle == 0){h2[i]->Fit("gaus","","",900,1200);}
	  else if(particle == 1){h2[i]->Fit("gaus","","",60, 120);}
	  h2[i]->GetFunction("gaus")->SetLineColor(4);
	  h2[i]->GetFunction("gaus")->SetLineWidth(5);
	//h2[i]->GetXaxis()->SetTitle();
	h2[i]->GetYaxis()->SetTitle("Toy MCs");
	h2[i]->Draw("same");
	}*/



	TCanvas* c_errors = new TCanvas("errors", "errors", 900, 800);
	gPad->SetLeftMargin(0.15);
	for(int i = 0; i < params_size; ++i){
		c_errors->cd();
		framesError.at(i)->GetYaxis()->SetTitle("Toy MCs");
		framesError.at(i)->GetYaxis()->SetTitleOffset(1.4);
		framesError.at(i)->Draw();
	}

	for(int i = 0; i < params_size; ++i){

		n_signal_init = params[i].getVal();
		n_signal_error_init = params[i].getError();

		c_pull->cd();
		h1[i]->SetTitle("");
		h1[i]->Draw();
		c_pull->Update();
		h1[i]->Fit("gaus","","",-5,5);
		h1[i]->GetFunction("gaus")->SetLineColor(4);
		h1[i]->GetFunction("gaus")->SetLineWidth(5);
		h1[i]->GetXaxis()->SetTitle(Form("%s Pull",XName[i].Data()));
		h1[i]->GetYaxis()->SetTitle("");
		h1[i]->Draw("APsame");
		c_errors->cd();

		h3[i]->SetTitle("");
		h3[i]->Draw();
		c_errors->Update();
		h3[i]->Fit("gaus","","",n_signal_error_init*0.7,n_signal_error_init*1.3);
		h3[i]->GetFunction("gaus")->SetLineColor(4);
		h3[i]->GetFunction("gaus")->SetLineWidth(5);
		h3[i]->GetXaxis()->SetTitle(Form("%s Error",XName[i].Data()));
		h3[i]->GetYaxis()->SetTitle("");
		h3[i]->Draw("APsame");

		c_params->cd();

		h2[i]->SetTitle("");
		h2[i]->Draw();
		c_params->Update();
		h2[i]->Fit("gaus","","",n_signal_init * 0.7,n_signal_init * 1.3);
		h2[i]->GetFunction("gaus")->SetLineColor(4);
		h2[i]->GetFunction("gaus")->SetLineWidth(5);
		h2[i]->GetXaxis()->SetTitle(Form("%s Mean",XName[i].Data()));
		h2[i]->GetYaxis()->SetTitle("");
		h2[i]->Draw("APsame");

		if(tree=="ntKp"){
			//c_pull->SaveAs("./mcstudy/pulls_poisson_Bu.pdf");
			c_pull->SaveAs(Form("new/pull_signal_full%d_%s_%d_%s.png",full,variable.Data(),_count,tree.Data()));
			c_params->SaveAs(Form("new/params_signal_full%d_%s_%d_%s.png",full,variable.Data(),_count,tree.Data()));
			c_errors->SaveAs(Form("new/error_signal_full%d_%s_%d_%s.png",full,variable.Data(),_count,tree.Data()));
			//c_errors->SaveAs("./results/BP/pulls/pulls_error.gif");
		}
		else if(tree=="ntphi"){
			c_pull->SaveAs(Form("newFIT/pull_signal_%s_%d_%d_%d_%d_%d_Bs.png",variable.Data(),CentMin,CentMax,ptMin,ptMax,i));
			//c_pull->SaveAs("./mcstudy/pulls_poisson_Bs.pdf");
			c_params->SaveAs(Form("newFIT/param_signal_%s_%d_%d_%d_%d_%d_Bs.png",variable.Data(),CentMin,CentMax,ptMin,ptMax,i));
			c_errors->SaveAs(Form("newFIT/error_signal_%s_%d_%d_%d_%d_%d_Bs.png",variable.Data(),CentMin,CentMax,ptMin,ptMax,i));
			//c_errors->SaveAs("./results/Bs/pulls/pulls_error.gif");
		}
	}
	/*
	   else if(tree=="ntphi"){
	   c_pull->SaveAs(Form("new/StepScan/pull/pull_signal_full%d_%s_%d_%d_Bs.png",full,variable.Data(),_count,NTrial));
//c_pull->SaveAs("./mcstudy/pulls_poisson_Bs.pdf");
c_params->SaveAs(Form("new/StepScan/param/param_signal_full%d_%s_%d_%d_Bs.png",full,variable.Data(),_count,NTrial));
c_errors->SaveAs(Form("new/StepScan/error/error_signal_full%d_%s_%d_%d_Bs.png",full,variable.Data(),_count,NTrial));
//c_errors->SaveAs("./results/Bs/pulls/pulls_error.gif");
}
*/
}


