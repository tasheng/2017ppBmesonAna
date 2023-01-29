#include "roofitB.h"
#include "CMS_lumi.C"
#include "parameter.h"                             
#include "parametersNew.h"
#include "TSystem.h"
#include <string>
#include <sstream>
#include <TGraph.h>
#include "TMultiGraph.h"
#include "TGraphErrors.h"

template<typename... Targs>
void plot_mcfit(RooWorkspace& w, RooAbsPdf* model, RooDataSet* ds, TString plotName,  Targs... options);
void read_samples(RooWorkspace& w, std::vector<TString>, TString fName, TString treeName, TString sample);

// PDF VARIATION FOR SYST STUDIES
int syst_study=1;

void roofitB(int doubly = 0, TString tree = "ntphi", int full = 0, int usePbPb = 0, int fitOnSaved = 0, TString inputdata = "", TString inputmc = "", TString varExp = "", TString trgselection = "",  TString cut = "", TString cutmcgen = "", int isMC = 0, Double_t luminosity = 1., int doweight = 1, TString outputfile = "", TString outplotf = "", TString npfit = "", int doDataCor = 0){ 

	//Create the Folders
	gSystem->mkdir("filesbp",true); 
	gSystem->mkdir("filesbs",true); 
	gSystem->mkdir(Form("%s", outplotf.Data()),true); 

	double MyBackground;
	double yield;
	int _nBins;
	TString seldata;
	TString selmc;
	TString selmcgen;

// Bpt Bpt Bpt Bpt	
	if(varExp == "Bpt"){ 
	if(full == 1){ _nBins = 1 ;}
	else if(full == 0) {
		if(tree=="ntphi"){ _nBins = nptBins;}
		//if(tree=="ntphi"){ _nBins = nptBins_test;}
		else if(tree=="ntKp"){ _nBins = nptBinsBP;}
		//else if(tree=="ntKp"){ _nBins = nptBinsBP_test;}
		}
	} else if(varExp == "By"){
		if(full == 1){_nBins = 1;}
		else if(full == 0){ cout << "bins?" << nyBins_both << endl;
			_nBins = nyBins_both;}
	}
		else if(varExp == "nMult"){
		if(full == 1){_nBins = 1;}
		else if(full == 0){_nBins = nmBins_both;}
	}

	cout << "number of bins: " << _nBins << endl;	
	double _ptBins[_nBins+1];

	if(varExp == "Bpt"){ 
	if(full == 1){
            if(tree=="ntphi"){for( int c=0; c<_nBins+1; c++){_ptBins[c]=ptBins_full[c];}}
			else if(tree=="ntKp"){for( int c=0; c<_nBins+1; c++){_ptBins[c]=ptBins_fullBP[c];}}
	}else if(full == 0) {
		if(tree=="ntphi"){for( int c=0; c<_nBins+1; c++){_ptBins[c]=ptbinsvec[c];}}
		//if(tree=="ntphi"){for( int c=0; c<_nBins+1; c++){_ptBins[c]=ptbinsvec_test[c];}}
		else if(tree=="ntKp"){for( int c=0; c<_nBins+1; c++){_ptBins[c]=ptbinsvecBP[c];}}
		//else if(tree=="ntKp"){for( int c=0; c<_nBins+1; c++){_ptBins[c]=ptbinsvecBP_test[c];}}
	                    }
	} else if(varExp == "By"){
		if(full == 0){
			for(int c=0; c<_nBins+1; c++){_ptBins[c]=ybinsvec[c];}
		}
		else if(full == 1){
			for(int c=0; c<_nBins+1; c++){_ptBins[c]=ybinsvec_full[c];}
		}
	}
		else if(varExp == "nMult"){
		if(full == 0){
			for(int c=0; c<_nBins+1; c++){_ptBins[c]=nmbinsvec[c];}
		}
		else if(full == 1){
			for(int c=0; c<_nBins+1; c++){_ptBins[c]=nmbinsvec_full[c];}
		}
	}

std::cout<<"Variable "<<varExp<<std::endl;
cout << tree << " BINS: ";
for(int t; t< sizeof(_ptBins)/sizeof(_ptBins[0]);t++){cout <<"__"<<  _ptBins[t]<<"__";}
cout << endl << endl;

/*	if(varExp=="nMult"){
		if(full==0) {
			nBins_check = nBins_Mult;
			_ptBins = MultBins;
			_nBins = nBins_Mult;
			_ptBins = MultBins;
		}
		else if(full==1){
			nBins_check = nBins_full;
			_ptBins = nMults_full;
		}
	}
	else if(varExp=="By"){
		if(full==0) {
			nBins_check = nBinsY;
			_ptBins = ptBinsY;
		}
		else if(full==1){
			nBins_check = nBins_full;
			_ptBins = hiBins_full;
		}
	}*/

	/*	if(varExp == "By"){
		_nBins = nBinsY;
		_ptBins = ptBinsY;
		}*/
	
	if (!(usePbPb==1||usePbPb==0)) std::cout<<"ERROR!!, you are using a non valid isPbPb option"<<std::endl;
	bool isPbPb=(bool)(usePbPb);

	if(!isPbPb)
	{
		std::cout<<"DEBUG 4 !PbPb"<<std::endl;
		seldata = Form("%s",trgselection.Data());
		std::cout<<"DEBUG 5"<<std::endl;
		selmc = Form("%s",cut.Data());
		std::cout<<"DEBUG 6"<<std::endl;
		selmcgen = Form("%s",cutmcgen.Data());
		std::cout<<"DEBUG 7"<<std::endl;
	}
	else
	{
		std::cout<<"DEBUG 4 else"<<std::endl;
		seldata = Form("%s&&%s",trgselection.Data(),cut.Data());
		std::cout<<"DEBUG 5"<<std::endl;
		selmc = Form("%s",cut.Data());
		std::cout<<"DEBUG 6"<<std::endl;
		selmcgen = Form("%s",cutmcgen.Data());
		std::cout<<"DEBUG 7"<<std::endl;
		std::cout<<"cut = "<<cut<<std::endl;
		std::cout<<"seldata= "<<seldata<<std::endl;
		std::cout<<"selmc= "<<selmc<<std::endl;
	}

	gStyle->SetTextSize(0.05);
	gStyle->SetTextFont(42);
	gStyle->SetPadRightMargin(cRightMargin);
	gStyle->SetPadLeftMargin(cLeftMargin);
	gStyle->SetPadTopMargin(cTopMargin);
	gStyle->SetPadBottomMargin(cBottomMargin);
	cout << "cBottomMargin = " << cBottomMargin << endl;
	gStyle->SetPadBottomMargin(0.45);
	gStyle->SetTitleX(.0f);

	TFile* inf = new TFile(inputdata.Data());
	TTree* skimtree_new = (TTree*)inf->Get(tree);
	TFile* infMC = new TFile(inputmc.Data());
	TTree* skimtreeMC_new = (TTree*)infMC->Get(tree);

	TH1D* h;
	TH1D* hMC;
	TH1D* hpull;

	TH1D* hPt = new TH1D("hPt","",_nBins,_ptBins);
	RooWorkspace* ws = new RooWorkspace("ws");
	TH1D* hPtMC = new TH1D("hPtMC","",_nBins,_ptBins);
	//TH1D* hPtGen = new TH1D("hPtGen","",_nBins,_ptBins);
	TH1D* hMean = new TH1D("hMean","",_nBins,_ptBins);                       

	std::cout<<"Histograms"<<std::endl;
	RooRealVar* Bgen = new RooRealVar("Bgen", "Bgen", 0, 30000);
	RooRealVar* mass = new RooRealVar("Bmass","Bmass",minhisto,maxhisto);
	RooRealVar* pt = new RooRealVar("Bpt","Bpt",0,100);
	RooRealVar* y = new RooRealVar("By","By",-2.4, 2.4);
	RooRealVar* nMult = new RooRealVar("nMult","nMult",0,100);
	RooRealVar BDT_pt_5_7("BDT_pt_5_7", "BDT_pt_5_7", -1, 1);
    RooRealVar BDT_pt_7_10("BDT_pt_7_10", "BDT_pt_7_10", -1, 1);
    RooRealVar BDT_pt_10_15("BDT_pt_10_15", "BDT_pt_10_15", -1, 1);
    RooRealVar BDT_pt_15_20("BDT_pt_15_20", "BDT_pt_15_20", -1, 1);
    RooRealVar BDT_pt_20_50("BDT_pt_20_50", "BDT_pt_20_50", -1, 1);
	RooRealVar BDT_pt_50_60("BDT_pt_50_60", "BDT_pt_50_60", -1, 1);

		ws->import(*mass);
		ws->import(*y);
		ws->import(*pt);
		ws->import(*Bgen);
		ws->import(BDT_pt_5_7);
    	ws->import(BDT_pt_7_10);
    	ws->import(BDT_pt_10_15);
    	ws->import(BDT_pt_15_20);
    	ws->import(BDT_pt_20_50);
		ws->import(BDT_pt_50_60);

	RooRealVar* trackSelection = new RooRealVar("track", "track", 0, 5);
	RooDataSet* ds = new RooDataSet();
	RooDataSet* dsMC = new RooDataSet();   
	std::cout<<"Created dataset"<<std::endl;
	RooDataHist* dh = new RooDataHist();   
	RooDataHist* dhMC = new RooDataHist();   
	std::cout<<"Created roodatahists"<<std::endl;
	RooPlot* frame = new RooPlot();
	RooHist* datahist = new RooHist();
	RooCurve* modelcurve = new RooCurve();

	//weightgen = weightgen_pp;
	weightmc  = weightmc_pp;
	if(usePbPb){
		weightgen = weightgen_PbPb;
		//weightmc = weightmc_PbPb;
		if(doweight == 0 ) weightmc = "1"; 
		if(doweight == 1) weightmc = "pthatweightNew";
	    weightmc = "1"; 
		
		std::cout<<"weights defined"<<std::endl;	
				}

	TString _prefix = "";
	TString _isMC = "data";
	if(isMC) _isMC = "mcAsData";
	TString _isPbPb = "pp";

	TString outputf;
	outputf = Form("%s",outputfile.Data());
	TFile* outf = new TFile(outputf.Data(),"recreate");
	outf->cd();
	dsMC = new RooDataSet(Form("dsMC%d",_count),"",skimtreeMC_new,RooArgSet(*mass, *pt, *y, *nMult, *trackSelection));
	ds = new RooDataSet(Form("ds%d",_count),"",skimtree_new,RooArgSet(*mass, *pt, *y, *nMult, *trackSelection));

	//MODELS for syst studies
	std::vector<std::string> background = {"1st", "2nd","mass_range"};
	std::vector<std::string> signal = {"3gauss", "fixed", "gauss_cb"};
	//MODELS for syst studies
	
	std::vector<std::vector<double>> background_syst;
	std::vector<std::vector<double>> signal_syst;
	std::vector<std::vector<double>> general_syst;
	std::vector<double> nominal_yields;
	std::vector<std::vector<double>> back_syst_rel_values;
	std::vector<std::vector<double>> sig_syst_rel_values;
	std::vector<std::vector<double>> stat_error;
	double yield_vec[_nBins];
	double yield_vec_err_low[_nBins];
	double yield_vec_err_high[_nBins];
	double scale_vec[_nBins];
	double scale_vec_err_low[_nBins];
	double scale_vec_err_high[_nBins];
	double resol_vec[_nBins];
	double resol_vec_err_low[_nBins];
	double resol_vec_err_high[_nBins];
	double yield_vec_systerr_low[_nBins];
	double yield_vec_systerr_high[_nBins];
	double var_mean_av[_nBins];
	double hori_av_low[_nBins];
	double hori_av_high[_nBins];

	// FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp
	//Fit the J/psi pi MC sample
   	//The shapes of J/psi pi peak is determined
	//The relative yield between signal and J/psi pi is fixed
	if(npfit != "1"){

		//PDF MODELS PDF MODELS PDF MODELS
		//inclusive MC signal Model
		RooRealVar* meannp = 0;
		RooRealVar* sigma1np = 0;
		RooProduct* sigma2np;
		RooRealVar* ratio_sigma12np = 0;
		RooRealVar* cofs_b_np = 0;
		meannp = new RooRealVar("meannp","meannp",5.281,5.24,5.32);
		sigma1np = new RooRealVar("sigma1np","sigma1np",0.02,0.01,0.03);
		ratio_sigma12np = new RooRealVar("ratio_sigma12np","ratio_sigma12np", 2.4, 0.1, 5);
		sigma2np = new RooProduct("sigma2np", "sigma2np", RooArgList(*sigma1np, *ratio_sigma12np));
		cofs_b_np = new RooRealVar("cofs_b_np", "cofs_b_np", 0.5, 0., 1.);
		RooGaussian* signal1_b_np = new RooGaussian("signal1_b_np","signal_gauss1_b_np",*mass,*meannp,*sigma1np);
		RooGaussian* signal2_b_np = new RooGaussian("signal2_b_np","signal_gauss2_b_np",*mass,*meannp,*sigma2np); 
		RooAddPdf* signalnp = new RooAddPdf("signalnp", "signalnp", RooArgList(*signal1_b_np,*signal2_b_np),*cofs_b_np);
		//inclusive MC signal Model

		//inclusive MC jpsipi Model
		RooRealVar* m_jpsipi_fraction2 = 0;
		RooRealVar* m_jpsipi_mean1 = 0;
		RooRealVar* m_jpsipi_sigma1l = 0;
		RooRealVar* m_jpsipi_sigma1r = 0;
		m_jpsipi_fraction2 = new RooRealVar("m_jpsipi_fraction2","m_jpsipi_fraction2",0.234646,0.0,0.8);
		m_jpsipi_mean1 = new RooRealVar("m_jpsipi_mean1","m_jpsipi_mean1",5.35, 5.3, 5.5);
		RooRealVar m_jpsipi_sigma2l("m_jpsipi_sigma2l","m_jpsipi_sigma2l",0.0994712,0.020,0.500);
		RooRealVar m_jpsipi_sigma2r("m_jpsipi_sigma2r","m_jpsipi_sigma2r",0.0994712,0.020,0.500);
		m_jpsipi_sigma1l = new RooRealVar("m_jpsipi_sigma1l","m_jpsipi_sigma1l",0.0290762,0.010,0.150);
		m_jpsipi_sigma1r = new RooRealVar("m_jpsipi_sigma1r","m_jpsipi_sigma1r",0.0652519,0.010,0.350);
		RooBifurGauss m_jpsipi_gaussian2("m_jpsipi_gaussian2", "m_jpsipi_gaussian2", *mass, *m_jpsipi_mean1, m_jpsipi_sigma2l, m_jpsipi_sigma2r);
		RooBifurGauss m_jpsipi_gaussian1("m_jpsipi_gaussian1", "m_jpsipi_gaussian1", *mass, *m_jpsipi_mean1, *m_jpsipi_sigma1l, *m_jpsipi_sigma1r);
		RooAddPdf* jpsipi = new RooAddPdf("jpsipi", "jpsipi", RooArgList(m_jpsipi_gaussian2, m_jpsipi_gaussian1), RooArgList(*m_jpsipi_fraction2));
		//inclusive MC jpsipi Model
		//PDF MODELS PDF MODELS PDF MODELS

		// PREPARE DATA SETS
		std::vector<TString> jpsi_vars = {"By", "Bpt", "Bgen","BDT_pt_5_7", "BDT_pt_7_10", "BDT_pt_10_15","BDT_pt_15_20", "BDT_pt_20_50"};
		//read_samples(*ws, jpsi_vars, "/data3/hlegoinha/data/jpsinp_inclusive.root", "ntnp", "jpsinp");
		read_samples(*ws, jpsi_vars, "/afs/cern.ch/user/t/tsheng/public/forHenrique/trk5/jpsinp_inclusive.root", "ntnp", "jpsinp");
		RooDataSet* full_data_MC = (RooDataSet*) ws->data("jpsinp");
		full_data_MC = (RooDataSet*)full_data_MC->reduce("(BDT_pt_5_7 > 0.08 && Bpt >= 5 && Bpt < 7) || (BDT_pt_7_10 > 0.07 && Bpt >= 7 && Bpt < 10) || (BDT_pt_10_15 > 0.0 && Bpt >= 10 && Bpt < 15) || (BDT_pt_15_20 > 0.02 && Bpt >= 15 && Bpt < 20) || (BDT_pt_20_50 > 0.04 && Bpt >= 20 && Bpt < 50) || (Bpt >= 20 && Bpt < 50) ");
		full_data_MC = (RooDataSet*)full_data_MC->reduce("(Bpt < 10 &&  abs(By) > 1.5 ) || (Bpt > 10)");
		
		// FORM INCLUSIVE SIGNAL AND PEAKING Background BINS
		RooDataSet* ds_sig = (RooDataSet*) full_data_MC->reduce("Bgen == 23333");
		RooDataSet* fullds_JPSI_shape_fix = (RooDataSet*)full_data_MC->reduce("Bgen == 23335");
		// FORM INCLUSIVE SIGNAL AND PEAKING Background BINS
		// PREPARE DATA SETS

		//[START] FIX SHAPE (J/Psi pi) 
		RooRealVar n_jpsipi_ext("n_jpsipi_ext", "n_jpsipi_ext", 1000 , 0., (fullds_JPSI_shape_fix->sumEntries())*2);
		RooExtendPdf jpsipi_ext("jpsipi_ext", "extended jpsipi", *jpsipi, n_jpsipi_ext);
		// FIT
		mass->setRange("bjpsipi", 5.22, 5.8);
		auto jpsipi_result = jpsipi_ext.fitTo(*fullds_JPSI_shape_fix, Range("bjpsipi"), Save(), Extended());
		// FIT
		plot_mcfit(*ws, &jpsipi_ext, fullds_JPSI_shape_fix, "./results/BP/InclusiveMC_JPsipi_fit.pdf", RooFit::Name("MCFit"), NormRange("bjpsipi"), LineColor(kRed), LineStyle(1), LineWidth(2));
		ws->import(*jpsipi);
		fix_parameters(*ws, "jpsipi" );
		//[END] FIX SHAPE (J/Psi pi) 

		//[START] FIT inclusiveMC SIGNAL 
  		RooRealVar n_signal_np("n_signal_np", "n_signal_np", 1000, 0., 150000); 
		RooExtendPdf signal_ext("signal_ext", "extended signal pdf", *signalnp, n_signal_np);
		// FIT
		mass->setRange("bmc", 5.18, 5.38);
		auto signal_result = signal_ext.fitTo(*ds_sig, Range("bmc"), Save(), Extended());
		// FIT
		plot_mcfit(*ws, &signal_ext, ds_sig, "./results/BP/InclusiveMC_Signal_fit.pdf", RooFit::Name("MCFit"), Range("bmc"), LineColor(kRed),LineStyle(1), LineWidth(2));
		ws->import(*signalnp);
		fix_parameters(*ws, "signalnp");
		//[END] FIT inclusiveMC SIGNAL 

		// Fix the ratio of jpsipi to signal
		RooRealVar jpsipi_to_signal_ratio("jpsipi_to_signal_ratio", "jpsipi_to_signal_ratio",0.05, 0, 1);
  		jpsipi_to_signal_ratio.setVal(n_jpsipi_ext.getVal() / n_signal_np.getVal());
  		cout << "jpsipi_to_signal_ratio " << jpsipi_to_signal_ratio.getVal()<< endl;
  		jpsipi_to_signal_ratio.setConstant();
		ws->import(jpsipi_to_signal_ratio);
		// Fix the ratio of jpsipi to signal
						}
	// FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp

	//BIN ANALYSIS START
	//BIN ANALYSIS START
	for(int i=0;i<_nBins;i++){
		_count++;
		TCanvas* c= new TCanvas(Form("c%d",_count),"",600,550);
		TCanvas* cMC= new TCanvas(Form("cMC%d",_count),"",600,600);
		
		RooDataSet* ds_cut;
		if(doubly==0) {
			if(varExp == "Bpt"){
				ds_cut = new RooDataSet(Form("ds_cut%d",_count),"",ds,RooArgSet(*mass, *pt, *y, *trackSelection),Form("(Bpt>%f && Bpt < %f)&&((Bpt < 10 &&  abs(By) > 1.5 ) || (Bpt > 10))",_ptBins[i] , _ptBins[i+1]));
				var_mean_av[i] = ds_cut->mean(*pt);}     	
			
			else if(varExp == "By"){
				ds_cut = new RooDataSet(Form("ds_cut%d", _count),"", ds, RooArgSet(*mass, *pt, *y, *trackSelection), Form("By>%f && By< %f", _ptBins[i], _ptBins[i+1]));
				var_mean_av[i] = ds_cut->mean(*pt);}

			else if(varExp == "nMult"){
				ds_cut = new RooDataSet(Form("ds_cut%d", _count),"", ds, RooArgSet(*mass, *pt, *y, *trackSelection), Form("nMult>%f && nMult< %f", _ptBins[i], _ptBins[i+1]));
				var_mean_av[i] = ds_cut->mean(*pt);}
  						}
		
		if(doubly==1)ds_cut = new RooDataSet(Form("ds_cut%d",_count),"",ds, RooArgSet(*mass, *pt, *y, *nMult), Form("%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto)); 
		if(doubly==2)ds_cut = new RooDataSet(Form("ds_cut%d",_count),"",ds, RooArgSet(*mass, *pt, *y), Form("%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto)); 

		RooDataSet* dsMC_cut;
		if(doubly==0) 	{		
			if(varExp == "Bpt"){dsMC_cut = new RooDataSet(Form("dsMC_cut%d",_count), "", dsMC,
            RooArgSet(*mass, *pt, *y, *trackSelection), Form("(%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f)&&((Bpt < 10 &&  abs(By) > 1.5 ) || (Bpt > 10))",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto));}
        	
			else if(varExp=="By" || varExp=="nMult"){
        	dsMC_cut = new RooDataSet(Form("dsMC_cut%d", _count),"", dsMC,  RooArgSet(*mass, *pt, *y, *trackSelection), Form("%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto));}
        
		/*	else if(varExp=="nMult"){
        	dsMC_cut = new RooDataSet(Form("dsMC_cut%d", _count),"", dsMC,  RooArgSet(*mass, *pt, *y, *trackSelection), Form("%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto));}	*/       
   					}
		
		if(doubly==1) dsMC_cut = new RooDataSet(Form("dsMC_cut%d",_count),"",dsMC, RooArgSet(*mass, *pt, *y, *nMult), Form("%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto), "1"); 
		if(doubly==2) dsMC_cut = new RooDataSet(Form("dsMC_cut%d",_count),"",dsMC, RooArgSet(*mass, *pt,  *y), Form("%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto), "1"); 


    // Apply track selection cut
    std::cout << "data entries: " << ds_cut->sumEntries() << "\n";
    std::cout << "MC entries: " << dsMC_cut->sumEntries() << "\n";
    ds_cut = (RooDataSet*) ds_cut->reduce(seldata);
    dsMC_cut = (RooDataSet*) dsMC_cut->reduce(selmc);
    std::cout << "data entries: " << ds_cut->sumEntries() << "\n";
    std::cout << "MC entries: " << dsMC_cut->sumEntries() << "\n";

		// create RooDataHist
		h = new TH1D(Form("h%d",_count),"",nbinsmasshisto,minhisto,maxhisto);
		hMC = new TH1D(Form("hMC%d",_count),"",nbinsmasshisto,minhisto,maxhisto);
	
		skimtreeMC_new->Project(Form("hMC%d",_count),"Bmass",Form("%s*(%s&&%s>%f&&%s<%f)",weightmc.Data(),Form("%s&&BgenNew==23333",selmc.Data()),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1]));
	
		dh = new RooDataHist(Form("dh%d",_count),"",*mass,Import(*h));
		dhMC = new RooDataHist(Form("dhMC%d",_count),"",*mass,Import(*hMC));
		
		h->SetAxisRange(0,h->GetMaximum()*1.4,"Y");
		outputw->import(*ds);
		outputw->import(*dsMC);
		outputw->import(*dh);
		outputw->import(*dhMC);

////////// FITFITFITFITFITFITFITFITFITFITFITFIT

		cout << "Starting the fiting fucntion" << endl;
		RooFitResult* f = fit("", "", tree, c, cMC, ds_cut, dsMC_cut, dh, mass, frame, _ptBins[i], _ptBins[i+1], isMC, npfit, *ws);

////////// FITFITFITFITFITFITFITFITFITFITFITFIT


	//	validate_fit(w_val, tree, varExp, full,centmin, centmax, _ptBins[i], _ptBins[i+1]);
		//scan_significance(w_val, tree, varExp, full,centmin, centmax, _ptBins[i], _ptBins[i+1]);
		/*for(int q= 0; q < 100; q++){
			validate_fit(w_val, tree, varExp, full,q);}*/
		//datahist = frame->getHist("ds");
		//TGraphAsymmErrors* datagraph = static_cast<TGraphAsymmErrors*>(datahist);

		RooRealVar* fitYield = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index(Form("nsig%d",_count))));
		modelcurve = frame->getCurve(Form("model%d",_count));
		yield = fitYield->getVal();
		RooRealVar* BackGround = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index(Form("nbkg%d",_count))));
		MyBackground = BackGround->getVal();
		RooRealVar* width_scale = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index("scale")));
		double Myscale = width_scale->getVal();
		double Myscale_err = width_scale->getError();
		scale_vec[i] = Myscale;
		scale_vec_err_low[i] = Myscale_err;
		scale_vec_err_high[i] = Myscale_err;
		nominal_yields.push_back(yield);
		double yieldErr = fitYield->getError();
		printf("yield: %f, yieldErr: %f\n", yield, yieldErr);
		double _ErrCor=1;
		yieldErr = yieldErr*_ErrCor;
		yield_vec[i]=yield;
		yield_vec_err_low[i]=yieldErr;
		yield_vec_err_high[i]=yieldErr;
		yield_vec[i]=yield_vec[i]/(_ptBins[i+1]-_ptBins[i]);
		yield_vec_err_low[i]=yield_vec_err_low[i]/(_ptBins[i+1]-_ptBins[i]);
		yield_vec_err_high[i]=yield_vec_err_high[i]/(_ptBins[i+1]-_ptBins[i]);
		hori_av_low[i] = var_mean_av[i]-_ptBins[i];
		hori_av_high[i] = _ptBins[i+1]-var_mean_av[i];

////////////////////////////	
//Resolution MC
		RooRealVar* sigma1 = static_cast<RooRealVar*>(f->constPars().at(f->constPars().index(Form("sigma1%d", _count))));
		double Mysigma1 = sigma1->getVal();
		double Mysigma1_err = sigma1->getError();
		RooRealVar* sigma2 = static_cast<RooRealVar*>(f->constPars().at(f->constPars().index(Form("sigma2%d", _count))));
		double Mysigma2 = sigma2->getVal();
		double Mysigma2_err = sigma2->getError();
		RooRealVar* weight = static_cast<RooRealVar*>(f->constPars().at(f->constPars().index(Form("sig1frac%d", _count))));
		double Myweight  = weight->getVal();
		double Myweight_err = weight->getError();
		double scale_err_rel = Myscale_err / Myscale;
		double resol = sqrt(Myweight * pow(Mysigma1, 2) + (1 - Myweight) * pow(Mysigma2, 2)) * Myscale ;
		double resol_err = scale_err_rel * resol;
		resol_vec[i] = resol;
		resol_vec_err_low[i] = resol_err;
		resol_vec_err_high[i] = resol_err;
//Resolution 
		std::vector<double> aa;
		double a=yield_vec_err_low[i]/yield_vec[i]*100;
		aa.push_back(a);

		if(fitOnSaved == 0){
			TH1D* htest = new TH1D(Form("htest%d",_count),"",nbinsmasshisto,minhisto,maxhisto);
			TString sideband = "(abs(Bmass-5.367)>0.2&&abs(Bmass-5.367)<0.3";
			std::cout<<"yield bkg sideband: "<<htest->GetEntries()<<std::endl;
		}
		if(varExp!="nMult"){
			hPt->SetBinContent(i+1,yield/(_ptBins[i+1]-_ptBins[i]));
			hPt->SetBinError(i+1,yieldErr/(_ptBins[i+1]-_ptBins[i]));
		}
		else{
			hPt->SetBinContent(i+1,yield/(_ptBins[1] - _ptBins[0]));
			hPt->SetBinError(i+1,yieldErr/(_ptBins[1] - _ptBins[0]));
		}
		if(f->floatParsFinal().index(Form("nsig%d",_count)) != -1){
			RooRealVar* fitMean = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index(Form("mean%d",_count))));
			hMean->SetBinContent(i+1,fitMean->getVal());
			hMean->SetBinError(i+1,fitMean->getError());  
		}
	////////////////////////////

		TLatex* tex_pt;
		TLatex* tex_nMult;
	  	TLatex* tex_y;
		TLatex* tex_y1;
		TLatex* tex_y11;
		TLatex* tex_y2;
		
	if(varExp=="Bpt"){
      //for the paper run these
      if (drawLegend) {
        tex_pt = new TLatex(0.65,0.4,Form("%d < p_{T} < %d GeV/c",(int)_ptBins[i],(int)_ptBins[i+1]));
        tex_y = new TLatex(0.65,0.34,"p_{T} > 10 GeV/c : |y| < 2.4");
        tex_y2 = new TLatex(0.65,0.28,"p_{T} < 10 GeV/c : 1.5 < |y| < 2.4");
        tex_y1 = new TLatex(0.65,0.34,"1.5 < |y| < 2.4");
        tex_y11 =new TLatex(0.65,0.34,"|y| < 2.4");
        tex_nMult = new TLatex(0.21,0.62,"0 < nTrks < 100");
      } else {
        //for the AN run these
        tex_pt = new TLatex(0.65,0.8,Form("%d < p_{T} < %d GeV/c",(int)_ptBins[i],(int)_ptBins[i+1]));
        tex_y = new TLatex(0.65,0.74,"p_{T} > 10 GeV/c : |y| < 2.4");
        tex_y2 = new TLatex(0.65,0.68,"p_{T} < 10 GeV/c : 1.5 < |y| < 2.4");
        tex_y1 = new TLatex(0.65,0.74,"1.5 < |y| < 2.4");
        tex_y11 = new TLatex(0.65,0.74,"|y| < 2.4");
        tex_nMult = new TLatex(0.21,0.62,"0 < nTrks < 100");
      }
		}

if(varExp=="By"){
      //for the paper run these
      if (drawLegend) {
        tex_pt = new TLatex(0.55,0.4,"0 < p_{T} < 100 GeV/c");
        tex_y = new TLatex(0.55,0.34,Form("%d < y < %d ",(int)_ptBins[i],(int)_ptBins[i+1]));
        tex_nMult = new TLatex(0.21,0.62,"0 < nTrks < 100");
      } else {
        //fr the AN run these
        tex_pt = new TLatex(0.55,0.8,"0 < p_{T} < 100 GeV/c");
        tex_y = new TLatex(0.55,0.74,Form("%d < y < %d ",(int)_ptBins[i],(int)_ptBins[i+1]));
        tex_nMult = new TLatex(0.21,0.62,"0 < nTrks < 100");
      }
		}
if(varExp=="nMult"){ 
	//for the paper run these
      if (drawLegend) {
      	tex_nMult = new TLatex(0.21,0.62,Form("%d < nTrks < %d",(int)_ptBins[i],(int)_ptBins[i+1]));
        tex_pt = new TLatex(0.55,0.4,"0 < p_{T} < 100 GeV/c");
        tex_y = new TLatex(0.55,0.34,"|y| < 2.4");
      } else {
        //fr the AN run these
        tex_nMult = new TLatex(0.21,0.62,Form("%d < nTrks < %d",(int)_ptBins[i],(int)_ptBins[i+1]));
        tex_pt = new TLatex(0.55,0.8,"0 < p_{T} < 100 GeV/c");
        tex_y = new TLatex(0.55,0.74,"|y| < 2.4");
      }
	}

		tex_pt->SetNDC();
		tex_pt->SetTextFont(42);
		tex_pt->SetTextSize(0.045);
		tex_pt->SetLineWidth(2);
		tex_pt->Draw();
		tex_nMult->SetNDC();
		tex_nMult->SetTextFont(42);
		tex_nMult->SetTextSize(0.045);
		tex_nMult->SetLineWidth(2);
		//tex_nMult->Draw();
		tex_y->SetNDC();
		tex_y->SetTextFont(42);
		tex_y->SetTextSize(0.045);
		tex_y->SetLineWidth(2);
	
	if (varExp=="Bpt"){
		tex_y1->SetNDC();
		tex_y1->SetTextFont(42);
		tex_y1->SetTextSize(0.045);
		tex_y1->SetLineWidth(2);
		tex_y11->SetNDC();
		tex_y11->SetTextFont(42);
		tex_y11->SetTextSize(0.045);
		tex_y11->SetLineWidth(2);
		tex_y2->SetNDC();
		tex_y2->SetTextFont(42);
		tex_y2->SetTextSize(0.045);
		tex_y2->SetLineWidth(2);

		if(_ptBins[i] >= 10){tex_y11->Draw();}
		else if(_ptBins[i+1]==50 || _ptBins[i+1]==60){
			tex_y->Draw();
			tex_y2->Draw();
			}
		else{tex_y1->Draw();}
	} 
	else{tex_y->Draw();}

	CMS_lumi(c,19011,0);
	c->Update();
	TLatex* texB;
	if(tree=="ntphi"){ texB = new TLatex(0.21,0.82, "B^{0}_{s}");}
	if(tree=="ntKp"){ texB = new TLatex(0.21,0.8, "B^{#pm}");}
	texB->SetNDC();
	texB->SetTextFont(62);
	texB->SetTextSize(0.065);
	texB->SetLineWidth(2);
	texB->Draw();

		/*TLatex *lat = new TLatex();
		lat->SetNDC();
		lat->SetTextSize(0.035);
		lat->DrawLatex(0.64,0.85,Form("S = %.1f", yield));
		lat->DrawLatex(0.64,0.80,Form("S_err = %.1f", yieldErr));		
		lat->DrawLatex(0.64,0.75,Form("B = %.1f", bkgd));
		//lat->DrawLatex(0.48,0.70,Form("Significance: S/#sqrt{S+B} = %.1f", Significance));
		lat->DrawLatex(0.64,0.70,Form("Significance: %.1f", real_significance));*/

		c->SaveAs(Form("%s%s/%s_%s_%s_%d_%d_cutY%d_",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),varExp.Data(),(int)_ptBins[i],(int)_ptBins[i+1], doubly)+tree+".pdf");
		cMC->SaveAs(Form("%s%s/%s_%s_%s_%d_%d_cutY%d_",outplotf.Data(),_prefix.Data(),"mc",_isPbPb.Data(),varExp.Data(), (int)_ptBins[i], (int)_ptBins[i+1], doubly)+tree+".pdf");

		RooCurve* modelcurve_back = new RooCurve();
		std::vector<double> back_variation; 
		std::vector<double> back_err;  
		RooCurve* modelcurve_signal = new RooCurve();
		std::vector<double> signal_variation; 
		std::vector<double> signal_err;
		std::vector<double> general_err;  
		double max_signal=0.;
		double max_back=0.;
		double full_err=0;
	
		if(syst_study==1){

				for(int j=0; j<background.size(); j++){
				RooFitResult* f_back = fit("background", background[j], tree, c, cMC, ds_cut, dsMC_cut, dh, mass, frame, _ptBins[i], _ptBins[i+1], isMC, npfit, *ws);
				
				texB->Draw();
				tex_pt->Draw();
				if (varExp=="Bpt"){
					if(_ptBins[i] >= 10){tex_y11->Draw();}
					else if(_ptBins[i+1]==50){
							tex_y->Draw();
							tex_y2->Draw();}
					else{tex_y1->Draw();}
				}else{tex_y->Draw();}
				CMS_lumi(c,19011,0);
				c->Update();
				c->SaveAs(Form("%s/%s_%s_%s_%d_%d_%s_cutY%d_", outplotf.Data(), _isMC.Data(), _isPbPb.Data(), varExp.Data(),(int)_ptBins[i],(int)_ptBins[i+1],background[j].c_str(), doubly)+tree+".pdf");
			
				modelcurve_back = frame->getCurve(Form("model%d",_count));
				RooRealVar* fitYield_back = static_cast<RooRealVar*>(f_back->floatParsFinal().at(f_back->floatParsFinal().index(Form("nsig%d",_count))));
				back_variation.push_back(fitYield_back->getVal());
				back_err.push_back(abs(((yield-fitYield_back->getVal())/yield)*100));
				if(abs(((yield-fitYield_back->getVal())/yield)*100)>max_back) max_back=abs(((yield-fitYield_back->getVal())/yield)*100);
			}

			general_err.push_back(max_back);

			for(int j=0; j<signal.size(); j++){
				RooFitResult* f_signal = fit("signal", signal[j], tree, c, cMC, ds_cut, dsMC_cut, dh, mass, frame, _ptBins[i], _ptBins[i+1], isMC, npfit, *ws);
				
				texB->Draw();
				tex_pt->Draw();
				if (varExp=="Bpt"){
					if(_ptBins[i] >= 10){tex_y11->Draw();}
					else if(_ptBins[i+1]==50){
							tex_y->Draw();
							tex_y2->Draw();}
					else{tex_y1->Draw();}
				} else{tex_y->Draw();}
				CMS_lumi(c,19011,0);
				c->Update();

				cMC->SaveAs(Form("%s%s/%s_%s_%s_%d_%d_%s_cutY%d_",outplotf.Data(),_prefix.Data(),"mc",_isPbPb.Data(),varExp.Data(), (int)_ptBins[i], (int)_ptBins[i+1],signal[j].c_str(), doubly)+tree+".pdf");
				c->SaveAs(Form("%s%s/%s_%s_%s_%d_%d_%s_cutY%d_",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),varExp.Data(),(int)_ptBins[i],(int)_ptBins[i+1],signal[j].c_str(), doubly)+tree+".pdf");
				modelcurve_signal = frame->getCurve(Form("model%d",_count));
				RooRealVar* fitYield_signal = static_cast<RooRealVar*>(f_signal->floatParsFinal().at(f_signal->floatParsFinal().index(Form("nsig%d",_count))));
				signal_variation.push_back(fitYield_signal->getVal());
				signal_err.push_back(abs(((yield-fitYield_signal->getVal())/yield)*100));
				if(abs(((yield-fitYield_signal->getVal())/yield)*100)>max_signal) max_signal=abs(((yield-fitYield_signal->getVal())/yield)*100);
			}
			general_err.push_back(max_signal);
			full_err=sqrt(max_back*max_back+max_signal*max_signal);
			general_err.push_back(full_err);
			stat_error.push_back(aa);
			background_syst.push_back(back_variation);
			signal_syst.push_back(signal_variation);
			back_syst_rel_values.push_back(back_err);
			sig_syst_rel_values.push_back(signal_err);
			general_syst.push_back(general_err);
			yield_vec_systerr_low[i] = general_err[2] / 100 * yield_vec[i];
			yield_vec_systerr_high[i] = general_err[2] / 100 * yield_vec[i];

		}

	//validate_fit(outputw,1,Form("model%d",_count),Form("nsig%d",_count));
	}

	hMean->Write();
	hPt->Write();
	outf->Close();	

	string Path;
	if(tree == "ntphi"){ Path = "./filesbs/";}
	else if (tree == "ntKp"){Path = "./filesbp/";}
	std::ofstream myfile;
	myfile.open (Path + "systematics_" + tree.Data() + ".txt");

	if(syst_study==1){ 
		for(int i=0; i<_nBins; i++){
			std::cout<<"pt bin = "<<_ptBins[i]<<" "<<_ptBins[i+1]<<std::endl;
			std::cout<<"nominal yield = "<<nominal_yields[i]<<std::endl;
			std::cout<<"count= "<<i<<std::endl;

			for(int j=0; j<background.size(); j++){
				std::cout<<" back sys in bin "<<i<<" ; with pdf "<< background[j] << " ="<<background_syst[i][j]<<std::endl;
				std::cout<<" back sys in bin "<<i<<" ; with pdf "<< background[j] << " ="<<back_syst_rel_values[i][j]<<" % "<<std::endl;
				myfile<<" back sys in bin "<<i<<" ; with pdf "<< background[j] << " ="<<back_syst_rel_values[i][j]<<" % "<<std::endl;
			}
			for(int j=0; j<signal.size(); j++){
				std::cout<<" signal sys in bin "<<i<<" ; with pdf "<< signal[j] << " ="<<signal_syst[i][j]<<std::endl;
				std::cout<<" signal sys in bin "<<i<<" ; with pdf "<< signal[j] << " ="<<sig_syst_rel_values[i][j]<<" % "<<std::endl;
				myfile<<" signal sys in bin "<<i<<" ; with pdf "<< signal[j] << " ="<<sig_syst_rel_values[i][j]<<" % "<<std::endl;
			}
		}

	}

	myfile.close();

	TCanvas* cPt =  new TCanvas("cPt","",600,600);
	cPt->SetLogy();
	hPt->SetXTitle("B_{s} p_{T} (GeV/c)");
	hPt->SetYTitle("Uncorrected dN(B_{s})/dp_{T}");
	hPt->Draw();

	std::vector<std::string> labels_back = {"Linear", "2nd Poly", "mass range" };
	std::vector<std::string> col_name_back;
	std::vector<std::string> labels_signal = {"Triple Gaussian", "Fixed Mean", "CB+Gaussian", "Double CB"};
	std::vector<std::string> labels_general = {"Background", "Signal", "Total"};
	std::vector<std::string> labels_general_stat = {"Statistical error"};
	std::vector<std::string> col_name_general;
	std::vector<std::string> col_name_signal;
	std::vector<std::string> col_name_general_stat;

	string name;
	int m=_nBins;
	
	col_name_back.push_back("Background Model");
	col_name_signal.push_back("Signal Model");
	col_name_general.push_back("Systematic Source");
	col_name_general_stat.push_back(" ");

	if(varExp=="Bpt"){name="$<p_T<$";} else if(varExp=="By"){name="$<y<$";} else if(varExp=="nMult"){name="$<nTrks<$";}
	for(int i=0;i<m;i++){
		std::ostringstream collabel;
		std::ostringstream siglabel;
		std::ostringstream genlabel;
		std::ostringstream statlabel;
		collabel<<_ptBins[i]<<name<<_ptBins[i+1];
		siglabel<<_ptBins[i]<<name<<_ptBins[i+1];
		genlabel<<_ptBins[i]<<name<<_ptBins[i+1];
		statlabel<<" "<<_ptBins[i]<<name<<_ptBins[i+1];
		std::string label1 = collabel.str();
		col_name_back.push_back(label1);
		col_name_signal.push_back(label1);
		col_name_general.push_back(label1);
		col_name_general_stat.push_back(label1);
	}
	
	if(syst_study==1 && full==0){
		latex_table(Path + "background_systematics_table_"+std::string (varExp.Data())+"_"+std::string (tree.Data()), _nBins+1,  (int)(1+background.size()),  col_name_back,labels_back,back_syst_rel_values, "Background PDF Systematic Errors");
		latex_table(Path + "signal_systematics_table_"+std::string (varExp.Data())+"_"+std::string (tree.Data()), _nBins+1, (int)(1+signal.size()),    col_name_signal, labels_signal,sig_syst_rel_values, "Signal PDF Systematic Errors");
		latex_table(Path + "general_systematics_table_"+std::string (varExp.Data())+"_"+std::string (tree.Data()),  _nBins+1, 4 , col_name_general, labels_general, general_syst, "Overall PDF Variation Systematic Errors");	
		latex_table(Path + "Statistical_error_table_"+std::string (varExp.Data())+"_"+std::string (tree.Data()),  _nBins+1, 2 , col_name_general_stat, labels_general_stat, stat_error, "Statistical error");	
	}

	cout << "Final Background = " << MyBackground << endl;
	cout << "Final Yield = " << yield << endl;

// Differential plot part starts
	gSystem->mkdir("./results/Graphs",true); 
	TFile *ratio_f= new TFile(Form("./results/%s_%s_ratio.root",tree.Data(),varExp.Data()),"recreate");
	ratio_f->cd();
	
	 TCanvas c_diff;
	 TMultiGraph* mg = new TMultiGraph();
	 TLegend *leg_d = new TLegend(0.7,0.7,0.9,0.9);

	 TGraphAsymmErrors* gr_staterr = new TGraphAsymmErrors(_nBins,var_mean_av,yield_vec,hori_av_low,hori_av_high,yield_vec_err_low,yield_vec_err_high);
	 gr_staterr->SetLineColor(1); 
	 mg->Add(gr_staterr);

	 if(syst_study==1){
		TGraphAsymmErrors* gr_systerr = new TGraphAsymmErrors(_nBins,var_mean_av,yield_vec,nullptr,nullptr,yield_vec_systerr_low,yield_vec_systerr_high);
		gr_systerr->SetLineColor(2);
		mg->Add(gr_systerr,"syst");
		leg_d->AddEntry(gr_systerr, "Systematic Uncertainty", "e");
	}
	 if(varExp == "By"){
		 mg->GetXaxis()->SetTitle("Rapidity (y)");
		 mg->GetYaxis()->SetTitle("dN_{S}/dy");
		 mg->GetXaxis()->SetLimits(-2.4 ,2.4);
	 }
	 if(varExp == "Bpt"){
		 mg->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		 mg->GetYaxis()->SetTitle("dN_{S}/dp_{T}");
		if (tree == "ntKp"){ mg->GetXaxis()->SetLimits(0 ,80); }
		if (tree == "ntphi"){ mg->GetXaxis()->SetLimits(0 ,60); }
	 }
	 if(varExp == "nMult"){
		 mg->GetXaxis()->SetTitle("Multiplicity (Mult)");
		 mg->GetYaxis()->SetTitle("dN_{S}/dMult");
		 mg->GetXaxis()->SetLimits(0, 110);
	 }

	 mg->Draw("ap");
	 //mg->SetTitle("Differential Signal Yield");  
	 
	 leg_d->AddEntry(gr_staterr, "Statistical Uncertainty", "e");
	 //leg_d->AddEntry(grs, "Systematic Uncertainty", "e");
	 leg_d->SetBorderSize(0);
	 leg_d->SetFillStyle(0);
	 leg_d->SetTextSize(0);
	 leg_d->Draw();

	 const char* pathc =Form("./results/Graphs/raw_yield_%s_%s.png",tree.Data(),varExp.Data());
	 c_diff.SaveAs(pathc);
	 ratio_f->Close();
// Differential plot part ends

// Parameters vs variables part starts
	double scale_max = 0;
	for(int i = 0; i < _nBins; i++){
		if(scale_vec[i] > scale_max){scale_max = scale_vec[i];}
									}

	 TCanvas c_par;
	 TMultiGraph* mg_par = new TMultiGraph();

	 TGraphAsymmErrors* gr_scale = new TGraphAsymmErrors(_nBins,var_mean_av,scale_vec,hori_av_low,hori_av_high,scale_vec_err_low,scale_vec_err_high);
	 gr_scale->SetLineColor(1); 
	
	 if(varExp == "By"){
		 mg_par->GetXaxis()->SetTitle("Rapidity (y)");
		 mg_par->GetYaxis()->SetTitle("Mass resolution scale factor");
		 mg_par->GetXaxis()->SetLimits(-2.4 ,2.4);
	 }
	 if(varExp == "Bpt"){
		 mg_par->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		 mg_par->GetYaxis()->SetTitle("Mass resolution scale factor");
		 if (tree == "ntKp"){ mg->GetXaxis()->SetLimits(0 ,80); }
		 if (tree == "ntphi"){ mg->GetXaxis()->SetLimits(0 ,60); }
	 }
	 if(varExp == "nMult"){
		 mg_par->GetXaxis()->SetTitle("Multiplicity (Mult)");
		 mg_par->GetYaxis()->SetTitle("Mass resolution scale factor");
		 mg_par->GetXaxis()->SetLimits(0, 110);
	 }
	 mg_par->Add(gr_scale);
	 mg_par->GetYaxis()->SetRangeUser(0,scale_max*1.4);
	 mg_par->Draw("ap");
/*	 
     TLegend *leg_par = new TLegend(0.7,0.7,0.9,0.9);
	 leg_par->AddEntry(gr_scale, "Scale", "e");
	 //leg_d->AddEntry(grs, "Systematic Uncertainty", "e");
	 leg_par->SetBorderSize(0);
	 leg_par->SetFillStyle(0);
	 leg_par->SetTextSize(0);
	 leg_par->Draw();
*/
	 const char* pathc_par =Form("./results/Graphs/parameters_variation_%s_%s.png",tree.Data(),varExp.Data()); 
	 c_par.SaveAs(pathc_par);
//Parameters vs variables part ends

//Resolution plot part starts
	 double resol_max = 0;
	 double resol_min = 100000;
	for(int i = 0; i < _nBins; i++){
		if(resol_vec[i] > resol_max){resol_max = resol_vec[i];}
		if(resol_vec[i] < resol_min){resol_min = resol_vec[i];}
									}
	 TCanvas c_resol;
	 TMultiGraph* mg_resol = new TMultiGraph();

	 TGraphAsymmErrors* gr_resol = new TGraphAsymmErrors(_nBins,var_mean_av,resol_vec,hori_av_low,hori_av_high,resol_vec_err_low,resol_vec_err_high);
	 gr_resol->SetLineColor(1); 
	
	 if(varExp == "By"){
		 mg_resol->GetXaxis()->SetTitle("Rapidity (y)");
		 mg_resol->GetYaxis()->SetTitle("Resolution");
		 mg_resol->GetXaxis()->SetLimits(-2.4 ,2.4);

	 }
	 if(varExp == "Bpt"){
		 mg_resol->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		 mg_resol->GetYaxis()->SetTitle("Resolution");
		 if (tree == "ntKp"){ mg_resol->GetXaxis()->SetLimits(0 ,80); }
		 if (tree == "ntphi"){ mg_resol->GetXaxis()->SetLimits(0 ,60); }
	 }
	 if(varExp == "nMult"){
		 mg_resol->GetXaxis()->SetTitle("Multiplicity (Mult)");
		 mg_resol->GetYaxis()->SetTitle("Resolution");
		 mg_resol->GetXaxis()->SetLimits(0, 110);
		 
	 }
	 mg_resol->GetYaxis()->SetRangeUser(resol_min*0.6, resol_max*1.4);
	 mg_resol->Add(gr_resol);
	 mg_resol->Draw("ap");

	 const char* pathc_resol =Form("./results/Graphs/resolution_%s_%s.png",tree.Data(),varExp.Data()); 
	 c_resol.SaveAs(pathc_resol);

//Resolution plot part ends

}


template<typename... Targs>
void plot_mcfit(RooWorkspace& w, RooAbsPdf* model, RooDataSet* ds, TString plotName, Targs... options) {
  
  TCanvas can_mc;
  can_mc.cd();
  TPad *p1 = new TPad("p1","p1",0.05,0.05,0.9,0.9);
  p1->SetBorderMode(1); 
  p1->SetFrameBorderMode(0); 
  p1->SetBorderSize(2);
  p1->SetBottomMargin(0.10);
  p1->Draw(); 
  p1->cd();

  RooRealVar Bmass = *(w.var("Bmass"));
  Bmass.setRange("bmass", 5.0, 6.0);
  RooPlot* massframe = Bmass.frame(Title(" "));
  ds->plotOn(massframe, MarkerSize(0.8));
  model->plotOn(massframe, options...);
  model->paramOn(massframe, Layout(0.55, 1, 0.9), "", Format("NEU", AutoPrecision(1)) ) ;
  massframe->GetXaxis()->SetRangeUser(5.15,5.45);
  if (plotName == "./results/BP/InclusiveMC_JPsipi_fit.pdf") {massframe->GetXaxis()->SetRangeUser(5.1,5.9);}
  massframe->Draw();
  can_mc.SaveAs(plotName);
}

void read_samples(RooWorkspace& w, std::vector<TString> label, TString fName, TString treeName, TString sample){
  TFile* fin = new TFile(fName);
  TTree* t1;
  t1 = (TTree*) fin->Get(treeName);

  RooArgList arg_list ("arg_list");
  // read the fitting variable
  arg_list.add(*(w.var("Bmass")));

  // read additional variables
  for(auto lab : label){arg_list.add(*(w.var(lab)));}

  RooDataSet* data_s = new RooDataSet(sample, sample, t1, arg_list);
  cout << "input filename = " << fName << "; entries: " << data_s->sumEntries() << endl;
  w.import(*data_s, Rename(sample));
}
