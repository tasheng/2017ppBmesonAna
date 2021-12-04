#ifndef _TMVACLASSIFICATION_H_
#define _TMVACLASSIFICATION_H_

#include <string>
#include <vector>
#include <TString.h>

#include "xjjcuti.h"
#include "ntuple.h"

namespace mytmva
{
  struct tmvavar
  {
    std::string varname;
    std::string vartex;
    std::string var;
    std::string cutsign;
    float varmin;
    float varmax;
    tmvavar(const std::string& varname_, const std::string& var_, const std::string& cutsign_, const std::string& vartex_, const float& varmin_, const float& varmax_) 
      : varname(varname_), var(var_), cutsign(cutsign_), vartex(vartex_), varmin(varmin_), varmax(varmax_) { ; }
  };

  const std::vector<mytmva::tmvavar> varlist = {
 
    mytmva::tmvavar("Btrk1Pt"  , "Btrk1Pt"                                                                                        , "FMax", "#pi_{1} p_{T} (GeV/c)"                , 0   , 10) , // 8
    mytmva::tmvavar("Btrk2Pt"  , "Btrk2Pt"                                                                                        , "FMax", "#pi_{2} p_{T} (GeV/c)"                , 0   , 10) , // 9	 
	mytmva::tmvavar("Trk1DCAz"  , "Trk1DCAz := abs(Btrk1Dz1/Btrk1DzError1)"                                                                                        , "FMax", "#Eta_{1} "                , -3   , 3) , // 8
    mytmva::tmvavar("Trk2DCAz"  , "Trk2DCAz := abs(Btrk2Dz1/Btrk2DzError1)"                                                                                        , "FMax", "#Eta_{2}"                , -3   , 3) , // 9	
    mytmva::tmvavar("Trk1DCAxy"   , "Trk1DCAxy := abs(Btrk1Dxy1/Btrk1DxyError1)",                                                                 "FMin", "#DeltaR(#pi_{1},J/#psi)"              , 0   , 0.5), // 1
    mytmva::tmvavar("Trk2DCAxy"   , "Trk2DCAxy := abs(Btrk2Dxy1/Btrk2DxyError1)"              ,                                    "FMin",      "#DeltaR(#pi_{2},J/#psi)"                     , 0   , 0.5), // 2
    mytmva::tmvavar("MassDis"   , "MassDis := abs(Btktkmass-1.019455)"                                                           , "FMin", "Q (GeV/c^{2})"                        , -0.2   , 0.2), // 3
	mytmva::tmvavar("dls"    , "dls := BsvpvDistance/BsvpvDisErr"                                                 , "FMax", "l_{xyz}/#sigma(l_{xyz})"              , 0   , 10) , // 6
	mytmva::tmvavar("Balpha"   , "Balpha"                                                                                         , "FMin", "#alpha"                               , 0   , 3.2), // 4
    mytmva::tmvavar("dls2D"    , "dls2D := Bd0"                                                                            , "FMax", "l_{xy}/#sigma(l_{xy})"                , 0   , 10) , // 7
	mytmva::tmvavar("cos(Bdtheta)" , "cos(Bdtheta)"                                                                , "FMax", "cos(#theta)"                          , -1  , 1)  , 
	mytmva::tmvavar("Bchi2cl"  , "Bchi2cl"                                                                                        , "FMax", "vertex #chi^{2} prob"                 , 0   , 1),   // 0
	mytmva::tmvavar("Btrk1Eta" , "Btrk1Eta"                                                                , "", "Trak1 Eta"                          , -2.4  , 2.4)  , 
	mytmva::tmvavar("Btrk2Eta"  , "Btrk2Eta"                                                                                        , "", "Trak 2 Eta"                 , 0   , 2.4),   // 0
	mytmva::tmvavar("Bmass"  , "Bmass"                                                                                        , "", "B+ mass"                 , 5.0   , 6.0),   // 0
	mytmva::tmvavar("BDT_pt_1_2"  , "BDT_pt_1_2"                                                                                        , "", "BDT 1 - 2 GeV/c"                 , -1.0   , 1.0),   // 0
	mytmva::tmvavar("BDT_pt_2_3"  , "BDT_pt_2_3"                                                                                        , "", "BDT 2 - 3 GeV/c"                 , -1.0   , 1.0),   // 0
	mytmva::tmvavar("BDT_pt_3_5"  , "BDT_pt_3_5"                                                                                        , "", "BDT 3 - 5 GeV/c"                 , -1.0   , 1.0),   // 0
	mytmva::tmvavar("BDT_pt_5_7"  , "BDT_pt_5_7"                                                                                        , "", "BDT 5 - 7 GeV/c"                 , -1.0   , 1.0),   // 0	
	mytmva::tmvavar("BDT_pt_7_10"  , "BDT_pt_7_10"                                                                                        , "", "BDT 7 - 10 GeV/c"                 , -1.0   , 1.0),   // 0
	mytmva::tmvavar("BDT_pt_10_15"  , "BDT_pt_10_15"                                                                                        , "", "BDT 10 - 15 GeV/c"                 , -1.0   , 1.0),   // 0
	mytmva::tmvavar("BDT_pt_15_20"  , "BDT_pt_15_20"                                                                                        , "", "BDT 15 - 20 GeV/c"                 , -1.0   , 1.0),   // 0
	mytmva::tmvavar("BDT_pt_20_30"  , "BDT_pt_20_30"                                                                                        , "", "BDT 20 - 30 GeV/c"                 , -1.0   , 1.0),   // 0	
        mytmva::tmvavar("BDT_pt_30_50"  , "BDT_pt_30_50"                                                                                        , "", "BDT 30 - 50 GeV/c"                 , -1.0   , 1.0),   // 0  
	mytmva::tmvavar("Btktkmass"  , "Btktkmass"                                                                                        , "", "BDT 20 - 50 GeV/c"                 , 1.004455   , 1.034455),   // 0	
	mytmva::tmvavar("2*Btktkmass"  , "2*Btktkmass"                                                                                        , "", "BDT 20 - 50 GeV/c"                 , 2.008910   , 2.068910),   // 0
	mytmva::tmvavar("-1*Btktkmass"  , "-1*Btktkmass"                                                                                        , "", "BDT 20 - 50 GeV/c"           , -2.068910       ,	-2.008910  ),  //0	
	mytmva::tmvavar("abs(Btktkmass-1.019455)"  , "abs(Btktkmass-1.019455)"                                                                                        , "", "BDT 20 - 50 GeV/c"                 , -0.015   , 0.015),   // 0	
	mytmva::tmvavar("2*abs(Btktkmass-1.019455)"  , "2*abs(Btktkmass-1.019455)"                                                                                        , "", "BDT 20 - 50 GeV/c"                 , -0.030   , 0.030)   // 0	
  };

  class varval
  {
  public:
    varval(mytmva::ntuple* nt) : fnt(nt), fvalid(true) { fval.clear(); fvalid = checkvarlist(); }
    varval(TTree* nttree) : fnt(0), fvalid(true) { fnt = new mytmva::ntuple(nttree); fval.clear(); fvalid = checkvarlist(); }
    float getval(std::string varname, int j) { refreshval(j); if(fval.find(varname) == fval.end()) { std::cout<<"==> "<<__FUNCTION__<<": invalid varname key "<<varname<<std::endl; return 0; } ; return fval[varname]; }
    mytmva::ntuple* getnt() { return fnt; }
    bool isvalid() { return fvalid; }

  private:
    bool fvalid;
    std::map<std::string, float> fval;
    mytmva::ntuple* fnt; //~
    void refreshval(int j)
    {
      fval["Bmass"]     = j<0?0:(float)fnt->Bmass[j];
      fval["dRtrk1"]    = j<0?0:(float)TMath::Sqrt(pow(TMath::ACos(TMath::Cos(fnt->Bujphi[j] - fnt->Btrk1Phi[j])), 2) + pow(fnt->Bujeta[j] - fnt->Btrk1Eta[j], 2));
      fval["dRtrk2"]    = j<0?0:(float)TMath::Sqrt(pow(TMath::ACos(TMath::Cos(fnt->Bujphi[j] - fnt->Btrk2Phi[j])), 2) + pow(fnt->Bujeta[j] - fnt->Btrk2Eta[j], 2));
      fval["MassDis"]    = j<0?0:(float)TMath::Abs(fnt->Btktkmass[j]-1.019455);
      fval["cos(Bdtheta)"]  = j<0?0:(float)TMath::Cos(fnt->Bdtheta[j]);
      fval["dls"]     = j<0?0:(float)fnt->BsvpvDistance[j]/fnt->BsvpvDisErr[j];
      fval["Bchi2cl"]   = j<0?0:(float)fnt->Bchi2cl[j];
      fval["Btrk1Pt"]   = j<0?0:(float)fnt->Btrk1Pt[j];
      fval["Btrk2Pt"]   = j<0?0:(float)fnt->Btrk2Pt[j];
      fval["Balpha"]    = j<0?0:(float)fnt->Balpha[j];
      fval["Trk1DCA"]  = j<0?0:(float)fnt->Btrk1Dxy[j]/fnt->Btrk1D0Err[j];
      fval["Trk2CA"]  = j<0?0:(float)fnt->Btrk2Dxy[j]/fnt->Btrk2D0Err[j];
      fval["dls2D"]     = j<0?0:(float)TMath::Abs(fnt->BsvpvDistance_2D[j]/fnt->BsvpvDisErr_2D[j]);
      fval["trkptimba"] = j<0?0:(float)TMath::Abs((fnt->Btrk1Pt[j]-fnt->Btrk2Pt[j]) / (fnt->Btrk1Pt[j]+fnt->Btrk2Pt[j]));
      fval["By"]        = j<0?0:(float)fnt->By[j];
      fval["Btrk1Eta"]  = j<0?0:(float)fnt->Btrk1Eta[j];
      fval["Btrk2Eta"]  = j<0?0:(float)fnt->Btrk2Eta[j];
	  fval["Bmass"]  = j<0?0:(float)fnt->Bmass[j];
	  fval["Btktkmass"]  = j<0?0:(float)fnt->Btktkmass[j];
	}
    bool checkvarlist() 
    {  
      refreshval(-1);
      for(auto& vn : varlist) 
        { if(fval.find(vn.varname) == fval.end()) { std::cout<<"==> "<<__FUNCTION__<<": invalid varname key "<<vn.varname<<std::endl; return false; } }
      return true;
    }
  };

}

#endif
