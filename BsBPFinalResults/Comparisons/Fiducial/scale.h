#include <vector>
#include "TString.h"
#include "ROOT/RDataFrame.hxx"

using dvec = std::vector<double>;
/** Compute the MC extrapolation from smaller rapidity region to a larger one
    input: filename, treename
    output: vector of scaling factors
*/
dvec scaleFactor(TString inFile, TString inTree, dvec ptbins) {

  TString ynarrow = "By < 2.4 && By > 1.5";
  TString yall = "By < 2.4";

  ROOT::RDataFrame df(inTree.Data(), inFile.Data());
  dvec factors;
  for (auto i = 0; i < ptbins.size() - 1; ++i) {
    double ptlow = ptbins[i];
    double pthigh = ptbins[i + 1];
    TString ptcut = TString::Format("Bpt > %f && Bpt < %f", ptlow, pthigh);
    auto dpt = df.Filter(ptcut.Data());
    double ent_all = dpt.Filter(yall.Data()).Count().GetValue();
    double ent_narrow = dpt.Filter(ynarrow.Data()).Count().GetValue();
    factors.push_back(ent_all / ent_narrow);
  }
  return factors;
}
