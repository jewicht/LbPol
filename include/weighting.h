#include "TVectorD.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooChebychev.h"
#include "RooProdPdf.h"
#include "RooCmdConfig.h"
#include "RooConstraintSum.h"
#include "RooNLLVar.h"
#include "RooProduct.h"
#include "RooAddition.h"
#include "RooMinimizer.h"

#include "calcacceptanceclass.h"

void computeamplitudes(RooRealVar&, RooRealVar&, RooRealVar&, RooRealVar&, double, double, double, double, const TMatrixDSym&);
TMatrixDSym getVCV(const TMatrixDSym&, const TMatrixDSym&);
TMatrixDSym makesymm(const TMatrixD&);


void deletecolumn(RooDataSet*&, const RooRealVar&);
void copyvar(RooDataSet*, const TString&, RooRealVar&, double scale=1.);

RooAbsPdf* factorizedacc(RooRealVar&, RooRealVar&, RooRealVar&, const TString&, RooArgSet*&, int, int, int, const TString& filename="");
void addacceptanceweight(RooDataSet*, RooRealVar&, bool, bool, calcacceptanceclass* mycalcacceptanceclass_ll, calcacceptanceclass* mycalcacceptanceclass_dd);
void addacceptanceweightfactpdf(RooDataSet*, RooRealVar&, bool, bool, RooAbsPdf* mycalcacceptanceclass_ll, RooAbsPdf* mycalcacceptanceclass_dd);
void multiplyweights(RooDataSet*, RooRealVar&, const TString&, const TString&, bool normalizetowmass, TString = "");

double calcPsample(RooDataSet*, RooRealVar&, RooAbsPdf*, RooAbsPdf*, double);
double getsumwsumw2(RooDataSet*, const TString& wname="");

RooFitResult* myfitto(RooDataSet*, RooAbsPdf*, const RooCmdArg&, int elevel=1, bool useminos=true);
void plotlikelihoodprofiles(const TString&, const TCanvas&, const RooArgList&, RooDataSet*, RooAbsPdf*, const RooCmdArg&);
void plotlikelihoodscans(const TString&, const TCanvas&, const RooArgList&, RooDataSet*, RooAbsPdf*, const RooCmdArg&);
void plotlikelihoodscans(const TString&, const TCanvas&, RooRealVar&, double, double, const RooArgList&, RooDataSet*, RooAbsPdf*, const RooCmdArg&);

