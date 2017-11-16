#include "RooRealVar.h"
#include "legendre3f3.h"
#include "RooLegendre3Pdfv2.h"
//include "RooChebychev3Pdf.h"
#include "RooFitResult.h"
#include "TFile.h"

#include "RooLbtoJpsiL0wAccPDF3.h"
#include "functions.h"
#include "functions-roofit.h"





//void writecoeff3(const TString&, RooRealVar*[ordermax3i][ordermax3j][ordermax3k]);

void readcoeff3_rrv(const TString&, RooRealVar*[ordermax3i][ordermax3j][ordermax3k]);

void setconstcoeff3(RooRealVar*[ordermax3i][ordermax3j][ordermax3k]);

void setconstcoeff3_maxorder(RooRealVar*[ordermax3i][ordermax3j][ordermax3k], int , int, int);
void setconstcoeff3_maxtotorder(RooRealVar*[ordermax3i][ordermax3j][ordermax3k], int);

void disableconstcoeff3_maxorder(RooRealVar*[ordermax3i][ordermax3j][ordermax3k], int, int, int);
void setmaxordercoeff3(RooRealVar*[ordermax3i][ordermax3j][ordermax3k], int, int, int);

void constornotconstcoeff3(RooRealVar*[ordermax3i][ordermax3j][ordermax3k]);


void createRRV3(bool, RooRealVar*[ordermax3i][ordermax3j][ordermax3k],TString, int, int, int, int, int, const TString& filename="");
void createRooLegendre3v2(RooRealVar&, RooRealVar&, RooRealVar&, TString, RooAbsPdf*&, RooRealVar*[ordermax3i][ordermax3j][ordermax3k]);

//void createRooChebychev3(RooRealVar&, RooRealVar&, RooRealVar&, TString, RooAbsPdf*&, RooRealVar*[ordermax3i][ordermax3j][ordermax3k], bool, int, int, int, int, int, const TString& filename="");

void createRooLbtoJpsiL0wAccPDF3(RooRealVar&, RooRealVar&, RooRealVar&, RooRealVar&, RooAbsReal&, RooRealVar&, RooAbsReal&, RooAbsReal&, TString, RooAbsPdf*&, RooRealVar*[ordermax3i][ordermax3j][ordermax3k], bool, const TString& filename="");

