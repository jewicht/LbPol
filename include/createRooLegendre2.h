#include "RooRealVar.h"
#include "RooLegendre2Pdf.h"
#include "RooFitResult.h"
#include "TFile.h"

#include "functions.h"
#include "functions-roofit.h"

#include "legendre2f2.h"



//void writecoeff3(const TString&, RooRealVar*[ordermax2i][ordermax2j]);

void readcoeff2_rrv(const TString&, RooRealVar*[ordermax2i][ordermax2j]);

void setconstcoeff2(RooRealVar*[ordermax2i][ordermax2j]);

void setconstcoeff2_maxorder(RooRealVar*[ordermax2i][ordermax2j], int , int, int);
void setconstcoeff2_maxtotorder(RooRealVar*[ordermax2i][ordermax2j], int);

void disableconstcoeff2_maxorder(RooRealVar*[ordermax2i][ordermax2j], int, int, int);
void setmaxordercoeff2(RooRealVar*[ordermax2i][ordermax2j], int, int, int);

void constornotconstcoeff2(RooRealVar*[ordermax2i][ordermax2j]);


void createRRV2(bool, RooRealVar*[ordermax2i][ordermax2j],TString, int, int, int, int, const TString& filename="");
void createRooLegendre2(RooRealVar&, RooRealVar&, TString, RooAbsPdf*&, RooRealVar*[ordermax2i][ordermax2j]);

/* void createRooChebychev3(RooRealVar&, RooRealVar&, RooRealVar&, TString, RooAbsPdf*&, RooRealVar*[ordermax2i][ordermax2j], bool, int, int, int, int, int, const TString& filename=""); */

/* void createRooLbtoJpsiL0wAccPDF3(RooRealVar&, RooRealVar&, RooRealVar&, RooRealVar&, RooAbsReal&, RooRealVar&, RooAbsReal&, RooAbsReal&, TString, RooAbsPdf*&, RooRealVar*[ordermax2i][ordermax2j], bool, const TString& filename=""); */

