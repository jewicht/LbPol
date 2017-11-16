#include "RooRealVar.h"
#include "RooLegendre5Pdfv2.h"
#include "legendre5f5.h"
#include <fstream>

void writecoeff5(const TString&, RooRealVar*[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m]);
void readcoeff5_rrv(const TString&, RooRealVar*[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m]);
void setconstcoeff5(RooRealVar*[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m]);
void setconstcoeff5_maxtotorder(RooRealVar*[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m], int);

void constornotconstcoeff5(RooRealVar*[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m]);

void createRooLegendre5v2(RooRealVar&, RooRealVar&, RooRealVar&, RooRealVar&, RooRealVar&, TString, RooLegendre5Pdfv2*&, RooRealVar*[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m], int, int, int ,int ,int ,int, int, const TString& filename="");
