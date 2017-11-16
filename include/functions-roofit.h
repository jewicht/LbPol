#ifndef FUNCTIONSROOFITH
#define FUNCTIONSROOFITH

#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCut.h"
#include "TTree.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TEventList.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TText.h"

#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "RooMCStudy.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooHist.h"
#include "RooGenericPdf.h"
#include "RooGaussian.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooPullVar.h"

#include "functions.h"
#include "rootlogon.h"
#include "progressbar.h"


bool checkweirdpulls(RooFitResult*, RooArgList*);

RooArgList* addvariation(RooDataSet*, const RooArgList&);

void th2ftolatextable(TH2F*, const TString&);

void plotprofile(const TCanvas&, const TString&, RooDataSet*, RooRealVar*, RooRealVar*);
void plotprofile(const TCanvas&, TTree*, const TString&, const TString&, const TCut& , const TString&);
void plotprofile(const TCanvas&, TTree*, const std::vector<TString>&, const std::vector<TString>&, const TCut& , const TString&);

bool isconverged(RooFitResult*, bool debug=false);

void myprintLatex(const RooArgList&, const char*, Int_t, const RooCmdArg*);


RooFitResult* getRooFitResult(const TString&, TString frname="");

void readvarsfromtfile(const RooArgList&, const TString&,  TString frname="", bool randomize=false);
Double_t findnextbin(RooDataSet*, Double_t, const TString&, Double_t, Double_t);

RooDataSet* filldataset(const TString&, const TString&, const RooArgList&, TTree*, const TCut&);
void filldataset_copytree(RooDataSet*&, const TString&, const TString&, const RooArgList&, TTree*, const TCut&, bool keep=false);

void plot(TCanvas&, const RooArgList&, const RooAbsPdf&, RooAbsData*, const TString&, bool plotpull=false, const RooArgList* components=NULL, const TString& mycut="", RooCategory* cat=NULL, bool addlhcb=false);

/* void plotwcomponents(TCanvas&, const RooRealVar&, const RooAbsPdf&, const RooArgList&, RooAbsData*, const TString&, bool plotpull=false); */
/* void plotwcomponents(TCanvas&, const RooArgList&, const RooAbsPdf&, const RooArgList&, RooAbsData*, const TString&, bool plotpull=false); */

//void plotwpull(TCanvas&, const RooRealVar&, const RooAbsPdf&, RooAbsData*, const TString&);
//void plotwpull(TCanvas&, const RooArgList&, const RooAbsPdf&, RooAbsData*, const TString&);

/* void plotwcut(const TCanvas&, const RooRealVar&, const RooAbsPdf&, RooAbsData*, const TString&, const TString&); */
/* void plotwcut(const TCanvas&, const RooArgList&, const RooAbsPdf&, RooAbsData*, const TString&, const TString&); */

/* void plotwcutwcomponents(const TCanvas&, const RooRealVar&, const RooAbsPdf&, const RooArgList&, RooDataSet*, const TString&, const TString&); */
/* void plotwcutwcomponents(const TCanvas&, const RooArgList&, const RooAbsPdf&, const RooArgList&, RooDataSet*, const TString&, const TString&); */

void plotpdf(const TCanvas&, const RooRealVar&, const RooAbsPdf&, const TString&);
void plotpdflist(const TCanvas&, const RooRealVar&, const RooArgList&, const TString&);

void plotdata(const TCanvas&, RooRealVar&, RooDataSet*, const TString&, TString mycut="", bool fitGauss=false, RooRealVar* mean=NULL, RooRealVar* sigma=NULL);
void plotdata(const TCanvas&, const RooArgList&, RooDataSet*, const TString&, TString mycut="", bool fitGauss=false, RooRealVar* mean=NULL, RooRealVar* sigma=NULL);

void plotdatavsdata2d(TCanvas&, RooRealVar*, RooRealVar*, RooDataSet*, RooDataSet*, const TString&, bool normbinning=false, bool writerootfile=true, int nbinx=0, int nbiny=0);
void plotdatavsdata(TCanvas&, const RooRealVar&, RooDataSet*, RooDataSet*, const TString&, bool normbinning=false, bool plotpull=true, bool writerootfile=true, double cutval=-666., bool gp=true);
void plotdatavsdata(TCanvas&, const RooArgList&, RooDataSet*, RooDataSet*, const TString&, bool normbinning=false, bool plotpull=true, bool writerootfile=true);

void plot_mcstudy(const TCanvas&, RooMCStudy&, const TString&, const RooRealVar&);
void plot_mcstudy(const TCanvas&, RooMCStudy&, const TString&, const RooArgList&);

void plotttree(TCanvas&, const RooArgList&, TTree*, const TCut&, const TCut&, const TString&, bool plotpull=false, TString leg1="", TString leg2="");

double sumofvarindata(RooDataSet*, const TString&);
double minofdataset(RooDataSet*, const RooRealVar&);
double maxofdataset(RooDataSet*, const RooRealVar&);

#endif
