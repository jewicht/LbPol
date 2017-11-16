#include "stdlib.h"
#include "TFile.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooFitResult.h"

#include <fstream>
#include "legendre3f3.h"
//#include "chebychev3f3.h"
#include "functions.h"
#include "math.h"
#include "functions-roofit.h"

#ifndef CALCACCEPTANCECLASSH
#define CALCACCEPTANCECLASSH

class calcacceptanceclass 
{
public:
  calcacceptanceclass(const TString&);
  ~calcacceptanceclass();

  double evaluate(double, double, double);
  void setcoeff(bool rand=false);

 private:
  bool debug;
  TString filename;
  
  //  RooFitResult* fr;
  
  static constexpr double epsilon=0.0001;
  static constexpr int idiv=19;
  static constexpr int jdiv=19;
  static constexpr int kdiv=19;

  double dcoeff3[ordermax3i][ordermax3j][ordermax3k];

  void findminimum(double&, double&, double&, double&, double, double, double, double, double, double);
  void findmaximum(double&, double&, double&, double&, double, double, double, double, double, double);
  //  double findmaximum(double, double, double, double, double, double, double);

  double integral(double, double, double, double, double, double);

  double func0_min,func0_max;

};


 #endif
