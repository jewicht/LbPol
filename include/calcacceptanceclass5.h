#include <fstream>
#include "TRandom3.h"

#include "legendre5f5.h"
#include "functions.h"

#include "math.h"

#ifndef CALCACCEPTANCECLASS5H
#define CALCACCEPTANCECLASS5H

class calcacceptanceclass5 
{
public:
  calcacceptanceclass5(const TString&);
  ~calcacceptanceclass5(){};

  double evaluate(double, double, double,double,double);

 private:

  double dcoeff5[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m];
  double findminimum(double, double, double, double, double, double, double, double, double, double, double);  
  double findmaximum(double, double, double, double, double, double, double, double, double, double, double);
  double func0_min,func0_max;

  static constexpr double epsilon=0.0001;
  static const int idiv=9;
  static const int jdiv=9;
  static const int kdiv=9;
  static const int ldiv=9;
  static const int mdiv=19;
  
};


 #endif
