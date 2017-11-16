#include "TF3.h"
#include "math.h"
#include "TClass.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"

double findminimum(const TF3&, double, double, double, double, double, double, double);
double findminimum2(TF3&);
double findmaximum(const TF3&, double, double, double, double, double, double, double);
double findmaximum2(TF3&);
