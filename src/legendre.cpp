#include "legendre.h"
#include <iostream>


double leg(int n, double x)
{
  const double x2=x*x;
  const double x3=x*x2;
  const double x4=x*x3;
  const double x5=x*x4;
  const double x6=x*x5;
  const double x7=x*x6;
  const double x8=x*x7;
  const double x9=x*x8;
  const double x10=x*x9;

  switch(n) {

  case  0  : return  1. ; break;
  case  1  : return  x  ; break;
  case  2  : return  1.50000000000000*x2 - 0.500000000000000  ; break;
  case  3  : return  2.50000000000000*x3 - 1.50000000000000*x  ; break;
  case  4  : return  4.37500000000000*x4 - 3.75000000000000*x2 + 0.375000000000000  ; break;
  case  5  : return  7.87500000000000*x5 - 8.75000000000000*x3 + 1.87500000000000*x  ; break;
  case  6  : return  14.4375000000000*x6 - 19.6875000000000*x4 + 6.56250000000000*x2 - 0.312500000000000  ; break;
  case  7  : return  26.8125000000000*x7 - 43.3125000000000*x5 + 19.6875000000000*x3 - 2.18750000000000*x  ; break;
  case  8  : return  50.2734375000000*x8 - 93.8437500000000*x6 + 54.1406250000000*x4 - 9.84375000000000*x2 + 0.273437500000000  ; break;
  case  9  : return  94.9609375000000*x9 - 201.093750000000*x7 + 140.765625000000*x5 - 36.0937500000000*x3 + 2.46093750000000*x  ; break;
  case  10  : return 180.425781250000*x10 - 427.324218750000*x8 + 351.914062500000*x6 - 117.304687500000*x4 + 13.5351562500000*x2 - 0.246093750000000  ; break;
    
  default: std::cerr << "error legendre n>10" << std::endl; return -1.; break;
  }
}

double legprim(int n, double x)
{
  const double x2=x*x;
  const double x3=x*x2;
  const double x4=x*x3;
  const double x5=x*x4;
  const double x6=x*x5;
  const double x7=x*x6;
  const double x8=x*x7;
  const double x9=x*x8;
  const double x10=x*x9;
  const double x11=x*x10;

  switch(n) {

  case  0  : return  x ; break;
  case  1  : return  0.5*x2  ; break;
  case  2  : return  0.5*x3 - 0.5*x  ; break;
  case  3  : return  0.625*x4 - 0.75*x2  ; break;
  case  4  : return  0.875*x5 - 1.25*x3 + 0.375*x  ; break;
  case  5  : return  1.3125*x6 - 2.1875*x4 + 0.9375*x2  ; break;
  case  6  : return  2.0625*x7 - 3.9375*x5 + 2.1875*x3 - 0.3125*x  ; break;
  case  7  : return  3.3515625*x8 - 7.21875*x6 + 4.921875*x4 - 1.09375*x2  ; break;
  case  8  : return  5.5859375*x9 - 13.40625*x7 + 10.828125*x5 - 3.28125*x3 + 0.2734375*x  ; break;
  case  9  : return  9.49609375*x10 - 25.13671875*x8 + 23.4609375*x6 - 9.0234375*x4 + 1.23046875*x2  ; break;
  case  10  : return  16.40234375*x11 - 47.48046875*x9 + 50.2734375*x7 - 23.4609375*x5 + 4.51171875*x3 - 0.24609375*x  ; break;

  default: std::cerr << "error legendre n>10" << std::endl; return -1.; break;
  }
}
