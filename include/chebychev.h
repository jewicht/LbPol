#include <iostream>

#ifndef CHEBYCHEVH
#define CHEBYCHEVH


inline double cheb(int n, double x)
{
  const double x2=x*x;
  const double x3=x*x2;
  const double x4=x*x3;
  const double x5=x*x4;
  const double x6=x*x5;
  const double x7=x*x6;
  const double x8=x*x7;
  const double x9=x*x8;
  //  const double x10=x*x9;

  switch(n) {

  case  0  : return       1. ; break;
  case  1  : return       x  ; break;
  case  2  : return    2.*x2 -   1. ; break;
  case  3  : return    4.*x3 -   3.*x; break;
  case  4  : return    8.*x4 -   8.*x2 +   1.; break;
  case  5  : return   16.*x5 -  20.*x3 +   5.*x; break;
  case  6  : return   32.*x6 -  48.*x4 +  18.*x2 -   1.; break;
  case  7  : return   64.*x7 - 112.*x5 +  56.*x3 -   7.*x;  break;
  case  8  : return  128.*x8 - 256.*x6 + 160.*x4 -  32.*x2 + 1.; break;
  case  9  : return  256.*x9 - 576.*x7 + 432.*x5 - 120.*x3 + 9.*x; break;
  default: std::cerr << "error cheb n>9" << std::endl; return -1.; break;
  }
}

inline double chebprim(int n, double x)
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
  //  const double x11=x*x10;

  switch(n) {

  case  0  : return       x; break;
  case  1  : return       x2/2.  ; break;
  case  2  : return    2.*x3/3.   -      x; break;
  case  3  : return    4.*x4/4.   -   3.*x2/2.; break;
  case  4  : return    8.*x5/5.   -   8.*x3/3. +   x; break;
  case  5  : return   16.*x6/6.   -  20.*x4/4. +   5.*x2/2.; break;
  case  6  : return   32.*x7/7.   -  48.*x5/5. +  18.*x3/3. -   x; break;
  case  7  : return   64.*x8/8.   - 112.*x6/6. +  56.*x4/4. -   7.*x2/2.;  break;
  case  8  : return  128.*x9/9.   - 256.*x7/7. + 160.*x5/5. -  32.*x3/3. + x; break;
  case  9  : return  256.*x10/10. - 576.*x8/8. + 432.*x6/6. - 120.*x4/4. + 9.*x2/2; break;

  default: std::cerr << "error cheb n>9" << std::endl; return -1.; break;
  }
}


#endif
