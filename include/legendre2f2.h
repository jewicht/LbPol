#include "legendre.h"

#ifndef LEGENDRE2F2H
#define LEGENDRE2F2H

const int ordermax2i=10;
const int ordermax2j=10;


inline double pdfval(double costh0, double costh1, double coeff[ordermax2i][ordermax2j])
{
   double value=0.;
   for (int i=0; i<ordermax2i;i++) {
     for (int j=0; j<ordermax2j;j++) {
       if (coeff[i][j]!=0.) value+=coeff[i][j]*leg(i,costh0)*leg(j,costh1);
     }
   }
   return value;
}

inline double primitiveval(double costh0, double costh1, double coeff[ordermax2i][ordermax2j])
{
   double value=0.;
   for (int i=0; i<ordermax2i;i++) {
     for (int j=0; j<ordermax2j;j++) {
       if (coeff[i][j]!=0.) value+=coeff[i][j]*legprim(i,costh0)*legprim(j,costh1);
     }
   }
   return value;
}



inline double primitivedphi1val(double costh0, double costh1, double coeff[ordermax2i][ordermax2j])
{
   double value=0.;
   for (int i=0; i<ordermax2i;i++) {
     for (int j=0; j<ordermax2j;j++) {
       if (coeff[i][j]!=0.) value+=coeff[i][j]*leg(i,costh0)*legprim(j,costh1);
     }
   }
   return value;
}

inline double primitivedphi2val(double costh0, double costh1, double coeff[ordermax2i][ordermax2j])
{
   double value=0.;
   for (int i=0; i<ordermax2i;i++) {
     for (int j=0; j<ordermax2j;j++) {
       if (coeff[i][j]!=0.) value+=coeff[i][j]*legprim(i,costh0)*leg(j,costh1);
     }
   }
   return value;
}

/* inline double primitivedcostheta2val(double costh0, double costh1, double costh2, double coeff[ordermax3i][ordermax3j][ordermax3k]) */
/* { */
/*    double value=0.; */
/*    for (int i=0; i<ordermax3i;i++) { */
/*      for (int j=0; j<ordermax3j;j++) { */
/*        for (int k=0; k<ordermax3k;k++) { */
/* 	 if (coeff[i][j][k]!=0.) value+=coeff[i][j][k]*legprim(i,costh0)*legprim(j,costh1)*leg(k,costh2); */
/*        }   */
/*      } */
/*    } */
/*    return value; */
/* } */


#endif
