#include "legendre.h"

#ifndef LEGENDRE3F3H
#define LEGENDRE3F3H

const int ordermax3i=7;
const int ordermax3j=5;
const int ordermax3k=15;

inline double pdfval(double costh0, double costh1, double costh2, double coeff[ordermax3i][ordermax3j][ordermax3k])
{
   double value=0.;
   for (int i=0; i<ordermax3i;i++) {
     for (int j=0; j<ordermax3j;j++) {
       for (int k=0; k<ordermax3k;k++) {
	 if (coeff[i][j][k]!=0.) value+=coeff[i][j][k]*leg(i,costh0)*leg(j,costh1)*leg(k,costh2);
       }  
     }
   }
   return value;
}

inline double primitiveval(double costh0, double costh1, double costh2, double coeff[ordermax3i][ordermax3j][ordermax3k])
{
   double value=0.;
   for (int i=0; i<ordermax3i;i++) {
     for (int j=0; j<ordermax3j;j++) {
       for (int k=0; k<ordermax3k;k++) {
	 if (coeff[i][j][k]!=0.) value+=coeff[i][j][k]*legprim(i,costh0)*legprim(j,costh1)*legprim(k,costh2);
       }  
     }
   }
   return value;
}



inline double primitivedcosthetaval(double costh0, double costh1, double costh2, double coeff[ordermax3i][ordermax3j][ordermax3k])
{
   double value=0.;
   for (int i=0; i<ordermax3i;i++) {
     for (int j=0; j<ordermax3j;j++) {
       for (int k=0; k<ordermax3k;k++) {
	 if (coeff[i][j][k]!=0.) value+=coeff[i][j][k]*leg(i,costh0)*legprim(j,costh1)*legprim(k,costh2);
       }  
     }
   }
   return value;
}

inline double primitivedcostheta1val(double costh0, double costh1, double costh2, double coeff[ordermax3i][ordermax3j][ordermax3k])
{
   double value=0.;
   for (int i=0; i<ordermax3i;i++) {
     for (int j=0; j<ordermax3j;j++) {
       for (int k=0; k<ordermax3k;k++) {
	 if (coeff[i][j][k]!=0.) value+=coeff[i][j][k]*legprim(i,costh0)*leg(j,costh1)*legprim(k,costh2);
       }  
     }
   }
   return value;
}

inline double primitivedcostheta2val(double costh0, double costh1, double costh2, double coeff[ordermax3i][ordermax3j][ordermax3k])
{
   double value=0.;
   for (int i=0; i<ordermax3i;i++) {
     for (int j=0; j<ordermax3j;j++) {
       for (int k=0; k<ordermax3k;k++) {
	 if (coeff[i][j][k]!=0.) value+=coeff[i][j][k]*legprim(i,costh0)*legprim(j,costh1)*leg(k,costh2);
       }  
     }
   }
   return value;
}


#endif
