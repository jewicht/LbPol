#include "chebychev.h"
#include "legendre3f3.h"

#ifndef CHEBYCHEV3F3H
#define CHEBYCHEV3F3H

/* const int ordermax3i=7; */
/* const int ordermax3j=7; */
/* const int ordermax3k=7; */

inline double chebychev3f3(double costh0, double costh1, double costh2, double coeff[ordermax3i][ordermax3j][ordermax3k])
{
   double value=0.;
   for (int i=0; i<ordermax3i;i++) {
     for (int j=0; j<ordermax3j;j++) {
       for (int k=0; k<ordermax3k;k++) {
	 if (coeff[i][j][k]!=0.) value+=coeff[i][j][k]*cheb(i,costh0)*cheb(j,costh1)*cheb(k,costh2);
       }  
     }
   }
   return value;
}

inline double chebychev3f3primitive(double costh0, double costh1, double costh2, double coeff[ordermax3i][ordermax3j][ordermax3k])
{
   double value=0.;
   for (int i=0; i<ordermax3i;i++) {
     for (int j=0; j<ordermax3j;j++) {
       for (int k=0; k<ordermax3k;k++) {
	 if (coeff[i][j][k]!=0.) value+=coeff[i][j][k]*chebprim(i,costh0)*chebprim(j,costh1)*chebprim(k,costh2);
       }  
     }
   }
   return value;
}



inline double chebychev3f3primitivedcostheta(double costh0, double costh1, double costh2, double coeff[ordermax3i][ordermax3j][ordermax3k])
{
   double value=0.;
   for (int i=0; i<ordermax3i;i++) {
     for (int j=0; j<ordermax3j;j++) {
       for (int k=0; k<ordermax3k;k++) {
	 if (coeff[i][j][k]!=0.) value+=coeff[i][j][k]*cheb(i,costh0)*chebprim(j,costh1)*chebprim(k,costh2);
       }  
     }
   }
   return value;
}

inline double chebychev3f3primitivedcostheta1(double costh0, double costh1, double costh2, double coeff[ordermax3i][ordermax3j][ordermax3k])
{
   double value=0.;
   for (int i=0; i<ordermax3i;i++) {
     for (int j=0; j<ordermax3j;j++) {
       for (int k=0; k<ordermax3k;k++) {
	 if (coeff[i][j][k]!=0.) value+=coeff[i][j][k]*chebprim(i,costh0)*cheb(j,costh1)*chebprim(k,costh2);
       }  
     }
   }
   return value;
}

inline double chebychev3f3primitivedcostheta2(double costh0, double costh1, double costh2, double coeff[ordermax3i][ordermax3j][ordermax3k])
{
   double value=0.;
   for (int i=0; i<ordermax3i;i++) {
     for (int j=0; j<ordermax3j;j++) {
       for (int k=0; k<ordermax3k;k++) {
	 if (coeff[i][j][k]!=0.) value+=coeff[i][j][k]*chebprim(i,costh0)*chebprim(j,costh1)*cheb(k,costh2);
       }  
     }
   }
   return value;
}


#endif
