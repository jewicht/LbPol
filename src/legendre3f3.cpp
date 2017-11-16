#include "legendre3f3.h"
#include "legendre.h"
#include <iostream>

double pdfval(double costh0, double costh1, double costh2, double coeff[4][4][7])
{
   double value=0.;
   for (int i=0; i<ordermaxi;i++) {
     for (int j=0; j<ordermaxj;j++) {
       for (int k=0; k<ordermaxk;k++) {
	 if (((i%2==0) or (i==1)) and 
	     ((j%2==0) or (j==1)) and 
	     ((k%2==0) or (k==1))) {
	   value+=coeff[i][j][k]*leg(i,costh0)*leg(j,costh1)*leg(k,costh2);
	 }
       }  
     }
   }
   return value;
}


double primitiveval(double costh0, double costh1, double costh2, double coeff[4][4][7])
{
   double value=0.;
   for (int i=0; i<ordermaxi;i++) {
     for (int j=0; j<ordermaxj;j++) {
       for (int k=0; k<ordermaxk;k++) {
	 if (((i%2==0) or (i==1)) and 
	     ((j%2==0) or (j==1)) and 
	     ((k%2==0) or (k==1))) {
	   value+=coeff[i][j][k]*legprim(i,costh0)*legprim(j,costh1)*legprim(k,costh2);
	 }
       }  
     }
   }
   return value;
}


