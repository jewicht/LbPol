#include "legendre5f5.h"
#include "legendre.h"
#include <iostream>

double pdfval5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[4][4][4][4][7])
{
   double value=0.;
   for (int i=0; i<ordermaxi;i++) {
     for (int j=0; j<ordermaxj;j++) {
       for (int k=0; k<ordermaxk;k++) {
	 for (int l=0; l<ordermaxl;l++) {
	   for (int m=0; m<ordermaxm;m++) {
	     if (((i%2==0) or (i==1)) and 
		 ((j%2==0) or (j==1)) and 
		 ((k%2==0) or (k==1)) and
		 ((l%2==0) or (l==1)) and
		 ((m%2==0) or (m==1))) {
	       //	       std::cout << " coeff = " << coeff[i][j][k][l][m];// << std::endl;
	       value+=coeff[i][j][k][l][m]*leg(i,costh0)*leg(j,costh1)*leg(k,costh2)*leg(l,phi1)*leg(m,phi2);
	     }
	   }
	 }
       }  
     }
   }
   //   std::cout << "value = " << value << std::endl;
   return value;
}


double primitiveval5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[4][4][4][4][7])
{
  double value=0.;
   for (int i=0; i<ordermaxi;i++) {
     for (int j=0; j<ordermaxj;j++) {
       for (int k=0; k<ordermaxk;k++) {
	 for (int l=0; l<ordermaxl;l++) {
	   for (int m=0; m<ordermaxm;m++) {
	     if (((i%2==0) or (i==1)) and 
		 ((j%2==0) or (j==1)) and 
		 ((k%2==0) or (k==1)) and
		 ((l%2==0) or (l==1)) and
		 ((m%2==0) or (m==1))) {
	       value+=
		 coeff[i][j][k][l][m]
		 *legprim(i,costh0)
		 *legprim(j,costh1)
		 *legprim(k,costh2)
		 *legprim(l,phi1)
		 *legprim(m,phi2);
	     }
	   }
	 }
       }  
     }
   }
   return value;
}



double primitivedcostheta0val5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[4][4][4][4][7])
{
  double value=0.;
   for (int i=0; i<ordermaxi;i++) {
     for (int j=0; j<ordermaxj;j++) {
       for (int k=0; k<ordermaxk;k++) {
	 for (int l=0; l<ordermaxl;l++) {
	   for (int m=0; m<ordermaxm;m++) {
	     if (((i%2==0) or (i==1)) and 
		 ((j%2==0) or (j==1)) and 
		 ((k%2==0) or (k==1)) and
		 ((l%2==0) or (l==1)) and
		 ((m%2==0) or (m==1))) {
	       value+=
		 coeff[i][j][k][l][m]
		 *leg(i,costh0)
		 *legprim(j,costh1)
		 *legprim(k,costh2)
		 *legprim(l,phi1)
		 *legprim(m,phi2);
	     }
	   }
	 }
       }  
     }
   }
   return value;
}

double primitivedcostheta1val5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[4][4][4][4][7])
{
  double value=0.;
   for (int i=0; i<ordermaxi;i++) {
     for (int j=0; j<ordermaxj;j++) {
       for (int k=0; k<ordermaxk;k++) {
	 for (int l=0; l<ordermaxl;l++) {
	   for (int m=0; m<ordermaxm;m++) {
	     if (((i%2==0) or (i==1)) and 
		 ((j%2==0) or (j==1)) and 
		 ((k%2==0) or (k==1)) and
		 ((l%2==0) or (l==1)) and
		 ((m%2==0) or (m==1))) {
	       value+=
		 coeff[i][j][k][l][m]
		 *legprim(i,costh0)
		 *leg(j,costh1)
		 *legprim(k,costh2)
		 *legprim(l,phi1)
		 *legprim(m,phi2);
	     }
	   }
	 }
       }  
     }
   }
   return value;
}

double primitivedcostheta2val5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[4][4][4][4][7])
{
  double value=0.;
   for (int i=0; i<ordermaxi;i++) {
     for (int j=0; j<ordermaxj;j++) {
       for (int k=0; k<ordermaxk;k++) {
	 for (int l=0; l<ordermaxl;l++) {
	   for (int m=0; m<ordermaxm;m++) {
	     if (((i%2==0) or (i==1)) and 
		 ((j%2==0) or (j==1)) and 
		 ((k%2==0) or (k==1)) and
		 ((l%2==0) or (l==1)) and
		 ((m%2==0) or (m==1))) {
	       value+=
		 coeff[i][j][k][l][m]
		 *legprim(i,costh0)
		 *legprim(j,costh1)
		 *leg(k,costh2)
		 *legprim(l,phi1)
		 *legprim(m,phi2);
	     }
	   }
	 }
       }  
     }
   }
   return value;
}

double primitivedphi1val5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[4][4][4][4][7])
{
  double value=0.;
   for (int i=0; i<ordermaxi;i++) {
     for (int j=0; j<ordermaxj;j++) {
       for (int k=0; k<ordermaxk;k++) {
	 for (int l=0; l<ordermaxl;l++) {
	   for (int m=0; m<ordermaxm;m++) {
	     if (((i%2==0) or (i==1)) and 
		 ((j%2==0) or (j==1)) and 
		 ((k%2==0) or (k==1)) and
		 ((l%2==0) or (l==1)) and
		 ((m%2==0) or (m==1))) {
	       value+=
		 coeff[i][j][k][l][m]
		 *legprim(i,costh0)
		 *legprim(j,costh1)
		 *legprim(k,costh2)
		 *leg(l,phi1)
		 *legprim(m,phi2);
	     }
	   }
	 }
       }  
     }
   }
   return value;
}



double primitivedphi2val5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[4][4][4][4][7])
{
  double value=0.;
   for (int i=0; i<ordermaxi;i++) {
     for (int j=0; j<ordermaxj;j++) {
       for (int k=0; k<ordermaxk;k++) {
	 for (int l=0; l<ordermaxl;l++) {
	   for (int m=0; m<ordermaxm;m++) {
	     if (((i%2==0) or (i==1)) and 
		 ((j%2==0) or (j==1)) and 
		 ((k%2==0) or (k==1)) and
		 ((l%2==0) or (l==1)) and
		 ((m%2==0) or (m==1))) {
	       value+=
		 coeff[i][j][k][l][m]
		 *legprim(i,costh0)
		 *legprim(j,costh1)
		 *legprim(k,costh2)
		 *legprim(l,phi1)
		 *leg(m,phi2);
	     }
	   }
	 }
       }  
     }
   }
   return value;
}


