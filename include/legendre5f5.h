#include "legendre.h"

#ifndef LEGENDRE5F5H
#define LEGENDRE5F5H

const int ordermax5i=7;
const int ordermax5j=7;
const int ordermax5k=7;
const int ordermax5l=9;
const int ordermax5m=9;

inline double pdfval5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m])
{
   double value=0.;
   for (int i=0; i<ordermax5i;i++) {
     for (int j=0; j<ordermax5j;j++) {
       for (int k=0; k<ordermax5k;k++) {
	 for (int l=0; l<ordermax5l;l++) {
	   for (int m=0; m<ordermax5m;m++) {
	     if (coeff[i][j][k][l][m]!=0.) value+=
	       coeff[i][j][k][l][m]*
	       leg(i,costh0)*
	       leg(j,costh1)*
	       leg(k,costh2)*
	       leg(l,phi1)*
	       leg(m,phi2);
	   }
	 }
       }  
     }
   }
   return value;
}


inline double primitiveval5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m])
{
  double value=0.;
  for (int i=0; i<ordermax5i;i++) {
    for (int j=0; j<ordermax5j;j++) {
      for (int k=0; k<ordermax5k;k++) {
	for (int l=0; l<ordermax5l;l++) {
	  for (int m=0; m<ordermax5m;m++) {
	    if (coeff[i][j][k][l][m]!=0.) value+=
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
  return value;
}



inline double primitivedcostheta0val5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m])
{
  double value=0.;
  for (int i=0; i<ordermax5i;i++) {
    for (int j=0; j<ordermax5j;j++) {
      for (int k=0; k<ordermax5k;k++) {
	for (int l=0; l<ordermax5l;l++) {
	  for (int m=0; m<ordermax5m;m++) {
	    if (coeff[i][j][k][l][m]!=0.) value+=
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
  return value;
}

inline double primitivedcostheta1val5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m])
{
  double value=0.;
  for (int i=0; i<ordermax5i;i++) {
    for (int j=0; j<ordermax5j;j++) {
      for (int k=0; k<ordermax5k;k++) {
	for (int l=0; l<ordermax5l;l++) {
	  for (int m=0; m<ordermax5m;m++) {
	    if (coeff[i][j][k][l][m]!=0.) value+=
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
  return value;
}

inline double primitivedcostheta2val5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m])
{
  double value=0.;
  for (int i=0; i<ordermax5i;i++) {
    for (int j=0; j<ordermax5j;j++) {
      for (int k=0; k<ordermax5k;k++) {
	for (int l=0; l<ordermax5l;l++) {
	  for (int m=0; m<ordermax5m;m++) {
	    if (coeff[i][j][k][l][m]!=0.) value+=
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
  return value;
}

inline double primitivedphi1val5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m])
{
  double value=0.;
  for (int i=0; i<ordermax5i;i++) {
    for (int j=0; j<ordermax5j;j++) {
      for (int k=0; k<ordermax5k;k++) {
	for (int l=0; l<ordermax5l;l++) {
	  for (int m=0; m<ordermax5m;m++) {
	    if (coeff[i][j][k][l][m]!=0.) value+=
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
  return value;
}



inline double primitivedphi2val5f5(double costh0, double costh1, double costh2, double phi1, double phi2, double coeff[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m])
{
  double value=0.;
  for (int i=0; i<ordermax5i;i++) {
    for (int j=0; j<ordermax5j;j++) {
      for (int k=0; k<ordermax5k;k++) {
	for (int l=0; l<ordermax5l;l++) {
	  for (int m=0; m<ordermax5m;m++) {
	    if (coeff[i][j][k][l][m]!=0.) value+=
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
  return value;
}

#endif
