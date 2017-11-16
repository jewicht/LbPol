#include "createRooLegendre5.h"


void writecoeff5(const TString& filename, RooRealVar* coeff5[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m])
{
 std::fstream filestr;
 filestr.open(filename, std::fstream::out);

 for (int i=0; i<ordermax5i;i++) {
   for (int j=0; j<ordermax5j;j++) {
     for (int k=0; k<ordermax5k;k++) {
       for (int l=0; l<ordermax5l;l++) {
	 for (int m=0; m<ordermax5m;m++) {
	   
	   if (coeff5[i][j][k][l][m]) {
	     filestr << i << " " << j << " " << k << " " << l << " " << m << " "
		     << coeff5[i][j][k][l][m]->getVal() << " " 
		     << coeff5[i][j][k][l][m]->getError() << " "
		     << coeff5[i][j][k][l][m]->isConstant() << std::endl;
	   } else {
	     filestr << i << " " 
		     << j << " " 
		     << k << " " 
		     << l << " " 
		     << m << " 0.0 0.0 0" << std::endl;
	   }
	 }
       }
     }
   }
 }
 filestr.close();
}

void readcoeff5_rrv(const TString& filename, RooRealVar* coeff5[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m])
{

  const bool debug=false;
  
  std::ifstream filestr;
  filestr.open(filename);
  if (!filestr.good()) {
    std::cerr << "Can't open " << filename <<  " . Bye bye..." << std::endl;
    return;
  }

  double value, error;
  bool isconstant;
  int ii,jj,kk,ll,mm;
  while (!filestr.eof()) {
    filestr >>  ii >> jj >> kk >> ll >> mm >> value >> error >> isconstant;// >> std::endl;
    if (debug) std::cout << ii << " " << jj << " " << kk << " " << ll << " "<< mm << " "<< value <<  " " << error << " " << isconstant << std::endl;
    if (coeff5[ii][jj][kk][ll][mm]) {
      if (debug) std::cout << "coeff5_" << ii << jj << kk << ll << mm << " is set" << std::endl; 
      coeff5[ii][jj][kk][ll][mm]->setVal(value);
      coeff5[ii][jj][kk][ll][mm]->setError(error);
      //	  coeff3[i][j][k][ll][mm]->setConstant(isconstant);
      //      coeff5[ii][jj][kk][ll][mm]->setConstant();
    }
  }

  //   for (int i=0; i<ordermax5i;i++) {
  //     for (int j=0; j<ordermax5j;j++) {
  //       for (int k=0; k<ordermax5k;k++) {
  // 	for (int l=0; l<ordermax5l;l++) {
  // 	  for (int m=0; m<ordermax5m;m++) {
  
  
  // 	  }
  // 	}
  //       }
  //     }
  //   }
  filestr.close();
}

void constornotconstcoeff5(RooRealVar* coeff5[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m])
{
  int nConst=0;
  for (int i=0; i<ordermax5i;i++) {
    for (int j=0; j<ordermax5j;j++) {
      for (int k=0; k<ordermax5k;k++) {
	for (int l=0; l<ordermax5l;l++) {
	  for (int m=0; m<ordermax5m;m++) {
	    if (coeff5[i][j][k][l][m]->isConstant()) nConst++;
	  }
	}
      }
    }
  }
  std::cout << nConst << " constant parameters out of " << ordermax5i*ordermax5j*ordermax5k*ordermax5l*ordermax5m << std::endl;
}

void setconstcoeff5(RooRealVar* coeff5[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m])
{
  for (int i=0; i<ordermax5i;i++) {
    for (int j=0; j<ordermax5j;j++) {
      for (int k=0; k<ordermax5k;k++) {
	for (int l=0; l<ordermax5l;l++) {
	  for (int m=0; m<ordermax5m;m++) {
	    if (((i%2==0) or (i==1)) and 
		((j%2==0) or (j==1)) and 
		((k%2==0) or (k==1)) and 
		((l%2==0) or (l==1)) and 
		((m%2==0) or (m==1))) {
	    } else {
	      coeff5[i][j][k][l][m]->setVal(0.);
	      coeff5[i][j][k][l][m]->setConstant();
	    }
	  }
	}
      }
    }
  }
}

void setconstcoeff5_maxtotorder(RooRealVar* coeff5[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m], int maxtotorder)
{
  for (int i=0; i<ordermax5i;i++) {
    for (int j=0; j<ordermax5j;j++) {
      for (int k=0; k<ordermax5k;k++) {
	for (int l=0; l<ordermax5l;l++) {
	  for (int m=0; m<ordermax5m;m++) {
	    if (i+j+k+l+m>=maxtotorder) {
	      coeff5[i][j][k][l][m]->setVal(0.);
	      coeff5[i][j][k][l][m]->setConstant();
	    }
	  }
	}
      }
    }
  }
}

void createRooLegendre5v2(RooRealVar& costheta, RooRealVar& costheta1, RooRealVar& costheta2, RooRealVar& phi1, RooRealVar& phi2, TString suffix, RooLegendre5Pdfv2*& pdf, RooRealVar* coeff5[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m], int maxorderi, int maxorderj, int maxorderk, int maxorderl, int maxorderm, int maxtotorder, int maxcorrel, const TString& filename)
{

  RooArgList varlist;
  const double dmin=-1.;
  const double dmax=1.;

  std::cout << "Adding "; 
  for (int i=0; i<ordermax5i;i++) {
    for (int j=0; j<ordermax5j;j++) {
      for (int k=0; k<ordermax5k;k++) {
	for (int l=0; l<ordermax5l;l++) {
	  for (int m=0; m<ordermax5m;m++) {
	    int correl=0;
	    if (i>0) correl++;
	    if (j>0) correl++;
	    if (k>0) correl++;
	    if (l>0) correl++;
	    if (m>0) correl++;

	    if (((i%2==0) or (i==1)) and 
		((j%2==0) or (j==1)) and 
		((k%2==0) or (k==1)) and 
		((l%2==0) or (l==1)) and 
		((m%2==0) or (m==1)) and
		(i<maxorderi) and
		(j<maxorderj) and
		(k<maxorderk) and
		(l<maxorderl) and
		(m<maxorderm) and
		(correl<=maxcorrel) and
		(i+j+k+l+m<=maxtotorder)) {
	      std::cout << "coeff" << i << j << k << l << m << "; "; 
	      coeff5[i][j][k][l][m] = new RooRealVar(suffix + "-c" + istr(i) + istr(j) + istr(k) + istr(l) + istr(m), suffix + "-c" + istr(i) + istr(j) + istr(k) + istr(l) + istr(m),0.,dmin,dmax);
	      varlist.add(*coeff5[i][j][k][l][m]);
	    } else {
	      coeff5[i][j][k][l][m] = NULL;
	    }
	  }
	}
      }
    }
  }
  std::cout << std::endl;
  coeff5[0][0][0][0][0]->setVal(1.);
  coeff5[0][0][0][0][0]->setConstant();
  if (filename!="") readcoeff5_rrv(filename, coeff5);

  std::cout << "Number of parameters = " << varlist.getSize() << std::endl;
  pdf = new RooLegendre5Pdfv2(suffix,suffix,costheta,costheta1,costheta2,phi1,phi2,varlist);//,ordermax5i,ordermax5j,ordermax5k,ordermax5l,ordermax5m);
}


