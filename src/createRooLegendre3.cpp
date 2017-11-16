#include "createRooLegendre3.h"

// void writecoeff3(const TString& filename, RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k])
// {
//  std::fstream filestr;
//  filestr.open(filename, std::fstream::out);

//  for (int i=0; i<ordermax3i;i++) {
//    for (int j=0; j<ordermax3j;j++) {
//      for (int k=0; k<ordermax3k;k++) {
       
//        if (coeff3[i][j][k]) {
// 	 filestr << i << " " 
// 		 << j << " " 
// 		 << k << " "
// 		 << coeff3[i][j][k]->getVal() << " " 
// 		 << coeff3[i][j][k]->getError() << " "
// 		 << coeff3[i][j][k]->isConstant() << std::endl;
//        } else {
// 	 filestr << i << " " 
// 		 << j << " " 
// 		 << k << " 0.0 0.0 0" << std::endl;
//        }
//      }
//    }
//  }
//  filestr.close();
// }


void coeff3_rrv2double(RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k], double dcoeff3[ordermax3i][ordermax3j][ordermax3k])
{
   for (int i=0; i<ordermax3i;i++) {
    for (int j=0; j<ordermax3j;j++) {
      for (int k=0; k<ordermax3k;k++) {
	if (coeff3[i][j][k]) {
	  dcoeff3[i][j][k]=coeff3[i][j][k]->getVal();
	} else {
	  dcoeff3[i][j][k]=0.;
	}
      }
    }
   } 
   
}

void readcoeff3_rrv(const TString& filename, RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k])
{
  RooFitResult* fr = getRooFitResult(filename);

  const RooArgList& floatPars=fr->floatParsFinal();

  for (int i=0; i<ordermax3i;i++) {
    for (int j=0; j<ordermax3j;j++) {
      for (int k=0; k<ordermax3k;k++) {
	if (i+j+k>0 and coeff3[i][j][k]) {
	  RooRealVar* var = (RooRealVar*)floatPars.find(coeff3[i][j][k]->GetName());
	  if (var) {
	    coeff3[i][j][k]->setVal(var->getVal());
	    coeff3[i][j][k]->setError(var->getError());
	    std::cerr << "readcoeff3_rrv: could find " <<  coeff3[i][j][k]->GetName() << std::endl;
	  } else {
	    std::cerr << "readcoeff3_rrv: can't find " <<  coeff3[i][j][k]->GetName() << std::endl;
	  }
	}
      }
    }
  }
  coeff3[0][0][0]->setVal(1.);
  coeff3[0][0][0]->setConstant();
  delete fr;
}

// void readcoeff3_rrv(const TString& filename, RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k])
// {
//   const bool debug=false;

//   std::ifstream filestr;
//   filestr.open("/afs/cern.ch/user/j/jwicht/fitLambdab2JpsiL0/" + filename);
//   if (!filestr.good()) {
//     std::cerr << "Can't open " << filename <<  " . Bye bye..." << std::endl;
//     exit(1);
//   }

//   double value, error;
//   bool isconstant;
//   int ii,jj,kk;
//   while (!filestr.eof()) {
//     filestr >>  ii >> jj >> kk >> value >> error >> isconstant;// >> std::endl;
//     if (debug) std::cout << ii << " " << jj << " " << kk << " " << value <<  " " << error << " " << isconstant << std::endl;
//     if (coeff3[ii][jj][kk]) {
//       if (debug) std::cout << "coeff3_" << ii << jj << kk << " is set" << std::endl; 
//       coeff3[ii][jj][kk]->setVal(value);
//       coeff3[ii][jj][kk]->setError(error);
//       //	  coeff3[i][j][k]->setConstant(isconstant);
//       //      coeff3[ii][jj][kk]->setConstant();
//     } else {
//       if (value!=0.) std::cerr << "readcoeff3_rrv: RooRealVar " << ii << jj << kk << " does not exist" << std::endl;
//     }
//   }
//   if (debug) mygetchar;

//   //   for (int i=0; i<ordermax3i;i++) {
//   //     for (int j=0; j<ordermax3j;j++) {
//   //       for (int k=0; k<ordermax3k;k++) {
  
//   //       }
//   //     }
//   //   }
//   filestr.close();
// }


// void constornotconstcoeff3(RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k])
// {
//   int nConst=0;
//   for (int i=0; i<ordermax3i;i++) {
//     for (int j=0; j<ordermax3j;j++) {
//       for (int k=0; k<ordermax3k;k++) {
// 	if (coeff3[i][j][k]->isConstant()) nConst++;
//       }
//     }
//   }
//   std::cout << nConst << " constant parameters out of " << ordermax3i*ordermax3j*ordermax3k << std::endl;
// }

// void setconstcoeff3(RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k])
// {
//   for (int i=0; i<ordermax3i;i++) {
//     for (int j=0; j<ordermax3j;j++) {
//       for (int k=0; k<ordermax3k;k++) {
// 	if (((i%2==0) or (i==1)) and 
// 	    ((j%2==0) or (j==1)) and 
// 	    ((k%2==0) or (k==1))) {
// 	} else {
// 	  coeff3[i][j][k]->setVal(0.);
// 	  coeff3[i][j][k]->setConstant();
// 	}
//       }
//     }
//   }
//   coeff3[0][0][0]->setVal(1.);
//   coeff3[0][0][0]->setConstant();
// }



// void setconstcoeff3_maxtotorder(RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k], int maxtotorder)
// {
//   for (int i=0; i<ordermax3i;i++) {
//     for (int j=0; j<ordermax3j;j++) {
//       for (int k=0; k<ordermax3k;k++) {
// 	if (i+j+k>=maxtotorder) {
// 	  coeff3[i][j][k]->setVal(0.);
// 	  coeff3[i][j][k]->setConstant();
// 	}
//       }
//     }
//   }
// }



// void setconstcoeff3_maxorder(RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k], int maxi, int maxj, int maxk)
// {
//   for (int i=0; i<ordermax3i;i++) {
//     for (int j=0; j<ordermax3j;j++) {
//       for (int k=0; k<ordermax3k;k++) {
// 	if (i>=maxi or j>=maxj or k>=maxk) {
// 	  coeff3[i][j][k]->setVal(0.);
// 	  coeff3[i][j][k]->setConstant();
// 	}
//       }
//     }
//   }
// }


// void disableconstcoeff3_maxorder(RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k], int maxi, int maxj, int maxk)
// {
//   for (int i=0; i<ordermax3i;i++) {
//     for (int j=0; j<ordermax3j;j++) {
//       for (int k=0; k<ordermax3k;k++) {
// 	if (i>=maxi or j>=maxj or k>=maxk) {
// 	  coeff3[i][j][k]->setVal(0.);
// 	  coeff3[i][j][k]->setConstant(false);
// 	}
//       }
//     }
//   }
//   setconstcoeff3(coeff3);
// }

// void setmaxordercoeff3(RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k], int maxi, int maxj, int maxk)
// {
//   for (int i=0; i<ordermax3i;i++) {
//     for (int j=0; j<ordermax3j;j++) {
//       for (int k=0; k<ordermax3k;k++) {
// 	if (i>=maxi or j>=maxj or k>=maxk) {
// 	  coeff3[i][j][k]->setVal(0.);
// 	  coeff3[i][j][k]->setConstant();
// 	}
//       }
//     }
//   }
// }




void createRRV3(bool doBd2JpsiKS, RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k],TString prefix, int maxorderi, int maxorderj, int maxorderk, int maxtotorder, int maxcorrel, const TString& filename)
{
  const double dmin=-1.;
  const double dmax=1.;
  for (int i=0; i<ordermax3i;i++) {
    for (int j=0; j<ordermax3j;j++) {
      for (int k=0; k<ordermax3k;k++) {
	int correl=0;
	if (i>0) correl++;
	if (j>0) correl++;
	if (k>0) correl++;
	coeff3[i][j][k]=NULL;
	if (doBd2JpsiKS) {
	  if (
	      ((i%2==0) or (i==1)) and 
	      ((j%2==0) or (j==1)) and 
	      ((k%2==0) or (k==1)) and 
	      //	      (i%2==0) and 
	      //	      (j%2==0) and 
	      //	      (k%2==0) and 
	      (i<maxorderi) and
	      (j<maxorderj) and
	      (k<maxorderk) and	
	      (correl<=maxcorrel) and
	      (i+j<=8) and
	      (i+j+k<=maxtotorder)) {
	    coeff3[i][j][k] = new RooRealVar(prefix + "-c" + istr(i,"02") + istr(j,"02") + istr(k,"02"), prefix + "-c" + istr(i,"02") + istr(j,"02") + istr(k,"02"),0.,dmin,dmax);
	  }
	} else {
	  if (
	      ((i%2==0) or (i==1)) and 
	      //	    ((j%2==0) or (j==1)) and 
	      ((k%2==0) or (k==1)) and 
	      //	      ((i%2==0)) and 
	      //	    ((j%2==1) or (j==0)) and 
	      //	    ((k%2==0)) and 
	      (i<maxorderi) and
	      (j<maxorderj) and
	      (k<maxorderk) and	
	      (correl<=maxcorrel) and
	      (i+j+k<=maxtotorder)) {
	    coeff3[i][j][k] = new RooRealVar(prefix + "-c" + istr(i,"02") + istr(j,"02") + istr(k,"02"), prefix + "-c" + istr(i,"02") + istr(j,"02") + istr(k,"02"),0.,dmin,dmax);
	  }
	}
      }
    }
  }
  //  std::cout << coeff3[0][0][0] << std::endl;
  //  mygetchar;
  //  coeff3[0][0][0]->setMax(1000.);
  coeff3[0][0][0]->setVal(1.);
  coeff3[0][0][0]->setConstant();
  //  mygetchar;
  if (filename!="") readcoeff3_rrv(filename,coeff3);
}

void createRooLegendre3v2(RooRealVar& costheta, RooRealVar& costheta1, RooRealVar& costheta2, TString prefix, RooAbsPdf*& pdf, RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k])
{
  RooArgList varlist;
  for (int i=0; i<ordermax3i;i++) {
    for (int j=0; j<ordermax3j;j++) {
      for (int k=0; k<ordermax3k;k++) {
	if (coeff3[i][j][k]) varlist.add(*coeff3[i][j][k]);	  
      }
    }
  }
  pdf = new RooLegendre3Pdfv2(prefix,prefix,costheta,costheta1,costheta2,varlist);//, ordermax3i, ordermax3j, ordermax3k);
}


// void createRooChebychev3(RooRealVar& costheta, RooRealVar& costheta1, RooRealVar& costheta2, TString suffix, RooAbsPdf*& pdf, RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k], bool createrrv, int maxorderi, int maxorderj, int maxorderk, int maxtotorder, int maxcorrel, const TString& filename)
// {

//   const double dmin=-1.;
//   const double dmax=1.;
//   RooArgList varlist;
//   for (int i=0; i<ordermax3i;i++) {
//     for (int j=0; j<ordermax3j;j++) {
//       for (int k=0; k<ordermax3k;k++) {
// 	int correl=0;
// 	if (i>0) correl++;
// 	if (j>0) correl++;
// 	if (k>0) correl++;
	
// 	if (
// 	    ((i%2==0) or (i==1)) and 
// 	    //	    ((j%2==0) or (j==1)) and 
// 	    ((k%2==0) or (k==1)) and 
// 	    (i<maxorderi) and
// 	    (j<maxorderj) and
// 	    (k<maxorderk) and	
// 	    (correl<=maxcorrel) and
// 	    (i+j+k<=maxtotorder)) {
// 	  if (createrrv) {
// 	    coeff3[i][j][k] = new RooRealVar(suffix + "-c" + istr(i) + istr(j) + istr(k), suffix + "-c" + istr(i) + istr(j) + istr(k),0.,dmin,dmax);
// 	    //	    std::cout << "Adding coeff"<<i<<j<<k<< " " << coeff3[i][j][k] << std::endl;
// 	  }
// 	  varlist.add(*coeff3[i][j][k]);
	  
// 	} else {
// 	  coeff3[i][j][k]=NULL;
// 	}
//       }
//     }
//   }
//   //  std::cout << coeff3[0][0][0] << std::endl;
//   //  mygetchar;
//   //  coeff3[0][0][0]->setMax(1000.);
//   coeff3[0][0][0]->setVal(1.);
//   coeff3[0][0][0]->setConstant();
//   //  mygetchar;
//   if (filename!="") readcoeff3_rrv(filename,coeff3);

//   pdf = new RooChebychev3Pdf(suffix,suffix,costheta,costheta1,costheta2,varlist);//, ordermax3i, ordermax3j, ordermax3k);
//   std::cout << "Creating RooChebychev3" << std::endl;
//   varlist.Print("v");
// }


void createRooLbtoJpsiL0wAccPDF3(RooRealVar& costheta, RooRealVar& costheta1, RooRealVar& costheta2, RooRealVar& P_b, RooAbsReal& alpha_b,RooRealVar& alpha_lambda, RooAbsReal& r_0, RooAbsReal& r_1, TString suffix,RooAbsPdf*& pdf, RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k], bool createrrv, const TString& filename)
{
  const double dmin=-1.;
  const double dmax=1.;
  
  if (createrrv) {
    for (int i=0; i<ordermax3i;i++) {
      for (int j=0; j<ordermax3j;j++) {
	for (int k=0; k<ordermax3k;k++) {
	  coeff3[i][j][k]=NULL;
	  coeff3[i][j][k]=new RooRealVar(suffix + "-c"+ istr(i,"02") + istr(j,"02") + istr(k,"02"), suffix + "-c"+ istr(i,"02") + istr(j,"02") + istr(k,"02"),0.,dmin,dmax);
	}
      }
    }
    
    coeff3[0][0][0]->setVal(1.);
    coeff3[0][0][0]->setConstant();
   
  }

  if (filename!="") readcoeff3_rrv(filename,coeff3);

  RooArgList varlist;
  for (int i=0; i<ordermax3i;i++) {
    for (int j=0; j<ordermax3j;j++) {
      for (int k=0; k<ordermax3k;k++) {
	if (coeff3[i][j][k]) {
	  coeff3[i][j][k]->setConstant();
	  if (coeff3[i][j][k]->getVal()!=0.) varlist.add(*coeff3[i][j][k]);
	}
      }
    }
  }
  

  pdf = new RooLbtoJpsiL0wAccPDF3(suffix,suffix,costheta,costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1,varlist);

}




