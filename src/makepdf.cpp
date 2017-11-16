#include <iostream>
#include <fstream>

#include "RooClassFactory.h"
#include "TString.h"

#include "RooRealVar.h"
#include "RooLegendre.h"
#include "RooProdPdf.h"

#include "TRandom3.h"

TString istr(int x,TString fill="", bool zerohack=false)
{
  if (zerohack and (x==0)) return "--";
  char buffer1[50];
  //  int n = 
  sprintf(buffer1, "%" + fill + "d", x);
  return TString(buffer1);
}


void legendre3()
{

  const int ordermaxi=4;
  const int ordermaxj=4;
  const int ordermaxk=7;

  TString varlist="costheta0,costheta1,costheta2";
  for (int i=0; i<ordermaxi;i++) {
    for (int j=0; j<ordermaxj;j++) {
      for (int k=0; k<ordermaxk;k++) {
	if (((i%2==0) or (i==1)) and 
	    ((j%2==0) or (j==1)) and 
	    ((k%2==0) or (k==1))) {
	  if (i+j+k<=12) varlist+=",c" + istr(i) + istr(j) + istr(k);
  	}
      }  
    }
  }
  //RooClassFactory::makePdf("RooLegendre3Pdf",varlist,0,"1.",true);
  

  std::cout << "static const int ordermaxi=" << ordermaxi << ";" << std::endl;
  std::cout << "static const int ordermaxj=" << ordermaxj << ";" << std::endl;
  std::cout << "static const int ordermaxk=" << ordermaxk << ";" << std::endl;
  std::cout << "double coeff3[ordermaxi][ordermaxj][ordermaxk];" << std::endl;
  
  for (int i=0; i<ordermaxi;i++) {
    for (int j=0; j<ordermaxj;j++) {
      for (int k=0; k<ordermaxk;k++) {
	if (((i%2==0) or (i==1)) and 
	    ((j%2==0) or (j==1)) and 
	    ((k%2==0) or (k==1))) {
	  if (i+j+k<=12) std::cout << "coeff3[" << i << "]["<< j << "]["<< k << "] = c" << i << j << k << ";" << std::endl; 
  	}
      }  
    }
  }

  TRandom3 rnd;

  for (int i=0; i<ordermaxi;i++) {
    for (int j=0; j<ordermaxj;j++) {
      for (int k=0; k<ordermaxk;k++) {
	if (((i%2==0) or (i==1)) and 
	    ((j%2==0) or (j==1)) and 
	    ((k%2==0) or (k==1))) {
	  if (i+j+k<=12) std::cout << "RooRealVar c" << i << j << k << "(\"c" <<  i << j << k << "\", \"c" << i << j << k << "\"," << 0. <<   ",-1.,1.);" << std::endl; 
	}
      }  
    }
  }


  std::cout << varlist << std::endl;


}



void legendre5()
{
  const int ordermaxi=4;
  const int ordermaxj=4;
  const int ordermaxk=4;
  const int ordermaxl=4;
  const int ordermaxm=7;

  TString varlist="costheta0,costheta1,costheta2,phi1,phi2";
  TString coefflist="";
  TString func5="";
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
	      if (i+j+k+l+m<=12) { 
		varlist+=",c" + istr(i) + istr(j) + istr(k) + istr(l) + istr(m);
		coefflist+=" c" + istr(i) + istr(j) + istr(k) + istr(l) + istr(m);
		func5+="+c"+ istr(i) + istr(j) + istr(k) + istr(l) + istr(m) + "*legendre_P(" + istr(i) + ",costheta0)*legendre_P(" + istr(j) + ",costheta1)*legendre_P(" + istr(k) + ",costheta2)*legendre_P(" + istr(l) + ",phi1)*legendre_P(" + istr(m) + ",phi2)";
	      }
	    }
  	  }
  	}
      }  
    }
  }


  std::fstream filestr;
  filestr.open("vars5.txt", std::fstream::out);
  filestr << coefflist << std::endl;
  filestr.close();
  filestr.open("func5.txt", std::fstream::out);
  filestr << func5 << std::endl;
  filestr.close();
  exit(1);


  //  RooClassFactory::makePdf("RooLegendre5Pdf",varlist,0,"1.",true);
  

  std::cout << "static const int ordermax5i=" << ordermaxi << ";" << std::endl;
  std::cout << "static const int ordermax5j=" << ordermaxj << ";" << std::endl;
  std::cout << "static const int ordermax5k=" << ordermaxk << ";" << std::endl;
  std::cout << "static const int ordermax5l=" << ordermaxl << ";" << std::endl;
  std::cout << "static const int ordermax5m=" << ordermaxm << ";" << std::endl;
  std::cout << "double coeff5[" << ordermaxi <<"]["<<ordermaxj<<"]["<<ordermaxk<<"]["<<ordermaxl<<"]["<<ordermaxm<<"];" << std::endl;
  
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
	       if (i+j+k+l+m<=12) std::cout << "coeff5[" << i << "]["<< j << "]["<< k << "]["<< l << "]["<< m << "] = c" << i << j << k << l << m << ";" << std::endl; 
	    }
	  }
  	}
      }  
    }
  }

  TRandom3 rnd;

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

	      //	      std::cout << "RooRealVar c" << i << j << k << l << m << "(\"" <<  i << j << k << l << m << "\", \"" << i << j << k << l << m << "\"," << rnd.Rndm()*2.-1. <<   ",-1.,1.);" << std::endl; 
	      if (i+j+k+l+m<=12) std::cout << "RooRealVar c" << i << j << k << l << m << "(\"c" <<  i << j << k << l << m << "\", \"c" << i << j << k << l << m << "\",0.,-100000.,100000.);" << std::endl; 
	    }
	  }
  	}
      }  
    }
  }


  std::cout << varlist << std::endl;


}



int main(int argc, char** argv)
{
  //  legendre3();exit(1);
  legendre5();exit(1);

  // static const double pi = 3.14159;

  // RooRealVar costheta("costheta",   "cos(#theta)",     -1., 1.);//RooArgList(theta));
  // RooRealVar costheta1("costheta1", "cos(#theta_{1})", -1., 1.);// RooArgList(theta1));
  // RooRealVar costheta2("costheta2", "cos(#theta_{2})", -1., 1.);// RooArgList(theta2));

  // RooRealVar phi1("phi1","#phi_{1}",-1.,1.);
  // RooRealVar phi2("phi2","#phi_{2}",-1.,1.);



  // RooLegendre leg1("leg1","leg1", costheta, 1);
  // RooLegendre leg2("leg2","leg2", costheta1, 1);
  // RooLegendre leg3("leg3","leg3", costheta2, 1);
  // RooLegendre leg4("leg4","leg4", phi1, 1);
  // RooLegendre leg5("leg5","leg5", phi2, 1);

  // RooProdPdf prod("prod","prod",RooArgList(leg1,leg2));
  //  RooProdPdf prod("prod","prod",RooArgList(leg1,leg2,leg3,leg4,leg5));


  // const int ordermaxi=5;
  // const int ordermaxj=5;
  // const int ordermaxk=5;
  // const int ordermaxl=5;
  // const int ordermaxm=5;
  // RooLegendre* legendrecostheta0[ordermaxi];
  // RooLegendre* legendrecostheta1[ordermaxj];
  // RooLegendre* legendrecostheta2[ordermaxk];
  // RooLegendre* legendrephi1[ordermaxl];
  // RooLegendre* legendrephi2[ordermaxm];

  // for (int i=0; i<ordermaxi;i++) legendrecostheta0[i] = new RooLegendre("legendrecostheta0"+istr(i),"legendrecostheta0"+istr(i), costheta,  i);
  // for (int j=0; j<ordermaxj;j++) legendrecostheta1[j] = new RooLegendre("legendrecostheta1"+istr(j),"legendrecostheta1"+istr(j), costheta1, j);
  // for (int k=0; k<ordermaxk;k++) legendrecostheta2[k] = new RooLegendre("legendrecostheta2"+istr(k),"legendrecostheta2"+istr(k), costheta2, k);
  // for (int l=0; l<ordermaxl;l++) legendrephi1[l]      = new RooLegendre("legendrephi1"+istr(l),     "legendrephi1"+istr(l),      phi1,      l);
  // for (int m=0; m<ordermaxm;m++) legendrephi2[m]      = new RooLegendre("legendrephi2"+istr(m),     "legendrephi2"+istr(m),      phi2,      m);

  // RooProdPdf* legendreprod[ordermaxi][ordermaxj][ordermaxk][ordermaxl][ordermaxm];
  // RooRealVar* coeff5[ordermaxi][ordermaxj][ordermaxk][ordermaxl][ordermaxm];
  // for (int i=0; i<ordermaxi;i++) {
  //   for (int j=0; j<ordermaxj;j++) {
  //     for (int k=0; k<ordermaxk;k++) {
  //  	for (int l=0; l<ordermaxl;l++) {
  //  	  for (int m=0; m<ordermaxm;m++) {
	    

  // 	    std::cout << "ijklm " << i << " "<< j << " "<< k << " "<< l << " "<< m << std::endl; 
  // 	    std::cout << legendrecostheta0[i]->getVal() << std::endl;
  // 	    //	    std::cout << leg1.getVal() << std::endl;
  // 	    std::cout << legendrecostheta1[j]->GetName() << std::endl;
  // 	    std::cout << legendrecostheta2[k]->GetName() << std::endl;
  // 	    std::cout << legendrephi1[l]->GetName() << std::endl;
  // 	    std::cout << legendrephi2[m]->GetName() << std::endl;

  // 	    coeff5[i][j][k][l][m] = new RooRealVar("coeff5" + istr(i)  + istr(j)  + istr(k)  + istr(l) + istr(m),
  // 						  "coeff5" + istr(i)  + istr(j)  + istr(k)  + istr(l) + istr(m),
  // 						  0.,-1.,1.);

  // 	    TString pdftitle = "legendreprod" + istr(i)  + istr(j)  + istr(k)  + istr(l) + istr(m);
  // 	    std::cout << pdftitle << std::endl;
	    
  // 	    legendreprod[i][j][k][l][m] = new RooProdPdf(pdftitle,pdftitle,
  // 	    						 RooArgList( *legendrecostheta0[i], 
  // 								     *legendrecostheta1[j],  
  // 								     *legendrecostheta2[k],  
  // 								     *legendrephi1[l], 
  // 								     *legendrephi2[m] 
  // 								     ) 
  // 							 );



  // 	    //	    legendreprod[i][j][k][l][m] = new RooProdPdf(pdftitle,pdftitle,
  // 	    //						 RooArgList( leg1,leg2,leg3,leg4,leg5) );


  // 	  }
  // 	}
  //     }  
  //   }
  //  } 

}
