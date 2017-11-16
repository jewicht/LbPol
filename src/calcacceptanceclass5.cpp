#include "calcacceptanceclass5.h"
#include "stdlib.h"

calcacceptanceclass5::calcacceptanceclass5(const TString& filename)
{
  std::ifstream filestr;
  filestr.open(TString("/scratch/z5/jwicht/fitLambdab2JpsiL0/") +  filename);

  if (!filestr.is_open()) {
    std::cerr << "calcacceptanceclass5::calcacceptanceclass5: file not found " << filename << std::endl;
    exit(1);
  }

  memset(dcoeff5, 0, ordermax5i*ordermax5j*ordermax5k*ordermax5l*ordermax5m*sizeof(double));
  //build TF3
  for (int i=0; i<ordermax5i;i++) {
    for (int j=0; j<ordermax5j;j++) {
      for (int k=0; k<ordermax5k;k++) {	
	for (int l=0; l<ordermax5l;l++) {	
	  for (int m=0; m<ordermax5m;m++) {	
	    double value, error;
	    bool isconstant;
	    int ii,jj,kk,ll,mm;
	    filestr >>  ii >> jj >> kk >> ll >> mm >> value >> error >> isconstant;// >> std::endl;
	    //	std::cout << ii << " " << jj << " " << kk << " " << value <<  " " << error << " " << isconstant << std::endl;
	    //	std::cout << "coeff3_" << ii << jj << kk << " is set" << std::endl; 
	    dcoeff5[ii][jj][kk][ll][mm] = value;
	    
	  }
	}
      }
    }
  }
  filestr.close();

  //  func0 = new TF3("func0",this, &calcacceptanceclass5::myfunction0, -1., 1., -1., 1., -1., 1., 0, "calcacceptanceclass5","myfunction0");
  //  func0 = new TF3("func0",myfunction0, -1., 1., -1., 1., -1., 1.);

  func0_min = findminimum(1000.,-1., 1., -1., 1., -1., 1.,-1.,1.,-1.,1.);
  //  double func0_min2 = findminimum2(*func0);
  //  std::cout << "min comp " << func0_min << " " << func0_min2 << std::endl;

  if (func0_min<0.) {
    std::cerr << "==========" << std::endl;
    std::cerr << "function goes below zero... file = " << filename << std::endl;
    std::cerr << "==========" << std::endl;
  }
  func0_max = findmaximum(-1000.,-1., 1., -1., 1., -1., 1.,-1.,1.,-1.,1.);

  
  TRandom3 rnd(0);
  double amax=-1000.;
  double amin=10000.;
  for (int i=0;i<100000;i++) {
    const double val=pdfval5f5(rnd.Rndm()*2.-1,
                               rnd.Rndm()*2.-1,
                               rnd.Rndm()*2.-1,
                               rnd.Rndm()*2.-1,
                               rnd.Rndm()*2.-1,
                               dcoeff5);
    if (val>amax) amax=val;
    if (val<amin) amin=val;
    if (val<func0_min) std::cerr << "minimum problem" << std::endl;
    if (val>func0_max) std::cerr << "maximum problem" << std::endl;
  }
  std::cout << "min - max = " << amin << " " << amax << std::endl;

  
  //  double func0_max2 = findmaximum2(*func0);
  //  std::cout << "max comp " << func0_max << " " << func0_max2 << std::endl;
  
}

// Double_t calcacceptanceclass5::myfunction0(double costheta, double costheta1, double costheta2, double phi1, double phi2)
// {
//   //  Double_t costheta  = x[0];
//   //  Double_t costheta1 = x[1];
//   //  Double_t costheta2 = x[2];
//   return pdfval5f5(costheta,costheta1,costheta2,phi1,phi2,dcoeff5);
// }


double calcacceptanceclass5::evaluate(double costheta, double costheta1, double costheta2, double phi1, double phi2)
{
  const double value = pdfval5f5(costheta,costheta1,costheta2,phi1,phi2,dcoeff5);
  //  return (value - func0_min) / (func0_max - func0_min);
  if (value<0.) return 99999999.;
  return value / func0_max;
}

double calcacceptanceclass5::findminimum(double currmin, double x0, double x1, double y0, double y1, double z0, double z1, double a0, double a1, double b0, double b1)
{
  std::cout << "(F5 min) searching in " << currmin <<" " << x0 << " " << x1 << " " << y0 << " " << y1 << "  " << z0 << " " << z1 << " " << a0 << " " << a1 << "  " << b0 << " " << b1 << std::endl;
  double thisminimum=currmin;
  int imin=0;
  int jmin=0;
  int kmin=0;
  int lmin=0;
  int mmin=0;
  for (int i=0;i<idiv;i++) {
    for (int j=0;j<jdiv;j++) {
      for (int k=0;k<kdiv;k++) {
	for (int l=0;l<ldiv;l++) {
	  for (int m=0;m<mdiv;m++) {
	    double value=pdfval5f5(x0+(x1-x0)*i/idiv,
				   y0+(y1-y0)*j/jdiv,
				   z0+(z1-z0)*k/kdiv,
				   a0+(a1-a0)*l/ldiv,
				   b0+(b1-b0)*m/mdiv,
				   dcoeff5
				   );
	    if (value<thisminimum) {
	      thisminimum=value;
	      imin=i;
	      jmin=j;
	      kmin=k;
	      lmin=l;
	      mmin=m;
	    }
	  }
	}
      }
    }
  }
  std::cout << "(F5 min) thisminimum " << thisminimum << " at  " << x0+(x1-x0)*imin/idiv  << "  "<< y0+(y1-y0)*jmin/jdiv  << "  " << z0+(z1-z0)*kmin/kdiv  << "  " << a0+(a1-a0)*lmin/ldiv  << "  " << b0+(b1-b0)*mmin/mdiv  << std::endl
;
  if (fabs(currmin-thisminimum)<epsilon) return thisminimum;
  return findminimum(thisminimum,
		     x0+(x1-x0)*(((imin==0)?0:imin-1))/idiv, 
		     x0+(x1-x0)*(((imin==idiv)?idiv:imin+1))/idiv,
		     y0+(y1-y0)*(((jmin==0)?0:jmin-1))/jdiv, 
		     y0+(y1-y0)*(((jmin==jdiv)?jdiv:jmin+1))/jdiv ,
		     z0+(z1-z0)*(((kmin==0)?0:kmin-1))/kdiv, 
		     z0+(z1-z0)*(((kmin==kdiv)?kdiv:kmin+1))/kdiv,
		     a0+(a1-a0)*(((lmin==0)?0:lmin-1))/ldiv, 
		     a0+(a1-a0)*(((lmin==ldiv)?ldiv:lmin+1))/ldiv, 
		     b0+(b1-b0)*(((mmin==0)?0:mmin-1))/mdiv, 
		     b0+(b1-b0)*(((mmin==mdiv)?mdiv:mmin+1))/mdiv 
		     );
  
}

double calcacceptanceclass5::findmaximum(double currmax, double x0, double x1, double y0, double y1, double z0, double z1, double a0, double a1, double b0, double b1)
{
  std::cout << "(F5 max) searching in " << currmax <<" " << x0 << " " << x1 << " " << y0 << " " << y1 << "  " << z0 << " " << z1 << " " << a0 << " " << a1 << "  " << b0 << " " << b1 << std::endl;
  double thismaximum=currmax;
  int imin=0;
  int jmin=0;
  int kmin=0;
  int lmin=0;
  int mmin=0;
  for (int i=0;i<idiv;i++) {
    for (int j=0;j<jdiv;j++) {
      for (int k=0;k<kdiv;k++) {
	for (int l=0;l<ldiv;l++) {
	  for (int m=0;m<mdiv;m++) {
	    double value=pdfval5f5(x0+(x1-x0)*i/idiv,
				   y0+(y1-y0)*j/jdiv,
				   z0+(z1-z0)*k/kdiv,
				   a0+(a1-a0)*l/ldiv,
				   b0+(b1-b0)*m/mdiv,
				   dcoeff5
				   );
	    if (value>thismaximum) {
	      thismaximum=value;
	      imin=i;
	      jmin=j;
	      kmin=k;
	      lmin=l;
	      mmin=m;
	    }
	  }
	}
      }
    }
  }
  std::cout << "(F5 max) thismaximum " << thismaximum << " at  " << x0+(x1-x0)*imin/idiv  << "  "<< y0+(y1-y0)*jmin/jdiv  << "  " << z0+(z1-z0)*kmin/kdiv  << "  " << a0+(a1-a0)*lmin/ldiv  << "  " << b0+(b1-b0)*mmin/mdiv  << std::endl
;
  if (fabs(currmax-thismaximum)<epsilon) return thismaximum;
  return findmaximum(thismaximum,
		     x0+(x1-x0)*(((imin==0)?0:imin-1))/idiv, 
		     x0+(x1-x0)*(((imin==idiv)?idiv:imin+1))/idiv,
		     y0+(y1-y0)*(((jmin==0)?0:jmin-1))/jdiv, 
		     y0+(y1-y0)*(((jmin==jdiv)?jdiv:jmin+1))/jdiv ,
		     z0+(z1-z0)*(((kmin==0)?0:kmin-1))/kdiv, 
		     z0+(z1-z0)*(((kmin==kdiv)?kdiv:kmin+1))/kdiv,
		     a0+(a1-a0)*(((lmin==0)?0:lmin-1))/ldiv, 
		     a0+(a1-a0)*(((lmin==ldiv)?ldiv:lmin+1))/ldiv, 
		     b0+(b1-b0)*(((mmin==0)?0:mmin-1))/mdiv, 
		     b0+(b1-b0)*(((mmin==mdiv)?mdiv:mmin+1))/mdiv 
		     );
  
}
