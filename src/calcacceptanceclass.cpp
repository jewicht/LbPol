#include "calcacceptanceclass.h"

calcacceptanceclass::calcacceptanceclass(const TString& _filename)
{
  debug=true;

  filename=_filename;
  memset(dcoeff3, 0, ordermax3i*ordermax3j*ordermax3k*sizeof(double));
  std::cout << "calcacceptanceclass: reading acceptance from " << filename << std::endl; 
  setcoeff(false);
}

calcacceptanceclass::~calcacceptanceclass()
{

}

void calcacceptanceclass::setcoeff(bool rand)
{
  RooFitResult* fr = getRooFitResult(filename);
  if (!fr) {
    std::cerr << "calcacceptanceclass: can't find RooFitResult in " << filename << std::endl;
    exit(1);
  }

  const RooArgList *floatPars,*constPars;
  if (rand) {
    floatPars=&fr->randomizePars();
  } else {
    floatPars=&fr->floatParsFinal();
  }
  constPars=&fr->constPars();

  const int len=strlen(floatPars->at(0)->GetName());
  
  for (int i=0;i<floatPars->getSize();i++) {
    const char* buf=floatPars->at(i)->GetName();
    const int vari=10*chtoint(buf[len-6])+chtoint(buf[len-5]);
    const int varj=10*chtoint(buf[len-4])+chtoint(buf[len-3]);
    const int vark=10*chtoint(buf[len-2])+chtoint(buf[len-1]);
    RooRealVar* var = (RooRealVar*)floatPars->at(i); 
    dcoeff3[vari][varj][vark]=var->getVal();
  }

  for (int i=0;i<constPars->getSize();i++) {
    const char* buf=constPars->at(i)->GetName();
    const int vari=10*chtoint(buf[len-6])+chtoint(buf[len-5]);
    const int varj=10*chtoint(buf[len-4])+chtoint(buf[len-3]);
    const int vark=10*chtoint(buf[len-2])+chtoint(buf[len-1]);
    RooRealVar* var = (RooRealVar*)constPars->at(i); 
    dcoeff3[vari][varj][vark]=var->getVal();
  }


  func0_min=1000.;
  double xmin,ymin,zmin;
  findminimum(func0_min, xmin, ymin,zmin,-1., 1., -1., 1., -1., 1.);

  if (func0_min<0.) {
    std::cerr << "==========" << std::endl;
    std::cerr << "function goes below zero... file = " << filename << std::endl;
    std::cout << "(F3 min)    minimum = " << func0_min << " at  " << xmin  << " , "<< ymin  << " , " << zmin  << std::endl;
    std::cerr << "==========" << std::endl;
  }
  func0_max=-1000.;
  double xmax,ymax,zmax;
  findmaximum(func0_max,xmax,ymax,zmax,-1., 1., -1., 1., -1., 1.);

  double integralval=integral(-1., 1., -1., 1., -1., 1.);
  std::cerr << "min = " << func0_min << " ; max = " << func0_max << " ; integral = " << integralval << std::endl;

  //  func0_max=8.;
}


double calcacceptanceclass::evaluate(double costheta, double costheta1, double costheta2)
{
  //  const double value = chebychev3f3(costheta,costheta1,costheta2,dcoeff3);
  const double value = pdfval(costheta,costheta1,costheta2,dcoeff3);
  //  return (value - func0_min) / (func0_max - func0_min);
  //  if (value<0.) {
    //    std::cerr << "Problem in " << filename << " at " << costheta << " " << costheta1 << " " << costheta2 << std::endl; 
    //    mygetchar;
  //    return NAN;
    //    return 0.;
  //  }
  //  return value / 8.;
  return value / func0_max;
}

void calcacceptanceclass::findminimum(double& currmin, double& xmin, double& ymin, double& zmin, double x0, double x1, double y0, double y1, double z0, double z1)
{
 if (debug)  std::cout << "(F3 min) searching in " << currmin <<" " << x0 << " " << x1 << " " << y0 << " " << y1 << "  " << z0 << " " << z1 << std::endl;
  double thisminimum=currmin;
  int imin=0;
  int jmin=0;
  int kmin=0;
  for (int i=0;i<idiv;i++) {
    for (int j=0;j<jdiv;j++) {
      for (int k=0;k<kdiv;k++) {
	//	double value=chebychev3f3(
	double value=pdfval(
			    x0+(x1-x0)*i/idiv,
			    y0+(y1-y0)*j/jdiv,
			    z0+(z1-z0)*k/kdiv,
			    dcoeff3
			    );
	if (value<thisminimum) {
	  thisminimum=value;
	  imin=i;
	  jmin=j;
	  kmin=k;
	  xmin=x0+(x1-x0)*imin/idiv;
	  ymin=y0+(y1-y0)*jmin/jdiv;
	  zmin=z0+(z1-z0)*kmin/kdiv;
	}
      }
    }
  }
  if (debug) std::cout << "(F3 min)    minimum = " << thisminimum << " at  " << xmin  << " , "<< ymin  << " , " << zmin  << std::endl;
  if (fabs(currmin-thisminimum)<epsilon) {
    currmin=thisminimum;
    return;// thisminimum;
  }
  currmin=thisminimum;
  findminimum(currmin,
	      xmin,ymin,zmin,
	      x0+(x1-x0)*(((imin==0)?0:imin-1))/idiv, 
	      x0+(x1-x0)*(((imin==idiv)?idiv:imin+1))/idiv,
	      y0+(y1-y0)*(((jmin==0)?0:jmin-1))/jdiv, 
	      y0+(y1-y0)*(((jmin==jdiv)?jdiv:jmin+1))/jdiv ,
	      z0+(z1-z0)*(((kmin==0)?0:kmin-1))/kdiv, 
	      z0+(z1-z0)*(((kmin==kdiv)?kdiv:kmin+1))/kdiv
	      );
}


double calcacceptanceclass::integral(double x0, double x1, double y0, double y1, double z0, double z1)
{
  const int idiv0=50;
  const int jdiv0=50;
  const int kdiv0=50;
  double value=0.;

  const double cellsize= 1.0
    *(x1-x0)/idiv0
    *(y1-y0)/jdiv0 
    *(z1-z0)/kdiv0;

  //  std::cout << "cellsize = " << cellsize << std::endl;

  for (int i=0;i<idiv0;i++) {
    for (int j=0;j<jdiv0;j++) {
      for (int k=0;k<kdiv0;k++) {
	value+=pdfval(
		      x0+(x1-x0)*(0.5+i)/idiv0,
		      y0+(y1-y0)*(0.5+j)/jdiv0,
		      z0+(z1-z0)*(0.5+k)/kdiv0,
		      dcoeff3
		      );
      }
    }
  }
  return value*cellsize;
}

void calcacceptanceclass::findmaximum(double& currmax, double& xmax, double& ymax, double& zmax, double x0, double x1, double y0, double y1, double z0, double z1)
{
  if (debug) std::cout << "(F3 max) searching in " << currmax <<" " << x0 << " " << x1 << " " << y0 << " " << y1 << "  " << z0 << " " << z1 << std::endl;
  double thismaximum=currmax;
  int imin=0;
  int jmin=0;
  int kmin=0;
  for (int i=0;i<idiv;i++) {
    for (int j=0;j<jdiv;j++) {
      for (int k=0;k<kdiv;k++) {
	//	double value=chebychev3f3(
	double value=pdfval(
			    x0+(x1-x0)*i/idiv,
			    y0+(y1-y0)*j/jdiv,
			    z0+(z1-z0)*k/kdiv,
			    dcoeff3
			    );
	if (value>thismaximum) {
	  thismaximum=value;
	  imin=i;
	  jmin=j;
	  kmin=k;
	  xmax=x0+(x1-x0)*imin/idiv;
	  ymax=y0+(y1-y0)*jmin/jdiv;
	  zmax=z0+(z1-z0)*kmin/kdiv;
	}
      }
    }
  }

  if (debug) std::cout << "(F3 max)     maximum = " << thismaximum << " at  " << xmax  << "  "<< ymax  << "  " << zmax  << std::endl;
  if (fabs(currmax-thismaximum)<epsilon) {
    currmax=thismaximum;
    return;
  }
  currmax=thismaximum;
  findmaximum(currmax,
	      xmax,ymax,zmax,
	      x0+(x1-x0)*(((imin==0)?0:imin-1))/idiv, 
	      x0+(x1-x0)*(((imin==idiv)?idiv:imin+1))/idiv,
	      y0+(y1-y0)*(((jmin==0)?0:jmin-1))/jdiv, 
	      y0+(y1-y0)*(((jmin==jdiv)?jdiv:jmin+1))/jdiv ,
	      z0+(z1-z0)*(((kmin==0)?0:kmin-1))/kdiv, 
	      z0+(z1-z0)*(((kmin==kdiv)?kdiv:kmin+1))/kdiv
	      );
  
}
// double calcacceptanceclass::findmaximum(double currmax, double x0, double x1, double y0, double y1, double z0, double z1)
// {
//   if (debug) std::cout << "(F3 max) searching in " << currmax <<" " << x0 << " " << x1 << " " << y0 << " " << y1 << "  " << z0 << " " << z1 << std::endl;
//   double thismaximum=currmax;
//   int imin=0;
//   int jmin=0;
//   int kmin=0;
//   for (int i=0;i<idiv;i++) {
//     for (int j=0;j<jdiv;j++) {
//       for (int k=0;k<kdiv;k++) {
// 	//	double value=chebychev3f3(
// 	double value=pdfval(
// 			    x0+(x1-x0)*i/idiv,
// 			    y0+(y1-y0)*j/jdiv,
// 			    z0+(z1-z0)*k/kdiv,
// 			    dcoeff3
// 			    );
// 	if (value>thismaximum) {
// 	  thismaximum=value;
// 	  imin=i;
// 	  jmin=j;
// 	  kmin=k;
// 	}
//       }
//     }
//   }

//   if (debug) std::cout << "(F3 max)     maximum = " << thismaximum << " at  " << x0+(x1-x0)*imin/idiv  << "  "<< y0+(y1-y0)*jmin/jdiv  << "  " << z0+(z1-z0)*kmin/kdiv  << "  " << std::endl;
//   if (fabs(currmax-thismaximum)<epsilon) return thismaximum;
//   return findmaximum(thismaximum,
// 		     x0+(x1-x0)*(((imin==0)?0:imin-1))/idiv, 
// 		     x0+(x1-x0)*(((imin==idiv)?idiv:imin+1))/idiv,
// 		     y0+(y1-y0)*(((jmin==0)?0:jmin-1))/jdiv, 
// 		     y0+(y1-y0)*(((jmin==jdiv)?jdiv:jmin+1))/jdiv ,
// 		     z0+(z1-z0)*(((kmin==0)?0:kmin-1))/kdiv, 
// 		     z0+(z1-z0)*(((kmin==kdiv)?kdiv:kmin+1))/kdiv
// 		     );
  
// }
