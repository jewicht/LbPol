#include "findminmaxtf3.h"

double findminimum(const TF3& f3, double currmin, double x0, double x1, double y0, double y1, double z0, double z1)
{
  std::cout << "(TF3 min) searching in " << currmin <<" " << x0 << " " << x1 << " " << y0 << " " << y1 << "  " << z0 << " " << z1 << std::endl;
  const double epsilon=0.0001;
  const int idiv=99;
  const int jdiv=99;
  const int kdiv=99;
  double thisminimum=currmin;
  int imin=0;
  int jmin=0;
  int kmin=0;
  for (int i=0;i<idiv;i++) {
    for (int j=0;j<jdiv;j++) {
      for (int k=0;k<kdiv;k++) {
	double value=f3.Eval(x0+(x1-x0)*i/idiv,y0+(y1-y0)*j/jdiv, z0+(z1-z0)*k/kdiv);
	if (value<thisminimum) {
	  thisminimum=value;
	  imin=i;
	  jmin=j;
	  kmin=k;
	}
      }
    }
  }
  std::cout << "(TF3 min) thisminimum " << thisminimum << " at  " << x0+(x1-x0)*imin/idiv  << "  "<< y0+(y1-y0)*jmin/jdiv  << "  " << z0+(z1-z0)*kmin/kdiv  << "  " << std::endl
;
  if (fabs(currmin-thisminimum)<epsilon) return thisminimum;
  return findminimum(f3,thisminimum,
		     x0+(x1-x0)*(((imin==0)?0:imin-1))/idiv, 
		     x0+(x1-x0)*(((imin==idiv)?idiv:imin+1))/idiv,
		     y0+(y1-y0)*(((jmin==0)?0:jmin-1))/jdiv, 
		     y0+(y1-y0)*(((jmin==jdiv)?jdiv:jmin+1))/jdiv ,
		     z0+(z1-z0)*(((kmin==0)?0:kmin-1))/kdiv, 
		     z0+(z1-z0)*(((kmin==kdiv)?kdiv:kmin+1))/kdiv );

}

double findminimum2(TF3& f3)//, double currmin, double x0, double x1, double y0, double y1, double z0, double z1)
{
  double x(0.),y(0.),z(0.);
  f3.GetMinimumXYZ(x,y,z);
  return f3.Eval(x,y,z);
}


void GetMaximumXYZ(TF3& f3, Double_t &x, Double_t &y, Double_t &z)
{
// Return the X, Y and Z values corresponding to the minimum value of the function
// on its range. To find the minimum on a subrange, use the SetRange() function first.
// Method:
//   First, a grid search is performed to find the initial estimate of the 
//   minimum location. The range of the function is divided 
//   into fNpx,fNpy and fNpz sub-ranges. If the function is "good" (or "bad"), 
//   these values can be changed by SetNpx(), SetNpy() and SetNpz() functions.
//   Then, Minuit minimization is used with starting values found by the grid search


   //First do a grid search with step size fNpx adn fNpy

  Double_t fXmin=-1.;
  Double_t fXmax=+1.;
  Double_t fYmin=-1.;
  Double_t fYmax=+1.;
  Double_t fZmin=-1.;
  Double_t fZmax=+1.;

  Int_t fNpx=10;
  Int_t fNpy=10;  
  Int_t fNpz=10;

   Double_t xx, yy, zz, tt;
   Double_t dx = (fXmax - fXmin)/fNpx;
   Double_t dy = (fYmax - fYmin)/fNpy;
   Double_t dz = (fZmax - fZmin)/fNpz;

   Double_t xxmin = fXmin;
   Double_t yymin = fYmin;
   Double_t zzmin = fZmin;
   Double_t ttmin = f3.Eval(xxmin, yymin, zzmin+dz);
   for (Int_t i=0; i<fNpx; i++){
      xx=fXmin + (i+0.5)*dx;
      for (Int_t j=0; j<fNpy; j++){
         yy=fYmin+(j+0.5)*dy;
         for (Int_t k=0; k<fNpz; k++){
            zz = fZmin+(k+0.5)*dz;
            tt = f3.Eval(xx, yy, zz);
            if (tt>ttmin) {xxmin = xx, yymin = yy; zzmin = zz; ttmin=tt;}
         }
      }
   }

   x = TMath::Min(fXmax, xxmin);
   y = TMath::Min(fYmax, yymin);
   z = TMath::Min(fZmax, zzmin);

   //go to minuit for the final minimization
   char f[]="TFitter";

   Int_t strdiff = 0;
   if (TVirtualFitter::GetFitter()){
      //If the fitter is already set and it's not minuit, delete it and 
      //create a minuit fitter
      strdiff = strcmp(TVirtualFitter::GetFitter()->IsA()->GetName(), f);
      if (strdiff!=0) 
         delete TVirtualFitter::GetFitter();
   }

   TVirtualFitter *minuit = TVirtualFitter::Fitter(&f3, 3);
   minuit->Clear();
   minuit->SetFitMethod("F3Maximizer");
   Double_t arglist[10];
   arglist[0]=-1;
   minuit->ExecuteCommand("SET PRINT", arglist, 1);

   minuit->SetParameter(0, "x", x, 0.1, 0, 0);
   minuit->SetParameter(1, "y", y, 0.1, 0, 0);
   minuit->SetParameter(2, "z", z, 0.1, 0, 0);
   arglist[0] = 5;
   arglist[1] = 1e-5;
   // minuit->ExecuteCommand("CALL FCN", arglist, 1);

   Int_t fitResult = minuit->ExecuteCommand("MIGRAD", arglist, 0);
   if (fitResult!=0){
      //migrad might have not converged
     //      Warning("GetMinimumXYZ", "Abnormal termination of minimization");
   }
   Double_t xtemp = minuit->GetParameter(0);
   Double_t ytemp = minuit->GetParameter(1);
   Double_t ztemp = minuit->GetParameter(2);
   if (xtemp>fXmax || xtemp<fXmin || ytemp>fYmax || ytemp<fYmin || ztemp>fZmax || ztemp<fZmin){
      //converged to something outside limits, redo with bounds 
      minuit->SetParameter(0, "x", x, 0.1, fXmin, fXmax);
      minuit->SetParameter(1, "y", y, 0.1, fYmin, fYmax);
      minuit->SetParameter(2, "z", z, 0.1, fZmin, fZmax);
      fitResult = minuit->ExecuteCommand("MIGRAD", arglist, 0);
      if (fitResult!=0){
         //migrad might have not converged
	//         Warning("GetMinimumXYZ", "Abnormal termination of minimization");
      }
   }
   x = minuit->GetParameter(0);
   y = minuit->GetParameter(1);
   z = minuit->GetParameter(2);
}

double findmaximum2(TF3& f3)
{
  double x(0.),y(0.),z(0.);
  GetMaximumXYZ(f3,x,y,z);
  return f3.Eval(x,y,z);

 // const int npar = 3;              // the number of parameters
 //  TMinuit minuit(npar);
 //  minuit.SetFCN(fcn);

 //  double par[npar];               // the start values
 //  double stepSize[npar];          // step sizes 
 //  double minVal[npar];            // minimum bound on parameter 
 //  double maxVal[npar];            // maximum bound on parameter
 //  std::string parName[npar];

 //  for (int i=0; i<npar; i++){
 //    minuit.DefineParameter(i, parName[i].c_str(), 
 // 			   par[i], stepSize[i], minVal[i], maxVal[i]);
 //  }
  
 //  // Do the minimization!
  
 //  minuit.Migrad();       // Minuit's best minimization algorithm
 //  double outpar[npar], err[npar];
 //  for (int i=0; i<npar; i++){
 //    minuit.GetParameter(i,outpar[i],err[i]);
 //  }
  

}

double findmaximum(const TF3& f3, double currmax, double x0, double x1, double y0, double y1, double z0, double z1)
{
  std::cout << "(TF3 max) searching in " << currmax <<" " << x0 << " " << x1 << " " << y0 << " " << y1 << "  " << z0 << " " << z1 << std::endl;
  const double epsilon=0.0001;
  const int idiv=99;
  const int jdiv=99;
  const int kdiv=99;
  double thismaximum=currmax;
  int imax=0;
  int jmax=0;
  int kmax=0;
  for (int i=0;i<=idiv;i++) {
    for (int j=0;j<=jdiv;j++) {
      for (int k=0;k<=kdiv;k++) {
	double value=f3.Eval(x0+(x1-x0)*i/idiv,
			     y0+(y1-y0)*j/jdiv, 
			     z0+(z1-z0)*k/kdiv);
	//	std::cout << value << std::endl;

	if (value>thismaximum) {
	  thismaximum=value;
	  imax=i;
	  jmax=j;
	  kmax=k;
	}
      }
    }
  }
  
  //  int imaxm1=imax-1;
  std::cout << "(TF3 max) thismaximum " << thismaximum << " at  " << x0+(x1-x0)*imax/idiv  << "  "<< y0+(y1-y0)*jmax/jdiv  << "  " << z0+(z1-z0)*kmax/kdiv  << "  " << std::endl
;
  if (fabs(currmax-thismaximum)<epsilon) return thismaximum;
  return findmaximum(f3,thismaximum,
		     x0+(x1-x0)*(((imax==0)?0:imax-1))/idiv, 
		     x0+(x1-x0)*(((imax==idiv)?idiv:imax+1))/idiv,
		     y0+(y1-y0)*(((jmax==0)?0:jmax-1))/jdiv, 
		     y0+(y1-y0)*(((jmax==jdiv)?jdiv:jmax+1))/jdiv ,
		     z0+(z1-z0)*(((kmax==0)?0:kmax-1))/kdiv, 
		     z0+(z1-z0)*(((kmax==kdiv)?kdiv:kmax+1))/kdiv );
}

