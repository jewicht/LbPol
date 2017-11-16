#include "computew3.h"


computew3::computew3()
{
}

bool computew3::checknan()
{
  return false;
  if (r_0+r_1<0.) return true;
  if (r_0-r_1<0.) return true;
  if (alpha_b+1.-r_0-r_1<0.) return true;
  if (1.-alpha_b-r_0+r_1<0.) return true;
  return false;
}

// TComplex computew3::a_plus(double _alpha_b, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi) const
// {
//   return sqrt(0.5*(_r_0+_r_1)) * TComplex::Exp(TComplex::I() * _alpha_plus);
// }

// TComplex computew3::a_minus(double _alpha_b, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi) const
// {
//   return sqrt(0.5*(_r_0-_r_1)) * TComplex::Exp(TComplex::I() * _alpha_minus);
// }

// TComplex computew3::b_plus(double _alpha_b, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi) const
// {
//   return sqrt(fabs(0.5*(_alpha_b+1.-_r_0-_r_1))) * TComplex::Exp(TComplex::I() * beta_plus);
// }

// TComplex computew3::b_minus(double _alpha_b, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi) const
// {
//   return sqrt(0.5*(1.-_alpha_b-_r_0+_r_1)) * TComplex::Exp(TComplex::I() * (_alpha_plus + _alpha_minus - beta_plus - _chi));
// }


double computew3::a_plus_sq(double _alpha_b, double _r_0, double _r_1) const
{
  return 0.5*(_r_0+_r_1);
}

double computew3::a_minus_sq(double _alpha_b, double _r_0, double _r_1) const
{
  return 0.5*(_r_0-_r_1);
}

double computew3::b_plus_sq(double _alpha_b, double _r_0, double _r_1) const
{
  return 0.5*(_alpha_b+1.-_r_0-_r_1);
}

double computew3::b_minus_sq(double _alpha_b, double _r_0, double _r_1) const
{
  return 0.5*(1.-_alpha_b-_r_0+_r_1);
}
//beta_minus = alphaplus + alphaminus - betaplus - chi
//chi = alphaplus + alphaminus - betaplus - betaminus
//alphaplus+alphaminus+betaplus+betaminus=2pi

double computew3::F(int i) const
{
  switch (i) {
    
  case  0 : return  1. ;break;
  case  1 : return  costheta0 ;break;
  case  2 : return  costheta1 ;break;
  case  3 : return  costheta0*costheta1 ;break;
  case  4 : return  1.5*costheta2pow2 - 0.5 ;break;
  case  5 : return  0.5*(3.*costheta2pow2 - 1.)*costheta0 ;break;
  case  6 : return  0.5*(3.*costheta2pow2 - 1.)*costheta1 ;break;
  case  7 : return  0.5*(3.*costheta2pow2 - 1.)*costheta0*costheta1 ;break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


double computew3::F_pi_costheta0(int i) const
{
  switch (i) {
  case  0 : return  costheta0 ;break;
  case  1 : return  0.5*costheta0pow2 ;break;
  case  2 : return  costheta0*costheta1 ;break;
  case  3 : return  0.5*costheta0pow2*costheta1 ;break;
  case  4 : return  0.5*(3.*costheta2pow2 - 1.)*costheta0 ;break;
  case  5 : return  0.25*(3.*costheta2pow2 - 1.)*costheta0pow2 ;break;
  case  6 : return  0.5*(3.*costheta2pow2 - 1.)*costheta0*costheta1 ;break;
  case  7 : return  0.25*(3.*costheta2pow2 - 1.)*costheta0pow2*costheta1 ;break;
    
  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}

double computew3::F_pi_costheta1(int i) const
{
  switch (i) {
  
  case  0 : return  costheta1 ;break;
  case  1 : return  costheta0*costheta1 ;break;
  case  2 : return  0.5*costheta1pow2 ;break;
  case  3 : return  0.5*costheta0*costheta1pow2 ;break;
  case  4 : return  0.5*(3.*costheta2pow2 - 1.)*costheta1 ;break;
  case  5 : return  0.5*(3.*costheta2pow2 - 1.)*costheta0*costheta1 ;break;
  case  6 : return  0.25*(3.*costheta2pow2 - 1.)*costheta1pow2 ;break;
  case  7 : return  0.25*(3.*costheta2pow2 - 1.)*costheta0*costheta1pow2 ;break;
  
  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


double computew3::F_pi_costheta2(int i) const
{
  switch (i) {

  case  0 : return  costheta2 ;break;
  case  1 : return  costheta0*costheta2 ;break;
  case  2 : return  costheta1*costheta2 ;break;
  case  3 : return  costheta0*costheta1*costheta2 ;break;
  case  4 : return  0.5*costheta2pow3 - 0.5*costheta2 ;break;
  case  5 : return  0.5*(costheta2pow3 - costheta2)*costheta0 ;break;
  case  6 : return  0.5*(costheta2pow3 - costheta2)*costheta1 ;break;
  case  7 : return  0.5*(costheta2pow3 - costheta2)*costheta0*costheta1 ;break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}

double computew3::F_pi_costheta0_costheta1(int i) const
{
  switch (i) {
    	
  case  0 : return  costheta0*costheta1 ;break;
  case  1 : return  0.5*costheta0pow2*costheta1 ;break;
  case  2 : return  0.5*costheta0*costheta1pow2 ;break;
  case  3 : return  0.25*costheta0pow2*costheta1pow2 ;break;
  case  4 : return  0.5*(3.*costheta2pow2 - 1.)*costheta0*costheta1 ;break;
  case  5 : return  0.25*(3.*costheta2pow2 - 1.)*costheta0pow2*costheta1 ;break;
  case  6 : return  0.25*(3.*costheta2pow2 - 1.)*costheta0*costheta1pow2 ;break;
  case  7 : return  0.125*(3.*costheta2pow2 - 1.)*costheta0pow2*costheta1pow2 ;break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}

double computew3::F_pi_costheta0_costheta2(int i) const
{
  switch (i) {

  case  0 : return  costheta0*costheta2 ;break;
  case  1 : return  0.5*costheta0pow2*costheta2 ;break;
  case  2 : return  costheta0*costheta1*costheta2 ;break;
  case  3 : return  0.5*costheta0pow2*costheta1*costheta2 ;break;
  case  4 : return  0.5*(costheta2pow3 - costheta2)*costheta0 ;break;
  case  5 : return  0.25*(costheta2pow3 - costheta2)*costheta0pow2 ;break;
  case  6 : return  0.5*(costheta2pow3 - costheta2)*costheta0*costheta1 ;break;
  case  7 : return  0.25*(costheta2pow3 - costheta2)*costheta0pow2*costheta1 ;break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


double computew3::F_pi_costheta1_costheta2(int i) const
{
  switch (i) {
    
  case  0 : return  costheta1*costheta2 ;break;
  case  1 : return  costheta0*costheta1*costheta2 ;break;
  case  2 : return  0.5*costheta1pow2*costheta2 ;break;
  case  3 : return  0.5*costheta0*costheta1pow2*costheta2 ;break;
  case  4 : return  0.5*(costheta2pow3 - costheta2)*costheta1 ;break;
  case  5 : return  0.5*(costheta2pow3 - costheta2)*costheta0*costheta1 ;break;
  case  6 : return  0.25*(costheta2pow3 - costheta2)*costheta1pow2 ;break;
  case  7 : return  0.25*(costheta2pow3 - costheta2)*costheta0*costheta1pow2 ;break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}



double computew3::F_primitive(int i) const
{
  switch (i) {
        	
  case  0 : return  costheta0*costheta1*costheta2 ;break;
  case  1 : return  0.5*costheta0pow2*costheta1*costheta2 ;break;
  case  2 : return  0.5*costheta0*costheta1pow2*costheta2 ;break;
  case  3 : return  0.25*costheta0pow2*costheta1pow2*costheta2 ;break;
  case  4 : return  0.5*(costheta2pow3 - costheta2)*costheta0*costheta1 ;break;
  case  5 : return  0.25*(costheta2pow3 - costheta2)*costheta0pow2*costheta1 ;break;
  case  6 : return  0.25*(costheta2pow3 - costheta2)*costheta0*costheta1pow2 ;break;
  case  7 : return  0.125*(costheta2pow3 - costheta2)*costheta0pow2*costheta1pow2 ;break;
  
  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


double computew3::f2(int i) const
{
  switch(i) {
  case  0: return 1.; break;
  case  1: return P_b; break;
  case  2: return alpha_lambda; break;
  case  3: return P_b*alpha_lambda; break;
  case  4: return 1.; break;
  case  5: return P_b; break;
  case  6: return alpha_lambda; break;
  case  7: return P_b*alpha_lambda; break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


double computew3::f2_pi_P_b(int i) const
{
  switch(i) {
  case  0: return P_b; break;
  case  1: return 0.5*pow(P_b,2); break;
  case  2: return alpha_lambda*P_b; break;
  case  3: return 0.5*pow(P_b,2)*alpha_lambda; break;
  case  4: return P_b; break;
  case  5: return 0.5*pow(P_b,2); break;
  case  6: return alpha_lambda*P_b; break;
  case  7: return 0.5*pow(P_b,2)*alpha_lambda; break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


double computew3::f2_pi_alpha_lambda(int i) const
{
  switch(i) {
  case  0: return alpha_lambda; break;
  case  1: return P_b*alpha_lambda; break;
  case  2: return 0.5*pow(alpha_lambda,2); break;
  case  3: return P_b*0.5*pow(alpha_lambda,2); break;
  case  4: return alpha_lambda; break;
  case  5: return P_b*alpha_lambda; break;
  case  6: return 0.5*pow(alpha_lambda,2); break;
  case  7: return P_b*0.5*pow(alpha_lambda,2); break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}



// double computew3::f1(int i) const
// {
//   switch (i) {
//   case  0: return  aplusrho2 + aminusrho2 + bplusrho2 + bminusrho2;break;
//   case  1: return  aplusrho2 - aminusrho2 + bplusrho2 - bminusrho2;break;
//   case  2: return  aplusrho2 - aminusrho2 - bplusrho2 + bminusrho2;break;
//   case  3: return  aplusrho2 + aminusrho2 - bplusrho2 - bminusrho2;break;

//   case  4: return -aplusrho2 - aminusrho2 + 0.5*bplusrho2 + 0.5*bminusrho2;break;
//   case  5: return -aplusrho2 + aminusrho2 + 0.5*bplusrho2 - 0.5*bminusrho2;break;
//   case  6: return -aplusrho2 + aminusrho2 - 0.5*bplusrho2 + 0.5*bminusrho2;break;
//   case  7: return -aplusrho2 - aminusrho2 - 0.5*bplusrho2 - 0.5*bminusrho2;break;

//   default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
//   }
// }

double computew3::f1(int i) const
{
  switch (i) {
  case  0: return  1.;break;
  case  1: return  alpha_b; break;
  case  2: return  2.*r_1-alpha_b;break;
  case  3: return  2.*r_0-1.;break;
  case  4: return 0.5-1.5*r_0;break;
  case  5: return 0.5*alpha_b-1.5*r_1;break;
  case  6: return -0.5*alpha_b-0.5*r_1;break;
  case  7: return -0.5*r_0-0.5;break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


void computew3::init(double _costheta0, double _costheta1, double _costheta2, double _P_b , double _alpha_b, double _alpha_lambda, double _r_0, double _r_1)
{
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;

  //speed-up
  // aplusrho2=a_plus_sq(alpha_b,r_0,r_1);
  // aminusrho2=a_minus_sq(alpha_b,r_0,r_1);
  // bplusrho2=b_plus_sq(alpha_b,r_0,r_1);
  // bminusrho2=b_minus_sq(alpha_b,r_0,r_1);

  //   aplusconj=TComplex::Conjugate(aplus);
  //   aminusconj=TComplex::Conjugate(aminus);
  //   bplusconj=TComplex::Conjugate(bplus);
  //   bminusconj=TComplex::Conjugate(bminus);
  
  //   aplusrho2=aplus.Rho2();
  //   aminusrho2=aminus.Rho2();
  //   bplusrho2=bplus.Rho2();
  //   bminusrho2=bminus.Rho2();

  costheta0pow2=costheta0*costheta0;
  costheta1pow2=costheta1*costheta1;
  costheta2pow2=costheta2*costheta2;
  costheta2pow3=costheta2*costheta2*costheta2;

}


double computew3::w(double _costheta0, double _costheta1, double _costheta2,
		    double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1/*, double xivector[8]*/)
{
  //  std::cerr << __LINE__ << std::endl;
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_r_0,_r_1);
  //  std::cerr << __LINE__ << std::endl;
  if (checknan()) return NAN;
  //  std::cerr << __LINE__ << std::endl;
  double val=0.;
  for (int i=0; i<=7; i++) val+=f1(i)*f2(i)*F(i);
  return val/16./pi;
}



double computew3::w_pi_costheta0(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_r_0,_r_1);
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1(i)*f2(i)*F_pi_costheta0(i);
  return val/16./pi;
}

double computew3::w_pi_costheta1(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_r_0,_r_1);
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1(i)*f2(i)*F_pi_costheta1(i);
  return val/16./pi;
}

double computew3::w_pi_costheta2(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_r_0,_r_1);
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1(i)*f2(i)*F_pi_costheta2(i);
  return val/16./pi;
}

double computew3::w_pi_costheta0_costheta1(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_r_0,_r_1);
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1(i)*f2(i)*F_pi_costheta0_costheta1(i);
  return val/16./pi;
}

double computew3::w_pi_costheta0_costheta2(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_r_0,_r_1);
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1(i)*f2(i)*F_pi_costheta0_costheta2(i);
  return val/16./pi;
}

double computew3::w_pi_costheta1_costheta2(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_r_0,_r_1);
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1(i)*f2(i)*F_pi_costheta1_costheta2(i);
  return val/16./pi;
}

double computew3::w_pi_P_b(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_r_0,_r_1);
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1(i)*f2_pi_P_b(i)*F(i);
  return val/16./pi;
}

double computew3::w_pi_alpha_lambda(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_r_0,_r_1);
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1(i)*f2_pi_alpha_lambda(i)*F(i);
  return val/16./pi;
}


double computew3::w_primitive(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1/*, double xivector[8]*/)
{ 
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_r_0,_r_1);
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1(i)*f2(i)*F_primitive(i);
  return val/16./pi;
}

