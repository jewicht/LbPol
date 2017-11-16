#include "computew3v2.h"


computew3v2::computew3v2()
{
}

bool computew3v2::checknan()
{
  //  std::cout << "|a+|^2 = " << a_plus_sq(alpha_b,beta,gamma) << std::endl;
  //  std::cout << "|a-|^2 = " << a_minus_sq(alpha_b,beta,gamma) << std::endl;
  //  std::cout << "|b+|^2 = " << b_plus_sq(alpha_b,beta,gamma) << std::endl;
  //  std::cout << "|b-|^2 = " << b_minus_sq(alpha_b,beta,gamma) << std::endl;

  // if (fabs(a_plus_sq(alpha_b,beta,gamma)-0.5)>0.5) return true;
  // if (fabs(a_minus_sq(alpha_b,beta,gamma)-0.5)>0.5) return true;
  // if (fabs(b_plus_sq(alpha_b,beta,gamma)-0.5)>0.5) return true;
  // if (fabs(b_minus_sq(alpha_b,beta,gamma)-0.5)>0.5) return true;
  
  return false;

}



double computew3v2::a_plus_sq(double _alpha_b, double _beta, double _gamma) const
{
  return 1./12.*(+3.*_alpha_b + 3.*_beta - 4.* _gamma + 2.  );
}

double computew3v2::a_minus_sq(double _alpha_b, double _beta, double _gamma) const
{
  return 1./12.*(-3.*_alpha_b - 3.*_beta - 4.* _gamma + 2.  );
}

double computew3v2::b_plus_sq(double _alpha_b, double _beta, double _gamma) const
{
  return 1./12.*(+3.*_alpha_b - 3.*_beta + 4.* _gamma + 4.  );
}

double computew3v2::b_minus_sq(double _alpha_b, double _beta, double _gamma) const
{
  return 1./12.*(-3.*_alpha_b + 3.*_beta + 4.* _gamma + 4.  );
}

double computew3v2::F(int i) const
{
  switch (i) {
    
  case  0 : return  1. ;break;
  case  1 : return  costheta0 ;break;
  case  2 : return  costheta1 ;break;
  case  3 : return  costheta0*costheta1 ;break;
  case  4 : return  0.5*(3.*costheta2pow2 - 1.);break;
  case  5 : return  0.5*(3.*costheta2pow2 - 1.)*costheta0 ;break;
  case  6 : return  0.5*(3.*costheta2pow2 - 1.)*costheta1 ;break;
  case  7 : return  0.5*(3.*costheta2pow2 - 1.)*costheta0*costheta1 ;break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


double computew3v2::F_pi_costheta0(int i) const
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

double computew3v2::F_pi_costheta1(int i) const
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


double computew3v2::F_pi_costheta2(int i) const
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

double computew3v2::F_pi_costheta0_costheta1(int i) const
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

double computew3v2::F_pi_costheta0_costheta2(int i) const
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


double computew3v2::F_pi_costheta1_costheta2(int i) const
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



double computew3v2::F_primitive(int i) const
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


void computew3v2::initf2()
{
  f2[0]=1.;
  f2[1]=P_b;
  f2[2]=alpha_lambda;
  f2[3]=P_b*alpha_lambda ;
  f2[4]=1. ;
  f2[5]=P_b ;
  f2[6]=alpha_lambda ;
  f2[7]=P_b*alpha_lambda;
}


void computew3v2::initf1()
{
  f1[0]=1.;
  f1[1]=alpha_b;
  f1[2]=beta;
  f1[3]=-1./3.*(1.+4*gamma);
  f1[4]=gamma;
  f1[5]=-0.25*(alpha_b+3.*beta);
  f1[6]=-0.25*(3.*alpha_b+beta);
  f1[7]=1./3.*(gamma-2.);
}


void computew3v2::init(double _costheta0, double _costheta1, double _costheta2, double _P_b , double _alpha_b, double _alpha_lambda, double _beta, double _gamma)
{
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  beta          =_beta;
  gamma         =_gamma;

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


double computew3v2::w(double _costheta0, double _costheta1, double _costheta2,
		    double _P_b, double _alpha_b, double _alpha_lambda, double _beta, double _gamma/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_beta,_gamma);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F(i);
  return val/16./pi;
}



double computew3v2::w_pi_costheta0(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _beta, double _gamma/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_beta,_gamma);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_pi_costheta0(i);
  return val/16./pi;
}

double computew3v2::w_pi_costheta1(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _beta, double _gamma/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_beta,_gamma);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_pi_costheta1(i);
  return val/16./pi;
}

double computew3v2::w_pi_costheta2(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _beta, double _gamma/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_beta,_gamma);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_pi_costheta2(i);
  return val/16./pi;
}

double computew3v2::w_pi_costheta0_costheta1(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _beta, double _gamma/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_beta,_gamma);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_pi_costheta0_costheta1(i);
  return val/16./pi;
}

double computew3v2::w_pi_costheta0_costheta2(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _beta, double _gamma/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_beta,_gamma);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_pi_costheta0_costheta2(i);
  return val/16./pi;
}

double computew3v2::w_pi_costheta1_costheta2(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _beta, double _gamma/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_beta,_gamma);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_pi_costheta1_costheta2(i);
  return val/16./pi;
}

double computew3v2::w_primitive(double _costheta0, double _costheta1, double _costheta2,
				double _P_b, double _alpha_b, double _alpha_lambda, double _beta, double _gamma/*, double xivector[8]*/)
{ 
  init(_costheta0,_costheta1,_costheta2,_P_b,_alpha_b,_alpha_lambda,_beta,_gamma);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_primitive(i);
  return val/16./pi;
}

