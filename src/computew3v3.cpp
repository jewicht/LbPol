#include "computew3v3.h"


computew3v3::computew3v3()
{
}

bool computew3v3::checknan()
{
  return false;
//   if (r_0+r_1<0.) return true;
//   if (r_0-r_1<0.) return true;
//   if (alpha_b+1.-r_0-r_1<0.) return true;
//   if (1.-alpha_b-r_0+r_1<0.) return true;
//   return false;
}



double computew3v3::a_plus_sq(double _ap2, double _am2, double _bp2) const
{
  return _ap2;
}

double computew3v3::a_minus_sq(double _ap2, double _am2, double _bp2) const
{
  return _am2;
}

double computew3v3::b_plus_sq(double _ap2, double _am2, double _bp2) const
{
  return _bp2;
}

double computew3v3::b_minus_sq(double _ap2, double _am2, double _bp2) const
{
  return 1.-_ap2-_am2-_bp2;
}

double computew3v3::F(int i) const
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


double computew3v3::F_pi_costheta0(int i) const
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

double computew3v3::F_pi_costheta1(int i) const
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


double computew3v3::F_pi_costheta2(int i) const
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

double computew3v3::F_pi_costheta0_costheta1(int i) const
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

double computew3v3::F_pi_costheta0_costheta2(int i) const
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


double computew3v3::F_pi_costheta1_costheta2(int i) const
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



double computew3v3::F_primitive(int i) const
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


void computew3v3::initf2()
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


void computew3v3::initf1()
{
  f1[0]=1.;
  f1[1]=+ap2 - am2 + bp2 - (1.-ap2-am2-bp2);
  f1[2]=+ap2 - am2 - bp2 + (1.-ap2-am2-bp2);
  f1[3]=+ap2 + am2 - bp2 - (1.-ap2-am2-bp2);
  f1[4]=-ap2 - am2 + 0.5*bp2 + 0.5*(1.-ap2-am2-bp2);
  f1[5]=-ap2 + am2 + 0.5*bp2 - 0.5*(1.-ap2-am2-bp2);
  f1[6]=-ap2 + am2 - 0.5*bp2 + 0.5*(1.-ap2-am2-bp2);
  f1[7]=-ap2 - am2 - 0.5*bp2 - 0.5*(1.-ap2-am2-bp2);
}


void computew3v3::init(double _costheta0, double _costheta1, double _costheta2, double _P_b , double _ap2, double _alpha_lambda, double _am2, double _bp2)
{
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;

  P_b          =_P_b;
  ap2      =_ap2;
  alpha_lambda =_alpha_lambda;
  am2          =_am2;
  bp2         =_bp2;

  //speed-up
  // aplusrho2=a_plus_sq(ap2,r_0,r_1);
  // aminusrho2=a_minus_sq(ap2,r_0,r_1);
  // bplusrho2=b_plus_sq(ap2,r_0,r_1);
  // bminusrho2=b_minus_sq(ap2,r_0,r_1);

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


double computew3v3::w(double _costheta0, double _costheta1, double _costheta2,
		    double _P_b, double _ap2, double _alpha_lambda, double _am2, double _bp2/*, double xivector[8]*/)
{
  //  std::cerr << __LINE__ << std::endl;
  init(_costheta0,_costheta1,_costheta2,_P_b,_ap2,_alpha_lambda,_am2,_bp2);
  initf1();
  initf2();
  //  std::cerr << __LINE__ << std::endl;
  if (checknan()) return NAN;
  //  std::cerr << __LINE__ << std::endl;
  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F(i);
  return val/16./pi;
}



double computew3v3::w_pi_costheta0(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _ap2, double _alpha_lambda, double _am2, double _bp2/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_ap2,_alpha_lambda,_am2,_bp2);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_pi_costheta0(i);
  return val/16./pi;
}

double computew3v3::w_pi_costheta1(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _ap2, double _alpha_lambda, double _am2, double _bp2/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_ap2,_alpha_lambda,_am2,_bp2);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_pi_costheta1(i);
  return val/16./pi;
}

double computew3v3::w_pi_costheta2(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _ap2, double _alpha_lambda, double _am2, double _bp2/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_ap2,_alpha_lambda,_am2,_bp2);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_pi_costheta2(i);
  return val/16./pi;
}

double computew3v3::w_pi_costheta0_costheta1(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _ap2, double _alpha_lambda, double _am2, double _bp2/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_ap2,_alpha_lambda,_am2,_bp2);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_pi_costheta0_costheta1(i);
  return val/16./pi;
}

double computew3v3::w_pi_costheta0_costheta2(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _ap2, double _alpha_lambda, double _am2, double _bp2/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_ap2,_alpha_lambda,_am2,_bp2);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_pi_costheta0_costheta2(i);
  return val/16./pi;
}

double computew3v3::w_pi_costheta1_costheta2(double _costheta0, double _costheta1, double _costheta2,
		   double _P_b, double _ap2, double _alpha_lambda, double _am2, double _bp2/*, double xivector[8]*/)
{
  init(_costheta0,_costheta1,_costheta2,_P_b,_ap2,_alpha_lambda,_am2,_bp2);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_pi_costheta1_costheta2(i);
  return val/16./pi;
}

double computew3v3::w_primitive(double _costheta0, double _costheta1, double _costheta2,
				double _P_b, double _ap2, double _alpha_lambda, double _am2, double _bp2/*, double xivector[8]*/)
{ 
  init(_costheta0,_costheta1,_costheta2,_P_b,_ap2,_alpha_lambda,_am2,_bp2);
  initf1();
  initf2();
  if (checknan()) return NAN;

  double val=0.;
  for (int i=0; i<=7; i++) val+=f1[i]*f2[i]*F_primitive(i);
  return val/16./pi;
}

