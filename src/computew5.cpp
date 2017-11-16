#include "computew5.h"


computew5::computew5()
{
  threedivsqrt2=3./sqrt(2.);
  mynan=sqrt(-1.);
  pi2=pi*pi;
}

bool computew5::checknan()
{
  //  if (r_0+r_1<0.) return true;
  //  if (r_0-r_1<0.) return true;
  //  if (alpha_b+1.-r_0-r_1<0.) return true;
  //  if (1.-alpha_b-r_0+r_1<0.) return true;
  return false;
}

inline TComplex mysqrt(double x)
{
  if (x>0.) return sqrt(x);
  return -sqrt(-x);
  //  return TComplex::Sqrt(x);
}

TComplex computew5::a_plus(double _alpha_b, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi) const
{
  return mysqrt(0.5*(_r_0+_r_1)) * TComplex::Exp(TComplex::I() * _alpha_plus);
}

TComplex computew5::a_minus(double _alpha_b, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi) const
{
  return mysqrt(0.5*(_r_0-_r_1)) * TComplex::Exp(TComplex::I() * _alpha_minus);
}

TComplex computew5::b_plus(double _alpha_b, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi) const
{
  return mysqrt(0.5*(_alpha_b+1.-_r_0-_r_1)) * TComplex::Exp(TComplex::I() * beta_plus);
}

TComplex computew5::b_minus(double _alpha_b, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi) const
{
  return mysqrt(0.5*(1.-_alpha_b-_r_0+_r_1)) * TComplex::Exp(TComplex::I() * (_alpha_plus + _alpha_minus - beta_plus - _chi));
}
//beta_minus = alphaplus + alphaminus - betaplus - chi
//chi = alphaplus + alphaminus - betaplus - betaminus
//alphaplus+alphaminus+betaplus+betaminus=2pi


double computew5::f1(int i) const
{
  switch (i) {
  case  0: return  aplusrho2 + aminusrho2 + bplusrho2 + bminusrho2;break;
  case  1: return  aplusrho2 - aminusrho2 + bplusrho2 - bminusrho2;break;
  case  2: return  aplusrho2 - aminusrho2 - bplusrho2 + bminusrho2;break;
  case  3: return  aplusrho2 + aminusrho2 - bplusrho2 - bminusrho2;break;

  case  4: return -aplusrho2 - aminusrho2 + 0.5*bplusrho2 + 0.5*bminusrho2;break;
  case  5: return -aplusrho2 + aminusrho2 + 0.5*bplusrho2 - 0.5*bminusrho2;break;
  case  6: return -aplusrho2 + aminusrho2 - 0.5*bplusrho2 + 0.5*bminusrho2;break;
  case  7: return -aplusrho2 - aminusrho2 - 0.5*bplusrho2 - 0.5*bminusrho2;break;

  case  8: return -3.  * (apxam).Re(); break;
  case  9: return  3.  * (apxam).Im(); break;
  case 10: return -1.5 * (bmxbp).Re();  break;
  case 11: return  1.5 * (bmxbp).Im();  break;

  case 12: return -threedivsqrt2 * (bmxap + amxbp).Re(); break;
  case 13: return  threedivsqrt2 * (bmxap + amxbp).Im(); break;
  case 14: return -threedivsqrt2 * (bmxam + apxbp).Re(); break;
  case 15: return  threedivsqrt2 * (bmxam + apxbp).Im(); break;
    
  case 16: return  threedivsqrt2 * (amxbp - bmxap).Re(); break;
  case 17: return -threedivsqrt2 * (amxbp - bmxap).Im(); break;
  case 18: return  threedivsqrt2 * (bmxam - apxbp).Re(); break;
  case 19: return -threedivsqrt2 * (bmxam - apxbp).Im(); break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}



double computew5::f2(int i) const
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
  case  8: return P_b*alpha_lambda; break;
  case  9: return P_b*alpha_lambda; break;
  case 10: return P_b*alpha_lambda; break;
  case 11: return P_b*alpha_lambda; break;
  case 12: return P_b*alpha_lambda; break;
  case 13: return P_b*alpha_lambda; break;
  case 14: return P_b*alpha_lambda; break;
  case 15: return P_b*alpha_lambda; break;
  case 16: return P_b; break;
  case 17: return P_b; break;
  case 18: return alpha_lambda; break;
  case 19: return alpha_lambda; break; 
  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


double computew5::F(int i) const
{
  switch (i) {

case  0 : return  1. ;break;
case  1 : return  costheta0 ;break;
case  2 : return  costheta1 ;break;
case  3 : return  costheta0*costheta1 ;break;
case  4 : return  3./2.*costheta2pow2 - 1./2. ;break;
case  5 : return  1./2.*(3.*costheta2pow2 - 1.)*costheta0 ;break;
case  6 : return  1./2.*(3.*costheta2pow2 - 1.)*costheta1 ;break;
case  7 : return  1./2.*(3.*costheta2pow2 - 1.)*costheta0*costheta1 ;break;
case  8 : return  -(costheta2pow2 - 1.)*sintheta1*sintheta0*cos(pi*phi1) ;break;
case  9 : return  -(costheta2pow2 - 1.)*sintheta1*sintheta0*sin(pi*phi1) ;break;
case  10 : return  -(costheta2pow2 - 1.)*sintheta1*sintheta0*cos(pi*phi1 + 2.*pi*phi2) ;break;
case  11 : return  -(costheta2pow2 - 1.)*sintheta1*sintheta0*sin(pi*phi1 + 2.*pi*phi2) ;break;
case  12 : return  sintheta2*sintheta0*costheta1*costheta2*cos(pi*phi2) ;break;
case  13 : return  sintheta2*sintheta0*costheta1*costheta2*sin(pi*phi2) ;break;
case  14 : return  sintheta2*sintheta1*costheta0*costheta2*cos(pi*phi1 + pi*phi2) ;break;
case  15 : return  sintheta2*sintheta1*costheta0*costheta2*sin(pi*phi1 + pi*phi2) ;break;
case  16 : return  sintheta2*sintheta0*costheta2*cos(pi*phi2) ;break;
case  17 : return  sintheta2*sintheta0*costheta2*sin(pi*phi2) ;break;
case  18 : return  sintheta2*sintheta1*costheta2*cos(pi*phi1 + pi*phi2) ;break;
case  19 : return  sintheta2*sintheta1*costheta2*sin(pi*phi1 + pi*phi2) ;break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


double computew5::F_pi_costheta0(int i) const
{
  switch (i) {

case  0 : return  costheta0 ;break;
case  1 : return  1./2.*costheta0pow2 ;break;
case  2 : return  costheta0*costheta1 ;break;
case  3 : return  1./2.*costheta0pow2*costheta1 ;break;
case  4 : return  1./2.*(3*costheta2pow2 - 1)*costheta0 ;break;
case  5 : return  1./4.*(3*costheta2pow2 - 1)*costheta0pow2 ;break;
case  6 : return  1./2.*(3*costheta2pow2 - 1)*costheta0*costheta1 ;break;
case  7 : return  1./4.*(3*costheta2pow2 - 1)*costheta0pow2*costheta1 ;break;
case  8 : return  -1./2.*(costheta2pow2 - 1.)*sintheta1*(sintheta0*costheta0 + asin(costheta0))*cos(pi*phi1) ;break;
case  9 : return  -1./2.*(costheta2pow2 - 1.)*sintheta1*(sintheta0*costheta0 + asin(costheta0))*sin(pi*phi1) ;break;
case  10 : return  -1./2.*(costheta2pow2 - 1.)*sintheta1*(sintheta0*costheta0 + asin(costheta0))*cos(pi*phi1 + 2.*pi*phi2) ;break;
case  11 : return  -1./2.*(costheta2pow2 - 1.)*sintheta1*(sintheta0*costheta0 + asin(costheta0))*sin(pi*phi1 + 2.*pi*phi2) ;break;
case  12 : return  1./2.*sintheta2*(sintheta0*costheta0 + asin(costheta0))*costheta1*costheta2*cos(pi*phi2) ;break;
case  13 : return  1./2.*sintheta2*(sintheta0*costheta0 + asin(costheta0))*costheta1*costheta2*sin(pi*phi2) ;break;
case  14 : return  1./2.*sintheta2*sintheta1*costheta0pow2*costheta2*cos(pi*phi1 + pi*phi2) ;break;
case  15 : return  1./2.*sintheta2*sintheta1*costheta0pow2*costheta2*sin(pi*phi1 + pi*phi2) ;break;
case  16 : return  1./2.*sintheta2*(sintheta0*costheta0 + asin(costheta0))*costheta2*cos(pi*phi2) ;break;
case  17 : return  1./2.*sintheta2*(sintheta0*costheta0 + asin(costheta0))*costheta2*sin(pi*phi2) ;break;
case  18 : return  sintheta2*sintheta1*costheta0*costheta2*cos(pi*phi1 + pi*phi2) ;break;
case  19 : return  sintheta2*sintheta1*costheta0*costheta2*sin(pi*phi1 + pi*phi2) ;break;

    
  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
    
  }
}

double computew5::F_pi_costheta1(int i) const
{
  switch (i) {

case  0 : return  costheta1 ;break;
case  1 : return  costheta0*costheta1 ;break;
case  2 : return  1./2.*costheta1pow2 ;break;
case  3 : return  1./2.*costheta0*costheta1pow2 ;break;
case  4 : return  1./2.*(3*costheta2pow2 - 1)*costheta1 ;break;
case  5 : return  1./2.*(3*costheta2pow2 - 1)*costheta0*costheta1 ;break;
case  6 : return  1./4.*(3*costheta2pow2 - 1)*costheta1pow2 ;break;
case  7 : return  1./4.*(3*costheta2pow2 - 1)*costheta0*costheta1pow2 ;break;
case  8 : return  -1./2.*(costheta2pow2 - 1.)*sintheta0*(sintheta1*costheta1 + asin(costheta1))*cos(pi*phi1) ;break;
case  9 : return  -1./2.*(costheta2pow2 - 1.)*sintheta0*(sintheta1*costheta1 + asin(costheta1))*sin(pi*phi1) ;break;
case  10 : return  -1./2.*(costheta2pow2 - 1.)*sintheta0*(sintheta1*costheta1 + asin(costheta1))*cos(pi*phi1 + 2.*pi*phi2) ;break;
case  11 : return  -1./2.*(costheta2pow2 - 1.)*sintheta0*(sintheta1*costheta1 + asin(costheta1))*sin(pi*phi1 + 2.*pi*phi2) ;break;
case  12 : return  1./2.*sintheta2*sintheta0*costheta1pow2*costheta2*cos(pi*phi2) ;break;
case  13 : return  1./2.*sintheta2*sintheta0*costheta1pow2*costheta2*sin(pi*phi2) ;break;
case  14 : return  1./2.*sintheta2*(sintheta1*costheta1 + asin(costheta1))*costheta0*costheta2*cos(pi*phi1 + pi*phi2) ;break;
case  15 : return  1./2.*sintheta2*(sintheta1*costheta1 + asin(costheta1))*costheta0*costheta2*sin(pi*phi1 + pi*phi2) ;break;
case  16 : return  sintheta2*sintheta0*costheta1*costheta2*cos(pi*phi2) ;break;
case  17 : return  sintheta2*sintheta0*costheta1*costheta2*sin(pi*phi2) ;break;
case  18 : return  1./2.*sintheta2*(sintheta1*costheta1 + asin(costheta1))*costheta2*cos(pi*phi1 + pi*phi2) ;break;
case  19 : return  1./2.*sintheta2*(sintheta1*costheta1 + asin(costheta1))*costheta2*sin(pi*phi1 + pi*phi2) ;break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);

  }
}


double computew5::F_pi_costheta2(int i) const
{
  switch (i) {

case  0 : return  costheta2 ;break;
case  1 : return  costheta0*costheta2 ;break;
case  2 : return  costheta1*costheta2 ;break;
case  3 : return  costheta0*costheta1*costheta2 ;break;
case  4 : return  1./2.*costheta2pow3 - 1./2.*costheta2 ;break;
case  5 : return  1./2.*(costheta2pow3 - costheta2)*costheta0 ;break;
case  6 : return  1./2.*(costheta2pow3 - costheta2)*costheta1 ;break;
case  7 : return  1./2.*(costheta2pow3 - costheta2)*costheta0*costheta1 ;break;
case  8 : return  -1./3.*sintheta1*sintheta0*(costheta2pow3 - 3.*costheta2)*cos(pi*phi1) ;break;
case  9 : return  -1./3.*sintheta1*sintheta0*(costheta2pow3 - 3.*costheta2)*sin(pi*phi1) ;break;
case  10 : return  -1./3.*sintheta1*sintheta0*(costheta2pow3 - 3.*costheta2)*cos(pi*phi1 + 2.*pi*phi2) ;break;
case  11 : return  -1./3.*sintheta1*sintheta0*(costheta2pow3 - 3.*costheta2)*sin(pi*phi1 + 2.*pi*phi2) ;break;
case  12 : return  -1./3.*pow1p5_of_1mcostheta2pow2*sintheta0*costheta1*cos(pi*phi2) ;break;
case  13 : return  -1./3.*pow1p5_of_1mcostheta2pow2*sintheta0*costheta1*sin(pi*phi2) ;break;
case  14 : return  -1./3.*pow1p5_of_1mcostheta2pow2*sintheta1*costheta0*cos(pi*phi1 + pi*phi2) ;break;
case  15 : return  -1./3.*pow1p5_of_1mcostheta2pow2*sintheta1*costheta0*sin(pi*phi1 + pi*phi2) ;break;
case  16 : return  -1./3.*pow1p5_of_1mcostheta2pow2*sintheta0*cos(pi*phi2) ;break;
case  17 : return  -1./3.*pow1p5_of_1mcostheta2pow2*sintheta0*sin(pi*phi2) ;break;
case  18 : return  -1./3.*pow1p5_of_1mcostheta2pow2*sintheta1*cos(pi*phi1 + pi*phi2) ;break;
case  19 : return  -1./3.*pow1p5_of_1mcostheta2pow2*sintheta1*sin(pi*phi1 + pi*phi2) ;break;


  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


double computew5::F_pi_phi1(int i) const
{
  switch (i) {

case  0 : return  phi1 ;break;
case  1 : return  costheta0*phi1 ;break;
case  2 : return  costheta1*phi1 ;break;
case  3 : return  costheta0*costheta1*phi1 ;break;
case  4 : return  1./2.*(3*costheta2pow2 - 1)*phi1 ;break;
case  5 : return  1./2.*(3*costheta2pow2 - 1)*costheta0*phi1 ;break;
case  6 : return  1./2.*(3*costheta2pow2 - 1)*costheta1*phi1 ;break;
case  7 : return  1./2.*(3*costheta2pow2 - 1)*costheta0*costheta1*phi1 ;break;
case  8 : return  -(costheta2pow2 - 1.)*sintheta1*sintheta0*sin(pi*phi1)/pi ;break;
case  9 : return  (costheta2pow2 - 1.)*sintheta1*sintheta0*cos(pi*phi1)/pi ;break;
case  10 : return  -(costheta2pow2 - 1.)*sintheta1*sintheta0*sin(pi*phi1 + 2.*pi*phi2)/pi ;break;
case  11 : return  (costheta2pow2 - 1.)*sintheta1*sintheta0*cos(pi*phi1 + 2.*pi*phi2)/pi ;break;
case  12 : return  sintheta2*sintheta0*costheta1*costheta2*phi1*cos(pi*phi2) ;break;
case  13 : return  sintheta2*sintheta0*costheta1*costheta2*phi1*sin(pi*phi2) ;break;
case  14 : return  sintheta2*sintheta1*costheta0*costheta2*sin(pi*phi1 + pi*phi2)/pi ;break;
case  15 : return  -sintheta2*sintheta1*costheta0*costheta2*cos(pi*phi1 + pi*phi2)/pi ;break;
case  16 : return  sintheta2*sintheta0*costheta2*phi1*cos(pi*phi2) ;break;
case  17 : return  sintheta2*sintheta0*costheta2*phi1*sin(pi*phi2) ;break;
case  18 : return  sintheta2*sintheta1*costheta2*sin(pi*phi1 + pi*phi2)/pi ;break;
case  19 : return  -sintheta2*sintheta1*costheta2*cos(pi*phi1 + pi*phi2)/pi ;break;


  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}

double computew5::F_pi_phi2(int i) const
{
  switch (i) {

case  0 : return  phi2 ;break;
case  1 : return  costheta0*phi2 ;break;
case  2 : return  costheta1*phi2 ;break;
case  3 : return  costheta0*costheta1*phi2 ;break;
case  4 : return  1./2.*(3*costheta2pow2 - 1)*phi2 ;break;
case  5 : return  1./2.*(3*costheta2pow2 - 1)*costheta0*phi2 ;break;
case  6 : return  1./2.*(3*costheta2pow2 - 1)*costheta1*phi2 ;break;
case  7 : return  1./2.*(3*costheta2pow2 - 1)*costheta0*costheta1*phi2 ;break;
case  8 : return  -(costheta2pow2 - 1.)*sintheta1*sintheta0*phi2*cos(pi*phi1) ;break;
case  9 : return  -(costheta2pow2 - 1.)*sintheta1*sintheta0*phi2*sin(pi*phi1) ;break;
case  10 : return  -1./2.*(costheta2pow2 - 1.)*sintheta1*sintheta0*sin(pi*phi1 + 2.*pi*phi2)/pi ;break;
case  11 : return  1./2.*(costheta2pow2 - 1.)*sintheta1*sintheta0*cos(pi*phi1 + 2.*pi*phi2)/pi ;break;
case  12 : return  sintheta2*sintheta0*costheta1*costheta2*sin(pi*phi2)/pi ;break;
case  13 : return  -sintheta2*sintheta0*costheta1*costheta2*cos(pi*phi2)/pi ;break;
case  14 : return  sintheta2*sintheta1*costheta0*costheta2*sin(pi*phi1 + pi*phi2)/pi ;break;
case  15 : return  -sintheta2*sintheta1*costheta0*costheta2*cos(pi*phi1 + pi*phi2)/pi ;break;
case  16 : return  sintheta2*sintheta0*costheta2*sin(pi*phi2)/pi ;break;
case  17 : return  -sintheta2*sintheta0*costheta2*cos(pi*phi2)/pi ;break;
case  18 : return  sintheta2*sintheta1*costheta2*sin(pi*phi1 + pi*phi2)/pi ;break;
case  19 : return  -sintheta2*sintheta1*costheta2*cos(pi*phi1 + pi*phi2)/pi ;break;


  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}

double computew5::F_pi_phi2_phi1(int i) const
{
  switch (i) {

case  0 : return  phi1*phi2 ;break;
case  1 : return  costheta0*phi1*phi2 ;break;
case  2 : return  costheta1*phi1*phi2 ;break;
case  3 : return  costheta0*costheta1*phi1*phi2 ;break;
case  4 : return  1./2.*(3*costheta2pow2 - 1)*phi1*phi2 ;break;
case  5 : return  1./2.*(3*costheta2pow2 - 1)*costheta0*phi1*phi2 ;break;
case  6 : return  1./2.*(3*costheta2pow2 - 1)*costheta1*phi1*phi2 ;break;
case  7 : return  1./2.*(3*costheta2pow2 - 1)*costheta0*costheta1*phi1*phi2 ;break;
case  8 : return  -(costheta2pow2 - 1.)*sintheta1*sintheta0*phi2*sin(pi*phi1)/pi ;break;
case  9 : return  (costheta2pow2 - 1.)*sintheta1*sintheta0*phi2*cos(pi*phi1)/pi ;break;
case  10 : return  1./2.*(costheta2pow2 - 1.)*sintheta1*sintheta0*cos(pi*phi1 + 2.*pi*phi2)/pi2 ;break;
case  11 : return  1./2.*(costheta2pow2 - 1.)*sintheta1*sintheta0*sin(pi*phi1 + 2.*pi*phi2)/pi2 ;break;
case  12 : return  sintheta2*sintheta0*costheta1*costheta2*phi1*sin(pi*phi2)/pi ;break;
case  13 : return  -sintheta2*sintheta0*costheta1*costheta2*phi1*cos(pi*phi2)/pi ;break;
case  14 : return  -sintheta2*sintheta1*costheta0*costheta2*cos(pi*phi1 + pi*phi2)/pi2 ;break;
case  15 : return  -sintheta2*sintheta1*costheta0*costheta2*sin(pi*phi1 + pi*phi2)/pi2 ;break;
case  16 : return  sintheta2*sintheta0*costheta2*phi1*sin(pi*phi2)/pi ;break;
case  17 : return  -sintheta2*sintheta0*costheta2*phi1*cos(pi*phi2)/pi ;break;
case  18 : return  -sintheta2*sintheta1*costheta2*cos(pi*phi1 + pi*phi2)/pi2 ;break;
case  19 : return  -sintheta2*sintheta1*costheta2*sin(pi*phi1 + pi*phi2)/pi2 ;break;


  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


double computew5::F_primitive(int i) const
{
  switch (i) {

case  0 : return  costheta0*costheta1*costheta2*phi1*phi2 ;break;
case  1 : return  1./2.*costheta0pow2*costheta1*costheta2*phi1*phi2 ;break;
case  2 : return  1./2.*costheta0*costheta1pow2*costheta2*phi1*phi2 ;break;
case  3 : return  1./4.*costheta0pow2*costheta1pow2*costheta2*phi1*phi2 ;break;
case  4 : return  1./2.*(costheta2pow3 - costheta2)*costheta0*costheta1*phi1*phi2 ;break;
case  5 : return  1./4.*(costheta2pow3 - costheta2)*costheta0pow2*costheta1*phi1*phi2 ;break;
case  6 : return  1./4.*(costheta2pow3 - costheta2)*costheta0*costheta1pow2*phi1*phi2 ;break;
case  7 : return  1./8.*(costheta2pow3 - costheta2)*costheta0pow2*costheta1pow2*phi1*phi2 ;break;
case  8 : return  -1./12.*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*phi2*sin(pi*phi1)/pi ;break;
case  9 : return  1./12.*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*phi2*cos(pi*phi1)/pi ;break;
case  10 : return  1./24.*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*cos(pi*phi1 + 2.*pi*phi2)/pi2 ;break;
case  11 : return  1./24.*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*sin(pi*phi1 + 2.*pi*phi2)/pi2 ;break;
case  12 : return  -1./12.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1pow2*phi1*sin(pi*phi2)/pi ;break;
case  13 : return  1./12.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1pow2*phi1*cos(pi*phi2)/pi ;break;
case  14 : return  1./12.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0pow2*cos(pi*phi1 + pi*phi2)/pi2 ;break;
case  15 : return  1./12.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0pow2*sin(pi*phi1 + pi*phi2)/pi2 ;break;
case  16 : return  -1./6.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1*phi1*sin(pi*phi2)/pi ;break;
case  17 : return  1./6.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1*phi1*cos(pi*phi2)/pi ;break;
case  18 : return  1./6.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0*cos(pi*phi1 + pi*phi2)/pi2 ;break;
case  19 : return  1./6.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0*sin(pi*phi1 + pi*phi2)/pi2 ;break;


  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}

double computew5::F_primitive_dcostheta0(int i) const
{
  switch (i) {

case  0 : return  costheta1*costheta2*phi1*phi2 ;break;
case  1 : return  costheta0*costheta1*costheta2*phi1*phi2 ;break;
case  2 : return  1./2.*costheta1pow2*costheta2*phi1*phi2 ;break;
case  3 : return  1./2.*costheta0*costheta1pow2*costheta2*phi1*phi2 ;break;
case  4 : return  1./2.*(costheta2pow3 - costheta2)*costheta1*phi1*phi2 ;break;
case  5 : return  1./2.*(costheta2pow3 - costheta2)*costheta0*costheta1*phi1*phi2 ;break;
case  6 : return  1./4.*(costheta2pow3 - costheta2)*costheta1pow2*phi1*phi2 ;break;
case  7 : return  1./4.*(costheta2pow3 - costheta2)*costheta0*costheta1pow2*phi1*phi2 ;break;
case  8 : return  -1./6.*sintheta0*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*phi2*sin(pi*phi1)/pi ;break;
case  9 : return  1./6.*sintheta0*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*phi2*cos(pi*phi1)/pi ;break;
case  10 : return  1./12.*sintheta0*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*cos(pi*phi1 + 2.*pi*phi2)/pi2 ;break;
case  11 : return  1./12.*sintheta0*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*sin(pi*phi1 + 2.*pi*phi2)/pi2 ;break;
case  12 : return  -1./6.*pow1p5_of_1mcostheta2pow2*sintheta0*costheta1pow2*phi1*sin(pi*phi2)/pi ;break;
case  13 : return  1./6.*pow1p5_of_1mcostheta2pow2*sintheta0*costheta1pow2*phi1*cos(pi*phi2)/pi ;break;
case  14 : return  1./6.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0*cos(pi*phi1 + pi*phi2)/pi2 ;break;
case  15 : return  1./6.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0*sin(pi*phi1 + pi*phi2)/pi2 ;break;
case  16 : return  -1./3.*pow1p5_of_1mcostheta2pow2*sintheta0*costheta1*phi1*sin(pi*phi2)/pi ;break;
case  17 : return  1./3.*pow1p5_of_1mcostheta2pow2*sintheta0*costheta1*phi1*cos(pi*phi2)/pi ;break;
case  18 : return  1./6.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*cos(pi*phi1 + pi*phi2)/pi2 ;break;
case  19 : return  1./6.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*sin(pi*phi1 + pi*phi2)/pi2 ;break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}

double computew5::F_primitive_dcostheta1(int i) const
{
  switch (i) {
case  0 : return  costheta0*costheta2*phi1*phi2 ;break;
case  1 : return  1./2.*costheta0pow2*costheta2*phi1*phi2 ;break;
case  2 : return  costheta0*costheta1*costheta2*phi1*phi2 ;break;
case  3 : return  1./2.*costheta0pow2*costheta1*costheta2*phi1*phi2 ;break;
case  4 : return  1./2.*(costheta2pow3 - costheta2)*costheta0*phi1*phi2 ;break;
case  5 : return  1./4.*(costheta2pow3 - costheta2)*costheta0pow2*phi1*phi2 ;break;
case  6 : return  1./2.*(costheta2pow3 - costheta2)*costheta0*costheta1*phi1*phi2 ;break;
case  7 : return  1./4.*(costheta2pow3 - costheta2)*costheta0pow2*costheta1*phi1*phi2 ;break;
case  8 : return  -1./6.*sintheta1*(sintheta0*costheta0 + asin(costheta0))*(costheta2pow3 - 3.*costheta2)*phi2*sin(pi*phi1)/pi ;break;
case  9 : return  1./6.*sintheta1*(sintheta0*costheta0 + asin(costheta0))*(costheta2pow3 - 3.*costheta2)*phi2*cos(pi*phi1)/pi ;break;
case  10 : return  1./12.*sintheta1*(sintheta0*costheta0 + asin(costheta0))*(costheta2pow3 - 3.*costheta2)*cos(pi*phi1 + 2.*pi*phi2)/pi2 ;break;
case  11 : return  1./12.*sintheta1*(sintheta0*costheta0 + asin(costheta0))*(costheta2pow3 - 3.*costheta2)*sin(pi*phi1 + 2.*pi*phi2)/pi2 ;break;
case  12 : return  -1./6.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1*phi1*sin(pi*phi2)/pi ;break;
case  13 : return  1./6.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1*phi1*cos(pi*phi2)/pi ;break;
case  14 : return  1./6.*pow1p5_of_1mcostheta2pow2*sintheta1*costheta0pow2*cos(pi*phi1 + pi*phi2)/pi2 ;break;
case  15 : return  1./6.*pow1p5_of_1mcostheta2pow2*sintheta1*costheta0pow2*sin(pi*phi1 + pi*phi2)/pi2 ;break;
case  16 : return  -1./6.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*phi1*sin(pi*phi2)/pi ;break;
case  17 : return  1./6.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*phi1*cos(pi*phi2)/pi ;break;
case  18 : return  1./3.*pow1p5_of_1mcostheta2pow2*sintheta1*costheta0*cos(pi*phi1 + pi*phi2)/pi2 ;break;
case  19 : return  1./3.*pow1p5_of_1mcostheta2pow2*sintheta1*costheta0*sin(pi*phi1 + pi*phi2)/pi2 ;break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}

double computew5::F_primitive_dcostheta2(int i) const
{
  switch (i) {
case  0 : return  costheta0*costheta1*phi1*phi2 ;break;
case  1 : return  1./2.*costheta0pow2*costheta1*phi1*phi2 ;break;
case  2 : return  1./2.*costheta0*costheta1pow2*phi1*phi2 ;break;
case  3 : return  1./4.*costheta0pow2*costheta1pow2*phi1*phi2 ;break;
case  4 : return  1./2.*(3*costheta2pow2 - 1)*costheta0*costheta1*phi1*phi2 ;break;
case  5 : return  1./4.*(3*costheta2pow2 - 1)*costheta0pow2*costheta1*phi1*phi2 ;break;
case  6 : return  1./4.*(3*costheta2pow2 - 1)*costheta0*costheta1pow2*phi1*phi2 ;break;
case  7 : return  1./8.*(3*costheta2pow2 - 1)*costheta0pow2*costheta1pow2*phi1*phi2 ;break;
case  8 : return  -1./4.*(costheta2pow2 - 1.)*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*phi2*sin(pi*phi1)/pi ;break;
case  9 : return  1./4.*(costheta2pow2 - 1.)*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*phi2*cos(pi*phi1)/pi ;break;
case  10 : return  1./8.*(costheta2pow2 - 1.)*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*cos(pi*phi1 + 2.*pi*phi2)/pi2 ;break;
case  11 : return  1./8.*(costheta2pow2 - 1.)*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*sin(pi*phi1 + 2.*pi*phi2)/pi2 ;break;
case  12 : return  1./4.*sintheta2*(sintheta0*costheta0 + asin(costheta0))*costheta1pow2*costheta2*phi1*sin(pi*phi2)/pi ;break;
case  13 : return  -1./4.*sintheta2*(sintheta0*costheta0 + asin(costheta0))*costheta1pow2*costheta2*phi1*cos(pi*phi2)/pi ;break;
case  14 : return  -1./4.*sintheta2*(sintheta1*costheta1 + asin(costheta1))*costheta0pow2*costheta2*cos(pi*phi1 + pi*phi2)/pi2 ;break;
case  15 : return  -1./4.*sintheta2*(sintheta1*costheta1 + asin(costheta1))*costheta0pow2*costheta2*sin(pi*phi1 + pi*phi2)/pi2 ;break;
case  16 : return  1./2.*sintheta2*(sintheta0*costheta0 + asin(costheta0))*costheta1*costheta2*phi1*sin(pi*phi2)/pi ;break;
case  17 : return  -1./2.*sintheta2*(sintheta0*costheta0 + asin(costheta0))*costheta1*costheta2*phi1*cos(pi*phi2)/pi ;break;
case  18 : return  -1./2.*sintheta2*(sintheta1*costheta1 + asin(costheta1))*costheta0*costheta2*cos(pi*phi1 + pi*phi2)/pi2 ;break;
case  19 : return  -1./2.*sintheta2*(sintheta1*costheta1 + asin(costheta1))*costheta0*costheta2*sin(pi*phi1 + pi*phi2)/pi2 ;break;

  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}

double computew5::F_primitive_dphi1(int i) const
{
  switch (i) {

case  0 : return  costheta0*costheta1*costheta2*phi2 ;break;
case  1 : return  1./2.*costheta0pow2*costheta1*costheta2*phi2 ;break;
case  2 : return  1./2.*costheta0*costheta1pow2*costheta2*phi2 ;break;
case  3 : return  1./4.*costheta0pow2*costheta1pow2*costheta2*phi2 ;break;
case  4 : return  1./2.*(costheta2pow3 - costheta2)*costheta0*costheta1*phi2 ;break;
case  5 : return  1./4.*(costheta2pow3 - costheta2)*costheta0pow2*costheta1*phi2 ;break;
case  6 : return  1./4.*(costheta2pow3 - costheta2)*costheta0*costheta1pow2*phi2 ;break;
case  7 : return  1./8.*(costheta2pow3 - costheta2)*costheta0pow2*costheta1pow2*phi2 ;break;
case  8 : return  -1./12.*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*phi2*cos(pi*phi1) ;break;
case  9 : return  -1./12.*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*phi2*sin(pi*phi1) ;break;
case  10 : return  -1./24.*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*sin(pi*phi1 + 2.*pi*phi2)/pi ;break;
case  11 : return  1./24.*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*cos(pi*phi1 + 2.*pi*phi2)/pi ;break;
case  12 : return  -1./12.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1pow2*sin(pi*phi2)/pi ;break;
case  13 : return  1./12.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1pow2*cos(pi*phi2)/pi ;break;
case  14 : return  -1./12.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0pow2*sin(pi*phi1 + pi*phi2)/pi ;break;
case  15 : return  1./12.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0pow2*cos(pi*phi1 + pi*phi2)/pi ;break;
case  16 : return  -1./6.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1*sin(pi*phi2)/pi ;break;
case  17 : return  1./6.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1*cos(pi*phi2)/pi ;break;
case  18 : return  -1./6.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0*sin(pi*phi1 + pi*phi2)/pi ;break;
case  19 : return  1./6.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0*cos(pi*phi1 + pi*phi2)/pi ;break;


  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}

double computew5::F_primitive_dphi2(int i) const
{
  switch (i) {
case  0 : return  costheta0*costheta1*costheta2*phi1 ;break;
case  1 : return  1./2.*costheta0pow2*costheta1*costheta2*phi1 ;break;
case  2 : return  1./2.*costheta0*costheta1pow2*costheta2*phi1 ;break;
case  3 : return  1./4.*costheta0pow2*costheta1pow2*costheta2*phi1 ;break;
case  4 : return  1./2.*(costheta2pow3 - costheta2)*costheta0*costheta1*phi1 ;break;
case  5 : return  1./4.*(costheta2pow3 - costheta2)*costheta0pow2*costheta1*phi1 ;break;
case  6 : return  1./4.*(costheta2pow3 - costheta2)*costheta0*costheta1pow2*phi1 ;break;
case  7 : return  1./8.*(costheta2pow3 - costheta2)*costheta0pow2*costheta1pow2*phi1 ;break;
case  8 : return  -1./12.*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*sin(pi*phi1)/pi ;break;
case  9 : return  1./12.*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*cos(pi*phi1)/pi ;break;
case  10 : return  -1./12.*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*sin(pi*phi1 + 2.*pi*phi2)/pi ;break;
case  11 : return  1./12.*(sintheta0*costheta0 + asin(costheta0))*(sintheta1*costheta1 + asin(costheta1))*(costheta2pow3 - 3.*costheta2)*cos(pi*phi1 + 2.*pi*phi2)/pi ;break;
case  12 : return  -1./12.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1pow2*phi1*cos(pi*phi2) ;break;
case  13 : return  -1./12.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1pow2*phi1*sin(pi*phi2) ;break;
case  14 : return  -1./12.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0pow2*sin(pi*phi1 + pi*phi2)/pi ;break;
case  15 : return  1./12.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0pow2*cos(pi*phi1 + pi*phi2)/pi ;break;
case  16 : return  -1./6.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1*phi1*cos(pi*phi2) ;break;
case  17 : return  -1./6.*pow1p5_of_1mcostheta2pow2*(sintheta0*costheta0 + asin(costheta0))*costheta1*phi1*sin(pi*phi2) ;break;
case  18 : return  -1./6.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0*sin(pi*phi1 + pi*phi2)/pi ;break;
case  19 : return  1./6.*pow1p5_of_1mcostheta2pow2*(sintheta1*costheta1 + asin(costheta1))*costheta0*cos(pi*phi1 + pi*phi2)/pi ;break;


  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}

double computew5::f2_pi_P_b(int i) const
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
  case  8: return 0.5*pow(P_b,2)*alpha_lambda; break;
  case  9: return 0.5*pow(P_b,2)*alpha_lambda; break;
  case 10: return 0.5*pow(P_b,2)*alpha_lambda; break;
  case 11: return 0.5*pow(P_b,2)*alpha_lambda; break;
  case 12: return 0.5*pow(P_b,2)*alpha_lambda; break;
  case 13: return 0.5*pow(P_b,2)*alpha_lambda; break;
  case 14: return 0.5*pow(P_b,2)*alpha_lambda; break;
  case 15: return 0.5*pow(P_b,2)*alpha_lambda; break;
  case 16: return 0.5*pow(P_b,2); break;
  case 17: return 0.5*pow(P_b,2); break;
  case 18: return alpha_lambda*P_b; break;
  case 19: return alpha_lambda*P_b; break; 
  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}


double computew5::f2_pi_alpha_lambda(int i) const
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
  case  8: return P_b*0.5*pow(alpha_lambda,2); break;
  case  9: return P_b*0.5*pow(alpha_lambda,2); break;
  case 10: return P_b*0.5*pow(alpha_lambda,2); break;
  case 11: return P_b*0.5*pow(alpha_lambda,2); break;
  case 12: return P_b*0.5*pow(alpha_lambda,2); break;
  case 13: return P_b*0.5*pow(alpha_lambda,2); break;
  case 14: return P_b*0.5*pow(alpha_lambda,2); break;
  case 15: return P_b*0.5*pow(alpha_lambda,2); break;
  case 16: return P_b*alpha_lambda; break;
  case 17: return P_b*alpha_lambda; break;
  case 18: return 0.5*pow(alpha_lambda,2); break;
  case 19: return 0.5*pow(alpha_lambda,2); break; 
  default: std::cerr << "bug " << __LINE__ << std::endl;exit(1);
  }
}



void computew5::speedup()
{
  //speed-up

  aplusrho2=0.5*(r_0+r_1);
  aminusrho2=0.5*(r_0-r_1);
  bplusrho2=0.5*(alpha_b+1.-r_0-r_1);
  bminusrho2=0.5*(1.-alpha_b-r_0+r_1);

  aplus=  a_plus( alpha_b, r_0, r_1, alpha_plus, alpha_minus, chi);
  aminus= a_minus(alpha_b, r_0, r_1, alpha_plus, alpha_minus, chi);
  bplus=  b_plus( alpha_b, r_0, r_1, alpha_plus, alpha_minus, chi);
  bminus= b_minus(alpha_b, r_0, r_1, alpha_plus, alpha_minus, chi);

  aplusconj=TComplex::Conjugate(aplus);
  aminusconj=TComplex::Conjugate(aminus);
  bplusconj=TComplex::Conjugate(bplus);
  bminusconj=TComplex::Conjugate(bminus);

  
  
  //  aplusrho2 =aplus.Rho2();
  //  aminusrho2=aminus.Rho2();
  //  bplusrho2 =bplus.Rho2();
  //  bminusrho2=bminus.Rho2();

  apxam=aplus  * aminusconj;
  amxbp=aminus * bplusconj;
  apxbp=aplus  * bplusconj;

  bmxbp=bminus * bplusconj;
  bmxap=bminus * aplusconj;
  bmxam=bminus * aminusconj;


  sintheta0=sqrt(1.-costheta0*costheta0);
  sintheta1=sqrt(1.-costheta1*costheta1);
  sintheta2=sqrt(1.-costheta2*costheta2);

  // cosphi1=cos(phi1);
  // sinphi1=sin(phi1);
  // cosphi2=cos(phi2);
  // sinphi2=sin(phi2);

  // cosphi1pphi2=cos(phi1+phi2);
  // sinphi1pphi2=sin(phi1+phi2);

  // cosphi1p2phi2=cos(phi1+2.*phi2);
  // sinphi1p2phi2=sin(phi1+2.*phi2);

  costheta0pow2=costheta0*costheta0;
  costheta1pow2=costheta1*costheta1;
  costheta2pow2=costheta2*costheta2;
  sintheta2pow2=sintheta2*sintheta2;
  costheta2pow3=costheta2*costheta2*costheta2;

  asincostheta0=asin(costheta0);
  asincostheta1=asin(costheta1);
  asincostheta2=asin(costheta2);

  pow1p5_of_1mcostheta2pow2 = pow(-costheta2pow2 + 1., 1.5);
  
}

void computew5::setparams(double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{
  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  //  if (checknan()) return mynan;
  speedup();
}

double computew5::w(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F(i);
  return val/pow(4.*pi,3);
}



double computew5::w_pi_costheta0(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F_pi_costheta0(i);
  return val/pow(4.*pi,3);
}

double computew5::w_pi_costheta1(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F_pi_costheta1(i);
  return val/pow(4.*pi,3);
}

double computew5::w_pi_costheta2(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F_pi_costheta2(i);
  return val/pow(4.*pi,3);
}

double computew5::w_pi_phi1(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F_pi_phi1(i);
  return val/pow(4.*pi,3);
}

double computew5::w_pi_phi2(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  speedup();
  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F_pi_phi2(i);
  return val/pow(4.*pi,3);
}

double computew5::w_pi_phi2_phi1(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  speedup();
  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F_pi_phi2_phi1(i);
  return val/pow(4.*pi,3);
}

double computew5::w_pi_P_b(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2_pi_P_b(i)*F(i);
  return val/pow(4.*pi,3);
}

double computew5::w_pi_alpha_lambda(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2_pi_alpha_lambda(i)*F(i);
  return val/pow(4.*pi,3);
}


double computew5::w_primitive(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{ 
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F_primitive(i);
  return val/pow(4.*pi,3);
}

double computew5::w_primitive_dcostheta0(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{ 
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;


  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F_primitive_dcostheta0(i);
  return val/pow(4.*pi,3);
}

double computew5::w_primitive_dcostheta1(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{ 
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F_primitive_dcostheta1(i);
  return val/pow(4.*pi,3);
}

double computew5::w_primitive_dcostheta2(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{ 
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F_primitive_dcostheta2(i);
  return val/pow(4.*pi,3);
}

double computew5::w_primitive_dphi1(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{ 
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;


  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F_primitive_dphi1(i);
  return val/pow(4.*pi,3);
}

double computew5::w_primitive_dphi2(double _costheta0, double _costheta1, double _costheta2, double _phi1, double _phi2,
		   double _P_b, double _alpha_b, double _alpha_lambda, double _r_0, double _r_1, double _alpha_plus, double _alpha_minus, double _chi)
{ 
  costheta0 =_costheta0;
  costheta1 =_costheta1;
  costheta2 =_costheta2;
  phi1   =_phi1;
  phi2   =_phi2;

  P_b          =_P_b;
  alpha_b      =_alpha_b;
  alpha_lambda =_alpha_lambda;
  r_0          =_r_0;
  r_1          =_r_1;
  alpha_plus   =_alpha_plus;
  alpha_minus  =_alpha_minus;
  chi          =_chi;

  if (checknan()) return mynan;
  speedup();

  double val=0.;
  for (int i=0; i<=19; i++) val+=f1(i)*f2(i)*F_primitive_dphi2(i);
  return val/pow(4.*pi,3);
}



