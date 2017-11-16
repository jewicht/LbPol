#ifndef COMPUTEW5H
#define COMPUTEW5H

#include <iostream>
#include "TComplex.h"

class computew5 
{
public:
  computew5();
  ~computew5(){};
  double w(double, double, double, double, double, double, double, double, double, double, double, double, double);

  double w_pi_costheta0(double, double, double, double, double, double, double, double, double, double, double, double, double);
  double w_pi_costheta1(double, double, double, double, double, double, double, double, double, double, double, double, double);
  double w_pi_costheta2(double, double, double, double, double, double, double, double, double, double, double, double, double);
  double w_pi_phi1(double, double, double, double, double, double, double, double, double, double, double, double, double);
  double w_pi_phi2(double, double, double, double, double, double, double, double, double, double, double, double, double);
  double w_pi_phi2_phi1(double, double, double, double, double, double, double, double, double, double, double, double, double);

  double w_pi_alpha_lambda(double, double, double, double, double, double, double, double, double, double, double, double, double); 
  double w_pi_P_b(double, double, double, double, double, double, double, double, double, double, double, double, double); 

  double w_primitive(double, double, double, double, double, double, double, double, double, double, double, double, double); 

  double w_primitive_dcostheta0(double, double, double, double, double, double, double, double, double, double, double, double, double);
  double w_primitive_dcostheta1(double, double, double, double, double, double, double, double, double, double, double, double, double);
  double w_primitive_dcostheta2(double, double, double, double, double, double, double, double, double, double, double, double, double);
  double w_primitive_dphi1(double, double, double, double, double, double, double, double, double, double, double, double, double);
  double w_primitive_dphi2(double, double, double, double, double, double, double, double, double, double, double, double, double);

  TComplex a_plus(double, double, double, double, double, double) const;
  TComplex a_minus(double, double, double, double, double, double) const;
  TComplex b_plus(double, double, double, double, double, double) const;
  TComplex b_minus(double, double, double, double, double, double) const;
  
  void setparams(double, double, double, double, double, double, double, double);
  double f1(int) const; 
  double f2(int) const; 

private:
  
  double alpha_b, alpha_lambda, P_b, r_0, r_1, alpha_plus, alpha_minus, chi;
  
  double mynan;

  static constexpr double beta_plus=0.;
  static constexpr double pi=3.141592654;
  double pi2;
  
  double costheta0, costheta1, costheta2, phi1, phi2;

  bool checknan();
  
  //  void init();
  
  double F(int) const;
  double F_primitive(int) const; 
  double F_primitive_dcostheta0(int) const;
  double F_primitive_dcostheta1(int) const;
  double F_primitive_dcostheta2(int) const;
  double F_primitive_dphi1(int) const;
  double F_primitive_dphi2(int) const;

  double F_pi_costheta0(int) const;
  double F_pi_costheta1(int) const;
  double F_pi_costheta2(int) const;
  double F_pi_phi1(int) const;
  double F_pi_phi2(int) const;
  double F_pi_phi2_phi1(int) const;

  double f2_pi_P_b(int) const; 
  double f2_pi_alpha_lambda(int) const; 

  void speedup();
  
  TComplex aplus,aminus,bplus,bminus;
  TComplex aplusconj,aminusconj,bplusconj,bminusconj;

  TComplex apxam, bmxbp, bmxap, bmxam, amxbp, apxbp;

  double aplusrho2, aminusrho2;
  double bplusrho2, bminusrho2;

  double sintheta0,sintheta1,sintheta2;

  /* double cosphi1, sinphi1; */
  /* double cosphi2, sinphi2; */
  /* double cosphi1pphi2, sinphi1pphi2; */
  /* double cosphi1p2phi2, sinphi1p2phi2; */
  
  double costheta0pow2, costheta1pow2, costheta2pow2,costheta2pow3;
  double sintheta2pow2;
  double threedivsqrt2;

  double asincostheta0, asincostheta1, asincostheta2;
  
  double  pow1p5_of_1mcostheta2pow2;
};

#endif
