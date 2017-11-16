
#include <iostream>
#include "TComplex.h"

class computew3 
{
public:
  computew3();
  ~computew3(){};
  double w(double, double, double, double, double, double, double, double);

  double F(int) const;
  double f1(int) const; 
  double f2(int) const; 

  double w_pi_costheta0(double, double, double, double, double, double, double, double);
  double w_pi_costheta1(double, double, double, double, double, double, double, double);
  double w_pi_costheta2(double, double, double, double, double, double, double, double);

  double w_pi_costheta0_costheta1(double, double, double, double, double, double, double, double);
  double w_pi_costheta0_costheta2(double, double, double, double, double, double, double, double);
  double w_pi_costheta1_costheta2(double, double, double, double, double, double, double, double);

  double w_pi_alpha_lambda(double, double, double, double, double, double, double, double);
  double w_pi_P_b(double, double, double, double, double, double, double, double);

  double w_primitive(double, double, double, double, double, double, double, double);

  double a_plus_sq(double, double, double) const;
  double a_minus_sq(double, double, double) const;
  double b_plus_sq(double, double, double) const;
  double b_minus_sq(double, double, double) const;

  void init(double, double, double, double, double, double, double, double);

private:
  
  double alpha_b, alpha_lambda, P_b, r_0, r_1;
  
  static constexpr double beta_plus=0.;
  static constexpr double pi=3.141592654;
  
  double costheta0, costheta1, costheta2;

  bool checknan();
  
  double F_primitive(int) const; 

  double F_pi_costheta0(int) const;
  double F_pi_costheta1(int) const;
  double F_pi_costheta2(int) const;

  double F_pi_costheta0_costheta1(int) const;
  double F_pi_costheta0_costheta2(int) const;
  double F_pi_costheta1_costheta2(int) const;

  double f2_pi_P_b(int) const; 
  double f2_pi_alpha_lambda(int) const; 

  //  TComplex aplus,aminus,bplus,bminus;
  //  TComplex aplusconj,aminusconj,bplusconj,bminusconj;

  /* double aplusrho2, aminusrho2; */
  /* double bplusrho2, bminusrho2; */

  double costheta0pow2, costheta1pow2, costheta2pow2,costheta2pow3;
  //  double threedivsqrt2;

};
