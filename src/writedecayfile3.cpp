#include "computew3.h"


void test()
{

  computew3 w3;
  
  double costheta=0.846551;
  double costheta1=-0.566535;
  double costheta2=0.785396;
  double Pb=0.;
  double alphab=0.;
  double alpha_lambda=0.;
  double r0=0.981656;
  double r1=0.019779;

  std::cout << w3.w(costheta, costheta1, costheta2, Pb,alphab,alpha_lambda, r0, r1) << std::endl; 
  std::cout << "a+ = " << w3.a_plus_sq(alphab, r0, r1) << std::endl;
  std::cout << "a- = " << w3.a_minus_sq(alphab, r0, r1) << std::endl;
  std::cout << "b+ = " << w3.b_plus_sq(alphab, r0, r1) << std::endl;
  std::cout << "b- = " << w3.b_minus_sq(alphab, r0, r1) << std::endl;
}

int main(int argc, char** argv)
{

  //  test();
  //  return 1;

  if (argc!=4 and argc!=5) {
    std::cout << argv[0] << " alpha_b r_0 r_1" << std::endl;
    std::cout << argv[0] << " bm2 ap2 am2 bp2" << std::endl;
    return 1;
  }

  if (argc==4) {
    const double alpha_b=atof(argv[1]);
    const double r_0=atof(argv[2]);
    const double r_1=atof(argv[3]);
    
    computew3 w3;
    const double aplus_sq= w3.a_plus_sq (alpha_b,r_0,r_1);
    const double aminus_sq=w3.a_minus_sq(alpha_b,r_0,r_1);
    const double bplus_sq= w3.b_plus_sq (alpha_b,r_0,r_1);
    const double bminus_sq=w3.b_minus_sq(alpha_b,r_0,r_1);
    
    std::cout << "bminus_sq = " << bminus_sq  << std::endl;
    std::cout << "aplus_sq  = " << aplus_sq  << std::endl;
    std::cout << "aminus_sq = " << aminus_sq  << std::endl;
    std::cout << "bplus_sq  = " << bplus_sq  << std::endl;

    std::cout << "bminus = " << sqrt(bminus_sq)  << std::endl;
    std::cout << "aplus  = " << sqrt(aplus_sq)  << std::endl;
    std::cout << "aminus = " << sqrt(aminus_sq)  << std::endl;
    std::cout << "bplus  = " << sqrt(bplus_sq)  << std::endl;

  }

  if (argc==5) {
    TComplex bminus(sqrt(atof(argv[1])),0.,true);
    TComplex aplus( sqrt(atof(argv[2])),0.,true);
    TComplex aminus(sqrt(atof(argv[3])),0.,true);
    TComplex bplus( sqrt(atof(argv[4])),0.,true);
    
    //normalization
    const double norm = sqrt(aplus.Rho2() + aminus.Rho2() + bplus.Rho2() + bminus.Rho2());
    //    const double norm = aplus.Rho2() + aminus.Rho2() + bplus.Rho2() + bminus.Rho2();

    aplus/=norm;
    aminus/=norm;
    bplus/=norm;
    bminus/=norm;

    //rotation so that beta_plus==0
    const TComplex rotation = TComplex::Exp(- TComplex::I() * bplus.Theta());
    
    aplus*=rotation;
    aminus*=rotation;
    bplus*=rotation;
    bminus*=rotation;

    std::cout 
      << bminus.Rho() << " " << bminus.Theta() << " "
      << aplus.Rho()  << " " << aplus.Theta()  << " "
      << aminus.Rho() << " " << aminus.Theta() << " "
      << bplus.Rho()  << " " << bplus.Theta()  << std::endl;

    std::cout << "sum         = " << aplus.Rho2() + aminus.Rho2() + bplus.Rho2() + bminus.Rho2()  << std::endl;
    std::cout << "alpha_b     = " << aplus.Rho2() - aminus.Rho2() + bplus.Rho2() - bminus.Rho2()  << std::endl;
    std::cout << "r_0         = " << aplus.Rho2() + aminus.Rho2() << std::endl;
    std::cout << "r_1         = " << aplus.Rho2() - aminus.Rho2() << std::endl;
    
    std::cout << "alpha_plus  = " << aplus.Theta()  << std::endl;
    std::cout << "alpha_minus = " << aminus.Theta() << std::endl;
    std::cout << "beta_plus   = " << bplus.Theta()  << std::endl;
    std::cout << "beta_minus  = " << bminus.Theta() << std::endl;
    std::cout << "chi         = " << aplus.Theta() + aminus.Theta() - bplus.Theta() - bminus.Theta() << std::endl;
  }
}
