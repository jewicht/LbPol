#include "computew5.h"
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

  if (argc!=7 and argc!=9) {
    std::cout << argv[0] << " alpha_b r_0 r_1 alpha_plus alpha_minus chi" << std::endl;
    std::cout << argv[0] << " bm(r th) ap(r th) am(r th) bp(r th)" << std::endl;
    return 1;
  }

  if (argc==7) {
    const double alpha_b=atof(argv[1]);
    const double r_0=atof(argv[2]);
    const double r_1=atof(argv[3]);
    
    const double alpha_plus=atof(argv[4]);
    const double alpha_minus=atof(argv[5]);
    const double chi=atof(argv[6]);
    
    computew5 w5;
    
    const TComplex aplus= w5.a_plus (alpha_b,r_0,r_1,alpha_plus,alpha_minus,chi);
    const TComplex aminus=w5.a_minus(alpha_b,r_0,r_1,alpha_plus,alpha_minus,chi);
    const TComplex bplus= w5.b_plus (alpha_b,r_0,r_1,alpha_plus,alpha_minus,chi);
    const TComplex bminus=w5.b_minus(alpha_b,r_0,r_1,alpha_plus,alpha_minus,chi);
    
    std::cout 
      << bminus.Rho() << " " << bminus.Theta() << " "
      << aplus.Rho()  << " " << aplus.Theta()  << " "
      << aminus.Rho() << " " << aminus.Theta() << " "
      << bplus.Rho()  << " " << bplus.Theta()  << std::endl;
  }

  if (argc==9) {
    TComplex bminus(atof(argv[1]),atof(argv[2]),true);
    TComplex aplus( atof(argv[3]),atof(argv[4]),true);
    TComplex aminus(atof(argv[5]),atof(argv[6]),true);
    TComplex bplus( atof(argv[7]),atof(argv[8]),true);
    
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
