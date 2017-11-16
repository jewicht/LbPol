#include "RooLbtoJpsiL0PDF3v2.h"

ClassImp(  RooLbtoJpsiL0PDF3v2)


//_____________________________________________________________________________
RooLbtoJpsiL0PDF3v2::RooLbtoJpsiL0PDF3v2(const char *name, const char *title
					 , RooAbsReal& _costheta0, RooAbsReal& _costheta1, RooAbsReal& _costheta2
					 , RooAbsReal& _P_b, RooAbsReal&  _alpha_b,  RooAbsReal& _alpha_lambda,  RooAbsReal& _beta,  RooAbsReal& _gamma
					 ) : RooAbsPdf(name,title)
				       , costheta0("costheta0","Costheta0",this,_costheta0)
				       , costheta1("costheta1","Costheta1",this,_costheta1)
				       , costheta2("costheta2","Costheta2",this,_costheta2)
				       , P_b("P_b","P_b",this,_P_b)
				       , alpha_b("alpha_b","alpha_b",this,_alpha_b)
				       , alpha_lambda("alpha_lambda","alpha_lambda",this,_alpha_lambda)
				       , beta("beta","beta",this,_beta)
				       , gamma("gamma","gamma",this,_gamma)
{
}

//_____________________________________________________________________________
RooLbtoJpsiL0PDF3v2::RooLbtoJpsiL0PDF3v2(const RooLbtoJpsiL0PDF3v2& other, const char* name) : 
  RooAbsPdf(other,name), 
  costheta0("costheta0",this,other.costheta0), 
  costheta1("costheta1",this,other.costheta1), 
  costheta2("costheta2",this,other.costheta2), 
  P_b("P_b",this,other.P_b),
  alpha_b("alpha_b",this,other.alpha_b),
  alpha_lambda("alpha_lambda",this,other.alpha_lambda),
  beta("beta",this,other.beta),
  gamma("gamma",this,other.gamma)//,
{
}

//_____________________________________________________________________________
Double_t RooLbtoJpsiL0PDF3v2::evaluate() const
{
  // return 
  //   1. + 
  //   alpha_b * P_b * costheta0 + 
  //   beta * alpha_lambda * costheta1 +
  //   -1./3.* (1.+4.*gamma) * P_b * alpha_lambda * costheta0 * costheta1 +
  //   gamma * 0.5 * (3.*costheta2*costheta2 - 1.) +
  //   -1./4. * (alpha_b+3.*beta) * P_b * 0.5 * (3.*costheta2*costheta2 - 1.) * costheta0 +
  //   -1./4. * (3.*alpha_b+beta) * alpha_lambda * 0.5 * (3.*costheta2*costheta2 - 1.) * costheta1 +
  //   1./3. * (gamma - 2.) * P_b * alpha_lambda * 0.5 * (3.*costheta2*costheta2 - 1.) * costheta0 * costheta1;
  
  

  computew3v2 w3v2;
  const Double_t tmp=w3v2.w(costheta0,costheta1,costheta2,P_b,alpha_b,alpha_lambda,beta,gamma);//,xivector);
  //  if (tmp<0.) return 0.000001;
  return tmp;
}


//_____________________________________________________________________________
Int_t RooLbtoJpsiL0PDF3v2::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  //3
  if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(),costheta1.arg(),costheta2.arg()) ) ) return 9;

  //2
  if ( matchArgs(allVars,analVars, RooArgSet(costheta1.arg(), costheta2.arg()))) return  6;
  if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(), costheta2.arg()))) return  7;
  if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(), costheta1.arg()))) return  8;

  //1
  if ( matchArgs(allVars,analVars,costheta0)    ) return 1 ;
  if ( matchArgs(allVars,analVars,costheta1)    ) return 2 ;
  if ( matchArgs(allVars,analVars,costheta2)    ) return 3 ;
  //  if ( matchArgs(allVars,analVars,P_b)          ) return 4 ;
  //  if ( matchArgs(allVars,analVars,alpha_lambda) ) return 5 ;
  
  return 0;
}



//_____________________________________________________________________________
Double_t RooLbtoJpsiL0PDF3v2::analyticalIntegral(Int_t code, const char* rangeName) const
{
  computew3v2 w3v2;
  Double_t result=0.;
  Double_t val[5][2];

  switch (code) {
  case 1: 
    return w3v2.w_pi_costheta0(costheta0.max(),costheta1,costheta2,P_b,alpha_b,alpha_lambda,beta,gamma)
      -    w3v2.w_pi_costheta0(costheta0.min(),costheta1,costheta2,P_b,alpha_b,alpha_lambda,beta,gamma);
    break;
  case 2: 
    return w3v2.w_pi_costheta1(costheta0,costheta1.max(),costheta2,P_b,alpha_b,alpha_lambda,beta,gamma)
      -    w3v2.w_pi_costheta1(costheta0,costheta1.min(),costheta2,P_b,alpha_b,alpha_lambda,beta,gamma);
    break;
  case 3:
    return w3v2.w_pi_costheta2(costheta0,costheta1,costheta2.max(),P_b,alpha_b,alpha_lambda,beta,gamma)
      -    w3v2.w_pi_costheta2(costheta0,costheta1,costheta2.min(),P_b,alpha_b,alpha_lambda,beta,gamma);
    break;
  // case 4:
  //   return w3v2.w_pi_P_b(costheta0,costheta1,costheta2,P_b.max(),alpha_b,alpha_lambda,beta,gamma)
  //     -    w3v2.w_pi_P_b(costheta0,costheta1,costheta2,P_b.min(),alpha_b,alpha_lambda,beta,gamma);
  //   break;
  // case 5:
  //   return w3v2.w_pi_alpha_lambda(costheta0,costheta1,costheta2,P_b,alpha_b,alpha_lambda.max(),beta,gamma)
  //     -    w3v2.w_pi_alpha_lambda(costheta0,costheta1,costheta2,P_b,alpha_b,alpha_lambda.min(),beta,gamma);
  //   break;

  case 9:
    val[0][0]=costheta0.min();
    val[0][1]=costheta0.max();
    val[1][0]=costheta1.min();
    val[1][1]=costheta1.max();
    val[2][0]=costheta2.min();
    val[2][1]=costheta2.max();
    
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
    	for (int c=0; c<=1; c++) {
  	  result+=pow(-1., a+b+c+1) * w3v2.w_primitive(val[0][a],val[1][b],val[2][c],P_b,alpha_b,alpha_lambda,beta,gamma);
    	}
      }
    }
    
    return result;
    //      +    w3v2.wprimitive(costheta.max(),costheta1.max(),costheta2.max(),phi1.max(),phi2.max(),P_b,alpha_b,alpha_lambda,beta,gamma);
    break;

    
  case 6://for costheta0
    val[0][0]=costheta1.min();
    val[0][1]=costheta1.max();
    val[1][0]=costheta2.min();
    val[1][1]=costheta2.max();
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
  	result+=pow(-1., a+b) * w3v2.w_pi_costheta1_costheta2(costheta0,val[0][a],val[1][b],P_b,alpha_b,alpha_lambda,beta,gamma);
      }
    }
    return result;
    break;

  case 7://for costheta1
    val[0][0]=costheta0.min();
    val[0][1]=costheta0.max();
    val[1][0]=costheta2.min();
    val[1][1]=costheta2.max();

    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
  	result+=pow(-1., a+b) * w3v2.w_pi_costheta0_costheta2(val[0][a],costheta1,val[1][b],P_b,alpha_b,alpha_lambda,beta,gamma);
      }
    }
    return result;
    break;


  case 8://for costheta2
    val[0][0]=costheta0.min();
    val[0][1]=costheta0.max();
    val[1][0]=costheta1.min();
    val[1][1]=costheta1.max();
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
  	result+=pow(-1., a+b) * w3v2.w_pi_costheta0_costheta1(val[0][a],val[1][b],costheta2,P_b,alpha_b,alpha_lambda,beta,gamma);
      }
    }
    return result;
    break;


  default:
    return 0.;
  }
}
