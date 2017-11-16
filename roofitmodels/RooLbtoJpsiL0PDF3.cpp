#include "RooLbtoJpsiL0PDF3.h"

ClassImp(  RooLbtoJpsiL0PDF3)

//_____________________________________________________________________________
RooLbtoJpsiL0PDF3::RooLbtoJpsiL0PDF3(const char *name, const char *title
				     , RooAbsReal& _costheta0, RooAbsReal& _costheta1, RooAbsReal& _costheta2
				     , RooAbsReal& _P_b, RooAbsReal&  _alpha_b,  RooAbsReal& _alpha_lambda,  RooAbsReal& _r_0,  RooAbsReal& _r_1
				     //				     , RooAbsReal& _xi0, RooAbsReal& _xi1, RooAbsReal& _xi2, RooAbsReal& _xi3
				     //				     ; RooAbsReal& _xi4, RooAbsReal& _xi5, RooAbsReal& _xi6, RooAbsReal& _xi7
				     ) : RooAbsPdf(name,title)
				       , costheta0("costheta0","Costheta0",this,_costheta0)
				       , costheta1("costheta1","Costheta1",this,_costheta1)
				       , costheta2("costheta2","Costheta2",this,_costheta2)
				       , P_b("P_b","P_b",this,_P_b)
				       , alpha_b("alpha_b","alpha_b",this,_alpha_b)
				       , alpha_lambda("alpha_lambda","alpha_lambda",this,_alpha_lambda)
				       , r_0("r_0","r_0",this,_r_0)
				       , r_1("r_1","r_1",this,_r_1)
				       // , xi0("xi0","xi0",this,_xi0)
				       // , xi1("xi1","xi1",this,_xi1)
				       // , xi2("xi2","xi2",this,_xi2)
				       // , xi3("xi3","xi3",this,_xi3)
				       // , xi4("xi4","xi4",this,_xi4)
				       // , xi5("xi5","xi5",this,_xi5)
				       // , xi6("xi6","xi6",this,_xi6)
				       // , xi7("xi7","xi7",this,_xi7)
{
}

//_____________________________________________________________________________
RooLbtoJpsiL0PDF3::RooLbtoJpsiL0PDF3(const RooLbtoJpsiL0PDF3& other, const char* name) : 
  RooAbsPdf(other,name), 
  costheta0("costheta0",this,other.costheta0), 
  costheta1("costheta1",this,other.costheta1), 
  costheta2("costheta2",this,other.costheta2), 
  P_b("P_b",this,other.P_b),
  alpha_b("alpha_b",this,other.alpha_b),
  alpha_lambda("alpha_lambda",this,other.alpha_lambda),
  r_0("r_0",this,other.r_0),
  r_1("r_1",this,other.r_1)//,
  // xi0("xi0",this,other.xi0),
  // xi1("xi1",this,other.xi1),
  // xi2("xi2",this,other.xi2),
  // xi3("xi3",this,other.xi3),
  // xi4("xi4",this,other.xi4),
  // xi5("xi5",this,other.xi5),
  // xi6("xi6",this,other.xi6),
  // xi7("xi7",this,other.xi7)
{
}

//_____________________________________________________________________________
Double_t RooLbtoJpsiL0PDF3::evaluate() const
{
  // Double_t xivector[8];
  // xivector[0]=xi0;
  // xivector[1]=xi1;
  // xivector[2]=xi2;
  // xivector[3]=xi3;
  // xivector[4]=xi4;
  // xivector[5]=xi5;
  // xivector[6]=xi6;
  // xivector[7]=xi7;

  // return 
  //   1. + 
  //   alpha_b * P_b * costheta0 + 
  //   (2.*r_1-alpha_b) * alpha_lambda * costheta1 +
  //   (2.*r_0 - 1.) * P_b * alpha_lambda * costheta0 * costheta1 +
  //   (0.5-1.5*r_0) * 0.5 * (3.*costheta2*costheta2 - 1.) +
  //   (0.5 * alpha_b  - 1.5 * r_1) * P_b * 0.5 * (3.*costheta2*costheta2 - 1.) * costheta0 +
  //   (-0.5 * alpha_b  - 0.5 * r_1) * alpha_lambda * 0.5 * (3.*costheta2*costheta2 - 1.) * costheta1 +
  //   (-0.5 * r_0  - 0.5) * P_b * alpha_lambda * 0.5 * (3.*costheta2*costheta2 - 1.) * costheta0 * costheta1;
  

  computew3 w3;
  const Double_t tmp=w3.w(costheta0,costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);//,xivector);
  //  if (tmp<0.) return 0.000001;
  return tmp;
}


//_____________________________________________________________________________
Int_t RooLbtoJpsiL0PDF3::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
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
  if ( matchArgs(allVars,analVars,P_b)          ) return 4 ;
  if ( matchArgs(allVars,analVars,alpha_lambda) ) return 5 ;
  
  return 0;
}



//_____________________________________________________________________________
Double_t RooLbtoJpsiL0PDF3::analyticalIntegral(Int_t code, const char* rangeName) const
{
  computew3 w3;
  Double_t result=0.;
  Double_t val[5][2];

  // Double_t xivector[8];
  // xivector[0]=xi0;
  // xivector[1]=xi1;
  // xivector[2]=xi2;
  // xivector[3]=xi3;
  // xivector[4]=xi4;
  // xivector[5]=xi5;
  // xivector[6]=xi6;
  // xivector[7]=xi7;

  switch (code) {
  case 1: 
    return w3.w_pi_costheta0(costheta0.max(),costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1)
      -    w3.w_pi_costheta0(costheta0.min(),costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);
    break;
  case 2: 
    return w3.w_pi_costheta1(costheta0,costheta1.max(),costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1)
      -    w3.w_pi_costheta1(costheta0,costheta1.min(),costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);
    break;
  case 3:
    return w3.w_pi_costheta2(costheta0,costheta1,costheta2.max(),P_b,alpha_b,alpha_lambda,r_0,r_1)
      -    w3.w_pi_costheta2(costheta0,costheta1,costheta2.min(),P_b,alpha_b,alpha_lambda,r_0,r_1);
    break;
  case 4:
    return w3.w_pi_P_b(costheta0,costheta1,costheta2,P_b.max(),alpha_b,alpha_lambda,r_0,r_1)
      -    w3.w_pi_P_b(costheta0,costheta1,costheta2,P_b.min(),alpha_b,alpha_lambda,r_0,r_1);
    break;
  case 5:
    return w3.w_pi_alpha_lambda(costheta0,costheta1,costheta2,P_b,alpha_b,alpha_lambda.max(),r_0,r_1)
      -    w3.w_pi_alpha_lambda(costheta0,costheta1,costheta2,P_b,alpha_b,alpha_lambda.min(),r_0,r_1);
    break;

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
	  result+=pow(-1., a+b+c+1) * w3.w_primitive(val[0][a],val[1][b],val[2][c],P_b,alpha_b,alpha_lambda,r_0,r_1);
    	}
      }
    }
    
    return result;
    //      +    w3.wprimitive(costheta.max(),costheta1.max(),costheta2.max(),phi1.max(),phi2.max(),P_b,alpha_b,alpha_lambda,r_0,r_1);
    break;

    
  case 6://for costheta0
    val[0][0]=costheta1.min();
    val[0][1]=costheta1.max();
    val[1][0]=costheta2.min();
    val[1][1]=costheta2.max();
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
	result+=pow(-1., a+b) * w3.w_pi_costheta1_costheta2(costheta0,val[0][a],val[1][b],P_b,alpha_b,alpha_lambda,r_0,r_1);
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
	result+=pow(-1., a+b) * w3.w_pi_costheta0_costheta2(val[0][a],costheta1,val[1][b],P_b,alpha_b,alpha_lambda,r_0,r_1);
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
	result+=pow(-1., a+b) * w3.w_pi_costheta0_costheta1(val[0][a],val[1][b],costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);
      }
    }
    return result;
    break;


  default:
    return 0.;
  }
}
