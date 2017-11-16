#include "RooLbtoJpsiL0PDF3v3.h"

ClassImp(  RooLbtoJpsiL0PDF3v3)

//computew3v3 w3v3;

//_____________________________________________________________________________
RooLbtoJpsiL0PDF3v3::RooLbtoJpsiL0PDF3v3(const char *name, const char *title
					 , RooAbsReal& _costheta0, RooAbsReal& _costheta1, RooAbsReal& _costheta2
					 , RooAbsReal& _P_b, RooAbsReal&  _ap2,  RooAbsReal& _alpha_lambda,  RooAbsReal& _am2,  RooAbsReal& _bp2
					 ) : RooAbsPdf(name,title)
				       , costheta0("costheta0","Costheta0",this,_costheta0)
				       , costheta1("costheta1","Costheta1",this,_costheta1)
				       , costheta2("costheta2","Costheta2",this,_costheta2)
				       , P_b("P_b","P_b",this,_P_b)
				       , ap2("ap2","ap2",this,_ap2)
				       , alpha_lambda("alpha_lambda","alpha_lambda",this,_alpha_lambda)
				       , am2("am2","am2",this,_am2)
				       , bp2("bp2","bp2",this,_bp2)
{
}

//_____________________________________________________________________________
RooLbtoJpsiL0PDF3v3::RooLbtoJpsiL0PDF3v3(const RooLbtoJpsiL0PDF3v3& other, const char* name) : 
  RooAbsPdf(other,name), 
  costheta0("costheta0",this,other.costheta0), 
  costheta1("costheta1",this,other.costheta1), 
  costheta2("costheta2",this,other.costheta2), 
  P_b("P_b",this,other.P_b),
  ap2("ap2",this,other.ap2),
  alpha_lambda("alpha_lambda",this,other.alpha_lambda),
  am2("am2",this,other.am2),
  bp2("bp2",this,other.bp2)//,
{
}

//_____________________________________________________________________________
Double_t RooLbtoJpsiL0PDF3v3::evaluate() const
{
  // return 
  //   1. + 
  //   ap2 * P_b * costheta0 + 
  //   am2 * alpha_lambda * costheta1 +
  //   -1./3.* (1.+4.*bp2) * P_b * alpha_lambda * costheta0 * costheta1 +
  //   bp2 * 0.5 * (3.*costheta2*costheta2 - 1.) +
  //   -1./4. * (ap2+3.*am2) * P_b * 0.5 * (3.*costheta2*costheta2 - 1.) * costheta0 +
  //   -1./4. * (3.*ap2+am2) * alpha_lambda * 0.5 * (3.*costheta2*costheta2 - 1.) * costheta1 +
  //   1./3. * (bp2 - 2.) * P_b * alpha_lambda * 0.5 * (3.*costheta2*costheta2 - 1.) * costheta0 * costheta1;
  
  

  computew3v3 w3v3;
  const Double_t tmp=w3v3.w(costheta0,costheta1,costheta2,P_b,ap2,alpha_lambda,am2,bp2);//,xivector);
  //  if (tmp<0.) return 0.000001;
  return tmp;
}


//_____________________________________________________________________________
Int_t RooLbtoJpsiL0PDF3v3::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  // //3
  // if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(),costheta1.arg(),costheta2.arg()) ) ) return 9;

  // //2
  // if ( matchArgs(allVars,analVars, RooArgSet(costheta1.arg(), costheta2.arg()))) return  6;
  // if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(), costheta2.arg()))) return  7;
  // if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(), costheta1.arg()))) return  8;

  // //1
  // if ( matchArgs(allVars,analVars,costheta0)    ) return 1 ;
  // if ( matchArgs(allVars,analVars,costheta1)    ) return 2 ;
  // if ( matchArgs(allVars,analVars,costheta2)    ) return 3 ;
  // //  if ( matchArgs(allVars,analVars,P_b)          ) return 4 ;
  // //  if ( matchArgs(allVars,analVars,alpha_lambda) ) return 5 ;
  
  return 0;
}



//_____________________________________________________________________________
Double_t RooLbtoJpsiL0PDF3v3::analyticalIntegral(Int_t code, const char* rangeName) const
{
  Double_t result=0.;
  Double_t val[5][2];
  computew3v3 w3v3;
  switch (code) {
  case 1: 
    return w3v3.w_pi_costheta0(costheta0.max(),costheta1,costheta2,P_b,ap2,alpha_lambda,am2,bp2)
      -    w3v3.w_pi_costheta0(costheta0.min(),costheta1,costheta2,P_b,ap2,alpha_lambda,am2,bp2);
    break;
  case 2: 
    return w3v3.w_pi_costheta1(costheta0,costheta1.max(),costheta2,P_b,ap2,alpha_lambda,am2,bp2)
      -    w3v3.w_pi_costheta1(costheta0,costheta1.min(),costheta2,P_b,ap2,alpha_lambda,am2,bp2);
    break;
  case 3:
    return w3v3.w_pi_costheta2(costheta0,costheta1,costheta2.max(),P_b,ap2,alpha_lambda,am2,bp2)
      -    w3v3.w_pi_costheta2(costheta0,costheta1,costheta2.min(),P_b,ap2,alpha_lambda,am2,bp2);
    break;
  // case 4:
  //   return w3v3.w_pi_P_b(costheta0,costheta1,costheta2,P_b.max(),ap2,alpha_lambda,am2,bp2)
  //     -    w3v3.w_pi_P_b(costheta0,costheta1,costheta2,P_b.min(),ap2,alpha_lambda,am2,bp2);
  //   break;
  // case 5:
  //   return w3v3.w_pi_alpha_lambda(costheta0,costheta1,costheta2,P_b,ap2,alpha_lambda.max(),am2,bp2)
  //     -    w3v3.w_pi_alpha_lambda(costheta0,costheta1,costheta2,P_b,ap2,alpha_lambda.min(),am2,bp2);
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
  	  result+=pow(-1., a+b+c+1) * w3v3.w_primitive(val[0][a],val[1][b],val[2][c],P_b,ap2,alpha_lambda,am2,bp2);
    	}
      }
    }
    
    return result;
    //      +    w3v3.wprimitive(costheta.max(),costheta1.max(),costheta2.max(),phi1.max(),phi2.max(),P_b,ap2,alpha_lambda,am2,bp2);
    break;

    
  case 6://for costheta0
    val[0][0]=costheta1.min();
    val[0][1]=costheta1.max();
    val[1][0]=costheta2.min();
    val[1][1]=costheta2.max();
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
  	result+=pow(-1., a+b) * w3v3.w_pi_costheta1_costheta2(costheta0,val[0][a],val[1][b],P_b,ap2,alpha_lambda,am2,bp2);
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
  	result+=pow(-1., a+b) * w3v3.w_pi_costheta0_costheta2(val[0][a],costheta1,val[1][b],P_b,ap2,alpha_lambda,am2,bp2);
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
  	result+=pow(-1., a+b) * w3v3.w_pi_costheta0_costheta1(val[0][a],val[1][b],costheta2,P_b,ap2,alpha_lambda,am2,bp2);
      }
    }
    return result;
    break;


  default:
    return 0.;
  }
}
