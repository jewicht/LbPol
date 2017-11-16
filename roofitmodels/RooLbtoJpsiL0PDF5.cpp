#include "RooLbtoJpsiL0PDF5.h"

ClassImp(  RooLbtoJpsiL0PDF5)

computew5 w2;

//_____________________________________________________________________________
RooLbtoJpsiL0PDF5::RooLbtoJpsiL0PDF5(const char *name, const char *title,
		   RooAbsReal& _costheta0, RooAbsReal& _costheta1, RooAbsReal& _costheta2, RooAbsReal& _phi1, RooAbsReal& _phi2, 
		   RooAbsReal& _P_b, RooAbsReal&  _alpha_b,  RooAbsReal& _alpha_lambda,  RooAbsReal& _r_0,  RooAbsReal& _r_1, RooAbsReal&  _alpha_plus,  RooAbsReal& _alpha_minus,  RooAbsReal& _chi) : RooAbsPdf(name,title),
  costheta0("costheta0","Costheta0",this,_costheta0),
  costheta1("costheta1","Costheta1",this,_costheta1),
  costheta2("costheta2","Costheta2",this,_costheta2),
  phi1("phi1","Phi1",this,_phi1),
  phi2("phi2","Phi2",this,_phi2),
  P_b("P_b","P_b",this,_P_b),
  alpha_b("alpha_b","alpha_b",this,_alpha_b),
  alpha_lambda("alpha_lambda","alpha_lambda",this,_alpha_lambda),
  r_0("r_0","r_0",this,_r_0),
  r_1("r_1","r_1",this,_r_1),
  alpha_plus("alpha_plus","alpha_plus",this,_alpha_plus), 
  alpha_minus("alpha_minus","alpha_minus",this,_alpha_minus),
  chi("chi","chi",this,_chi)
{
}

//_____________________________________________________________________________
RooLbtoJpsiL0PDF5::RooLbtoJpsiL0PDF5(const RooLbtoJpsiL0PDF5& other, const char* name) : 
  RooAbsPdf(other,name), 
  costheta0("costheta0",this,other.costheta0), 
  costheta1("costheta1",this,other.costheta1), 
  costheta2("costheta2",this,other.costheta2), 
  phi1("phi1",this,other.phi1),
  phi2("phi2",this,other.phi2),
  P_b("P_b",this,other.P_b),
  alpha_b("alpha_b",this,other.alpha_b),
  alpha_lambda("alpha_lambda",this,other.alpha_lambda),
  r_0("r_0",this,other.r_0),
  r_1("r_1",this,other.r_1),
  alpha_plus("alpha_plus",this,other.alpha_plus), 
  alpha_minus("alpha_minus",this,other.alpha_minus),
  chi("chi",this,other.chi)
{
}

//_____________________________________________________________________________
Double_t RooLbtoJpsiL0PDF5::evaluate() const
{
  return w2.w(costheta0,costheta1,costheta2,phi1,phi2,P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
}


//_____________________________________________________________________________
Int_t RooLbtoJpsiL0PDF5::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  //5
  if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(),costheta1.arg(),costheta2.arg(),phi1.arg(),phi2.arg()) ) ) return 8;

  //4
  if ( matchArgs(allVars,analVars, RooArgSet(costheta1.arg(), costheta2.arg(), phi1.arg(),   phi2.arg()))) return  9;
  if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(), costheta2.arg(), phi1.arg(),   phi2.arg()))) return 10;
  if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(), costheta1.arg(), phi1.arg(),   phi2.arg()))) return 11;
  if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(), costheta1.arg(), costheta2.arg(), phi2.arg()))) return 12;
  if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(), costheta1.arg(), costheta2.arg(), phi1.arg()))) return 13;

  //3
  if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(), costheta1.arg(), costheta2.arg()))) return 14;
  
  //1
  if ( matchArgs(allVars,analVars,costheta0)    ) return 1 ;
  if ( matchArgs(allVars,analVars,costheta1)    ) return 2 ;
  if ( matchArgs(allVars,analVars,costheta2)    ) return 3 ;
  if ( matchArgs(allVars,analVars,phi1)         ) return 4 ;
  if ( matchArgs(allVars,analVars,phi2)         ) return 5 ;
  if ( matchArgs(allVars,analVars,P_b)          ) return 6 ;
  if ( matchArgs(allVars,analVars,alpha_lambda) ) return 7 ;
  
  return 0;
}



//_____________________________________________________________________________
Double_t RooLbtoJpsiL0PDF5::analyticalIntegral(Int_t code, const char* rangeName) const
{

  // Double_t x1[5];
  // Double_t par1[8];
  // Double_t x0[5];
  // Double_t par0[8];
  // computew w;

  //  x0[0]=theta;
  // x0[1]=theta_1;
  // x0[2]=theta_2;
  // x0[3]=phi_1;
  // x0[4]=phi2;
  
  // par0[0]=P_b;
  // par0[1]=alpha_b;
  // par0[2]=alpha_lambda;
  // par0[3]=r_0;
  // par0[4]=r_1;
  // par0[5]=alpha_plus;
  // par0[6]=alpha_minus;
  // par0[7]=chi;
  
  // x1[0]=theta;
  // x1[1]=theta_1;
  // x1[2]=theta_2;
  // x1[3]=phi_1;
  // x1[4]=phi2;
  
  // par1[0]=P_b;
  // par1[1]=alpha_b;
  // par1[2]=alpha_lambda;
  // par1[3]=r_0;
  // par1[4]=r_1;
  // par1[5]=alpha_plus;
  // par1[6]=alpha_minus;
  // par1[7]=chi;
  
  //  computew w;

  Double_t result;
  result=0.;
  Double_t val[5][2];


  switch (code) {
  case 1: 
    return w2.w_pi_costheta0(costheta0.max(),costheta1,costheta2,phi1,phi2,P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi)
      -    w2.w_pi_costheta0(costheta0.min(),costheta1,costheta2,phi1,phi2,P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    break;
  case 2: 
    return w2.w_pi_costheta1(costheta0,costheta1.max(),costheta2,phi1,phi2,P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi)
      -    w2.w_pi_costheta1(costheta0,costheta1.min(),costheta2,phi1,phi2,P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    break;
  case 3:
    return w2.w_pi_costheta2(costheta0,costheta1,costheta2.max(),phi1,phi2,P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi)
      -    w2.w_pi_costheta2(costheta0,costheta1,costheta2.min(),phi1,phi2,P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    break;
  case 4:
    return w2.w_pi_phi1(costheta0,costheta1,costheta2,phi1.max(),phi2,P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi)
      -    w2.w_pi_phi1(costheta0,costheta1,costheta2,phi1.min(),phi2,P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    break;
  case 5:
    return w2.w_pi_phi2(costheta0,costheta1,costheta2,phi1,phi2.max(),P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi)
      -    w2.w_pi_phi2(costheta0,costheta1,costheta2,phi1,phi2.min(),P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    break;
  case 6:
    return w2.w_pi_P_b(costheta0,costheta1,costheta2,phi1,phi2,P_b.max(),alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi)
      -    w2.w_pi_P_b(costheta0,costheta1,costheta2,phi1,phi2,P_b.min(),alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    break;
  case 7:
    return w2.w_pi_alpha_lambda(costheta0,costheta1,costheta2,phi1,phi2,P_b,alpha_b,alpha_lambda.max(),r_0,r_1,alpha_plus,alpha_minus,chi)
      -    w2.w_pi_alpha_lambda(costheta0,costheta1,costheta2,phi1,phi2,P_b,alpha_b,alpha_lambda.min(),r_0,r_1,alpha_plus,alpha_minus,chi);
    break;

  case 8:
    val[0][0]=costheta0.min();
    val[0][1]=costheta0.max();
    val[1][0]=costheta1.min();
    val[1][1]=costheta1.max();
    val[2][0]=costheta2.min();
    val[2][1]=costheta2.max();
    val[3][0]=phi1.min();
    val[3][1]=phi1.max();
    val[4][0]=phi2.min();
    val[4][1]=phi2.max();
    
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
    	for (int c=0; c<=1; c++) {
    	  for (int d=0; d<=1; d++) {
    	    for (int e=0; e<=1; e++) {
    	      result= result + pow(-1., a+b+c+d+e+1) * w2.w_primitive(val[0][a],val[1][b],val[2][c],val[3][d],val[4][e],P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    	    }
    	  }
    	}
      }
    }
    
    return result;
    //      +    w2.wprimitive(costheta.max(),costheta1.max(),costheta2.max(),phi1.max(),phi2.max(),P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    break;

    
  case 9://for costheta0
    val[0][0]=costheta1.min();
    val[0][1]=costheta1.max();
    val[1][0]=costheta2.min();
    val[1][1]=costheta2.max();
    val[2][0]=phi1.min();
    val[2][1]=phi1.max();
    val[3][0]=phi2.min();
    val[3][1]=phi2.max();
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
    	for (int c=0; c<=1; c++) {
    	  for (int d=0; d<=1; d++) {
	    result= result + pow(-1., a+b+c+d) * w2.w_primitive_dcostheta0(costheta0,val[0][a],val[1][b],val[2][c],val[3][d],P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    	  }
    	}
      }
    }
    return result;
    break;

  case 10://for costheta1
    val[0][0]=costheta0.min();
    val[0][1]=costheta0.max();
    val[1][0]=costheta2.min();
    val[1][1]=costheta2.max();
    val[2][0]=phi1.min();
    val[2][1]=phi1.max();
    val[3][0]=phi2.min();
    val[3][1]=phi2.max();
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
    	for (int c=0; c<=1; c++) {
    	  for (int d=0; d<=1; d++) {
	    result= result + pow(-1., a+b+c+d) * w2.w_primitive_dcostheta1(val[0][a],costheta1,val[1][b],val[2][c],val[3][d],P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    	  }
    	}
      }
    }
    return result;
    break;


  case 11://for costheta2
    val[0][0]=costheta0.min();
    val[0][1]=costheta0.max();
    val[1][0]=costheta1.min();
    val[1][1]=costheta1.max();
    val[2][0]=phi1.min();
    val[2][1]=phi1.max();
    val[3][0]=phi2.min();
    val[3][1]=phi2.max();
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
    	for (int c=0; c<=1; c++) {
    	  for (int d=0; d<=1; d++) {
	    result= result + pow(-1., a+b+c+d) * w2.w_primitive_dcostheta2(val[0][a],val[1][b],costheta2,val[2][c],val[3][d],P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    	  }
    	}
      }
    }
    return result;
    break;


  case 12://for phi1
    val[0][0]=costheta0.min();
    val[0][1]=costheta0.max();
    val[1][0]=costheta1.min();
    val[1][1]=costheta1.max();
    val[2][0]=costheta2.min();
    val[2][1]=costheta2.max();
    val[3][0]=phi2.min();
    val[3][1]=phi2.max();
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
    	for (int c=0; c<=1; c++) {
    	  for (int d=0; d<=1; d++) {
	    result= result + pow(-1., a+b+c+d) * w2.w_primitive_dphi1(val[0][a],val[1][b],val[2][c],phi1,val[3][d],P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    	  }
    	}
      }
    }
    return result;
    break;


  case 13://for phi2
    val[0][0]=costheta0.min();
    val[0][1]=costheta0.max();
    val[1][0]=costheta1.min();
    val[1][1]=costheta1.max();
    val[2][0]=costheta2.min();
    val[2][1]=costheta2.max();
    val[3][0]=phi1.min();
    val[3][1]=phi1.max();
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
    	for (int c=0; c<=1; c++) {
    	  for (int d=0; d<=1; d++) {
	    result= result + pow(-1., a+b+c+d) * w2.w_primitive_dphi2(val[0][a],val[1][b],val[2][c],val[3][d],phi2,P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    	  }
    	}
      }
    }
    return result;
    break;

  case 14://for phi1, phi2
    val[0][0]=costheta0.min();
    val[0][1]=costheta0.max();
    val[1][0]=costheta1.min();
    val[1][1]=costheta1.max();
    val[2][0]=costheta2.min();
    val[2][1]=costheta2.max();
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
    	for (int c=0; c<=1; c++) {
	  result= result + pow(-1., 1+a+b+c) * w2.w_pi_phi2_phi1(val[0][a],val[1][b],val[2][c],phi1,phi2,P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);
    	}
      }
    }
    return result;
    break;


  default:
    return 0.;
  }
}
