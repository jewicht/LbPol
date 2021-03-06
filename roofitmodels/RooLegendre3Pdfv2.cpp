/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooLegendre3Pdfv2.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

#include "legendre3f3.h"

ClassImp(RooLegendre3Pdfv2) 

RooLegendre3Pdfv2::RooLegendre3Pdfv2(const char *name, const char *title, 
				     RooAbsReal& _costheta0,
				     RooAbsReal& _costheta1,
				     RooAbsReal& _costheta2,
				     RooArgList& _varlist) :
  RooAbsPdf(name,title), 
  costheta0("costheta0","costheta0",this,_costheta0),
  costheta1("costheta1","costheta1",this,_costheta1),
  costheta2("costheta2","costheta2",this,_costheta2),
  varlist("varlist","varlist",this)
{ 
  
  TIterator* coefIter = _varlist.createIterator() ;
  RooAbsArg* coef ;
  
  while((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      std::cout << "RooLegendre3Pdfv2::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() 
		<< " is not of type RooAbsReal" << std::endl ;
      assert(0) ;
    }
    varlist.add(*coef) ;
  }
  delete coefIter ;
  
  
  const int varlen=strlen(varlist[0].GetName());
  
  for (int i=0;i<varlist.getSize();i++) {
    
    const char* buf=varlist[i].GetName();
    const int len=strlen(buf);
    const int vari=10*chtoint(buf[len-6])+chtoint(buf[len-5]);
    const int varj=10*chtoint(buf[len-4])+chtoint(buf[len-3]);
    const int vark=10*chtoint(buf[len-2])+chtoint(buf[len-1]);
    
    if (len!=varlen or 
	vari<0 or vari>=ordermax3i or
	varj<0 or varj>=ordermax3j or
	vark<0 or vark>=ordermax3k) assert(i);
  }
} 


RooLegendre3Pdfv2::RooLegendre3Pdfv2(const RooLegendre3Pdfv2& other, const char* name) :  
  RooAbsPdf(other,name), 
  costheta0("costheta0",this,other.costheta0),
  costheta1("costheta1",this,other.costheta1),
  costheta2("costheta2",this,other.costheta2),
  varlist("varlist",this,other.varlist)
{ 
} 



Double_t RooLegendre3Pdfv2::evaluate() const 
{ 
  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
  

  double coeff3[ordermax3i][ordermax3j][ordermax3k];
  memset(coeff3,0,ordermax3i*ordermax3j*ordermax3k*sizeof(double));
  
  const Int_t len=strlen(varlist[0].GetName());
  for (int i=0;i<varlist.getSize();i++) {
    const char* buf=varlist[i].GetName();
    const int vari=10*chtoint(buf[len-6])+chtoint(buf[len-5]);
    const int varj=10*chtoint(buf[len-4])+chtoint(buf[len-3]);
    const int vark=10*chtoint(buf[len-2])+chtoint(buf[len-1]);
    coeff3[vari][varj][vark]=((RooAbsReal&)varlist[i]).getVal(); 
  }
  
  const double value=pdfval(costheta0,costheta1,costheta2,coeff3) ; 
  //  if (costheta0<-0.999 and costheta1>0.999 and costheta2<-0.999) std::cerr << value << std::endl;
  
  return value;
  //  return (value>0.)?value:0.;
} 



Int_t RooLegendre3Pdfv2::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const  
{ 
  // LIST HERE OVER WHICH VARIABLES ANALYTICAL INTEGRATION IS SUPPORTED, 
  // ASSIGN A NUMERIC CODE FOR EACH SUPPORTED (SET OF) PARAMETERS 
  // THE EXAMPLE BELOW ASSIGNS CODE 1 TO INTEGRATION OVER VARIABLE X
  // YOU CAN ALSO IMPLEMENT MORE THAN ONE ANALYTICAL INTEGRAL BY REPEATING THE matchArgs 
  // EXPRESSION MULTIPLE TIMES
  
  // if (matchArgs(allVars,analVars,x)) return 1 ;    
  if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(),costheta1.arg(),costheta2.arg()) ) ) return 1;
  if ( matchArgs(allVars,analVars, RooArgSet(costheta1.arg(),costheta2.arg()) ) ) return 2;
  if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(),costheta2.arg()) ) ) return 3;
  if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(),costheta1.arg()) ) ) return 4;
  
  return 0;
} 



Double_t RooLegendre3Pdfv2::analyticalIntegral(Int_t code, const char* rangeName) const  
{ 
  // RETURN ANALYTICAL INTEGRAL DEFINED BY RETURN CODE ASSIGNED BY getAnalyticalIntegral
  // THE MEMBER FUNCTION x.min(rangeName) AND x.max(rangeName) WILL RETURN THE INTEGRATION
  // BOUNDARIES FOR EACH OBSERVABLE x

  // assert(code==1) ; 
  // return (x.max(rangeName)-x.min(rangeName)) ; 
  
  double coeff3[ordermax3i][ordermax3j][ordermax3k];
  memset(coeff3,0,ordermax3i*ordermax3j*ordermax3k*sizeof(double));
  
  const Int_t len=strlen(varlist[0].GetName());
  for (int i=0;i<varlist.getSize();i++) {
    const char* buf=varlist[i].GetName();
    const int vari=10*chtoint(buf[len-6])+chtoint(buf[len-5]);
    const int varj=10*chtoint(buf[len-4])+chtoint(buf[len-3]);
    const int vark=10*chtoint(buf[len-2])+chtoint(buf[len-1]);
    coeff3[vari][varj][vark]=((RooAbsReal&)varlist[i]).getVal(); 
  }

  Double_t val[3][2];
  memset(val, 0, 3*2*sizeof(Double_t));
  val[0][0]=costheta0.min(rangeName);
  val[0][1]=costheta0.max(rangeName);
  val[1][0]=costheta1.min(rangeName);
  val[1][1]=costheta1.max(rangeName);
  val[2][0]=costheta2.min(rangeName);
  val[2][1]=costheta2.max(rangeName);
  
  Double_t result=0.;
  switch(code) {
  case 1: 
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
	for (int c=0; c<=1; c++) {
	  result+= pow(-1., a+b+c+1) * primitiveval(val[0][a],val[1][b],val[2][c],coeff3);
	}
      }
    }
    return result;
    break;
  case 2://for costheta0
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
	result+=pow(-1., a+b) * primitivedcosthetaval(costheta0,val[1][a],val[2][b],coeff3);
      }
    }
    return result;
    break;
    
  case 3://for costheta1
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
	result+=pow(-1., a+b) * primitivedcostheta1val(val[0][a],costheta1,val[2][b],coeff3);
      }
    }
    return result;
    break;
    
  case 4://for costheta2
    for (int a=0; a<=1; a++) {
      for (int b=0; b<=1; b++) {
	result+=pow(-1., a+b) * primitivedcostheta2val(val[0][a],val[1][b],costheta2,coeff3);
      }
    }
    return result;
    break;
    
  }

  return 0 ; 
} 



