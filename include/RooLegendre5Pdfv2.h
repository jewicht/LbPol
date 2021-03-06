/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOLEGENDRE5PDFV2
#define ROOLEGENDRE5PDFV2

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "functions.h"
 
class RooLegendre5Pdfv2 : public RooAbsPdf {
public:
  RooLegendre5Pdfv2() {} ; 
  RooLegendre5Pdfv2(const char *name, const char *title,
		    RooAbsReal& _costheta0,
		    RooAbsReal& _costheta1,
		    RooAbsReal& _costheta2,
		    RooAbsReal& _phi1,
		    RooAbsReal& _phi2,
		    RooArgList& _varlist);

  RooLegendre5Pdfv2(const RooLegendre5Pdfv2& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooLegendre5Pdfv2(*this,newname); }
  inline virtual ~RooLegendre5Pdfv2() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:
  Double_t evaluate() const ;

  RooRealProxy costheta0 ;
  RooRealProxy costheta1 ;
  RooRealProxy costheta2 ;
  RooRealProxy phi1 ;
  RooRealProxy phi2 ;
  RooListProxy varlist;
  
private:
  ClassDef(RooLegendre5Pdfv2,1) // Your description goes here...
};
 
#endif
