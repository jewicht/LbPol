
#ifndef ROO_LBTOJPSIL0PDF3V4
#define ROO_LBTOJPSIL0PDF3V4

#include "RooFit.h"
#include "Riostream.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooArgList.h"
#include "legendre2f2.h"
#include "functions.h"
#include "computew5.h"

class RooLbtoJpsiL0PDF3v4 : public RooAbsPdf {
public:
  RooLbtoJpsiL0PDF3v4() {} ;
  RooLbtoJpsiL0PDF3v4(const char *name, const char *title,
		      RooAbsReal& _costheta0, RooAbsReal& _costheta1, RooAbsReal& _costheta2,
		      RooAbsReal& _P_b, RooAbsReal&  _alpha_b,  RooAbsReal& _alpha_lambda,  RooAbsReal& _r_0,  RooAbsReal& _r_1,
		      RooAbsReal& _alpha_plus, RooAbsReal& _alpha_minus, RooAbsReal& _chi,
		      RooArgList& _varlist
		    );
  RooLbtoJpsiL0PDF3v4(const RooLbtoJpsiL0PDF3v4& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooLbtoJpsiL0PDF3v4(*this,newname); }
  inline virtual ~RooLbtoJpsiL0PDF3v4() { }

  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  //  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  //  void generateEvent(Int_t code);

  //

  //private:
  /* Double_t x[5]; */
  /* Double_t par[8]; */
  /* Double_t x0[5]; */
  /* Double_t par0[8]; */
  /* Double_t x1[5]; */
  /* Double_t par1[8]; */


protected:

  //var
  RooRealProxy costheta0, costheta1, costheta2;
  //param
  RooRealProxy P_b, alpha_b, alpha_lambda, r_0, r_1;
  RooRealProxy alpha_plus, alpha_minus, chi;
  RooListProxy varlist;
  //  RooListProxy fXlist;

  Double_t evaluate() const ;

private:
  ClassDef(RooLbtoJpsiL0PDF3v4,1) // Gaussian PDF
};

#endif
