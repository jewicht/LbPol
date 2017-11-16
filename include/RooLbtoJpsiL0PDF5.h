
#ifndef ROO_LBTOJPSIL0PDF2
#define ROO_LBTOJPSIL0PDF2

#include "RooFit.h"
#include "Riostream.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"

#include "computew5.h"

class RooLbtoJpsiL0PDF5 : public RooAbsPdf {
public:
  RooLbtoJpsiL0PDF5() {} ;
  RooLbtoJpsiL0PDF5(const char *name, const char *title,
		   RooAbsReal& _costheta0, RooAbsReal& _costheta1, RooAbsReal& _costheta2, RooAbsReal& _phi1, RooAbsReal& _phi2, 
		   RooAbsReal& P_b, RooAbsReal&  alpha_b,  RooAbsReal& alpha_lambda,  RooAbsReal& r_0,  RooAbsReal& r_1, RooAbsReal&  alpha_plus,  RooAbsReal& alpha_minus,  RooAbsReal& chi);
  RooLbtoJpsiL0PDF5(const RooLbtoJpsiL0PDF5& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooLbtoJpsiL0PDF5(*this,newname); }
  inline virtual ~RooLbtoJpsiL0PDF5() { }

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
  RooRealProxy costheta0, costheta1, costheta2, phi1, phi2 ;
  //param
  RooRealProxy P_b, alpha_b, alpha_lambda, r_0, r_1, alpha_plus, alpha_minus, chi;


  Double_t evaluate() const ;

private:
  ClassDef(RooLbtoJpsiL0PDF5,1) // Gaussian PDF
};

#endif
