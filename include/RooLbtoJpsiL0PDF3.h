
#ifndef ROO_LBTOJPSIL0PDF3
#define ROO_LBTOJPSIL0PDF3

#include "RooFit.h"
#include "Riostream.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"

#include "computew3.h"

class RooLbtoJpsiL0PDF3 : public RooAbsPdf {
public:
  RooLbtoJpsiL0PDF3() {} ;
  RooLbtoJpsiL0PDF3(const char *name, const char *title,
		    RooAbsReal& _costheta0, RooAbsReal& _costheta1, RooAbsReal& _costheta2,
		    RooAbsReal& P_b, RooAbsReal&  alpha_b,  RooAbsReal& alpha_lambda,  RooAbsReal& r_0,  RooAbsReal& r_1//,
		    //		    RooAbsReal& xi0, RooAbsReal& xi1, RooAbsReal& xi2, RooAbsReal& xi3, 
		    //		    RooAbsReal& xi4, RooAbsReal& xi5, RooAbsReal& xi6, RooAbsReal& xi7
		    );
  RooLbtoJpsiL0PDF3(const RooLbtoJpsiL0PDF3& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooLbtoJpsiL0PDF3(*this,newname); }
  inline virtual ~RooLbtoJpsiL0PDF3() { }

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
  //  RooRealProxy xi0,xi1,xi2,xi3,xi4,xi5,xi6,xi7;
  
  Double_t evaluate() const ;

private:
  ClassDef(RooLbtoJpsiL0PDF3,1) // Gaussian PDF
};

#endif
