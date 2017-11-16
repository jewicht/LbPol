
#ifndef ROO_LBTOJPSIL0PDF3V3
#define ROO_LBTOJPSIL0PDF3V3

#include "RooFit.h"
#include "Riostream.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"

#include "computew3v3.h"

class RooLbtoJpsiL0PDF3v3 : public RooAbsPdf {
public:
  RooLbtoJpsiL0PDF3v3() {} ;
  RooLbtoJpsiL0PDF3v3(const char *name, const char *title,
		      RooAbsReal& _costheta0, RooAbsReal& _costheta1, RooAbsReal& _costheta2,
		      RooAbsReal& P_b, RooAbsReal&  alpha_b,  RooAbsReal& alpha_lambda,  RooAbsReal& beta,  RooAbsReal& gamma
		      );
  RooLbtoJpsiL0PDF3v3(const RooLbtoJpsiL0PDF3v3& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooLbtoJpsiL0PDF3v3(*this,newname); }
  inline virtual ~RooLbtoJpsiL0PDF3v3() { }

  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  //  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  //  void generateEvent(Int_t code);

protected:

  //var
  RooRealProxy costheta0, costheta1, costheta2;
  //param
  RooRealProxy P_b,  ap2,alpha_lambda, am2, bp2;
  
  Double_t evaluate() const ;

private:

    ClassDef(RooLbtoJpsiL0PDF3v3,1) // Gaussian PDF

};

#endif
