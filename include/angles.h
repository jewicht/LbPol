#include "TLorentzVector.h"
#include "TVector3.h"
#include "functions.h"

TLorentzVector myboost(const TLorentzVector&, const TVector3&);
void calculateangles(const TLorentzVector&, 
		     const TLorentzVector&, const TLorentzVector&, const TLorentzVector&, 
		     const TLorentzVector&, const TLorentzVector&, const TLorentzVector&,
		     Float_t&, Float_t&, Float_t&, Float_t&, Float_t&, const TVector3&);

void calculateangles(const TLorentzVector&, 
		     const TLorentzVector&, const TLorentzVector&, const TLorentzVector&, 
		     const TLorentzVector&, const TLorentzVector&, const TLorentzVector&,
		     Float_t&, Float_t&, Float_t&, Float_t&, Float_t&);

TString lvcout(const TLorentzVector&);
TString lvcout(const TVector3&);
