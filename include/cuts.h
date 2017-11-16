#include "TCut.h"
#include "TString.h"

#include "RooRealVar.h"
#include "functions.h"


const double Jpsipdgmass=3096.916;

double getLambda0mass(bool);

const TCut Lambda0_LL("piminus_TRACK_Type==3");
const TCut Lambda0_DD("piminus_TRACK_Type==5");

const TCut MagnetUp("Polarity==1");
const TCut MagnetDown("Polarity==-1");


void varnames(bool, TString&, TString&, TString&);
void simplecut(bool, TCut&, TCut&);
void simplecut2(bool, TCut&, TCut&);
void tmvarectcut(bool, TCut&, TCut&);
void bdtsel(bool, TCut&, TCut&);
//TCut pidcut(bool);

TCut truemcsel(bool);
TCut masscut(bool, const RooRealVar&);

TCut triggercut(bool);
