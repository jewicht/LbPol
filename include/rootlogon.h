#include "TROOT.h"
#include "TStyle.h"
#include <iostream>

#ifndef ROOTLOGONH
#define ROOTLOGONH

  // Use times new roman, precision 2 
//  Int_t lhcbFont        = 132;  // Old LHCb style: 62;
static const Int_t  lhcbFont = 42;
  // Line thickness
static const  Double_t lhcbWidth    = 3.00; // Old LHCb style: 3.00;
  // Text size
static const Double_t lhcbTSize    = 0.06; 

void rootlogon();

void lhcbstyle();

#endif
