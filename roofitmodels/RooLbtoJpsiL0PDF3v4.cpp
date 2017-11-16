#include "RooLbtoJpsiL0PDF3v4.h"

ClassImp(  RooLbtoJpsiL0PDF3v4)

//_____________________________________________________________________________
RooLbtoJpsiL0PDF3v4::RooLbtoJpsiL0PDF3v4(const char *name, const char *title
					 , RooAbsReal& _costheta0, RooAbsReal& _costheta1, RooAbsReal& _costheta2
					 , RooAbsReal& _P_b, RooAbsReal&  _alpha_b,  RooAbsReal& _alpha_lambda,  RooAbsReal& _r_0,  RooAbsReal& _r_1
					 , RooAbsReal& _alpha_plus, RooAbsReal& _alpha_minus, RooAbsReal& _chi
					 , RooArgList& _varlist
					 //					 , RooArgList& _fXlist
				     ) : RooAbsPdf(name,title)
				       , costheta0("costheta0","Costheta0",this,_costheta0)
				       , costheta1("costheta1","Costheta1",this,_costheta1)
				       , costheta2("costheta2","Costheta2",this,_costheta2)
				       , P_b("P_b","P_b",this,_P_b)
				       , alpha_b("alpha_b","alpha_b",this,_alpha_b)
				       , alpha_lambda("alpha_lambda","alpha_lambda",this,_alpha_lambda)
				       , r_0("r_0","r_0",this,_r_0)
				       , r_1("r_1","r_1",this,_r_1)
  , alpha_plus("alpha_plus","alpha_plus",this,_alpha_plus)
  , alpha_minus("alpha_minus","alpha_minus",this,_alpha_minus)
  , chi("chi","chi",this,_chi)
                                       , varlist("varlist","varlist",this)  
  //                                       , fXlist("fXlist","fXlist",this)  
{
  TIterator* coefIter = _varlist.createIterator() ;
  RooAbsArg* coef ;
  
  while((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      std::cout << "RooLbtoJpsiL0PDF3v4::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() 
		<< " is not of type RooAbsReal" << std::endl ;
      assert(0) ;
    }
    varlist.add(*coef) ;
  }
  delete coefIter ;


  // coefIter = _fXlist.createIterator() ;
  
  // while((coef = (RooAbsArg*)coefIter->Next())) {
  //   if (!dynamic_cast<RooAbsReal*>(coef)) {
  //     cout << "RooLbtoJpsiL0PDF3v4::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() 
  // 	   << " is not of type RooAbsReal" << endl ;
  //     assert(0) ;
  //   }
  //   fXlist.add(*coef) ;
  // }
  // delete coefIter ;
  


  
  const int varlen=strlen(varlist[0].GetName());
  
  for (int i=0;i<varlist.getSize();i++) {
    
    const char* buf=varlist[i].GetName();
    const int len=strlen(buf);
    const int vari=10*chtoint(buf[len-4])+chtoint(buf[len-3]);
    const int varj=10*chtoint(buf[len-2])+chtoint(buf[len-1]);
    //    const int vark=10*chtoint(buf[len-2])+chtoint(buf[len-1]);
    
    if (len!=varlen or 
	vari<0 or vari>=ordermax2i or
	varj<0 or varj>=ordermax2j)  assert(i);
  }
}

//_____________________________________________________________________________
RooLbtoJpsiL0PDF3v4::RooLbtoJpsiL0PDF3v4(const RooLbtoJpsiL0PDF3v4& other, const char* name) : 
  RooAbsPdf(other,name), 
  costheta0("costheta0",this,other.costheta0), 
  costheta1("costheta1",this,other.costheta1), 
  costheta2("costheta2",this,other.costheta2), 
  P_b("P_b",this,other.P_b),
  alpha_b("alpha_b",this,other.alpha_b),
  alpha_lambda("alpha_lambda",this,other.alpha_lambda),
  r_0("r_0",this,other.r_0),
  r_1("r_1",this,other.r_1),
  alpha_plus("alpha_plus",this,other.alpha_plus),
  alpha_minus("alpha_minus",this,other.alpha_minus),
  chi("chi",this,other.chi),
  varlist("varlist",this,other.varlist)
  //  fXlist("fXlist",this,other.fXlist)
{
}

//_____________________________________________________________________________
Double_t RooLbtoJpsiL0PDF3v4::evaluate() const
{
  double c[ordermax2i][ordermax2j];
  memset(c,0,ordermax2i*ordermax2j*sizeof(double));
  
  const Int_t len=strlen(varlist[0].GetName());
  for (int i=0;i<varlist.getSize();i++) {
    const char* buf=varlist[i].GetName();
    const int vari=10*chtoint(buf[len-4])+chtoint(buf[len-3]);
    const int varj=10*chtoint(buf[len-2])+chtoint(buf[len-1]);
    //    const int vark=10*chtoint(buf[len-2])+chtoint(buf[len-1]);

    //    std::cout << varlist[i].GetName() << " " << vari << " " << varj << std::endl;
    c[vari][varj]=((RooAbsReal&)varlist[i]).getVal(); 
  }


  Double_t f[20], g[20];

  computew5 w5;
  w5.setparams(P_b, alpha_b, alpha_lambda, r_0, r_1, alpha_plus, alpha_minus, chi);
  for (int i=0;i<=19;i++) f[i]=w5.f1(i)*w5.f2(i);

  g[0]=1.;
  g[1]=alpha_b*P_b;
  g[2]=(2.*r_1-alpha_b)*alpha_lambda;
  g[3]=(2.*r_0-1.)*P_b*alpha_lambda;
  g[4]=0.5*(1.-3.*r_0);
  g[5]=0.5*(alpha_b-3.*r_1)*P_b;
  g[6]=-0.5*(alpha_b+r_1)*alpha_lambda;
  g[7]=-0.5*(1.+r_0)*P_b*alpha_lambda;

  // if (f[0]!=1.) {
  //   std::cout << f[0] << std::endl;
  //   std::cout << "alpha_b = " << alpha_b << std::endl;
  //   std::cout << "r_0 = " << r_0 << std::endl;
  //   std::cout << "r_1 = " << r_1 << std::endl;
  // }

  //  for (int i=0; i<8; i++) std::cout << "i=" << i << " " << f[i] << " vs " << g[i] << std::endl;



  // // for (int i=0;i<fXlist.getSize();i++) {
  // //   f[8+i]=((RooAbsReal&)fXlist[i]).getVal(); 
  // // }
  
  // for (int i=8;i<=15;i++) f[i]*=P_b*alpha_lambda;
  // f[16]*=P_b;
  // f[17]*=P_b;
  // f[18]*=alpha_lambda;
  // f[19]*=alpha_lambda;

  const Double_t pi=4.*atan(1.);
  const Double_t pi2=pi*pi;
  const Double_t pi3=pow(pi,3);
  const Double_t pi4=pow(pi,4);
  const Double_t pi5=pow(pi,5);
  const Double_t pi6=pow(pi,6);
  const Double_t pi7=pow(pi,7);
  const Double_t pi8=pow(pi,8);
  const Double_t pi9=pow(pi,9);
  const Double_t pi10=pow(pi,10);

  const Double_t costheta0sq=pow(costheta0,2);
  const Double_t costheta1sq=pow(costheta1,2);
  const Double_t costheta2sq=pow(costheta2,2);




  return 1./512.*(512.*pi10*c[0][0]*costheta0*costheta1*f[3] + 512.*pi10*c[0][0]*costheta0*f[1] + 512.*pi10*c[0][0]*costheta1*f[2] + 512.*pi10*c[0][0]*f[0] - sqrt(-costheta1sq + 1.)*sqrt(-costheta0sq + 1.)*(512.*(3.*pi8*c[2][0] - (3.*pi8*c[2][0] + 5.*(2.*pi8 - 21.*pi6)*c[4][0] + 21.*(pi8 - 60.*pi6 + 495.*pi4)*c[6][0] + 9*(4.*pi8 - 770.*pi6 + 30030.*pi4 - 225225.*pi2)*c[8][0])*costheta2sq + 5.*(2.*pi8 - 21.*pi6)*c[4][0] + 21.*(pi8 - 60.*pi6 + 495.*pi4)*c[6][0] + 9*(4.*pi8 - 770.*pi6 + 30030.*pi4 - 225225.*pi2)*c[8][0])*f[8] - 512.*(pi9*c[1][0] - (pi9*c[1][0] + (pi9 - 15.*pi7)*c[3][0] + (pi9 - 105.*pi7 + 945.*pi5)*c[5][0] + (pi9 - 378.*pi7 + 17325.*pi5 - 135135.*pi3)*c[7][0] + (34459425.*pi + pi9 - 990.*pi7 + 135135.*pi5 - 4729725.*pi3)*c[9][0])*costheta2sq + (pi9 - 15.*pi7)*c[3][0] + (pi9 - 105.*pi7 + 945.*pi5)*c[5][0] + (pi9 - 378.*pi7 + 17325.*pi5 - 135135.*pi3)*c[7][0] + (34459425.*pi + pi9 - 990.*pi7 + 135135.*pi5 - 4729725.*pi3)*c[9][0])*f[9] - 2.*(192.*pi7*c[1][2] + 384.*pi7*c[2][1] - (192.*pi7*c[1][2] + 384.*pi7*c[2][1] + 192.*(pi7 - 15.*pi5)*c[3][2] + 640.*(2.*pi7 - 21.*pi5)*c[4][1] + 96*(4.*pi7 - 15.*pi5)*c[2][3] + 80.*(8.*pi7 - 21.*pi5)*c[1][4] + 192.*(pi7 - 105.*pi5 + 945.*pi3)*c[5][2] + 2688.*(pi7 - 60.*pi5 + 495.*pi3)*c[6][1] + 80.*(8.*pi7 - 141.*pi5 + 315.*pi3)*c[3][4] + 160.*(8.*pi7 - 114.*pi5 + 315.*pi3)*c[4][3] + 24.*(16*pi7 - 420.*pi5 + 945.*pi3)*c[2][5] + 84.*(16*pi7 - 240.*pi5 + 495.*pi3)*c[1][6] - 84.*(7425.*pi - 16*pi7 + 480.*pi5 - 4095.*pi3)*c[3][6] - 672.*(7425.*pi - 4.*pi7 + 255.*pi5 - 2880.*pi3)*c[6][3] - 40.*(19845.*pi - 32.*pi7 + 1176*pi5 - 10710.*pi3)*c[4][5] - 80.*(19845.*pi - 8.*pi7 + 861.*pi5 - 9765.*pi3)*c[5][4] - 6*(135135.*pi - 64.*pi7 + 6048.*pi5 - 69300.*pi3)*c[2][7] - 192.*(135135.*pi - pi7 + 378.*pi5 - 17325.*pi3)*c[7][2] - 9*(225225.*pi - 256*pi7 + 12320.*pi5 - 120120.*pi3)*c[1][8] - 1152.*(225225.*pi - 4.*pi7 + 770.*pi5 - 30030.*pi3)*c[8][1])*costheta2sq + 192.*(pi7 - 15.*pi5)*c[3][2] + 640.*(2.*pi7 - 21.*pi5)*c[4][1] + 96*(4.*pi7 - 15.*pi5)*c[2][3] + 80.*(8.*pi7 - 21.*pi5)*c[1][4] + 192.*(pi7 - 105.*pi5 + 945.*pi3)*c[5][2] + 2688.*(pi7 - 60.*pi5 + 495.*pi3)*c[6][1] + 80.*(8.*pi7 - 141.*pi5 + 315.*pi3)*c[3][4] + 160.*(8.*pi7 - 114.*pi5 + 315.*pi3)*c[4][3] + 24.*(16*pi7 - 420.*pi5 + 945.*pi3)*c[2][5] + 84.*(16*pi7 - 240.*pi5 + 495.*pi3)*c[1][6] - 84.*(7425.*pi - 16*pi7 + 480.*pi5 - 4095.*pi3)*c[3][6] - 672.*(7425.*pi - 4.*pi7 + 255.*pi5 - 2880.*pi3)*c[6][3] - 40.*(19845.*pi - 32.*pi7 + 1176*pi5 - 10710.*pi3)*c[4][5] - 80.*(19845.*pi - 8.*pi7 + 861.*pi5 - 9765.*pi3)*c[5][4] - 6*(135135.*pi - 64.*pi7 + 6048.*pi5 - 69300.*pi3)*c[2][7] - 192.*(135135.*pi - pi7 + 378.*pi5 - 17325.*pi3)*c[7][2] - 9*(225225.*pi - 256*pi7 + 12320.*pi5 - 120120.*pi3)*c[1][8] - 1152.*(225225.*pi - 4.*pi7 + 770.*pi5 - 30030.*pi3)*c[8][1])*f[11] - (256*pi8*c[1][1] - 1152.*pi6*c[2][2] - (256*pi8*c[1][1] - 1152.*pi6*c[2][2] - 1920.*(2.*pi6 - 21.*pi4)*c[4][2] - 480.*(8.*pi6 - 21.*pi4)*c[2][4] + 256*(pi8 - 15.*pi6)*c[3][1] + 64.*(4.*pi8 - 15.*pi6)*c[1][3] - 8064.*(pi6 - 60.*pi4 + 495.*pi2)*c[6][2] - 3456*(4.*pi6 - 770.*pi4 + 30030.*pi2 - 225225)*c[8][2] - 3360.*(8.*pi6 - 501.*pi4 + 5220.*pi2 - 10395)*c[6][4] - 504.*(16*pi6 - 240.*pi4 + 495.*pi2)*c[2][6] - 800.*(16*pi6 - 210.*pi4 + 441.*pi2)*c[4][4] - 840.*(32.*pi6 - 816*pi4 + 6030.*pi2 - 10395)*c[4][6] - 54.*(256*pi6 - 12320.*pi4 + 120120.*pi2 - 225225)*c[2][8] + 256*(pi8 - 105.*pi6 + 945.*pi4)*c[5][1] + 64.*(4.*pi8 - 75.*pi6 + 225.*pi4)*c[3][3] + 16*(16*pi8 - 420.*pi6 + 945.*pi4)*c[1][5] + 256.*(pi8 - 990.*pi6 + 135135.*pi4 - 4729725.*pi2 + 34459425)*c[9][1] + 256*(pi8 - 378.*pi6 + 17325.*pi4 - 135135.*pi2)*c[7][1] + 64.*(4.*pi8 - 1527.*pi6 + 74970.*pi4 - 800415.*pi2 + 2027025)*c[7][3] + 64.*(4.*pi8 - 435.*pi6 + 5355.*pi4 - 14175.*pi2)*c[5][3] + 16*(16*pi8 - 2100.*pi6 + 60165.*pi4 - 496125.*pi2 + 893025)*c[5][5] + 16*(16*pi8 - 660.*pi6 + 7245.*pi4 - 14175.*pi2)*c[3][5] + 4.*(64.*pi8 - 7008.*pi6 + 160020.*pi4 - 1174635.*pi2 + 2027025)*c[3][7] + 4.*(64.*pi8 - 6048.*pi6 + 69300.*pi4 - 135135.*pi2)*c[1][7] + (256*pi8 - 63360.*pi6 + 2162160.*pi4 - 18918900.*pi2 + 34459425)*c[1][9])*costheta2sq - 1920.*(2.*pi6 - 21.*pi4)*c[4][2] - 480.*(8.*pi6 - 21.*pi4)*c[2][4] + 256*(pi8 - 15.*pi6)*c[3][1] + 64.*(4.*pi8 - 15.*pi6)*c[1][3] - 8064.*(pi6 - 60.*pi4 + 495.*pi2)*c[6][2] - 3456*(4.*pi6 - 770.*pi4 + 30030.*pi2 - 225225)*c[8][2] - 3360.*(8.*pi6 - 501.*pi4 + 5220.*pi2 - 10395)*c[6][4] - 504.*(16*pi6 - 240.*pi4 + 495.*pi2)*c[2][6] - 800.*(16*pi6 - 210.*pi4 + 441.*pi2)*c[4][4] - 840.*(32.*pi6 - 816*pi4 + 6030.*pi2 - 10395)*c[4][6] - 54.*(256*pi6 - 12320.*pi4 + 120120.*pi2 - 225225)*c[2][8] + 256*(pi8 - 105.*pi6 + 945.*pi4)*c[5][1] + 64.*(4.*pi8 - 75.*pi6 + 225.*pi4)*c[3][3] + 16*(16*pi8 - 420.*pi6 + 945.*pi4)*c[1][5] + 256*(pi8 - 990.*pi6 + 135135.*pi4 - 4729725.*pi2 + 34459425)*c[9][1] + 256*(pi8 - 378.*pi6 + 17325.*pi4 - 135135.*pi2)*c[7][1] + 64.*(4.*pi8 - 1527.*pi6 + 74970.*pi4 - 800415.*pi2 + 2027025)*c[7][3] + 64.*(4.*pi8 - 435.*pi6 + 5355.*pi4 - 14175.*pi2)*c[5][3] + 16*(16*pi8 - 2100.*pi6 + 60165.*pi4 - 496125.*pi2 + 893025)*c[5][5] + 16*(16*pi8 - 660.*pi6 + 7245.*pi4 - 14175.*pi2)*c[3][5] + 4.*(64.*pi8 - 7008.*pi6 + 160020.*pi4 - 1174635.*pi2 + 2027025)*c[3][7] + 4.*(64.*pi8 - 6048.*pi6 + 69300.*pi4 - 135135.*pi2)*c[1][7] + (256.*pi8 - 63360.*pi6 + 2162160.*pi4 - 18918900.*pi2 + 34459425)*c[1][9])*f[10]) + 256*(3.*pi10*c[0][0]*costheta2sq - pi10*c[0][0])*f[4] + 256*(3.*pi10*c[0][0]*costheta1*costheta2sq - pi10*c[0][0]*costheta1)*f[6] + 256*(3.*pi10*c[0][0]*costheta0*costheta2sq - pi10*c[0][0]*costheta0)*f[5] + 256*(3.*pi10*c[0][0]*costheta0*costheta1*costheta2sq - pi10*c[0][0]*costheta0*costheta1)*f[7] - 512.*sqrt(-costheta2sq + 1.)*(sqrt(-costheta1sq + 1.)*((3.*pi7*c[1][2] + 3.*pi7*c[2][1] + 3.*(pi7 - 15.*pi5)*c[2][3] + 3.*(pi7 - 15.*pi5)*c[3][2] + 5.*(2.*pi7 - 21.*pi5)*c[1][4] + 5.*(2.*pi7 - 21.*pi5)*c[4][1] + 3.*(pi7 - 105.*pi5 + 945.*pi3)*c[2][5] + 3.*(pi7 - 105.*pi5 + 945.*pi3)*c[5][2] + 21.*(pi7 - 60.*pi5 + 495.*pi3)*c[1][6] + 21.*(pi7 - 60.*pi5 + 495.*pi3)*c[6][1] + 5.*(2.*pi7 - 51.*pi5 + 315.*pi3)*c[3][4] + 5.*(2.*pi7 - 51.*pi5 + 315.*pi3)*c[4][3] - 21.*(7425.*pi - pi7 + 75.*pi5 - 1395.*pi3)*c[3][6] - 21.*(7425.*pi - pi7 + 75.*pi5 - 1395.*pi3)*c[6][3] - 5.*(19845.*pi - 2.*pi7 + 231.*pi5 - 4095.*pi3)*c[4][5] - 5.*(19845.*pi - 2.*pi7 + 231.*pi5 - 4095.*pi3)*c[5][4] - 3.*(135135.*pi - pi7 + 378.*pi5 - 17325.*pi3)*c[2][7] - 3.*(135135.*pi - pi7 + 378.*pi5 - 17325.*pi3)*c[7][2] - 9*(225225.*pi - 4.*pi7 + 770.*pi5 - 30030.*pi3)*c[1][8] - 9*(225225.*pi - 4.*pi7 + 770.*pi5 - 30030.*pi3)*c[8][1])*costheta0*costheta2*f[15] + (pi8*c[1][1] - 9*pi6*c[2][2] - 15.*(2.*pi6 - 21.*pi4)*c[2][4] - 15.*(2.*pi6 - 21.*pi4)*c[4][2] + (pi8 - 15.*pi6)*c[1][3] + (pi8 - 15.*pi6)*c[3][1] - 63.*(pi6 - 60.*pi4 + 495.*pi2)*c[2][6] - 63.*(pi6 - 60.*pi4 + 495.*pi2)*c[6][2] - 105.*(2.*pi6 - 141.*pi4 + 2250.*pi2 - 10395)*c[4][6] - 105.*(2.*pi6 - 141.*pi4 + 2250.*pi2 - 10395)*c[6][4] - 27.*(4.*pi6 - 770.*pi4 + 30030.*pi2 - 225225)*c[2][8] - 27.*(4.*pi6 - 770.*pi4 + 30030.*pi2 - 225225)*c[8][2] - 25.*(4.*pi6 - 84.*pi4 + 441.*pi2)*c[4][4] + (pi8 - 105.*pi6 + 945.*pi4)*c[1][5] + (pi8 - 105.*pi6 + 945.*pi4)*c[5][1] + (pi8 - 30.*pi6 + 225.*pi4)*c[3][3] + (pi8 - 990.*pi6 + 135135.*pi4 - 4729725.*pi2 + 34459425)*c[1][9] + (pi8 - 990.*pi6 + 135135.*pi4 - 4729725.*pi2 + 34459425)*c[9][1] + (pi8 - 393.*pi6 + 22995.*pi4 - 395010.*pi2 + 2027025)*c[3][7] + (pi8 - 393.*pi6 + 22995.*pi4 - 395010.*pi2 + 2027025)*c[7][3] + (pi8 - 378.*pi6 + 17325.*pi4 - 135135.*pi2)*c[1][7] + (pi8 - 378.*pi6 + 17325.*pi4 - 135135.*pi2)*c[7][1] + (pi8 - 210.*pi6 + 12915.*pi4 - 198450.*pi2 + 893025)*c[5][5] + (pi8 - 120.*pi6 + 2520.*pi4 - 14175.*pi2)*c[3][5] + (pi8 - 120.*pi6 + 2520.*pi4 - 14175.*pi2)*c[5][3])*costheta0*costheta2*f[14] + (3.*pi7*c[1][2] + 3.*pi7*c[2][1] + 3.*(pi7 - 15.*pi5)*c[2][3] + 3.*(pi7 - 15.*pi5)*c[3][2] + 5.*(2.*pi7 - 21.*pi5)*c[1][4] + 5.*(2.*pi7 - 21.*pi5)*c[4][1] + 3.*(pi7 - 105.*pi5 + 945.*pi3)*c[2][5] + 3.*(pi7 - 105.*pi5 + 945.*pi3)*c[5][2] + 21.*(pi7 - 60.*pi5 + 495.*pi3)*c[1][6] + 21.*(pi7 - 60.*pi5 + 495.*pi3)*c[6][1] + 5.*(2.*pi7 - 51.*pi5 + 315.*pi3)*c[3][4] + 5.*(2.*pi7 - 51.*pi5 + 315.*pi3)*c[4][3] - 21.*(7425.*pi - pi7 + 75.*pi5 - 1395.*pi3)*c[3][6] - 21.*(7425.*pi - pi7 + 75.*pi5 - 1395.*pi3)*c[6][3] - 5.*(19845.*pi - 2.*pi7 + 231.*pi5 - 4095.*pi3)*c[4][5] - 5.*(19845.*pi - 2.*pi7 + 231.*pi5 - 4095.*pi3)*c[5][4] - 3.*(135135.*pi - pi7 + 378.*pi5 - 17325.*pi3)*c[2][7] - 3.*(135135.*pi - pi7 + 378.*pi5 - 17325.*pi3)*c[7][2] - 9*(225225.*pi - 4.*pi7 + 770.*pi5 - 30030.*pi3)*c[1][8] - 9*(225225.*pi - 4.*pi7 + 770.*pi5 - 30030.*pi3)*c[8][1])*costheta2*f[19] + (pi8*c[1][1] - 9*pi6*c[2][2] - 15.*(2.*pi6 - 21.*pi4)*c[2][4] - 15.*(2.*pi6 - 21.*pi4)*c[4][2] + (pi8 - 15.*pi6)*c[1][3] + (pi8 - 15.*pi6)*c[3][1] - 63.*(pi6 - 60.*pi4 + 495.*pi2)*c[2][6] - 63.*(pi6 - 60.*pi4 + 495.*pi2)*c[6][2] - 105.*(2.*pi6 - 141.*pi4 + 2250.*pi2 - 10395)*c[4][6] - 105.*(2.*pi6 - 141.*pi4 + 2250.*pi2 - 10395)*c[6][4] - 27.*(4.*pi6 - 770.*pi4 + 30030.*pi2 - 225225)*c[2][8] - 27.*(4.*pi6 - 770.*pi4 + 30030.*pi2 - 225225)*c[8][2] - 25.*(4.*pi6 - 84.*pi4 + 441.*pi2)*c[4][4] + (pi8 - 105.*pi6 + 945.*pi4)*c[1][5] + (pi8 - 105.*pi6 + 945.*pi4)*c[5][1] + (pi8 - 30.*pi6 + 225.*pi4)*c[3][3] + (pi8 - 990.*pi6 + 135135.*pi4 - 4729725.*pi2 + 34459425)*c[1][9] + (pi8 - 990.*pi6 + 135135.*pi4 - 4729725.*pi2 + 34459425)*c[9][1] + (pi8 - 393.*pi6 + 22995.*pi4 - 395010.*pi2 + 2027025)*c[3][7] + (pi8 - 393.*pi6 + 22995.*pi4 - 395010.*pi2 + 2027025)*c[7][3] + (pi8 - 378.*pi6 + 17325.*pi4 - 135135.*pi2)*c[1][7] + (pi8 - 378.*pi6 + 17325.*pi4 - 135135.*pi2)*c[7][1] + (pi8 - 210.*pi6 + 12915.*pi4 - 198450.*pi2 + 893025)*c[5][5] + (pi8 - 120.*pi6 + 2520.*pi4 - 14175.*pi2)*c[3][5] + (pi8 - 120.*pi6 + 2520.*pi4 - 14175.*pi2)*c[5][3])*costheta2*f[18]) + sqrt(-costheta0sq + 1.)*((3.*pi8*c[0][2] + 5.*(2.*pi8 - 21.*pi6)*c[0][4] + 21.*(pi8 - 60.*pi6 + 495.*pi4)*c[0][6] + 9*(4.*pi8 - 770.*pi6 + 30030.*pi4 - 225225.*pi2)*c[0][8])*costheta1*costheta2*f[12] - (pi9*c[0][1] + (pi9 - 15.*pi7)*c[0][3] + (pi9 - 105.*pi7 + 945.*pi5)*c[0][5] + (pi9 - 378.*pi7 + 17325.*pi5 - 135135.*pi3)*c[0][7] + (34459425.*pi + pi9 - 990.*pi7 + 135135.*pi5 - 4729725.*pi3)*c[0][9])*costheta1*costheta2*f[13] + (3.*pi8*c[0][2] + 5.*(2.*pi8 - 21.*pi6)*c[0][4] + 21.*(pi8 - 60.*pi6 + 495.*pi4)*c[0][6] + 9*(4.*pi8 - 770.*pi6 + 30030.*pi4 - 225225.*pi2)*c[0][8])*costheta2*f[16] - (pi9*c[0][1] + (pi9 - 15.*pi7)*c[0][3] + (pi9 - 105.*pi7 + 945.*pi5)*c[0][5] + (pi9 - 378.*pi7 + 17325.*pi5 - 135135.*pi3)*c[0][7] + (34459425.*pi + pi9 - 990.*pi7 + 135135.*pi5 - 4729725.*pi3)*c[0][9])*costheta2*f[17])))/(pi10*c[0][0]);





}

//_____________________________________________________________________________
Int_t RooLbtoJpsiL0PDF3v4::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  // //3
  // if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(),costheta1.arg(),costheta2.arg()) ) ) return 9;

  // //2
  // if ( matchArgs(allVars,analVars, RooArgSet(costheta1.arg(), costheta2.arg()))) return  6;
  // if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(), costheta2.arg()))) return  7;
  // if ( matchArgs(allVars,analVars, RooArgSet(costheta0.arg(), costheta1.arg()))) return  8;

  // //1
  // if ( matchArgs(allVars,analVars,costheta0)    ) return 1 ;
  // if ( matchArgs(allVars,analVars,costheta1)    ) return 2 ;
  // if ( matchArgs(allVars,analVars,costheta2)    ) return 3 ;
  // if ( matchArgs(allVars,analVars,P_b)          ) return 4 ;
  // if ( matchArgs(allVars,analVars,alpha_lambda) ) return 5 ;
  
  return 0;
}



//_____________________________________________________________________________
Double_t RooLbtoJpsiL0PDF3v4::analyticalIntegral(Int_t code, const char* rangeName) const
{
  // computew3 w3;
  // Double_t result=0.;
  // Double_t val[5][2];

  // // Double_t xivector[8];
  // // xivector[0]=xi0;
  // // xivector[1]=xi1;
  // // xivector[2]=xi2;
  // // xivector[3]=xi3;
  // // xivector[4]=xi4;
  // // xivector[5]=xi5;
  // // xivector[6]=xi6;
  // // xivector[7]=xi7;

  // switch (code) {
  // case 1: 
  //   return w3.w_pi_costheta0(costheta0.max(),costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1)
  //     -    w3.w_pi_costheta0(costheta0.min(),costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);
  //   break;
  // case 2: 
  //   return w3.w_pi_costheta1(costheta0,costheta1.max(),costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1)
  //     -    w3.w_pi_costheta1(costheta0,costheta1.min(),costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);
  //   break;
  // case 3:
  //   return w3.w_pi_costheta2(costheta0,costheta1,costheta2.max(),P_b,alpha_b,alpha_lambda,r_0,r_1)
  //     -    w3.w_pi_costheta2(costheta0,costheta1,costheta2.min(),P_b,alpha_b,alpha_lambda,r_0,r_1);
  //   break;
  // case 4:
  //   return w3.w_pi_P_b(costheta0,costheta1,costheta2,P_b.max(),alpha_b,alpha_lambda,r_0,r_1)
  //     -    w3.w_pi_P_b(costheta0,costheta1,costheta2,P_b.min(),alpha_b,alpha_lambda,r_0,r_1);
  //   break;
  // case 5:
  //   return w3.w_pi_alpha_lambda(costheta0,costheta1,costheta2,P_b,alpha_b,alpha_lambda.max(),r_0,r_1)
  //     -    w3.w_pi_alpha_lambda(costheta0,costheta1,costheta2,P_b,alpha_b,alpha_lambda.min(),r_0,r_1);
  //   break;

  // case 9:
  //   val[0][0]=costheta0.min();
  //   val[0][1]=costheta0.max();
  //   val[1][0]=costheta1.min();
  //   val[1][1]=costheta1.max();
  //   val[2][0]=costheta2.min();
  //   val[2][1]=costheta2.max();
    
  //   for (int a=0; a<=1; a++) {
  //     for (int b=0; b<=1; b++) {
  //   	for (int c=0; c<=1; c++) {
  // 	  result+=pow(-1., a+b+c+1) * w3.w_primitive(val[0][a],val[1][b],val[2][c],P_b,alpha_b,alpha_lambda,r_0,r_1);
  //   	}
  //     }
  //   }
    
  //   return result;
  //   //      +    w3.wprimitive(costheta.max(),costheta1.max(),costheta2.max(),phi1.max(),phi2.max(),P_b,alpha_b,alpha_lambda,r_0,r_1);
  //   break;

    
  // case 6://for costheta0
  //   val[0][0]=costheta1.min();
  //   val[0][1]=costheta1.max();
  //   val[1][0]=costheta2.min();
  //   val[1][1]=costheta2.max();
  //   for (int a=0; a<=1; a++) {
  //     for (int b=0; b<=1; b++) {
  // 	result+=pow(-1., a+b) * w3.w_pi_costheta1_costheta2(costheta0,val[0][a],val[1][b],P_b,alpha_b,alpha_lambda,r_0,r_1);
  //     }
  //   }
  //   return result;
  //   break;

  // case 7://for costheta1
  //   val[0][0]=costheta0.min();
  //   val[0][1]=costheta0.max();
  //   val[1][0]=costheta2.min();
  //   val[1][1]=costheta2.max();

  //   for (int a=0; a<=1; a++) {
  //     for (int b=0; b<=1; b++) {
  // 	result+=pow(-1., a+b) * w3.w_pi_costheta0_costheta2(val[0][a],costheta1,val[1][b],P_b,alpha_b,alpha_lambda,r_0,r_1);
  //     }
  //   }
  //   return result;
  //   break;


  // case 8://for costheta2
  //   val[0][0]=costheta0.min();
  //   val[0][1]=costheta0.max();
  //   val[1][0]=costheta1.min();
  //   val[1][1]=costheta1.max();
  //   for (int a=0; a<=1; a++) {
  //     for (int b=0; b<=1; b++) {
  // 	result+=pow(-1., a+b) * w3.w_pi_costheta0_costheta1(val[0][a],val[1][b],costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);
  //     }
  //   }
  //   return result;
  //   break;


  // default:
  //   return 0.;
  // }

  return 0.;
}
