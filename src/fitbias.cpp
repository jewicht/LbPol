#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TRandom3.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

#include "rootlogon.h"
#include "functions-roofit.h"
#include "weighting.h"



TMatrixDSym* computeCorrMatrix(const TMatrixD& _VM)
{
// Now construct correlation matrix from it
  TMatrixDSym* _CM = (TMatrixDSym*) _VM.Clone() ;
  for (Int_t i=0 ; i<_CM->GetNrows() ; i++) {
    for (Int_t j=0 ; j<_CM->GetNcols() ; j++) {
      if (i!=j) {
	(*_CM)(i,j) = (*_CM)(i,j) / sqrt((*_CM)(i,i)*(*_CM)(j,j)) ;
      }
    }
  }
  for (Int_t i=0 ; i<_CM->GetNrows() ; i++) {
    (*_CM)(i,i) = 1.0 ;
  }
  return _CM;
}



void printrrv(const RooRealVar& var)
{
  std::cout << var.GetName() << " = " << var.getVal() << " +/- " << var.getError() << std::endl;
}

int main(int argc, char** argv)
{

  lhcbstyle();



  TCanvas c("c","c",100,100);

  RooRealVar P_b("P_b","P_{b}",0.);
  RooRealVar alpha_b("alpha_b","#alpha_{b}",0.);//0.457
  RooRealVar r_0("r_0","r_{0}",0.);//0.5
  RooRealVar r_1("r_1","r_{1}",0.);//0.1

  RooRealVar ap2("ap2","|a+|^{2}",0.);
  RooRealVar am2("am2","|a-|^{2}",0.);
  RooRealVar bp2("bp2","|b+|^{2}",0.);
  RooRealVar bm2("bm2","|b-|^{2}",0.);

  RooRealVar ap2corr("ap2corr","|a+|^{2}",0.);
  RooRealVar am2corr("am2corr","|a-|^{2}",0.);
  RooRealVar bp2corr("bp2corr","|b+|^{2}",0.);
  RooRealVar bm2corr("bm2corr","|b-|^{2}",0.);
  RooRealVar Pbcorr("Pbcorr","P_{b}",0.);
  RooRealVar abcorr("abcorr","#alpha_{b}",0.);
  RooRealVar r0corr("r0corr","r_{0}",0.);
  RooRealVar r1corr("r1corr","r_{1}",0.);

  const TString filename("datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetAll-sim.root");
  TFile file(filename);
  RooFitResult* fr = (RooFitResult*)file.Get("AngFit");


  TMatrixD M(4,4);
  for (int i=0;i<4;i++) for (int j=0;j<4;j++) M(i,j)=0.;
  
  M(0,2)=+0.5;
  M(0,3)=+0.5;
  M(1,2)=+0.5;
  M(1,3)=-0.5;
  
  M(2,1)=+0.5;
  M(2,2)=-0.5;
  M(2,3)=-0.5;
  
  M(3,1)=-0.5;
  M(3,2)=-0.5;
  M(3,3)=+0.5;
  
  TMatrixD MT(TMatrixD::kTransposed,M);

  TRandom3 rnd(0);

  RooDataSet parscorrected("parscorrected","parscorrected",RooArgList(Pbcorr,abcorr,r0corr,r1corr,ap2corr,am2corr,bp2corr,bm2corr));

  double Pbcorrnom=0.;
  double abcorrnom=0.;
  double r0corrnom=0.;
  double r1corrnom=0.;
  double ap2corrnom=0.;
  double am2corrnom=0.;
  double bp2corrnom=0.;
  double bm2corrnom=0.;

  const int maxtoy=500;

  TFile corrfile("toys/dataparam0/corr.root");

  for (int itoy=0;itoy<maxtoy;itoy++) {

    // TMatrixDSym Ccorr(4);
    // for (int i=0;i<4;i++) for (int j=0;j<4;j++) Ccorr(i,j)=0.;
    // Ccorr(0,0)=pow(itoy==0?1.070:rnd.Gaus(1.070,0.037),2.);
    // Ccorr(1,1)=pow(itoy==0?1.685:rnd.Gaus(1.685,0.059),2.);
    // Ccorr(2,2)=pow(itoy==0?1.037:rnd.Gaus(1.037,0.036),2.);
    // Ccorr(3,3)=pow(itoy==0?1.705:rnd.Gaus(1.705,0.080),2.);

    // Ccorr.Print("v");

    //    std::cout << corrfile.GetListOfKeys()->At(0)->GetName() << std::endl;
    //    std::cout << corrfile.GetListOfKeys()->At(1)->GetName() << std::endl;

    TVectorD meansvec(*(TVectorD*)corrfile.Get("meansvector"));
    TMatrixDSym Ccorr(*(TMatrixDSym*)corrfile.Get("covmatrix"));

    for (int i=0; i<4;i++) for (int j=0;j<4;j++) if (i!=j) Ccorr(i,j)=0.;
    

    TMatrixD Covcorr(fr->covarianceMatrix(),TMatrixD::kMult,Ccorr);
    //    TMatrixD& Covcorrsym = Covcorr;
    //    TMatrixD Covcorr(fr->covarianceMatrix(),TMatrixD::kMult,TMatrixD(*Ccorr,TMatrixD::kMult,fr->covarianceMatrix()));
    TMatrixDSym Covcorrsym = makesymm(Covcorr);

    //    TMatrixDSym Covcorrsym = getVCV(fr->covarianceMatrix(), *Ccorr);

    if (itoy==0) {
      fr->Print("v");
      meansvec.Print("v");

      std::cout << "Cov matrix" << std::endl;
      fr->covarianceMatrix().Print("v");

      std::cout << "Ccorr matrix" << std::endl;
      Ccorr.Print("v");

      std::cout << "Covcorr matrix" << std::endl;
      Covcorr.Print("v");

      std::cout << "CovCorrSym matrix" << std::endl;
      Covcorrsym.Print("v");

      TMatrixDSym* CM = computeCorrMatrix(Covcorrsym);
      std::cout << "Correlation matrix" << std::endl;
      CM->Print("v");
      mygetchar;
    }

    P_b=((RooRealVar*)fr->floatParsFinal().find("P_b"))->getVal();
    alpha_b=((RooRealVar*)fr->floatParsFinal().find("alpha_b"))->getVal();
    r_0=((RooRealVar*)fr->floatParsFinal().find("r_0"))->getVal();
    r_1=((RooRealVar*)fr->floatParsFinal().find("r_1"))->getVal();

    Pbcorr.setError(sqrt(Covcorrsym(0,0)));
    abcorr.setError(sqrt(Covcorrsym(1,1)));
    r0corr.setError(sqrt(Covcorrsym(2,2)));
    r1corr.setError(sqrt(Covcorrsym(3,3)));

    // TVectorD musq(4);
    // musq(0)=pow((itoy==0)?meansvec(0):rnd.Gaus(meansvec(0),0.053),2.);
    // musq(1)=pow((itoy==0)?meansvec(1):rnd.Gaus(meansvec(1),0.25),2.);
    // musq(2)=pow((itoy==0)?meansvec(2):rnd.Gaus(meansvec(2),0.051),2.);
    // musq(3)=pow((itoy==0)?meansvec(3):rnd.Gaus(meansvec(3),0.25),2.);

    
    
    const bool bug=false;
    if (bug) {
      Pbcorr=P_b.getVal()     + Pbcorr.getError() * ((itoy==0)?meansvec(0):rnd.Gaus(meansvec(0),0.053));
      abcorr=alpha_b.getVal() + abcorr.getError() * ((itoy==0)?meansvec(1):rnd.Gaus(meansvec(1),0.25));
      r0corr=r_0.getVal()     + r0corr.getError() * ((itoy==0)?meansvec(2):rnd.Gaus(meansvec(2),0.051));
      r1corr=r_1.getVal()     + r1corr.getError() * ((itoy==0)?meansvec(3):rnd.Gaus(meansvec(3),0.25));
    } else {
      Pbcorr=P_b.getVal()     - Pbcorr.getError() * ((itoy==0)?meansvec(0):rnd.Gaus(meansvec(0),0.053));
      abcorr=alpha_b.getVal() - abcorr.getError() * ((itoy==0)?meansvec(1):rnd.Gaus(meansvec(1),0.25));
      r0corr=r_0.getVal()     - r0corr.getError() * ((itoy==0)?meansvec(2):rnd.Gaus(meansvec(2),0.051));
      r1corr=r_1.getVal()     - r1corr.getError() * ((itoy==0)?meansvec(3):rnd.Gaus(meansvec(3),0.25));
    }
    //    Pbcorr=P_b.getVal()+ Pbcorr.getError() * ((itoy==0)?-0.0767:rnd.Gaus(-0.0767,0.053));
    //    abcorr=alpha_b.getVal()+ abcorr.getError() * ((itoy==0)?-0.2660:rnd.Gaus(-0.2660,0.084));
    //    r0corr=r_0.getVal()+ r0corr.getError() * ((itoy==0)?-0.1943:rnd.Gaus(-0.1943,0.051));
    //    r1corr=r_1.getVal()+r1corr.getError() * ((itoy==0)?-0.1682:rnd.Gaus(-0.1682,0.085));


    {
      TVectorD vec3(4);
      vec3(0)=P_b.getVal();
      vec3(1)=alpha_b.getVal();
      vec3(2)=r_0.getVal();
      vec3(3)=r_1.getVal();
      //	  vec(0)=P_b.getVal();
      //	  vec(1)=alpha_b.getVal();
      //	  vec(2)=r_0.getVal();
      //	  vec(3)=r_1.getVal();
      
      TVectorD vec4=M*vec3;
      
      ap2=vec4(0);
      am2=vec4(1);
      bp2=0.5+vec4(2);
      bm2=0.5+vec4(3);
    }


    TVectorD vec(4);
    vec(0)=Pbcorr.getVal();
    vec(1)=abcorr.getVal();
    vec(2)=r0corr.getVal();
    vec(3)=r1corr.getVal();
    //	  vec(0)=P_b.getVal();
    //	  vec(1)=alpha_b.getVal();
    //	  vec(2)=r_0.getVal();
    //	  vec(3)=r_1.getVal();
    
    TVectorD vec2=M*vec;
    
    ap2corr=vec2(0);
    am2corr=vec2(1);
    bp2corr=0.5+vec2(2);
    bm2corr=0.5+vec2(3);




    TMatrixD Covampl(M,TMatrixD::kMult,TMatrixD(Covcorr,TMatrixD::kMult,MT)) ;
    
    ap2corr.setError(sqrt(Covampl(0,0)));
    am2corr.setError(sqrt(Covampl(1,1)));
    bp2corr.setError(sqrt(Covampl(2,2)));
    bm2corr.setError(sqrt(Covampl(3,3)));	  
    
    Pbcorr=Pbcorr.getVal()-Pbcorrnom;
    abcorr=abcorr.getVal()-abcorrnom;
    r0corr=r0corr.getVal()-r0corrnom;
    r1corr=r1corr.getVal()-r1corrnom;
    ap2corr=ap2corr.getVal()-ap2corrnom;
    am2corr=am2corr.getVal()-am2corrnom;
    bp2corr=bp2corr.getVal()-bp2corrnom;
    bm2corr=bm2corr.getVal()-bm2corrnom;

    if (itoy==0) {
      printrrv(P_b);
      printrrv(alpha_b);
      printrrv(r_0);
      printrrv(r_1);
      printrrv(Pbcorr);
      printrrv(abcorr);
      printrrv(r0corr);
      printrrv(r1corr);
      printrrv(ap2corr);
      printrrv(am2corr);
      printrrv(bp2corr);
      printrrv(bm2corr);

      Pbcorrnom=Pbcorr.getVal();
      abcorrnom=abcorr.getVal();
      r0corrnom=r0corr.getVal();
      r1corrnom=r1corr.getVal();
      ap2corrnom=ap2corr.getVal();
      am2corrnom=am2corr.getVal();
      bp2corrnom=bp2corr.getVal();
      bm2corrnom=bm2corr.getVal();
      mygetchar;
    
    } else {
      parscorrected.add(RooArgSet(Pbcorr,abcorr,r0corr,r1corr,ap2corr,am2corr,bp2corr,bm2corr));
    }

  }
  plotdata(c, RooArgSet(Pbcorr,abcorr,r0corr,r1corr,ap2corr,am2corr,bp2corr,bm2corr), &parscorrected, "fitbias", "", true);

}


