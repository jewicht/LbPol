//g++ -o testroominuit testroominuit.C `root-config --cflags --libs` -lRooFit

#include <iostream>

#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"

#include "RooRealVar.h"
#include "RooLbtoJpsiL0PDF3.h"
#include "createRooLegendre3.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooAddition.h"


#include "rootlogon.h"

const double pi=atan(1.)*4.


int main(int argc, char** argv)
{

  rootlogon();
  
  TCanvas c("c","c",100,100);

  // RooRealVar theta("theta",  "#theta",    0.,pi);
  // RooRealVar theta1("theta1","#theta_{1}",0.,pi);
  // RooRealVar theta2("theta2","#theta_{2}",0.,pi);

  //  RooFormulaVar costheta("costheta",   "cos(theta)",  RooArgList(theta));
  //  RooFormulaVar costheta1("costheta1", "cos(theta1)", RooArgList(theta1));
  //  RooFormulaVar costheta2("costheta2", "cos(theta2)", RooArgList(theta2));

  RooRealVar costheta("costheta",   "cos(#theta)",     -1., 1.);//RooArgList(theta));
  RooRealVar costheta1("costheta1", "cos(#theta_{1})", -1., 1.);// RooArgList(theta1));
  RooRealVar costheta2("costheta2", "cos(#theta_{2})", -1., 1.);// RooArgList(theta2));

  RooRealVar phi1("phi1","#phi_{1}",-pi,pi);
  RooRealVar phi2("phi2","#phi_{2}",-pi,pi);

  RooRealVar keep("keep","keep",0.,1.);
  RooRealVar acceptance("w","w",0.,1.);

  RooRealVar lb_id("lb_id","lb_id",-6000.,6000.);

  RooRealVar P_b("P_b","P_{b}",1.00, -1.5, 1.5);//0.50
  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}",0.642,0.,1.);
  alpha_lambda.setConstant();

  RooRealVar alpha_b("alpha_b","#alpha_{b}",-0.45767,-1.1,1.1);//0.457
  RooRealVar r_0("r_0","r_{0}",0.251733,-0.3,2.3);//0.5
  RooRealVar r_1("r_1","r_{1}",0.116484,-1.3,1.3);//0.1

  RooRealVar* xi[8];
  for (int i=0; i<8; i++) xi[i] = new RooRealVar("xi" + istr(i),"xi" + istr(i),1.);

  RooLbtoJpsiL0PDF3 w3("w3","w3",costheta,costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);//,*xi[0],*xi[1],*xi[2],*xi[3],*xi[4],*xi[5],*xi[6],*xi[7]);


  //create acceptance

  RooLegendre3Pdfv2* accpdf = NULL;

  std::ifstream filestr;
  filestr.open("coeff3.txt");
  
  //build TF3
  double dcoeff3[ordermax3i][ordermax3j][ordermax3k];
  RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k];
  createRooLegendre3v2(costheta, costheta1, costheta2, "accpdf", accpdf, coeff3, true,4,4,7,10);
  accpdf->forceNumInt();
  for (int i=0; i<ordermax3i;i++) {
    for (int j=0; j<ordermax3j;j++) {
      for (int k=0; k<ordermax3k;k++) {
	double error;
	bool isconstant;
        filestr >>  dcoeff3[i][j][k] >> error >> isconstant;// >> std::endl;
	std::cout << dcoeff3[i][j][k] <<  " " << error << " " << isconstant << std::endl;
	if (coeff3[i][j][k]) {
	  coeff3[i][j][k]->setVal(dcoeff3[i][j][k]);
	  coeff3[i][j][k]->setError(error);
	  //	  coeff3[i][j][k]->setConstant(isconstant);
	  coeff3[i][j][k]->setConstant();
	}
      }
    }
  }
  filestr.close();


  int iSample=atoi(argv[1]);
  char buffer[50];
  sprintf( buffer, "%04d", iSample );
  
  RooDataSet* data = RooDataSet::read("toys/toy_" + TString(buffer) + ".dat", RooArgSet(costheta,costheta1,costheta2,phi1,phi2));
  
  RooAbsReal* nll_acc = accpdf->createNLL(*data);
  RooAbsReal* nll_model = w3.createNLL(*data);

  RooAddition nll("nll","nll",RooArgSet(*nll_acc,*nll_model));

  
  RooMinuit m(*nll_model);
  m.setVerbose(kTRUE) ;
  m.migrad();
}
