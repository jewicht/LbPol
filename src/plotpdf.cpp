#include "RooLbtoJpsiL0PDF3.h"

#include "functions-roofit.h"

int main(int argc, char** argv)
{
  RooRealVar costheta("costheta",   "cos#theta",     -1., 1.);
  RooRealVar costheta1("costheta1", "cos#theta_{1}", -1., 1.);
  RooRealVar costheta2("costheta2", "cos#theta_{2}", -1., 1.);

  lhcbstyle();
  //rootlogon();
  TCanvas c("c","c",100,100);
  RooPlot* plot=costheta.frame();
  plot->SetMinimum(0.);
  plot->Draw();
  c.SaveAs("aaa.eps");
  return 1;




  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}",0.642);

  RooRealVar Pb1("Pb1","P_{b}",-0.5);
  RooRealVar ab1("ab1","#alpha_{b}",-0.7);
  RooRealVar r01("r01","r_{0}",0.3);
  RooRealVar r1("r1","r_{1}", -0.5);

  RooRealVar Pb2("Pb2","P_{b}",0.);
  RooRealVar ab2("ab2","#alpha_{b}",-0.7);
  RooRealVar r02("r02","r_{0}",0.5);
  //  RooRealVar r12("r12","r_{1}", -0.5);

  RooRealVar Pb3("Pb3","P_{b}",0.5);
  RooRealVar ab3("ab3","#alpha_{b}",-0.7);
  RooRealVar r03("r03","r_{0}",0.5);
  //  RooRealVar r13("r13","r_{1}", -0.5);

  RooLbtoJpsiL0PDF3 sig_ang_model1("sig_ang_model1","sig_ang_model1",costheta,costheta1,costheta2, Pb1, ab1, alpha_lambda, r01, r1);
  RooLbtoJpsiL0PDF3 sig_ang_model2("sig_ang_model2","sig_ang_model2",costheta,costheta1,costheta2, Pb2, ab2, alpha_lambda, r02, r1);
  RooLbtoJpsiL0PDF3 sig_ang_model3("sig_ang_model3","sig_ang_model3",costheta,costheta1,costheta2, Pb3, ab3, alpha_lambda, r03, r1);



  plotpdflist(c,costheta,RooArgList(sig_ang_model1,sig_ang_model2,sig_ang_model3), "plotpdf-Pbvar-costheta0.eps");
  plotpdflist(c,costheta1,RooArgList(sig_ang_model1,sig_ang_model2,sig_ang_model3),"plotpdf-Pbvar-costheta1.eps");  
  plotpdflist(c,costheta2,RooArgList(sig_ang_model1,sig_ang_model2,sig_ang_model3),"plotpdf-Pbvar-costheta2.eps");

  Pb1=0.5;
  Pb2=0.5;
  Pb3=0.5;
  ab1=-0.5;
  ab2=0.;
  ab3=+0.5;

  plotpdflist(c,costheta,RooArgList(sig_ang_model1,sig_ang_model2,sig_ang_model3), "plotpdf-abvar-costheta0.eps");
  plotpdflist(c,costheta1,RooArgList(sig_ang_model1,sig_ang_model2,sig_ang_model3),"plotpdf-abvar-costheta1.eps");  
  plotpdflist(c,costheta2,RooArgList(sig_ang_model1,sig_ang_model2,sig_ang_model3),"plotpdf-abvar-costheta2.eps");

  Pb1=0.5;
  Pb2=0.5;
  Pb3=0.5;
  ab1=0.;
  ab2=0.;
  ab3=0.;

  r01=0.1;
  r02=0.;
  r03=-0.1;

  r1=0.4;

  plotpdflist(c,costheta,RooArgList(sig_ang_model1,sig_ang_model2,sig_ang_model3), "plotpdf-r0var-costheta0.eps");
  plotpdflist(c,costheta1,RooArgList(sig_ang_model1,sig_ang_model2,sig_ang_model3),"plotpdf-r0var-costheta1.eps");  
  plotpdflist(c,costheta2,RooArgList(sig_ang_model1,sig_ang_model2,sig_ang_model3),"plotpdf-r0var-costheta2.eps");
}
