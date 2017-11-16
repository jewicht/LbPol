//g++ -o test test.C `root-config --cflags --libs` -lRooFit

#include <iostream>

#include "TROOT.h"
#include "TComplex.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

#include "RooFit.h"
#include "Riostream.h"

#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"


#include "RooRealProxy.h"
#include "RooStats/SPlot.h"



//include "computew.h"
#include "RooLbtoJpsiL0PDF3.h"
#include "RooLbtoJpsiL0wAccPDF3.h"
#include "RooLbtoJpsiL0PDF5.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooDataHist.h"

#include "RooPlot.h"

#include "RooMCStudy.h"
#include "RooMinuit.h"

#include "RooCategory.h"
#include "RooEfficiency.h"

#include "createRooLegendre3.h"
#include "createRooLegendre5.h"


#include "calcacceptanceclass.h"
#include "functions-roofit.h"
#include "rootlogon.h"

const double pi=atan(1.)*4.

using namespace RooFit;


void polaire(const TComplex& c)
{
  //  std::cout << c << std::endl;
  std::cout << c.Rho() << " * exp(i * " << c.Theta() << " )"<< std::endl;
  //  std::cout << c.Rho() * TComplex::Exp(TComplex::I()*c.Theta())<< std::endl;

}

void test3()
{
  TComplex ap = -0.0176 - 0.4290 * TComplex::I();
  TComplex am =  0.0867 + 0.2454 * TComplex::I();
  TComplex bp = -0.0810 - 0.2837 * TComplex::I();
  TComplex bm =  0.0296 + 0.8124 * TComplex::I();
  

  std::cout << ap.Rho2() + am.Rho2() + bp.Rho2() + bm.Rho2() << std::endl;
  std::cout << "alpha_b = " <<  ap.Rho2() - am.Rho2() + bp.Rho2() - bm.Rho2() << std::endl;
  polaire(ap);
  polaire(am);
  polaire(bp);
  polaire(bm);

  std::cout << "sum theta = " << ap.Theta()+am.Theta()+bp.Theta()+bm.Theta() << std::endl;
}

void test4()
{
  TComplex a=0.813 + 1.534 * TComplex::I();
  TComplex b=0.429 - 1.612 * TComplex::I();
  TComplex c=0.260 + 1.231 * TComplex::I();
  TComplex d=0.295 - 1.849 * TComplex::I();
  polaire(a);
  polaire(b);
  polaire(c);
  polaire(d);

  std::cout << "sum theta = " << a.Theta()+b.Theta()+c.Theta()+d.Theta() << std::endl;

  std::cout << a << std::endl;
  std::cout << a.Rho() * TComplex::Exp( a.Theta() * TComplex::I()) << std::endl;

}

void addacceptanceweight(RooDataSet* data, RooRealVar& w)
{
  RooDataSet dataww6("dataww6","dataww6",RooArgList(w));
  
  calcacceptanceclass mycalcacceptanceclass("coeff3.txt");
  
  for (int i=0;i<data->sumEntries();i++) {
    const RooArgSet* blu = data->get(i);
    
    w.setVal(1./mycalcacceptanceclass.evaluate(blu->getRealValue("costheta"),
					       blu->getRealValue("costheta1"),
					       blu->getRealValue("costheta2")));
    //       std::cout << "cos0 = " << costheta.getVal() << std::endl;
    //       std::cout << "cos1 = " << costheta1.getVal() << std::endl;
    //       std::cout << "cos2 = " << costheta2.getVal() << std::endl;
    //    std::cout << "w    = " << w.getVal() << std::endl;
    dataww6.add(RooArgSet(w));
    //    const RooArgSet blu2(costheta,costheta1,costheta2,phi1,phi2);
    //    blu2.Print("v");
    //    std::cout << "accpdf3              = " << accpdf3->getVal() << std::endl;
    //    std::cout << "data_acc5_pdf        = " << data_acc5_pdf.getVal() << std::endl;
    //    std::cout << "data_acc5_pdf (blu)  = " << data_acc5_pdf.getVal(blu) << std::endl;
    //    std::cout << "data_acc5_pdf (blu2) = " << data_acc5_pdf.getVal(&blu2) << std::endl;
  }
  data->merge(&dataww6);
}

void calculatew1w2(RooDataSet* data, const TString& w1str, const TString& w2str, RooRealVar& w)
{
    RooDataSet dataww6("dataww6","dataww6",RooArgList(w));
    
    Float_t w1,w2;

    for (int i=0;i<data->sumEntries();i++) {
      const RooArgSet* blu = data->get(i);
      //    blu->Print("v");
      //	dcostheta=costheta.getVal();
      //	dcostheta1=costheta1.getVal();
      //	dcostheta2=costheta2.getVal();
      w1=blu->getRealValue(w1str);
      w2=blu->getRealValue(w2str);
      w.setVal( w1*w2);
      //      std::cout << "w1   = " << w1 << std::endl;
      //      std::cout << "w2   = " << w2 << std::endl;
      //      std::cout << "w    = " << w.getVal() << std::endl;
      dataww6.add(RooArgSet(w));
      //    const RooArgSet blu2(costheta,costheta1,costheta2,phi1,phi2);
      //    blu2.Print("v");
      //    std::cout << "accpdf3              = " << accpdf3->getVal() << std::endl;
      //    std::cout << "data_acc5_pdf        = " << data_acc5_pdf.getVal() << std::endl;
    //    std::cout << "data_acc5_pdf (blu)  = " << data_acc5_pdf.getVal(blu) << std::endl;
    //    std::cout << "data_acc5_pdf (blu2) = " << data_acc5_pdf.getVal(&blu2) << std::endl;
    }
    data->merge(&dataww6);
}


void usage(char* argv0)
{
  std::cout << argv0 << " gen nVariables nSamples nEventspersample" << std::endl;
  std::cout << argv0 << " fit nVariables iSample" << std::endl;
  std::cout << argv0 << " fitwmass nVariables iSample" << std::endl;
  std::cout << argv0 << " merge nVariables nSamples" << std::endl;
  std::cout << argv0 << " all nVariables nSamples nEventspersample" << std::endl;
  exit(1);
}

int main(int argc, char** argv)
{

  //  test3();exit(0);
  //  test1();return 0;
  //  test3();test4();return 0;

  if (argc<1) usage(argv[0]);

  bool dogen=false;
  bool domerge=false;
  bool doall=false;
  bool dofit=false;
  bool dofitmc=false;
  bool dotest=false;
  bool dofitwithmass=false;
  
  
  int nSamples=0; 
  int nEventspersample=0;
  int iSample=0;


  const TString opmode(argv[1]);
  TString filename;
  if (opmode=="gen") {
    if (argc!=5) usage(argv[0]);
    dogen=true;
    nSamples=atoi(argv[3]);
    nEventspersample=atoi(argv[4]);
  } else if (opmode=="fit") {
    if (argc!=4) usage(argv[0]);
    dofit=true;
    iSample=atoi(argv[3]);
  } else if (opmode=="fitwmass") {
    if (argc!=4) usage(argv[0]);
    dofitwithmass=true;
    iSample=atoi(argv[3]);
  } else if (opmode=="merge") {
    if (argc!=4) usage(argv[0]);
    domerge=true;
    nSamples=atoi(argv[3]);
  } else if (opmode=="all") {
    if (argc!=5) usage(argv[0]);
    doall=true;
    nSamples=atoi(argv[3]);
    nEventspersample=atoi(argv[4]);    
  } else if (opmode=="fitmc") {
    dofitmc=true;
    filename=argv[3];
  } else if (opmode=="test") {
    dotest=true;
  } else {
    usage(argv[0]);
  }
  const int nVariables=atoi(argv[2]);

  rootlogon();
  
  TCanvas c("c","c",100,100);

  // RooRealVar theta("theta",  "#theta",    0.,pi);
  // RooRealVar theta1("theta1","#theta_{1}",0.,pi);
  // RooRealVar theta2("theta2","#theta_{2}",0.,pi);

  //  RooFormulaVar costheta("costheta",   "cos(theta)",  RooArgList(theta));
  //  RooFormulaVar costheta1("costheta1", "cos(theta1)", RooArgList(theta1));
  //  RooFormulaVar costheta2("costheta2", "cos(theta2)", RooArgList(theta2));

  RooRealVar mass("mass","M(#LambdaJ/#psi)",5500.,5750.,"MeV/c^{2}");

  RooRealVar costheta("costheta",   "cos(#theta)",     -1., 1.);//RooArgList(theta));
  RooRealVar costheta1("costheta1", "cos(#theta_{1})", -1., 1.);// RooArgList(theta1));
  RooRealVar costheta2("costheta2", "cos(#theta_{2})", -1., 1.);// RooArgList(theta2));

  RooRealVar phi1("phi1","#phi_{1}/#pi",-1.,1.);
  RooRealVar phi2("phi2","#phi_{2}/#pi",-1.,1.);

  RooRealVar keep("keep","keep",0.,1.);
  //  RooRealVar acceptance("w","w",0.,1.);

  RooRealVar wacc("wacc","wacc",0.,5000.);
  RooRealVar w("w","w",-100000.,100000.);

  //  RooRealVar w2mu("w2mu","w2mu",1.,100000000.);//,10000.);

  RooRealVar lb_id("lb_id","lb_id",-6000.,6000.);

  RooRealVar P_b("P_b","P_{b}",0.40, -1.5, 1.5);//0.50
  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}",0.642, 0., 1.);
  //  alpha_lambda.setConstant();
  RooGaussian alpha_lambda_constrain("alpha_lambda_constrain","alpha_lambda_constrain",alpha_lambda,RooConst(0.642),RooConst(0.013));
  
  
  RooRealVar alpha_b("alpha_b","#alpha_{b}",-0.45767,-1.1,1.1);//0.457
  RooRealVar r_0("r_0","r_{0}",0.251733,-0.3,2.3);//0.5
  RooRealVar r_1("r_1","r_{1}",0.116484,-1.3,1.3);//0.1


//   P_b=0.;
//   alpha_b=0.;
//   r_0=0.5;
//   r_1=0.;

  RooRealVar f1("f1","f1",0.5,-0.1,1.1);
  RooRealVar f2("f2","f2",0.5,-0.1,1.1);
  RooRealVar f3("f3","f3",0.5,-0.1,1.1);

  RooFormulaVar aplus("aplus","f1",RooArgList(f1));
  RooFormulaVar aminus("aminus","(1.-f1)*f2",RooArgList(f1,f2));
  RooFormulaVar bplus("bplus","(1.-f1)*(1.-f2)*f3",RooArgList(f1,f2,f3));
  RooFormulaVar bminus("bminus","(1.-f1)*(1.-f2)*(1.-f3)",RooArgList(f1,f2,f3));

  RooFormulaVar alpha_b_form("alpha_b_form","aplus*aplus-aminus*aminus+bplus*bplus-bminus*bminus",RooArgList(aplus,aminus,bplus,bminus));
  RooFormulaVar r_0_form("r_0_form","aplus*aplus+aminus*aminus",RooArgList(aplus,aminus));
  RooFormulaVar r_1_form("r_1_form","aplus*aplus-aminus*aminus",RooArgList(aplus,aminus));
  
  RooRealVar alpha_plus("alpha_plus",  "#alpha_{+}", 0.237, -0.5*pi,2.5*pi);
  RooRealVar alpha_minus("alpha_minus","#alpha_{-}", 3.08,  -0.5*pi,2.5*pi);
  RooRealVar chi("chi","#chi",6.21719,0.,5.*pi);

  RooRealVar frac("frac","frac", 0.46, 0., 1.);

  RooRealVar mean1("mean1","mean1",5619.4,5500.,5700.);
  RooRealVar sigma1("sigma1","sigma1",12.440,0.,100.);
  RooGaussian gauss1("gauss1","gauss1",mass,mean1,sigma1);
  
  //  RooRealVar mean2("mean2","mean2",5620.,5500.,5700.);
  RooRealVar sigma2("sigma2","sigma2",5.8403,0.,100.);
  RooGaussian gauss2("gauss2","gauss2",mass,mean1,sigma2);
  
  RooAddPdf sigmassmodel("sigmassmodel","sigmassmodel",RooArgList(gauss1,gauss2),RooArgList(frac));

  RooRealVar c0("c0","c0",-0.13,-1.1,1.1);
  //  c0.setConstant();
  RooRealVar c1("c1","c1",-0.0415,-1.1,1.1);
  
  RooChebychev bkgmassmodel("bkgmassmodel","bkgmassmodel",mass,RooArgList(c0,c1));

  RooRealVar c2("c2","c2",0.,-1.1,1.1);
  RooRealVar c3("c3","c3",0.,-1.1,1.1);
  RooRealVar c4("c4","c4",0.,-1.1,1.1);

  RooChebychev polc2("polc2","polc2",costheta,RooArgList(c2));
  RooChebychev polc3("polc3","polc3",costheta1,RooArgList(c3));
  RooChebychev polc4("polc4","polc4",costheta2,RooArgList(c4));

  RooProdPdf bkgangmodel("bkgangmodel","bkgangmodel",RooArgList(polc2,polc3,polc4));

  RooProdPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(bkgmassmodel,bkgangmodel));

  RooRealVar nSig("nSig","nSig",2000.,0.,10000.);
  RooRealVar nBkg("nBkg","nBkg",3000.,0.,8000.);
  RooAddPdf massmodel("massmodel","massmodel",RooArgList(sigmassmodel,bkgmassmodel),RooArgList(nSig,nBkg));

  //  RooMCStudy massmodelstudy(massmodel,RooArgList(mass),Extended(true));
  //  massmodelstudy.generateAndFit(300, 10000) ;
  //  plot_mcstudy(massmodelstudy,c,RooArgList(nSig,nBkg,frac,mean1,sigma1,sigma2,c0,c1),99);
  

  // RooRealVar* xi[8];
  // double xival[8];
  // xival[0] = 278215;
  // xival[1] = 2107.85;
  // xival[2] = -70508.6;
  // xival[3] = -1.56197;
  // xival[4] = 1946.35;
  // xival[5] = 132.977;
  // xival[6] = -101.589;
  // xival[7] = -126.39;
  // for (int i=0; i<8; i++) xival[i] = 1.;
  // for (int i=0; i<8; i++) xi[i] = new RooRealVar("xi" + istr(i),"xi" + istr(i),1.);


  //create acceptance

  //  RooLegendre3Pdf* accpdf = NULL;
  RooLbtoJpsiL0wAccPDF3* w3wacc = NULL;

  //build TF3

  //  double dcoeff3[ordermax3i][ordermax3j][ordermax3k];
  RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k];
  RooRealVar* coeff5[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m];
  RooRealVar* coeff32[ordermax3i][ordermax3j][ordermax3k];
  //  createRooLegendre3(costheta, costheta1, costheta2, "accpdf", accpdf, coeff3);
  //  createRooLegendre3(costheta, costheta1, costheta2, "accpdf", accpdf, coeff3);
  createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b, alpha_lambda, r_0, r_1, "w3wacc", w3wacc, coeff3, true, "coeff3.txt");
  //  readcoeff3_rrv("coeff3.txt",coeff3);

  RooAbsPdf* accpdf3 = NULL;
  RooLegendre5Pdfv2* accpdf5 = NULL;
  createRooLegendre3v2(costheta, costheta1, costheta2,             "accpdf3", accpdf3, coeff32, true, 4,4,7,15,5, "coeff3.txt");
  createRooLegendre5v2(costheta, costheta1, costheta2, phi1, phi2, "accpdf5", accpdf5, coeff5, 4,4,4,4,7,15,2, "coeff5-nofit.txt");
  //  readcoeff3_rrv("coeff3.txt",coeff32);

  //  std::cerr << __FILE__ << ":" << __LINE__ << std::endl;

  //RooLbtoJpsiL0PDF3 w3("w3","w3",costheta,costheta1,costheta2,P_b,alpha_b_form,alpha_lambda,r_0_form,r_1_form);
  RooLbtoJpsiL0PDF3 w3("w3","w3",costheta,costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);//,*xi[0],*xi[1],*xi[2],*xi[3],*xi[4],*xi[5],*xi[6],*xi[7]);

  //  RooLbtoJpsiL0PDF5 w5("w5","w5",costheta,costheta1,costheta2,phi1,phi2,P_b,alpha_b_form,alpha_lambda,r_0_form,r_1_form,alpha_plus,alpha_minus,chi);
  RooLbtoJpsiL0PDF5 w5("w5","w5",costheta,costheta1,costheta2,phi1,phi2,P_b,alpha_b,alpha_lambda,r_0,r_1,alpha_plus,alpha_minus,chi);

  RooDataSet* testdata = w5.generate(RooArgList(costheta,costheta1,costheta2,phi1,phi2), 10000);
  RooDataSet* testdata2 = (RooDataSet*)testdata->reduce("phi1>0.");
  plotdatavsdata(c,RooArgList(costheta,costheta1,costheta2),testdata,testdata2,"phi1cut");
  return 1;

  RooCategory cut("cut","cutr");
  cut.defineType("accept",1);
  cut.defineType("reject",0);
  
  //  RooFormulaVar rfv("rfv","rfv","costheta*costheta*costheta1*costheta1*costheta2*costheta2",RooArgList(costheta,costheta1,costheta2));
  RooEfficiency efficiency3("eff3","eff3",*accpdf3,cut,"accept");
  RooProdPdf w3xacc("w3xacc","w3xacc",w3,Conditional(efficiency3,cut));

  //  plotpdflist(c,costheta,RooArgList(w3xacc,*w3wacc),"comp-pdf-costheta.eps");
  //  plotpdflist(c,costheta1,RooArgList(w3xacc,*w3wacc),"comp-pdf-costheta1.eps");
  //  plotpdflist(c,costheta2,RooArgList(w3xacc,*w3wacc),"comp-pdf-costheta2.eps");
  //  return 0;
  
  if (false) {

    RooDataSet* dataww = w3xacc.generate(RooArgSet(costheta,costheta1,costheta2,cut),100000);
    dataww->write("dataww.dat");
    plotdata(c, costheta,  dataww, "dataww-costheta.eps" , "cut==1");
    plotdata(c, costheta1, dataww, "dataww-costheta1.eps", "cut==1");
    plotdata(c, costheta2, dataww, "dataww-costheta2.eps", "cut==1");
    
    costheta.setBins(4);
    costheta1.setBins(4);
    costheta2.setBins(4);
    phi1.setBins(4);
    phi2.setBins(10);
    
    RooDataSet* data_acc5 = RooDataSet::read("acc_histo5.dat", RooArgList(costheta,costheta1,costheta2,phi1,phi2));
    RooDataHist data_acc5_binned("data_acc5_binned","data_acc5_binned",RooArgList(costheta,costheta1,costheta2,phi1,phi2),*data_acc5);
    RooHistPdf data_acc5_pdf("data_acc5_pdf","data_acc5_pdf",RooArgList(costheta,costheta1,costheta2,phi1,phi2),data_acc5_binned);
    
    RooEfficiency efficiency5("eff5","eff5",data_acc5_pdf,cut,"accept");
    
    RooProdPdf w5wacc("w5wacc","w5wacc",w5,Conditional(efficiency5,cut));
    RooDataSet* dataww5 = w5wacc.generate(RooArgSet(costheta,costheta1,costheta2,phi1,phi2,cut),100000);

    

    RooDataSet dataww6("dataww6","dataww6",RooArgList(w));
    
    
    for (int i=0;i<dataww5->sumEntries();i++) {
      const RooArgSet* blu = dataww5->get(i);
      //    blu->Print("v");
      //	dcostheta=costheta.getVal();
      //	dcostheta1=costheta1.getVal();
      //	dcostheta2=costheta2.getVal();
      costheta.setVal(blu->getRealValue("costheta"));
      costheta1.setVal(blu->getRealValue("costheta1"));
      costheta2.setVal(blu->getRealValue("costheta2"));
      phi1.setVal(blu->getRealValue("phi1"));
      phi2.setVal(blu->getRealValue("phi2"));
      w.setVal(data_acc5_pdf.getVal());
      dataww6.add(RooArgSet(w));
      //    const RooArgSet blu2(costheta,costheta1,costheta2,phi1,phi2);
      //    blu2.Print("v");
      //    std::cout << "accpdf3              = " << accpdf3->getVal() << std::endl;
      //    std::cout << "data_acc5_pdf        = " << data_acc5_pdf.getVal() << std::endl;
    //    std::cout << "data_acc5_pdf (blu)  = " << data_acc5_pdf.getVal(blu) << std::endl;
    //    std::cout << "data_acc5_pdf (blu2) = " << data_acc5_pdf.getVal(&blu2) << std::endl;
    }
    dataww5->merge(&dataww6);
    
    RooRealVar* rrv = (RooRealVar*)dataww5->addColumn(data_acc5_pdf);
    RooFormulaVar rfv("rfv","rfv","1./"+TString(rrv->GetName()),RooArgList(*rrv));
    
    dataww5->addColumn(rfv);
    
    dataww5->Print("v");
    dataww5->write("dataww5.dat");
    plotdata(c, costheta,  dataww5, "dataww5-costheta.eps" , "cut==1");
    plotdata(c, costheta1, dataww5, "dataww5-costheta1.eps", "cut==1");
    plotdata(c, costheta2, dataww5, "dataww5-costheta2.eps", "cut==1");
    plotdata(c, phi1,      dataww5, "dataww5-phi1.eps",      "cut==1");
    plotdata(c, phi2,      dataww5, "dataww5-phi2.eps",      "cut==1");
    
    return 1;
  }


  // TRandom3 rnd;
  // for (int i=0;i < 10000; i++) {
  //   costheta.setVal(-1.+2.*rnd.Rndm());
  //   costheta1.setVal(-1.+2.*rnd.Rndm());
  //   costheta2.setVal(-1.+2.*rnd.Rndm());
  //   std::cout << w3xacc->getVal() << std::endl;
  //   //    std::cout << w3.getVal() << std::endl;

  // }


  // P_b.setVal(0.);
  // alpha_b.setVal(0.);

  // plotpdf(c,costheta,w3, "w3-acceptancetest-costheta-before.eps");
  // plotpdf(c,costheta1,w3,"w3-acceptancetest-costheta1-before.eps");
  // plotpdf(c,costheta2,w3,"w3-acceptancetest-costheta2-before.eps");

  // for (int i=0; i<8; i++) xi[i]->setVal(xival[i]);

  // plotpdf(c,costheta, w3,"w3-acceptancetest-costheta-after.eps");
  // plotpdf(c,costheta1,w3,"w3-acceptancetest-costheta1-after.eps");
  // plotpdf(c,costheta2,w3,"w3-acceptancetest-costheta2-after.eps");

  // //  return 0;




  RooAbsPdf* wn;
  RooArgSet* vars;
  if (nVariables==3) {
    //    vars = new RooArgSet(costheta,costheta1,costheta2,phi1,phi2,w);
    vars = new RooArgSet(costheta,costheta1,costheta2);
    wn=&w3;
    //    wn=w3wacc;
    //    w=w3wacc;
  } else {
    vars = new RooArgSet(costheta,costheta1,costheta2,phi1,phi2);
    wn=&w5;
  }

  //  if (dogen) w.forceNumInt();

  //  RooMCStudy mgr(w,RooArgSet(theta,theta_1,theta_2,phi_1,phi_2)) ;

  P_b=0.4;
  alpha_b=-0.45767;
  r_0=0.251733;
  r_1=0.116484;
  RooMCStudy mgr(*wn,*vars) ;

  
  // Setup PDF
  //  RooRealVar x("x","x",-5,15) ;
  //  RooRealVar mean("mean","mean of gaussian",-1.,-3.,3.) ;
  //  RooRealVar sigma("sigma","width of gaussian",4.,2.,6.) ;
  //  RooGaussian gauss("gauss","gaussian PDF",x,mean,sigma) ;
  //  RooMCStudy mgr(gauss,RooArgSet(x));
  // Create manager


  // Generate and fit 1000 experiments of 100 events each
  //  mgr.generateAndFit(1000,100) ;

  if (dogen) {
    if (nSamples>1) {
      mgr.generate(nSamples, nEventspersample, kTRUE, "toys/toy_%04d.dat") ;
    } else { 
      //      RooDataSet* data = bkgangmodel.generate(*vars,nEventspersample);
      RooDataSet* data = wn->generate(*vars,nEventspersample);
      
      const TString output("sample.root");
      TFile *f = new TFile(output,"recreate");
      f->mkdir("BtoqllGeneratedDistributions");
      f->cd("BtoqllGeneratedDistributions");
      //      TDirectory *direc = new TDirectory("BtoqllGeneratedDistributions","BtoqllGeneratedDistributions");
      TTree *T = new TTree("Lb2JpsiL0","test friend trees");
      //      direc->Add(T);
      Int_t lb_id;
      Float_t dcostheta,dcostheta1,dcostheta2;
      T->Branch("lb_id",    &lb_id,    "lb_id/I"); 
      T->Branch("costheta",    &dcostheta,    "costheta/F"); 
      T->Branch("costheta1",   &dcostheta1,   "costheta1/F"); 
      T->Branch("costheta2",   &dcostheta2,   "costheta2/F"); 
      //      T->Branch("keep", &keep, "keep/I"); 
      //      T->Branch("w",    &w,    "w/F"); 

      for (int i=0;i<data->sumEntries();i++) {
	const RooArgSet* blu = data->get(i);
	//	dcostheta=costheta.getVal();
	//	dcostheta1=costheta1.getVal();
	//	dcostheta2=costheta2.getVal();
	dcostheta=blu->getRealValue("costheta");
	dcostheta1=blu->getRealValue("costheta1");
	dcostheta2=blu->getRealValue("costheta2");
	lb_id=5122;
	T->Fill();
      }
      
      T->Write();
      f->Close();
      //      data->write("sample.dat");
    }
    return 0;
  }

  if (domerge) {
    for (int i=0; i<nSamples;i++) {
      TFile* frfile = new TFile("toy_" + istr(i,"04") + "_" + istr(nVariables) + "variables_result.root");
      if (frfile and !frfile->IsZombie()) {
	//        RooFitResult* fr = (RooFitResult*)frfile->Get("fitresult_w"+istr(nVariables)+"_databinned");
	//        RooFitResult* fr = (RooFitResult*)frfile->Get("fitresult_w"+istr(nVariables)+"_data");
        const TString frname = frfile->GetListOfKeys()->At(0)->GetName();
        RooFitResult* fr = (RooFitResult*)frfile->Get(frname);
	//	std::cout << fr->GetName() << std::endl;
	if (fr) {
	  if (fr->status()==0 and fr->covQual()>-3)  {
	    std::cerr << "CONVERGED     for toy_" + istr(i,"04") + "_" + istr(nVariables) + "variables_result.root (status=" << fr->status() << " : covQual=" << fr->covQual() << ")" << std::endl;
	    mgr.addFitResult(*fr);
	  } else {
	    std::cerr << "NOT CONVERGED for toy_" + istr(i,"04") + "_" + istr(nVariables) + "variables_result.root (status=" << fr->status() << " : covQual=" << fr->covQual() << ")" << std::endl;
	  }
	} else {
	  std::cerr << "fr does not exist in toy_" + istr(i,"04") + "_" + istr(nVariables) + "variables_result.root" << std::endl;
	}
      }
      
      frfile->Close();
    }
    mygetchar;
    plot_mcstudy(c,mgr,RooArgList(P_b,alpha_b,r_0,r_1),nVariables);

    if (nVariables==5) {
      plot_mcstudy(c,mgr,RooArgList(alpha_plus,alpha_minus,chi),nVariables);
    }
    return 0;
  }

  if (doall) {
    mgr.generateAndFit(nSamples, nEventspersample) ;
    plot_mcstudy(c,mgr,RooArgList(P_b,alpha_b,r_0,r_1),nVariables);
    if (nVariables==5) {
      plot_mcstudy(c,mgr,RooArgList(alpha_plus,alpha_minus,chi),nVariables);
    }
  }


  if (dofit) {

    const bool weightedfit=false;

    RooDataSet* dataunbinned = RooDataSet::read("toys/toy_" + istr(iSample,"04") + ".dat", RooArgList(costheta,costheta1,costheta2));
    
    // TH3F datahisto("datahisto","datahisto",
    // 		   10,-1.,1.,
    // 		   10,-1.,1.,
    // 		   10,-1.,1.);
    // dataunbinned->fillHistogram(&datahisto,RooArgList(costheta,costheta1,costheta2));
    
    // TFile accfile("acc_histo.root");
    // TH3F* acc_histo = (TH3F*)accfile.Get("acc_histo");
    // datahisto.Divide(acc_histo);
    // RooDataHist databinned("databinned","databinned",RooArgList(costheta,costheta1,costheta2),&datahisto);
   
    //    vars->add(w);
    //    RooDataSet* data = new RooDataSet("data","data",bkgdata,RooArgList(costheta,costheta1,costheta2,w),0,w.GetName());
    //    RooDataSet* data = new RooDataSet("data","data",bkgdata,RooArgList(costheta,costheta1,costheta2,w),0,"nSig_sw");

    //    RooDataSet* data = dataunbinned;
    RooDataSet* data = dataunbinned;
    if (weightedfit) data = new RooDataSet("data","data",dataunbinned,*dataunbinned->get(),0,wacc.GetName());


    RooAbsPdf* fitfunction = w3wacc;

    data->Print("v");
    mygetchar;
    //    RooAbsData* data = dataunbinned;

    
    
    //    RooFitResult* frorig = w->fitTo(*data, Save(true), SumW2Error(false), Minos(true));

    RooFitResult* fr = fitfunction->fitTo(*data, Save(true), SumW2Error(weightedfit));//, Minos(true));
    //    RooFitResult* fr = w->fitTo(*data, Save(true));
    //    frorig->Print("v");
    //    std::cout << "Fit result: status=" << frorig->status() << " : covQual=" << frorig->covQual() << std::endl;
    //    mygetchar;
    fr->Print("v");
    std::cout << "Fit result: status=" << fr->status() << " : covQual=" << fr->covQual() << std::endl;
    //    mygetchar;
    //    mgr.fit(iSample, "output/toy_%04d.dat");
    //    const RooFitResult* fr = mgr.fitResult();//iSample);

    plot(c, RooArgList(costheta,costheta1,costheta2),  *fitfunction, data, "toy_" + istr(iSample,"04") + "_" + istr(nVariables) + "variables");
    if (nVariables==5) {
      plot(c, RooArgList(phi1,phi2),      *fitfunction, data, "toy_" + istr(iSample,"04") + "_" + istr(nVariables) + "variables");
    }

    fr->SaveAs(TString("toy_" + istr(iSample,"04") + "_" + istr(nVariables) + "variables_result.root"));
    return 0;
  }
  
  if (dofitwithmass) {
    RooDataSet* dataunbinned = RooDataSet::read("toys/toy_" + istr(iSample,"04") + ".dat", RooArgList(costheta,costheta1,costheta2));//,wacc));
    
    //add mass for signal
    RooDataSet* sigdata = sigmassmodel.generate(RooArgSet(mass), *dataunbinned, dataunbinned->sumEntries());
    
    //add background
    RooDataSet* bkgdata = bkgmodel.generate(RooArgSet(mass,costheta,costheta1,costheta2),10000);
    addacceptanceweight(bkgdata, wacc);
    
    //    RooCategory* blindState = new RooCategory("blindState","Blinding State") ;
    //    bkgdata->addColumn(*blindState) ;

    sigdata->Print("v");
    bkgdata->Print("v");

    bkgdata->append(*sigdata);
    //    sigdata->merge(bkgdata);
    mygetchar;
    RooFitResult* massfr = massmodel.fitTo(*bkgdata,Extended(true),Save(true));

    mass.setBins(100);
    plot(c, mass, massmodel, bkgdata,  "toy_" + istr(iSample,"04") + "_" + istr(nVariables) + "variables_mass.eps",true);


    massfr->Print("v");

    frac.setConstant();
    c0.setConstant();
    c1.setConstant();
    mean1.setConstant();
    sigma1.setConstant();
    sigma2.setConstant();
    RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
					       *bkgdata, &massmodel, RooArgList(nSig,nBkg) );

    //    bkgdata->Print("v");mygetchar;
    //    calculatew1w2(bkgdata, "nSig_sw", "wacc", w);
    //    bkgdata->Print("v");mygetchar;

    //    RooRealVar* nSig_sw
    

    //    dataunbinned->Print("v");
    //    calculateweight3(dataunbinned, w, *accpdf3, costheta, costheta1, costheta2);
    //    dataunbinned->Print("v");
    //    RooRealVar* rrv = (RooRealVar*)dataunbinned->addColumn(*accpdf3);
    //    RooFormulaVar rfv("rfv","rfv","1./"+TString(rrv->GetName()),RooArgList(*rrv));
    //    dataunbinned->addColumn(rfv);

    // TH3F datahisto("datahisto","datahisto",
    // 		   10,-1.,1.,
    // 		   10,-1.,1.,
    // 		   10,-1.,1.);
    // dataunbinned->fillHistogram(&datahisto,RooArgList(costheta,costheta1,costheta2));
    
    // TFile accfile("acc_histo.root");
    // TH3F* acc_histo = (TH3F*)accfile.Get("acc_histo");
    // datahisto.Divide(acc_histo);
    // RooDataHist databinned("databinned","databinned",RooArgList(costheta,costheta1,costheta2),&datahisto);
   
    //    vars->add(w);
    //    RooDataSet* data = new RooDataSet("data","data",bkgdata,RooArgList(costheta,costheta1,costheta2,w),0,w.GetName());
    //    RooDataSet* data = new RooDataSet("data","data",bkgdata,RooArgList(costheta,costheta1,costheta2,w),0,"nSig_sw");

    RooDataSet* data = new RooDataSet("data","data",bkgdata,*bkgdata->get(),0,"nSig_sw");

    RooAbsPdf* fitfunction = w3wacc;

    data->Print("v");
    mygetchar;
    //    RooAbsData* data = dataunbinned;

    
    
    //    RooFitResult* frorig = w->fitTo(*data, Save(true), SumW2Error(false), Minos(true));

    RooFitResult* fr = fitfunction->fitTo(*data, Save(true), SumW2Error(true));//, Minos(true));
    //    RooFitResult* fr = w->fitTo(*data, Save(true));
    //    frorig->Print("v");
    //    std::cout << "Fit result: status=" << frorig->status() << " : covQual=" << frorig->covQual() << std::endl;
    //    mygetchar;
    fr->Print("v");
    std::cout << "Fit result: status=" << fr->status() << " : covQual=" << fr->covQual() << std::endl;
    //    mygetchar;
    //    mgr.fit(iSample, "output/toy_%04d.dat");
    //    const RooFitResult* fr = mgr.fitResult();//iSample);

    costheta.setBins(25);
    costheta1.setBins(25);
    costheta2.setBins(25);

    plot(c, RooArgList(costheta,costheta1,costheta2),  *fitfunction, data, "toy_" + istr(iSample,"04") + "_" + istr(nVariables) + "variables");
    if (nVariables==5) {
      plot(c, RooArgList(phi1,phi2),      *fitfunction, data, "toy_" + istr(iSample,"04") + "_" + istr(nVariables) + "variables");
    }

    fr->SaveAs(TString("toy_" + istr(iSample,"04") + "_" + istr(nVariables) + "variables_result.root"));
    return 0;
  }



  // if (dofitmc) {
  //   std::cout << __LINE__ << std::endl;
    
  //   TFile myfile(filename);
  //   TDirectory* mydirectory = (TDirectory*)myfile.Get("BtoqllGeneratedDistributions");
  //   TTree* mytree = (TTree*)mydirectory->Get("Lb2JpsiL0");
    
  //   TString friendfile(filename);
  //   friendfile.ReplaceAll(".root",".acceptance.root");
  //   mytree->AddFriend("Lb2JpsiL0",friendfile);


  //   //    RooFormulaVar unmw("unmw","unmw","1.-w",RooArgList(acceptance));    
  //   RooDataSet origdata("data","data",mytree, RooArgSet(costheta,costheta1,costheta2,phi1,phi2,lb_id,w,keep) , "lb_id>0&&keep==1","w");
  //   //    RooDataSet* reddata = (RooDataSet*)origdata.reduce(RooFit::EventRange(1,5000));

  //   //    RooDataSet& data = *reddata;
  //   RooDataSet& data = origdata;

  //   data.Print("v");
  //   std::cout << "here" << std::endl;
  //   mygetchar;

  //   std::cout << __LINE__ << std::endl;
    
  //   TString output( TString(filename).ReplaceAll(".root",""));
    
  //   //     P_b.setVal(1.);
  //   //     alpha_b.setVal(-1.);
  //   //     r_0.setVal(0.);
  //   //     r_1.setVal(0.);
    
  //   // P_b.setVal(0.);
  //   // alpha_b.setVal(0.);
  //   // alpha_lambda.setVal(0.);
  //   // r_0.setVal(0.5);
  //   // r_1.setVal(0.);
    
  //   // alpha_plus.setVal(0.); 
  //   // alpha_minus.setVal(0.);
  //   // chi.setVal(0.); 
    
  //   //     alpha_plus.setVal(0.); alpha_plus.setConstant();
  //   //     alpha_minus.setVal(0.); alpha_minus.setConstant();
  //   //     chi.setVal(0.); chi.setConstant();
    
  //   //     w->forceNumInt();
    
    
  //   RooFitResult* fr = wn->fitTo(data, Save(true), Minos(true));
  //   fr->Print("v");
    
  //   plot(c, RooArgList(costheta,costheta1,costheta2),  *wn, &data, output);
    
  //   if (nVariables==5) {
  //     plot(c, RooArgList(phi1,phi2),   *wn, &data, output);
  //   }
  // }


 if (dofitmc) {
    std::cout << __LINE__ << std::endl;
    
    TFile myfile(filename);
    TDirectory* mydirectory = (TDirectory*)myfile.Get("BtoqllGeneratedDistributions");
    TTree* mytree = (TTree*)mydirectory->Get("Lb2JpsiL0");
    
    P_b=1.;
    alpha_b=-1.;
    r_0=0;
    r_1=0;
    //    alpha_lambda=-0.642;
    
    //    RooFormulaVar unmw("unmw","unmw","1.-w",RooArgList(acceptance));    
    RooDataSet origdata("data","data",mytree, RooArgSet(costheta,costheta1,costheta2,phi1,phi2,lb_id,w,keep) , "lb_id>0");
    //    RooDataSet* reddata = (RooDataSet*)origdata.reduce(RooFit::EventRange(1,5000));

    //    RooDataSet& data = *reddata;
    RooDataSet& data = origdata;

    data.Print("v");
    std::cout << "here" << std::endl;
    mygetchar;

    std::cout << __LINE__ << std::endl;
    
    TString output( TString(filename).ReplaceAll(".root",""));
    
    //     P_b.setVal(1.);
    //     alpha_b.setVal(-1.);
    //     r_0.setVal(0.);
    //     r_1.setVal(0.);
    
    // P_b.setVal(0.);
    // alpha_b.setVal(0.);
    // alpha_lambda.setVal(0.);
    // r_0.setVal(0.5);
    // r_1.setVal(0.);
    
    // alpha_plus.setVal(0.); 
    // alpha_minus.setVal(0.);
    // chi.setVal(0.); 
    
    //     alpha_plus.setVal(0.); alpha_plus.setConstant();
    //     alpha_minus.setVal(0.); alpha_minus.setConstant();
    //     chi.setVal(0.); chi.setConstant();
    
    //     w->forceNumInt();
    
    P_b=0.4;
    w3.fitTo(data, Save(false), Minos(false), ExternalConstraints(alpha_lambda_constrain));
    RooFitResult* fr = wn->fitTo(data, Save(true), Minos(true), ExternalConstraints(alpha_lambda_constrain));
    fr->Print("v");
    
    plot(c, RooArgList(costheta,costheta1,costheta2),  *wn, &data, output);
    
    if (nVariables==5) {
      plot(c, RooArgList(phi1,phi2),   *wn, &data, output);
    }
  }

  if (dotest) {
    const bool donotfit=true;

    if (donotfit) {
      P_b.setConstant();
      alpha_b.setConstant();
      r_0.setConstant();
      r_1.setConstant();
      alpha_plus.setConstant();
      alpha_minus.setConstant();
      chi.setConstant();
    }
   
    //    w->forceNumInt();
    std::cout << "generate toy data" << std::endl;

    alpha_plus.setVal(0.);
    alpha_minus.setVal(0.);
    chi.setVal(0.);


    P_b.setVal(1.);

    r_0.setVal(1.);
    r_1.setVal(1.);
    alpha_b.setVal(1.);
    RooDataSet* toydata_ap = wn->generate(RooArgSet(costheta,costheta1,costheta2,phi1,phi2), 10000);

    //    mygetchar;

    r_0.setVal(1.);
    r_1.setVal(-1.);
    alpha_b.setVal(-1.);
    RooDataSet* toydata_am = wn->generate(RooArgSet(costheta,costheta1,costheta2,phi1,phi2), 10000);

    //    mygetchar;


    r_0.setVal(0.);
    r_1.setVal(0.);
    alpha_b.setVal(1.);
    RooDataSet* toydata_bp = wn->generate(RooArgSet(costheta,costheta1,costheta2,phi1,phi2), 10000);

    //    mygetchar;

    r_0.setVal(0.);
    r_1.setVal(0.);
    alpha_b.setVal(-1.);
    RooDataSet* toydata_bm = wn->generate(RooArgSet(costheta,costheta1,costheta2,phi1,phi2), 10000);

    //    mygetchar;

    plotdata(c,costheta,  toydata_ap, "toy_ap_costheta.eps");
    plotdata(c,costheta1, toydata_ap, "toy_ap_costheta1.eps");
    plotdata(c,costheta2, toydata_ap, "toy_ap_costheta2.eps");
    plotdata(c,phi1, toydata_ap, "toy_ap_phi1.eps");
    plotdata(c,phi2, toydata_ap, "toy_ap_phi2.eps");
   


    plotdata(c,costheta,  toydata_am, "toy_am_costheta.eps");
    plotdata(c,costheta1, toydata_am, "toy_am_costheta1.eps");
    plotdata(c,costheta2, toydata_am, "toy_am_costheta2.eps");
    plotdata(c,phi1, toydata_am, "toy_am_phi1.eps");
    plotdata(c,phi2, toydata_am, "toy_am_phi2.eps");
   


    plotdata(c,costheta,  toydata_bp, "toy_bp_costheta.eps");
    plotdata(c,costheta1, toydata_bp, "toy_bp_costheta1.eps");
    plotdata(c,costheta2, toydata_bp, "toy_bp_costheta2.eps");
    plotdata(c,phi1, toydata_bp, "toy_bp_phi1.eps");
    plotdata(c,phi2, toydata_bp, "toy_bp_phi2.eps");
   


    plotdata(c,costheta,  toydata_bm, "toy_bm_costheta.eps");
    plotdata(c,costheta1, toydata_bm, "toy_bm_costheta1.eps");
    plotdata(c,costheta2, toydata_bm, "toy_bm_costheta2.eps");
    plotdata(c,phi1, toydata_bm, "toy_bm_phi1.eps");
    plotdata(c,phi2, toydata_bm, "toy_bm_phi2.eps");
   

    //    toydata->addColumn(theta);
    //    toydata->addColumn(theta1);
    //    toydata->addColumn(theta2);
    
    //    toydata->SaveAs("toy.root");

    // std::cout << "fit" << std::endl;
    // //    RooFitResult* result = w.fitTo(*toydata, Save(true));
    // //    result->Print("v");
    
    // std::cout << "plot" << std::endl;
    
    // plot(c, costheta,  *w, toydata, "costheta.eps");
    // plot(c, costheta1, *w, toydata, "costheta1.eps");
    // plot(c, costheta2, *w, toydata, "costheta2.eps");
    // plot(c, phi1,      *w, toydata, "phi1.eps");
    // plot(c, phi2,      *w, toydata, "phi2.eps");
    
    
    
  }

  const bool docompareroofitandevtgen=false;
  if (docompareroofitandevtgen) {

    std::cout << "docompareroofitandevtgen" << std::endl;
    
    TFile datafile("Tuple-Lb2L0Jpsi-Polb-ATLAS.root");
    TDirectory* datadirectory = (TDirectory*)datafile.Get("BtoqllGeneratedDistributions");
    TTree* datatree = (TTree*)datadirectory->Get("Lb2JpsiL0");
    
    const int nbin=50;
    TH1F evtgen_costheta("evtgen_costheta",  "evtgen_costheta", nbin,-1.,1.); 
    TH1F evtgen_costheta1("evtgen_costheta1","evtgen_costheta1", nbin,-1.,1.); 
    TH1F evtgen_costheta2("evtgen_costheta2","evtgen_costheta2", nbin,-1.,1.); 
    TH1F evtgen_phi1("evtgen_phi1","evtgen_phi1", nbin,-pi,pi); 
    TH1F evtgen_phi2("evtgen_phi2","evtgen_phi2", nbin,-pi,pi); 
    
    const int nevents=100000;
    datatree->Draw("costheta>>evtgen_costheta","lb_id>0","",nevents);
    datatree->Draw("costheta1>>evtgen_costheta1","lb_id>0","",nevents);
    datatree->Draw("costheta2>>evtgen_costheta2","lb_id>0","",nevents);
    datatree->Draw("phi1>>evtgen_phi1","lb_id>0","",nevents);
    datatree->Draw("phi2>>evtgen_phi2","lb_id>0","",nevents);
    
    TH1F roofit_costheta("roofit_costheta",  "roofit_costheta", 100,-1.,1.); 
    TH1F roofit_costheta1("roofit_costheta1","roofit_costheta1", 100,-1.,1.); 
    TH1F roofit_costheta2("roofit_costheta2","roofit_costheta2", 100,-1.,1.); 
    TH1F roofit_phi1("roofit_phi1","roofit_phi1", 100,-pi,pi); 
    TH1F roofit_phi2("roofit_phi2","roofit_phi2", 100,-pi,pi); 
    
    RooDataSet* data  = w5.generate(RooArgSet(costheta,costheta1,costheta2,phi1,phi2), nevents);

    data->fillHistogram(&roofit_costheta,  RooArgList(costheta));
    data->fillHistogram(&roofit_costheta1, RooArgList(costheta1));
    data->fillHistogram(&roofit_costheta2, RooArgList(costheta2));
    data->fillHistogram(&roofit_phi1, RooArgList(phi1));
    data->fillHistogram(&roofit_phi2, RooArgList(phi2));


    evtgen_costheta.GetXaxis()->SetTitle(costheta.GetTitle());
    evtgen_costheta.SetMinimum(0.);
    evtgen_costheta.SetLineColor(kRed);
    evtgen_costheta.Draw();
    roofit_costheta.Draw("SAME");
    c.SaveAs("comproofitevtgen_costheta.eps");


    evtgen_costheta1.GetXaxis()->SetTitle(costheta1.GetTitle());
    evtgen_costheta1.SetMinimum(0.);
    evtgen_costheta1.SetLineColor(kRed);
    evtgen_costheta1.Draw();
    roofit_costheta1.Draw("SAME");
    c.SaveAs("comproofitevtgen_costheta1.eps");


    evtgen_costheta2.GetXaxis()->SetTitle(costheta2.GetTitle());
    evtgen_costheta2.SetMinimum(0.);
    evtgen_costheta2.SetLineColor(kRed);
    evtgen_costheta2.Draw();
    roofit_costheta2.Draw("SAME");
    c.SaveAs("comproofitevtgen_costheta2.eps");


    evtgen_phi1.GetXaxis()->SetTitle(phi1.GetTitle());
    evtgen_phi1.SetMinimum(0.);
    evtgen_phi1.SetLineColor(kRed);
    evtgen_phi1.Draw();
    roofit_phi1.Draw("SAME");
    c.SaveAs("comproofitevtgen_phi1.eps");


    evtgen_phi2.GetXaxis()->SetTitle(phi2.GetTitle());
    evtgen_phi2.SetMinimum(0.);
    evtgen_phi2.SetLineColor(kRed);
    evtgen_phi2.Draw();
    roofit_phi2.Draw("SAME");
    c.SaveAs("comproofitevtgen_phi2.eps");


  }

  return 0;


}
