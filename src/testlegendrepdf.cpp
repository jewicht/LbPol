#include <fstream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TF3.h"
#include "TH3F.h"
#include "TCut.h"

#include "rootlogon.h"

#include "RooRealVar.h"
#include "RooLegendre3Pdfv2.h"
#include "RooLegendre5Pdfv2.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooChebychev.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooLbtoJpsiL0PDF3.h"
#include "createRooLegendre3.h"
#include "createRooLegendre5.h"

#include "legendre3f3.h"
#include "findminmaxtf3.h"

#include "functions-roofit.h"

double dcoeff3[ordermax3i][ordermax3j][ordermax3k];

Double_t myfunction0(Double_t *x, Double_t *par)
{
  const Double_t costheta  = x[0];
  const Double_t costheta1 = x[1];
  const Double_t costheta2 = x[2];
  return pdfval(costheta,costheta1,costheta2,dcoeff3);
}
Double_t myfunction1(Double_t *x, Double_t *par)
{
  const Double_t costheta  = x[0];
  const Double_t costheta1 = x[1];
  const Double_t costheta2 = x[2];
  return costheta*pdfval(costheta,costheta1,costheta2,dcoeff3);
}
Double_t myfunction2(Double_t *x, Double_t *par)
{
  const Double_t costheta  = x[0];
  const Double_t costheta1 = x[1];
  const Double_t costheta2 = x[2];
  return costheta1*pdfval(costheta,costheta1,costheta2,dcoeff3);
}
Double_t myfunction3(Double_t *x, Double_t *par)
{
  const Double_t costheta  = x[0];
  const Double_t costheta1 = x[1];
  const Double_t costheta2 = x[2];
  return costheta*costheta1*pdfval(costheta,costheta1,costheta2,dcoeff3);
}
Double_t myfunction4(Double_t *x, Double_t *par)
{
  const Double_t costheta  = x[0];
  const Double_t costheta1 = x[1];
  const Double_t costheta2 = x[2];
  return 0.5*(3.*costheta2*costheta2-1.)*pdfval(costheta,costheta1,costheta2,dcoeff3);
}
Double_t myfunction5(Double_t *x, Double_t *par)
{
  const Double_t costheta  = x[0];
  const Double_t costheta1 = x[1];
  const Double_t costheta2 = x[2];
  return 0.5*(3.*costheta2*costheta2-1.)*costheta*pdfval(costheta,costheta1,costheta2,dcoeff3);
}
Double_t myfunction6(Double_t *x, Double_t *par)
{
  const Double_t costheta  = x[0];
  const Double_t costheta1 = x[1];
  const Double_t costheta2 = x[2];
  return 0.5*(3.*costheta2*costheta2-1.)*costheta1*pdfval(costheta,costheta1,costheta2,dcoeff3);
}
Double_t myfunction7(Double_t *x, Double_t *par)
{
  const Double_t costheta  = x[0];
  const Double_t costheta1 = x[1];
  const Double_t costheta2 = x[2];
  return 0.5*(3.*costheta2*costheta2-1.)*costheta*costheta1*pdfval(costheta,costheta1,costheta2,dcoeff3);
}


void testintegral()
{
  RooRealVar x("x","x",-10.,10.);
  x.setRange("a",-1.,1.);
  x.setRange("b",-10.,10.);
  RooRealVar mean1("mean1","mean1",1.,0.,10.);
  RooRealVar sigma1("sigma1","sigma1",2.,0.,10.);
  RooRealVar mean2("mean2","mean2",1.,0.,10.);
  RooRealVar sigma2("sigma2","sigma2",1.,0.,10.);
  
  RooGaussian gauss1("gauss1","gauss1",x,mean1,sigma1);
  RooGaussian gauss2("gauss2","gauss2",x,mean2,sigma2);
  
  RooGenericPdf test("test","test","exp(-0.5*(x-mean1)*(x-mean1)/(sigma1*sigma1))*exp(-0.5*(x-mean2)*(x-mean2)/(sigma2*sigma2))",RooArgList(x,mean1,sigma1,mean2,sigma2));

  RooProdPdf mult("mult","mult",RooArgList(gauss1,gauss2));

  std::cout << "Integral = " << mult.createIntegral(RooArgSet(x))->getVal() << std::endl;
  std::cout << "Integral = " << mult.createIntegral(x,x,"a")->getVal() << std::endl;
  std::cout << "Integral = " << mult.createIntegral(x,x,"b")->getVal() << std::endl;
  std::cout << "Integral = " << mult.createRunningIntegral(RooArgSet(x))->getVal() << std::endl;
  std::cout << "Integral = " << test.createRunningIntegral(RooArgSet(x))->getVal() << std::endl;
  std::cout << "Integral = " << test.createRunningIntegral(RooArgSet(x),"a")->getVal() << std::endl;
 
}

void initdcoeff3(RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k])
{
  memset(dcoeff3,0,ordermax3i*ordermax3j*ordermax3k*sizeof(double));
  for (int i=0; i<ordermax3i;i++) {
    for (int j=0; j<ordermax3j;j++) {
      for (int k=0; k<ordermax3k;k++) {
	dcoeff3[i][j][k]=0.;
	if (coeff3[i][j][k]) dcoeff3[i][j][k]= coeff3[i][j][k]->getVal() / 8.;
      }
    }
  }
}

void sexmaxorder(RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k], int maxorder)
{
  for (int i=0; i<ordermax3i;i++) {
    for (int j=0; j<ordermax3j;j++) {
      for (int k=0; k<ordermax3k;k++) {
	if (coeff3[i][j][k] and (i>maxorder or j>maxorder or k>maxorder)) {
	  coeff3[i][j][k]->setVal(0.);
	  coeff3[i][j][k]->setConstant();
	}
      }
    }
  }
}





void setcoeff5(RooRealVar* coeff5[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m], RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k], TString coeff3name, int ii, int jj, int kk)
{
   for (int i=0; i<ordermax3i;i++) {
     for (int j=0; j<ordermax3j;j++) {
       for (int k=0; k<ordermax3k;k++) {
	 int i5[5];
	 i5[0]=0;
	 i5[1]=0;
	 i5[2]=0;
	 i5[3]=0;
	 i5[4]=0;
	 i5[ii]=i;
	 i5[jj]=j;
	 i5[kk]=k;
	 if (coeff5[i5[0]][i5[1]][i5[2]][i5[3]][i5[4]] && coeff3[i][j][k]) {
	   std::cout << "Setting coeff5[" <<i5[0] << "][" <<i5[1]<<"]["<<i5[2]<<"]["<< i5[3] << "][" << i5[4]<< "]  from " << coeff3name << " " << i << " " << j << " " << k << std::endl;
	   coeff5[i5[0]][i5[1]][i5[2]][i5[3]][i5[4]]->setVal(coeff3[i][j][k]->getVal());
	   if (fabs(coeff3[i][j][k]->getVal()) < fabs(coeff3[i][j][k]->getError())) {
	     coeff5[i5[0]][i5[1]][i5[2]][i5[3]][i5[4]]->setVal(0.);
	     coeff5[i5[0]][i5[1]][i5[2]][i5[3]][i5[4]]->setConstant(0.);
	     std::cout << " to 0. and constant" << std::endl;

	   }
	   
	   //	   coeff5[i5[0]][i5[1]][i5[2]][i5[3]][i5[4]]->setConstant();
	 }
       }
     }
   }
}


int main()
{  

  //  rootlogon();
  lhcbstyle();
  TCanvas c("c","c",100,100);
  
 // testintegral();return 0;

  static const double pi=atan(1.)*4.

  RooRealVar costheta("costheta",   "cos(#theta)",     -1., 1.);//RooArgList(theta));
  RooRealVar costheta1("costheta1", "cos(#theta_{1})", -1., 1.);// RooArgList(theta1));
  RooRealVar costheta2("costheta2", "cos(#theta_{2})", -1., 1.);// RooArgList(theta2));

  costheta.setRange("costhetapos",0.,1.);
  costheta1.setRange("costheta1pos",0.,1.);
  costheta2.setRange("costheta2pos",0.,1.);

  costheta.setRange("costhetaneg",-1.,0.);
  costheta1.setRange("costheta1neg",-1.,0.);
  costheta2.setRange("costheta2neg",-1.,0.);

  //  RooRealVar phi1("LambdabAngles_phi1","#phi_{1}",-pi,pi);
  //  RooRealVar phi2("LambdabAngles_phi2","#phi_{2}",-pi,pi);
  RooRealVar phi1("phi1","#phi_{1}/#pi",-1.,1.);
  RooRealVar phi2("phi2","#phi_{2}/#pi",-1.,1.);

  phi1.setRange("phi1pos",0.,1.);
  phi1.setRange("phi1neg",-1.,0.);
  phi2.setRange("phi2pos",0.,1.);
  phi2.setRange("phi2neg",-1.,0.);


  RooRealVar wacc("wacc","wacc",0.,5000.);


  RooRealVar P_b("P_b","P_{b}",0., -1.5, 1.5);//0.50
  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}",0.,0.,1.);
  alpha_lambda.setConstant();

  RooRealVar alpha_b("alpha_b","#alpha_{b}",0.,-1.1,1.1);//0.457
  RooRealVar r_0("r_0","r_{0}",0.,-0.3,2.3);//0.5
  RooRealVar r_1("r_1","r_{1}",0.,-1.3,1.3);//0.1

  const TString lambdabname = "lambda_b0";
  const TString lambdaname = "lambda0";
  const TString jpsiname = "jpsi1s";


  // RooFormulaVar phi1divpi("phi1divpi","LambdabAngles_phi1/3.14159",RooArgList(phi1));
  // RooFormulaVar phi2divpi("phi2divpi","LambdabAngles_phi2/3.14159",RooArgList(phi2));


  const TString sigfilename("DVTuples-Lb2JpsiL0-MC11a.root");
  const TString datafilename("DVTuples_data_stripping17.root");
  const TString dirname("JpsiL0Tuple_betas");

  TString sigbdtfilename(sigfilename);
  sigbdtfilename.ReplaceAll(".root",".bdt.root");

  TString siganglesfilename(sigfilename);
  siganglesfilename.ReplaceAll(".root", ".angles.root");

  TFile sigfile(sigfilename);
  TDirectory* sigdir = (TDirectory*)sigfile.Get(dirname);
  TTree* sigtree = (TTree*)sigdir->Get("DecayTree");
  sigtree->AddFriend(dirname, sigbdtfilename);
  sigtree->AddFriend(dirname, siganglesfilename);


   
  TString databdtfilename(datafilename);
  databdtfilename.ReplaceAll(".root",".bdt.root");

  TString dataanglesfilename(datafilename);
  dataanglesfilename.ReplaceAll(".root", ".angles.root");

  TFile datafile(datafilename);
  TDirectory* datadir = (TDirectory*)datafile.Get(dirname);
  TTree* datatree = (TTree*)datadir->Get("DecayTree");
  datatree->AddFriend(dirname, databdtfilename);
  datatree->AddFriend(dirname, dataanglesfilename);

  const TCut Triggercut("Hlt1DiMuonHighMassDecision==1&&Hlt2DiMuonDetachedJPsiDecision==1");
  const TCut bdtcut_LL("bdt>0.");
  const TCut bdtcut_DD("bdt>0.1");

  const TCut Lambda0_LL("piminus_TRACK_Type==3");
  const TCut Lambda0_DD("piminus_TRACK_Type==5");

  const TCut TrueLambdab("abs(lambda_b0_TRUEID)==5122");

  plotttree(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), sigtree, TrueLambdab + Lambda0_LL, TrueLambdab + Lambda0_DD, "acc_llvsdd" );

  plotttree(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), sigtree, TrueLambdab + Lambda0_LL + TCut("bdt<0."), TrueLambdab + Lambda0_LL + TCut("bdt>0."), "acc_ll_bdt" );
  plotttree(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), sigtree, TrueLambdab + Lambda0_DD + TCut("bdt<0.1"), TrueLambdab + Lambda0_DD + TCut("bdt>0.1"), "acc_dd_bdt" );

  plotttree(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), sigtree, TrueLambdab + TCut("Hlt2DiMuonDetachedDecision==1"), TrueLambdab + TCut("Hlt2DiMuonDetachedDecision==0"), "acc_hlt2detached" );

  plotttree(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), sigtree, TrueLambdab + TCut("Hlt2DiMuonJPsiDecision==1"), TrueLambdab + TCut("Hlt2DiMuonJPsiDecision==0"), "acc_hlt2jpsi" );

  //  plotttree(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), datatree, Lambda0_LL + TCut("bdt>0.") + TCut("lambda_b0_M<5550."),  Lambda0_LL + TCut("bdt>0.") + TCut("lambda_b0_M>5700."), "data_ll_mass" );
  //  plotttree(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), datatree, Lambda0_DD + TCut("bdt>0.1") + TCut("lambda_b0_M<5550."),  Lambda0_DD + TCut("bdt>0.1") + TCut("lambda_b0_M>5700."), "data_dd_mass" );

  //  plotttree(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), datatree, Lambda0_DD + TCut("lambda_b0_M<5550."),  Lambda0_DD + TCut("lambda_b0_M>5700."), "data_dd_masswoutbdt" );


  //  return 1;


  //  TTree* subtree_ll = sigtree->CopyTree(TrueLambdab + Lambda0_LL + bdtcut_LL + Triggercut);
  //  TTree* subtree_dd = sigtree->CopyTree(TrueLambdab + Lambda0_DD + bdtcut_DD + Triggercut);

  //  std::cout << "subtree_ll candidates = " << subtree_ll->GetEntries() << std::endl;  
  //  std::cout << "subtree_dd candidates = " << subtree_dd->GetEntries() << std::endl;
  //  return 1;



  //  TTree* subtree = sigtree->CopyTree(TrueLambdab + bdtcut_LL + Triggercut);
  //  std::cout << "Ncorr rec =" << subtree->GetEntries() << std::endl;

//   TH1F acc_mc_costheta_histo("acc_mc_costheta_histo","acc_mc_costheta_histo",25,-1.,1.);
//   TH1F acc_mc_costheta1_histo("acc_mc_costheta1_histo","acc_mc_costheta1_histo",25,-1.,1.);
//   TH1F acc_mc_costheta2_histo("acc_mc_costheta2_histo","acc_mc_costheta2_histo",25,-1.,1.);

//   TH1F acc_toy_costheta_histo("acc_toy_costheta_histo","acc_toy_costheta_histo",25,-1.,1.);
//   TH1F acc_toy_costheta1_histo("acc_toy_costheta1_histo","acc_toy_costheta1_histo",25,-1.,1.);
//   TH1F acc_toy_costheta2_histo("acc_toy_costheta2_histo","acc_toy_costheta2_histo",25,-1.,1.);

//   RooDataSet* toydata = RooDataSet::read("toys/toy_0000.dat", RooArgList(costheta,costheta1,costheta2,wacc));
//   toydata->fillHistogram(&acc_toy_costheta_histo,RooArgList(costheta));
//   toydata->fillHistogram(&acc_toy_costheta1_histo,RooArgList(costheta1));
//   toydata->fillHistogram(&acc_toy_costheta2_histo,RooArgList(costheta2));
//   subtree->Draw("costheta>>acc_mc_costheta_histo");
//   subtree->Draw("costheta1>>acc_mc_costheta1_histo");
//   subtree->Draw("costheta2>>acc_mc_costheta2_histo");

//   acc_mc_costheta_histo.Sumw2();
//   acc_mc_costheta1_histo.Sumw2();
//   acc_mc_costheta2_histo.Sumw2();

//   acc_toy_costheta_histo.Sumw2();
//   acc_toy_costheta1_histo.Sumw2();
//   acc_toy_costheta2_histo.Sumw2();

//   acc_mc_costheta_histo.GetXaxis()->SetTitle(costheta.GetTitle());
//   acc_mc_costheta1_histo.GetXaxis()->SetTitle(costheta1.GetTitle());
//   acc_mc_costheta2_histo.GetXaxis()->SetTitle(costheta2.GetTitle());
//   acc_mc_costheta_histo.SetMinimum(0.);
//   acc_mc_costheta1_histo.SetMinimum(0.);
//   acc_mc_costheta2_histo.SetMinimum(0.);

//   acc_toy_costheta_histo.SetLineColor(kRed);
//   acc_toy_costheta1_histo.SetLineColor(kRed);
//   acc_toy_costheta2_histo.SetLineColor(kRed);

//   acc_mc_costheta_histo.DrawNormalized("E");
//   acc_toy_costheta_histo.DrawNormalized("SAME");
//   c.SaveAs("acc_mcvstoy_costheta_histo.eps");
  
//   acc_mc_costheta1_histo.DrawNormalized("E");
//   acc_toy_costheta1_histo.DrawNormalized("SAME");
//   c.SaveAs("acc_mcvstoy_costheta1_histo.eps");
  
//   acc_mc_costheta2_histo.DrawNormalized("E");
//   acc_toy_costheta2_histo.DrawNormalized("SAME");
//   c.SaveAs("acc_mcvstoy_costheta2_histo.eps");

//   mygetchar;

//   TH3F acc_histo("acc_histo","acc_histo",
// 		 10,-1.,1.,
// 		 10,-1.,1.,
// 		 10,-1.,1.);

//   subtree->Draw("costheta:costheta1:costheta2>>acc_histo");
//   acc_histo.SaveAs("acc_histo.root");

  //  RooRealVar piminus_TRACK_Type("piminus_TRACK_Type","piminus_TRACK_Type",0.,10.);

  //  RooDataSet data("data","data",subtree,RooArgSet(costheta,costheta1,costheta2,phi1,phi2));
  RooDataSet* datatmp = filldataset("data","data",RooArgList(costheta,costheta1,costheta2,phi1,phi2),sigtree,TrueLambdab + Lambda0_LL);
  RooDataSet& data = *datatmp;

  //add fake point at 0.0412852  0.999992  -1
  RooRealVar weight("weight","weight",1.);
  data.addColumn(weight);

  costheta.setVal(0.0412852);
  costheta1.setVal(0.99999);
  costheta2.setVal(-1.);
  costheta.setVal(-1.);
  costheta1.setVal(1.);
  costheta2.setVal(-1.);
  phi1.setVal(0.);
  phi2.setVal(0.);
  weight.setVal(1.);
  data.add(RooArgSet(costheta,costheta1,costheta2,phi1,phi2,weight));

  costheta.setVal(-1.);
  costheta1.setVal(1.);
  costheta2.setVal(-.5);
  phi1.setVal(0.);
  phi2.setVal(0.);
  weight.setVal(1.);
  data.add(RooArgSet(costheta,costheta1,costheta2,phi1,phi2,weight));

  //  data.write("bbaba.txt");
  
  RooDataSet wdata("wdata","wdata",&data,*data.get(),0,weight.GetName());
  //  costheta.setBins(6);
  //  costheta1.setBins(6);
  //  costheta2.setBins(6);
  //  phi1.setBins(6);
  //  phi2.setBins(10);
  //  RooDataHist* databinned = new RooDataHist("databinned","databinned",RooArgList(costheta,costheta1,costheta2,phi1,phi2), data);

  RooDataSet data2("data2","data2",&data,RooArgSet(costheta,costheta1,costheta2,phi1,phi2));
  data2.write("acc_histo5.dat");
  //  TFile *f = new TFile("acc_histo5.root","recreate");
  //  databinned->Write("databinned");
  //  f->Close();
    

  //  mygetchar;

  // RooAbsPdf* prodpdf[8];

  // RooGenericPdf* Fi[8];
  // //  Fi[0] = new RooGenericPdf("F0","F0","1.");
  // Fi[1] = new RooGenericPdf("F1","F1","costheta",RooArgList(costheta));
  // Fi[2] = new RooGenericPdf("F2","F2","costheta1",RooArgList(costheta1));
  // Fi[3] = new RooGenericPdf("F3","F3","costheta*costheta1",RooArgList(costheta,costheta1));
  // Fi[4] = new RooGenericPdf("F4","F4","0.5*(3.*costheta2*costheta2-1.)",RooArgList(costheta2));
  // Fi[5] = new RooGenericPdf("F5","F5","0.5*(3.*costheta2*costheta2-1.)*costheta",RooArgList(costheta, costheta2));
  // Fi[6] = new RooGenericPdf("F6","F6","0.5*(3.*costheta2*costheta2-1.)*costheta1",RooArgList(costheta1, costheta2));
  // Fi[7] = new RooGenericPdf("F7","F7","0.5*(3.*costheta2*costheta2-1.)*costheta*costheta1",RooArgList(costheta,costheta1, costheta2));

  
 // RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k];
 // RooRealVar* doeff3[ordermax3i][ordermax3j][ordermax3k];
 // RooRealVar* eoeff3[ordermax3i][ordermax3j][ordermax3k];
 // RooRealVar* foeff3[ordermax3i][ordermax3j][ordermax3k];

 
// RooLegendre3Pdf test3("test3","test3",costheta,costheta1,costheta2,c000,c001,c002,c010,c011,c012,c020,c021,c022,c100,c101,c102,c110,c111,c112,c120,c121,c122,c200,c201,c202,c210,c211,c212,c220,c221,c222);
// RooLegendre3Pdf test3("test3","test3",costheta,costheta1,costheta2,c000,c001,c002,c004,c006,c010,c011,c012,c014,c016,c020,c021,c022,c024,c026,c100,c101,c102,c104,c106,c110,c111,c112,c114,c116,c120,c121,c122,c124,c126,c200,c201,c202,c204,c206,c210,c211,c212,c214,c216,c220,c221,c222,c224,c226);

 // RooLegendre3Pdf* test3 = NULL;
 // RooLegendre3Pdf* test4 = NULL;
 // RooLegendre3Pdf* test5 = NULL;
 // RooLegendre3Pdf* test6 = NULL;

 // std::cout << __LINE__ << " " << test3 << std::endl;
 // createRooLegendre3(costheta, costheta1, costheta2, "test3", test3, coeff3);
 // createRooLegendre3(costheta, phi1,      phi2,      "test4", test4, doeff3);
 // createRooLegendre3(costheta, costheta1, phi2,      "test5", test5, eoeff3);
 // createRooLegendre3(costheta, costheta2, phi2,      "test6", test6, foeff3);
 // std::cout << __LINE__ << " " << test3 << std::endl;
 // // mygetchar;

 // prodpdf[0] = test3; 
 // for (int i=1; i<8; i++ ) prodpdf[i] = new RooProdPdf("prodpdf" + istr(i),"prodpdf" + istr(i),RooArgList(*test3,*Fi[i]));

 // RooLegendre3Pdf test4("test4","test4",costheta,phi1,phi2,d000,d001,d002,d004,d006,d010,d011,d012,d014,d016,d020,d021,d022,d024,d026,d100,d101,d102,d104,d106,d110,d111,d112,d114,d116,d120,d121,d122,d124,d126,d200,d201,d202,d204,d206,d210,d211,d212,d214,d216,d220,d221,d222,d224,d226);


 // test3.forceNumInt(true);



 RooRealVar *cc0c1c2[ordermax3i][ordermax3j][ordermax3k], *cc0c1p1[ordermax3i][ordermax3j][ordermax3k], *cc0c1p2[ordermax3i][ordermax3j][ordermax3k];
 RooRealVar *cc0c2p1[ordermax3i][ordermax3j][ordermax3k], *cc0c2p2[ordermax3i][ordermax3j][ordermax3k];
 RooRealVar *cc0p1p2[ordermax3i][ordermax3j][ordermax3k];
 RooRealVar *cc1c2p1[ordermax3i][ordermax3j][ordermax3k], *cc1c2p2[ordermax3i][ordermax3j][ordermax3k];
 RooRealVar *cc1p1p2[ordermax3i][ordermax3j][ordermax3k];
 RooRealVar *cc2p1p2[ordermax3i][ordermax3j][ordermax3k];

 RooAbsPdf *c0c1c2, *c0c1p1, *c0c1p2;
 RooAbsPdf *c0c2p1, *c0c2p2;
 RooAbsPdf *c0p1p2;
 RooAbsPdf *c1c2p1, *c1c2p2;
 RooAbsPdf *c1p1p2;
 RooAbsPdf *c2p1p2;

 createRooLegendre3v2(costheta, costheta1, costheta2, "c0c1c2", c0c1c2, cc0c1c2, true, 4, 4, 7, 10,5); 
 createRooLegendre3v2(costheta, costheta1,      phi1, "c0c1p1", c0c1p1, cc0c1p1, true, 4, 4, 7, 10,5); 
 createRooLegendre3v2(costheta, costheta1,      phi2, "c0c1p2", c0c1p2, cc0c1p2, true, 4, 4, 7, 10,5); 

 createRooLegendre3v2(costheta, costheta2,      phi1, "c0c2p1", c0c2p1, cc0c2p1, true, 4, 4, 7, 10,5); 
 createRooLegendre3v2(costheta, costheta2,      phi2, "c0c2p2", c0c2p2, cc0c2p2, true, 4, 4, 7, 10,5); 

 createRooLegendre3v2(costheta,      phi1,      phi2, "c0p1p2", c0p1p2, cc0p1p2, true, 4, 4, 7, 10,5); 

 createRooLegendre3v2(costheta1, costheta2,      phi1, "c1c2p1", c1c2p1, cc1c2p1, true, 4, 4, 7, 10,5); 
 createRooLegendre3v2(costheta1, costheta2,      phi2, "c1c2p2", c1c2p2, cc1c2p2, true, 4, 4, 7, 10,5); 

 createRooLegendre3v2(costheta1,      phi1,      phi2, "c1p1p2", c1p1p2, cc1p1p2, true, 4, 4, 7, 10,5); 

 createRooLegendre3v2(costheta2,      phi1,      phi2, "c2p1p2", c2p1p2, cc2p1p2, true, 4, 4, 7, 10,5); 


 RooRealVar* coeff5[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m];

 // std::cout << data.numEntries() << std::endl;
 // // mygetchar;
 // for (int i=2; i<ordermaxi;i++) {
 //   for (int j=2; j<ordermaxj;j++) {
 //     for (int k=2; k<ordermaxk;k++) {
 //       for (int l=2; l<ordermaxl;l++) {
 // 	 for (int m=2; m<ordermaxm;m++) {
 // 	   if (((i%2==0) or (i==1)) and 
 // 	       ((j%2==0) or (j==1)) and 
 // 	       ((k%2==0) or (k==1)) and 
 // 	       ((l%2==0) or (l==1)) and 
 // 	       ((m%2==0) or (m==1))) {

 // 	     coeff5[i][j][k][l][m]->setVal(0.);
 // 	     coeff5[i][j][k][l][m]->setConstant();
 // 	   }
 // 	 }
 //       }
 //     }  
 //   }
 // }


 RooLegendre5Pdfv2* test = NULL;
 createRooLegendre5v2(costheta,costheta1,costheta2,phi1,phi2,"test",test, coeff5, 4,4,4,4,7,8,2);
 

 // for (int i=0; i<ordermaxi;i++) {
 //   for (int j=0; j<ordermaxj;j++) {
 //     for (int k=0; k<ordermaxk;k++) {
 //       if (coeff3[i][j][k] && k>=4) coeff3[i][j][k]->setConstant(true);
 //     }  
 //   }
 // }

 // RooFitResult* test4fr = test4.fitTo(data, RooFit::Save(true));
 // plot(c, costheta, test4, &data, "acceptance-test4-costheta.eps");
 // plot(c, phi1,     test4, &data, "acceptance-test4-phi1.eps");
 // plot(c, phi2,     test4, &data, "acceptance-test4-phi2.eps");

 // RooArgSet tmp(test4fr->floatParsFinal());
 // tmp.writeToFile("test4.txt");

 // mygetchar;
 
 // for (int i=0; i<ordermaxi;i++) {
 //   for (int j=0; j<ordermaxj;j++) {
 //     for (int k=0; k<ordermaxk;k++) {
 //       if (coeff3[i][j][k]) {
 // 	 coeff3[i][j][k]->setMin(-1.1);
 // 	 coeff3[i][j][k]->setMax(+1.1);
 //       }
 //     }  
 //   }
 // }
 // coeff3[0][0][0]->setVal(1.);
 // coeff3[0][0][0]->setConstant();
 // for (int i=0; i<ordermaxi;i++) {
 //   for (int j=0; j<ordermaxj;j++) {
 //     for (int k=0; k<ordermaxk;k++) {
 //       if (coeff3[i][j][k]) {
 // 	 if ((i>4) or (j>4) or (k>4)) {
 // 	   coeff3[i][j][k]->setVal(0.);
 // 	   coeff3[i][j][k]->setConstant();
 // 	 }
 //       }
 //     }  
 //   }
 // }

 // RooArgList params3("params3");
 //  for (int i=0; i<ordermaxi;i++) {
 //   for (int j=0; j<ordermaxj;j++) {
 //     for (int k=0; k<ordermaxk;k++) {
 //       if (coeff3[i][j][k]) params3.add(*coeff3[i][j][k]);
 //     }  
 //   }
 // }

 const bool do3dfit=false;

 TF3 func0("func0",myfunction0, -1., 1., -1., 1., -1., 1.);
 TF3 func1("func1",myfunction1, -1., 1., -1., 1., -1., 1.);
 TF3 func2("func2",myfunction2, -1., 1., -1., 1., -1., 1.);
 TF3 func3("func3",myfunction3, -1., 1., -1., 1., -1., 1.);
 TF3 func4("func4",myfunction4, -1., 1., -1., 1., -1., 1.);
 TF3 func5("func5",myfunction5, -1., 1., -1., 1., -1., 1.);
 TF3 func6("func6",myfunction6, -1., 1., -1., 1., -1., 1.);
 TF3 func7("func7",myfunction7, -1., 1., -1., 1., -1., 1.);


 if (do3dfit) {

  costheta.setBins(25);
  costheta1.setBins(25);
  costheta2.setBins(25);


   sexmaxorder(cc0c1c2,3);
   RooFitResult* frc0c1c2 = c0c1c2->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
   frc0c1c2->SaveAs("frc0c1c2.root");
   frc0c1c2->Print("v");//floatParsFinal().printMultiline(std::cout);
   writecoeff3("cc0c1c2.txt", cc0c1c2);

   initdcoeff3(cc0c1c2);
   
   RooRealVar* xi[8];
   xi[0]= new RooRealVar("xi0","xi0",func0.Integral(-1.,1.,-1.,1.,-1.,1.));
   xi[1]= new RooRealVar("xi1","xi1",func1.Integral(-1.,1.,-1.,1.,-1.,1.));
   xi[2]= new RooRealVar("xi2","xi2",func2.Integral(-1.,1.,-1.,1.,-1.,1.));
   xi[3]= new RooRealVar("xi3","xi3",func3.Integral(-1.,1.,-1.,1.,-1.,1.));
   xi[4]= new RooRealVar("xi4","xi4",func4.Integral(-1.,1.,-1.,1.,-1.,1.));
   xi[5]= new RooRealVar("xi5","xi5",func5.Integral(-1.,1.,-1.,1.,-1.,1.));
   xi[6]= new RooRealVar("xi6","xi6",func6.Integral(-1.,1.,-1.,1.,-1.,1.));
   xi[7]= new RooRealVar("xi7","xi7",func7.Integral(-1.,1.,-1.,1.,-1.,1.));
   for (int i=0; i<8; i++) std::cout << "xi[" << i << "] = " << xi[i]->getVal() << std::endl;

   RooLbtoJpsiL0PDF3 w3("w3","w3",costheta,costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);//,*xi[0],*xi[1],*xi[2],*xi[3],*xi[4],*xi[5],*xi[6],*xi[7]);
   
   RooLbtoJpsiL0wAccPDF3* fullpdf = NULL;
   
   createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b, alpha_lambda, r_0, r_1, "wxacc", fullpdf, cc0c1c2, false);

   plotpdf(c,costheta,w3,"w3mod-costheta.eps");
   plotpdf(c,costheta1,w3,"w3mod-costheta1.eps");
   plotpdf(c,costheta2,w3,"w3mod-costheta2.eps");

   plotpdf(c,costheta,*fullpdf,"w3wacc-costheta.eps");
   plotpdf(c,costheta1,*fullpdf,"w3wacc-costheta1.eps");
   plotpdf(c,costheta2,*fullpdf,"w3wacc-costheta2.eps");

   for (int i=0; i<8; i++) xi[i]->setVal(1.);
   
   plotpdf(c,costheta,w3,"w3-costheta.eps");
   plotpdf(c,costheta1,w3,"w3-costheta1.eps");
   plotpdf(c,costheta2,w3,"w3-costheta2.eps");

   plot(c, RooArgList(costheta,costheta1,costheta2), *c0c1c2, &data, "acceptance-c0c1c2");

   mygetchar;


   RooFitResult* frc0c1p1 = c0c1p1->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
   frc0c1p1->Print("v");
   writecoeff3("cc0c1p1.txt", cc0c1p1);
   mygetchar;
   plot(c, RooArgList(costheta,costheta1,phi1), *c0c1p1, &data, "acceptance-c0c1p1");
   RooFitResult* frc0c1p2 = c0c1p2->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
   frc0c1p2->Print("v");
   writecoeff3("cc0c1p2.txt", cc0c1p2);
   plot(c, RooArgList(costheta,costheta1,phi2), *c0c1p2, &data, "acceptance-c0c1p2");
   mygetchar;

   RooFitResult* frc0c2p1 = c0c2p1->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
   frc0c2p1->Print("v");
   writecoeff3("cc0c2p1.txt", cc0c2p1);
   mygetchar;
   plot(c, RooArgList(costheta,costheta2,phi1), *c0c2p1, &data, "acceptance-c0c2p1");
   RooFitResult* frc0c2p2 = c0c2p2->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
   frc0c2p2->Print("v");
   writecoeff3("cc0c2p2.txt", cc0c2p2);
   mygetchar;
   plot(c, RooArgList(costheta,costheta2,phi2), *c0c2p2, &data, "acceptance-c0c2p2");

   RooFitResult* frc0p1p2 = c0p1p2->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
   frc0p1p2->Print("v");
   writecoeff3("cc0p1p2.txt", cc0p1p2);
   mygetchar;
   plot(c, RooArgList(costheta,phi1,phi2), *c0p1p2, &data, "acceptance-c0p1p2");

   RooFitResult* frc1c2p1 = c1c2p1->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
   frc1c2p1->Print("v");
   writecoeff3("cc1c2p1.txt", cc1c2p1);
   mygetchar;
   plot(c, RooArgList(costheta1,costheta2,phi1), *c1c2p1, &data, "acceptance-c1c2p1");
   RooFitResult* frc1c2p2 = c1c2p2->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
   frc1c2p2->Print("v");
   writecoeff3("cc1c2p2.txt", cc1c2p2);
   mygetchar;
   plot(c, RooArgList(costheta1,costheta2,phi2), *c1c2p2, &data, "acceptance-c1c2p2");

   RooFitResult* frc1p1p2 = c1p1p2->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
   frc1p1p2->Print("v");
   writecoeff3("cc1p1p2.txt", cc1p1p2);
   mygetchar;
   plot(c, RooArgList(costheta1,phi1,phi2), *c1p1p2, &data, "acceptance-c1p1p2");

   RooFitResult* frc2p1p2 = c2p1p2->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
   frc2p1p2->Print("v");
   writecoeff3("cc2p1p2.txt", cc2p1p2);
   mygetchar;
   plot(c, RooArgList(costheta2,phi1,phi2), *c2p1p2, &data, "acceptance-c2p1p2");

 } else {

   readcoeff3_rrv("cc0c1c2.txt", 
		   cc0c1c2);
   readcoeff3_rrv("cc0c1p1.txt", 
		   cc0c1p1);
   readcoeff3_rrv("cc0c1p2.txt", 
		   cc0c1p2);

   readcoeff3_rrv("cc0c2p1.txt", 
		   cc0c2p1);
   readcoeff3_rrv("cc0c2p2.txt", 
		   cc0c2p2);


   readcoeff3_rrv("cc0p1p2.txt", 
		   cc0p1p2);

   readcoeff3_rrv("cc1c2p1.txt", 
		   cc1c2p1);
   readcoeff3_rrv("cc1c2p2.txt", 
		   cc1c2p2);

   readcoeff3_rrv("cc1p1p2.txt", 
		   cc1p1p2);
   readcoeff3_rrv("cc2p1p2.txt", 
		   cc2p1p2);

 }
 // exit(1);


 // RooFitResult* test3fr = test3->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
 // plot(c, costheta,  *test3, &data, "acceptance-test3-costheta.eps");
 // plot(c, costheta1, *test3, &data, "acceptance-test3-costheta1.eps");
 // plot(c, costheta2, *test3, &data, "acceptance-test3-costheta2.eps");
 // RooFitResult* test4fr = test4->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
 // plot(c, costheta, *test4, &data, "acceptance-test4-costheta.eps");
 // plot(c,     phi1, *test4, &data, "acceptance-test4-phi1.eps");
 // plot(c,     phi2, *test4, &data, "acceptance-test4-phi2.eps");
 // RooFitResult* test5fr = test5->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
 // plot(c,  costheta, *test5, &data, "acceptance-test5-costheta.eps");
 // plot(c, costheta1, *test5, &data, "acceptance-test5-costheta1.eps");
 // plot(c,      phi2, *test5, &data, "acceptance-test5-phi2.eps");
 // RooFitResult* test6fr = test6->fitTo(data, RooFit::Save(true), RooFit::Minos(false));
 // plot(c,  costheta, *test6, &data, "acceptance-test6-costheta.eps");
 // plot(c, costheta2, *test6, &data, "acceptance-test6-costheta2.eps");
 // plot(c,      phi2, *test6, &data, "acceptance-test6-phi2.eps");

 std::cout << "fit 3D done" << std::endl;
 mygetchar;


 setcoeff5(coeff5,cc0c1c2,"cc0c1c2",0,1,2);
 setcoeff5(coeff5,cc0c1p1,"cc0c1p1",0,1,3);
 setcoeff5(coeff5,cc0c1p2,"cc0c1p2",0,1,4);

 setcoeff5(coeff5,cc0c2p1,"cc0c2p1",0,2,3);
 setcoeff5(coeff5,cc0c2p2,"cc0c2p2",0,2,4);

 setcoeff5(coeff5,cc0p1p2,"cc0p1p2",0,3,4);

 setcoeff5(coeff5,cc1c2p1,"cc1c2p1",1,2,3);
 setcoeff5(coeff5,cc1c2p2,"cc1c2p2",1,2,4);

 setcoeff5(coeff5,cc1p1p2,"cc1p1p2",1,3,4);

 setcoeff5(coeff5,cc2p1p2,"cc2p1p2",2,3,4);

 // for (int i=0; i<ordermax3i;i++) {
 //   for (int j=0; j<ordermax3j;j++) {
 //     for (int k=0; k<ordermax3k;k++) {
 // 	   //	   if (((i%2==0) or (i==1)) and 
 // 	   //	       ((j%2==0) or (j==1)) and 
 // 	   //	       ((k%2==0) or (k==1)) and 
 // 	   //	       ((l%2==0) or (l==1)) and 
 // 	   //	       ((m%2==0) or (m==1))) {

 //       if (coeff5[i][j][k][0][0] && coeff3[i][j][k])  {
 // 	 std::cout << "Setting coeff5 from coeff3 " << i << " " << j << " " << k << std::endl;
 // 	 coeff5[i][j][k][0][0]->setVal(coeff3[i][j][k]->getVal());
 //       }

 //       if (coeff5[i][0][0][j][k] && doeff3[i][j][k])  {
 // 	 std::cout << "Setting coeff5 from doeff3 " << i << " " << j << " " << k << std::endl;
 // 	 coeff5[i][0][0][j][k]->setVal(doeff3[i][j][k]->getVal());
 //       }


 //       if (coeff5[i][j][0][0][k] && eoeff3[i][j][k])  {
 // 	 std::cout << "Setting coeff5 from doeff3 " << i << " " << j << " " << k << std::endl;
 // 	 coeff5[i][j][0][0][k]->setVal(eoeff3[i][j][k]->getVal());
 //       }


 //       if (coeff5[i][0][j][0][k] && foeff3[i][j][k])  {
 // 	 std::cout << "Setting coeff5 from doeff3 " << i << " " << j << " " << k << std::endl;
 // 	 coeff5[i][0][j][0][k]->setVal(foeff3[i][j][k]->getVal());
 //       }

 // 	     //	     coeff5[i][j][k][l][m]->setConstant();
 // 	     //	   }
 //     }  
 //   }
 // }


 for (int i=0; i<ordermax5i;i++) {
   for (int j=0; j<ordermax5j;j++) {
     for (int k=0; k<ordermax5k;k++) {
       for (int l=0; l<ordermax5l;l++) {
	 for (int m=0; m<ordermax5m;m++) {
	   if (coeff5[i][j][k][l][m]) {
	     if (coeff5[i][j][k][l][m]->getVal()==0.) {
	       //std::cout << "coeff5:" << i << j << k << l << m << " set constant" << std::endl;
	       //	       coeff5[i][j][k][l][m]->setConstant(true);
	     } else {
	       //	       coeff5[i][j][k][l][m]->setConstant(false);
	     }
	     //	     if (i+j+k+l+m>4)  coeff5[i][j][k][l][m]->setConstant(true);
	   }
	 }
       }
     }
   }  
 }
 coeff5[0][0][0][0][0]->setConstant(true);

 int paramconstant=0;
 int paramnonconstant=0;
 
 for (int i=0; i<ordermax5i;i++) {
   for (int j=0; j<ordermax5j;j++) {
     for (int k=0; k<ordermax5k;k++) {
       for (int l=0; l<ordermax5l;l++) {
	 for (int m=0; m<ordermax5m;m++) {
	   const RooRealVar* rrv = coeff5[i][j][k][l][m];
	   if (rrv) {
	     if (rrv->isConstant()) {paramconstant++;} else {paramnonconstant++;}
	     std::cout << "coeff5:" << i << j << k << l << m << " ("<< rrv->isConstant() << ") = " << rrv->getVal() << " +- " << rrv->getError()  << std::endl;
	   }
	 }
       }
     }
   }  
 }
 std::cout << "Param constants     = " << paramconstant << std::endl;
 std::cout << "Param non constants = " << paramnonconstant << std::endl;
 

 std::cout << "fit 5D starting" << std::endl;
 mygetchar;

 writecoeff5("coeff5-nofit.txt", coeff5);

 // plot(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), *test, databinned, "acceptance-testbinned-nofit");
 plot(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), *test, &data, "acceptance-test-nofit");

 plot(c, RooArgList(costheta1,costheta2,phi1,phi2), *test, &data,  "acceptance-test-nofit",false,NULL,"costhetaneg"); 
 plot(c, RooArgList(costheta1,costheta2,phi1,phi2), *test, &data,  "acceptance-test-nofit" ,false,NULL,"costhetapos"); 
 plot(c, RooArgList(costheta,costheta2,phi1,phi2), *test, &data,  "acceptance-test-nofit" ,false,NULL,"costheta1neg"); 
 plot(c, RooArgList(costheta,costheta2,phi1,phi2), *test, &data,  "acceptance-test-nofit" ,false,NULL,"costheta1pos"); 
 plot(c, RooArgList(costheta,costheta1,phi1,phi2), *test, &data,  "acceptance-test-nofit" ,false,NULL,"costheta2neg"); 
 plot(c, RooArgList(costheta,costheta1,phi1,phi2), *test, &data, "acceptance-test-nofit" ,false,NULL,"costheta2pos"); 

 // std::cout << "binned fit" << std::endl;

 // test->fitTo(*databinned,RooFit::NumCPU(8));
 // plot(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), *test, databinned, "acceptance-testbinnned");
 // std::cout << "binned fit done" << std::endl;
 // mygetchar;
  
 // plot(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), *test, databinned, "acceptance-testbinned");
 
 // test->forceNumInt();

 test->fitTo(data);//,RooFit::NumCPU(8));
 plot(c, RooArgList(costheta,costheta1,costheta2,phi1,phi2), *test, &data, "acceptance-test");

 writecoeff5("coeff5.txt", coeff5);

 std::cout << "fit 5D done" << std::endl;
 mygetchar;



 // params3.writeToFile("params3.txt");


 // plot(c, costheta,  *test3, &data, "acceptance-test3-costheta.eps");
 // plot(c, costheta1, *test3, &data, "acceptance-test3-costheta1.eps");
 // plot(c, costheta2, *test3, &data, "acceptance-test3-costheta2.eps");
 // plotwcut(c, costheta,  *test3, &data, "costheta1neg", "acceptance-test3-costheta-costheta1neg.eps");
 // plotwcut(c, costheta,  *test3, &data, "costheta1pos", "acceptance-test3-costheta-costheta1pos.eps");
 // plotwcut(c, costheta,  *test3, &data, "costheta2neg", "acceptance-test3-costheta-costheta2neg.eps");
 // plotwcut(c, costheta,  *test3, &data, "costheta2pos", "acceptance-test3-costheta-costheta2pos.eps");

 // plotwcut(c, costheta1,  *test3, &data, "costhetaneg", "acceptance-test3-costheta1-costhetaneg.eps");
 // plotwcut(c, costheta1,  *test3, &data, "costhetapos", "acceptance-test3-costheta1-costhetapos.eps");
 // plotwcut(c, costheta1,  *test3, &data, "costheta2neg", "acceptance-test3-costheta1-costheta2neg.eps");
 // plotwcut(c, costheta1,  *test3, &data, "costheta2pos", "acceptance-test3-costheta1-costheta2pos.eps");

 // plotwcut(c, costheta2,  *test3, &data, "costheta1neg", "acceptance-test3-costheta2-costheta1neg.eps");
 // plotwcut(c, costheta2,  *test3, &data, "costheta1pos", "acceptance-test3-costheta2-costheta1pos.eps");
 // plotwcut(c, costheta2,  *test3, &data, "costhetaneg", "acceptance-test3-costheta2-costhetaneg.eps");
 // plotwcut(c, costheta2,  *test3, &data, "costhetapos", "acceptance-test3-costheta2-costhetapos.eps");


  std::cout << " new function min " << findminimum(func0, 10., -1., 1., -1., 1., -1., 1.) << std::endl;
 std::cout << " new function max " << findmaximum(func0, 0., -1., 1., -1., 1., -1., 1.) << std::endl;
 mygetchar;


 // exit(1);
 
  // RooRealVar P_b("P_b","P_{b}",1.00, -1.5, 1.5);//0.50
  // RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}",0.642,0.,1.);
  // alpha_lambda.setConstant();

  // RooRealVar alpha_b("alpha_b","#alpha_{b}",-0.45767,-1.1,1.1);//0.457
  // RooRealVar r_0("r_0","r_{0}",0.251733,-0.3,2.3);//0.5
  // RooRealVar r_1("r_1","r_{1}",0.116484,-1.3,1.3);//0.1

  // RooRealVar* xirrv[8];
  // for (int i=0; i<8; i++) xirrv[i] = new RooRealVar("xi" + istr(i),"xi" + istr(i),1.);

  // RooLbtoJpsiL0PDF3 w3("w3","w3",costheta,costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1,*xirrv[0],*xirrv[1],*xirrv[2],*xirrv[3],*xirrv[4],*xirrv[5],*xirrv[6],*xirrv[7]);

  // RooProdPdf testpdf("testpdf","testpdf",RooArgList(w3,*test3));

  // plotpdf(c,costheta,*test3,"jw-costheta-test3.eps");
  // plotpdf(c,costheta1,*test3,"jw-costheta1-test3.eps");
  // plotpdf(c,costheta2,*test3,"jw-costheta2-test3.eps");

  // plotpdf(c,costheta,testpdf,"jw-costheta-testpdf.eps");
  // plotpdf(c,costheta1,testpdf,"jw-costheta1-testpdf.eps");
  // plotpdf(c,costheta2,testpdf,"jw-costheta2-testpdf.eps");




 // TString datafilename="Tuple-Lb2L0Jpsi-Polb-ATLAS.root";
 // TFile datafile(datafilename);
 // TDirectory* datadirectory = (TDirectory*)datafile.Get("BtoqllGeneratedDistributions");
 // TTree* datatree = (TTree*)datadirectory->Get("Lb2JpsiL0");
 
 RooRealVar keep("keep","keep",0.,1.);
 RooRealVar acceptance("w","w",0.,1.);
 RooRealVar lb_id("lb_id","lb_id",-6000.,6000.);
  
 RooDataSet origdata("realdata","realdata",datatree, RooArgSet(costheta,costheta1,costheta2,phi1,phi2,lb_id,keep,acceptance) , "lb_id>0&&keep==1");
 RooDataSet* reddata = (RooDataSet*)origdata.reduce(RooFit::EventRange(1,1000));
 
 RooDataSet& realdata = *reddata;


 //fix acceptance
  for (int i=0; i<ordermax3i;i++) {
   for (int j=0; j<ordermax3j;j++) {
     for (int k=0; k<ordermax3k;k++) {
       if (cc0c1c2[i][j][k]) cc0c1c2[i][j][k]->setConstant();
     }
   }
 }

  // plot(c, costheta,  testpdf, &realdata, "blu-costheta.eps");
  // plot(c, costheta1, testpdf, &realdata, "blu-costheta1.eps");
  // plot(c, costheta2, testpdf, &realdata, "blu-costheta2.eps");


  // RooFitResult* rdfr = testpdf.fitTo(realdata,RooFit::Save(true)); 


  // rdfr->Print("v");
 mygetchar;

 double xi[8];
 xi[0]=func0.Integral(-1.,1.,-1.,1.,-1.,1.);
 xi[1]=func1.Integral(-1.,1.,-1.,1.,-1.,1.);
 xi[2]=func2.Integral(-1.,1.,-1.,1.,-1.,1.);
 xi[3]=func3.Integral(-1.,1.,-1.,1.,-1.,1.);
 xi[4]=func4.Integral(-1.,1.,-1.,1.,-1.,1.);
 xi[5]=func5.Integral(-1.,1.,-1.,1.,-1.,1.);
 xi[6]=func6.Integral(-1.,1.,-1.,1.,-1.,1.);
 xi[7]=func7.Integral(-1.,1.,-1.,1.,-1.,1.);

   // for (int i=0; i<8; i++) {
   //   xi[i]=prodpdf[i]->createIntegral(RooArgSet(costheta,costheta1,costheta2))->getVal();
   // }
 for (int i=0; i<8; i++) std::cout << "xi[" << i << "] = " << xi[i] << std::endl;

 std::cout << "func0 max = " << func0.GetMaximum() << std::endl;
 mygetchar;

 // RooArgSet tmp2(test3fr->floatParsFinal());
 // tmp2.writeToFile("test3.txt");

 
 

 // RooArgSet tmp;tmp.readFromFile("test4.txt");
 // RooArgSet tmp2;tmp2.readFromFile("test4.txt");


 // return 1;


 // TString mycut3="costhetaneg";
 // plotwcut(c, costheta,  test3, &data, mycut3, "acceptance-test3-wcut1-costheta.eps");
 // plotwcut(c, costheta1, test3, &data, mycut3, "acceptance-test3-wcut1-costheta1.eps");
 // plotwcut(c, costheta2, test3, &data, mycut3, "acceptance-test3-wcut1-costheta2.eps");

 // mycut3="costhetapos";
 // plotwcut(c, costheta,  test3, &data, mycut3, "acceptance-test3-wcut2-costheta.eps");
 // plotwcut(c, costheta1, test3, &data, mycut3, "acceptance-test3-wcut2-costheta1.eps");
 // plotwcut(c, costheta2, test3, &data, mycut3, "acceptance-test3-wcut2-costheta2.eps");



 
  return 0;



 // return 0;


 //no correlation

 RooRealVar costhc0("costhc0","costhc0",0.,-1.5,1.5);
 RooRealVar costhc1("costhc1","costhc1",0.,-1.5,1.5);
 RooRealVar costhc2("costhc2","costhc2",0.,-1.5,1.5);
 RooRealVar costhc3("costhc3","costhc3",0.,-1.5,1.5);
 RooRealVar costh1c0("costh1c0","costh1c0",0.,-1.5,1.5);
 RooRealVar costh1c1("costh1c1","costh1c1",0.,-1.5,1.5);
 RooRealVar costh1c2("costh1c2","costh1c2",0.,-1.5,1.5);
 RooRealVar costh2c0("costh2c0","costh2c0",0.,-1.5,1.5);
 RooRealVar costh2c1("costh2c1","costh2c1",0.,-1.5,1.5);
 RooRealVar costh2c2("costh2c2","costh2c2",0.,-1.5,1.5);
 RooChebychev costheta0pdf("costheta0pdf","costheta0pdf",costheta, RooArgList(costhc0,costhc1,costhc2,costhc3));
 RooChebychev costheta1pdf("costheta1pdf","costheta1pdf",costheta1,RooArgList(costh1c0,costh1c1,costh1c2));
 RooChebychev costheta2pdf("costheta2pdf","costheta2pdf",costheta2,RooArgList(costh2c0,costh2c1,costh2c2));

 RooRealVar ph1c0("ph1c0","ph1c0",0.,-1.5,1.5);
 RooRealVar ph1c1("ph1c1","ph1c1",0.,-1.5,1.5);
 RooRealVar ph1c2("ph1c2","ph1c2",0.,-1.5,1.5);
 RooRealVar ph2c0("ph2c0","ph2c0",0.,-1.5,1.5);
 RooRealVar ph2c1("ph2c1","ph2c1",0.,-1.5,1.5);
 RooRealVar ph2c2("ph2c2","ph2c2",0.,-1.5,1.5);
 RooRealVar ph2c3("ph2c3","ph2c3",0.,-1.5,1.5);
 RooRealVar ph2c4("ph2c4","ph2c4",0.,-1.5,1.5);
 RooRealVar ph2c5("ph2c5","ph2c5",0.,-1.5,1.5);
 RooChebychev phi1pdf("phi1pdf","phi1pdf",phi1,RooArgList(ph1c0,ph1c1,ph1c2));
 RooChebychev phi2pdf("phi2pdf","phi2pdf",phi2,RooArgList(ph2c0,ph2c1,ph2c2,ph2c3,ph2c4,ph2c5));
 

 RooProdPdf ncpdf("ncpdf","ncpdf",RooArgList(costheta0pdf,costheta1pdf,costheta2pdf,phi1pdf,phi2pdf));

 // RooRealVar costhc4("costhc4","costhc4",0.,-1.5,1.5);
 // RooRealVar costhc5("costhc5","costhc5",0.,-1.5,1.5);
 // RooChebychev costheta0pdf2("costheta0pdf2","costheta0pdf2",costheta, RooArgList(costhc4,costhc5));


 // RooRealVar frac("frac","frac",0.5,-0.1,1.1);
 // RooAddPdf model("model","model",RooArgList(ncpdf,costheta0pdf2),RooArgList(frac));

 RooAbsPdf& model = ncpdf;


 // RooDataSet* data = test.generate(RooArgList(costheta,costheta1,costheta2,phi1,phi2),10);

 //

 model.fitTo(data);


 plot(c, costheta,  model, &data, "acceptance-costheta.eps");
 plot(c, costheta1, model, &data, "acceptance-costheta1.eps");
 plot(c, costheta2, model, &data, "acceptance-costheta2.eps");
 plot(c, phi1,      model, &data, "acceptance-phi1.eps");
 plot(c, phi2,      model, &data, "acceptance-phi2.eps");

 TString mycut="costhetaneg";
 plot(c, costheta,  model, &data, "acceptance-wcut1-costheta.eps",false,NULL,mycut);
 plot(c, costheta1, model, &data, "acceptance-wcut1-costheta1.eps",false,NULL,mycut);
 plot(c, costheta2, model, &data, "acceptance-wcut1-costheta2.eps",false,NULL,mycut);
 plot(c, phi1,      model, &data, "acceptance-wcut1-phi1.eps",false,NULL,mycut);
 plot(c, phi2,      model, &data, "acceptance-wcut1-phi2.eps",false,NULL,mycut);

 mycut="costhetapos";
 plot(c, costheta,  model, &data, "acceptance-wcut2-costheta.eps",false,NULL,mycut);
 plot(c, costheta1, model, &data, "acceptance-wcut2-costheta1.eps",false,NULL,mycut);
 plot(c, costheta2, model, &data, "acceptance-wcut2-costheta2.eps",false,NULL,mycut);
 plot(c, phi1,      model, &data, "acceptance-wcut2-phi1.eps",false,NULL,mycut);
 plot(c, phi2,      model, &data, "acceptance-wcut2-phi2.eps",false,NULL,mycut);

 



}
