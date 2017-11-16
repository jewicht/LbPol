#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TF3.h"
#include "TH3F.h"
#include "TCut.h"
#include "TProfile.h"

#include "RooRealVar.h"
#include "RooNDKeysPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooChebychev.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"

#include "rootlogon.h"
#include "functions-roofit.h"
#include "createRooLegendre5.h"
#include "cuts.h"
#include "calcacceptanceclass5.h"

RooRealVar costheta("costheta",   "cos(#theta)",     -1., 1.);
RooRealVar costheta1("costheta1", "cos(#theta_{1})", -1., 1.);
RooRealVar costheta2("costheta2", "cos(#theta_{2})", -1., 1.);
RooRealVar phi1("phi1", "#phi_{1}/#pi", -1., 1.);
RooRealVar phi2("phi2", "#phi_{2}/#pi", -1., 1.);



void plotprofile(const TCanvas& c, TTree* sigtree, const RooArgList& varlist, const TCut& selection , const TString& plotname)
{
  for (int i=0; i<varlist.getSize(); i++) {
    for (int j=i+1; j<varlist.getSize(); j++) {
      RooAbsArg* var1 = varlist.at(i);
      RooAbsArg* var2 = varlist.at(j);  
      TProfile prof("prof","prof",25,-1.,1.);
      sigtree->Draw(TString(var1->GetName())+":"+var2->GetName() + ">>prof", selection);
      prof.Fit("pol1");
      prof.GetXaxis()->SetTitle(TString(var1->GetTitle())+" vs "+var2->GetTitle());
      prof.Draw();
      c.SaveAs(plotname + "-" + TString(var1->GetName()) + "-vs-" + var2->GetName() +".eps");
    }
  }
}

RooAbsPdf* factorizedacc(int a, int b, int c, int d, int e) 
{
  RooArgList costhetacoeff;
  RooArgList costheta1coeff;
  RooArgList costheta2coeff;
  RooArgList phi1coeff;
  RooArgList phi2coeff;
  for (int i=0;i<a;i++) costhetacoeff.add(*(RooRealVar*)new RooRealVar("coeffa"+istr(i),"coeffa"+istr(i), 0., -1., 1. ));
  for (int i=0;i<b;i++) costheta1coeff.add(*(RooRealVar*)new RooRealVar("coeffb"+istr(i),"coeffb"+istr(i), 0., -1., 1. ));
  for (int i=0;i<c;i++) costheta2coeff.add(*(RooRealVar*)new RooRealVar("coeffc"+istr(i),"coeffc"+istr(i), 0., -1., 1. ));
  for (int i=0;i<d;i++) phi1coeff.add(*(RooRealVar*)new RooRealVar("coeffd"+istr(i),"coeffd"+istr(i), 0., -1., 1. ));
  for (int i=0;i<e;i++) phi2coeff.add(*(RooRealVar*)new RooRealVar("coeffe"+istr(i),"coeffe"+istr(i), 0., -1., 1. ));

  return (RooAbsPdf*)new RooProdPdf("model","model",RooArgList(*(RooChebychev*)new RooChebychev("polc0","polc0",costheta, costhetacoeff),
							       *(RooChebychev*)new RooChebychev("polc1","polc1",costheta1, costheta1coeff),
							       *(RooChebychev*)new RooChebychev("polc2","polc2",costheta2, costheta2coeff),
							       *(RooChebychev*)new RooChebychev("polc3","polc3",phi1, phi1coeff),
							       *(RooChebychev*)new RooChebychev("polc4","polc4",phi2, phi2coeff))
				    );
}


int main(int argc, char** argv)
{

  
  //  calcacceptanceclass5 acc("Lb-cc0c1c2p1p2-ll.txt");
  //  return 1;

  if (argc!=3) {
    std::cerr << argv[0] << " Lb/B0 LL/DD" << std::endl;
    exit(2);
  }

  const TString mode(argv[1]);
  const bool doBd2JpsiKS0=mode.Contains("B0");
  
  const TString llordd(argv[2]);
  const bool doLL=llordd.Contains("LL");

  TString Lambdabname, Lambda0name, Jpsiname;

  varnames(doBd2JpsiKS0, Lambdabname, Lambda0name, Jpsiname);

  lhcbstyle();
  TCanvas c("c","c",100,100);

  costheta.setRange("costhetapos",0.,1.);
  costheta.setRange("costhetaneg",-1.,0.);
  costheta1.setRange("costheta1pos",0.,1.);
  costheta1.setRange("costheta1neg",-1.,0.);
  costheta2.setRange("costheta2pos",0.,1.);
  costheta2.setRange("costheta2neg",-1.,0.);
  phi1.setRange("phi1pos",0.,1.);
  phi1.setRange("phi1neg",-1.,0.);
  phi2.setRange("phi2pos",0.,1.);
  phi2.setRange("phi2neg",-1.,0.);

  TString sigfilename("DVTuples-Lb2JpsiL0-MC11a.root");
  TString dirname("Tuple_JpsiL0_betas");
  if (doBd2JpsiKS0) {
    sigfilename="DVTuples-Bd2JpsiKS-MC11a.root";
    dirname="Tuple_JpsiKs_detached_betas";
  }

  TFile sigfile(sigfilename);
  TDirectory* sigdir = (TDirectory*)sigfile.Get(dirname);
  TTree* sigtree = (TTree*)sigdir->Get("DecayTree");

  TString siganglesfilename(sigfilename);
  sigtree->AddFriend(dirname, TString(sigfilename).ReplaceAll(".root", ".angles.root"));


  TString coeff;
  TString plotname("fitacceptance5-");

  if (doBd2JpsiKS0) {
    coeff+="B0-";
    plotname+="B0-";
  } else {
    coeff+="Lb-";
    plotname+="Lb-";
  }
  coeff+="cc0c1c2p1p2-";
  if (doLL) {
    coeff+="ll.txt";
    plotname+="ll";
  } else {
    coeff+="dd.txt";
    plotname+="dd";
  }


  //cut based selection
  TCut LL_selection;
  TCut DD_selection;

  simplecut(doBd2JpsiKS0, LL_selection, DD_selection);
  
  TCut TrueLambdab=truemcsel(doBd2JpsiKS0);
  TCut TriggerSel=triggercut(doBd2JpsiKS0);
  LL_selection+=TrueLambdab;
  DD_selection+=TrueLambdab;
  LL_selection+=TriggerSel;
  DD_selection+=TriggerSel;

  RooArgList varlist(costheta,costheta1,costheta2,phi1,phi2);
  if (doLL) {
    plotprofile(c,sigtree,varlist, DD_selection, plotname + "-profile");
  } else {
    plotprofile(c,sigtree,varlist, DD_selection, plotname + "-profile");
  }

  RooDataSet* data; 
  if (doLL) {
    data = filldataset("mcll","mcll", RooArgList(costheta,costheta1,costheta2,phi1,phi2), sigtree, LL_selection);
  } else {
    data = filldataset("mcdd","mcdd", RooArgList(costheta,costheta1,costheta2,phi1,phi2), sigtree, DD_selection);
  }

  costheta.setBins(25);
  costheta1.setBins(25);
  costheta2.setBins(25);
  phi1.setBins(25);
  phi2.setBins(25);

  RooAbsPdf* model = factorizedacc(4,4,4,4,7);

  model->fitTo(*data);
  plot(c,RooArgList(costheta,costheta1,costheta2,phi1,phi2),*model,data, plotname + "-factorpdf");

  //  return 1;


  RooLegendre5Pdfv2* corrpdf=NULL;
  RooRealVar* coeff5[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m];

  //  7,5,3,
  if (doBd2JpsiKS0) {
    createRooLegendre5v2( costheta,costheta1,costheta2,phi1,phi2, "test", corrpdf,coeff5, 7,5,3,4,7,8,2);
  } else {
    createRooLegendre5v2( costheta,costheta1,costheta2,phi1,phi2, "test", corrpdf,coeff5, 4,4,4,9,9,10,2);
  }
  //  RooDataHist mcll_binned("mcll_binned","mcll_binned",RooArgList(costheta,costheta1,costheta2,phi1,phi2),*mcll);
  //  setconstcoeff5(coeff5);
  //  setconstcoeff5_maxtotorder(coeff5,7);
  
  readcoeff5_rrv(coeff,coeff5);
  std::cout << data->numEntries() << std::endl;
  //   RooDataSet* mcllred = (RooDataSet*)mcll->reduce(RooFit::EventRange(1,1000));

  mygetchar;
  RooFitResult* frll = corrpdf->fitTo(*data, RooFit::Save(true),RooFit::Minos(false), RooFit::Hesse(false));//, RooFit::NumCPU(4));
  writecoeff5(coeff, coeff5);
  mygetchar;
  plot(c,RooArgList(costheta,costheta1,costheta2,phi1,phi2),*corrpdf,data,plotname + "-corrpdf");
  plot(c, RooArgList(costheta1,costheta2,phi1,phi2), *corrpdf, data, plotname + "-corrpdf" + "-proj",false,NULL,"costhetaneg");
  plot(c, RooArgList(costheta1,costheta2,phi1,phi2), *corrpdf, data, plotname + "-corrpdf" + "-proj",false,NULL,"costhetaneg");
  plot(c, RooArgList(costheta,costheta2,phi1,phi2), *corrpdf, data, plotname + "-corrpdf" + "-proj",false,NULL,"costheta1neg");
  plot(c, RooArgList(costheta,costheta2,phi1,phi2), *corrpdf, data, plotname + "-corrpdf" + "-proj",false,NULL,"costheta1pos");
  plot(c, RooArgList(costheta,costheta1,phi1,phi2), *corrpdf, data, plotname + "-corrpdf" + "-proj",false,NULL,"costheta2neg");
  plot(c, RooArgList(costheta,costheta1,phi1,phi2), *corrpdf, data, plotname + "-corrpdf" + "-proj",false,NULL,"costheta2pos");
  plot(c, RooArgList(costheta,costheta1,costheta2,phi2), *corrpdf, data, plotname + "-corrpdf" + "-proj",false,NULL,"phi1neg");
  plot(c, RooArgList(costheta,costheta1,costheta2,phi2), *corrpdf, data, plotname + "-corrpdf" + "-proj",false,NULL,"phi1pos");
  plot(c, RooArgList(costheta,costheta1,costheta1,phi1), *corrpdf, data, plotname + "-corrpdf" + "-proj",false,NULL,"phi2neg");
  plot(c, RooArgList(costheta,costheta1,costheta1,phi1), *corrpdf, data, plotname + "-corrpdf" + "-proj",false,NULL,"phi2pos");
  frll->Print("v");
  mygetchar;

  return 0;

  //  RooHistPdf  histpdf_ll("histpdf_ll","histpdf_ll",RooArgList(costheta,costheta1,costheta2,phi1,phi2),mcll_binned,3);
  //  plot(c,RooArgList(costheta,costheta1,costheta2,phi1,phi2),histpdf_ll,mcll,"fitacceptance5-ll-histpdf");
  

//   std::cout << "Making ll keys pdf out of Nevents = " << mcll->numEntries() << std::endl;
//   RooNDKeysPdf keyspdf3_ll("keyspdf3_ll","keyspdf3_ll", RooArgList(costheta,costheta1,costheta2), *mcll);
//   plot(c,RooArgList(costheta,costheta1,costheta2),keyspdf3_ll,mcll,"fitacceptance5-ll-keys3");
//   RooNDKeysPdf keyspdf_ll("keyspdf_ll","keyspdf_ll", RooArgList(costheta,costheta1,costheta2,phi1,phi2), *mcll);
//   plot(c,RooArgList(costheta,costheta1,costheta2,phi1,phi2),keyspdf_ll,mcll,"fitacceptance-ll-keys5");

//   std::cout << "Making dd keys pdf out of Nevents = " << mcdd->numEntries() << std::endl;
//   RooNDKeysPdf keyspdf_dd("keyspdf_dd","keyspdf_dd", RooArgList(costheta,costheta1,costheta2,phi1,phi2), *mcdd);
//   plot(c,RooArgList(costheta,costheta1,costheta2,phi1,phi2),keyspdf_dd,mcdd,"fitacceptance5-dd-keys");


//   return 0;
  
}
