#include <iostream>
#include "math.h"

#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooCBShape.h"

#include "functions.h"
#include "functions-roofit.h"

#include "rootlogon.h"

#include "cuts.h"

int main(int argc, char** argv)
{

  lhcbstyle();

  const TString mode(argv[1]);
  if (mode!="B0" and mode!="Lb") exit(1);//usage(argv[0]);
  const bool doBd2JpsiKS0=mode.Contains("B0");

  TString dirname("Tuple_JpsiL0_betas");
  if (doBd2JpsiKS0) {
    dirname="Tuple_JpsiKs_detached_betas";
  }    

  //  const TString sigfilename("DVTuples-Lb2JpsiL0-MC11a.root");
  //   TFile sigfile(sigfilename);
  //   TDirectory* sigdirectory = (TDirectory*)sigfile.Get(dirname);
  //   TTree* sigtree = (TTree*)sigdirectory->Get("DecayTree");
  // //  sigtree->AddFriend(dirname, TString(sigfilename).ReplaceAll(".root",".bdt.root"));
  
  const TString bkgfilename("DVTuples_data_stripping17.reduced.root");
  TFile bkgfile(bkgfilename);
  TDirectory* bkgdirectory = (TDirectory*)bkgfile.Get(dirname);
  TTree* bkgtree = (TTree*)bkgdirectory->Get("DecayTree");
//  bkgtree->AddFriend(dirname, TString(bkgfilename).ReplaceAll(".root",".bdt.root"));


  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0, Lambdabname, Lambda0name, Jpsiname);
  
  const TCut TriggerCut = triggercut(doBd2JpsiKS0);  
  //  const TCut pid        = pidcut(false);
  const TCut cutall(TriggerCut + TCut(Lambdabname + "_DIRA_OWNPV>0.99&&" + Lambdabname + "_TAU>0.&&" + Lambdabname + "_IPCHI2_OWNPV<30."));
  const TCut cutsignal = truemcsel(doBd2JpsiKS0);//("abs(1.*lambda_b0_TRUEID)==5122");
  const TCut cutbkg(TString(Lambdabname + "_M>5800.&&" + Lambdabname + "_M<6000."));

  TCanvas c("c","c",100,100);


  RooRealVar mass(Lambdabname + "_M","M(#LambdaJ/#psi)",5550.,5700.,"MeV/c^{2}");

  if (doBd2JpsiKS0) {
    mass.SetTitle("M(J/#psiK_{s}^{0})");
    mass.setMin(5230.);    
    mass.setMax(5450.);
  }
  const TCut masscut(TString(Lambdabname + "_M>" + str(mass.getMin()) + "&&" + Lambdabname + "_M<"+ str(mass.getMax())));

  mass.setBins(50);
  mass.setRange("signalregion", 5600.,5640.);

  RooRealVar nSig("nSig","nSig",2000.,0.,300000.);
  //  RooRealVar nSig2("nSig2","nSig2",2000.,0.,5000.);
  RooRealVar nBkg("nBkg","nBkg",10000.,0.,1000000.);

  RooRealVar sig_dd_mean("sig_dd_mean","sig_dd_mean",5620.,5500.,5700.);
  RooRealVar sig_ll_mean("sig_ll_mean","sig_ll_mean",5620.,5500.,5700.);
  // if (doBd2JpsiKS0) {
  //   sig_dd_mean.setVal(5280.);
  //   sig_dd_mean.setMin(5270.);
  //   sig_dd_mean.setMax(5290.);
  //   sig_ll_mean.setVal(5280.);
  //   sig_ll_mean.setMin(5270.);
  //   sig_ll_mean.setMax(5290.);
  // }

  if (doBd2JpsiKS0) {
    sig_dd_mean=5280.;
    sig_dd_mean.setMin(5260.);
    sig_dd_mean.setMax(5300.);
    sig_ll_mean=5280.;
    sig_ll_mean.setMin(5260.);
    sig_ll_mean.setMax(5300.);
  }

  RooRealVar sig_dd_sigma("sig_dd_sigma","sig_dd_sigma",7.5,0.,100.);
  RooRealVar sig_ll_sigma("sig_ll_sigma","sig_ll_sigma",7.5,0.,100.);

  RooRealVar sig_dd_frac("sig_dd_frac","sig_dd_frac", 0.5, -0.1, 1.1);
  RooRealVar sig_ll_frac("sig_ll_frac","sig_ll_frac", 0.5, -0.1, 1.1);

  RooRealVar sig_dd_alpha1("sig_dd_alpha1","sig_dd_alpha1", 1.,0.,10.);
  RooRealVar sig_dd_alpha2("sig_dd_alpha2","sig_dd_alpha2",-1.,-10.,0.);

  RooRealVar sig_dd_n1("sig_dd_n1","sig_dd_n1",10.,0.,200.);
  RooRealVar sig_dd_n2("sig_dd_n2","sig_dd_n2",10.,0.,200.);

  RooRealVar sig_ll_alpha1("sig_ll_alpha1","sig_ll_alpha1", 1.,0.,10.);
  RooRealVar sig_ll_alpha2("sig_ll_alpha2","sig_ll_alpha2",-1.,-10.,0.);

  RooRealVar sig_ll_n1("sig_ll_n1","sig_ll_n1",10.,0.,200.);
  RooRealVar sig_ll_n2("sig_ll_n2","sig_ll_n2",10.,0.,200.);
  
  RooCBShape sig_dd_cbshape1("sig_dd_cbshape1", "sig_dd_cbshape1", mass, sig_dd_mean, sig_dd_sigma, sig_dd_alpha1, sig_dd_n1);//, RooAbsReal& _n)
  RooCBShape sig_dd_cbshape2("sig_dd_cbshape2", "sig_dd_cbshape2", mass, sig_dd_mean, sig_dd_sigma, sig_dd_alpha2, sig_dd_n2);//, RooAbsReal& _n)
  RooAddPdf sig_dd_twocb("sig_dd_twocb","sig_dd_twocb",RooArgList(sig_dd_cbshape1,sig_dd_cbshape2),RooArgList(sig_dd_frac));
  RooCBShape sig_ll_cbshape1("sig_ll_cbshape1", "sig_ll_cbshape1", mass, sig_ll_mean, sig_ll_sigma, sig_ll_alpha1, sig_ll_n1);//, RooAbsReal& _n)
  RooCBShape sig_ll_cbshape2("sig_ll_cbshape2", "sig_ll_cbshape2", mass, sig_ll_mean, sig_ll_sigma, sig_ll_alpha2, sig_ll_n2);//, RooAbsReal& _n)
  RooAddPdf sig_ll_twocb("sig_ll_twocb","sig_ll_twocb",RooArgList(sig_ll_cbshape1,sig_ll_cbshape2),RooArgList(sig_ll_frac));

  const RooArgList sig_dd_vars(sig_dd_frac,sig_dd_mean,sig_dd_sigma,sig_dd_alpha1,sig_dd_alpha2,sig_dd_n1,sig_dd_n2);
  const RooArgList sig_dd_varsnottoconst(sig_dd_mean,sig_dd_sigma);
  //  const RooArgList mass_dd_pdfs(sig_dd_mass_model,bkg_dd_mass_model,bkg_dd_Lbbar_mass_model);
  
  const RooArgList sig_ll_vars(sig_ll_frac,sig_ll_mean,sig_ll_sigma,sig_ll_alpha1,sig_ll_alpha2,sig_ll_n1,sig_ll_n2);
  const RooArgList sig_ll_varsnottoconst(sig_ll_mean,sig_ll_sigma);
  //  const RooArgList mass_ll_pdfs(sig_ll_mass_model,bkg_ll_mass_model,bkg_ll_Lbbar_mass_model);
  
  
  if (doBd2JpsiKS0) {
  readvarsfromtfile(sig_ll_vars, TString("fitmcmass-Bd2JpsiKS0-") + dirname + "-mc-ll-fr.root", "", false);
  readvarsfromtfile(sig_dd_vars, TString("fitmcmass-Bd2JpsiKS0-") + dirname + "-mc-dd-fr.root", "", false);//Krandsigparams and dorand!=0);
} else {
  readvarsfromtfile(sig_dd_vars, TString("fitmcmass-Lb2JpsiL0-") + dirname + "-mc-dd-fr.root", "", false);//randsigparams and dorand!=0);
  readvarsfromtfile(sig_ll_vars, TString("fitmcmass-Lb2JpsiL0-") + dirname + "-mc-ll-fr.root", "", false);//randsigparams and dorand!=0);
  }

  sig_dd_frac.setConstant();
  sig_ll_frac.setConstant();

  sig_dd_alpha1.setConstant();
  sig_dd_alpha2.setConstant();

  sig_dd_n1.setConstant();
  sig_dd_n2.setConstant();

  sig_ll_alpha1.setConstant();
  sig_ll_alpha2.setConstant();

  sig_ll_n1.setConstant();
  sig_ll_n2.setConstant();
  
  RooAbsPdf& sig_dd_mass_model = sig_dd_twocb;
  RooAbsPdf& sig_ll_mass_model = sig_ll_twocb;


  RooRealVar c0("c0","c0",0.,-1.1,1.1);
  RooChebychev bkg_mass_model("bkg_mass_model","bkg_mass_model",mass,RooArgList(c0));
  RooAddPdf tot_dd_mass_model("tot_dd_mass_model","tot_dd_mass_model",RooArgList(sig_dd_mass_model,bkg_mass_model),RooArgList(nSig,nBkg));
  RooAddPdf tot_ll_mass_model("tot_ll_mass_model","tot_ll_mass_model",RooArgList(sig_ll_mass_model,bkg_mass_model),RooArgList(nSig,nBkg));


  RooDataSet* datadd = filldataset("datadd","datadd",RooArgList(mass),bkgtree,cutall + masscut + Lambda0_DD);

  RooFitResult* frdd = tot_dd_mass_model.fitTo(*datadd,RooFit::Save(true));
  RooArgSet floatpars(frdd->floatParsFinal()); floatpars.writeToFile("optbdt-dd-float.txt");
  RooArgSet constpars(frdd->constPars()); constpars.writeToFile("optbdt-dd-const.txt");
  RooAbsReal* bkgddint    = bkg_mass_model.createIntegral(RooArgSet(mass));
  RooAbsReal* bkgddint_sr = bkg_mass_model.createIntegral(RooArgSet(mass),"signalregion");
  RooAbsReal* sigddint = sig_dd_mass_model.createIntegral(RooArgSet(mass));
  RooAbsReal* sigddint_sr = sig_dd_mass_model.createIntegral(RooArgSet(mass),"signalregion");
  std::cout << "bkgddint    = " << bkgddint->getVal() << std::endl;
  std::cout << "bkgddint_sr = " << bkgddint_sr->getVal() << std::endl;
  std::cout << "nSig      = " << nSig.getVal() << std::endl;
  std::cout << "nBkg      = " << nBkg.getVal() << std::endl;
  std::cout << "nSig_sr   = " << nSig.getVal()*sigddint_sr->getVal()/sigddint->getVal() << std::endl;
  std::cout << "nBkg_sr   = " << nBkg.getVal()*bkgddint_sr->getVal()/bkgddint->getVal() << std::endl;

  plot(c,mass,tot_dd_mass_model,datadd,"optbdt-dd-"+mode);
  
  RooDataSet* datall = filldataset("datall","datall",RooArgList(mass),bkgtree,cutall + masscut + Lambda0_LL);
  RooFitResult* frll = tot_ll_mass_model.fitTo(*datall,RooFit::Save(true));
  floatpars = (frll->floatParsFinal()); floatpars.writeToFile("optbdt-ll-float.txt");
  constpars = (frll->constPars()); constpars.writeToFile("optbdt-ll-const.txt");
  RooAbsReal* bkgllint = bkg_mass_model.createIntegral(RooArgSet(mass));
  RooAbsReal* bkgllint_sr = bkg_mass_model.createIntegral(RooArgSet(mass),"signalregion");
  RooAbsReal* sigllint = sig_ll_mass_model.createIntegral(RooArgSet(mass));
  RooAbsReal* sigllint_sr = sig_ll_mass_model.createIntegral(RooArgSet(mass),"signalregion");
  std::cout << "bkgllint    = " << bkgllint->getVal() << std::endl;
  std::cout << "bkgllint_sr = " << bkgllint_sr->getVal() << std::endl;
  std::cout << "nSig      = " << nSig.getVal() << std::endl;
  std::cout << "nBkg      = " << nBkg.getVal() << std::endl;
  std::cout << "nSig_sr   = " << nSig.getVal()*sigllint_sr->getVal()/sigllint->getVal() << std::endl;
  std::cout << "nBkg_sr   = " << nBkg.getVal()*bkgllint_sr->getVal()/bkgllint->getVal() << std::endl;
  plot(c,mass,tot_ll_mass_model,datall,"optbdt-ll-"+mode);
  
  std::cout << "DD fit" << std::endl;
  frdd->Print("v");

  std::cout << "LL fit" << std::endl;
  frll->Print("v");

//LL: 8x1600 (12800) bkg in signal region, 2314 sig
//best cut -0.0281
//DD: 8x6000 bkg, 6950 sig
//best cut: -0.0202
  return 0;


  // bkgtree->Draw(Lambdabname + "_M" , cutall + cutbkg);
  // c.SaveAs("c.eps");

  // const double luminosity = 1.1;//fb-1
  // const double bbcs       = 0.45445E12 ; //fb

  // const double decprodcuteff = 0.16;
  // const int nSigTot    = 532997 + 528298;
  // const double bf = 4.7E-05 * 0.0593 * 0.639;
  // const double f_Lb = 0.20;

  // const double N_Lambdab = luminosity * bbcs * f_Lb * 2.;

  // std::cout << "N(Lambdab) = " << N_Lambdab << std::endl; 
  // std::cout << "N(Lambdab to Jpsi L0) = " << N_Lambdab * bf << std::endl;
  // std::cout << "N(Lambdab to Jpsi L0) in acceptance = " << N_Lambdab * bf * decprodcuteff << std::endl;

  // std::cout << "N(Lambdab to Jpsi L0) reconstructable = " << N_Lambdab * bf * decprodcuteff * sigtree->GetEntries(cutall) / nSigTot << std::endl;

  
  // const double bdtmin=-0.5;
  // const double bdtmax=0.5;
  // const int nstep=20;

  // const double bkgnorm=1.25*40./200.;

  // for (int i=0; i<nstep; i++) {
  //   const double bdtcutval = bdtmin + (bdtmax-bdtmin) * i/nstep;
  //   const TCut bdtcut("bdt>" + str(bdtcutval));
  //   const double nsig = N_Lambdab * bf * decprodcuteff * sigtree->GetEntries(cutall + bdtcut + Lambda0_LL) / nSigTot;
  //   const double nbkg = bkgtree->GetEntries(cutall + cutbkg + bdtcut + Lambda0_LL) * bkgnorm;

  //   const double signif = nsig / sqrt(nsig+nbkg);

  //   std::cout << "bdt > " << bdtcutval << std::endl;
  //   std::cout << "nsig = " << nsig << std::endl;
  //   std::cout << "nbkg = " << nbkg << std::endl;
  //   std::cout << "signif = " << signif << std::endl;
  // }


  // for (int i=0; i<10; i++) std::cout << std::endl;

  // for (int i=0; i<nstep; i++) {
  //   const double bdtcutval = bdtmin + (bdtmax-bdtmin) * i/nstep;
  //   const TCut bdtcut("bdt>" + str(bdtcutval));
  //   const double nsig = N_Lambdab * bf * decprodcuteff * sigtree->GetEntries(cutall + bdtcut + Lambda0_DD) / nSigTot;
  //   const double nbkg = bkgtree->GetEntries(cutall + cutbkg + bdtcut + Lambda0_DD) * bkgnorm;

  //   const double signif = nsig / sqrt(nsig+nbkg);

  //   std::cout << "bdt > " << bdtcutval << std::endl;
  //   std::cout << "nsig = " << nsig << std::endl;
  //   std::cout << "nbkg = " << nbkg << std::endl;
  //   std::cout << "signif = " << signif << std::endl;
  // }

}
