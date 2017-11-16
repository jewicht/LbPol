#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"

#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooProdPdf.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"

#include "rootlogon.h"
#include "functions-roofit.h"
#include "cuts.h"

RooRealVar mass("M","M(#LambdaJ/#psi)",5500.,5750.,"MeV/c^{2}");
RooRealVar ksmass("ksmass","M(#Lambda)",1000.,1200.,"MeV/c^{2}");
//RooLbtoJpsiL0wAccPDF3* w3wacc = NULL;

RooRealVar nSig("nSig","nSig",2000.,0.,300000.);
//  RooRealVar nSig2("nSig2","nSig2",2000.,0.,5000.);
RooRealVar nBkg("nBkg","nBkg",10000.,0.,1000000.);

RooCmdArg latexformat = RooFit::Format("NE",RooFit::AutoPrecision(2),RooFit::VerbatimName(true));

const bool usebdt=true;
const bool writeresult=false;


void simplemassfit(TCanvas& c, RooAbsPdf& model, RooDataSet* data, const TString& plotname)
{
  RooFitResult* fr = model.fitTo(*data,RooFit::Minos(false),RooFit::Save(true));
  std::cout << plotname << " mass fit (status=" << fr->status() << "; covQual=" << fr->covQual() << "; edm=" << fr->edm() << ")" << std::endl;
  if (fr->status()!=0 or fr->covQual()!=3) {
    std::cout << plotname << " mass fit - not successful (status=" << fr->status() << "; covQual=" << fr->covQual() << "; edm=" << fr->edm() << ")" << std::endl;
    mygetchar;
  }
  isconverged(fr,true);
  fr->Print("v");
  plot(c,mass,model,data, plotname + "-mass");
  plot(c,mass,model,data, plotname + "-masswpull",true);

  if (writeresult) {
    fr->SaveAs(plotname + "-fr.root");
    
    RooArgSet floatpars(fr->floatParsFinal());
    floatpars.writeToFile(plotname + ".txt");
    fr->floatParsFinal().printLatex( latexformat, RooFit::OutputFile(plotname + "-float.tex" ), RooFit::Columns(1) );
    //    myprintLatex(fr->floatParsFinal(), plotname + "-float.tex" , 2, &latexformat);
  }
  delete fr;
  mygetchar;
}

int main(int argc, char** argv)
{

  if (argc!=3) {
    std::cerr << argv[0] << " Lb/B0 LL/DD" << std::endl;
    exit(2);
  }
  const TString mode(argv[1]);
  const bool doBd2JpsiKS0=mode.Contains("B0");

  TString llordd(argv[2]);
  llordd.ToLower();
  const bool doLL=llordd.Contains("ll");

  lhcbstyle();
  TCanvas c("c","c",100,100);
  c.SetLogy();

  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0,Lambdabname,Lambda0name,Jpsiname);
  
  mass.SetName(Lambdabname + "_M");
  mass.setUnit("MeV/c^{2}");

  ksmass.SetName(Lambda0name + "_M");
  ksmass.setUnit("MeV/c^{2}");

  if (doBd2JpsiKS0) {
    mass.setMin(5150.);
    mass.setMax(5400.);
    mass.SetTitle("M(J/#psiK_{S}^{0})");
    ksmass.setMin(480.);
    ksmass.setMax(520.);
    ksmass.SetTitle("M_{#pi#pi}");
  } else {
    mass.setMin(5500.);
    mass.setMax(5750.);
    ksmass.setMin(1100.);
    ksmass.setMax(1130.);
    ksmass.SetTitle("M_{p#pi}");
  }


  RooRealVar piminus_TRACK_Type("piminus_TRACK_Type","piminus_TRACK_Type",0.,10.);


  std::vector<std::string> directorylist;
  if (doBd2JpsiKS0) {
    directorylist.push_back("Tuple_JpsiKs_detached_betas");
    //    directorylist.push_back("Tuple_JpsiKs_prescaled_betas");
    //    directorylist.push_back("Tuple_JpsiKs_detached_betas_refitted");
    //    directorylist.push_back("Tuple_JpsiKs_B2XMuMu");
    //    directorylist.push_back("Tuple_JpsiKs_B2XMuMu_refitted");
  } else {
    directorylist.push_back("Tuple_JpsiL0_betas");
    //    directorylist.push_back("Tuple_JpsiL0_betas_refitted");
    //    directorylist.push_back("Tuple_JpsiL0_B2XMuMu");
    //    directorylist.push_back("Tuple_JpsiL0_B2XMuMu_refitted");
  }

  TTree* mctree[directorylist.size()];
  memset(mctree,0,directorylist.size()*sizeof(TTree*));

  //  TString mcfilename("DVTuples-Lb2JpsiL0-MC11a-reco12a.reduced.root");
  TString mcfilename("DVTuples-Lb2JpsiL0-MC11a.reduced.root");
  if (doBd2JpsiKS0) {
    mcfilename = "DVTuples-Bd2JpsiKS-MC11a.reduced.root";
  }
  TFile mcfile(mcfilename);
  for (unsigned int i=0;i<directorylist.size();i++) {
    const TString& dirname = directorylist[i];
    TDirectory* mcdir = (TDirectory*)mcfile.Get(dirname);
    if (!mcdir) continue;
    mctree[i] = (TTree*)mcdir->Get("DecayTree");
    if (usebdt) mctree[i]->AddFriend(dirname, TString(mcfilename).ReplaceAll(".root",".bdt.root"));
    //    mctree[i]->AddFriend(dirname, mcanglesfilename);
    //    if (userwmc) mctree[i]->AddFriend(dirname, mcrwfilename);
  }

  TCut LL_selection;
  TCut DD_selection;
  simplecut2(doBd2JpsiKS0, LL_selection,DD_selection);

  //bdt selection
  if (usebdt) {
    bdtsel(doBd2JpsiKS0, LL_selection,DD_selection);
  }

  const TCut mcut=masscut(doBd2JpsiKS0, mass);
  const TCut TrueMC=truemcsel(doBd2JpsiKS0);
  const TCut TriggerCut=triggercut(doBd2JpsiKS0);

  std::cout << "DD trigger eff: " << 100.*mctree[0]->GetEntries(DD_selection + TrueMC + mcut + TriggerCut)/mctree[0]->GetEntries(DD_selection + TrueMC + mcut) << std::endl;
  std::cout << "LL trigger eff: " << 100.*mctree[0]->GetEntries(LL_selection + TrueMC + mcut + TriggerCut)/mctree[0]->GetEntries(LL_selection + TrueMC + mcut) << std::endl;

  //  return 1;

  LL_selection+=mcut;
  DD_selection+=mcut;
  LL_selection+=TrueMC;
  DD_selection+=TrueMC;
  LL_selection+=TriggerCut;
  DD_selection+=TriggerCut;

  const TCut selection = LL_selection or DD_selection;

  RooDataSet* mcdata[directorylist.size()];
  RooDataSet* mcdd[directorylist.size()];
  RooDataSet* mcll[directorylist.size()];

  RooArgList vars(mass,ksmass,piminus_TRACK_Type);
  for (unsigned int i=0; i<directorylist.size(); i++) {
    mcdata[i]=filldataset("mcdata-" + directorylist[i],"mcdata-" + directorylist[i], vars, mctree[i], selection);
    mcdata[i]->Print("v");
  }

  for (unsigned int i=0; i<directorylist.size(); i++) {
    if (!mcdata[i]) continue;
    mcdd[i] = (RooDataSet*)mcdata[i]->reduce(Lambda0_DD); 
    mcll[i] = (RooDataSet*)mcdata[i]->reduce(Lambda0_LL); 
  }

  RooRealVar mean("sig_" + llordd + "_" + "mean","sig_" + llordd + "_" +"mean",5620.,5500.,5700.);
  //  RooRealVar mean3("mean3","mean3",5620.,5500.,5700.);
  //  RooRealVar mean1("mean1","mean1",5250.,5000.,5700.);
  if (doBd2JpsiKS0) {
    mean.setVal(5280.);
    mean.setMin(5260.);
    mean.setMax(5300.);
    //    mean3.setVal(5290.);
    //    mean3.setMin(5260.);
    //    mean3.setMax(5320.);
  }

  RooRealVar frac("sig_" + llordd + "_" +"frac","sig_" + llordd + "_" +"frac", 0.5, 0., 1.);
//   RooRealVar frac2("frac2","frac2", 0.9, 0., 1.);
  
  RooRealVar sigma("sig_" + llordd + "_" +"sigma","sig_" + llordd + "_" +"sigma",5.,0.,20.);
//   RooRealVar sigma2("sigma2","sigma2",20.,0.,100.);
//   RooRealVar sigma3("sigma3","sigma3",30.,0.,100.);
  
//   RooRealVar dm2("dm2","dm2",0.,-60.,60);
//   RooRealVar dm3("dm3","dm3",0.,-60.,60);
//   RooFormulaVar mean2("mean2","mean1+dm2",RooArgList(mean1,dm2));
//   RooFormulaVar mean3("mean3","mean1+dm3",RooArgList(mean1,dm3));

//   RooGaussian gauss1("gauss1","gauss1",mass,mean1,sigma1);
//   RooGaussian gauss2("gauss2","gauss2",mass,mean2,sigma2);
//   RooGaussian gauss3("gauss3","gauss3",mass,mean3,sigma3);

  RooRealVar alpha1("sig_" + llordd + "_" +"alpha1","sig_" + llordd + "_" +"alpha1", 1.,-1.,10.);
  RooRealVar alpha2("sig_" + llordd + "_" +"alpha2","sig_" + llordd + "_" +"alpha2",-1.,-10.,1.);

  RooRealVar n1("sig_" + llordd + "_" +"n1","sig_" + llordd + "_" +"n1",10.,0.,1000.);
  RooRealVar n2("sig_" + llordd + "_" +"n2","sig_" + llordd + "_" +"n2",10.,0.,1000.);

  RooCBShape cbshape1("cbshape1", "cbshape1", mass, mean, sigma, alpha1,n1);//, RooAbsReal& _n)
  RooCBShape cbshape2("cbshape2", "cbshape2", mass, mean, sigma, alpha2,n2);//, RooAbsReal& _n)


  //  RooAddPdf twogauss("twogauss","twogauss",RooArgList(gauss1,gauss2),RooArgList(frac1));
  //  RooAddPdf threegauss("threegauss","threegauss",RooArgList(gauss1,gauss2,gauss3),RooArgList(frac1,frac2),true);

  RooAddPdf twocb("twocb","twocb",RooArgList(cbshape1,cbshape2),RooArgList(frac));

  //  RooAddPdf gausspcbshape("gausspcbshape","gausspcbshape",RooArgList(gauss1,cbshape1),RooArgList(frac1));

  //  RooAbsPdf& sig_mass_model = gausspcbshape;
  //  RooAbsPdf& sig_mass_model = twogauss;
  //  RooAbsPdf& sig_mass_model = cbshape;
  
  //  RooAbsPdf* sig_mass_model = &twogauss;
  //  if (doBd2JpsiKS0) 
  //  RooAbsPdf* sig_mass_model = &threegauss;
  RooAbsPdf* sig_mass_model = &twocb;

  RooRealVar meanbw("meanbw","meanbw",500.,400.,600.);
  if (!doBd2JpsiKS0) {
    meanbw.setMin(1050.);   
    meanbw.setMax(1150.);
    meanbw.setVal(1115.);
  }
  RooRealVar sigmabw("sigmabw","sigmabw", 3., 0., 100.);

  //  RooBreitWigner bw("bw","bw",ksmass,meanbw,sigmabw);
  RooGaussian bw("bw","bw",ksmass,meanbw,sigmabw);

  //  bw.fitTo(*mcll[0]);
  //  plot(c,ksmass,bw,mcll[0], "bw-mass");
  // return 0;
  
  for (unsigned int i=0; i<directorylist.size(); i++) {
    TString plotprefix("fitmcmass-");;
    if (doBd2JpsiKS0) {plotprefix+="Bd2JpsiKS0";}else {plotprefix+="Lb2JpsiL0";}
    plotprefix+="-"+directorylist[i];
    c.SetLogy();    
    if (doLL) {
      simplemassfit( c,*sig_mass_model, mcll[i], plotprefix + "-mc-ll");
    } else {
      simplemassfit( c,*sig_mass_model, mcdd[i], plotprefix + "-mc-dd");
    }
    
  }

}

