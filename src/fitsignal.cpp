//g++ -o fitsignal fitsignal.cxx -O2 `root-config --cflags --libs` -lRooFit

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TEventList.h"
#include "TH1F.h"
#include "TH3F.h"

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

#include "RooStats/SPlot.h"

#include "rootlogon.h"

#include "functions-roofit.h"

#include "RooLbtoJpsiL0wAccPDF3.h"
#include "RooLbtoJpsiL0PDF3.h"
#include "createRooLegendre3.h"

#include "calcacceptanceclass.h"

RooRealVar* mass;
RooLbtoJpsiL0wAccPDF3* w3wacc = NULL;

RooRealVar nSig("nSig","nSig",2000.,0.,300000.);
//  RooRealVar nSig2("nSig2","nSig2",2000.,0.,5000.);
RooRealVar nBkg("nBkg","nBkg",10000.,0.,1000000.);

RooRealVar costheta("costheta",   "cos(#theta)",     -1., 1.);
RooRealVar costheta1("costheta1", "cos(#theta_{1})", -1., 1.);
RooRealVar costheta2("costheta2", "cos(#theta_{2})", -1., 1.);

//TH3F* acc_histo;


void simplemassfit(TCanvas& c, RooAbsPdf& model, RooDataSet* data, const TString& plotname)
{
  RooFitResult* frll = model.fitTo(*data,RooFit::Minos(true),RooFit::Save(true));
  std::cout << plotname << " mass fit (status=" << frll->status() << "; covQual=" << frll->covQual() << "; edm=" << frll->edm() << ")" << std::endl;
  if (frll->status()!=0 or frll->covQual()!=3) {
    std::cout << plotname << " mass fit - not successful (status=" << frll->status() << "; covQual=" << frll->covQual() << "; edm=" << frll->edm() << ")" << std::endl;
    mygetchar;
  }
  frll->Print("v");

  RooArgSet floatpars(frll->floatParsFinal());
  floatpars.writeToFile(plotname + ".txt");
  delete frll;
  mygetchar;

  plot(c,*mass,model,data, plotname + "-mass.eps");
  plotwpull(c,*mass,model,data, plotname + "-masswpull.eps");
}

void dosplot(TCanvas& c, RooAbsPdf& model, RooDataSet* data, const TString& plotname,bool dosplotplots, const RooArgList& splotvars, const RooArgList& splotvarstoplot)
{
  RooFitResult* frll = model.fitTo(*data,RooFit::Minos(true),RooFit::Save(true),RooFit::Extended(true));
  std::cout << plotname << " mass fit (status=" << frll->status() << "; covQual=" << frll->covQual() << "; edm=" << frll->edm() << ")" << std::endl;
  if (frll->status()!=0 or frll->covQual()!=3) {
    std::cout << plotname << " mass fit - not successful (status=" << frll->status() << "; covQual=" << frll->covQual() << "; edm=" << frll->edm() << ")" << std::endl;
    mygetchar;
  }
  frll->Print("v");
  RooArgSet floatpars(frll->floatParsFinal());
  floatpars.writeToFile(plotname + "-float.txt");
  RooArgSet constpars(frll->constPars());
  constpars.writeToFile(plotname + "-const.txt");
  delete frll;
  mygetchar;

  plot(c,*mass,model,data, plotname + "-mass.eps");
  plotwpull(c,*mass,model,data, plotname + "-masswpull.eps");

  if (!dosplotplots) return;

  for (int i=0; i<splotvars.getSize(); i++) {
    RooRealVar* myvar = (RooRealVar*)splotvars.at(i);
    myvar->setConstant();
  }

  RooStats::SPlot* sData = new RooStats::SPlot("sData-" + plotname,"An SPlot",
					       *data, &model, RooArgList(nSig,nBkg) );

  sData->Print("v");

  RooDataSet data_sig("data_sig","data_sig",data,*data->get(),0,"nSig_sw") ;
  RooDataSet data_bkg("data_bkg","data_bkg",data,*data->get(),0,"nBkg_sw") ;
  
  data_sig.get()->Print("v");

  plotdata(c,splotvarstoplot, &data_sig, plotname + "-splot-signal");
  mygetchar;
  plotdata(c,splotvarstoplot, &data_bkg, plotname + "-splot-bkg");


  for (int i=0; i<splotvars.getSize(); i++) {
    RooRealVar* myvar = (RooRealVar*)splotvars.at(i);
    myvar->setConstant(false);
  }
}

void compsigsplot(TCanvas& c, const RooArgList& varstoplot, RooDataSet* data1, RooDataSet* data2, const TString& plotname)
{
  
  RooDataSet data1_sig("data1_sig","data1_sig",data1,*data1->get(),0,"nSig_sw") ;
  RooDataSet data2_sig("data2_sig","data2_sig",data2,*data2->get(),0,"nSig_sw") ;

  plotdatavsdata(c,varstoplot, &data1_sig, &data2_sig, plotname);
}

void compsigsplotmc(TCanvas& c, const RooArgList& varstoplot, RooDataSet* data, RooDataSet* mc, const TString& plotname)
{
  
  RooDataSet data_sig("data_sig","data_sig",data,*data->get(),0,"nSig_sw") ;
  //  RooDataSet data2_sig("data2_sig","data2_sig",data2,*data2->get(),0,"nSig_sw") ;

  plotdatavsdata(c,varstoplot, &data_sig, mc, plotname);
}


void fit(TCanvas& c, RooAbsPdf& model, RooDataSet* data, RooDataSet* mcdata, const TString& plotname, bool plotsplot, const RooArgList& splotvars, const RooArgList& splotvarstoplot, bool doangfit, const RooArgList& angfitvars)
{

  RooFitResult* frll = model.fitTo(*data,RooFit::Minos(true),RooFit::Save(true),RooFit::Extended(true));
  if (frll->status()!=0 or frll->covQual()!=3) {
  std::cout << plotname << " mass fit - not successful (status=" << frll->status() << "; covQual=" << frll->covQual() << ")" << std::endl;
    mygetchar;
  }
  frll->Print("v");
  delete frll;

  plot(c,*mass,model,data, plotname + "-mass.eps");
  plotwpull(c,*mass,model,data, plotname + "-masswpull.eps");
  //  plotvar(c,costheta,&data,"signalregion", plotname + "-costheta.eps");
  //  plotvar(c,costheta1,&data,"signalregion", plotname + "-costheta1.eps");
  //  plotvar(c,costheta2,&data,"signalregion", plotname + "-costheta2.eps");

  if (!plotsplot) return;

  //  int i=0;
  for (int i=0; i<splotvars.getSize(); i++) {
    RooRealVar* myvar = (RooRealVar*)splotvars.at(i);
    myvar->setConstant();
  }

  // frac.setConstant();
  // c0.setConstant();
  // c1.setConstant();
  // mean1.setConstant();
  // //  mean2.setConstant();
  // sigma1.setConstant();
  // sigma2.setConstant();
  RooStats::SPlot* sData = new RooStats::SPlot("sData-" + plotname,"An SPlot",
					       *data, &model, RooArgList(nSig,nBkg) );

  sData->Print("v");
  //  mygetchar();

  RooDataSet data_sig("data_sig","data_sig",data,*data->get(),0,"nSig_sw") ;
  RooDataSet data_bkg("data_bkg","data_bkg",data,*data->get(),0,"nBkg_sw") ;
  
  data_sig.get()->Print("v");
  //  mygetchar();
  //  data_sig.write(plotname + "-data.dat");

  //  RooArgList vars(costheta,costheta1,costheta2,phi1,phi2);
  //  vars.add(RooArgList(Lambdab_BDIRA,Lambdab_IPCHI2,Lambdab_FDCHI2,Lambdab_PT,Lambda0_PT,Lambdab_TAU,Lambda0_TAU));
  //  vars.add(RooArgList(*ratio1,*ratio2,*ratio3));

  //  TH1F costhetahisto("costhetahisto","costhetahisto",25,-1.,1.);


  const int nbin=10;
  TH3F costhetahisto("costhetahisto","costhetahisto",nbin,-1.,1.,nbin,-1.,1.,nbin,-1.,1.);
  TH3F mccosthetahisto("mccosthetahisto","mccosthetahisto",nbin,-1.,1.,nbin,-1.,1.,nbin,-1.,1.);
  data_sig.fillHistogram(&costhetahisto,RooArgList(costheta,costheta1,costheta2));
  mcdata->fillHistogram(&mccosthetahisto,RooArgList(costheta,costheta1,costheta2));
  std::cout << "mccosthetahisto entries = " << mccosthetahisto.GetEntries() << std::endl;
  //  mygetchar;
  costhetahisto.Divide(&mccosthetahisto);

  TH1D* his0=costhetahisto.ProjectionX();
  his0->SetMinimum(0.);
  his0->GetXaxis()->SetTitle(costheta.GetTitle());
  his0->Draw();
  c.SaveAs(plotname + "-test-costheta.eps");

  TH1D* his1=costhetahisto.ProjectionY();
  his1->SetMinimum(0.);
  his1->GetXaxis()->SetTitle(costheta1.GetTitle());
  his1->Draw();
  c.SaveAs(plotname + "-test-costheta1.eps");

  TH1D* his2=costhetahisto.ProjectionZ();
  his2->SetMinimum(0.);
  his2->GetXaxis()->SetTitle(costheta2.GetTitle());
  his2->Draw();
  c.SaveAs(plotname + "-test-costheta2.eps");




  plotdata(c,splotvarstoplot, &data_sig, plotname + "-splot-signal");
  if (mcdata) plotdatavsdata(c,splotvarstoplot, &data_sig, mcdata, plotname + "-splot-datavsmc");
  mygetchar;
  //  plotdata(c,costheta1,&data_sig, plotname + "-splot-signal-costheta1.eps");
  //  plotdata(c,costheta2,&data_sig, plotname + "-splot-signal-costheta2.eps");
  //  plotdata(c,phi1,&data_sig, plotname + "-splot-signal-phi1.eps");
  //  plotdata(c,phi2,&data_sig, plotname + "-splot-signal-phi2.eps");

  plotdata(c,splotvarstoplot, &data_bkg, plotname + "-splot-bkg");
  //  plotdata(c,costheta1,&data_bkg, plotname + "-splot-bkg-costheta1.eps");
  //  plotdata(c,costheta2,&data_bkg, plotname + "-splot-bkg-costheta2.eps");
  //  plotdata(c,phi1,&data_bkg, plotname + "-splot-bkg-phi1.eps");
  //  plotdata(c,phi2,&data_bkg, plotname + "-splot-bkg-phi2.eps");

  for (int i=0; i<splotvars.getSize(); i++) {
    RooRealVar* myvar = (RooRealVar*)splotvars.at(i);
    myvar->setConstant(false);
  }

  // frac.setConstant(false);
  // c0.setConstant(false);
  // c1.setConstant(false);
  // mean1.setConstant(false);
  // //  mean2.setConstant(false);
  // sigma1.setConstant(false);
  // sigma2.setConstant(false);


  if (!doangfit) return;


  RooFitResult* angfit = w3wacc->fitTo(data_sig,RooFit::Save(true), RooFit::SumW2Error(true));//,RooFit::Minos(true));
  angfit->Print("v");
  mygetchar;

  plot(c, angfitvars, *w3wacc, &data_sig, plotname + "-fitdata-ang");


}


void fitbinbybin(TCanvas& c, RooRealVar& var, RooAbsPdf& model, RooDataSet* datadd_Lambdab0, const TString& title, const RooArgList& vars)
{
  //  frac.setConstant();
  //  sigma1.setConstant();
  //  sigma2.setConstant();
  //  mean1.setConstant();
  
  for (int i=0; i<vars.getSize(); i++) {
    RooRealVar* myvar = (RooRealVar*)vars.at(i);
    myvar->setConstant();
  }


  const double varmin=var.getMin();
  const double varmax=var.getMax();
  const TString varname(var.GetName());
  const int ndiv=25;

  // RooCategory masscat("masscat","masscat");
  // for (int i=0;i<ndiv;i++) {
  //   //    TString varcut(varname + ">" + str(varmin+(varmax-varmin)/ndiv*i) + "&&" + varname + "<" + str(varmin+(varmax-varmin)/ndiv*(i+1)));
  //   TString rangename=varname + "-reg" + istr(i);
  //   var.setRange(rangename,varmin+(varmax-varmin)/ndiv*i, varmin+(varmax-varmin)/ndiv*(i+1));
  //   masscat.setRange(rangename,rangename);
  // }
  // mygetchar;

  TH1F histo("histo","histo",ndiv,varmin,varmax);
  for (int i=0;i<ndiv;i++) {
    TString varcut(varname + ">" + str(varmin+(varmax-varmin)/ndiv*i) + "&&" + 
                   varname + "<" + str(varmin+(varmax-varmin)/ndiv*(i+1)));
    std::cout << "costhetacut = " << varcut << std::endl;
    //    filldataset(datadd_Lambdab0,"datadd_Lambdab0","data", RooArgList(mass,costheta,costheta1,costheta2,phi1,phi2), sigtree,DD_selection + isLambdab0 + costhetacut);
    RooDataSet* datadd_Lambdab0_red = (RooDataSet*)datadd_Lambdab0->reduce(varcut);
    fit(c,model,datadd_Lambdab0_red, NULL,   title + "-reg" + istr(i),false,vars,vars,false,vars);
    std::cout << "nSig = " << nSig.getVal() << " +- " << nSig.getError() << std::endl;
    histo.SetBinContent(i+1,nSig.getVal());
    histo.SetBinError(i+1,nSig.getError());
    delete datadd_Lambdab0_red;
  }
  histo.SetMinimum(0.);
  histo.GetXaxis()->SetTitle(var.GetTitle());
  histo.Draw("E");
  c.SaveAs(title + ".eps");

  
  for (int i=0; i<vars.getSize(); i++) {
    RooRealVar* myvar = (RooRealVar*)vars.at(i);
    myvar->setConstant(false);
  }

  //  frac.setConstant(false);
  //  sigma1.setConstant(false);
  //  sigma2.setConstant(false);
  //  mean1.setConstant(false);

}

void addacceptanceweight(RooDataSet* data, RooRealVar& w, bool normalizeweights=false)
{
  RooDataSet dataww6("dataww6","dataww6",RooArgList(w));
  
  calcacceptanceclass mycalcacceptanceclass_ll("cc0c1c2-ll.txt");
  calcacceptanceclass mycalcacceptanceclass_dd("cc0c1c2-dd.txt");
  

  //find sum of weight
  double sumofweights=0.;

  const int nEntries=data->numEntries();
  double weights[nEntries];

  for (int i=0;i<nEntries;i++) {  
    const RooArgSet* blu = data->get(i);
    if (blu->getRealValue("piminus_TRACK_Type")==3.) {
      weights[i]=1./mycalcacceptanceclass_ll.evaluate(blu->getRealValue("costheta"),
						      blu->getRealValue("costheta1"),
							blu->getRealValue("costheta2"));
    } else {
      weights[i]=1./mycalcacceptanceclass_dd.evaluate(blu->getRealValue("costheta"),
						      blu->getRealValue("costheta1"),
						      blu->getRealValue("costheta2"));
    }
    sumofweights+=weights[i];
    
  }
  
  if (!normalizeweights) {
    sumofweights=nEntries;
  }
  
  
  for (int i=0;i<nEntries;i++) {
    w.setVal(1.*nEntries/sumofweights * weights[i]);
    dataww6.add(RooArgSet(w));
  }
  data->merge(&dataww6);

  //   const bool debug=true;
  //   if (debug) {
  //     double newsumofweights=0.;
  //     for (int i=0;i<nEntries;i++) {  
  //       const RooArgSet* blu = data->get(i);
  //       newsumofweights+=blu->getRealValue(w.GetName());
  //     }
  //     std::cout << "New sum of weights = " << newsumofweights << std::endl;
  //     std::cout << "nEntries           = " << nEntries << std::endl;
  //     getchar();
  //   }
}


void  fit4dim(TCanvas& c, RooDataSet* datadd, RooAbsPdf& tot_mass_model, RooAbsPdf& tot_model, const RooArgList& totmodelpdfs, 
              const TString& plotprefix)
{      
  RooDataSet newdata1("newdata1","newdata1", datadd,*datadd->get());

  RooRealVar wacc("wacc","wacc",0.);
  addacceptanceweight(&newdata1, wacc,true);
      
  RooArgList varswacc(*newdata1.get());
  varswacc.add(RooArgList(wacc));

  RooDataSet newdata2("newdata2","newdata2", &newdata1, varswacc, 0, "wacc");

  tot_mass_model.fitTo(newdata2, RooFit::SumW2Error(true), RooFit::Extended(true));
  //  plotwcomponents(c,RooArgList(*mass), tot_mass_model, massmodelpdfs, datadd, plotname + "-massfit");

  tot_model.fitTo(newdata2, RooFit::SumW2Error(true), RooFit::Extended(true));//, RooFit::Minos(true));
  mygetchar;
  
  plotwcomponents(c,RooArgList(*mass,costheta,costheta1,costheta2), tot_model, totmodelpdfs, &newdata2, plotprefix + "-totfit");
  plotwcutwcomponents(c,RooArgList(costheta,costheta1,costheta2),   tot_model, totmodelpdfs, &newdata2, "signalregion", plotprefix + "-totfit");
  plotwcutwcomponents(c,RooArgList(costheta,costheta1,costheta2),   tot_model, totmodelpdfs, &newdata2, "lowersideband",plotprefix + "-totfit");
  plotwcutwcomponents(c,RooArgList(costheta,costheta1,costheta2),   tot_model, totmodelpdfs, &newdata2, "uppersideband",plotprefix + "-totfit");
}



int main(int argc, char** argv)
{
  //  const bool doBd2JpsiKS0=true;
  //const bool doBd2JpsiKS0=false;
  const bool fillacceptanceweight=false;
  const bool userwmc=false;
  const bool usebdt=false;
  const bool usesimplecuts=false;
  const bool antitrigger=false;
  
  if (argc!=2) {
    std::cerr << argv[0] << " Lb/B0" << std::endl;
    exit(2);
  }

  //  rootlogon();
  lhcbstyle();


  const TString mode(argv[1]);
  const bool doBd2JpsiKS0=mode.Contains("B0");
  
  std::vector<std::string> directorylist;
  if (doBd2JpsiKS0) {
    //    directorylist.push_back("JpsiKSTuple_B2XMuMu");
    directorylist.push_back("JpsiKSTuple_detached_betas");
    //    directorylist.push_back("JpsiKSTuple_prescaled_betas");
  } else {
    directorylist.push_back("JpsiL0Tuple_B2XMuMu");
    directorylist.push_back("JpsiL0Tuple_betas");
    directorylist.push_back("JpsiL0Tuple_B2XMuMu_refitted");
    directorylist.push_back("JpsiL0Tuple_betas_refitted");
  }

  
  //  TString dirname("JpsiL0Tuple_B2XMuMu");
  //  TString dirname("JpsiL0Tuple_betas");
  //  if (doBd2JpsiKS0) dirname = "JpsiKSTuple_prescaled_betas";
  //  if (doBd2JpsiKS0) dirname = "JpsiKSTuple_B2XMuMu";
  //  if (doBd2JpsiKS0) dirname = "JpsiKSTuple_detached_betas";
  

  //  TFile accfile("acc_histo.root");
  //  acc_histo = (TH3F*)accfile.Get("acc_histo");


  //Lambda_b0 to J/psi Lambda0
  TString Lambdabname("lambda_b0");
  TString Lambda0name("lambda0");
  TString Jpsiname("jpsi1s");
  
  //B0 to J/psi KS0
  if (doBd2JpsiKS0) {
    Lambdabname="b0";
    Lambda0name="ks0";
  }
  
  mass = new RooRealVar(Lambdabname + "_M","M_{#Lambda^{0}J/#psi(1S)}",5500.,5750.,"MeV/c^{2}");

  if (doBd2JpsiKS0) {
    mass->setMin(5230.);
    mass->setMax(5450.);
    mass->SetTitle("M_{J/#psi(1S)K_{S}^{0}}");
  }

  //  RooRealVar frac("frac","frac", 0.5, 0., 1.);
  RooRealVar frac1("frac1","frac1", 0.9, -0.1, 1.1);
  RooRealVar frac2("frac2","frac2", 0.5, -0.1, 1.1);
  
  //  const double pi=atan(1.)*4.;//3.14159;
  RooRealVar phi1("phi1", "#phi_{1}/#pi", -1., 1.);
  RooRealVar phi2("phi2", "#phi_{2}/#pi", -1., 1.);

  RooRealVar mean1("mean1","mean1",5620.,5500.,5700.);
  //  RooRealVar mean1("mean1","mean1",5250.,5000.,5700.);
  if (doBd2JpsiKS0) {
    mean1.setVal(5290.);
    mean1.setMin(5260.);
    mean1.setMax(5320.);
  }
  
  RooRealVar sigma1("sigma1","sigma1",5.,0.,100.);
  RooGaussian gauss1("gauss1","gauss1",*mass,mean1,sigma1);
  
  //RooRealVar mean2("mean2","mean2",5620.,5500.,5700.);
  RooRealVar sigma2("sigma2","sigma2",20.,0.,100.);
  RooRealVar sigma3("sigma3","sigma3",5.,0.,100.);
  RooGaussian gauss2("gauss2","gauss2",*mass,mean1,sigma2);
  RooGaussian gauss3("gauss3","gauss3",*mass,mean1,sigma3);
  
  RooRealVar c0("c0","c0",0.,-1.1,1.1);
  //  c0.setConstant();
  RooRealVar c1("c1","c1",0.,-1.1,1.1);
  
  
  //  RooRealVar Lambdab_ENDVERTEX_CHI2(Lambdabname + "_ENDVERTEX_CHI2", Lambdabname + "_ENDVERTEX_CHI2", 0., 100.);
  //  RooRealVar Lambdab_ENDVERTEX_NDOF(Lambdabname + "_ENDVERTEX_NDOF", Lambdabname + "_ENDVERTEX_NDOF", 0., 100.);
  RooRealVar Lambdab_ENDVERTEX_CHI2(Lambdabname + "_ENDVERTEX_CHI2", Lambdabname + "_ENDVERTEX_CHI2", 0.);
  RooRealVar Lambdab_ENDVERTEX_NDOF(Lambdabname + "_ENDVERTEX_NDOF", Lambdabname + "_ENDVERTEX_NDOF", 0.);
  RooFormulaVar Lambdab_ENDVERTEX_CHI2NDOF(Lambdabname + "_ENDVERTEX_CHI2NDOF", 
                                           Lambdabname + "_ENDVERTEX_CHI2NDOF",
                                           Lambdabname + "_ENDVERTEX_CHI2/" + Lambdabname + "_ENDVERTEX_NDOF", 
                                           RooArgList(Lambdab_ENDVERTEX_CHI2,Lambdab_ENDVERTEX_NDOF));
  
  
  //  RooRealVar Jpsi_ENDVERTEX_CHI2(Jpsiname + "_ENDVERTEX_CHI2", Jpsiname + "_ENDVERTEX_CHI2", 0., 100.);
  //  RooRealVar Jpsi_ENDVERTEX_NDOF(Jpsiname + "_ENDVERTEX_NDOF", Jpsiname + "_ENDVERTEX_NDOF", 0., 100.);
  RooRealVar Jpsi_ENDVERTEX_CHI2(Jpsiname + "_ENDVERTEX_CHI2", Jpsiname + "_ENDVERTEX_CHI2", 0.);
  RooRealVar Jpsi_ENDVERTEX_NDOF(Jpsiname + "_ENDVERTEX_NDOF", Jpsiname + "_ENDVERTEX_NDOF", 0.);
  RooFormulaVar Jpsi_ENDVERTEX_CHI2NDOF(Jpsiname + "_ENDVERTEX_CHI2NDOF", 
                                        Jpsiname + "_ENDVERTEX_CHI2NDOF",
                                        Jpsiname + "_ENDVERTEX_CHI2/" + Jpsiname + "_ENDVERTEX_NDOF", 
                                        RooArgList(Jpsi_ENDVERTEX_CHI2,Jpsi_ENDVERTEX_NDOF));
  
  //  RooRealVar Lambda0_ENDVERTEX_CHI2(Lambda0name + "_ENDVERTEX_CHI2", Lambda0name + "_ENDVERTEX_CHI2", 0., 100.);
  //  RooRealVar Lambda0_ENDVERTEX_NDOF(Lambda0name + "_ENDVERTEX_NDOF", Lambda0name + "_ENDVERTEX_NDOF", 0., 100.);
  
  
  //  RooRealVar Lambdab_BDIRA(Lambdabname + "_DIRA_OWNPV",Lambdabname + "_DIRA_OWNPV",0.9999,1.);
  //  RooRealVar Lambdab_IPCHI2(Lambdabname + "_IPCHI2_OWNPV",Lambdabname + "_IPCHI2_OWNPV",0.,30.);
  //  RooRealVar Lambdab_FDCHI2(Lambdabname + "_FDCHI2_OWNPV",Lambdabname + "_FDCHI2_OWNPV",100.,100000.);
  //  RooRealVar Lambdab_PT(Lambdabname + "_PT",Lambdabname + "_PT",0.,50000.);
  //  RooRealVar Lambda0_PT(Lambda0name + "_PT",Lambda0name + "_PT",0.,20000.);

  RooRealVar Lambdab_BDIRA(Lambdabname + "_DIRA_OWNPV",Lambdabname + "_DIRA_OWNPV",0.);
  RooRealVar Lambdab_IPCHI2(Lambdabname + "_IPCHI2_OWNPV",Lambdabname + "_IPCHI2_OWNPV",0.);
  RooRealVar Lambdab_FDCHI2(Lambdabname + "_FDCHI2_OWNPV",Lambdabname + "_FDCHI2_OWNPV",0.);
  RooRealVar Lambdab_PT(Lambdabname + "_PT",Lambdabname + "_PT",0.);
  RooRealVar Lambdab_P(Lambdabname + "_P",Lambdabname + "_P",0.);
  RooRealVar Lambdab_eta(Lambdabname + "_eta",Lambdabname + "_eta",0.);

  
  RooRealVar Lambdab_TAU(Lambdabname + "_TAU",Lambdabname + "_TAU",0.);
  //  RooRealVar Lambdab_TAU(Lambdabname + "_TAU",Lambdabname + "_TAU",0.,0.01);
  //  RooRealVar Lambda0_TAU(Lambda0name + "_TAU",Lambda0name + "_TAU",0.,1.);


  RooRealVar Lambda0_ENDVERTEX_CHI2(Lambda0name + "_ENDVERTEX_CHI2", Lambda0name + "_ENDVERTEX_CHI2", 0.);
  RooRealVar Lambda0_ENDVERTEX_NDOF(Lambda0name + "_ENDVERTEX_NDOF", Lambda0name + "_ENDVERTEX_NDOF", 0.);
  RooFormulaVar Lambda0_ENDVERTEX_CHI2NDOF(Lambda0name + "_ENDVERTEX_CHI2NDOF", 
                                           Lambda0name + "_ENDVERTEX_CHI2NDOF",
                                           Lambda0name + "_ENDVERTEX_CHI2/" + Lambda0name + "_ENDVERTEX_NDOF", 
                                           RooArgList(Lambda0_ENDVERTEX_CHI2,Lambda0_ENDVERTEX_NDOF));
  RooRealVar Lambda0_P(Lambda0name + "_P",Lambda0name + "_P",0.);
  RooRealVar Lambda0_PT(Lambda0name + "_PT",Lambda0name + "_PT",0.);
  RooRealVar Lambda0_TAU(Lambda0name + "_TAU",Lambda0name + "_TAU",0.);
  RooRealVar Lambda0_phi(Lambda0name + "_phi",Lambda0name + "_phi",0.);
  RooRealVar Lambda0_ENDVERTEX_Z(Lambda0name + "_ENDVERTEX_Z", Lambda0name + "_ENDVERTEX_Z", 0.);  

  RooRealVar Lambdab_ID(Lambdabname + "_ID",Lambdabname + "_ID",-10000.,10000.);
  RooRealVar piminus_TRACK_Type("piminus_TRACK_Type","piminus_TRACK_Type",0.,10.);
  RooRealVar Polarity("Polarity","Polarity",-10.,10.);
  
  RooRealVar bdt("bdt","BDT",0.);

  RooRealVar Jpsi_M(Jpsiname  + "_M",Jpsiname  + "_M",0. );

  //build TF3
  //  double dcoeff3[ordermax3i][ordermax3j][ordermax3k];
  RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k];
  //  RooRealVar* coeff5[ordermax5i][ordermax5j][ordermax5k][ordermax5l][ordermax5m];
  //  RooRealVar* coeff32[ordermax3i][ordermax3j][ordermax3k];
  //  createRooLegendre3(costheta, costheta1, costheta2, "accpdf", accpdf, coeff3);
  //  createRooLegendre3(costheta, costheta1, costheta2, "accpdf", accpdf, coeff3);


  RooRealVar P_b("P_b","P_{b}",0.0, -1.5, 1.5);//0.50
  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}",0.642,0.,1.);
  alpha_lambda.setConstant();

  //  P_b=0.;r_0=0.5;r_1=0.;

  RooRealVar alpha_b("alpha_b","#alpha_{b}",0.,-1.1,1.1);//0.457
  RooRealVar r_0("r_0","r_{0}",0.5,-0.3,2.3);//0.5
  RooRealVar r_1("r_1","r_{1}",0.,-1.3,1.3);//0.1

  createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b, alpha_lambda, r_0, r_1, "w3wacc", w3wacc, coeff3, true, "coeff3.txt");

  costheta.setBins(20);
  costheta1.setBins(20);
  costheta2.setBins(20);
  phi1.setBins(20);
  phi2.setBins(20);

  mass->setBins(50);
  mass->setRange("signalregion", 5600.,5640.);
  mass->setRange("lowersideband",5500.,5570.);
  mass->setRange("uppersideband",5670.,5750.);

  //  c1.setConstant();
  //  RooChebychev pol("pol","pol",mass,RooArgList(c0,c1));
  RooChebychev bkg_mass_model("bkg_mass_model","bkg_mass_model",*mass,RooArgList(c0));
  RooAddPdf twogauss("twogauss","twogauss",RooArgList(gauss1,gauss2),RooArgList(frac1));
  RooAddPdf threegauss("threegauss","threegauss",RooArgList(gauss1,gauss2,gauss3),RooArgList(frac1,frac2),true);

  RooAbsPdf& sig_mass_model = twogauss;
  //  RooAbsPdf& sig_mass_model = gauss1;
  RooAddPdf tot_mass_model("tot_mass_model","tot_mass_model",RooArgList(sig_mass_model,bkg_mass_model),RooArgList(nSig,nBkg));
  //  RooAddPdf model("model","model",RooArgList(gauss,pol),RooArgList(frac));

  RooRealVar bkg_costheta_c0("bkg_costheta_c0","bkg_costheta_c0",0.,-1.1,1.1);
  RooRealVar bkg_costheta_c1("bkg_costheta_c1","bkg_costheta_c1",0.,-1.1,1.1);
  RooRealVar bkg_costheta_c2("bkg_costheta_c2","bkg_costheta_c2",0.,-1.1,1.1);
  RooRealVar bkg_costheta1_c0("bkg_costheta1_c0","bkg_costheta1_c0",0.,-1.1,1.1);
  RooRealVar bkg_costheta1_c1("bkg_costheta1_c1","bkg_costheta1_c1",0.,-1.1,1.1);
  RooRealVar bkg_costheta1_c2("bkg_costheta1_c2","bkg_costheta1_c2",0.,-1.1,1.1);
  RooRealVar bkg_costheta2_c0("bkg_costheta2_c0","bkg_costheta2_c0",0.,-1.1,1.1);
  RooRealVar bkg_costheta2_c1("bkg_costheta2_c1","bkg_costheta2_c1",0.,-1.1,1.1);
  RooRealVar bkg_costheta2_c2("bkg_costheta2_c2","bkg_costheta2_c2",0.,-1.1,1.1);
  RooChebychev bkg_costheta_model("bkg_costheta_model",  "bkg_costheta_model", costheta, RooArgList(bkg_costheta_c0, bkg_costheta_c1, bkg_costheta_c2));
  RooChebychev bkg_costheta1_model("bkg_costheta1_model","bkg_costheta1_model",costheta1,RooArgList(bkg_costheta1_c0,bkg_costheta1_c1,bkg_costheta1_c2));
  RooChebychev bkg_costheta2_model("bkg_costheta2_model","bkg_costheta2_model",costheta2,RooArgList(bkg_costheta2_c0,bkg_costheta2_c1,bkg_costheta2_c2));

  RooProdPdf bkg_ang_model("bkg_ang_model","bkg_ang_model",RooArgList(bkg_costheta_model,bkg_costheta1_model,bkg_costheta2_model));

  RooLbtoJpsiL0PDF3 sig_ang_model("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1);


  RooProdPdf sig_model("sig_model","sig_model",RooArgList(sig_mass_model,sig_ang_model));
  RooProdPdf bkg_model("bkg_model","bkg_model",RooArgList(bkg_mass_model,bkg_ang_model));

  RooAddPdf tot_model("tot_model","tot_model",RooArgList(sig_model,bkg_model), RooArgList(nSig,nBkg));
  
  //  sigtree->SetEventList(&llsellist);

  TCanvas c("c","c",100,100);


  //  TFile sigfile(argv[1]);
  //  TFile sigfile("MergedLambdabTuple.reduced.root");

  TString mcfilename("DVTuples-Lb2JpsiL0-MC11a.root");
  //  const TString mcfilename("DVTuples-Bd2JpsiKS-MC11a.root");
  if (doBd2JpsiKS0) {
    mcfilename = "DVTuples-Bd2JpsiKS-MC11a.root";
  }
  //  if (doBd2JpsiKS0) mcfilename = "small.root";


  //  const TString sigfilename("DVTuples_data_stripping17.reduced.root");
  const TString sigfilename("DVTuples_data_stripping17.root");


  //  std::cout << "directory name = " << dirname << std::endl;
  //  plotprefix+="-"+dirname;

  TString sigbdtfilename(sigfilename);
  sigbdtfilename.ReplaceAll(".root",".bdt.root");
  TString siganglesfilename(sigfilename);
  siganglesfilename.ReplaceAll(".root", ".angles.root");
  TString sigksfilename(sigfilename);
  sigksfilename.ReplaceAll(".root", ".ksmass.root");

  TString mcbdtfilename(mcfilename);
  mcbdtfilename.ReplaceAll(".root",".bdt.root");
  TString mcanglesfilename(mcfilename);
  mcanglesfilename.ReplaceAll(".root", ".angles.root");
  TString mcksfilename(mcfilename);
  mcksfilename.ReplaceAll(".root", ".ksmass.root");

  TString mcrwfilename(mcfilename);
  mcrwfilename.ReplaceAll(".root",".reweighted.root");



  std::cout << "data file      = " << sigfilename << std::endl;
  std::cout << "mc file        = " << mcfilename << std::endl;
  std::cout << "mc angles file = " << mcanglesfilename << std::endl;
  std::cout << "directory list = ";
  for (unsigned int i=0;i<directorylist.size();i++) std::cout << directorylist[i] << " ";
  std::cout << std::endl;

  TTree* sigtree[directorylist.size()];
  TTree* mctree[directorylist.size()];

  memset(sigtree,0,directorylist.size()*sizeof(TTree*));
  memset(mctree,0,directorylist.size()*sizeof(TTree*));

  TFile sigfile(sigfilename);
  for (unsigned int i=0;i<directorylist.size();i++) {
    const TString& dirname = directorylist[i];
    TDirectory* sigdir = (TDirectory*)sigfile.Get(dirname);
    if (!sigdir) continue;
    sigtree[i] = (TTree*)sigdir->Get("DecayTree");
    if (usebdt) sigtree[i]->AddFriend(dirname, sigbdtfilename);
    sigtree[i]->AddFriend(dirname, siganglesfilename);
    sigtree[i]->AddFriend(dirname, TString(sigfilename).ReplaceAll(".root",".missvar.root"));
  }
  //  sigtree->AddFriend(dirname, ksfilename);


  //  RooDataSet* bkgdata = filldataset("data-","data", RooArgList(costheta,costheta1,costheta2), sigtree[0], "lambda_b0_M>5670." );
  //  bkgdata->write("bkgdata.txt");
  //  return 0;

  TFile mcfile(mcfilename);
  for (unsigned int i=0;i<directorylist.size();i++) {
    const TString& dirname = directorylist[i];
    TDirectory* mcdir = (TDirectory*)mcfile.Get(dirname);
    if (!mcdir) continue;
    mctree[i] = (TTree*)mcdir->Get("DecayTree");
    if (usebdt) mctree[i]->AddFriend(dirname, mcbdtfilename);
    mctree[i]->AddFriend(dirname, mcanglesfilename);
    if (userwmc) mctree[i]->AddFriend(dirname, mcrwfilename);
  }

  const TCut Lambda0_LL("piminus_TRACK_Type==3");
  const TCut Lambda0_DD("piminus_TRACK_Type==5");

  //  const TCut Lambdab_BDIRA_LL(Lambdabname + "_DIRA_OWNPV>0.999968");
  //  const TCut Lambdab_BDIRA_DD(Lambdabname + "_DIRA_OWNPV>0.999978");

  const TCut Lambdab_BDIRA_LL(Lambdabname + "_DIRA_OWNPV>0.9999");
  const TCut Lambdab_BDIRA_DD(Lambdabname + "_DIRA_OWNPV>0.9999");

  const TCut Lambdab_IPCHI2_LL(Lambdabname + "_IPCHI2_OWNPV<8.");
  const TCut Lambdab_IPCHI2_DD(Lambdabname + "_IPCHI2_OWNPV<7.");

  const TCut Lambdab_FDCHI2_LL(Lambdabname + "_FDCHI2_OWNPV>130.");
  const TCut Lambdab_FDCHI2_DD(Lambdabname + "_FDCHI2_OWNPV>180.");

  const TCut Lambdab_VCHI2DOF_LL(Lambdabname + "_ENDVERTEX_CHI2/" + Lambdabname + "_ENDVERTEX_NDOF<5.");
  const TCut Lambdab_VCHI2DOF_DD(Lambdabname + "_ENDVERTEX_CHI2/" + Lambdabname + "_ENDVERTEX_NDOF<4.");

  TCut Lambdab_M(Lambdabname + "_M>5500.&&" + Lambdabname + "_M<5750.");
  TCut Lambda0_M_LL("abs(" + Lambda0name + "_M-1115.683)<20.");
  TCut Lambda0_M_DD("abs(" + Lambda0name + "_M-1115.683)<15.");
  if (doBd2JpsiKS0) {
    Lambdab_M = TString(Lambdabname + "_M>5230.&&" + Lambdabname + "_M<5450.");
    Lambda0_M_LL = TString("abs(" + Lambda0name + "_M-497.614)<21.");
    Lambda0_M_DD = TString("abs(" + Lambda0name + "_M-497.614)<12.");
  }

  const TCut Lambda0_PT_LL(Lambda0name + "_PT>800.");
  const TCut Lambda0_PT_DD(Lambda0name + "_PT>1000.");

  const TCut Jpsi_VCHI2DOF_LL(Jpsiname + "_ENDVERTEX_CHI2/" + Jpsiname + "_ENDVERTEX_NDOF<9.");
  const TCut Jpsi_VCHI2DOF_DD(Jpsiname + "_ENDVERTEX_CHI2/" + Jpsiname + "_ENDVERTEX_NDOF<8.");

  const TCut MagnetUp("Polarity==1");
  const TCut MagnetDown("Polarity==-1");

  const TCut TriggerCut("L0DiMuonDecision==1&&Hlt1DiMuonHighMassDecision==1&&Hlt2DiMuonJPsiDecision==1.&&jpsi1sL0Global_TOS==1&&jpsi1sHlt1Global_TOS==1&&jpsi1sHlt2Global_TOS==1" );

  const TCut Lambda0_ENDVERTEX_Z_DD(Lambda0name + "_ENDVERTEX_Z>1500.");

  //  const TCut Lambdab_eta_cut(Lambdabname + "_eta>0.");

  //cut based selection
  //  TCut LL_selection = TriggerCut + Lambda0_LL + Lambdab_BDIRA_LL +  Lambdab_IPCHI2_LL + Lambdab_FDCHI2_LL + Lambdab_VCHI2DOF_LL  + Lambda0_M_LL  +  Lambda0_PT_LL + Jpsi_VCHI2DOF_LL;
  //  TCut DD_selection = TriggerCut + Lambda0_DD + Lambdab_BDIRA_DD +  Lambdab_IPCHI2_DD + Lambdab_FDCHI2_DD + Lambdab_VCHI2DOF_DD  + Lambda0_M_DD  +  Lambda0_PT_DD + Jpsi_VCHI2DOF_DD;

  TCut LL_selection =  Lambda0_LL + Lambdab_M + Lambdab_BDIRA_LL +  Lambdab_IPCHI2_LL + Lambdab_FDCHI2_LL + Lambdab_VCHI2DOF_LL  + Lambda0_M_LL  +  Lambda0_PT_LL + Jpsi_VCHI2DOF_LL;
  TCut DD_selection = Lambda0_DD + Lambdab_M + Lambdab_BDIRA_DD +  Lambdab_IPCHI2_DD + Lambdab_FDCHI2_DD + Lambdab_VCHI2DOF_DD  + Lambda0_M_DD  +  Lambda0_PT_DD + Jpsi_VCHI2DOF_DD;


  //  TCut LL_selection = Lambda0_LL + Lambdab_M + Lambdab_BDIRA_LL;// +  Lambdab_IPCHI2_LL + Lambdab_FDCHI2_LL + Lambdab_VCHI2DOF_LL  + Lambda0_M_LL  +  Lambda0_PT_LL + Jpsi_VCHI2DOF_LL;
  //  TCut DD_selection = Lambda0_DD + Lambdab_M + Lambdab_BDIRA_DD;// +  Lambdab_IPCHI2_DD + Lambdab_FDCHI2_DD + Lambdab_VCHI2DOF_DD  + Lambda0_M_DD  +  Lambda0_PT_DD + Jpsi_VCHI2DOF_DD;
  //  TCut DD_selection = Lambda0_DD + Lambda0_ENDVERTEX_Z_DD + Lambdab_M + Lambdab_BDIRA_DD;// +  Lambdab_IPCHI2_DD + Lambdab_FDCHI2_DD + Lambdab_VCHI2DOF_DD  + Lambda0_M_DD  +  Lambda0_PT_DD + Jpsi_VCHI2DOF_DD;

  //  TCut LL_selection = Lambda0_LL + Lambdab_M;// + Lambdab_eta_cut + Lambdab_BDIRA_LL;// +  Lambdab_IPCHI2_LL + Lambdab_FDCHI2_LL + Lambdab_VCHI2DOF_LL  + Lambda0_M_LL  +  Lambda0_PT_LL + Jpsi_VCHI2DOF_LL;
  //  TCut DD_selection = Lambda0_DD + Lambdab_M;// + Lambdab_eta_cut + Lambdab_BDIRA_DD;// +  Lambdab_IPCHI2_DD + Lambdab_FDCHI2_DD + Lambdab_VCHI2DOF_DD  + Lambda0_M_DD  +  Lambda0_PT_DD + Jpsi_VCHI2DOF_DD;



  //bdt selection
  if (usebdt) {
    //    LL_selection = Lambda0_LL + TriggerCut + Lambdab_M + TCut("bdt>0.0170");
    //    DD_selection = Lambda0_DD + TriggerCut + Lambdab_M + TCut("bdt>0.0432");
    const TCut BDTPRECUT=TCut(Lambdabname + "_DIRA_OWNPV>0.7&&" + Lambdabname + "_TAU>0.&&" + Lambdabname + "_IPCHI2_OWNPV<30.&&" + Lambda0name + "_TAU>0.");
    LL_selection = Lambda0_LL + Lambdab_M + BDTPRECUT + TCut("bdt>-0.02");
    DD_selection = Lambda0_DD + Lambdab_M + BDTPRECUT + TCut("bdt>-0.02");
  }

  if (usesimplecuts) {
    TCut simplecut_cutall(TriggerCut + Lambdab_M + TCut(Lambdabname + "_DIRA_OWNPV>0.7&&" + Lambdabname + "_TAU>0.&&" + Lambdabname + "_IPCHI2_OWNPV<30.&&" + Lambda0name + "_TAU>0."));
    if (antitrigger) simplecut_cutall = !TriggerCut + Lambdab_M + TCut(Lambdabname + "_DIRA_OWNPV>0.7&&" + Lambdabname + "_TAU>0.&&" + Lambdabname + "_IPCHI2_OWNPV<30.&&" + Lambda0name + "_TAU>0.");
    simplecut_cutall = Lambdab_M + TCut(Lambdabname + "_DIRA_OWNPV>0.7&&" + Lambdabname + "_TAU>0.&&" + Lambdabname + "_IPCHI2_OWNPV<30.&&" + Lambda0name + "_TAU>0.");
    LL_selection = Lambda0_LL + simplecut_cutall;
    DD_selection = Lambda0_DD + simplecut_cutall;
  }

  TString Lambdab_PDGID("5122");
  if (doBd2JpsiKS0) Lambdab_PDGID = "511";
  const TCut trueMC("abs(" + Lambdabname + "_TRUEID)==" + Lambdab_PDGID);

  TCut additionalmcsel = trueMC;// + TCut("L0DiMuonDecision==1");
  if (userwmc) additionalmcsel+=TCut("rw==1");

  const TCut selection = LL_selection or DD_selection;

  
  //  plotttree(c,RooArgList(costheta,costheta1,costheta2), mctree[0], trueMC + Lambda0_LL + TCut(Lambda0name +"_ENDVERTEX_Z<1500."), trueMC + Lambda0_LL + TCut(Lambda0name +"_ENDVERTEX_Z>1500."), "acceptance-ENDVERTEX_Z");
  //  plotttree(c,RooArgList(costheta,costheta1,costheta2), mctree[0], trueMC + Lambda0_DD + TCut(Lambda0name +"_ENDVERTEX_Z<1500."), trueMC + Lambda0_DD + TCut(Lambda0name +"_ENDVERTEX_Z>1500."), "acceptance-dd-ENDVERTEX_Z");
  //  plotttree(c,RooArgList(costheta,costheta1,costheta2), mctree[0], trueMC + Lambda0_LL + TCut(Lambdabname +"_PT<5000."), trueMC + Lambda0_LL + TCut(Lambdabname +"_PT>5000."), "acceptance-ll-b0_PT");
  //  plotttree(c,RooArgList(costheta,costheta1,costheta2), mctree[0], trueMC + Lambda0_DD + TCut(Lambdabname +"_PT<5000."), trueMC + Lambda0_DD + TCut(Lambdabname +"_PT>5000."), "acceptance-dd-b0_PT");
  //  return 0;
  

  //   for (unsigned int i=0; i<directorylist.size(); i++) {
  //     std::cout << "Directory: " << directorylist[i]  << std::endl;
  //     std::cout << "MC: All entries          : " << mctree[i]->GetEntries() << std::endl;
  //     std::cout << "MC: All entries (trigger
  
  //     std::cout << "data: All entries          : " << sigtree[i]->GetEntries() << std::endl;
  //     std::cout << "data: All entries (trigger): " << sigtree[i]->GetEntries(TriggerCut) << std::endl;
  //   }
  //   return 0;
  
  std::cout << "LL selection = "<< LL_selection << std::endl;
  std::cout << "DD selection = "<< DD_selection << std::endl;
  std::cout << "selection = "   << selection << std::endl;
  std::cout << "  MC sel  = "   << additionalmcsel << std::endl;
  mygetchar;


  //  sigtree->Draw("Lambdab_M_pmisid",TCut("bdt>0.") + Lambdab_M + TCut("Lambdab_M_pmisid>5100."));
  //  c.SaveAs("Lambdab_M_pmisid.eps");

  //  sigtree->Draw("Lambda0_M_pmisid",TCut("bdt>0.") + Lambdab_M);
  //  c.SaveAs("Lambda0_M_pmisid.eps");

  //  return 1;


  
  // TEventList llsellist("llsellist","llsellist");
  // TEventList ddsellist("ddsellist","ddsellist");
  
  // sigtree->Draw(">>llsellist",LL_selection);
  // sigtree->Draw(">>ddsellist",DD_selection);


  //  const TCut isLambdab0(Lambdabname + "_ID==5122.");
  //  const TCut isLambdab0bar(Lambdabname + "_ID==-5122.");

  const TCut isLambdab0(Lambdabname + "_ID>0.");
  const TCut isLambdab0bar(Lambdabname + "_ID<0.");

  //  TFile* outputfile = TFile::Open("/tmp/tmpfile.root","RECREATE");


  //  TTree* LLtree_isLambdab0 = sigtree->CopyTree( LL_selection + isLambdab0 );
  //  RooDataSet datall_Lambdab0("datall_Lambdab0","data",   LLtree_isLambdab0,  RooArgList(mass,costheta,costheta1,costheta2,phi1,phi2));//, RooFit::Import(*LLtree_isLambdab0));//, RooFit::Cut(LL_selection + isLambdab0));


  

  RooArgList vars(*mass,costheta,costheta1,costheta2,phi1,phi2);
  vars.add(RooArgList(Lambdab_ID,piminus_TRACK_Type,Polarity));
  //Lambdab stuff
  vars.add(RooArgList(Lambdab_BDIRA,Lambdab_IPCHI2,Lambdab_FDCHI2,Lambdab_PT,Lambdab_P,Lambdab_TAU,Lambdab_eta));
  vars.add(RooArgList(Lambdab_ENDVERTEX_CHI2,Lambdab_ENDVERTEX_NDOF));
  //Lambda0 stuff
  vars.add(RooArgList(Lambda0_P,Lambda0_PT,Lambda0_TAU,Lambda0_ENDVERTEX_CHI2,Lambda0_ENDVERTEX_NDOF,Lambda0_ENDVERTEX_Z,Lambda0_phi));
  //Jpsi stuff
  vars.add(RooArgList(Jpsi_ENDVERTEX_CHI2,Jpsi_ENDVERTEX_NDOF));
  vars.add(RooArgList(Jpsi_M));

  RooArgList varstoplot(costheta,costheta1,costheta2,phi1,phi2);
  varstoplot.add(RooArgList(Lambdab_BDIRA,Lambdab_IPCHI2,Lambdab_FDCHI2,Lambdab_PT,Lambdab_P,Lambdab_eta));
  varstoplot.add(RooArgList(Lambda0_PT,Lambdab_TAU,Lambda0_TAU,Lambda0_ENDVERTEX_Z,Lambda0_P,Lambda0_phi));
  varstoplot.add(RooArgList(Jpsi_M));
  
  if (usebdt) {
    vars.add(bdt);
    varstoplot.add(RooArgList(bdt));
  }
  //  varstoplot.add(RooArgList(Lambdab_ENDVERTEX_CHI2,Lambdab_ENDVERTEX_NDOF));
  //  varstoplot.add(RooArgList(Jpsi_ENDVERTEX_CHI2,Jpsi_ENDVERTEX_NDOF));
  //  varstoplot.add(RooArgList(Lambda0_ENDVERTEX_CHI2,Lambda0_ENDVERTEX_NDOF));

  

  RooArgList angfitvars(costheta,costheta1,costheta2);

  RooDataSet* data[directorylist.size()];
  RooDataSet* mcdata[directorylist.size()];
//   RooDataSet* dataddup[directorylist.size()];
//   RooDataSet* datallup[directorylist.size()];
//   RooDataSet* datadddown[directorylist.size()];
//   RooDataSet* datalldown[directorylist.size()];

  RooDataSet* datadd[directorylist.size()];
  RooDataSet* mcdd[directorylist.size()];
  RooDataSet* datall[directorylist.size()];
  RooDataSet* mcll[directorylist.size()];

  RooDataSet* datadd_Lambdab0[directorylist.size()];
  RooDataSet* mcdd_Lambdab0[directorylist.size()];
  RooDataSet* datall_Lambdab0[directorylist.size()];
  RooDataSet* mcll_Lambdab0[directorylist.size()];
  RooDataSet* datadd_Lambdab0bar[directorylist.size()];
  RooDataSet* mcdd_Lambdab0bar[directorylist.size()];
  RooDataSet* datall_Lambdab0bar[directorylist.size()];
  RooDataSet* mcll_Lambdab0bar[directorylist.size()];

//   RooDataSet* data_uppersideband_dd[directorylist.size()];
//   RooDataSet* data_lowersideband_dd[directorylist.size()];
  
//   RooDataSet* data_uppersideband_ll[directorylist.size()];
//   RooDataSet* data_lowersideband_ll[directorylist.size()];


  memset(mcdata,0,sizeof(RooDataSet*)*directorylist.size());
  memset(data,0,sizeof(RooDataSet*)*directorylist.size());


  RooRealVar* Lambdab_ENDVERTEX_CHI2NDOF_rrv = NULL;
  RooRealVar* Jpsi_ENDVERTEX_CHI2NDOF_rrv = NULL;
  RooRealVar* Lambda0_ENDVERTEX_CHI2NDOF_rrv = NULL;
  for (unsigned int i=0; i<directorylist.size(); i++) {
    if (!sigtree[i]) continue;
    RooDataSet* dataorig = filldataset("data-" + directorylist[i],"data-"+ directorylist[i], vars, sigtree[i], selection );
    
    RooRealVar wacc("wacc","wacc",0.);
    //    RooDataSet* data = NULL;
    if (fillacceptanceweight) {
      addacceptanceweight(dataorig, wacc);
      dataorig->Print("v");
      
      RooArgList varswacc(vars);
      varswacc.add(RooArgList(wacc));
      data[i] = new RooDataSet(TString("datawacc-") + directorylist[i], TString("datawacc-") + directorylist[i], dataorig, varswacc, 0, "wacc");
    } else {
      data[i] = dataorig;
    }
    
    data[i]->Print("v");
 
  
  //  mygetchar;

    Lambdab_ENDVERTEX_CHI2NDOF_rrv = (RooRealVar*)data[i]->addColumn(Lambdab_ENDVERTEX_CHI2NDOF);
    Jpsi_ENDVERTEX_CHI2NDOF_rrv = (RooRealVar*)data[i]->addColumn(Jpsi_ENDVERTEX_CHI2NDOF);
    Lambda0_ENDVERTEX_CHI2NDOF_rrv = (RooRealVar*)data[i]->addColumn(Lambda0_ENDVERTEX_CHI2NDOF);
    Lambdab_ENDVERTEX_CHI2NDOF_rrv->setMin(0.);
    Lambdab_ENDVERTEX_CHI2NDOF_rrv->setMax(10.);
    Lambda0_ENDVERTEX_CHI2NDOF_rrv->setMin(0.);
    Lambda0_ENDVERTEX_CHI2NDOF_rrv->setMax(10.);
    Jpsi_ENDVERTEX_CHI2NDOF_rrv->setMin(0.);
    Jpsi_ENDVERTEX_CHI2NDOF_rrv->setMax(10.);
    Lambdab_ENDVERTEX_CHI2NDOF_rrv->setBins(25);
    Lambda0_ENDVERTEX_CHI2NDOF_rrv->setBins(25);
    Jpsi_ENDVERTEX_CHI2NDOF_rrv->setBins(25);
    
    if (i==0) varstoplot.add(RooArgList(*Lambdab_ENDVERTEX_CHI2NDOF_rrv,*Jpsi_ENDVERTEX_CHI2NDOF_rrv,*Lambda0_ENDVERTEX_CHI2NDOF_rrv));
  }

  for (unsigned int i=0; i<directorylist.size(); i++) {
    mcdata[i]=filldataset("mcdata-" + directorylist[i],"mcdata-" + directorylist[i], vars, mctree[i], selection + additionalmcsel);
    mcdata[i]->addColumn(Lambdab_ENDVERTEX_CHI2NDOF);
    mcdata[i]->addColumn(Jpsi_ENDVERTEX_CHI2NDOF);
    mcdata[i]->addColumn(Lambda0_ENDVERTEX_CHI2NDOF);
    mcdata[i]->Print("v");
  }


  Jpsi_M.setMin(0.);
  Jpsi_M.setMax(5500.);
  Jpsi_M.setBins(100.);
  

  Lambdab_BDIRA.setMin(0.9999);
  Lambdab_BDIRA.setMax(1.0);
  Lambdab_IPCHI2.setMin(0.);
  Lambdab_IPCHI2.setMax(30.);
  Lambdab_FDCHI2.setMin(100.);
  Lambdab_FDCHI2.setMax(100000.);
  Lambdab_PT.setMin(0.);
  Lambdab_PT.setMax(50000.);
  if (doBd2JpsiKS0)   Lambdab_PT.setMax(30000.);
  Lambdab_P.setMin(0.);
  Lambdab_P.setMax(300000.);
  Lambdab_TAU.setMin(0.);
  Lambdab_TAU.setMax(0.01);
  Lambdab_eta.setMin(0.);
  Lambdab_eta.setMax(10.);

  Lambda0_PT.setMin(0.);
  Lambda0_PT.setMax(20000.);
  Lambda0_TAU.setMin(0.);
  Lambda0_TAU.setMax(1.);
  if (doBd2JpsiKS0) Lambda0_TAU.setMax(0.1);
  Lambda0_P.setMin(0.);
  Lambda0_P.setMax(200000.);
  Lambda0_ENDVERTEX_Z.setMin(-400.);
  Lambda0_ENDVERTEX_Z.setMax(3000.);
  Lambda0_phi.setMin(-1.);
  Lambda0_phi.setMax(1.);

  bdt.setMin(0.);
  bdt.setMax(1.);
  
  Lambdab_BDIRA.setBins(25);
  Lambdab_IPCHI2.setBins(25);
  Lambdab_FDCHI2.setBins(25);
  Lambdab_PT.setBins(25);
  Lambdab_P.setBins(25);
  Lambdab_TAU.setBins(25);
  Lambdab_eta.setBins(25);

  Lambda0_PT.setBins(25);
  Lambda0_TAU.setBins(25);
  Lambda0_ENDVERTEX_Z.setBins(25);
  Lambda0_P.setBins(25);
  Lambda0_phi.setBins(25);

  bdt.setBins(25);


  // RooDataSet* dataddup = NULL;
  // filldataset(dataddup,"dataddup","data", vars, sigtree,DD_selection + MagnetUp );
  // RooDataSet* datadddown = NULL;
  // filldataset(datadddown,"datadddown","data", vars, sigtree,DD_selection + MagnetDown );
  // RooDataSet* datallup = NULL;
  // filldataset(datallup,"datallup","data", vars, sigtree,LL_selection + MagnetUp );
  // RooDataSet* datalldown = NULL;
  // filldataset(datalldown,"datalldown","data", vars, sigtree,LL_selection + MagnetDown );

  // RooDataSet* datadd_Lambdab0 = NULL;
  // filldataset(datadd_Lambdab0,"datadd_Lambdab0","data", vars, sigtree,DD_selection + isLambdab0 );
  // RooDataSet* mcdd_Lambdab0 = NULL;
  // filldataset(mcdd_Lambdab0,"mcdd_Lambdab0","data", vars, mctree, DD_selection + isLambdab0 );

  // RooDataSet* datadd_Lambdab0bar=NULL;
  // filldataset(datadd_Lambdab0bar,"datadd_Lambdab0bar","data", vars,sigtree,DD_selection + isLambdab0bar);
  // RooDataSet* mcdd_Lambdab0bar=NULL;
  // filldataset(mcdd_Lambdab0bar,"mcdd_Lambdab0bar","data", vars,mctree,DD_selection + isLambdab0bar);

  // RooDataSet* datall_Lambdab0 = NULL;
  // filldataset(datall_Lambdab0,"datall_Lambdab0","data", vars, sigtree,LL_selection + isLambdab0); 
  // RooDataSet* mcll_Lambdab0 = NULL;
  // filldataset(mcll_Lambdab0,"mcll_Lambdab0","data", vars, mctree,LL_selection + isLambdab0); 

  // RooDataSet* datall_Lambdab0bar = NULL;
  // filldataset(datall_Lambdab0bar,"datall_Lambdab0bar","data",vars,sigtree,LL_selection + isLambdab0bar );
  // RooDataSet* mcll_Lambdab0bar = NULL;
  // filldataset(mcll_Lambdab0bar,"mcll_Lambdab0bar","data",vars,mctree,LL_selection + isLambdab0bar );


  //  RooArgList splotconstvars(frac1,frac2,c0,c1,mean1,sigma1,sigma2);
  RooArgList splotconstvars(c0,c1,mean1,sigma1);


  for (unsigned int i=0; i<directorylist.size(); i++) {
    if (!mcdata[i]) continue;
    mcdd[i] = (RooDataSet*)mcdata[i]->reduce(Lambda0_DD); 
    mcll[i] = (RooDataSet*)mcdata[i]->reduce(Lambda0_LL); 
    mcdd_Lambdab0[i] = (RooDataSet*)mcdata[i]->reduce(isLambdab0 + Lambda0_DD); 
    mcll_Lambdab0[i] = (RooDataSet*)mcdata[i]->reduce(isLambdab0 + Lambda0_LL); 
    mcdd_Lambdab0bar[i] = (RooDataSet*)mcdata[i]->reduce(isLambdab0bar + Lambda0_DD); 
    mcll_Lambdab0bar[i] = (RooDataSet*)mcdata[i]->reduce(isLambdab0bar + Lambda0_LL); 
  }
    
    
  for (unsigned int i=0; i<directorylist.size(); i++) {
//     dataddup[i] = (RooDataSet*)data[i]->reduce(MagnetUp + Lambda0_DD );
//     datallup[i] = (RooDataSet*)data[i]->reduce(MagnetUp + Lambda0_LL );
//     datadddown[i] = (RooDataSet*)data[i]->reduce(MagnetDown + Lambda0_DD );
//     datalldown[i] = (RooDataSet*)data[i]->reduce(MagnetDown + Lambda0_LL );
    if (!data[i]) continue;
    datadd[i] = (RooDataSet*)data[i]->reduce(Lambda0_DD); 
    datall[i] = (RooDataSet*)data[i]->reduce(Lambda0_LL); 
    
    datadd_Lambdab0[i] = (RooDataSet*)data[i]->reduce(isLambdab0 + Lambda0_DD); 
    datall_Lambdab0[i] = (RooDataSet*)data[i]->reduce(isLambdab0 + Lambda0_LL); 
    datadd_Lambdab0bar[i] = (RooDataSet*)data[i]->reduce(isLambdab0bar + Lambda0_DD); 
    datall_Lambdab0bar[i] = (RooDataSet*)data[i]->reduce(isLambdab0bar + Lambda0_LL); 
    
    
    
//     data_uppersideband_dd[i] = (RooDataSet*)data[i]->reduce(Lambda0_DD + TCut(Lambdabname + "_M>5670.")); 
//     data_lowersideband_dd[i] = (RooDataSet*)data[i]->reduce(Lambda0_DD + TCut(Lambdabname + "_M<5570.")); 
    
//     data_uppersideband_ll[i] = (RooDataSet*)data[i]->reduce(Lambda0_LL + TCut(Lambdabname + "_M>5670.")); 
//     data_lowersideband_ll[i] = (RooDataSet*)data[i]->reduce(Lambda0_LL + TCut(Lambdabname + "_M<5570.")); 
  }


  
//   for (int i=0;i<mcll[0]->sumEntries();i++) {
//     const RooArgSet* blu = mcll[0]->get(i);
//     std::cout << " (MC) b0_eta = " <<  blu->getRealValue("b0_eta") << std::endl;
//     std::cout << " (MC) piminus_TRACK_Type = " <<  blu->getRealValue("piminus_TRACK_Type") << std::endl;
//   }
//   for (int i=0;i<datall[0]->sumEntries();i++) {
//     const RooArgSet* blu = datall[0]->get(i);
//     std::cout << " (DATA) b0_eta = " <<  blu->getRealValue("b0_eta") << std::endl;
//     std::cout << " (DATA) piminus_TRACK_Type = " <<  blu->getRealValue("piminus_TRACK_Type") << std::endl;
//   }
//   exit(1);
	

  // //DD
  // bkg_ang_model.fitTo(*data_uppersideband_dd, RooFit::SumW2Error(true));
  // plot(c, RooArgList(costheta,costheta1,costheta2), bkg_ang_model, data_uppersideband_dd, plotprefix + "-blind-dd-bkgangfit-uppersideband");  
  // bkg_ang_model.fitTo(*data_lowersideband_dd, RooFit::SumW2Error(true));
  // plot(c, RooArgList(costheta,costheta1,costheta2), bkg_ang_model, data_lowersideband_dd, plotprefix + "-blind-dd-bkgangfit-lowersideband");  



  // //LL
  // bkg_ang_model.fitTo(*data_uppersideband_ll, RooFit::SumW2Error(true));
  // plot(c, RooArgList(costheta,costheta1,costheta2), bkg_ang_model, data_uppersideband_ll, plotprefix + "-blind-ll-bkgangfit-uppersideband");  
  // bkg_ang_model.fitTo(*data_lowersideband_ll, RooFit::SumW2Error(true));
  // plot(c, RooArgList(costheta,costheta1,costheta2), bkg_ang_model, data_lowersideband_ll, plotprefix + "-blind-ll-bkgangfit-lowersideband");  

  // tot_mass_model.fitTo(*datall, RooFit::SumW2Error(true));
  // plotwcomponents(c,RooArgList(*mass), tot_mass_model, RooArgList(sig_mass_model,bkg_mass_model), datall, plotprefix + "-blind-ll-massfit");
  // tot_model.fitTo(*datall, RooFit::SumW2Error(true));

  // plotwcomponents(c,RooArgList(*mass,costheta,costheta1,costheta2), tot_model, RooArgList(sig_model,bkg_model), datall, plotprefix + "-blind-ll-totfit");
  // plotwcutwcomponents(c,RooArgList(costheta,costheta1,costheta2), tot_model, RooArgList(sig_model,bkg_model), datall, "signalregion", plotprefix + "-blind-ll-totfit");
  // plotwcutwcomponents(c,RooArgList(costheta,costheta1,costheta2), tot_model, RooArgList(sig_model,bkg_model), datall, "lowersideband",plotprefix + "-blind-ll-totfit");
  // plotwcutwcomponents(c,RooArgList(costheta,costheta1,costheta2), tot_model, RooArgList(sig_model,bkg_model), datall, "uppersideband",plotprefix + "-blind-ll-totfit");



  //  return 0;
  for (unsigned int i=0; i<directorylist.size(); i++) {
    TString plotprefix;
    if (doBd2JpsiKS0) {plotprefix="Bd2JpsiKS0";}else {plotprefix="Lb2JpsiL0";}
    plotprefix+="-"+directorylist[i];
    
    if (!doBd2JpsiKS0) {
      simplemassfit( c,sig_mass_model, mcll[i], plotprefix + "-mc-ll");
      simplemassfit( c,sig_mass_model, mcdd[i], plotprefix + "-mc-dd");
      continue;
      RooArgSet sigparamsll(frac1,sigma1,mean1,sigma2);sigparamsll.readFromFile(plotprefix + "-mc-ll.txt");
      RooArgSet sigparamsdd(frac1,sigma1,mean1,sigma2);sigparamsdd.readFromFile(plotprefix + "-mc-dd.txt");

      sigma2.setConstant();
      frac1.setConstant();


      sigma2.setConstant();
      frac1.setConstant();
      sigma2.setVal(sigparamsll.getRealValue("sigma2"));
      frac1.setVal(sigparamsll.getRealValue("frac1"));
      //      std::cout << sigma2.getVal() << std::endl;
      //      mygetchar;
      dosplot( c,tot_mass_model, datall[i], plotprefix + "-ll", true, splotconstvars,varstoplot);
      compsigsplotmc(c, varstoplot, datall[i], mcll[i], plotprefix + "-ll-datavsmc");
      
      sigma2.setVal(sigparamsdd.getRealValue("sigma2"));
      frac1.setVal(sigparamsdd.getRealValue("frac1"));
      dosplot( c,tot_mass_model, datadd[i], plotprefix + "-dd", true, splotconstvars,varstoplot);
      compsigsplotmc(c, varstoplot, datadd[i], mcdd[i], plotprefix + "-dd-datavsmc");

      // if (usesimplecuts) {
      // 	sigma2.setVal(sigparamsdd.getRealValue("sigma2"));
      // 	frac1.setVal(sigparamsdd.getRealValue("frac1"));
      // 	fit(        c,tot_mass_model, datadd[i], mcdd[i], plotprefix + "-dd", true, splotconstvars,varstoplot, false, angfitvars);
      // 	sigma2.setVal(sigparamsll.getRealValue("sigma2"));
      // 	frac1.setVal(sigparamsll.getRealValue("frac1"));
      // 	fit(        c,tot_mass_model, datall[i], mcll[i], plotprefix + "-ll", true, splotconstvars,varstoplot, false, angfitvars);
	
      // } else {
	
      // 	P_b=0.;r_0=0.5;r_1=0.;
      // 	mygetchar;
      // 	sigma2.setVal(sigparamsdd.getRealValue("sigma2"));
      // 	frac1.setVal(sigparamsdd.getRealValue("frac1"));
      // 	fit(        c,tot_mass_model, datadd_Lambdab0[i], mcdd_Lambdab0[i], plotprefix + "-dd-Lambdab0", true, splotconstvars,varstoplot, false, angfitvars);
      // 	//    fitbinbybin(c,costheta,tot_mass_model, datadd_Lambdab0, plotprefix +"-dd-Lambdab0-costhetabinbybin", splotconstvars);
	
	
      // 	P_b=0.;r_0=0.5;r_1=0.;
      // 	sigma2.setVal(sigparamsdd.getRealValue("sigma2"));
      // 	frac1.setVal(sigparamsdd.getRealValue("frac1"));
      // 	fit(        c,tot_mass_model, datadd_Lambdab0bar[i], mcdd_Lambdab0bar[i], plotprefix + "-dd-Lambdab0bar", true,splotconstvars, varstoplot, false, angfitvars);
      // 	//    fitbinbybin(c,costheta,tot_mass_model, datadd_Lambdab0bar, plotprefix + "-dd-Lambdab0bar-costhetabinbybin", splotconstvars);
	
      // 	P_b=0.;r_0=0.5;r_1=0.;
      // 	sigma2.setVal(sigparamsll.getRealValue("sigma2"));
      // 	frac1.setVal(sigparamsll.getRealValue("frac1"));
      // 	fit(c,        tot_mass_model, datall_Lambdab0[i], mcll_Lambdab0[i],plotprefix + "-ll-Lambdab0", true, splotconstvars, varstoplot, false, angfitvars);
      // 	//    fitbinbybin(c,costheta,tot_mass_model, datall_Lambdab0, plotprefix + "-ll-Lambdab0-costhetabinbybin", splotconstvars);
	
      // 	P_b=0.;r_0=0.5;r_1=0.;
      // 	sigma2.setVal(sigparamsll.getRealValue("sigma2"));
      // 	frac1.setVal(sigparamsll.getRealValue("frac1"));
      // 	fit(        c, tot_mass_model, datall_Lambdab0bar[i], mcll_Lambdab0bar[i], plotprefix + "-ll-Lambdab0bar",true, splotconstvars,varstoplot,  false, angfitvars);
      // 	//    fitbinbybin(c, costheta, tot_mass_model, datall_Lambdab0bar, plotprefix + "-ll-Lambdab0bar-costhetabinbybin", splotconstvars);
	
	
      // 	compsigsplot(c,varstoplot,datadd_Lambdab0[i], datadd_Lambdab0bar[i], plotprefix + "-dd-splot-Lambdab0vsLambdab0bar");
      // 	compsigsplot(c,varstoplot,datall_Lambdab0[i], datall_Lambdab0bar[i], plotprefix + "-ll-splot-Lambdab0vsLambdab0bar");
      //    }
      //      P_b=0.;alpha_b=0.;r_0=0.5;r_1=0.;
      //      fit4dim(c, datadd_Lambdab0[i], tot_mass_model, tot_model, RooArgList(sig_model,bkg_model), plotprefix + "-blind-dd-Lambdab0");
      //      P_b=0.;alpha_b=0.;r_0=0.5;r_1=0.;
      //      fit4dim(c, datall_Lambdab0[i], tot_mass_model, tot_model, RooArgList(sig_model,bkg_model), plotprefix + "-blind-ll-Lambdab0");
  } else {
      //      simplemassfit( c,sig_mass_model, mcll[i], plotprefix + "-mc-ll");
      //      simplemassfit( c,sig_mass_model, mcdd[i], plotprefix + "-mc-dd");

      RooArgSet sigparamsll(frac1,sigma1,mean1,sigma2);sigparamsll.readFromFile(plotprefix + "-mc-ll.txt");
      RooArgSet sigparamsdd(frac1,sigma1,mean1,sigma2);sigparamsdd.readFromFile(plotprefix + "-mc-dd.txt");

      sigma2.setConstant();
      frac1.setConstant();
      sigma2.setVal(sigparamsll.getRealValue("sigma2"));
      frac1.setVal(sigparamsll.getRealValue("frac1"));
      //      std::cout << sigma2.getVal() << std::endl;
      //      mygetchar;
      dosplot( c,tot_mass_model, datall[i], plotprefix + "-ll", true, splotconstvars,varstoplot);
      compsigsplotmc(c, varstoplot, datall[i], mcll[i], plotprefix + "-ll-datavsmc");

      sigma2.setVal(sigparamsdd.getRealValue("sigma2"));
      frac1.setVal(sigparamsdd.getRealValue("frac1"));
      dosplot( c,tot_mass_model, datadd[i], plotprefix + "-dd", true, splotconstvars,varstoplot);
      compsigsplotmc(c, varstoplot, datadd[i], mcdd[i], plotprefix + "-dd-datavsmc");

      //      fit(        c,tot_mass_model, datadd[i], mcdd[i], plotprefix + "-dd", true, splotconstvars,varstoplot, false, angfitvars);
      //      fit(        c,tot_mass_model, datall[i], mcll[i], plotprefix + "-ll", true, splotconstvars,varstoplot, false, angfitvars);
    }

    
    
    //    return 0;
    
    
    //     fit(c,tot_mass_model, data[i], mcdata[i], plotprefix + "-alldata", true,   RooArgList(frac,c0,c1,mean1,sigma1,sigma2), 
    // 	varstoplot, false, angfitvars);
    
    //     fit(        c,tot_mass_model, datadd[i], mcdd[i], plotprefix + "-dd", true, splotconstvars,varstoplot, true, angfitvars);
    //     fit(        c,tot_mass_model, datall[i], mcll[i], plotprefix + "-ll", true, splotconstvars,varstoplot, true, angfitvars);
    
    
    //     if (!doBd2JpsiKS0) {
    //       fit(        c,tot_mass_model, dataddup[i], NULL, plotprefix + "-dd-MagUp",     true, 
    // 		  splotconstvars, varstoplot, false, angfitvars);
    //       fit(        c,tot_mass_model, datadddown[i], NULL, plotprefix + "-dd-MagDown", true, 
    // 		  splotconstvars, varstoplot, false, angfitvars);
    //       fit(        c,tot_mass_model, datallup[i], NULL, plotprefix + "-ll-MagUp",     true, splotconstvars,varstoplot, false, angfitvars);
    //       fit(        c,tot_mass_model, datalldown[i], NULL, plotprefix + "-ll-MagDown", true, splotconstvars,varstoplot, false, angfitvars);
    //     }
  }
  //  return 0;

  if (directorylist.size()>=2) {
    compsigsplot(c,varstoplot,datadd[0], datadd[1], "Bd2JpsiKS0-dd-splot-" + directorylist[0] + "_vs_" + directorylist[1]);
    compsigsplot(c,varstoplot,datall[0], datall[1], "Bd2JpsiKS0-ll-splot-" + directorylist[0] + "_vs_" + directorylist[1]);
  }

  //  TTree* DDtree_isLambdab0 = sigtree->CopyTree( DD_selection + isLambdab0 );

  // datall_Lambdab0->Print("v");
  // datall_Lambdab0bar->Print("v");

  // datadd_Lambdab0->Print("v");
  // datadd_Lambdab0bar->Print("v");
  // mygetchar(__FILE__,__LINE__);



  //  outputfile->Close();
  //  remove("/tmp/tmpfile.root");
  //  TTree* DDtree = sigtree->CopyTree( DD_selection );

  //  RooDataSet datall_Lambdab0("datall_Lambdab0","data",      LLtree_isLambdab0,    RooArgList(mass,costheta,costheta1,costheta2,phi1,phi2));




  //  RooDataHist* databinned = data.binnedClone();





  //  RooPlot* c = costheta.frame() ; 
  //  data_sig->plotOn(frame2, DataError(RooAbsData::SumW2) ) ; 



  // RooFitResult* frdd = tot_mass_model.fitTo(datadd,RooFit::Minos(true),RooFit::Save(true));
  // frdd->Print("v");
  // std::cout << "DD" << std::endl;
  // mygetchar;
  // plot(c,mass,tot_mass_model,&datadd,"massdd.eps");



}
