//g++ -o toy toy.cxx `root-config --lbs --cflags` -lRooFit

#include "TCanvas.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TVector.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooMinuit.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooStats/SPlot.h"
#include "RooMCStudy.h"
#include "RooCBShape.h"
#include "RooMinimizer.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooMultiVarGaussian.h"
#include "RooPullVar.h"
#include "RooBifurGauss.h"
#include "RooProduct.h"
#include "RooAddition.h"
#include "RooCmdConfig.h"
#include "RooConstraintSum.h"


#include "RooLbtoJpsiL0PDF3.h"
#include "RooLbtoJpsiL0PDF3v2.h"
#include "RooLbtoJpsiL0PDF3v3.h"
#include "calcacceptanceclass.h"
#include "rootlogon.h"

#include "functions.h"
#include "functions-roofit.h"
#include "ConfigFile.h"
#include "weighting.h"
#include "cuts.h"

const bool normalizetowmass=true;
const bool usebdt=true;
const bool userwmc=true;
const bool usetmvarectcut=false;
const bool doBd2JpsiKS0=false;

struct angfitstr {
  TString fitname;
  RooAbsPdf* model;
  RooAbsPdf* model4mc;
  RooCmdArg sumw2arg;
  RooCmdArg minosarg;
  RooFitResult* fr;
  RooMCStudy* mcstudy;
  RooArgList* vars;
};

void getBdbkg(int itoy, const RooArgList& vars, RooDataSet*& datadd, int ndd, RooDataSet*& datall, int nll)
{

  const TString sigfilename("/afs/cern.ch/user/j/jwicht/fitLambdab2JpsiL0/DVTuples-Bd2JpsiKS-MC11a.reduced.root");
  const TString dirname("Tuple_JpsiL0_betas");
  TFile sigfile(sigfilename);

  TDirectory* sigdir = (TDirectory*)sigfile.Get(dirname);
  TTree* sigtree = (TTree*)sigdir->Get("DecayTree");
  sigtree->AddFriend(dirname, TString(sigfilename).ReplaceAll(".root", ".angles.root"));
  sigtree->AddFriend(dirname, TString(sigfilename).ReplaceAll(".root", ".bdt.root"));
  
  TCut LL_selection(Lambda0_LL);
  TCut DD_selection(Lambda0_DD);

  //  bdtsel(doBd2JpsiKS0, LL_selection,DD_selection);
  //  LL_selection+=TCut("lambda_b0_M>5500.&&lambda_b0_M<5700.");
  //  DD_selection+=TCut("lambda_b0_M>5500.&&lambda_b0_M<5700.");

  RooDataSet* accdatall_full = filldataset("accdatall","accdatall",vars, sigtree, LL_selection); 
  RooDataSet* accdatadd_full = filldataset("accdatadd","accdatadd",vars, sigtree, DD_selection);

  std::cout << "Bd bkg: ll = " << accdatall_full->numEntries() << std::endl;
  std::cout << "Bd bkg: dd = " << accdatadd_full->numEntries() << std::endl;
  
  std::cout << "ndd = " << ndd << std::endl;
  std::cout << "nll = " << nll << std::endl;
  //  mygetchar;

  RooDataSet* dataddred = (RooDataSet*) accdatadd_full->reduce(RooFit::EventRange( itoy * ndd, (itoy+1) * ndd ) );
  RooDataSet* datallred = (RooDataSet*) accdatall_full->reduce(RooFit::EventRange( itoy * nll, (itoy+1) * nll ) );

  datadd->merge(dataddred);
  datall->merge(datallred);

  delete accdatall_full;
  delete accdatadd_full;
  delete dataddred;
  delete datallred;
}

void addpull(RooDataSet* data, RooRealVar* var, const RooRealVar& varval)
{
  TString name(var->GetName()), title(var->GetTitle()) ;
  name.Append("pull") ; title.Append(" Pull") ;
  RooRealVar pvar(name,title,0.) ;
  RooDataSet dataww6("dataww6","dataww6",RooArgList(pvar));
  const int nEntries=data->numEntries();

  for (int i=0;i<nEntries;i++) {     
    const RooArgSet* blu = data->get(i);
    RooRealVar* tmpvar = (RooRealVar*)blu->find(var->GetName());
    pvar=(tmpvar->getVal()  -varval.getVal())/tmpvar->getError();
    dataww6.add(RooArgSet(pvar));

    //    std::cout << i << " " << tmpvar->getVal() << " " << varval << " " << tmpvar->getError() << std::endl;

  }
  data->merge(&dataww6);
  //  mygetchar;
}


void plotPull(TCanvas& c, TString outputname, const RooRealVar& param, RooDataSet* data, Double_t lo, Double_t hi, Int_t nbins, Bool_t fitGauss, RooRealVar* mean=NULL, RooRealVar* sigma=NULL) 
{
  // Create a RooPlot of the pull distribution for the given
  // parameter.  The range lo-hi is plotted in nbins.  If fitGauss is
  // set, an unbinned ML fit of the distribution to a Gaussian p.d.f
  // is performed. The fit result is overlaid on the returned RooPlot
  // and a box with the fitted mean and sigma is added.

  TString name(param.GetName()), title(param.GetTitle()) ;
  name.Append("pull") ; title.Append(" Pull") ;
  RooRealVar pvar(name,title,lo,hi) ;
  pvar.setBins(nbins) ;

  RooPlot* frame = pvar.frame() ;
  data->plotOn(frame) ;

  if (fitGauss) {
    RooRealVar pullMean("pullMean","Mean of pull",0.,lo,hi) ;
    RooRealVar pullSigma("pullSigma","Width of pull",1.,0.,5.) ;
    RooGaussian* pullGauss = NULL;
    if (mean and sigma) {
      pullGauss = new RooGaussian("pullGauss","Gaussian of pull",
				  pvar,*mean,*sigma) ;
    } else {
      pullGauss = new RooGaussian("pullGauss","Gaussian of pull",
				  pvar,pullMean,pullSigma) ;
    }
    pullGauss->fitTo(*data,RooFit::Minos(false)) ;
    pullGauss->plotOn(frame) ;
    pullGauss->paramOn(frame,data) ;
    delete pullGauss;
  }

  frame->Draw();
  c.SaveAs(outputname + "-" + pvar.GetName() + ".eps");
  //  mygetchar;
  delete frame ;
}

bool hascategory(std::map<std::string,RooDataSet*>* datamapping, const std::string& category)
{
  return ( datamapping->find(category) != datamapping->end() );
}

void setminmax(RooArgList& list)
{
  for (int i=0;i<list.getSize();i++) {
    RooRealVar* var = (RooRealVar*)list.at(i);

    std::cout << "variable " << var->GetName() << " = " << var->getVal() << " +/- " << var->getError() << std::endl;
    var->setMin( var->getVal() - 3.* var->getError());
    var->setMax( var->getVal() + 3.* var->getError());
  }
}


RooDataSet* generatetoy(const TString& name, const TString& title, RooAbsPdf& model, const RooArgSet& vars, const int nevents)
{
  if (nevents==0) return new RooDataSet(name,title,vars);
  RooDataSet* newdata=model.generate(vars,nevents);
  newdata->SetName(name);
  newdata->SetTitle(title);
  return newdata;
}

bool isphysical(double alpha_b, double r_0, double r_1)
{
  if (r_0+r_1<0.) return false;
  if (r_0-r_1<0.) return false;
  if (alpha_b+1.-r_0-r_1<0.) return false;
  if (1.-alpha_b-r_0+r_1<0.) return false;
  return true;
}

int countphysicalsolutions(const RooDataSet& data)
{
  int retval=0;
  for (int i=0; i<data.numEntries(); i++) {
    const RooArgSet* blu = data.get(i);
    double alpha_b=blu->getRealValue("alpha_b");
    double r_0=blu->getRealValue("r_0");
    double r_1=blu->getRealValue("r_1");
    if (isphysical(alpha_b,r_0,r_1)) {
      std::cout << "alpha_b = " << alpha_b << " ; r_0 = " << r_0 << " ; r_1 = " << r_1 << std::endl;
      retval++;
    }
  }
  return retval;
}



RooDataSet* generatetoy(const TString& name, const TString& title, RooAbsPdf& model, const RooArgSet& vars, const int nevents, calcacceptanceclass* acceptance,TRandom3* rnd)
{
  progressbar pbar(std::string("Generate toy " + name));
  int ievents=0;
  RooDataSet* finaltoy = new RooDataSet(name,title,vars);
  if (nevents==0) return finaltoy;
  
  while (true) {
    RooDataSet* toy = model.generate(vars,1000);
    
    for (int i=0; i<toy->numEntries(); i++) {
      const RooArgSet* blu = toy->get(i);

      const double normvalue = acceptance->evaluate(blu->getRealValue("costheta"),
						    blu->getRealValue("costheta1"),
						    blu->getRealValue("costheta2"));
      const double rndvalue = rnd->Rndm();
      if (normvalue > rndvalue) {
	finaltoy->add(*blu);
	ievents++;
	pbar.print(100.*ievents/nevents);	
	if (ievents==nevents) {
	  delete toy;
	  pbar.finish();
	  return finaltoy;
	}
      }
    }
    delete toy;
  }
}
RooRealVar* findvar(RooDataSet* data, TString varname)
{
  const RooArgSet* blu = data->get(0);
  return (RooRealVar*)blu->find(varname);
}


int main(int argc, char** argv)
{
  const TString mode(argv[1]);

  if (argc!=3 or (mode!="fit" and mode!="merge")) {
    std::cerr << argv[0] << " fit iSample" << std::endl;
    std::cerr << argv[0] << " merge nSamples" << std::endl;
    exit(1);
  }
  const int iSample=atoi(argv[2]);

  lhcbstyle();
  TCanvas c("c","c",100,100);
  TRandom3* rnd;
  if (iSample==0) {
    rnd = new TRandom3(666);
    RooRandom::randomGenerator()->SetSeed(666);
  } else {
    rnd = new TRandom3(0);
    RooRandom::randomGenerator()->SetSeed(0);
  }
  
  RooCmdArg latexformat = RooFit::Format("E",RooFit::AutoPrecision(2),RooFit::VerbatimName(true));
 
  RooRealVar costheta("costheta",   "cos#theta",     -1., 1.);//RooArgList(theta));
  RooRealVar costheta1("costheta1", "cos#theta_{1}", -1., 1.);// RooArgList(theta1));
  RooRealVar costheta2("costheta2", "cos#theta_{2}", -1., 1.);// RooArgList(theta2));

  costheta.setBins(25);
  costheta1.setBins(25);
  costheta2.setBins(25);

  RooRealVar piminus_TRACK_Type("piminus_TRACK_Type","piminus_TRACK_Type",0.);

  RooRealVar P_b("P_b","P_{b}",0.40, -1.5, 1.5);//0.50

  const double al_val=0.642;
  const double al_err=0.013;
  // const double albar_val=0.71;
  // const double albar_err=0.08;

  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}",
			  al_val,
			  al_val-7.*al_err,  
			  al_val+7.*al_err);

  alpha_lambda.setMin(-1.);  
  alpha_lambda.setMax(+1.);

  // RooRealVar alpha_lambdabar("alpha_lambdabar","#alpha_{#bar{#lambda}}",
  // 			     albar_val, 
  // 			     albar_val-7.*albar_err,  
  // 			     albar_val+7.*albar_err);
  
  RooGaussian alpha_lambda_constrain("alpha_lambda_constrain",
				     "alpha_lambda_constrain",
				     alpha_lambda,
				     RooFit::RooConst(al_val),
				     RooFit::RooConst(al_err));
  // RooGaussian alpha_lambdabar_constrain("alpha_lambdabar_constrain",
  // 					"alpha_lambdabar_constrain",
  // 					alpha_lambdabar,
  // 					RooFit::RooConst(albar_val),
  // 					RooFit::RooConst(albar_err));
  
  RooRealVar alpha_b("alpha_b","#alpha_{b}",-0.45767,-1.1,1.1);//0.457
  RooRealVar r_0("r_0","r_{0}",0.251733,-0.3,2.3);//0.5
  RooRealVar r_1("r_1","r_{1}",0.116484,-1.3,1.3);//0.1

  // RooFormulaVar ap2_rfv("ap2_rrv","|a+|^{2}","0.5*(@0+@1)",RooArgList(r_0,r_1));
  // RooFormulaVar am2_rfv("am2_rrv","|a-|^{2}","0.5*(@0-@1)",RooArgList(r_0,r_1));
  // RooFormulaVar bp2_rfv("bp2_rrv","|b+|^{2}","0.5*(@0+1.-@1-@2)",RooArgList(alpha_b,r_0,r_1));
  // RooFormulaVar bm2_rfv("bm2_rrv","|b-|^{2}","0.5*(1.-@0-@1+@2)",RooArgList(alpha_b,r_0,r_1));
  RooFormulaVar* ap2_rfv = NULL;
  RooFormulaVar* am2_rfv = NULL;
  RooFormulaVar* bp2_rfv = NULL;
  RooFormulaVar* bm2_rfv = NULL;

  RooRealVar ap2("ap2","|a+|^{2}",0.);
  RooRealVar am2("am2","|a-|^{2}",0.);
  RooRealVar bp2("bp2","|b+|^{2}",0.);
  RooRealVar bm2_rrv("bm2_rrv","|b-|^{2}",0.);
  RooFormulaVar bm2("bm2","1.0-ap2-am2-bp2",RooArgList(ap2,am2,bp2));
  //  RooFormulaVar alpha_b_rfv("alpha_b_rfv","ap2-am2+bp2-bm2",RooArgList(ap2,am2,bp2,bm2));
  //  RooFormulaVar r_0_rfv("r_0_rfv","ap2+am2",RooArgList(ap2,am2));
  //  RooFormulaVar r_1_rfv("r_1_rfv","ap2-am2",RooArgList(ap2,am2));


  RooRealVar ap2corr("ap2corr","|a+|^{2}",0.);
  RooRealVar am2corr("am2corr","|a-|^{2}",0.);
  RooRealVar bp2corr("bp2corr","|b+|^{2}",0.);
  RooRealVar bm2corr("bm2corr","|b-|^{2}",0.);

  RooRealVar beta("beta","#beta",0.0,-5.0,1.3);//0.5
  RooRealVar gammarrv("gamma","#gamma", 0.,-5.0,5.0);//0.1
  
  RooRealVar mass("lambda_b0_M","M(#LambdaJ/#psi)",5550.,5700.,"MeV/c^{2}");  
  mass.setRange("signalregion", 5600.,5640.);

  RooRealVar sig_dd_mean("sig_dd_mean","sig_dd_mean",5620.,5500.,5700.);
  RooRealVar sig_ll_mean("sig_ll_mean","sig_ll_mean",5620.,5500.,5700.);

  RooRealVar sig_dd_frac("sig_dd_frac","sig_dd_frac", 0.5, -0.1, 1.1);
  RooRealVar sig_ll_frac("sig_ll_frac","sig_ll_frac", 0.5, -0.1, 1.1);

  RooRealVar sig_dd_sigma("sig_dd_sigma","sig_dd_sigma",7.5,0.,100.);
  RooRealVar sig_ll_sigma("sig_ll_sigma","sig_ll_sigma",7.5,0.,100.);


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
  RooCBShape sig_ll_cbshape1("sig_ll_cbshape1", "sig_ll_cbshape1", mass, sig_ll_mean, sig_ll_sigma, sig_ll_alpha1, sig_ll_n1);//, RooAbsReal& _n)
  RooCBShape sig_ll_cbshape2("sig_ll_cbshape2", "sig_ll_cbshape2", mass, sig_ll_mean, sig_ll_sigma, sig_ll_alpha2, sig_ll_n2);//, RooAbsReal& _n)

  RooAddPdf sig_dd_mass_model("sig_dd_mass_model","sig_dd_mass_model",RooArgList(sig_dd_cbshape1,sig_dd_cbshape2),RooArgList(sig_dd_frac));
  RooAddPdf sig_ll_mass_model("sig_ll_mass_model","sig_ll_mass_model",RooArgList(sig_ll_cbshape1,sig_ll_cbshape2),RooArgList(sig_ll_frac));


  RooRealVar bkg_dd_c0("bkg_dd_c0","bkg_dd_c0",0.,-1.1,1.1);
  RooRealVar bkg_ll_c0("bkg_ll_c0","bkg_ll_c0",0.,-1.1,1.1);

  RooRealVar bkg_dd_Lbbar_c0("bkg_dd_Lbbar_c0","bkg_dd_Lbbar_c0",0.,-1.1,1.1);
  RooRealVar bkg_ll_Lbbar_c0("bkg_ll_Lbbar_c0","bkg_ll_Lbbar_c0",0.,-1.1,1.1);


  RooChebychev bkg_dd_mass_model("bkg_dd_mass_model","bkg_dd_mass_model",mass,RooArgList(bkg_dd_c0));
  RooChebychev bkg_ll_mass_model("bkg_ll_mass_model","bkg_ll_mass_model",mass,RooArgList(bkg_ll_c0));
  RooChebychev bkg_dd_Lbbar_mass_model("bkg_dd_Lbbar_mass_model","bkg_dd_Lbbar_mass_model",mass,RooArgList(bkg_dd_Lbbar_c0));
  RooChebychev bkg_ll_Lbbar_mass_model("bkg_ll_Lbbar_mass_model","bkg_ll_Lbbar_mass_model",mass,RooArgList(bkg_ll_Lbbar_c0));

  RooRealVar nSig_dd("nSig_dd","nSig_dd",2000.,-2000.,300000.);
  RooRealVar nSig_ll("nSig_ll","nSig_ll",2000.,-2000.,300000.);
  RooRealVar nSig_dd_Lbbar("nSig_dd_Lbbar","nSig_dd_Lbbar",2000.,-2000.,300000.);
  RooRealVar nSig_ll_Lbbar("nSig_ll_Lbbar","nSig_ll_Lbbar",2000.,-2000.,300000.);

  RooRealVar nBkg_dd("nBkg_dd","nBkg_dd",10000.,-2000.,1000000.);
  RooRealVar nBkg_ll("nBkg_ll","nBkg_ll",10000.,-2000.,1000000.);
  RooRealVar nBkg_dd_Lbbar("nBkg_dd_Lbbar","nBkg_dd_Lbbar",10000.,-2000.,1000000.);
  RooRealVar nBkg_ll_Lbbar("nBkg_ll_Lbbar","nBkg_ll_Lbbar",10000.,-2000.,1000000.);

  RooAddPdf tot_dd_mass_model("tot_dd_mass_model","tot_dd_mass_model",RooArgList(sig_dd_mass_model,bkg_dd_mass_model),RooArgList(nSig_dd,nBkg_dd));
  RooAddPdf tot_ll_mass_model("tot_ll_mass_model","tot_ll_mass_model",RooArgList(sig_ll_mass_model,bkg_ll_mass_model),RooArgList(nSig_ll,nBkg_ll));

  RooAddPdf tot_dd_Lbbar_mass_model("tot_dd_Lbbar_mass_model","tot_dd_Lbbar_mass_model",RooArgList(sig_dd_mass_model,bkg_dd_Lbbar_mass_model),RooArgList(nSig_dd_Lbbar,nBkg_dd_Lbbar));
  RooAddPdf tot_ll_Lbbar_mass_model("tot_ll_Lbbar_mass_model","tot_ll_Lbbar_mass_model",RooArgList(sig_ll_mass_model,bkg_ll_Lbbar_mass_model),RooArgList(nSig_ll_Lbbar,nBkg_ll_Lbbar));

  ConfigFile cf("toyconfig.txt");
  
  const int nSigEvents_dd=atoi(cf.Value ("main","nSig_dd","6000").c_str());
  const int nSigEvents_dd_Lbbar=atoi(cf.Value ("main","nSig_dd_Lbbar","6000").c_str());
  const int nSigEvents_ll=atoi(cf.Value ("main","nSig_ll","6000").c_str());
  const int nSigEvents_ll_Lbbar=atoi(cf.Value ("main","nSig_ll_Lbbar","6000").c_str());
  const int nBkgEvents_dd=atoi(cf.Value ("main","nBkg_dd","6000").c_str());
  const int nBkgEvents_dd_Lbbar=atoi(cf.Value ("main","nBkg_dd_Lbbar","6000").c_str());
  const int nBkgEvents_ll=atoi(cf.Value ("main","nBkg_ll","6000").c_str());
  const int nBkgEvents_ll_Lbbar=atoi(cf.Value ("main","nBkg_ll_Lbbar","6000").c_str());
  
  nSig_dd=nSigEvents_dd;
  nSig_dd_Lbbar=nSigEvents_dd_Lbbar;
  nSig_ll=nSigEvents_ll;
  nSig_ll_Lbbar=nSigEvents_ll_Lbbar;

  nBkg_dd=nBkgEvents_dd;
  nBkg_dd_Lbbar=nBkgEvents_dd_Lbbar;
  nBkg_ll=nBkgEvents_ll;
  nBkg_ll_Lbbar=nBkgEvents_ll_Lbbar;

  const bool debug=atob(cf.Value("main","debug","false"));
  const bool doplot=atob(cf.Value("main","doplot","true"));
  //  const bool useformparam=atob(cf.Value("main","useformparam","false"));
  const bool inclacc=atob(cf.Value("main","inclacc","true"));
  const bool inclbkg=(nBkgEvents_dd+nBkgEvents_ll+nBkgEvents_dd_Lbbar+nBkgEvents_ll_Lbbar!=0);
  const bool floatalphalambda=atob(cf.Value("main","floatalphalambda","true"));
  const bool floatSigParams=atob(cf.Value("main","floatsigparams","true"));
  const int useparam=atoi(cf.Value("main","useparam","2").c_str());

  P_b=atof(cf.Value ("amplitudes","P_b","0.2").c_str());
  alpha_b=atof(cf.Value ("amplitudes","alpha_b","-0.45767").c_str());
  r_0=atof(cf.Value ("amplitudes","r_0","0.251733").c_str());
  r_1=atof(cf.Value ("amplitudes","r_1","0.116484").c_str());
  beta=atof(cf.Value ("amplitudes","beta","0.251733").c_str());
  gammarrv=atof(cf.Value ("amplitudes","gamma","0.116484").c_str());

  bkg_dd_c0=atof(cf.Value("bkgshape","bkg_dd_c0","0.").c_str());
  bkg_dd_Lbbar_c0=atof(cf.Value("bkgshape","bkg_dd_Lbbar_c0","0.").c_str());
  bkg_ll_c0=atof(cf.Value("bkgshape","bkg_ll_c0","0.").c_str());
  bkg_ll_Lbbar_c0=atof(cf.Value("bkgshape","bkg_ll_Lbbar_c0","0.").c_str());


  //  computew3 w3;
  //  ap2=w3.a_plus_sq(alpha_b.getVal(), r_0.getVal(), r_1.getVal());
  //  am2=w3.a_minus_sq(alpha_b.getVal(), r_0.getVal(), r_1.getVal());
  //  bp2=w3.b_plus_sq(alpha_b.getVal(), r_0.getVal(), r_1.getVal());

  const RooArgList sig_dd_vars(sig_dd_frac,sig_dd_mean,sig_dd_sigma,sig_dd_alpha1,sig_dd_alpha2,sig_dd_n1,sig_dd_n2);
  const RooArgList mass_dd_pdfs(sig_dd_mass_model,bkg_dd_mass_model);
  
  const RooArgList sig_ll_vars(sig_ll_frac,sig_ll_mean,sig_ll_sigma,sig_ll_alpha1,sig_ll_alpha2,sig_ll_n1,sig_ll_n2);
  const RooArgList mass_ll_pdfs(sig_ll_mass_model,bkg_ll_mass_model);

  const TString directorylist = "Tuple_JpsiL0_betas";
  readvarsfromtfile(sig_dd_vars,TString("fitmcmass-Lb2JpsiL0-") + directorylist + "-mc-dd-fr.root");
  readvarsfromtfile(sig_ll_vars,TString("fitmcmass-Lb2JpsiL0-") + directorylist + "-mc-ll-fr.root");

  
  RooArgList testdd(sig_dd_frac,sig_dd_alpha1,sig_dd_alpha2,sig_dd_n1,sig_dd_n2);
  RooArgList testll(sig_ll_frac,sig_ll_alpha1,sig_ll_alpha2,sig_ll_n1,sig_ll_n2);
  if (floatSigParams) {
    setminmax(testdd);
    setminmax(testll);
  }

  RooFitResult* frdd=NULL;
  RooFitResult* frll=NULL;
  if (doBd2JpsiKS0) {
    frdd=getRooFitResult(TString("fitmcmass-Bd2JpsiKS0-") + directorylist + "-mc-dd-fr.root");
    frll=getRooFitResult(TString("fitmcmass-Bd2JpsiKS0-") + directorylist + "-mc-ll-fr.root");
  } else {
    frdd=getRooFitResult(TString("fitmcmass-Lb2JpsiL0-") + directorylist + "-mc-dd-fr.root");
    frll=getRooFitResult(TString("fitmcmass-Lb2JpsiL0-") + directorylist + "-mc-ll-fr.root");
  }
  if (!frdd or !frll) {
    std::cerr << "Can't read signal shape for RooMultiVarGaussian" << std::endl;
    return 1;
  }
  RooMultiVarGaussian sig_dd_constrain("sig_dd_constrain","sig_dd_constrain",testdd,  *frdd);
  RooMultiVarGaussian sig_ll_constrain("sig_ll_constrain","sig_ll_constrain",testll,  *frll);
  delete frdd;
  delete frll;



  sig_dd_alpha1.setConstant(!floatSigParams);
  sig_dd_alpha2.setConstant(!floatSigParams);
  sig_dd_n1.setConstant(!floatSigParams);
  sig_dd_n2.setConstant(!floatSigParams);
  sig_dd_frac.setConstant(!floatSigParams);
  sig_ll_alpha1.setConstant(!floatSigParams);
  sig_ll_alpha2.setConstant(!floatSigParams);
  sig_ll_n1.setConstant(!floatSigParams);
  sig_ll_n2.setConstant(!floatSigParams);
  sig_ll_frac.setConstant(!floatSigParams);

  RooCmdArg massfitconstrains = RooCmdArg::none();
  if (floatSigParams) massfitconstrains = RooFit::ExternalConstraints(RooArgSet(sig_dd_constrain,sig_ll_constrain));

  RooCmdArg angfitconstrains = RooCmdArg::none();
  if (floatalphalambda) {
    angfitconstrains = RooFit::ExternalConstraints(RooArgSet(alpha_lambda_constrain));//,alpha_lambdabar_constrain));
  } else {
    //    alpha_lambda.setConstant();
    //    alpha_lambdabar.setConstant();
  }
  

  RooAbsPdf* sig_ang_model       = NULL;
  RooAbsPdf* sig_ang_model_Lbbar = NULL;

  if (useparam==0) {
    sig_ang_model = new RooLbtoJpsiL0PDF3v3("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, ap2, alpha_lambda, am2, bp2);
    sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF3v3("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,P_b, ap2, alpha_lambda, am2, bp2);

    ap2.setMin(-0.3);
    ap2.setMax(1.3);
    am2.setMin(-0.3);
    am2.setMax(1.3);
    bp2.setMin(-0.3);
    bp2.setMax(1.3);

    ap2_rfv = new RooFormulaVar("ap2_rrv","|a+|^{2}","@0",RooArgList(ap2));
    am2_rfv = new RooFormulaVar("am2_rrv","|a-|^{2}","@0",RooArgList(am2));
    bp2_rfv = new RooFormulaVar("bp2_rrv","|b+|^{2}","@0",RooArgList(bp2));
    bm2_rfv = new RooFormulaVar("bm2_rrv","|b-|^{2}","1.-@0-@1-@2",RooArgList(ap2,am2,bp2));


  } 
  if (useparam==1)  {
    
    //    sig_ang_model = new RooLbtoJpsiL0PDF5("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,phi1,phi2,P_b, alpha_b, alpha_lambda, r_0, r_1, alpha_plus, alpha_minus, chi);
    //    sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF5("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,phi1,phi2,P_b, alpha_b, alpha_lambdabar, r_0, r_1, alpha_plus, alpha_minus, chi);

    sig_ang_model = new RooLbtoJpsiL0PDF3("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1);
    sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF3("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1);
    
    ap2_rfv = new RooFormulaVar("ap2_rrv","|a+|^{2}","0.5*(@0+@1)",RooArgList(r_0,r_1));
    am2_rfv = new RooFormulaVar("am2_rrv","|a-|^{2}","0.5*(@0-@1)",RooArgList(r_0,r_1));
    bp2_rfv = new RooFormulaVar("bp2_rrv","|b+|^{2}","0.5*(1.+@0-@1-@2)",RooArgList(alpha_b,r_0,r_1));
    bm2_rfv = new RooFormulaVar("bm2_rrv","|b-|^{2}","0.5*(1.-@0-@1+@2)",RooArgList(alpha_b,r_0,r_1));
    
  } 
  if (useparam==2) {
    sig_ang_model = new RooLbtoJpsiL0PDF3v2("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, beta, gammarrv);
    sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF3v2("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, beta, gammarrv);

    ap2_rfv = new RooFormulaVar("ap2_rrv","|a+|^{2}","1./12.*(+3.*@0+3.*@1-4.*@2+2.)",RooArgList(alpha_b,beta,gammarrv));
    am2_rfv = new RooFormulaVar("am2_rrv","|a-|^{2}","1./12.*(-3.*@0-3.*@1-4.*@2+2.)",RooArgList(alpha_b,beta,gammarrv));
    bp2_rfv = new RooFormulaVar("bp2_rrv","|b+|^{2}","1./12.*(+3.*@0-3.*@1+4.*@2+4.)",RooArgList(alpha_b,beta,gammarrv));
    bm2_rfv = new RooFormulaVar("bm2_rrv","|b-|^{2}","1./12.*(-3.*@0+3.*@1+4.*@2+4.)",RooArgList(alpha_b,beta,gammarrv));

  }



  if (iSample==0 or mode=="merge") {
    std::cout << "Config" << std::endl;
    std::cout << "  nSigEvents_dd       = " << nSigEvents_dd << std::endl;
    std::cout << "  nSigEvents_dd_Lbbar = " << nSigEvents_dd_Lbbar << std::endl;
    std::cout << "  nSigEvents_ll       = " << nSigEvents_ll << std::endl;
    std::cout << "  nSigEvents_ll_Lbbar = " << nSigEvents_ll_Lbbar << std::endl;
    std::cout << "  nBkgEvents_dd       = " << nBkgEvents_dd << std::endl;
    std::cout << "  nBkgEvents_dd_Lbbar = " << nBkgEvents_dd_Lbbar << std::endl;
    std::cout << "  nBkgEvents_ll       = " << nBkgEvents_ll << std::endl;
    std::cout << "  nBkgEvents_ll_Lbbar = " << nBkgEvents_ll_Lbbar << std::endl;

    std::cout << "  debug        = " << debug << std::endl;  
    std::cout << "  doplot       = " << doplot << std::endl;
    //    std::cout << "  useformparam = " << useformparam << std::endl;
    std::cout << "  inclacc      = " << inclacc << std::endl;
    std::cout << "  inclbkg      = " << inclbkg << std::endl;
    std::cout << "  P_b          = " << P_b.getVal() << std::endl;
    if (useparam==0) {
      std::cout << "  ap2          = " << ap2.getVal() << std::endl;
      std::cout << "  am2          = " << am2.getVal() << std::endl;
      std::cout << "  bp2          = " << bp2.getVal() << std::endl;
    }
    //      std::cout << "  bm2          = " << bm2.getVal() << std::endl;
    //      std::cout << "  sum          = " << ap2.getVal()+am2.getVal()+bp2.getVal()+bm2.getVal() << std::endl;
      //      std::cout << "  alpha_b      = " << alpha_b_rfv.getVal() << std::endl;
      //      std::cout << "  r_0          = " << r_0_rfv.getVal() << std::endl;
      //      std::cout << "  r_1          = " << r_1_rfv.getVal() << std::endl;
      //    } else {
    if (useparam==1) {
      std::cout << "  alpha_b      = " << alpha_b.getVal() << std::endl;
      std::cout << "  r_0          = " << r_0.getVal() << std::endl;
      std::cout << "  r_1          = " << r_1.getVal() << std::endl;
    } 
    if (useparam==2) {
      std::cout << "  alpha_b      = " << alpha_b.getVal() << std::endl;
      std::cout << "  beta         = " << beta.getVal() << std::endl;
      std::cout << "  gamma        = " << gammarrv.getVal() << std::endl;
    }
    if (useparam==1 or useparam==2) {
      std::cout << "    ap2        = " << ap2_rfv->getVal() << std::endl;
      std::cout << "    am2        = " << am2_rfv->getVal() << std::endl;
      std::cout << "    bp2        = " << bp2_rfv->getVal() << std::endl;
      std::cout << "    bm2        = " << bm2_rfv->getVal() << std::endl;
    }
//    }
    std::cout << "  sig_dd_mean         = " << sig_dd_mean.getVal() << std::endl;
    std::cout << "  sig_dd_sigma        = " << sig_dd_sigma.getVal() << std::endl;
    std::cout << "  sig_dd_alpha1       = " << sig_dd_alpha1.getVal() << " +/- " << sig_dd_alpha1.getError() << std::endl;
    std::cout << "  sig_dd_alpha2       = " << sig_dd_alpha2.getVal() << " +/- " << sig_dd_alpha2.getError() << std::endl;
    std::cout << "  sig_dd_n1           = " << sig_dd_n1.getVal() << " +/- " << sig_dd_n1.getError() << std::endl;
    std::cout << "  sig_dd_n2           = " << sig_dd_n2.getVal() << " +/- " << sig_dd_n2.getError() << std::endl;
    std::cout << "  sig_dd_frac         = " << sig_dd_frac.getVal() << " +/- " << sig_dd_frac.getError() << std::endl;
    std::cout << "  sig_ll_mean         = " << sig_ll_mean.getVal() << std::endl;
    std::cout << "  sig_ll_sigma        = " << sig_ll_sigma.getVal() << std::endl;
    std::cout << "  sig_ll_alpha1       = " << sig_ll_alpha1.getVal() << " +/- " << sig_ll_alpha1.getError() << std::endl;
    std::cout << "  sig_ll_alpha2       = " << sig_ll_alpha2.getVal() << " +/- " << sig_ll_alpha2.getError() << std::endl;
    std::cout << "  sig_ll_n1           = " << sig_ll_n1.getVal() << " +/- " << sig_ll_n1.getError() << std::endl;
    std::cout << "  sig_ll_n2           = " << sig_ll_n2.getVal() << " +/- " << sig_ll_n2.getError() << std::endl;
    std::cout << "  sig_ll_frac         = " << sig_ll_frac.getVal() << " +/- " << sig_ll_frac.getError() << std::endl;

    std::cout << "  bkg_dd_c0           = " << bkg_dd_c0.getVal() << std::endl;
    std::cout << "  bkg_dd_Lbbar_c0     = " << bkg_dd_Lbbar_c0.getVal() << std::endl;
    std::cout << "  bkg_ll_c0           = " << bkg_ll_c0.getVal() << std::endl;
    std::cout << "  bkg_ll_Lbbar_c0     = " << bkg_ll_Lbbar_c0.getVal() << std::endl;
  }

  RooRealVar ap2_initval("ap2_initval","ap2_initval",ap2_rfv->getVal());
  RooRealVar am2_initval("am2_initval","am2_initval",am2_rfv->getVal());
  RooRealVar bp2_initval("bp2_initval","bp2_initval",bp2_rfv->getVal());
  RooRealVar bm2_initval("bm2_initval","bm2_initval",bm2_rfv->getVal());
  //  const double ap2_initval=ap2_rfv->getVal();
  //  const double am2_initval=am2_rfv->getVal();
  //  const double bp2_initval=bp2_rfv->getVal();
  //  const double bm2_initval=bm2_rfv->getVal();


  // const double Pb_initval=P_b.getVal();
  // const double ab_initval=alpha_b.getVal();
  // const double r0_initval=r_0.getVal();
  // const double r1_initval=r_1.getVal();

  RooRealVar Pb_initval("Pb_initval","Pb_initval",P_b.getVal());
  RooRealVar ab_initval("ab_initval","ab_initval",alpha_b.getVal());
  RooRealVar r0_initval("r0_initval","r0_initval",r_0.getVal());
  RooRealVar r1_initval("r1_initval","r1_initval",r_1.getVal());

  if (nSigEvents_dd==0.) nSig_dd.setConstant();
  if (nSigEvents_dd_Lbbar==0.) nSig_dd_Lbbar.setConstant();
  if (nSigEvents_dd+nSigEvents_dd_Lbbar==0.) {
    sig_dd_alpha1.setConstant();
    sig_dd_alpha2.setConstant();
    sig_dd_frac.setConstant();
    sig_dd_mean.setConstant();
    sig_dd_n1.setConstant();
    sig_dd_n2.setConstant();
    sig_dd_sigma.setConstant();
  }
  
  if (nSigEvents_ll==0.) nSig_ll.setConstant();
  if (nSigEvents_ll_Lbbar==0.) nSig_ll_Lbbar.setConstant();
  if (nSigEvents_ll+nSigEvents_ll_Lbbar==0.) {
    sig_ll_alpha1.setConstant();
    sig_ll_alpha2.setConstant();
    sig_ll_frac.setConstant();
    sig_ll_mean.setConstant();
    sig_ll_n1.setConstant();
    sig_ll_n2.setConstant();
    sig_ll_sigma.setConstant();
  }

  if (nBkgEvents_dd==0.) {
    nBkg_dd.setConstant();
    bkg_dd_c0.setConstant();
  }
  if (nBkgEvents_dd_Lbbar==0.) {
    nBkg_dd_Lbbar.setConstant();
    bkg_dd_Lbbar_c0.setConstant();
  }
  if (nBkgEvents_ll==0.) {
    nBkg_ll.setConstant();
    bkg_ll_c0.setConstant();
  }
  if (nBkgEvents_ll_Lbbar==0.) {
    nBkg_ll_Lbbar.setConstant();
    bkg_ll_Lbbar_c0.setConstant();
  }



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

  //  RooLbtoJpsiL0PDF3* sig_ang_model = new RooLbtoJpsiL0PDF3("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1);
  //  RooLbtoJpsiL0PDF3* sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF3("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1);
  // RooLbtoJpsiL0PDF3* sig_ang_model_fp = new RooLbtoJpsiL0PDF3("sig_ang_model_fp","sig_ang_model_fp",costheta,costheta1,costheta2,P_b, alpha_b_rfv, alpha_lambda, r_0_rfv, r_1_rfv);
  // RooLbtoJpsiL0PDF3* sig_ang_model_Lbbar_fp = new RooLbtoJpsiL0PDF3("sig_ang_model_Lbbar_fp","sig_ang_model_Lbbar_fp",costheta,costheta1,costheta2,P_b, alpha_b_rfv, alpha_lambda, r_0_rfv, r_1_rfv);

  // RooLbtoJpsiL0PDF3v2* sig_ang_model = new RooLbtoJpsiL0PDF3v2("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, beta, gammarrv);
  // RooLbtoJpsiL0PDF3v2* sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF3v2("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, beta, gammarrv);


  RooProdPdf bkg_dd_model("bkg_dd_model","bkg_dd_model",RooArgList(bkg_dd_mass_model,bkg_ang_model));
  RooProdPdf bkg_ll_model("bkg_ll_model","bkg_ll_model",RooArgList(bkg_ll_mass_model,bkg_ang_model));
  RooProdPdf bkg_dd_Lbbar_model("bkg_dd_Lbbar_model","bkg_dd_Lbbar_model",RooArgList(bkg_dd_Lbbar_mass_model,bkg_ang_model));
  RooProdPdf bkg_ll_Lbbar_model("bkg_ll_Lbbar_model","bkg_ll_Lbbar_model",RooArgList(bkg_ll_Lbbar_mass_model,bkg_ang_model));


  RooCategory cat("cat","cat");  
  std::vector<TString> catlist;
  if (nSigEvents_dd + nBkgEvents_dd != 0) catlist.push_back("ddLb"   );
  if (nSigEvents_ll + nBkgEvents_ll != 0) catlist.push_back("llLb"   );
  if (nSigEvents_dd_Lbbar + nBkgEvents_dd_Lbbar != 0) catlist.push_back("ddLbbar");
  if (nSigEvents_ll_Lbbar + nBkgEvents_ll_Lbbar != 0) catlist.push_back("llLbbar");
  for (unsigned int i=0;i<catlist.size();i++) cat.defineType(catlist[i],i+1);
  

  std::map<std::string,RooAbsPdf*> masspdfmapping;
  if (nSigEvents_dd + nBkgEvents_dd != 0) masspdfmapping[ "ddLb"    ] = &tot_dd_mass_model;
  if (nSigEvents_ll + nBkgEvents_ll != 0) masspdfmapping[ "llLb"    ] = &tot_ll_mass_model;
  if (nSigEvents_dd_Lbbar + nBkgEvents_dd_Lbbar != 0) masspdfmapping[ "ddLbbar" ] = &tot_dd_Lbbar_mass_model;
  if (nSigEvents_ll_Lbbar + nBkgEvents_ll_Lbbar != 0) masspdfmapping[ "llLbbar" ] = &tot_ll_Lbbar_mass_model;
  
  
  std::map<std::string,RooAbsPdf*> angpdfmapping;
  if (nSigEvents_dd + nBkgEvents_dd != 0) angpdfmapping[ "ddLb"    ] = sig_ang_model;
  if (nSigEvents_ll + nBkgEvents_ll != 0) angpdfmapping[ "llLb"    ] = sig_ang_model;
  if (nSigEvents_dd_Lbbar + nBkgEvents_dd_Lbbar != 0) angpdfmapping[ "ddLbbar" ] = sig_ang_model_Lbbar;
  if (nSigEvents_ll_Lbbar + nBkgEvents_ll_Lbbar != 0) angpdfmapping[ "llLbbar" ] = sig_ang_model_Lbbar;

  // std::map<std::string,RooAbsPdf*> angpdfmapping_fp;
  // if (nSigEvents_dd + nBkgEvents_dd != 0)  angpdfmapping_fp[ "ddLb"    ] = sig_ang_model_fp;
  // if (nSigEvents_ll + nBkgEvents_ll != 0) angpdfmapping_fp[ "llLb"    ] = sig_ang_model_fp;
  // if (nSigEvents_dd_Lbbar + nBkgEvents_dd_Lbbar != 0) angpdfmapping_fp[ "ddLbbar" ] = sig_ang_model_Lbbar_fp;
  // if (nSigEvents_ll_Lbbar + nBkgEvents_ll_Lbbar != 0) angpdfmapping_fp[ "llLbbar" ] = sig_ang_model_Lbbar_fp;
  
  RooSimultaneous masssimpdf("masssimpdf","masssimpdf",masspdfmapping,cat);
  RooSimultaneous angsimpdf("angsimpdf","angsimpdf",angpdfmapping,cat);
  //  RooSimultaneous angsimpdf_fp("angsimpdf_fp","angsimpdf_fp",angpdfmapping_fp,cat);

  RooAbsPdf&  tot_mass_model = masssimpdf;
  RooAbsPdf&  ang_model      = angsimpdf;
  //  RooAbsPdf&  ang_model_fp   = angsimpdf_fp;


  RooRealVar wmass("wmass","wmass",0.);
  RooRealVar wacc("wacc","wacc",0.);
  RooRealVar wtot("wtot","wtot",0.);

  //  RooArgList vars(P_b,alpha_b,r_0,r_1);
  RooArgList vars(P_b);
  if (useparam==0) vars.add(RooArgList(ap2,am2,bp2));
  if (useparam==1) vars.add(RooArgList(alpha_b,r_0,r_1));
  if (useparam==2) vars.add(RooArgList(alpha_b,beta,gammarrv));
  RooArgList varsfp(P_b,ap2,am2,bp2);

  if (floatalphalambda) {
    vars.add(RooArgList(alpha_lambda));
    varsfp.add(RooArgList(alpha_lambda));
  }
  
  std::vector<angfitstr> angfitv;
  angfitstr angfit;

  angfit.fitname="fr_angfit_woutsumw2";
  angfit.model=&angsimpdf;
  angfit.model4mc=sig_ang_model;
  angfit.vars=&vars;
  angfit.sumw2arg=RooFit::SumW2Error(false);
  angfit.minosarg=RooFit::Minos(true);
  //  angfitv.push_back(angfit);

  angfit.fitname="fr_angfit_woutsumw2_woutminos";
  angfit.sumw2arg=RooFit::SumW2Error(false);
  angfit.minosarg=RooFit::Minos(false);
  //  angfitv.push_back(angfit);

  // angfit.fitname="fr_angfit_wsumw2";
  // angfit.sumw2arg=RooFit::SumW2Error(inclbkg or inclacc);
  // angfit.minosarg=RooFit::Minos(false);
  // angfitv.push_back(angfit);

  // angfit.fitname="fr_angfit_woutsumw2_fp";
  // angfit.model=&angsimpdf_fp;
  // angfit.model4mc=sig_ang_model_fp;
  // angfit.vars=&varsfp;
  // angfit.sumw2arg=RooFit::SumW2Error(false);
  // angfit.minosarg=RooFit::Minos(true);
  // angfitv.push_back(angfit);

  // angfit.fitname="fr_angfit_woutsumw2_woutminos_fp";
  // angfit.sumw2arg=RooFit::SumW2Error(false);
  // angfit.minosarg=RooFit::Minos(false);
  // angfitv.push_back(angfit);
  
  // angfit.fitname="fr_angfit_wsumw2_fp";
  // angfit.sumw2arg=RooFit::SumW2Error(inclbkg or inclacc);
  // angfit.minosarg=RooFit::Minos(false);
  // angfitv.push_back(angfit);
  
  
  if (mode=="merge") {

    RooRealVar* means[vars.getSize()];
    RooRealVar* sigmas[vars.getSize()];
    //    RooGaussian* gauss[vars.getSize()];
    for (int i=0;i<vars.getSize();i++) {
      RooRealVar* var = (RooRealVar*)vars.at(i);
      const TString varname(var->GetName());
      means[i] = new RooRealVar(varname + "mean", varname + "mean",0.,-5.,5.); 
      sigmas[i] = new RooRealVar(varname + "sigma", varname + "sigma",1.,0.,3.); 
      //      gauss[i] = new RooGaussian(varname + "gauss", varname + "gauss",*var,*means[i],*sigmas[i]);
    }


    RooRealVar Pbcorr("Pbcorr","Pbcorr",0.);
    RooRealVar abcorr("abcorr","abcorr",0.);
    RooRealVar r0corr("r0corr","r0corr",0.);
    RooRealVar r1corr("r1corr","r1corr",0.);

    //    RooFormulaVar* ap2corr_rfv = new RooFormulaVar("ap2corr_rrv","|a+|^{2}","0.5*(@0+@1)",RooArgList(r0corr,r1corr));
    //    RooFormulaVar* am2corr_rfv = new RooFormulaVar("am2corr_rrv","|a-|^{2}","0.5*(@0-@1)",RooArgList(r0corr,r1corr));
    //    RooFormulaVar* bp2corr_rfv = new RooFormulaVar("bp2corr_rrv","|b+|^{2}","0.5*(1.+@0-@1-@2)",RooArgList(abcorr,r0corr,r1corr));
    //    RooFormulaVar* bm2corr_rfv = new RooFormulaVar("bm2corr_rrv","|b-|^{2}","0.5*(1.-@0-@1+@2)",RooArgList(abcorr,r0corr,r1corr));


    RooDataSet parscorrected("parscorrected","parscorrected",RooArgList(Pbcorr,abcorr,r0corr,r1corr,ap2corr,am2corr,bp2corr,bm2corr));

    RooDataSet amplitudespulls("amplitudespulls","amplitudespulls",RooArgList(ap2,am2,bp2,bm2_rrv));

    TMatrixDSym Ccorr(4);
    for (int i=0;i<4;i++) for (int j=0;j<4;j++) Ccorr(i,j)=0.;
    Ccorr(0,0)=1.;
    Ccorr(1,1)=pow(1.7,2.);
    Ccorr(2,2)=1.;
    Ccorr(3,3)=pow(1.7,2.);


    /*
      |      0    |      1    |      2    |      3    |
      ---------------------------------------------------------
      0 |      1.144      0.3669      0.1493      0.3579 
      1 |     0.3669       2.837      0.5749       2.695 
      2 |     0.1493      0.5749       1.074       0.575 
      3 |     0.3579       2.695       0.575       2.905 
    */


    //     //      |      0    |      1    |      2    |      3    |
    //     //      ---------------------------------------------------------
    //     //      0 |     
    //     Ccorr(0,0)=1.144;      Ccorr(0,1)=0.3669;      Ccorr(0,2)=0.1493;      Ccorr(0,3)=0.3579; 
    //     //      1 |     
    //     Ccorr(1,0)=0.3669;       Ccorr(1,1)=2.837;      Ccorr(1,2)=0.5749;       Ccorr(1,3)=2.695; 
    //     //      2 |     
    //     Ccorr(2,0)=0.1493;      Ccorr(2,1)=0.5749;       Ccorr(2,2)=1.074;       Ccorr(2,3)=0.575; 
    //     //      3 |     
    //     Ccorr(3,0)=0.3579;       Ccorr(3,1)=2.695;       Ccorr(3,2)=0.575;       Ccorr(3,3)=2.905; 


    Ccorr.Print("v");
    mygetchar;
    

 //    Ccorr(0,0)=1.144;
//     Ccorr(1,1)=2.837;
//     Ccorr(2,2)=1.074;
//     Ccorr(3,3)=2.905;
//     Ccorr(1,0)=0.3669;
//     Ccorr(2,0)=0.1493;
//     Ccorr(3,0)=0.3579;
//     Ccorr(2,1)=0.5749;
//     Ccorr(3,1)=2.695;
//     Ccorr(3,2)=0.575;


    TMatrixD Ccorrinv(Ccorr);
    Double_t det=0.;
    Ccorrinv.Invert(&det) ;

    TMatrixD M(4,4);
    for (int i=0;i<4;i++) for (int j=0;j<4;j++) M(i,j)=0.;
    // M(2,0)=+0.5;
    // M(3,0)=+0.5;
    // M(2,1)=+0.5;
    // M(3,1)=-0.5;
    // M(1,2)=+0.5;
    // M(2,2)=-0.5;
    // M(3,2)=+0.5;
    // M(1,3)=-0.5;
    // M(2,3)=-0.5;
    // M(3,3)=+0.5;

    //    mygetchar;

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

    RooArgList massvars;
    if (nSigEvents_dd + nBkgEvents_dd != 0) massvars.add(RooArgList(nBkg_dd, bkg_dd_c0, nSig_dd));
    if (nSigEvents_ll + nBkgEvents_ll != 0) massvars.add(RooArgList(nBkg_ll, bkg_ll_c0, nSig_ll));
    if (nSigEvents_dd_Lbbar + nBkgEvents_dd_Lbbar != 0) massvars.add(RooArgList(nBkg_dd_Lbbar, nSig_dd_Lbbar, bkg_dd_Lbbar_c0));
    if (nSigEvents_ll_Lbbar + nBkgEvents_ll_Lbbar != 0) massvars.add(RooArgList(nBkg_ll_Lbbar, nSig_ll_Lbbar, bkg_ll_Lbbar_c0));
    if (nSigEvents_dd + nBkgEvents_dd + nSigEvents_dd_Lbbar + nBkgEvents_dd_Lbbar!= 0) massvars.add(RooArgList(sig_dd_mean,sig_dd_sigma));
    if (nSigEvents_ll + nBkgEvents_ll + nSigEvents_ll_Lbbar + nBkgEvents_ll_Lbbar!= 0) massvars.add(RooArgList(sig_ll_mean,sig_ll_sigma));
    if (floatSigParams) {
      massvars.add(RooArgList(sig_dd_alpha1, sig_dd_alpha2, sig_dd_frac, sig_dd_n1, sig_dd_n2));
      massvars.add(RooArgList(sig_ll_alpha1, sig_ll_alpha2, sig_ll_frac, sig_ll_n1, sig_ll_n2));
    }

    RooArgList fitvars(costheta,costheta1,costheta2);
    RooMCStudy massmcstudy(masssimpdf,RooArgList(mass));
    for (unsigned int i=0; i<angfitv.size(); i++) angfitv[i].mcstudy = new RooMCStudy(*angfitv[i].model, fitvars);
    RooMCStudy angmcstudy(angsimpdf,fitvars);

    int ntot[angfitv.size()];
    int nok[angfitv.size()];
    memset(ntot,0,angfitv.size()*sizeof(int));
    memset(nok, 0,angfitv.size()*sizeof(int));
    for (int i=0; i<iSample;i++) {
      const TString filename="toy-" + istr(i,"07") + "-result.root";
      TFile file(filename);
      if (file.IsZombie()) continue;
      RooFitResult* fr_massfit = (RooFitResult*)file.Get("fr_massfit");
      
      //      checkweirdpulls(fr_massfit, &massvars);
      if (isconverged(fr_massfit)) {
	massmcstudy.addFitResult(*fr_massfit);
	for (unsigned int j=0; j<angfitv.size(); j++) {
	  RooFitResult* fr = (RooFitResult*)file.Get(angfitv[j].fitname);//getRooFitResult(filename, angfitv[j].fitname);
	  if (fr) {
	    //	    std::cout << filename << ":" <<  angfitv[j].fitname << " : status=" << status << " ; covQual=" << covQual << std::endl;
	    ntot[j]++;
	    if (isconverged(fr) and checkweirdpulls(fr, angfitv[j].vars)) { // status==0 and (covQual==3 or covQual==-1) 
	      nok[j]++;
	      angfitv[j].mcstudy->addFitResult(*fr);
	    }
	  }
	}
	RooFitResult* fr = (RooFitResult*)file.Get("fr_angfit_myfitto");
	
	if (fr and isconverged(fr) and checkweirdpulls(fr, &vars)) {
	  angmcstudy.addFitResult(*fr);

	  ap2=ap2_rfv->getVal();
	  am2=am2_rfv->getVal();
	  bp2=bp2_rfv->getVal();
	  bm2_rrv=bm2_rfv->getVal();

	  ap2.setError(ap2_rfv->getPropagatedError(*fr));
	  am2.setError(am2_rfv->getPropagatedError(*fr));
	  bp2.setError(bp2_rfv->getPropagatedError(*fr));
	  bm2_rrv.setError(bm2_rfv->getPropagatedError(*fr));

	  amplitudespulls.add(RooArgSet(ap2,am2,bp2,bm2_rrv));

	  //	  TMatrixDSym Covampl(4);
	  // Calculate corrected covariance matrix = V C-1 V
	  //	  Covampl = M * Ccorr;// * MT;

	  TMatrixD Covcorr(fr->covarianceMatrix(),TMatrixD::kMult,Ccorr);
	  //	  TMatrixD Covcorr(fr->covarianceMatrix(),TMatrixD::kMult,TMatrixD(Ccorr,TMatrixD::kMult,fr->covarianceMatrix()));

	  Int_t n =Covcorr.GetNrows() ;
	  TMatrixDSym Covcorrsym(n);
	  for (Int_t i=0 ; i<n ; i++) {
	    for (Int_t j=i ; j<n ; j++) {
	      if (i==j) {
		Covcorrsym(i,j) = Covcorr(i,j) ;
	      }
	      if (i!=j) {
		Double_t deltaRel = (Covcorr(i,j)-Covcorr(j,i))/sqrt(Covcorr(i,i)*Covcorr(j,j)) ;
		if (fabs(deltaRel)>1e-3) {
		  //		  coutW(Fitting) << "RooAbsPdf::fitTo(" << GetName() << ") WARNING: Corrected covariance matrix is not (completely) symmetric: V[" << i << "," << j << "] = " 
		  //				 << Covcorr(i,j) << " V[" << j << "," << i << "] = " << Covcorr(j,i) << " explicitly restoring symmetry by inserting average value" << endl ;
		}
		Covcorrsym(i,j) = (Covcorr(i,j)+Covcorr(j,i))/2 ;
	      }
	    }
	  }
	  Pbcorr.setError(sqrt(Covcorrsym(0,0)));
	  abcorr.setError(sqrt(Covcorrsym(1,1)));
	  r0corr.setError(sqrt(Covcorrsym(2,2)));
	  r1corr.setError(sqrt(Covcorrsym(3,3)));

	  //	  Pbcorr.setError(P_b.getError());
	  //	  abcorr.setError(alpha_b.getError()*1.7);
	  //	  r0corr.setError(r_0.getError());
	  //	  r1corr.setError(r_1.getError()*1.7);

	  Pbcorr=P_b.getVal()+0.04*Pbcorr.getError();
	  abcorr=alpha_b.getVal()+0.20*abcorr.getError();
	  r0corr=r_0.getVal()+0.21*r0corr.getError();
	  r1corr=r_1.getVal()+0.13*r1corr.getError();

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

	  //	  ap2corr=ap2corr_rfv->getVal();
	  //	  am2corr=am2corr_rfv->getVal();
	  //	  bp2corr=bp2corr_rfv->getVal();
	  //	  bm2corr=bm2corr_rfv->getVal();
	  
	  //	  Covampl = M * fr->covarianceMatrix() * Ccorr * M.T();
	  TMatrixD Covampl(M,TMatrixD::kMult,TMatrixD(Covcorr,TMatrixD::kMult,MT)) ;
	  
	  ap2corr.setError(sqrt(Covampl(0,0)));
	  am2corr.setError(sqrt(Covampl(1,1)));
	  bp2corr.setError(sqrt(Covampl(2,2)));
	  bm2corr.setError(sqrt(Covampl(3,3)));
	  //	  am2corr.setError(am2corr_rfv->getPropagatedError(*fr));
	  //	  bp2corr.setError(bp2corr_rfv->getPropagatedError(*fr));
	  //	  bm2corr.setError(bm2corr_rfv->getPropagatedError(*fr));

	  parscorrected.add(RooArgSet(Pbcorr,abcorr,r0corr,r1corr,ap2corr,am2corr,bp2corr,bm2corr));
	}
      } else {
	std::cout << filename << ": mass fit failed!!!" << std::endl; 
      }
      file.Close();
    }
    //    mygetchar;      

    plot_mcstudy(c,massmcstudy, "mcstudy-massfit", massvars);
    plot_mcstudy(c,angmcstudy, "mcstudy-angfit_myfitto", vars);
    {
      const RooDataSet& data = angmcstudy.fitParDataSet();
      RooDataSet data2(data,"data2");

      data2.Print("v");
      data2.write("data2.txt");

      RooArgList blu(*findvar(&data2,"P_bpull"),*findvar(&data2,"alpha_bpull"),*findvar(&data2,"r_0pull"), *findvar(&data2,"r_1pull"));
      blu.Print("v");

      TMatrixDSym* cov = data2.covarianceMatrix(blu);
      std::cout << "Covariance between pulls" << std::endl;
      cov->Print("v");
      mygetchar;

      TMatrixDSym* corrvars = data2.correlationMatrix(vars);
      std::cout << "Correlation between vars" << std::endl;
      corrvars->Print("v");
      mygetchar;

      TMatrixDSym* corr = data2.correlationMatrix(blu);
      std::cout << "Correlation between pulls" << std::endl;
      corr->Print("v");
      mygetchar;



      RooRealVar* ap2_rrv = (RooRealVar*)data2.addColumn(*ap2_rfv);
      RooRealVar* am2_rrv = (RooRealVar*)data2.addColumn(*am2_rfv);
      RooRealVar* bp2_rrv = (RooRealVar*)data2.addColumn(*bp2_rfv);
      RooRealVar* bm2_rrv2 = (RooRealVar*)data2.addColumn(*bm2_rfv);

      addpull(&amplitudespulls,&ap2, ap2_initval);
      addpull(&amplitudespulls,&am2, am2_initval);
      addpull(&amplitudespulls,&bp2, bp2_initval);
      addpull(&amplitudespulls,&bm2_rrv, bm2_initval);
      
      plotdata(c, RooArgList(*ap2_rrv, *am2_rrv, *bp2_rrv, *bm2_rrv2), &data2, "mcstudy-angfit_myfitto-ampl");

      for (int i=0;i<vars.getSize();i++) {
	RooRealVar* var = (RooRealVar*)vars.at(i);
	plotPull(c, "mcstudy-angfit_myfitto", *var,&data2, -7., 7., 25, true,means[i],sigmas[i]);
      }
      
      TVectorD meansvec(4);
      for (int i=0;i<vars.getSize();i++) meansvec[i]=means[i]->getVal();
      TFile afile("corr.root","recreate");
      //      meansvec.SetName("meansvector");
      //      cov->SetName("covmatrix");
      meansvec.Write("meansvector");
      cov->Write("covmatrix");
      afile.Close();

      RooArgSet tmp(Pb_initval, ab_initval, r0_initval, r1_initval);
      for (int i=0;i<vars.getSize();i++) tmp.add(RooArgSet(*means[i],*sigmas[i]));
      tmp.writeToFile("mcstudy-angfit_myfitto-pulls.txt");
      tmp.printLatex(latexformat, RooFit::OutputFile("mcstudy-angfit_myfitto-pulls.tex"), RooFit::Columns(1));

      plotPull(c,"mcstudy-angfit_myfitto-ampl", ap2, &amplitudespulls, -7., 7., 25, true);
      plotPull(c,"mcstudy-angfit_myfitto-ampl", am2, &amplitudespulls, -7., 7., 25, true);
      plotPull(c,"mcstudy-angfit_myfitto-ampl", bp2, &amplitudespulls, -7., 7., 25, true);
      plotPull(c,"mcstudy-angfit_myfitto-ampl", bm2_rrv, &amplitudespulls, -7., 7., 25, true);


      addpull(&parscorrected,&Pbcorr, Pb_initval);
      addpull(&parscorrected,&abcorr, ab_initval);
      addpull(&parscorrected,&r0corr, r0_initval);
      addpull(&parscorrected,&r1corr, r1_initval);
      addpull(&parscorrected,&ap2corr, ap2_initval);
      addpull(&parscorrected,&am2corr, am2_initval);
      addpull(&parscorrected,&bp2corr, bp2_initval);
      addpull(&parscorrected,&bm2corr, bm2_initval);
      

      plotdata(c, RooArgList(ap2corr, am2corr, bp2corr, bm2corr), &parscorrected, "mcstudy-angfit_myfitto-corrected-ampl");

      plotPull(c,"mcstudy-angfit_myfitto-corrected-init", Pbcorr, &parscorrected, -7., 7., 25, true);
      plotPull(c,"mcstudy-angfit_myfitto-corrected-init", abcorr, &parscorrected, -7., 7., 25, true);
      plotPull(c,"mcstudy-angfit_myfitto-corrected-init", r0corr, &parscorrected, -7., 7., 25, true);
      plotPull(c,"mcstudy-angfit_myfitto-corrected-init", r1corr, &parscorrected, -7., 7., 25, true);

      plotPull(c,"mcstudy-angfit_myfitto-corrected-ampl", ap2corr, &parscorrected, -7., 7., 25, true);
      plotPull(c,"mcstudy-angfit_myfitto-corrected-ampl", am2corr, &parscorrected, -7., 7., 25, true);
      plotPull(c,"mcstudy-angfit_myfitto-corrected-ampl", bp2corr, &parscorrected, -7., 7., 25, true);
      plotPull(c,"mcstudy-angfit_myfitto-corrected-ampl", bm2corr, &parscorrected, -7., 7., 25, true);

    }


    for (unsigned int j=0; j<angfitv.size(); j++) {
      const RooDataSet& data = angfitv[j].mcstudy->fitParDataSet();
      RooDataSet data2(data,"data2");
      data2.write("mcstudy-" + angfitv[j].fitname + ".txt");
      plot_mcstudy(c,*angfitv[j].mcstudy, "mcstudy-" + angfitv[j].fitname, *angfitv[j].vars);
      
      data2.get()->Print("v");
      RooRealVar* P_bpull = (RooRealVar*)data2.get()->find("P_bpull");
      RooRealVar mean("mean","mean",0., -1., 1.);
      RooRealVar sigmal("sigmal","#sigma_{L}",1.6, 0.1, 2.);
      RooRealVar sigmar("sigmar","#sigma_{R}",1.6, 0.1, 2.);
      RooBifurGauss bfgauss("bfgauss","bfgauss",*P_bpull, mean, sigmal, sigmar);
      RooFitResult* fr = bfgauss.fitTo(data2,RooFit::Save(true));
      fr->Print("v");
      P_bpull->setMin(-7.);
      P_bpull->setMax(7.);
      //      plot(c,RooArgList(*P_bpull),bfgauss,&data2, "mcstudy-bfgauss-" + angfitv[j].fitname);
      
      RooPlot* frame = P_bpull->frame() ;
      data2.plotOn(frame) ;
      bfgauss.plotOn(frame) ;
      bfgauss.paramOn(frame,&data2) ;
      frame->Draw();
      c.SaveAs(TString("mcstudy-bfgauss-P_bpull-") + angfitv[j].fitname + ".eps");
      

      //      RooDataSet* datared = (RooDataSet*)data2.reduce(TString(bm2_rrv->GetName()) + ">-0.091");


      
    }
    
    for (unsigned int j=0; j<angfitv.size(); j++) {
      std::cout << angfitv[j].fitname << ": "<<  nok[j] << "/" << ntot[j] << std::endl;
    }

    {
      const RooDataSet& datainit = angmcstudy.fitParDataSet();
      RooDataSet data(datainit,"data2");
      RooRealVar* ap2_rrv = (RooRealVar*)data.addColumn(*ap2_rfv);
      RooRealVar* am2_rrv = (RooRealVar*)data.addColumn(*am2_rfv);
      RooRealVar* bp2_rrv = (RooRealVar*)data.addColumn(*bp2_rfv);
      RooRealVar* bm2_rrv = (RooRealVar*)data.addColumn(*bm2_rfv);

      plotdata(c, RooArgList(*ap2_rrv, *am2_rrv, *bp2_rrv, *bm2_rrv), &data, "mcstudy-fuck");
      
      //      std::cout << "N(bm2<0.091) = " << 100.*data.reduce(TString(bm2_rrv->GetName()) + "<-0.091")->numEntries()/data.numEntries() << " %" << std::endl;
      std::cout << "Ntoys            = " << data.numEntries() << std::endl;
      std::cout << "Ntoys (physical) = " << countphysicalsolutions(data) << std::endl;

      std::cout << "N(Pb>20%) = " << 100.*data.reduce("P_b>0.20")->numEntries()/data.numEntries() << " %" << std::endl;
      std::cout << "N(ab>49%) = " << 100.*data.reduce("alpha_b>0.49")->numEntries()/data.numEntries() << " %" << std::endl;
      std::cout << "N(ab>77%) = " << 100.*data.reduce("alpha_b>0.77")->numEntries()/data.numEntries() << " %" << std::endl;
      std::cout << "N(bm2<0.078) = " << 100.*data.reduce(TString(bm2_rrv->GetName()) + "<-0.078")->numEntries()/data.numEntries() << " %" << std::endl;

    }

    return 0;
  }



  TString acceptancefile="fitacceptance3";
  if (userwmc) {
    acceptancefile+="-mcrw";
  } else {
    acceptancefile+="-mcnonrw";
  }
  if (usebdt) acceptancefile+="-wbdt";
  if (usetmvarectcut) acceptancefile+="-wtmvarectcut";
  if (doBd2JpsiKS0) {
    acceptancefile+="-B0-";
  } else {
    acceptancefile+="-Lb-";
  }
  acceptancefile+=directorylist;
  calcacceptanceclass* acceptancedd=NULL;
  calcacceptanceclass* acceptancell=NULL;
  if (inclacc) {
    acceptancedd = new calcacceptanceclass(acceptancefile + "-dd.root"); 
    acceptancell = new calcacceptanceclass(acceptancefile + "-ll.root");
  }

  //  calcacceptanceclass* acceptancebddd= new calcacceptanceclass("fitacceptance3-mcrw-wbdt-B0-Tuple_JpsiKs_detached_betas-dd.root");
  //  calcacceptanceclass* acceptancebdll= new calcacceptanceclass("fitacceptance3-mcrw-wbdt-B0-Tuple_JpsiKs_detached_betas-ll.root");

  RooDataSet* sig_dd_ang_data=NULL;
  RooDataSet* sig_dd_Lbbar_ang_data=NULL;
  RooDataSet* sig_ll_ang_data=NULL;
  RooDataSet* sig_ll_Lbbar_ang_data=NULL;
  if (inclacc) {
    sig_dd_ang_data = generatetoy("sig_dd_ang_data","sig_dd_ang_data",*sig_ang_model,RooArgSet(costheta,costheta1,costheta2), nSigEvents_dd, acceptancedd, rnd);
    sig_dd_Lbbar_ang_data = generatetoy("sig_dd_Lbbar_ang_data","sig_dd_Lbbar_ang_data",*sig_ang_model_Lbbar,RooArgSet(costheta,costheta1,costheta2), nSigEvents_dd_Lbbar, acceptancedd, rnd);
    sig_ll_ang_data = generatetoy("sig_ll_ang_data","sig_ll_ang_data",*sig_ang_model,RooArgSet(costheta,costheta1,costheta2), nSigEvents_ll, acceptancell, rnd);
    sig_ll_Lbbar_ang_data = generatetoy("sig_ll_Lbbar_ang_data","sig_ll_Lbbar_ang_data",*sig_ang_model_Lbbar,RooArgSet(costheta,costheta1,costheta2), nSigEvents_ll_Lbbar, acceptancell, rnd);

    // piminus_TRACK_Type=5.;
    // sig_dd_ang_data->addColumn(piminus_TRACK_Type);
    // addacceptanceweight(sig_dd_ang_data, wacc,false,false,acceptancell,acceptancedd);
    // RooDataSet test("test","test", sig_dd_ang_data, *sig_dd_ang_data->get(), "", wacc.GetName());
    // plotdata(c, RooArgList(costheta,costheta1,costheta2), &test, "test-ll");
    // return 1;

  } else {
    sig_dd_ang_data = generatetoy("sig_dd_ang_data","sig_dd_ang_data",*sig_ang_model,RooArgSet(costheta,costheta1,costheta2), nSigEvents_dd);
    sig_dd_Lbbar_ang_data = generatetoy("sig_dd_Lbbar_ang_data","sig_dd_Lbbar_ang_data",*sig_ang_model_Lbbar,RooArgSet(costheta,costheta1,costheta2), nSigEvents_dd_Lbbar);
    sig_ll_ang_data = generatetoy("sig_ll_ang_data","sig_ll_ang_data",*sig_ang_model,RooArgSet(costheta,costheta1,costheta2), nSigEvents_ll);
    sig_ll_Lbbar_ang_data = generatetoy("sig_ll_Lbbar_ang_data","sig_ll_Lbbar_ang_data",*sig_ang_model_Lbbar,RooArgSet(costheta,costheta1,costheta2), nSigEvents_ll_Lbbar);
  }

  RooDataSet* sig_dd_mass_data       = generatetoy("sig_dd_mass_data",      "sig_dd_mass_data",      sig_dd_mass_model,RooArgSet(mass),nSigEvents_dd);
  RooDataSet* sig_dd_Lbbar_mass_data = generatetoy("sig_dd_Lbbar_mass_data","sig_dd_Lbbar_mass_data",sig_dd_mass_model,RooArgSet(mass),nSigEvents_dd_Lbbar);
  RooDataSet* sig_ll_mass_data       = generatetoy("sig_ll_mass_data",      "sig_ll_mass_data",      sig_ll_mass_model,RooArgSet(mass),nSigEvents_ll);
  RooDataSet* sig_ll_Lbbar_mass_data = generatetoy("sig_ll_Lbbar_mass_data","sig_ll_Lbbar_mass_data",sig_ll_mass_model,RooArgSet(mass),nSigEvents_ll_Lbbar);

  //  if (debug) {
  //    sig_ang_data->Print("v");
  //    mygetchar;
  //    sig_mass_data->Print("v");
  //    mygetchar;
  //  }
  sig_dd_ang_data->merge(sig_dd_mass_data);
  sig_dd_Lbbar_ang_data->merge(sig_dd_Lbbar_mass_data);
  sig_ll_ang_data->merge(sig_ll_mass_data);
  sig_ll_Lbbar_ang_data->merge(sig_ll_Lbbar_mass_data);
  if (debug) {
    //    sig_ang_data->Print("v");
    mygetchar;
  }

  if (inclbkg) {
    bkg_costheta_c0.setVal(rnd->Rndm()/2.-0.25);
    bkg_costheta_c1.setVal(rnd->Rndm()/2.-0.25);
    bkg_costheta_c2.setVal(rnd->Rndm()/2.-0.25);
    bkg_costheta1_c0.setVal(rnd->Rndm()/2.-0.25);
    bkg_costheta1_c1.setVal(rnd->Rndm()/2.-0.25);
    bkg_costheta1_c2.setVal(rnd->Rndm()/2.-0.25);
    bkg_costheta2_c0.setVal(rnd->Rndm()/2.-0.25);
    bkg_costheta2_c1.setVal(rnd->Rndm()/2.-0.25);
    bkg_costheta2_c2.setVal(rnd->Rndm()/2.-0.25);

    
    const double fracBdll=0.;
    const double fracBddd=0.1;

    RooDataSet* bkg_dd_data       = generatetoy("bkg_dd_data",      "bkg_dd_data",      bkg_dd_model,      RooArgSet(mass,costheta,costheta1,costheta2),nBkgEvents_dd * (1.-fracBddd));
    RooDataSet* bkg_dd_Lbbar_data = generatetoy("bkg_dd_Lbbar_data","bkg_dd_Lbbar_data",bkg_dd_Lbbar_model,RooArgSet(mass,costheta,costheta1,costheta2),nBkgEvents_dd_Lbbar* (1.-fracBddd));
    RooDataSet* bkg_ll_data       = generatetoy("bkg_ll_data",      "bkg_ll_data",      bkg_ll_model,      RooArgSet(mass,costheta,costheta1,costheta2),nBkgEvents_ll* (1.-fracBdll));
    RooDataSet* bkg_ll_Lbbar_data = generatetoy("bkg_ll_Lbbar_data","bkg_ll_Lbbar_data",bkg_ll_Lbbar_model,RooArgSet(mass,costheta,costheta1,costheta2),nBkgEvents_ll_Lbbar* (1.-fracBdll));
    

    RooDataSet* bkgbs_dd_data=generatetoy("bkgbs_dd_data",      "bkgbs_dd_data",      bkg_dd_model,      RooArgSet(mass),nBkgEvents_dd * fracBddd);
    RooDataSet* bkgbs_ll_data=generatetoy("bkgbs_ll_data",      "bkgbs_ll_data",      bkg_ll_model,      RooArgSet(mass),nBkgEvents_ll * fracBdll);

    //    getBdbkg(iSample,RooArgList(costheta,costheta1,costheta2), bkgbs_dd_data, nBkgEvents_dd * fracBddd, bkgbs_ll_data, nBkgEvents_ll * fracBdll);



    bkg_dd_data->append(*bkgbs_dd_data);
    bkg_ll_data->append(*bkgbs_ll_data);

    sig_dd_ang_data->append(*bkg_dd_data);
    sig_dd_Lbbar_ang_data->append(*bkg_dd_Lbbar_data);
    sig_ll_ang_data->append(*bkg_ll_data);
    sig_ll_Lbbar_ang_data->append(*bkg_ll_Lbbar_data);
  }


  RooDataSet* tot_dd_data = sig_dd_ang_data;
  RooDataSet* tot_dd_Lbbar_data = sig_dd_Lbbar_ang_data;
  RooDataSet* tot_ll_data = sig_ll_ang_data;
  RooDataSet* tot_ll_Lbbar_data = sig_ll_Lbbar_ang_data;

  if (inclacc) {
    piminus_TRACK_Type=5.;
    tot_dd_data->addColumn(piminus_TRACK_Type);
    tot_dd_Lbbar_data->addColumn(piminus_TRACK_Type);
    addacceptanceweight(tot_dd_data,       wacc,false,false,acceptancell,acceptancedd);
    addacceptanceweight(tot_dd_Lbbar_data, wacc,false,false,acceptancell,acceptancedd);

    piminus_TRACK_Type=3.;
    tot_ll_data->addColumn(piminus_TRACK_Type);
    tot_ll_Lbbar_data->addColumn(piminus_TRACK_Type);
    addacceptanceweight(tot_ll_data,       wacc,false,false,acceptancell,acceptancedd);
    addacceptanceweight(tot_ll_Lbbar_data, wacc,false,false,acceptancell,acceptancedd);
  }

  std::map<std::string,RooDataSet*> mapping;
  if (nSigEvents_dd + nBkgEvents_dd != 0)             mapping[ "ddLb"    ] = tot_dd_data;
  if (nSigEvents_ll + nBkgEvents_ll != 0)             mapping[ "llLb"    ] = tot_ll_data;
  if (nSigEvents_dd_Lbbar + nBkgEvents_dd_Lbbar != 0) mapping[ "ddLbbar" ] = tot_dd_Lbbar_data;
  if (nSigEvents_ll_Lbbar + nBkgEvents_ll_Lbbar != 0) mapping[ "llLbbar" ] = tot_ll_Lbbar_data;
  
  RooFitResult* fr_mass=NULL;
  
  RooDataSet simdata("simdata", "simdata", RooArgList(mass,costheta,costheta1,costheta2,wacc,cat), RooFit::Index(cat),RooFit::Import(mapping));
  simdata.Print("v");
  //  simdata.write("simdata.txt");
  RooDataSet* tot_data       = &simdata;

  if (inclbkg) {
    fr_mass = tot_mass_model.fitTo(*tot_data, RooFit::Save(true), RooFit::Minos(false),massfitconstrains);
    
    if ((iSample<10) and (doplot)) {
      const RooArgList pdflist(sig_dd_mass_model,sig_ll_mass_model,bkg_dd_mass_model,bkg_ll_mass_model,bkg_dd_Lbbar_mass_model,bkg_ll_Lbbar_mass_model);
      plot(c,RooArgList(mass), tot_mass_model, tot_data, "toy-" + istr(iSample,"07") + "-fr",false,&pdflist,"", &cat);
    }
    
    RooArgSet* blu = tot_mass_model.getVariables();
    TIterator* iter = blu->createIterator();
    //  for (int i=0;blu.getSize();i++) std::cout << blu.at(i)->GetName() << std::endl;
    RooRealVar* arg;
    //const everything except nSig, nBkg and _M
    //    std::vector<RooRealVar*> varstounconst;
    RooArgList yieldlist_dd;
    RooArgList yieldlist_dd_Lbbar;
    RooArgList yieldlist_ll;
    RooArgList yieldlist_ll_Lbbar;
    //    std::cout << "Making constants: ";
    while((arg=(RooRealVar*)iter->Next())) {
      TString varname(arg->GetName());
      // if (!varname.Contains("_M") and !varname.Contains("mass") and !varname.Contains("cat") and !varname.Contains("nSig") and !varname.Contains("nBkg") and (!arg->isConstant())) {
      // 	std::cout << varname << " ";
      // 	arg->setConstant();
      // 	varstounconst.push_back(arg);
      // }
      if (varname.Contains("nSig") or varname.Contains("nBkg")) {
	if (varname.Contains("ll")) {
	  if (varname.Contains("Lbbar")) {yieldlist_ll_Lbbar.add(*arg);} else {yieldlist_ll.add(*arg);}
	} else {
	  if (varname.Contains("Lbbar")) {yieldlist_dd_Lbbar.add(*arg);} else {yieldlist_dd.add(*arg);}
	}
      }
      
    }
    //    std::cout << std::endl;
    //      RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
    //						   *tot_data, &tot_mass_model, yieldlist );
    
    
    //      double Psample_dd      =calcPsample(tot_dd_data,       mass, &sig_dd_mass_model, &bkg_dd_mass_model,       nSig_dd.getVal()/nBkg_dd.getVal());
    //      double Psample_dd_Lbbar=calcPsample(tot_dd_Lbbar_data, mass, &sig_dd_mass_model, &bkg_dd_Lbbar_mass_model, nSig_dd_Lbbar.getVal()/nBkg_dd_Lbbar.getVal());
    //      double Psample_ll      =calcPsample(tot_ll_data,       mass, &sig_ll_mass_model, &bkg_ll_mass_model,       nSig_ll.getVal()/nBkg_ll.getVal());
    //      double Psample_ll_Lbbar=calcPsample(tot_ll_Lbbar_data, mass, &sig_ll_mass_model, &bkg_ll_Lbbar_mass_model, nSig_ll_Lbbar.getVal()/nBkg_ll_Lbbar.getVal());
    //    double Psample_dd      =1.;
    //    double Psample_dd_Lbbar=1.;
    //    double Psample_ll      =1.;
    //    double Psample_ll_Lbbar=1.;
    if (debug) mygetchar;
    
    if (tot_dd_data->numEntries()>0.) {
      RooStats::SPlot sData("sData","An SPlot", *tot_dd_data, &tot_dd_mass_model, yieldlist_dd );
      copyvar(tot_dd_data,"nSig_dd_sw",wmass);
      //copyvar(tot_dd_data,"L_nSig_dd",wmass,Psample_dd);
      if (debug) {
	sData.Print("v");
	mygetchar;
      }
    }
    if (tot_dd_Lbbar_data->numEntries()>0.) {
      RooStats::SPlot sData("sData","An SPlot",*tot_dd_Lbbar_data, &tot_dd_Lbbar_mass_model, yieldlist_dd_Lbbar );
      copyvar(tot_dd_Lbbar_data,"nSig_dd_Lbbar_sw",wmass);
      //copyvar(tot_dd_Lbbar_data,"L_nSig_dd_Lbbar",wmass,Psample_dd_Lbbar);
    }
    if (tot_ll_data->numEntries()>0.) {
      RooStats::SPlot sData("sData","An SPlot",*tot_ll_data, &tot_ll_mass_model, yieldlist_ll);
      copyvar(tot_ll_data,"nSig_ll_sw",wmass);
      //copyvar(tot_ll_data,"L_nSig_ll",wmass,Psample_ll);
    }
    if (tot_ll_Lbbar_data->numEntries()>0.) {
      RooStats::SPlot sData("sData","An SPlot", *tot_ll_Lbbar_data, &tot_ll_Lbbar_mass_model, yieldlist_ll_Lbbar );
      copyvar(tot_ll_Lbbar_data,"nSig_ll_Lbbar_sw",wmass);
      //copyvar(tot_ll_Lbbar_data,"L_nSig_ll_Lbbar",wmass,Psample_ll_Lbbar);
    }
    std::cout << "Sum of wmass in tot_dd_data       = " << sumofvarindata(tot_dd_data,"wmass") << std::endl;
    std::cout << "Sum of wmass in tot_dd_Lbbar_data = " << sumofvarindata(tot_dd_Lbbar_data,"wmass") << std::endl;
    std::cout << "Sum of wmass in tot_ll_data       = " << sumofvarindata(tot_ll_data,"wmass") << std::endl;
    std::cout << "Sum of wmass in tot_ll_Lbbar_data = " << sumofvarindata(tot_ll_Lbbar_data,"wmass") << std::endl;
  }
  if (inclacc) {
    //    if (tot_dd_data->numEntries()>0.) addacceptanceweight(tot_dd_data,      wacc,false,false,acceptancedd,acceptancedd);
    //    if (tot_dd_Lbbar_data->numEntries()>0.) addacceptanceweight(tot_dd_Lbbar_data,wacc,false,false,acceptancedd,acceptancedd);
    //    if (tot_ll_data->numEntries()>0.) addacceptanceweight(tot_ll_data,      wacc,false,false,acceptancell,acceptancell);
    //    if (tot_ll_Lbbar_data->numEntries()>0.) addacceptanceweight(tot_ll_Lbbar_data,wacc,false,false,acceptancell,acceptancell);
    std::cout << "Sum of wacc in tot_dd_data        = " << sumofvarindata(tot_dd_data,"wacc") << std::endl;
    std::cout << "Sum of wacc in tot_dd_Lbbar_data  = " << sumofvarindata(tot_dd_Lbbar_data,"wacc") << std::endl;
    std::cout << "Sum of wacc in tot_ll_data        = " << sumofvarindata(tot_ll_data,"wacc") << std::endl;
    std::cout << "Sum of wacc in tot_ll_Lbbar_data  = " << sumofvarindata(tot_ll_Lbbar_data,"wacc") << std::endl;
  }
  if (inclacc and inclbkg) {
    if (tot_dd_data->numEntries()>0.)       multiplyweights(tot_dd_data,      wtot,wacc.GetName(),wmass.GetName(), normalizetowmass);
    if (tot_dd_Lbbar_data->numEntries()>0.) multiplyweights(tot_dd_Lbbar_data,wtot,wacc.GetName(),wmass.GetName(), normalizetowmass);
    if (tot_ll_data->numEntries()>0.)       multiplyweights(tot_ll_data,      wtot,wacc.GetName(),wmass.GetName(), normalizetowmass);
    if (tot_ll_Lbbar_data->numEntries()>0.) multiplyweights(tot_ll_Lbbar_data,wtot,wacc.GetName(),wmass.GetName(), normalizetowmass);
    std::cout << "Sum of wtot in tot_dd_data        = " << sumofvarindata(tot_dd_data,"wtot") << std::endl;
    std::cout << "Sum of wtot in tot_dd_Lbbar_data  = " << sumofvarindata(tot_dd_Lbbar_data,"wtot") << std::endl;
    std::cout << "Sum of wtot in tot_ll_data        = " << sumofvarindata(tot_ll_data,"wtot") << std::endl;
    std::cout << "Sum of wtot in tot_ll_Lbbar_data  = " << sumofvarindata(tot_ll_Lbbar_data,"wtot") << std::endl;
  }
  
  RooDataSet simdata2("simdata2", "simdata2", RooArgList(costheta,costheta1,costheta2,wtot), RooFit::Index(cat),RooFit::Import(mapping));
  simdata2.Print("v");
  //    simdata2.write("simdata2.txt");
  
  tot_data=&simdata2;
  
  //    RooDataSet data_sig("data_sig","data_sig",tot_data,*tot_data->get());//,0,"nSig_sw") ;
  if (debug) {
    tot_data->Print("v");
    mygetchar;
  }
  RooDataSet* tot_data_w=NULL;
  if (inclbkg  and inclacc)  tot_data_w = new RooDataSet("dataw","dataw",tot_data,*tot_data->get(),0,wtot.GetName());
  if (inclbkg  and !inclacc) tot_data_w = new RooDataSet("dataw","dataw",tot_data,*tot_data->get(),0,wmass.GetName());
  if (!inclbkg and inclacc)  tot_data_w = new RooDataSet("dataw","dataw",tot_data,*tot_data->get(),0,wacc.GetName());
  if (!inclbkg and !inclacc) tot_data_w = tot_data;
  
  if (debug) {
    tot_data_w->Print("v");
    mygetchar;
  }

  RooCmdArg minimizer=RooCmdArg::none();
  //      RooCmdArg minimizer=RooFit::Minimizer("Minuit2","minimize");
  
  for (unsigned j=0; j<angfitv.size();j++) {
    angfitv[j].fr = angfitv[j].model->fitTo(*tot_data_w, RooFit::Save(true), angfitv[j].minosarg, angfitv[j].sumw2arg, minimizer, angfitconstrains);
    //    if (doplot and j==0 and iSample<10) {
    //      plot(c,RooArgList(costheta,costheta1,costheta2), ang_model, tot_data_w, "toy-" + istr(iSample,"07") + "-fr",false,NULL,"",&cat);
    //    }
  }
  RooFitResult* fr_myfitto = myfitto(tot_data_w,&angsimpdf, angfitconstrains,-1);
  if (doplot and iSample<10) plot(c,RooArgList(costheta,costheta1,costheta2), ang_model, tot_data_w, "toy-" + istr(iSample,"07") + "-fr",false,NULL,"",&cat);

  TFile outputfile(TString("toy-" + istr(iSample,"07") + "-result.root"), "recreate");
  
  if (fr_mass) {
    std::cout << "Mass fit : status=" << fr_mass->status() << " ; covqual=" << fr_mass->covQual() << std::endl; 
    fr_mass->SetName("fr_massfit");
    fr_mass->Print("v");
    fr_mass->Write();
  }
  for (unsigned j=0; j<angfitv.size();j++) {
    std::cout << angfitv[j].fitname << " : status=" << angfitv[j].fr->status()  << " ; covqual=" << angfitv[j].fr->covQual() << std::endl;   
    angfitv[j].fr->SetName(angfitv[j].fitname);
    angfitv[j].fr->Print("v");
    angfitv[j].fr->Write();
  }

  fr_myfitto->SetName("fr_angfit_myfitto");
  fr_myfitto->Print("v");
  fr_myfitto->Write();

  outputfile.Close();

  // Extract covariance and correlation matrix as TMatrixDSym
  //   const TMatrixDSym& cor = fr->correlationMatrix() ;
  //   const TMatrixDSym& cov = fr->covarianceMatrix() ;
  
  //   // Print correlation, covariance matrix
  //   cout << "correlation matrix" << endl ;
  //   cor.Print() ;
  //   cout << "covariance matrix" << endl ;
  //   cov.Print();
  //   fr->SaveAs(TString("toy_" + istr(iSample,"04") + "_" + istr(nVariables) + "variables_result.root"));
  
  //   return 0;
  
  // TCanvas c("c","c",100,100);
  
  // RooPlot* plot = costheta.frame(50);
  // testw->plotOn(plot);
  // w3.plotOn(plot);
  // plot->Draw();
  // c.SaveAs("toy.eps");

}



// void myfitto(RooDataSet* tot_data, RooAbsPdf* sig_ang_model)
// {
  
//   RooNLLVar nll("nll","-log(L)",*sig_ang_model,*tot_data);//,NumCPU(nCPU));
//   //  nll.applyWeightSquared(true);
//   nll.Print("v");
//   //    getchar();
  
//   //      RooMinuit m1(nll);
//   RooMinimizer m1(nll) ;
//   //  m1.setVerbose(kFALSE) ;
//   //      m1.setProfile(1);
//   //  m1.setPrintLevel(3);
//   //  m1.setEps(1e-12);
//   //  m1.setStrategy(1);
//   //    m1.save()->Print("v");
//   //      m1.simplex();
//   //    getchar();

//   //  m1.hesse();
//   //  m1.migrad();
//   m1.minimize("Minuit","migrad") ;
//   m1.hesse();
//   //      m1.migrad();
//   //      m1.hesse();
//   //    RooFitResult* r1 = m1.save();
//   //     if(false){
//   //       m1.hesse();
//   //       if (m1.save()->covQual() != 3) {
//   // 	m1.migrad();
//   // 	m1.hesse();
//   //       }
//   //     }
  
//   RooFitResult* fr = m1.save();
//   std::cout << "First fit: status=" << fr->status() << " ; covqual=" << fr->covQual() << std::endl; 
//   if (fr->status()==0 and fr->covQual()==3) {
    
//     // Make list of RooNLLVar components of FCN
//     list<RooNLLVar*> nllComponents ;
//     RooArgSet* comps = nll.getComponents() ;
//     RooAbsArg* arg ;
//     TIterator* citer = comps->createIterator() ;
//     while((arg=(RooAbsArg*)citer->Next())) {
//       RooNLLVar* nllComp = dynamic_cast<RooNLLVar*>(arg) ;
//       if (nllComp) {
// 	nllComponents.push_back(nllComp) ;
//       }
//     }
//     delete citer ;
//     delete comps ;  
    
//     // Calculated corrected errors for weighted likelihood fits
//     RooFitResult* rw = m1.save() ;
//     for (list<RooNLLVar*>::iterator iter1=nllComponents.begin() ; iter1!=nllComponents.end() ; iter1++) {
//       (*iter1)->applyWeightSquared(kTRUE) ;
//     }
//     std::cout << "RooAbsPdf::fitTo Calculating sum-of-weights-squared correction matrix for covariance matrix" << endl ;
//     m1.hesse() ;
//     RooFitResult* rw2 = m1.save() ;
//     for (list<RooNLLVar*>::iterator iter2=nllComponents.begin() ; iter2!=nllComponents.end() ; iter2++) {
//       (*iter2)->applyWeightSquared(kFALSE) ;
//     }
    
//     // Apply correction matrix
//     const TMatrixDSym& V = rw->covarianceMatrix() ;
//     TMatrixDSym  C = rw2->covarianceMatrix() ;
    
//     // Invert C
//     Double_t det(0) ;
//     C.Invert(&det) ;
//     if (det==0) {
//       std::cerr << "RooAbsPdf::fitTo ERROR: Cannot apply sum-of-weights correction to covariance matrix: correction matrix calculated with weight-squared is singular" <<endl ;
//     } else {
      
// 	  // Calculate corrected covariance matrix = V C-1 V
//       TMatrixD VCV(V,TMatrixD::kMult,TMatrixD(C,TMatrixD::kMult,V)) ; 
	  
//       // Make matrix explicitly symmetric
//       Int_t n = VCV.GetNrows() ;
//       TMatrixDSym VCVsym(n) ;
//       for (Int_t i=0 ; i<n ; i++) {
// 	for (Int_t j=i ; j<n ; j++) {
// 	  if (i==j) {
// 	    VCVsym(i,j) = VCV(i,j) ;
// 	  }
// 	  if (i!=j) {
// 	    Double_t deltaRel = (VCV(i,j)-VCV(j,i))/sqrt(VCV(i,i)*VCV(j,j)) ;
// 	    if (fabs(deltaRel)>1e-3) {
// 	      std::cout << "RooAbsPdf::fitTo WARNING: Corrected covariance matrix is not (completely) symmetric: V[" << i << "," << j << "] = " 
// 			<< VCV(i,j) << " V[" << j << "," << i << "] = " << VCV(j,i) << " explicitly restoring symmetry by inserting average value" << endl ;
// 		}
// 	    VCVsym(i,j) = (VCV(i,j)+VCV(j,i))/2 ;
// 	  }
// 	}
//       }
      
//       // Propagate corrected errors to parameters objects
//       //      m1.applyCovarianceMatrix(VCVsym) ;
//       const RooArgList* _floatParamList = &rw2->floatParsFinal();//new RooArgList(sig_ang_model->getDependents());
//       Int_t _nDim=_floatParamList->getSize();
//       for (Int_t i=0 ; i<_nDim ; i++) {
// 	// Skip fixed parameters
// 	if (_floatParamList->at(i)->isConstant()) {
// 	  continue ;
// 	}
// 	//SetPdfParamErr(i, sqrt(V(i,i))) ;
// 	RooRealVar* myvar = ((RooRealVar*)_floatParamList->at(i));
// 	myvar->setError(sqrt(VCVsym(i,i)));
// 	std::cout << myvar->GetName() << " = "  << myvar->getVal() << " +/- " << myvar->getError() << std::endl;
//       }
      

//       std::cout << "new cov matrix" << std::endl;
//       VCVsym.Print("v");
//     }

//     std::cout << "rw status=" << rw->status() << " ; covqual=" << rw->covQual() << std::endl; 
//     rw->Print("v");
//     delete rw ;
//     std::cout << "rw2 status=" << rw2->status() << " ; covqual=" << rw2->covQual() << std::endl; 
//     rw2->Print("v");

//     rw2->covarianceMatrix().Print("v");
    
//     delete rw2 ;
//   }
  
  
//   //  fr = m1.save();

  
//   //    r1->Print("v") ;
//   //  fr->Print("v") ;
  
// }

