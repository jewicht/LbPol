//g++ -o toy toy.cxx `root-config --lbs --cflags` -lRooFit

#include "TCanvas.h"
#include "TRandom3.h"
#include "TMatrixD.h"

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

#include "RooLbtoJpsiL0PDF3.h"
#include "calcacceptanceclass.h"
#include "rootlogon.h"

#include "functions.h"
#include "functions-roofit.h"
#include "ConfigFile.h"
#include "weighting.h"
#include "createRooLegendre3.h"

// void removecolumns(RooDataSet*& data, const RooArgSet& vars)
// {
//   RooDataSet* newdata = (RooDataSet*)data->reduce(vars);
// }

class angfitstr {
public:
  TString fitname;
  RooAbsPdf* model;
  RooAbsPdf* model4mc;
  RooCmdArg sumw2arg;
  RooCmdArg minosarg;
  RooFitResult* fr;
  RooMCStudy* mcstudy;
  RooArgList* vars;
};


RooDataSet* generatetoy(const TString& name, const TString& title, RooAbsPdf& model, const RooArgSet& vars, const int nevents)
{
  if (nevents==0) return new RooDataSet(name,title,vars);
  RooDataSet* newdata=model.generate(vars,nevents);
  newdata->SetName(name);
  newdata->SetTitle(title);
  return newdata;
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

void myfitto(RooDataSet* tot_data, RooAbsPdf* sig_ang_model)
{
  
  RooNLLVar nll("nll","-log(L)",*sig_ang_model,*tot_data);//,NumCPU(nCPU));
  //  nll.applyWeightSquared(true);
  nll.Print("v");
  //    getchar();
  
  //      RooMinuit m1(nll);
  RooMinimizer m1(nll) ;
  //  m1.setVerbose(kFALSE) ;
  //      m1.setProfile(1);
  //  m1.setPrintLevel(3);
  //  m1.setEps(1e-12);
  //  m1.setStrategy(1);
  //    m1.save()->Print("v");
  //      m1.simplex();
  //    getchar();

  //  m1.hesse();
  //  m1.migrad();
  m1.minimize("Minuit","migrad") ;
  m1.hesse();
  //      m1.migrad();
  //      m1.hesse();
  //    RooFitResult* r1 = m1.save();
  //     if(false){
  //       m1.hesse();
  //       if (m1.save()->covQual() != 3) {
  // 	m1.migrad();
  // 	m1.hesse();
  //       }
  //     }
  
  RooFitResult* fr = m1.save();
  std::cout << "First fit: status=" << fr->status() << " ; covqual=" << fr->covQual() << std::endl; 
  if (fr->status()==0 and fr->covQual()==3) {
    
    // Make list of RooNLLVar components of FCN
    std::list<RooNLLVar*> nllComponents ;
    RooArgSet* comps = nll.getComponents() ;
    RooAbsArg* arg ;
    TIterator* citer = comps->createIterator() ;
    while((arg=(RooAbsArg*)citer->Next())) {
      RooNLLVar* nllComp = dynamic_cast<RooNLLVar*>(arg) ;
      if (nllComp) {
	nllComponents.push_back(nllComp) ;
      }
    }
    delete citer ;
    delete comps ;  
    
    // Calculated corrected errors for weighted likelihood fits
    RooFitResult* rw = m1.save() ;
    for (std::list<RooNLLVar*>::iterator iter1=nllComponents.begin() ; iter1!=nllComponents.end() ; iter1++) {
      (*iter1)->applyWeightSquared(kTRUE) ;
    }
    std::cout << "RooAbsPdf::fitTo Calculating sum-of-weights-squared correction matrix for covariance matrix" << std::endl ;
    m1.hesse() ;
    RooFitResult* rw2 = m1.save() ;
    for (std::list<RooNLLVar*>::iterator iter2=nllComponents.begin() ; iter2!=nllComponents.end() ; iter2++) {
      (*iter2)->applyWeightSquared(kFALSE) ;
    }
    
    // Apply correction matrix
    const TMatrixDSym& V = rw->covarianceMatrix() ;
    TMatrixDSym  C = rw2->covarianceMatrix() ;
    
    // Invert C
    Double_t det(0) ;
    C.Invert(&det) ;
    if (det==0) {
      std::cerr << "RooAbsPdf::fitTo ERROR: Cannot apply sum-of-weights correction to covariance matrix: correction matrix calculated with weight-squared is singular" << std::endl ;
    } else {
      
	  // Calculate corrected covariance matrix = V C-1 V
      TMatrixD VCV(V,TMatrixD::kMult,TMatrixD(C,TMatrixD::kMult,V)) ; 
	  
      // Make matrix explicitly symmetric
      Int_t n = VCV.GetNrows() ;
      TMatrixDSym VCVsym(n) ;
      for (Int_t i=0 ; i<n ; i++) {
	for (Int_t j=i ; j<n ; j++) {
	  if (i==j) {
	    VCVsym(i,j) = VCV(i,j) ;
	  }
	  if (i!=j) {
	    Double_t deltaRel = (VCV(i,j)-VCV(j,i))/sqrt(VCV(i,i)*VCV(j,j)) ;
	    if (fabs(deltaRel)>1e-3) {
	      std::cout << "RooAbsPdf::fitTo WARNING: Corrected covariance matrix is not (completely) symmetric: V[" << i << "," << j << "] = " 
			<< VCV(i,j) << " V[" << j << "," << i << "] = " << VCV(j,i) << " explicitly restoring symmetry by inserting average value" << std::endl ;
		}
	    VCVsym(i,j) = (VCV(i,j)+VCV(j,i))/2 ;
	  }
	}
      }
      
      // Propagate corrected errors to parameters objects
      //      m1.applyCovarianceMatrix(VCVsym) ;
      const RooArgList* _floatParamList = &rw2->floatParsFinal();//new RooArgList(sig_ang_model->getDependents());
      Int_t _nDim=_floatParamList->getSize();
      for (Int_t i=0 ; i<_nDim ; i++) {
	// Skip fixed parameters
	if (_floatParamList->at(i)->isConstant()) {
	  continue ;
	}
	//SetPdfParamErr(i, sqrt(V(i,i))) ;
	RooRealVar* myvar = ((RooRealVar*)_floatParamList->at(i));
	myvar->setError(sqrt(VCVsym(i,i)));
	std::cout << myvar->GetName() << " = "  << myvar->getVal() << " +/- " << myvar->getError() << std::endl;
      }
      

      std::cout << "new cov matrix" << std::endl;
      VCVsym.Print("v");
    }

    std::cout << "rw status=" << rw->status() << " ; covqual=" << rw->covQual() << std::endl; 
    rw->Print("v");
    delete rw ;
    std::cout << "rw2 status=" << rw2->status() << " ; covqual=" << rw2->covQual() << std::endl; 
    rw2->Print("v");

    rw2->covarianceMatrix().Print("v");
    
    delete rw2 ;
  }
  
  
  //  fr = m1.save();

  
  //    r1->Print("v") ;
  //  fr->Print("v") ;
  
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
 
  RooRealVar costheta("costheta",   "cos(#theta)",     -1., 1.);//RooArgList(theta));
  RooRealVar costheta1("costheta1", "cos(#theta_{1})", -1., 1.);// RooArgList(theta1));
  RooRealVar costheta2("costheta2", "cos(#theta_{2})", -1., 1.);// RooArgList(theta2));

  costheta.setBins(50);
  costheta1.setBins(50);
  costheta2.setBins(50);

  RooRealVar piminus_TRACK_Type("piminus_TRACK_Type","piminus_TRACK_Type",0.);

  RooRealVar P_b("P_b","P_{b}",0.40, -1.5, 1.5);//0.50
  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}",0.642, 0., 1.);
  //  alpha_lambda.setError(0.013);
  //  alpha_lambda.setConstant();
  RooRealVar alpha_lambda_Lbbar("alpha_lambda_Lbbar","#alpha_{#lambda}",0.71,0.,1.);//0.71

 RooGaussian alpha_lambda_constrain("alpha_lambda_constrain",
				    "alpha_lambda_constrain",
				    alpha_lambda,
				    RooFit::RooConst(0.642),
				    RooFit::RooConst(0.013));
 RooGaussian alpha_lambda_Lbbar_constrain("alpha_lambda_Lbbar_constrain",
					  "alpha_lambda_Lbbar_constrain",
					  alpha_lambda_Lbbar,
					  RooFit::RooConst(0.71),
					  RooFit::RooConst(0.08));

  RooRealVar alpha_b("alpha_b","#alpha_{b}",-0.45767,-1.1,1.1);//0.457
  RooRealVar r_0("r_0","r_{0}",0.251733,-0.3,2.3);//0.5
  RooRealVar r_1("r_1","r_{1}",0.116484,-1.3,1.3);//0.1

  RooRealVar ap2("ap2","|a+|^{2}",0.25,-0.1,1.1);
  RooRealVar am2("am2","|a-|^{2}",0.25,-0.1,1.1);
  RooRealVar bp2("bp2","|b+|^{2}",0.25,-0.1,1.1);
  RooFormulaVar bm2("bm2","1.0-ap2-am2-bp2",RooArgList(ap2,am2,bp2));
  RooFormulaVar alpha_b_rfv("alpha_b_rfv","ap2-am2+bp2-bm2",RooArgList(ap2,am2,bp2,bm2));
  RooFormulaVar r_0_rfv("r_0_rfv","ap2+am2",RooArgList(ap2,am2));
  RooFormulaVar r_1_rfv("r_1_rfv","ap2-am2",RooArgList(ap2,am2));
 
  RooRealVar mass("mass","M(#LambdaJ/#psi)",5500.,5750.,"MeV/c^{2}");  
  mass.setRange("signalregion", 5600.,5640.);

  RooRealVar sig_dd_mean1("sig_dd_mean1","sig_dd_mean1",5620.,5500.,5700.);
  RooRealVar sig_ll_mean1("sig_ll_mean1","sig_ll_mean1",5620.,5500.,5700.);

  RooRealVar sig_dd_frac1("sig_dd_frac1","sig_dd_frac1", 0.5, -0.1, 1.1);
  RooRealVar sig_ll_frac1("sig_ll_frac1","sig_ll_frac1", 0.5, -0.1, 1.1);

  RooRealVar sig_dd_sigma1("sig_dd_sigma1","sig_dd_sigma1",7.5,0.,100.);
  RooRealVar sig_ll_sigma1("sig_ll_sigma1","sig_ll_sigma1",7.5,0.,100.);


  RooRealVar sig_dd_alpha1("sig_dd_alpha1","sig_dd_alpha1", 1.,0.,10.);
  RooRealVar sig_dd_alpha2("sig_dd_alpha2","sig_dd_alpha2",-1.,-10.,0.);

  RooRealVar sig_dd_n1("sig_dd_n1","sig_dd_n1",10.,0.,200.);
  RooRealVar sig_dd_n2("sig_dd_n2","sig_dd_n2",10.,0.,200.);

  RooRealVar sig_ll_alpha1("sig_ll_alpha1","sig_ll_alpha1", 1.,0.,10.);
  RooRealVar sig_ll_alpha2("sig_ll_alpha2","sig_ll_alpha2",-1.,-10.,0.);

  RooRealVar sig_ll_n1("sig_ll_n1","sig_ll_n1",10.,0.,200.);
  RooRealVar sig_ll_n2("sig_ll_n2","sig_ll_n2",10.,0.,200.);
  
  RooCBShape sig_dd_cbshape1("sig_dd_cbshape1", "sig_dd_cbshape1", mass, sig_dd_mean1, sig_dd_sigma1, sig_dd_alpha1, sig_dd_n1);//, RooAbsReal& _n)
  RooCBShape sig_dd_cbshape2("sig_dd_cbshape2", "sig_dd_cbshape2", mass, sig_dd_mean1, sig_dd_sigma1, sig_dd_alpha2, sig_dd_n2);//, RooAbsReal& _n)
  RooCBShape sig_ll_cbshape1("sig_ll_cbshape1", "sig_ll_cbshape1", mass, sig_ll_mean1, sig_ll_sigma1, sig_ll_alpha1, sig_ll_n1);//, RooAbsReal& _n)
  RooCBShape sig_ll_cbshape2("sig_ll_cbshape2", "sig_ll_cbshape2", mass, sig_ll_mean1, sig_ll_sigma1, sig_ll_alpha2, sig_ll_n2);//, RooAbsReal& _n)

  RooAddPdf sig_dd_mass_model("sig_dd_mass_model","sig_dd_mass_model",RooArgList(sig_dd_cbshape1,sig_dd_cbshape2),RooArgList(sig_dd_frac1));
  RooAddPdf sig_ll_mass_model("sig_ll_mass_model","sig_ll_mass_model",RooArgList(sig_ll_cbshape1,sig_ll_cbshape2),RooArgList(sig_ll_frac1));


  RooRealVar bkg_dd_c0("bkg_dd_c0","bkg_dd_c0",0.,-1.1,1.1);
  RooRealVar bkg_ll_c0("bkg_ll_c0","bkg_ll_c0",0.,-1.1,1.1);

  RooRealVar bkg_dd_Lbbar_c0("bkg_dd_Lbbar_c0","bkg_dd_Lbbar_c0",0.,-1.1,1.1);
  RooRealVar bkg_ll_Lbbar_c0("bkg_ll_Lbbar_c0","bkg_ll_Lbbar_c0",0.,-1.1,1.1);


  RooChebychev bkg_dd_mass_model("bkg_dd_mass_model","bkg_dd_mass_model",mass,RooArgList(bkg_dd_c0));
  RooChebychev bkg_ll_mass_model("bkg_ll_mass_model","bkg_ll_mass_model",mass,RooArgList(bkg_ll_c0));
  RooChebychev bkg_dd_Lbbar_mass_model("bkg_dd_Lbbar_mass_model","bkg_dd_Lbbar_mass_model",mass,RooArgList(bkg_dd_Lbbar_c0));
  RooChebychev bkg_ll_Lbbar_mass_model("bkg_ll_Lbbar_mass_model","bkg_ll_Lbbar_mass_model",mass,RooArgList(bkg_ll_Lbbar_c0));

  RooRealVar nSig_dd("nSig_dd","nSig_dd",2000.,0.,300000.);
  RooRealVar nBkg_dd("nBkg_dd","nBkg_dd",10000.,0.,1000000.);
  RooRealVar nSig_ll("nSig_ll","nSig_ll",2000.,0.,300000.);
  RooRealVar nBkg_ll("nBkg_ll","nBkg_ll",10000.,0.,1000000.);
  
  RooRealVar nSig_dd_Lbbar("nSig_dd_Lbbar","nSig_dd_Lbbar",2000.,0.,300000.);
  RooRealVar nBkg_dd_Lbbar("nBkg_dd_Lbbar","nBkg_dd_Lbbar",10000.,0.,1000000.);
  RooRealVar nSig_ll_Lbbar("nSig_ll_Lbbar","nSig_ll_Lbbar",2000.,0.,300000.);
  RooRealVar nBkg_ll_Lbbar("nBkg_ll_Lbbar","nBkg_ll_Lbbar",10000.,0.,1000000.);

  RooAddPdf tot_dd_mass_model("tot_dd_mass_model","tot_dd_mass_model",RooArgList(sig_dd_mass_model,bkg_dd_mass_model),RooArgList(nSig_dd,nBkg_dd));
  RooAddPdf tot_ll_mass_model("tot_ll_mass_model","tot_ll_mass_model",RooArgList(sig_ll_mass_model,bkg_ll_mass_model),RooArgList(nSig_ll,nBkg_ll));

  RooAddPdf tot_dd_Lbbar_mass_model("tot_dd_Lbbar_mass_model","tot_dd_Lbbar_mass_model",RooArgList(sig_dd_mass_model,bkg_dd_Lbbar_mass_model),RooArgList(nSig_dd_Lbbar,nBkg_dd_Lbbar));
  RooAddPdf tot_ll_Lbbar_mass_model("tot_ll_Lbbar_mass_model","tot_ll_Lbbar_mass_model",RooArgList(sig_ll_mass_model,bkg_ll_Lbbar_mass_model),RooArgList(nSig_ll_Lbbar,nBkg_ll_Lbbar));


  ConfigFile cf("toyconfig.txt");
  
  const int nSigEvents_dd=atoi(cf.Value ("main","nSigEvents_dd","6000").c_str());
  const int nSigEvents_dd_Lbbar=atoi(cf.Value ("main","nSigEvents_dd_Lbbar","6000").c_str());
  const int nSigEvents_ll=atoi(cf.Value ("main","nSigEvents_ll","6000").c_str());
  const int nSigEvents_ll_Lbbar=atoi(cf.Value ("main","nSigEvents_ll_Lbbar","6000").c_str());
  const int nBkgEvents_dd=atoi(cf.Value ("main","nBkgEvents_dd","6000").c_str());
  const int nBkgEvents_dd_Lbbar=atoi(cf.Value ("main","nBkgEvents_dd_Lbbar","6000").c_str());
  const int nBkgEvents_ll=atoi(cf.Value ("main","nBkgEvents_ll","6000").c_str());
  const int nBkgEvents_ll_Lbbar=atoi(cf.Value ("main","nBkgEvents_ll_Lbbar","6000").c_str());

  const bool debug=atob(cf.Value("main","debug","false"));
  const bool simplefit=atob(cf.Value("main","simplefit","false"));
  const bool doplot=atob(cf.Value("main","doplot","true"));
  //  const bool useformparam=atob(cf.Value("main","useformparam","false"));
  const bool inclacc=atob(cf.Value("main","inclacc","true"));
  const bool inclbkg=(nBkgEvents_dd+nBkgEvents_ll+nBkgEvents_dd_Lbbar+nBkgEvents_ll_Lbbar!=0);
  const bool usesbforbkg=atob(cf.Value("main","usesbforbkg","true"));

  P_b=atof(cf.Value ("amplitudes","P_b","0.2").c_str());
  alpha_b=atof(cf.Value ("amplitudes","alpha_b","-0.45767").c_str());
  r_0=atof(cf.Value ("amplitudes","r_0","0.251733").c_str());
  r_1=atof(cf.Value ("amplitudes","r_1","0.116484").c_str());

  bkg_dd_c0=atof(cf.Value("bkgshape","c0","0.").c_str());
  bkg_dd_Lbbar_c0=bkg_dd_c0.getVal();
  bkg_ll_c0=bkg_dd_c0.getVal();
  bkg_ll_Lbbar_c0=bkg_dd_c0.getVal();

  computew3 w3;
  ap2=w3.a_plus_sq(alpha_b.getVal(), r_0.getVal(), r_1.getVal());
  am2=w3.a_minus_sq(alpha_b.getVal(), r_0.getVal(), r_1.getVal());
  bp2=w3.b_plus_sq(alpha_b.getVal(), r_0.getVal(), r_1.getVal());

  const RooArgList sig_dd_vars(sig_dd_frac1,sig_dd_mean1,sig_dd_sigma1,sig_dd_alpha1,sig_dd_alpha2,sig_dd_n1,sig_dd_n2);
  const RooArgList mass_dd_pdfs(sig_dd_mass_model,bkg_dd_mass_model);
  
  const RooArgList sig_ll_vars(sig_ll_frac1,sig_ll_mean1,sig_ll_sigma1,sig_ll_alpha1,sig_ll_alpha2,sig_ll_n1,sig_ll_n2);
  const RooArgList mass_ll_pdfs(sig_ll_mass_model,bkg_ll_mass_model);

  const TString directorylist = "Tuple_JpsiL0_betas";
  readvarsfromtfile(sig_dd_vars,TString("fitmcmass-Lb2JpsiL0-") + directorylist + "-mc-dd-fr.root");
  readvarsfromtfile(sig_ll_vars,TString("fitmcmass-Lb2JpsiL0-") + directorylist + "-mc-ll-fr.root");

  RooGaussian sig_dd_frac1_constrain("sig_dd_frac1_constrain", "sig_dd_frac1_constrain", sig_dd_frac1, RooFit::RooConst(sig_dd_frac1.getVal()), RooFit::RooConst(sig_dd_frac1.getError()));
  RooGaussian sig_ll_frac1_constrain("sig_ll_frac1_constrain", "sig_ll_frac1_constrain", sig_ll_frac1, RooFit::RooConst(sig_ll_frac1.getVal()), RooFit::RooConst(sig_ll_frac1.getError()));
  RooGaussian sig_dd_alpha1_constrain("sig_dd_alpha1_constrain", "sig_dd_alpha1_constrain", sig_dd_alpha1, RooFit::RooConst(sig_dd_alpha1.getVal()), RooFit::RooConst(sig_dd_alpha1.getError()));
  RooGaussian sig_ll_alpha1_constrain("sig_ll_alpha1_constrain", "sig_ll_alpha1_constrain", sig_ll_alpha1, RooFit::RooConst(sig_ll_alpha1.getVal()), RooFit::RooConst(sig_ll_alpha1.getError()));
  RooGaussian sig_dd_alpha2_constrain("sig_dd_alpha2_constrain", "sig_dd_alpha2_constrain", sig_dd_alpha2, RooFit::RooConst(sig_dd_alpha2.getVal()), RooFit::RooConst(sig_dd_alpha2.getError()));
  RooGaussian sig_ll_alpha2_constrain("sig_ll_alpha2_constrain", "sig_ll_alpha2_constrain", sig_ll_alpha2, RooFit::RooConst(sig_ll_alpha2.getVal()), RooFit::RooConst(sig_ll_alpha2.getError()));
  RooGaussian sig_dd_n1_constrain("sig_dd_n1_constrain", "sig_dd_n1_constrain", sig_dd_n1, RooFit::RooConst(sig_dd_n1.getVal()), RooFit::RooConst(sig_dd_n1.getError()));
  RooGaussian sig_ll_n1_constrain("sig_ll_n1_constrain", "sig_ll_n1_constrain", sig_ll_n1, RooFit::RooConst(sig_ll_n1.getVal()), RooFit::RooConst(sig_ll_n1.getError()));
  RooGaussian sig_dd_n2_constrain("sig_dd_n2_constrain", "sig_dd_n2_constrain", sig_dd_n2, RooFit::RooConst(sig_dd_n2.getVal()), RooFit::RooConst(sig_dd_n2.getError()));
  RooGaussian sig_ll_n2_constrain("sig_ll_n2_constrain", "sig_ll_n2_constrain", sig_ll_n2, RooFit::RooConst(sig_ll_n2.getVal()), RooFit::RooConst(sig_ll_n2.getError()));
  
  
  RooArgSet massfitblu(sig_dd_frac1_constrain, sig_dd_alpha1_constrain, sig_dd_alpha2_constrain, sig_dd_n1_constrain, sig_dd_n2_constrain);
  massfitblu.add(RooArgSet(sig_ll_frac1_constrain, sig_ll_alpha1_constrain, sig_ll_alpha2_constrain, sig_ll_n1_constrain, sig_ll_n2_constrain));
  RooCmdArg massfitconstrains =   RooFit::ExternalConstraints(massfitblu);

  RooCmdArg angfitconstrains = RooFit::ExternalConstraints(RooArgSet(alpha_lambda_constrain,alpha_lambda_Lbbar_constrain));
  

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
    std::cout << "  usesbforbkg  = " << usesbforbkg << std::endl;
    std::cout << "  P_b          = " << P_b.getVal() << std::endl;
    //    if (useformparam) {
      std::cout << "  ap2          = " << ap2.getVal() << std::endl;
      std::cout << "  am2          = " << am2.getVal() << std::endl;
      std::cout << "  bp2          = " << bp2.getVal() << std::endl;
      std::cout << "  bm2          = " << bm2.getVal() << std::endl;
      std::cout << "  sum          = " << ap2.getVal()+am2.getVal()+bp2.getVal()+bm2.getVal() << std::endl;
      std::cout << "  alpha_b      = " << alpha_b_rfv.getVal() << std::endl;
      std::cout << "  r_0          = " << r_0_rfv.getVal() << std::endl;
      std::cout << "  r_1          = " << r_1_rfv.getVal() << std::endl;
      //    } else {
      std::cout << "  alpha_b      = " << alpha_b.getVal() << std::endl;
      std::cout << "  r_0          = " << r_0.getVal() << std::endl;
      std::cout << "  r_1          = " << r_1.getVal() << std::endl;
      //    }
    std::cout << "  sig_dd_mean         = " << sig_dd_mean1.getVal() << std::endl;
    std::cout << "  sig_dd_sigma1       = " << sig_dd_sigma1.getVal() << std::endl;
    std::cout << "  sig_dd_alpha1       = " << sig_dd_alpha1.getVal() << " +/- " << sig_dd_alpha1.getError() << std::endl;
    std::cout << "  sig_dd_alpha2       = " << sig_dd_alpha2.getVal() << " +/- " << sig_dd_alpha2.getError() << std::endl;
    std::cout << "  sig_dd_n1           = " << sig_dd_n1.getVal() << " +/- " << sig_dd_n1.getError() << std::endl;
    std::cout << "  sig_dd_n2           = " << sig_dd_n2.getVal() << " +/- " << sig_dd_n2.getError() << std::endl;
    std::cout << "  sig_dd_frac1        = " << sig_dd_frac1.getVal() << " +/- " << sig_dd_frac1.getError() << std::endl;
    std::cout << "  sig_ll_mean         = " << sig_ll_mean1.getVal() << std::endl;
    std::cout << "  sig_ll_sigma1       = " << sig_ll_sigma1.getVal() << std::endl;
    std::cout << "  sig_ll_alpha1       = " << sig_ll_alpha1.getVal() << " +/- " << sig_ll_alpha1.getError() << std::endl;
    std::cout << "  sig_ll_alpha2       = " << sig_ll_alpha2.getVal() << " +/- " << sig_ll_alpha2.getError() << std::endl;
    std::cout << "  sig_ll_n1           = " << sig_ll_n1.getVal() << " +/- " << sig_ll_n1.getError() << std::endl;
    std::cout << "  sig_ll_n2           = " << sig_ll_n2.getVal() << " +/- " << sig_ll_n2.getError() << std::endl;
    std::cout << "  sig_ll_frac1        = " << sig_ll_frac1.getVal() << " +/- " << sig_ll_frac1.getError() << std::endl;

    std::cout << "  bkg_dd_c0           = " << bkg_dd_c0.getVal() << std::endl;
    std::cout << "  bkg_dd_Lbbar_c0     = " << bkg_dd_Lbbar_c0.getVal() << std::endl;
    std::cout << "  bkg_ll_c0           = " << bkg_ll_c0.getVal() << std::endl;
    std::cout << "  bkg_ll_Lbbar_c0     = " << bkg_ll_Lbbar_c0.getVal() << std::endl;
  }


  if (nSigEvents_dd==0.) nSig_dd.setConstant();
  if (nSigEvents_dd_Lbbar==0.) nSig_dd_Lbbar.setConstant();
  if (nSigEvents_dd+nSigEvents_dd_Lbbar==0.) {
    sig_dd_alpha1.setConstant();
    sig_dd_alpha2.setConstant();
    sig_dd_frac1.setConstant();
    sig_dd_mean1.setConstant();
    sig_dd_n1.setConstant();
    sig_dd_n2.setConstant();
    sig_dd_sigma1.setConstant();
  }
  
  if (nSigEvents_ll==0.) nSig_ll.setConstant();
  if (nSigEvents_ll_Lbbar==0.) nSig_ll_Lbbar.setConstant();
  if (nSigEvents_ll+nSigEvents_ll_Lbbar==0.) {
    sig_ll_alpha1.setConstant();
    sig_ll_alpha2.setConstant();
    sig_ll_frac1.setConstant();
    sig_ll_mean1.setConstant();
    sig_ll_n1.setConstant();
    sig_ll_n2.setConstant();
    sig_ll_sigma1.setConstant();
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

  RooLbtoJpsiL0PDF3* sig_ang_model = new RooLbtoJpsiL0PDF3("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1);
  RooLbtoJpsiL0PDF3* sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF3("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda_Lbbar, r_0, r_1);
  RooLbtoJpsiL0PDF3* sig_ang_model_fp = new RooLbtoJpsiL0PDF3("sig_ang_model_fp","sig_ang_model_fp",costheta,costheta1,costheta2,P_b, alpha_b_rfv, alpha_lambda, r_0_rfv, r_1_rfv);
  RooLbtoJpsiL0PDF3* sig_ang_model_Lbbar_fp = new RooLbtoJpsiL0PDF3("sig_ang_model_Lbbar_fp","sig_ang_model_Lbbar_fp",costheta,costheta1,costheta2,P_b, alpha_b_rfv, alpha_lambda_Lbbar, r_0_rfv, r_1_rfv);


  const bool usebdt=false;
  const bool userwmc=true;
  const bool usetmvarectcut=true;
  const bool doBd2JpsiKS0=false;
  TString acceptancefile="fitacceptance3";
  if (userwmc) {
    acceptancefile+="-mcrw";
  } else {
    acceptancefile+="-mcnonrw";
  }
  if (usebdt) acceptancefile+="-wbdt";
  if (usetmvarectcut) acceptancefile+="-wtmvarectcut";
  if (doBd2JpsiKS0) {
    acceptancefile+="-B0-cc0c1c2";
  } else {
    acceptancefile+="-Lb-cc0c1c2";
  }
  

  RooRealVar *acccoeffll[ordermax3i][ordermax3j][ordermax3k];  
  RooRealVar *acccoeffdd[ordermax3i][ordermax3j][ordermax3k];  
//  RooAbsPdf *accll, *accdd;
//  createRooLegendre3v2(costheta, costheta1, costheta2, "c0c1c2ll", accll, acccoeffll, true, 7,5,7,8,2, acceptancefile + "-ll.txt");
//  createRooLegendre3v2(costheta, costheta1, costheta2, "c0c1c2dd", accdd, acccoeffdd, true, 7,5,7,8,2, acceptancefile + "-dd.txt");

   for (int i=0;i<ordermax3i;i++) {
     for (int j=0;j<ordermax3j;j++) {
       for (int k=0;k<ordermax3k;k++) {
	 acccoeffll[i][j][k]= new RooRealVar("acccoeffll" + istr(i) + istr(j) + istr(k), "acccoeffll" + istr(i) + istr(j) + istr(k), 0.);
	 acccoeffdd[i][j][k]= new RooRealVar("acccoeffdd" + istr(i) + istr(j) + istr(k), "acccoeffdd" + istr(i) + istr(j) + istr(k), 0.);
	   //	 if (acccoeffll[i][j][k]) acccoeffll[i][j][k]->setConstant();
	   //	 if (acccoeffdd[i][j][k]) acccoeffdd[i][j][k]->setConstant();
       }
     }
   }

  RooAbsPdf *tot_ang_model_ll, *tot_ang_model_Lbbar_ll, *tot_ang_model_dd,  *tot_ang_model_Lbbar_dd;
  RooAbsPdf *tot_ang_model_fp_ll, *tot_ang_model_Lbbar_fp_ll, *tot_ang_model_fp_dd,  *tot_ang_model_Lbbar_fp_dd;

if (inclacc){
  readcoeff3_rrv(acceptancefile + "-ll.txt", acccoeffll);
  readcoeff3_rrv(acceptancefile + "-dd.txt", acccoeffdd);

  createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b,alpha_lambda, r_0, r_1, "tot_ang_model_ll", tot_ang_model_ll, acccoeffll, false);
  createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b,alpha_lambda_Lbbar, r_0, r_1, "tot_ang_model_Lbbar_ll", tot_ang_model_Lbbar_ll, acccoeffll, false);
  createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b,alpha_lambda, r_0, r_1, "tot_ang_model_dd", tot_ang_model_dd, acccoeffdd, false);
  createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b,alpha_lambda_Lbbar, r_0, r_1, "tot_ang_model_Lbbar_dd", tot_ang_model_Lbbar_dd, acccoeffdd, false);

  createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b_rfv,alpha_lambda, r_0_rfv, r_1_rfv, "tot_ang_model_fp_ll", tot_ang_model_fp_ll, acccoeffll, false);
  createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b_rfv,alpha_lambda_Lbbar, r_0_rfv, r_1_rfv, "tot_ang_model_Lbbar_fp_ll", tot_ang_model_Lbbar_fp_ll, acccoeffll, false);
  createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b_rfv,alpha_lambda, r_0_rfv, r_1_rfv, "tot_ang_model_fp_dd", tot_ang_model_fp_dd, acccoeffdd, false);
  createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b_rfv,alpha_lambda_Lbbar, r_0_rfv, r_1_rfv, "tot_ang_model_Lbbar_fp_dd", tot_ang_model_Lbbar_fp_dd, acccoeffdd, false);


 } else {
  tot_ang_model_ll = sig_ang_model;
  tot_ang_model_Lbbar_ll = sig_ang_model_Lbbar;
  tot_ang_model_dd =  sig_ang_model;
  tot_ang_model_Lbbar_dd =  sig_ang_model_Lbbar;

  tot_ang_model_fp_ll = sig_ang_model_fp;
  tot_ang_model_Lbbar_fp_ll = sig_ang_model_Lbbar_fp;
  tot_ang_model_fp_dd = sig_ang_model_fp;
  tot_ang_model_Lbbar_fp_dd = sig_ang_model_Lbbar_fp;


}
  RooProdPdf bkg_dd_model("bkg_dd_model","bkg_dd_model",RooArgList(bkg_dd_mass_model,bkg_ang_model));
  RooProdPdf bkg_ll_model("bkg_ll_model","bkg_ll_model",RooArgList(bkg_ll_mass_model,bkg_ang_model));
  RooProdPdf bkg_dd_Lbbar_model("bkg_dd_Lbbar_model","bkg_dd_Lbbar_model",RooArgList(bkg_dd_Lbbar_mass_model,bkg_ang_model));
  RooProdPdf bkg_ll_Lbbar_model("bkg_ll_Lbbar_model","bkg_ll_Lbbar_model",RooArgList(bkg_ll_Lbbar_mass_model,bkg_ang_model));


  RooCategory cat("cat","cat");  
  std::vector<TString> catlist;
  catlist.push_back("llLb"   );
  catlist.push_back("llLbbar");
  catlist.push_back("ddLb"   );
  catlist.push_back("ddLbbar");
  for (unsigned int i=0;i<catlist.size();i++) cat.defineType(catlist[i],i+1);
  

  std::map<std::string,RooAbsPdf*> masspdfmapping;
  masspdfmapping[ "ddLb"    ] = &tot_dd_mass_model;
  masspdfmapping[ "ddLbbar" ] = &tot_dd_Lbbar_mass_model;
  masspdfmapping[ "llLb"    ] = &tot_ll_mass_model;
  masspdfmapping[ "llLbbar" ] = &tot_ll_Lbbar_mass_model;
  
  
  std::map<std::string,RooAbsPdf*> angpdfmapping;
  angpdfmapping[ "ddLb"    ] = tot_ang_model_dd;
  angpdfmapping[ "ddLbbar" ] = tot_ang_model_Lbbar_dd;
  angpdfmapping[ "llLb"    ] = tot_ang_model_ll;
  angpdfmapping[ "llLbbar" ] = tot_ang_model_Lbbar_ll;

  std::map<std::string,RooAbsPdf*> angpdfmapping_fp;
  angpdfmapping_fp[ "ddLb"    ] = tot_ang_model_fp_dd;
  angpdfmapping_fp[ "ddLbbar" ] = tot_ang_model_Lbbar_fp_dd;
  angpdfmapping_fp[ "llLb"    ] = tot_ang_model_fp_ll;
  angpdfmapping_fp[ "llLbbar" ] = tot_ang_model_Lbbar_fp_ll;
  
  RooSimultaneous masssimpdf("masssimpdf","masssimpdf",masspdfmapping,cat);
  RooSimultaneous angsimpdf("angsimpdf","angsimpdf",angpdfmapping,cat);
  RooSimultaneous angsimpdf_fp("angsimpdf_fp","angsimpdf_fp",angpdfmapping_fp,cat);

  RooAbsPdf&  tot_mass_model = masssimpdf;
  RooAbsPdf&  ang_model      = angsimpdf;
  //  RooAbsPdf&  ang_model_fp   = angsimpdf_fp;


  RooRealVar wmass("wmass","wmass",0.);
  RooRealVar wacc("wacc","wacc",0.);
  RooRealVar wtot("wtot","wtot",0.);

  RooArgList vars(P_b,alpha_b,r_0,r_1);
  RooArgList varsfp(P_b,ap2,am2,bp2);
  
  std::vector<angfitstr> angfitv;
  angfitstr angfit;

  angfit.fitname="fr_angfit_woutsumw2";
  angfit.model=&angsimpdf;
  angfit.model4mc=sig_ang_model;
  angfit.vars=&vars;
  angfit.sumw2arg=RooFit::SumW2Error(false);
  angfit.minosarg=RooFit::Minos(true);
  angfitv.push_back(angfit);

  angfit.fitname="fr_angfit_woutsumw2_woutminos";
  angfit.sumw2arg=RooFit::SumW2Error(false);
  angfit.minosarg=RooFit::Minos(false);
  angfitv.push_back(angfit);

  angfit.fitname="fr_angfit_wsumw2";
  angfit.sumw2arg=RooFit::SumW2Error(inclbkg or inclacc);
  angfit.minosarg=RooFit::Minos(false);
  angfitv.push_back(angfit);


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
    RooArgList fitvars(costheta,costheta1,costheta2);

    for (unsigned int i=0; i<angfitv.size(); i++) angfitv[i].mcstudy = new RooMCStudy(*angfitv[i].model4mc, fitvars);

    TH1F histo1("histo1","histo1",50,-5.,5.);
    TH1F histo2("histo2","histo2",50,-5.,5.);
    TH1F histo3("histo3","histo3",50,-5.,5.);
    
    for (int i=0; i<iSample;i++) {
      const TString filename="toy_" + istr(i,"04") + "_result.root";
      TFile* frfile = new TFile(filename);
      if (frfile and !frfile->IsZombie()) {

	for (unsigned int j=0; j<angfitv.size(); j++) {
	  RooFitResult* fr = (RooFitResult*)frfile->Get(angfitv[j].fitname);
	  if (fr) angfitv[j].mcstudy->addFitResult(*fr);
	}

	//        RooFitResult* fr_wsumw2 = (RooFitResult*)frfile->Get(frname_wsumw2);
	//        RooFitResult* fr_woutsumw2 = (RooFitResult*)frfile->Get(frname_woutsumw2);
	//        RooFitResult* fr_woutsumw2_woutminos = (RooFitResult*)frfile->Get(frname_woutsumw2_woutminos);
	//	if (fr_wsumw2) mgr_wsumw2.addFitResult(*fr_wsumw2);
	//	if (fr_woutsumw2) mgr_woutsumw2.addFitResult(*fr_woutsumw2);
	//	if (fr_woutsumw2_woutminos) mgr_woutsumw2_woutminos.addFitResult(*fr_woutsumw2_woutminos);
	//	std::cout << fr->GetName() << std::endl;
	// if (fr) {
	//   if (fr->status()==0 and fr->covQual()>-3)  {
	//     std::cerr << filename << ": CONVERGED:    status=" << fr->status() << " ; covQual=" << fr->covQual() << std::endl;
	//     mgr.addFitResult(*fr);
	//     mgr.addFitResult(*fr);
	//     if (useformparam) {
	//       histo1.Fill( ( alpha_b_rfv.getVal(fr->floatParsFinal()) - alpha_b.getVal()  ) / alpha_b_rfv.getPropagatedError(*fr) );
	//       histo2.Fill( ( r_0_rfv.getVal(fr->floatParsFinal()) - r_0.getVal()  ) / r_0_rfv.getPropagatedError(*fr) );
	//       histo3.Fill( ( r_1_rfv.getVal(fr->floatParsFinal()) - r_1.getVal()  ) / r_1_rfv.getPropagatedError(*fr) );
	//     }
	//   } else {
	//     std::cerr << filename << ": NOT CONVERGED: status=" << fr->status() << " ; covQual=" << fr->covQual() << std::endl;
	//   }
      } else {
	std::cerr << filename << ": fit result not found" << std::endl;
      }
      frfile->Close();
    }
    //    mygetchar;
    //    TCanvas c("c","c",100,100);
    //    plot_mcstudy(c,mgr_wsumw2,"mcstudy-wsumw2",vars,nVariables);
    //  plot_mcstudy(c,mgr_woutsumw2,"mcstudy-woutsumw2",vars,nVariables);
    //  plot_mcstudy(c,mgr_woutsumw2_woutminos,"mcstudy-woutsumw2-woutminos",vars,nVariables);
    for (unsigned int j=0; j<angfitv.size(); j++) {
      plot_mcstudy(c,*angfitv[j].mcstudy, "mcstudy-" + angfitv[j].fitname, *angfitv[j].vars);
    }

    // if (useformparam) { 
    //   histo1.SetStats(1);
    //   histo2.SetStats(1);
    //   histo3.SetStats(1);
    //   gStyle->SetOptFit(1);
    //   histo1.GetXaxis()->SetTitle(alpha_b.GetTitle());
    //   histo1.Fit("gaus");
    //   histo1.Draw("E");
    //   c.SaveAs("histo-alpha_b.eps");
    //   histo2.GetXaxis()->SetTitle(r_0.GetTitle());
    //   histo2.Fit("gaus");
    //   histo2.Draw("E");
    //   c.SaveAs("histo-r_0.eps");
    //   histo3.GetXaxis()->SetTitle(r_1.GetTitle());
    //   histo3.Fit("gaus");
    //   histo3.Draw("E");
    //   c.SaveAs("histo-r_1.eps");
    // }
    return 0;
    

  }




  RooDataSet* sig_dd_ang_data=NULL;
  RooDataSet* sig_dd_Lbbar_ang_data=NULL;
  RooDataSet* sig_ll_ang_data=NULL;
  RooDataSet* sig_ll_Lbbar_ang_data=NULL;

  sig_dd_ang_data = generatetoy("sig_dd_ang_data","sig_dd_ang_data",*tot_ang_model_dd,RooArgSet(costheta,costheta1,costheta2), nSigEvents_dd);
  sig_dd_Lbbar_ang_data = generatetoy("sig_dd_Lbbar_ang_data","sig_dd_Lbbar_ang_data",*tot_ang_model_Lbbar_dd,RooArgSet(costheta,costheta1,costheta2), nSigEvents_dd_Lbbar);
  sig_ll_ang_data = generatetoy("sig_ll_ang_data","sig_ll_ang_data",*tot_ang_model_ll,RooArgSet(costheta,costheta1,costheta2), nSigEvents_ll);
  sig_ll_Lbbar_ang_data = generatetoy("sig_ll_Lbbar_ang_data","sig_ll_Lbbar_ang_data",*tot_ang_model_Lbbar_ll,RooArgSet(costheta,costheta1,costheta2), nSigEvents_ll_Lbbar);
  
  // sig_dd_ang_data->Print("v");
  // deletecolumn( sig_dd_ang_data, costheta);
  // sig_dd_ang_data->Print("v");
  // mygetchar;

  RooDataSet* sig_dd_mass_data       = generatetoy("sig_dd_mass_data",      "sig_dd_mass_data",      sig_dd_mass_model,RooArgSet(mass),nSigEvents_dd);
  RooDataSet* sig_dd_Lbbar_mass_data = generatetoy("sig_dd_Lbbar_mass_data","sig_dd_Lbbar_mass_data",sig_dd_mass_model,RooArgSet(mass),nSigEvents_dd_Lbbar);
  RooDataSet* sig_ll_mass_data       = generatetoy("sig_ll_mass_data",      "sig_ll_mass_data",      sig_ll_mass_model,RooArgSet(mass),nSigEvents_ll);
  RooDataSet* sig_ll_Lbbar_mass_data = generatetoy("sig_ll_Lbbar_mass_data","sig_ll_Lbbar_mass_data",sig_ll_mass_model,RooArgSet(mass),nSigEvents_ll_Lbbar);

  if (debug) {
    //    sig_ang_data->Print("v");
    mygetchar;
    //    sig_mass_data->Print("v");
    mygetchar;
  }
  sig_dd_ang_data->merge(sig_dd_mass_data);
  sig_dd_Lbbar_ang_data->merge(sig_dd_Lbbar_mass_data);
  sig_ll_ang_data->merge(sig_ll_mass_data);
  sig_ll_Lbbar_ang_data->merge(sig_ll_Lbbar_mass_data);
  if (debug) {
    //    sig_ang_data->Print("v");
    mygetchar;
  }

  if (inclbkg) {
    //    RooDataSet* bkg_ang_data=NULL;
    if (usesbforbkg) {
      
      // RooDataSet* bkg_ang_data_full = RooDataSet::read("bkgdata.txt", RooArgList(costheta,costheta1,costheta2));
      
      // if ((iSample+1)*nBkgEvents> bkg_ang_data_full->numEntries()) {
      // 	std::cerr << "Not enough bkg events" << std::endl;
      // 	exit(1);
      // }
    
      // RooDataSet* bkg_ang_data = (RooDataSet*)bkg_ang_data_full->reduce(RooFit::EventRange(iSample*nBkgEvents, (iSample+1)*nBkgEvents));
      // if (debug) {
      // 	bkg_ang_data->Print("v");
      // 	mygetchar;
      // }
      // RooDataSet* bkg_mass_data = bkg_mass_model.generate(RooArgList(mass),nBkgEvents);
      // if (debug) {
      // 	bkg_mass_data->Print("v");
      // 	mygetchar;
      // }
      // bkg_ang_data->merge(bkg_mass_data);
      // if (debug) {
      // 	bkg_ang_data->Print("v");
      // 	mygetchar;
      // }
      // sig_ang_data->append(*bkg_ang_data);
      // if (debug) {
      // 	sig_ang_data->Print("v");
      // 	mygetchar;
      // }
    } else {
      bkg_costheta_c0.setVal(rnd->Rndm()/2.-0.25);
      bkg_costheta_c1.setVal(rnd->Rndm()/2.-0.25);
      bkg_costheta_c2.setVal(rnd->Rndm()/2.-0.25);
      bkg_costheta1_c0.setVal(rnd->Rndm()/2.-0.25);
      bkg_costheta1_c1.setVal(rnd->Rndm()/2.-0.25);
      bkg_costheta1_c2.setVal(rnd->Rndm()/2.-0.25);
      bkg_costheta2_c0.setVal(rnd->Rndm()/2.-0.25);
      bkg_costheta2_c1.setVal(rnd->Rndm()/2.-0.25);
      bkg_costheta2_c2.setVal(rnd->Rndm()/2.-0.25);
      RooDataSet* bkg_dd_data       = generatetoy("bkg_dd_data",      "bkg_dd_data",      bkg_dd_model,      RooArgSet(mass,costheta,costheta1,costheta2),nBkgEvents_dd);
      RooDataSet* bkg_dd_Lbbar_data = generatetoy("bkg_dd_Lbbar_data","bkg_dd_Lbbar_data",bkg_dd_Lbbar_model,RooArgSet(mass,costheta,costheta1,costheta2),nBkgEvents_dd_Lbbar);
      RooDataSet* bkg_ll_data       = generatetoy("bkg_ll_data",      "bkg_ll_data",      bkg_ll_model,      RooArgSet(mass,costheta,costheta1,costheta2),nBkgEvents_ll);
      RooDataSet* bkg_ll_Lbbar_data = generatetoy("bkg_ll_Lbbar_data","bkg_ll_Lbbar_data",bkg_ll_Lbbar_model,RooArgSet(mass,costheta,costheta1,costheta2),nBkgEvents_ll_Lbbar);

      sig_dd_ang_data->append(*bkg_dd_data);
      sig_dd_Lbbar_ang_data->append(*bkg_dd_Lbbar_data);
      sig_ll_ang_data->append(*bkg_ll_data);
      sig_ll_Lbbar_ang_data->append(*bkg_ll_Lbbar_data);
    }
  } else {
    // nBkg=0.;
    // nBkg.setConstant();
    // c0.setConstant();
    // bkg_costheta_c0.setConstant();
    // bkg_costheta_c1.setConstant();
    // bkg_costheta_c2.setConstant();
    // bkg_costheta1_c0.setConstant();
    // bkg_costheta1_c1.setConstant();
    // bkg_costheta1_c2.setConstant();
    // bkg_costheta2_c0.setConstant();
    // bkg_costheta2_c1.setConstant();
    // bkg_costheta2_c2.setConstant();
  }

  piminus_TRACK_Type.setVal(5.);
  sig_dd_ang_data->addColumn(piminus_TRACK_Type);
  sig_dd_Lbbar_ang_data->addColumn(piminus_TRACK_Type);
  piminus_TRACK_Type.setVal(3.);
  sig_ll_ang_data->addColumn(piminus_TRACK_Type);
  sig_ll_Lbbar_ang_data->addColumn(piminus_TRACK_Type);


  RooDataSet* tot_dd_data = sig_dd_ang_data;
  RooDataSet* tot_dd_Lbbar_data = sig_dd_Lbbar_ang_data;
  RooDataSet* tot_ll_data = sig_ll_ang_data;
  RooDataSet* tot_ll_Lbbar_data = sig_ll_Lbbar_ang_data;

  std::map<std::string,RooDataSet*> mapping;
  mapping[ "ddLb"    ] = tot_dd_data;
  mapping[ "ddLbbar" ] = tot_dd_Lbbar_data;
  mapping[ "llLb"    ] = tot_ll_data;
  mapping[ "llLbbar" ] = tot_ll_Lbbar_data;

  RooFitResult* fr_mass=NULL;
  
  RooDataSet simdata("simdata", "simdata", *tot_dd_data->get(), RooFit::Index(cat),RooFit::Import(mapping));
  simdata.Print("v");
  //  simdata.write("simdata.txt");
  RooDataSet* tot_data       = &simdata;
  if (simplefit) {
    // RooRealVar w("w","w",0.);
    // if (inclacc) addacceptanceweight(tot_data,w,false,true,acceptance,acceptance);
    // RooDataSet* tot_data_w = new RooDataSet("dataw","dataw",sig_ang_data,RooArgList(mass,costheta,costheta1,costheta2,w),0,"w");
    // if (debug) {
    //   tot_data_w->Print("v");
    //   mygetchar;
    // }
    
    // tot_mass_model.fitTo(*tot_data_w,RooFit::SumW2Error(inclacc), RooFit::Minos(false));
    // fr = tot_model.fitTo(*tot_data_w,RooFit::SumW2Error(inclacc), RooFit::Save(true), RooFit::Minos(false));
    // fr->Print("v");
    // std::cout << "status=" << fr->status() << " ; covqual=" << fr->covQual() << std::endl; 
    
    // if (doplot) {
    //   TCanvas c("c","c",100,100);
    //   //      plotwpull(c,RooArgList(mass,costheta,costheta1,costheta2), tot_model, dataw, "toy-" + istr(iSample,"04") );
    //   const RooArgList pdflist(sig_model,bkg_model);
    //   plot(c,RooArgList(mass,costheta,costheta1,costheta2), tot_model, tot_data_w, "toy-" + istr(iSample,"04") + "-fr",false, &pdflist);
    //   plot(c,RooArgList(mass,costheta,costheta1,costheta2), tot_model, tot_data_w, "toy-" + istr(iSample,"04") + "-sr", false, &pdflist,"signalregion");
    // }

  } else {
    if (inclbkg) {
      fr_mass = tot_mass_model.fitTo(*tot_data, RooFit::Save(true), RooFit::Minos(false),massfitconstrains);
      
      if ((iSample<10) and (doplot)) {
	const RooArgList pdflist(sig_dd_mass_model,sig_ll_mass_model,bkg_dd_mass_model,bkg_ll_mass_model,bkg_dd_Lbbar_mass_model,bkg_ll_Lbbar_mass_model);
	plot(c,RooArgList(mass), tot_mass_model, tot_data, "toy-" + istr(iSample,"04") + "-fr",false,&pdflist,"", &cat);
      }

      RooArgSet* blu = tot_mass_model.getVariables();
      TIterator* iter = blu->createIterator();
      //  for (int i=0;blu.getSize();i++) std::cout << blu.at(i)->GetName() << std::endl;
      RooRealVar* arg;
      //const everything except nSig, nBkg and _M
      std::vector<RooRealVar*> varstounconst;
      RooArgList yieldlist_dd;
      RooArgList yieldlist_dd_Lbbar;
      RooArgList yieldlist_ll;
      RooArgList yieldlist_ll_Lbbar;
      std::cout << "Making constants: ";
      while((arg=(RooRealVar*)iter->Next())) {
	TString varname(arg->GetName());
	if (!varname.Contains("_M") and !varname.Contains("mass") and !varname.Contains("cat") and !varname.Contains("nSig") and !varname.Contains("nBkg") and (!arg->isConstant())) {
	  std::cout << varname << " ";
	  arg->setConstant();
	  varstounconst.push_back(arg);
	}
	if (varname.Contains("nSig") or varname.Contains("nBkg")) {
	  if (varname.Contains("ll")) {
	    if (varname.Contains("Lbbar")) {yieldlist_ll_Lbbar.add(*arg);} else {yieldlist_ll.add(*arg);}
	  } else {
	    if (varname.Contains("Lbbar")) {yieldlist_dd_Lbbar.add(*arg);} else {yieldlist_dd.add(*arg);}
	  }
	}
	
      }
      std::cout << std::endl;
      //      RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
      //						   *tot_data, &tot_mass_model, yieldlist );

      
      double Psample_dd      =calcPsample(tot_dd_data,       mass, &sig_dd_mass_model, &bkg_dd_mass_model,       nSig_dd.getVal()/nBkg_dd.getVal());
      double Psample_dd_Lbbar=calcPsample(tot_dd_Lbbar_data, mass, &sig_dd_mass_model, &bkg_dd_Lbbar_mass_model, nSig_dd_Lbbar.getVal()/nBkg_dd_Lbbar.getVal());
      double Psample_ll      =calcPsample(tot_ll_data,       mass, &sig_ll_mass_model, &bkg_ll_mass_model,       nSig_ll.getVal()/nBkg_ll.getVal());
      double Psample_ll_Lbbar=calcPsample(tot_ll_Lbbar_data, mass, &sig_ll_mass_model, &bkg_ll_Lbbar_mass_model, nSig_ll_Lbbar.getVal()/nBkg_ll_Lbbar.getVal());

      Psample_dd=1.;
      Psample_dd_Lbbar=1.;
      Psample_ll=1.;
      Psample_ll_Lbbar=1.;
      mygetchar;
      
      if (tot_dd_data->numEntries()>0.) {
	RooStats::SPlot sData1("sData1","An SPlot",
			       *tot_dd_data, &tot_dd_mass_model, yieldlist_dd );
	copyvar(tot_dd_data,"nSig_dd_sw",wmass,Psample_dd);
	if (debug) {
	  sData1.Print("v");
	  mygetchar;
	}
      }
      if (tot_dd_Lbbar_data->numEntries()>0.) {
	RooStats::SPlot sData2("sData2","An SPlot",
			       *tot_dd_Lbbar_data, &tot_dd_Lbbar_mass_model, yieldlist_dd_Lbbar );
	copyvar(tot_dd_Lbbar_data,"nSig_dd_Lbbar_sw",wmass,Psample_dd_Lbbar);
      }
      if (tot_ll_data->numEntries()>0.) {
	RooStats::SPlot sData3("sData3","An SPlot",
			       *tot_ll_data, &tot_ll_mass_model, yieldlist_ll );
	copyvar(tot_ll_data,"nSig_ll_sw",wmass,Psample_ll);
      }
      if (tot_ll_Lbbar_data->numEntries()>0.) {
	RooStats::SPlot sData4("sData4","An SPlot",
			       *tot_ll_Lbbar_data, &tot_ll_Lbbar_mass_model, yieldlist_ll_Lbbar );
	copyvar(tot_ll_Lbbar_data,"nSig_ll_Lbbar_sw",wmass,Psample_ll_Lbbar);
      }
      std::cout << "Sum of wmass in tot_dd_data       = " << sumofvarindata(tot_dd_data,"wmass") << std::endl;
      std::cout << "Sum of wmass in tot_dd_Lbbar_data = " << sumofvarindata(tot_dd_Lbbar_data,"wmass") << std::endl;
      std::cout << "Sum of wmass in tot_ll_data       = " << sumofvarindata(tot_ll_data,"wmass") << std::endl;
      std::cout << "Sum of wmass in tot_ll_Lbbar_data = " << sumofvarindata(tot_ll_Lbbar_data,"wmass") << std::endl;


    }
    // if (inclacc) {
    //   addacceptanceweight(tot_dd_data,      wacc,false,false,acceptancell,acceptancedd);
    //   addacceptanceweight(tot_dd_Lbbar_data,wacc,false,false,acceptancell,acceptancedd);
    //   addacceptanceweight(tot_ll_data,      wacc,false,false,acceptancell,acceptancedd);
    //   addacceptanceweight(tot_ll_Lbbar_data,wacc,false,false,acceptancell,acceptancedd);
    //   std::cout << "Sum of wacc in tot_dd_data        = " << sumofvarindata(tot_dd_data,"wacc") << std::endl;
    //   std::cout << "Sum of wacc in tot_dd_Lbbar_data  = " << sumofvarindata(tot_dd_Lbbar_data,"wacc") << std::endl;
    //   std::cout << "Sum of wacc in tot_ll_data        = " << sumofvarindata(tot_ll_data,"wacc") << std::endl;
    //   std::cout << "Sum of wacc in tot_ll_Lbbar_data  = " << sumofvarindata(tot_ll_Lbbar_data,"wacc") << std::endl;
    // }
    // if (inclacc and inclbkg) {
    //   multiplyweights(tot_dd_data,      wtot,wacc.GetName(),wmass.GetName());
    //   multiplyweights(tot_dd_Lbbar_data,wtot,wacc.GetName(),wmass.GetName());
    //   multiplyweights(tot_ll_data,      wtot,wacc.GetName(),wmass.GetName());
    //   multiplyweights(tot_ll_Lbbar_data,wtot,wacc.GetName(),wmass.GetName());
    //   std::cout << "Sum of wtot in tot_dd_data        = " << sumofvarindata(tot_dd_data,"wtot") << std::endl;
    //   std::cout << "Sum of wtot in tot_dd_Lbbar_data  = " << sumofvarindata(tot_dd_Lbbar_data,"wtot") << std::endl;
    //   std::cout << "Sum of wtot in tot_ll_data        = " << sumofvarindata(tot_ll_data,"wtot") << std::endl;
    //   std::cout << "Sum of wtot in tot_ll_Lbbar_data  = " << sumofvarindata(tot_ll_Lbbar_data,"wtot") << std::endl;
    // }
    
    RooDataSet simdata2("simdata2", "simdata2", *tot_dd_data->get(), RooFit::Index(cat),RooFit::Import(mapping));
    simdata2.Print("v");
    //    simdata2.write("simdata2.txt");
    
    tot_data=&simdata2;

    //    RooDataSet data_sig("data_sig","data_sig",tot_data,*tot_data->get());//,0,"nSig_sw") ;
    if (debug) {
      tot_data->Print("v");
      mygetchar;
    }
    RooDataSet* tot_data_w=NULL;
    if (inclbkg  and inclacc)  tot_data_w = new RooDataSet("dataw","dataw",tot_data,*tot_data->get(),0,wmass.GetName());
    if (inclbkg  and !inclacc) tot_data_w = new RooDataSet("dataw","dataw",tot_data,*tot_data->get(),0,wmass.GetName());
    if (!inclbkg and inclacc)  tot_data_w = tot_data;//new RooDataSet("dataw","dataw",tot_data,*tot_data->get(),0,wacc.GetName());
    if (!inclbkg and !inclacc) tot_data_w = tot_data;

    if (debug) {
      tot_data_w->Print("v");
      mygetchar;
    }
    bool useroominuit=false;
    if (!useroominuit) {  
      RooCmdArg minimizer=RooCmdArg::none();
      //      RooCmdArg minimizer=RooFit::Minimizer("Minuit2","minimize");
      //      RooFitResult* fr_test = ang_model.fitTo(*tot_data_w, RooFit::SumW2Error(true), RooFit::Save(true), RooFit::Minos(false), minimizer, angfitconstrains);
      //      fr_ang_nosumw2_nominos = ang_model.fitTo(*tot_data_w, RooFit::Save(true), RooFit::Minos(false), RooFit::SumW2Error(false), minimizer, angfitconstrains);
      //      fr_ang_nosumw2 = ang_model.fitTo(*tot_data_w, RooFit::Save(true), RooFit::Minos(true), RooFit::SumW2Error(false), minimizer, angfitconstrains);
      //      fr_ang_sumw2   = ang_model.fitTo(*tot_data_w, RooFit::Save(true), RooFit::Minos(false), RooFit::SumW2Error(inclbkg or inclacc), minimizer, angfitconstrains);

      for (unsigned j=0; j<angfitv.size();j++) angfitv[j].fr = angfitv[j].model->fitTo(*tot_data_w, RooFit::Save(true), angfitv[j].minosarg, angfitv[j].sumw2arg, minimizer, angfitconstrains);

      
      //      std::cout << " false "  << std::endl;

      //      std::cout << " true "  << std::endl;
      //      fr_test->Print("v");
      //      mygetchar;
      //      std::cout << " *** my fitTo *** " << std::endl;
      //      myfitto(tot_data_w, ang_model);

    } else {
      ang_model.fitTo(*tot_data_w,RooFit::SumW2Error(inclbkg or inclacc), RooFit::Minos(false));
      RooNLLVar nll("nll","-log(L)",ang_model,*tot_data);//,NumCPU(nCPU));
      //      nll.applyWeightSquared(true);
      nll.Print("v");
      //    getchar();
      
      //      RooMinuit m1(nll);
      RooMinimizer m1(nll) ;
      m1.setVerbose(kFALSE) ;
      //      m1.setProfile(1);
      m1.setPrintLevel(3);
      m1.setEps(1e-12);
      m1.setStrategy(1);
      //    m1.save()->Print("v");
      //      m1.simplex();
      //    getchar();
      //      m1.migrad();
      m1.hesse();
      m1.minimize("OldMinuit","minuit") ;
      m1.hesse();
      //      m1.migrad();
      //      m1.hesse();
      //    RooFitResult* r1 = m1.save();
      //     if(false){
      //       m1.hesse();
      //       if (m1.save()->covQual() != 3) {
      // 	m1.migrad();
      // 	m1.hesse();
      //       }
      //     }
    }

    if (iSample<10) {
      if (doplot) {
	//      plotwpull(c,RooArgList(mass,costheta,costheta1,costheta2), tot_model, dataw, "toy-" + istr(iSample,"04") );
	plot(c,RooArgList(costheta,costheta1,costheta2), ang_model, tot_data_w, "toy-" + istr(iSample,"04") + "-fr",false,NULL,"",&cat);
	//      plotwcutwcomponents(c,RooArgList(mass,costheta,costheta1,costheta2), tot_model, RooArgList(sig_model,bkg_model),tot_data_w, "signalregion", "toy-" + istr(iSample,"04") + "-sr" );
	tot_data_w->Print("v");
	tot_data_w->write("toy-" + istr(iSample,"04") + ".txt");
      }
    }
  }

  TFile outputfile(TString("toy_" + istr(iSample,"04") + "_result.root"), "recreate");

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
