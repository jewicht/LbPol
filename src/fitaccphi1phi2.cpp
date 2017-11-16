#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TF3.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TCut.h"
#include "TProfile.h"
#include "TRandom3.h"

#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooProdPdf.h"
#include "RooNDKeysPdf.h"
#include "RooChi2Var.h"

#include "rootlogon.h"
#include "createRooLegendre2.h"
#include "functions-roofit.h"
//#include "calcacceptanceclass.h"
#include "cuts.h"
#include "legendre2f2.h"
//include "weighting.h"
#include "ConfigFile.h"


const bool usetmvarectcut=false;
bool usebdt=true;
bool usesimplecut=false;
bool userwmc=true;
bool runongrid=true;
const bool plotonly=false;

const TString rootdir("/afs/cern.ch/user/j/jwicht/fitLambdab2JpsiL0/");

void addtodata(RooRealVar& costheta, RooRealVar& costheta1, RooRealVar& costheta2,RooDataSet* data, double a, double b, double c)
{
  costheta.setVal(a);
  costheta1.setVal(b);
  costheta2.setVal(c);
  data->add(RooArgSet(costheta,costheta1,costheta2));
  std::cout << "Adding costheta=" << a << " ; costheta1= " << b << " ; costheta2 = " << c << std::endl;
}


void addweight(RooDataSet* data,RooRealVar& w)
{  
  RooDataSet dataww6("dataww6","dataww6",RooArgList(w));
  const int nEntries=data->numEntries();
 
  //  double weights[nEntries];
  double costheta2;
  for (int i=0;i<nEntries;i++) {  
    const RooArgSet* blu = data->get(i);
    costheta2=blu->getRealValue("costheta2");
    if (fabs(costheta2)>0.99) {
      w.setVal(1.);
    } else {
      w.setVal(1./(1.-costheta2*costheta2));
    }
    //    w.setVal((1.-costheta2*costheta2));
    dataww6.add(RooArgSet(w));
  }
  data->merge(&dataww6);
}

void usage(char* argv0)
{
  std::cerr << argv0 << " Lb/B0 LL/DD itoy" << std::endl;
  exit(1);
}

int main(int argc, char** argv)
{

  if (argc!=4) usage(argv[0]);

  const TString mode(argv[1]);
  if (mode!="B0" and mode!="Lb") usage(argv[0]);
  const bool doBd2JpsiKS0=mode.Contains("B0");

  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0,Lambdabname,Lambda0name,Jpsiname);

  TString llordd(argv[2]);
  llordd.ToLower();
  if (llordd!="dd" and llordd!="ll") usage(argv[0]);
  const bool doLL=llordd.Contains("ll");

  const int itoy=atoi(argv[3]);

  ConfigFile cf("fitaccphi1phi2config.txt");
  userwmc=atob(cf.Value("main","userwmc","true").c_str());
  runongrid=atob(cf.Value("main","runongrid","false").c_str());
  usebdt=atob(cf.Value("main","usebdt","true").c_str());
  usesimplecut=atob(cf.Value("main","usesimplecut","false").c_str());


  const std::string confsec(mode+"-"+llordd);
  const int Lb_i=atoi(cf.Value(confsec,"i","0").c_str());
  const int Lb_j=atoi(cf.Value(confsec,"j","0").c_str());
  const int Lb_corr=atoi(cf.Value(confsec,"corr","0").c_str());
  const int Lb_tot=atoi(cf.Value(confsec,"tot","0").c_str());

  const int B0_i=atoi(cf.Value(confsec,"i","0").c_str());
  const int B0_j=atoi(cf.Value(confsec,"j","0").c_str());
  const int B0_corr=atoi(cf.Value(confsec,"corr","0").c_str());
  const int B0_tot=atoi(cf.Value(confsec,"tot","0").c_str());

  if (!userwmc and itoy>0) usage(argv[0]);


  //  rootlogon();
  lhcbstyle();
  TCanvas c("c","c",100,100);

  RooRealVar phi1("phi1", "#phi_{1}/#pi", -1., 1.);// RooArgList(theta1));
  RooRealVar phi2("phi2", "#phi_{2}/#pi", -1., 1.);// RooArgList(theta2));
  //  RooRealVar fakevar("fakevar", "fakevar",  1.);// RooArgList(theta2));
  
  phi1.setRange("phi1pos",0.,1.);
  phi1.setRange("phi1neg",-1.,0.);
  phi2.setRange("phi2pos",0.,1.);
  phi2.setRange("phi2neg",-1.,0.);

  //  phi2.setRange("fitrange",-0.5,0.5);

  RooRealVar Lambdab_ID(Lambdabname + "_ID",Lambdabname + "_ID",-10000.,10000.);
  RooRealVar Polarity("Polarity","Polarity",-10.,10.);

  RooRealVar *cc0c1c2[ordermax2i][ordermax2j];
  RooAbsPdf *c0c1c2;


  if (doBd2JpsiKS0) {
    createRRV2(doBd2JpsiKS0,cc0c1c2, "B0"+llordd+"c0c1c2", B0_i, B0_j,  B0_tot, B0_corr);  //7,5,7,8,2);
    createRooLegendre2(phi1, phi2, "B0"+llordd+"c0c1c2", c0c1c2, cc0c1c2);
  } else {
    createRRV2(doBd2JpsiKS0,cc0c1c2, "Lb"+llordd+"c0c1c2",  Lb_i, Lb_j, Lb_tot, Lb_corr);  //7,5,7,8,2); 
    createRooLegendre2(phi1, phi2, "Lb"+llordd+"c0c1c2", c0c1c2, cc0c1c2);
  }

  RooAbsPdf* model=c0c1c2;

//   const TString filename("Tuple-Bd2JpsiKS-NoPol.root");
//   const TString dirname2("BtoqllGeneratedDistributions");
  
//   TFile file(filename);
//   TDirectory* dir2 = (TDirectory*)file.Get(dirname2);
//   TTree* sigtree2 = (TTree*)dir2->Get("B02JpsiKS");
  
//   RooDataSet* data  = filldataset("data","data", RooArgList(phi,phi1,phi2), sigtree2, "1==1" );

//   RooRealVar c0("c0","c0",0.,-1.,1.);
//   RooRealVar c1("c1","c1",0.,-1.,1.);
//   RooRealVar c2("c2","c2",0.,-1.,1.);
//   //  RooChebychev pol2("pol2","pol2",phi2,RooArgList(c0,c1,c2));
//   RooPolynomial pol2("pol2","pol2",phi2,RooArgList(c0,c1,c2));
  
//   RooGenericPdf genpdf("genpdf","genpdf","1.-phi2*phi2",RooArgList(phi2));

//   RooFitResult* blu = pol2.fitTo(*data,RooFit::Save(true));
//   blu->Print("v");
//   plot(c,phi2,genpdf,data,"blublu-phi2.eps");
//   return 0;



//   RooRealVar xvar("xvar","xvar",-5.,5.);  
//   RooRealVar mean("mean","mean",0.,-1.,1.);  
//   RooRealVar width("width","width",1.,0.,2.);
  
//   RooGaussian gauss("gauss","gauss",xvar,mean,width);
  
//   RooDataSet* data = gauss.generate(RooArgSet(xvar),10000);

//   gauss.fitTo(*data);

//   xvar.setBins(100);
//   plotwpull(c,xvar,gauss,data,"gauss.eps");

//   return 1;



  TString sigfilename(rootdir + "DVTuples-Lb2JpsiL0-MC11a-reco12a.reduced.root");

  TString dirname("Tuple_JpsiL0_betas");
  if (doBd2JpsiKS0) {
    sigfilename=rootdir + "DVTuples-Bd2JpsiKS-MC11a.reduced.root";
    dirname="Tuple_JpsiKs_detached_betas";
    //    dirname="Tuple_JpsiKs_prescaled_betas";
  }
  //  const TString datafilename("DVTuples_data_stripping17.root");

  if (runongrid) sigfilename="rfio:/castor/cern.ch/grid/lhcb/user/j/jwicht/" + sigfilename;

  TFile* sigfile = TFile::Open(sigfilename);
  //  TFile* sigfile = new TFile(sigfilename);
  TDirectory* sigdir = (TDirectory*)sigfile->Get(dirname);
  TTree* sigtree = (TTree*)sigdir->Get("DecayTree");

  TString siganglesfilename(sigfilename);
  sigtree->AddFriend(dirname, TString(sigfilename).ReplaceAll(".root", ".angles.root"));
  if (userwmc) sigtree->AddFriend(dirname, TString(sigfilename).ReplaceAll(".root", ".reweighted.root"));
  if (usebdt) sigtree->AddFriend(dirname, TString(sigfilename).ReplaceAll(".root", ".bdt.root"));



  //cut based selection
  TCut LL_selection, DD_selection;
  simplecut(doBd2JpsiKS0, LL_selection, DD_selection);



  TCut TrueLambdab=truemcsel(doBd2JpsiKS0);

  TCut TriggerSel=triggercut(doBd2JpsiKS0);
  
  if (usebdt and usetmvarectcut) {
    std::cerr << "Both bdt and tmvarectcut are set ... " << std::endl;
    return  1;
  }
  if (usebdt) bdtsel(doBd2JpsiKS0, LL_selection,DD_selection);
  if (usetmvarectcut) tmvarectcut(doBd2JpsiKS0, LL_selection,DD_selection);
  if (usesimplecut) simplecut2(doBd2JpsiKS0, LL_selection,DD_selection);

  LL_selection+=TrueLambdab;
  DD_selection+=TrueLambdab;

  LL_selection+=TriggerSel;
  DD_selection+=TriggerSel;


  if (userwmc) {
    TCut rwsel("rw==1");

    if (itoy>0) {
      rwsel=TCut("rw"+istr(itoy)+"==1");
    }

    LL_selection+=rwsel;
    DD_selection+=rwsel;
  }

  std::cout << " selection (LL) = " << LL_selection << std::endl;


  TString outputname("fitaccphi1phi2-");
  
  if (userwmc) {
    outputname+="mcrw-";
  } else {
    outputname+="mcnonrw-";
  }

  if (usebdt) {
    outputname+="wbdt-";
  } 
  if (usetmvarectcut) {
    outputname+="wtmvarectcut-";
  }
  if (usesimplecut) {
    outputname+="wsimplecut-";
  }
  if (doBd2JpsiKS0) {
    outputname+="B0-";
  } else {
    outputname+="Lb-";
  }

  outputname+=dirname + "-";

  if (doLL) {
    outputname+="ll";
  } else {
    outputname+="dd";
  }

  if (itoy>0) {      
    readcoeff2_rrv(outputname+".root", cc0c1c2);
    outputname+="-toy" + istr(itoy,"06");
  }

  //  return 1;
  phi1.setBins(25);
  phi2.setBins(25);

  RooDataSet* datau; 
  if (doLL) {
    datau = filldataset("datall","datall", RooArgList(phi1,phi2, Polarity, Lambdab_ID), sigtree, LL_selection);
  } else {
    datau = filldataset("datadd","datadd", RooArgList(phi1,phi2, Polarity, Lambdab_ID), sigtree, DD_selection);
  }

  datau=(RooDataSet*)datau->reduce(RooFit::EventRange(0,10000));
  sigfile->Close();

  //  datau->addColumn(fakevar);

  RooRealVar w("w","w",1.);

  //  RooAbsData* data= & datab;
  RooAbsData* data= datau;

  if (doBd2JpsiKS0) {
    addweight(datau,w);
    data=new RooDataSet("dataw","dataw",datau, *datau->get(),0,"w");
    //    plotdata(c,RooArgList(costheta,costheta1,costheta2),(RooDataSet*)data,"blu");
  } else {
    //    data=new RooDataSet("dataw","dataw",datau, *datau->get(),0,"fakevar");
  }

  //  plotdata(c, RooArgList(costheta,costheta1,costheta2), (RooDataSet*)data, "whatsgoingon");
  //  RooRealVar w("w","w",0.);
  //  if (doBd2JpsiKS0) {
    //    addaccweight(datall,w);
    //    addaccweight(datadd,w);
  //  }



//   RooDataSet* datall_antisel = filldataset("datall_antisel","datall_antisel", RooArgList(costheta,costheta1,costheta2), sigtree, LL_antisel );
//   RooDataSet* datadd_antisel = filldataset("datadd_antisel","datadd_antisel", RooArgList(costheta,costheta1,costheta2), sigtree, DD_antisel );
 

//   plotdatavsdata(c, RooArgList(costheta,costheta1,costheta2), datall, datadd, "fitaccphi1phi2-llvsdd");

//   plotdatavsdata(c, RooArgList(costheta,costheta1,costheta2), datall, datall_antisel, "fitaccphi1phi2-ll-bdtcheck");
//   plotdatavsdata(c, RooArgList(costheta,costheta1,costheta2), datadd, datadd_antisel, "fitaccphi1phi2-dd-bdtcheck");



  
  //  mygetchar;
  if (!plotonly) {    
    //    if (!factorizedpdf) readcoeff3_rrv(coeff, cc0c1c2);
    
    //    RooCmdArg fitrange = RooCmdArg::none();
    //    if (doBd2JpsiKS0) fitrange = RooFit::Range("fitrange");
    
    RooFitResult* frc0c1c2 = model->fitTo(*data, RooFit::PrintLevel(-1), RooFit::Save(true), RooFit::Hesse(itoy==0), RooFit::Minos(false),RooFit::SumW2Error(false));//, RooFit::Minimizer("Minuit2","minimize"));


    if (doBd2JpsiKS0) {
      frc0c1c2->SetName("accB0"+llordd);
    } else {
      frc0c1c2->SetName("accLb"+llordd);
    }
    frc0c1c2->Print("v");
    TFile outputfile(outputname+".root", "recreate");
    frc0c1c2->Write();
    outputfile.Close();
    //    mygetchar;
    
  } else {
    readcoeff2_rrv(outputname+".root", cc0c1c2);
  }
  
  
  if (itoy==0) {
    plot(c,RooArgList(phi1,phi2),*model,data, outputname, true);
  }


}
