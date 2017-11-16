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
#include "createRooLegendre3.h"
#include "functions-roofit.h"
#include "calcacceptanceclass.h"
#include "cuts.h"
#include "legendre3f3.h"
#include "weighting.h"
#include "ConfigFile.h"


const bool usetmvarectcut=false;
bool usebdt=true;
bool usesimplecut=false;
bool userwmc=true;
bool runongrid=true;
const bool plotonly=true;
const bool factorizedpdf=false;

const bool plotpull=false;
const bool plotlhcb=true;


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

  ConfigFile cf("fitaccconfig.txt");
  userwmc=atob(cf.Value("main","userwmc","true").c_str());
  runongrid=atob(cf.Value("main","runongrid","false").c_str());
  usebdt=atob(cf.Value("main","usebdt","true").c_str());
  usesimplecut=atob(cf.Value("main","usesimplecut","false").c_str());


  const std::string confsec(mode+"-"+llordd);
  const int Lb_i=atoi(cf.Value(confsec,"i","0").c_str());
  const int Lb_j=atoi(cf.Value(confsec,"j","0").c_str());
  const int Lb_k=atoi(cf.Value(confsec,"k","0").c_str());
  const int Lb_corr=atoi(cf.Value(confsec,"corr","0").c_str());
  const int Lb_tot=atoi(cf.Value(confsec,"tot","0").c_str());

  const int B0_i=atoi(cf.Value(confsec,"i","0").c_str());
  const int B0_j=atoi(cf.Value(confsec,"j","0").c_str());
  const int B0_k=atoi(cf.Value(confsec,"k","0").c_str());
  const int B0_corr=atoi(cf.Value(confsec,"corr","0").c_str());
  const int B0_tot=atoi(cf.Value(confsec,"tot","0").c_str());

  if (!userwmc and itoy>0) usage(argv[0]);


  //  rootlogon();
  lhcbstyle();
  TCanvas c("c","c",100,100);

  RooRealVar costheta("costheta",   "cos#theta",     -1., 1.);//RooArgList(theta));
  RooRealVar costheta1("costheta1", "cos#theta_{1}", -1., 1.);// RooArgList(theta1));
  RooRealVar costheta2("costheta2", "cos#theta_{2}", -1., 1.);// RooArgList(theta2));
  //  RooRealVar fakevar("fakevar", "fakevar",  1.);// RooArgList(theta2));
  
  costheta.setRange("costhetapos",0.,1.);
  costheta.setRange("costhetaneg",-1.,0.);
  costheta1.setRange("costheta1pos",0.,1.);
  costheta1.setRange("costheta1neg",-1.,0.);
  costheta2.setRange("costheta2pos",0.,1.);
  costheta2.setRange("costheta2neg",-1.,0.);

  //  costheta2.setRange("fitrange",-0.5,0.5);

  RooRealVar Lambdab_ID(Lambdabname + "_ID",Lambdabname + "_ID",-10000.,10000.);
  RooRealVar Polarity("Polarity","Polarity",-10.,10.);

  RooRealVar *cc0c1c2[ordermax3i][ordermax3j][ordermax3k];
  RooAbsPdf *c0c1c2;
  RooAbsPdf* facpdf = NULL;
  RooArgSet* facvars=NULL;


  if (doBd2JpsiKS0) {
    createRRV3(doBd2JpsiKS0,cc0c1c2, "B0"+llordd+"c0c1c2", B0_i, B0_j, B0_k, B0_tot, B0_corr);  //7,5,7,8,2);
    createRooLegendre3v2(costheta, costheta1, costheta2, "B0"+llordd+"c0c1c2", c0c1c2, cc0c1c2);
    facpdf=factorizedacc(costheta,costheta1,costheta2,"B0"+llordd, facvars, 7,5,3);
  } else {
    createRRV3(doBd2JpsiKS0,cc0c1c2, "Lb"+llordd+"c0c1c2",  Lb_i, Lb_j, Lb_k, Lb_tot, Lb_corr);  //7,5,7,8,2); 
    createRooLegendre3v2(costheta, costheta1, costheta2, "Lb"+llordd+"c0c1c2", c0c1c2, cc0c1c2);
      
    facpdf=factorizedacc(costheta,costheta1,costheta2,"Lb"+llordd, facvars,5,4,7);
  }

  RooAbsPdf* model=c0c1c2;
  if (factorizedpdf) model=facpdf;

//   const TString filename("Tuple-Bd2JpsiKS-NoPol.root");
//   const TString dirname2("BtoqllGeneratedDistributions");
  
//   TFile file(filename);
//   TDirectory* dir2 = (TDirectory*)file.Get(dirname2);
//   TTree* sigtree2 = (TTree*)dir2->Get("B02JpsiKS");
  
//   RooDataSet* data  = filldataset("data","data", RooArgList(costheta,costheta1,costheta2), sigtree2, "1==1" );

//   RooRealVar c0("c0","c0",0.,-1.,1.);
//   RooRealVar c1("c1","c1",0.,-1.,1.);
//   RooRealVar c2("c2","c2",0.,-1.,1.);
//   //  RooChebychev pol2("pol2","pol2",costheta2,RooArgList(c0,c1,c2));
//   RooPolynomial pol2("pol2","pol2",costheta2,RooArgList(c0,c1,c2));
  
//   RooGenericPdf genpdf("genpdf","genpdf","1.-costheta2*costheta2",RooArgList(costheta2));

//   RooFitResult* blu = pol2.fitTo(*data,RooFit::Save(true));
//   blu->Print("v");
//   plot(c,costheta2,genpdf,data,"blublu-costheta2.eps");
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


  TString outputname("fitacceptance3-");
  
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
  if (factorizedpdf) outputname+="facpdf-";

  if (doLL) {
    outputname+="ll";
  } else {
    outputname+="dd";
  }

  if (itoy>0) {      
    readcoeff3_rrv(outputname+".root", cc0c1c2);
    outputname+="-toy" + istr(itoy,"06");
  }

  RooArgList varlist(costheta,costheta1,costheta2);
  std::vector<TString> varlist1,varlist2;
  varlist1.push_back("costheta");
  varlist1.push_back("costheta1");
  varlist1.push_back("costheta2");
  varlist2.push_back(Lambdabname + "_DIRA_OWNPV");
  varlist2.push_back(Lambdabname + "_FDCHI2_OWNPV");
  varlist2.push_back(Lambdabname + "_TAU");
  varlist2.push_back(Lambdabname + "_PT");  
  varlist2.push_back(Lambdabname + "_P");
  varlist2.push_back(Lambdabname + "_IPCHI2_OWNPV");
  varlist2.push_back(Lambdabname + "_ENDVERTEX_CHI2/" + Lambdabname + "_ENDVERTEX_NDOF" );
  varlist2.push_back(Lambda0name + "_PT");
  varlist2.push_back(Lambda0name + "_TAU");
  varlist2.push_back(Lambda0name + "_ENDVERTEX_Z");

  if (!doBd2JpsiKS0) {
    varlist2.push_back("piminus_ProbNNk-piminus_ProbNNpi");
    varlist2.push_back("piminus_ProbNNp-piminus_ProbNNpi");
    varlist2.push_back("pplus_ProbNNk-pplus_ProbNNp");
    varlist2.push_back("pplus_ProbNNpi-pplus_ProbNNp");
    varlist2.push_back("piminus_PIDK");
    varlist2.push_back("piminus_PIDp");
    varlist2.push_back("pplus_PIDK");
    varlist2.push_back("pplus_PIDp");
    
    varlist2.push_back("piminus_P");
    varlist2.push_back("pplus_P");
    
    varlist2.push_back("piminus_PT");
    varlist2.push_back("pplus_PT");
    
    varlist2.push_back("bdt");
  }

  TCut profsel = DD_selection;
  if (doLL) profsel = LL_selection;
  
  if (!plotonly) {
    if (itoy==0) plotprofile(c,sigtree,varlist1, varlist1, profsel,outputname + "-profile");
    if (itoy==0) plotprofile(c,sigtree,varlist2, varlist1, profsel,outputname + "-profile");
  }
  //  return 1;
  costheta.setBins(25);
  costheta1.setBins(25);
  costheta2.setBins(25);

  RooDataSet* datau; 
  if (doLL) {
    datau = filldataset("datall","datall", RooArgList(costheta,costheta1,costheta2, Polarity, Lambdab_ID), sigtree, LL_selection);
  } else {
    datau = filldataset("datadd","datadd", RooArgList(costheta,costheta1,costheta2, Polarity, Lambdab_ID), sigtree, DD_selection);
  }
  sigfile->Close();

  if (!doBd2JpsiKS0 and itoy==0) {
    RooDataSet* data_up = (RooDataSet*)datau->reduce(MagnetUp);
    RooDataSet* data_down = (RooDataSet*)datau->reduce(MagnetDown);
    plotdatavsdata(c, RooArgList(costheta,costheta1,costheta2), data_up, data_down, "fitacceptance3-" + llordd + "-MagnetUpvsDown");
    
    const TCut isLambdab0(Lambdabname + "_ID>0.");
    const TCut isLambdab0bar(Lambdabname + "_ID<0.");
    RooDataSet* data_Lb    = (RooDataSet*)datau->reduce(isLambdab0);
    RooDataSet* data_Lbbar = (RooDataSet*)datau->reduce(isLambdab0bar);
    plotdatavsdata(c, RooArgList(costheta,costheta1,costheta2), data_Lb, data_Lbbar, "fitacceptance3-" + llordd + "-LbvsLbbar");
  }
  //  return 0;

  //  datau->addColumn(fakevar);


  //   RooNDKeysPdf keyspdf("keyspdf","keyspdf",RooArgList(costheta,costheta1,costheta2),*datau);
  //   RooDataSet* toydata = keyspdf.generate(RooArgList(costheta,costheta1,costheta2),10000);
  //   plotdata(c,RooArgList(costheta,costheta1,costheta2),toydata, "keyspdf-toy");
  //   plot(c,RooArgList(costheta,costheta1,costheta2),keyspdf,datau, "keyspdf");
  //   return 1;


//   double numEntries = datau->numEntries();
//   int target=50;
//   int nBins=pow(numEntries/target,1./3);
//   //  std::cout << "nEntries = " << sigtree-->GetEntries()
//   std::cout << "numEntries = " << numEntries << std::endl;
//   std::cout << "nBins = " << nBins << std::endl;
//   double entriesperbin=numEntries/nBins;

//   Double_t xbin[nBins+1];
//   xbin[0]=costheta.getMin();
//   xbin[nBins]=costheta.getMax();
//   Double_t ybin[nBins+1];
//   ybin[0]=costheta1.getMin();
//   ybin[nBins]=costheta1.getMax();
//   Double_t zbin[nBins+1];
//   zbin[0]=costheta2.getMin();
//   zbin[nBins]=costheta2.getMax();


//   for (int i=1; i<nBins; i++) {
//     //    std::cout << "Find bin " << i << std::endl;
//     xbin[i]=findnextbin(datau,entriesperbin,"costheta",xbin[i-1],xbin[nBins]);
//     ybin[i]=findnextbin(datau,entriesperbin,"costheta1",ybin[i-1],ybin[nBins]);
//     zbin[i]=findnextbin(datau,entriesperbin,"costheta2",zbin[i-1],zbin[nBins]);
//   }
  
//   TH3F histo("histo","histo",nBins,xbin,nBins,ybin,nBins,zbin);
//   datau->fillHistogram(&histo, RooArgList(costheta,costheta1,costheta2));
  
//   histo.Draw("box");
//   c.SaveAs("box.eps");

//   for (int i=1;i<=nBins;i++) {
//     for (int j=1;j<=nBins;j++) {
//       for (int k=1;k<=nBins;k++) {
// 	std::cout << i << " " << j << " " << k << " " << histo.GetBinContent(i,j,k) << std::endl;
//       }
//     }
//   }

  //  RooDataHist datab("datab","datab",RooArgList(costheta,costheta1,costheta2), RooFit::Import(histo));
  //  return 1;


  //  RooDataHist datallbinned("datallbinned","datallbinned", RooArgList(costheta,costheta1,costheta2), *datall);
  //   RooHistPdf  histpdfll("histpdfll","histpdfll",RooArgList(costheta,costheta1,costheta2),datallbinned,3);
  
  //   plot(c,RooArgList(costheta,costheta1,costheta2),histpdfll,datall,"lolololll");
  //   return 0;
  
  const bool addfakepoints=true;
  if (addfakepoints) {
    if (doBd2JpsiKS0) {
      if (doLL) {

	addtodata(costheta,costheta1,costheta2,datau, -0.624782,  -1.,  0.999987);
	addtodata(costheta,costheta1,costheta2,datau, -0.718326,  0.999999,  -1);
	addtodata(costheta,costheta1,costheta2,datau, 0.000145794,  -0.522963,  -1);
	addtodata(costheta,costheta1,costheta2,datau,-0.615979,  -0.720076,  -1  ); 
	addtodata(costheta,costheta1,costheta2,datau,0.474559 , 0.637265 , -1);
	addtodata(costheta,costheta1,costheta2,datau,0.877533 , -1 , -1);
	addtodata(costheta,costheta1,costheta2,datau,0.000145794 , 0.15527 , -1);
	addtodata(costheta,costheta1,costheta2,datau,-0.387957 , -0.483307 , -1);
	addtodata(costheta,costheta1,costheta2,datau,0.000145794 , 0.397871 , -1);
	addtodata(costheta,costheta1,costheta2,datau,0.000145794 , -0.348885 , -1);
	addtodata(costheta,costheta1,costheta2,datau, -1 , -1 , 0.426739);
	addtodata(costheta,costheta1,costheta2,datau, 0.000145794 , -0.732322 , -1);
	addtodata(costheta,costheta1,costheta2,datau,-1 , -1 , -0.257618);
	addtodata(costheta,costheta1,costheta2,datau, -0.329057 , -1 , -1);
	addtodata(costheta,costheta1,costheta2,datau, -0.473684 , -1 , 0.999999);
	addtodata(costheta,costheta1,costheta2,datau,-0.578947 , -1 , -1);
	addtodata(costheta,costheta1,costheta2,datau,-1 , -1 , -1);
	addtodata(costheta,costheta1,costheta2,datau,-1 , 1 , -1);
	addtodata(costheta,costheta1,costheta2,datau,0.684211 , -1 , -1);
	addtodata(costheta,costheta1,costheta2,datau,0.80345 , -1 , 0.999999);
	addtodata(costheta,costheta1,costheta2,datau,-0.717743 , -1 , 0.43432);
	//	TRandom3 rnd(0);
	//	for (int i=0;i<100;i++) addtodata(costheta,costheta1,costheta2,datau,rnd.Rndm()*2.-1.,rnd.Rndm()*2.-1.,costheta2val);
	//	std::cout << "numEntries = " << datau->numEntries() << std::endl;
	//	datau->write("blu.txt");

      } else {
	addtodata(costheta,costheta1,costheta2,datau,-0.000145794 , 0.34247 , -1);
	addtodata(costheta,costheta1,costheta2,datau,0.886281 , -1 , -1 );
	addtodata(costheta,costheta1,costheta2,datau,0.000145794 , -0.369879 , -1);
	addtodata(costheta,costheta1,costheta2,datau,0.314769 , 0.506634 , -1);
	addtodata(costheta,costheta1,costheta2,datau,-6.37675e-08 , 0.999999 , -1);
	addtodata(costheta,costheta1,costheta2,datau,0.000145794 , -1 , -1);
      }
    } else {
      if (doLL) {
	//	addtodata(costheta,costheta1,costheta2,datau,-0.0263023 , 0.999987 , -1);
	addtodata(costheta,costheta1,costheta2,datau,-1 , 0.999999 , -0.668769);
	addtodata(costheta,costheta1,costheta2,datau,-1 , 0.999999 , 0.669608);
	addtodata(costheta,costheta1,costheta2,datau,0.0414911 , 0.999999 , -0.800555);
	addtodata(costheta,costheta1,costheta2,datau,-0.0526316 , 0.999999 , -1);

	addtodata(costheta,costheta1,costheta2,datau,-1 , 0.999987 , -0.0447224);
	addtodata(costheta,costheta1,costheta2,datau,0.690873 , 0.999999 , 0.0279254);
	addtodata(costheta,costheta1,costheta2,datau,0.0733811 , 0.999999 , 0.763762);
	addtodata(costheta,costheta1,costheta2,datau,0.578946 , 0.999999 , -1);
	addtodata(costheta,costheta1,costheta2,datau,-1 , 0.999999 , -1);
	addtodata(costheta,costheta1,costheta2,datau,-1 , 0.999999 , 0.999999);
	addtodata(costheta,costheta1,costheta2,datau,-0.684211 , 0.999999 , -1);
      } else {

	// addtodata(costheta,costheta1,costheta2,datau,-0.185656,  0.999999,  -1);
	// addtodata(costheta,costheta1,costheta2,datau,0.418158 , 0.999999 , -1);	

      }
    }
  }

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
 

//   plotdatavsdata(c, RooArgList(costheta,costheta1,costheta2), datall, datadd, "fitacceptance3-llvsdd");

//   plotdatavsdata(c, RooArgList(costheta,costheta1,costheta2), datall, datall_antisel, "fitacceptance3-ll-bdtcheck");
//   plotdatavsdata(c, RooArgList(costheta,costheta1,costheta2), datadd, datadd_antisel, "fitacceptance3-dd-bdtcheck");



  
  //  mygetchar;
  if (!plotonly) {    
    //    if (!factorizedpdf) readcoeff3_rrv(coeff, cc0c1c2);
    
    //    RooCmdArg fitrange = RooCmdArg::none();
    //    if (doBd2JpsiKS0) fitrange = RooFit::Range("fitrange");
    
    RooFitResult* frc0c1c2 = model->fitTo(*data, RooFit::Save(true), RooFit::Hesse(itoy==0), RooFit::Minos(false),RooFit::SumW2Error(false));//, RooFit::Minimizer("Minuit2","minimize"));


    costheta.setBins(10);
    costheta1.setBins(10);
    costheta2.setBins(10);
    RooDataHist binneddata("binneddata","binneddata",RooArgList(costheta,costheta1,costheta2), *data);
    RooChi2Var chi2("chi2","chi2",*model,binneddata);
    const double dof = 10*10*10 - frc0c1c2->floatParsFinal().getSize();


    std::ofstream file;
    file.open("chi2dof.txt");
    file << "chi2      = " << chi2.getVal() << std::endl;
    file << "dof       = " << dof << std::endl;
    file << "chi2/dof  = " << chi2.getVal()/dof << std::endl;
    file.close();


    //for plots
    costheta.setBins(25);
    costheta1.setBins(25);
    costheta2.setBins(25);

    

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

    if (!factorizedpdf) {
      //      writecoeff3(coeff+".txt", cc0c1c2);
    } else {
      facvars->writeToFile(outputname+".txt");
    }
  } else {
    if (!factorizedpdf) { 
      readcoeff3_rrv(outputname+".root", cc0c1c2);
    } else {
      facvars->readFromFile(outputname);
    }
    //     costheta=-1;
    //     costheta1=0.829022;
    //     costheta2=-0.792282;
    //     std::cout << "getval = " << model->getVal(RooArgSet(costheta,costheta1,costheta2)) << std::endl;
    //     costheta=1.;
    //     costheta1=1.;
    //     costheta2=-0.;
    //     std::cout << "getval = " << model->getVal(RooArgSet(costheta,costheta1,costheta2)) << std::endl;
    //     mygetchar;
  }
  

  if (!factorizedpdf) calcacceptanceclass mycalcacceptanceclass(outputname+".root");
  if (itoy==0) {
    plot(c,RooArgList(costheta,costheta1,costheta2),*model,data, outputname, plotpull, NULL, "", NULL, plotlhcb);
    if (!plotonly) {
      plot(c, RooArgList(costheta1,costheta2), *model, data,  outputname + "-proj",false,NULL,"costhetaneg");
      plot(c, RooArgList(costheta1,costheta2), *model, data,  outputname + "-proj",false,NULL,"costhetapos");
      plot(c, RooArgList(costheta,costheta2), *model, data,  outputname + "-proj",false,NULL,"costheta1neg");
      plot(c, RooArgList(costheta,costheta2), *model, data,  outputname + "-proj",false,NULL,"costheta1pos");
      plot(c, RooArgList(costheta,costheta1), *model, data, outputname + "-proj",false,NULL,"costheta2neg");
      plot(c, RooArgList(costheta,costheta1), *model, data,  outputname + "-proj",false,NULL,"costheta2pos");
    }
  }

  return 1;
  TRandom3 rnd(0);
  for (int i=0;i<1000000;i++) {
    costheta=rnd.Rndm()*2.-1.;
    costheta1=rnd.Rndm()*2.-1.;
    costheta2=rnd.Rndm()*2.-1.;
    //    std::cout << "getval = " << model->getVal(RooArgSet(costheta,costheta1,costheta2)) << std::endl;
    double val=model->getVal(RooArgSet(costheta,costheta1,costheta2));
    if (val<0. or val>1. or i<10) std::cout << costheta.getVal() << " " << costheta1.getVal() << " " << costheta2.getVal() << " " << val << std::endl;
    
  }


}
