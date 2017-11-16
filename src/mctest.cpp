#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TRandom3.h"
#include "TH3F.h"
#include "TLorentzVector.h"

#include "RooGaussian.h"
#include "RooCategory.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooCBShape.h"
#include "RooSimultaneous.h"
#include "RooPolynomial.h"
#include "RooConstVar.h"
#include "RooRandom.h"
#include "RooMultiVarGaussian.h"
#include "RooUnblindPrecision.h"
#include "RooNDKeysPdf.h"

#include "RooStats/SPlot.h"


#include "rootlogon.h"
#include "functions-roofit.h"
#include "RooLbtoJpsiL0PDF3.h"
#include "RooLbtoJpsiL0PDF5.h"
#include "RooLbtoJpsiL0PDF3v2.h"
#include "RooLbtoJpsiL0PDF3v3.h"
#include "RooLbtoJpsiL0PDF3v4.h"
#include "createRooLegendre2.h"
#include "createRooLegendre3.h"
#include "calcacceptanceclass.h"
#include "weighting.h"
#include "cuts.h"
#include "ConfigFile.h"
#include "angles.h"
#include "legendre2f2.h"
#include "legendre3f3.h"

const double pi=atan(1.)*4.;
const TString rootdir("/afs/cern.ch/user/j/jwicht/fitLambdab2JpsiL0/");

bool userwmc=true;
bool usesimplecut=false;
const bool usetmvarectcut=false;
bool usebdt=true;

const bool normalizetowmass=true;
const bool correctangdistrib=false;
const bool usePID=false;

int dorand;
const bool blindPb=false;

bool extendedmassfit=true;

bool separateLbLbar=false;
bool dontfitLb=false;
bool dontfitLbbar=false;
bool dontfitdd=false;
bool dontfitll=false;

bool floatSigParams=false;
bool floatalphalambda=false;
bool assumecp=true;
int useparam=0;

bool randacceptance=false;
bool randmccalib=false;
bool randsfit=false;
bool randsigparams=false;
bool randalphalambda=false;
bool randbeam=false;
bool randphases=false;


bool donotsubtractll=false;

bool useminos=true;

bool usekeyspdf=true;
RooAbsPdf* accpdfll=NULL;
RooAbsPdf* accpdfdd=NULL;

RooArgList* angvars;

RooRealVar mass("M","M(#LambdaJ/#psi)",5500.,5750.,"MeV/c^{2}");
RooRealVar costheta("costheta",   "cos#theta",     -1., 1.);
RooRealVar costheta1("costheta1", "cos#theta_{1}", -1., 1.);
RooRealVar costheta2("costheta2", "cos#theta_{2}", -1., 1.);
RooRealVar phi1("phi1", "#phi_{1}", -1., 1.);
RooRealVar phi2("phi2", "#phi_{2}", -1., 1.);
RooRealVar wmass("wmass","wmass",0.);
RooRealVar wacc("wacc","wacc",0.);
RooRealVar wcorr("wcorr","wcorr",0.);
RooRealVar wtot("wtot","wtot",0.);
RooCategory cat("cat","cat");


RooRealVar P_b("P_b","P_{b}",0.0, -1.5, 1.5);//0.50
RooRealVar alpha_b("alpha_b","#alpha_{b}",0.0,-1.1,1.1);//0.457
RooRealVar r_0("r_0","r_{0}",0.0,-0.1,1.2);//0.5
RooRealVar r_1("r_1","r_{1}", 0.0,-1.2,1.2);//0.1
RooRealVar alpha_plus("alpha_plus",  "#alpha_{+}", 0.237, -0.5*pi,2.5*pi);
RooRealVar alpha_minus("alpha_minus","#alpha_{-}", 3.08,  -0.5*pi,2.5*pi);
RooRealVar chi("chi","#chi",6.21719,0.,5.*pi);

RooRealVar beta("beta","#beta",0.0,-5.0,1.3);//0.5
RooRealVar gammarrv("gamma","#gamma", 0.,-5.0,5.0);//0.1


RooRealVar ap2("ap2","|a+|^{2}",0.25);
RooRealVar am2("am2","|a-|^{2}",0.25);
RooRealVar bp2_rrv("bp2_rrv","|b+|^{2}",0.25);
RooRealVar bm2("bm2","|b-|^{2}",0.25);
RooFormulaVar bp2("bp2","1.0-ap2-am2-bm2",RooArgList(ap2,am2,bm2));


// RooFormulaVar ap2_rfv("ap2_rrv","|a+|^{2}","0.5*(@0+@1)",RooArgList(r_0,r_1));
// RooFormulaVar am2_rfv("am2_rrv","|a-|^{2}","0.5*(@0-@1)",RooArgList(r_0,r_1));
// RooFormulaVar bp2_rfv("bp2_rrv","|b+|^{2}","0.5*(1.+@0-@1-@2)",RooArgList(alpha_b,r_0,r_1));
// RooFormulaVar bm2_rfv("bm2_rrv","|b-|^{2}","0.5*(1.-@0-@1+@2)",RooArgList(alpha_b,r_0,r_1));


RooFormulaVar* ap2_rfv;
RooFormulaVar* am2_rfv;
RooFormulaVar* bp2_rfv;
RooFormulaVar* bm2_rfv;
 
// RooFormulaVar h1("h1","h1","@0",RooArgList(costheta));
// RooFormulaVar h2("h2","h2","@0",RooArgList(costheta1));
// RooFormulaVar h3("h3","h3","@0*@1",RooArgList(costheta,costheta1));
// RooFormulaVar h4("h4","h4","0.5*(3.*@0*@0-1.)",RooArgList(costheta2));
// RooFormulaVar h5("h5","h5","0.5*(3.*@0*@0-1.)*@1",RooArgList(costheta2,costheta));
// RooFormulaVar h6("h6","h6","0.5*(3.*@0*@0-1.)*@1",RooArgList(costheta2,costheta1));
// RooFormulaVar h7("h7","h7","0.5*(3.*@0*@0-1.)*@1*@2",RooArgList(costheta2,costheta,costheta1));

//RooRealVar* hrrv[8];

RooAbsPdf *corrpdfdd, *corrpdfll;

calcacceptanceclass* acceptancell;
calcacceptanceclass* acceptancedd;

RooCmdArg latexformat = RooFit::Format("NE",RooFit::AutoPrecision(2),RooFit::VerbatimName(true));


RooRealVar lambdab_PE("lambda_b0_PE","lambda_b0_PE",0.);
RooRealVar lambdab_PX("lambda_b0_PX","lambda_b0_PX",0.);
RooRealVar lambdab_PY("lambda_b0_PY","lambda_b0_PY",0.);
RooRealVar lambdab_PZ("lambda_b0_PZ","lambda_b0_PZ",0.);

RooRealVar lambda0_PE("lambda0_PE","lambda0_PE",0.);
RooRealVar lambda0_PX("lambda0_PX","lambda0_PX",0.);
RooRealVar lambda0_PY("lambda0_PY","lambda0_PY",0.);
RooRealVar lambda0_PZ("lambda0_PZ","lambda0_PZ",0.);

RooRealVar piminus_PE("piminus_PE","piminus_PE",0.);
RooRealVar piminus_PX("piminus_PX","piminus_PX",0.);
RooRealVar piminus_PY("piminus_PY","piminus_PY",0.);
RooRealVar piminus_PZ("piminus_PZ","piminus_PZ",0.);

RooRealVar pplus_PE("pplus_PE","pplus_PE",0.);
RooRealVar pplus_PX("pplus_PX","pplus_PX",0.);
RooRealVar pplus_PY("pplus_PY","pplus_PY",0.);
RooRealVar pplus_PZ("pplus_PZ","pplus_PZ",0.);

RooRealVar jpsi_PE("jpsi1s_PE","jpsi1s_PE",0.);
RooRealVar jpsi_PX("jpsi1s_PX","jpsi1s_PX",0.);
RooRealVar jpsi_PY("jpsi1s_PY","jpsi1s_PY",0.);
RooRealVar jpsi_PZ("jpsi1s_PZ","jpsi1s_PZ",0.);

RooRealVar muminus_PE("muminus_PE","muminus_PE",0.);
RooRealVar muminus_PX("muminus_PX","muminus_PX",0.);
RooRealVar muminus_PY("muminus_PY","muminus_PY",0.);
RooRealVar muminus_PZ("muminus_PZ","muminus_PZ",0.);

RooRealVar muplus_PE("muplus_PE","muplus_PE",0.);
RooRealVar muplus_PX("muplus_PX","muplus_PX",0.);
RooRealVar muplus_PY("muplus_PY","muplus_PY",0.);
RooRealVar muplus_PZ("muplus_PZ","muplus_PZ",0.);


// RooArgList* addvariation(RooDataSet* data, RooArgList& varstoplot)
// {

//   RooArgList* variationlist = new RooArgList();

//   double_t val0[varstoplot.getSize()];

//   for (int i=0;i<varstoplot.getSize();i++) {
//     RooRealVar* var = (RooRealVar*)varstoplot.at(i);
//     val0[i]=data->get(0)->getRealValue(var->GetName());
//     RooRealVar* var2 = new RooRealVar(TString(var->GetName()) + "_var", TString(var->GetTitle()) + " variation",0.);
//     variationlist->add(RooArgList(*var2));
//   }
  
//   mygetchar;

//   RooDataSet dataww6("dataww6","dataww6",*variationlist);

//   const int nEntries=data->numEntries();
//   for (int i=0;i<nEntries;i++) {  
//     const RooArgSet* blu = data->get(i);
//     for (int j=0;j<varstoplot.getSize();j++) {
//       RooRealVar* var = (RooRealVar*)variationlist->at(j);
//       var->setVal( blu->getRealValue(varstoplot.at(j)->GetName()) - val0[j] );
//       //      var->setVal( 100. * (val0[j]-blu->getRealValue(varstoplot.at(j)->GetName()))/val0[j]  );
//     }
//     dataww6.add(*variationlist);
//   }
//  data->merge(&dataww6);

//  return variationlist;
// }



void addangles(RooDataSet* data, TRandom3& rnd)
{
  RooDataSet dataww6("dataww6","dataww6",RooArgList(costheta,costheta1,costheta2,phi1,phi2));
  
  TVector3 proddir(0.,0.,1.);
  if (dorand>0) {
    const double cs=rnd.Gaus(0.,0.001);
    proddir=TVector3( sin(cs), 0., cos(cs));
  }

  const int nEntries=data->numEntries();

  Float_t costhetaval;
  Float_t costheta1val;
  Float_t costheta2val;
  Float_t phi1val;
  Float_t phi2val;

  for (int i=0;i<nEntries;i++) {  
    const RooArgSet* blu = data->get(i);
    TLorentzVector lambdab_LV(blu->getRealValue(lambdab_PX.GetName()),
			      blu->getRealValue(lambdab_PY.GetName()),
			      blu->getRealValue(lambdab_PZ.GetName()),
			      blu->getRealValue(lambdab_PE.GetName()));
    TLorentzVector lambda0_LV(blu->getRealValue(lambda0_PX.GetName()),
			      blu->getRealValue(lambda0_PY.GetName()),
			      blu->getRealValue(lambda0_PZ.GetName()),
			      blu->getRealValue(lambda0_PE.GetName()));
    TLorentzVector piminus_LV(blu->getRealValue(piminus_PX.GetName()),
			      blu->getRealValue(piminus_PY.GetName()),
			      blu->getRealValue(piminus_PZ.GetName()),
			      blu->getRealValue(piminus_PE.GetName()));
    TLorentzVector pplus_LV(blu->getRealValue(pplus_PX.GetName()),
			    blu->getRealValue(pplus_PY.GetName()),
			    blu->getRealValue(pplus_PZ.GetName()),
			    blu->getRealValue(pplus_PE.GetName()));
    TLorentzVector jpsi_LV(blu->getRealValue(jpsi_PX.GetName()),
			   blu->getRealValue(jpsi_PY.GetName()),
			   blu->getRealValue(jpsi_PZ.GetName()),
			   blu->getRealValue(jpsi_PE.GetName()));
    TLorentzVector muplus_LV(blu->getRealValue(muplus_PX.GetName()),
			     blu->getRealValue(muplus_PY.GetName()),
			     blu->getRealValue(muplus_PZ.GetName()),
			     blu->getRealValue(muplus_PE.GetName()));
    TLorentzVector muminus_LV(blu->getRealValue(muminus_PX.GetName()),
			      blu->getRealValue(muminus_PY.GetName()),
			      blu->getRealValue(muminus_PZ.GetName()),
			      blu->getRealValue(muminus_PE.GetName()));
    
    calculateangles(lambdab_LV, 
		    lambda0_LV, pplus_LV, piminus_LV,
		    jpsi_LV, muplus_LV, muminus_LV,
		    costhetaval, costheta1val, costheta2val, phi1val, phi2val, proddir);

    costheta=costhetaval;
    costheta1=costheta1val;
    costheta2=costheta2val;
    phi1=phi1val;
    phi2=phi2val;
    dataww6.add(RooArgSet(costheta,costheta1,costheta2,phi1,phi2));
  }
  data->merge(&dataww6);
}


void fillwmassll(RooDataSet* data, RooRealVar& w)
{
  RooDataSet dataww6("dataww6","dataww6",RooArgList(w));
  const int nEntries=data->numEntries();
  double weights[nEntries];
  
  double massval;
  for (int i=0;i<nEntries;i++) {  
    const RooArgSet* blu = data->get(i);
    massval=blu->getRealValue(mass.GetName());
    weights[i]=0.;
    if (massval>5600. and massval<5640.) weights[i]=1.;
    //    std::cout << massval << " " << weights[i] << std::endl;
  }
  
  for (int i=0;i<nEntries;i++) {
    w.setVal(weights[i]);
    dataww6.add(RooArgSet(w));
  }
  data->merge(&dataww6);

}



bool hascategory(std::map<std::string,RooDataSet*>* datamapping, const std::string& category)
{
  return ( datamapping->find(category) != datamapping->end() );
}


void printamplitudes(RooFitResult* fr, TString plotname)
{
  std::cout << "ap2 = " << ap2_rfv->getVal() << " +/- " << ap2_rfv->getPropagatedError(*fr) << std::endl;
  std::cout << "am2 = " << am2_rfv->getVal() << " +/- " << am2_rfv->getPropagatedError(*fr) << std::endl;
  std::cout << "bp2 = " << bp2_rfv->getVal() << " +/- " << bp2_rfv->getPropagatedError(*fr) << std::endl;
  std::cout << "bm2 = " << bm2_rfv->getVal() << " +/- " << bm2_rfv->getPropagatedError(*fr) << std::endl;
  
  ap2.setVal(  ap2_rfv->getVal());
  ap2.setError(ap2_rfv->getPropagatedError(*fr));
  am2.setVal(  am2_rfv->getVal());
  am2.setError(am2_rfv->getPropagatedError(*fr));
  bp2_rrv.setVal(  bp2_rfv->getVal());
  bp2_rrv.setError(bp2_rfv->getPropagatedError(*fr));
  bm2.setVal(  bm2_rfv->getVal());
  bm2.setError(bm2_rfv->getPropagatedError(*fr));
  
  RooArgList tmp(P_b,ap2,am2,bp2_rrv,bm2);
  tmp.printLatex( latexformat, RooFit::OutputFile(plotname + "-angfit-amplitudes.tex" ), RooFit::Columns(1) );
  
  RooArgSet tmp2(tmp);
  tmp2.writeToFile(plotname + "-angfit-amplitudes.txt");
}

void fitwsplot(TCanvas& c, bool doBd2JpsiKS0, RooDataSet* tot_data, RooAbsPdf& tot_mass_model, const RooArgList& masspdfs, RooArgList& angpdflist, TString plotname,
	       const RooCmdArg& extconstmfit=RooCmdArg::none(), const RooCmdArg& extconstangfit=RooCmdArg::none(),
	       std::map<std::string,RooDataSet*>* datamapping=NULL, std::map<std::string,RooAbsPdf*>* masspdfsmapping=NULL
	       )
{
  bool debug=false;
  bool doplot=true;

  RooFitResult* fr_mass = tot_mass_model.fitTo(*tot_data, RooFit::Extended(extendedmassfit), RooFit::Minos(false), RooFit::Save(true), extconstmfit);
  fr_mass->Print("v");
  std::cout << plotname << ": status = " << fr_mass->status() << " ; covQual = " << fr_mass->covQual() <<  std::endl;
  if (dorand==0) {
    RooArgSet constpars(fr_mass->constPars());constpars.writeToFile(plotname + "-massfit-const.txt");
    RooArgSet floatpars(fr_mass->floatParsFinal());floatpars.writeToFile(plotname + "-massfit-float.txt");
    fr_mass->constPars().printLatex( latexformat, RooFit::OutputFile(plotname + "-massfit-const.tex" ), RooFit::Columns(2) );
    fr_mass->floatParsFinal().printLatex( latexformat, RooFit::OutputFile(plotname + "-massfit-float.tex" ), RooFit::Columns(2) );
  }
  isconverged(fr_mass);

  //  myprintLatex(fr_mass->constPars(), plotname + "-massfit-const.tex" , 2, &latexformat);
  //  myprintLatex(fr_mass->floatParsFinal(), plotname + "-massfit-float.tex" , 2, &latexformat);


  if (doplot and dorand==0) {
    if (doBd2JpsiKS0) c.SetLogy();
    plot(c,RooArgList(mass), tot_mass_model, tot_data, plotname, true, &masspdfs, "", &cat); 
    c.SetLogy(0);
    mygetchar;
  }
  

  
  if (randsfit and dorand>0) {
    const RooArgList& rndPars=fr_mass->randomizePars();
    RooArgSet* blu = tot_mass_model.getVariables();
    TIterator* iter = blu->createIterator();
    RooRealVar* arg;
    while((arg=(RooRealVar*)iter->Next())) {
      TString varname(arg->GetName());
      RooRealVar* rrv = (RooRealVar*)rndPars.find(varname);
      if (rrv) {
	std::cout << varname << " rrv = " << rrv->getVal() << std::endl;
	std::cout << varname << " arg = " << arg->getVal() << std::endl;
	arg->setVal(rrv->getVal());
      }
    }
    delete iter;
  }

  RooArgSet* blu = tot_mass_model.getVariables();
  TIterator* iter = blu->createIterator();
  //  for (int i=0;blu.getSize();i++) std::cout << blu.at(i)->GetName() << std::endl;
  RooRealVar* arg;
  //const everything except nSig, nBkg and _M
  std::vector<RooRealVar*> varstounconst;
  RooArgList yieldlist;
  RooArgList yieldlist_dd;
  RooArgList yieldlist_dd_Lbbar;
  RooArgList yieldlist_ll;
  RooArgList yieldlist_ll_Lbbar;
  //  std::cout << "Making these variables constant: ";
  while((arg=(RooRealVar*)iter->Next())) {
    TString varname(arg->GetName());
    if ((!varname.Contains("_M")) and !varname.Contains("cat") and !varname.Contains("nSig") and !varname.Contains("nBkg") and (!arg->isConstant())) {
    //      std::cout << varname << " ";
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
//  std::cout << std::endl;
  yieldlist.add(yieldlist_ll);
  yieldlist.add(yieldlist_ll_Lbbar);
  yieldlist.add(yieldlist_dd);
  yieldlist.add(yieldlist_dd_Lbbar);
  delete iter;
  //  for (int i=0; i<splotvars.getSize(); i++) ((RooRealVar*)splotvars.at(i))->setConstant();
  //  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
  //					       *tot_data, &tot_mass_model, RooArgList(nSig,nBkg) );
  
  
      //      RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
      //						   *tot_data, &tot_mass_model, yieldlist );

  if (masspdfsmapping) {
    if (doBd2JpsiKS0) {
      if (!dontfitdd) {
	RooStats::SPlot sData1("sData1","An SPlot",
			       *(*datamapping)["dd"], (*masspdfsmapping)["dd"], yieldlist_dd );
	copyvar((*datamapping)["dd"],"nSig_dd_sw",wmass);
	if (debug) {
	  sData1.Print("v");
	  mygetchar;
	}
      }
      if (!dontfitll) {
	RooStats::SPlot sData3("sData1","An SPlot",
			       *(*datamapping)["ll"], (*masspdfsmapping)["ll"], yieldlist_ll );
	copyvar((*datamapping)["ll"],"nSig_ll_sw",wmass);
      }
      if (!dontfitdd) {
	multiplyweights((*datamapping)["dd"],     wtot,wacc.GetName(),wmass.GetName(), normalizetowmass);
	std::cout << "Sum of wtot in tot_dd_data        = " << sumofvarindata((*datamapping)["dd"],"wtot") << std::endl;
      }
      if (!dontfitll) {	
	multiplyweights((*datamapping)["ll"],     wtot,wacc.GetName(),wmass.GetName(), normalizetowmass);
	std::cout << "Sum of wtot in tot_ll_data        = " << sumofvarindata((*datamapping)["ll"],"wtot") << std::endl;      
      }

    } else {
      if (separateLbLbar) {
	if (!dontfitLb) {
	  if (!dontfitdd) {
	    RooStats::SPlot sData1("sData1","An SPlot",
				   *(*datamapping)["ddLb"], (*masspdfsmapping)["ddLb"], yieldlist_dd );
	    copyvar((*datamapping)["ddLb"],"nSig_dd_sw",wmass);
	    if (debug) {
	      sData1.Print("v");
	      mygetchar;
	    }
	  }
	  if (!dontfitll) {
	    if (donotsubtractll) {
	      fillwmassll((*datamapping)["llLb"],wmass);
	    } else {
	      RooStats::SPlot sData3("sData1","An SPlot",
				     *(*datamapping)["llLb"], (*masspdfsmapping)["llLb"], yieldlist_ll );
	      copyvar((*datamapping)["llLb"],"nSig_ll_sw",wmass);
	    }
	  }
	}
	
	if (!dontfitLbbar) {
	  if (!dontfitdd) {
	    RooStats::SPlot sData2("sData1","An SPlot",
				   *(*datamapping)["ddLbbar"], (*masspdfsmapping)["ddLbbar"], yieldlist_dd_Lbbar );
	    copyvar((*datamapping)["ddLbbar"],"nSig_dd_Lbbar_sw",wmass);
	  }
	  if (!dontfitll) {
	    if (donotsubtractll) {
	      fillwmassll((*datamapping)["llLbbar"],wmass);
	    } else {
	      RooStats::SPlot sData4("sData1","An SPlot",
				     *(*datamapping)["llLbbar"], (*masspdfsmapping)["llLbbar"], yieldlist_ll_Lbbar );
	      copyvar((*datamapping)["llLbbar"],"nSig_ll_Lbbar_sw",wmass);
	    }
	  }
	}
	
	if (!dontfitLbbar) {
	  if (!dontfitdd) {
	    multiplyweights((*datamapping)["ddLbbar"],  wtot,wacc.GetName(),wmass.GetName(), normalizetowmass, wcorr.GetName());
	    std::cout << "Sum of wtot in tot_dd_Lbbar_data  = " << sumofvarindata((*datamapping)["ddLbbar"],"wtot") << std::endl;
	  }
	  if (!dontfitll) {
	    multiplyweights((*datamapping)["llLbbar"],  wtot,wacc.GetName(),wmass.GetName(), normalizetowmass, wcorr.GetName());
	    std::cout << "Sum of wtot in tot_ll_Lbbar_data  = " << sumofvarindata((*datamapping)["llLbbar"],"wtot") << std::endl;
	  }
	}
	if (!dontfitLb) {
	  if (!dontfitdd) {
	    multiplyweights((*datamapping)["ddLb"],     wtot,wacc.GetName(),wmass.GetName(), normalizetowmass, wcorr.GetName());
	    std::cout << "Sum of wtot in tot_dd_data        = " << sumofvarindata((*datamapping)["ddLb"],"wtot") << std::endl;
	  }
	  if (!dontfitll) {
	    multiplyweights((*datamapping)["llLb"],     wtot,wacc.GetName(),wmass.GetName(), normalizetowmass, wcorr.GetName());
	    std::cout << "Sum of wtot in tot_ll_data        = " << sumofvarindata((*datamapping)["llLb"],"wtot") << std::endl;
	  }
	}
      } else { 
	if (!dontfitdd) {
	  RooStats::SPlot sData1("sData1","An SPlot",
				 *(*datamapping)["dd"], (*masspdfsmapping)["dd"], yieldlist_dd );
	  copyvar((*datamapping)["dd"],"nSig_dd_sw",wmass);
	  //	  copyvar((*datamapping)["dd"],"L_nSig_dd",wmass);
	  if (dorand==0) {
	    (*datamapping)["dd"]->Print("v");
	    mygetchar;
	  }

	  if (debug) {
	    sData1.Print("v");
	    mygetchar;
	  }
	}
	if (!dontfitll) {
	  if (donotsubtractll) {
	    fillwmassll((*datamapping)["ll"],wmass);
	  } else {
	    RooStats::SPlot sData3("sData1","An SPlot",
				   *(*datamapping)["ll"], (*masspdfsmapping)["ll"], yieldlist_ll );
	    copyvar((*datamapping)["ll"],"nSig_ll_sw",wmass);
	    //	    copyvar((*datamapping)["ll"],"L_nSig_ll",wmass);
	  }
	}
	//  RooStats::SPlot sData("sData","An SPlot",
	//  			*tot_data, &tot_mass_model, yieldlist );
	//    RooDataSet data_sig("data_sig","data_sig",tot_data,*tot_data->get());//,0,"nSig_sw") ;
	
	if (!dontfitdd) {
	  multiplyweights((*datamapping)["dd"],     wtot,wacc.GetName(),wmass.GetName(), normalizetowmass);
	  std::cout << "Sum of wtot in tot_dd_data        = " << sumofvarindata((*datamapping)["dd"],"wtot") << std::endl;
	}
	if (!dontfitll) {
	  multiplyweights((*datamapping)["ll"],     wtot,wacc.GetName(),wmass.GetName(), normalizetowmass);
	  std::cout << "Sum of wtot in tot_ll_data        = " << sumofvarindata((*datamapping)["ll"],"wtot") << std::endl;      
	}
      }
    } 
  } else {
    RooStats::SPlot sData("sData","An SPlot",
			  *tot_data, &tot_mass_model, yieldlist );
    
    copyvar(tot_data,"nSig_ll_Lbbar_sw",wmass);
    addacceptanceweight(tot_data, wacc,doBd2JpsiKS0,false,acceptancell,acceptancedd);
    multiplyweights(tot_data,  wtot,wacc.GetName(),wmass.GetName(), normalizetowmass);
    std::cout << "Sum of wtot in tot_data        = " << sumofvarindata(tot_data,"wtot") << std::endl;
  }
  RooDataSet* tot_data_w=NULL;
  if (masspdfsmapping) {
    //    RooArgList myvars(*(*datamapping)["dd"]->get());
    //    myvars.
    RooArgSet myvars(cat,mass,wmass,wacc,wtot);
    myvars.add(*angvars);
    //    myvars.add(RooArgList(*hrrv[1],*hrrv[2],*hrrv[3],*hrrv[4],*hrrv[5],*hrrv[6],*hrrv[7]));
    tot_data_w = new RooDataSet("simdata2", "simdata2", myvars, RooFit::Index(cat),RooFit::Import(*datamapping), RooFit::WeightVar(wtot));
  } else {
    tot_data_w = new RooDataSet("dataw","dataw",tot_data,*tot_data->get(),0,"wtot");
  }

  //  tot_data_w->Print("v");
  //  if (debug) {
  //  RooDataHist* rdh;
  // if (dorand==0) {
  //   TH3F histo("histo","histo",
  // 	       25,-1.,1.,
  // 	       25,-1.,1.,
  // 	       25,-1.,1.);
    
  //   tot_data_w->fillHistogram(&histo,RooArgList(costheta,costheta1,costheta2));
  //   histo.SaveAs(plotname + "-thf3.root");
  //   //    rdh=new RooDataHist("rdh","rdh",RooArgList(costheta,costheta1,costheta2),&histo);
  //   tot_data_w->write(plotname + "-data.txt");
  // }
    //    tot_data->Print("v");
  //  }
  //  RooDataSet* tot_data_w = new RooDataSet("dataw","dataw",tot_data,*tot_data->get(),0,"wtot");
  
  if (debug) {
    tot_data_w->Print("v");
    mygetchar;
  }
  
  RooCmdArg minimizer=RooCmdArg::none();
  //  RooCmdArg minimizer=RooFit::Minimizer("Minuit2","minimize");
  // does not converge  RooCmdArg minimizer=RooFit::Minimizer("GSLSimAn","-");
  //  RooCmdArg minimizer=RooFit::Minimizer("GSLMultiMin","conjugatefr");
  //  RooCmdArg minimizer=RooFit::Minimizer("GSLMultiMin","conjugatepr");

  //   RooFitResult* fr_woutw = sig_ang_model.fitTo(*tot_data_w,RooFit::SumW2Error(false), RooFit::Save(true), RooFit::Minos(false), minimizer, extconst);
  //   if (fr_woutw->status()!=0 or fr_woutw->covQual()!=3) {
  //     std::cout << plotname << ": status = " << fr_woutw->status() << " ; covQual = " << fr_woutw->covQual() <<  std::endl;
  //     std::cout << "The real fit will fail" << std::endl;
  //     mygetchar;
  //   } else {
  //     std::cout << "The real fit will succeed" << std::endl;
  //     mygetchar;
  //   }

  //  RooFitResult* fr = sig_ang_model.fitTo(*tot_data_w,RooFit::SumW2Error(false), RooFit::Save(true), RooFit::Minos(true), minimizer, extconstangfit);

  const TString plotnameorig(plotname);
  for (int i=0; i<angpdflist.getSize(); i++) {

    RooAbsPdf* sig_ang_model=(RooAbsPdf*)angpdflist.at(i);
    const RooArgSet* obs = sig_ang_model->getObservables(tot_data_w);
    obs->Print("v");
    mygetchar;

    plotname=plotnameorig;
    if (obs->getSize()<3) {
      plotname+="-loop" + istr(i);
    }

    const bool doprefit=false;
    if (doprefit) {
      P_b=0.;
      P_b.setConstant();
      alpha_b=0.;
      alpha_b.setConstant();
      RooFitResult* frtmp1 = myfitto(tot_data_w,sig_ang_model,extconstangfit);
      frtmp1->Print("v");
      if (!doBd2JpsiKS0 and (dorand==0)) printamplitudes(frtmp1,plotname);
      delete frtmp1;
      mygetchar;
      
      alpha_b=-0.50;
      alpha_b.setConstant(false);
      RooFitResult* frtmp2 = myfitto(tot_data_w,sig_ang_model,extconstangfit);
      if (doplot and dorand==0) plot(c,*obs, *sig_ang_model, tot_data_w, plotname, true, NULL, "", &cat);
      frtmp2->Print("v");
      if (!doBd2JpsiKS0 and (dorand==0)) printamplitudes(frtmp2,plotname);
      delete frtmp2;
      mygetchar;
      
      //  alpha_b=0.;
      //  alpha_b.setConstant();
      P_b.setConstant(false);
      alpha_plus.setConstant(false);
      alpha_minus.setConstant(false);
    }
    
    //    alpha_b=0.77;
    //    alpha_b.setConstant();
    //    P_b=0.20;
    //    P_b.setConstant();


    //  RooDataHist rdh("rdh","rdh",RooArgList(costheta,costheta1,costheta2), *tot_data_w);
    //    P_b.setConstant();
    RooFitResult* fr = myfitto(tot_data_w,sig_ang_model,extconstangfit, 1, useminos);
    //  RooFitResult* fr = sig_ang_model.chi2FitTo(*tot_data_w, RooFit::Save(true));
    fr->Print("v");
    std::cout << plotname << ": status = " << fr->status() << " ; covQual = " << fr->covQual() <<  std::endl;
    isconverged(fr);
    if (doplot and dorand==0) {
      //    TCanvas c("c","c",100,100);
      //      plotwpull(c,RooArgList(mass,costheta,costheta1,costheta2), tot_model, dataw, "toy-" + istr(iSample,"04") );
      plot(c,*obs, *sig_ang_model, tot_data_w, plotname, true, NULL, "", &cat);
      
      //      plotlikelihoodprofiles(plotname, c,RooArgList(P_b,alpha_b,r_0,r_1),tot_data_w,sig_ang_model,extconstangfit);
      
      // if (useparam==0) {
      // 	//	plotlikelihoodscans(plotname,c,RooArgList(ap2,bm2),tot_data_w,sig_ang_model,extconstangfit);

      // 	P_b=0.;
      // 	am2=0.;
      // 	bm2=0.;
      // 	plotlikelihoodscans(plotname,c,ap2,0.,0.2,RooArgList(P_b,am2,bm2),tot_data_w,sig_ang_model,extconstangfit);

      // 	P_b=0.;
      // 	ap2=0.;
      // 	bm2=0.;
      // 	plotlikelihoodscans(plotname,c,bm2,0.,0.2,RooArgList(P_b,ap2,am2),tot_data_w,sig_ang_model,extconstangfit);

      // }
      // if (useparam==1) {
      // 	alpha_b=0.;
      // 	r_0=0.5;
      // 	r_1=0.0;
      // 	//	plotlikelihoodscans(plotname,c,P_b, -0.5, 0.5, RooArgList(alpha_b, r_0, r_1), tot_data_w,sig_ang_model,extconstangfit);
      // 	P_b=0.;
      // 	r_0=0.5;
      // 	r_1=0.0;
      // 	plotlikelihoodscans(plotname,c,alpha_b, -0.6, 0.8, RooArgList(P_b, r_0, r_1), tot_data_w,sig_ang_model,extconstangfit);
      // }

      std::cout << "NLL = " << fr->minNll() << std::endl;

      std::cout << "Correlation Matrix" << std::endl;
      fr->correlationMatrix().Print(".2f");

      std::cout << "Covariance Matrix" << std::endl;
      fr->covarianceMatrix().Print(".2f");

      //      plotwcutwcomponents(c,RooArgList(mass,costheta,costheta1,costheta2), tot_model, RooArgList(sig_model,bkg_model),tot_data_w, "signalregion", "toy-" + istr(iSample,"04") + "-sr" );
      
      //    for (int i=1;i<=7;i++) plotdata(c,*hrrv[i], tot_data_w, plotname + "-" + hrrv[i]->GetName() + ".eps");
      
    }
    if (!doBd2JpsiKS0 and (dorand==0)) printamplitudes(fr,plotname);
    
    if (dorand==0) {
      RooArgSet constpars2(fr->constPars());constpars2.writeToFile(plotname + "-angfit-const.txt");
      RooArgSet floatpars2(fr->floatParsFinal());floatpars2.writeToFile(plotname + "-angfit-float.txt");
      fr->constPars().printLatex( latexformat, RooFit::OutputFile(plotname + "-angfit-const.tex" ), RooFit::Columns(1) );
      fr->floatParsFinal().printLatex( latexformat, RooFit::OutputFile(plotname + "-angfit-float.tex" ), RooFit::Columns(1) );
    }
    //  myprintLatex(fr->constPars(), plotname + "-angfit-const.tex" , 2, &latexformat);
    //  myprintLatex(fr->floatParsFinal(), plotname + "-angfit-float.tex" , 2, &latexformat);
    
    if (!doBd2JpsiKS0 and dorand==0) {
      computeamplitudes(ap2, am2, bp2_rrv, bm2, P_b.getVal(), alpha_b.getVal(), r_0.getVal(), r_1.getVal(), fr->covarianceMatrix());
      RooArgSet test(ap2,am2,bp2_rrv,bm2);
      test.writeToFile(plotname + "-angfit-amplitudestest.txt");
      test.printLatex( latexformat, RooFit::OutputFile(plotname + "-angfit-amplitudestest.tex" ), RooFit::Columns(1) );
    }
    
    if (dorand==0) mygetchar;
    
    fr->SetName("AngFit");
    fr_mass->SetName("MassFit");
    
    TFile outputfile(plotname + ".root","recreate");
    fr->Write();
    fr_mass->Write();
    outputfile.Close();
    delete fr;
  }
  //  for (int i=0; i<splotvars.getSize(); i++) ((RooRealVar*)splotvars.at(i))->setConstant(false);
  for (unsigned int i=0; i<varstounconst.size();i++) varstounconst[i]->setConstant(false);
  delete tot_data_w;
  delete fr_mass;
}

// void resettoinitialvalues(const RooArgSet& initialvalues)
// {
//   TIterator* iter = initialvalues.createIterator();
//   RooRealVar* arg;
//   while((arg=(RooRealVar*)iter->Next())) {
//     if (!arg->isConstant()) {
//       const double currvalue = arg->getVal();
//       const double initvalue = initialvalues.getRealValue(arg->GetName());
//       std::cout << arg->GetName() << " " << currvalue << " " << initvalue << std::endl;
//       arg->setVal(initialvalues.getRealValue(arg->GetName()));
//     }
//   }
//   mygetchar;
// }
  

void usage(char* argv0)
{
  std::cerr << argv0 << " fit/merge Lb/B0 dorand(0...1000)" << std::endl;
  exit(2);
}

int main(int argc, char** argv)
{

  // RooRealVar x("x","x",-10.,10.);
  // RooRealVar mean("mean","mean",0.);
  // RooRealVar sigma("sigma","sigma",1.);
  
  // RooGaussian gaus("gaus","gaus",x,mean,sigma);
  // RooArgSet* blu = gaus.getVariables();
  // RooRealVar* xcopy = (RooRealVar*)blu->find("x");
  // RooDataSet* datagen = gaus.generate(RooArgSet(x),20);
  
  // for (int i=0;i<datagen->numEntries();i++) {
  //   const RooArgSet* blu = datagen->get(i);
  //   double xval = blu->getRealValue("x");
  //   xcopy->setVal(xval);
  //   std::cout << "x     = " << xval << std::endl;
  //   std::cout << "gauss = " << gaus.getVal() << std::endl; 
  // }

  // return 1;

  if (argc!=4) {
    usage(argv[0]);
  } 
  const TString oper(argv[1]);
  if (oper!="fit" and oper!="merge") usage(argv[0]);
  const TString mode(argv[2]);
  if (mode!="B0" and mode!="Lb") usage(argv[0]);
  const bool doBd2JpsiKS0=mode.Contains("B0");
   dorand=atoi(argv[3]);
   if (dorand<0) usage(argv[0]);

  // const TString oper("fit");
  // if (oper!="fit" and oper!="merge") usage(argv[0]);
  // const TString mode("Lb");
  // if (mode!="B0" and mode!="Lb") usage(argv[0]);
  // const bool doBd2JpsiKS0=mode.Contains("B0");
  // dorand=0;

  lhcbstyle();
  TCanvas c("c","c",100,100);

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
  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0,Lambdabname,Lambda0name,Jpsiname);


  ConfigFile cf("mctestconfig.txt");
  userwmc=atob(cf.Value("main","userwmc","true"));
  usesimplecut=atob(cf.Value("main","usesimplecut","false"));
  usebdt=atob(cf.Value("main","usebdt","true"));
  floatalphalambda=atob(cf.Value("main","floatalphalambda","true"));
  floatSigParams=atob(cf.Value("main","floatsigparams","true"));
  dontfitLb=atob(cf.Value("main","dontfitLb","false"));
  dontfitLbbar=atob(cf.Value("main","dontfitLbbar","false"));
  dontfitdd=atob(cf.Value("main","dontfitdd","false"));
  dontfitll=atob(cf.Value("main","dontfitll","false"));
  separateLbLbar=atob(cf.Value("main","separateLbLbar","false"));
  useparam=atoi(cf.Value("main","useparam","0").c_str());
  randacceptance=atob(cf.Value("main","randacceptance","false"));
  randmccalib=atob(cf.Value("main","randmccalib","false"));
  randsfit=atob(cf.Value("main","randsfit","false"));
  randsigparams=atob(cf.Value("main","randsigparams","false"));
  randalphalambda=atob(cf.Value("main","randalphalambda","false"));
  randbeam=atob(cf.Value("main","randbeam","false"));
  randphases=atob(cf.Value("main","randphases","false"));
  usekeyspdf=atob(cf.Value("main","usekeyspdf","false"));
  donotsubtractll=atob(cf.Value("main","donotsubtractll","false"));
  useminos=atob(cf.Value("main","useminos","true"));

  const int selectLb=atoi(cf.Value("main","selectLb","0").c_str());
  // 0 = all, 1=up only, 2=down only
  const int fitPolarity=atoi(cf.Value("main","Polarity","0").c_str());

  std::cout << "Config" << std::endl;
  std::cout << "userwmc          = " << btoa(userwmc) << std::endl;
  std::cout << "usesimplecut     = " << btoa(usesimplecut) << std::endl;
  std::cout << "usebdt           = " << btoa(usebdt) << std::endl;
  std::cout << "floatalphalambda = " << btoa(floatalphalambda) << std::endl;
  std::cout << "floatSigParams   = " << btoa(floatSigParams) << std::endl;
  std::cout << "dontfitLb        = " << btoa(dontfitLb) << std::endl;
  std::cout << "dontfitLbbar     = " << btoa(dontfitLbbar) << std::endl;
  std::cout << "dontfitdd        = " << btoa(dontfitdd) << std::endl;
  std::cout << "dontfitll        = " << btoa(dontfitll) << std::endl;
  std::cout << "separateLbLbar   = " << btoa(separateLbLbar) << std::endl;
  std::cout << "useparam         = " << useparam << std::endl;
  std::cout << "randacceptance   = " << btoa(randacceptance) << std::endl;
  std::cout << "randmccalib      = " << btoa(randmccalib) << std::endl;
  std::cout << "randsfit         = " << btoa(randsfit) << std::endl;
  std::cout << "randsigparams    = " << btoa(randsigparams) << std::endl;
  std::cout << "randalphalambda  = " << btoa(randalphalambda) << std::endl;
  std::cout << "randbeam         = " << btoa(randbeam) << std::endl;
  std::cout << "usekeyspdf       = " << btoa(usekeyspdf) << std::endl;
  std::cout << "selectLb         = " << selectLb << std::endl;
  std::cout << "fitPolarity      = " << fitPolarity << std::endl;


  mass.SetName(Lambdabname + "_M");
  mass.setUnit("MeV/c^{2}");

  if (doBd2JpsiKS0) {
    mass.setMin(5230.);
    //    mass.setMin(5150.);
    mass.setMax(5450.);
    mass.SetTitle("M(J/#psiK_{S}^{0})");
  } else {
    //    mass.setMin(5500.);
    //    mass.setMax(5750.);
    mass.setMin(5550.);
    mass.setMax(5700.);

  }

  //  angvars=new RooArgList(costheta,costheta1,costheta2,phi1,phi2);
  angvars=new RooArgList(costheta,costheta1,costheta2);

  // P_b = -0.0120926 ;// +/-  (-0.019131, 0.019327) L(-0.5 - 0.5) 
  // alpha_b = -0.0292658 ;//  +/-  (-0.076503, 0.075983) L(-2 - 2) 
  // r_0 = -0.486078  ;// +/-  (-0.032089, 0.032402) L(-2 - 2) 
  // r_1 = -0.615183  ;// +/-  (-0.041266, 0.041439) L(-2 - 2) 
  //  P_b =  0.063513 ;//  +/- 0.058984 L(-1.5 - 1.5) 
  //  alpha_b = -0.00747469 ;//  +/- 0.098906 L(-1.1 - 1.1) 
  //  r_0 =  0.57440 ;//  +/- 0.022483 L(-0.3 - 2) 
  //  r_1 = -0.603807 ;//  +/- 0.054908 L(-1.3 - 1.3) 

  // P_b.setConstant();
  // alpha_b.setConstant();
  // r_0.setConstant();
  // r_1.setConstant();


  const TCut isLambdab0(Lambdabname + "_ID>0.");
  const TCut isLambdab0bar(Lambdabname + "_ID<0.");

  TCut LL_selection;
  TCut DD_selection;
  if (usesimplecut) simplecut2(doBd2JpsiKS0, LL_selection,DD_selection);

  if (usebdt and usetmvarectcut) {
    std::cerr << "Both bdt and tmvarectcut are set ... " << std::endl;
    return  1;
  }
  if (usebdt) bdtsel(doBd2JpsiKS0, LL_selection,DD_selection);
  if (usetmvarectcut) tmvarectcut(doBd2JpsiKS0, LL_selection,DD_selection);


  LL_selection+=masscut(doBd2JpsiKS0, mass);
  DD_selection+=masscut(doBd2JpsiKS0, mass);

  LL_selection+=triggercut(doBd2JpsiKS0);
  DD_selection+=triggercut(doBd2JpsiKS0);

  //  if (usePID) {
  //    LL_selection+=pidcut(doBd2JpsiKS0);
  //    DD_selection+=pidcut(doBd2JpsiKS0);
  //  }

  //  LL_selection+=TCut("lambda0_PT>3000.");
  //  DD_selection+=TCut("lambda0_PT>3000.");

  if (fitPolarity==1) {
    LL_selection+=MagnetUp;
    DD_selection+=MagnetUp;
  }
  if (fitPolarity==2) {
    LL_selection+=MagnetDown;
    DD_selection+=MagnetDown;
  }

  if (selectLb==1) {
    LL_selection+=isLambdab0;
    DD_selection+=isLambdab0;
  }
  if (selectLb==2) {
    LL_selection+=isLambdab0bar;
    DD_selection+=isLambdab0bar;
  }

  std::cout << "LL_selection: " << LL_selection << std::endl;
  std::cout << "DD_selection: " << DD_selection << std::endl;
  const TCut selection = LL_selection or DD_selection;
  if (dorand==0) mygetchar;

  std::vector<TString> catlist;
  if (doBd2JpsiKS0) {
    if (!dontfitdd) catlist.push_back("dd");
    if (!dontfitll) catlist.push_back("ll");
    for (unsigned int i=0;i<catlist.size();i++) cat.defineType(catlist[i],i);
  } else {
    if (separateLbLbar) {
      if (!dontfitLb    and !dontfitdd) catlist.push_back("ddLb"   );
      if (!dontfitLbbar and !dontfitdd) catlist.push_back("ddLbbar");
      if (!dontfitLb    and !dontfitll) catlist.push_back("llLb"   );
      if (!dontfitLbbar and !dontfitll) catlist.push_back("llLbbar");
    } else {
      if (!dontfitdd) catlist.push_back("dd");
      if (!dontfitll) catlist.push_back("ll");
    }
    for (unsigned int i=0;i<catlist.size();i++) cat.defineType(catlist[i],i);
  }


  RooRealVar bkg_dd_c0("bkg_dd_c0","bkg_dd_c0",0.,-1.1,1.1);
  RooRealVar bkg_ll_c0("bkg_ll_c0","bkg_ll_c0",0.,-1.1,1.1);

  RooRealVar bkg_dd_Lbbar_c0("bkg_dd_Lbbar_c0","bkg_dd_Lbbar_c0",0.,-1.1,1.1);
  RooRealVar bkg_ll_Lbbar_c0("bkg_ll_Lbbar_c0","bkg_ll_Lbbar_c0",0.,-1.1,1.1);

  const bool linearbkg=true;
  RooAbsPdf *bkg_dd_mass_model, *bkg_ll_mass_model, *bkg_dd_Lbbar_mass_model, *bkg_ll_Lbbar_mass_model;

  if (linearbkg) {
    bkg_dd_mass_model = new RooChebychev("bkg_dd_mass_model","bkg_dd_mass_model",mass,RooArgList(bkg_dd_c0));
    bkg_ll_mass_model = new RooChebychev("bkg_ll_mass_model","bkg_ll_mass_model",mass,RooArgList(bkg_ll_c0));
    bkg_dd_Lbbar_mass_model = new RooChebychev("bkg_dd_Lbbar_mass_model","bkg_dd_Lbbar_mass_model",mass,RooArgList(bkg_dd_Lbbar_c0));
    bkg_ll_Lbbar_mass_model = new RooChebychev("bkg_ll_Lbbar_mass_model","bkg_ll_Lbbar_mass_model",mass,RooArgList(bkg_ll_Lbbar_c0));
  } else {
    bkg_dd_mass_model = new RooExponential("bkg_dd_mass_model","bkg_dd_mass_model",mass,bkg_dd_c0);
    bkg_ll_mass_model = new RooExponential("bkg_ll_mass_model","bkg_ll_mass_model",mass,bkg_ll_c0);
    bkg_dd_Lbbar_mass_model = new RooExponential("bkg_dd_Lbbar_mass_model","bkg_dd_Lbbar_mass_model",mass,bkg_dd_Lbbar_c0);
    bkg_ll_Lbbar_mass_model = new RooExponential("bkg_ll_Lbbar_mass_model","bkg_ll_Lbbar_mass_model",mass,bkg_ll_Lbbar_c0);

    RooArgList tmp(bkg_dd_c0, bkg_ll_c0, bkg_dd_Lbbar_c0, bkg_ll_Lbbar_c0);
    for (int i=0;i<tmp.getSize();i++) {
      RooRealVar* var = (RooRealVar*)tmp.at(i);
      var->setVal(1.);
      var->setMin(0.);
      var->setMax(100000.);
    }

  }

  RooRealVar sig_dd_mean("sig_dd_mean","sig_dd_mean",5620.,5500.,5700.);
  RooRealVar sig_ll_mean("sig_ll_mean","sig_ll_mean",5621.,5500.,5700.);
  if (doBd2JpsiKS0) {
    sig_dd_mean.setVal(5280.);
    sig_dd_mean.setMin(5270.);
    sig_dd_mean.setMax(5290.);
    sig_ll_mean.setVal(5280.);
    sig_ll_mean.setMin(5270.);
    sig_ll_mean.setMax(5290.);
  }
  RooRealVar sig_dd_frac("sig_dd_frac","sig_dd_frac", 0.5, -0.1, 1.1);
  RooRealVar sig_ll_frac("sig_ll_frac","sig_ll_frac", 0.5, -0.1, 1.1);


  RooRealVar sig_dd_sigma("sig_dd_sigma","sig_dd_sigma",7.5,0.,100.);
  RooRealVar sig_ll_sigma("sig_ll_sigma","sig_ll_sigma",7.2597,0.,100.);


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

  
  RooAbsPdf& sig_dd_mass_model = sig_dd_twocb;
  RooAbsPdf& sig_ll_mass_model = sig_ll_twocb;

  RooFormulaVar bkgbs_dd_mean("bkgbs_dd_mean","bkgbs_dd_mean","@0+87.35",RooArgList(sig_dd_mean));
  RooFormulaVar bkgbs_ll_mean("bkgbs_ll_mean","bkgbs_ll_mean","@0+87.35",RooArgList(sig_ll_mean));

  RooCBShape bkgbs_dd_cbshape1("bkgbs_dd_cbshape1", "bkgbs_dd_cbshape1", mass, bkgbs_dd_mean, sig_dd_sigma, sig_dd_alpha1, sig_dd_n1);//, RooAbsReal& _n)
  RooCBShape bkgbs_dd_cbshape2("bkgbs_dd_cbshape2", "bkgbs_dd_cbshape2", mass, bkgbs_dd_mean, sig_dd_sigma, sig_dd_alpha2, sig_dd_n2);//, RooAbsReal& _n)
  RooAddPdf bkgbs_dd_twocb("bkgbs_dd_twocb","bkgbs_dd_twocb",RooArgList(bkgbs_dd_cbshape1,bkgbs_dd_cbshape2),RooArgList(sig_dd_frac));
  RooCBShape bkgbs_ll_cbshape1("bkgbs_ll_cbshape1", "bkgbs_ll_cbshape1", mass, bkgbs_ll_mean, sig_ll_sigma, sig_ll_alpha1, sig_ll_n1);//, RooAbsReal& _n)
  RooCBShape bkgbs_ll_cbshape2("bkgbs_ll_cbshape2", "bkgbs_ll_cbshape2", mass, bkgbs_ll_mean, sig_ll_sigma, sig_ll_alpha2, sig_ll_n2);//, RooAbsReal& _n)
  RooAddPdf bkgbs_ll_twocb("bkgbs_ll_twocb","bkgbs_ll_twocb",RooArgList(bkgbs_ll_cbshape1,bkgbs_ll_cbshape2),RooArgList(sig_ll_frac));

  RooAbsPdf& bkgbs_dd_mass_model = bkgbs_dd_twocb;
  RooAbsPdf& bkgbs_ll_mass_model = bkgbs_ll_twocb;


  RooRealVar nSig_dd("nSig_dd","nSig_dd",2000.,0.,300000.);
  RooRealVar nBkg_dd("nBkg_dd","nBkg_dd",10000.,0.,1000000.);
  RooRealVar nSig_ll("nSig_ll","nSig_ll",2000.,0.,300000.);
  RooRealVar nBkg_ll("nBkg_ll","nBkg_ll",10000.,0.,1000000.);
  RooRealVar nSig_dd_Lbbar("nSig_dd_Lbbar","nSig_dd_Lbbar",2000.,0.,300000.);
  RooRealVar nBkg_dd_Lbbar("nBkg_dd_Lbbar","nBkg_dd_Lbbar",10000.,0.,1000000.);
  RooRealVar nSig_ll_Lbbar("nSig_ll_Lbbar","nSig_ll_Lbbar",2000.,0.,300000.);
  RooRealVar nBkg_ll_Lbbar("nBkg_ll_Lbbar","nBkg_ll_Lbbar",10000.,0.,1000000.);
  
  RooRealVar nBkgBs_ll("nBkgBs_ll","nBkgBs_ll",100.,0.,100000.);
  RooRealVar nBkgBs_dd("nBkgBs_dd","nBkgBs_dd",100.,0.,100000.);


  if (doBd2JpsiKS0) {
    
    bkg_dd_c0 = -0.133033;
    bkg_ll_c0 = -0.165032;
    nBkgBs_dd =  94.604;
    nBkgBs_ll =  96.681;
    nBkg_dd =  15721;
    nBkg_ll =  2507.1;
    nSig_dd =  24671;
    nSig_ll =  10353;
  } else {
    bkg_dd_Lbbar_c0 = -0.309009;
    bkg_dd_c0 = -0.245135;
    bkg_ll_Lbbar_c0 = -0.161986;
    bkg_ll_c0 = -0.21779;
    nBkg_dd =  4383.8;
    nBkg_dd_Lbbar =  4111.5;
    nBkg_ll =  526.73;
    nBkg_ll_Lbbar =  476.85;


    nBkg_dd = 0.;
    nBkg_dd_Lbbar = 0.;
    nBkg_ll =  0.;
    nBkg_ll_Lbbar = 0.;

    nSig_dd =  2791.2;
    nSig_dd_Lbbar =  2607.5;
    nSig_ll =  1101.3;
    nSig_ll_Lbbar =  949.16;
  }

  
  RooAbsPdf* tot_dd_mass_model;
  RooAbsPdf* tot_ll_mass_model;
  if (doBd2JpsiKS0) {
    tot_dd_mass_model = new RooAddPdf("tot_dd_mass_model","tot_dd_mass_model",RooArgList(sig_dd_mass_model,*bkg_dd_mass_model,bkgbs_dd_mass_model),RooArgList(nSig_dd,nBkg_dd,nBkgBs_dd));
    tot_ll_mass_model = new RooAddPdf("tot_ll_mass_model","tot_ll_mass_model",RooArgList(sig_ll_mass_model,*bkg_ll_mass_model,bkgbs_ll_mass_model),RooArgList(nSig_ll,nBkg_ll,nBkgBs_ll));
  } else {
    if (extendedmassfit) {
      tot_dd_mass_model = new RooAddPdf("tot_dd_mass_model","tot_dd_mass_model",RooArgList(sig_dd_mass_model,*bkg_dd_mass_model),RooArgList(nSig_dd,nBkg_dd));
      tot_ll_mass_model = new RooAddPdf("tot_ll_mass_model","tot_ll_mass_model",RooArgList(sig_ll_mass_model,*bkg_ll_mass_model),RooArgList(nSig_ll,nBkg_ll));
    } else {
      tot_dd_mass_model = new RooAddPdf("tot_dd_mass_model","tot_dd_mass_model",RooArgList(sig_dd_mass_model,*bkg_dd_mass_model),RooArgList(nSig_dd));
      tot_ll_mass_model = new RooAddPdf("tot_ll_mass_model","tot_ll_mass_model",RooArgList(sig_ll_mass_model,*bkg_ll_mass_model),RooArgList(nSig_ll));
      nSig_dd=0.5;
      nSig_ll=0.5;
    }
  }

  RooAddPdf* tot_dd_Lbbar_mass_model = new RooAddPdf("tot_dd_Lbbar_mass_model","tot_dd_Lbbar_mass_model",RooArgList(sig_dd_mass_model,*bkg_dd_Lbbar_mass_model),RooArgList(nSig_dd_Lbbar,nBkg_dd_Lbbar));
  RooAddPdf* tot_ll_Lbbar_mass_model = new RooAddPdf("tot_ll_Lbbar_mass_model","tot_ll_Lbbar_mass_model",RooArgList(sig_ll_mass_model,*bkg_ll_Lbbar_mass_model),RooArgList(nSig_ll_Lbbar,nBkg_ll_Lbbar));

//   RooArgSet* blu = tot_mass_model.getVariables();
//   TIterator* iter = blu->createIterator();
//   //  for (int i=0;blu.getSize();i++) std::cout << blu.at(i)->GetName() << std::endl;
//   RooAbsArg* arg;
//   while((arg=(RooAbsArg*)iter->Next())) {
//     std::cout << arg->GetName() << std::endl;
//   }
//   return 1;
  

  // RooRealVar& P_b_unblind = P_b;
  // if (blindPb) {
  //   P_b_unblind = *new RooUnblindPrecision("P_b_unblind","P_{b} (unblind)","TheBlindingString",-1.5,1.5, P_b) ;
  // }
  //  P_b=0.;
  //  P_b.setConstant();

  //  alpha_b=0.;
  // alpha_b.setConstant();

  const double al_val=0.0;
  const double al_err=0.013;
  const double albar_val=0.71;
  const double albar_err=0.08;

  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}",
			  al_val,
			  al_val-7.*al_err,  
			  al_val+7.*al_err);
  TRandom3 rnd(0);
  RooRandom::randomGenerator()->SetSeed(0);
  if (randalphalambda and dorand>0) alpha_lambda=rnd.Gaus(al_val,al_err);

  RooRealVar& alpha_lambdabar=alpha_lambda;
  if (!assumecp) {
    alpha_lambdabar = *new RooRealVar("alpha_lambdabar","#alpha_{#bar{#lambda}}",
				      albar_val, 
				      albar_val - 7.*albar_err,  
				      albar_val + 7.*albar_err);
  }
  RooGaussian alpha_lambda_constrain("alpha_lambda_constrain",
				     "alpha_lambda_constrain",
				     alpha_lambda,
				     RooFit::RooConst(al_val),
				     RooFit::RooConst(al_err));
  RooGaussian alpha_lambdabar_constrain("alpha_lambdabar_constrain",
					"alpha_lambdabar_constrain",
					alpha_lambdabar,
					RooFit::RooConst(albar_val),
					RooFit::RooConst(albar_err));
  
  RooFormulaVar alpha_b_rfv("alpha_b_rfv","ap2-am2+bp2-bm2",RooArgList(ap2,am2,bp2,bm2));
  RooFormulaVar r_0_rfv("r_0_rfv","ap2+am2",RooArgList(ap2,am2));
  RooFormulaVar r_1_rfv("r_1_rfv","ap2-am2",RooArgList(ap2,am2));


  computew3 w3;
  if (isnan(w3.w(0.,0.,0.,P_b.getVal(),alpha_b.getVal(),alpha_lambda.getVal(),r_0.getVal(),r_1.getVal()))) {
    std::cerr << "Ca va pas marcher" << std::endl;
    return 1;
  }



  if (doBd2JpsiKS0) {
    alpha_lambda=0.;
    alpha_lambdabar=0.;
    alpha_lambda.setConstant();
    alpha_lambdabar.setConstant();
  }

  RooAbsPdf* sig_ang_model       = NULL;
  RooAbsPdf* sig_ang_model_Lbbar = NULL;
  //  RooLbtoJpsiL0PDF3* sig_ang_model       = NULL;
  //  RooLbtoJpsiL0PDF3* sig_ang_model_Lbbar = NULL;
  RooFormulaVar r_1_alpha_b("r_1_alpha_b","alpha_b/3.",RooArgSet(alpha_b));
  //  if (doBd2JpsiKS0) {
  //    sig_ang_model       = new RooLbtoJpsiL0PDF3("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1_alpha_b);
  //    sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF3("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambdabar, r_0, r_1_alpha_b);
  //  } else {
  //  sig_ang_model       = new RooLbtoJpsiL0PDF3("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1);
  //  sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF3("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambdabar, r_0, r_1);
  //  }
  
  //  RooRealVar *acccoeffll[ordermax3i][ordermax3j][ordermax3k];  
  //  RooRealVar *acccoeffdd[ordermax3i][ordermax3j][ordermax3k]; 
  
  if (useparam==0) {
    sig_ang_model = new RooLbtoJpsiL0PDF3v3("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, ap2, alpha_lambda, am2, bp2);
    sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF3v3("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,P_b, ap2, alpha_lambdabar, am2, bp2);

    ap2.setMin(-0.3);
    ap2.setMax(1.3);
    am2.setMin(-0.3);
    am2.setMax(1.3);
    bm2.setMin(-0.3);
    bm2.setMax(1.3);

    ap2.setConstant(false);
    am2.setConstant(false);
    bm2.setConstant(false);

    ap2_rfv = new RooFormulaVar("ap2_rrv","|a+|^{2}","@0",RooArgList(ap2));
    am2_rfv = new RooFormulaVar("am2_rrv","|a-|^{2}","@0",RooArgList(am2));
    bp2_rfv = new RooFormulaVar("bp2_rrv","|b+|^{2}","@0",RooArgList(bp2));
    bm2_rfv = new RooFormulaVar("bm2_rrv","|b-|^{2}","1.-@0-@1-@2",RooArgList(ap2,am2,bp2));


  } 
  if (useparam==1)  {
    
    //    sig_ang_model = new RooLbtoJpsiL0PDF5("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,phi1,phi2,P_b, alpha_b, alpha_lambda, r_0, r_1, alpha_plus, alpha_minus, chi);
    //    sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF5("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,phi1,phi2,P_b, alpha_b, alpha_lambdabar, r_0, r_1, alpha_plus, alpha_minus, chi);

    if (randphases) {

      RooRealVar *coeff2[ordermax2i][ordermax2j];
      createRRV2(doBd2JpsiKS0,coeff2, "Lbllc0c1c2",  9, 9, 10, 1);
      //      readcoeff2_rrv("fitaccphi1phi2-mcrw-wbdt-Lb-Tuple_JpsiL0_betas-ll.root", coeff2);

      RooArgList* varlist = new RooArgList();
      for (int i=0; i<ordermax2i;i++) {
	for (int j=0; j<ordermax2j;j++) {
	  if (coeff2[i][j]) {
	    coeff2[i][j]->setConstant();
	    varlist->add(*coeff2[i][j]);
	  }
	}
      }

      readvarsfromtfile(*varlist, "fitaccphi1phi2-mcrw-wbdt-Lb-Tuple_JpsiL0_betas-ll.root", "", randphases and dorand!=0);
      for (int i=0;i<varlist->getSize();i++) ((RooRealVar*)varlist->at(i))->setConstant();


      alpha_plus=2.*pi*rnd.Rndm();
      alpha_minus=2.*pi*rnd.Rndm();
      chi=2.*pi*rnd.Rndm();

      // alpha_plus=0.;
      // alpha_minus=0.;
      // chi=0.;
      alpha_plus.setConstant();
      alpha_minus.setConstant();
      chi.setConstant();

      sig_ang_model = new RooLbtoJpsiL0PDF3v4("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1, alpha_plus, alpha_minus, chi, *varlist);
      sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF3v4("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambdabar, r_0, r_1, alpha_plus, alpha_minus, chi, *varlist);


    } else {
      sig_ang_model = new RooLbtoJpsiL0PDF3("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1);
      sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF3("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambdabar, r_0, r_1);
    }

    ap2_rfv = new RooFormulaVar("ap2_rrv","|a+|^{2}","0.5*(@0+@1)",RooArgList(r_0,r_1));
    am2_rfv = new RooFormulaVar("am2_rrv","|a-|^{2}","0.5*(@0-@1)",RooArgList(r_0,r_1));
    bp2_rfv = new RooFormulaVar("bp2_rrv","|b+|^{2}","0.5*(1.+@0-@1-@2)",RooArgList(alpha_b,r_0,r_1));
    bm2_rfv = new RooFormulaVar("bm2_rrv","|b-|^{2}","0.5*(1.-@0-@1+@2)",RooArgList(alpha_b,r_0,r_1));
    
  } 
  if (useparam==2) {
    sig_ang_model = new RooLbtoJpsiL0PDF3v2("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, beta, gammarrv);
    sig_ang_model_Lbbar = new RooLbtoJpsiL0PDF3v2("sig_ang_model_Lbbar","sig_ang_model_Lbbar",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambdabar, beta, gammarrv);

    ap2_rfv = new RooFormulaVar("ap2_rrv","|a+|^{2}","1./12.*(+3.*@0+3.*@1-4.*@2+2.)",RooArgList(alpha_b,beta,gammarrv));
    am2_rfv = new RooFormulaVar("am2_rrv","|a-|^{2}","1./12.*(-3.*@0-3.*@1-4.*@2+2.)",RooArgList(alpha_b,beta,gammarrv));
    bp2_rfv = new RooFormulaVar("bp2_rrv","|b+|^{2}","1./12.*(+3.*@0-3.*@1+4.*@2+4.)",RooArgList(alpha_b,beta,gammarrv));
    bm2_rfv = new RooFormulaVar("bm2_rrv","|b-|^{2}","1./12.*(-3.*@0+3.*@1+4.*@2+4.)",RooArgList(alpha_b,beta,gammarrv));

  }
  
    //    createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b,alpha_lambda, r_0, r_1, "Lbllc0c1c2", sig_ang_model, acccoeffll, true,acceptancefile + "-ll.root");


    //  0  1. + 
    //  1  alpha_b * P_b * costheta + 
    //  2  (2.*r_1-alpha_b) * alpha_lambda * costheta1 +
    //  3  (2.*r_0 - 1.) * P_b * alpha_lambda * costheta * costheta1 +
    //  4  (0.5-1.5*r_0) * 0.5 * (3.*costheta2*costheta2 - 1.) +
    //  5  (0.5 * alpha_b  - 1.5 * r_1) * P_b * 0.5 * (3.*costheta2*costheta2 - 1.) * costheta +
    //  6  (-0.5 * alpha_b  - 0.5 * r_1) * alpha_lambda * 0.5 * (3.*costheta2*costheta2 - 1.) * costheta1 +
    //  7  (-0.5 * r_0  - 0.5) * P_b * alpha_lambda * 0.5 * (3.*costheta2*costheta2 - 1.) * costheta * costheta1

    //  0  1. + 
    //  1  @4                      * @3      * @0 + 
    //  2  (2.*@7-@4)              * @5      * @1 +
    //  3  (2.*@6 - 1.)            * @3 * @5 * @0 * @1 +
    //  4  (0.5-1.5*@6)            * 1.0     * 0.5 * (3.*@2*@2 - 1.) +
    //  5  (0.5 * @4  - 1.5 * @7)  * @3      * 0.5 * (3.*@2*@2 - 1.) * @0 +
    //  6  (-0.5 * @4  - 0.5 * @7) * @5      * 0.5 * (3.*@2*@2 - 1.) * @1 +
    //  7  (-0.5 * @6  - 0.5)      * @3 * @5 * 0.5 * (3.*@2*@2 - 1.) * @0 * @1

    // sig_ang_model = new RooGenericPdf("sig_ang_model","sig_ang_model", " 1. + @4 * @3 * @0 + (2.*@7-@4) * @5 * @1 + (2.*@6 - 1.) * @3 * @5 * @0 * @1 + (0.5-1.5*@6) * 0.5 * (3.*@2*@2 - 1.) + (0.5 * @4  - 1.5 * @7) * @3 * 0.5 * (3.*@2*@2 - 1.) * @0 + (-0.5 * @4  - 0.5 * @7) * @5 * 0.5 * (3.*@2*@2 - 1.) * @1 + (-0.5 * @6  - 0.5) * @3 * @5 * 0.5 * (3.*@2*@2 - 1.) * @0 * @1" ,RooArgList(costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1));

    // sig_ang_model_Lbbar = new RooGenericPdf("sig_ang_model_Lbbar","sig_ang_model_Lbbar", " 1. + @4 * @3 * @0 + (2.*@7-@4) * @5 * @1 + (2.*@6 - 1.) * @3 * @5 * @0 * @1 + (0.5-1.5*@6) * 0.5 * (3.*@2*@2 - 1.) + (0.5 * @4  - 1.5 * @7) * @3 * 0.5 * (3.*@2*@2 - 1.) * @0 + (-0.5 * @4  - 0.5 * @7) * @5 * 0.5 * (3.*@2*@2 - 1.) * @1 + (-0.5 * @6  - 0.5) * @3 * @5 * 0.5 * (3.*@2*@2 - 1.) * @0 * @1" ,RooArgList(costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambdabar, r_0, r_1));

  
  const double maxchebval=2.0;
  RooRealVar bkg_dd_costheta_c0("bkg_dd_costheta0_c0","bkg_dd_costheta0_c0",0.,-maxchebval,maxchebval);
  //  bkg_dd_costheta_c0.setConstant();
  RooRealVar bkg_dd_costheta1_c0("bkg_dd_costheta1_c0","bkg_dd_costheta1_c0",0.,-maxchebval,maxchebval);
  //  bkg_dd_costheta1_c0.setConstant();
  RooRealVar bkg_dd_costheta2_c0("bkg_dd_costheta2_c0","bkg_dd_costheta2_c0",0.,-maxchebval,maxchebval);
  RooRealVar bkg_dd_costheta2_c1("bkg_dd_costheta2_c1","bkg_dd_costheta2_c1",-0.9,-maxchebval,maxchebval);
  RooRealVar bkg_dd_costheta2_c2("bkg_dd_costheta2_c2","bkg_dd_costheta2_c2",0.,-maxchebval,maxchebval);
  //  bkg_dd_costheta2_c1=0.; bkg_dd_costheta2_c1.setConstant();
  //  bkg_dd_costheta2_c0.setConstant();
  bkg_dd_costheta2_c2.setConstant();

  RooRealVar bkg_ll_costheta_c0("bkg_ll_costheta0_c0","bkg_ll_costheta0_c0",0.,-maxchebval,maxchebval);
  //  bkg_ll_costheta_c0.setConstant();
  RooRealVar bkg_ll_costheta1_c0("bkg_ll_costheta1_c0","bkg_ll_costheta1_c0",0.,-maxchebval,maxchebval);
  //  bkg_ll_costheta1_c0.setConstant();
  RooRealVar bkg_ll_costheta2_c0("bkg_ll_costheta2_c0","bkg_ll_costheta2_c0",0.,-maxchebval,maxchebval);
  RooRealVar bkg_ll_costheta2_c1("bkg_ll_costheta2_c1","bkg_ll_costheta2_c1",-0.9,-maxchebval,maxchebval);
  RooRealVar bkg_ll_costheta2_c2("bkg_ll_costheta2_c2","bkg_ll_costheta2_c2",0.,-maxchebval,maxchebval);
  //  bkg_ll_costheta2_c1=0.; bkg_ll_costheta2_c1.setConstant();
  //  bkg_ll_costheta2_c0.setConstant();
  bkg_ll_costheta2_c2.setConstant();

  RooAbsPdf *bkg_dd_costheta_model,*bkg_dd_costheta1_model,*bkg_dd_costheta2_model;
  RooAbsPdf *bkg_ll_costheta_model,*bkg_ll_costheta1_model,*bkg_ll_costheta2_model;

  const bool usecheb=true;
  if (usecheb) {
    bkg_dd_costheta_model = new RooChebychev("bkg_dd_costheta_model",  "bkg_dd_costheta_model", costheta, RooArgList(bkg_dd_costheta_c0));//, bkg_dd_costheta_c1, bkg_dd_costheta_c2));
    bkg_dd_costheta1_model = new RooChebychev("bkg_dd_costheta1_model","bkg_dd_costheta1_model",costheta1,RooArgList(bkg_dd_costheta1_c0));//,bkg_dd_costheta1_c1,bkg_dd_costheta1_c2));

    bkg_dd_costheta2_model = new RooChebychev("bkg_dd_costheta2_model","bkg_dd_costheta2_model",costheta2,RooArgList(bkg_dd_costheta2_c0,bkg_dd_costheta2_c1,bkg_dd_costheta2_c2));

    //    bkg_dd_costheta2_model = new RooGenericPdf("bkg_dd_costheta2_model","bkg_dd_costheta2_model","1-@0*@0*@1",RooArgList(costheta2,bkg_dd_costheta2_c1));

    bkg_ll_costheta_model = new RooChebychev("bkg_ll_costheta_model",  "bkg_ll_costheta_model", costheta, RooArgList(bkg_ll_costheta_c0));//, bkg_ll_costheta_c1, bkg_ll_costheta_c2));
    bkg_ll_costheta1_model = new RooChebychev("bkg_ll_costheta1_model","bkg_ll_costheta1_model",costheta1,RooArgList(bkg_ll_costheta1_c0));//,bkg_ll_costheta1_c1,bkg_ll_costheta1_c2));
    bkg_ll_costheta2_model = new RooChebychev("bkg_ll_costheta2_model","bkg_ll_costheta2_model",costheta2,RooArgList(bkg_ll_costheta2_c0,bkg_ll_costheta2_c1,bkg_ll_costheta2_c2));

    //    bkg_ll_costheta2_model = new RooGenericPdf("bkg_ll_costheta2_model","bkg_ll_costheta2_model","1-@0*@0*@1",RooArgList(costheta2,bkg_ll_costheta2_c1));
  } else {
    bkg_dd_costheta_model = new RooPolynomial("bkg_dd_costheta_model",  "bkg_dd_costheta_model", costheta, RooArgList(bkg_dd_costheta_c0));//, bkg_dd_costheta_c1, bkg_dd_costheta_c2));
    bkg_dd_costheta1_model = new RooPolynomial("bkg_dd_costheta1_model","bkg_dd_costheta1_model",costheta1,RooArgList(bkg_dd_costheta1_c0));//,bkg_dd_costheta1_c1,bkg_dd_costheta1_c2));
    bkg_dd_costheta2_model = new RooPolynomial("bkg_dd_costheta2_model","bkg_dd_costheta2_model",costheta2,RooArgList(bkg_dd_costheta2_c0,bkg_dd_costheta2_c1,bkg_dd_costheta2_c2));
    bkg_ll_costheta_model = new RooPolynomial("bkg_ll_costheta_model",  "bkg_ll_costheta_model", costheta, RooArgList(bkg_ll_costheta_c0));//, bkg_ll_costheta_c1, bkg_ll_costheta_c2));
    bkg_ll_costheta1_model = new RooPolynomial("bkg_ll_costheta1_model","bkg_ll_costheta1_model",costheta1,RooArgList(bkg_ll_costheta1_c0));//,bkg_ll_costheta1_c1,bkg_ll_costheta1_c2));
    bkg_ll_costheta2_model = new RooPolynomial("bkg_ll_costheta2_model","bkg_ll_costheta2_model",costheta2,RooArgList(bkg_ll_costheta2_c0,bkg_ll_costheta2_c1,bkg_ll_costheta2_c2));
  }

  RooProdPdf bkg_dd_ang_model("bkg_dd_ang_model","bkg_dd_ang_model",RooArgList(*bkg_dd_costheta_model,*bkg_dd_costheta1_model,*bkg_dd_costheta2_model));
  RooProdPdf bkg_ll_ang_model("bkg_ll_ang_model","bkg_ll_ang_model",RooArgList(*bkg_ll_costheta_model,*bkg_ll_costheta1_model,*bkg_ll_costheta2_model));

  //  const RooArgSet* obs = bkg_dd_ang_model.getObservables();
  //  const RooArgList& pdfList = bkg_dd_ang_model.pdfList();
  //  for (int i=0;i<pdfList.getSize();i++) {
  //    std::cout << pdfList.at(i)->GetName() << std::endl;
  //  }
  //  mygetchar;

  //  RooAbsPdf&  bkg_dd_ang_model = *bkg_dd_costheta2_model;
  //  RooAbsPdf&  bkg_ll_ang_model = *bkg_ll_costheta2_model;

  corrpdfdd = &bkg_dd_ang_model;
  corrpdfll = &bkg_ll_ang_model;

  //  RooProdPdf bkg_ang_model("bkg_ang_model","bkg_ang_model",RooArgList(*bkg_costheta_model,*bkg_costheta1_model));
  //  RooAbsPdf& bkg_ang_model = *bkg_costheta2_model;
  //   r_0=1.;
  //   plotpdf(c,costheta2,sig_ang_model,"w3-costheta2-r0eq1.eps");
  //   bkg_costheta2_c1=-1.;
  //   plotpdf(c,costheta2,bkg_costheta2_model,"bkgmodel-costheta2.eps");
  //   r_0=0.;
  //   plotpdf(c,costheta2,sig_ang_model,"w3-costheta2-r0eq0.eps");
  //   return 1;



  //  RooAbsPdf* ang_model = sig_ang_model;
  //  if (doBd2JpsiKS0) ang_model=&bkg_dd_ang_model;


  RooRealVar Lambdab_ID(Lambdabname + "_ID",Lambdabname + "_ID",-10000.,10000.);
  RooRealVar piminus_TRACK_Type("piminus_TRACK_Type","piminus_TRACK_Type",0.,10.);
  RooRealVar Polarity("Polarity","Polarity",-10.,10.);


  if (oper=="merge") {


    RooArgList vars(P_b);
    if (useparam==1) vars.add(RooArgList(alpha_b,r_0,r_1));
    if (useparam==2) vars.add(RooArgList(alpha_b,beta,gammarrv));
    RooArgList fitvars(costheta,costheta1,costheta2);

    TString plotprefix("mctest-");
    if (userwmc) {plotprefix+="mcrw-";}else{plotprefix+="mcnonrw-";}
    if (doBd2JpsiKS0) {plotprefix+="Bd2JpsiKS0";}else{plotprefix+="Lb2JpsiL0";}
    const TString& thisdirectory = directorylist[0];
    plotprefix+="-"+thisdirectory; 
    switch (fitPolarity) {
    case 0: plotprefix+="-MagnetAll";break;
    case 1: plotprefix+="-MagnetUp";break;
    case 2: plotprefix+="-MagnetDown";break;
    }
    TString plotname=plotprefix + "-sim";

    
    RooFitResult* frnom = getRooFitResult(plotname + ".root" , "AngFit");
    const RooArgList& floatPars=frnom->floatParsFinal();
    for (int i=0;i<vars.getSize();i++) {
      RooRealVar* arg = (RooRealVar*)(vars.at(i));
      RooRealVar* var = (RooRealVar*)floatPars.find(arg->GetName());
      if (var) arg->setVal(var->getVal());
    }
    delete frnom;
    RooMCStudy mcstudy(*sig_ang_model,fitvars);
    for (int i=1;i<dorand;i++) {
      TString plotname=plotprefix + "-sim-toy"+istr(i,"07");    
      RooFitResult* fr_ang = getRooFitResult(plotname + ".root" , "AngFit");
      RooFitResult* fr_mass = getRooFitResult(plotname + ".root" , "MassFit");
      if (fr_mass and fr_ang) {
	std::cout << plotname << ": mass status = " << fr_mass->status() << " ; covQual = " << fr_mass->covQual() <<  std::endl;
	std::cout << plotname << ": ang  status = " << fr_ang->status() << " ; covQual = " << fr_ang->covQual() <<  std::endl;
	mcstudy.addFitResult(*fr_ang);
	delete fr_ang;
	delete fr_mass;
      }
    }
    //    plot_mcstudy(c,mcstudy, "datafit-mcstudy", vars);

    const RooDataSet& data = mcstudy.fitParDataSet();
    RooDataSet data2(data,"data2");

    TMatrixDSym* correl = data2.correlationMatrix(RooArgList(P_b,alpha_b,r_0,r_1));
    std::cout << "Correlation matrix ->" << std::endl;
    correl->Print("v");
    

    RooRealVar* ap2_rrv = (RooRealVar*)data2.addColumn(*ap2_rfv);
    RooRealVar* am2_rrv = (RooRealVar*)data2.addColumn(*am2_rfv);
    RooRealVar* bp2_rrv = (RooRealVar*)data2.addColumn(*bp2_rfv);
    RooRealVar* bm2_rrv = (RooRealVar*)data2.addColumn(*bm2_rfv);

    RooArgList varstoplot(vars);
    varstoplot.add(RooArgList(*ap2_rrv, *am2_rrv, *bp2_rrv, *bm2_rrv));

    data2.Print("v");
    data2.write("blu.txt");



    RooArgList* variationlist = addvariation(&data2,varstoplot);

    //redefine boundaries
    for (int i=0;i<varstoplot.getSize();i++) {
      RooRealVar* var = (RooRealVar*)varstoplot.at(i);
      var->setMin(minofdataset(&data2,*var));
      var->setMax(maxofdataset(&data2,*var));
      var->setBins(25);
    }

    //redefine boundaries
    for (int i=0;i<variationlist->getSize();i++) {
      RooRealVar* var = (RooRealVar*)variationlist->at(i);
      const double varmin=minofdataset(&data2,*var);
      const double varmax=maxofdataset(&data2,*var);
      if (varmax>-varmin) {
	var->setMin(-varmax);
	var->setMax(varmax);
      } else {
	var->setMin(varmin);
	var->setMax(-varmin);
      }
      var->setBins(25);
    }

    TString outputname("mctest-syst-");
    if (randacceptance) outputname+="randacceptance";
    if (randmccalib) outputname+="randmccalib";
    if (randsfit) outputname+="randsfit";
    if (randsigparams) outputname+="randsigparams";
    if (randalphalambda) outputname+="randalphalambda";
    if (randbeam) outputname+="randbeam";
    if (randphases) outputname+="randphases";


    plotdata(c, varstoplot, &data2, outputname, "", true);
    plotdata(c, *variationlist, &data2, outputname, "", true);
   
    return 0;
  }


  //read acceptance
  TString acceptancefile="fitacceptance3";
  if (userwmc) {
    acceptancefile+="-mcrw";
  } else {
    acceptancefile+="-mcnonrw";
  }
  if (usebdt) acceptancefile+="-wbdt";
  if (usetmvarectcut) acceptancefile+="-wtmvarectcut";
  if (usesimplecut) acceptancefile+="-wsimplecut";
  if (doBd2JpsiKS0) {
    acceptancefile+="-B0-";
  } else {
    acceptancefile+="-Lb-";
  }
  acceptancefile+=directorylist[0];

  RooRandom::randomGenerator()->SetSeed(0);

  if (usekeyspdf) {
    
    //    TString sigfilename("DVTuples-Lb2JpsiL0-MC11a-reco12a.reduced.root");
    TString sigfilename(rootdir + "DVTuples-Lb2JpsiL0-MC11a-reco12a.reduced.root");
    TString dirname("Tuple_JpsiL0_betas");
    if (doBd2JpsiKS0) {
      sigfilename=rootdir + "DVTuples-Bd2JpsiKS-MC11a.reduced.root";
      dirname="Tuple_JpsiKs_detached_betas";
      //    dirname="Tuple_JpsiKs_prescaled_betas";
    }
    //  const TString datafilename("DVTuples_data_stripping17.root");
    
    
    TFile sigfile(sigfilename);
    TDirectory* sigdir = (TDirectory*)sigfile.Get(dirname);
    TTree* sigtree = (TTree*)sigdir->Get("DecayTree");

    
    if (randbeam and dorand>0) {
      std::cout << "kikoo" << std::endl;
      sigtree->AddFriend(dirname + "-" + istr(dorand), TString(sigfilename).ReplaceAll(".root", ".angles.root"));
    } else {
      sigtree->AddFriend(dirname, TString(sigfilename).ReplaceAll(".root", ".angles.root"));
    }
    if (userwmc) sigtree->AddFriend(dirname, TString(sigfilename).ReplaceAll(".root", ".reweighted.root"));
    if (usebdt) sigtree->AddFriend(dirname, TString(sigfilename).ReplaceAll(".root", ".bdt.root"));

    const int maxevents=30000;

    TCut selrwmc("1==1");
    if (userwmc) {
      if (randmccalib and dorand>0) {
	selrwmc=TCut("rw" + istr(dorand) + "==1");
      } else {
	selrwmc=TCut("rw==1");
      }
    }

    RooDataSet* accdatall_full = filldataset("accdatall","accdatall",*angvars, sigtree, LL_selection + truemcsel(doBd2JpsiKS0) + selrwmc); 
    RooDataSet* accdatall = (RooDataSet*)accdatall_full->reduce(RooFit::EventRange(0,maxevents));
    std::cout << "Making accpdfll out of " << accdatall->numEntries() << " events" << std::endl;
    accpdfll = new RooNDKeysPdf("accpdfll","accpdfll",*angvars,*accdatall,"am");

    RooDataSet* accdatadd_full = filldataset("accdatadd","accdatadd",*angvars, sigtree, DD_selection + truemcsel(doBd2JpsiKS0) + selrwmc); 
    RooDataSet* accdatadd = (RooDataSet*)accdatadd_full->reduce(RooFit::EventRange(0,maxevents));
    std::cout << "Making accpdfdd out of " << accdatadd->numEntries() << " events" << std::endl;
    accpdfdd = new RooNDKeysPdf("accpdfdd","accpdfdd",*angvars,*accdatadd,"am");
    
  } else {

    if (randmccalib and dorand>0) {
      acceptancell = new calcacceptanceclass(acceptancefile + "-ll-toy" + istr(dorand,"06") + ".root");
      acceptancedd = new calcacceptanceclass(acceptancefile + "-dd-toy" + istr(dorand,"06") + ".root");
    } else {
      acceptancell = new calcacceptanceclass(acceptancefile + "-ll.root");
      acceptancedd = new calcacceptanceclass(acceptancefile + "-dd.root");
    }
    if (randacceptance and dorand>0) {
      acceptancell->setcoeff(true);
      acceptancedd->setcoeff(true);
    }
  }
  
  // if (doBd2JpsiKS0) {
  //   if (usekeyspdf) {
  //     RooArgSet *llvars*ddvars;A
  //     accpdfll = factorizedacc(costheta,costheta1,costheta2,"B0ll",llvars,7,5,4,"fitacceptance3-mcrw-B0-cc0c1c2-facpdf-ll.txt");
  //     accpdfdd = factorizedacc(costheta,costheta1,costheta2,"B0dd",ddvars,7,5,4,"fitacceptance3-mcrw-B0-cc0c1c2-facpdf-dd.txt");
  //   //    plotpdf(c,costheta,*accpdfll,"tmpjw-accpdfll-costheta.eps");
  //   //    plotpdf(c,costheta1,*accpdfll,"tmpjw-accpdfll-costheta1.eps");
  //   //    plotpdf(c,costheta2,*accpdfll,"tmpjw-accpdfll-costehta2.eps");
      
  //     // // // TRandom3 rnd(0);
  //     // // // double maxval=0.;
  //     // // // for (int i=0;i<1000000;i++) {
  //     // // // 	costheta=rnd.Rndm()*2.-1.;
  //     // // // 	costheta1=rnd.Rndm()*2.-1.;
  //     // // // 	costheta2=rnd.Rndm()*2.-1.;
  //     // // // 	//    std::cout << "getval = " << model->getVal(RooArgSet(costheta,costheta1,costheta2)) << std::endl;
  //     // // // 	double val=accpdfll->getVal(RooArgSet(costheta,costheta1,costheta2));
  //     // // // 	if (val<0. or val>1. or i<10) std::cout << costheta.getVal() << " " << costheta1.getVal() << " " << costheta2.getVal() << " " << val << std::endl;
  //     // // // 	if (val>maxval) maxval=val;
  //     // // // }
  //     // // // std::cout << "maxval = " << maxval << std::endl;
  //     //      mygetchar;
  //   }
    
  //   //    exit(1);
  //   //    llvars->readFromFile("fitacceptance3-mcrw-B0-cc0c1c2-facpdf-ll.txt",0,0,true);
  //   //    ddvars->readFromFile("fitacceptance3-mcrw-B0-cc0c1c2-facpdf-dd.txt",0,0,true);
  // } else {
  //   if (usekeyspdf) {
  //     RooArgSet *llvars,*ddvars;
  //     accpdfll = factorizedacc(costheta,costheta1,costheta2,"Lbll",llvars,5,4,7,"fitacceptance3-mcrw-Lb-cc0c1c2-facpdf-ll.txt");
  //     accpdfdd = factorizedacc(costheta,costheta1,costheta2,"Lbdd",ddvars,5,4,7,"fitacceptance3-mcrw-Lb-cc0c1c2-facpdf-dd.txt");
  //   }
  //   //    llvars->readFromFile("fitacceptance3-mcrw-Lb-cc0c1c2-facpdf-ll.txt",0,0,true);
  //   //    ddvars->readFromFile("fitacceptance3-mcrw-Lb-cc0c1c2-facpdf-dd.txt",0,0,true);
  // }




  const TString datafilename(rootdir + "DVTuples-Lb2JpsiL0-MC11a-reco12a.reduced.root");
  //  const TString datafilename("DVTuples_data_stripping19a.reduced.root");

  TTree* datatree[directorylist.size()];
  memset(datatree,0,directorylist.size()*sizeof(TTree*));


  TFile datafile(datafilename);
  for (unsigned int i=0;i<directorylist.size();i++) {
    const TString& dirname = directorylist[i];
    TDirectory* datadir = (TDirectory*)datafile.Get(dirname);
    if (!datadir) continue;
    datatree[i] = (TTree*)datadir->Get("DecayTree");
    if (usebdt) datatree[i]->AddFriend(dirname, TString(datafilename).ReplaceAll(".root",".bdt.root"));
    //    datatree[i]->AddFriend(dirname, TString(datafilename).ReplaceAll(".root",".ksmass.root"));
    if (randbeam and dorand>0) {
      //      datatree[i]->AddFriend(dirname + "-" + istr(dorand), TString(datafilename).ReplaceAll(".root",".angles.root"));
    } else {
      datatree[i]->AddFriend(dirname, TString(datafilename).ReplaceAll(".root",".angles.root"));
    }
    // std::vector<TString> varlist1,varlist2;
    // varlist1.push_back("costheta");
    // varlist1.push_back("costheta1");
    // varlist1.push_back("costheta2");
    // varlist2.push_back(mass.GetName());
    // plotprofile(c,datatree[i],varlist2, varlist1, TCut(TString(mass.GetName()) + "<5600.") , "datafit-testprofile");

    //    sigtree[i]->AddFriend(dirname, TString(datafilename).ReplaceAll(".root",".missvar.root"));
  }

  RooArgList vars(mass);
  if (randbeam and dorand>0) {
    vars.add(RooArgList(lambdab_PE,lambdab_PX,lambdab_PY,lambdab_PZ));
    vars.add(RooArgList(lambda0_PE,lambda0_PX,lambda0_PY,lambda0_PZ));
    vars.add(RooArgList(piminus_PE,piminus_PX,piminus_PY,piminus_PZ));
    vars.add(RooArgList(pplus_PE,pplus_PX,pplus_PY,pplus_PZ));
    vars.add(RooArgList(jpsi_PE,jpsi_PX,jpsi_PY,jpsi_PZ));
    vars.add(RooArgList(muplus_PE,muplus_PX,muplus_PY,muplus_PZ));
    vars.add(RooArgList(muminus_PE,muminus_PX,muminus_PY,muminus_PZ));
  } else {
    vars.add(*angvars);
  }
  vars.add(RooArgList(Lambdab_ID,piminus_TRACK_Type,Polarity));

  RooDataSet* data[directorylist.size()];
  //  RooDataSet* dataddup[directorylist.size()];
  //  RooDataSet* datallup[directorylist.size()];
  //  RooDataSet* datadddown[directorylist.size()];
  //  RooDataSet* datalldown[directorylist.size()];

  RooDataSet* datadd[directorylist.size()];
  RooDataSet* datall[directorylist.size()];

  RooDataSet* datadd_Lambdab0[directorylist.size()];
  RooDataSet* datall_Lambdab0[directorylist.size()];
  RooDataSet* datadd_Lambdab0bar[directorylist.size()];
  RooDataSet* datall_Lambdab0bar[directorylist.size()];

  if (doBd2JpsiKS0) {
    costheta.setBins(50);
    costheta1.setBins(50);
    costheta2.setBins(50);
  } else {
    costheta.setBins(25);
    costheta1.setBins(25);
    costheta2.setBins(25);
  }


  for (unsigned int i=0; i<directorylist.size(); i++) {
    if (!datatree[i]) continue;
    data[i] = filldataset("data-" + directorylist[i],"data-"+ directorylist[i], vars, datatree[i], selection );

    // mass=5621.;
    // costheta=0.5;
    // costheta1=1.;
    // costheta2=1.;
    // Lambdab_ID=5122.;
    // piminus_TRACK_Type=3.;
    // Polarity=1.;
    // data[i]->add(RooArgList(mass,costheta,costheta1,costheta2,Lambdab_ID,piminus_TRACK_Type,Polarity));


    // hrrv[1] = (RooRealVar*)data[i]->addColumn(h1);
    // hrrv[2] = (RooRealVar*)data[i]->addColumn(h2);
    // hrrv[3] = (RooRealVar*)data[i]->addColumn(h3);
    // hrrv[4] = (RooRealVar*)data[i]->addColumn(h4);
    // hrrv[5] = (RooRealVar*)data[i]->addColumn(h5);
    // hrrv[6] = (RooRealVar*)data[i]->addColumn(h6);
    // hrrv[7] = (RooRealVar*)data[i]->addColumn(h7);

    // for (int j=1;j<=7;j++) {
    //   hrrv[j]->setMin(-1.);
    //   hrrv[j]->setMax(1.);
    // }

    if (randbeam and dorand>0) {
      addangles(data[i], rnd);
    }

    if (usekeyspdf) {
      addacceptanceweightfactpdf(data[i], wacc, doBd2JpsiKS0, false, accpdfll, accpdfdd);
    } else {
      addacceptanceweight(data[i],    wacc,doBd2JpsiKS0,false,acceptancell,acceptancedd);
    }

    if (correctangdistrib) {
      addacceptanceweightfactpdf(data[i], wcorr, doBd2JpsiKS0, false, corrpdfll, corrpdfdd);
    }

    data[i]->Print("v");
    if (dorand==0) mygetchar;

  }
  datafile.Close();

  for (unsigned int i=0; i<directorylist.size(); i++) {
    if (!data[i]) continue;
    datadd[i] = (RooDataSet*)data[i]->reduce(Lambda0_DD); 
    datall[i] = (RooDataSet*)data[i]->reduce(Lambda0_LL); 

    datadd[i]->SetName("datadd"+istr(i));
    datall[i]->SetName("datall"+istr(i));

    datadd[i] = (RooDataSet*)datadd[i]->reduce(RooFit::EventRange(1,1800));
    datall[i] = (RooDataSet*)datall[i]->reduce(RooFit::EventRange(1,5500));

    //    dataddup[i] = (RooDataSet*)data[i]->reduce(MagnetUp + Lambda0_DD );
    //    datallup[i] = (RooDataSet*)data[i]->reduce(MagnetUp + Lambda0_LL );
    //    datadddown[i] = (RooDataSet*)data[i]->reduce(MagnetDown + Lambda0_DD );
    //    datalldown[i] = (RooDataSet*)data[i]->reduce(MagnetDown + Lambda0_LL );
    
    datadd_Lambdab0[i] = (RooDataSet*)data[i]->reduce(isLambdab0 + Lambda0_DD); 
    datall_Lambdab0[i] = (RooDataSet*)data[i]->reduce(isLambdab0 + Lambda0_LL); 
    datadd_Lambdab0bar[i] = (RooDataSet*)data[i]->reduce(isLambdab0bar + Lambda0_DD); 
    datall_Lambdab0bar[i] = (RooDataSet*)data[i]->reduce(isLambdab0bar + Lambda0_LL); 

    datadd_Lambdab0[i]->SetName("datadd_Lambdab0" + istr(i));
    datall_Lambdab0[i]->SetName("datall_Lambdab0" + istr(i));
    datadd_Lambdab0bar[i]->SetName("datadd_Lambdab0bar" + istr(i));
    datall_Lambdab0bar[i]->SetName("datall_Lambdab0bar" + istr(i));
    
  }


  const RooArgList sig_dd_vars(sig_dd_frac,sig_dd_mean,sig_dd_sigma,sig_dd_alpha1,sig_dd_alpha2,sig_dd_n1,sig_dd_n2);
  
  const RooArgList sig_ll_vars(sig_ll_frac,sig_ll_mean,sig_ll_sigma,sig_ll_alpha1,sig_ll_alpha2,sig_ll_n1,sig_ll_n2);

  const RooArgList mass_dd_pdfs(sig_dd_mass_model,*bkg_dd_mass_model,*bkg_dd_Lbbar_mass_model);
  const RooArgList mass_ll_pdfs(sig_ll_mass_model,*bkg_ll_mass_model,*bkg_ll_Lbbar_mass_model);
  
  if (doBd2JpsiKS0) {
    readvarsfromtfile(sig_ll_vars, TString("fitmcmass-Bd2JpsiKS0-") + directorylist[0] + "-mc-ll-fr.root", "", randsigparams and dorand!=0);
    readvarsfromtfile(sig_dd_vars, TString("fitmcmass-Bd2JpsiKS0-") + directorylist[0] + "-mc-dd-fr.root", "", randsigparams and dorand!=0);
  } else {
    readvarsfromtfile(sig_dd_vars, TString("fitmcmass-Lb2JpsiL0-") + directorylist[0] + "-mc-dd-fr.root", "", randsigparams and dorand!=0);
    readvarsfromtfile(sig_ll_vars, TString("fitmcmass-Lb2JpsiL0-") + directorylist[0] + "-mc-ll-fr.root", "", randsigparams and dorand!=0);
  }

  if (!doBd2JpsiKS0 and correctangdistrib) {
    RooArgSet corrvars(bkg_dd_costheta_c0,bkg_dd_costheta1_c0,bkg_dd_costheta2_c0);
    corrvars.add(RooArgSet(bkg_ll_costheta_c0,bkg_ll_costheta1_c0,bkg_ll_costheta2_c0));
    readvarsfromtfile(corrvars, "datafit-mcnonrw-Bd2JpsiKS0-Tuple_JpsiKs_detached_betas-sim.root", "AngFit", dorand!=0);
  }

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

  //  RooArgSet initialparams(bkg_costheta_c0,bkg_costheta1_c0,bkg_costheta2_c0,bkg_costheta2_c1,bkg_costheta2_c2,P_b,alpha_b,r_0,r_1);
  //  initialparams.add(RooArgSet(sig_dd_mean,sig_dd_sigma,sig_dd_frac));
  //  initialparams.add(RooArgSet(sig_ll_mean,sig_ll_sigma,sig_ll_frac));
  //  initialparams.writeToFile("initpars.txt");

  //  RooGaussian sig_dd_alpha1_constrain("sig_dd_alpha1_constrain", "sig_dd_alpha1_constrain", sig_dd_alpha1, RooFit::RooConst(sig_dd_alpha1.getVal()), RooFit::RooConst(sig_dd_alpha1.getError()));
  //  RooGaussian sig_dd_alpha2_constrain("sig_dd_alpha2_constrain", "sig_dd_alpha2_constrain", sig_dd_alpha2, RooFit::RooConst(sig_dd_alpha2.getVal()), RooFit::RooConst(sig_dd_alpha2.getError()));
  //  RooGaussian sig_dd_n1_constrain("sig_dd_n1_constrain", "sig_dd_n1_constrain", sig_dd_n1, RooFit::RooConst(sig_dd_n1.getVal()), RooFit::RooConst(sig_dd_n1.getError()));
  //  RooGaussian sig_dd_n2_constrain("sig_dd_n2_constrain", "sig_dd_n2_constrain", sig_dd_n2, RooFit::RooConst(sig_dd_n2.getVal()), RooFit::RooConst(sig_dd_n2.getError()));
  //  RooGaussian sig_dd_frac_constrain("sig_dd_frac_constrain", "sig_dd_frac_constrain", sig_dd_frac, RooFit::RooConst(sig_dd_frac.getVal()), RooFit::RooConst(sig_dd_frac.getError()));
  //  RooGaussian sig_ll_alpha1_constrain("sig_ll_alpha1_constrain", "sig_ll_alpha1_constrain", sig_ll_alpha1, RooFit::RooConst(sig_ll_alpha1.getVal()), RooFit::RooConst(sig_ll_alpha1.getError()));
  //  RooGaussian sig_ll_alpha2_constrain("sig_ll_alpha2_constrain", "sig_ll_alpha2_constrain", sig_ll_alpha2, RooFit::RooConst(sig_ll_alpha2.getVal()), RooFit::RooConst(sig_ll_alpha2.getError()));
  //  RooGaussian sig_ll_n1_constrain("sig_ll_n1_constrain", "sig_ll_n1_constrain", sig_ll_n1, RooFit::RooConst(sig_ll_n1.getVal()), RooFit::RooConst(sig_ll_n1.getError()));
  //  RooGaussian sig_ll_n2_constrain("sig_ll_n2_constrain", "sig_ll_n2_constrain", sig_ll_n2, RooFit::RooConst(sig_ll_n2.getVal()), RooFit::RooConst(sig_ll_n2.getError()));
  //  RooGaussian sig_ll_frac_constrain("sig_ll_frac_constrain", "sig_ll_frac_constrain", sig_ll_frac, RooFit::RooConst(sig_ll_frac.getVal()), RooFit::RooConst(sig_ll_frac.getError()));

  RooArgList testdd(sig_dd_frac,sig_dd_alpha1,sig_dd_alpha2,sig_dd_n1,sig_dd_n2);
  RooArgList testll(sig_ll_frac,sig_ll_alpha1,sig_ll_alpha2,sig_ll_n1,sig_ll_n2);
  RooFitResult* frdd=NULL;
  RooFitResult* frll=NULL;
  if (doBd2JpsiKS0) {
    frdd=getRooFitResult(TString("fitmcmass-Bd2JpsiKS0-") + directorylist[0] + "-mc-dd-fr.root");
    frll=getRooFitResult(TString("fitmcmass-Bd2JpsiKS0-") + directorylist[0] + "-mc-ll-fr.root");
  } else {
    frdd=getRooFitResult(TString("fitmcmass-Lb2JpsiL0-") + directorylist[0] + "-mc-dd-fr.root");
    frll=getRooFitResult(TString("fitmcmass-Lb2JpsiL0-") + directorylist[0] + "-mc-ll-fr.root");
  }
  if (!frdd or !frll) {
    std::cerr << "Can't read signal shape for RooMultiVarGaussian" << std::endl;
    return 1;
  }
  RooMultiVarGaussian sig_dd_constrain("sig_dd_constrain","sig_dd_constrain",testdd,  *frdd);
  RooMultiVarGaussian sig_ll_constrain("sig_ll_constrain","sig_ll_constrain",testll,  *frll);
  delete frdd;
  delete frll;
  //  RooArgSet massfitblu(    sig_dd_frac_constrain, sig_dd_alpha1_constrain, sig_dd_alpha2_constrain, sig_dd_n1_constrain, sig_dd_n2_constrain);
  //  massfitblu.add(RooArgSet(sig_ll_frac_constrain, sig_ll_alpha1_constrain, sig_ll_alpha2_constrain, sig_ll_n1_constrain, sig_ll_n2_constrain));
  //  RooCmdArg massfitconstrains =   RooFit::ExternalConstraints(massfitblu);

  RooCmdArg massfitconstrains =   RooCmdArg::none();
  if (floatSigParams) {
    massfitconstrains =   RooFit::ExternalConstraints(RooArgSet(sig_dd_constrain,sig_ll_constrain));
  }

  
  RooCmdArg angfitconstrains = RooCmdArg::none();
  if (floatalphalambda) {
    if (assumecp) {
      angfitconstrains = RooFit::ExternalConstraints(RooArgSet(alpha_lambda_constrain));
    } else {
      angfitconstrains = RooFit::ExternalConstraints(RooArgSet(alpha_lambda_constrain,alpha_lambdabar_constrain));
    }
  } else {
    alpha_lambda.setConstant();
    alpha_lambdabar.setConstant();
  }

//   if (!doBd2JpsiKS0) {
//     if (dontfitLb) {
//       nBkg_dd.setVal(0.);
//       nBkg_ll.setVal(0.);
//       nSig_dd.setVal(0.);
//       nSig_ll.setVal(0.);
//       nBkg_dd.setConstant();
//       nBkg_ll.setConstant();
//       nSig_dd.setConstant();
//       nSig_ll.setConstant();
//     }
    
//     if (dontfitLbbar) {
//       nBkg_dd_Lbbar.setVal(0.);
//       nBkg_ll_Lbbar.setVal(0.);
//       nSig_dd_Lbbar.setVal(0.);
//       nSig_ll_Lbbar.setVal(0.);
//       nBkg_dd_Lbbar.setConstant();
//       nBkg_ll_Lbbar.setConstant();
//       nSig_dd_Lbbar.setConstant();
//       nSig_ll_Lbbar.setConstant();
//     }
//   }
 
  if (dontfitdd) {	
    nBkg_dd.setVal(0.);
    nSig_dd.setVal(0.);
    nBkg_dd.setConstant();
    nSig_dd.setConstant();
    nBkg_dd_Lbbar.setVal(0.);
    nSig_dd_Lbbar.setVal(0.);
    nBkg_dd_Lbbar.setConstant();
    nSig_dd_Lbbar.setConstant();
    sig_dd_alpha1.setConstant();
    sig_dd_alpha2.setConstant();
    sig_dd_frac.setConstant();
    sig_dd_mean.setConstant();
    sig_dd_n1.setConstant();
    sig_dd_n2.setConstant();
    sig_dd_sigma.setConstant();
  }
      
  if (dontfitll) {
    nBkg_ll_Lbbar.setVal(0.);
    nSig_ll_Lbbar.setVal(0.);
    nBkg_ll_Lbbar.setConstant();
    nSig_ll_Lbbar.setConstant();
    nBkg_ll.setVal(0.);
    nSig_ll.setVal(0.);
    nBkg_ll.setConstant();
    nSig_ll.setConstant();
    sig_ll_alpha1.setConstant();
    sig_ll_alpha2.setConstant();
    sig_ll_frac.setConstant();
    sig_ll_mean.setConstant();
    sig_ll_n1.setConstant();
    sig_ll_n2.setConstant();
    sig_ll_sigma.setConstant();
  }


  


  for (unsigned int i=0; i<directorylist.size(); i++) {
    TString plotprefix("mctest-");
    if (userwmc) {
      plotprefix+="mcrw-";
    } else {
      plotprefix+="mcnonrw-";
    }
    if (doBd2JpsiKS0) {plotprefix+="Bd2JpsiKS0";}else {plotprefix+="Lb2JpsiL0";}
    const TString& thisdirectory = directorylist[i];
    plotprefix+="-"+thisdirectory;

    switch (fitPolarity) {
    case 0: plotprefix+="-MagnetAll";break;
    case 1: plotprefix+="-MagnetUp";break;
    case 2: plotprefix+="-MagnetDown";break;
    }
    
    if (doBd2JpsiKS0) {

      //      alpha_b=1.;
      //      alpha_b.setConstant();

      std::map<std::string,RooDataSet*> datamapping;
      if (!dontfitdd) datamapping[ "dd" ]    = datadd[i];
      if (!dontfitll) datamapping[ "ll" ]    = datall[i];


      RooDataSet simdata("simdata", "simdata", RooArgList(mass,costheta,costheta1,costheta2,wacc,Polarity,piminus_TRACK_Type,Lambdab_ID,cat), RooFit::Index(cat),RooFit::Import(datamapping));
      simdata.Print("v");

 
      std::map<std::string,RooAbsPdf*> masspdfmapping;
      if (!dontfitdd) masspdfmapping[ "dd" ] = tot_dd_mass_model;
      if (!dontfitll) masspdfmapping[ "ll" ] = tot_ll_mass_model;

      RooSimultaneous masssimpdf("masssimpdf","masssimpdf",masspdfmapping,cat);

      std::map<std::string,RooAbsPdf*> angpdfmapping;
      if (!dontfitdd) angpdfmapping[ "dd" ] = &bkg_dd_ang_model;
      if (!dontfitll) angpdfmapping[ "ll" ] = &bkg_ll_ang_model;
      RooSimultaneous angsimpdf("angsimpdf","angsimpdf",angpdfmapping,cat);

      std::map<std::string,RooAbsPdf*> angcos0pdfmapping;
      if (!dontfitdd) angcos0pdfmapping[ "dd" ] = bkg_dd_costheta_model;
      if (!dontfitll) angcos0pdfmapping[ "ll" ] = bkg_ll_costheta_model;
      RooSimultaneous angcos0simpdf("angcos0simpdf","angcos0simpdf",angcos0pdfmapping,cat);
      std::map<std::string,RooAbsPdf*> angcos1pdfmapping;
      if (!dontfitdd) angcos1pdfmapping[ "dd" ] = bkg_dd_costheta1_model;
      if (!dontfitll) angcos1pdfmapping[ "ll" ] = bkg_ll_costheta1_model;
      RooSimultaneous angcos1simpdf("angcos1simpdf","angcos1simpdf",angcos1pdfmapping,cat);
      std::map<std::string,RooAbsPdf*> angcos2pdfmapping;
      if (!dontfitdd) angcos2pdfmapping[ "dd" ] = bkg_dd_costheta2_model;
      if (!dontfitll) angcos2pdfmapping[ "ll" ] = bkg_ll_costheta2_model;
      RooSimultaneous angcos2simpdf("angcos2simpdf","angcos2simpdf",angcos2pdfmapping,cat);


      RooArgList masspdfs;
      masspdfs.add(mass_dd_pdfs);
      masspdfs.add(mass_ll_pdfs);

      RooArgList angpdflist(angcos0simpdf,angcos1simpdf,angcos2simpdf);

      TString plotname=plotprefix + "-sim";
      if (dorand!=0) plotname+="-toy"+istr(dorand,"07");
      fitwsplot(c,doBd2JpsiKS0, &simdata, masssimpdf, masspdfs, angpdflist,  plotname, massfitconstrains, RooCmdArg::none(), &datamapping, &masspdfmapping);

      // initialparams.readFromFile("initpars.txt");
      // fitwsplot(c, doBd2JpsiKS0, datall[i], tot_ll_mass_model, mass_ll_pdfs, *ang_model,  plotprefix + "-ll");

      // initialparams.readFromFile("initpars.txt");
      // fitwsplot(c, doBd2JpsiKS0, datadd[i], tot_dd_mass_model, mass_dd_pdfs, *ang_model,  plotprefix + "-dd");

      // initialparams.readFromFile("initpars.txt");
      // fitwsplot(c, doBd2JpsiKS0, datallup[i], tot_ll_mass_model, mass_ll_pdfs, *ang_model,  plotprefix + "-llup");

      // initialparams.readFromFile("initpars.txt");
      // fitwsplot(c, doBd2JpsiKS0, datalldown[i], tot_ll_mass_model, mass_ll_pdfs, *ang_model,  plotprefix + "-lldown");
      

      // initialparams.readFromFile("initpars.txt");
      // fitwsplot(c, doBd2JpsiKS0, dataddup[i], tot_dd_mass_model, mass_dd_pdfs, *ang_model,  plotprefix + "-ddup");
      // initialparams.readFromFile("initpars.txt");
      // fitwsplot(c, doBd2JpsiKS0, datadddown[i], tot_dd_mass_model, mass_dd_pdfs, *ang_model,  plotprefix + "-dddown");
      
    } else {

      //       initialparams.readFromFile("initpars.txt");
      //       fitwsplot(c, doBd2JpsiKS0, datadd[i], tot_dd_mass_model, mass_dd_pdfs, sig_ang_model,  plotprefix + "-dd-tot", RooFit::ExternalConstraints(RooArgSet(alpha_lambda_constrain)) );      
      //       initialparams.readFromFile("initpars.txt");
      //       fitwsplot(c, doBd2JpsiKS0, datadd_Lambdab0[i], tot_dd_mass_model, mass_dd_pdfs, sig_ang_model,  plotprefix + "-dd-Lb", RooFit::ExternalConstraints(RooArgSet(alpha_lambda_constrain)));      
      //       initialparams.readFromFile("initpars.txt");
      //       fitwsplot(c, doBd2JpsiKS0, datadd_Lambdab0bar[i], tot_dd_mass_model, mass_dd_pdfs, sig_ang_model_Lbbar,  plotprefix + "-dd-Lbbar", RooFit::ExternalConstraints(RooArgSet(alpha_lambdabar_constrain)));      
      
      
      //       initialparams.readFromFile("initpars.txt");
      //       fitwsplot(c, doBd2JpsiKS0, datall[i], tot_ll_mass_model, mass_ll_pdfs, sig_ang_model,  plotprefix + "-ll-tot", RooFit::ExternalConstraints(RooArgSet(alpha_lambda_constrain)));
      //       initialparams.readFromFile("initpars.txt");
      //       fitwsplot(c, doBd2JpsiKS0, datall_Lambdab0[i], tot_ll_mass_model, mass_ll_pdfs, sig_ang_model,  plotprefix + "-ll-Lb", RooFit::ExternalConstraints(RooArgSet(alpha_lambda_constrain)));
      //       initialparams.readFromFile("initpars.txt");
      //       fitwsplot(c, doBd2JpsiKS0, datall_Lambdab0bar[i], tot_ll_mass_model, mass_ll_pdfs, sig_ang_model_Lbbar,  plotprefix + "-ll-Lbbar", RooFit::ExternalConstraints(RooArgSet(alpha_lambdabar_constrain)));
      

      
      std::map<std::string,RooDataSet*> datamapping; 
      if (separateLbLbar) {
	if (!dontfitLb    and !dontfitdd) datamapping[ "ddLb" ]    = datadd_Lambdab0[i];
	if (!dontfitLbbar and !dontfitdd) datamapping[ "ddLbbar" ] = datadd_Lambdab0bar[i];
	if (!dontfitLb    and !dontfitll) datamapping[ "llLb" ]    = datall_Lambdab0[i];
	if (!dontfitLbbar and !dontfitll) datamapping[ "llLbbar" ] = datall_Lambdab0bar[i];
      } else {
	if (!dontfitdd) datamapping[ "dd" ]    = datadd[i];
	if (!dontfitll) datamapping[ "ll" ]    = datall[i];
      }

      RooDataSet simdata("simdata", "simdata", RooArgList(mass,costheta,costheta1,costheta2,wacc), RooFit::Index(cat),RooFit::Import(datamapping));
      simdata.Print("v");

 
      std::map<std::string,RooAbsPdf*> masspdfmapping;
      if (separateLbLbar) {
	if (!dontfitLb    and !dontfitdd) masspdfmapping[ "ddLb"    ] = tot_dd_mass_model;
	if (!dontfitLbbar and !dontfitdd) masspdfmapping[ "ddLbbar" ] = tot_dd_Lbbar_mass_model;
	if (!dontfitLb    and !dontfitll) masspdfmapping[ "llLb"    ] = tot_ll_mass_model;
	if (!dontfitLbbar and !dontfitll) masspdfmapping[ "llLbbar" ] = tot_ll_Lbbar_mass_model;
      } else {
	if (!dontfitdd) masspdfmapping[ "dd" ] = tot_dd_mass_model;
	if (!dontfitll) masspdfmapping[ "ll" ] = tot_ll_mass_model;
      }

      std::map<std::string,RooAbsPdf*> angpdfmapping;
      if (separateLbLbar) {
	if (!dontfitLb    and !dontfitdd) angpdfmapping[ "ddLb"    ] = sig_ang_model;
	if (!dontfitLbbar and !dontfitdd) angpdfmapping[ "ddLbbar" ] = sig_ang_model_Lbbar;
	if (!dontfitLb    and !dontfitll) angpdfmapping[ "llLb"    ] = sig_ang_model;
	if (!dontfitLbbar and !dontfitll) angpdfmapping[ "llLbbar" ] = sig_ang_model_Lbbar;
      } else {
	if (!dontfitdd) angpdfmapping[ "dd" ] = sig_ang_model;
	if (!dontfitll) angpdfmapping[ "ll" ] = sig_ang_model;
      }
      RooSimultaneous masssimpdf("masssimpdf","masssimpdf",masspdfmapping,cat);
      RooSimultaneous angsimpdf("angsimpdf","angsimpdf",angpdfmapping,cat);

      RooArgList masspdfs;
      masspdfs.add(mass_dd_pdfs);
      masspdfs.add(mass_ll_pdfs);

      TString plotname=plotprefix + "-sim";
      if (dorand!=0) plotname+="-toy"+istr(dorand,"07");

      RooArgList angpdflist(angsimpdf);

      fitwsplot(c,doBd2JpsiKS0, &simdata, masssimpdf, masspdfs, angpdflist,  plotname, massfitconstrains, angfitconstrains, &datamapping, &masspdfmapping);


      //      fitwsplot(c,doBd2JpsiKS0, &simdata, masssimpdf, masspdfs, angsimpdf,  plotname, massfitconstrains, angfitconstrains, &datamapping, &masspdfmapping);
      //      RooFitResult* fr_mass = masssimpdf.fitTo(simdata, RooFit::Save(true), RooFit::Minos(true));
      


      //      RooFitResult* fr = angsimpdf.fitTo(simdata, RooFit::SumW2Error(true), RooFit::Save(true), RooFit::Minos(false));      

      //      fr->Print("v");

      //      plot(c,RooArgList(costheta,costheta1,costheta2), simpdf, &simdata, plotprefix + "-simpdf", true);


    }
  }

  for (unsigned int i=0; i<directorylist.size(); i++) {
    delete data[i];
    delete datadd[i];
    delete datall[i];

    delete datadd_Lambdab0[i];
    delete datall_Lambdab0[i];
    delete datadd_Lambdab0bar[i];
    delete datall_Lambdab0bar[i];
  }

}
