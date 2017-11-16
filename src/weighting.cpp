#include "weighting.h"


TMatrixDSym makesymm(const TMatrixD& VCV)
{
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
  return VCVsym;
}


TMatrixDSym getVCV(const TMatrixDSym& V, const TMatrixDSym& C0)
{
  // Invert C
  TMatrixDSym C(C0);
  Double_t det(0) ;
  C.Invert(&det) ;
  if (det==0) {
    std::cerr << "getVCV: det==0 : correction matrix is singular" << std::endl ;
    TMatrixDSym VCVsym(V.GetNrows());
    return VCVsym;
  } else {
    
    // Calculate corrected covariance matrix = V C-1 V
    TMatrixD VCV(V,TMatrixD::kMult,TMatrixD(C,TMatrixD::kMult,V)) ; 
    return makesymm(VCV);
  }
}     

void computeamplitudes(RooRealVar& ap2, RooRealVar& am2, RooRealVar& bp2, RooRealVar& bm2, double Pb, double ab, double r0, double r1, const TMatrixDSym& Cov)
{
  TMatrixD M(4,4);
  for (int i=0;i<4;i++) for (int j=0;j<4;j++) M(i,j)=0.;
  
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
  
  TVectorD vec(4);
  vec(0)=Pb;//corr.getVal();
  vec(1)=ab;//corr.getVal();
  vec(2)=r0;//corr.getVal();
  vec(3)=r1;//corr.getVal();
  
  TVectorD vec2=M*vec;
  
  ap2=vec2(0);
  am2=vec2(1);
  bp2=0.5+vec2(2);
  bm2=0.5+vec2(3);
  
  TMatrixD Covampl(M,TMatrixD::kMult,TMatrixD(Cov,TMatrixD::kMult,MT)) ;
  ap2.setError(sqrt(Covampl(0,0)));
  am2.setError(sqrt(Covampl(1,1)));
  bp2.setError(sqrt(Covampl(2,2)));
  bm2.setError(sqrt(Covampl(3,3)));
  std::cout << "computeamplitudes::ap2 = " << ap2.getVal() << " +/- " << ap2.getError() << std::endl;
  std::cout << "computeamplitudes::am2 = " << am2.getVal() << " +/- " << am2.getError() << std::endl;
  std::cout << "computeamplitudes::bp2 = " << bp2.getVal() << " +/- " << bp2.getError() << std::endl;
  std::cout << "computeamplitudes::bm2 = " << bm2.getVal() << " +/- " << bm2.getError() << std::endl;
  
}


void deletecolumn(RooDataSet*& data, const RooRealVar& var)
{
  RooArgSet vars(*data->get());
  vars.remove(var,true,true);
  //  vars.Print("v");
  RooDataSet* reddata=(RooDataSet*)data->reduce(vars);
  //  reddata->Print("v");
  delete data;
  data=reddata;
}

//bool 
//{
//  
//}

struct point
{
  double costheta,costheta1,costheta2;
};
typedef std::vector<point> points;

void addpoint(points& list, double a, double b, double c)
{
  point tmp;
  tmp.costheta=a;
  tmp.costheta1=b;
  tmp.costheta2=c;
  list.push_back(tmp);
}

bool removesomepoints(double icostheta, double icostheta1, double icostheta2, double wmass, double wacc) 
{
  const double epsilon=0.01;
  
  points pointslist;

  // if (wacc>30.) {
  //   std::cerr << "Weird for c0=" << icostheta << " ; c1=" << icostheta1 << " ; c2=" << icostheta2 << " ; wmass=" << wmass << " ; wacc = " << wacc << std::endl;
  //   mygetchar;
  //   return true;
  // }

  // addpoint(pointslist,0.54623,0.987738,-0.0312308);
  // addpoint(pointslist,-0.0497229,0.914422,-0.990999);
  // addpoint(pointslist,0.54623,0.987738,-0.0312308);
  // addpoint(pointslist,0.505994,0.988784,0.167523);
  // addpoint(pointslist,-0.35751, 0.844934, -0.919114);
  // addpoint(pointslist,-0.81424, 0.931943, 0.987428);
  // addpoint(pointslist, -0.717746, 0.954329, -0.820874);
  // addpoint(pointslist, -0.974342, 0.487171, 0.99);
  // addpoint(pointslist,-0.75, 0.974342, -0.974342);
  // addpoint(pointslist,0.0877695, 0.729683, 0.902522);
  // addpoint(pointslist,-0.37406, 0.968159, -0.343404);

  for (unsigned int i=0;i<pointslist.size();i++) {
    if ((fabs(icostheta  - pointslist[i].costheta)<epsilon)  and 
	(fabs(icostheta1 - pointslist[i].costheta1)<epsilon) and 
	(fabs(icostheta2 - pointslist[i].costheta2)<epsilon)) {
      std::cerr << "Setting wacc=0. for c0=" << icostheta << " ; c1=" << icostheta1 << " ; c2=" << icostheta2 << " ; wmass=" << wmass << " ; wacc = " << wacc << std::endl;
      mygetchar;
      return true;
    }
  }

  return false;
}



RooAbsPdf* factorizedacc(RooRealVar& costheta, RooRealVar& costheta1, RooRealVar& costheta2, const TString& prefix, RooArgSet*& vars, int a, int b, int c, const TString& filename ) 
{
  RooArgList costhetacoeff;
  RooArgList costheta1coeff;
  RooArgList costheta2coeff;
  for (int i=0;i<a;i++) costhetacoeff.add(*(RooRealVar*)new RooRealVar(prefix + "-a"+istr(i),prefix + "-a"+istr(i), 0., -1., 1. ));
  for (int i=0;i<b;i++) costheta1coeff.add(*(RooRealVar*)new RooRealVar(prefix + "-b"+istr(i),prefix + "-b"+istr(i), 0., -1., 1. ));
  for (int i=0;i<c;i++) costheta2coeff.add(*(RooRealVar*)new RooRealVar(prefix + "-c"+istr(i),prefix + "-c"+istr(i), 0., -1., 1. ));

  vars = new RooArgSet();
  vars->add(costhetacoeff);
  vars->add(costheta1coeff);
  vars->add(costheta2coeff);
  
  if (filename!="") vars->readFromFile(filename,0,0,true);

  return (RooAbsPdf*)new RooProdPdf(prefix + "-facacc", prefix + "-facacc",RooArgList(*(RooChebychev*)new RooChebychev(prefix + "pol0",prefix + "pol0",costheta, costhetacoeff),
										      *(RooChebychev*)new RooChebychev(prefix + "pol1",prefix + "pol1",costheta1, costheta1coeff),
										      *(RooChebychev*)new RooChebychev(prefix + "pol2",prefix + "pol2",costheta2, costheta2coeff)));
}



void copyvar(RooDataSet* data, const TString& origvar, RooRealVar& target, double scale)
{
  RooDataSet dataww6("dataww6","dataww6",RooArgList(target));
  const int nEntries=data->numEntries();
  if (nEntries==0) {
    std::cerr << "copyvar: dataset is empty: " << data->GetName() << std::endl;
    mygetchar;
    return;
  }
  if (!data->get(0)->find(origvar)) {
    std::cerr << "copyvar: can't find " << origvar << " in " << data->GetName() << std::endl;
    mygetchar;
    return;
  }
  if (data->get(0)->find(target.GetName())) {
    std::cerr << "copyvar: " <<  target.GetName() << " already exists in " << data->GetName() << std::endl;
    mygetchar;
    return;
  }
  for (int i=0;i<nEntries;i++) {     
    const RooArgSet* blu = data->get(i);
    target.setVal(scale*blu->getRealValue(origvar));    
    dataww6.add(RooArgSet(target));
  }
  data->merge(&dataww6);
}

double calcPsample(RooDataSet* data, RooRealVar& mass, RooAbsPdf* sigpdf, RooAbsPdf* bkgpdf, double Fs)
{
  double Psample=0.;
  const int nEntries=data->numEntries();
  for (int i=0;i<nEntries;i++) {  
    const RooArgSet* blu = data->get(i);
    mass.setVal(blu->getRealValue(mass.GetName()));
    if (i<4) {
      std::cout << "sigpdf = " << sigpdf->getVal() << std::endl;
      std::cout << "bkgpdf = " << bkgpdf->getVal() << std::endl;
    }
    Psample+=pow( sigpdf->getVal()/bkgpdf->getVal() ,2);
  }
  std::cout << "Psample : " << data->GetName() << " " << Psample << " " << Fs << std::endl;
  return Psample*Fs;
}

double getsumwsumw2(RooDataSet* data, const TString& wname)
{
  const int nEntries=data->numEntries();

  double sumw=0.;
  double sumw2=0.;
  double w;
  for (int i=0;i<nEntries;i++) {  
    const RooArgSet* blu = data->get(i);
    if (wname=="") {
      w=data->weight();
    } else {
      w=blu->getRealValue(wname);
    }

    sumw+=w;
    sumw2+=w*w;
  }    
  return sumw/sumw2;
}


void addacceptanceweight(RooDataSet* data, RooRealVar& w,bool doBd2JpsiKS0, bool normalizeweights, calcacceptanceclass* mycalcacceptanceclass_ll, calcacceptanceclass* mycalcacceptanceclass_dd)
{
  const int nEntries=data->numEntries();
  if (nEntries==0) return;

  RooDataSet dataww6("dataww6","dataww6",RooArgList(w));

  //find max of weights
  double maxweight=-1000.;

  //  const double sumEntries=data->sumEntries();
  double weights[nEntries];

  //  bool doBd2JpsiKS0=true;
  
  double costheta,costheta1,costheta2,piminus_ttype,wmass;
  calcacceptanceclass* mycalcacceptanceclass;
  for (int i=0;i<nEntries;i++) {  
    const RooArgSet* blu = data->get(i);
    costheta=blu->getRealValue("costheta");
    costheta1=blu->getRealValue("costheta1");
    costheta2=blu->getRealValue("costheta2");
    wmass=blu->getRealValue("wmass");
    piminus_ttype = blu->getRealValue("piminus_TRACK_Type",0.);

    //    std::cout << "piminus_TRACK_Type == " << piminus_ttype << std::endl;
    
    if (piminus_ttype==3.) {
      mycalcacceptanceclass=mycalcacceptanceclass_ll;
    } else if (piminus_ttype==5.) {
      mycalcacceptanceclass=mycalcacceptanceclass_dd;
    } else {
      std::cout << "No piminus_TRACK_Type in dataset" << std::endl;
      exit(1);
    }

    const double value = mycalcacceptanceclass->evaluate(costheta,costheta1,costheta2);
    if (value<0.) {
      weights[i]=0.;
    } else if (value>=1.) {
      std::cerr << "costheta=" << costheta << " ; costheta1=" << costheta1 << " ; costheta2= " << costheta2 << " ; acc=" << value << std::endl;
      weights[i]=1./value;
    } else {
      weights[i]=1./value;
    }
    //    if (doBd2JpsiKS0) weights[i]/=1.-costheta2*costheta2;
    //    if (doBd2JpsiKS0) weights[i]*=1.-costheta2*costheta2;
    if (doBd2JpsiKS0 and fabs(costheta2)>0.98) weights[i]=0.;
    if (removesomepoints(costheta,costheta1,costheta2,wmass,weights[i])) weights[i]=0.;
    if (weights[i]>maxweight) maxweight=weights[i];

    //    std::cout << costheta << " " << costheta1 << " " << costheta2 << " " << weights[i] << std::endl;

  }
  
  if (!normalizeweights) {
    maxweight=1.;
  }
  
  
  for (int i=0;i<nEntries;i++) {
    w.setVal(weights[i] / maxweight);
    dataww6.add(RooArgSet(w));
  }
  data->merge(&dataww6);

}


void addacceptanceweightfactpdf(RooDataSet* data, RooRealVar& w, bool doBd2JpsiKS0, bool normalizeweights, RooAbsPdf* mycalcacceptanceclass_ll, RooAbsPdf* mycalcacceptanceclass_dd)
{
  const int nEntries=data->numEntries();
  if (nEntries==0) return;

  RooDataSet dataww6("dataww6","dataww6",RooArgList(w));

  //find max of weights
  double maxweight=-1000.;

  //  const double sumEntries=data->sumEntries();
  double weights[nEntries];

  RooArgSet* blublu = mycalcacceptanceclass_dd->getVariables();
  RooRealVar* costheta= (RooRealVar*)blublu->find("costheta");
  RooRealVar* costheta1= (RooRealVar*)blublu->find("costheta1");
  RooRealVar* costheta2= (RooRealVar*)blublu->find("costheta2");

  RooRealVar* phi1=(RooRealVar*)blublu->find("phi1");
  RooRealVar* phi2=(RooRealVar*)blublu->find("phi2");

  //  double costheta,costheta1,costheta2;
  RooAbsPdf* mycalcacceptanceclass;

  progressbar pbar(std::string("accacceptanceweight in " + TString(data->GetName()) ));
  for (int i=0;i<nEntries;i++) {  
    pbar.print(100.*i/nEntries);
    const RooArgSet* blu = data->get(i);

    const double costhetaval=blu->getRealValue("costheta");
    const double costheta1val=blu->getRealValue("costheta1");
    const double costheta2val=blu->getRealValue("costheta2");
    const double wmassval=blu->getRealValue("wmass");

    costheta->setVal(costhetaval);
    costheta1->setVal(costheta1val);
    costheta2->setVal(costheta2val);
    
    double phi1val=-666.;
    double phi2val=-666.;
    if (phi1 and phi2) {
      phi1val=blu->getRealValue("phi1");
      phi2val=blu->getRealValue("phi2");
      phi1->setVal(phi1val);
      phi2->setVal(phi2val);
    }


    if (blu->getRealValue("piminus_TRACK_Type")==3.) {
      mycalcacceptanceclass=mycalcacceptanceclass_ll;
    } else {
      mycalcacceptanceclass=mycalcacceptanceclass_dd;
    }

    //    const double value = mycalcacceptanceclass->evaluate(costheta,costheta1,costheta2);
    const double value = mycalcacceptanceclass->getVal();//RooArgSet(*costheta,*costheta1,*costheta2));
    //    std::cout << "costheta  = " << costhetaval << std::endl;
    //    std::cout << "costheta1 = " << costheta1val << std::endl;
    //    std::cout << "costheta2 = " << costheta2val << std::endl;
    //    std::cout << "phi1      = " << phi1val << std::endl;
    //    std::cout << "phi2      = " << phi2val << std::endl;
    //    std::cout << "value     = " << value << std::endl;
    // if (value>1.) std::cerr << "fuck" << std::endl;
    if (value<0.) { // or (doBd2JpsiKS0 and fabs(costheta2)>0.999)) {
      weights[i]=0.;
    } else {
      weights[i]=1./value;
      //      weights[i]=1.-value;
    }

    if (doBd2JpsiKS0) {
      weights[i]*=(1.-costheta2val*costheta2val);
    }
    if (removesomepoints(costhetaval,costheta1val,costheta2val,wmassval,weights[i])) weights[i]=0.;
    //    sumofweights+=weights[i];
    if (weights[i]>maxweight) maxweight=weights[i];
  }
  
  if (!normalizeweights) {
    maxweight=1.;
  }
  
  for (int i=0;i<nEntries;i++) {
    w.setVal(weights[i] / maxweight);
    dataww6.add(RooArgSet(w));
  }
  data->merge(&dataww6);
  pbar.finish();
}






void multiplyweights(RooDataSet* data, RooRealVar& wtot, const TString& waccname, const TString& wmassname, bool normalizetowmass, TString wcorrname)
{
  const int nEntries=data->numEntries();
  if (nEntries<=0) return;
  RooDataSet dataww6("dataww6","dataww6",RooArgList(wtot));
  
  const double defaultval=-666.;
  const double defaultwcorrval=-666.;

  double sumofwmass=0.;
  double sumofwacc=0.;
  double sumofwtot=0.;
  //  double sumofweights=0.;

  const bool debug=false;

  double wtotvec[nEntries];
  double wmassval,waccval;
  double wcorrval=1.;
  
  bool haswcorr=false;
  if (wcorrname!="") haswcorr=true;

  const RooArgSet* blu = data->get(0);
  wmassval=blu->getRealValue(wmassname, defaultval);;
  waccval=blu->getRealValue(waccname,defaultval);
  if (wmassval==defaultval) {
    std::cout << "No " << wmassname << " in dataset" << std::endl;
    exit(1);
  }
  if (waccval==defaultval) {
    std::cout << "No " << waccname << " in dataset" << std::endl;
    exit(1);
  }
  // if (haswcorr) {
  //   wcorrval=blu->getRealValue(wcorrname,defaultval);
  //   if (wcorrval==defaultval) {
  //     std::cout << "No " << wcorrname << " in dataset" << std::endl;
  //     exit(1);
  //   }
  // }
  for (int i=0;i<nEntries;i++) {  
    const RooArgSet* blu = data->get(i);
    wmassval=blu->getRealValue(wmassname, defaultval);;
    waccval=blu->getRealValue(waccname,defaultval);
    if (haswcorr) wcorrval=blu->getRealValue(wcorrname,defaultwcorrval);
    if (debug) std::cout << "wacc = " << waccval << " ; wmass = " << wmassval << " ; wcorr = " << wcorrval <<std::endl;

    wtotvec[i]= waccval*wmassval*wcorrval;

    sumofwmass+=wmassval;
    sumofwacc+=waccval;
    sumofwtot+=wtotvec[i];
  }

  //  if (debug) mygetchar;
  //rw weights
  // make that sum (wtot) == sum(wmass)
  for (int i=0;i<nEntries;i++) { 
    if (normalizetowmass) wtotvec[i]*=sumofwmass/sumofwtot;
    wtot.setVal(wtotvec[i]);
    dataww6.add(RooArgSet(wtot));
  }

  data->merge(&dataww6);
  
}

RooAbsReal* makenll(RooDataSet* tot_data, RooAbsPdf* sig_ang_model, const RooCmdArg& angfitconstrains)
{
  RooCmdConfig pc("blue");
  pc.defineSet("extCons","ExternalConstraints",0,0) ;
  RooLinkedList cmdList;
  cmdList.Add((TObject*)&angfitconstrains);
  pc.process(cmdList) ;
  

  RooNLLVar* nll = new RooNLLVar("nll","-log(L)",*sig_ang_model,*tot_data);
  
  //  double sumwsumw2=getsumwsumw2(tot_data,"wmass");
  //  double sumwsumw2=getsumwsumw2(tot_data,"wacc");
  const double sumwsumw2=getsumwsumw2(tot_data)*2.0;
  //  const double sumwsumw2=getsumwsumw2(tot_data);
 
  //  std::cout << "sumwsumw2 = " << sumwsumw2 << std::endl;

  RooRealVar* offset = new RooRealVar("offset","offset",sumwsumw2);
  
  RooProduct* prodnll = new RooProduct("prodnll","#alpha log(L)",RooArgSet(*nll,*offset));

  const RooArgSet* extCons = pc.getSet("extCons") ;

  RooAbsReal* finalnll;
  RooConstraintSum* nllCons=NULL;
  if (extCons) {
    //    allConstraints.add(*extCons) ;
    RooArgSet* cPars = sig_ang_model->getParameters(tot_data,kFALSE);
    nllCons = new RooConstraintSum("nllCons","nllCons",*extCons, *cPars);
    finalnll = new RooAddition("finalnll","finalnll",RooArgSet(*prodnll,*nllCons)) ;
    std::cout << "constrains added ...." << std::endl;
    extCons->Print("v");
  } else {
    finalnll = prodnll;
  }

  return finalnll;

}

RooFitResult* myfitto(RooDataSet* tot_data, RooAbsPdf* sig_ang_model, const RooCmdArg& angfitconstrains, int elevel, bool useminos)
{

  RooAbsReal* finalnll = makenll(tot_data,sig_ang_model,angfitconstrains);

  //  std::cout << elevel << std::endl;
  //  mygetchar;

  RooMinimizer m1(*finalnll);
  m1.setErrorLevel(elevel);
  m1.setPrintEvalErrors(elevel);
  m1.setEvalErrorWall(false);

  //  m1.minimize("Minuit","migrad") ;
  m1.hesse();
  m1.minimize("Minuit","migrad") ;

  //  m1.migrad();
  m1.hesse();
  if (useminos) m1.minos();
  RooFitResult* r1 = m1.save();
  //  if (extCons) {
  //    delete nllCons;
  //    delete finalnll;
  //  }
  delete finalnll;
  return r1;
}

void plotlikelihoodscans(const TString& plotprefix, const TCanvas& c, const RooArgList& vars, RooDataSet* tot_data, RooAbsPdf* sig_ang_model, const RooCmdArg& angfitconstrains)
{

  const int maxstep=100;
  const double nsigma=2.;

  RooFitResult* fr0 = myfitto(tot_data, sig_ang_model, angfitconstrains, -1, true);
  const double nll0=fr0->minNll();
  delete fr0;

  double varval0[vars.getSize()];
  double errval[vars.getSize()];

  for (int i=0;i<vars.getSize();i++) {
    RooRealVar* var = (RooRealVar*)vars.at(i);
    varval0[i]=var->getVal();
    errval[i]=var->getError();
  }

  for (int i=0;i<vars.getSize();i++) {
    RooRealVar* var = (RooRealVar*)vars.at(i);
    const TString varname(var->GetName());
    
    var->setConstant();
    TGraph graph(maxstep);
    for (int istep=0;istep<maxstep;istep++) {
      const double varval=varval0[i] - nsigma*errval[i] + 2.*nsigma*errval[i] * istep / maxstep;
      var->setVal( varval );
      
      RooFitResult* fr = myfitto(tot_data, sig_ang_model, angfitconstrains, -1, true);
      graph.SetPoint(istep, varval, fr->minNll()-nll0);
      delete fr;
    }

    graph.Draw("AC");
    c.SaveAs(plotprefix + "-likelihoodscan-" + varname + ".eps");

    var->setVal(varval0[i]);
    var->setConstant(false);
  }
}



void plotlikelihoodprofiles(const TString& plotprefix, const TCanvas& c, const RooArgList& vars, RooDataSet* tot_data, RooAbsPdf* sig_ang_model, const RooCmdArg& angfitconstrains)
{
  RooAbsReal* nll = makenll(tot_data,sig_ang_model,angfitconstrains);

  for (int i=0;i<vars.getSize();i++) {
    RooRealVar* var = (RooRealVar*)vars.at(i);
    const TString varname(var->GetName());
    //    const double nsigma=(varname.Contains("P_b"))?1:0.5;
    const double nsigma=0.5;
    std::cout << var->GetName() << ": " << var->getVal()-nsigma*var->getError() << " to " << var->getVal()+nsigma*var->getError() << std::endl;
    mygetchar;
    RooPlot* frame = var->frame(var->getVal()-nsigma*var->getError(), var->getVal()+nsigma*var->getError());
    nll->plotOn(frame,RooFit::ShiftToZero());
    frame->Draw();
    c.SaveAs(plotprefix + "-likelihoodprofile-" + varname + ".eps");
    delete frame;
  }

  delete nll;
}


void plotlikelihoodscans(const TString& plotprefix, const TCanvas& c, RooRealVar& varscan, double min, double max, const RooArgList& varlist, RooDataSet* tot_data, RooAbsPdf* sig_ang_model, const RooCmdArg& angfitconstrains)
{

  double othervar[varlist.getSize()];
  for (int i=0;i<varlist.getSize();i++) othervar[i]=((RooRealVar*)varlist.at(i))->getVal();


  RooFitResult* fr0 = myfitto(tot_data, sig_ang_model, angfitconstrains, -1, true);

  fr0->Print("v");
  mygetchar;

  const int nstep=50;

  const double var0=varscan.getVal();

  TGraph graph(nstep);

  varscan.setConstant();

  //from  minimum to maximum
  for (int i=0;i<=nstep;i++) {
    for (int j=0;j<varlist.getSize();j++) ((RooRealVar*)varlist.at(j))->setVal(othervar[j]);
    const double varval = min + (max-min) * i  / (nstep-1);
    varscan.setVal(varval);
    RooFitResult* fr = myfitto(tot_data, sig_ang_model, angfitconstrains, -1, false);
    graph.SetPoint(i, varval, fr->minNll() - fr0->minNll());
    delete fr;
  }

  // //from fit minimum to maximum
  // for (int i=0;i<=nstep;i++) {
  //   //    for (int j=0;j<varlist.getSize();j++) ((RooRealVar*)varlist.at(j))->setVal(othervar[j]);
  //   const double varval = min + (max-min) * i  / (nstep-1);
  //   if (varval>var0) {
  //     varscan.setVal(varval);
  //     RooFitResult* fr = myfitto(tot_data, sig_ang_model, angfitconstrains, -1, false);
  //     graph.SetPoint(i, varval, fr->minNll() - fr0->minNll());
  //     delete fr;
  //   }
  // }
  // //from fit minimum to minimum
  // for (int i=nstep;i>=0;i--) {
  //   //    for (int j=0;j<varlist.getSize();j++) ((RooRealVar*)varlist.at(j))->setVal(othervar[j]);
  //   const double varval = min + (max-min) * i  / (nstep-1);
  //   if (varval<var0) {
  //     varscan.setVal(varval);
  //     RooFitResult* fr = myfitto(tot_data, sig_ang_model, angfitconstrains, -1, false);
  //     graph.SetPoint(i, varval, fr->minNll() - fr0->minNll());
  //     delete fr;
  //   }
  // }

  varscan.setConstant(false);

  graph.GetXaxis()->SetTitle(varscan.GetTitle());
  graph.GetYaxis()->SetTitle("NLL");

  graph.Draw("A*C");
  c.SaveAs(plotprefix + "-graph-" + varscan.GetName() + ".eps");

  TFile outputfile(plotprefix + "-graph-" + varscan.GetName() + ".root", "recreate");
  graph.Write();
  outputfile.Close();


  varscan.setVal(var0);

  delete fr0;

}
