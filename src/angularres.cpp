#include "TRandom3.h"
#include "TGraphErrors.h"

#include "RooDataSet.h"


#include "rootlogon.h"
#include "functions-roofit.h"
#include "RooLbtoJpsiL0PDF3.h"
#include "cuts.h"

RooRealVar costheta("costheta",   "cos#theta",     -1., 1.);
RooRealVar costheta1("costheta1", "cos#theta_{1}", -1., 1.);
RooRealVar costheta2("costheta2", "cos#theta_{2}", -1., 1.);

RooRealVar true_costheta("true_costheta",   "cos#theta",     -1., 1.);
RooRealVar true_costheta1("true_costheta1", "cos#theta_{1}", -1., 1.);
RooRealVar true_costheta2("true_costheta2", "cos#theta_{2}", -1., 1.);

const int maxbin=10;
const double pi=atan(1.)*4.;

RooDataSet* randomdataset(RooDataSet* data, TRandom3* rnd, TGraphErrors* resol[3])
{
  RooDataSet* newdata = new RooDataSet("dataww6","dataww6",RooArgList(costheta,costheta1,costheta2));
  const int nEntries=data->numEntries();

  for (int i=0;i<nEntries;i++) {     
    const RooArgSet* blu = data->get(i);
    const double theta=acos(blu->getRealValue("costheta"));
    const double theta1=acos(blu->getRealValue("costheta1"));
    const double theta2=acos(blu->getRealValue("costheta2"));

    double thetaresol,theta1resol,theta2resol, tmp;

    resol[0]->GetPoint(theta/pi*maxbin, tmp, thetaresol);
    resol[1]->GetPoint(theta1/pi*maxbin, tmp, theta1resol);
    resol[2]->GetPoint(theta2/pi*maxbin, tmp, theta2resol);
    //    std::cout << theta2/pi << " " << tmp << " " << theta2resol << std::endl;
    //    mygetchar;


    costheta=cos(theta*(1.+rnd->Gaus(0.,thetaresol)));
    costheta1=cos(theta1*(1.+rnd->Gaus(0.,theta1resol)));
    costheta2=cos(theta2*(1.+rnd->Gaus(0.,theta2resol)));

    //    target.setVal(scale*blu->getRealValue(origvar));    
    newdata->add(RooArgSet(costheta,costheta1,costheta2));
  }
  return newdata;
}

void usage(char* argv0)
{
  std::cerr << argv0 << " ll/dd " << std::endl;
  exit(1);
}

int main(int argc, char** argv)
{
  if (argc!=2) usage(argv[0]);

  TString llordd(argv[1]);
  llordd.ToLower();
  if (llordd!="dd" and llordd!="ll") usage(argv[0]);
  const bool doLL=llordd.Contains("ll");

  lhcbstyle();
  TCanvas c("c","c",100,100);

  const TString sigfilename("DVTuples-Lb2JpsiL0-MC11a-reco12a-nodecprodcut.root");  
  const TString dirname("Tuple_JpsiL0_betas");

  TFile sigfile(sigfilename);
  TDirectory* sigdir = (TDirectory*)sigfile.Get(dirname);
  TTree* sigtree = (TTree*)sigdir->Get("DecayTree");
  sigtree->AddFriend(dirname, TString(sigfilename).ReplaceAll(".root", ".angles.root"));
  
  RooFormulaVar thetares("thetares","(#theta^{meas}-#theta^{true})/#theta^{true}","(acos(@0)-acos(@1))/acos(@1)",RooArgList(costheta,true_costheta));
  RooFormulaVar theta1res("theta1res","(#theta_{1}^{meas}-#theta_{1}^{true})/#theta_{1}^{true}","(acos(@0)-acos(@1))/acos(@1)",RooArgList(costheta1,true_costheta1));
  RooFormulaVar theta2res("theta2res","(#theta_{2}^{meas}-#theta_{2}^{true})/#theta_{2}^{true}","(acos(@0)-acos(@1))/acos(@1)",RooArgList(costheta2,true_costheta2));


  RooDataSet* sigdata_orig;
  TString outputname("angularres");
  if (doLL) {
    outputname+="-ll";
    sigdata_orig = filldataset("sigdata_orig","sigdata_orig",RooArgList(costheta,costheta1,costheta2,true_costheta,true_costheta1,true_costheta2), sigtree, Lambda0_LL + truemcsel(false)); 
  } else {
    outputname+="-dd";
    sigdata_orig = filldataset("sigdata_orig","sigdata_orig",RooArgList(costheta,costheta1,costheta2,true_costheta,true_costheta1,true_costheta2), sigtree, Lambda0_DD + truemcsel(false)); 
  }
  
  RooRealVar* thetares_rrv = (RooRealVar*)sigdata_orig->addColumn(thetares);
  RooRealVar* theta1res_rrv = (RooRealVar*)sigdata_orig->addColumn(theta1res);
  RooRealVar* theta2res_rrv = (RooRealVar*)sigdata_orig->addColumn(theta2res);
  
  RooDataSet* sigdata = (RooDataSet*)sigdata_orig->reduce("abs(thetares)<0.03&&abs(theta2res)<0.03");
  //  RooDataSet* sigdata = sigdata_orig;

  //  thetares_rrv->setMin(-0.03);thetares_rrv->setMax(+0.03);
  //  theta1res_rrv->setMin(-0.2);theta1res_rrv->setMax(+0.2);
  //  theta2res_rrv->setMin(-0.03);theta2res_rrv->setMax(+0.03);

  TGraphErrors* resol[3];
  for (int i=0;i<3;i++) resol[i] = new TGraphErrors(10);

  for (int ivar=0;ivar<3;ivar++) {
    RooRealVar *var,*dvar;
    if (ivar==0) {
      var=&true_costheta;
      dvar=thetares_rrv;
    } else if (ivar==1) {
      var=&true_costheta1;
      dvar=theta1res_rrv;
    } else {
      var=&true_costheta2;
      dvar=theta2res_rrv;
    } 

    const TString varname(var->GetName());
    const TString vartitle(var->GetTitle());


    for (int ibin=0;ibin<maxbin;ibin++) {
      TString mycut("acos(" + varname + ")/3.14159>" + str(1.*ibin/maxbin) + "&&acos(" + varname + ")/3.14159<" + str(1.*(ibin+1)/maxbin));
      std::cout << "mycut = " << mycut << std::endl;
      RooDataSet* sigdata0 = (RooDataSet*)sigdata->reduce(mycut);

      RooRealVar mean("mean","mean",0.0,-1.,1.);
      RooRealVar sigma("sigma","sigma",0.01,0.,1.);
      plotdata(c,RooArgList(*dvar), sigdata0, outputname + "-" +istr(ibin), "", true, &mean, &sigma);
      resol[ivar]->SetPoint(ibin, 1./2./maxbin*(2.*ibin+1.), sigma.getVal());
      resol[ivar]->SetPointError(ibin, 1./maxbin/2., sigma.getError());

      delete sigdata0;
    }

    //    resol[ivar]->Fit("pol3");
    if (ivar==0) resol[ivar]->GetXaxis()->SetTitle("#theta^{true}/#pi");
    if (ivar==1) resol[ivar]->GetXaxis()->SetTitle("#theta_{1}^{true}/#pi");
    if (ivar==2) resol[ivar]->GetXaxis()->SetTitle("#theta_{2}^{true}/#pi");
    resol[ivar]->GetYaxis()->SetTitle("#sigma(" + TString(dvar->GetTitle()) + ")");
    resol[ivar]->SetMinimum(0.);
    resol[ivar]->Draw("ap");
    c.SaveAs(outputname + "-graph-" + varname + ".eps");
  }


  //  exit(1);    
  
  

  RooRealVar P_b("P_b","P_{b}",0.06, -1.5, 1.5);//0.50
  RooRealVar alpha_b("alpha_b","#alpha_{b}",-0.,-1.1,1.1);//0.457
  RooRealVar r_0("r_0","r_{0}",0.5,-0.1,1.2);//0.5
  RooRealVar r_1("r_1","r_{1}", -0.5,-1.2,1.2);//0.1
  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}", 0.642);

  RooFormulaVar ap2_rfv("ap2_rrv","|a+|^{2}","0.5*(@0+@1)",RooArgList(r_0,r_1));
  RooFormulaVar am2_rfv("am2_rrv","|a-|^{2}","0.5*(@0-@1)",RooArgList(r_0,r_1));
  RooFormulaVar bp2_rfv("bp2_rrv","|b+|^{2}","0.5*(1.+@0-@1-@2)",RooArgList(alpha_b,r_0,r_1));
  RooFormulaVar bm2_rfv("bm2_rrv","|b-|^{2}","0.5*(1.-@0-@1+@2)",RooArgList(alpha_b,r_0,r_1));


  RooLbtoJpsiL0PDF3* sig_ang_model       = new RooLbtoJpsiL0PDF3("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1);


  RooDataSet* dataorig = sig_ang_model->generate(RooArgList(costheta,costheta1,costheta2),5000);

  TRandom3 rnd(0);

  RooMCStudy mcstudy(*sig_ang_model,RooArgSet(costheta,costheta1,costheta2));

  
  RooFitResult* fr0 = sig_ang_model->fitTo(*dataorig, RooFit::Save(true), RooFit::PrintLevel(-1));
  if (isconverged(fr0)) {
    plot(c,RooArgList(costheta,costheta1,costheta2),*sig_ang_model,dataorig,outputname + "-firstfit");
    mcstudy.addFitResult(*fr0);
  } else {
    std::cout << "horrible mistake" << std::endl;
    exit(1);
  }
  delete fr0;

  const bool debug=false;

  progressbar pbar("running toys");
  const int maxtoy=500;
  for (int itoy=0; itoy<maxtoy; itoy++) {
    RooDataSet* newdata = randomdataset(dataorig,&rnd,resol);

    RooFitResult* fr = sig_ang_model->fitTo(*newdata,RooFit::Save(true),RooFit::PrintLevel(-1));

    if (0==itoy) {
    }

    if (debug) {
      fr->Print("v");
      mygetchar;
    }
    if (isconverged(fr)) mcstudy.addFitResult(*fr);
    
    delete fr;
    delete newdata;    
    pbar.print(100.*itoy/maxtoy);
  }
  pbar.finish();

  plot_mcstudy(c,mcstudy,outputname , RooArgList(P_b,alpha_b,r_0,r_1));

  RooDataSet data2(mcstudy.fitParDataSet(),"data2"); 

  RooRealVar* ap2_rrv = (RooRealVar*)data2.addColumn(ap2_rfv);
  RooRealVar* am2_rrv = (RooRealVar*)data2.addColumn(am2_rfv);
  RooRealVar* bp2_rrv = (RooRealVar*)data2.addColumn(bp2_rfv);
  RooRealVar* bm2_rrv = (RooRealVar*)data2.addColumn(bm2_rfv);

  RooArgList* variationlist = addvariation(&data2,RooArgList(P_b,alpha_b,r_0,r_1,*ap2_rrv,*am2_rrv,*bp2_rrv,*bm2_rrv));

  plotdata(c,*variationlist, &data2, outputname, "", true);

}
