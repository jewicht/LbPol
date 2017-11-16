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
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooProdPdf.h"
#include "RooStats/SPlot.h"

#include "rootlogon.h"
#include "functions-roofit.h"
#include "cuts.h"


const bool usebdt=true;
const bool applybdtcut=true;
const bool applytmvarectcut=false;
const bool applysimplecut=false;

  //  const bool usesimplecuts=false;
  //  const bool antitrigger=false;
  

struct var12 {
  RooRealVar* var1;
  RooRealVar* var2;
};

RooRealVar mass("M","M(#LambdaJ/#psi)",5500.,5750.,"MeV/c^{2}");
//RooLbtoJpsiL0wAccPDF3* w3wacc = NULL;

RooRealVar nSig("nSig","nSig",2000.,0.,300000.);
//  RooRealVar nSig2("nSig2","nSig2",2000.,0.,5000.);
RooRealVar nBkg("nBkg","nBkg",10000.,0.,1000000.);

RooRealVar costheta("costheta",   "cos(#theta)",     -1., 1.);
RooRealVar costheta1("costheta1", "cos(#theta_{1})", -1., 1.);
RooRealVar costheta2("costheta2", "cos(#theta_{2})", -1., 1.);
RooRealVar phi1("phi1", "#phi_{1}/#pi", -1., 1.);
RooRealVar phi2("phi2", "#phi_{2}/#pi", -1., 1.);

// RooRealVar F1("F1","F1",-1.,1.);
// RooRealVar F2("F2","F2",-1.,1.);
// RooRealVar F3("F3","F3",-1.,1.);
// RooRealVar F4("F4","F4",-1.,1.);
// RooRealVar F5("F5","F5",-1.,1.);
// RooRealVar F6("F6","F6",-1.,1.);
// RooRealVar F7("F7","F7",-1.,1.);

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

  plot(c,mass,model,data, plotname + "-mass.eps");
  plot(c,mass,model,data, plotname + "-masswpull.eps", true);
}

void fitmassandfillsweights(TCanvas& c, RooAbsPdf& model, RooDataSet* data, const TString& plotname,bool dosplotplots, const RooArgList& splotvarstoplot)
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

  //  c.SetLogy();
  plot(c,mass,model,data, plotname, true);

  //  c.SetLogy(0);
  if (!dosplotplots) return;

  RooArgSet* blu = model.getVariables();
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
  while((arg=(RooRealVar*)iter->Next())) {
    TString varname(arg->GetName());
    if ((!varname.Contains("_M")) and !varname.Contains("cat") and !varname.Contains("nSig") and !varname.Contains("nBkg") and (!arg->isConstant())) {
      std::cout << "Making " << varname << " constant" << std::endl;
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
  yieldlist.add(yieldlist_ll);
  yieldlist.add(yieldlist_ll_Lbbar);
  yieldlist.add(yieldlist_dd);
  yieldlist.add(yieldlist_dd_Lbbar);

  RooStats::SPlot* sData = new RooStats::SPlot("sData-" + plotname,"An SPlot",
					       *data, &model, yieldlist );

  sData->Print("v");

  for (unsigned int i=0; i<varstounconst.size();i++) varstounconst[i]->setConstant(false);

  //   RooDataSet data_sig("data_sig","data_sig",data,*data->get(),0,"nSig_sw") ;
  //   RooDataSet data_bkg("data_bkg","data_bkg",data,*data->get(),0,"nBkg_sw") ;
  
  //   data_sig.get()->Print("v");
  
  //   plotdata(c,splotvarstoplot, &data_sig, plotname + "-splot-signal");
  //   mygetchar;
  //   plotdata(c,splotvarstoplot, &data_bkg, plotname + "-splot-bkg");


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
  plotdatavsdata(c,varstoplot, &data_sig, mc, plotname, true);
}

void compsigsplotmc2d(TCanvas& c, RooRealVar* var1, RooRealVar* var2, RooDataSet* data, RooDataSet* mc, const TString& plotname, int nbinx, int nbiny)
{
  
  RooDataSet data_sig("data_sig","data_sig",data,*data->get(),0,"nSig_sw") ;
  plotdatavsdata2d(c,var1,var2, &data_sig, mc, plotname + "-" +var1->GetName() + "_vs_" + var2->GetName(), true, true, nbinx, nbiny);
}



int usage(char* argv0)
{
  std::cerr << argv0 << " Lb/B0 ll/dd rwmc(0/1)" << std::endl;
  return 2;
}

int main(int argc, char** argv)
{

  //  rootlogon();
  lhcbstyle();

  if (argc!=4) return usage(argv[0]);

  const TString mode(argv[1]);
  const bool doBd2JpsiKS0=mode.Contains("B0");

  TString llordd(argv[2]);
  llordd.ToLower();
  if (llordd!="dd" and llordd!="ll") return usage(argv[0]);
  const bool doLL=llordd.Contains("ll");

  const int www=atoi(argv[3]);
  bool userwmc;
  if (www==1) {
    userwmc=true;
  } else if (www==0) {
    userwmc=false;
  } else {
    return usage(argv[0]);
  }

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

  
  //  TString dirname("JpsiL0Tuple_B2XMuMu");
  //  TString dirname("JpsiL0Tuple_betas");
  //  if (doBd2JpsiKS0) dirname = "JpsiKSTuple_prescaled_betas";
  //  if (doBd2JpsiKS0) dirname = "JpsiKSTuple_B2XMuMu";
  //  if (doBd2JpsiKS0) dirname = "JpsiKSTuple_detached_betas";
  

  //  TFile accfile("acc_histo.root");
  //  acc_histo = (TH3F*)accfile.Get("acc_histo");


  //Lambda_b0 to J/psi Lambda0

  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0,Lambdabname,Lambda0name,Jpsiname);
  
  mass.SetName(Lambdabname + "_M");
  mass.setUnit("MeV/c^{2}");

  if (doBd2JpsiKS0) {
    mass.setMin(5230.);
    mass.setMax(5450.);
    mass.SetTitle("M(J/#psiK_{S}^{0})");
  } else {
    mass.setMin(5500.);
    mass.setMax(5750.);
  }

  TString LbLatex="#Lambda_{b}^{0}";
  if (doBd2JpsiKS0) LbLatex="B^{0}";
  TString L0Latex="#Lambda";
  if (doBd2JpsiKS0) L0Latex="K_{S}^{0}";
  TString JpsiLatex="J/#psi";

  RooRealVar Lambdab_ENDVERTEX_CHI2(Lambdabname + "_ENDVERTEX_CHI2", LbLatex + ": vtx #chi^{2}", 0.);
  RooRealVar Lambdab_ENDVERTEX_NDOF(Lambdabname + "_ENDVERTEX_NDOF", LbLatex + ": vtx N_{dof.}", 0.);
  RooFormulaVar Lambdab_ENDVERTEX_CHI2NDOF(Lambdabname + "_ENDVERTEX_CHI2NDOF", 
                                           LbLatex + ": vtx #chi^{2}/N_{dof.}",
					   "@0/@1",
					   RooArgList(Lambdab_ENDVERTEX_CHI2,Lambdab_ENDVERTEX_NDOF));
  
  
  RooRealVar Jpsi_ENDVERTEX_CHI2(Jpsiname + "_ENDVERTEX_CHI2", JpsiLatex + ": vtx #chi^{2}", 0.);
  RooRealVar Jpsi_ENDVERTEX_NDOF(Jpsiname + "_ENDVERTEX_NDOF", JpsiLatex + ": vtx N_{dof.}", 0.);
  RooFormulaVar Jpsi_ENDVERTEX_CHI2NDOF(Jpsiname + "_ENDVERTEX_CHI2NDOF", 
                                        JpsiLatex + ": vtx #chi^{2}/N_{dof.}",
					"@0/@1",
                                        RooArgList(Jpsi_ENDVERTEX_CHI2,Jpsi_ENDVERTEX_NDOF));
  
  RooRealVar Lambdab_BDIRA(Lambdabname + "_DIRA_OWNPV", LbLatex + ": DIRA",0.);
  RooFormulaVar Lambdab_BDIRALOG(Lambdabname + "_BDIRALOG", 
				 LbLatex + ": log_{10}(1-DIRA)",
				 "log10(1.0000001-@0)",
				 RooArgList(Lambdab_BDIRA));
  RooRealVar Lambdab_IPCHI2(Lambdabname + "_IPCHI2_OWNPV",LbLatex + ": IP #chi^{2}",0.);
  RooRealVar Lambdab_FDCHI2(Lambdabname + "_FDCHI2_OWNPV",LbLatex + ": FD #chi^{2}",0.);
  RooFormulaVar Lambdab_FDCHI2LOG(Lambdabname + "_FDCHI2LOG", 
				  LbLatex + ": FD log_{10}(#chi^{2})",
				 "log10(@0)",
				  RooArgList(Lambdab_FDCHI2));
  RooRealVar Lambdab_PT(Lambdabname + "_PT",LbLatex + ": p_{T}",0.,"MeV");
  RooRealVar Lambdab_PZ(Lambdabname + "_PZ",LbLatex + ": p_{z}",0.,"MeV");
  RooRealVar Lambdab_PE(Lambdabname + "_PE",  LbLatex + ": E",0.,"MeV");
  RooRealVar Lambdab_P(Lambdabname + "_P",  LbLatex + ": p",0.,"MeV");

  RooRealVar Lambdab_eta(Lambdabname + "_eta",LbLatex + ": #eta",0.);

  RooFormulaVar Lambdab_PT_rfv(Lambdabname + "_PT_rfv",LbLatex + ": p_{T}","@0/1000.",RooArgList(Lambdab_PT));
  RooFormulaVar Lambdab_P_rfv(Lambdabname + "_P_rfv",  LbLatex + ": p", "@0/1000."   ,RooArgList(Lambdab_P));
  Lambdab_PT_rfv.setUnit("GeV");
  Lambdab_P_rfv.setUnit("GeV");


  //  RooFormulaVar Lambdab_y(Lambdabname + "_y", LbLatex + ": y", "0.5*log((@0+@1)/(@0-@1))", RooArgList(Lambdab_PE,Lambdab_PZ));
  RooRealVar Lambdab_y(Lambdabname + "_y", LbLatex + ": y",0.);

  RooRealVar Lambdab_TAU(Lambdabname + "_TAU", LbLatex + ": #tau",0.);
  //  RooRealVar Lambdab_TAU(Lambdabname + "_TAU",Lambdabname + "_TAU",0.,0.01);
  //  RooRealVar Lambda0_TAU(Lambda0name + "_TAU",Lambda0name + "_TAU",0.,1.);


  RooRealVar Lambda0_ENDVERTEX_CHI2(Lambda0name + "_ENDVERTEX_CHI2", L0Latex + ": vtx #chi^{2}", 0.);
  RooRealVar Lambda0_ENDVERTEX_NDOF(Lambda0name + "_ENDVERTEX_NDOF", L0Latex + ": vtx N_{dof.}", 0.);
  RooFormulaVar Lambda0_ENDVERTEX_CHI2NDOF(Lambda0name + "_ENDVERTEX_CHI2NDOF", 
                                           L0Latex + ": vtx #chi^{2}/N_{dof.}",
					   "@0/@1",
                                           RooArgList(Lambda0_ENDVERTEX_CHI2,Lambda0_ENDVERTEX_NDOF));
  RooRealVar Lambda0_P(Lambda0name + "_P",L0Latex + ": p",0.,"MeV");
  RooRealVar Lambda0_PE(Lambda0name + "_PE",L0Latex + ": E",0.,"MeV");
  RooRealVar Lambda0_PT(Lambda0name + "_PT",L0Latex + ": p_{T}",0.,"MeV");
  RooRealVar Lambda0_PZ(Lambda0name + "_PZ",L0Latex + ": p_{z}",0.,"MeV");
  RooRealVar Lambda0_M(Lambda0name + "_M",L0Latex + ": M",0.,"MeV/c^{2}");
  RooRealVar Lambda0_TAU(Lambda0name + "_TAU",L0Latex + ": #tau",0.);
  RooRealVar Lambda0_phi(Lambda0name + "_phi",L0Latex + ": #phi",0.);
  RooRealVar Lambda0_ENDVERTEX_Z(Lambda0name + "_ENDVERTEX_Z", L0Latex + ": vtx z", 0.);  

  //  RooFormulaVar Lambda0_y(Lambda0name + "_y", L0Latex + ": y", "0.5*log((@0+@1)/(@0-@1))", RooArgList(Lambda0_PE,Lambda0_PZ));
  RooRealVar Lambda0_y(Lambda0name + "_y", L0Latex + ": y", 0.);

  RooFormulaVar Lambda0_DM(Lambda0name + "_DM", L0Latex + ": #DeltaM","abs(@0-" + str(getLambda0mass(doBd2JpsiKS0)) + ")",RooArgList(Lambda0_M));
  Lambda0_DM.setUnit("MeV/c^{2}");

  RooRealVar Lambda0_eta(Lambda0name + "_eta",L0Latex + ": #eta",0.);


  RooRealVar Lambdab_ID(Lambdabname + "_ID",Lambdabname + "_ID",-10000.,10000.);
  RooRealVar piminus_TRACK_Type("piminus_TRACK_Type","piminus_TRACK_Type",0.,10.);
  RooRealVar Polarity("Polarity","Polarity",-10.,10.);
  
  RooRealVar piminus_P("piminus_P","#pi: p_{tot}",0.,"MeV");
  RooRealVar pplus_P("pplus_P","p: p_{tot}",0.,"MeV");

  RooRealVar piminus_PT("piminus_PT","#pi: p_{T}",0.,"MeV");
  RooRealVar piminus_PZ("piminus_PZ","#pi: p_{z}",0.,"MeV");
  RooRealVar piminus_PE("piminus_PE","#pi: E",0.,"MeV");
  //  RooFormulaVar piminus_y("piminus_y", "#pi: y", "0.5*log((@0+@1)/(@0-@1))", RooArgList(piminus_PE,piminus_PZ));
  RooRealVar piminus_y("piminus_y", "#pi: y", 0.);


  RooRealVar pplus_PT("pplus_PT","p: p_{T}",0.,"MeV");
  RooRealVar pplus_y("pplus_y", "p: y", 0.);  
  if (doBd2JpsiKS0) {
    pplus_P.SetName("piplus_P");
    pplus_P.SetTitle("#pi^{+}: p_{tot}");
    piminus_P.SetTitle("#pi^{-}: p_{tot}");
    pplus_PT.SetName("piplus_PT");
    pplus_PT.SetTitle("#pi^{+}: p_{T}");
    piminus_PT.SetTitle("#pi^{-}: p_{T}");

    //    pplus_PZ.SetName("piplus_PZ");
    //    pplus_PZ.SetTitle("#pi^{+}: p_{z}");
    piminus_PZ.SetTitle("#pi^{-}: p_{z}");

    pplus_y.SetName("piplus_y");
    pplus_y.SetTitle("#pi^{+}: y");

  }

  RooRealVar bdt("bdt","BDT",0.);

  RooRealVar Jpsi_M(Jpsiname  + "_MM",JpsiLatex  + ": MM",0. ,"MeV/c^{2}");
  RooFormulaVar Jpsi_DMM(Jpsiname + "_DMM", JpsiLatex + ": #DeltaMM", "abs(@0-" + str(Jpsipdgmass) + ")",RooArgList(Jpsi_M));
  Jpsi_DMM.setUnit("MeV/c^{2}");

  //pid
  RooRealVar piminus_ProbNNk ("piminus_ProbNNk", "piminus_ProbNNk",0.);
  RooRealVar piminus_ProbNNp ("piminus_ProbNNp", "piminus_ProbNNp",0.);
  RooRealVar piminus_ProbNNpi("piminus_ProbNNpi","piminus_ProbNNpi",0.);
  RooRealVar piminus_PIDK    ("piminus_PIDK", "piminus_PIDK",0.);
  RooRealVar piminus_PIDp    ("piminus_PIDp", "piminus_PIDp",0.);

  RooFormulaVar piminus_kpi  ("piminus_kpi","#pi: P_{K}-P_{#pi}","@0-@1",RooArgList(piminus_ProbNNk,piminus_ProbNNpi));;
  RooFormulaVar piminus_ppi  ("piminus_ppi","#pi: P_{p}-P_{#pi}","@0-@1",RooArgList(piminus_ProbNNp,piminus_ProbNNpi));;

  TString pplusname="pplus";
  if (doBd2JpsiKS0) pplusname="piplus";
  
  RooRealVar pplus_ProbNNk (pplusname + "_ProbNNk" ,pplusname + "_ProbNNk",0.);
  RooRealVar pplus_ProbNNp (pplusname + "_ProbNNp" ,pplusname + "_ProbNNp",0.);
  RooRealVar pplus_ProbNNpi(pplusname + "_ProbNNpi",pplusname + "_ProbNNpi",0.);
  RooRealVar pplus_PIDK    (pplusname + "_PIDK", pplusname + "_PIDK",0.);
  RooRealVar pplus_PIDp    (pplusname + "_PIDp", pplusname + "_PIDp",0.);
  RooFormulaVar pplus_pk   (pplusname + "_kp","p: P_{p}-P_{K}","@0-@1",RooArgList(pplus_PIDp,pplus_PIDK));;
  //  RooFormulaVar pplus_pip  (pplusname + "_pip","p: P_{#pi}-P_{p}","@0-@1",RooArgList(pplus_ProbNNpi,pplus_ProbNNp));;




  RooRealVar P_b("P_b","P_{b}",0.0, -1.5, 1.5);//0.50
  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#Lambda}",0.642,0.,1.);
  alpha_lambda.setConstant();

  if (doBd2JpsiKS0) alpha_lambda=0.;

  RooRealVar alpha_b("alpha_b","#alpha_{b}",0.,-1.1,1.1);//0.457
  RooRealVar r_0("r_0","r_{0}",0.5,-0.3,2.3);//0.5
  RooRealVar r_1("r_1","r_{1}",0.,-1.3,1.3);//0.1

  costheta.setBins(20);
  costheta1.setBins(20);
  costheta2.setBins(20);
  phi1.setBins(20);
  phi2.setBins(20);

  mass.setBins(50);
  mass.setRange("signalregion", 5600.,5640.);
  mass.setRange("lowersideband",5500.,5570.);
  mass.setRange("uppersideband",5670.,5750.);

  //  c1.setConstant();
  //  RooChebychev pol("pol","pol",mass,RooArgList(c0,c1));
  RooRealVar c0("c0","c0",0.,-1.1,1.1);
  //  c0.setConstant();
  RooRealVar c1("c1","c1",0.,-1.1,1.1);
  RooChebychev bkg_mass_model("bkg_mass_model","bkg_mass_model",mass,RooArgList(c0));

  //  RooLbtoJpsiL0PDF3 sig_ang_model("sig_ang_model","sig_ang_model",costheta,costheta1,costheta2,P_b, alpha_b, alpha_lambda, r_0, r_1);
  RooRealVar sig_dd_mean("sig_dd_mean","sig_dd_mean",5620.,5500.,5700.);
  RooRealVar sig_ll_mean("sig_ll_mean","sig_ll_mean",5620.,5500.,5700.);
  if (doBd2JpsiKS0) {
    sig_dd_mean.setVal(5280.);
    sig_dd_mean.setMin(5270.);
    sig_dd_mean.setMax(5290.);
    sig_ll_mean.setVal(5280.);
    sig_ll_mean.setMin(5270.);
    sig_ll_mean.setMax(5290.);
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

  
  RooAbsPdf& sig_dd_mass_model = sig_dd_twocb;
  RooAbsPdf& sig_ll_mass_model = sig_ll_twocb;

  RooRealVar bkg_dd_c0("bkg_dd_c0","bkg_dd_c0",0.,-1.1,1.1);
  RooRealVar bkg_ll_c0("bkg_ll_c0","bkg_ll_c0",0.,-1.1,1.1);

  RooChebychev bkg_dd_mass_model("bkg_dd_mass_model","bkg_dd_mass_model",mass,RooArgList(bkg_dd_c0));
  RooChebychev bkg_ll_mass_model("bkg_ll_mass_model","bkg_ll_mass_model",mass,RooArgList(bkg_ll_c0));

  RooAddPdf tot_dd_mass_model("tot_dd_mass_model","tot_dd_mass_model",RooArgList(sig_dd_mass_model,bkg_dd_mass_model),RooArgList(nSig,nBkg));
  RooAddPdf tot_ll_mass_model("tot_ll_mass_model","tot_ll_mass_model",RooArgList(sig_ll_mass_model,bkg_ll_mass_model),RooArgList(nSig,nBkg));

  const RooArgList sig_dd_vars(sig_dd_frac,sig_dd_mean,sig_dd_sigma,sig_dd_alpha1,sig_dd_alpha2,sig_dd_n1,sig_dd_n2);
  const RooArgList sig_dd_varsnottoconst(sig_dd_mean,sig_dd_sigma);
  //  const RooArgList mass_dd_pdfs(sig_dd_mass_model,bkg_dd_mass_model,bkg_dd_Lbbar_mass_model);
  
  const RooArgList sig_ll_vars(sig_ll_frac,sig_ll_mean,sig_ll_sigma,sig_ll_alpha1,sig_ll_alpha2,sig_ll_n1,sig_ll_n2);
  const RooArgList sig_ll_varsnottoconst(sig_ll_mean,sig_ll_sigma);
  //  const RooArgList mass_ll_pdfs(sig_ll_mass_model,bkg_ll_mass_model,bkg_ll_Lbbar_mass_model);
  
  if (doBd2JpsiKS0) {
    readvarsfromtfile(sig_dd_vars,TString("fitmcmass-Bd2JpsiKS0-") + directorylist[0] + "-mc-dd-fr.root");
    readvarsfromtfile(sig_ll_vars,TString("fitmcmass-Bd2JpsiKS0-") + directorylist[0] + "-mc-ll-fr.root");
  } else {
    readvarsfromtfile(sig_dd_vars,TString("fitmcmass-Lb2JpsiL0-") + directorylist[0] + "-mc-dd-fr.root");
    readvarsfromtfile(sig_ll_vars,TString("fitmcmass-Lb2JpsiL0-") + directorylist[0] + "-mc-ll-fr.root");
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


  TCanvas c("c","c",100,100);


  TString mcfilename("DVTuples-Lb2JpsiL0-MC11a-reco12a.reduced.root");
  if (doBd2JpsiKS0) {
    mcfilename = "DVTuples-Bd2JpsiKS-MC11a.reduced.root";
  }
  //  if (doBd2JpsiKS0) mcfilename = "small.root";

  const TString datafilename("DVTuples_data_stripping17.reduced.root");


  //  std::cout << "directory name = " << dirname << std::endl;
  //  plotprefix+="-"+dirname;

  std::cout << "data file      = " << datafilename << std::endl;
  std::cout << "mc file        = " << mcfilename << std::endl;
  std::cout << "directory list = ";
  for (unsigned int i=0;i<directorylist.size();i++) std::cout << directorylist[i] << " ";
  std::cout << std::endl;

  TTree* datatree[directorylist.size()];
  TTree* mctree[directorylist.size()];

  memset(datatree,0,directorylist.size()*sizeof(TTree*));
  memset(mctree,0,directorylist.size()*sizeof(TTree*));

  TFile datafile(datafilename);
  for (unsigned int i=0;i<directorylist.size();i++) {
    const TString& dirname = directorylist[i];
    TDirectory* datadir = (TDirectory*)datafile.Get(dirname);
    if (!datadir) continue;
    datatree[i] = (TTree*)datadir->Get("DecayTree");
    if (usebdt) datatree[i]->AddFriend(dirname, TString(datafilename).ReplaceAll(".root",".bdt.root"));
    datatree[i]->AddFriend(dirname, TString(datafilename).ReplaceAll(".root",".angles.root"));
    //    sigtree[i]->AddFriend(dirname, TString(datafilename).ReplaceAll(".root",".missvar.root"));
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
    if (usebdt) mctree[i]->AddFriend(dirname, TString(mcfilename).ReplaceAll(".root",".bdt.root"));
    mctree[i]->AddFriend(dirname, TString(mcfilename).ReplaceAll(".root",".angles.root"));
    if (userwmc) mctree[i]->AddFriend(dirname, TString(mcfilename).ReplaceAll(".root",".reweighted.root"));
  }


  TCut LL_selection;
  TCut DD_selection;
  //  simplecut(doBd2JpsiKS0, LL_selection,DD_selection)
  //bdt selection
  if (applybdtcut) {
    bdtsel(doBd2JpsiKS0, LL_selection,DD_selection);
  }
  if (applytmvarectcut) tmvarectcut(doBd2JpsiKS0, LL_selection,DD_selection);
  if (applysimplecut) simplecut2(doBd2JpsiKS0, LL_selection,DD_selection);

  const TCut TriggerCut=triggercut(doBd2JpsiKS0);
  const TCut trueMC = truemcsel(doBd2JpsiKS0);
  const TCut MassCut = masscut(doBd2JpsiKS0, mass); 

  LL_selection+=MassCut;
  DD_selection+=MassCut;

  LL_selection+=TriggerCut;
  DD_selection+=TriggerCut;
  
  TCut additionalmcsel = trueMC;// + TCut("L0DiMuonDecision==1");
  if (userwmc) additionalmcsel+=TCut("rw==1");
  const TCut selection = LL_selection or DD_selection;

  LL_selection+=Lambda0_LL;
  DD_selection+=Lambda0_DD;
  
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
  //  std::cout << "selection = "   << selection << std::endl;
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


  

  RooArgList vars(mass,costheta,costheta1,costheta2,phi1,phi2);
  //  vars.add(RooArgList(F1,F2,F3,F4,F5,F6,F7));
  vars.add(RooArgList(Lambdab_ID,piminus_TRACK_Type,Polarity));
  //Lambdab stuff
  vars.add(RooArgList(Lambdab_BDIRA,Lambdab_IPCHI2,Lambdab_FDCHI2,Lambdab_PT,Lambdab_P,Lambdab_TAU,Lambdab_eta,Lambdab_PZ,Lambdab_PE));
  vars.add(RooArgList(Lambdab_ENDVERTEX_CHI2,Lambdab_ENDVERTEX_NDOF));
  //Lambda0 stuff
  vars.add(RooArgList(Lambda0_P,Lambda0_M,Lambda0_PT,Lambda0_TAU,Lambda0_ENDVERTEX_CHI2,Lambda0_ENDVERTEX_NDOF,Lambda0_ENDVERTEX_Z,Lambda0_phi,Lambda0_eta));
  vars.add(RooArgList(Lambda0_PZ,Lambda0_PE));
  //Jpsi stuff
  vars.add(RooArgList(Jpsi_ENDVERTEX_CHI2,Jpsi_ENDVERTEX_NDOF));
  vars.add(RooArgList(Jpsi_M));

  //pid
  vars.add(RooArgList(piminus_ProbNNk,piminus_ProbNNp,piminus_ProbNNpi));
  vars.add(RooArgList(pplus_ProbNNk,pplus_ProbNNp,pplus_ProbNNpi));
  vars.add(RooArgList(piminus_PIDK,piminus_PIDp));
  vars.add(RooArgList(pplus_PIDK,pplus_PIDp));


  vars.add(RooArgList(pplus_P,piminus_P));
  vars.add(RooArgList(pplus_PT,piminus_PT));
  vars.add(RooArgList(piminus_PZ,piminus_PE));

  vars.add(RooArgList(Lambdab_y,Lambda0_y,piminus_y,pplus_y));

  RooArgList varstoplot(costheta,costheta1,costheta2,phi1,phi2);
  //  varstoplot.add(RooArgList(F1,F2,F3,F4,F5,F6,F7));
  varstoplot.add(RooArgList(Lambdab_BDIRA,Lambdab_IPCHI2,Lambdab_FDCHI2,Lambdab_eta, Lambdab_PT));// Lambdab_PT,Lambdab_P
  varstoplot.add(RooArgList(Lambda0_PT,Lambdab_TAU,Lambda0_TAU,Lambda0_ENDVERTEX_Z,Lambda0_P,Lambda0_phi,Lambda0_eta));
  varstoplot.add(RooArgList(Jpsi_M));
  varstoplot.add(RooArgList(piminus_PIDK,piminus_PIDp));
  varstoplot.add(RooArgList(pplus_PIDK,pplus_PIDp));

  varstoplot.add(RooArgList(pplus_P,piminus_P));
  varstoplot.add(RooArgList(pplus_PT,piminus_PT));

  varstoplot.add(RooArgList(Lambdab_y,Lambda0_y,piminus_y,pplus_y));


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

  // RooDataSet* dataddup[directorylist.size()];
  // RooDataSet* datallup[directorylist.size()];
  // RooDataSet* datadddown[directorylist.size()];
  // RooDataSet* datalldown[directorylist.size()];

  // RooDataSet* mcddup[directorylist.size()];
  // RooDataSet* mcllup[directorylist.size()];
  // RooDataSet* mcdddown[directorylist.size()];
  // RooDataSet* mclldown[directorylist.size()];

  RooDataSet* datadd[directorylist.size()];
  RooDataSet* mcdd[directorylist.size()];
  RooDataSet* datall[directorylist.size()];
  RooDataSet* mcll[directorylist.size()];

  // RooDataSet* datadd_Lambdab0[directorylist.size()];
  // RooDataSet* mcdd_Lambdab0[directorylist.size()];
  // RooDataSet* datall_Lambdab0[directorylist.size()];
  // RooDataSet* mcll_Lambdab0[directorylist.size()];
  // RooDataSet* datadd_Lambdab0bar[directorylist.size()];
  // RooDataSet* mcdd_Lambdab0bar[directorylist.size()];
  // RooDataSet* datall_Lambdab0bar[directorylist.size()];
  // RooDataSet* mcll_Lambdab0bar[directorylist.size()];

  RooDataSet* data_uppersideband_dd[directorylist.size()];
  RooDataSet* data_uppersideband_ll[directorylist.size()];

//   RooDataSet* data_lowersideband_dd[directorylist.size()];
//   RooDataSet* data_lowersideband_ll[directorylist.size()];


  memset(mcdata,0,sizeof(RooDataSet*)*directorylist.size());
  memset(data,0,sizeof(RooDataSet*)*directorylist.size());


  RooRealVar* Lambdab_ENDVERTEX_CHI2NDOF_rrv = NULL;
  RooRealVar* Jpsi_ENDVERTEX_CHI2NDOF_rrv = NULL;
  RooRealVar* Lambda0_ENDVERTEX_CHI2NDOF_rrv = NULL;
  RooRealVar* Lambda0_DM_rrv = NULL;
  RooRealVar* Jpsi_DMM_rrv = NULL;
  RooRealVar* Lambdab_BDIRALOG_rrv = NULL;
  RooRealVar* Lambdab_FDCHI2LOG_rrv = NULL;

  RooRealVar* piminus_kpi_rrv = NULL;
  RooRealVar* piminus_ppi_rrv = NULL;

  RooRealVar* pplus_pk_rrv = NULL;
  //  RooRealVar* pplus_pip_rrv = NULL;
  
  RooRealVar* Lambdab_PT_rrv = NULL;
  RooRealVar* Lambdab_P_rrv = NULL;

  const bool plotalldata=false;
  if (plotalldata) {
    // TCut simplecut_ll,simplecut_dd;
    // TCut bdtcut_ll,bdtcut_dd;

    // simplecut(doBd2JpsiKS0, simplecut_ll, simplecut_dd);
    // bdtsel(doBd2JpsiKS0,bdtcut_ll, bdtcut_dd);
    // RooDataSet* mcdd_bdt = filldataset("mcdd_bdt" , "mcdd_bdt", RooArgList(costheta,costheta1,costheta2), mctree[0], bdtcut_dd + trueMC + TriggerCut + Lambda0_DD);
    // RooDataSet* mcdd_simple = filldataset("mcdd_simple" , "mcdd_simple", RooArgList(costheta,costheta1,costheta2), mctree[0], simplecut_dd + trueMC + TriggerCut + Lambda0_DD);

    // plotdatavsdata(c,RooArgList(costheta,costheta1,costheta2), mcdd_simple,mcdd_bdt,TString("compmc-dd-bdtvssimple-"), false, false, false);
    // return 0;
    
    TCut sbmasscut(Lambdabname + "_M>5670.&&" + Lambdabname + "_M<5750.");
    if (doBd2JpsiKS0) sbmasscut=TCut(Lambdabname + "_M>5400.");

    data_uppersideband_dd[0] = filldataset("datausdd" , "datausdd", vars, datatree[0], Lambda0_DD + TriggerCut + sbmasscut); 
    data_uppersideband_ll[0] = filldataset("datausll" , "datausll", vars, datatree[0], Lambda0_LL + TriggerCut + sbmasscut);
    
    RooDataSet* mcnocutll = filldataset("mcnocutll" , "mcnocutll", vars, mctree[0], trueMC + TriggerCut + Lambda0_LL);// + TCut(Lambdabname + "_M>5670.&&" + Lambdabname + "_M<5750.")); 
    RooDataSet* mcnocutdd = filldataset("mcnocutdd" , "mcnocutdd", vars, mctree[0], trueMC + TriggerCut + Lambda0_DD);// + TCut(Lambdabname + "_M>5670.&&" + Lambdabname + "_M<5750.")); 

    Lambdab_ENDVERTEX_CHI2NDOF_rrv = (RooRealVar*)data_uppersideband_dd[0]->addColumn(Lambdab_ENDVERTEX_CHI2NDOF);
    Jpsi_ENDVERTEX_CHI2NDOF_rrv = (RooRealVar*)data_uppersideband_dd[0]->addColumn(Jpsi_ENDVERTEX_CHI2NDOF);
    Lambda0_ENDVERTEX_CHI2NDOF_rrv = (RooRealVar*)data_uppersideband_dd[0]->addColumn(Lambda0_ENDVERTEX_CHI2NDOF);
    Lambda0_DM_rrv = (RooRealVar*)data_uppersideband_dd[0]->addColumn(Lambda0_DM);
    Jpsi_DMM_rrv = (RooRealVar*)data_uppersideband_dd[0]->addColumn(Jpsi_DMM);
    Lambdab_BDIRALOG_rrv = (RooRealVar*)data_uppersideband_dd[0]->addColumn(Lambdab_BDIRALOG);
    Lambdab_FDCHI2LOG_rrv = (RooRealVar*)data_uppersideband_dd[0]->addColumn(Lambdab_FDCHI2LOG);
    piminus_kpi_rrv = (RooRealVar*)data_uppersideband_dd[0]->addColumn(piminus_kpi);
    piminus_ppi_rrv = (RooRealVar*)data_uppersideband_dd[0]->addColumn(piminus_ppi);
    pplus_pk_rrv = (RooRealVar*)data_uppersideband_dd[0]->addColumn(pplus_pk);
    //    pplus_pip_rrv = (RooRealVar*)data_uppersideband_dd[0]->addColumn(pplus_pip);

    data_uppersideband_ll[0]->addColumn(Lambdab_ENDVERTEX_CHI2NDOF);
    data_uppersideband_ll[0]->addColumn(Jpsi_ENDVERTEX_CHI2NDOF);
    data_uppersideband_ll[0]->addColumn(Lambda0_ENDVERTEX_CHI2NDOF);
    data_uppersideband_ll[0]->addColumn(Lambda0_DM);
    data_uppersideband_ll[0]->addColumn(Jpsi_DMM);
    data_uppersideband_ll[0]->addColumn(Lambdab_BDIRALOG);
    data_uppersideband_ll[0]->addColumn(Lambdab_FDCHI2LOG);
    data_uppersideband_ll[0]->addColumn(piminus_kpi);
    data_uppersideband_ll[0]->addColumn(piminus_ppi);
    data_uppersideband_ll[0]->addColumn(pplus_pk);
    //    data_uppersideband_ll[0]->addColumn(pplus_pip);
    data_uppersideband_ll[0]->addColumn(Lambdab_y);
    data_uppersideband_ll[0]->addColumn(Lambda0_y);

    mcnocutll->addColumn(Lambdab_ENDVERTEX_CHI2NDOF);
    mcnocutll->addColumn(Jpsi_ENDVERTEX_CHI2NDOF);
    mcnocutll->addColumn(Lambda0_ENDVERTEX_CHI2NDOF);
    mcnocutll->addColumn(Lambda0_DM);
    mcnocutll->addColumn(Jpsi_DMM);
    mcnocutll->addColumn(Lambdab_BDIRALOG);
    mcnocutll->addColumn(Lambdab_FDCHI2LOG);
    mcnocutll->addColumn(piminus_kpi);
    mcnocutll->addColumn(piminus_ppi);
    mcnocutll->addColumn(pplus_pk);
    //    mcnocutll->addColumn(pplus_pip);
    mcnocutll->addColumn(Lambdab_y);
    mcnocutll->addColumn(Lambda0_y);

    mcnocutdd->addColumn(Lambdab_ENDVERTEX_CHI2NDOF);
    mcnocutdd->addColumn(Jpsi_ENDVERTEX_CHI2NDOF);
    mcnocutdd->addColumn(Lambda0_ENDVERTEX_CHI2NDOF);
    mcnocutdd->addColumn(Lambda0_DM);
    mcnocutdd->addColumn(Jpsi_DMM);
    mcnocutdd->addColumn(Lambdab_BDIRALOG);
    mcnocutdd->addColumn(Lambdab_FDCHI2LOG);
    mcnocutdd->addColumn(piminus_kpi);
    mcnocutdd->addColumn(piminus_ppi);
    mcnocutdd->addColumn(pplus_pk);
    //    mcnocutdd->addColumn(pplus_pip);
    mcnocutdd->addColumn(Lambdab_y);
    mcnocutdd->addColumn(Lambda0_y);

    Lambdab_ENDVERTEX_CHI2NDOF_rrv->setMin(0.);
    Lambdab_ENDVERTEX_CHI2NDOF_rrv->setMax(10.);
    Lambda0_ENDVERTEX_CHI2NDOF_rrv->setMin(0.);
    Lambda0_ENDVERTEX_CHI2NDOF_rrv->setMax(10.);
    Jpsi_ENDVERTEX_CHI2NDOF_rrv->setMin(0.);
    Jpsi_ENDVERTEX_CHI2NDOF_rrv->setMax(10.);
    Lambdab_ENDVERTEX_CHI2NDOF_rrv->setBins(50);
    Lambda0_ENDVERTEX_CHI2NDOF_rrv->setBins(50);
    Jpsi_ENDVERTEX_CHI2NDOF_rrv->setBins(50);

    Lambda0_DM_rrv->setMin(0.);
    Lambda0_DM_rrv->setMax(15.);
    if (doBd2JpsiKS0) Lambda0_DM_rrv->setMax(30.);
    Lambda0_DM_rrv->setBins(50);

    Jpsi_DMM_rrv->setMin(0.);
    Jpsi_DMM_rrv->setMax(100.);
    Jpsi_DMM_rrv->setBins(50);

    varstoplot.add(RooArgList(*Lambdab_ENDVERTEX_CHI2NDOF_rrv,*Jpsi_ENDVERTEX_CHI2NDOF_rrv,*Lambda0_ENDVERTEX_CHI2NDOF_rrv));
    varstoplot.add(*Lambda0_DM_rrv);
    varstoplot.add(*Jpsi_DMM_rrv);
    varstoplot.add(*Lambdab_BDIRALOG_rrv);
    varstoplot.add(*Lambdab_FDCHI2LOG_rrv);

    varstoplot.add(*piminus_kpi_rrv);
    varstoplot.add(*piminus_ppi_rrv);
    varstoplot.add(*pplus_pk_rrv);
    //    varstoplot.add(*pplus_pip_rrv);

    std::cout << "mcnocutll->numEntries() = " << mcnocutll->numEntries() << std::endl;
    std::cout << "mcnocutdd->numEntries() = " << mcnocutdd->numEntries() << std::endl;


    Jpsi_M.setMin(3000.);
    Jpsi_M.setMax(3200.);
    Jpsi_M.setBins(100.);
  

    Lambdab_BDIRA.setMin(0.999);
    Lambdab_BDIRA.setMax(1.0);
    Lambdab_IPCHI2.setMin(0.);
    Lambdab_IPCHI2.setMax(15.);
    Lambdab_FDCHI2.setMin(0.);
    Lambdab_FDCHI2.setMax(10000.);
    Lambdab_PT.setMin(0.);
    Lambdab_PT.setMax(25000.);
    if (doBd2JpsiKS0)   Lambdab_PT.setMax(30000.);
    Lambdab_P.setMin(0.);
    Lambdab_P.setMax(300000.);
    Lambdab_TAU.setMin(0.);
    Lambdab_TAU.setMax(0.005);
    Lambdab_eta.setMin(2.);
    Lambdab_eta.setMax(7.);


    Lambda0_PT.setMin(500.);
    Lambda0_PT.setMax(8000.);
    Lambda0_TAU.setMin(0.);
    Lambda0_TAU.setMax(0.5);
    if (doBd2JpsiKS0) Lambda0_TAU.setMax(0.1);
    Lambda0_P.setMin(0.);
    Lambda0_P.setMax(200000.);
    Lambda0_ENDVERTEX_Z.setMin(-400.);
    Lambda0_ENDVERTEX_Z.setMax(3000.);
    Lambda0_phi.setMin(-1.);
    Lambda0_phi.setMax(1.);
    Lambda0_eta.setMin(1.);
    Lambda0_eta.setMax(6.);

    bdt.setMin(0.);
    bdt.setMax(1.);
  
    Lambdab_BDIRA.setBins(50);
    Lambdab_IPCHI2.setBins(50);
    Lambdab_FDCHI2.setBins(100);
    Lambdab_PT.setBins(50);
    Lambdab_P.setBins(50);
    Lambdab_TAU.setBins(50);
    Lambdab_eta.setBins(50);

    Lambda0_PT.setBins(50);
    Lambda0_TAU.setBins(50);
    Lambda0_ENDVERTEX_Z.setBins(50);
    Lambda0_P.setBins(50);
    Lambda0_phi.setBins(50);
    Lambda0_eta.setBins(50);

    bdt.setBins(50);

    Lambdab_BDIRALOG_rrv->setMin(-8.);
    Lambdab_BDIRALOG_rrv->setMax(0.);
    Lambdab_BDIRALOG_rrv->setBins(50);

    Lambdab_FDCHI2LOG_rrv->setMin(0.);
    Lambdab_FDCHI2LOG_rrv->setMax(6.);
    Lambdab_FDCHI2LOG_rrv->setBins(50);

    piminus_kpi_rrv->setMin(-1.);
    piminus_kpi_rrv->setMax(1.);
    piminus_kpi_rrv->setBins(50);
    piminus_ppi_rrv->setMin(-1.);
    piminus_ppi_rrv->setMax(1.);
    piminus_ppi_rrv->setBins(50);

    // pplus_pk_rrv->setMin(-1.);
    // pplus_pk_rrv->setMax(1.);
    // pplus_pk_rrv->setBins(50);
    //    pplus_pip_rrv->setMin(-1.);
    //    pplus_pip_rrv->setMax(1.);
    //    pplus_pip_rrv->setBins(50);

    pplus_PIDK.setMin(-50.); pplus_PIDK.setMax(50.); pplus_PIDK.setBins(50);
    pplus_PIDp.setMin(-50.); pplus_PIDp.setMax(50.); pplus_PIDp.setBins(50);
    piminus_PIDK.setMin(-50.); piminus_PIDK.setMax(50.); piminus_PIDK.setBins(50);
    piminus_PIDp.setMin(-50.); piminus_PIDp.setMax(50.); piminus_PIDp.setBins(50);

    Lambdab_y.setMin(2.);
    Lambdab_y.setMax(7.);
    Lambdab_y.setBins(50);

    Lambda0_y.setMin(2.);
    Lambda0_y.setMax(7.);
    Lambda0_y.setBins(50);

    std::cout << "mcnocutll->numEntries() = " << mcnocutll->numEntries() << std::endl;
    std::cout << "mcnocutdd->numEntries() = " << mcnocutdd->numEntries() << std::endl;

    std::map<RooAbsReal*, double> cutsdd;
    cutsdd[&Lambdab_BDIRA]=0.9999;
    cutsdd[Lambdab_BDIRALOG_rrv]=log10(1.0000001-0.9999);
    cutsdd[Lambdab_FDCHI2LOG_rrv]=log10(180.);
    cutsdd[&Lambdab_IPCHI2]=7.;
    cutsdd[&Lambdab_FDCHI2]=180.;
    cutsdd[Lambdab_ENDVERTEX_CHI2NDOF_rrv]=4.;
    cutsdd[Lambda0_DM_rrv]=15.;
    cutsdd[&Lambda0_PT]=1000.;
    cutsdd[Jpsi_ENDVERTEX_CHI2NDOF_rrv]=8.;
    cutsdd[&Lambdab_TAU]=0.0002;
    cutsdd[Jpsi_DMM_rrv]=40.;

    std::map<RooAbsReal*, double> cutsll;
    cutsll[&Lambdab_BDIRA]=0.9999;
    cutsll[Lambdab_BDIRALOG_rrv]=log10(1.0000001-0.9999);
    cutsll[Lambdab_FDCHI2LOG_rrv]=log10(130.);
    cutsll[&Lambdab_IPCHI2]=8.;
    cutsll[&Lambdab_FDCHI2]=130.;
    cutsll[Lambdab_ENDVERTEX_CHI2NDOF_rrv]=5.;
    cutsll[Lambda0_DM_rrv]=20.;
    cutsll[&Lambda0_PT]=800.;
    cutsll[Jpsi_ENDVERTEX_CHI2NDOF_rrv]=9.;
    cutsll[&Lambdab_TAU]=0.0002;
    cutsll[Jpsi_DMM_rrv]=40.;

    TString plotprefix("compmcdata-");
    if (doBd2JpsiKS0) {plotprefix+="Bd2JpsiKS0";}else {plotprefix+="Lb2JpsiL0";}
    plotprefix+="-signalvsbkg";
      
    for (int i=0;i<varstoplot.getSize();i++){
      RooRealVar* var = (RooRealVar*)varstoplot.at(i);
      const TString varname=var->GetName();
      if (varname.Contains("_PT") or varname.Contains("_FDCHI2")) c.SetLogx();
      if (varname.Contains("_FDCHI2LOG")) c.SetLogx(0);
      double cutvaldd = -666.;
      if (cutsdd.find(var) != cutsdd.end()) cutvaldd = cutsdd[var];
      double cutvalll = -666.;
      if (cutsll.find(var) != cutsll.end()) cutvalll = cutsll[var];

      bool gp=true;
      if (varname.Contains("DIRALOG") or varname.Contains("_DM") or  varname.Contains("CHI2NDOF") or varname.Contains("IPCHI2")) gp=false;
      
      plotdatavsdata(c,*var, mcnocutdd, data_uppersideband_dd[0], plotprefix + "-dd-" + var->GetName() + ".eps", false, false, false,cutvaldd,gp);
      plotdatavsdata(c,*var, mcnocutll, data_uppersideband_ll[0], plotprefix + "-ll-"+ var->GetName() + ".eps", false, false, false,cutvalll,gp);
      c.SetLogx(0);
    }
    //    plotdatavsdata(c,varstoplot, mcnocutdd, data_uppersideband_dd[0],"compmcdata-dd-signalvsbkg", false, false, false);
    //    plotdatavsdata(c,varstoplot, mcnocutll, data_uppersideband_ll[0],"compmcdata-ll-signalvsbkg", false, false, false);
    
    //     c.SetLogx();
    //     RooArgList varstoplotlogx(Lambda0_PT,Lambda0_P,Lambdab_FDCHI2,Lambdab_P,Lambdab_PT);
    //     plotdatavsdata(c,varstoplotlogx, mcnocutdd, data_uppersideband_dd[0],"compmcdata-dd-signalvsbkg", false, false, false);
    //     plotdatavsdata(c,varstoplotlogx, mcnocutll, data_uppersideband_ll[0],"compmcdata-ll-signalvsbkg", false, false, false);
    
    return 0;
  }



  for (unsigned int i=0; i<directorylist.size(); i++) {
    if (!datatree[i]) continue;
    RooDataSet* dataorig = filldataset("data-" + directorylist[i],"data-"+ directorylist[i], vars, datatree[i], selection );
    
    data[i] = dataorig;
    data[i]->Print("v");
 
    //  mygetchar;

    Lambdab_ENDVERTEX_CHI2NDOF_rrv = (RooRealVar*)data[i]->addColumn(Lambdab_ENDVERTEX_CHI2NDOF);
    Jpsi_ENDVERTEX_CHI2NDOF_rrv = (RooRealVar*)data[i]->addColumn(Jpsi_ENDVERTEX_CHI2NDOF);
    Lambda0_ENDVERTEX_CHI2NDOF_rrv = (RooRealVar*)data[i]->addColumn(Lambda0_ENDVERTEX_CHI2NDOF);
    Lambdab_FDCHI2LOG_rrv = (RooRealVar*)data[i]->addColumn(Lambdab_FDCHI2LOG);
    Lambdab_PT_rrv = (RooRealVar*)data[i]->addColumn(Lambdab_PT_rfv);
    Lambdab_P_rrv = (RooRealVar*)data[i]->addColumn(Lambdab_P_rfv);
    pplus_pk_rrv =  (RooRealVar*)data[i]->addColumn(pplus_pk);


    Lambdab_ENDVERTEX_CHI2NDOF_rrv->setMin(0.);
    Lambdab_ENDVERTEX_CHI2NDOF_rrv->setMax(10.);
    Lambda0_ENDVERTEX_CHI2NDOF_rrv->setMin(0.);
    Lambda0_ENDVERTEX_CHI2NDOF_rrv->setMax(10.);
    Jpsi_ENDVERTEX_CHI2NDOF_rrv->setMin(0.);
    Jpsi_ENDVERTEX_CHI2NDOF_rrv->setMax(10.);
    Lambdab_ENDVERTEX_CHI2NDOF_rrv->setBins(25);
    Lambda0_ENDVERTEX_CHI2NDOF_rrv->setBins(25);
    Jpsi_ENDVERTEX_CHI2NDOF_rrv->setBins(25);
    
    Lambdab_FDCHI2LOG_rrv->setMin(0.);
    Lambdab_FDCHI2LOG_rrv->setMax(6.);
    Lambdab_FDCHI2LOG_rrv->setBins(50);

    Lambdab_PT_rrv->setMin(0.);
    Lambdab_PT_rrv->setMax(25.);
    Lambdab_PT_rrv->setBins(25);

    Lambdab_P_rrv->setMin(0.);
    Lambdab_P_rrv->setMax(300.);
    Lambdab_P_rrv->setBins(25);



    if (i==0) {
      varstoplot.add(RooArgList(*Lambdab_ENDVERTEX_CHI2NDOF_rrv,*Jpsi_ENDVERTEX_CHI2NDOF_rrv,*Lambda0_ENDVERTEX_CHI2NDOF_rrv,*Lambdab_FDCHI2LOG_rrv,*Lambdab_PT_rrv,*Lambdab_P_rrv));
      varstoplot.add(*pplus_pk_rrv);
    }
  }

  for (unsigned int i=0; i<directorylist.size(); i++) {
    mcdata[i]=filldataset("mcdata-" + directorylist[i],"mcdata-" + directorylist[i], vars, mctree[i], selection + additionalmcsel);
    mcdata[i]->addColumn(Lambdab_ENDVERTEX_CHI2NDOF);
    mcdata[i]->addColumn(Jpsi_ENDVERTEX_CHI2NDOF);
    mcdata[i]->addColumn(Lambda0_ENDVERTEX_CHI2NDOF);
    mcdata[i]->addColumn(Lambdab_FDCHI2LOG);   
    mcdata[i]->addColumn(Lambdab_PT_rfv);
    mcdata[i]->addColumn(Lambdab_P_rfv);
    mcdata[i]->addColumn(pplus_pk);
    mcdata[i]->Print("v");
  }

  //  std::cout << "Avant: " << data[0]->numEntries() << std::endl;

  Jpsi_M.setMin(0.);
  Jpsi_M.setMax(5500.);
  Jpsi_M.setBins(100.);
  

  Lambdab_BDIRA.setMin(0.9999);
  Lambdab_BDIRA.setMax(1.0);
  Lambdab_IPCHI2.setMin(0.);
  Lambdab_IPCHI2.setMax(15.);
  Lambdab_FDCHI2.setMin(100.);
  Lambdab_FDCHI2.setMax(100000.);
  Lambdab_PT.setMin(0.);
  Lambdab_PT.setMax(25000.);
  if (doBd2JpsiKS0)   Lambdab_PT.setMax(30000.);
  Lambdab_P.setMin(0.);
  Lambdab_P.setMax(300000.);
  Lambdab_TAU.setMin(0.);
  Lambdab_TAU.setMax(0.01);
  Lambdab_eta.setMin(2.);
  Lambdab_eta.setMax(7.);

  Lambda0_PT.setMin(0.);
  Lambda0_PT.setMax(8000.);
  Lambda0_TAU.setMin(0.);
  Lambda0_TAU.setMax(1.);
  if (doBd2JpsiKS0) Lambda0_TAU.setMax(0.1);
  Lambda0_P.setMin(0.);
  Lambda0_P.setMax(200000.);
  Lambda0_ENDVERTEX_Z.setMin(-400.);
  Lambda0_ENDVERTEX_Z.setMax(3000.);
  Lambda0_phi.setMin(-1.);
  Lambda0_phi.setMax(1.);
  Lambda0_eta.setMin(1.);
  Lambda0_eta.setMax(6.);

  bdt.setMin(0.);
  bdt.setMax(1.);
  
  piminus_P.setMin(2000.);
  piminus_P.setMax(20000.);
  pplus_P.setMin(2000.);
  pplus_P.setMax(100000.);
  if (doBd2JpsiKS0) piminus_P.setMax(100000.);


  piminus_PT.setMin(0.);
  piminus_PT.setMax(2000.);
  pplus_PT.setMin(0.);
  pplus_PT.setMax(7000.);

  if (doBd2JpsiKS0) {
    piminus_PT.setMin(0.);
    piminus_PT.setMax(4000.);
    pplus_PT.setMin(0.);
    pplus_PT.setMax(4000.);
  }

  Lambdab_BDIRA.setBins(50);
  Lambdab_IPCHI2.setBins(50);
  Lambdab_FDCHI2.setBins(50);
  Lambdab_PT.setBins(50);
  Lambdab_P.setBins(50);
  Lambdab_TAU.setBins(50);
  Lambdab_eta.setBins(50);

  Lambda0_PT.setBins(50);
  Lambda0_TAU.setBins(50);
  Lambda0_ENDVERTEX_Z.setBins(50);
  Lambda0_P.setBins(50);
  Lambda0_phi.setBins(50);
  Lambda0_eta.setBins(50);

  bdt.setBins(50);

  piminus_P.setBins(50);
  pplus_P.setBins(50);

  
  Lambdab_y.setMin(1.5);
  Lambdab_y.setMax(5.);
  Lambdab_y.setBins(50);
  
  Lambda0_y.setMin(1.5);
  Lambda0_y.setMax(5.);
  Lambda0_y.setBins(50);
  
  piminus_y.setMin(1.5);
  piminus_y.setMax(5.);
  piminus_y.setBins(50);

  pplus_y.setMin(1.5);
  pplus_y.setMax(5.);
  pplus_y.setBins(50);
  
  pplus_PIDK.setMin(-50.); pplus_PIDK.setMax(50.); pplus_PIDK.setBins(50);
  pplus_PIDp.setMin(-10.); pplus_PIDp.setMax(100.); pplus_PIDp.setBins(50);
  piminus_PIDK.setMin(-50.); piminus_PIDK.setMax(50.); piminus_PIDK.setBins(50);
  piminus_PIDp.setMin(-50.); piminus_PIDp.setMax(50.); piminus_PIDp.setBins(50);

  pplus_pk_rrv->setMin(-50.); pplus_pk_rrv->setMax(100.); pplus_pk_rrv->setBins(50);

  //  std::cout << "After: " << data[0]->numEntries() << std::endl;
  //  return 1;


  for (unsigned int i=0; i<directorylist.size(); i++) {
    if (!mcdata[i]) continue;
    mcdd[i] = (RooDataSet*)mcdata[i]->reduce(Lambda0_DD); 
    mcll[i] = (RooDataSet*)mcdata[i]->reduce(Lambda0_LL); 

    //    mcddup[i] = (RooDataSet*)mcdata[i]->reduce(MagnetUp + Lambda0_DD );
    //    mcllup[i] = (RooDataSet*)mcdata[i]->reduce(MagnetUp + Lambda0_LL );
    //    mcdddown[i] = (RooDataSet*)mcdata[i]->reduce(MagnetDown + Lambda0_DD );
    //    mclldown[i] = (RooDataSet*)mcdata[i]->reduce(MagnetDown + Lambda0_LL );
    
    //    mcdd_Lambdab0[i] = (RooDataSet*)mcdata[i]->reduce(isLambdab0 + Lambda0_DD); 
    //    mcll_Lambdab0[i] = (RooDataSet*)mcdata[i]->reduce(isLambdab0 + Lambda0_LL); 
    //    mcdd_Lambdab0bar[i] = (RooDataSet*)mcdata[i]->reduce(isLambdab0bar + Lambda0_DD); 
    //    mcll_Lambdab0bar[i] = (RooDataSet*)mcdata[i]->reduce(isLambdab0bar + Lambda0_LL); 
  }
    
    
  for (unsigned int i=0; i<directorylist.size(); i++) {
    if (!data[i]) continue;
    datadd[i] = (RooDataSet*)data[i]->reduce(Lambda0_DD); 
    datall[i] = (RooDataSet*)data[i]->reduce(Lambda0_LL); 

    // dataddup[i] = (RooDataSet*)data[i]->reduce(MagnetUp + Lambda0_DD );
    // datallup[i] = (RooDataSet*)data[i]->reduce(MagnetUp + Lambda0_LL );
    // datadddown[i] = (RooDataSet*)data[i]->reduce(MagnetDown + Lambda0_DD );
    // datalldown[i] = (RooDataSet*)data[i]->reduce(MagnetDown + Lambda0_LL );
     
    // datadd_Lambdab0[i] = (RooDataSet*)data[i]->reduce(isLambdab0 + Lambda0_DD); 
    // datall_Lambdab0[i] = (RooDataSet*)data[i]->reduce(isLambdab0 + Lambda0_LL); 
    // datadd_Lambdab0bar[i] = (RooDataSet*)data[i]->reduce(isLambdab0bar + Lambda0_DD); 
    // datall_Lambdab0bar[i] = (RooDataSet*)data[i]->reduce(isLambdab0bar + Lambda0_LL); 
    
    
    
//     data_lowersideband_dd[i] = (RooDataSet*)data[i]->reduce(Lambda0_DD + TCut(Lambdabname + "_M<5570.")); 
    
//     data_lowersideband_ll[i] = (RooDataSet*)data[i]->reduce(Lambda0_LL + TCut(Lambdabname + "_M<5570.")); 
  }


  //  RooArgList varstoplot();  

  
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
  // plotwcomponents(c,RooArgList(mass), tot_mass_model, RooArgList(sig_mass_model,bkg_mass_model), datall, plotprefix + "-blind-ll-massfit");
  // tot_model.fitTo(*datall, RooFit::SumW2Error(true));

  // plotwcomponents(c,RooArgList(mass,costheta,costheta1,costheta2), tot_model, RooArgList(sig_model,bkg_model), datall, plotprefix + "-blind-ll-totfit");
  // plotwcutwcomponents(c,RooArgList(costheta,costheta1,costheta2), tot_model, RooArgList(sig_model,bkg_model), datall, "signalregion", plotprefix + "-blind-ll-totfit");
  // plotwcutwcomponents(c,RooArgList(costheta,costheta1,costheta2), tot_model, RooArgList(sig_model,bkg_model), datall, "lowersideband",plotprefix + "-blind-ll-totfit");
  // plotwcutwcomponents(c,RooArgList(costheta,costheta1,costheta2), tot_model, RooArgList(sig_model,bkg_model), datall, "uppersideband",plotprefix + "-blind-ll-totfit");



  std::vector<var12> varstoplot2d;

  var12 avar12;
  avar12.var1=&Lambdab_PT;
  avar12.var2=&Lambdab_y;
  varstoplot2d.push_back(avar12);
  avar12.var1=&Lambda0_PT;
  avar12.var2=&Lambda0_y;
  varstoplot2d.push_back(avar12);
  avar12.var1=&piminus_PT;
  avar12.var2=&piminus_y;
  varstoplot2d.push_back(avar12);
  avar12.var1=&pplus_PT;
  avar12.var2=&pplus_y;
  varstoplot2d.push_back(avar12);


  int nbinx,nbiny;
  if (doBd2JpsiKS0) {
    nbinx=15;
    nbiny=4;
  } else {
    if (doLL) {
      nbinx=5;
      nbiny=2;
    } else {
      nbinx=10;
      nbiny=3;
    }
  }
  

  //  return 0;
  for (unsigned int i=0; i<directorylist.size(); i++) {
    TString plotprefix("compmcdata-");
    if (userwmc) {plotprefix+="mcrw-"; } else {plotprefix+="mcnonrw-";}
    if (doBd2JpsiKS0) {plotprefix+="Bd2JpsiKS0";}else {plotprefix+="Lb2JpsiL0";}
    const TString& thisdirectory = directorylist[i];
    plotprefix+="-"+thisdirectory;

    if (applybdtcut) {
      plotprefix+="-bdtapplied";
    }

    if (applytmvarectcut) {
      plotprefix+="-tmvaapplied";
    }

    if (applysimplecut) {
      plotprefix+="-simplecutapplied";
    }

    if (doLL) {
      fitmassandfillsweights( c,tot_ll_mass_model, datall[i], plotprefix + "-ll", true, varstoplot);

      for (unsigned int j=0;j<varstoplot2d.size();j++) {
	compsigsplotmc2d(c, varstoplot2d[j].var1, varstoplot2d[j].var2, datall[i], mcll[i], plotprefix + "-ll-datavsmc",nbinx,nbiny);
      }

      compsigsplotmc(c, varstoplot, datall[i], mcll[i], plotprefix + "-ll-datavsmc");

    } else {
      //      fitwsplot(c, doBd2JpsiKS0, datadd[i], tot_mass_model, RooArgList(sig_mass_model,bkg_mass_model), sig_ang_model, RooArgList(sigma,mean,c0), "fitwsplot-dd");
      //      simplemassfit( c,tot_mass_model, datadd[i], plotprefix + "-data-dd");
      fitmassandfillsweights( c,tot_dd_mass_model, datadd[i], plotprefix + "-dd", true, varstoplot);

      for (unsigned int j=0;j<varstoplot2d.size();j++) {
	compsigsplotmc2d(c, varstoplot2d[j].var1, varstoplot2d[j].var2, datadd[i], mcdd[i], plotprefix + "-dd-datavsmc",nbinx,nbiny);
      }

      compsigsplotmc(c, varstoplot, datadd[i], mcdd[i], plotprefix + "-dd-datavsmc");

    }
  }
}    






























      // if (usesimplecuts) {
      // 	sigma2.setVal(sigparamsdd.getRealValue("sigma2"));
      // 	frac.setVal(sigparamsdd.getRealValue("frac"));
      // 	fit(        c,tot_mass_model, datadd[i], mcdd[i], plotprefix + "-dd", true, splotconstvars,varstoplot, false, angfitvars);
      // 	sigma2.setVal(sigparamsll.getRealValue("sigma2"));
      // 	frac.setVal(sigparamsll.getRealValue("frac"));
      // 	fit(        c,tot_mass_model, datall[i], mcll[i], plotprefix + "-ll", true, splotconstvars,varstoplot, false, angfitvars);
	
      // } else {
	
      // 	P_b=0.;r_0=0.5;r_1=0.;
      // 	mygetchar;
      // 	sigma2.setVal(sigparamsdd.getRealValue("sigma2"));
      // 	frac.setVal(sigparamsdd.getRealValue("frac"));
      // 	fit(        c,tot_mass_model, datadd_Lambdab0[i], mcdd_Lambdab0[i], plotprefix + "-dd-Lambdab0", true, splotconstvars,varstoplot, false, angfitvars);
      // 	//    fitbinbybin(c,costheta,tot_mass_model, datadd_Lambdab0, plotprefix +"-dd-Lambdab0-costhetabinbybin", splotconstvars);
	
	
      // 	P_b=0.;r_0=0.5;r_1=0.;
      // 	sigma2.setVal(sigparamsdd.getRealValue("sigma2"));
      // 	frac.setVal(sigparamsdd.getRealValue("frac"));
      // 	fit(        c,tot_mass_model, datadd_Lambdab0bar[i], mcdd_Lambdab0bar[i], plotprefix + "-dd-Lambdab0bar", true,splotconstvars, varstoplot, false, angfitvars);
      // 	//    fitbinbybin(c,costheta,tot_mass_model, datadd_Lambdab0bar, plotprefix + "-dd-Lambdab0bar-costhetabinbybin", splotconstvars);
	
      // 	P_b=0.;r_0=0.5;r_1=0.;
      // 	sigma2.setVal(sigparamsll.getRealValue("sigma2"));
      // 	frac.setVal(sigparamsll.getRealValue("frac"));
      // 	fit(c,        tot_mass_model, datall_Lambdab0[i], mcll_Lambdab0[i],plotprefix + "-ll-Lambdab0", true, splotconstvars, varstoplot, false, angfitvars);
      // 	//    fitbinbybin(c,costheta,tot_mass_model, datall_Lambdab0, plotprefix + "-ll-Lambdab0-costhetabinbybin", splotconstvars);
	
      // 	P_b=0.;r_0=0.5;r_1=0.;
      // 	sigma2.setVal(sigparamsll.getRealValue("sigma2"));
      // 	frac.setVal(sigparamsll.getRealValue("frac"));
      // 	fit(        c, tot_mass_model, datall_Lambdab0bar[i], mcll_Lambdab0bar[i], plotprefix + "-ll-Lambdab0bar",true, splotconstvars,varstoplot,  false, angfitvars);
      // 	//    fitbinbybin(c, costheta, tot_mass_model, datall_Lambdab0bar, plotprefix + "-ll-Lambdab0bar-costhetabinbybin", splotconstvars);
	
	
      // 	compsigsplot(c,varstoplot,datadd_Lambdab0[i], datadd_Lambdab0bar[i], plotprefix + "-dd-splot-Lambdab0vsLambdab0bar");
      // 	compsigsplot(c,varstoplot,datall_Lambdab0[i], datall_Lambdab0bar[i], plotprefix + "-ll-splot-Lambdab0vsLambdab0bar");
      //    }
      //      P_b=0.;alpha_b=0.;r_0=0.5;r_1=0.;
      //      fit4dim(c, datadd_Lambdab0[i], tot_mass_model, tot_model, RooArgList(sig_model,bkg_model), plotprefix + "-blind-dd-Lambdab0");
      //      P_b=0.;alpha_b=0.;r_0=0.5;r_1=0.;
      //      fit4dim(c, datall_Lambdab0[i], tot_mass_model, tot_model, RooArgList(sig_model,bkg_model), plotprefix + "-blind-ll-Lambdab0");

    
    
    //    return 0;
    
    
    //     fit(c,tot_mass_model, data[i], mcdata[i], plotprefix + "-alldata", true,   RooArgList(frac,c0,c1,mean,sigma,sigma2), 
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
  //  return 0;

  //   if (directorylist.size()>=2) {
  //     compsigsplot(c,varstoplot,datadd[0], datadd[1], "Bd2JpsiKS0-dd-splot-" + directorylist[0] + "_vs_" + directorylist[1]);
  //     compsigsplot(c,varstoplot,datall[0], datall[1], "Bd2JpsiKS0-ll-splot-" + directorylist[0] + "_vs_" + directorylist[1]);
  //   }
  
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







// void compdataandmcwithsplot(TCanvas& c, RooAbsPdf& model, RooDataSet* data, RooDataSet* mcdata, const TString& plotname, bool plotsplot, const RooArgList& splotvars, const RooArgList& splotvarstoplot, bool doangfit, const RooArgList& angfitvars)
// {

//   RooFitResult* frll = model.fitTo(*data,RooFit::Minos(true),RooFit::Save(true),RooFit::Extended(true));
//   if (frll->status()!=0 or frll->covQual()!=3) {
//   std::cout << plotname << " mass fit - not successful (status=" << frll->status() << "; covQual=" << frll->covQual() << ")" << std::endl;
//     mygetchar;
//   }
//   frll->Print("v");
//   delete frll;

//   plot(c,mass,model,data, plotname + "-mass.eps");
//   plotwpull(c,mass,model,data, plotname + "-masswpull.eps");
//   //  plotvar(c,costheta,&data,"signalregion", plotname + "-costheta.eps");
//   //  plotvar(c,costheta1,&data,"signalregion", plotname + "-costheta1.eps");
//   //  plotvar(c,costheta2,&data,"signalregion", plotname + "-costheta2.eps");

//   if (!plotsplot) return;

//   //  int i=0;
//   for (int i=0; i<splotvars.getSize(); i++) {
//     RooRealVar* myvar = (RooRealVar*)splotvars.at(i);
//     myvar->setConstant();
//   }

//   // frac.setConstant();
//   // c0.setConstant();
//   // c1.setConstant();
//   // mean.setConstant();
//   // //  mean2.setConstant();
//   // sigma.setConstant();
//   // sigma2.setConstant();
//   RooStats::SPlot* sData = new RooStats::SPlot("sData-" + plotname,"An SPlot",
// 					       *data, &model, RooArgList(nSig,nBkg) );

//   sData->Print("v");
//   //  mygetchar();

//   RooDataSet data_sig("data_sig","data_sig",data,*data->get(),0,"nSig_sw") ;
//   RooDataSet data_bkg("data_bkg","data_bkg",data,*data->get(),0,"nBkg_sw") ;
  
//   data_sig.get()->Print("v");
//   //  mygetchar();
//   //  data_sig.write(plotname + "-data.dat");

//   //  RooArgList vars(costheta,costheta1,costheta2,phi1,phi2);
//   //  vars.add(RooArgList(Lambdab_BDIRA,Lambdab_IPCHI2,Lambdab_FDCHI2,Lambdab_PT,Lambda0_PT,Lambdab_TAU,Lambda0_TAU));
//   //  vars.add(RooArgList(*ratio1,*ratio2,*ratio3));

//   //  TH1F costhetahisto("costhetahisto","costhetahisto",25,-1.,1.);


// //   const int nbin=10;
// //   TH3F costhetahisto("costhetahisto","costhetahisto",nbin,-1.,1.,nbin,-1.,1.,nbin,-1.,1.);
// //   TH3F mccosthetahisto("mccosthetahisto","mccosthetahisto",nbin,-1.,1.,nbin,-1.,1.,nbin,-1.,1.);
// //   data_sig.fillHistogram(&costhetahisto,RooArgList(costheta,costheta1,costheta2));
// //   mcdata->fillHistogram(&mccosthetahisto,RooArgList(costheta,costheta1,costheta2));
// //   std::cout << "mccosthetahisto entries = " << mccosthetahisto.GetEntries() << std::endl;
// //   //  mygetchar;
// //   costhetahisto.Divide(&mccosthetahisto);

// //   TH1D* his0=costhetahisto.ProjectionX();
// //   his0->SetMinimum(0.);
// //   his0->GetXaxis()->SetTitle(costheta.GetTitle());
// //   his0->Draw();
// //   c.SaveAs(plotname + "-test-costheta.eps");

// //   TH1D* his1=costhetahisto.ProjectionY();
// //   his1->SetMinimum(0.);
// //   his1->GetXaxis()->SetTitle(costheta1.GetTitle());
// //   his1->Draw();
// //   c.SaveAs(plotname + "-test-costheta1.eps");

// //   TH1D* his2=costhetahisto.ProjectionZ();
// //   his2->SetMinimum(0.);
// //   his2->GetXaxis()->SetTitle(costheta2.GetTitle());
// //   his2->Draw();
// //   c.SaveAs(plotname + "-test-costheta2.eps");


//   plotdata(c,splotvarstoplot, &data_sig, plotname + "-splot-signal");
//   if (mcdata) plotdatavsdata(c,splotvarstoplot, &data_sig, mcdata, plotname + "-splot-datavsmc");
//   mygetchar;
//   //  plotdata(c,costheta1,&data_sig, plotname + "-splot-signal-costheta1.eps");
//   //  plotdata(c,costheta2,&data_sig, plotname + "-splot-signal-costheta2.eps");
//   //  plotdata(c,phi1,&data_sig, plotname + "-splot-signal-phi1.eps");
//   //  plotdata(c,phi2,&data_sig, plotname + "-splot-signal-phi2.eps");

//   plotdata(c,splotvarstoplot, &data_bkg, plotname + "-splot-bkg");
//   //  plotdata(c,costheta1,&data_bkg, plotname + "-splot-bkg-costheta1.eps");
//   //  plotdata(c,costheta2,&data_bkg, plotname + "-splot-bkg-costheta2.eps");
//   //  plotdata(c,phi1,&data_bkg, plotname + "-splot-bkg-phi1.eps");
//   //  plotdata(c,phi2,&data_bkg, plotname + "-splot-bkg-phi2.eps");

//   for (int i=0; i<splotvars.getSize(); i++) {
//     RooRealVar* myvar = (RooRealVar*)splotvars.at(i);
//     myvar->setConstant(false);
//   }

//   // frac.setConstant(false);
//   // c0.setConstant(false);
//   // c1.setConstant(false);
//   // mean.setConstant(false);
//   // //  mean2.setConstant(false);
//   // sigma.setConstant(false);
//   // sigma2.setConstant(false);


//   if (!doangfit) return;


//   //  RooFitResult* angfit = w3wacc->fitTo(data_sig,RooFit::Save(true), RooFit::SumW2Error(true));//,RooFit::Minos(true));
//   //  angfit->Print("v");
//   //  mygetchar;

//   //  plot(c, angfitvars, *w3wacc, &data_sig, plotname + "-fitdata-ang");


// }


// void fitbinbybin(TCanvas& c, RooRealVar& var, RooAbsPdf& model, RooDataSet* datadd_Lambdab0, const TString& title, const RooArgList& vars)
// {
//   //  frac.setConstant();
//   //  sigma.setConstant();
//   //  sigma2.setConstant();
//   //  mean.setConstant();
  
//   for (int i=0; i<vars.getSize(); i++) {
//     RooRealVar* myvar = (RooRealVar*)vars.at(i);
//     myvar->setConstant();
//   }


//   const double varmin=var.getMin();
//   const double varmax=var.getMax();
//   const TString varname(var.GetName());
//   const int ndiv=25;

//   // RooCategory masscat("masscat","masscat");
//   // for (int i=0;i<ndiv;i++) {
//   //   //    TString varcut(varname + ">" + str(varmin+(varmax-varmin)/ndiv*i) + "&&" + varname + "<" + str(varmin+(varmax-varmin)/ndiv*(i+1)));
//   //   TString rangename=varname + "-reg" + istr(i);
//   //   var.setRange(rangename,varmin+(varmax-varmin)/ndiv*i, varmin+(varmax-varmin)/ndiv*(i+1));
//   //   masscat.setRange(rangename,rangename);
//   // }
//   // mygetchar;

//   TH1F histo("histo","histo",ndiv,varmin,varmax);
//   for (int i=0;i<ndiv;i++) {
//     TString varcut(varname + ">" + str(varmin+(varmax-varmin)/ndiv*i) + "&&" + 
//                    varname + "<" + str(varmin+(varmax-varmin)/ndiv*(i+1)));
//     std::cout << "costhetacut = " << varcut << std::endl;
//     //    filldataset(datadd_Lambdab0,"datadd_Lambdab0","data", RooArgList(mass,costheta,costheta1,costheta2,phi1,phi2), sigtree,DD_selection + isLambdab0 + costhetacut);
//     RooDataSet* datadd_Lambdab0_red = (RooDataSet*)datadd_Lambdab0->reduce(varcut);
//     fit(c,model,datadd_Lambdab0_red, NULL,   title + "-reg" + istr(i),false,vars,vars,false,vars);
//     std::cout << "nSig = " << nSig.getVal() << " +- " << nSig.getError() << std::endl;
//     histo.SetBinContent(i+1,nSig.getVal());
//     histo.SetBinError(i+1,nSig.getError());
//     delete datadd_Lambdab0_red;
//   }
//   histo.SetMinimum(0.);
//   histo.GetXaxis()->SetTitle(var.GetTitle());
//   histo.Draw("E");
//   c.SaveAs(title + ".eps");

  
//   for (int i=0; i<vars.getSize(); i++) {
//     RooRealVar* myvar = (RooRealVar*)vars.at(i);
//     myvar->setConstant(false);
//   }

//   //  frac.setConstant(false);
//   //  sigma.setConstant(false);
//   //  sigma2.setConstant(false);
//   //  mean.setConstant(false);

// }



// void  fit4dim(TCanvas& c, RooDataSet* datadd, RooAbsPdf& tot_mass_model, RooAbsPdf& tot_model, const RooArgList& totmodelpdfs, 
//               const TString& plotprefix, bool doBd2JpsiKS0)
// {      
//   RooDataSet newdata1("newdata1","newdata1", datadd,*datadd->get());

//   RooRealVar wacc("wacc","wacc",0.);
//   addacceptanceweight(&newdata1, wacc,doBd2JpsiKS0, true);
      
//   RooArgList varswacc(*newdata1.get());
//   varswacc.add(RooArgList(wacc));

//   RooDataSet newdata2("newdata2","newdata2", &newdata1, varswacc, 0, "wacc");

//   tot_mass_model.fitTo(newdata2, RooFit::SumW2Error(true), RooFit::Extended(true));
//   //  plotwcomponents(c,RooArgList(mass), tot_mass_model, massmodelpdfs, datadd, plotname + "-massfit");

//   tot_model.fitTo(newdata2, RooFit::SumW2Error(true), RooFit::Extended(true));//, RooFit::Minos(true));
//   mygetchar;
  
//   plotwcomponents(c,RooArgList(mass,costheta,costheta1,costheta2), tot_model, totmodelpdfs, &newdata2, plotprefix + "-totfit");
//   plotwcutwcomponents(c,RooArgList(costheta,costheta1,costheta2),   tot_model, totmodelpdfs, &newdata2, "signalregion", plotprefix + "-totfit");
//   plotwcutwcomponents(c,RooArgList(costheta,costheta1,costheta2),   tot_model, totmodelpdfs, &newdata2, "lowersideband",plotprefix + "-totfit");
//   plotwcutwcomponents(c,RooArgList(costheta,costheta1,costheta2),   tot_model, totmodelpdfs, &newdata2, "uppersideband",plotprefix + "-totfit");
// }


// void fitwsplot(TCanvas& c, bool doBd2JpsiKS0, RooDataSet* tot_data, RooAbsPdf& tot_mass_model, const RooArgList& masspdfs, RooAbsPdf& sig_ang_model, const RooArgList& splotvars, const TString& plotname)
// {
//   bool inclacc=true;
//   bool inclbkg=true;
//   bool debug=false;
//   bool doplot=true;

//   RooFitResult* fr_mass = tot_mass_model.fitTo(*tot_data, RooFit::Minos(true), RooFit::Save(true));
//   fr_mass->Print("v");

//   if (doplot) plotwcomponents(c,RooArgList(mass), tot_mass_model, masspdfs,tot_data, plotname); 
   
//   for (int i=0; i<splotvars.getSize(); i++) ((RooRealVar*)splotvars.at(i))->setConstant();
//   RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
// 					       *tot_data, &tot_mass_model, RooArgList(nSig,nBkg) );
//   if (debug) {
//     sData->Print("v");
//     mygetchar;
//   }
//   //    RooDataSet data_sig("data_sig","data_sig",tot_data,*tot_data->get());//,0,"nSig_sw") ;
//   RooRealVar wacc("wacc","wacc",0.);
//   RooRealVar wtot("wtot","wtot",0.);
//   if (inclacc) {
//     //    addacceptanceweight(tot_data,wacc,*acceptance,true);  
//     addacceptanceweight(tot_data, wacc,doBd2JpsiKS0, false);
//     multiplyweights(tot_data,wtot,wacc.GetName(),"nSig_sw", nSig.getVal());
//   }
//   if (debug) {
//     tot_data->Print("v");
//     mygetchar;
//   }
//   RooDataSet* tot_data_w=NULL;
//   if (inclbkg and inclacc) tot_data_w = new RooDataSet("dataw","dataw",tot_data,*tot_data->get(),0,"wtot");
//   if (inclbkg and !inclacc) tot_data_w = new RooDataSet("dataw","dataw",tot_data,*tot_data->get(),0,"nSig_sw");
//   if (!inclbkg and inclacc) tot_data_w = new RooDataSet("dataw","dataw",tot_data,*tot_data->get(),0,"wacc");
//   if (!inclbkg and !inclacc) tot_data_w = tot_data;
  
//   if (debug) {
//     tot_data_w->Print("v");
//     mygetchar;
//   }
  
//   RooFitResult* fr = sig_ang_model.fitTo(*tot_data_w,RooFit::SumW2Error(inclbkg or inclacc), RooFit::Save(true), RooFit::Minos(false));
//   fr->Print("v");
//   std::cout << "status=" << fr->status() << " ; covqual=" << fr->covQual() << std::endl; 
  
//   if (doplot) {
//     //    TCanvas c("c","c",100,100);
//     //      plotwpull(c,RooArgList(mass,costheta,costheta1,costheta2), tot_model, dataw, "toy-" + istr(iSample,"04") );
//     plotwpull(c,RooArgList(costheta,costheta1,costheta2), sig_ang_model, tot_data_w, plotname);
//       //      plotwcutwcomponents(c,RooArgList(mass,costheta,costheta1,costheta2), tot_model, RooArgList(sig_model,bkg_model),tot_data_w, "signalregion", "toy-" + istr(iSample,"04") + "-sr" );
//   }

//   for (int i=0; i<splotvars.getSize(); i++) ((RooRealVar*)splotvars.at(i))->setConstant(false);
// }


