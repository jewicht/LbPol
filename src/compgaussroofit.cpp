#include "TFile.h"
#include "RooLbtoJpsiL0PDF3.h"
#include "functions-roofit.h"
#include "rootlogon.h"

//./writedecayfile -0.15 0.1 0.05 0. 0. 0.
//0.74162 0 0.273861 0 0.158114 0 0.591608 0


struct astr {
  double pb;
  double ap;
  double am;
  double bp;
  double bm;
};

int main(int argc, char** argv)
{
  lhcbstyle();

  RooRealVar lb_pt("lb_pt",   "#Lambda_{b}^{0}: p_{T}", 0.,60000.,"MeV/c");
  RooRealVar l0_pt("l0_pt",   "#Lambda: p_{T}",         0.,25000.,"MeV/c");
  RooRealVar pion_pt("pion_pt",   "#pi: p_{T}",         0.,5000.,"MeV/c");

  RooRealVar costheta("costheta",   "cos#theta",     -1., 1.);
  RooRealVar costheta1("costheta1", "cos#theta_{1}", -1., 1.);
  RooRealVar costheta2("costheta2", "cos#theta_{2}", -1., 1.);
  
  RooRealVar P_b("P_b","P_{b}",0.40, -1.5, 1.5);//0.50
  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}",0.642);
  
  RooRealVar alpha_b("alpha_b","#alpha_{b}",-0.15,-1.1,1.1);//0.457
  RooRealVar r_0("r_0","r_{0}",0.1,-0.3,2.3);//0.5
  RooRealVar r_1("r_1","r_{1}",0.05,-1.3,1.3);//0.1
  
  TCanvas c("c","c",100,100);


  std::vector<TString> files;
  files.push_back("Tuple-Lb2L0Jpsi-Polb-HELAMP10000000.root");
  //  files.push_back("Tuple-Lb2L0Jpsi-Polb-HELAMP00000010.root");
  //  files.push_back("Tuple-Lb2L0Jpsi-Polb-HELAMP00001000.root");
  //  files.push_back("Tuple-Lb2L0Jpsi-Polb-HELAMP00100000.root");
 
  //  files.push_back("Tuple-Lb2L0Jpsi-Polb-HELAMP00000010-Pb00.root");
  //  files.push_back("Tuple-Lb2L0Jpsi-Polb-HELAMP00000010-Pb20.root");
  //  files.push_back("Tuple-Lb2L0Jpsi-Polb-HELAMP00000010-Pb40.root");
  //  files.push_back("Tuple-Lb2L0Jpsi-Polb-HELAMP00000010-Pb60.root");
  files.push_back("Tuple-Lb2L0Jpsi-Polb-HELAMP00000010-Pb80.root");
  {
    
    RooDataSet* data[files.size()];
    RooDataSet* data_lbptcutp[files.size()];
    RooDataSet* data_lbptcutm[files.size()];
    RooDataSet* data_l0ptcutp[files.size()];
    RooDataSet* data_l0ptcutm[files.size()];
    for (unsigned int i=0;i<files.size();i++) {
      TFile inputfile(files[i]);
      TDirectory* mydirectory = (TDirectory*)inputfile.Get("BtoqllGeneratedDistributions");
      TTree* mytree = (TTree*)mydirectory->Get("Lb2JpsiL0");
      data[i] =  filldataset("data" + istr(i), "data", RooArgList(costheta,costheta1,costheta2,l0_pt,lb_pt,pion_pt), mytree, "lb_id>0.");
      data_lbptcutm[i] = (RooDataSet*)data[i]->reduce("lb_pt<20000.");
      data_lbptcutp[i] = (RooDataSet*)data[i]->reduce("lb_pt>20000.");
      data_l0ptcutm[i] = (RooDataSet*)data[i]->reduce("l0_pt<5000.");
      data_l0ptcutp[i] = (RooDataSet*)data[i]->reduce("l0_pt>5000.");

      //      std::cout << data[i]->numEntries() << std::endl;
      //      std::cout << data_ptcut[i]->numEntries() << std::endl;

      plotdatavsdata(c, RooArgList(costheta,costheta1,costheta2), data_lbptcutp[i], data_lbptcutm[i], "compgaussroofit-lbptcut" + istr(i), false, true, false);
      plotdatavsdata(c, RooArgList(costheta,costheta1,costheta2), data_l0ptcutp[i], data_l0ptcutm[i], "compgaussroofit-l0ptcut" + istr(i), false, true, false);

    }
    plotprofile(c, "compgaussroofit-profile", data[0], &costheta, &lb_pt);
    plotprofile(c, "compgaussroofit-profile", data[0], &costheta1, &lb_pt);
    plotprofile(c, "compgaussroofit-profile", data[0], &costheta, &l0_pt);
    plotprofile(c, "compgaussroofit-profile", data[0], &costheta1, &l0_pt);

    plotdatavsdata(c, RooArgList(l0_pt,lb_pt,pion_pt), data[0], data[1], "compgaussroofit-ptcomp", false, true, false);
    exit(1);

  }

  RooLbtoJpsiL0PDF3 w3("w3","w3",costheta,costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);
  
  if (1==2) {
    TFile inputfile("Tuple-Lb2L0Jpsi-Polb-PRD65074030.root");
    TDirectory* mydirectory = (TDirectory*)inputfile.Get("BtoqllGeneratedDistributions");
    TTree* mytree = (TTree*)mydirectory->Get("Lb2JpsiL0");
    RooDataSet* data =  filldataset("data", "data", RooArgList(costheta,costheta1,costheta2), mytree, "lb_id>0.");
    
  
    
    RooDataSet* data2 = w3.generate(RooArgList(costheta,costheta1,costheta2), data->numEntries());
    
    plotdatavsdata(c,RooArgList(costheta,costheta1,costheta2),data,data2,"compgaussroofit");
    
    delete data;
    delete data2;
    inputfile.Close();
  }
  
  std::vector<astr> astrv;

  astr a;  
  a.pb=+0.4; a.bm=1.; a.ap=0.; a.am=0.; a.bp=0.; astrv.push_back(a);
  a.pb=+0.2; a.bm=0.; a.ap=1.; a.am=0.; a.bp=0.; astrv.push_back(a);
  a.pb=-0.2; a.bm=0.; a.ap=0.; a.am=1.; a.bp=0.; astrv.push_back(a);
  a.pb=-0.4; a.bm=0.; a.ap=0.; a.am=0.; a.bp=1.; astrv.push_back(a);

  //ap, bm, bp, am 
  for (unsigned int i=0;i<astrv.size();i++) {
    TString helampname="HELAMP" 
      + istr(astrv[i].bp) + "0" 
      + istr(astrv[i].am) + "0"
      + istr(astrv[i].ap) + "0" 
      + istr(astrv[i].bm) + "0"; 

    TComplex bminus(astrv[i].bm,0.,true);
    TComplex aplus( astrv[i].ap,0.,true);
    TComplex aminus(astrv[i].am,0.,true);
    TComplex bplus(astrv[i].bp,0. ,true);

    TFile inputfile("Tuple-Lb2L0Jpsi-Polb-" + helampname +  ".root");
    TDirectory* mydirectory = (TDirectory*)inputfile.Get("BtoqllGeneratedDistributions");
    TTree* mytree = (TTree*)mydirectory->Get("Lb2JpsiL0");
    RooDataSet* data =  filldataset("data", "data", RooArgList(costheta,costheta1,costheta2), mytree, "lb_id>0.");
    
    P_b=0.;
    alpha_b=0.;
    r_0=0.5;
    r_1=0.;

    //    w3.fitTo(*data);
    //    std::cout << helampname << std::endl;
    //    mygetchar;

    P_b=astrv[i].pb;
    alpha_b=  aplus.Rho2() - aminus.Rho2() + bplus.Rho2() - bminus.Rho2();
    r_0 =  aplus.Rho2() + aminus.Rho2();
    r_1 = aplus.Rho2() - aminus.Rho2();
    RooDataSet* data2 = w3.generate(RooArgList(costheta,costheta1,costheta2), data->numEntries());
    
    std::cout << helampname << std::endl;
    std::cout << "Pb = " << P_b.getVal() << std::endl;
    std::cout << "ab = " << alpha_b.getVal() << std::endl;
    std::cout << "r0 = " << r_0.getVal() << std::endl;
    std::cout << "r1 = " << r_1.getVal() << std::endl;
    plotdatavsdata(c,RooArgList(costheta,costheta1,costheta2),data,data2,"compgaussroofit-" + helampname);
    delete data;
    delete data2;
    inputfile.Close();
  }
  
}
