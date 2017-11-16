#include "TFile.h"
#include "RooLbtoJpsiL0PDF3.h"
#include "functions-roofit.h"
#include "rootlogon.h"




int main(int argc, char** argv)
{
  lhcbstyle();
  
  RooRealVar costheta("costheta",   "cos#theta",     -1., 1.);
  RooRealVar costheta1("costheta1", "cos#theta_{1}", -1., 1.);
  RooRealVar costheta2("costheta2", "cos#theta_{2}", -1., 1.);
  RooRealVar P_b("P_b","P_{b}",0.40, -1.5, 1.5);//0.50
  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}",0.642);
  
  RooRealVar alpha_b("alpha_b","#alpha_{b}",-0.15,-1.1,1.1);//0.457
  RooRealVar r_0("r_0","r_{0}",0.1,-0.3,2.3);//0.5
  RooRealVar r_1("r_1","r_{1}",0.05,-1.3,1.3);//0.1
  
  RooFormulaVar ap2_rfv("ap2_rrv","|a+|^{2}","0.5*(@0+@1)",RooArgList(r_0,r_1));
  RooFormulaVar am2_rfv("am2_rrv","|a-|^{2}","0.5*(@0-@1)",RooArgList(r_0,r_1));
  RooFormulaVar bp2_rfv("bp2_rrv","|b+|^{2}","0.5*(1.+@0-@1-@2)",RooArgList(alpha_b,r_0,r_1));
  RooFormulaVar bm2_rfv("bm2_rrv","|b-|^{2}","0.5*(1.-@0-@1+@2)",RooArgList(alpha_b,r_0,r_1));
  
  RooLbtoJpsiL0PDF3 w3("w3","w3",costheta,costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);
  
  TCanvas c("c","c",100,100);
  
  for (int i=1;i<argc;i++) {
    TString input(argv[i]);
    TString output(input);
    output.ReplaceAll(".root","");

    TFile inputfile(input);
    TDirectory* mydirectory = (TDirectory*)inputfile.Get("BtoqllGeneratedDistributions");
    TTree* mytree = (TTree*)mydirectory->Get("Lb2JpsiL0");
    RooDataSet* data =  filldataset("data", "data", RooArgList(costheta,costheta1,costheta2), mytree, "lb_id>0.");
    
    alpha_b=0.;
    r_0=0.5;
    r_1=0.;

    RooFitResult* fr=w3.fitTo(*data,RooFit::Save(true));
    plot(c,RooArgList(costheta,costheta1,costheta2),w3,data, output ,true);
    
    fr->Print("v");
  
    std::cout << output << std::endl;
    std::cout << "ap2 = " << ap2_rfv.getVal() << " +/- " << ap2_rfv.getPropagatedError(*fr) << std::endl;
    std::cout << "am2 = " << am2_rfv.getVal() << " +/- " << am2_rfv.getPropagatedError(*fr) << std::endl;
    std::cout << "bp2 = " << bp2_rfv.getVal() << " +/- " << bp2_rfv.getPropagatedError(*fr) << std::endl;
    std::cout << "bm2 = " << bm2_rfv.getVal() << " +/- " << bm2_rfv.getPropagatedError(*fr) << std::endl;
    std::cout << std::endl;
    mygetchar;
    delete data;
    delete fr;
  }
}

