#include "TRandom3.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooLbtoJpsiL0PDF3.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"

#include "rootlogon.h"
#include "functions-roofit.h"

int main()
{
  rootlogon();

  RooRealVar costheta("costheta",   "cos(#theta)",     -1., 1.);
  RooRealVar costheta1("costheta1", "cos(#theta_{1})", -1., 1.);
  RooRealVar costheta2("costheta2", "cos(#theta_{2})", -1., 1.);

  const double pi = 4.* atan(1.);
  
  RooRealVar theta("theta",   "#theta",     0., pi);
  RooRealVar theta1("theta1", "#theta_{1}", 0., pi);
  RooRealVar theta2("theta2", "#theta_{2}", 0., pi);

  RooFormulaVar costhetaf("costhetaf","cos(@0)",theta);
  RooFormulaVar costheta1f("costheta1f","cos(@0)",theta1);
  RooFormulaVar costheta2f("costheta2f","cos(@0)",theta2);

  RooFormulaVar thetaf( "thetaf", "acos(costheta)", costheta);
  RooFormulaVar theta1f("theta1f","acos(costheta1)",costheta1);
  RooFormulaVar theta2f("theta2f","acos(costheta2)",costheta2);


  // Define resolution R(psi)
  RooRealVar gbias("gbias","gbias",0.,0.,1) ;
  RooRealVar greso("greso","greso",0.005,0.1,1.0) ;
  RooGaussian Rtheta("Rtheta","Rtheta",thetaf,gbias,greso) ;



  RooRealVar P_b("P_b","P_{b}",0.40, -1.5, 1.5);//0.50
  RooRealVar alpha_lambda("alpha_lambda","#alpha_{#lambda}",0.642, 0., 1.);
  alpha_lambda.setConstant();
  //  RooGaussian alpha_lambda_constrain("alpha_lambda_constrain","alpha_lambda_constrain",alpha_lambda,RooConst(0.642),RooConst(0.013));
  
  
  RooRealVar alpha_b("alpha_b","#alpha_{b}",-0.45767,-1.1,1.1);//0.457
  RooRealVar r_0("r_0","r_{0}",0.251733,-0.3,2.3);//0.5
  RooRealVar r_1("r_1","r_{1}",0.116484,-1.3,1.3);//0.1

  //  RooLbtoJpsiL0PDF3 w3("w3","w3",costheta,costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);
  RooLbtoJpsiL0PDF3 w3("w3","w3",costheta,costheta1,costheta2,P_b,alpha_b,alpha_lambda,r_0,r_1);

  RooDataSet* data = w3.generate(RooArgSet(costheta,costheta1,costheta2),10000);

  RooDataSet* data2 = new RooDataSet("data2","data2",RooArgList(costheta,costheta1,costheta2));
    
  TRandom3 rnd;
  for (int i=0;i<data->sumEntries();i++) {
      const RooArgSet* blu = data->get(i);
      //    blu->Print("v");
      //	dcostheta=costheta.getVal();
      //	dcostheta1=costheta1.getVal();
      //	dcostheta2=costheta2.getVal();
      costheta.setVal( cos(acos(blu->getRealValue(costheta.GetName()))   + rnd.Gaus(gbias.getVal(),greso.getVal())));
      costheta1.setVal(cos(acos(blu->getRealValue(costheta1.GetName()))  + rnd.Gaus(gbias.getVal(),greso.getVal())));
      costheta2.setVal(cos(acos(blu->getRealValue(costheta2.GetName()))  + rnd.Gaus(gbias.getVal(),greso.getVal())));

      data2->add(RooArgSet(costheta,costheta1,costheta2));
	      //      w2=ablu->getRealValue(w2str);
  }
  
  RooFFTConvPdf w3conv("Mf","Mf",thetaf,costheta,w3,Rtheta) ;
  w3conv.setBufferFraction(0.);

  TCanvas c("c","c",100,100);

  RooPlot* costhetaplot = costheta.frame(25);
  //  w3.plotOn(costhetaplot,RooFit::LineColor(kRed));
  //  w3conv.plotOn(costhetaplot,RooFit::LineColor(kBlue));

  data->plotOn(costhetaplot,RooFit::LineColor(kRed));
  data2->plotOn(costhetaplot,RooFit::LineColor(kBlue));
  costhetaplot->Draw();

  c.SaveAs("w3conv.eps");


}
