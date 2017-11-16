#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSpline.h"
#include "TFile.h"
#include "TAxis.h"
#include "TMatrixDSym.h"

#include "rootlogon.h"
#include "functions.h"

TF1 *nll_tf1, *lik_tf1;//, *syst_tf1, *liktimessyst_tf1, *nll_conv_tf1, *lik_conv_tf1;
//TF1 *nlltimessyst_tf1;
TSpline3* nllspline;


double nll_func(double *x, double *par)
{
  return nllspline->Eval(x[0])*par[0];
}

double lik_func(double *x, double *par)
{
  return exp(-nll_tf1->Eval(x[0]));//    nllspline->Eval(x[0]));
}


// double syst_func(double *x, double *par)
// {
//   return TMath::Gaus(x[0],par[0],par[1],true);
// }


// double liktimessyst_func( double *x, double *par )
// {
//   return lik_tf1->Eval(x[0])*syst_tf1->Eval( par[0] - x[0] );
//   //  return lik_tf1->Eval(par[0] - x[0])*syst_tf1->Eval( x[0] );
// }

// double lik_conv_func( double *x, double *par )
// {
//   liktimessyst_tf1->SetParameter(0, x[0]);
//   return liktimessyst_tf1->Integral(par[0], par[1]);
// }

// // double nlltimessyst_func( double *x, double *par )
// // {
// //   return nll_tf1->Eval(x[0])*syst_tf1->Eval( par[0] - x[0] );
// // }

// double nll_conv_func( double *x, double *par )
// {
//   return -log(lik_conv_tf1->Eval(x[0]));
//   //  nlltimessyst_tf1->SetParameter(0, x[0]);
//   //  return nlltimessyst_tf1->Integral(par[0], par[1]);
// }


// double int_lik_conv_func(double* x, double* par)
// {
//   return lik_conv_tf1->Integral(par[0], x[0])/par[2];
// }

double int_lik_func(double* x, double* par)
{
  return lik_tf1->Integral(par[0], x[0])/par[2];
}


int main(int argc, char** argv)
{
  lhcbstyle();

  if (argc!=2) return 1;

  const TString varname(argv[1]);
  TString nicename;

  TFile inputfile("datafit-mcrw-Lb2JpsiL0-Tuple_JpsiL0_betas-MagnetAll-sim-graph-" + varname + ".root");
  TGraph* graph = (TGraph*)inputfile.Get("Graph");
  nllspline = new TSpline3("nllspline", graph);

  double valerr(0.), min(0.), max(1.), likelihoodscale(1.);//,  systsigma(1.);
  //  double systmean(0.);

  std::vector<double> significances;
  std::vector<double> limits;

  min=graph->GetXaxis()->GetXmin();
  max=graph->GetXaxis()->GetXmax();

  TFile corrfile("/afs/cern.ch/user/j/jwicht/fitLambdab2JpsiL0/toys/dataparam0/corr.root");
  //  TVectorD meansvec(*(TVectorD*)corrfile.Get("meansvector"));
  TMatrixDSym Ccorr(*(TMatrixDSym*)corrfile.Get("covmatrix"));
  corrfile.Close();

  if (varname=="alpha_b") {
    nicename="#alpha_{b}";
    likelihoodscale=1./Ccorr(1,1);

    valerr=0.17;

    // systmean=0.0;
    // systsigma=0.07;

    significances.push_back(0.49);
    significances.push_back(0.77);
    significances.push_back(0.77-0.09);
  }

  if (varname=="P_b") {
    nicename="P_{b}";
    likelihoodscale=1./Ccorr(0,0);
    //    likelihoodscale=1.;
    
    valerr=0.07;

    // systmean=0.0;
    // systsigma=0.02;

    significances.push_back(-0.20);
    significances.push_back(+0.20);
  }


  if (varname=="ap2") {
    nicename="|a+|^{2}";
    likelihoodscale=1./2.4;
    
    // systmean=0.;
    // systsigma=0.01;

    limits.push_back(0.9);
    limits.push_back(0.95);
  }

  if (varname=="bm2") {
    nicename="|b-|^{2}";
    likelihoodscale=1./2.4;
    
    // systmean=0.;
    // systsigma=0.01;

    limits.push_back(0.9);
    limits.push_back(0.95);
  }


  nll_tf1 = new TF1("nll_tf1", nll_func, min, max, 1);
  //  nll_tf1->SetNpx(1000);
  nll_tf1->SetParameter(0, likelihoodscale);

  lik_tf1 = new TF1("lik_tf1", lik_func, min, max, 0);
  //  lik_tf1->SetNpx(1000);
  //  lik_tf1->SetParameter(0, likelihoodscale);
  //  lik_tf1->SetParameter(0, 1.);


  //  systmean=nll_tf1->GetMinimumX(min,max);
  //  std::cout << systmean << std::endl;
  //  const double systmin=systmean-7.*systsigma;
  //  const double systmax=systmean+7.*systsigma;

  // syst_tf1 = new TF1("syst_tf1", syst_func, min, max, 2);
  // syst_tf1->SetNpx(10000);
  // syst_tf1->SetParameter(0, systmean);
  // syst_tf1->SetParameter(1, systsigma);

  // liktimessyst_tf1 = new TF1("liktimessyst_tf1", liktimessyst_func, min, max, 1);
  // liktimessyst_tf1->SetNpx(1000);
  // liktimessyst_tf1->SetParameter(0, 0.);

  // nlltimessyst_tf1 = new TF1("nlltimessyst_tf1", nlltimessyst_func, min, max, 1);
  // nlltimessyst_tf1->SetNpx(1000);
  // nlltimessyst_tf1->SetParameter(0, 0.);

  // lik_conv_tf1 = new TF1("lik_conv_tf1", lik_conv_func, min, max, 2 );
  // //  lik_conv_tf1->SetParameter(0,-1.);
  // //  lik_conv_tf1->SetParameter(1, 1.);
  // lik_conv_tf1->SetParameter(0,min);
  // lik_conv_tf1->SetParameter(1,max);
  // lik_conv_tf1->SetNpx(1000);

  // nll_conv_tf1 = new TF1("nll_conv_tf1", nll_conv_func, min, max, 0 );
  // nll_conv_tf1->SetNpx(1000);

  TF1* int_lik_tf1 = new TF1("int_lik_tf1", int_lik_func, min, max, 3);
  int_lik_tf1->SetParameter(0,0.);
  int_lik_tf1->SetParameter(1,max);
  int_lik_tf1->SetParameter(2,lik_tf1->Integral(0.,max));  

  // TF1* int_lik_conv_tf1 = new TF1("int_lik_conv_tf1", int_lik_conv_func, 0., max, 3);
  // int_lik_conv_tf1->SetParameter(0,0.);
  // int_lik_conv_tf1->SetParameter(1,max);
  // int_lik_conv_tf1->SetParameter(2,lik_conv_tf1->Integral(0.,max));  



  // std::cout << nll_tf1->GetMinimumX(0.5, -0.2, 0.4) << std::endl;
  // std::cout << nll_tf1->GetX(0.5, -0.2, 0.2) << std::endl;
  // std::cout << nll_tf1->GetX(0.5, 0.2, 0.4) << std::endl;
  //  return 1;

  TCanvas *c1 = new TCanvas("c1", "c1", 100, 100);
  //  c1->Divide(2, 3);
 
  //  c1->cd(1);
  nll_tf1->GetXaxis()->SetTitle(nicename);
  nll_tf1->GetYaxis()->SetTitle("NLL");
  nll_tf1->Draw();
  c1->SaveAs("convandintegrate-nll-"+ varname+".eps");

  // c1->cd(2);
  lik_tf1->GetXaxis()->SetTitle(nicename);
  lik_tf1->GetYaxis()->SetTitle("L");
  lik_tf1->Draw();
  c1->SaveAs("convandintegrate-lik-"+ varname+".eps");

  // // c1->cd(3);
  // nll_conv_tf1->GetXaxis()->SetTitle(nicename);
  // nll_conv_tf1->GetYaxis()->SetTitle("NLL-conv");
  // nll_conv_tf1->Draw();
  // c1->SaveAs("convandintegrate-nll_conv-"+ varname+".eps");

  // //  c1->cd(4);
  // lik_conv_tf1->GetXaxis()->SetTitle(nicename);
  // lik_conv_tf1->GetYaxis()->SetTitle("L-conv");
  // lik_conv_tf1->Draw();
  // c1->SaveAs("convandintegrate-lik_conv-"+ varname+".eps");

  //  c1->cd(5);
  if (limits.size()>0) {
    // int_lik_conv_tf1->GetXaxis()->SetTitle(nicename);
    // int_lik_conv_tf1->GetYaxis()->SetTitle("int L-conv");
    // int_lik_conv_tf1->Draw();
    // c1->SaveAs("convandintegrate-int_lik-"+ varname+".eps");
  }

  // //  c1->cd(6);
  // syst_tf1->GetXaxis()->SetTitle(nicename);
  // syst_tf1->GetYaxis()->SetTitle("syst");
  // syst_tf1->Draw();
  // c1->SaveAs("convandintegrate-syst-"+ varname+".eps");

  
  if (valerr>0.) {
    const Double_t valmin=nll_tf1->GetMinimumX();
    std::cout << "wout syst: NLL(" << valmin-valerr  << ") = " << nll_tf1->Eval( valmin-valerr ) << std::endl;
    std::cout << "wout syst: NLL(" << valmin+valerr  << ") = " << nll_tf1->Eval( valmin+valerr ) << std::endl;
    //    std::cout << "Significance at valmin+valerr =" << valmin+valerr  << " wout syst : " << sqrt(2.* (nll_tf1->Eval( valmin+valerr ) - nll_tf1->GetMinimum())) << std::endl;  
  }

  for (unsigned int i=0; i<significances.size(); i++) {

    std::cout << "Significance at " << significances[i] << " wout syst : " << sqrt(2.* (nll_tf1->Eval( significances[i] ) - nll_tf1->GetMinimum())) << std::endl;
    std::cout << "Significance at " << significances[i] << " with syst : " << sqrt(2.* (nll_conv_tf1->Eval( significances[i] ) - nll_conv_tf1->GetMinimum())) << std::endl;
  }

  for (unsigned int i=0; i<limits.size(); i++) {
    std::cout << "wout syst: " << varname << " < " << int_lik_tf1->GetX(limits[i],min,max) << " at " << istr(limits[i]*100.) << "%" << std::endl;
    //    std::cout << "with syst: " << varname << " < " << int_lik_conv_tf1->GetX(limits[i],min,max) << " at " << istr(limits[i]*100.) << "%" << std::endl;
  }
}
