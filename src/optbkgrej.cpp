//g++ -o optbkgrej optbkgrej.cxx  -O2 `root-config --cflags --libs` 
 
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

#include <iostream>

#include "functions.h"


int main(int argc, char** argv)
{

  const TString dirname="JpsiL0Tuple_betas";

  TFile sigfile("DVTuples-Lb2JpsiL0-MC11a.truemc.root");
  TDirectory* sigdirectory = (TDirectory*)sigfile.Get(dirname);
  TTree* sigtree = (TTree*)sigdirectory->Get("DecayTree");

  sigtree->AddFriend(dirname, "DVTuples-Lb2JpsiL0-MC11a.truemc.bdt.root");

  TFile bkgfile("DVTuples_data_stripping17.root");
  TDirectory* bkgdirectory = (TDirectory*)bkgfile.Get(dirname);
  TTree* bkgtree = (TTree*)bkgdirectory->Get("DecayTree");
  
  bkgtree->AddFriend(dirname, "DVTuples_data_stripping17.bdt.root");


  const TString Lambdabname("lambda_b0");
  const TString Lambda0name("lambda0");
  const TString Jpsiname("jpsi1s");

  const TCut Lambda0_LL("piminus_TRACK_Type==3");
  const TCut Lambda0_DD("piminus_TRACK_Type==5");

  const TCut cutall(Lambdabname + "_DIRA_OWNPV>0.7&&" + Lambdabname + "_TAU>0.&&" + Lambdabname + "_IPCHI2_OWNPV<30.&&" + Lambda0name + "_TAU>0.");
  //  const TCut cutsignal("1==1");
  //  const TCut cutsignal("abs(1.*lambda_b0_TRUEID)==5122");
  const TCut cutbkg(TString(Lambdabname + "_M>5800.&&" + Lambdabname + "_M<6000."));

  TCanvas c("c","c",100,100);
  bkgtree->Draw(Lambdabname + "_M" , cutall + cutbkg);
  c.SaveAs("c.eps");

  const double luminosity = 1.1;//fb-1
  const double bbcs       = 0.45445E12 ; //fb

  const double decprodcuteff = 0.16;
  const int nSigTot    = 532997 + 528298;
  const double bf = 4.7E-05 * 0.0593 * 0.639;
  const double f_Lb = 0.20;

  const double N_Lambdab = luminosity * bbcs * f_Lb * 2.;

  std::cout << "N(Lambdab) = " << N_Lambdab << std::endl; 
  std::cout << "N(Lambdab to Jpsi L0) = " << N_Lambdab * bf << std::endl;
  std::cout << "N(Lambdab to Jpsi L0) in acceptance = " << N_Lambdab * bf * decprodcuteff << std::endl;

  std::cout << "N(Lambdab to Jpsi L0) reconstructable = " << N_Lambdab * bf * decprodcuteff * sigtree->GetEntries(cutall) / nSigTot << std::endl;

  
  const double bdtmin=-0.5;
  const double bdtmax=0.5;
  const int nstep=20;

  const double bkgnorm=1.25*40./200.;

  for (int i=0; i<nstep; i++) {
    const double bdtcutval = bdtmin + (bdtmax-bdtmin) * i/nstep;
    const TCut bdtcut("bdt>" + str(bdtcutval));
    const double nsig = N_Lambdab * bf * decprodcuteff * sigtree->GetEntries(cutall + bdtcut + Lambda0_LL) / nSigTot;
    const double nbkg = bkgtree->GetEntries(cutall + cutbkg + bdtcut + Lambda0_LL) * bkgnorm;

    const double signif = nsig / sqrt(nsig+nbkg);

    std::cout << "bdt > " << bdtcutval << std::endl;
    std::cout << "nsig = " << nsig << std::endl;
    std::cout << "nbkg = " << nbkg << std::endl;
    std::cout << "signif = " << signif << std::endl;
  }


  for (int i=0; i<10; i++) std::cout << std::endl;

  for (int i=0; i<nstep; i++) {
    const double bdtcutval = bdtmin + (bdtmax-bdtmin) * i/nstep;
    const TCut bdtcut("bdt>" + str(bdtcutval));
    const double nsig = N_Lambdab * bf * decprodcuteff * sigtree->GetEntries(cutall + bdtcut + Lambda0_DD) / nSigTot;
    const double nbkg = bkgtree->GetEntries(cutall + cutbkg + bdtcut + Lambda0_DD) * bkgnorm;

    const double signif = nsig / sqrt(nsig+nbkg);

    std::cout << "bdt > " << bdtcutval << std::endl;
    std::cout << "nsig = " << nsig << std::endl;
    std::cout << "nbkg = " << nbkg << std::endl;
    std::cout << "signif = " << signif << std::endl;
  }

}
