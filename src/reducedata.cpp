//g++ -o reducedata reducedata.C -O2 `root-config --cflags --libs` 

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TCut.h"
#include "cuts.h"


int main(int argc, char** argv)
{
  if (argc!=2)  {
    std::cout << argv[0] << " file.root" << std::endl;
    return 1;
  }

  std::vector<std::string> directorylist;
  directorylist.push_back("Tuple_JpsiKs_detached_betas");
  directorylist.push_back("Tuple_JpsiKs_prescaled_betas");
  //  directorylist.push_back("Tuple_JpsiKs_detached_betas_refitted");
  //    directorylist.push_back("Tuple_JpsiKs_B2XMuMu");
  //    directorylist.push_back("Tuple_JpsiKs_B2XMuMu_refitted");
  directorylist.push_back("Tuple_JpsiL0_betas");
  //  directorylist.push_back("Tuple_JpsiL0_betas_refitted");
  //    directorylist.push_back("Tuple_JpsiL0_B2XMuMu");
  //    directorylist.push_back("Tuple_JpsiL0_B2XMuMu_refitted");


  const TString filename(argv[1]);
  TFile myfile(filename);

  TString output(filename);
  output.ReplaceAll(".root", ".reduced.root");
  TFile *f = new TFile(output,"recreate");

  for (unsigned int idl=0; idl<directorylist.size(); idl++) {

    //lambda_b0_DIRA_OWNPV>0.7&&lambda_b0_TAU>0.&&lambda_b0_IPCHI2_OWNPV<30.&&lambda0_TAU>0.

    const TString& directoryname = directorylist[idl];

    const bool doBd2JpsiKS0=directoryname.Contains("JpsiKs");
    TString Lambdabname, Lambda0name, Jpsiname;
    varnames(doBd2JpsiKS0,Lambdabname, Lambda0name, Jpsiname);
    
    std::vector<TString> varlist;
    varlist.push_back("Polarity");
    varlist.push_back(Lambdabname + "_M");
    varlist.push_back(Lambdabname + "_P");
    varlist.push_back(Lambdabname + "_PT");
    varlist.push_back(Lambdabname + "_TAU");
    varlist.push_back(Lambdabname + "_DIRA_OWNPV");
    varlist.push_back(Lambdabname + "_IPCHI2_OWNPV");
    varlist.push_back(Lambdabname + "_FDCHI2_OWNPV");
    varlist.push_back(Lambdabname + "_ENDVERTEX_CHI2");
    varlist.push_back(Lambdabname + "_ENDVERTEX_NDOF");
    if (filename.Contains("MC")) varlist.push_back(Lambdabname + "_TRUEID");
    varlist.push_back(Lambdabname + "_ID");
    varlist.push_back(Lambdabname + "_PX");
    varlist.push_back(Lambdabname + "_PY");
    varlist.push_back(Lambdabname + "_PZ");
    varlist.push_back(Lambdabname + "_PE");

    varlist.push_back(Lambda0name + "_ENDVERTEX_CHI2");
    varlist.push_back(Lambda0name + "_ENDVERTEX_NDOF");
    varlist.push_back(Lambda0name + "_M");
    varlist.push_back(Lambda0name + "_PT");
    varlist.push_back(Lambda0name + "_P");
    varlist.push_back(Lambda0name + "_TAU");
    varlist.push_back(Lambda0name + "_ENDVERTEX_Z");
    varlist.push_back(Lambda0name + "_ID");
    varlist.push_back(Lambda0name + "_PX");
    varlist.push_back(Lambda0name + "_PY");
    varlist.push_back(Lambda0name + "_PZ");
    varlist.push_back(Lambda0name + "_PE");
    
    varlist.push_back("piminus_ID");
    varlist.push_back("piminus_TRACK_Type");
    varlist.push_back("piminus_PX");
    varlist.push_back("piminus_PY");
    varlist.push_back("piminus_PZ");
    varlist.push_back("piminus_PE");
    varlist.push_back("piminus_P");
    varlist.push_back("piminus_PT");

    if (doBd2JpsiKS0) {
      varlist.push_back("piplus_ID");
      varlist.push_back("piplus_PX");
      varlist.push_back("piplus_PY");
      varlist.push_back("piplus_PZ");
      varlist.push_back("piplus_PE");
      varlist.push_back("piplus_P");
      varlist.push_back("piplus_PT");
    } else {
      varlist.push_back("pplus_ID");
      varlist.push_back("pplus_PX");
      varlist.push_back("pplus_PY");
      varlist.push_back("pplus_PZ");
      varlist.push_back("pplus_PE");
      varlist.push_back("pplus_P");
      varlist.push_back("pplus_PT");
    }
    if (!doBd2JpsiKS0) {
      varlist.push_back("piminus_ProbNNk");
      varlist.push_back("piminus_ProbNNp");
      varlist.push_back("piminus_ProbNNpi");
      varlist.push_back("pplus_ProbNNk");
      varlist.push_back("pplus_ProbNNp");
      varlist.push_back("pplus_ProbNNpi");

      varlist.push_back("piminus_PIDK");
      varlist.push_back("piminus_PIDp");
      varlist.push_back("pplus_PIDK");
      varlist.push_back("pplus_PIDp");

    } else {

      varlist.push_back("piminus_ProbNNk");
      varlist.push_back("piminus_ProbNNp");
      varlist.push_back("piminus_ProbNNpi");
      varlist.push_back("piplus_ProbNNk");
      varlist.push_back("piplus_ProbNNp");
      varlist.push_back("piplus_ProbNNpi");

      varlist.push_back("piminus_PIDK");
      varlist.push_back("piminus_PIDp");
      varlist.push_back("piplus_PIDK");
      varlist.push_back("piplus_PIDp");

    }
    
    varlist.push_back(Jpsiname + "_ENDVERTEX_CHI2");
    varlist.push_back(Jpsiname + "_ENDVERTEX_NDOF");
    varlist.push_back(Jpsiname + "_M");
    varlist.push_back(Jpsiname + "_MM");
    varlist.push_back(Jpsiname + "L0MuonDecision_TOS");
    varlist.push_back(Jpsiname + "L0DiMuonDecision_TOS");
    varlist.push_back(Jpsiname + "Hlt1TrackMuonDecision_TOS");
    varlist.push_back(Jpsiname + "Hlt1TrackAllL0Decision_TOS");
    varlist.push_back(Jpsiname + "Hlt1DiMuonHighMassDecision_TOS");
    varlist.push_back(Jpsiname + "Hlt2DiMuonDetachedJPsiDecision_TOS");
    varlist.push_back(Jpsiname + "_ID");
    varlist.push_back(Jpsiname + "_PX");
    varlist.push_back(Jpsiname + "_PY");
    varlist.push_back(Jpsiname + "_PZ");
    varlist.push_back(Jpsiname + "_PE");

    varlist.push_back("muplus_ID");
    varlist.push_back("muplus_PX");
    varlist.push_back("muplus_PY");
    varlist.push_back("muplus_PZ");
    varlist.push_back("muplus_PE");
    varlist.push_back("muplus_P");
    varlist.push_back("muplus_PT");

    varlist.push_back("muminus_ID");
    varlist.push_back("muminus_PX");
    varlist.push_back("muminus_PY");
    varlist.push_back("muminus_PZ");
    varlist.push_back("muminus_PE");
    varlist.push_back("muminus_P");
    varlist.push_back("muminus_PT");

    TDirectory* mydirectory = (TDirectory*)myfile.Get(directoryname);
    TTree* mytree = (TTree*)mydirectory->Get("DecayTree");

    mytree->SetBranchStatus("*",0);  // disable all branches
    for (unsigned int i=0;i<varlist.size();i++) mytree->SetBranchStatus(varlist[i],1);  // activate branchname

    f->mkdir(directoryname);
    f->cd(directoryname);
    
    //    TString mycut="lambda_b0_M>5500.&&lambda_b0_M<5750.&&lambda_b0_DIRA_OWNPV>0.7&&lambda_b0_TAU>0.&&lambda_b0_IPCHI2_OWNPV<30.&&lambda0_TAU>0.";
    //    if (directoryname.Contains("JpsiKSTuple")) mycut="b0_M>5000.&&b0_M<5500.&&b0_DIRA_OWNPV>0.7&&b0_TAU>0.&&b0_IPCHI2_OWNPV<30.&&ks0_TAU>0.";

    //    TString mycut="lambda_b0_M>5500.&&lambda_b0_M<5750.&&lambda_b0_DIRA_OWNPV>0.7";
    //    if (directoryname.Contains("JpsiKSTuple")) mycut="b0_M>5000.&&b0_M<5500.&&b0_DIRA_OWNPV>0.7";

    TCut mycut = triggercut(false);

    std::cout << "Reducing " << filename << ":" << directoryname <<std::endl;
    TTree *T = mytree->CopyTree(mycut);
    //    T->Write();
    T->AutoSave();

  }
  f->Close();
  myfile.Close();
}



// int main(int argc, char** argv)
// {
//   if (argc==2) {
//     const TString filename(argv[1]);
//     TFile myfile(filename);
//     TTree* mytree = (TTree*)myfile.Get("DecayTree");
    
//     TString output(filename);
//     output.ReplaceAll(".root",".reduced.root");
    
//     TFile *f = new TFile(output,"recreate");
//     TTree *T = mytree->CopyTree("Lambda_b_M>5400.&&Lambda_b_M<5800.&&Lambda_b_DIRA_OWNPV>0.9999&&Lambda_b_FDCHI2_OWNPV>100.&&Lambda_b_IPCHI2_OWNPV<9.&&Lambda_b_ENDVERTEX_CHI2/Lambda_b_ENDVERTEX_NDOF<8.&&Lambda0_PT>500.&&abs(Lambda0_M-1115.)<30.&&J_psi_1S_ENDVERTEX_CHI2/J_psi_1S_ENDVERTEX_NDOF<9.&&Lambda_b_LifetimeFit_Lambda0_ctau/0.2999>2.");
    
//     T->Write();
//     f->Close();
//   } else if (argc==3) {
//     const TString filename(argv[1]);
//     const TString directoryname(argv[2]);
   
//     TFile myfile(filename);

//     TDirectory* mydirectory = (TDirectory*)myfile.Get(argv[2]);
//     TTree* mytree = (TTree*)mydirectory->Get("DecayTree");
    
//     TString output(filename);
//     output.ReplaceAll(".root", "." + directoryname + ".reduced.root");
    
//     TFile *f = new TFile(output,"recreate");
//     f->mkdir(directoryname);
//     f->cd(directoryname);
//     //    TTree *T = mytree->CopyTree("lambda_b0_M>5400.&&lambda_b0_M<5800.&&lambda_b0_DIRA_OWNPV>0.9999&&lambda_b0_FDCHI2_OWNPV>100.&&lambda_b0_IPCHI2_OWNPV<9.&&lambda_b0_ENDVERTEX_CHI2/lambda_b0_ENDVERTEX_NDOF<8.&&lambda0_PT>500.&&abs(lambda0_M-1115.)<30.&&jpsi1s_ENDVERTEX_CHI2/jpsi1s_ENDVERTEX_NDOF<9.&&lambda_b0_TAU>0.002");
//     TTree *T = mytree->CopyTree("lambda_b0_M>5400.&&lambda_b0_M<5800.&&lambda_b0_DIRA_OWNPV>0.999");
    
//     T->Write();
//     f->Close();
    
    
//   }
// }
