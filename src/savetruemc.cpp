//g++ -o savetruemc savetruemc.cxx -O2 `root-config --cflags --libs` 

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TCut.h"


int main(int argc, char** argv)
{
  if (argc!=2) {
    std::cerr << argv[0] << " file.root"  << std::endl;
    return 1;
  }

  const TString filename(argv[1]);
  TFile myfile(filename);
  
  TString output(filename);
  output.ReplaceAll(".root", ".truemc.root");
  
  TFile *f = new TFile(output,"recreate");
  
  std::vector<std::string> directorylist;
  
  directorylist.push_back("JpsiL0Tuple_betas");
  directorylist.push_back("JpsiL0Tuple_B2XMuMu");
  directorylist.push_back("JpsiKSTuple_prescaled_betas");
  directorylist.push_back("JpsiKSTuple_detached_betas");
  directorylist.push_back("JpsiKSTuple_B2XMuMu");
  
  for (unsigned int i=0; i<directorylist.size(); i++) {
    const TString& directoryname = directorylist[i];

    const bool doBd2JpsiKS0=directoryname.Contains("KSTuple");

    TDirectory* mydirectory = (TDirectory*)myfile.Get(directoryname);
    if (!mydirectory) continue;

    f->mkdir(directoryname);
    f->cd(directoryname);

    TTree* mytree = (TTree*)mydirectory->Get("DecayTree");

    TCut mycut="abs(1.*lambda_b0_TRUEID)==5122";//&&(L0DiMuonDecision==1)&&Hlt1DiMuonHighMassDecision==1&&Hlt2DiMuonJPsiDecision==1.&&jpsi1sL0Global_TOS==1&&jpsi1sHlt1Global_TOS==1&&jpsi1sHlt2Global_TOS==1";
    if (doBd2JpsiKS0) mycut="abs(1.*b0_TRUEID)==511";//&&(L0DiMuonDecision==1)&&Hlt1DiMuonHighMassDecision==1&&Hlt2DiMuonJPsiDecision==1.&&jpsi1sL0Global_TOS==1&&jpsi1sHlt1Global_TOS==1&&jpsi1sHlt2Global_TOS==1";

    TTree *T = mytree->CopyTree(mycut);
    T->Write();
  }
  f->Close();
  
}
