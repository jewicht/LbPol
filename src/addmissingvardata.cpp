#include <iostream>

#include "TFile.h"
#include "TTree.h"

int main(int argc, char** argv)
{

  if (argc!=2) {
    std::cerr << argv[0] << " file" << std::endl;
    return 1;
  }

  const bool debug=false;

  const TString filename(argv[1]);
  TFile myfile(filename);
  

  std::vector<TString> directorylist;
  //  directorylist.push_back("JpsiKSTuple_B2XMuMu");
  //  directorylist.push_back("JpsiKSTuple_prescaled_betas");
  //  directorylist.push_back("JpsiKSTuple_detached_betas");
  //  directorylist.push_back("JpsiL0Tuple_B2XMuMu");
  directorylist.push_back("Tuple_JpsiL0_betas");
  directorylist.push_back("Tuple_JpsiKs_detached_betas");

  TString output(filename);
  output.ReplaceAll(".root", ".missvar.root");
  TFile *f = NULL;
  if (!debug) f = new TFile(output,"recreate");

  //  Float_t L0DiMuonDecision=1.;
  Float_t lambdab0_TRUEID=0.;
  Float_t rw=0.;
  for (unsigned int idl=0; idl<directorylist.size(); idl++) {
    
    TDirectory* mydirectory = (TDirectory*)myfile.Get(directorylist[idl]);//.c_str());
    if (!mydirectory) continue;
    TTree* mytree = (TTree*)mydirectory->Get("DecayTree");

    const bool isJpsiKs=directorylist[idl].Contains("JpsiKs");

    TTree *T =NULL;
    if (!debug) {
      T = new TTree(directorylist[idl],"test friend trees");
      //      T->Branch("L0DiMuonDecision", &L0DiMuonDecision, "L0DiMuonDecision/F");

      T->Branch("rw",&rw,"rw/F");
      if (isJpsiKs) {
	T->Branch("b0_TRUEID",&lambdab0_TRUEID,"b0_TRUEID/F");
      } else {
	T->Branch("lambda_b0_TRUEID",&lambdab0_TRUEID,"lambda_b0_TRUEID/F");
      }
    }

    const Long64_t nentries=debug?4:mytree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) {
      if (!debug) T->Fill();
    }
    
    //    if (!debug) T->Write()AutoSave();
    if (!debug) T->AutoSave();
  }
  if (!debug) f->Close();
  myfile.Close();
}
