//g++ -o rewritephi rewritephi.cxx -O2 `root-config --cflags --libs` 

#include <iostream>

#include "TFile.h"
#include "TTree.h"



int main(int argc, char** argv)
{
  const TString filename(argv[1]);

  TFile myfile(filename);

    TDirectory* mydirectory = (TDirectory*)myfile.Get("NewTuple");
  //  std::cout << mydirectory << std::endl;
  TTree* mytree = (TTree*)mydirectory->Get("DecayTree");

  Float_t  
    costheta,costheta1,costheta2,phi1,phi2;

  mytree->SetBranchAddress("LambdabAngles_costheta",  &costheta);
  mytree->SetBranchAddress("LambdabAngles_costheta1", &costheta1);
  mytree->SetBranchAddress("LambdabAngles_costheta2", &costheta2);
  mytree->SetBranchAddress("LambdabAngles_phi1", &phi1);
  mytree->SetBranchAddress("LambdabAngles_phi2", &phi2);

  Float_t phi1norm, phi2norm;

  TString output(filename);
  output.ReplaceAll(".root","");
  output+=".phinorm.root";

  TFile *f = new TFile(output,"recreate");
  TTree *T = new TTree("DecayTree","test friend trees");
  T->Branch("costheta", &costheta, "costheta/F");
  T->Branch("costheta1",&costheta1,"costheta1/F");
  T->Branch("costheta2",&costheta2,"costheta2/F");
  T->Branch("phi1norm",&phi1norm,"phi1norm/F");
  T->Branch("phi2norm",&phi2norm,"phi2norm/F");
  
  
  for (Long64_t i=0; i<mytree->GetEntries(); i++) {
  //  for (Int_t i=0; i<10; i++) {
    mytree->GetEvent(i);
    
    phi1norm=phi1/3.14159;
    if (phi1norm>1.) phi1norm=1.;
    if (phi1norm<-1.) phi1norm=-1.;
    phi2norm=phi2/3.14159;
    if (phi2norm>1.) phi2norm=1.;
    if (phi2norm<-1.) phi2norm=-1.;

    T->Fill();
  }

  T->Write();
  f->Close();

}
