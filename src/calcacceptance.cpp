//g++ -o calcacceptance calcacceptance.cxx -O2 `root-config --cflags --libs` 

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TRandom3.h"

#include "calcacceptanceclass.h"


int main(int argc, char** argv)
{
  if (argc!=2) {
    std::cerr << argv[0] << " filename" << std::endl;
    return 1;
  }


  TRandom3* rnd = new TRandom3();
  calcacceptanceclass mycalcacceptanceclass("coeff3.txt");

  //acceptance histogram
  TFile accfile("acc_histo.root");
  TH3F* acc_histo = (TH3F*)accfile.Get("acc_histo");
  const Double_t acc_histo_min=acc_histo->GetMinimum();
  const Double_t acc_histo_max=acc_histo->GetMaximum();

  std::cout << "acc_histo_min = " << acc_histo_min << std::endl;
  std::cout << "acc_histo_max = " << acc_histo_max << std::endl;
  

  const TString filename(argv[1]);
  TFile myfile(filename);

    TDirectory* mydirectory = (TDirectory*)myfile.Get("BtoqllGeneratedDistributions");
  //  std::cout << mydirectory << std::endl;
  TTree* mytree = (TTree*)mydirectory->Get("Lb2JpsiL0");

  Float_t costheta,costheta1,costheta2;
  Float_t phi1,phi2;

  mytree->SetBranchAddress("costheta",  &costheta);
  mytree->SetBranchAddress("costheta1", &costheta1);
  mytree->SetBranchAddress("costheta2", &costheta2);
  mytree->SetBranchAddress("phi1",      &phi1);
  mytree->SetBranchAddress("phi2",      &phi2);

  Int_t keep;
  Float_t w;

  TString output(filename);
  output.ReplaceAll(".root",".acceptance.root");

  TFile *f = new TFile(output,"recreate");
  TTree *T = new TTree("Lb2JpsiL0","test friend trees");
  T->Branch("keep", &keep, "keep/I"); 
  T->Branch("w",    &w,    "w/F"); 
  //  T->Branch("costheta", &costheta, "costheta/F");
  //  T->Branch("costheta1",&costheta1,"costheta1/F");
  //  T->Branch("costheta2",&costheta2,"costheta2/F");
  //  T->Branch("phi1norm",&phi1norm,"phi1norm/F");
  //  T->Branch("phi2norm",&phi2norm,"phi2norm/F");
  
  
  for (Long64_t i=0; i<mytree->GetEntries(); i++) {
  //  for (Int_t i=0; i<10; i++) {
    mytree->GetEvent(i);
    

    const Int_t xbin=1+int((costheta+1.)/2.*10.);
    const Int_t ybin=1+int((costheta1+1.)/2.*10.);
    const Int_t zbin=1+int((costheta2+1.)/2.*10.);

    const Double_t normvalue_hist = (acc_histo->GetBinContent(xbin,ybin,zbin)-acc_histo_min)/(acc_histo_max-acc_histo_min);

    const Double_t normvalue = mycalcacceptanceclass.evaluate(costheta,costheta1,costheta2);

    const Double_t rndvalue = rnd->Rndm();
    
    //    std::cout << "normvalue   = " << normvalue << std::endl;
    //    std::cout << "normvaluemu = " << normvaluemu << std::endl;

    keep = (normvalue > rndvalue);
    w= 1./normvalue;
    //    w2mu=1./normvalue;

    // phi1norm=phi1/3.14159;
    // if (phi1norm>1.) phi1norm=1.;
    // if (phi1norm<-1.) phi1norm=-1.;
    // phi2norm=phi2/3.14159;
    // if (phi2norm>1.) phi2norm=1.;
    // if (phi2norm<-1.) phi2norm=-1.;

    T->Fill();
  }

  T->Write();
  f->Close();

}
