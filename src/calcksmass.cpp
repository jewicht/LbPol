//g++ -o calcksmass calcksmass.cxx -O2 `root-config --cflags --libs` 

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"





int main(int argc, char** argv)
{

  if (argc!=2) {
    std::cerr << argv[0] << " file" << std::endl;
    return 1;
  }
  const double pi=atan(1.)*4.;

  const TString filename(argv[1]);

  TFile myfile(filename);

  Float_t  
    muplus_PE,
    muplus_PX,
    muplus_PY,
    muplus_PZ,
    
    muminus_PE,
    muminus_PX,
    muminus_PY,
    muminus_PZ,
    
    piminus_PE,
    piminus_PX,
    piminus_PY,
    piminus_PZ,
    
    pplus_PE,
    pplus_PX,
    pplus_PY,
    pplus_PZ;

  Int_t
    muplus_ID,
    muminus_ID,
    piminus_ID,
    pplus_ID;

  TString output(filename);
  output.ReplaceAll(".root", ".ksmass.root");
  TFile *f = new TFile(output,"recreate");


  std::vector<TString> directorylist;
  directorylist.push_back("Tuple_JpsiL0_betas");
  directorylist.push_back("Tuple_JpsiL0_betas_refitted");
  directorylist.push_back("Tuple_JpsiL0_B2XMuMu");
  directorylist.push_back("Tuple_JpsiL0_B2XMuMu_refitted");
  directorylist.push_back("Tuple_JpsiKs_detached_betas");
  directorylist.push_back("Tuple_JpsiKs_prescaled_betas");
  directorylist.push_back("Tuple_JpsiKs_detached_betas_refitted");
  directorylist.push_back("Tuple_JpsiKs_B2XMuMu");
  directorylist.push_back("Tuple_JpsiKs_B2XMuMu_refitted");

  for (unsigned int idl=0; idl<directorylist.size(); idl++) {

    TDirectory* mydirectory = (TDirectory*)myfile.Get(directorylist[idl]);
    if (!mydirectory) continue;

    TTree* mytree = (TTree*)mydirectory->Get("DecayTree");
    
    mytree->SetBranchAddress("muplus_ID", &muplus_ID);
    mytree->SetBranchAddress("muplus_PX", &muplus_PX);
    mytree->SetBranchAddress("muplus_PY", &muplus_PY);
    mytree->SetBranchAddress("muplus_PZ", &muplus_PZ);
    mytree->SetBranchAddress("muplus_PE", &muplus_PE);
    
    mytree->SetBranchAddress("muminus_ID", &muminus_ID);
    mytree->SetBranchAddress("muminus_PX", &muminus_PX);
    mytree->SetBranchAddress("muminus_PY", &muminus_PY);
    mytree->SetBranchAddress("muminus_PZ", &muminus_PZ);
    mytree->SetBranchAddress("muminus_PE", &muminus_PE);
    
    if (directorylist[idl].Contains("JpsiKs")) {
      mytree->SetBranchAddress("piplus_ID", &pplus_ID);
      mytree->SetBranchAddress("piplus_PX", &pplus_PX);
      mytree->SetBranchAddress("piplus_PY", &pplus_PY);
      mytree->SetBranchAddress("piplus_PZ", &pplus_PZ);
      mytree->SetBranchAddress("piplus_PE", &pplus_PE);
    } else {
      mytree->SetBranchAddress("pplus_ID", &pplus_ID);
      mytree->SetBranchAddress("pplus_PX", &pplus_PX);
      mytree->SetBranchAddress("pplus_PY", &pplus_PY);
      mytree->SetBranchAddress("pplus_PZ", &pplus_PZ);
      mytree->SetBranchAddress("pplus_PE", &pplus_PE);
    }


    mytree->SetBranchAddress("piminus_ID", &piminus_ID);
    mytree->SetBranchAddress("piminus_PX", &piminus_PX);
    mytree->SetBranchAddress("piminus_PY", &piminus_PY);
    mytree->SetBranchAddress("piminus_PZ", &piminus_PZ);
    mytree->SetBranchAddress("piminus_PE", &piminus_PE);
    
    Float_t Lambda0_M_pmisid;
    Float_t Lambdab_M_pmisid;
    
    TTree *T = new TTree(directorylist[idl],"test friend trees");
    T->Branch("Lambda0_M_pmisid",&Lambda0_M_pmisid,"Lambda0_M_pmisid/F");
    T->Branch("Lambdab_M_pmisid",&Lambdab_M_pmisid,"Lambdab_M_pmisid/F");
    
    const double pionmass= 139.57018;
    const double lambda0mass = 1115.683;
    
    for (Long64_t i=0; i<mytree->GetEntries(); i++) {
      //  for (Int_t i=0; i<10; i++) {
      mytree->GetEvent(i);
      
      const TLorentzVector lv_muplus(muplus_PX, muplus_PY, muplus_PZ, muplus_PE);
      const TLorentzVector lv_muminus(muminus_PX, muminus_PY, muminus_PZ, muminus_PE);
      
      //E2=P2+M2
      
      const TVector3 vec3_pionmisid(pplus_PX, pplus_PY, pplus_PZ);
      const TLorentzVector lv_pionmisid(vec3_pionmisid, sqrt(pionmass*pionmass+vec3_pionmisid.Mag2()));
      const TLorentzVector lv_pion(piminus_PX, piminus_PY, piminus_PZ, piminus_PE);
      
      // const TVector3 vec3_test(pplus_PX, pplus_PY, pplus_PZ);
      // const TLorentzVector lv_test(vec3_test, sqrt(lambda0mass*lambda0mass+vec3_test.Mag2()));
      
      // std::cout << pplus_PE << " vs " << lv_test.E() << std::endl;
      
      Lambda0_M_pmisid=(lv_pion+lv_pionmisid).M();
      Lambdab_M_pmisid=(lv_pion+lv_pionmisid+lv_muplus+lv_muminus).M();
      
      T->Fill();
    }

    T->Write();
  }
  f->Close();
}
