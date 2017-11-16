#include "stdlib.h"
#include <iostream>
#include <fstream>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

void usage(char* argv0)
{
    std::cout << argv0 << " filename nVar neventsperfile" << std::endl;
    exit(1);
}

int main(int argc, char** argv)
{  
  
  if (argc!=4) usage(argv[0]);
  TString filename(argv[1]);
  const int nvar=atoi(argv[2]);
  if (nvar!=3 and nvar!=5) usage(argv[0]);
  const int neventsperfile=atoi(argv[3]);

  TFile myfile(filename);
  TDirectory* mydirectory = (TDirectory*)myfile.Get("BtoqllGeneratedDistributions");
  std::cout << mydirectory << std::endl;
  TTree* mytree = (TTree*)mydirectory->Get("Lb2JpsiL0");

  TString friendfile(filename);
  friendfile.ReplaceAll(".root",".acceptance.root");

  mytree->AddFriend("Lb2JpsiL0",friendfile);

  Int_t lb_id,keep;
  Float_t costheta, costheta1, costheta2;
  Float_t phi1, phi2;
  Float_t w;
  mytree->SetBranchAddress("lb_id", &lb_id);
  mytree->SetBranchAddress("keep",  &keep);
  mytree->SetBranchAddress("costheta",  &costheta);
  mytree->SetBranchAddress("costheta1", &costheta1);
  mytree->SetBranchAddress("costheta2", &costheta2);
  if (nvar==5) {
    mytree->SetBranchAddress("phi1", &phi1);
    mytree->SetBranchAddress("phi2", &phi2);
  }
  mytree->SetBranchAddress("w",   &w);
  //  mytree->SetBranchAddress("w2mu", &w2mu);

  std::fstream filestr;

  int ievent=0;
  int ifile=0;

  for (int i=0; i< mytree->GetEntries(); i++) {
    mytree->GetEvent(i);
    
    if (lb_id>0 && keep==1) {
      //if (lb_id>0) {
   
      if (ievent%neventsperfile == 0) {
	char buffer[50];
	sprintf( buffer, "%04d", ifile );
	TString output = "toys/toy_" +TString(buffer) +".dat";
	std::cout << "opening " << output;
	filestr.open(output, std::fstream::out);
      }
      if (nvar==5) {
	filestr << costheta << " " << costheta1 << " " << costheta2 << " " << phi1 << " " << phi2 << " " << w << std::endl;
      } else {
	filestr << costheta << " " << costheta1 << " " << costheta2 << " " << w << std::endl;
      }
      //      std::cout << costheta << " " << costheta1 << " " << costheta2 << " " << phi1 << " " << phi2 << std::endl;
      
      if (ievent%neventsperfile == neventsperfile-1) {
	std::cout << " - closing" << std::endl;
	filestr.close();
	ifile++;
      }
      ievent++;
    }
    //    std::cout << __LINE__ << std::endl;
  }


}
