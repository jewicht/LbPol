#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TRandom3.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TSpline.h"

#include "functions.h"


int main()
{
  std::vector<std::string> directorylist;
  directorylist.push_back("JpsiKSTuple_B2XMuMu");


  std::vector<TString> varstocorr;
  varstocorr.push_back("b0_PT");


  const TString inputfile="DVTuples-Bd2JpsiKS-MC11a.root";
  TFile myfile(inputfile);

  
  TCanvas c("c","c",100,100);

  TString outputname=inputfile;
  outputname.ReplaceAll(".root",".reweighted.root");
  std::cout << "outputname = " << outputname << std::endl;
  TFile* rwfile = new TFile(outputname);
  std::cerr << __LINE__ << std::endl;

  for (unsigned int idir=0; idir<directorylist.size(); idir++) {
    const TString dirname = directorylist[idir];
    TDirectory* mydir = (TDirectory*)myfile.Get(dirname);
    std::cerr << __LINE__ << std::endl;
    TTree* mytree = (TTree*)mydir->Get("DecayTree");
    std::cerr << __LINE__ << std::endl;
    TTree* rwtree = (TTree*)rwfile->Get(dirname);
    std::cerr << __LINE__ << std::endl;    

    for (unsigned int j=0; j<varstocorr.size(); j++) {

      rwtree->Draw(varstocorr[j]);
      std::cerr << __LINE__ << std::endl;
      //      mytree->Draw(varstocorr[j],"","SAME");
      //      std::cerr << __LINE__ << std::endl;
      c.SaveAs("rw_" + varstocorr[j] + ".eps");
    }
  }
}
1
