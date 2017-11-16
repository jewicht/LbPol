#include "TFile.h"
#include "TTree.h"
#include "math.h"
#include "TMVA/Reader.h"
#include "cuts.h"
#include "functions-roofit.h"

TString Lambdabname,Lambda0name,Jpsiname,pplusname;

//for stripping17
struct fd
{
  Float_t f;
  Float_t d;
};

//For stripping19a
// struct fd
// {
//   Float_t f;
//   Double_t d;
// };


fd Lb_DIRA_OWNPV, Lb_IPCHI2_OWNPV, Lb_FDCHI2_OWNPV, Lb_TAU;

fd Lb_ENDVERTEX_CHI2;
fd L0_M;
fd Jpsi_ENDVERTEX_CHI2, Jpsi_MM;
fd piminus_ProbNNk, piminus_ProbNNp, piminus_ProbNNpi;
fd pplus_ProbNNk, pplus_ProbNNp, pplus_ProbNNpi;

Int_t Lb_ENDVERTEX_NDOF;
Int_t Jpsi_ENDVERTEX_NDOF;
Int_t piminus_TRACK_Type;


Float_t Lb_FDCHI2_LOG,Lb_DIRA_LOG, Lb_ENDVERTEX_ratio, Jpsi_ENDVERTEX_ratio, Jpsi_DMM,L0_DM;
Float_t pion_DLL_Kpi,pion_DLL_ppi, proton_DLL_Kp, proton_DLL_pip;
Float_t bdt;


const bool usepid=false;


void addvariables(TMVA::Reader* readerll)
{
  readerll->AddVariable("log10(1.0000001-" + Lambdabname + "_DIRA_OWNPV)", &Lb_DIRA_LOG);
  readerll->AddVariable(Lambdabname + "_IPCHI2_OWNPV", &Lb_IPCHI2_OWNPV.f);
  readerll->AddVariable("log10(" + Lambdabname + "_FDCHI2_OWNPV)", &Lb_FDCHI2_LOG);
  readerll->AddVariable(Lambdabname + "_ENDVERTEX_CHI2/" + Lambdabname + "_ENDVERTEX_NDOF", &Lb_ENDVERTEX_ratio);
  readerll->AddVariable(Lambdabname +  "_TAU", &Lb_TAU.f);
 
  readerll->AddVariable("abs(" + Lambda0name + "_M-" + str(getLambda0mass(false)) + ")", &L0_DM); 
  //  readerll->AddVariable(Lambda0name + "_PT", &L0_PT);
  //  readerll->AddVariable(Lambda0name + "_TAU", &L0_TAU);
  
  readerll->AddVariable("abs(" + Jpsiname + "_MM-" + str(Jpsipdgmass) + ")", &Jpsi_DMM);
  readerll->AddVariable(Jpsiname + "_ENDVERTEX_CHI2/" + Jpsiname + "_ENDVERTEX_NDOF", &Jpsi_ENDVERTEX_ratio);

  if (usepid) {
    readerll->AddVariable("piminus_ProbNNk-piminus_ProbNNpi", &pion_DLL_Kpi);
    readerll->AddVariable("piminus_ProbNNp-piminus_ProbNNpi", &pion_DLL_ppi);
    
    readerll->AddVariable("pplus_ProbNNk-pplus_ProbNNp",  &proton_DLL_Kp);
    readerll->AddVariable("pplus_ProbNNpi-pplus_ProbNNp", &proton_DLL_pip);
  }
  //  readerll->AddSpectator(Lambdabname + "_M", &Lb_M);
}


int main(int argc, char** argv)
{
  if (argc!=2) {
    std::cerr << argv[0] << " filename" << std::endl;
    return 1;
  }

  const TString filename(argv[1]);
  TFile myfile(filename);

  std::vector<TString> directorylist;
  directorylist.push_back("Tuple_JpsiKs_detached_betas");
  directorylist.push_back("Tuple_JpsiKs_prescaled_betas");
  directorylist.push_back("Tuple_JpsiL0_betas");
  //    directorylist.push_back("Tuple_JpsiKs_detached_betas_refitted");
  //    directorylist.push_back("Tuple_JpsiKs_B2XMuMu");
  //    directorylist.push_back("Tuple_JpsiKs_B2XMuMu_refitted");
  //} else {
  // directorylist.push_back("JpsiL0Tuple_betas");
  //  directorylist.push_back("JpsiL0Tuple_B2XMuMu");
  //  directorylist.push_back("JpsiKSTuple_prescaled_betas");
  //  directorylist.push_back("JpsiKSTuple_detached_betas");
  //  directorylist.push_back("JpsiKSTuple_B2XMuMu");

  varnames(false, Lambdabname,Lambda0name,Jpsiname);

  //  const TString readerddxml("weights_bkgrej-dd-20120621-1225/TMVA_BDT.weights.xml");
  //  const TString readerllxml("weights_bkgrej-ll-20120621-1225/TMVA_BDT.weights.xml");

  const TString readerddxml("weights_bkgrej-dd-20120719-1115/TMVA_BDT.weights.xml");
  const TString readerllxml("weights_bkgrej-ll-20120719-1115/TMVA_BDT.weights.xml");

  TMVA::Reader* readerll = new TMVA::Reader( "!Color:!Silent" );
  addvariables(readerll);
  readerll->BookMVA( "BDT", readerllxml ); 
  

  TMVA::Reader* readerdd = new TMVA::Reader( "!Color:!Silent" );    
  addvariables(readerdd);
  readerdd->BookMVA( "BDT", readerddxml ); 


  TString output(filename);
  output.ReplaceAll(".root",".bdt.root");
  TFile *f = new TFile(output,"recreate");
  
  for (unsigned int idl=0; idl<directorylist.size(); idl++) {

    TDirectory* mydirectory = (TDirectory*)myfile.Get(directorylist[idl]);
    if (!mydirectory) continue;
    //  std::cout << mydirectory << std::endl;
    TTree* mytree = (TTree*)mydirectory->Get("DecayTree");
  
    const bool isJpsiKS=directorylist[idl].Contains("JpsiKs");  
    varnames(isJpsiKS, Lambdabname,Lambda0name,Jpsiname);

    if (isJpsiKS) {
      pplusname="piplus";
    } else {
      pplusname="pplus";
    }
 
    mytree->SetBranchAddress(Lambdabname + "_DIRA_OWNPV",&Lb_DIRA_OWNPV.d); 
    mytree->SetBranchAddress(Lambdabname + "_IPCHI2_OWNPV",&Lb_IPCHI2_OWNPV.d);
    mytree->SetBranchAddress(Lambdabname + "_FDCHI2_OWNPV", &Lb_FDCHI2_OWNPV.d);
    
    mytree->SetBranchAddress(Lambdabname + "_TAU",&Lb_TAU.d); 
    //    mytree->SetBranchAddress(Lambdabname + "_M", &Lb_M.d);
    mytree->SetBranchAddress(Lambdabname + "_ENDVERTEX_CHI2",&Lb_ENDVERTEX_CHI2.d);
    mytree->SetBranchAddress(Lambdabname + "_ENDVERTEX_NDOF", &Lb_ENDVERTEX_NDOF);
    
    
    mytree->SetBranchAddress(Lambda0name + "_M", &L0_M.d);
    //    mytree->SetBranchAddress(Lambda0name + "_TAU", &L0_TAU.d);
    //    mytree->SetBranchAddress(Lambda0name + "_PT", &L0_PT.d);
    mytree->SetBranchAddress(Jpsiname + "_ENDVERTEX_CHI2", &Jpsi_ENDVERTEX_CHI2.d);
    mytree->SetBranchAddress(Jpsiname + "_ENDVERTEX_NDOF", &Jpsi_ENDVERTEX_NDOF);
    mytree->SetBranchAddress(Jpsiname + "_MM", &Jpsi_MM.d);

    mytree->SetBranchAddress("piminus_ProbNNk",  &piminus_ProbNNk.d);
    mytree->SetBranchAddress("piminus_ProbNNp",  &piminus_ProbNNp.d);
    mytree->SetBranchAddress("piminus_ProbNNpi", &piminus_ProbNNpi.d);
    mytree->SetBranchAddress("piminus_TRACK_Type", &piminus_TRACK_Type);

    mytree->SetBranchAddress(pplusname + "_ProbNNk",  &pplus_ProbNNk.d);
    mytree->SetBranchAddress(pplusname + "_ProbNNp",  &pplus_ProbNNp.d);
    mytree->SetBranchAddress(pplusname + "_ProbNNpi", &pplus_ProbNNpi.d);



    TTree *T = new TTree(directorylist[idl],"test friend trees");
    T->Branch("bdt", &bdt, "bdt/F"); 
    
    const Long64_t nEntries=mytree->GetEntries();
    progressbar pbar(std::string("Calculating BDT for " + filename + ":" + directorylist[idl]));    
    for (Long64_t i=0; i<nEntries; i++) {
      mytree->GetEvent(i);
      Lb_ENDVERTEX_ratio=Lb_ENDVERTEX_CHI2.d/Lb_ENDVERTEX_NDOF;
      Jpsi_ENDVERTEX_ratio=Jpsi_ENDVERTEX_CHI2.d/Jpsi_ENDVERTEX_NDOF;

      Lb_FDCHI2_LOG=log10(Lb_FDCHI2_OWNPV.d);
      Lb_DIRA_LOG=log10(1.0000001-Lb_DIRA_OWNPV.d);

      L0_DM=fabs(L0_M.d-getLambda0mass(isJpsiKS));
      Jpsi_DMM=fabs(Jpsi_MM.d-Jpsipdgmass);
      pion_DLL_Kpi=piminus_ProbNNk.d-piminus_ProbNNpi.d;
      pion_DLL_ppi=piminus_ProbNNp.d-piminus_ProbNNpi.d;
      proton_DLL_Kp=pplus_ProbNNk.d-pplus_ProbNNp.d;
      proton_DLL_pip=pplus_ProbNNpi.d-pplus_ProbNNp.d;
      
      //      if (isJpsiKS) {
      //	L0_DM/=2.;
      //	pion_DLL_Kpi=-1.;
      //	pion_DLL_ppi=-1.;
      //	proton_DLL_Kp=-0.6;
      //      proton_DLL_pip=-0.6;
      //      }

      
      //transform double to float...
      Lb_IPCHI2_OWNPV.f=Lb_IPCHI2_OWNPV.d;
      Lb_TAU.f=Lb_TAU.d;


      if (piminus_TRACK_Type==3) {
	bdt=readerll->EvaluateMVA("BDT");
      } else if (piminus_TRACK_Type==5) {
	bdt=readerdd->EvaluateMVA("BDT");
      } else {
	std::cerr << "piminus_TRACK_Type???" << std::endl;
	exit(1);
      }

      pbar.print(100.*i/nEntries);
      T->Fill();
    }
    T->Write();
    pbar.finish();
  }
  f->Close();
}

