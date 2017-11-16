#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TEventList.h"
#include "TMVA/Factory.h"
#include "TMVA/Config.h"
#include "functions.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/MethodBDT.h"

#include "cuts.h"


void usage(char* argv0)
{
  std::cerr << argv0 << " Lb/B0 LL/DD" << std::endl;
  exit(1);
}

int main(int argc, char** argv)
{
  bool userwmc=false;
  bool usepid=true;
  bool noothervars=true;

  if (argc!=3) usage(argv[0]);

  const TString mode(argv[1]);
  if (mode!="B0" and mode!="Lb") usage(argv[0]);
  const bool doBd2JpsiKS0=mode.Contains("B0");

  TString LbLatex="#Lambda_{b}^{0}";
  if (doBd2JpsiKS0) LbLatex="B^{0}";
  TString L0Latex="#Lambda^{0}";
  if (doBd2JpsiKS0) L0Latex="K_{S}^{0}";
  TString JpsiLatex="J/#psi";
  

  TString llordd(argv[2]);
  llordd.ToLower();
  if (llordd!="dd" and llordd!="ll") usage(argv[0]);
  //  const bool doLL=llordd.Contains("ll");

  int doLL;
  if (llordd=="ll") {
    doLL=0;
  } else if (llordd=="dd") {
    doLL=1;
  } else {
    doLL=2;
  }

  //  std::cout << nicedate() << std::endl;
  //  return 0;

  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0, Lambdabname, Lambda0name, Jpsiname);
  const bool debug=true;

  TString sigfilename("DVTuples-Lb2JpsiL0-MC11a.reduced.root");
  const TString bkgfilename("DVTuples_data_stripping17.reduced.root");
  TString dirname("Tuple_JpsiL0_betas");
  if (doBd2JpsiKS0) {
    sigfilename="DVTuples-Bd2JpsiKS-MC11a.reduced.root";
    dirname="Tuple_JpsiKs_detached_betas";
  }


  TFile sigfile(sigfilename);
  TDirectory* sigdirectory = (TDirectory*)sigfile.Get(dirname);
  if (!sigdirectory) exit(1);
  TTree* sigtree = (TTree*)sigdirectory->Get("DecayTree");
  sigtree->AddFriend(dirname,TString(sigfilename).ReplaceAll(".root",".angles.root"));
  if (userwmc) sigtree->AddFriend(dirname,TString(sigfilename).ReplaceAll(".root",".reweighted.root"));

  if (debug) std::cerr << __FILE__ << " " <<  __LINE__ << std::endl;

  TFile bkgfile(bkgfilename);
  TDirectory* bkgdirectory = (TDirectory*)bkgfile.Get(dirname);
  TTree* bkgtree = (TTree*)bkgdirectory->Get("DecayTree");
  bkgtree->AddFriend(dirname,TString(bkgfilename).ReplaceAll(".root",".angles.root"));
  bkgtree->AddFriend(dirname,TString(bkgfilename).ReplaceAll(".root",".missvar.root"));

  // bkgfile.GetListOfKeys()->Print();
  // bkgdirectory->GetListOfKeys()->Print();
  // TTree* bkgtree1 = (TTree*)bkgdirectory->Get("DecayTree;1");
  // TTree* bkgtree2 = (TTree*)bkgdirectory->Get("DecayTree;2");

  // std::cout << bkgtree1 << " : "  <<bkgtree2 << std::endl;


  // TTree* sigtree1 = (TTree*)sigdirectory->Get("DecayTree;1");
  // TTree* sigtree2 = (TTree*)sigdirectory->Get("DecayTree;2");

  // std::cout << sigtree1 << " : "  <<sigtree2 << std::endl;


  //  bkgtree1->Print();
  //  bkgtree2->Print();

  if (debug) std::cerr << __FILE__ << " " <<  __LINE__ << std::endl;

  const TCut TriggerCut = triggercut(doBd2JpsiKS0);
  const TCut pid        = TCut("1==1");
  //  const TCut pid        = pidcut(doBd2JpsiKS0);
  TCut cutall(TriggerCut + pid + TCut(Lambdabname + "_DIRA_OWNPV>0.99&&" + Lambdabname + "_TAU>0.&&" + Lambdabname + "_IPCHI2_OWNPV<30."));
  if (doLL==0) cutall+=Lambda0_LL;
  if (doLL==1) cutall+=Lambda0_DD;


  
  //  const TCut cutsignal("1==1");
  TCut cutsignal=truemcsel(doBd2JpsiKS0);
  if (userwmc) cutsignal+=TCut("rw==1");
  TCut cutbkg(TString(Lambdabname + "_M>5800."));
  if (doBd2JpsiKS0) cutbkg=TString(Lambdabname + "_M>5400.");
  
  // //  TTree* sigtreered = sigtree->CopyTree(cutsignal);
  // TEventList sigmctrue("sigmctrue","sigmctrue");
  // sigtree->Draw(">>sigmctrue",cutsignal);
  // std::cout << "sigtree nevents = " << sigtree->GetEntries() << std::endl;
  // sigtree->SetEventList(&sigmctrue);
  // std::cout << "sigtree nevents = " << sigtree->GetEntries() << std::endl;
  // //  TCanvas c("c","c",100,100);
  // //  sigtree->Draw("abs(1.*lambda_b0_TRUEID)");
  // //  c.SaveAs("fuck.eps");


  std::cout << "cutall = " << cutall << std::endl;
  std::cout << "cutsig = " << cutsignal << std::endl;
  std::cout << "cutbkg = " << cutbkg << std::endl;
  //  return 1;


  TString outputname("bkgrej-");
  if (doLL==0) outputname+="ll-";
  if (doLL==1) outputname+="dd-";
  if (doLL==2) outputname+="all-";
  outputname+=nicedate();
  const TString outfileName("TMVA-" + outputname + ".root");
  TFile* outputFile = TFile::Open(outfileName,"RECREATE");
  if (outputFile->IsZombie() or (!outputFile->IsOpen())) {
    std::cerr << outfileName << " already exists (?!)" << std::endl;
    return 1;
  }

  TMVA::Factory* factory = new TMVA::Factory("TMVA", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
					     //Form("!V:%sColor", gROOT->IsBatch()?"!":"") );

  (TMVA::gConfig().GetIONames()).fWeightFileDir = "weights_" + outputname;

  if (debug) std::cerr << __FILE__ << " " <<  __LINE__ << std::endl;

  const double signalWeight=1.0;
  const double backgroundWeight=1.0;

  factory->AddSignalTree    (sigtree,      signalWeight);
  factory->AddBackgroundTree(bkgtree, backgroundWeight);

  factory->PrepareTrainingAndTestTree(cutall + cutsignal, cutall + cutbkg,
				      "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

  //				      "nTrain_Signal=10000:nTrain_Background=20000:nTest_Signal=9000:nTest_Background=40000:SplitMode=Random:NormMode=NumEvents:!V");

  if (!noothervars) {

    //  factory->AddVariable(Lambdabname + "_DIRA_OWNPV", 'F');
    factory->AddVariable("log10(1.0000001-" + Lambdabname + "_DIRA_OWNPV)", LbLatex + ": log_{10}(1-DIRA)", "", 'F'); 
    factory->AddVariable(Lambdabname + "_IPCHI2_OWNPV", LbLatex + ": IP #chi^{2}", "", 'F');
    factory->AddVariable("log10(" + Lambdabname + "_FDCHI2_OWNPV)",LbLatex + ": FD log_{10}(#chi^{2})", "", 'F');
    factory->AddVariable(Lambdabname + "_ENDVERTEX_CHI2/" + Lambdabname + "_ENDVERTEX_NDOF",  LbLatex + ": vtx #chi^{2}/N_{dof.}", "", 'F');
    factory->AddVariable(Lambdabname +  "_TAU", LbLatex + ": #tau", "", 'F');
    //  factory->AddVariable(Lambdabname +  "_PT", 'F');
    
    factory->AddVariable("abs(" + Lambda0name + "_M-" + str(getLambda0mass(doBd2JpsiKS0)) + ")",  L0Latex + ": |M(p#pi)-m(#Lambda)|","",'F');
    //  factory->AddVariable(Lambda0name + "_PT", 'F');
    //  factory->AddVariable(Lambda0name + "_TAU", 'F');
    
    factory->AddVariable("abs(" + Jpsiname + "_MM-" + str(Jpsipdgmass) + ")", JpsiLatex + ": |M(#mu#mu) - m(J/#psi)|", "", 'F');
    factory->AddVariable(Jpsiname + "_ENDVERTEX_CHI2/" + Jpsiname + "_ENDVERTEX_NDOF", JpsiLatex + ": vtx #chi^{2}/N_{dof.}", "", 'F');
    
    //  factory->AddVariable("piminus_ProbNNk-piminus_ProbNNpi");
    //  factory->AddVariable("piminus_ProbNNp-piminus_ProbNNpi");
    
    //  factory->AddVariable("pplus_ProbNNk-pplus_ProbNNpi");
    //  factory->AddVariable("pplus_ProbNNp-pplus_ProbNNpi");
  }

  if (usepid) {
    //    factory->AddVariable("log(piminus_ProbNNk)-log(piminus_ProbNNpi)");
    //    factory->AddVariable("log(piminus_ProbNNp)-log(piminus_ProbNNpi)");
    
    //    factory->AddVariable("log(pplus_ProbNNk)-log(pplus_ProbNNp)");
    //    factory->AddVariable("log(pplus_ProbNNpi)-log(pplus_ProbNNp)");

    factory->AddVariable("piminus_PIDK","#pi: PID(K)","",'F');
    factory->AddVariable("piminus_PIDp","#pi: PID(p)","",'F');

    factory->AddVariable("pplus_PIDK","p: PID(K)","",'F');
    factory->AddVariable("pplus_PIDp","p: PID(p)","",'F');

    //    factory->AddVariable("piminus_ProbNNk-piminus_ProbNNpi");
    //    factory->AddVariable("piminus_ProbNNp-piminus_ProbNNpi");
    
    //    factory->AddVariable("pplus_ProbNNk-pplus_ProbNNp");
    //    factory->AddVariable("pplus_ProbNNpi-pplus_ProbNNp");
  }

  //  factory->AddSpectator(Lambdabname + "_M", 'F');

  // factory->AddVariable("costheta", 'F');
  // factory->AddVariable("costheta1", 'F');
  // factory->AddVariable("costheta2", 'F');
  // factory->AddVariable("phi1", 'F');
  // factory->AddVariable("phi2", 'F');
  //  factory->AddSpectator("costheta", 'F');
  //  factory->AddSpectator("costheta1", 'F');
  //  factory->AddSpectator("costheta2", 'F');
  //  factory->AddSpectator("phi1", 'F');
  //  factory->AddSpectator("phi2", 'F');

  if (doLL==0) {
    factory->BookMethod( TMVA::Types::kBDT, "BDT",
			 "!H:!V:NTrees=50:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=10:PruneMethod=NoPruning" );
  } else if (doLL==1) {
    factory->BookMethod( TMVA::Types::kBDT, "BDT",
			 "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=10:PruneMethod=NoPruning" );
  } else {
    factory->BookMethod( TMVA::Types::kBDT, "BDT",
			 "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
  }
  //  factory->BookMethod( TMVA::Types::kCuts, "Cuts",
  //  		       "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
  //  factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
  //		       "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale:VarProp=FSmart" );


  //  factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
  //    		       "H:!V:FitMethod=GA:EffSel:Steps=40:Cycles=3:PopSize=1000:SC_steps=10:SC_rate=5:SC_factor=0.95");//:VarProp=FSmart");//:EffSMin=0.7:EffSMax=1.0" );  
  


  //  factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
  factory->TrainAllMethods();
  
  if (debug) std::cerr << __FILE__ << " " <<  __LINE__ << std::endl;
  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();
   
  if (debug) std::cerr << __FILE__ << " " <<  __LINE__ << std::endl;
  
  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();    
  
  if (debug) std::cerr << __FILE__ << " " <<  __LINE__ << std::endl;
  // --------------------------------------------------------------

  TMVA::MethodCuts* cuts = dynamic_cast<TMVA::MethodCuts*>(factory->GetMethod("CutsGA"));
  TMVA::MethodBDT* bdt = dynamic_cast<TMVA::MethodBDT*>(factory->GetMethod("BDT"));
  if (cuts) {
    for (Double_t eff=0.88; eff<=0.97; eff += 0.01) {
      cuts->PrintCuts( eff );
    }
  }

  if (bdt) {
    
  }

  // Save the output
  outputFile->Close();
  
  if (debug) std::cerr << __FILE__ << " " <<  __LINE__ << std::endl;  
  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAnalysis is done!" << std::endl;      
  
  delete factory;

}
