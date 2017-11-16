#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TEventList.h"
#include "TMVA/Factory.h"
#include "TMVA/Config.h"


int main(int argc, char** argv)
{

  const TString Lambdabname("lambda_b0");
  const TString Lambda0name("lambda0");
  const TString Jpsiname("jpsi1s");
  
  const bool debug=true;

  TFile sigfile("DVTuples-Lb2JpsiL0-MC11a.truemc.root");

  const TString dirname("JpsiL0Tuple_betas");
  
  TDirectory* sigdirectory = (TDirectory*)sigfile.Get(dirname);
  TTree* sigtree = (TTree*)sigdirectory->Get("DecayTree");

  sigtree->AddFriend(dirname,"DVTuples-Lb2JpsiL0-MC11a.truemc.angles.root");

  if (debug) std::cerr << __FILE__ << " " <<  __LINE__ << std::endl;

  TFile bkgfile("DVTuples_data_stripping17.root");
  TDirectory* bkgdirectory = (TDirectory*)bkgfile.Get(dirname);
  TTree* bkgtree = (TTree*)bkgdirectory->Get("DecayTree");

  bkgtree->AddFriend(dirname,"DVTuples_data_stripping17.angles.root");

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

  //  const TCut cutsignal(TString("abs(" + Lambdabname + "_TRUEID)==5122"));


  const TCut cutall(Lambdabname + "_DIRA_OWNPV>0.7&&" + Lambdabname + "_TAU>0.&&" + Lambdabname + "_IPCHI2_OWNPV<30.&&" + Lambda0name + "_TAU>0.");
  //  const TCut cutsignal("1==1");
  //  const TCut cutsignal("abs(1.*lambda_b0_TRUEID)==5122");
  const TCut cutbkg(TString(Lambdabname + "_M>5800."));

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
  std::cout << "cutbkg = " << cutbkg << std::endl;
  //  return 1;


  const TString outputname("bkgrej-test");
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

  factory->PrepareTrainingAndTestTree(cutall, cutall + cutbkg,
				      "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

  //				      "nTrain_Signal=10000:nTrain_Background=20000:nTest_Signal=9000:nTest_Background=40000:SplitMode=Random:NormMode=NumEvents:!V");

  factory->AddVariable(Lambdabname + "_DIRA_OWNPV", 'F');
  factory->AddVariable(Lambdabname + "_IPCHI2_OWNPV", 'F');
  factory->AddVariable(Lambdabname + "_FDCHI2_OWNPV", 'F');
  factory->AddVariable(Lambdabname + "_ENDVERTEX_CHI2/" + Lambdabname + "_ENDVERTEX_NDOF", 'F');
  factory->AddVariable(Lambdabname +  "_TAU", 'F');
  factory->AddVariable("abs(" + Lambda0name + "_M-1115.683)", 'F');

  factory->AddVariable(Lambda0name + "_PT", 'F');
  factory->AddVariable(Lambda0name + "_TAU", 'F');

  factory->AddVariable(Jpsiname + "_ENDVERTEX_CHI2/" + Jpsiname + "_ENDVERTEX_NDOF", 'F');

  //  factory->AddVariable("piminus_ProbNNk-piminus_ProbNNpi");
  //  factory->AddVariable("piminus_ProbNNp-piminus_ProbNNpi");

  //  factory->AddVariable("pplus_ProbNNk-pplus_ProbNNpi");
  //  factory->AddVariable("pplus_ProbNNp-pplus_ProbNNpi");

  factory->AddVariable("log(piminus_ProbNNk)-log(piminus_ProbNNpi)");
  factory->AddVariable("log(piminus_ProbNNp)-log(piminus_ProbNNpi)");

  factory->AddVariable("log(pplus_ProbNNk)-log(pplus_ProbNNp)");
  factory->AddVariable("log(pplus_ProbNNpi)-log(pplus_ProbNNp)");

  factory->AddSpectator(Lambdabname + "_M", 'F');

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

  factory->BookMethod( TMVA::Types::kBDT, "BDT",
		       "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
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
  
  // Save the output
  outputFile->Close();
  
  if (debug) std::cerr << __FILE__ << " " <<  __LINE__ << std::endl;  
  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAnalysis is done!" << std::endl;      
  
  delete factory;

}
