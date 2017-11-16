#include <iostream>
#include "stdlib.h"

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TRandom3.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TSpline.h"
#include "TEventList.h"
#include "TLatex.h"
#include "TH2Poly.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "Hoption.h"
#include "Hparam.h"
#include "TPainter3dAlgorithms.h"
#include "TCutG.h"
#include "THistPainter.h"

#include "functions.h"
#include "cuts.h"
#include "rootlogon.h"

#include "functions-roofit.h"
#include "progressbar.h"
#include "ConfigFile.h"

const bool debug=false;
const bool applybdtcut=true;
const bool applytmvarectcut=false;
const bool applysimplecut=false;

const bool userandom=false;


struct var12 {
  TString var1;
  TString var2;
  TString nvar1;
  TString nvar2;
  bool mode;
  TString dirname;
};




TH1* copyandrandom(TH1* corrhist, TRandom3* rnd, int i)
{
  if (debug) std::cout << "Making " << i << "th copy of " << corrhist->GetName() << std::endl; 

  TH1* result = (TH1*)corrhist->Clone(corrhist->GetName() + istr(i));
  for (int ibin=1; ibin<=result->GetBin(10000,10000,10000); ibin++) {
    double histoval=result->GetBinContent(ibin);
    double histoerr=result->GetBinError(ibin);
    result->SetBinContent(ibin, rnd->Gaus(histoval,histoerr));
  }
  const double maximum=result->GetMaximum();
  //  std::cout << "Max = " << maximum << std::endl;
  result->Scale(1./maximum);
  return result;
}


TH1F* loadcorr(bool doBd2JpsiKS0, const TString& varname, const TString& dirname, const TString& llordd)
{
  TCanvas c("c","c",100,100);

  TString filename="compmcdata-mcnonrw-";
  if (doBd2JpsiKS0) {
    filename+="Bd2JpsiKS0-";
  } else {
    filename+="Lb2JpsiL0-";
  }
  filename+=dirname;

  if (applybdtcut) {
    filename+="-bdtapplied";
  }
  if (applytmvarectcut) {
    filename+="-tmvaapplied";
  }
  if (applysimplecut) {
    filename+="-simplecutapplied";
  }
  filename+="-";
  filename+=llordd;
  filename+="-datavsmc-";
  filename+=varname;
  filename+=".root";

  TFile afile(filename);
  if (afile.IsZombie()) {
    std::cerr << "loadcorr: can't open " << filename << std::endl;
    exit(1);
  }
  if (!afile.Get("histo1nb") or !afile.Get("histo2nb")) {
    std::cerr << "loadcorr: can't find histograms in " << filename << std::endl;
    exit(2);
  }
  TH1F histo1(*(TH1F*)afile.Get("histo1nb"));
  TH1F histo2(*(TH1F*)afile.Get("histo2nb"));
  afile.Close();

  // histo1.Draw();
  // c.SaveAs("tmp-histo1-" + varname + "-" + llordd + ".eps");
  // histo2.Draw();
  // c.SaveAs("tmp-histo2-" + varname + "-" + llordd + ".eps");


  //   histo1.Smooth(5);
  //   histo2.Smooth(5);


  const bool manually=false;
  TH1F* corr;
  if (manually) {
    
    corr = new TH1F("corr-"+varname,"corr-"+varname, histo1.GetNbinsX(),histo1.GetXaxis()->GetXbins()->GetArray());
    //  corr.Divide(histo2);
    
    //  TSpline5 spline(&corr);
    
    double corrval[corr->GetNbinsX()];
    double maximum=-1000.;
    for (int i=1; i<=corr->GetNbinsX(); i++) {
      corrval[i-1]=histo1.GetBinContent(i)/histo2.GetBinContent(i);
      if (corrval[i-1]>maximum) maximum=corrval[i-1];
    }
    
    //corr.Scale(1./maximum);
    std::cout << "maximum = " << maximum << std::endl;
    
    for (int i=1; i<=corr->GetNbinsX(); i++) {
      corr->SetBinContent(i,corrval[i-1]/maximum);
      std::cout << i << " " << histo1.GetBinContent(i) << " " <<  histo2.GetBinContent(i) << " " << corr->GetBinContent(i) << std::endl;
    }
  } else {
    //    corr = new TH1F(*histo2);
    corr = (TH1F*)histo1.Clone("ratio-"+varname);
    corr->Divide(&histo2);
    const double maximum=corr->GetMaximum();
    std::cout << "Max = " << maximum << std::endl;
    corr->Scale(1./maximum);
    // for (int i=1; i<=corr->GetNbinsX(); i++) {
    //   corr->SetBinContent(i,corr->GetBinContent(i)/maximum);  
    //   std::cout << i << "/" << corr->GetNbinsX() << " " << histo1.GetBinContent(i) << " " <<  histo2.GetBinContent(i) << " " << corr->GetBinContent(i) << std::endl;
    // }
    // corr->Draw();
    // c.SaveAs("tmp-corr-" + varname + "-" + llordd + ".eps");
  }


  filename.ReplaceAll("compmcdata-mcnonrw-","");
  filename.ReplaceAll(".root","");
  corr->Smooth(1);
  corr->Draw("HIST");
  c.SaveAs("reweightmc-" + filename + ".eps");
  
  return corr;
}


TH2F* loadcorr2d(const var12& avar12, const TString& llordd)
{
  TCanvas c("c","c",100,100);

  TString filename="compmcdata-mcnonrw-";
  if (avar12.mode) {
    filename+="Bd2JpsiKS0-";
  } else {
    filename+="Lb2JpsiL0-";
  }
  filename+=avar12.dirname;

  if (applybdtcut) {
    filename+="-bdtapplied";
  }
  if (applytmvarectcut) {
    filename+="-tmvaapplied";
  }
  if (applysimplecut) {
    filename+="-simplecutapplied";
  }
  filename+="-";
  filename+=llordd;
  filename+="-datavsmc-";
  filename+=avar12.var1 + "_vs_" + avar12.var2;
  filename+=".root";

  TFile afile(filename);
  if (afile.IsZombie()) {
    std::cerr << "loadcorr2d: can't open " << filename << std::endl;
    exit(1);
  }
  if (!afile.Get("histo1nb") or !afile.Get("histo2nb")) {
    std::cerr << "loadcorr2d: can't find histograms in " << filename << std::endl;
    exit(2);
  }
  TH2F* h1src = (TH2F*)afile.Get("histo1nb");
  TH2F* h2src = (TH2F*)afile.Get("histo2nb");

  TH2F histo1(*(TH2F*)h1src->DrawNormalized());
  TH2F histo2(*(TH2F*)h2src->DrawNormalized());
  afile.Close();

  // histo1.Draw();
  // c.SaveAs("tmp-histo1-" + varname + "-" + llordd + ".eps");
  // histo2.Draw();
  // c.SaveAs("tmp-histo2-" + varname + "-" + llordd + ".eps");


  //   histo1.Smooth(5);
  //   histo2.Smooth(5);


  const bool manually=false;
  TH2F* corr;
  if (manually) {
    
    // corr = new TH1F("corr-"+varname,"corr-"+varname, histo1.GetNbinsX(),histo1.GetXaxis()->GetXbins()->GetArray());
    // //  corr.Divide(histo2);
    
    // //  TSpline5 spline(&corr);
    
    // double corrval[corr->GetNbinsX()];
    // double maximum=-1000.;
    // for (int i=1; i<=corr->GetNbinsX(); i++) {
    //   corrval[i-1]=histo1.GetBinContent(i)/histo2.GetBinContent(i);
    //   if (corrval[i-1]>maximum) maximum=corrval[i-1];
    // }
    
    // //corr.Scale(1./maximum);
    // std::cout << "maximum = " << maximum << std::endl;
    
    // for (int i=1; i<=corr->GetNbinsX(); i++) {
    //   corr->SetBinContent(i,corrval[i-1]/maximum);
    //   std::cout << i << " " << histo1.GetBinContent(i) << " " <<  histo2.GetBinContent(i) << " " << corr->GetBinContent(i) << std::endl;
    // }
  } else {
    //    corr = new TH1F(*histo2);
    corr = (TH2F*)histo1.Clone("ratio-"+avar12.var1 + "_vs_" + avar12.var2);
    corr->Divide(&histo2);
    const double maximum=corr->GetMaximum();
    std::cout << "Max = " << maximum << std::endl;
    corr->Scale(1./maximum);
    std::cout << "Max = " << corr->GetMaximum() << std::endl;
    // for (int i=1; i<=corr->GetNbinsX(); i++) {
    //   corr->SetBinContent(i,corr->GetBinContent(i)/maximum);  
    //   std::cout << i << "/" << corr->GetNbinsX() << " " << histo1.GetBinContent(i) << " " <<  histo2.GetBinContent(i) << " " << corr->GetBinContent(i) << std::endl;
    // }
    // corr->Draw();
    // c.SaveAs("tmp-corr-" + varname + "-" + llordd + ".eps");
  }


  filename.ReplaceAll("compmcdata-mcnonrw-","");
  filename.ReplaceAll(".root","");
  //  corr->Smooth(1);

  corr->GetXaxis()->SetTitle(avar12.nvar1);
  corr->GetYaxis()->SetTitle(avar12.nvar2);
  //  corr->SetMinimum(0.);
  corr->Draw("E");
  c.SaveAs("reweightmc-" + filename + ".eps");

  //  Int_t colors[] = {19,18,17,16,15,14,13,12};//, 5, 6}; // #colors >= #levels - 1
  //  gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
  gStyle->SetPalette(52,0);
  // #levels <= #colors + 1 (notes: +-1.79e308 = +-DBL_MAX; +1.17e-38 = +FLT_MIN)
  //  Double_t levels[] = {0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  //  Double_t levels[] = {0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0};
  //  corr->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
  
  gStyle->SetPaintTextFormat(".2f");
  corr->SetMarkerSize(.5);
  //  corr->Draw("AXIG");
  corr->Draw("COLZ");
  //  gPad->RedrawAxis("g");
  c.SaveAs("reweightmc-" + filename + "-COLZ.eps");
  corr->Draw("TEXT90E");
  //  myDraw(corr, "TEXT90E");
  c.SaveAs("reweightmc-" + filename + "-TEXT.eps");

  th2ftolatextable(corr, "reweightmc-" + filename + ".tex");
  

  return corr;
}


int main(int argc, char** argv)
{

  //  test();
  //  return 1;
  if (argc!=2) {
    std::cerr << argv[0] << " Lb/B0" << std::endl;
    exit(2);
  }
  const TString mode(argv[1]);
  const bool doBd2JpsiKS0=mode.Contains("B0");

  lhcbstyle();
  gStyle->SetPadRightMargin(0.15);

  std::vector<std::string> directorylist;
  if (doBd2JpsiKS0) {
    directorylist.push_back("Tuple_JpsiKs_detached_betas");
    //    directorylist.push_back("Tuple_JpsiKs_prescaled_betas");
  } else {
    directorylist.push_back("Tuple_JpsiL0_betas");
  }


  TString mcfilename("DVTuples-Lb2JpsiL0-MC11a-reco12a.reduced.root");
  if (doBd2JpsiKS0) {
    mcfilename = "DVTuples-Bd2JpsiKS-MC11a.reduced.root";
  }
  //  if (
  TFile myfile(mcfilename);

  TString anglefile(mcfilename);
  anglefile.ReplaceAll(".root",".angles.root");

  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0,Lambdabname,Lambda0name,Jpsiname);

  std::vector<TString> varstocorr;
  //  varstocorr.push_back(Lambdabname + "_PT");
  //  varstocorr.push_back(Lambda0name + "_PT");

  //  varstocorr.push_back(Lambdabname + "_IPCHI2_OWNPV");
  //  varstocorr.push_back(Lambda0name + "_ENDVERTEX_Z");
  //  if (!doBd2JpsiKS0) {
  //    varstocorr.push_back("piminus_y");  
  //  }
  //  varstocorr.push_back("pplus_PT");

  const TString Lambdablatex("#Lambda_{b}^{0}");  
  const TString Lambda0latex("#Lambda");
  const TString pilatex("#pi");

  std::vector<var12> varstocorr2d;
  var12 avar12;
  avar12.var1=Lambdabname + "_PT";
  avar12.var2=Lambdabname + "_y";
  avar12.nvar1=Lambdablatex + ": p_{T} (MeV)";
  avar12.nvar2=Lambdablatex + ": y";
  avar12.mode=doBd2JpsiKS0;
  avar12.dirname=directorylist[0];
  varstocorr2d.push_back(avar12);

  avar12.var1=Lambda0name + "_PT";
  avar12.var2=Lambda0name + "_y";
  avar12.nvar1=Lambda0latex + ": p_{T} (MeV)";
  avar12.nvar2=Lambda0latex + ": y";
  avar12.mode=doBd2JpsiKS0;
  avar12.dirname=directorylist[0];
  varstocorr2d.push_back(avar12);

  //  if (doBd2JpsiKS0) {
    avar12.var1="piminus_PT";
    avar12.var2="piminus_y";
    avar12.nvar1=pilatex + ": p_{T} (MeV)";
    avar12.nvar2=pilatex + ": y";
    avar12.mode=true;
    avar12.dirname="Tuple_JpsiKs_detached_betas";
    varstocorr2d.push_back(avar12);
    //  }




  TString outputname=mcfilename;
  outputname.ReplaceAll(".root",".reweighted.root");
  TFile *newfile = NULL;
  if (!debug) newfile = new TFile(outputname,"recreate");
  TRandom3 rnd(0);

  
  //  TCut LL_selection, DD_selection;
  //  simplecut(doBd2JpsiKS0, LL_selection, DD_selection);
  //  TCut selection = LL_selection or DD_selection;

  TCut TriggerCut = triggercut(doBd2JpsiKS0);
  TCut trueMC = truemcsel(doBd2JpsiKS0);

  const int maxtoy=200;

  for (unsigned int idir=0; idir<directorylist.size(); idir++) {
    const TString dirname = directorylist[idir];


    TTree* newtree=NULL;
    //    bool keep=true;

    bool keep[maxtoy];
    memset(keep,true,maxtoy*sizeof(bool));

    bool keepvec[varstocorr.size()];
    //    int rw=0;
    if (!debug) {
      //      newfile->mkdir(dirname);
      //      newfile->cd(dirname);
      //      newtree = mytree->CloneTree(0);
      newtree=new TTree(dirname,"test friend trees");
      newtree->Branch("rw",    &keep[0],    "rw/b");      
      for (int itoy=1;itoy<maxtoy;itoy++) newtree->Branch("rw" + istr(itoy),    &keep[itoy],    "rw" + istr(itoy) + "/b"); 
 
      for (unsigned int i=0; i<varstocorr.size();i++) {
	newtree->Branch("rw-" + varstocorr[i],    &keepvec[i],    "rw-" + varstocorr[i] + "/b"); 
      }
      //      newfile->Add(newtree);
    }


    TDirectory* mydir = (TDirectory*)myfile.Get(dirname);
    TTree* mytree = (TTree*)mydir->Get("DecayTree");
    mytree->AddFriend(dirname,anglefile);


    Float_t vars[varstocorr.size()];
    Float_t vars2d[varstocorr2d.size()][2];
    Int_t piminus_TRACK_Type;
    TH1F* corrhist[maxtoy][2][varstocorr.size()];
    TH2F* corrhist2d[maxtoy][2][varstocorr2d.size()];
    mytree->SetBranchAddress("piminus_TRACK_Type", &piminus_TRACK_Type);
    for (unsigned int i=0; i<varstocorr.size();i++) {
      mytree->SetBranchAddress(varstocorr[i], &vars[i]);
      corrhist[0][0][i]=loadcorr(doBd2JpsiKS0, varstocorr[i],dirname,"ll");
      corrhist[0][1][i]=loadcorr(doBd2JpsiKS0, varstocorr[i],dirname,"dd");

      for (int itoy=1;itoy<maxtoy;itoy++) {
	corrhist[itoy][0][i] = (TH1F*)copyandrandom(corrhist[0][0][i],&rnd,itoy);
	corrhist[itoy][1][i] = (TH1F*)copyandrandom(corrhist[0][1][i],&rnd,itoy);
      }

    }
    for (unsigned int i=0; i<varstocorr2d.size();i++) {
      mytree->SetBranchAddress(varstocorr2d[i].var1, &vars2d[i][0]);
      mytree->SetBranchAddress(varstocorr2d[i].var2, &vars2d[i][1]);
      corrhist2d[0][0][i]=loadcorr2d(varstocorr2d[i],"ll");
      corrhist2d[0][1][i]=loadcorr2d(varstocorr2d[i],"dd");

      for (int itoy=1;itoy<maxtoy;itoy++) {
	corrhist2d[itoy][0][i] = (TH2F*)copyandrandom(corrhist2d[0][0][i],&rnd,itoy);
	corrhist2d[itoy][1][i] = (TH2F*)copyandrandom(corrhist2d[0][1][i],&rnd,itoy);
      }

    }
    //    exit(1);


    TEventList myeventlist("myeventlist","myeventlist");
    mytree->Draw(">>myeventlist", TriggerCut + trueMC);

    Long64_t kept=0;
    const Long64_t nEvents=mytree->GetEntries();
    const Long64_t nEvents_triggered=myeventlist.GetN();
    progressbar pbar(std::string("Writing rw for " + dirname));
    for (Long64_t ievent=0; ievent<nEvents; ievent++) {
      if (myeventlist.Contains(ievent)) {
	mytree->GetEvent(ievent);
	//	keep=true;    
	memset(keep,true,maxtoy*sizeof(bool));

	for (unsigned int j=0; j<varstocorr.size(); j++) keepvec[j]=true;
	int index=1;
	if (piminus_TRACK_Type==3) index=0;
	
	for (int itoy=0;itoy<maxtoy;itoy++) {

	  //1D
	  for (unsigned int j=0; j<varstocorr.size(); j++) {
	    const double varsj=vars[j];
	    const double histoval=corrhist[itoy][index][j]->GetBinContent(corrhist[itoy][index][j]->FindBin(varsj));
	    const double histoerr=corrhist[itoy][index][j]->GetBinError(corrhist[itoy][index][j]->FindBin(varsj));
	    const double random=rnd.Rndm();
	    const double histornd=userandom?rnd.Gaus(histoval,histoerr):histoval;
	    //	  std::cout << histornd <<  " " << histoval << std::endl;
	    //	  std::cout << histoval << " " << histoerr << " " << histornd << std::endl;
	    //	    if (itoy==0)
	    keepvec[j]=(histornd>random);
	    keep[itoy]=keep[itoy] and keepvec[j];
	    //	  keep=keep and (histoval>random);
	    //	  if (debug) std::cout << "b0_eta = " << varsj << " ; histoval = " << histoval << " ; random = " << random << " ; keep = " << keep << std::endl;
	  }
	  
	  //2D
	  for (unsigned int j=0; j<varstocorr2d.size(); j++) {
	    const double varsj1=vars2d[j][0];
	    const double varsj2=vars2d[j][1];
	    //	  std::cout << varstocorr2d[j].var1 << " = " << vars2d[j][0] << std::endl;
	    //	  std::cout << varstocorr2d[j].var2 << " = " << vars2d[j][1] << std::endl;
	    const double histoval=corrhist2d[itoy][index][j]->GetBinContent(corrhist2d[itoy][index][j]->FindBin(varsj1,varsj2));
	    const double histoerr=corrhist2d[itoy][index][j]->GetBinError(corrhist2d[itoy][index][j]->FindBin(varsj1,varsj2));
	    const double random=rnd.Rndm();
	    const double histornd=userandom?rnd.Gaus(histoval,histoerr):histoval;
	    //	  std::cout << histornd <<  " " << histoval << std::endl;
	    //	  std::cout << histoval << " " << histoerr << " " << histornd << std::endl;
	    
	    //	  keepvec[j]=(histornd>random);
	    //	  keep=keep and keepvec[j];
	    
	    keep[itoy]=keep[itoy] and (histornd>random);
	    //	  keep=keep and (histoval>random);
	    //	  if (debug) std::cout << "b0_eta = " << varsj << " ; histoval = " << histoval << " ; random = " << random << " ; keep = " << keep << std::endl;
	  }
	}

      } else {
	//	keep=false;
	memset(keep,false,maxtoy*sizeof(bool));
      }

      if (keep[0]) kept++;
      if (!debug) {
	newtree->Fill();
      } else {
	// if (keep) {
	//   std::cout << "keep" << std::endl;
	// } else {
	//   std::cout << "throw away" << std::endl;
	// }
      }

      pbar.print(100.*ievent/nEvents);
    }
    if (!debug) {
      //      newtree->Write();
      //      delete newtree;
    }
    pbar.finish();
    std::cerr << "nEvents           = " << nEvents << std::endl;
    std::cerr << "nEvents_triggered = " << nEvents_triggered << std::endl;
    std::cerr << "kept              = " << kept << std::endl;
  }

  if (!debug) {
    newfile->Write();
    newfile->Close();
  }
  //  delete newfile;




  //compare costheta etc... 
  RooRealVar costheta("costheta",   "cos(#theta)",     -1., 1.);
  RooRealVar costheta1("costheta1", "cos(#theta_{1})", -1., 1.);
  RooRealVar costheta2("costheta2", "cos(#theta_{2})", -1., 1.);
  costheta.setBins(25);
  costheta1.setBins(25);
  costheta2.setBins(25);
  for (unsigned int idir=0; idir<directorylist.size(); idir++) {
    const TString dirname = directorylist[idir];


    TTree* newtree=NULL;
    bool keep[maxtoy];
    memset(keep,true,maxtoy*sizeof(bool));
    //    int rw=0;
    if (!debug) {
      //      newfile->mkdir(dirname);
      //      newfile->cd(dirname);
      //      newtree = mytree->CloneTree(0);
      newtree=new TTree(dirname,"test friend trees");
      newtree->Branch("rw",    &keep[0],    "rw/b"); 
      //      for (int i=1;i<maxtoy;i++) newtree->Branch("rw" + istr(i),    &keep[i],    "rw" + istr(i) + "/b"); 
      //      newfile->Add(newtree);
    }


    TDirectory* mydir = (TDirectory*)myfile.Get(dirname);
    if (!mydir) continue;
    TTree* mytree = (TTree*)mydir->Get("DecayTree");
    //    TString rwfile=mcfilename;
    //    outputname.ReplaceAll(".root",".reweighted.root");


    TString outputname="reweightmc-";
    if (doBd2JpsiKS0) {
      outputname+="Bd2JpsiKS0-";
    } else {
      outputname+="Lb2JpsiL0-";
    }
    outputname+=dirname;
   
    if (applybdtcut) {
      outputname+="-bdtapplied";
    }
    if (applytmvarectcut) {
      outputname+="-tmvaapplied";
    }
    
    if (applysimplecut) {
      outputname+="-simplecutapplied";
    }
    mytree->AddFriend(dirname,TString(mcfilename).ReplaceAll(".root",".reweighted.root"));
    mytree->AddFriend(dirname,TString(mcfilename).ReplaceAll(".root",".angles.root"));
    TCanvas c("c","c",1000,1000);  
    plotttree(c,RooArgList(costheta,costheta1,costheta2), mytree,
	      Lambda0_LL + trueMC + TriggerCut + TCut("rw!=1"), 
	      Lambda0_LL + trueMC + TriggerCut + TCut("rw==1"), 
	      outputname + "-ll-rwcheck", false, "Rejected", "Accepted");
    plotttree(c,RooArgList(costheta,costheta1,costheta2), mytree,
	      Lambda0_DD + trueMC + TriggerCut + TCut("rw!=1"), 
	      Lambda0_DD + trueMC + TriggerCut + TCut("rw==1"), 
	      outputname + "-dd-rwcheck", false, "Rejected", "Accepted");

  }

  
  myfile.Close();
  return 0;

//   //compare




//   //  std::cout << "outputname = " << outputname << std::endl;
//   //  TFile* rwfile = new TFile(outputname);
//   //  std::cerr << __LINE__ << std::endl;


//   //PT
//   const double varmin=0.;
//   const double varmax=30000.;

//   //eta
//   //  const double varmin=0.;
//   //  const double varmax=10.;


//   TFile myfile2(inputfile);
  
//   for (unsigned int idir=0; idir<directorylist.size(); idir++) {
//     const TString dirname = directorylist[idir];
//     TDirectory* mydir = (TDirectory*)myfile2.Get(dirname);
//     //    std::cerr << __LINE__ << std::endl;
//     TTree* mytree = (TTree*)mydir->Get("DecayTree");
//     //    std::cerr << __LINE__ << std::endl;
//     mytree->AddFriend(dirname,outputname);
//     //    mytree->AddFriend(dirname,anglefile);

//     //    Int_t piminus_TRACK_Type;
//     //    mytree->SetBranchAddress("piminus_TRACK_Type", &piminus_TRACK_Type);
    
//     for (unsigned int j=0; j<varstocorr.size(); j++) {


//       TH1F fhisto1("fhisto1","fhisto1",25,varmin,varmax);
//       TH1F fhisto2("fhisto2","fhisto2",25,varmin,varmax);

//       mytree->Draw(varstocorr[j]+">>fhisto1","piminus_TRACK_Type==3");
//       mytree->Draw(varstocorr[j]+">>fhisto2","piminus_TRACK_Type==3&&rw==1");

//       TFile histofile("Bd2JpsiKS0-" + dirname + "-ll-datavsmc-" + varstocorr[j] + ".root");
//       TH1F histo3(*(TH1F*)histofile.Get("histo1"));
//       TH1F histo4(*(TH1F*)histofile.Get("histo2"));

//       histo3.SetLineColor(kBlue);
//       histo4.SetLineColor(kRed);

//       //      fhisto2.SetMaximum(0.2);
//       fhisto1.SetLineColor(kRed);
//       fhisto2.SetLineColor(kBlue);


//       fhisto2.DrawNormalized();
//       fhisto1.DrawNormalized("SAME");
//       histo3.DrawNormalized("SAME");
//       histo4.DrawNormalized("SAME");
      
      
//       //      std::cerr << __LINE__ << std::endl;
//       //      mytree->Draw(varstocorr[j],"","SAME");
//       //      std::cerr << __LINE__ << std::endl;
//       c.SaveAs("rw_" + varstocorr[j] + ".eps");
//     }
//   }
}
