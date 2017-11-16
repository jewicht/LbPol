#include "functions-roofit.h"

bool checkweirdpulls(RooFitResult* fr, RooArgList* vars)
{
  for (int i=0;i<vars->getSize();i++) {
    RooRealVar* vartrue = (RooRealVar*)vars->at(i);
    const RooArgList& fpars = fr->floatParsFinal();
    RooRealVar* varmeas = (RooRealVar*)fpars.find(vars->at(i)->GetName());

    if (vartrue and varmeas) {
      RooPullVar pull("pull","pull",*varmeas,*vartrue);
      
      const double pullval = pull.getVal();
      
      if ((fabs(pullval)>7.) or (pullval!=pullval) or isinf(pullval) or isnan(pullval) )  {
	std::cout << "Weird pull for " << vars->at(i)->GetName() << ": pull=" << pullval << std::endl;
	fr->Print("v");
	mygetchar;
	return false;
      }
    }
  }
  return true;
}


inline std::string finoupas(int i, int j)
{
  if (i==j) return " \\\\ ";
  return " & ";
}

void th2ftolatextable(TH2F* histo, const TString& outputname)
{

  const Int_t dimx=histo->GetNbinsX();
  const Int_t dimy=histo->GetNbinsY();

  std::fstream fileout;
  fileout.open(outputname, std::fstream::out);

  //write header
  fileout << "\\begin{tabular}{l";
  for (Int_t iy=1; iy<=dimy; iy++) fileout << "c";
  fileout << "}" << std::endl;
  fileout << " & ";
  for (Int_t iy=1; iy<=dimy; iy++)	{
    fileout << "$ " << histo->GetYaxis()->GetBinLowEdge(iy) << " - " << histo->GetYaxis()->GetBinLowEdge(iy+1) << " $" << finoupas(iy,dimy);
  }
  fileout << " \\hline" << std::endl;

  //write table
  for (Int_t ix=1; ix<=dimx; ix++) {
    //    fileout << "$ " <<histo->GetYaxis()->GetBinLowEdge(histo->GetBin(1,iy)) << " - " << histo->GetYaxis()->GetBinLowEdge(histo->GetBin(1,iy+1)) << " $ & ";
    fileout << "$ " <<histo->GetXaxis()->GetBinLowEdge(ix) << " - " << histo->GetXaxis()->GetBinLowEdge(ix+1) << " $ & ";
    for (Int_t iy=1; iy<=dimy; iy++) {  
      char buffer1[50];
      sprintf(buffer1, "$ %.3f \\pm %.3f $", histo->GetBinContent(histo->GetBin(ix,iy)), histo->GetBinError(histo->GetBin(ix,iy)));
      fileout << TString(buffer1) << finoupas(iy,dimy);
    }
    if (ix==dimx) fileout << " \\hline";
    fileout << std::endl;
  }  

  //write header
  fileout << "\\end{tabular}" << std::endl;
  fileout.close();
}


RooArgList* addvariation(RooDataSet* data, const RooArgList& varstoplot)
{

  RooArgList* variationlist = new RooArgList();

  double_t val0[varstoplot.getSize()];

  for (int i=0;i<varstoplot.getSize();i++) {
    RooRealVar* var = (RooRealVar*)varstoplot.at(i);
    val0[i]=data->get(0)->getRealValue(var->GetName());
    RooRealVar* var2 = new RooRealVar(TString(var->GetName()) + "_var", "#delta(" + TString(var->GetTitle()) + ")",0.);
    variationlist->add(RooArgList(*var2));
  }
  
  mygetchar;

  RooDataSet dataww6("dataww6","dataww6",*variationlist);

  const int nEntries=data->numEntries();
  for (int i=0;i<nEntries;i++) {  
    const RooArgSet* blu = data->get(i);
    for (int j=0;j<varstoplot.getSize();j++) {
      RooRealVar* var = (RooRealVar*)variationlist->at(j);
      var->setVal( blu->getRealValue(varstoplot.at(j)->GetName()) - val0[j] );
      //      var->setVal( 100. * (val0[j]-blu->getRealValue(varstoplot.at(j)->GetName()))/val0[j]  );
    }
    dataww6.add(*variationlist);
  }
 data->merge(&dataww6);

 return variationlist;
}

void plotprofile(const TCanvas& c, const TString& output, RooDataSet* data, RooRealVar* var1, RooRealVar* var2) 
{
  TH2F histo("histo","histo",
	     50, var1->getMin(), var1->getMax(),
	     50, var2->getMin(), var2->getMax()
	     );

  data->fillHistogram(&histo,RooArgList(*var1,*var2));
  
  TProfile* prof=histo.ProfileX();
  
  prof->GetXaxis()->SetTitle(TString(var1->GetTitle()) + " vs " + var2->GetTitle() );

  prof->Draw();
  c.SaveAs(output + "-" + var1->GetName() + "vs" + var2->GetName() + ".eps");

  delete prof;
}

void plotprofile(const TCanvas& c, TTree* sigtree, const TString& var1, const TString& var2, const TCut& selection , const TString& plotname)
{
   TProfile prof("prof","prof",25,-1.,1.);
   sigtree->Draw(var1 + ":" + var2 + ">>prof", selection);
   //   prof.Fit("pol1");
   prof.GetXaxis()->SetTitle(var1 + " vs " + var2 );
   prof.Draw();

   TString var1tmp(var1);
   TString var2tmp(var2);
   var1tmp.ReplaceAll("/","div");
   var2tmp.ReplaceAll("/","div");
   c.SaveAs(plotname + "-" + var1tmp + "-vs-" + var2tmp + ".eps");
}



void plotprofile(const TCanvas& c, TTree* sigtree, const std::vector<TString>& varlist1, const std::vector<TString>& varlist2, const TCut& selection , const TString& plotname)
{
  for (unsigned int i=0; i<varlist1.size(); i++) {
    for (unsigned int j=0; j<varlist2.size(); j++) {
      //      RooAbsArg* var1 = varlist.at(i);
      //      RooAbsArg* var2 = varlist.at(j); 
      if (varlist1[i]!=varlist2[j]) plotprofile(c,sigtree, varlist1[i], varlist2[j], selection, plotname);
      //      TProfile prof("prof","prof",25,-1.,1.);
      //      sigtree->Draw(TString(var1->GetName())+":"+var2->GetName() + ">>prof", selection);
      //      prof.Fit("pol1");
      //      prof.GetXaxis()->SetTitle(TString(var1->GetTitle())+" vs "+var2->GetTitle());
      //      prof.Draw();
      //      c.SaveAs(plotname + "-" + TString(var1->GetName()) + "-vs-" + var2->GetName() +".eps");
    }
  }
}

bool isconverged(RooFitResult* fr, bool debug)
{
  bool result=true;
  if (debug) std::cout << "status = " << fr->status() << " ; covQual = " << fr->covQual() <<  std::endl;
  for (unsigned int i=0; i<fr->numStatusHistory(); i++) {
    if (debug) std::cout << "i=" << i << " ; label = " << fr->statusLabelHistory(i) << " ; status = " << fr->statusCodeHistory(i) << std::endl;
    result = result and (fr->statusCodeHistory(i)==0);
  }
  result = result and (fr->covQual()==3);
  return result;
}


void myprintLatex(const RooArgList& mylist, const char* filename, Int_t ncol, const RooCmdArg* formatCmd)
{
  // Internal implementation function of printLatex

  std::fstream ofs(filename);

  // Count number of rows to print
  Int_t nrow = (Int_t) ((mylist.getSize()) / ncol) ;
  //  if (mylist.getSize()%2==1) nrow++;
  //  std::cout << "mylist.getSize() = " << mylist.getSize() << std::endl;
  //  std::cout << "ncol = " << ncol << std::endl;
  //  std::cout << "nrow = " << nrow << std::endl;

  Int_t i,j,k ;

  // Sibling list do not need to print their name as it is supposed to be the same  
  //  RooCmdArg sibFormatCmd = RooFit::Format("E",RooFit::AutoPrecision(2),RooFit::VerbatimName(true));
  
  //  sibFormatCmd = *formatCmd ;
  // TString tmp = formatCmd->_s[0] ;
  // tmp.ReplaceAll("N","") ;    
  // tmp.ReplaceAll("n","") ;    
  // static char buf[100] ;
  // strlcpy(buf,tmp.Data(),100) ;
  // sibFormatCmd._s[0] = buf ;


  // Make list of lists ;
  RooLinkedList listList ;
  listList.Add((RooAbsArg*)&mylist) ;
  //  RooFIter sIter = siblingList.fwdIterator() ;
  RooAbsCollection* col ;
  //  while((col=(RooAbsCollection*)sIter.next())) {
  //    listList.Add(col) ;
  //  }

  RooLinkedList listListRRV ;

  // Make list of RRV-only components
  RooFIter lIter = listList.fwdIterator() ;
  RooArgList* prevList = 0 ;
  while((col=(RooAbsCollection*)lIter.next())) {
    RooArgList* list = new RooArgList ;
    RooFIter iter = col->fwdIterator() ;
    RooAbsArg* arg ;
    while((arg=iter.next())) {    
      
      RooRealVar* rrv = dynamic_cast<RooRealVar*>(arg) ;
      if (rrv) {
	list->add(*rrv) ;
      } else {
	//	coutW(InputArguments) << "RooAbsCollection::printLatex: can only print RooRealVar in LateX, skipping non-RooRealVar object named "
	//	     << arg->GetName() << endl ;      
      }
      if (prevList && TString(rrv->GetName()).CompareTo(prevList->at(list->getSize()-1)->GetName())) {
	//	coutW(InputArguments) << "RooAbsCollection::printLatex: WARNING: naming and/or ordering of sibling list is different" << endl ;
      }
    }
    listListRRV.Add(list) ;
    if (prevList && list->getSize() != prevList->getSize()) {
      //      coutW(InputArguments) << "RooAbsCollection::printLatex: ERROR: sibling list(s) must have same length as self" << endl ;
      delete list ;
      listListRRV.Delete() ;
      return ;
    }
    prevList = list ;
  }

  // Construct table header
  Int_t nlist = listListRRV.GetSize() ;
  TString subheader = "l" ;
  for (k=0 ; k<nlist ; k++) subheader += "c" ;

  TString header = "\\begin{tabular}{" ;
  for (j=0 ; j<ncol ; j++) {
    if (j>0) header += "|" ;
    header += subheader ;
  }
  header += "}" ;
  ofs << header << std::endl ;


  // Print contents, delegating actual printing to RooRealVar::format()
  for (i=0 ; i<nrow ; i++) {
    for (j=0 ; j<ncol ; j++) {
      for (k=0 ; k<nlist ; k++) {
	RooRealVar* par = (RooRealVar*) ((RooArgList*)listListRRV.At(k))->at(i+j*nrow) ;
	if (par) {
	  TString* tmp = par->format(*formatCmd) ;
	  ofs << *tmp ;
	  delete tmp ;
	}
	if (!(j==ncol-1 && k==nlist-1)) {
	  ofs << " & " ;
	}
      }
    }
    ofs << "\\\\" << std::endl ;
  }
  
  ofs << "\\end{tabular}" << std::endl ;
  listListRRV.Delete() ;
}


RooFitResult* getRooFitResult(const TString& filename, TString frname)
{
  //  TFile inputfile("/afs/cern.ch/user/j/jwicht/fitLambdab2JpsiL0/" + filename);

  std::ifstream ifs;
  ifs.open(filename);

  TString realfilename=filename;
  if (!ifs.good()) {
    ifs.close();
    ifs.open("/afs/cern.ch/user/j/jwicht/fitLambdab2JpsiL0/" + filename);
    if (!ifs.good()) {
      std::cerr << "getRooFitResult: can't find " << filename <<  " in curr or work dir" << std::endl;
      ifs.close();
      return NULL;  
    }
    realfilename="/afs/cern.ch/user/j/jwicht/fitLambdab2JpsiL0/" + filename;
    std::cout << "getRooFitResult: opening " << filename << " from work dir" << std::endl; 
  } else {
    std::cout << "getRooFitResult: opening " << filename << " from curr dir" << std::endl; 
  }
  ifs.close();


  TFile* inputfile = new TFile(realfilename);
  if (frname=="") frname = inputfile->GetListOfKeys()->At(0)->GetName();
  RooFitResult* fr = (RooFitResult*)inputfile->Get(frname); 
  if (!fr) {
    std::cerr << "getRooFitResult: can't find " << frname << std::endl;
    return NULL;
  }
  RooFitResult* frclone = new RooFitResult(*fr);
  inputfile->Close();
  delete inputfile;
  return frclone;
}


void readvarsfromtfile(const RooArgList& varlist, const TString& filename, TString frname, bool randomize)
{
  RooFitResult* fr = getRooFitResult(filename, frname);
  if (!fr) {
    std::cerr << "readvarsfromtfile: can't find RooFitResult in " << filename << ":" << frname << std::endl;
    return;
  }
  const RooArgList* floatPars;
  if (randomize) {
    floatPars = &fr->randomizePars();
  } else {
    floatPars = &fr->floatParsFinal();
  }
  //  floatPars.Print("v");
  for (int i=0; i<varlist.getSize(); i++) {
    RooRealVar* myvar = (RooRealVar*)varlist.at(i);
    RooRealVar* myvar2 = (RooRealVar*)floatPars->find(myvar->GetName());
    if (myvar2) {
      myvar->setVal(myvar2->getVal());
      myvar->setError(myvar2->getError());
    } else {
      std::cerr << "readvarsfromtfile: can't find " << myvar->GetName() << " in " << filename << ":" << frname << std::endl;
      //      exit(1);
    }
  }
  delete fr;
}

void filldataset_copytree(RooDataSet*& data, const TString& name, const TString& title, const RooArgList& vars, TTree* mytree, const TCut& mycut, bool keep)
{  
  const TString outputfilename("/tmp/" + name + ".root");
  TFile* outputfile = TFile::Open(outputfilename,"RECREATE");
  TTree* subtree = mytree->CopyTree(mycut);
  data = new RooDataSet(name,title,subtree,vars);
  //  std::cerr << __LINE__ << std::endl;
  outputfile->Write();
  //  std::cerr << __LINE__ << std::endl;
  outputfile->Close();
  //  std::cerr << __LINE__ << std::endl;
  //  delete subtree;
  //  std::cerr << __LINE__ << std::endl;
  if (!keep) remove(outputfilename);
}

RooDataSet* filldataset(const TString& name, const TString& title, const RooArgList& vars, TTree* mytree, const TCut& mycut)
{  

  progressbar pbar(std::string("Filling " + name));

  TEventList myeventlist("myeventlist","myeventlist");
  mytree->Draw(">>myeventlist",mycut);

  const int nvars = vars.getSize();
  TLeaf* leafs[nvars];

  //  mytree->SetBranchStatus("*",0);
  for (int i=0; i<nvars; i++) {
    RooRealVar* var = (RooRealVar*)vars.at(i);
    leafs[i] = mytree->GetLeaf(var->GetName());
    if (!leafs[i]) {
      std::cerr << "ERROR in filldataset: no " << var->GetName() << " variable found in TTree." << std::endl;
      return NULL;
    }
    //    mytree->SetBranchStatus(var->GetName(),1);
  }

  RooDataSet* data = new RooDataSet(name,title,vars);
  
  const Long64_t nevents=myeventlist.GetN();
  for (Long64_t ie=0; ie<nevents ; ie++) {
    //    std::cout << myeventlist.GetEntry(i) << std::endl;
    mytree->GetEvent(myeventlist.GetEntry(ie));

    RooArgSet newrow; 
    for (int i=0; i<nvars; i++) {
      RooRealVar* var = (RooRealVar*)vars.at(i);
      var->setVal(leafs[i]->GetValue());
      newrow.add(*var);
    }
    data->add(newrow);

    pbar.print(100.*ie/(nevents-1));
  }
  pbar.finish();
  return data;
}


void plot(TCanvas& c, const RooRealVar& var, const RooAbsPdf& mypdf, RooAbsData* data, const TString& plotname, bool plotpull, const RooArgList* components, const TString& mycut, const TString& mycat, RooCategory* cat, bool addlhcb)
{  
  const Int_t logy = c.GetLogy();
  // int colors[10];
  // colors[0]=kBlue;
  // colors[1]=kRed;
  // colors[2]=kMagenta;
  // colors[3]=kRed;
  // colors[4]=kMagenta;
  // colors[5]=kRed;
  // colors[6]=kMagenta;
  // colors[7]=kRed;
  // colors[8]=kMagenta;
  // colors[9]=kRed;

  RooCmdArg projwdataarg=RooCmdArg::none();
  RooCmdArg cutarg=RooCmdArg::none();
  RooCmdArg slicearg=RooCmdArg::none();
  if (mycat!="") {
    cutarg=RooFit::Cut(TString(cat->GetName()) + "==" + cat->GetName() + "::" + mycat);
    //    cat->setLabel(mycat);
    projwdataarg=RooFit::ProjWData(*cat,*data);
    slicearg=RooFit::Slice(*cat,mycat);
  } else {
    projwdataarg=RooFit::ProjWData(var,*data);
  }

  RooCmdArg cutrangearg=RooCmdArg::none();
  RooCmdArg projrangearg=RooCmdArg::none();
  if (mycut!="") {
    cutrangearg = RooFit::CutRange(mycut);
    projrangearg = RooFit::ProjectionRange(mycut);
  }
  RooCmdArg colorw=RooFit::LineColor(kBlack);
  if (!components) colorw=RooFit::LineColor(kBlue);

  TText lhcbLabel(0.75, 0.88, "LHCb");
  lhcbLabel.SetTextFont(lhcbFont);
  lhcbLabel.SetTextColor(1);
  lhcbLabel.SetTextSize(lhcbTSize*1.2);
  lhcbLabel.SetTextAlign(12);
  lhcbLabel.SetNDC();

  RooPlot* varplot = var.frame();  
  //  varplot->SetDrawOption("0");

  data->plotOn(varplot,cutrangearg,cutarg);
  
  const TString varname(var.GetName());
  
  std::cout << varname << std::endl;

  //official plots
  if (addlhcb) {

    if (varname == "lambda_b0_M") {

      varplot->SetLabelSize(1.1*lhcbTSize,"x");
      varplot->SetLabelSize(1.1*lhcbTSize,"y");

      //for mass plots
      varplot->SetLabelOffset(0.04,"X");
      //  varplot->SetLabelOffset(0.04,"Y");
      
      varplot->SetTitleOffset(1.3,"X");
      varplot->SetTitleOffset(1.3,"Y");


      varplot->GetYaxis()->SetNdivisions(505, true);
    }


    if (varname.Contains("costheta")) {


      //      if (varname == "costheta")
      varplot->SetMaximum(varplot->GetMaximum()*1.1);

      varplot->SetLabelSize(1.5*lhcbTSize,"x");
      varplot->SetLabelSize(1.5*lhcbTSize,"y");

      //      if (varname.Contains("costheta1")) varplot->GetYaxis()->SetNdivisions(508, true);

      varplot->SetTitleSize(1.5*lhcbTSize,"x");
      varplot->SetTitleSize(1.5*lhcbTSize,"y");

      //for costheta plots
      //      varplot->SetLabelOffset(0.025,"X");
      //  varplot->SetLabelOffset(0.04,"Y");
      
      //      varplot->SetTitleOffset(1.3,"X");      
      varplot->SetLabelOffset(0.04,"X");
      varplot->SetTitleOffset(1.05,"Y"); 

      lhcbLabel.SetTextSize(lhcbTSize*1.5);



      //acceptance
      varplot->GetYaxis()->SetNdivisions(505, true);
      //      varplot->SetTitleOffset(1.4,"Y"); 


    }



  }  


  const double binsize= (var.getMax() - var.getMin()) / var.getBins();
  
  char buffer1[80];
  if (varname.Contains("costheta")) {
    sprintf(buffer1, "%.2f", binsize);
  } else {
    sprintf(buffer1, "%.1f", binsize);
  }
  if ( TString(var.getUnit()) == "") {
    varplot->GetYaxis()->SetTitle("Candidates / " + TString(buffer1)  );
    varplot->GetXaxis()->SetTitle(var.GetTitle());
  } else {
    varplot->GetYaxis()->SetTitle("Candidates / (" + TString(buffer1) + " " +  TString(var.getUnit()) + ")");
    varplot->GetXaxis()->SetTitle(TString(var.GetTitle()) + " [" + var.getUnit() + "]");
  }
  mypdf.plotOn(varplot,colorw,projrangearg,slicearg,projwdataarg);  
  RooHist* hpull =NULL;
  if (plotpull) hpull = varplot->pullHist();
  if (components) {
    RooArgList signalpdfs,bkgpdfs,otherpdfs;
    for (int i=0; i<components->getSize(); i++) {
      RooAbsPdf* mycomppdf = (RooAbsPdf*)components->at(i);
      TString pdfname(mycomppdf->GetName());
      pdfname.ToLower();
      std::cout << "pdfname = " << pdfname << std::endl;
      if (pdfname.Contains("sig")) {
	signalpdfs.add(*mycomppdf);
      } else if (pdfname.Contains("bkg")) {
	bkgpdfs.add(*mycomppdf);
      } else {
	otherpdfs.add(*mycomppdf);
      }
    }
    if (signalpdfs.getSize()>0) mypdf.plotOn(varplot,RooFit::Components(signalpdfs),RooFit::LineColor(kBlue),projrangearg,projwdataarg,slicearg);
    if (bkgpdfs.getSize()>0) mypdf.plotOn(varplot,RooFit::Components(bkgpdfs),RooFit::LineColor(kRed), RooFit::LineStyle(kDashed),projrangearg,projwdataarg,slicearg);
    if (otherpdfs.getSize()>0) mypdf.plotOn(varplot,RooFit::Components(otherpdfs),RooFit::LineColor(kMagenta),projrangearg,projwdataarg,slicearg);
  }

  if (plotpull) {  
    RooPlot* varplot2 = var.frame() ;
    varplot2->addPlotable(hpull,"P") ;
    varplot2->GetXaxis()->SetTitle("");
    varplot2->GetYaxis()->SetTitle("");
    //    TPad *p1 = new TPad("i1", "i1", 0.0, 0.85, 1.0, 1.0);
    //    TPad *p2 = new TPad("i2", "i2", 0.0, 0.0, 1.0, 0.85);

    TPad *p1 = new TPad("i1", "i1", 0.0, 0.80, 1.0, 0.95);
    TPad *p2 = new TPad("i2", "i2", 0.0, 0.0, 1.0, 0.80);
    p1->SetBottomMargin(0.);
    p2->SetTopMargin(0.);
    varplot2->GetYaxis()->SetNdivisions(10);

    p2->SetLogy(logy);
    p1->Draw();p1->cd();varplot2->Draw();
    c.cd();
    p2->Draw();p2->cd();varplot->Draw();   
    if (addlhcb) lhcbLabel.Draw();
    c.SaveAs(plotname + ".eps");
    c.SaveAs(plotname + ".C");
    c.SaveAs(plotname + "-plot.root");
    delete varplot;
    delete varplot2;
    c.Clear();
    c.SetLogy(logy);
  } else {
//    varplot->SetMinimum(0.0001);
    varplot->Draw("0");
    if (addlhcb) lhcbLabel.Draw();
    c.SaveAs(plotname + ".eps");
    //    c.SaveAs(plotname + ".png");
    c.SaveAs(plotname + ".C");
    c.SaveAs(plotname + "-plot.root");
    delete varplot;
  }
  
}


void plot(TCanvas& c, const RooArgList& vars, const RooAbsPdf& w, RooAbsData* data, const TString& plotname, bool plotpull, const RooArgList* components, const TString& mycut, RooCategory* cat, bool addlhcb)
{
  
  for (int i=0; i<vars.getSize(); i++) {
    RooRealVar* var = (RooRealVar*)vars.at(i);
    if (var==(RooRealVar*)cat) continue;
    std::cerr << var->GetName() << std::endl;
    if (mycut!="") {
      plot(c,*var,w,data,plotname + "-cut_" +  mycut + "-" + var->GetName() ,plotpull, components, mycut, "", NULL, addlhcb);
      if (cat) {
	TIterator* tIter = cat->typeIterator() ;
	//	RooCatType* type = 0;
	while (RooCatType* type=(RooCatType*)tIter->Next()) {
	  plot(c,*var,w,data,plotname + "-cut_" +  mycut + "-cat_" + type->GetName() + "-" + var->GetName() ,plotpull, components, mycut, type->GetName(), cat, addlhcb);
	}
      }
    } else {
      plot(c,*var,w,data,plotname + "-" + var->GetName() ,plotpull, components, mycut, "", NULL, addlhcb);
      if (cat){
	TIterator* tIter = cat->typeIterator() ;
	//	RooCatType* type = 0;
	while ( RooCatType* type=(RooCatType*)tIter->Next()) {
	  plot(c,*var,w,data,plotname + "-cat_" + type->GetName() + "-" + var->GetName() ,plotpull, components, mycut, type->GetName(), cat, addlhcb);
	}
      }
    }
  }
}

void plotpdf(const TCanvas& c, const RooRealVar& var, const RooAbsPdf& w, const TString& plotname)
{
  RooPlot* varplot = var.frame();
  w.plotOn(varplot);
  varplot->Draw();
  c.SaveAs(plotname);
}

void plotpdflist(const TCanvas& c, const RooRealVar& var, const RooArgList& pdflist, const TString& plotname)
{
  Int_t colors[5] = {kBlack,kRed,kBlue,kMagenta,kYellow};

  RooPlot* varplot = var.frame();
  for (int i=0; i<pdflist.getSize(); i++) {
    RooAbsPdf* w = (RooAbsPdf*)pdflist.at(i);
    w->plotOn(varplot, RooFit::LineColor(colors[i]));    
  }
  varplot->GetYaxis()->SetTitle("");
  varplot->Draw();
  c.SaveAs(plotname);
}

double minofdataset(RooDataSet* data, const RooRealVar& var)
{
  double minval=1000000.;
  for (int i=0;i<data->sumEntries();i++) {
    const RooArgSet* blu = data->get(i);
    const double thisval=blu->getRealValue(var.GetName());
    if (thisval<minval) minval=thisval;
  }
  return minval;
}

double maxofdataset(RooDataSet* data, const RooRealVar& var)
{
  double maxval=-1000000.;
  for (int i=0;i<data->sumEntries();i++) {
    const RooArgSet* blu = data->get(i);
    const double thisval=blu->getRealValue(var.GetName());
    if (thisval>maxval) maxval=thisval;
  }
  return maxval;
}

void plotdata(const TCanvas& c, RooRealVar& var, RooDataSet* data, const TString& plotname, TString mycut, bool fitGauss, RooRealVar* mean, RooRealVar* sigma)
{

  RooPlot* varplot = NULL;
  double lo=0.;
  double hi=0.;
  if (var.hasMin() and var.hasMax()) {
    lo=var.getMin();
    hi=var.getMax();
  } else {
    lo=minofdataset(data,var);
    hi=maxofdataset(data,var);
  }
  
  RooRealVar pullMean("Mean","Mean of pull",lo+(hi-lo)/2.,lo,hi) ;
  RooRealVar pullSigma("Sigma","Width of pull",(hi-lo)/4.,0.,(hi-lo)) ;

  RooGaussian* pullGauss = NULL;

  //  bool fitGauss=true;
  if (fitGauss) {
    // std::cout << var.GetName() << std::endl;
    // std::cout << "lo = " << lo << std::endl;
    // std::cout << "hi = " << hi << std::endl;
    // mygetchar;
    //    RooGenericPdf pullGauss("pullGauss","Gaussian of pull",
    //			    "exp(-0.5*(@0-@1)*(@0-@1)/(@2*@2))",
    //			    RooArgSet(var,pullMean,pullSigma)) ;
    if (mean and sigma) {
      mean->setVal(lo+(hi-lo)/2.);
      mean->setMin(lo);
      mean->setMax(hi);
      sigma->setVal((hi-lo)/4.);
      sigma->setMin(0.);
      sigma->setMax((hi-lo));
      pullGauss = new RooGaussian("pullGauss","Gaussian of pull",
				  var,*mean,*sigma) ;
    } else {
      pullGauss = new RooGaussian("pullGauss","Gaussian of pull",
				  var,pullMean,pullSigma) ;
    }

    pullGauss->fitTo(*data,"mh") ;

    // if (mean and sigma) {
    //   lo=mean->getVal()-5.*sigma->getVal();
    //   hi=mean->getVal()+5.*sigma->getVal();
    // } else {
    //   lo=pullMean.getVal()-5.*pullSigma.getVal();
    //   hi=pullMean.getVal()+5.*pullSigma.getVal();
    // }
    
  }
  varplot = var.frame();
    

  data->plotOn(varplot, RooFit::Cut(mycut));
  if (pullGauss)  {
    pullGauss->plotOn(varplot) ;
    pullGauss->paramOn(varplot,data) ;
  }
  varplot->GetYaxis()->SetTitle("");

  varplot->Draw();
  c.SaveAs(plotname);
  delete varplot;
  if (pullGauss) delete pullGauss;
}

void plotdata(const TCanvas& c, const RooArgList& vars, RooDataSet* data, const TString& plotname, TString mycut, bool fitGauss, RooRealVar* mean, RooRealVar* sigma)
{
  for (int i=0; i<vars.getSize(); i++) {
    RooRealVar* var = (RooRealVar*)vars.at(i);
    plotdata(c,*var,data, plotname + "-" +TString(var->GetName()) + ".eps", mycut, fitGauss,mean,sigma);
  } 
}

TH1F pulltwohistos(const TH1F& histo1norm, const TH1F& histo2norm)
{
  TH1F historatio(histo1norm);
  historatio.Divide(&histo2norm);
  historatio.GetXaxis()->SetTitle("");
  //  historatio.SetMinimum(0.);

  TH1F histopull(*(TH1F*)historatio.Clone("histopull"));// ("histopull","histopull",historatio.GetNbinsX(),historatio.GetXaxis()->GetXbins()->GetArray());
  //  TH1F histopull("histopull","histopull",historatio.GetNbinsX(),historatio.GetXaxis()->GetXmin(),historatio.GetXaxis()->GetXmax());
  for (int i=1;i<=historatio.GetNbinsX();i++) {
    if (historatio.GetBinError(i)==0) {
      histopull.SetBinContent(i,0.);
    } else {
      histopull.SetBinContent(i,(historatio.GetBinContent(i)-1.)/historatio.GetBinError(i));
    }
  }

  if (histopull.GetMaximum()>-histopull.GetMinimum()) {
    histopull.SetMinimum(-histopull.GetMaximum());
  } else {
    histopull.SetMaximum(-histopull.GetMinimum());
  } 
  histopull.GetXaxis()->SetLabelSize(0);
  return histopull;
}


Double_t findnextbin(RooDataSet* data, Double_t entriesperbin, const TString& varname, Double_t smin, Double_t smax)
{
  //  bool converged=false;
  Double_t thismin=smin;
  Double_t thismax=smax;
  const Double_t maxdiff=5;
  
  int iiter=0;
  const int maxiter=50;
  while (iiter<maxiter) {
    const Double_t thiscut=(thismax+thismin)/2.;
    const double sumEntries = data->sumEntries(varname + ">=" + str(smin) + "&&" + varname + "<=" + str(thiscut));
    
    //     std::cout << "thismin    = " << thismin << std::endl;
    //     std::cout << "thismax    = " << thismax << std::endl;
    //     std::cout << "thiscut    = " << thiscut << std::endl;
    //     std::cout << "sumEntries = " << sumEntries << std::endl;
    
    if (fabs(sumEntries-entriesperbin)<maxdiff) return thiscut;
    if (sumEntries>entriesperbin) {
      thismax=thiscut;
    } else {
      thismin=thiscut;
    }
    iiter++;
  }
  return (thismin+thismax)/2;
}

void BinLogX(TH1*h)
{

  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for (int i = 0; i <= bins; i++) {
    //    new_bins[i] = TMath::Power(10, from + i * width);
    new_bins[i] = pow(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete new_bins;
} 

void plotdatavsdata2d(TCanvas& c, RooRealVar* var1, RooRealVar* var2, RooDataSet* data1, RooDataSet* data2, const TString& plotname, bool normbinning, bool writerootfile, int nbinx, int nbiny)
{
  //  std::cout << "data1 " << data1->numEntries() << " data2 " << data2->numEntries() << std::endl;
  // RooPlot* varplot = var.frame(25);

  // data1->plotOn(varplot, RooFit::Cut(mycut), RooFit::Rescale(1./data1->numEntries()));
  // data2->plotOn(varplot, RooFit::Cut(mycut), RooFit::Rescale(1./data2->numEntries()), RooFit::LineColor(kRed));
  // varplot->Draw();

  //find optbinning

  TH2F* histo1 = NULL;
  TH2F* histo2 = NULL;
  if (nbinx==0 and nbiny==0) {
    histo1=new TH2F("histo1","histo1",var1->getBins(),var1->getMin(),var1->getMax(),var2->getBins(),var2->getMin(),var2->getMax());
    histo2=new TH2F("histo2","histo2",var1->getBins(),var1->getMin(),var1->getMax(),var2->getBins(),var2->getMin(),var2->getMax());
  } else {
    histo1=new TH2F("histo1","histo1",nbinx,var1->getMin(),var1->getMax(),nbiny,var2->getMin(),var2->getMax());
    histo2=new TH2F("histo2","histo2",nbinx,var1->getMin(),var1->getMax(),nbiny,var2->getMin(),var2->getMax());
  }
  
  //  TH2F* histo2 = new TH1F("histo2","histo2",var1getBins(),var.getMin(),var.getMax());
  
  TH2F* histo1nb = NULL;  
  TH2F* histo2nb = NULL;
  
  if (normbinning) {

    RooDataSet* smalldata = NULL;
    std::cout << data1->GetName() << "->sumEntries() = " << data1->sumEntries() << std::endl;
    std::cout << data2->GetName() << "->sumEntries() = " << data2->sumEntries() << std::endl;
    if (data1->sumEntries()<data2->sumEntries()) {
      smalldata = (RooDataSet*)data1->reduce(RooArgSet(*var1,var2));
    } else {
      smalldata = (RooDataSet*)data2->reduce(RooArgSet(*var1,var2));
    }
    const TString var1name(var1->GetName());
    const TString var2name(var2->GetName());
    const double entriesperbin_goal=200;
    const double sumEntries = smalldata->sumEntries(var1name + ">=" + str(var1->getMin()) + "&&" + var1name + "<=" + str(var1->getMax()) + "&&" + var2name + ">=" + str(var2->getMin()) + "&&" + var2name + "<=" + str(var2->getMax()) );
    if (nbinx==0 and nbiny==0) {
      nbinx=int(sqrt(sumEntries/entriesperbin_goal));
      nbiny=nbinx;
      //      if (nbin<4) nbin=4;
    }

    std::cout << "sumEntries = " << sumEntries << std::endl;
    std::cout << "nbinx      = " << nbinx << std::endl;
    std::cout << "nbiny      = " << nbiny << std::endl;
    const double entriesperbin=sumEntries/nbinx/nbiny;
    Double_t xbin[nbinx+1];
    xbin[0]=var1->getMin();
    xbin[nbinx]=var1->getMax();
    Double_t ybin[nbiny+1];
    ybin[0]=var2->getMin();
    ybin[nbiny]=var2->getMax();
    for (int i=1; i<nbinx; i++) {
      //    std::cout << "Find bin " << i << std::endl;
      xbin[i]=findnextbin(smalldata,entriesperbin*nbiny,var1name,xbin[i-1],xbin[nbinx]);
    }
    for (int i=1; i<nbiny; i++) {
      ybin[i]=findnextbin(smalldata,entriesperbin*nbinx,var2name,ybin[i-1],ybin[nbiny]);
    }
    std::cout << var1->GetName() << " : xbin = ";
    for (int i=0; i<=nbinx; i++) std::cout << " " << xbin[i];
    std::cout << std::endl;
    std::cout << var2->GetName() << " : ybin = ";
    for (int i=0; i<=nbiny; i++) std::cout << " " << ybin[i];
    std::cout << std::endl;
    //    mygetchar;
    //   std::cout << "varname = " << varname << std::endl;
    //   for (int i=0;i<nbin;i++) std::cout << "Bin " << i << " : " << xbin[i] << " - " << xbin[i+1] << " : sumEntries =" << data1->sumEntries(varname + ">=" + str(xbin[i]) + "&&" + varname + "<=" + str(xbin[i+1])) << std::endl;
    //   mygetchar;
    
    histo1nb = new TH2F("histo1nb","histo1nb",nbinx,xbin,nbiny,ybin);
    histo2nb = new TH2F("histo2nb","histo2nb",nbinx,xbin,nbiny,ybin);
    delete smalldata;
  } 


  data1->fillHistogram(histo1,RooArgList(*var1,*var2));
  data2->fillHistogram(histo2,RooArgList(*var1,*var2));
  if (normbinning) {
    data1->fillHistogram(histo1nb,RooArgList(*var1,*var2));
    data2->fillHistogram(histo2nb,RooArgList(*var1,*var2));
  }
  
  TH2F histo1norm((TH2F&)*histo1->DrawNormalized());
  TH2F histo2norm((TH2F&)*histo2->DrawNormalized());

  if (writerootfile) {
    TString output(plotname);
    //    output.ReplaceAll(".eps",".root");
    output+=".root";

    TFile* f = new TFile(output,"recreate");
    histo1norm.Write();
    histo2norm.Write();
    if (normbinning) {
      histo1nb->Write();
      histo2nb->Write();
    }
    //  histo1.Write();
    //  histo2.Write();
    f->Close();
  }
  delete histo1;
  delete histo2;
  if (normbinning) {
    delete histo1nb;
    delete histo2nb;
  }
}



void plotdatavsdata(TCanvas& c, const RooRealVar& var, RooDataSet* data1, RooDataSet* data2, const TString& plotname, bool normbinning, bool plotpull, bool writerootfile, double cutval, bool gp)
{
  std::cout << "plotdatavsdata: var=" << var.GetName() << std::endl;

  //  std::cout << "data1 " << data1->numEntries() << " data2 " << data2->numEntries() << std::endl;
  // RooPlot* varplot = var.frame(25);

  // data1->plotOn(varplot, RooFit::Cut(mycut), RooFit::Rescale(1./data1->numEntries()));
  // data2->plotOn(varplot, RooFit::Cut(mycut), RooFit::Rescale(1./data2->numEntries()), RooFit::LineColor(kRed));
  // varplot->Draw();

  //find optbinning

  TH1F* histo1 = new TH1F("histo1","histo1",var.getBins(),var.getMin(),var.getMax());
  TH1F* histo2 = new TH1F("histo2","histo2",var.getBins(),var.getMin(),var.getMax());
  
  TH1F* histo1nb = NULL;  
  TH1F* histo2nb = NULL;
  
  if (normbinning) {

    RooDataSet* smalldata = NULL;
    std::cout << data1->GetName() << "->sumEntries() = " << data1->sumEntries() << std::endl;
    std::cout << data2->GetName() << "->sumEntries() = " << data2->sumEntries() << std::endl;
    if (data1->sumEntries()<data2->sumEntries()) {
      smalldata = (RooDataSet*)data1->reduce(RooArgSet(var));
    } else {
      smalldata = (RooDataSet*)data2->reduce(RooArgSet(var));
    }
    const TString varname(var.GetName());
    const double entriesperbin_goal=200;
    const double sumEntries = smalldata->sumEntries(varname + ">=" + str(var.getMin()) + "&&" + varname + "<=" + str(var.getMax()));
    const int nbin=int(sumEntries/entriesperbin_goal);
    std::cout << "sumEntries = " << sumEntries << std::endl;
    std::cout << "nbin       = " << nbin << std::endl;
    const double entriesperbin=sumEntries/nbin;
    Double_t xbin[nbin+1];
    xbin[0]=var.getMin();
    xbin[nbin]=var.getMax();
    for (int i=1; i<nbin; i++) {
      //    std::cout << "Find bin " << i << std::endl;
      xbin[i]=findnextbin(smalldata,entriesperbin,varname,xbin[i-1],xbin[nbin]);
    }
    
    //   std::cout << "varname = " << varname << std::endl;
    //   for (int i=0;i<nbin;i++) std::cout << "Bin " << i << " : " << xbin[i] << " - " << xbin[i+1] << " : sumEntries =" << data1->sumEntries(varname + ">=" + str(xbin[i]) + "&&" + varname + "<=" + str(xbin[i+1])) << std::endl;
    //   mygetchar;
    
    histo1nb = new TH1F("histo1nb","histo1nb",nbin,xbin);
    histo2nb = new TH1F("histo2nb","histo2nb",nbin,xbin);
    delete smalldata;
  } 


  data1->fillHistogram(histo1,var);
  data2->fillHistogram(histo2,var);
  if (normbinning) {
    data1->fillHistogram(histo1nb,var);
    data2->fillHistogram(histo2nb,var);
  }

  TH1F histo1norm((TH1F&)*histo1->DrawNormalized());
  TH1F histo2norm((TH1F&)*histo2->DrawNormalized());

  TH1F histopull=pulltwohistos(histo1norm, histo2norm);
    
  
  TLine zeroline(var.getMin(),0.,var.getMax(),0.);
  zeroline.SetLineColor(kBlack);
  TLine threepline(var.getMin(),3.,var.getMax(),3.);
  TLine threemline(var.getMin(),-3.,var.getMax(),-3.);
  threepline.SetLineColor(kGreen);
  threemline.SetLineColor(kGreen);
  threepline.SetLineStyle(kDashed);
  threemline.SetLineStyle(kDashed);
  TLine fivepline(var.getMin(),5.,var.getMax(),5.);
  TLine fivemline(var.getMin(),-5.,var.getMax(),-5.);
  fivepline.SetLineColor(kRed);
  fivemline.SetLineColor(kRed);
  fivepline.SetLineStyle(kDashed);
  fivemline.SetLineStyle(kDashed);

  histo1norm.SetMinimum(0.);
  histo2norm.SetMinimum(0.);
  histo1norm.GetXaxis()->SetTitle(var.getTitle(true));
  histo2norm.GetXaxis()->SetTitle(var.getTitle(true));
  histo1norm.SetLineColor(kBlue);
  histo1norm.SetFillColor(kBlue);

  histo2norm.SetLineColor(kRed);
  histo2norm.SetFillColor(kRed);
  histo2norm.SetFillStyle(3345);

  if (histo2norm.GetMaximum()>histo1norm.GetMaximum()) histo1norm.SetMaximum(histo2norm.GetMaximum()*1.1);

  TLine cutline(cutval,0.,cutval,0.5*(histo1norm.GetMaximum()+histo2norm.GetMaximum()));
  double arrowlength=0.1;
  if (!gp) arrowlength*=-1.;
  TArrow arrow(cutval,0.5*(histo1norm.GetMaximum()+histo2norm.GetMaximum()),cutval+ arrowlength*(var.getMax()-var.getMin()) , 0.5*(histo1norm.GetMaximum()+histo2norm.GetMaximum()));
  cutline.SetLineColor(kBlack);
  arrow.SetLineColor(kBlack);

  c.Clear();
  if (plotpull) {    
    c.cd();
    TPad *p2 = new TPad("i2", "i2", 
			0.0, 0.05, 
			1.0, 0.80);
    TPad *p1 = new TPad("i1", "i1", 
			0.0, 0.80, 
			1.0, 0.95);

    p1->SetBottomMargin(0.);
    p2->SetTopMargin(0.);

    p1->Draw();
    p2->Draw();

    histopull.GetYaxis()->SetNdivisions(10);
    p1->cd();
    histopull.Draw("HIST");
    zeroline.Draw();
    threepline.Draw();
    threemline.Draw();
    fivepline.Draw();
    fivemline.Draw();
    
    p2->cd();
    histo1norm.Draw("HIST");
    histo2norm.Draw("HIST SAME");
    if (cutval!=-666.) {
      cutline.Draw();
      arrow.Draw();
    }
    c.SaveAs(plotname);
    c.Clear();
  } else {
    histo1norm.Draw("HIST");
    histo2norm.Draw("HIST SAME");
    if (cutval!=-666.) {
      cutline.Draw();
      arrow.Draw();
    }
    c.SaveAs(plotname);
  }
  


  if (writerootfile) {
    TString output(plotname);
    output.ReplaceAll(".eps",".root");
    
    TFile* f = new TFile(output,"recreate");
    histo1norm.Write();
    histo2norm.Write();
    if (normbinning) {
      histo1nb->Write();
      histo2nb->Write();
    }
    //  histo1.Write();
    //  histo2.Write();
    f->Close();
  }
  delete histo1;
  delete histo2;
  if (normbinning) {
    delete histo1nb;
    delete histo2nb;
  }
}

void plotdatavsdata(TCanvas& c, const RooArgList& vars, RooDataSet* data1, RooDataSet* data2, const TString& plotname, bool normbinning, bool plotpull, bool writerootfile)
{
  for (int i=0; i<vars.getSize(); i++) {
    RooRealVar* var = (RooRealVar*)vars.at(i);
    plotdatavsdata(c,*var,data1,data2, plotname + "-" +TString(var->GetName()) + ".eps", normbinning, plotpull, writerootfile);
  } 
}



void plotttree(TCanvas& c, const RooArgList& vars, TTree* tree, const TCut& cut1, const TCut& cut2, const TString& output, bool plotpull, TString leg1, TString leg2)
{
  for (int i=0;i<vars.getSize();i++) {
    RooRealVar* var = (RooRealVar*)vars.at(i);
    TH1F histo1("histo1","histo1",25,var->getMin(),var->getMax());
    TH1F histo2("histo2","histo2",25,var->getMin(),var->getMax());
    tree->Draw(TString(var->GetName()) + ">>histo1",cut1);
    tree->Draw(TString(var->GetName()) + ">>histo2",cut2);
    
    //    std::cout << "histo1 = " << histo1.GetEntries() << " ; histo2 = " << histo2.GetEntries() << std::endl;

    histo1.Sumw2();
    histo2.Sumw2();

    histo1.GetXaxis()->SetTitle(var->getTitle(true));
    histo1.SetMinimum(0.); 
    histo2.SetLineColor(kRed);
    
    TH1F histo1norm((TH1F&)*histo1.DrawNormalized());
    TH1F histo2norm((TH1F&)*histo2.DrawNormalized());
    histo1norm.SetMinimum(0.);
    histo2norm.SetMinimum(0.);
    if (histo2norm.GetMaximum()>histo1norm.GetMaximum()) histo1norm.SetMaximum(histo2norm.GetMaximum()*1.1);

    TH1F histopull=pulltwohistos(histo1norm, histo2norm);

    TLine redline(var->getMin(),0.,var->getMax(),0.);
    redline.SetLineColor(kRed);

    TLine zeroline(var->getMin(),0.,var->getMax(),0.);
    zeroline.SetLineColor(kBlack);
    TLine threepline(var->getMin(),3.,var->getMax(),3.);
    TLine threemline(var->getMin(),-3.,var->getMax(),-3.);
    threepline.SetLineColor(kGreen);
    threemline.SetLineColor(kGreen);
    threepline.SetLineStyle(kDashed);
    threemline.SetLineStyle(kDashed);
    TLine fivepline(var->getMin(),5.,var->getMax(),5.);
    TLine fivemline(var->getMin(),-5.,var->getMax(),-5.);
    fivepline.SetLineColor(kRed);
    fivemline.SetLineColor(kRed);
    fivepline.SetLineStyle(kDashed);
    fivemline.SetLineStyle(kDashed);

    //    histo1norm.DrawNormalized("E");
    //    histo2norm.DrawNormalized("SAME");

    c.Clear();

    const bool doleg=((leg1!="") and (leg2!=""));
    TLegend* leg = NULL;
    if (doleg) {
      leg = new TLegend(0.5,0.7,0.7,0.9);
      leg->AddEntry(&histo1norm, leg1);
      leg->AddEntry(&histo2norm, leg2);
      leg->SetFillStyle(0);
    }
    if (plotpull) {
      TPad *p1 = new TPad("i1", "i1", 0.0, 0.80, 1.0, 0.95);
      TPad *p2 = new TPad("i2", "i2", 0.0, 0.05, 1.0, 0.80);
      p1->SetBottomMargin(0.);
      p2->SetTopMargin(0.);
      histopull.GetYaxis()->SetNdivisions(10);

      p1->Draw();p1->cd();histopull.Draw("HIST");redline.Draw();
      c.cd();
      p2->Draw();p2->cd();
      histo1norm.Draw("E");
      histo2norm.Draw("SAME");
    } else {
      histo1norm.Draw("E");
      histo2norm.Draw("SAME");
      if (doleg) leg->Draw();
    }
    //  c.SaveAs(plotname);
    c.SaveAs(output + "-" +TString(var->GetName()) + ".eps");
    c.Clear();
    if (doleg) delete leg;
  }
} 

void plot_mcstudy(const TCanvas& c, RooMCStudy& mgr, const TString& prefix, const RooRealVar& var)
{
  const TString outputvalue(prefix + "-value-" + TString(var.GetName()) + ".eps");
  const TString outputerror(prefix + "-error-" + TString(var.GetName()) + ".eps");
  const TString outputpull(prefix + "-pull-" + TString(var.GetName()) + ".eps");


  RooPlot* mframe = mgr.plotParam(var);
  mframe->Draw();
  c.SaveAs(outputvalue);

  // Plot the distribution of the error
  RooPlot* meframe = mgr.plotError(var);//,0.,0.1) ;
  meframe->Draw() ;
  c.SaveAs(outputerror);
  
  //plot pull
  RooPlot* pbframe = mgr.plotPull(var,-7.0, 7.0, 50, kTRUE) ;
  pbframe->Draw() ;
  c.SaveAs(outputpull);
}

void plot_mcstudy(const TCanvas& c, RooMCStudy& mgr,const TString& prefix, const RooArgList& varlist)
{
  for (int i=0; i<varlist.getSize(); i++) {
    plot_mcstudy(c,mgr,prefix,*(RooRealVar*)varlist.at(i));
  }
}



// void plotwcomponents(TCanvas& c, const RooRealVar& var, const RooAbsPdf& w, const RooArgList& components, RooAbsData* toydata, const TString& plotname,  bool plotpull)
// { 
//   const Int_t logy = c.GetLogy();

//   int colors[10];
//   colors[0]=kBlue;
//   colors[1]=kRed;
//   colors[2]=kMagenta;
  
//   RooPlot* varplot = var.frame();
//   toydata->plotOn(varplot);
//   w.plotOn(varplot,RooFit::LineColor(kBlack));
//   RooHist* hpull = varplot->pullHist() ;
//   for (int i=0; i<components.getSize(); i++) {
//     RooAbsPdf* mypdf = (RooAbsPdf*)components.at(i);
//     w.plotOn(varplot,RooFit::Components(*mypdf),RooFit::LineColor(colors[i]));
//   }

//   if (plotpull) {
//     //    RooPlot* varplot3 = var.frame();
//     //    toydata->plotOn(varplot3);
//     //    w.plotOn(varplot3,RooFit::LineColor(kBlack));
//     RooPlot* varplot2 = var.frame() ;
//     //  varplot2->addPlotable(hpull);//,"F");
//     //  varplot2->Set
//     //,RooFit::DataError(RooAbsData::None),RooFit::XErrorSize(0),RooFit::FillColor(kGray));
//     varplot2->addPlotable(hpull,"P") ;
//     varplot2->GetXaxis()->SetTitle("");
//     varplot2->GetYaxis()->SetTitle("");
//     TPad *p1 = new TPad("i1", "i1", 0.0, 0.85, 1.0, 1.0);
//     TPad *p2 = new TPad("i2", "i2", 0.0, 0.0, 1.0, 0.85);
//     p2->SetLogy(logy);
    
//     p1->Draw();p1->cd();varplot2->Draw();
//     c.cd();
//     p2->Draw();p2->cd();varplot->Draw();
//     c.SaveAs(plotname);
//     delete varplot;
//     delete varplot2;
//     c.Clear();
//     c.SetLogy(logy);
//   } else {
//     varplot->Draw();
//     c.SaveAs(plotname);
//     delete varplot;
//   }
// }


// void plotwcomponents(TCanvas& c, const RooArgList& vars, const RooAbsPdf& w, const RooArgList& components, RooAbsData* toydata, const TString& plotname, bool plotpull)
// {
//   for (int i=0; i<vars.getSize(); i++) {
//     RooRealVar* var = (RooRealVar*)vars.at(i);
//     plotwcomponents(c,*var,w,components,toydata,plotname + "-" + var->GetName() + ".eps", plotpull);
//   }
// }


// void plotwcut(const TCanvas& c, const RooRealVar& var, const RooAbsPdf& w, RooAbsData* toydata, const TString& mycut, const TString& plotname)
// {
//   RooPlot* varplot = var.frame();
//   toydata->plotOn(varplot,RooFit::CutRange(mycut));
//   w.plotOn(varplot,RooFit::ProjectionRange(mycut));
//   varplot->Draw();
//   c.SaveAs(plotname);
// }

// void plotwcut(const TCanvas& c, const RooArgList& vars, const RooAbsPdf& w, RooAbsData* toydata, const TString& mycut, const TString& plotname)
// {
//   for (int i=0; i<vars.getSize(); i++) {
//     RooRealVar* var = (RooRealVar*)vars.at(i);
//     plotwcut(c,*var, w, toydata, mycut, plotname+ "-" +  mycut + "-" + var->GetName() + ".eps");
// 	     //    RooPlot* varplot = var->frame(25);
// 	     //    toydata->plotOn(varplot,RooFit::CutRange(mycut));
// 	     //    w.plotOn(varplot,RooFit::ProjectionRange(mycut));
// 	     //    varplot->Draw();
// 	     //    c.SaveAs(plotname+ "-" +  mycut + "-" + var->GetName() + ".eps");
//   }
// }

// void plotwcutwcomponents(const TCanvas& c, const RooRealVar& var, const RooAbsPdf& w, const RooArgList& components, RooDataSet* toydata, const TString& mycut, const TString& plotname)
// {
//   int colors[10];
//   colors[0]=kBlue;
//   colors[1]=kRed;
//   colors[2]=kMagenta;

//   RooPlot* varplot = var.frame();
//   toydata->plotOn(varplot,RooFit::CutRange(mycut));
//   w.plotOn(varplot,RooFit::ProjectionRange(mycut),RooFit::LineColor(kBlack));
//   for (int i=0; i<components.getSize(); i++) {
//     RooAbsPdf* mypdf = (RooAbsPdf*)components.at(i);
//     w.plotOn(varplot,RooFit::ProjectionRange(mycut),RooFit::Components(*mypdf),RooFit::LineColor(colors[i]));
//   }
//   varplot->Draw();
//   c.SaveAs(plotname);
// }

// void plotwcutwcomponents(const TCanvas& c, const RooArgList& vars, const RooAbsPdf& w, const RooArgList& components, RooDataSet* toydata, const TString& mycut, const TString& plotname)
// {
//   for (int i=0; i<vars.getSize(); i++) {
//     RooRealVar* var = (RooRealVar*)vars.at(i);
//     plotwcutwcomponents(c,*var, w, components, toydata, mycut, plotname+ "-" +  mycut + "-" + var->GetName() + ".eps");
//     //    RooPlot* varplot = var->frame(25);
//     //    toydata->plotOn(varplot,RooFit::CutRange(mycut));
//     //    w.plotOn(varplot,RooFit::ProjectionRange(mycut));
//     //    varplot->Draw();
//     //    c.SaveAs(plotname+ "-" +  mycut + "-" + var->GetName() + ".eps");
//   }
// }




// void plotwpull(TCanvas& c, const RooRealVar& var, const RooAbsPdf& w, RooAbsData* toydata, const TString& plotname)
// {

//   const Int_t logy = c.GetLogy();

//   RooPlot* varplot = var.frame();
//   toydata->plotOn(varplot);//, RooFit::XErrorSize(0.), RooFit::YErrorSize(0.));
//   w.plotOn(varplot);

//   RooHist* hpull = varplot->pullHist() ;

//   RooPlot* varplot2 = var.frame() ;
//   //  varplot2->addPlotable(hpull);//,"F");
//   //  varplot2->Set
//   //,RooFit::DataError(RooAbsData::None),RooFit::XErrorSize(0),RooFit::FillColor(kGray));
//   varplot2->addPlotable(hpull,"P") ;
//   varplot2->GetXaxis()->SetTitle("");
//   varplot2->GetYaxis()->SetTitle("");
//   c.Clear();  

//   TPad *p1 = new TPad("i1", "i1", 0.0, 0.85, 1.0, 1.0);
//   TPad *p2 = new TPad("i2", "i2", 0.0, 0.0, 1.0, 0.85);
//   p2->SetLogy(logy);

//   p1->Draw();p1->cd();varplot2->Draw();
//   c.cd();
//   p2->Draw();p2->cd();varplot->Draw();
//   c.SaveAs(plotname);
  
//   //  c.Divide(1,2);
//   //  c.cd(1); varplot->Draw();
//   //  c.cd(2); varplot2->Draw();
//   //  c.SaveAs(plotname);


//   c.Clear();
//   c.SetLogy(logy);
//   delete varplot;
//   delete varplot2;
// }

// void plotwpull(TCanvas& c, const RooArgList& vars, const RooAbsPdf& w, RooAbsData* toydata, const TString& plotname)
// {
//   for (int i=0; i<vars.getSize(); i++) {
//     RooRealVar* var = (RooRealVar*)vars.at(i);
//     plotwpull(c,*var,w,toydata,plotname + "-" + var->GetName() + ".eps");
//   }
// }




double sumofvarindata(RooDataSet* data, const TString& varname)
{
  double result=0.;
  for (int i=0; i<data->numEntries(); i++) result+=data->get(i)->getRealValue(varname);
  return result;
}

