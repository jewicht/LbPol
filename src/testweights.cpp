//g++ testweights.cxx -o testweights `root-config --libs --cflags` -lRooFit

#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooChebychev.h"

#include "rootlogon.C"

void plotdata(const TCanvas& c, const RooRealVar& var, RooDataSet* toydata, TString plotname)
{
  RooPlot* varplot = var.frame(25);
  //  toydata->plotOn(varplot, RooFit::DataError(RooAbsData::SumW2) );
  toydata->plotOn(varplot);
  varplot->Draw();
  c.SaveAs(plotname);
}

int main(int argc, char** argv)
{

  rootlogon();

  RooRealVar x("x","x",0.,1.);
  RooRealVar w("w","w",0.,1.);

  RooFormulaVar wform("wform","x*x",RooArgList(x));
  RooRealVar c("c","c",0.);

  RooChebychev p1("p1","p1",x,RooArgList(c));

  RooDataSet* data = p1.generate(RooArgSet(x),10000);

  data->addColumn(wform);
  data->write("tmpxw.txt");

  data->Print();

  //  RooDataSet dataw("dataw","dataw",RooArgList(x,w),"w");
  //  dataw.read("tmpxw.txt",RooArgList(x,w));

  RooDataSet* dataw2 = RooDataSet::read("tmpxw.txt",RooArgList(x,w));
  RooDataSet* dataw = new RooDataSet("dataw","dataw",dataw2,RooArgSet(x,w),0,"w");

  dataw->Print("v");
  
  TCanvas canvas("canvas","canvas",100,100);

  plotdata(canvas,x,data,"test-x.eps");
  plotdata(canvas,x,dataw,"test-xw.eps");

}

